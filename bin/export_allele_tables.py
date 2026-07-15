#!/usr/bin/env python3
"""
Export per-sample allele comparison CSV tables from harmonized benchmark rows.

Each output row is one (sample, modality, gene) combination.
Columns: sample metadata + ground-truth alleles + one group of columns per
single tool (alphabetical) + WeightedConsensus columns + MajorityVote columns.

Usage:
    python3 bin/export_allele_tables.py --benchmark-dir analysis/benchmark_wgs_44samples
    python3 bin/export_allele_tables.py --all
"""

import argparse
import os
import sys
import pandas as pd

BASE = "/scratch/project_2008084/pihla-publish"

BENCHMARKS = [
    "analysis/1000g_realdata",
    "analysis/1000g_realdata/benchmark_wgs_wave1",
    "analysis/benchmark_wgs_44samples",
    "analysis/benchmark_wes_50samples",
    "analysis/benchmark_rna_50samples",
    "analysis/benchmark_trimodal_50samples",
    "analysis/full_cohort_benchmark",
]


def find_cohort_manifest(benchmark_dir):
    for candidate in [
        os.path.join(benchmark_dir, "manifests", "cohort_manifest.tsv"),
        os.path.join(benchmark_dir, "tables", "cohort_manifest.tsv"),
    ]:
        if os.path.exists(candidate):
            return candidate
    return None


def merge_consensus_methods(wide, benchmark_dir):
    """Merge WeightedConsensus and MajorityVote columns into the wide table."""
    tables_dir = os.path.join(benchmark_dir, "tables")

    specs = [
        (
            "weighted_consensus_calls.tsv",
            "WeightedConsensus",
            ["allele1", "allele2", "is_callable", "is_correct",
             "support_fraction", "total_weight", "discordance_tag"],
        ),
        (
            "majority_vote_baseline.tsv",
            "MajorityVote",
            ["allele1", "allele2", "is_callable", "is_correct",
             "support_fraction", "discordance_tag"],
        ),
    ]

    for filename, prefix, cols in specs:
        path = os.path.join(tables_dir, filename)
        if not os.path.exists(path):
            continue
        ens = pd.read_csv(path, sep="\t", dtype=str)
        available = [c for c in cols if c in ens.columns]
        ens = ens[["sample", "modality", "gene"] + available].copy()
        rename = {c: f"{prefix}_{c}" for c in available}
        rename["is_callable"] = f"{prefix}_callable"
        rename["is_correct"] = f"{prefix}_correct"
        ens = ens.rename(columns=rename)
        wide = wide.merge(ens, on=["sample", "modality", "gene"], how="left")

    return wide


def export_benchmark(benchmark_dir):
    rows_path = os.path.join(benchmark_dir, "tables", "harmonized_benchmark_rows.tsv")
    if not os.path.exists(rows_path):
        print(f"  SKIP (no harmonized_benchmark_rows.tsv): {benchmark_dir}")
        return

    df = pd.read_csv(rows_path, sep="\t", dtype=str)

    # Keep only the columns we need
    keep = [
        "sample", "modality", "tool", "gene",
        "truth_allele1", "truth_allele2",
        "allele1", "allele2",
        "is_callable", "is_correct",
        "confidence_score",
    ]
    missing = [c for c in keep if c not in df.columns]
    if missing:
        print(f"  SKIP (missing columns {missing}): {benchmark_dir}")
        return

    df = df[keep].copy()

    # Merge cohort metadata (population, split)
    manifest_path = find_cohort_manifest(benchmark_dir)
    meta_cols = ["sample"]
    if manifest_path:
        meta = pd.read_csv(manifest_path, sep="\t", dtype=str)
        for col in ["population", "split"]:
            if col in meta.columns:
                meta_cols.append(col)
        df = df.merge(meta[meta_cols].drop_duplicates("sample"), on="sample", how="left")
    else:
        df["population"] = ""
        df["split"] = ""

    # Truth alleles are the same for every tool on a given (sample, modality, gene).
    # Extract them once, then build the wide pivot.
    truth = (
        df[["sample", "modality", "gene", "truth_allele1", "truth_allele2"]]
        .drop_duplicates(["sample", "modality", "gene"])
    )

    # Build per-tool columns via pivot
    tool_cols = {}
    for col, suffix in [
        ("allele1",         "allele1"),
        ("allele2",         "allele2"),
        ("is_callable",     "callable"),
        ("is_correct",      "correct"),
        ("confidence_score","confidence"),
    ]:
        pivoted = df.pivot_table(
            index=["sample", "modality", "gene"],
            columns="tool",
            values=col,
            aggfunc="first",
        )
        pivoted.columns = [f"{t}_{suffix}" for t in pivoted.columns]
        tool_cols[suffix] = pivoted

    # Combine all tool columns in alphabetical tool order
    tools_sorted = sorted(df["tool"].dropna().unique())
    wide = pd.concat(tool_cols.values(), axis=1).reset_index()

    # Reorder columns: for each tool, put its 5 columns together
    tool_column_order = []
    for t in tools_sorted:
        for suffix in ["allele1", "allele2", "callable", "correct", "confidence"]:
            col = f"{t}_{suffix}"
            if col in wide.columns:
                tool_column_order.append(col)

    # Merge truth alleles into wide table
    wide = wide.merge(truth, on=["sample", "modality", "gene"], how="left")

    # Merge metadata
    if "population" in df.columns and "split" in df.columns:
        meta_wide = (
            df[["sample", "modality", "gene", "population", "split"]]
            .drop_duplicates(["sample", "modality", "gene"])
        )
        wide = wide.merge(meta_wide, on=["sample", "modality", "gene"], how="left")
    else:
        wide["population"] = ""
        wide["split"] = ""

    # Final column order
    front_cols = ["sample", "population", "split", "modality", "gene",
                  "truth_allele1", "truth_allele2"]
    final_cols = front_cols + [c for c in tool_column_order if c in wide.columns]
    # Include any remaining tool columns not yet in order
    for c in wide.columns:
        if c not in final_cols:
            final_cols.append(c)
    wide = wide[final_cols]

    # Merge ensemble method columns
    wide = merge_consensus_methods(wide, benchmark_dir)

    # Sort rows
    wide = wide.sort_values(["sample", "modality", "gene"]).reset_index(drop=True)

    out_path = os.path.join(benchmark_dir, "tables", "allele_comparison_table.csv")
    wide.to_csv(out_path, index=False)

    n_samples = wide["sample"].nunique()
    n_rows = len(wide)
    n_tools = len(tools_sorted)
    ensemble_cols = [c for c in wide.columns if c.startswith("WeightedConsensus_") or c.startswith("MajorityVote_")]
    print(f"  OK  {out_path}")
    print(f"      {n_rows} rows  |  {n_samples} samples  |  {n_tools} single tools + ensemble ({len(ensemble_cols)} cols)")


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--benchmark-dir", help="Path to a single benchmark directory")
    parser.add_argument("--all", action="store_true",
                        help="Run for all standard benchmark directories")
    args = parser.parse_args()

    if not args.benchmark_dir and not args.all:
        parser.print_help()
        sys.exit(1)

    dirs = []
    if args.all:
        dirs = [os.path.join(BASE, b) for b in BENCHMARKS]
    else:
        d = args.benchmark_dir
        if not os.path.isabs(d):
            d = os.path.join(BASE, d)
        dirs = [d]

    for d in dirs:
        print(f"\n[{os.path.relpath(d, BASE)}]")
        export_benchmark(d)

    print("\nDone.")


if __name__ == "__main__":
    main()
