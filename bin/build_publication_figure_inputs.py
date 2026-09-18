#!/usr/bin/env python3
"""Build a staged publication figure input root from authoritative benchmark runs."""

import argparse
import csv
import json
import shutil
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path


TABLES_TO_CONCAT = [
    "method_comparison.tsv",
    "method_per_gene.tsv",
    "per_gene_gain.tsv",
    "confidence_bin_summary.tsv",
    "confidence_calibration_summary.tsv",
    "tool_confidence_weights.tsv",
    "abstention_tradeoff.tsv",
    "summary_per_gene.tsv",
    "ambiguity_summary.tsv",
]

SOURCE_FILES_TO_COPY = [
    "truth_manifest.tsv",
    "sequencing_manifest.tsv",
    "cohort_manifest.tsv",
    "truth_long.tsv",
]


def parse_args():
    repo_root = Path(__file__).resolve().parent.parent
    parser = argparse.ArgumentParser(
        description="Stage publication figure inputs from authoritative WGS/WES/RNA benchmark roots."
    )
    parser.add_argument(
        "--wgs-root",
        type=Path,
        default=repo_root / "analysis/1000g_realdata/benchmark_wgs_wave2",
        help="Authoritative WGS benchmark root.",
    )
    parser.add_argument(
        "--wes-root",
        type=Path,
        default=repo_root / "analysis/1000g_realdata/benchmark_wes_truthbacked/run",
        help="Authoritative WES benchmark root.",
    )
    parser.add_argument(
        "--rna-root",
        type=Path,
        default=repo_root / "analysis/1000g_realdata/benchmark_rna_truthbacked/run",
        help="Authoritative RNA benchmark root.",
    )
    parser.add_argument(
        "--trimodal-root",
        type=Path,
        default=repo_root / "analysis/benchmark_trimodal_all_samples",
        help="Secondary trimodal benchmark root for publication figures 8/9.",
    )
    parser.add_argument(
        "--output-root",
        type=Path,
        default=repo_root / "analysis/publication_figure_inputs_v2",
        help="Output root for staged publication figure inputs.",
    )
    return parser.parse_args()


def read_tsv(path):
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path, rows, fieldnames):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def concat_tables(paths):
    merged = []
    fieldnames = []
    for path in paths:
        rows = read_tsv(path)
        if not rows:
            continue
        for key in rows[0].keys():
            if key not in fieldnames:
                fieldnames.append(key)
        merged.extend(rows)
    normalized = []
    for row in merged:
        normalized.append(dict((key, row.get(key, "")) for key in fieldnames))
    return normalized, fieldnames


def append_wgs_ensemble_rows(wgs_tables, merged_rows, fieldnames):
    extra_files = [
        "champion_challenger_method_comparison.tsv",
        "gated_method_comparison.tsv",
    ]
    rows = list(merged_rows)
    for fname in extra_files:
        path = wgs_tables / fname
        if not path.exists():
            continue
        extra_rows = read_tsv(path)
        if not extra_rows:
            continue
        for key in extra_rows[0].keys():
            if key not in fieldnames:
                fieldnames.append(key)
        rows.extend(extra_rows)
    normalized = []
    for row in rows:
        normalized.append(dict((key, row.get(key, "")) for key in fieldnames))
    return normalized


def append_wgs_ensemble_gene_rows(wgs_tables, merged_rows, fieldnames):
    extra_files = [
        "champion_challenger_method_per_gene.tsv",
        "gated_method_per_gene.tsv",
    ]
    rows = list(merged_rows)
    for fname in extra_files:
        path = wgs_tables / fname
        if not path.exists():
            continue
        extra_rows = read_tsv(path)
        if not extra_rows:
            continue
        for key in extra_rows[0].keys():
            if key not in fieldnames:
                fieldnames.append(key)
        rows.extend(extra_rows)
    normalized = []
    for row in rows:
        normalized.append(dict((key, row.get(key, "")) for key in fieldnames))
    return normalized


def merge_discordance(paths):
    counts = defaultdict(int)
    for path in paths:
        for row in read_tsv(path):
            counts[(row.get("scope", ""), row.get("tag", ""))] += int(row.get("n_events", "0") or 0)
    rows = []
    for key in sorted(counts):
        scope, tag = key
        rows.append({"scope": scope, "tag": tag, "n_events": str(counts[key])})
    return rows, ["scope", "tag", "n_events"]


def build_benchmark_metadata(wgs_root, wes_root, rna_root):
    wgs_meta = json.loads((wgs_root / "tables" / "benchmark_metadata.json").read_text(encoding="utf-8"))
    return {
        "publication_input_mode": "multi_root_staged_v2",
        "truth_source": wgs_meta.get("truth_source"),
        "truth_acquisition_date": wgs_meta.get("truth_acquisition_date"),
        "imgt_hla_version": wgs_meta.get("imgt_hla_version"),
        "supported_loci": wgs_meta.get("supported_loci", ["A", "B", "C"]),
        "tuned_consensus_min_support": "modality_specific",
        "wgs_authoritative_root": str(wgs_root),
        "wes_authoritative_root": str(wes_root),
        "rna_authoritative_root": str(rna_root),
        "wgs_split_counts": wgs_meta.get("split_counts", {}),
        "wgs_population_counts": wgs_meta.get("population_counts", {}),
        "notes": [
            "Primary WGS source is split-aware wave2 holdout benchmark.",
            "Secondary WES and RNA sources are truth-backed full-cohort modality benchmarks.",
            "Trimodal and bimodal publication figures remain sourced separately.",
        ],
    }


def copy_source_files(wgs_root, output_root):
    copied_any = False
    for fname in SOURCE_FILES_TO_COPY:
        source = wgs_root / fname
        if source.exists():
            shutil.copy2(str(source), str(output_root / fname))
            copied_any = True
    if copied_any:
        return
    alt_root = wgs_root.parent / "wgs_wave2_inputs"
    for fname in SOURCE_FILES_TO_COPY:
        source = alt_root / fname
        if source.exists():
            shutil.copy2(str(source), str(output_root / fname))


def write_caption_note(output_root):
    note = """# Publication Figure Inputs V2

This staging root merges the current authoritative modality-specific benchmarks for manuscript-facing figure generation.

- WGS authoritative source: `analysis/1000g_realdata/benchmark_wgs_wave2`
- WES authoritative source: `analysis/1000g_realdata/benchmark_wes_truthbacked/run`
- RNA authoritative source: `analysis/1000g_realdata/benchmark_rna_truthbacked/run`
- Trimodal secondary source: `analysis/benchmark_trimodal_all_samples`
"""
    (output_root / "captions" / "README.md").write_text(note, encoding="utf-8")


def main():
    args = parse_args()
    output_root = args.output_root.resolve()
    tables_dir = output_root / "tables"
    figures_dir = output_root / "figures"
    metadata_dir = output_root / "metadata"
    captions_dir = output_root / "captions"
    for directory in (tables_dir, figures_dir, metadata_dir, captions_dir):
        directory.mkdir(parents=True, exist_ok=True)

    roots = {
        "wgs": args.wgs_root.resolve(),
        "wes": args.wes_root.resolve(),
        "rna": args.rna_root.resolve(),
    }
    table_roots = {}
    for key, root in roots.items():
        table_roots[key] = root / "tables"

    for fname in TABLES_TO_CONCAT:
        paths = [table_roots["wgs"] / fname, table_roots["wes"] / fname, table_roots["rna"] / fname]
        merged_rows, fieldnames = concat_tables(paths)
        if fname == "method_comparison.tsv":
            merged_rows = append_wgs_ensemble_rows(table_roots["wgs"], merged_rows, fieldnames)
        elif fname == "method_per_gene.tsv":
            merged_rows = append_wgs_ensemble_gene_rows(table_roots["wgs"], merged_rows, fieldnames)
        write_tsv(tables_dir / fname, merged_rows, fieldnames)

    discord_rows, discord_fields = merge_discordance([
        table_roots["wgs"] / "discordance_summary.tsv",
        table_roots["wes"] / "discordance_summary.tsv",
        table_roots["rna"] / "discordance_summary.tsv",
    ])
    write_tsv(tables_dir / "discordance_summary.tsv", discord_rows, discord_fields)

    benchmark_metadata = build_benchmark_metadata(roots["wgs"], roots["wes"], roots["rna"])
    (tables_dir / "benchmark_metadata.json").write_text(
        json.dumps(benchmark_metadata, indent=2) + "\n",
        encoding="utf-8",
    )

    source_manifest = {
        "generated_from": datetime.now(timezone.utc).isoformat(),
        "wgs_root": str(args.wgs_root.resolve()),
        "wes_root": str(args.wes_root.resolve()),
        "rna_root": str(args.rna_root.resolve()),
        "trimodal_root": str(args.trimodal_root.resolve()),
        "wgs_authoritative_root": "analysis/1000g_realdata/benchmark_wgs_wave2",
        "wes_authoritative_root": "analysis/1000g_realdata/benchmark_wes_truthbacked/run",
        "rna_authoritative_root": "analysis/1000g_realdata/benchmark_rna_truthbacked/run",
        "notes": [
            "Main publication figures 1-7 should be generated from the staged tables root.",
            "Figures 8-9 remain sourced from the trimodal benchmark root.",
            "WGS wave2 is the authoritative manuscript-facing WGS benchmark.",
        ],
    }
    (metadata_dir / "source_manifest.json").write_text(
        json.dumps(source_manifest, indent=2) + "\n",
        encoding="utf-8",
    )

    copy_source_files(roots["wgs"], output_root)
    write_caption_note(output_root)

    print("Staged publication figure inputs written to %s" % output_root)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
