#!/usr/bin/env python3
"""Build phase-gated real-data inputs from existing 1000G run outputs."""

import argparse
import csv
from collections import defaultdict
from pathlib import Path

import hla_benchmark as hb

REQUIRED_MODALITIES = ["wgs", "wes", "rnaseq"]
DEFAULT_SAMPLES = ["NA06985", "NA06986", "NA06994"]
DEFAULT_GENES = ["A", "B", "C", "DRB1", "DQB1"]


def parse_args():
    parser = argparse.ArgumentParser(description="Build phase-gated 1000G truth/sequencing inputs from finished run outputs.")
    parser.add_argument("--truth-csv", required=True, help="Ground-truth CSV/TSV file in Gourraud-like wide format.")
    parser.add_argument("--output-dir", required=True, help="Output directory for generated truth and sequencing tables.")
    parser.add_argument("--wgs-results", default="", help="WGS result root (contains sample folders and/or by_tool).")
    parser.add_argument("--wes-results", default="", help="WES result root (contains sample folders and/or by_tool).")
    parser.add_argument("--rnaseq-results", default="", help="RNA-seq result root (contains sample folders and/or by_tool).")
    parser.add_argument("--samples", default=",".join(DEFAULT_SAMPLES), help="Comma-separated samples to include in this phase-gated release.")
    parser.add_argument("--supported-loci", default=",".join(DEFAULT_GENES), help="Comma-separated truth loci to export.")
    return parser.parse_args()


def parse_sample_list(value):
    samples = [hb.clean_token(token) for token in value.split(",")] if value else []
    return [sample for sample in samples if sample]


def detect_sample_tool_outputs(results_root, tool_suffixes):
    root = Path(results_root)
    sample_tools = defaultdict(set)
    if not results_root or not root.exists():
        return sample_tools

    for sample_dir in sorted([d for d in root.iterdir() if d.is_dir() and d.name not in {"by_tool", "pipeline_info"}]):
        for tool_dir in [d for d in sample_dir.iterdir() if d.is_dir()]:
            has_files = any(path.is_file() for path in tool_dir.rglob("*"))
            if has_files:
                sample_id = sample_dir.name.split('.')[0]
                sample_tools[sample_id].add(tool_dir.name.lower())

    by_tool = root / "by_tool"
    if by_tool.exists():
        for tool_dir in sorted([d for d in by_tool.iterdir() if d.is_dir()]):
            suffix = tool_suffixes.get(tool_dir.name.lower(), "_" + tool_dir.name.lower())
            for file_path in tool_dir.iterdir():
                if not file_path.is_file():
                    continue
                name = file_path.name
                if name.endswith(".txt") and suffix in name:
                    sample = hb.clean_token(name.split(suffix, 1)[0])
                    if sample:
                        sample_tools[sample].add(tool_dir.name.lower())
    return sample_tools


def modality_expected_tools(modality):
    if modality == "wgs":
        return ["spechla", "hlahd", "arcashla", "optitype", "polysolver", "kourami", "t1k"]
    if modality == "wes":
        return ["spechla", "hlahd", "optitype", "polysolver", "kourami", "t1k", "arcashla"]
    return ["arcashla", "optitype", "seq2hla", "t1k", "spechla", "hlahd"]


def write_truth_long(truth_path, output_path, target_samples, target_genes):
    rows = hb.read_table(Path(truth_path))
    if rows and len(rows[0].keys()) == 1 and "\t" in next(iter(rows[0].keys())):
        rows = hb.read_table(Path(truth_path), delimiter="\t")
    if not rows:
        raise ValueError("No rows found in truth file: %s" % truth_path)
    key_map = {hb.normalize_column_name(name): name for name in rows[0].keys()}
    sample_key = hb.first_existing(key_map, ["sampleid", "sample", "id"])
    if not sample_key:
        raise ValueError("Could not detect sample column in truth table: %s" % truth_path)
    sample_col = key_map[sample_key]

    output_rows = []
    for row in rows:
        sample = hb.clean_token(row.get(sample_col, ""))
        if sample not in target_samples:
            continue
        for gene in target_genes:
            key1 = next((col for col in row.keys() if hb.normalize_column_name(col) == hb.normalize_column_name("HLA-%s 1" % gene)), "")
            key2 = next((col for col in row.keys() if hb.normalize_column_name(col) == hb.normalize_column_name("HLA-%s 2" % gene)), "")
            a1 = hb.clean_token(row.get(key1, "")).split("/")[0] if key1 else ""
            a2 = hb.clean_token(row.get(key2, "")).split("/")[0] if key2 else ""
            if not a1 and not a2:
                continue
            output_rows.append({"sample": sample, "gene": gene, "allele1": a1, "allele2": a2})

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["sample", "gene", "allele1", "allele2"], delimiter="\t")
        writer.writeheader()
        for row in sorted(output_rows, key=lambda r: (r["sample"], hb.gene_sort_key(r["gene"]))):
            writer.writerow(row)


def main():
    args = parse_args()
    outdir = Path(args.output_dir)
    outdir.mkdir(parents=True, exist_ok=True)
    target_samples = parse_sample_list(args.samples)
    target_genes = hb.parse_gene_list(args.supported_loci)
    truth_long_path = outdir / "truth_long.tsv"
    write_truth_long(args.truth_csv, truth_long_path, target_samples, target_genes)

    tool_suffixes = {
        "spechla": "_spechla",
        "hlahd": "_hlahd",
        "arcashla": "_arcashla",
        "optitype": "_optitype",
        "polysolver": "_polysolver",
        "kourami": "_kourami",
        "t1k": "_t1k",
        "seq2hla": "_seq2hla",
    }
    per_modality_detected = {
        "wgs": detect_sample_tool_outputs(args.wgs_results, tool_suffixes),
        "wes": detect_sample_tool_outputs(args.wes_results, tool_suffixes),
        "rnaseq": detect_sample_tool_outputs(args.rnaseq_results, tool_suffixes),
    }

    sequencing_rows = []
    tool_rows = []
    roots = {"wgs": args.wgs_results, "wes": args.wes_results, "rnaseq": args.rnaseq_results}
    for sample in target_samples:
        for modality in REQUIRED_MODALITIES:
            detected_tools = sorted(per_modality_detected[modality].get(sample, set()))
            available = "1" if detected_tools else "0"
            sequencing_rows.append(
                {
                    "sample": sample,
                    "population": "",
                    "modality": modality,
                    "data_locator": roots.get(modality, ""),
                    "available": available,
                }
            )
            for tool in modality_expected_tools(modality):
                status = "available" if tool in detected_tools else "not_available"
                tool_rows.append({"sample": sample, "modality": modality, "tool": tool, "status": status})

    with (outdir / "sequencing_source.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["sample", "population", "modality", "data_locator", "available"], delimiter="\t")
        writer.writeheader()
        for row in sequencing_rows:
            writer.writerow(row)

    with (outdir / "tool_availability_seed.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["sample", "modality", "tool", "status"], delimiter="\t")
        writer.writeheader()
        for row in tool_rows:
            writer.writerow(row)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
