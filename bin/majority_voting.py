#!/usr/bin/env python3
"""Weighted consensus caller for harmonized HLA tool outputs."""

import argparse
import csv
import json
import sys
from collections import defaultdict
from pathlib import Path

MISSING_TOKENS = {"", "na", "n/a", "none", "null", "nan", "-", ".", "failed", "fail", "no_call", "no result"}


def parse_args():
    parser = argparse.ArgumentParser(description="Compute weighted consensus calls from harmonized tool outputs.")
    parser.add_argument("--calls", required=True, help="TSV with per-tool HLA calls.")
    parser.add_argument("--weights", required=True, help="JSON weights produced by hla_benchmark.py.")
    parser.add_argument("--output", required=True, help="Output TSV path.")
    parser.add_argument("--min-weight", type=float, default=0.0, help="Minimum winning total weight required to emit a call.")
    parser.add_argument("--use-gene-specific", action="store_true", help="Prefer gene-specific weights when available.")
    return parser.parse_args()


def clean_token(value):
    if value is None:
        return ""
    return str(value).strip()


def is_missing(value):
    return clean_token(value).lower() in MISSING_TOKENS


def normalize_pair(allele1, allele2):
    alleles = [clean_token(allele1), clean_token(allele2)]
    alleles = [allele for allele in alleles if not is_missing(allele)]
    if len(alleles) < 2:
        return ("", "")
    alleles = sorted(alleles)
    return (alleles[0], alleles[1])


def read_calls(path):
    with Path(path).open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return list(reader)


def load_weights(path):
    with Path(path).open("r", encoding="utf-8") as handle:
        return json.load(handle)


def is_callable(row):
    if "is_callable" in row:
        return clean_token(row.get("is_callable")) == "1"
    pair = normalize_pair(row.get("allele1"), row.get("allele2"))
    return all(pair)


def lookup_weight(weights, row, use_gene_specific):
    tool = clean_token(row.get("tool"))
    modality = clean_token(row.get("modality")).lower()
    gene = clean_token(row.get("gene"))
    if use_gene_specific:
        gene_weight = weights.get("gene_weights", {}).get(tool, {}).get(modality, {}).get(gene, {})
        value = gene_weight.get("final_weight")
        if value is not None:
            return float(value)
    tool_weight = weights.get("tool_weights", {}).get(tool, {}).get(modality, {})
    value = tool_weight.get("final_weight")
    if value is not None:
        return float(value)
    return 0.0


def build_consensus_rows(call_rows, weights, min_weight, use_gene_specific):
    grouped = defaultdict(list)
    for row in call_rows:
        if not is_callable(row):
            continue
        pair = normalize_pair(row.get("allele1"), row.get("allele2"))
        if not all(pair):
            continue
        grouped[(clean_token(row.get("sample")), clean_token(row.get("gene")))].append(row)

    output_rows = []
    for key, rows in sorted(grouped.items()):
        pair_scores = defaultdict(lambda: {"total_weight": 0.0, "agreeing_tools": 0, "contributors": []})
        total_contributing_tools = 0
        for row in rows:
            weight = lookup_weight(weights, row, use_gene_specific)
            if weight <= 0:
                continue
            total_contributing_tools += 1
            pair = normalize_pair(row.get("allele1"), row.get("allele2"))
            entry = pair_scores[pair]
            entry["total_weight"] += weight
            entry["agreeing_tools"] += 1
            entry["contributors"].append(clean_token(row.get("tool")))
        if not pair_scores:
            output_rows.append({
                "sample": key[0],
                "gene": key[1],
                "allele1": "",
                "allele2": "",
                "total_weight": "0.0000",
                "agreeing_tools": "0",
                "contributing_tools": "0",
                "consensus_status": "no_call",
                "chosen_by": "no_nonzero_weight",
            })
            continue
        winning_pair, winning_meta = max(
            pair_scores.items(),
            key=lambda item: (item[1]["total_weight"], item[1]["agreeing_tools"], "%s|%s" % item[0]),
        )
        if winning_meta["total_weight"] <= min_weight:
            output_rows.append({
                "sample": key[0],
                "gene": key[1],
                "allele1": "",
                "allele2": "",
                "total_weight": "%.4f" % winning_meta["total_weight"],
                "agreeing_tools": str(winning_meta["agreeing_tools"]),
                "contributing_tools": str(total_contributing_tools),
                "consensus_status": "no_call",
                "chosen_by": "below_min_weight",
            })
            continue
        output_rows.append({
            "sample": key[0],
            "gene": key[1],
            "allele1": winning_pair[0],
            "allele2": winning_pair[1],
            "total_weight": "%.4f" % winning_meta["total_weight"],
            "agreeing_tools": str(winning_meta["agreeing_tools"]),
            "contributing_tools": str(total_contributing_tools),
            "consensus_status": "called",
            "chosen_by": ",".join(sorted(winning_meta["contributors"])),
        })
    return output_rows


def write_tsv(path, rows, fieldnames):
    output_path = Path(path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow(dict((field, row.get(field, "")) for field in fieldnames))


def main():
    args = parse_args()
    call_rows = read_calls(args.calls)
    weights = load_weights(args.weights)
    output_rows = build_consensus_rows(call_rows, weights, args.min_weight, args.use_gene_specific)
    write_tsv(args.output, output_rows, ["sample", "gene", "allele1", "allele2", "total_weight", "agreeing_tools", "contributing_tools", "consensus_status", "chosen_by"])
    return 0


if __name__ == "__main__":
    sys.exit(main())
