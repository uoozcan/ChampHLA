#!/usr/bin/env python3
"""Weighted consensus caller for harmonized HLA tool outputs."""

import argparse
import csv
import json
import sys
from collections import defaultdict
from pathlib import Path

MISSING_TOKENS = {"", "na", "n/a", "none", "null", "nan", "-", ".", "failed", "fail", "no_call", "no result"}
TRUTH_SUPPORTED_GENES = ("A", "B", "C", "DRB1", "DQB1")


def parse_args():
    parser = argparse.ArgumentParser(description="Compute weighted consensus calls from harmonized tool outputs.")
    parser.add_argument("--calls", required=True, help="TSV with per-tool HLA calls.")
    parser.add_argument("--weights", required=True, help="JSON weights file (ignored when --weighting equal).")
    parser.add_argument("--output", required=True, help="Output TSV path.")
    parser.add_argument("--weighting", default="calibrated", choices=["equal", "calibrated"],
                        help="Weighting strategy: 'equal' (all tools weight 1.0) or 'calibrated' (use JSON weights file). Default: calibrated.")
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

def clip_01(value):
    try:
        parsed = float(value)
    except (TypeError, ValueError):
        return 0.0
    if parsed < 0.0:
        return 0.0
    if parsed > 1.0:
        return 1.0
    return parsed


def mean_tool_weight_from_raw_accuracy(raw_accuracy, tool):
    gene_scores = raw_accuracy.get(tool, {}) if isinstance(raw_accuracy, dict) else {}
    values = [clip_01(gene_scores.get(gene)) for gene in TRUTH_SUPPORTED_GENES]
    if not values:
        return 0.0
    return sum(values) / float(len(values))



def is_callable(row):
    if "is_callable" in row:
        return clean_token(row.get("is_callable")) == "1"
    pair = normalize_pair(row.get("allele1"), row.get("allele2"))
    return all(pair)


def lookup_weight(weights, row, use_gene_specific, equal_mode=False):
    if equal_mode:
        return 1.0
    tool = clean_token(row.get("tool"))
    modality = clean_token(row.get("modality")).lower()
    gene = clean_token(row.get("gene"))

    # Native runtime schema support.
    if "tool_weights" in weights or "gene_weights" in weights:
        if use_gene_specific:
            gene_weight = weights.get("gene_weights", {}).get(tool, {}).get(modality, {}).get(gene, {})
            value = gene_weight.get("final_weight")
            if value is not None:
                return clip_01(value)
        tool_weight = weights.get("tool_weights", {}).get(tool, {}).get(modality, {})
        value = tool_weight.get("final_weight")
        if value is not None:
            return clip_01(value)
        return 0.0

    # Backward compatibility for imported calibration JSON files.
    # Schema: {"raw_accuracy": {"tool": {"gene": score}}}
    raw_accuracy = weights.get("raw_accuracy", {})
    if isinstance(raw_accuracy, dict):
        if use_gene_specific:
            value = raw_accuracy.get(tool, {}).get(gene)
            if value is not None:
                return clip_01(value)
            return 0.0
        return mean_tool_weight_from_raw_accuracy(raw_accuracy, tool)

    return 0.0


def build_consensus_rows(call_rows, weights, min_weight, use_gene_specific, equal_mode=False):
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
            weight = lookup_weight(weights, row, use_gene_specific, equal_mode=equal_mode)
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


CLINICAL_WARNING_DESCRIPTIONS = {
    "NO_CALL":               "No consensus reached for this locus — report as failed/unresolved",
    "SINGLE_TOOL":           "Only one tool contributed — result unconfirmed, consider retyping",
    "TOOL_DISCORD":          "Fewer than half of tools agree — active disagreement, manual review required",
    "PUTATIVE_HOMOZYGOSITY": "Both alleles identical — verify true homozygosity vs. allele dropout",
    "LOW_CONFIDENCE":        "Total weight < 0.30 — low-confidence call, treat with caution",
    "LOW_RESOLUTION":        "Allele at lower than 2-field resolution — insufficient for transplant matching",
    "ALLELE_AMBIGUITY":      "Ambiguous allele notation (contains '/' or 'g' suffix) — multiple alleles possible",
}

LOW_CONFIDENCE_THRESHOLD = 0.30


def generate_warnings(row):
    """Return a list of warning codes for a consensus row."""
    warnings = []
    status       = clean_token(row.get("consensus_status"))
    agreeing     = int(row.get("agreeing_tools", 0) or 0)
    contributing = int(row.get("contributing_tools", 0) or 0)
    total_weight = float(row.get("total_weight", 0) or 0)
    allele1      = clean_token(row.get("allele1"))
    allele2      = clean_token(row.get("allele2"))

    if status == "no_call":
        warnings.append("NO_CALL")
        return warnings  # downstream checks not meaningful for no-calls

    if contributing == 1:
        warnings.append("SINGLE_TOOL")
    elif contributing > 1 and agreeing < contributing / 2:
        warnings.append("TOOL_DISCORD")

    if allele1 and allele2 and allele1 == allele2:
        warnings.append("PUTATIVE_HOMOZYGOSITY")

    if total_weight < LOW_CONFIDENCE_THRESHOLD:
        warnings.append("LOW_CONFIDENCE")

    for allele in (allele1, allele2):
        if not allele:
            continue
        if "*" in allele and ":" not in allele:
            warnings.append("LOW_RESOLUTION")
            break
        if "/" in allele or allele.lower().endswith("g") or ":g" in allele.lower():
            warnings.append("ALLELE_AMBIGUITY")
            break

    return warnings


def write_tsv(path, rows, fieldnames):
    output_path = Path(path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow(dict((field, row.get(field, "")) for field in fieldnames))


def write_clinical_warnings(path, rows):
    """Write a human-readable warnings report — only rows that have at least one flag."""
    flagged = []
    for row in rows:
        codes = generate_warnings(row)
        if codes:
            for code in codes:
                flagged.append({
                    "sample":      row.get("sample", ""),
                    "gene":        row.get("gene", ""),
                    "allele1":     row.get("allele1", ""),
                    "allele2":     row.get("allele2", ""),
                    "warning":     code,
                    "description": CLINICAL_WARNING_DESCRIPTIONS.get(code, ""),
                    "details":     "tools=%s/%s weight=%.4f status=%s" % (
                        row.get("agreeing_tools", "?"),
                        row.get("contributing_tools", "?"),
                        float(row.get("total_weight", 0) or 0),
                        row.get("consensus_status", "?"),
                    ),
                })
    write_tsv(path, flagged,
              ["sample", "gene", "allele1", "allele2", "warning", "description", "details"])
    return len(flagged)


def main():
    args = parse_args()
    call_rows = read_calls(args.calls)
    weights = load_weights(args.weights)
    equal_mode = (args.weighting == "equal")
    output_rows = build_consensus_rows(call_rows, weights, args.min_weight, args.use_gene_specific, equal_mode=equal_mode)

    # Attach warnings column to consensus output
    for row in output_rows:
        codes = generate_warnings(row)
        row["warnings"] = ",".join(codes)

    write_tsv(args.output, output_rows,
              ["sample", "gene", "allele1", "allele2", "total_weight",
               "agreeing_tools", "contributing_tools", "consensus_status", "chosen_by", "warnings"])

    # Write separate clinical warnings report
    warnings_path = str(args.output).replace(".tsv", "_clinical_warnings.tsv")
    n_flags = write_clinical_warnings(warnings_path, output_rows)
    if n_flags:
        print("[majority_voting] %d clinical warning(s) written to %s" % (n_flags, warnings_path),
              file=sys.stderr)

    return 0


if __name__ == "__main__":
    sys.exit(main())
