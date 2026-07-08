#!/usr/bin/env python3
"""
annotate_hla_calls.py — Interpretive annotation layer over CHAMPHLA consensus /
Champion-Challenger calls.

Combines two sources of interpretive value identified for CHAMPHLA:

1. Immuannot-style functional annotation (IPD-IMGT/HLA expression status):
   flags null (N), low (L), secreted (S), cytoplasm (C), aberrant (A) and
   questionable (Q) expression alleles, derived from the allele-name suffix
   and/or an optional IPD-IMGT/HLA allele-status table. These are clinically
   critical (e.g. a null allele means no surface expression) yet invisible in a
   bare consensus call.

2. Locityper cross-track validation: joins the independent, depth-aware
   sequence-level genotype (and its Phred GQ / novelty score) by sample+gene to
       - flag NOVEL_ALLELE_CANDIDATE where no DB haplotype fits well,
       - flag CROSS_TRACK_DISCORD where the sequence-level call contradicts the
         consensus at high confidence,
       - flag ALLELE_DROPOUT_SUSPECTED where a homozygous consensus call is
         contradicted by a heterozygous sequence-level call (true homozygosity
         vs. allele dropout — hardening the heuristic PUTATIVE_HOMOZYGOSITY).

Warning codes reuse majority_voting.CLINICAL_WARNING_DESCRIPTIONS (single source
of truth). Output mirrors the input calls TSV with added annotation columns and
an augmented `warnings` column.
"""

import argparse
import csv
import glob
import os
import re
import sys
from pathlib import Path

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from majority_voting import CLINICAL_WARNING_DESCRIPTIONS  # noqa: E402

# IPD-IMGT/HLA expression suffixes. NOTE: G/P are group designators, NOT
# expression, and are deliberately excluded.
EXPRESSION_LABELS = {
    "N": ("null", "NULL_ALLELE"),
    "L": ("low", "LOW_EXPRESSION"),
    "S": ("secreted", "LOW_EXPRESSION"),
    "C": ("cytoplasm", "LOW_EXPRESSION"),
    "A": ("aberrant", "LOW_EXPRESSION"),
    "Q": ("questionable", "LOW_EXPRESSION"),
}
EXPRESSION_SUFFIX_RE = re.compile(r"([NLSCAQ])$")

MISSING = {"", "-", "NA", "None", "none", "."}


def parse_args():
    p = argparse.ArgumentParser(description="Annotate CHAMPHLA calls with expression status and Locityper cross-validation")
    p.add_argument("--calls", required=True, help="Consensus / Champion-Challenger calls TSV")
    p.add_argument("--output", required=True, help="Annotated output TSV")
    p.add_argument("--locityper-calls", nargs="*", default=[],
                   help="Locityper standard TSV(s)/glob(s): <sample>_locityper.txt (for pair comparison)")
    p.add_argument("--locityper-confidence", nargs="*", default=[],
                   help="Locityper confidence side-channel TSV(s)/glob(s): <sample>_locityper.confidence.tsv")
    p.add_argument("--imgt-allele-status", default=None,
                   help="Optional IPD-IMGT/HLA allele-status TSV (allele<TAB>status)")
    p.add_argument("--novelty-threshold", type=float, default=0.5,
                   help="novelty_score at/above which to flag NOVEL_ALLELE_CANDIDATE (default 0.5 ~ GQ<=3)")
    p.add_argument("--gq-discord-min", type=float, default=10.0,
                   help="Minimum Locityper GQ to trust a disagreement enough to flag discord (default 10)")
    p.add_argument("--clinical-warnings-output", default=None,
                   help="Optional human-readable report of flagged rows only")
    return p.parse_args()


def clean(value):
    return (value or "").strip()


def is_missing(value):
    return clean(value) in MISSING


def expand_paths(patterns):
    out = []
    for pattern in patterns:
        matches = glob.glob(pattern, recursive=True)
        out.extend(matches if matches else ([pattern] if Path(pattern).exists() else []))
    return out


def sample_from_filename(path, suffix):
    name = Path(path).name
    if name.endswith(suffix):
        return name[: -len(suffix)]
    return Path(path).stem


def norm_gene(gene):
    return clean(gene).replace("HLA-", "").upper()


def norm_pair(a1, a2, resolution=2):
    """Two-field, order-independent allele pair for comparison."""
    def trunc(allele):
        token = clean(allele).replace("HLA-", "").replace("HLA_", "")
        if "*" not in token:
            return ""
        gene, fields = token.split("*", 1)
        fields = EXPRESSION_SUFFIX_RE.sub("", fields)  # drop expression suffix for matching
        parts = [p for p in fields.split(":") if p][:resolution]
        return "%s*%s" % (gene.upper(), ":".join(parts)) if parts else ""
    pair = sorted(p for p in (trunc(a1), trunc(a2)) if p)
    return tuple(pair)


def load_imgt_status(path):
    mapping = {}
    if not path:
        return mapping
    with open(path, "r", encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) >= 2 and parts[1].strip():
                mapping[parts[0].strip()] = parts[1].strip()[:1].upper()
    return mapping


def expression_letter(allele, imgt_status):
    """Return the expression suffix letter for an allele, or '' if normal/unknown."""
    token = clean(allele)
    if not token or "*" not in token:
        return ""
    # Direct suffix on the allele name.
    fields = token.split("*", 1)[1]
    match = EXPRESSION_SUFFIX_RE.search(fields)
    if match:
        return match.group(1)
    # Fall back to the IMGT status table (exact, then 2-field key).
    if token in imgt_status:
        return imgt_status[token]
    key2 = norm_pair(token, token)
    if key2 and key2[0] in imgt_status:
        return imgt_status[key2[0]]
    return ""


def load_locityper_calls(paths):
    """sample -> gene -> (allele1, allele2) from Locityper standard TSVs."""
    index = {}
    for path in paths:
        sample = sample_from_filename(path, "_locityper.txt")
        with open(path, "r", encoding="utf-8") as fh:
            for raw in fh:
                line = raw.strip()
                if not line or line.startswith("#"):
                    continue
                parts = line.split("\t")
                if len(parts) < 3 or parts[0].lower() == "gene":
                    continue
                gene = norm_gene(parts[0])
                index.setdefault(sample, {})[gene] = (parts[1].strip(), parts[2].strip())
    return index


def load_locityper_confidence(paths):
    """sample -> gene -> {gq, novelty_score, mean_depth} from confidence side-channels."""
    index = {}
    for path in paths:
        fallback_sample = sample_from_filename(path, "_locityper.confidence.tsv")
        with open(path, "r", encoding="utf-8") as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for row in reader:
                sample = clean(row.get("sample")) or fallback_sample
                gene = norm_gene(row.get("gene"))
                if not gene:
                    continue
                index.setdefault(sample, {})[gene] = {
                    "gq": _to_float(row.get("gq")),
                    "novelty_score": _to_float(row.get("novelty_score")),
                    "mean_depth": _to_float(row.get("mean_depth")),
                }
    return index


def _to_float(value):
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def annotate_row(row, loc_calls, loc_conf, imgt_status, args):
    """Return (annotation_dict, list_of_added_warning_codes)."""
    sample = clean(row.get("sample"))
    gene = norm_gene(row.get("gene"))
    a1 = clean(row.get("allele1"))
    a2 = clean(row.get("allele2"))
    added = []

    # --- Expression / functional annotation -------------------------------
    letters = [expression_letter(a, imgt_status) for a in (a1, a2)]
    letters = [l for l in letters if l]
    expression_flag = ",".join(sorted(set(letters)))
    is_null = any(l == "N" for l in letters)
    notes = []
    for l in letters:
        label, code = EXPRESSION_LABELS[l]
        notes.append(label)
        if code not in added:
            added.append(code)
    functional_note = ",".join(notes)

    # --- Locityper cross-track validation ---------------------------------
    loc_pair = loc_calls.get(sample, {}).get(gene)
    conf = loc_conf.get(sample, {}).get(gene, {})
    gq = conf.get("gq")
    novelty = conf.get("novelty_score")

    if novelty is not None and novelty >= args.novelty_threshold:
        added.append("NOVEL_ALLELE_CANDIDATE")

    consensus_pair = norm_pair(a1, a2)
    if loc_pair and all(consensus_pair):
        loc_norm = norm_pair(loc_pair[0], loc_pair[1])
        trust = (gq is None) or (gq >= args.gq_discord_min)
        if all(loc_norm) and loc_norm != consensus_pair and trust:
            added.append("CROSS_TRACK_DISCORD")
            # Homozygous consensus contradicted by heterozygous sequence-level call.
            if consensus_pair[0] == consensus_pair[1] and loc_norm[0] != loc_norm[1]:
                added.append("ALLELE_DROPOUT_SUSPECTED")

    annotation = {
        "expression_flag": expression_flag,
        "is_null_allele": "1" if is_null else "0",
        "functional_note": functional_note,
        "locityper_allele1": loc_pair[0] if loc_pair else "",
        "locityper_allele2": loc_pair[1] if loc_pair else "",
        "locityper_gq": "" if gq is None else gq,
        "locityper_novelty_score": "" if novelty is None else novelty,
    }
    return annotation, added


def merge_warnings(existing, added):
    codes = [c for c in (existing or "").replace(";", ",").split(",") if c.strip()]
    for code in added:
        if code not in codes:
            codes.append(code)
    return ",".join(codes)


def main():
    args = parse_args()
    loc_calls = load_locityper_calls(expand_paths(args.locityper_calls))
    loc_conf = load_locityper_confidence(expand_paths(args.locityper_confidence))
    imgt_status = load_imgt_status(args.imgt_allele_status)

    with open(args.calls, "r", encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        rows = list(reader)
        in_fields = list(reader.fieldnames or [])

    annotation_fields = ["expression_flag", "is_null_allele", "functional_note",
                         "locityper_allele1", "locityper_allele2",
                         "locityper_gq", "locityper_novelty_score"]
    out_fields = list(in_fields)
    for f in annotation_fields:
        if f not in out_fields:
            out_fields.append(f)
    if "warnings" not in out_fields:
        out_fields.append("warnings")

    flagged = []
    out_rows = []
    for row in rows:
        annotation, added = annotate_row(row, loc_calls, loc_conf, imgt_status, args)
        row.update(annotation)
        row["warnings"] = merge_warnings(row.get("warnings"), added)
        out_rows.append(row)
        if row["warnings"]:
            flagged.append(row)

    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=out_fields, delimiter="\t")
        writer.writeheader()
        for row in out_rows:
            writer.writerow({f: row.get(f, "") for f in out_fields})

    if args.clinical_warnings_output:
        with open(args.clinical_warnings_output, "w", encoding="utf-8") as handle:
            handle.write("# Clinical warnings — flagged calls only\n")
            handle.write("sample\tmodality\tgene\tallele1\tallele2\twarnings\tdescriptions\n")
            for row in flagged:
                codes = [c for c in row["warnings"].split(",") if c]
                descs = "; ".join(CLINICAL_WARNING_DESCRIPTIONS.get(c, c) for c in codes)
                handle.write("\t".join([
                    clean(row.get("sample")), clean(row.get("modality")), clean(row.get("gene")),
                    clean(row.get("allele1")), clean(row.get("allele2")), row["warnings"], descs,
                ]) + "\n")

    print(f"Annotation complete: {len(out_rows)} calls, {len(flagged)} flagged")


if __name__ == "__main__":
    main()
