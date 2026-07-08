#!/usr/bin/env python3
"""
parse_locityper_results.py — Convert Locityper genotyping output to standard
pipeline format, plus a sequence-level confidence/novelty side-channel.

Locityper (Prodanov & Bansal) genotypes complex loci (HLA, KIR, ...) from
short- or long-read WGS by aligning recruited reads against a *locus pangenome*
of full-length haplotype sequences and selecting the haplotype pair that
maximises a likelihood combining alignment quality and read-depth profile.

Output layout consumed here (per `locityper genotype -o <out>`):
    <out>/loci/<locus>/res.json[.gz]    # predicted genotype + quality per locus

The per-locus JSON is parsed flexibly (schema has varied across releases). We
look for the predicted genotype (a comma-separated haplotype pair), a Phred-like
genotype quality (GQ), the weighted alignment distance, and a mean read depth.

Two outputs are written:

1. Standard pipeline TSV (matches Kourami/T1K format) — feeds aggregation/voting:
       # Locityper results for <sample>
       Gene    Allele1     Allele2
       A       A*02:01     A*11:01

2. Confidence/novelty side-channel TSV (consumed by the benchmark
   `long_confidence_table` parser and by annotate_hla_calls.py):
       sample  gene  confidence_score  read_support  gq  weighted_dist  mean_depth  novelty_score

   - confidence_score: GQ mapped to a probability via 1 - 10^(-GQ/10)
     (the natural reading of a Phred-scaled genotype quality).
   - novelty_score:    10^(-GQ/10), high when no DB haplotype fits well —
     the signal used to flag putative novel / structurally divergent alleles.

Haplotype-name -> IMGT allele crosswalk:
   When the Locityper DB is built from IPD-IMGT/HLA genomic (`*_gen`) sequences
   the haplotype names already carry the allele id, so a trailing assembly/DB
   tag (e.g. ".hap1", ":HG002") is stripped and the value passes straight to the
   downstream `normalize_allele()`. When the DB is HPRC-derived, supply a
   two-column crosswalk TSV (haplotype<TAB>imgt_allele) via --crosswalk.
"""

import argparse
import csv
import gzip
import json
import re
import sys
from pathlib import Path


def parse_args():
    parser = argparse.ArgumentParser(description="Parse Locityper output to standard pipeline format")
    parser.add_argument("--input", required=True,
                        help="Locityper genotyping output directory (expects <input>/loci/<locus>/res.json[.gz])")
    parser.add_argument("--sample", required=True, help="Sample ID")
    parser.add_argument("--output", required=True, help="Standard pipeline TSV (Gene/Allele1/Allele2)")
    parser.add_argument("--confidence-output", default=None,
                        help="Optional confidence/novelty side-channel TSV")
    parser.add_argument("--crosswalk", default=None,
                        help="Optional haplotype->IMGT allele crosswalk TSV (for HPRC-derived DBs)")
    parser.add_argument("--loci", default=None,
                        help="Optional comma-separated locus allow-list (e.g. A,B,C). Default: all loci found.")
    return parser.parse_args()


def open_maybe_gzip(path):
    """Open a path transparently whether or not it is gzip-compressed."""
    path = Path(path)
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open("r", encoding="utf-8")


def load_crosswalk(path):
    if not path:
        return {}
    mapping = {}
    with open(path, "r", encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) >= 2:
                mapping[parts[0].strip()] = parts[1].strip()
    return mapping


def coerce_float(value):
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def first_present(payload, keys):
    """Return the first present, non-null value among keys in a dict."""
    if not isinstance(payload, dict):
        return None
    for key in keys:
        if key in payload and payload[key] is not None:
            return payload[key]
    return None


def clean_haplotype_name(raw, crosswalk):
    """Map a Locityper haplotype name to an IMGT-style allele token."""
    token = str(raw).strip()
    if not token or token == "*":
        return ""
    if token in crosswalk:
        return crosswalk[token]
    # Strip a trailing assembly/DB tag: "A*02:01:01.hap1" or "A*02:01:01:HG002".
    token = re.split(r"[.](?:hap\d+|\d+)$", token)[0]
    token = token.split(" ")[0]
    return token


def extract_genotype_record(payload):
    """Pull the predicted genotype + quality fields from a per-locus JSON payload.

    Locityper (>=1.3) writes the predicted genotype + quality at the TOP LEVEL of
    res.json (`genotype`, `quality`, `weight_dist`, `total_reads`), with a ranked
    `options` list whose entries hold per-genotype likelihoods (no `quality`). We
    therefore read the headline fields from the top level, falling back to
    options[0] only for the genotype string if the top-level is absent.
    """
    if not isinstance(payload, dict):
        return None

    genotype = first_present(payload, ["genotype", "pred", "best", "alleles"])
    if genotype is None:
        options = first_present(payload, ["options", "genotypes", "predictions", "results"])
        if isinstance(options, list) and options and isinstance(options[0], dict):
            genotype = first_present(options[0], ["genotype", "pred", "best", "alleles"])

    # Locityper `quality` is Phred-scaled: quality = -10*log10(1 - prob_of_best).
    gq = coerce_float(first_present(payload, ["quality", "gq", "GQ", "qual"]))
    weighted_dist = coerce_float(first_present(payload, ["weight_dist", "weighted_dist", "dist", "distance"]))
    mean_depth = coerce_float(first_present(payload, ["mean_depth", "read_depth", "depth", "coverage", "total_reads"]))
    return {
        "genotype": genotype,
        "gq": gq,
        "weighted_dist": weighted_dist,
        "mean_depth": mean_depth,
    }


def split_genotype(genotype):
    """Split a Locityper genotype value into two haplotype names."""
    if genotype is None:
        return []
    if isinstance(genotype, (list, tuple)):
        tokens = list(genotype)
    else:
        tokens = re.split(r"[,/]", str(genotype))
    return [t.strip() for t in tokens if str(t).strip()]


def gene_from_allele(allele):
    """Extract the gene symbol from an IMGT allele token (e.g. A*02:01 -> A)."""
    token = allele.replace("HLA-", "").replace("HLA_", "")
    if "*" in token:
        return token.split("*")[0].upper()
    return token.upper()


def phred_to_probability(gq):
    """Phred-scaled genotype quality -> success probability (1 - 10^(-GQ/10))."""
    if gq is None:
        return None
    if gq < 0:
        gq = 0.0
    return round(1.0 - (10.0 ** (-gq / 10.0)), 6)


def novelty_from_quality(gq):
    """Complement of the Phred probability: high when no DB haplotype fits well."""
    if gq is None:
        return None
    if gq < 0:
        gq = 0.0
    return round(10.0 ** (-gq / 10.0), 6)


def collect_locus_results(input_dir):
    """Yield (locus_dir_name, payload) for each per-locus result JSON found."""
    loci_root = Path(input_dir) / "loci"
    search_root = loci_root if loci_root.is_dir() else Path(input_dir)
    for res_path in sorted(search_root.glob("*/res.json*")):
        if res_path.suffix not in (".json", ".gz"):
            continue
        try:
            with open_maybe_gzip(res_path) as fh:
                payload = json.load(fh)
        except (OSError, ValueError) as exc:
            print(f"Warning: could not read {res_path}: {exc}", file=sys.stderr)
            continue
        yield res_path.parent.name, payload


def main():
    args = parse_args()
    crosswalk = load_crosswalk(args.crosswalk)
    locus_allow = None
    if args.loci:
        locus_allow = {g.strip().replace("HLA-", "").upper() for g in args.loci.split(",") if g.strip()}

    # gene -> {allele1, allele2, gq, weighted_dist, mean_depth, novelty_score}
    calls = {}

    for locus_name, payload in collect_locus_results(args.input):
        record = extract_genotype_record(payload)
        if record is None:
            continue
        haplotypes = [clean_haplotype_name(h, crosswalk) for h in split_genotype(record["genotype"])]
        haplotypes = [h for h in haplotypes if h]
        if not haplotypes:
            continue

        gene = gene_from_allele(haplotypes[0])
        if locus_allow is not None and gene not in locus_allow:
            continue

        allele1 = haplotypes[0]
        allele2 = haplotypes[1] if len(haplotypes) > 1 else haplotypes[0]  # homozygous fallback
        calls[gene] = {
            "allele1": allele1,
            "allele2": allele2,
            "gq": record["gq"],
            "weighted_dist": record["weighted_dist"],
            "mean_depth": record["mean_depth"],
            "confidence_score": phred_to_probability(record["gq"]),
            "novelty_score": novelty_from_quality(record["gq"]),
        }

    # 1) Standard pipeline TSV
    out_path = Path(args.output)
    with out_path.open("w", encoding="utf-8") as out:
        out.write(f"# Locityper results for {args.sample}\n")
        out.write("# Tool: Locityper -- locus-pangenome genotyping (sequence-level)\n")
        out.write("Gene\tAllele1\tAllele2\n")
        if not calls:
            out.write("# No alleles called\n")
            print(f"Warning: no Locityper genotypes parsed from {args.input}", file=sys.stderr)
        for gene in sorted(calls):
            c = calls[gene]
            out.write(f"{gene}\t{c['allele1']}\t{c['allele2']}\n")

    # 2) Confidence / novelty side-channel TSV
    if args.confidence_output:
        fields = ["sample", "gene", "confidence_score", "read_support",
                  "gq", "weighted_dist", "mean_depth", "novelty_score"]
        with Path(args.confidence_output).open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
            writer.writeheader()
            for gene in sorted(calls):
                c = calls[gene]
                writer.writerow({
                    "sample": args.sample,
                    "gene": gene,
                    "confidence_score": "" if c["confidence_score"] is None else c["confidence_score"],
                    "read_support": "" if c["mean_depth"] is None else c["mean_depth"],
                    "gq": "" if c["gq"] is None else c["gq"],
                    "weighted_dist": "" if c["weighted_dist"] is None else c["weighted_dist"],
                    "mean_depth": "" if c["mean_depth"] is None else c["mean_depth"],
                    "novelty_score": "" if c["novelty_score"] is None else c["novelty_score"],
                })

    print(f"Locityper parsing complete: {len(calls)} loci typed for {args.sample}")


if __name__ == "__main__":
    main()
