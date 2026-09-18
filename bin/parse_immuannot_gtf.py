#!/usr/bin/env python3
"""
parse_immuannot_gtf.py — Convert an Immuannot GTF annotation into the standard
pipeline TSV plus a haplotype->IMGT allele crosswalk.

Immuannot annotates assembled contigs against IPD-IMGT/HLA gene features and
emits a (gzipped) GTF. Each gene feature carries attributes including the
assigned allele, the source contig, and a structural completeness status. We
extract:

1. Standard pipeline TSV (Gene/Allele1/Allele2) — up to two alleles per HLA gene,
   so an Immuannot annotation of a diploid assembly can enter the ensemble.
2. A crosswalk TSV (contig/haplotype <TAB> imgt_allele) usable as
   parse_locityper_results.py --crosswalk when the Locityper DB is HPRC-derived.

GTF attribute names vary across Immuannot releases; we parse the attribute
column flexibly, looking for allele/gene/contig/status-like keys.
"""

import argparse
import gzip
import re
import sys
from collections import defaultdict
from pathlib import Path

CLASSICAL = {"A", "B", "C", "DRB1", "DQA1", "DQB1", "DPA1", "DPB1"}
ATTR_RE = re.compile(r'(\w+)\s+"([^"]*)"')


def parse_args():
    p = argparse.ArgumentParser(description="Parse Immuannot GTF to standard format + crosswalk")
    p.add_argument("--input", required=True, help="Immuannot GTF (.gtf or .gtf.gz)")
    p.add_argument("--sample", required=True, help="Sample ID")
    p.add_argument("--output", required=True, help="Standard pipeline TSV")
    p.add_argument("--crosswalk-output", default=None, help="Optional haplotype->allele crosswalk TSV")
    return p.parse_args()


def open_maybe_gzip(path):
    path = Path(path)
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open("r", encoding="utf-8")


def parse_attributes(field):
    return dict(ATTR_RE.findall(field))


def gene_of(allele_or_gene):
    token = allele_or_gene.replace("HLA-", "").replace("HLA_", "")
    if "*" in token:
        return token.split("*")[0].upper()
    return token.upper()


def main():
    args = parse_args()
    # gene -> list of (allele, contig), preserving first-seen order
    gene_alleles = defaultdict(list)
    crosswalk = {}

    try:
        with open_maybe_gzip(args.input) as fh:
            for line in fh:
                if not line.strip() or line.startswith("#"):
                    continue
                cols = line.rstrip("\n").split("\t")
                if len(cols) < 9:
                    continue
                feature = cols[2].lower()
                # Annotate at the gene/transcript level, not every exon.
                if feature not in ("gene", "transcript", "mrna"):
                    continue
                attrs = parse_attributes(cols[8])
                allele = (attrs.get("allele") or attrs.get("Allele")
                          or attrs.get("gene_id") or attrs.get("transcript_id") or "").strip()
                gene = (attrs.get("gene") or attrs.get("gene_name") or "").strip()
                contig = cols[0].strip()
                if not allele:
                    continue
                gene_sym = gene_of(gene or allele)
                if gene_sym not in CLASSICAL:
                    continue
                pair = (allele, contig)
                if pair not in gene_alleles[gene_sym]:
                    gene_alleles[gene_sym].append(pair)
                if contig and allele:
                    crosswalk[contig] = allele
    except (OSError, ValueError) as exc:
        print(f"Warning: could not read {args.input}: {exc}", file=sys.stderr)

    with Path(args.output).open("w", encoding="utf-8") as out:
        out.write(f"# Immuannot results for {args.sample}\n")
        out.write("# Tool: Immuannot -- assembly gene-structure annotation\n")
        out.write("Gene\tAllele1\tAllele2\n")
        if not gene_alleles:
            out.write("# No alleles annotated\n")
        for gene in sorted(gene_alleles):
            alleles = gene_alleles[gene]
            a1 = alleles[0][0] if alleles else "NA"
            a2 = alleles[1][0] if len(alleles) > 1 else a1  # homozygous fallback
            out.write(f"{gene}\t{a1}\t{a2}\n")

    if args.crosswalk_output:
        with Path(args.crosswalk_output).open("w", encoding="utf-8") as out:
            out.write("# haplotype/contig\timgt_allele\n")
            for contig in sorted(crosswalk):
                out.write(f"{contig}\t{crosswalk[contig]}\n")

    print(f"Immuannot parsing complete: {len(gene_alleles)} genes annotated for {args.sample}")


if __name__ == "__main__":
    main()
