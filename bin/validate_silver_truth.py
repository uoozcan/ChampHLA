#!/usr/bin/env python3
"""
validate_silver_truth.py — Validate a generated silver-standard truth table
against a gold truth table, per locus.

Silver truth from generate_silver_truth.py is only trustworthy once we know how
well it agrees with gold truth on samples that have both. This compares two
`truth_long.tsv` files (same schema: sample, gene, allele1, allele2) at 2-field
resolution and reports per-locus concordance + coverage, so the generator's own
error/abstention rate is known before silver truth is used on truth-less cohorts.

Concordance is computed only over (sample, gene) pairs present in BOTH tables
(silver may legitimately abstain). We also report the abstention footprint
(gold-covered loci the silver generator did not emit).

Output `silver_truth_validation.tsv`:
    locus  n_compared  n_concordant  concordance_2field  n_silver_only  n_gold_only  n_abstained
plus an ALL summary row.
"""

import argparse
import csv
import importlib.util
import os
from collections import defaultdict
from pathlib import Path

# Reuse the canonical allele normalisation from the benchmark harness.
_HB_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), "hla_benchmark.py")
_spec = importlib.util.spec_from_file_location("hla_benchmark", _HB_PATH)
hla_benchmark = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(hla_benchmark)


def parse_args():
    p = argparse.ArgumentParser(description="Validate generated silver truth against gold truth")
    p.add_argument("--silver", required=True, help="Generated silver truth_long.tsv")
    p.add_argument("--gold", required=True, help="Gold truth_long.tsv")
    p.add_argument("--output", required=True, help="Per-locus validation TSV")
    p.add_argument("--resolution", type=int, default=2, help="Comparison resolution (default 2-field)")
    return p.parse_args()


def load_truth_pairs(path, resolution):
    """(sample, gene) -> normalised, order-independent allele pair."""
    pairs = {}
    with open(path, "r", encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        norm = {hla_benchmark.normalize_column_name(name): name for name in (reader.fieldnames or [])}
        sample_key = hla_benchmark.first_existing(norm, ["sample", "sampleid", "id"])
        gene_key = hla_benchmark.first_existing(norm, ["gene", "locus"])
        a1_key = hla_benchmark.first_existing(norm, ["allele1", "truthallele1", "typedallele1"])
        a2_key = hla_benchmark.first_existing(norm, ["allele2", "truthallele2", "typedallele2"])
        if not (sample_key and gene_key and a1_key and a2_key):
            raise SystemExit(f"Could not find sample/gene/allele columns in {path}")
        for row in reader:
            sample = hla_benchmark.clean_token(row[norm[sample_key]])
            gene = hla_benchmark.normalize_gene(row[norm[gene_key]])
            if not sample or not gene:
                continue
            pair = hla_benchmark.sort_alleles([row[norm[a1_key]], row[norm[a2_key]]], resolution=resolution)
            if all(pair):
                pairs[(sample, gene)] = pair
    return pairs


def main():
    args = parse_args()
    silver = load_truth_pairs(args.silver, args.resolution)
    gold = load_truth_pairs(args.gold, args.resolution)

    loci = sorted({g for (_, g) in silver} | {g for (_, g) in gold}, key=hla_benchmark.gene_sort_key)
    stats = defaultdict(lambda: {"n_compared": 0, "n_concordant": 0,
                                 "n_silver_only": 0, "n_gold_only": 0, "n_abstained": 0})

    all_keys = set(silver) | set(gold)
    for key in all_keys:
        sample, gene = key
        s = stats[gene]
        in_silver = key in silver
        in_gold = key in gold
        if in_silver and in_gold:
            s["n_compared"] += 1
            if silver[key] == gold[key]:
                s["n_concordant"] += 1
        elif in_silver and not in_gold:
            s["n_silver_only"] += 1
        elif in_gold and not in_silver:
            s["n_gold_only"] += 1
            s["n_abstained"] += 1  # gold-covered locus the generator did not emit

    fields = ["locus", "n_compared", "n_concordant", "concordance_2field",
              "n_silver_only", "n_gold_only", "n_abstained"]
    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    totals = {"n_compared": 0, "n_concordant": 0, "n_silver_only": 0, "n_gold_only": 0, "n_abstained": 0}
    with out_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for locus in loci:
            s = stats[locus]
            conc = round(s["n_concordant"] / s["n_compared"], 4) if s["n_compared"] else ""
            writer.writerow({"locus": locus, "concordance_2field": conc, **{k: s[k] for k in totals}})
            for k in totals:
                totals[k] += s[k]
        all_conc = round(totals["n_concordant"] / totals["n_compared"], 4) if totals["n_compared"] else ""
        writer.writerow({"locus": "ALL", "concordance_2field": all_conc, **totals})

    print(f"Silver-vs-gold validation: {totals['n_compared']} loci compared, "
          f"overall 2-field concordance "
          f"{(totals['n_concordant'] / totals['n_compared']) if totals['n_compared'] else 0:.4f}, "
          f"{totals['n_abstained']} gold loci abstained by generator")


if __name__ == "__main__":
    main()
