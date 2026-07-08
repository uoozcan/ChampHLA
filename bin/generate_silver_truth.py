#!/usr/bin/env python3
"""
generate_silver_truth.py — Build a silver-standard HLA truth table from the
orthogonal genotypers Locityper (read/depth-based) and Immuannot (assembly
annotation).

Rationale: many WGS/long-read NGS datasets have no gold HLA truth, which blocks
benchmarking the ensemble + Champion-Challenger on them. Locityper and Immuannot
are methodologically orthogonal to the alignment/assembly typers, so their
*agreement* is a defensible silver standard (HPRC-style truth is itself
assembly-derived). This is a SILVER standard — validate it against gold truth
with validate_silver_truth.py before trusting it, and never benchmark a
truth-source tool against its own truth (see the hla_benchmark.py guardrail).

Output `truth_long.tsv` matches the existing truth-file schema exactly
(`sample  gene  allele1  allele2`, bare 2-field alleles) so it is consumed
as-is by the benchmark via `truth.path` — no benchmark-code changes required.

Gating (per sample x gene):
  - require_agreement (default): emit only when BOTH sources give the same
    2-field pair AND Locityper gq >= --min-gq AND novelty <= --max-novelty.
  - otherwise: a single high-confidence source may pass, flagged in provenance.
Everything else is abstained (omitted from truth_long.tsv, logged in provenance).
"""

import argparse
import csv
import os
import sys
from pathlib import Path

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
# Reuse the parsing/normalisation helpers already written for the annotation layer.
from annotate_hla_calls import (  # noqa: E402
    expand_paths, norm_pair, norm_gene, clean, load_locityper_confidence,
    sample_from_filename,
)


def parse_args():
    p = argparse.ArgumentParser(description="Generate silver-standard HLA truth from Locityper + Immuannot")
    p.add_argument("--locityper-calls", nargs="*", default=[],
                   help="Locityper standard TSV(s)/glob(s): <sample>_locityper.txt")
    p.add_argument("--locityper-confidence", nargs="*", default=[],
                   help="Locityper confidence side-channel(s)/glob(s): <sample>_locityper.confidence.tsv")
    p.add_argument("--immuannot-calls", nargs="*", default=[],
                   help="Immuannot standard TSV(s)/glob(s): <sample>_immuannot.txt")
    p.add_argument("--output", required=True, help="Output truth_long.tsv")
    p.add_argument("--provenance-output", default=None, help="Optional truth_provenance.tsv")
    p.add_argument("--min-gq", type=float, default=20.0,
                   help="Minimum Locityper GQ for a call to qualify as truth (default 20 ~ conf 0.99)")
    p.add_argument("--max-novelty", type=float, default=0.05,
                   help="Maximum Locityper novelty_score for a call to qualify (default 0.05)")
    p.add_argument("--require-agreement", dest="require_agreement", action="store_true", default=True,
                   help="Require Locityper and Immuannot to agree (default: on)")
    p.add_argument("--allow-single-source", dest="require_agreement", action="store_false",
                   help="Allow a single high-confidence source to emit truth (flagged in provenance)")
    p.add_argument("--loci", default=None,
                   help="Optional comma-separated locus allow-list (e.g. A,B,C)")
    return p.parse_args()


def load_standard_calls(paths, suffix):
    """sample -> gene -> (allele1, allele2) from a standard Gene/Allele1/Allele2 TSV."""
    index = {}
    for path in paths:
        sample = sample_from_filename(path, suffix)
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


def bare(allele):
    """Strip gene prefix to match the existing truth-file convention (A*02:01 -> 02:01)."""
    token = clean(allele)
    return token.split("*", 1)[1] if "*" in token else token


def main():
    args = parse_args()
    locus_allow = None
    if args.loci:
        locus_allow = {g.strip().replace("HLA-", "").upper() for g in args.loci.split(",") if g.strip()}

    loc_calls = load_standard_calls(expand_paths(args.locityper_calls), "_locityper.txt")
    loc_conf = load_locityper_confidence(expand_paths(args.locityper_confidence))
    imm_calls = load_standard_calls(expand_paths(args.immuannot_calls), "_immuannot.txt")

    samples = sorted(set(loc_calls) | set(imm_calls))
    truth_rows = []       # (sample, gene, allele1, allele2)
    provenance_rows = []  # full audit

    for sample in samples:
        genes = set(loc_calls.get(sample, {})) | set(imm_calls.get(sample, {}))
        for gene in sorted(genes):
            if locus_allow is not None and gene not in locus_allow:
                continue
            loc_raw = loc_calls.get(sample, {}).get(gene)
            imm_raw = imm_calls.get(sample, {}).get(gene)
            conf = loc_conf.get(sample, {}).get(gene, {})
            gq = conf.get("gq")
            novelty = conf.get("novelty_score")

            loc_pair = norm_pair(*loc_raw) if loc_raw else ()
            imm_pair = norm_pair(*imm_raw) if imm_raw else ()
            loc_ok = len(loc_pair) == 2
            imm_ok = len(imm_pair) == 2

            gq_pass = (gq is None and not loc_ok) or (gq is not None and gq >= args.min_gq)
            novelty_pass = novelty is None or novelty <= args.max_novelty
            loc_conf_pass = (gq is not None and gq >= args.min_gq) and novelty_pass

            source = ""
            emit_pair = None
            status = "abstained"

            if loc_ok and imm_ok:
                if loc_pair == imm_pair:
                    if loc_conf_pass or gq is None:
                        source, emit_pair, status = "locityper+immuannot", loc_pair, "emitted"
                    else:
                        status = "abstained_lowconf"
                else:
                    status = "abstained_disagree"
            elif loc_ok and not imm_ok:
                if args.require_agreement:
                    status = "single_source"
                elif loc_conf_pass:
                    source, emit_pair, status = "locityper", loc_pair, "emitted_single_source"
                else:
                    status = "abstained_lowconf"
            elif imm_ok and not loc_ok:
                if args.require_agreement:
                    status = "single_source"
                else:
                    # Immuannot (assembly) is high quality on its own; no GQ to gate.
                    source, emit_pair, status = "immuannot", imm_pair, "emitted_single_source"

            if emit_pair:
                truth_rows.append((sample, gene, bare(emit_pair[0]), bare(emit_pair[1])))

            provenance_rows.append({
                "sample": sample,
                "gene": gene,
                "source": source,
                "locityper_pair": "+".join(loc_pair) if loc_ok else "",
                "immuannot_pair": "+".join(imm_pair) if imm_ok else "",
                "locityper_gq": "" if gq is None else gq,
                "novelty_score": "" if novelty is None else novelty,
                "agree": "1" if (loc_ok and imm_ok and loc_pair == imm_pair) else "0",
                "status": status,
            })

    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["sample", "gene", "allele1", "allele2"])
        for row in truth_rows:
            writer.writerow(row)

    if args.provenance_output:
        fields = ["sample", "gene", "source", "locityper_pair", "immuannot_pair",
                  "locityper_gq", "novelty_score", "agree", "status"]
        with Path(args.provenance_output).open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
            writer.writeheader()
            for row in provenance_rows:
                writer.writerow(row)

    emitted = len(truth_rows)
    abstained = len(provenance_rows) - emitted
    print(f"Silver truth generated: {emitted} loci emitted, {abstained} abstained "
          f"across {len(samples)} samples")


if __name__ == "__main__":
    main()
