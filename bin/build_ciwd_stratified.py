#!/usr/bin/env python3
"""Build CIWD-stratified concordance + plausibility summaries for the champHLA report.

Reuses the vendored CIWD 3.0.0 catalogue (bin/ciwd.py, assets/ciwd_3.0.0.tsv) to annotate
*already-harmonized* benchmark rows — no HLA-typing tool is re-run. It aggregates two outputs the
July-2026 progress report embeds:

  1. summary_ciwd_stratified.tsv  — two-field concordance by (method, modality, CIWD stratum), for
     the three consensus methods (ChampHLA/ChampionChallenger, WeightedConsensus, MajorityVote).
     Answers the reviewer question "does the consensus win only on common alleles?"
  2. summary_ciwd_plausibility.tsv — per (tool, modality) QC: how often a callable call is
     biologically implausible (not-CIWD / novel), and whether such calls are error-enriched.

Sources are the coherent full-cohort 1000G truth-backed runs (RNA, WES) plus the complete WGS
wave-2 run — each carries champion_challenger_calls / weighted_consensus_calls /
majority_vote_baseline / harmonized_benchmark_rows in one directory, so method-level and per-tool
views come from the same samples.

CIWD is an allele-classification catalogue, NOT a ground-truth cohort (assets/README_CIWD.md); this
is an annotation layer over the existing ground-truth benchmark.
"""
import sys
from collections import defaultdict
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "bin"))
import ciwd  # noqa: E402

STRATA_ORDER = ["common", "intermediate", "well_documented", "not_ciwd", "unknown"]
_STRAT_RANK = {s: i for i, s in enumerate(STRATA_ORDER)}
_MOD_RANK = {"rnaseq": 0, "rna": 0, "wes": 1, "wgs": 2}

# Consensus call table -> canonical method label the figure expects.
METHOD_FILES = {
    "champion_challenger_calls.tsv": "ChampionChallenger",
    "weighted_consensus_calls.tsv": "WeightedConsensus",
    "majority_vote_baseline.tsv": "MajorityVote",
}

# Coherent full-cohort sources (one dir per modality with all method + harmonized tables).
SOURCE_DIRS = {
    "rna": REPO / "analysis/1000g_realdata/benchmark_rna_truthbacked/run/tables",
    "wes": REPO / "analysis/1000g_realdata/benchmark_wes_truthbacked/run/tables",
    "wgs": REPO / "analysis/1000g_realdata/benchmark_wgs_wave2/tables",
}

OUT_DIR = REPO / "analysis/ciwd_stratified/tables"


def read_tsv(path):
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            yield dict(zip(header, line.rstrip("\n").split("\t")))


def write_tsv(path, rows, columns):
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as fh:
        fh.write("\t".join(columns) + "\n")
        for r in rows:
            fh.write("\t".join(str(r.get(c, "")) for c in columns) + "\n")


def build_stratified(cat):
    """(method, modality, stratum) -> gene_rows / callable_rate / overall_correct_call_rate."""
    agg = defaultdict(lambda: {"gene_rows": 0, "callable": 0, "correct": 0})
    for _mod, tables in SOURCE_DIRS.items():
        for fname, method in METHOD_FILES.items():
            path = tables / fname
            if not path.exists():
                continue
            for row in read_tsv(path):
                modality = row.get("modality", "")
                stratum = cat.genotype_stratum(row.get("truth_allele1", ""), row.get("truth_allele2", ""))
                b = agg[(method, modality, stratum)]
                b["gene_rows"] += 1
                b["callable"] += 1 if row.get("is_callable", "0") == "1" else 0
                b["correct"] += 1 if row.get("is_correct", "0") == "1" else 0
    out = []
    for (method, modality, stratum), b in agg.items():
        n = b["gene_rows"]
        out.append({
            "method": method,
            "modality": modality,
            "ciwd_stratum": stratum,
            "gene_rows": n,
            "callable_rate": round(b["callable"] / n, 4) if n else "",
            "overall_correct_call_rate": round(b["correct"] / n, 4) if n else "",
        })
    out.sort(key=lambda r: (r["method"], _MOD_RANK.get(r["modality"], 9),
                            _STRAT_RANK.get(r["ciwd_stratum"], 99)))
    return out


def build_plausibility(cat):
    """(tool, modality) -> implausible-call rate + error-enrichment, from harmonized per-tool rows."""
    agg = defaultdict(lambda: {"callable": 0, "implausible": 0, "implausible_incorrect": 0})
    for tables in SOURCE_DIRS.values():
        path = tables / "harmonized_benchmark_rows.tsv"
        if not path.exists():
            continue
        for row in read_tsv(path):
            if row.get("is_callable", "0") != "1":
                continue
            b = agg[(row.get("tool", ""), row.get("modality", ""))]
            b["callable"] += 1
            implausible = cat.is_implausible(row.get("allele1", "")) or cat.is_implausible(row.get("allele2", ""))
            if implausible:
                b["implausible"] += 1
                if row.get("is_correct", "0") != "1":
                    b["implausible_incorrect"] += 1
    out = []
    for (tool, modality), b in agg.items():
        n, imp = b["callable"], b["implausible"]
        out.append({
            "tool": tool,
            "modality": modality,
            "callable_calls": n,
            "implausible_calls": imp,
            "implausible_rate": round(imp / n, 4) if n else "",
            "implausible_incorrect": b["implausible_incorrect"],
            "implausible_error_rate": round(b["implausible_incorrect"] / imp, 4) if imp else "",
        })
    out.sort(key=lambda r: (_MOD_RANK.get(r["modality"], 9), r["tool"]))
    return out


def main():
    cat = ciwd.load_ciwd()
    if len(cat) == 0:
        raise SystemExit("CIWD catalogue is empty — check assets/ciwd_3.0.0.tsv")

    strat = build_stratified(cat)
    plaus = build_plausibility(cat)

    strat_path = OUT_DIR / "summary_ciwd_stratified.tsv"
    plaus_path = OUT_DIR / "summary_ciwd_plausibility.tsv"
    write_tsv(strat_path, strat,
              ["method", "modality", "ciwd_stratum", "gene_rows", "callable_rate", "overall_correct_call_rate"])
    write_tsv(plaus_path, plaus,
              ["tool", "modality", "callable_calls", "implausible_calls", "implausible_rate",
               "implausible_incorrect", "implausible_error_rate"])

    print(f"CIWD catalogue: {len(cat)} two-field alleles")
    print(f"Wrote {strat_path} ({len(strat)} rows)")
    print(f"Wrote {plaus_path} ({len(plaus)} rows)")
    mods = sorted({r["modality"] for r in strat}, key=lambda m: _MOD_RANK.get(m, 9))
    methods = sorted({r["method"] for r in strat})
    print(f"  modalities: {mods}")
    print(f"  methods:    {methods}")


if __name__ == "__main__":
    main()
