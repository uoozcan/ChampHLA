#!/usr/bin/env python3
"""Populate the empty champion_challenger_* stubs in the *_cv_recalibrated dirs
from the nested-CV outputs, so the figure generators and downstream consumers
read the honest held-out ChampionChallenger numbers (satisfies plan item A2).

The nested-CV method name 'ChampionChallenger_nestedCV' is written back as
'ChampionChallenger' to match the schema the generators expect.
"""
import csv
from pathlib import Path

ROOT = Path("/scratch/project_2008084/pihla-publish/analysis")
NCV = ROOT / "nested_cv_champion_challenger"
MODS = {"wgs": "wgs", "wes": "wes", "rna": "rnaseq"}  # dir label -> modality token


def read_tsv(p):
    with open(p) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def write_tsv(p, rows, cols):
    with open(p, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t")
        w.writeheader()
        for r in rows:
            w.writerow({c: r.get(c, "") for c in cols})


for lbl, mod in MODS.items():
    ncv_dir = NCV / lbl
    dst = ROOT / f"benchmark_{lbl}_cv_recalibrated" / "tables"
    if not (ncv_dir / "nested_cv_method_comparison.tsv").exists() or not dst.exists():
        print(f"skip {lbl}: missing inputs")
        continue

    comp = read_tsv(ncv_dir / "nested_cv_method_comparison.tsv")
    cc = next(r for r in comp if r["method"].startswith("ChampionChallenger"))
    comp_cols = ["method", "method_type", "modality", "sample_count", "gene_rows",
                 "callable_rate", "accuracy_among_callable", "overall_correct_call_rate",
                 "overall_correct_call_rate_ci_lo", "overall_correct_call_rate_ci_hi"]
    comp_row = {
        "method": "ChampionChallenger", "method_type": "ensemble", "modality": mod,
        "sample_count": "", "gene_rows": cc["n_holdout_gene_rows"],
        "callable_rate": cc["callable_rate"], "accuracy_among_callable": cc["accuracy_among_callable"],
        "overall_correct_call_rate": cc["overall_correct_call_rate"],
        "overall_correct_call_rate_ci_lo": cc["ci_lo"], "overall_correct_call_rate_ci_hi": cc["ci_hi"],
    }
    write_tsv(dst / "champion_challenger_method_comparison.tsv", [comp_row], comp_cols)

    # per-gene
    pg = read_tsv(ncv_dir / "nested_cv_per_gene.tsv")
    pg_cols = ["method", "method_type", "modality", "gene", "sample_count", "gene_rows",
               "callable_rate", "accuracy_among_callable", "overall_correct_call_rate"]
    pg_rows = []
    for r in pg:
        if not r["method"].startswith("ChampionChallenger"):
            continue
        pg_rows.append({
            "method": "ChampionChallenger", "method_type": "ensemble", "modality": mod,
            "gene": r["gene"], "sample_count": "", "gene_rows": r["n_gene_rows"],
            "callable_rate": r["callable_rate"], "accuracy_among_callable": r["accuracy_among_callable"],
            "overall_correct_call_rate": r["overall_correct_call_rate"],
        })
    write_tsv(dst / "champion_challenger_method_per_gene.tsv", pg_rows, pg_cols)
    print(f"populated {lbl}: CC overall={comp_row['overall_correct_call_rate']} + {len(pg_rows)} per-gene rows")
