#!/usr/bin/env python3.11
"""
fig_fimm_loh_wes.py  [figure hub copy]

Stacked bar: LOH status distribution per HLA gene (A, B, C) for FIMM AML/MDS
samples using SpecHLA WES allele frequency data.

Original: bin/analysis/analyze_fimm_loh_wes.py  (classification + TSV output)
This hub copy focuses on the figure only; run the original for TSV/JSON outputs.

Hub standards: Arial font, legend outside axes, 300 DPI.

CAVEAT: WES allele frequencies reflect capture-panel read depth, not
genome-wide allele balance. Results are exploratory; WGS confirmation required.

Usage:
    python3.11 figure_hub/scripts/fig_fimm_loh_wes.py [--candidates PATH] [--outdir PATH]
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

# ── Hub style ──────────────────────────────────────────────────────────────
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from plot_style import apply_style, legend_outside, save_fig, WONG

apply_style()

import matplotlib.pyplot as plt  # after apply_style

# ── Paths ──────────────────────────────────────────────────────────────────
_LOCAL      = Path("/scratch/project_2008084/pihla_local")
DEFAULT_IN  = _LOCAL / "analysis" / "loh_analysis" / "fimm_loh_candidates_wes.tsv"
DEFAULT_IN_DEIDENT = _LOCAL / "analysis" / "loh_analysis" / "fimm_loh_candidates_wes_deident.tsv"
DEFAULT_OUT = _LOCAL / "analysis" / "figures_final"

GENE_ORDER  = ["A", "B", "C"]
STATUS_ORDER = ["candidate_loh", "allelic_imbalance_review",
                "balanced_heterozygous", "insufficient_evidence"]
STATUS_COLORS = {
    "candidate_loh":            WONG["vermil"],
    "allelic_imbalance_review": WONG["orange"],
    "balanced_heterozygous":    WONG["green"],
    "insufficient_evidence":    "#aaaaaa",
}
STATUS_LABELS = {
    "candidate_loh":            "Candidate LOH",
    "allelic_imbalance_review": "Allelic Imbalance (Review)",
    "balanced_heterozygous":    "Balanced Heterozygous",
    "insufficient_evidence":    "Insufficient Evidence",
}


def make_figure(df: pd.DataFrame, outdir: Path) -> None:
    counts = {g: {s: 0 for s in STATUS_ORDER} for g in GENE_ORDER}
    for _, row in df.iterrows():
        gene   = str(row.get("gene", ""))
        status = str(row.get("loh_status_wes", "insufficient_evidence"))
        if gene in counts and status in counts[gene]:
            counts[gene][status] += 1

    fig, ax = plt.subplots(figsize=(6, 4.5))
    x      = np.arange(len(GENE_ORDER))
    bottom = np.zeros(len(GENE_ORDER))

    for status in STATUS_ORDER:
        vals = np.array([counts[g][status] for g in GENE_ORDER], dtype=float)
        ax.bar(x, vals, bottom=bottom,
               color=STATUS_COLORS[status],
               label=STATUS_LABELS[status],
               edgecolor="white", linewidth=0.5)
        bottom += vals

    ax.set_xticks(x)
    ax.set_xticklabels([f"HLA-{g}" for g in GENE_ORDER])
    ax.set_ylabel("Number of samples")
    ax.set_title("FIMM WES — LOH Status Distribution per Gene\n"
                 "(SpecHLA WES allele frequencies, exploratory)",
                 loc="left", fontweight="bold")

    # Legend outside — right of axes
    legend_outside(ax, loc="upper left", anchor=(1.02, 1.0))

    ax.text(0.01, -0.22,
            "WES allele frequencies less reliable for LOH than WGS.\n"
            "WGS data required for confirmation of candidate LOH events.",
            transform=ax.transAxes, fontsize=6.5, color=WONG["vermil"], style="italic")

    fig.suptitle("Exploratory — not for clinical use without WGS validation",
                 fontsize=7.5, color="#888888", y=1.01)

    outdir.mkdir(parents=True, exist_ok=True)
    save_fig(fig, outdir / "figure_fimm_loh_wes")
    print(f"  Saved figure_fimm_loh_wes.{{pdf,svg,png}} → {outdir}")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--candidates", type=Path, default=DEFAULT_IN_DEIDENT if DEFAULT_IN_DEIDENT.exists() else DEFAULT_IN,
                    help="fimm_loh_candidates_wes.tsv (or _deident.tsv)")
    ap.add_argument("--outdir", type=Path, default=DEFAULT_OUT)
    args = ap.parse_args()
    if not args.candidates.exists():
        raise FileNotFoundError(f"Input not found: {args.candidates}")
    df = pd.read_csv(args.candidates, sep="\t", dtype=str)
    print(f"Loaded {len(df)} rows from {args.candidates.name}")
    make_figure(df, args.outdir)


if __name__ == "__main__":
    main()
