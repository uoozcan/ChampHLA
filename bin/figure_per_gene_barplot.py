#!/usr/bin/env python3.11
"""
figure_per_gene_barplot.py — Per-locus accuracy barplot for PIHLA benchmark.

Three subplots (HLA-A, HLA-B, HLA-C), each showing overall_correct_call_rate
per tool with grouped bars for WES, WGS, and RNAseq.

Usage:
    python3.11 bin/figure_per_gene_barplot.py
    python3.11 bin/figure_per_gene_barplot.py --tables_dir /path/to/tables --out_dir /path/to/out
"""

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

TABLES_DIR_DEFAULT = (
    Path(__file__).parent.parent
    / "analysis/benchmark_trimodal_50samples/tables"
)
OUT_DIR_DEFAULT = Path(__file__).parent.parent / "analysis/figures_final"

TOOL_ORDER = [
    "ArcasHLA", "HLA-HD", "Kourami", "OptiType", "POLYSOLVER",
    "Seq2HLA", "SpecHLA", "T1K", "MajorityVote", "WeightedConsensus",
]

MODALITIES = ["wes", "wgs", "rnaseq"]
MODALITY_LABELS = {"wes": "WES", "wgs": "WGS", "rnaseq": "RNAseq"}

# Wong colorblind-safe palette
MODALITY_COLORS = {
    "wes":    "#E69F00",
    "wgs":    "#56B4E9",
    "rnaseq": "#009E73",
}

GENES = ["A", "B", "C"]
GENE_TITLES = {"A": "HLA-A", "B": "HLA-B", "C": "HLA-C"}

BAR_WIDTH = 0.22
DPI = 300
WATERMARK = "PIHLA v2.0 · IMGT/HLA 3.59.0"

# ---------------------------------------------------------------------------
# Style
# ---------------------------------------------------------------------------

def apply_style():
    plt.rcParams.update({
        "font.family":        "DejaVu Sans",
        "font.size":          9,
        "axes.titlesize":     12,
        "axes.titleweight":   "bold",
        "axes.labelsize":     9,
        "xtick.labelsize":    8,
        "ytick.labelsize":    8,
        "axes.spines.top":    False,
        "axes.spines.right":  False,
        "axes.grid":          True,
        "axes.grid.axis":     "y",
        "grid.color":         "#e5e7eb",
        "grid.linewidth":     0.6,
        "legend.fontsize":    8,
        "legend.frameon":     False,
        "figure.dpi":         DPI,
        "savefig.dpi":        DPI,
        "savefig.bbox":       "tight",
        "savefig.pad_inches": 0.05,
        "pdf.fonttype":       42,
        "ps.fonttype":        42,
    })

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--tables_dir", type=Path, default=TABLES_DIR_DEFAULT)
    p.add_argument("--out_dir", type=Path, default=OUT_DIR_DEFAULT)
    return p.parse_args()


def load_data(tables_dir: Path) -> pd.DataFrame:
    fpath = tables_dir / "method_per_gene.tsv"
    df = pd.read_csv(fpath, sep="\t")
    df["overall_correct_call_rate"] = pd.to_numeric(
        df["overall_correct_call_rate"], errors="coerce"
    )
    return df


def draw_subplot(ax, acc: dict, counts: dict, gene: str, is_leftmost: bool):
    """Draw grouped barplot for one gene on ax.

    acc:    dict keyed (tool, modality) -> accuracy value
    counts: dict keyed (tool, modality) -> sample_count
    """
    n_tools = len(TOOL_ORDER)
    n_mod = len(MODALITIES)
    offsets = np.array([(i - (n_mod - 1) / 2) * BAR_WIDTH for i in range(n_mod)])
    x_centers = np.arange(n_tools)

    for mod_idx, mod in enumerate(MODALITIES):
        for tool_idx, tool in enumerate(TOOL_ORDER):
            val = acc.get((tool, mod))
            if val is None:
                continue
            xpos = x_centers[tool_idx] + offsets[mod_idx]
            ax.bar(
                xpos, val,
                width=BAR_WIDTH * 0.92,
                color=MODALITY_COLORS[mod],
                edgecolor="white",
                linewidth=0.4,
                zorder=3,
            )
            n = counts.get((tool, mod))
            if n is not None:
                ax.text(
                    xpos, val + 0.02, f"n={int(n)}",
                    ha="center", va="bottom",
                    fontsize=5.5, color="#374151",
                    rotation=90, zorder=4,
                )

    # Separator line between single tools and ensemble/baseline
    # MajorityVote is at index 8 in TOOL_ORDER (0-based)
    sep_idx = TOOL_ORDER.index("MajorityVote")
    ax.axvline(
        x=sep_idx - 0.5,
        color="#9ca3af",
        linewidth=0.8,
        linestyle="--",
        zorder=2,
    )

    ax.set_title(GENE_TITLES[gene])
    ax.set_xticks(x_centers)
    ax.set_xticklabels(TOOL_ORDER, rotation=45, ha="right")
    ax.set_xlim(-0.6, n_tools - 0.4)
    ax.set_ylim(0, 1.22)
    ax.set_yticks([0.0, 0.2, 0.4, 0.6, 0.8, 1.0])

    if is_leftmost:
        ax.set_ylabel("Accuracy (correct call rate)")

    # Light grey shading for ensemble group
    ax.axvspan(sep_idx - 0.5, n_tools - 0.4, color="#f3f4f6", zorder=1)


def make_figure(df: pd.DataFrame, out_dir: Path):
    apply_style()

    # Build lookups per gene: (tool, modality) -> accuracy / sample_count
    acc_lookup: dict[str, dict[tuple, float]] = {g: {} for g in GENES}
    cnt_lookup: dict[str, dict[tuple, int]]   = {g: {} for g in GENES}
    for _, row in df.iterrows():
        tool = row["method"]
        mod  = row["modality"]
        gene = row["gene"]
        val  = row["overall_correct_call_rate"]
        n    = row["sample_count"]
        if tool in TOOL_ORDER and mod in MODALITIES and gene in GENES and pd.notna(val):
            acc_lookup[gene][(tool, mod)] = float(val)
            cnt_lookup[gene][(tool, mod)] = int(n)

    fig, axes = plt.subplots(1, 3, figsize=(16, 4.5))

    for i, gene in enumerate(GENES):
        draw_subplot(axes[i], acc_lookup[gene], cnt_lookup[gene], gene, is_leftmost=(i == 0))

    # Shared legend below subplots
    legend_handles = [
        mpatches.Patch(facecolor=MODALITY_COLORS[m], label=MODALITY_LABELS[m])
        for m in MODALITIES
    ]
    fig.legend(
        handles=legend_handles,
        loc="lower center",
        ncol=3,
        bbox_to_anchor=(0.5, -0.02),
        frameon=False,
        fontsize=9,
    )

    fig.text(
        0.99, 0.005, WATERMARK,
        ha="right", va="bottom",
        fontsize=7, color="#9ca3af",
        transform=fig.transFigure,
    )

    fig.tight_layout(rect=[0, 0.06, 1, 1])

    out_dir.mkdir(parents=True, exist_ok=True)
    stem = "figure_per_gene_barplot"
    fig.savefig(str(out_dir / f"{stem}.png"), dpi=DPI)
    fig.savefig(str(out_dir / f"{stem}.pdf"))
    plt.close(fig)
    print(f"Saved: {out_dir / stem}.png / .pdf")


def main():
    args = parse_args()
    df = load_data(args.tables_dir)
    make_figure(df, args.out_dir)


if __name__ == "__main__":
    main()
