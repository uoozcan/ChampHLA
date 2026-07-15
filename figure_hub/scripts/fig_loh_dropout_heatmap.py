#!/usr/bin/env python3.11
"""
fig_loh_dropout_heatmap.py  [figure hub copy]

Two-panel publication figure from benchmark_homozygosity_error_summary.tsv:
  Panel A — False Duplicate Rate heatmap (tool × modality-gene)
  Panel B — Allele dropout count by modality (stacked bar, coloured by tool)

Original: bin/analysis/figure_loh_allele_dropout_heatmap.py
Hub standards: Arial font, legend outside axes, 300 DPI.

Usage:
    python3.11 figure_hub/scripts/fig_loh_dropout_heatmap.py [--outdir PATH]
"""

import argparse
import sys
from pathlib import Path

import matplotlib.patches as mpatches
import numpy as np
import pandas as pd

# ── Hub style ──────────────────────────────────────────────────────────────
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from plot_style import apply_style, legend_outside, rotate_xticklabels, save_fig, tight_with_legend, TOOL_COLORS

apply_style()

import matplotlib.pyplot as plt  # after apply_style

# ── Paths ──────────────────────────────────────────────────────────────────
_SCRATCH    = Path("/scratch/project_2008084/pihla-publish/analysis")
DEFAULT_INPUT  = _SCRATCH / "loh_analysis" / "benchmark_homozygosity_error_summary.tsv"
DEFAULT_OUTDIR = _SCRATCH / "figures_final"

TOOL_ORDER     = ["ArcasHLA", "HLA-HD", "Kourami", "OptiType", "POLYSOLVER",
                  "Seq2HLA", "SpecHLA", "T1K"]
MOD_LABELS     = {"wgs": "WGS", "wes": "WES", "rnaseq": "RNA"}
GENE_ORDER     = ["A", "B", "C"]
MODALITY_ORDER = ["wgs", "wes", "rnaseq"]


def load_data(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t")
    df["false_duplicate_rate"]  = pd.to_numeric(df["false_duplicate_rate"],  errors="coerce")
    df["possible_allele_dropout"] = pd.to_numeric(df["possible_allele_dropout"], errors="coerce").fillna(0)
    return df


def make_figure(df: pd.DataFrame, outdir: Path) -> None:
    col_labels = [f"{MOD_LABELS.get(m, m)}\n{g}"
                  for m in MODALITY_ORDER for g in GENE_ORDER]
    col_keys   = [(m, g) for m in MODALITY_ORDER for g in GENE_ORDER]
    tools_present = [t for t in TOOL_ORDER if t in df["tool"].values]

    # Heatmap matrix
    heatmap = np.full((len(tools_present), len(col_keys)), np.nan)
    for i, tool in enumerate(tools_present):
        for j, (mod, gene) in enumerate(col_keys):
            row = df[(df["tool"] == tool) & (df["modality"] == mod) & (df["gene"] == gene)]
            if not row.empty:
                heatmap[i, j] = row["false_duplicate_rate"].iloc[0]

    # Stacked bar counts
    bar_data = {}
    for tool in tools_present:
        counts = []
        for mod in MODALITY_ORDER:
            sub = df[(df["tool"] == tool) & (df["modality"] == mod)]
            counts.append(sub["possible_allele_dropout"].sum())
        bar_data[tool] = counts

    # Layout — extra right margin for legend outside
    fig = plt.figure(figsize=(15, 5.5))
    gs  = fig.add_gridspec(1, 2, width_ratios=[2.6, 1], wspace=0.38)
    ax_heat = fig.add_subplot(gs[0])
    ax_bar  = fig.add_subplot(gs[1])

    # Panel A — Heatmap
    cmap = plt.cm.RdYlGn_r
    cmap.set_bad(color="#e8e8e8")
    im = ax_heat.imshow(heatmap, cmap=cmap, vmin=0.0, vmax=1.0,
                        aspect="auto", interpolation="nearest")
    for x in np.arange(-0.5, len(col_keys), 1):
        ax_heat.axvline(x, color="white", lw=0.5)
    for y in np.arange(-0.5, len(tools_present), 1):
        ax_heat.axhline(y, color="white", lw=0.5)
    for sep in [2.5, 5.5]:
        ax_heat.axvline(sep, color="#555555", lw=1.2)
    for k, mod in enumerate(MODALITY_ORDER):
        centre = k * 3 + 1
        ax_heat.text(centre, -1.1, MOD_LABELS.get(mod, mod),
                     ha="center", va="bottom", fontsize=9,
                     fontweight="bold", transform=ax_heat.transData)
    for i in range(len(tools_present)):
        for j in range(len(col_keys)):
            val = heatmap[i, j]
            if np.isnan(val):
                continue
            txt   = f"{val:.2f}" if val > 0 else "0"
            color = "white" if val > 0.6 else "black"
            ax_heat.text(j, i, txt, ha="center", va="center", fontsize=6.5, color=color)
    ax_heat.set_xticks(range(len(col_keys)))
    ax_heat.set_xticklabels(col_labels, fontsize=7.5)
    ax_heat.set_yticks(range(len(tools_present)))
    ax_heat.set_yticklabels(tools_present, fontsize=8.5)
    ax_heat.set_title("A  False Duplicate Rate\n(proxy for allele dropout burden)",
                      fontsize=9.5, loc="left", fontweight="bold")
    cbar = fig.colorbar(im, ax=ax_heat, fraction=0.025, pad=0.02)
    cbar.set_label("False duplicate rate", fontsize=8)
    cbar.ax.tick_params(labelsize=7)
    cbar.set_ticks([0.0, 0.25, 0.50, 0.75, 1.0])

    # Panel B — Stacked bar + legend OUTSIDE
    x      = np.arange(len(MODALITY_ORDER))
    bottom = np.zeros(len(MODALITY_ORDER))
    bar_handles = []
    for tool in tools_present:
        vals = np.array(bar_data[tool])
        ax_bar.bar(x, vals, bottom=bottom,
                   color=TOOL_COLORS.get(tool, "#aaaaaa"),
                   label=tool, edgecolor="white", linewidth=0.4)
        bar_handles.append(mpatches.Patch(color=TOOL_COLORS.get(tool, "#aaaaaa"), label=tool))
        bottom += vals
    ax_bar.set_xticks(x)
    ax_bar.set_xticklabels([MOD_LABELS.get(m, m) for m in MODALITY_ORDER])
    ax_bar.set_ylabel("Possible allele dropout calls\n(sum across A + B + C genes)")
    ax_bar.set_title("B  Possible Allele Dropout\nby Modality",
                     fontsize=9.5, loc="left", fontweight="bold")
    # Legend outside — right of Panel B
    legend_outside(ax_bar, loc="upper left", anchor=(1.02, 1.0),
                   handles=bar_handles, fontsize=7)

    fig.suptitle(
        "HLA Allele Dropout Characterisation Across Tools and Modalities\n"
        "Benchmark QC using 1000G truth genotypes (WGS n=131, WES n=129, RNA n=107)",
        fontsize=10, y=1.03, fontweight="bold",
    )
    fig.text(
        0.01, -0.04,
        "Note: Direct allele frequency measurement is available for one VENEX sample (VX_92_2_D1), "
        "which shows balanced heterozygosity at all three HLA class I loci (no LOH detected). "
        "For 79/80 samples, input WGS data is unavailable; LOH assessment is limited to "
        "allele dropout pattern characterisation shown here.",
        fontsize=7, color="#555555", wrap=True,
    )

    rotate_xticklabels(ax_heat, angle=0, ha="center")
    tight_with_legend(fig, right=0.78, top=0.90, bottom=0.16, left=0.06)
    save_fig(fig, outdir / "figure_loh_allele_dropout_heatmap")
    print(f"  Saved figure_loh_allele_dropout_heatmap.{{pdf,svg,png}} → {outdir}")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--input",  type=Path, default=DEFAULT_INPUT)
    ap.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    args = ap.parse_args()
    if not args.input.exists():
        raise FileNotFoundError(f"Input not found: {args.input}")
    args.outdir.mkdir(parents=True, exist_ok=True)
    df = load_data(args.input)
    make_figure(df, args.outdir)


if __name__ == "__main__":
    main()
