#!/usr/bin/env python3.11
"""
fig_fimm_concordance.py  [figure hub copy]

Three publication figures from fimm_hla_calls.tsv:
  1. Per-gene concordance heatmap (tools × genes)
  2. Cross-modality comparison plot
  3. scRNA concordance heatmap

Original: bin/analyze_fimm_hla.py
Hub standards: Arial font, legends outside axes, 300 DPI.
Supports de-identified input (fimm_hla_calls_deident.tsv).

Usage:
    python3.11 figure_hub/scripts/fig_fimm_concordance.py [--calls PATH] [--outdir PATH]
"""

import argparse
import math
import sys
import warnings
warnings.filterwarnings("ignore")
from pathlib import Path

import numpy as np
import pandas as pd

# ── Hub style ──────────────────────────────────────────────────────────────
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from plot_style import apply_style, rotate_xticklabels, save_fig, tight_with_legend, WONG

apply_style()

import matplotlib.pyplot as plt          # after apply_style
import matplotlib.patches as mpatches
import seaborn as sns

# ── Paths ──────────────────────────────────────────────────────────────────
_LOCAL       = Path("/scratch/project_2008084/pihla_local")
DEFAULT_CALLS_DEIDENT = _LOCAL / "fimm_results" / "fimm_hla_calls_deident.tsv"
DEFAULT_CALLS = _LOCAL / "fimm_results" / "fimm_hla_calls.tsv"
DEFAULT_OUT   = _LOCAL / "analysis" / "figures_final"

GENES = ["A", "B", "C"]

# Modality column prefixes present in fimm_hla_calls.tsv
MODALITIES = {
    "scrna":   ["scrna_arcashla", "scrna_optitype"],
    "bulkrna": ["bulkrna_OptiType", "bulkrna_arcasHLA", "bulkrna_SpecHLA"],
    "wes":     ["wes_OptiType", "wes_arcasHLA", "wes_SpecHLA"],
}


def _allele_cols(df: pd.DataFrame, gene: str) -> list[str]:
    """Return all column names that contain calls for a given gene."""
    return [c for c in df.columns if c.endswith(f"_{gene}")]


def _pairwise_concordance(df: pd.DataFrame, col_a: str, col_b: str) -> float:
    """Fraction of rows where col_a == col_b (ignoring NaN rows)."""
    if col_a == col_b:
        valid = df[[col_a]].dropna()
        return float("nan") if valid.empty else 1.0
    valid = df[[col_a, col_b]].dropna()
    if valid.empty:
        return float("nan")
    return (valid[col_a] == valid[col_b]).mean()


def make_concordance_heatmap(df: pd.DataFrame, outdir: Path) -> None:
    """Figure 1: tool-pair concordance heatmap per gene."""
    tool_cols = [c for c in df.columns if c != "sample_id"]

    for gene in GENES:
        gene_cols = _allele_cols(df, gene)
        if len(gene_cols) < 2:
            continue

        n  = len(gene_cols)
        mat = np.full((n, n), np.nan)
        for i, ca in enumerate(gene_cols):
            for j, cb in enumerate(gene_cols):
                if i <= j:
                    v = _pairwise_concordance(df, ca, cb)
                    mat[i, j] = mat[j, i] = v

        fig, ax = plt.subplots(figsize=(max(6, n * 0.9), max(5, n * 0.85)))
        short_labels = [c.replace(f"_{gene}", "").replace("_", "\n") for c in gene_cols]
        sns.heatmap(mat, annot=True, fmt=".2f", cmap="YlOrRd",
                    xticklabels=short_labels, yticklabels=short_labels,
                    vmin=0, vmax=1, ax=ax, cbar_kws={"shrink": 0.7})
        ax.set_title(f"HLA-{gene}  Tool-Pair Concordance (FIMM cohort)",
                     fontweight="bold", pad=12)
        # Rotate x labels to avoid overlap
        rotate_xticklabels(ax, angle=45)
        ax.set_yticklabels(ax.get_yticklabels(), rotation=0, fontsize=7)

        save_fig(fig, outdir / f"fig_concordance_by_gene_{gene}")
        print(f"  Saved fig_concordance_by_gene_{gene}.{{pdf,svg,png}}")


def make_crossmodal_plot(df: pd.DataFrame, outdir: Path) -> None:
    """Figure 2: cross-modality concordance per gene."""
    mod_pairs = [
        ("scrna", "wes", "scRNA vs WES"),
        ("scrna", "bulkrna", "scRNA vs BulkRNA"),
        ("bulkrna", "wes", "BulkRNA vs WES"),
    ]
    fig, axes = plt.subplots(1, len(GENES), figsize=(12, 4.5), sharey=True)

    for gi, gene in enumerate(GENES):
        ax = axes[gi]
        xs, ys, labels, colors_list = [], [], [], []
        cmap_list = [WONG["blue"], WONG["orange"], WONG["green"]]

        for pi, (ma, mb, label) in enumerate(mod_pairs):
            # pick first available column for each modality × gene
            col_a = next((c for c in df.columns if c.startswith(ma) and c.endswith(f"_{gene}")), None)
            col_b = next((c for c in df.columns if c.startswith(mb) and c.endswith(f"_{gene}")), None)
            if col_a and col_b:
                conc = _pairwise_concordance(df, col_a, col_b)
                xs.append(pi)
                ys.append(conc)
                labels.append(label)
                colors_list.append(cmap_list[pi])

        ax.bar(range(len(ys)), ys, color=colors_list, edgecolor="white")
        ax.set_xticks(range(len(labels)))
        ax.set_xticklabels(labels, rotation=30, ha="right", fontsize=7)
        ax.set_ylim(0, 1.05)
        ax.set_title(f"HLA-{gene}", fontweight="bold")
        if gi == 0:
            ax.set_ylabel("Concordance rate")

    fig.suptitle("FIMM — Cross-Modality HLA Concordance",
                 fontsize=11, fontweight="bold", y=1.02)
    tight_with_legend(fig, right=0.98)
    save_fig(fig, outdir / "fig_crossmodal_comparison")
    print("  Saved fig_crossmodal_comparison.{pdf,svg,png}")


def make_scrna_heatmap(df: pd.DataFrame, outdir: Path) -> None:
    """Figure 3: per-sample scRNA concordance with WES/BulkRNA."""
    scrna_cols = [c for c in df.columns if c.startswith("scrna_")]
    wes_cols   = [c for c in df.columns if c.startswith("wes_")]
    if not scrna_cols or not wes_cols:
        print("  SKIP fig_scrna_concordance_heatmap — no scRNA or WES columns found")
        return

    rows = []
    for _, row in df.iterrows():
        for gene in GENES:
            sc = next((c for c in scrna_cols if c.endswith(f"_{gene}")), None)
            wc = next((c for c in wes_cols   if c.endswith(f"_{gene}")), None)
            if sc and wc and pd.notna(row.get(sc)) and pd.notna(row.get(wc)):
                rows.append({
                    "sample": row["sample_id"],
                    "gene": gene,
                    "match": int(row[sc] == row[wc]),
                })

    if not rows:
        print("  SKIP fig_scrna_concordance_heatmap — no matched rows")
        return

    mat_df = pd.DataFrame(rows).pivot_table(index="sample", columns="gene",
                                             values="match", aggfunc="mean")
    fig, ax = plt.subplots(figsize=(5, max(4, len(mat_df) * 0.35)))
    sns.heatmap(mat_df, annot=True, fmt=".0f", cmap="RdYlGn",
                vmin=0, vmax=1, ax=ax, cbar_kws={"shrink": 0.6})
    ax.set_title("scRNA vs WES Concordance per Sample & Gene",
                 fontweight="bold", pad=10)
    ax.set_xlabel("HLA Gene")
    ax.set_ylabel("")
    ax.set_yticklabels(ax.get_yticklabels(), rotation=0, fontsize=7)

    save_fig(fig, outdir / "fig_scrna_concordance_heatmap")
    print("  Saved fig_scrna_concordance_heatmap.{pdf,svg,png}")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--calls",  type=Path, default=DEFAULT_CALLS_DEIDENT if DEFAULT_CALLS_DEIDENT.exists() else DEFAULT_CALLS,
                    help="fimm_hla_calls.tsv or _deident.tsv")
    ap.add_argument("--outdir", type=Path, default=DEFAULT_OUT)
    args = ap.parse_args()

    if not args.calls.exists():
        raise FileNotFoundError(f"Input not found: {args.calls}")
    args.outdir.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(args.calls, sep="\t", dtype=str)
    print(f"Loaded {len(df)} rows, {len(df.columns)} columns from {args.calls.name}")

    make_concordance_heatmap(df, args.outdir)
    make_crossmodal_plot(df, args.outdir)
    make_scrna_heatmap(df, args.outdir)
    print(f"\nAll figures saved to {args.outdir}")


if __name__ == "__main__":
    main()
