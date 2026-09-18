#!/usr/bin/env python3.11
"""
generate_figures_pub.py — Secondary publication figures for cross-modality analyses.

  Figure 9: WES+RNA bimodal joint consensus per-gene accuracy (HLA-A/B/C)
            6 bars per gene: WES-MV, WES-WC, RNA-MV, RNA-WC, Bimodal-MV, Bimodal-WC

  Figure 10: Trimodal vs bimodal — two-panel
            Panel A: accuracy-among-callable for all 4 methods
            Panel B: callable-rate comparison showing coverage cost

Usage:
    python3 bin/generate_figures_pub.py \
        --trimodal-tables analysis/benchmark_trimodal_all_samples/tables \
        --out-dir analysis/figures_final
"""

import argparse, math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# ── Style constants ───────────────────────────────────────────────────────────
GRID_COLOR  = "#E5E7EB"
NEUTRAL_MID = "#9CA3AF"
DPI         = 300
WATERMARK   = "MVHLA · IMGT/HLA 3.59.0 · HLA-A, -B, -C"


def apply_style():
    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "font.size": 10,
        "axes.titlesize": 12, "axes.titleweight": "bold",
        "axes.labelsize": 10, "xtick.labelsize": 9, "ytick.labelsize": 9,
        "axes.spines.top": False, "axes.spines.right": False,
        "axes.grid": True, "axes.grid.axis": "y",
        "grid.color": GRID_COLOR, "grid.linewidth": 0.6,
        "legend.fontsize": 9, "legend.framealpha": 0.9,
        "figure.dpi": DPI, "savefig.dpi": DPI, "savefig.bbox": "tight",
        "pdf.fonttype": 42, "ps.fonttype": 42,
    })


def wm(fig):
    fig.text(0.99, 0.005, WATERMARK, ha="right", va="bottom",
             fontsize=7, color="#B8C0CC", transform=fig.transFigure)


def save(fig, out_dir, stem):
    out_dir.mkdir(parents=True, exist_ok=True)
    for ext in ("png", "pdf", "svg"):
        fig.savefig(str(out_dir / f"{stem}.{ext}"),
                    dpi=DPI if ext == "png" else None)
    plt.close(fig)


def wilson_ci(k, n, z=1.96):
    if n == 0:
        return 0.0, 0.0
    p = k / n
    d = 1 + z ** 2 / n
    c = (p + z ** 2 / (2 * n)) / d
    m = z * math.sqrt(p * (1 - p) / n + z ** 2 / (4 * n ** 2)) / d
    return max(0.0, c - m), min(1.0, c + m)


def _read(path):
    df = pd.read_csv(str(path), sep="\t", dtype=str)
    for col in df.columns:
        converted = pd.to_numeric(df[col], errors="coerce")
        if not converted.isna().all():
            df[col] = converted
    return df


# ── Figure 9: bimodal per-gene ────────────────────────────────────────────────

def _per_gene_acc(df, gene, modality_filter, method_filter, bimodal_samples):
    sub = df[
        (df["gene"] == gene) &
        (df["modality"] == modality_filter) &
        (df["method"] == method_filter) &
        (df["sample"].isin(bimodal_samples)) &
        (df["is_callable"].astype(float) == 1)
    ]
    if sub.empty:
        return np.nan, np.nan, np.nan
    n = len(sub)
    k = int(sub["is_correct"].astype(float).sum())
    acc = k / n
    lo, hi = wilson_ci(k, n)
    return acc, lo, hi


def _per_gene_acc_bimodal(df, gene, method):
    sub = df[
        (df["gene"] == gene) &
        (df["method"] == method) &
        (df["is_callable"].astype(float) == 1)
    ]
    if sub.empty:
        return np.nan, np.nan, np.nan
    n = len(sub)
    k = int(sub["is_correct"].astype(float).sum())
    acc = k / n
    lo, hi = wilson_ci(k, n)
    return acc, lo, hi


def fig9_bimodal_per_gene(bimodal_df, mv_df, wc_df, out_dir):
    """Figure 9: per-gene accuracy for WES-MV, WES-WC, RNA-MV, RNA-WC, Bimodal-MV, Bimodal-WC."""
    genes = ["A", "B", "C"]
    bimodal_samples = set(bimodal_df["sample"].unique())

    for df in [bimodal_df, mv_df, wc_df]:
        for col in ["is_callable", "is_correct"]:
            if col in df.columns:
                df[col] = pd.to_numeric(df[col], errors="coerce")

    comparisons = [
        ("WES_MajorityVote",          "WES Majority",        "#4C78A8"),
        ("WES_WeightedConsensus",     "WES Weighted",        "#2A6099"),
        ("RNA_MajorityVote",          "RNA Majority",        "#E69F00"),
        ("RNA_WeightedConsensus",     "RNA Weighted",        "#B87A00"),
        ("Bimodal_MajorityVote",      "WES+RNA Majority",    "#C44E52"),
        ("Bimodal_WeightedConsensus", "WES+RNA Weighted",    "#8B1A1A"),
    ]

    per_gene = {}
    for gene in genes:
        per_gene[gene] = {
            "WES_MajorityVote":          _per_gene_acc(mv_df,      gene, "wes",    "MajorityVote",          bimodal_samples),
            "WES_WeightedConsensus":     _per_gene_acc(wc_df,      gene, "wes",    "WeightedConsensus",     bimodal_samples),
            "RNA_MajorityVote":          _per_gene_acc(mv_df,      gene, "rnaseq", "MajorityVote",          bimodal_samples),
            "RNA_WeightedConsensus":     _per_gene_acc(wc_df,      gene, "rnaseq", "WeightedConsensus",     bimodal_samples),
            "Bimodal_MajorityVote":      _per_gene_acc_bimodal(bimodal_df, gene, "BimodalMajorityVote"),
            "Bimodal_WeightedConsensus": _per_gene_acc_bimodal(bimodal_df, gene, "BimodalWeightedConsensus"),
        }

    n_groups = len(comparisons)
    width    = 0.12
    x        = np.arange(len(genes))
    offsets  = np.linspace(-(n_groups - 1) * width / 2,
                            (n_groups - 1) * width / 2, n_groups)

    fig, ax = plt.subplots(figsize=(12, 6))

    for i, (comp_key, label, color) in enumerate(comparisons):
        vals, yerr_lo, yerr_hi = [], [], []
        for gene in genes:
            acc, ci_lo, ci_hi = per_gene[gene].get(comp_key, (np.nan, np.nan, np.nan))
            vals.append(acc)
            yerr_lo.append(0 if np.isnan(acc) else acc - ci_lo)
            yerr_hi.append(0 if np.isnan(acc) else ci_hi - acc)

        positions  = x + offsets[i]
        is_bimodal = "Bimodal" in comp_key
        ax.bar(positions, np.nan_to_num(vals, nan=0.0),
               width=width * 0.9,
               color=color,
               label=label,
               edgecolor="#111827" if is_bimodal else "white",
               linewidth=1.6 if is_bimodal else 0.4,
               zorder=3 if is_bimodal else 2)

        for pos, acc, lo, hi in zip(positions, vals, yerr_lo, yerr_hi):
            if np.isnan(acc):
                continue
            ax.errorbar(pos, acc, yerr=[[lo], [hi]], fmt="none",
                        ecolor="#374151", elinewidth=1.0, capsize=2.5, capthick=1.0,
                        zorder=4)

    for gap in [0.5, 1.5]:
        ax.axvline(gap, color=GRID_COLOR, lw=1.2, zorder=1)

    ax.set_xticks(x)
    ax.set_xticklabels([f"HLA-{g}" for g in genes], fontsize=11, fontweight="bold")
    ax.set_ylabel("Correct-call rate (callable loci)", fontsize=10)
    ax.set_ylim(0, 1.08)
    ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda v, _: f"{v:.0%}"))

    n_samples = len(bimodal_df["sample"].unique())
    ax.text(0.99, 0.02,
            f"n = {n_samples} samples (WES + RNA-seq); error bars = 95% Wilson CI",
            transform=ax.transAxes, ha="right", va="bottom",
            fontsize=8, color=NEUTRAL_MID)

    ax.legend(ncol=3, fontsize=8, loc="upper left", framealpha=0.9,
              title="Method (bold outline = bimodal joint)", title_fontsize=8)
    ax.set_title(
        "Figure 9. Secondary cross-modality analysis: WES + RNA-seq bimodal HLA consensus\n"
        "(HLA-A, -B, -C; same n per gene; bold outlines = bimodal joint)",
        fontsize=12, fontweight="bold", pad=10)

    wm(fig)
    save(fig, out_dir, "figure_09_bimodal_per_gene")
    print("  Figure 9 saved.")


# ── Figure 10: trimodal vs bimodal ───────────────────────────────────────────

def _overall_stats(df, method):
    """Return (accuracy_among_callable, callable_rate, n_total, acc_lo, acc_hi)."""
    sub = df[df["method"] == method].copy()
    sub["is_callable"] = pd.to_numeric(sub["is_callable"], errors="coerce")
    sub["is_correct"]  = pd.to_numeric(sub["is_correct"],  errors="coerce")
    n_total    = len(sub)
    n_callable = int(sub["is_callable"].sum())
    callable_sub = sub[sub["is_callable"] == 1]
    if callable_sub.empty or n_total == 0:
        return np.nan, np.nan, n_total, np.nan, np.nan
    k   = int(callable_sub["is_correct"].sum())
    acc = k / n_callable
    lo, hi = wilson_ci(k, n_callable)
    return acc, n_callable / n_total, n_total, lo, hi


def fig10_trimodal_comparison(bimodal_df, trimodal_df, out_dir):
    """Figure 10: two-panel — accuracy gain and callable-rate cost of adding WGS."""

    method_specs = [
        ("BimodalMajorityVote",       "Bimodal\nMajority Vote",       "#4C78A8", "solid"),
        ("BimodalWeightedConsensus",  "Bimodal\nWeighted Consensus",  "#2A6099", "solid"),
        ("TrimodalMajorityVote",      "Trimodal\nMajority Vote",      "#C44E52", "dashed"),
        ("TrimodalWeightedConsensus", "Trimodal\nWeighted Consensus", "#8B1A1A", "dashed"),
    ]

    # Restrict bimodal to class I only (genes A/B/C) to match trimodal scope
    bimodal_abc = bimodal_df[bimodal_df["gene"].isin(["A", "B", "C"])].copy()

    stats = {}
    for meth_key, _, _, _ in method_specs:
        src = bimodal_abc if meth_key.startswith("Bimodal") else trimodal_df
        stats[meth_key] = _overall_stats(src, meth_key)

    fig, (ax_acc, ax_call) = plt.subplots(1, 2, figsize=(12, 6))
    fig.suptitle(
        "Figure 10. Secondary cross-modality analysis: trimodal (WGS+WES+RNA) vs. bimodal (WES+RNA)\n"
        "(HLA-A, -B, -C; n=370 trimodal assessments, n≤391 bimodal)",
        fontsize=12, fontweight="bold", y=1.01)

    x     = np.arange(len(method_specs))
    width = 0.55

    # Panel A: accuracy among callable
    accs      = []
    errs_lo   = []
    errs_hi   = []
    colors    = []
    patterns  = []
    labels    = []
    for meth_key, label, color, ls in method_specs:
        acc, cr, nt, lo, hi = stats[meth_key]
        accs.append(acc if not np.isnan(acc) else 0.0)
        errs_lo.append(0.0 if np.isnan(acc) else acc - lo)
        errs_hi.append(0.0 if np.isnan(acc) else hi - acc)
        colors.append(color)
        patterns.append("///" if ls == "dashed" else "")
        labels.append(label)

    bars_a = ax_acc.bar(x, accs, width=width, color=colors,
                        hatch=patterns, edgecolor="white", linewidth=0.5)
    ax_acc.errorbar(x, accs, yerr=[errs_lo, errs_hi], fmt="none",
                    ecolor="#374151", elinewidth=1.2, capsize=4, capthick=1.2, zorder=5)

    # Annotate bars with exact values
    for xi, (acc, lo_e, hi_e) in enumerate(zip(accs, errs_lo, errs_hi)):
        if acc > 0:
            ax_acc.text(xi, acc + hi_e + 0.003, f"{acc:.1%}",
                        ha="center", va="bottom", fontsize=8, fontweight="bold")

    acc_vals = [v for v in accs if v > 0]
    y_min = min(acc_vals) - 0.03 if acc_vals else 0.90
    ax_acc.set_ylim(max(0.88, y_min), 1.00)
    ax_acc.set_xticks(x)
    ax_acc.set_xticklabels(labels, fontsize=9)
    ax_acc.set_ylabel("Accuracy among callable loci", fontsize=10)
    ax_acc.yaxis.set_major_formatter(plt.FuncFormatter(lambda v, _: f"{v:.1%}"))
    ax_acc.set_title("Panel A: Accuracy (callable loci)", fontsize=11, fontweight="bold")

    # Add hatching legend
    from matplotlib.patches import Patch
    legend_handles = [
        Patch(facecolor="#9CA3AF", label="Bimodal (WES+RNA)"),
        Patch(facecolor="#9CA3AF", hatch="///", label="Trimodal (WGS+WES+RNA)"),
    ]
    ax_acc.legend(handles=legend_handles, fontsize=8, loc="lower right", framealpha=0.9)

    # Panel B: callable rate
    call_rates = []
    n_totals   = []
    for meth_key, _, _, _ in method_specs:
        acc, cr, nt, lo, hi = stats[meth_key]
        call_rates.append(cr if not np.isnan(cr) else 0.0)
        n_totals.append(nt)

    bars_b = ax_call.bar(x, call_rates, width=width, color=colors,
                         hatch=patterns, edgecolor="white", linewidth=0.5)
    for xi, (cr, nt) in enumerate(zip(call_rates, n_totals)):
        if cr > 0:
            ax_call.text(xi, cr + 0.002, f"{cr:.1%}\n(n={nt})",
                         ha="center", va="bottom", fontsize=8, fontweight="bold")

    cr_vals = [v for v in call_rates if v > 0]
    y_min_c = min(cr_vals) - 0.06 if cr_vals else 0.88
    ax_call.set_ylim(max(0.88, y_min_c), 1.02)
    ax_call.set_xticks(x)
    ax_call.set_xticklabels(labels, fontsize=9)
    ax_call.set_ylabel("Callable rate (fraction of loci with a call)", fontsize=10)
    ax_call.yaxis.set_major_formatter(plt.FuncFormatter(lambda v, _: f"{v:.1%}"))
    ax_call.set_title("Panel B: Callable rate (coverage cost)", fontsize=11, fontweight="bold")
    ax_call.legend(handles=legend_handles, fontsize=8, loc="lower right", framealpha=0.9)

    fig.tight_layout()
    wm(fig)
    save(fig, out_dir, "figure_10_trimodal_comparison")
    print("  Figure 10 saved.")


# ── CLI ───────────────────────────────────────────────────────────────────────

def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--trimodal-tables", required=True,
                   help="Path to the trimodal benchmark tables used for secondary bimodal/trimodal publication figures.")
    p.add_argument("--out-dir", required=True,
                   help="Output directory for figures")
    return p.parse_args()


def main():
    args   = parse_args()
    td     = Path(args.trimodal_tables)
    out    = Path(args.out_dir)
    apply_style()
    print(f"Reading tables from: {td}")
    print(f"Writing figures to:  {out}\n")

    bimodal_df   = _read(td / "bimodal_wes_rna_consensus.tsv")
    trimodal_df  = _read(td / "trimodal_wgs_wes_rna_consensus.tsv")
    mv_df        = _read(td / "majority_vote_baseline.tsv")
    wc_df        = _read(td / "weighted_consensus_calls.tsv")

    print("Figure 9 — bimodal per-gene accuracy ...")
    fig9_bimodal_per_gene(bimodal_df, mv_df, wc_df, out)

    print("Figure 10 — trimodal vs bimodal comparison ...")
    fig10_trimodal_comparison(bimodal_df, trimodal_df, out)

    print(f"\nDone. Figures written to {out}/")


if __name__ == "__main__":
    main()
