#!/usr/bin/env python3.11
"""
fix_figure1.py — Regenerate figure_01_accuracy_overview.png as a clean,
                  readable Cleveland dot plot.

Bug in original: individual-tool stems drew gap (acc→1.0) while ensemble
stems drew fill (0→acc), making the panels unreadable.

Fix: all entries use the same idiom — stem from 0 to dot, large dot at
accuracy, thin CI span through the dot, value label at right.

Usage:
    PYTHONPATH=/users/ozcanumu/.local/lib/python3.11/site-packages \
    python3.11 bin/fix_figure1.py \
        --tables-dir analysis/benchmark_trimodal_all_samples/tables \
        --out-dir    analysis/report_v10_figures
"""
import argparse, math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import numpy as np
import pandas as pd

# ── palette ─────────────────────────────────────────────────────────────────
TOOL_COLORS = {
    "ArcasHLA":           "#7A7A7A",
    "HLA-HD":             "#4C78A8",
    "Kourami":            "#8F63B8",
    "OptiType":           "#2A9D8F",
    "POLYSOLVER":         "#E69F00",
    "Seq2HLA":            "#A6A57A",
    "SpecHLA":            "#E07B54",
    "T1K":                "#56B4E9",
    "MajorityVote":       "#C44E52",
    "WeightedConsensus":  "#7B3F9E",
}
ENSEMBLE_COLORS = {
    "MajorityVote":      "#C44E52",
    "WeightedConsensus": "#7B3F9E",
}
GRID_COLOR = "#E5E7EB"
DPI = 300
WATERMARK = "PIHLA · IMGT/HLA 3.59.0 · 2026-05-11"

MOD_LABELS = {
    "wgs":    "WGS\n(n≤138 samples)",
    "wes":    "WES\n(n≤130 samples)",
    "rnaseq": "RNA-seq\n(n≤107 samples)",
}
MOD_COLORS = {
    "wgs":    "#2A9D8F",
    "wes":    "#E69F00",
    "rnaseq": "#C44E52",
}
MODALITY_ORDER = ["wgs", "wes", "rnaseq"]

# Off-label tools get a grey dot outline to mark "not designed for this input"
OFF_LABEL = {
    "wgs":    {"ArcasHLA", "Seq2HLA"},   # RNA tools run on WGS
    "wes":    {"ArcasHLA", "Seq2HLA"},   # RNA tools run on WES
    "rnaseq": {"SpecHLA", "Kourami", "POLYSOLVER"},  # DNA-only on RNA
}


def wilson_ci(k, n, z=1.96):
    if n == 0:
        return 0.0, 0.0
    p = k / n
    d = 1 + z**2 / n
    c = (p + z**2 / (2*n)) / d
    m = z * math.sqrt(p*(1-p)/n + z**2/(4*n**2)) / d
    return max(0.0, c - m), min(1.0, c + m)


def load(tables_dir):
    def read(name):
        p = tables_dir / name
        if not p.exists():
            return pd.DataFrame()
        df = pd.read_csv(str(p), sep="\t", dtype=str)
        for col in df.columns:
            try:
                df[col] = pd.to_numeric(df[col])
            except (ValueError, TypeError):
                pass
        return df
    return read("summary_full_cohort.tsv"), read("method_comparison.tsv")


def make_figure(summ, meth, out_dir):
    fig, axes = plt.subplots(1, 3, figsize=(17, 7), sharey=False)
    plt.rcParams.update({
        "font.family": "DejaVu Sans", "font.size": 10,
        "axes.spines.top": False, "axes.spines.right": False,
        "axes.spines.left": False,
        "figure.dpi": DPI,
    })

    for ax, mod in zip(axes, MODALITY_ORDER):
        sub  = summ[summ["modality"] == mod].copy() if not summ.empty else pd.DataFrame()
        emeth = meth[meth["modality"] == mod].copy() if not meth.empty else pd.DataFrame()

        # ── Collect rows ────────────────────────────────────────────────────
        rows = []   # (label, acc, ci_lo, ci_hi, color, is_ensemble, off_label)

        # Individual tools, sorted descending by accuracy
        for _, r in sub.sort_values("overall_correct_call_rate", ascending=False).iterrows():
            tool = str(r["tool"])
            acc  = float(r["overall_correct_call_rate"])
            ci_lo = float(r.get("overall_correct_call_rate_ci_lo", acc))
            ci_hi = float(r.get("overall_correct_call_rate_ci_hi", acc))
            color = TOOL_COLORS.get(tool, "#9CA3AF")
            ol    = tool in OFF_LABEL.get(mod, set())
            rows.append((tool, acc, ci_lo, ci_hi, color, False, ol))

        # Separator
        rows.append(("", None, None, None, "white", False, False))

        # Ensemble methods: MajorityVote then WeightedConsensus
        for ens in ("MajorityVote", "WeightedConsensus"):
            r = emeth[emeth["method"] == ens]
            if r.empty:
                continue
            acc   = float(r["overall_correct_call_rate"].iloc[0])
            ci_lo = float(r["overall_correct_call_rate_ci_lo"].iloc[0])
            ci_hi = float(r["overall_correct_call_rate_ci_hi"].iloc[0])
            color = ENSEMBLE_COLORS[ens]
            rows.append((ens, acc, ci_lo, ci_hi, color, True, False))

        # ── Plot ────────────────────────────────────────────────────────────
        n = len(rows)
        y_positions = np.arange(n)

        ax.set_xlim(-0.02, 1.15)
        ax.set_ylim(-0.8, n - 0.2)
        ax.set_yticks([])
        ax.xaxis.set_major_formatter(plt.FuncFormatter(lambda v, _: f"{v*100:.0f}%"))
        ax.set_xticks([0, 0.25, 0.5, 0.75, 1.0])
        ax.tick_params(axis="x", labelsize=8.5)
        ax.axvline(1.0, color="#CBD5E1", lw=1, ls="--", zorder=0)
        ax.axvline(0,   color="#CBD5E1", lw=0.6, zorder=0)
        ax.set_facecolor("#FAFBFC")
        ax.grid(axis="x", color=GRID_COLOR, linewidth=0.5, zorder=0)

        # Alternating row shading
        for yi in range(n):
            if yi % 2 == 0:
                ax.axhspan(yi - 0.48, yi + 0.48, color="#F1F5F9", zorder=0, alpha=0.6)

        for yi, (label, acc, ci_lo, ci_hi, color, is_ens, off_lbl) in \
                enumerate(rows):
            if acc is None:
                continue

            # ── Stem from 0 to ci_lo ────────────────────────────────────
            ax.plot([0, ci_lo], [yi, yi],
                    color="#CBD5E1", lw=1.6, solid_capstyle="round", zorder=1)

            # ── CI span (ci_lo to ci_hi) ────────────────────────────────
            ci_lw = 2.5 if is_ens else 1.8
            ax.plot([ci_lo, ci_hi], [yi, yi],
                    color=color, lw=ci_lw, alpha=0.55, zorder=2)

            # ── CI tick marks ────────────────────────────────────────────
            tick_h = 0.22
            for xv in (ci_lo, ci_hi):
                ax.plot([xv, xv], [yi - tick_h, yi + tick_h],
                        color=color, lw=1.4, zorder=3)

            # ── Main dot ─────────────────────────────────────────────────
            dot_size  = 120 if is_ens else 80
            edge_color = "#555" if off_lbl else color
            edge_lw    = 1.8  if off_lbl else 0.8
            linestyle  = "--" if off_lbl else "-"
            ax.scatter([acc], [yi],
                       s=dot_size, color=color,
                       edgecolors=edge_color, linewidths=edge_lw,
                       zorder=4)
            if off_lbl:
                # Dashed ring to mark off-label
                ax.scatter([acc], [yi], s=dot_size * 2.4,
                           facecolors="none", edgecolors=edge_color,
                           linewidths=1.0, linestyles="--", zorder=3)

            # ── Tool label (left) ─────────────────────────────────────
            fw = "bold" if is_ens else "normal"
            ax.text(-0.02, yi, label, ha="right", va="center",
                    fontsize=8.5 if is_ens else 8,
                    fontweight=fw, color=color,
                    transform=ax.get_yaxis_transform())

            # ── Value label (right of dot) ────────────────────────────
            lbl_x = min(ci_hi + 0.012, 1.13)
            ax.text(lbl_x, yi, f"{acc*100:.1f}%",
                    ha="left", va="center",
                    fontsize=8 if is_ens else 7.5,
                    fontweight=fw, color=color, zorder=5)

        # Title
        mod_col = MOD_COLORS[mod]
        ax.set_title(MOD_LABELS[mod], fontsize=11, fontweight="bold",
                     color=mod_col, pad=10)
        ax.set_xlabel("Overall Correct-Call Rate", fontsize=9, labelpad=6)

    # ── Legend ──────────────────────────────────────────────────────────────
    legend_handles = [
        mlines.Line2D([], [], marker="o", color="w",
                      markerfacecolor=ENSEMBLE_COLORS["MajorityVote"],
                      markersize=9, label="MajorityVote (ensemble)"),
        mlines.Line2D([], [], marker="o", color="w",
                      markerfacecolor=ENSEMBLE_COLORS["WeightedConsensus"],
                      markersize=9, label="WeightedConsensus (ensemble)"),
        mlines.Line2D([], [], marker="o", color="w",
                      markerfacecolor="#9CA3AF", markersize=8,
                      label="Individual tool"),
        mlines.Line2D([], [], marker="o", color="w",
                      markerfacecolor="#9CA3AF", markersize=8,
                      markeredgecolor="#555", markeredgewidth=1.8,
                      label="Off-label tool *"),
    ]
    fig.legend(handles=legend_handles, loc="lower center",
               ncol=4, fontsize=8.5, framealpha=0.95,
               bbox_to_anchor=(0.5, -0.02))

    fig.suptitle(
        "HLA Typing Tool Accuracy: Single Tools and Ensemble Methods\n"
        "by Sequencing Modality",
        fontsize=13, fontweight="bold", y=1.01)

    fig.text(0.5, -0.07,
             "* Off-label: tool evaluated on a sequencing modality it was not designed for "
             "(dashed ring). Dots = overall correct-call rate; horizontal bars = 95% Wilson CI.",
             ha="center", fontsize=7.5, color="#64748B")

    fig.text(0.99, 0.0, WATERMARK,
             ha="right", va="bottom", fontsize=7, color="#B8C0CC")

    fig.tight_layout(rect=[0.09, 0.04, 1.0, 1.0])
    out_dir.mkdir(parents=True, exist_ok=True)
    for ext in ("png", "pdf"):
        fig.savefig(str(out_dir / f"figure_01_accuracy_overview.{ext}"),
                    dpi=DPI if ext == "png" else None, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved figure_01_accuracy_overview.png/pdf → {out_dir}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tables-dir", required=True, type=Path)
    ap.add_argument("--out-dir",    required=True, type=Path)
    args = ap.parse_args()

    summ, meth = load(args.tables_dir)
    make_figure(summ, meth, args.out_dir)


if __name__ == "__main__":
    main()
