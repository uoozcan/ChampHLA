#!/usr/bin/env python3.11
"""
generate_figures_v3.py — Streamlined publication figures focused on majority voting.

Five focused figures that tell a clear story:
  1. Accuracy overview: all tools + MajorityVote per modality (main comparison)
  2. MajorityVote advantage: gain over individual tools and average
  3. Tool agreement → accuracy by HLA-A/B/C gene: why majority voting works
  4. Per-gene accuracy: HLA-A, B, C breakdown with MajorityVote
  5. Cross-modality landscape: which tools work where

Usage:
    python3.11 bin/generate_figures_v3.py \
        --tables-dir analysis/benchmark_trimodal_all_samples/tables \
        --out-dir    analysis/figures_v3
"""

import argparse, csv, math, re, glob as globmod
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np
import pandas as pd

# ── Style constants ───────────────────────────────────────────────────────────
TOOL_COLORS = {
    # Muted, colorblind-safe Okabe-Ito-inspired palette.
    "ArcasHLA":          "#7A7A7A",
    "HLA-HD":            "#4C78A8",
    "Kourami":           "#8F63B8",
    "OptiType":          "#2A9D8F",
    "POLYSOLVER":        "#E69F00",
    "Seq2HLA":           "#A6A57A",
    "SpecHLA":           "#6C7A89",
    "T1K":               "#56B4E9",
    "MajorityVote":      "#C44E52",
}
MODALITY_ORDER  = ["wgs", "wes", "rnaseq"]
MODALITY_LABELS = {"wgs": "WGS\n(n=131)", "wes": "WES\n(n=75)", "rnaseq": "RNA-seq\n(n=92)"}
GENE_ORDER = ["A", "B", "C"]
GENE_COLORS = {"A": "#4C78A8", "B": "#E69F00", "C": "#2A9D8F"}
POSITIVE_COLOR = "#4C9F70"
NEGATIVE_COLOR = "#B279A2"
NEUTRAL_DARK = "#4B5563"
NEUTRAL_MID = "#9CA3AF"
GRID_COLOR = "#E5E7EB"
MISSING_CELL = "#ECEFF3"
DNA_MARKER = "#AFC6E9"
RNA_MARKER = "#E8C2D4"
HEATMAP_CMAP = plt.cm.cividis
WATERMARK = "MVHLA v3.0 · IMGT/HLA 3.59.0 · HLA-A, -B, -C · majority voting"
DPI = 300

TOOL_ORDER = ["OptiType","POLYSOLVER","HLA-HD","T1K","SpecHLA","Kourami","Seq2HLA","ArcasHLA"]
TOOL_DESIGN = {
    "ArcasHLA": "RNA", "HLA-HD": "Both", "Kourami": "DNA",
    "OptiType": "DNA/RNA", "POLYSOLVER": "DNA(WES)", "Seq2HLA": "RNA",
    "SpecHLA": "DNA", "T1K": "Both",
}

def apply_style():
    plt.rcParams.update({
        "font.family": "DejaVu Sans", "font.size": 10,
        "axes.titlesize": 12, "axes.titleweight": "bold",
        "axes.labelsize": 10, "xtick.labelsize": 9, "ytick.labelsize": 9,
        "axes.spines.top": False, "axes.spines.right": False,
        "axes.grid": True, "axes.grid.axis": "y",
        "grid.color": GRID_COLOR, "grid.linewidth": 0.6,
        "legend.fontsize": 9, "legend.framealpha": 0.9,
        "figure.dpi": DPI, "savefig.dpi": DPI, "savefig.bbox": "tight",
    })

def wm(fig):
    fig.text(0.99, 0.005, WATERMARK, ha="right", va="bottom",
             fontsize=7, color="#B8C0CC", transform=fig.transFigure)

def save(fig, out_dir, stem):
    out_dir.mkdir(parents=True, exist_ok=True)
    for ext in ("png", "pdf"):
        fig.savefig(str(out_dir / f"{stem}.{ext}"), dpi=DPI if ext == "png" else None)
    plt.close(fig)

def wilson_ci(k, n, z=1.96):
    if n == 0: return 0.0, 0.0
    p = k / n
    d = 1 + z**2 / n
    c = (p + z**2 / (2*n)) / d
    m = z * math.sqrt(p*(1-p)/n + z**2/(4*n**2)) / d
    return max(0, c-m), min(1, c+m)

def load_tables(tables_dir):
    files = {
        "summary":    "summary_full_cohort.tsv",
        "method":     "method_comparison.tsv",
        "per_gene":   "method_per_gene.tsv",
        "mv_calls":   "majority_vote_baseline.tsv",
        "harmonized": "harmonized_benchmark_rows.tsv",
        "bimodal":         "bimodal_accuracy_comparison.tsv",
        "bimodal_detail":  "bimodal_wes_rna_consensus.tsv",
        "wc_calls":        "weighted_consensus_calls.tsv",
    }
    tables = {}
    for key, fname in files.items():
        p = tables_dir / fname
        if not p.exists():
            tables[key] = pd.DataFrame()
            continue
        df = pd.read_csv(str(p), sep="\t", dtype=str)
        for col in df.columns:
            try: df[col] = pd.to_numeric(df[col], errors="ignore")
            except: pass
        tables[key] = df
    return tables

# ── Figure 1: Main accuracy overview ─────────────────────────────────────────
def fig1_accuracy_overview(tables, out_dir):
    """All tools + MajorityVote per modality — the headline figure."""
    summ  = tables["summary"].copy()
    meth  = tables["method"].copy()
    if summ.empty: return

    fig, axes = plt.subplots(1, 3, figsize=(16, 6), sharey=False)

    for ax, mod in zip(axes, MODALITY_ORDER):
        sub = summ[summ["modality"] == mod].copy()
        mv  = meth[(meth["method"] == "MajorityVote") & (meth["modality"] == mod)]

        # Build rows: (tool, acc, ci_lo, ci_hi, color, n_loci)
        rows = []
        for tool in TOOL_ORDER:
            r = sub[sub["tool"] == tool]
            if r.empty: continue
            acc = float(r["overall_correct_call_rate"].iloc[0])
            n   = int(r["gene_rows"].iloc[0]) if "gene_rows" in r.columns else 0
            ci_lo, ci_hi = wilson_ci(round(acc*n), n)
            rows.append((tool, acc, ci_lo, ci_hi, TOOL_COLORS.get(tool, NEUTRAL_MID)))

        # Add MajorityVote at bottom after gap
        if not mv.empty:
            acc_mv = float(mv["overall_correct_call_rate"].iloc[0])
            n_mv   = int(mv["gene_rows"].iloc[0]) if "gene_rows" in mv.columns else 0
            ci_lo_mv, ci_hi_mv = wilson_ci(round(acc_mv*n_mv), n_mv)
            rows.append(("_sep", 0, 0, 0, "white"))
            rows.append(("MajorityVote", acc_mv, ci_lo_mv, ci_hi_mv, TOOL_COLORS["MajorityVote"]))

        tools   = [r[0] for r in rows]
        accs    = [r[1] for r in rows]
        ci_los  = [r[2] for r in rows]
        ci_his  = [r[3] for r in rows]
        colors  = [r[4] for r in rows]

        x = np.arange(len(tools))
        bars = ax.bar(x, accs, color=colors, edgecolor="white", linewidth=0.4, width=0.75)

        # CI error bars (skip separator)
        for xi, (lo, hi, acc) in enumerate(zip(ci_los, ci_his, accs)):
            if acc == 0: continue
            ax.plot([xi, xi], [lo, hi], color=NEUTRAL_DARK, lw=1.4)
            ax.plot([xi-0.1, xi+0.1], [lo, lo], color=NEUTRAL_DARK, lw=1.4)
            ax.plot([xi-0.1, xi+0.1], [hi, hi], color=NEUTRAL_DARK, lw=1.4)

        # MajorityVote reference line across full width
        if not mv.empty:
            ax.axhline(acc_mv, color=TOOL_COLORS["MajorityVote"],
                       lw=1.8, ls="--", alpha=0.7, zorder=0)

        # Accuracy labels on bars
        for xi, (tool, acc) in enumerate(zip(tools, accs)):
            if acc == 0 or tool == "_sep": continue
            ax.text(xi, acc + 0.01, f"{acc:.3f}", ha="center", va="bottom",
                    fontsize=7.5,
                    fontweight="bold" if tool == "MajorityVote" else "normal",
                    color=TOOL_COLORS["MajorityVote"] if tool == "MajorityVote" else NEUTRAL_DARK)

        labels = [t if t != "_sep" else "" for t in tools]
        ax.set_xticks(x)
        ax.set_xticklabels(labels, rotation=35, ha="right", fontsize=8.5)
        for tick, lbl in zip(ax.get_xticklabels(), labels):
            tick.set_fontweight("bold" if lbl == "MajorityVote" else "normal")
            if lbl == "MajorityVote":
                tick.set_color(TOOL_COLORS["MajorityVote"])
        ax.set_ylim(0, 1.12)
        ax.set_ylabel("Overall correct-call rate", fontsize=9)
        ax.set_title(MODALITY_LABELS[mod], fontsize=11)

    fig.suptitle("HLA typing accuracy: individual tools vs majority voting\n"
                 "across WGS, WES, and RNA-seq modalities",
                 fontsize=13, fontweight="bold")
    fig.tight_layout()
    wm(fig)
    save(fig, out_dir, "figure_01_accuracy_overview")
    print("  Figure 1 saved.")

# ── Figure 2: MajorityVote advantage ─────────────────────────────────────────
def fig2_majority_vote_advantage(tables, out_dir):
    """Gain of MajorityVote over each tool and over best/average individual tool."""
    summ = tables["summary"].copy()
    meth = tables["method"].copy()
    if summ.empty or meth.empty: return

    fig, axes = plt.subplots(1, 3, figsize=(15, 5), sharey=False)

    for ax, mod in zip(axes, MODALITY_ORDER):
        sub = summ[summ["modality"] == mod].copy()
        sub["overall_correct_call_rate"] = pd.to_numeric(sub["overall_correct_call_rate"], errors="coerce")
        mv  = meth[(meth["method"] == "MajorityVote") & (meth["modality"] == mod)]
        if mv.empty or sub.empty: continue

        acc_mv = float(mv["overall_correct_call_rate"].iloc[0])

        gains, labels, colors = [], [], []
        for tool in TOOL_ORDER:
            r = sub[sub["tool"] == tool]
            if r.empty: continue
            acc = float(r["overall_correct_call_rate"].iloc[0])
            gain = acc_mv - acc
            gains.append(gain)
            labels.append(tool)
            colors.append(POSITIVE_COLOR if gain > 0 else NEGATIVE_COLOR)

        # Add average and best reference lines
        avg_acc  = float(sub["overall_correct_call_rate"].mean())
        best_acc = float(sub["overall_correct_call_rate"].max())
        best_tool = sub.loc[sub["overall_correct_call_rate"].idxmax(), "tool"]

        y = np.arange(len(labels))
        bars = ax.barh(y, gains, color=colors, alpha=0.85, edgecolor="white", height=0.65)

        # Zero line
        ax.axvline(0, color=NEUTRAL_DARK, lw=1.5)

        # Gain labels
        for yi, (g, lab) in enumerate(zip(gains, labels)):
            side = 0.003 if g >= 0 else -0.003
            ha = "left" if g >= 0 else "right"
            ax.text(g + side, yi, f"{g:+.3f}", va="center", ha=ha, fontsize=8)

        ax.set_yticks(y)
        ax.set_yticklabels(labels, fontsize=9)
        ax.set_xlabel("Gain (MajorityVote − individual tool)", fontsize=9)
        ax.set_title(MODALITY_LABELS[mod], fontsize=11)

        # Annotate best tool
        ax.text(0.02, -0.12, f"MajorityVote: {acc_mv:.3f}\nBest single ({best_tool}): {best_acc:.3f}",
                transform=ax.transAxes, fontsize=8, color=NEUTRAL_DARK,
                bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="#e5e7eb"))

    fig.suptitle("MajorityVote gain over individual HLA tools\n"
                 "(positive = MajorityVote better; negative = tool better)",
                 fontsize=13, fontweight="bold")
    fig.tight_layout()
    wm(fig)
    save(fig, out_dir, "figure_02_majority_vote_advantage")
    print("  Figure 2 saved.")

# ── Figure 3: Tool agreement → accuracy ──────────────────────────────────────
def fig3_agreement_accuracy(tables, out_dir):
    """When more tools agree, accuracy is higher — shown separately for HLA-A/B/C."""
    mv = tables["mv_calls"].copy()
    if mv.empty: return

    for col in ["is_correct","is_callable","agreeing_tools","contributing_tools"]:
        if col in mv.columns:
            mv[col] = pd.to_numeric(mv[col], errors="coerce")

    fig, axes = plt.subplots(3, 3, figsize=(15, 11), sharex=False, sharey=True)

    for row_i, gene in enumerate(GENE_ORDER):
        for col_i, mod in enumerate(MODALITY_ORDER):
            ax = axes[row_i, col_i]
            sub = mv[
                (mv["modality"] == mod) &
                (mv["gene"] == gene) &
                (mv["is_callable"] == 1)
            ].copy()
            if sub.empty:
                ax.set_axis_off()
                continue

            agreement_vals = sorted(sub["agreeing_tools"].dropna().unique())
            acc_by_agree = {}
            n_by_agree = {}
            for n_agree in agreement_vals:
                grp = sub[sub["agreeing_tools"] == n_agree]
                if len(grp) < 3:
                    continue
                acc_by_agree[n_agree] = grp["is_correct"].mean()
                n_by_agree[n_agree] = len(grp)

            if not acc_by_agree:
                ax.set_axis_off()
                continue

            xs = sorted(acc_by_agree.keys())
            ys = [acc_by_agree[x] for x in xs]
            ns = [n_by_agree[x] for x in xs]

            bars = ax.bar(
                xs, ys, color=GENE_COLORS[gene], alpha=0.9,
                edgecolor="white", width=0.72
            )

            for x, y, n in zip(xs, ys, ns):
                ax.text(x, y + 0.03, f"n={n}", ha="center", va="bottom",
                        fontsize=7, color=NEUTRAL_DARK)
                ax.text(x, max(y * 0.5, 0.06), f"{y:.0%}", ha="center", va="center",
                        fontsize=7.5, fontweight="bold", color="white")

            if len(xs) >= 3:
                z = np.polyfit(xs, ys, 1)
                p = np.poly1d(z)
                xfit = np.linspace(min(xs), max(xs), 50)
                ax.plot(xfit, p(xfit), "--", color=NEUTRAL_DARK, lw=1.3, alpha=0.45)

            ax.set_ylim(0, 1.20)
            ax.set_xticks(xs)
            ax.set_xticklabels([str(int(x)) for x in xs])
            ax.grid(True, axis="y", alpha=0.6)

            if row_i == 0:
                ax.set_title(MODALITY_LABELS[mod], fontsize=11)
            if col_i == 0:
                ax.set_ylabel(f"HLA-{gene}\nAccuracy of majority-vote call", fontsize=9)
            if row_i == len(GENE_ORDER) - 1:
                ax.set_xlabel("Tools agreeing on the majority call", fontsize=9)

    fig.suptitle(
        "Why majority voting works: accuracy by tool agreement level and HLA gene\n"
        "(HLA-A, HLA-B, and HLA-C shown separately within each modality)",
        fontsize=13, fontweight="bold"
    )
    fig.tight_layout()
    wm(fig)
    save(fig, out_dir, "figure_03_agreement_vs_accuracy")
    print("  Figure 3 saved.")

# ── Figure 4: Per-gene accuracy ───────────────────────────────────────────────
def fig4_per_gene(tables, out_dir):
    """Per-gene (HLA-A, B, C) accuracy for tools + MajorityVote."""
    pg = tables["per_gene"].copy()
    if pg.empty: return

    pg["overall_correct_call_rate"] = pd.to_numeric(
        pg.get("overall_correct_call_rate", pd.Series()), errors="coerce")

    genes = ["A","B","C"]
    fig, axes = plt.subplots(1, 3, figsize=(16, 5), sharey=False)

    for ax, mod in zip(axes, MODALITY_ORDER):
        sub = pg[pg["modality"] == mod].copy()
        if sub.empty: continue

        x = np.arange(len(genes))
        width = 0.08
        n_tools = len(TOOL_ORDER) + 1  # tools + MajorityVote

        all_methods = TOOL_ORDER + ["MajorityVote"]
        offsets = np.linspace(-(n_tools-1)*width/2, (n_tools-1)*width/2, n_tools)

        for i, method in enumerate(all_methods):
            vals = []
            for gene in genes:
                r = sub[(sub["method"] == method) & (sub["gene"] == gene)]
                if r.empty:
                    vals.append(0)
                else:
                    vals.append(float(r["overall_correct_call_rate"].iloc[0]))
            color = TOOL_COLORS.get(method, NEUTRAL_MID)
            lw = 2.5 if method == "MajorityVote" else 0.3
            ax.bar(x + offsets[i], vals, width=width*0.9, color=color,
                   edgecolor="white" if method != "MajorityVote" else "#8F3A3D",
                   linewidth=lw, alpha=1.0 if method == "MajorityVote" else 0.75,
                   zorder=3 if method == "MajorityVote" else 2,
                   label=method if mod == "wgs" else "")

        ax.set_xticks(x)
        ax.set_xticklabels([f"HLA-{g}" for g in genes], fontsize=10, fontweight="bold")
        ax.set_ylim(0, 1.12)
        ax.set_ylabel("Overall correct-call rate", fontsize=9)
        ax.set_title(MODALITY_LABELS[mod], fontsize=11)

    # Single legend for all panels
    axes[0].legend(loc="upper left", fontsize=7, ncol=1,
                   title="Tool / Method", title_fontsize=8)

    fig.suptitle("Per-locus accuracy: HLA-A, HLA-B, HLA-C\nall tools and majority voting",
                 fontsize=13, fontweight="bold")
    fig.tight_layout()
    wm(fig)
    save(fig, out_dir, "figure_04_per_gene_accuracy")
    print("  Figure 4 saved.")

# ── Figure 5: Cross-modality accuracy landscape ───────────────────────────────
def fig5_crossmodal_landscape(tables, out_dir):
    """Which tools work where — complete tool × modality heatmap."""
    summ = tables["summary"].copy()
    meth = tables["method"].copy()
    if summ.empty: return

    summ["overall_correct_call_rate"] = pd.to_numeric(
        summ.get("overall_correct_call_rate", pd.Series()), errors="coerce")
    meth["overall_correct_call_rate"] = pd.to_numeric(
        meth.get("overall_correct_call_rate", pd.Series()), errors="coerce")

    all_methods = TOOL_ORDER + ["MajorityVote"]
    mods = MODALITY_ORDER

    # Build matrix
    mat = np.full((len(all_methods), len(mods)), np.nan)
    for i, method in enumerate(all_methods):
        for j, mod in enumerate(mods):
            if method == "MajorityVote":
                r = meth[(meth["method"] == "MajorityVote") & (meth["modality"] == mod)]
            else:
                r = summ[(summ["tool"] == method) & (summ["modality"] == mod)]
            if not r.empty:
                mat[i, j] = float(r["overall_correct_call_rate"].iloc[0])

    fig, ax = plt.subplots(figsize=(10, 7))
    cmap = HEATMAP_CMAP
    cmap.set_bad(color=MISSING_CELL)
    im = ax.imshow(mat, cmap=cmap, vmin=0, vmax=1, aspect="auto")

    # Annotate cells
    for i in range(len(all_methods)):
        for j in range(len(mods)):
            val = mat[i, j]
            if np.isnan(val):
                ax.text(j, i, "—", ha="center", va="center", fontsize=10, color=NEUTRAL_MID)
            else:
                txt_color = "white" if val < 0.22 or val > 0.63 else NEUTRAL_DARK
                weight = "bold" if all_methods[i] == "MajorityVote" else "normal"
                ax.text(j, i, f"{val:.3f}", ha="center", va="center",
                        fontsize=10, fontweight=weight, color=txt_color)

    # Separator before MajorityVote
    ax.axhline(len(TOOL_ORDER) - 0.5, color=NEUTRAL_DARK, lw=2)

    ax.set_xticks(range(len(mods)))
    ax.set_xticklabels([MODALITY_LABELS[m].replace("\n", " ") for m in mods],
                       fontsize=11, fontweight="bold")
    ax.set_yticks(range(len(all_methods)))
    ylab = [f"{m}  [{TOOL_DESIGN.get(m,'—')}]" if m != "MajorityVote" else "★ MajorityVote"
            for m in all_methods]
    ax.set_yticklabels(ylab, fontsize=9)
    for tick, m in zip(ax.get_yticklabels(), all_methods):
        tick.set_fontweight("bold" if m == "MajorityVote" else "normal")
        if m == "MajorityVote":
            tick.set_color(TOOL_COLORS["MajorityVote"])

    cb = fig.colorbar(im, ax=ax, fraction=0.025, pad=0.03)
    cb.set_label("Overall correct-call rate", fontsize=9)

    # Tool design annotation
    for i, method in enumerate(all_methods):
        if method == "MajorityVote": continue
        design = TOOL_DESIGN.get(method, "")
        if design == "RNA":
            ax.add_patch(mpatches.FancyBboxPatch(
                (-0.45, i-0.45), 0.1, 0.9, boxstyle="round,pad=0.05",
                fc=RNA_MARKER, ec="none", alpha=0.45, transform=ax.transData))
        elif design == "DNA":
            ax.add_patch(mpatches.FancyBboxPatch(
                (-0.45, i-0.45), 0.1, 0.9, boxstyle="round,pad=0.05",
                fc=DNA_MARKER, ec="none", alpha=0.45, transform=ax.transData))

    ax.set_title("Cross-modality HLA typing accuracy landscape\n"
                 "(★ MajorityVote shown separately; tool design: DNA=blue, RNA=red strip)",
                 fontsize=12, fontweight="bold", pad=12)
    wm(fig)
    save(fig, out_dir, "figure_05_crossmodal_landscape")
    print("  Figure 5 saved.")

# ── Figure 6: Best-performing tool by gene and modality ──────────────────────
def fig6_best_tool_by_gene_modality(tables, out_dir):
    """Winner-focused per-gene comparison across modalities."""
    pg = tables["per_gene"].copy()
    if pg.empty:
        return

    pg["overall_correct_call_rate"] = pd.to_numeric(
        pg.get("overall_correct_call_rate", pd.Series()), errors="coerce")
    pg = pg[pg["method"] != "WeightedConsensus"].copy()

    genes = ["A", "B", "C"]
    fig, axes = plt.subplots(1, 3, figsize=(17, 6), sharey=True)

    for ax, mod in zip(axes, MODALITY_ORDER):
        sub = pg[pg["modality"] == mod].copy()
        if sub.empty:
            continue

        methods_present = [
            method for method in TOOL_ORDER + ["MajorityVote"]
            if not sub[sub["method"] == method].empty
        ]
        if not methods_present:
            continue

        n_methods = len(methods_present)
        width = min(0.12, 0.82 / max(n_methods, 1))
        x = np.arange(len(genes))
        offsets = np.linspace(-(n_methods - 1) * width / 2, (n_methods - 1) * width / 2, n_methods)

        winners = {}
        for gene in genes:
            gene_sub = sub[sub["gene"] == gene].copy()
            if gene_sub.empty:
                continue
            winner_idx = gene_sub["overall_correct_call_rate"].idxmax()
            winner_row = gene_sub.loc[winner_idx]
            winners[gene] = (
                winner_row["method"],
                float(winner_row["overall_correct_call_rate"]),
            )

        for i, method in enumerate(methods_present):
            vals = []
            alphas = []
            edgecolors = []
            linew = []
            for gene in genes:
                r = sub[(sub["method"] == method) & (sub["gene"] == gene)]
                if r.empty:
                    vals.append(np.nan)
                    alphas.append(0.0)
                    edgecolors.append("white")
                    linew.append(0.0)
                    continue
                acc = float(r["overall_correct_call_rate"].iloc[0])
                vals.append(acc)
                is_winner = winners.get(gene, (None, None))[0] == method
                alphas.append(1.0 if is_winner else 0.42)
                edgecolors.append("#111827" if is_winner else "white")
                linew.append(1.8 if is_winner else 0.4)

            bar_positions = x + offsets[i]
            bars = ax.bar(
                bar_positions,
                np.nan_to_num(vals, nan=0.0),
                width=width * 0.92,
                color=TOOL_COLORS.get(method, NEUTRAL_MID),
                edgecolor=edgecolors,
                linewidth=0.4,
                zorder=3 if method == "MajorityVote" else 2,
                label=method if mod == "wgs" else "",
            )

            for bar, alpha, lw, ec, val in zip(bars, alphas, linew, edgecolors, vals):
                if np.isnan(val):
                    bar.set_visible(False)
                    continue
                bar.set_alpha(alpha)
                bar.set_linewidth(lw)
                bar.set_edgecolor(ec)

        for j, gene in enumerate(genes):
            winner = winners.get(gene)
            if winner is None:
                continue
            winner_method, winner_acc = winner
            ax.text(
                x[j],
                min(winner_acc + 0.065, 1.11),
                f"Best: {winner_method}\n{winner_acc:.3f}",
                ha="center",
                va="bottom",
                fontsize=7.8,
                fontweight="bold",
                color=NEUTRAL_DARK,
                bbox=dict(boxstyle="round,pad=0.22", fc="white", ec="#D1D5DB", alpha=0.95),
            )

        ax.set_xticks(x)
        ax.set_xticklabels([f"HLA-{g}" for g in genes], fontsize=10, fontweight="bold")
        ax.set_ylim(0, 1.14)
        ax.set_title(MODALITY_LABELS[mod], fontsize=11)
        ax.set_ylabel("Overall correct-call rate", fontsize=9)
        ax.grid(True, axis="y", alpha=0.55)

    axes[0].legend(
        loc="upper left",
        fontsize=7,
        ncol=1,
        title="Tool / Method",
        title_fontsize=8,
    )

    fig.suptitle(
        "Best-performing HLA typing methods by gene and sequencing modality\n"
        "(individual tools plus MajorityVote; winning method highlighted in each HLA-A/B/C cluster)",
        fontsize=13,
        fontweight="bold",
    )
    fig.tight_layout()
    wm(fig)
    save(fig, out_dir, "figure_06_best_tool_by_gene_modality")
    print("  Figure 6 saved.")

# ── Main ─────────────────────────────────────────────────────────────────────
def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--tables-dir", required=True)
    p.add_argument("--out-dir", required=True)
    return p.parse_args()

# ── Figure 7: WES+RNA bimodal joint consensus comparison ─────────────────────
def _bimodal_per_gene_acc(detail_df, mv_df, wc_df, genes):
    """Compute per-gene accuracy for 6 comparison methods, restricted to bimodal-eligible samples."""
    bimodal_samples = set(detail_df["sample"].unique())

    def _acc(df, gene, modality_filter, method_filter):
        sub = df[
            (df["gene"] == gene) &
            (df["modality"] == modality_filter) &
            (df["method"] == method_filter) &
            (df["sample"].isin(bimodal_samples)) &
            (df["is_callable"].astype(str) == "1")
        ]
        if sub.empty:
            return np.nan, np.nan, np.nan
        n = len(sub)
        k = sub["is_correct"].astype(int).sum()
        acc = k / n
        # Wilson interval
        z = 1.96
        denom = 1 + z**2 / n
        centre = (k + z**2 / 2) / (n * denom)
        margin = z * np.sqrt(k * (n - k) / n + z**2 / 4) / (n * denom)
        return acc, centre - margin, centre + margin

    def _acc_bimodal(df, gene, method):
        sub = df[
            (df["gene"] == gene) &
            (df["method"] == method) &
            (df["is_callable"].astype(str) == "1")
        ]
        if sub.empty:
            return np.nan, np.nan, np.nan
        n = len(sub)
        k = sub["is_correct"].astype(int).sum()
        acc = k / n
        z = 1.96
        denom = 1 + z**2 / n
        centre = (k + z**2 / 2) / (n * denom)
        margin = z * np.sqrt(k * (n - k) / n + z**2 / 4) / (n * denom)
        return acc, centre - margin, centre + margin

    result = {}
    for gene in genes:
        result[gene] = {
            "WES_MajorityVote":      _acc(mv_df,  gene, "wes",    "MajorityVote"),
            "WES_WeightedConsensus": _acc(wc_df,  gene, "wes",    "WeightedConsensus"),
            "RNA_MajorityVote":      _acc(mv_df,  gene, "rnaseq", "MajorityVote"),
            "RNA_WeightedConsensus": _acc(wc_df,  gene, "rnaseq", "WeightedConsensus"),
            "Bimodal_MajorityVote":      _acc_bimodal(detail_df, gene, "BimodalMajorityVote"),
            "Bimodal_WeightedConsensus": _acc_bimodal(detail_df, gene, "BimodalWeightedConsensus"),
        }
    return result


def fig7_bimodal_comparison(tables, out_dir):
    """Grouped bar chart: per-gene (A/B/C) accuracy for WES-only, RNA-only, and bimodal methods."""
    detail = tables.get("bimodal_detail", pd.DataFrame())
    mv_df  = tables.get("mv_calls",      pd.DataFrame())
    wc_df  = tables.get("wc_calls",      pd.DataFrame())
    agg_df = tables.get("bimodal",       pd.DataFrame())

    if detail is None or detail.empty:
        print("  [skip] bimodal_wes_rna_consensus.tsv not found.")
        return

    for df in [detail, mv_df, wc_df]:
        for col in ["is_callable", "is_correct"]:
            if col in df.columns:
                df[col] = pd.to_numeric(df[col], errors="coerce")

    genes = ["A", "B", "C"]
    comparisons = [
        ("WES_MajorityVote",          "WES\nMajority",          "#4C78A8"),
        ("WES_WeightedConsensus",     "WES\nWeighted",          "#2A6099"),
        ("RNA_MajorityVote",          "RNA\nMajority",          "#E69F00"),
        ("RNA_WeightedConsensus",     "RNA\nWeighted",          "#B87A00"),
        ("Bimodal_MajorityVote",      "WES+RNA\nMajority",      "#C44E52"),
        ("Bimodal_WeightedConsensus", "WES+RNA\nWeighted",      "#8B1A1A"),
    ]

    per_gene = _bimodal_per_gene_acc(detail, mv_df, wc_df, genes)

    x = np.arange(len(genes))
    n_groups = len(comparisons)
    width = 0.12
    offsets = np.linspace(-(n_groups - 1) * width / 2, (n_groups - 1) * width / 2, n_groups)

    fig, ax = plt.subplots(figsize=(12, 6))

    for i, (comp_key, label, color) in enumerate(comparisons):
        vals, yerr_lo, yerr_hi = [], [], []
        for gene in genes:
            acc, ci_lo, ci_hi = per_gene[gene].get(comp_key, (np.nan, np.nan, np.nan))
            vals.append(acc)
            yerr_lo.append(0 if np.isnan(acc) else acc - ci_lo)
            yerr_hi.append(0 if np.isnan(acc) else ci_hi - acc)

        positions = x + offsets[i]
        is_bimodal = "Bimodal" in comp_key
        ax.bar(
            positions,
            np.nan_to_num(vals, nan=0.0),
            width=width * 0.9,
            color=color,
            label=label.replace("\n", " "),
            edgecolor="#111827" if is_bimodal else "white",
            linewidth=1.6 if is_bimodal else 0.4,
            zorder=3 if is_bimodal else 2,
        )
        for pos, acc, lo, hi in zip(positions, vals, yerr_lo, yerr_hi):
            if np.isnan(acc):
                continue
            ax.errorbar(pos, acc, yerr=[[lo], [hi]], fmt="none",
                        ecolor="#374151", elinewidth=1.0, capsize=2.5, capthick=1.0, zorder=4)

    for gap in [0.5, 1.5]:
        ax.axvline(gap, color=GRID_COLOR, lw=1.2, zorder=1)

    ax.set_xticks(x)
    ax.set_xticklabels([f"HLA-{g}" for g in genes], fontsize=11, fontweight="bold")
    ax.set_ylabel("Correct-call rate (callable loci)", fontsize=10)
    ax.set_ylim(0, 1.08)
    ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda v, _: f"{v:.0%}"))

    # Sample count from aggregate table
    n_samples = len(detail["sample"].unique())
    if not agg_df.empty:
        n_row = agg_df[agg_df["comparison"] == "Bimodal_MajorityVote"]
        if not n_row.empty:
            n_samples = int(n_row["sample_count"].iloc[0])
    ax.text(0.99, 0.02, f"n = {n_samples} samples (WES + RNA-seq); DRB1/DQB1 in Table S1",
            transform=ax.transAxes, ha="right", va="bottom",
            fontsize=8, color=NEUTRAL_MID)

    handles, labels_leg = ax.get_legend_handles_labels()
    ax.legend(handles, labels_leg, ncol=3, fontsize=8, loc="upper left",
              framealpha=0.9, title="Method", title_fontsize=8)

    ax.set_title("WES + RNA-seq bimodal joint HLA consensus vs. single-modality baselines\n"
                 "(HLA-A, -B, -C; same n per gene; bold outlines = bimodal joint)",
                 fontsize=12, fontweight="bold", pad=10)
    wm(fig)
    save(fig, out_dir, "figure_07_bimodal_comparison")
    print("  Figure 7 saved.")


def main():
    args = parse_args()
    td = Path(args.tables_dir)
    od = Path(args.out_dir)
    apply_style()
    print(f"Reading: {td}\nWriting: {od}\n")
    tables = load_tables(td)

    print("Figure 1 — Accuracy overview...")
    fig1_accuracy_overview(tables, od)

    print("Figure 2 — MajorityVote advantage...")
    fig2_majority_vote_advantage(tables, od)

    print("Figure 3 — Agreement → accuracy...")
    fig3_agreement_accuracy(tables, od)

    print("Figure 4 — Per-gene accuracy...")
    fig4_per_gene(tables, od)

    print("Figure 5 — Cross-modality landscape...")
    fig5_crossmodal_landscape(tables, od)

    print("Figure 6 — Best tool by gene and modality...")
    fig6_best_tool_by_gene_modality(tables, od)

    print("Figure 7 — WES+RNA bimodal joint consensus...")
    fig7_bimodal_comparison(tables, od)

    print(f"\nAll 7 figures written to {od}/")

if __name__ == "__main__":
    main()
