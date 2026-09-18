#!/usr/bin/env python3.11
"""
generate_new_figures.py — New publication-ready scientific figures for PIHLA.

Five new figures + two summary tables that expose findings not covered by
generate_figures_v2.py:
  A. Cross-modality accuracy landscape (tool × modality bubble heatmap)
  B. Calibration reliability diagrams (predicted vs observed accuracy)
  C. Resolution profile (2-field / 3-field / G-group / P-group)
  D. Voting dynamics (WGS outvoting mechanism)
  E. Per-modality abstention efficiency curves

Plus:
  Table A. Comprehensive accuracy with 95% Wilson CIs
  Table B. Resolution-stratified performance

Usage:
    python3.11 bin/generate_new_figures.py \
        --tables-dir analysis/benchmark_trimodal_all_samples/tables \
        --figures-dir analysis/figures_new
"""

import argparse
import csv
import math
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np
import pandas as pd

# ---------------------------------------------------------------------------
# Constants — same palette as generate_figures_v2.py
# ---------------------------------------------------------------------------

TOOL_COLORS = {
    "ArcasHLA":          "#6b7280",
    "HLA-HD":            "#3b82f6",
    "Kourami":           "#8b5cf6",
    "OptiType":          "#06b6d4",
    "POLYSOLVER":        "#f97316",
    "Seq2HLA":           "#84cc16",
    "SpecHLA":           "#64748b",
    "T1K":               "#0ea5e9",
    "MajorityVote":      "#d97706",
    "WeightedConsensus": "#0f766e",
}

GUARDRAIL_COLORS = {
    "applied":          "#16a34a",
    "poor_calibration": "#dc2626",
    "no_confidence":    "#9ca3af",
}

MODALITY_LABELS = {"wgs": "WGS", "wes": "WES", "rnaseq": "RNA-seq"}
MODALITY_ORDER  = ["wgs", "wes", "rnaseq"]

TOOL_DESIGN = {
    "ArcasHLA":  "RNA",
    "HLA-HD":    "Both",
    "Kourami":   "DNA",
    "OptiType":  "DNA/RNA",
    "POLYSOLVER":"DNA",
    "Seq2HLA":   "RNA",
    "SpecHLA":   "DNA",
    "T1K":       "Both",
}
DESIGN_COLORS = {"DNA": "#3b82f6", "RNA": "#ef4444", "Both": "#16a34a", "DNA/RNA": "#8b5cf6"}

TOOL_NOTES = {
    ("ArcasHLA", "wgs"):   "RNA design; off-label WGS",
    ("ArcasHLA", "wes"):   "RNA design; off-label WES",
    ("SpecHLA",  "rnaseq"):"DNA design; N-masking on RNA",
    ("Seq2HLA",  "rnaseq"):"n=12 callable (path coverage)",
}

WATERMARK = "PIHLA v2.0 · IMGT/HLA 3.59.0 · n=99 WGS / n=51 WES / n=50 RNA"
DPI = 300

# ---------------------------------------------------------------------------
# Style + save helpers
# ---------------------------------------------------------------------------

def apply_style():
    plt.rcParams.update({
        "font.family":       "DejaVu Sans",
        "font.size":         10,
        "axes.titlesize":    11,
        "axes.titleweight":  "bold",
        "axes.labelsize":    10,
        "xtick.labelsize":   9,
        "ytick.labelsize":   9,
        "axes.spines.top":   False,
        "axes.spines.right": False,
        "axes.grid":         True,
        "grid.color":        "#e5e7eb",
        "grid.linewidth":    0.6,
        "legend.fontsize":   9,
        "legend.framealpha": 0.9,
        "figure.dpi":        DPI,
        "savefig.dpi":       DPI,
        "savefig.bbox":      "tight",
    })


def watermark(fig):
    fig.text(0.99, 0.005, WATERMARK, ha="right", va="bottom",
             fontsize=7, color="#9ca3af", transform=fig.transFigure)


def save_fig(fig, out_dir: Path, stem: str):
    out_dir.mkdir(parents=True, exist_ok=True)
    for ext in ("png", "pdf"):
        fig.savefig(str(out_dir / f"{stem}.{ext}"), dpi=DPI if ext == "png" else None)
    plt.close(fig)


def wilson_ci(k, n, z=1.96):
    if n == 0:
        return 0.0, 0.0
    p = k / n
    denom = 1 + z**2 / n
    centre = (p + z**2 / (2 * n)) / denom
    margin = z * math.sqrt(p * (1 - p) / n + z**2 / (4 * n**2)) / denom
    return max(0.0, centre - margin), min(1.0, centre + margin)

# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------

def load_tables(tables_dir: Path) -> dict:
    files = {
        "summary":          "summary_full_cohort.tsv",
        "method_comp":      "method_comparison.tsv",
        "ambiguity":        "ambiguity_summary.tsv",
        "conf_bins":        "confidence_bin_summary.tsv",
        "conf_cal":         "confidence_calibration_summary.tsv",
        "tool_weights":     "tool_confidence_weights.tsv",
        "harmonized":       "harmonized_benchmark_rows.tsv",
        "wc_calls":         "weighted_consensus_calls.tsv",
        "mv_calls":         "majority_vote_baseline.tsv",
    }
    tables = {}
    for key, fname in files.items():
        p = tables_dir / fname
        if not p.exists():
            tables[key] = pd.DataFrame()
            continue
        df = pd.read_csv(str(p), sep="\t", dtype=str)
        for col in df.columns:
            try:
                df[col] = pd.to_numeric(df[col], errors="ignore")
            except Exception:
                pass
        tables[key] = df
    return tables

# ---------------------------------------------------------------------------
# Figure A — Cross-modality accuracy landscape
# ---------------------------------------------------------------------------

def fig_A_accuracy_landscape(tables: dict, out_dir: Path):
    df = tables["summary"].copy()
    if df.empty:
        return
    df["overall_correct_call_rate"] = pd.to_numeric(df["overall_correct_call_rate"], errors="coerce").fillna(0)
    df["callable_rate"] = pd.to_numeric(df["callable_rate"], errors="coerce").fillna(0)

    tools_present = [t for t in ["ArcasHLA","HLA-HD","Kourami","OptiType","POLYSOLVER","Seq2HLA","SpecHLA","T1K"] if t in df["tool"].values]
    mods = [m for m in MODALITY_ORDER if m in df["modality"].values]

    # Sort tools by mean accuracy across present modalities
    mean_acc = {t: df[df["tool"]==t]["overall_correct_call_rate"].mean() for t in tools_present}
    tools_sorted = sorted(tools_present, key=lambda t: mean_acc[t], reverse=True)

    fig, axes = plt.subplots(1, 1, figsize=(10, 7))
    ax = axes

    cmap = plt.cm.RdYlGn
    norm = mcolors.Normalize(vmin=0, vmax=1)

    for xi, mod in enumerate(mods):
        for yi, tool in enumerate(tools_sorted):
            row = df[(df["tool"]==tool) & (df["modality"]==mod)]
            if row.empty:
                ax.scatter(xi, yi, s=30, color="#e5e7eb", zorder=2)
                continue
            acc = float(row["overall_correct_call_rate"].iloc[0])
            cr  = float(row["callable_rate"].iloc[0])
            size = max(40, acc * 1800)
            color = cmap(norm(cr))
            ax.scatter(xi, yi, s=size, color=color, edgecolors="#374151", linewidths=0.5, zorder=3)
            ax.text(xi, yi, f"{acc:.2f}", ha="center", va="center", fontsize=8,
                    fontweight="bold", color="white" if cr > 0.4 else "#374151", zorder=4)

    # Design scope strip on right
    for yi, tool in enumerate(tools_sorted):
        design = TOOL_DESIGN.get(tool, "")
        dc = DESIGN_COLORS.get(design, "#9ca3af")
        ax.annotate(f" {design}", xy=(len(mods) - 0.05, yi), xycoords="data",
                    fontsize=8, color=dc, fontweight="bold", va="center")

    ax.set_xticks(range(len(mods)))
    ax.set_xticklabels([MODALITY_LABELS[m] for m in mods], fontsize=10, fontweight="bold")
    ax.set_yticks(range(len(tools_sorted)))
    ax.set_yticklabels(tools_sorted, fontsize=9)
    ax.set_xlim(-0.6, len(mods) + 0.8)
    ax.set_ylim(-0.7, len(tools_sorted) - 0.3)
    ax.set_title("Cross-modality HLA typing accuracy landscape", pad=12)
    ax.set_xlabel("Sequencing modality", labelpad=8)
    ax.grid(False)

    # Colorbar for callable_rate
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax, fraction=0.025, pad=0.02)
    cbar.set_label("Callable rate", fontsize=9)

    # Bubble size legend
    for acc_val, label in [(0.2, "0.20"), (0.5, "0.50"), (0.9, "0.90")]:
        ax.scatter([], [], s=max(40, acc_val * 1800), color="#9ca3af",
                   edgecolors="#374151", linewidths=0.5, label=f"Accuracy {label}")
    ax.legend(title="Bubble size = overall accuracy", loc="lower right", framealpha=0.9, fontsize=8)

    watermark(fig)
    save_fig(fig, out_dir, "figure_A_accuracy_landscape")
    print("  Figure A saved.")

# ---------------------------------------------------------------------------
# Figure B — Calibration reliability diagrams
# ---------------------------------------------------------------------------

def fig_B_calibration_reliability(tables: dict, out_dir: Path):
    bins_df = tables["conf_bins"].copy()
    cal_df  = tables["conf_cal"].copy()
    wt_df   = tables["tool_weights"].copy()
    if bins_df.empty:
        return

    # Tools with WGS confidence data
    wgs_bins = bins_df[bins_df["modality"] == "wgs"].copy()
    tools_with_data = sorted(wgs_bins["tool"].unique())
    n = len(tools_with_data)
    if n == 0:
        return

    ncols = min(3, n)
    nrows = math.ceil(n / ncols)
    fig, axes = plt.subplots(nrows, ncols, figsize=(4.5 * ncols, 4 * nrows), squeeze=False)

    for idx, tool in enumerate(tools_with_data):
        ax = axes[idx // ncols][idx % ncols]
        tdata = wgs_bins[wgs_bins["tool"] == tool].copy()
        tdata["mean_confidence"] = pd.to_numeric(tdata["mean_confidence"], errors="coerce")
        tdata["observed_accuracy"] = pd.to_numeric(tdata["observed_accuracy"], errors="coerce")
        tdata = tdata.dropna(subset=["mean_confidence", "observed_accuracy"])
        tdata = tdata.sort_values("mean_confidence")

        # Perfect calibration diagonal
        ax.plot([0, 1], [0, 1], "k--", lw=1, alpha=0.4, label="Perfect calibration")

        if not tdata.empty:
            color = TOOL_COLORS.get(tool, "#6b7280")
            ax.scatter(tdata["mean_confidence"], tdata["observed_accuracy"],
                       s=tdata["n_rows"].astype(float).fillna(10).clip(upper=300),
                       color=color, alpha=0.85, edgecolors="#374151", linewidths=0.4, zorder=3)
            ax.plot(tdata["mean_confidence"], tdata["observed_accuracy"],
                    color=color, lw=1.5, alpha=0.7, zorder=2)

        # Annotations from calibration summary
        cal_row = cal_df[(cal_df["tool"] == tool) & (cal_df["modality"] == "wgs")]
        wt_row  = wt_df[(wt_df["tool"] == tool) & (wt_df["modality"] == "wgs")]
        guard_status = wt_row["guardrail_status"].iloc[0] if not wt_row.empty and "guardrail_status" in wt_row.columns else ""
        guard_color = GUARDRAIL_COLORS.get(guard_status, "#9ca3af")

        if not cal_row.empty:
            brier = pd.to_numeric(cal_row["brier_score"].iloc[0], errors="coerce")
            ece   = pd.to_numeric(cal_row["expected_calibration_error"].iloc[0], errors="coerce")
            stats_txt = f"Brier={brier:.3f}\nECE={ece:.3f}"
            ax.text(0.03, 0.97, stats_txt, transform=ax.transAxes,
                    fontsize=8, va="top", color="#374151",
                    bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="#e5e7eb", alpha=0.85))

        guard_label = guard_status.replace("_", " ").title() if guard_status else "No data"
        ax.set_title(f"{tool}\n[{guard_label}]", color=guard_color, fontsize=10)
        ax.set_xlim(-0.02, 1.02)
        ax.set_ylim(-0.02, 1.02)
        ax.set_xlabel("Mean confidence score", fontsize=9)
        ax.set_ylabel("Observed accuracy", fontsize=9)
        ax.set_aspect("equal")

    # Hide unused subplots
    for idx in range(n, nrows * ncols):
        axes[idx // ncols][idx % ncols].set_visible(False)

    fig.suptitle("Confidence calibration reliability diagrams (WGS)", fontsize=13, fontweight="bold", y=1.01)
    watermark(fig)
    save_fig(fig, out_dir, "figure_B_calibration_reliability")
    print("  Figure B saved.")

# ---------------------------------------------------------------------------
# Figure C — Resolution profile (2-field / 3-field / G/P-group)
# ---------------------------------------------------------------------------

def fig_C_resolution_profile(tables: dict, out_dir: Path):
    df = tables["ambiguity"].copy()
    if df.empty:
        return

    for col in ["exact_2field_rate", "exact_3field_rate", "g_group_match_rate", "p_group_match_rate"]:
        df[col] = pd.to_numeric(df.get(col, pd.Series(dtype=float)), errors="coerce").fillna(0)

    mods = [m for m in MODALITY_ORDER if m in df["modality"].values]
    fig, axes = plt.subplots(1, len(mods), figsize=(5.5 * len(mods), 7), sharey=False)
    if len(mods) == 1:
        axes = [axes]

    for ax, mod in zip(axes, mods):
        sub = df[df["modality"] == mod].copy()
        # Sort by 2-field rate
        sub = sub.sort_values("exact_2field_rate", ascending=True)
        tools = sub["tool"].tolist()
        y = np.arange(len(tools))
        h = 0.18

        bars_2field = sub["exact_2field_rate"].values
        bars_3field = sub["exact_3field_rate"].values
        bars_g      = sub["g_group_match_rate"].values
        bars_p      = sub["p_group_match_rate"].values

        ax.barh(y + 1.5*h, bars_2field, h*1.4, color="#3b82f6", alpha=0.85, label="2-field exact")
        ax.barh(y + 0.5*h, bars_3field, h*1.4, color="#06b6d4", alpha=0.85, label="3-field exact")
        ax.barh(y - 0.5*h, bars_g,      h*1.4, color="#8b5cf6", alpha=0.7,  label="G-group")
        ax.barh(y - 1.5*h, bars_p,      h*1.4, color="#d97706", alpha=0.7,  label="P-group")

        # Annotate 3-field bars for tools that achieve it
        for yi, (tool, v3, v2) in enumerate(zip(tools, bars_3field, bars_2field)):
            if v3 > 0.05:
                ax.text(v3 + 0.01, yi + 0.5*h, f"{v3:.2f}", va="center", fontsize=7.5,
                        color="#0e7490", fontweight="bold")

        ax.set_yticks(y)
        ax.set_yticklabels(tools, fontsize=9)
        ax.set_xlim(0, 1.12)
        ax.set_xlabel("Rate", fontsize=9)
        ax.set_title(MODALITY_LABELS[mod], fontsize=11, fontweight="bold")
        ax.axvline(x=0, color="#374151", lw=0.5)
        if ax == axes[0]:
            ax.legend(loc="lower right", fontsize=8)

    fig.suptitle("Allele resolution profile across tools and modalities", fontsize=13, fontweight="bold")
    fig.tight_layout()
    watermark(fig)
    save_fig(fig, out_dir, "figure_C_resolution_profile")
    print("  Figure C saved.")

# ---------------------------------------------------------------------------
# Figure D — Voting dynamics (WGS)
# ---------------------------------------------------------------------------

def fig_D_voting_dynamics(tables: dict, out_dir: Path):
    harm = tables["harmonized"].copy()
    wc   = tables["wc_calls"].copy()
    if harm.empty or wc.empty:
        return

    harm_wgs = harm[harm["modality"] == "wgs"].copy()
    harm_wgs["is_correct"] = pd.to_numeric(harm_wgs["is_correct"], errors="coerce").fillna(0)
    wc_wgs   = wc[wc["modality"] == "wgs"].copy()
    wc_wgs["is_correct"] = pd.to_numeric(wc_wgs["is_correct"], errors="coerce").fillna(0)

    # For each (sample, gene): count how many tools called correctly
    vote_counts = harm_wgs.groupby(["sample", "gene"])["is_correct"].agg(
        n_correct="sum", n_tools="count"
    ).reset_index()

    # Join with WeightedConsensus outcome
    wc_outcome = wc_wgs[["sample", "gene", "is_correct", "call_status"]].rename(
        columns={"is_correct": "wc_correct", "call_status": "wc_status"})
    merged = vote_counts.merge(wc_outcome, on=["sample", "gene"], how="left")
    merged["wc_label"] = merged.apply(
        lambda r: "WC correct" if r.get("wc_correct") == 1
        else ("WC abstained" if str(r.get("wc_status", "")).startswith("low") else "WC incorrect"),
        axis=1
    )

    # Panel L: histogram of n_correct votes, stacked by WC outcome
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5))

    max_votes = int(merged["n_tools"].max()) if not merged.empty else 6
    bins_x = np.arange(0, max_votes + 1)
    label_order = ["WC correct", "WC incorrect", "WC abstained"]
    colors_wc = {"WC correct": "#16a34a", "WC incorrect": "#dc2626", "WC abstained": "#9ca3af"}

    bottom = np.zeros(max_votes + 1)
    for lbl in label_order:
        sub = merged[merged["wc_label"] == lbl]
        counts = sub["n_correct"].value_counts().reindex(bins_x, fill_value=0).values
        ax1.bar(bins_x, counts, bottom=bottom, color=colors_wc[lbl], alpha=0.85,
                label=lbl, edgecolor="white", linewidth=0.4)
        bottom += counts

    ax1.set_xlabel("Number of tools calling truth allele correctly", fontsize=10)
    ax1.set_ylabel("Number of WGS loci", fontsize=10)
    ax1.set_title("(A) Tool agreement with truth per WGS locus", fontsize=11, fontweight="bold")
    ax1.set_xticks(bins_x)
    ax1.legend(fontsize=9)
    ax1.grid(axis="y")

    # Panel R: for OptiType-correct loci, how many OTHER tools agree?
    optitype_correct = harm_wgs[
        (harm_wgs["tool"] == "OptiType") & (harm_wgs["is_correct"] == 1)
    ][["sample", "gene"]].copy()
    optitype_loci = set(zip(optitype_correct["sample"], optitype_correct["gene"]))

    other_tools = harm_wgs[harm_wgs["tool"] != "OptiType"]
    other_counts = other_tools.groupby(["sample", "gene"])["is_correct"].agg(
        n_other_correct="sum", n_other_tools="count"
    ).reset_index()
    other_counts["optitype_correct_locus"] = other_counts.apply(
        lambda r: (r["sample"], r["gene"]) in optitype_loci, axis=1
    )

    oc_sub = other_counts[other_counts["optitype_correct_locus"]]
    if not oc_sub.empty:
        n_other = int(oc_sub["n_other_tools"].max())
        vote_bins = np.arange(0, n_other + 1)
        counts2 = oc_sub["n_other_correct"].value_counts().reindex(vote_bins, fill_value=0).values
        colors2 = ["#dc2626" if x < n_other / 2 else "#16a34a" for x in vote_bins]
        bars = ax2.bar(vote_bins, counts2, color=colors2, alpha=0.85, edgecolor="white", linewidth=0.4)
        majority_threshold = n_other / 2
        ax2.axvline(x=majority_threshold, color="#374151", ls="--", lw=1.5,
                    label=f"Majority threshold ({majority_threshold:.1f})")
        # Annotate: fraction where majority agrees with OptiType
        n_majority = int(oc_sub[oc_sub["n_other_correct"] > majority_threshold].shape[0])
        n_total = len(oc_sub)
        pct = 100 * n_majority / n_total if n_total > 0 else 0
        ax2.text(0.97, 0.97,
                 f"Other tools agree with\nOptiType majority: {pct:.0f}%\n({n_majority}/{n_total} loci)",
                 transform=ax2.transAxes, ha="right", va="top", fontsize=9,
                 bbox=dict(boxstyle="round,pad=0.4", fc="white", ec="#e5e7eb"))

    ax2.set_xlabel("Other tools agreeing with truth (OptiType-correct loci)", fontsize=10)
    ax2.set_ylabel("Number of WGS loci", fontsize=10)
    ax2.set_title("(B) When OptiType is correct:\nhow many other tools agree?", fontsize=11, fontweight="bold")
    ax2.legend(fontsize=9)
    ax2.grid(axis="y")

    fig.suptitle("WGS voting dynamics: the ensemble outvoting problem", fontsize=13, fontweight="bold")
    fig.tight_layout()
    watermark(fig)
    save_fig(fig, out_dir, "figure_D_voting_dynamics")
    print("  Figure D saved.")

# ---------------------------------------------------------------------------
# Figure E — Per-modality abstention efficiency curves
# ---------------------------------------------------------------------------

def fig_E_abstention_efficiency(tables: dict, out_dir: Path):
    wc_all = tables["wc_calls"].copy()
    mv_all = tables["mv_calls"].copy()
    summ   = tables["summary"].copy()
    mc     = tables["method_comp"].copy()
    if wc_all.empty:
        return

    thresholds = np.arange(0.30, 0.95, 0.05)
    mods = [m for m in MODALITY_ORDER if m in wc_all["modality"].values]

    fig, axes = plt.subplots(1, len(mods), figsize=(5 * len(mods), 5), sharey=False)
    if len(mods) == 1:
        axes = [axes]

    for ax, mod in zip(axes, mods):
        wc = wc_all[wc_all["modality"] == mod].copy()
        wc["support_fraction"] = pd.to_numeric(wc.get("support_fraction", pd.Series()), errors="coerce")
        wc["is_correct"] = pd.to_numeric(wc["is_correct"], errors="coerce").fillna(0)
        wc["is_callable"] = pd.to_numeric(wc["is_callable"], errors="coerce").fillna(0)

        curve_x, curve_y = [], []
        for thr in thresholds:
            called = wc[wc["support_fraction"] >= thr] if "support_fraction" in wc.columns else wc[wc["is_callable"] == 1]
            total_loci = len(wc)
            n_called = len(called)
            call_rate = n_called / total_loci if total_loci > 0 else 0
            acc_among = called["is_correct"].mean() if n_called > 0 else 0
            curve_x.append(call_rate)
            curve_y.append(acc_among)

        # Main curve
        ax.plot(curve_x, curve_y, "o-", color="#0f766e", lw=2, ms=5, label="WeightedConsensus (varying θ)", zorder=3)

        # Annotate current default (θ=0.55)
        default_idx = np.argmin(np.abs(thresholds - 0.55))
        if default_idx < len(curve_x):
            ax.scatter([curve_x[default_idx]], [curve_y[default_idx]],
                       s=120, color="#0f766e", zorder=5, edgecolors="white", linewidths=1.5)
            ax.annotate(f"  θ=0.55\n  ({curve_x[default_idx]:.2f}, {curve_y[default_idx]:.2f})",
                        xy=(curve_x[default_idx], curve_y[default_idx]),
                        fontsize=7.5, color="#0f766e")

        # MajorityVote reference
        mv_row = mc[(mc["method"] == "MajorityVote") & (mc["modality"] == mod)]
        if not mv_row.empty:
            mv_acc = pd.to_numeric(mv_row["accuracy_among_callable"].iloc[0], errors="coerce")
            mv_cr  = pd.to_numeric(mv_row["callable_rate"].iloc[0], errors="coerce")
            ax.scatter([mv_cr], [mv_acc], s=120, color="#d97706", zorder=5,
                       marker="D", edgecolors="white", linewidths=1)
            ax.annotate(f"  MajorityVote\n  ({mv_cr:.2f}, {mv_acc:.2f})",
                        xy=(mv_cr, mv_acc), fontsize=7.5, color="#d97706")

        # Best single tool reference
        summ_mod = summ[summ["modality"] == mod].copy()
        summ_mod["overall_correct_call_rate"] = pd.to_numeric(summ_mod["overall_correct_call_rate"], errors="coerce")
        summ_mod["accuracy_among_callable"] = pd.to_numeric(summ_mod["accuracy_among_callable"], errors="coerce")
        summ_mod["callable_rate"] = pd.to_numeric(summ_mod["callable_rate"], errors="coerce")
        if not summ_mod.empty:
            best_row = summ_mod.loc[summ_mod["overall_correct_call_rate"].idxmax()]
            best_tool = best_row["tool"]
            best_acc  = float(best_row["accuracy_among_callable"])
            best_cr   = float(best_row["callable_rate"])
            ax.scatter([best_cr], [best_acc], s=120, marker="^", zorder=5,
                       color=TOOL_COLORS.get(best_tool, "#06b6d4"), edgecolors="white", linewidths=1)
            ax.annotate(f"  {best_tool}\n  ({best_cr:.2f}, {best_acc:.2f})",
                        xy=(best_cr, best_acc), fontsize=7.5, color=TOOL_COLORS.get(best_tool, "#06b6d4"))

        ax.set_xlabel("Call rate (1 − abstention rate)", fontsize=10)
        ax.set_ylabel("Accuracy among called loci", fontsize=10)
        ax.set_title(MODALITY_LABELS[mod], fontsize=11, fontweight="bold")
        ax.set_xlim(-0.02, 1.12)
        ax.set_ylim(-0.02, 1.12)
        if mod == "wgs":
            ax.legend(fontsize=8, loc="upper left")
        ax.annotate("← more selective", xy=(0.02, 0.04), xycoords="axes fraction", fontsize=8, color="#6b7280")
        ax.annotate("more calls →", xy=(0.72, 0.04), xycoords="axes fraction", fontsize=8, color="#6b7280")

    fig.suptitle("Per-modality abstention efficiency: coverage vs. accuracy tradeoff", fontsize=13, fontweight="bold")
    fig.tight_layout()
    watermark(fig)
    save_fig(fig, out_dir, "figure_E_abstention_efficiency")
    print("  Figure E saved.")

# ---------------------------------------------------------------------------
# Table A — Comprehensive accuracy with Wilson CIs
# ---------------------------------------------------------------------------

def table_A_comprehensive_accuracy(tables: dict, out_dir: Path):
    df = tables["summary"].copy()
    if df.empty:
        return

    for col in ["overall_correct_call_rate", "accuracy_among_callable", "callable_rate",
                "gene_rows", "sample_count", "median_runtime_hours", "median_max_ram_gb"]:
        df[col] = pd.to_numeric(df.get(col, pd.Series(dtype=float)), errors="coerce")

    rows = []
    for _, r in df.iterrows():
        tool = str(r.get("tool", ""))
        mod  = str(r.get("modality", ""))
        n_loci = int(r.get("gene_rows", 0) or 0)
        ocr    = float(r.get("overall_correct_call_rate", 0) or 0)
        aac    = float(r.get("accuracy_among_callable", 0) or 0)
        cr     = float(r.get("callable_rate", 0) or 0)
        n_correct = round(ocr * n_loci)
        ci_lo, ci_hi = wilson_ci(n_correct, n_loci)
        note = TOOL_NOTES.get((tool, mod), "")
        rows.append({
            "tool": tool, "modality": mod,
            "n_samples": int(r.get("sample_count", 0) or 0),
            "n_loci": n_loci,
            "callable_rate": f"{cr:.4f}",
            "accuracy_among_callable": f"{aac:.4f}",
            "overall_correct_call_rate": f"{ocr:.4f}",
            "ci_lo_95": f"{ci_lo:.4f}",
            "ci_hi_95": f"{ci_hi:.4f}",
            "median_runtime_h": f"{r.get('median_runtime_hours', ''):.3f}" if pd.notna(r.get("median_runtime_hours")) else "",
            "notes": note,
        })

    out_dir.mkdir(parents=True, exist_ok=True)
    fieldnames = ["tool","modality","n_samples","n_loci","callable_rate",
                  "accuracy_among_callable","overall_correct_call_rate",
                  "ci_lo_95","ci_hi_95","median_runtime_h","notes"]
    with open(out_dir / "table_A_comprehensive_accuracy.tsv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    print("  Table A saved.")

# ---------------------------------------------------------------------------
# Table B — Resolution-stratified performance
# ---------------------------------------------------------------------------

def table_B_resolution_stratified(tables: dict, out_dir: Path):
    df = tables["ambiguity"].copy()
    if df.empty:
        return

    for col in ["exact_2field_rate","exact_3field_rate","g_group_match_rate","p_group_match_rate",
                "g_group_comparable_rows","p_group_comparable_rows","gene_rows"]:
        df[col] = pd.to_numeric(df.get(col, pd.Series(dtype=float)), errors="coerce").fillna(0)

    rows = []
    for _, r in df.iterrows():
        tool = str(r["tool"])
        mod  = str(r["modality"])
        v2   = float(r["exact_2field_rate"])
        v3   = float(r["exact_3field_rate"])
        flag = "3-field capable" if v3 >= v2 * 0.8 and v3 > 0.1 else ""
        rows.append({
            "tool": tool, "modality": mod,
            "n_loci": int(r.get("gene_rows", 0)),
            "exact_2field_rate": f"{v2:.4f}",
            "exact_3field_rate": f"{v3:.4f}",
            "g_group_rows": int(r.get("g_group_comparable_rows", 0)),
            "g_group_match_rate": f"{float(r.get('g_group_match_rate', 0)):.4f}",
            "p_group_rows": int(r.get("p_group_comparable_rows", 0)),
            "p_group_match_rate": f"{float(r.get('p_group_match_rate', 0)):.4f}",
            "notes": flag,
        })

    out_dir.mkdir(parents=True, exist_ok=True)
    fieldnames = ["tool","modality","n_loci","exact_2field_rate","exact_3field_rate",
                  "g_group_rows","g_group_match_rate","p_group_rows","p_group_match_rate","notes"]
    with open(out_dir / "table_B_resolution_stratified.tsv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    print("  Table B saved.")

# ---------------------------------------------------------------------------
# Figure F — Graphical accuracy table with inline bars
# ---------------------------------------------------------------------------

def fig_F_accuracy_table_bars(tables: dict, out_dir: Path):
    summ = tables["summary"].copy()
    mc   = tables["method_comp"].copy()
    if summ.empty:
        return

    for col in ["overall_correct_call_rate", "callable_rate",
                "overall_correct_call_rate_ci_lo", "overall_correct_call_rate_ci_hi"]:
        for df in [summ, mc]:
            if col in df.columns:
                df[col] = pd.to_numeric(df[col], errors="coerce")

    mods = ["wgs", "wes", "rnaseq"]
    mod_labels = {"wgs": "WGS", "wes": "WES", "rnaseq": "RNA-seq"}

    # Tool order: single tools sorted by WGS accuracy, then ensemble at bottom
    single_order = ["OptiType","T1K","HLA-HD","POLYSOLVER","SpecHLA","Kourami","Seq2HLA","ArcasHLA"]
    ensemble_order = ["MajorityVote","WeightedConsensus"]

    def color_for(acc):
        if acc >= 0.90: return "#16a34a"
        if acc >= 0.70: return "#d97706"
        return "#dc2626"

    fig, axes = plt.subplots(1, 3, figsize=(16, 7))

    for ax, mod in zip(axes, mods):
        sub_single = summ[summ["modality"] == mod].copy()
        sub_ens    = mc[(mc["modality"] == mod) & (mc["method_type"].isin(["baseline","ensemble"]))].copy()

        # Build rows: (label, accuracy, ci_lo, ci_hi, color, is_ensemble)
        rows_data = []
        for tool in single_order:
            r = sub_single[sub_single["tool"] == tool]
            if r.empty: continue
            acc  = float(r["overall_correct_call_rate"].iloc[0])
            n    = int(r.get("gene_rows", r.get("sample_count", pd.Series([0]))).iloc[0])
            ocr_n = round(acc * n) if n > 0 else 0
            ci_lo, ci_hi = wilson_ci(ocr_n, n)
            cr   = float(r["callable_rate"].iloc[0])
            label = f"{tool}" + (f"\n  (callable {cr:.0%})" if cr < 0.99 else "")
            rows_data.append((label, acc, ci_lo, ci_hi, color_for(acc), False))

        # Separator
        rows_data.append(("", 0, 0, 0, "white", False))

        for method in ensemble_order:
            r = sub_ens[sub_ens["method"] == method]
            if r.empty: continue
            acc  = float(r["overall_correct_call_rate"].iloc[0])
            cr   = float(r["callable_rate"].iloc[0])
            n    = int(r["gene_rows"].iloc[0]) if "gene_rows" in r.columns else 0
            ocr_n = round(acc * n) if n > 0 else 0
            ci_lo, ci_hi = wilson_ci(ocr_n, n)
            label = f"{method}" + (f"\n  (callable {cr:.0%})" if cr < 0.99 else "")
            rows_data.append((label, acc, ci_lo, ci_hi, "#0f766e", True))

        y_pos = list(range(len(rows_data)))
        y_pos.reverse()  # top to bottom

        for idx, (label, acc, ci_lo, ci_hi, col, is_ens) in enumerate(rows_data):
            y = y_pos[idx]
            if not label:  # separator
                ax.axhline(y + 0.5, color="#e5e7eb", lw=1.5, ls="--")
                continue
            # Bar
            lw = 2.0 if is_ens else 1.0
            ax.barh(y, acc, height=0.65, color=col, alpha=0.85 if not is_ens else 1.0,
                    edgecolor="white", linewidth=0.3)
            # CI whisker
            if ci_hi > ci_lo:
                ax.plot([ci_lo, ci_hi], [y, y], color="#374151", lw=1.5, solid_capstyle="round")
                ax.plot([ci_lo, ci_lo], [y-0.15, y+0.15], color="#374151", lw=1.5)
                ax.plot([ci_hi, ci_hi], [y-0.15, y+0.15], color="#374151", lw=1.5)
            # Accuracy text inside/outside bar
            txt_x = max(acc - 0.03, 0.02)
            ax.text(txt_x, y, f"{acc:.3f}", va="center", ha="right",
                    fontsize=8, fontweight="bold" if is_ens else "normal",
                    color="white" if acc > 0.15 else "#374151")

        # Y-axis labels
        ax.set_yticks(y_pos)
        ax.set_yticklabels([r[0] for r in rows_data], fontsize=8.5)
        ax.set_xlim(0, 1.08)
        ax.set_xlabel("Overall correct-call rate", fontsize=9)
        ax.set_title(f"{mod_labels[mod]}\n(n={int(sub_single['sample_count'].max()) if not sub_single.empty else 0})",
                     fontsize=11, fontweight="bold")
        ax.axvline(x=0, color="#374151", lw=0.5)

        # Reference line at best single-tool accuracy
        if not sub_single.empty:
            best = float(sub_single["overall_correct_call_rate"].max())
            ax.axvline(x=best, color="#3b82f6", lw=1, ls=":", alpha=0.7)

        # Color legend on first panel
        if mod == "wgs":
            from matplotlib.patches import Patch
            patches = [Patch(color="#16a34a", label="≥ 90%"),
                       Patch(color="#d97706", label="70–90%"),
                       Patch(color="#dc2626", label="< 70%"),
                       Patch(color="#0f766e", label="Ensemble")]
            ax.legend(handles=patches, fontsize=8, loc="lower right", framealpha=0.9)

    fig.suptitle("HLA typing accuracy across tools, modalities, and ensemble methods\n"
                 "(bars = overall correct-call rate; whiskers = 95% Wilson CI; blue dotted = best single tool)",
                 fontsize=12, fontweight="bold")
    fig.tight_layout()
    watermark(fig)
    save_fig(fig, out_dir, "figure_F_accuracy_table_bars")
    print("  Figure F saved.")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def parse_args():
    p = argparse.ArgumentParser(description="Generate new publication figures for PIHLA.")
    p.add_argument("--tables-dir", required=True, help="Path to benchmark tables directory.")
    p.add_argument("--figures-dir", required=True, help="Output directory for figures and tables.")
    return p.parse_args()


def main():
    args = parse_args()
    tables_dir = Path(args.tables_dir)
    out_dir    = Path(args.figures_dir)

    apply_style()
    print(f"Reading tables from: {tables_dir}")
    print(f"Writing figures to:  {out_dir}")
    print()

    tables = load_tables(tables_dir)

    print("Generating Figure A (accuracy landscape)...")
    fig_A_accuracy_landscape(tables, out_dir)

    print("Generating Figure B (calibration reliability)...")
    fig_B_calibration_reliability(tables, out_dir)

    print("Generating Figure C (resolution profile)...")
    fig_C_resolution_profile(tables, out_dir)

    print("Generating Figure D (voting dynamics)...")
    fig_D_voting_dynamics(tables, out_dir)

    print("Generating Figure E (abstention efficiency)...")
    fig_E_abstention_efficiency(tables, out_dir)

    print("Generating Table A (comprehensive accuracy)...")
    table_A_comprehensive_accuracy(tables, out_dir)

    print("Generating Table B (resolution stratified)...")
    table_B_resolution_stratified(tables, out_dir)

    print("Generating Figure F (graphical accuracy table with bars)...")
    fig_F_accuracy_table_bars(tables, out_dir)

    print()
    print("─" * 50)
    print(f"All outputs written to {out_dir}/")


if __name__ == "__main__":
    main()
