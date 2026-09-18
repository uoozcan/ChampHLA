#!/usr/bin/env python3
"""
generate_figures_analysis.py  —  MVHLA supplementary analysis figures
Produces five new publication figures addressing critical peer-review weaknesses.

Figures generated:
  figure_a1_guardrail_scatter.{png,pdf}     W3/W2: Calibration guardrail decision map
  figure_a2_tool_ceiling.{png,pdf}          W1:    Tool-limitation ceiling bar chart
  figure_a3_weight_decomposition.{png,pdf}  W1:    Weight decomposition by tool/modality
  figure_a8_calibration_ci.{png,pdf}        W6:    Reliability diagram with CI bands
  figure_a9_weight_sensitivity.{png,pdf}    W4:    Weight sensitivity (alpha sweep)

Usage:
    python3.11 bin/generate_figures_analysis.py [--outdir PATH] [--verify]
"""

import argparse
import csv
import json
import math
import sys
from pathlib import Path

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches
    from matplotlib.lines import Line2D
    import numpy as np
except ImportError:
    sys.exit("matplotlib/numpy not installed. Run: python3.11 -m pip install --user matplotlib numpy")

# ─── Paths ────────────────────────────────────────────────────────────────────

_SCRATCH = Path("/scratch/project_2008084/pihla-publish")

BENCH = {
    "wgs":    _SCRATCH / "analysis/benchmark_wgs_all_samples/tables",
    "wes":    _SCRATCH / "analysis/benchmark_wes_all_samples/tables",
    "rna":    _SCRATCH / "analysis/benchmark_rna_all_samples/tables",
    "trimodal": _SCRATCH / "analysis/benchmark_trimodal_all_samples/tables",
}
WEIGHT_SENS = _SCRATCH / "analysis/weight_sensitivity"
DEFAULT_OUTDIR = _SCRATCH / "analysis/figures_final"

# ─── Design tokens ────────────────────────────────────────────────────────────
# Wong colorblind-safe palette
C_WGS  = "#E69F00"   # amber
C_WES  = "#0072B2"   # dark blue
C_RNA  = "#009E73"   # teal-green
C_GREY = "#999999"
C_DARK = "#333333"
C_CEILING = "#D55E00"  # vermillion for ceiling line

MODALITY_COLOR = {"wgs": C_WGS, "wes": C_WES, "rnaseq": C_RNA}
MODALITY_LABEL = {"wgs": "WGS", "wes": "WES", "rnaseq": "RNA-seq"}
MODALITY_ORDER = ["wgs", "wes", "rnaseq"]
# Map canonical modality names to benchmark directory keys
MODALITY_BENCH = {"wgs": "wgs", "wes": "wes", "rnaseq": "rna"}

METHOD_COLOR = {
    "single_tool": "#BBBBBB",
    "baseline":    C_WES,
    "ensemble":    C_RNA,
}

DPI = 300
FONT_MAIN  = 9
FONT_SMALL = 7.5
FONT_TITLE = 10

plt.rcParams.update({
    "font.size":        FONT_MAIN,
    "axes.labelsize":   FONT_MAIN,
    "axes.titlesize":   FONT_TITLE,
    "xtick.labelsize":  FONT_SMALL,
    "ytick.labelsize":  FONT_SMALL,
    "legend.fontsize":  FONT_SMALL,
    "figure.dpi":       DPI,
    "axes.spines.top":  False,
    "axes.spines.right":False,
    "figure.facecolor": "white",
})


# ─── Helpers ──────────────────────────────────────────────────────────────────

def read_tsv(path):
    with open(path) as f:
        return list(csv.DictReader(f, delimiter="\t"))


def f(x, decimals=4):
    """Coerce to float, return None if empty/missing."""
    if x is None or str(x).strip() == "":
        return None
    try:
        return round(float(x), decimals)
    except ValueError:
        return None


def wilson_ci(k, n, z=1.96):
    """Wilson binomial confidence interval. Returns (lo, hi)."""
    if n == 0:
        return (0.0, 1.0)
    p = k / n
    denom = 1 + z**2 / n
    centre = (p + z**2 / (2 * n)) / denom
    half = (z / denom) * math.sqrt(p * (1 - p) / n + z**2 / (4 * n**2))
    return (max(0.0, centre - half), min(1.0, centre + half))


def save_figure(fig, outdir, stem):
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    for ext in ("png", "pdf"):
        path = outdir / f"{stem}.{ext}"
        fig.savefig(path, dpi=DPI, bbox_inches="tight")
    print(f"  Saved {stem}.png/pdf")


# ─── Figure A1: Calibration guardrail scatter ─────────────────────────────────

def fig_a1_guardrail_scatter(outdir):
    """
    ECE vs Brier score for every (tool, modality) combination.
    Filled markers = blocked by guardrail; open = retained.
    Decision boundary lines at 0.35.
    """
    rows = []
    for mod in MODALITY_ORDER:
        path = BENCH[MODALITY_BENCH[mod]] / "tool_confidence_weights.tsv"
        if not path.exists():
            continue
        for r in read_tsv(path):
            brier = f(r.get("brier_score"))
            ece   = f(r.get("expected_calibration_error"))
            if brier is None or ece is None:
                # no_confidence tools — plot at edge
                status = r.get("guardrail_status", "")
                if status == "no_confidence":
                    rows.append(dict(tool=r["tool"], modality=r["modality"],
                                     brier=None, ece=None,
                                     status="no_confidence",
                                     final_weight=f(r.get("final_weight"))))
                continue
            rows.append(dict(
                tool=r["tool"], modality=r["modality"],
                brier=brier, ece=ece,
                status=r.get("guardrail_status", ""),
                final_weight=f(r.get("final_weight"))
            ))

    fig, ax = plt.subplots(figsize=(6.5, 5.5))

    THRESHOLD = 0.35
    ax.axhline(THRESHOLD, color=C_CEILING, lw=1.2, ls="--", alpha=0.8, zorder=1)
    ax.axvline(THRESHOLD, color=C_CEILING, lw=1.2, ls="--", alpha=0.8, zorder=1)
    ax.fill_between([0, THRESHOLD], [THRESHOLD, THRESHOLD], [0, 0],
                    color="#009E7320", zorder=0)
    ax.text(0.02, 0.02, "Retained\n(guardrail passed)",
            transform=ax.transAxes, fontsize=FONT_SMALL,
            color="#009E73", va="bottom")
    ax.fill_betweenx([THRESHOLD, 1.1], [0, 0], [THRESHOLD, THRESHOLD],
                     color="#D55E0015", zorder=0)
    ax.text(0.98, 0.98, "Blocked\n(poor calibration)",
            transform=ax.transAxes, fontsize=FONT_SMALL,
            color=C_CEILING, va="top", ha="right")

    marker_by_status = {
        "applied":          ("o",  True),
        "poor_calibration": ("o",  False),
        "no_confidence":    ("x",  False),
        "coverage_too_low": ("^",  False),
    }

    jitter = 0.008
    np.random.seed(42)
    label_placed = {}

    for row in rows:
        if row["brier"] is None:
            continue
        mod    = row["modality"]
        status = row["status"]
        color  = MODALITY_COLOR.get(mod, C_GREY)
        shape, filled = marker_by_status.get(status, ("o", False))
        jx = row["ece"]   + np.random.uniform(-jitter, jitter)
        jy = row["brier"] + np.random.uniform(-jitter, jitter)

        ax.scatter(jx, jy,
                   c=color if filled else "none",
                   edgecolors=color,
                   linewidths=1.4,
                   marker=shape, s=80, zorder=4)

        # Label each point with tool name (abbreviated)
        label_key = (row["tool"], mod)
        if label_key not in label_placed:
            tool_short = row["tool"].replace("POLYSOLVER", "POLY").replace("ArcasHLA", "ArcasHLA")
            ax.annotate(tool_short,
                        xy=(jx, jy),
                        xytext=(4, 3), textcoords="offset points",
                        fontsize=6.5, color=color, zorder=5)
            label_placed[label_key] = True

    # Legend
    legend_handles = [
        Line2D([0],[0], marker="o", color="w", markerfacecolor=MODALITY_COLOR["wgs"],
               markeredgecolor=MODALITY_COLOR["wgs"], ms=7, label="WGS"),
        Line2D([0],[0], marker="o", color="w", markerfacecolor=MODALITY_COLOR["wes"],
               markeredgecolor=MODALITY_COLOR["wes"], ms=7, label="WES"),
        Line2D([0],[0], marker="o", color="w", markerfacecolor=MODALITY_COLOR["rnaseq"],
               markeredgecolor=MODALITY_COLOR["rnaseq"], ms=7, label="RNA-seq"),
        Line2D([0],[0], linestyle="none"),
        Line2D([0],[0], marker="o", color="w", markerfacecolor="#666666",
               markeredgecolor="#666666", ms=7, label="Retained (filled)"),
        Line2D([0],[0], marker="o", color="w", markerfacecolor="none",
               markeredgecolor="#666666", ms=7, label="Blocked (open)"),
        Line2D([0],[0], marker="x", color="#666666", linestyle="none",
               ms=7, label="No confidence data"),
        Line2D([0],[0], color=C_CEILING, lw=1.2, ls="--",
               label="Guardrail threshold (0.35)"),
    ]
    ax.legend(handles=legend_handles, loc="upper left", framealpha=0.9, fontsize=FONT_SMALL)

    ax.set_xlabel("Expected Calibration Error (ECE)")
    ax.set_ylabel("Brier Score")
    ax.set_title("Fig A1 — Calibration guardrail decision map\n"
                 "Tools above/right of threshold are blocked from confidence contribution",
                 fontsize=FONT_TITLE)
    ax.set_xlim(-0.02, 1.0)
    ax.set_ylim(-0.02, 1.0)
    ax.plot([0, 1], [0, 1], color=C_GREY, lw=0.6, ls=":", zorder=0)  # diagonal

    fig.tight_layout()
    save_figure(fig, outdir, "figure_a1_guardrail_scatter")
    plt.close(fig)


# ─── Figure A2: Tool-limitation ceiling bar chart ─────────────────────────────

def fig_a2_tool_ceiling(outdir):
    """
    Horizontal bar chart of all methods per modality, with CI bars.
    Dashed vertical line marks the best single-tool accuracy (ceiling).
    """
    fig, axes = plt.subplots(1, 3, figsize=(14, 5.5), sharey=False)

    MODALITY_BENCH_KEY = {"wgs": "wgs", "wes": "wes", "rna": "rna"}
    modalities = [("wgs", "WGS"), ("wes", "WES"), ("rna", "RNA-seq")]

    for ax, (mod_key, mod_label) in zip(axes, modalities):
        bench_key = mod_key if mod_key != "rna" else "rna"
        # RNA stored as 'rnaseq' in modality column but files under benchmark_rna_all_samples
        bench_dir = BENCH["wgs"] if mod_key == "wgs" else (
            BENCH["wes"] if mod_key == "wes" else BENCH["rna"])
        rows = read_tsv(bench_dir / "method_comparison.tsv")

        # Filter to correct modality column value
        mod_filter = {"wgs": "wgs", "wes": "wes", "rna": "rnaseq"}[mod_key]
        rows = [r for r in rows if r.get("modality") == mod_filter or r.get("modality", "") == ""]

        # Find ceiling (best single tool)
        single_rates = [f(r["overall_correct_call_rate"]) for r in rows
                        if r.get("method_type") == "single_tool"
                        and f(r["overall_correct_call_rate"]) is not None]
        ceiling = max(single_rates) if single_rates else None

        # Sort by overall_correct_call_rate descending
        rows_sorted = sorted(rows,
            key=lambda r: f(r.get("overall_correct_call_rate")) or 0,
            reverse=False)  # horizontal bar starts at bottom

        labels, vals, lo_errs, hi_errs, colors = [], [], [], [], []
        for r in rows_sorted:
            rate   = f(r.get("overall_correct_call_rate"))
            ci_lo  = f(r.get("overall_correct_call_rate_ci_lo"))
            ci_hi  = f(r.get("overall_correct_call_rate_ci_hi"))
            mtype  = r.get("method_type", "single_tool")
            if rate is None:
                continue
            labels.append(r["method"])
            vals.append(rate)
            lo_errs.append(rate - ci_lo if ci_lo is not None else 0)
            hi_errs.append(ci_hi - rate if ci_hi is not None else 0)
            colors.append(METHOD_COLOR.get(mtype, C_GREY))

        y_pos = np.arange(len(labels))
        bars = ax.barh(y_pos, vals, color=colors, height=0.65,
                       alpha=0.85, edgecolor="white", linewidth=0.5, zorder=3)
        ax.errorbar(vals, y_pos, xerr=[lo_errs, hi_errs],
                    fmt="none", color=C_DARK, capsize=3, linewidth=1.2, zorder=4)

        if ceiling is not None:
            ax.axvline(ceiling, color=C_CEILING, lw=1.5, ls="--", zorder=5,
                       label=f"Tool ceiling ({ceiling:.3f})")
            ax.text(ceiling + 0.005, len(labels) - 0.5,
                    f"Tool ceiling\n{ceiling:.3f}",
                    fontsize=FONT_SMALL, color=C_CEILING, va="top")

        ax.set_yticks(y_pos)
        ax.set_yticklabels(labels, fontsize=FONT_SMALL)
        ax.set_xlabel("Overall correct-call rate")
        ax.set_title(f"{mod_label}\n(N={rows_sorted[0].get('sample_count','?') if rows_sorted else '?'})",
                     fontsize=FONT_TITLE)
        ax.set_xlim(0, 1.05)

        # Legend for method types
        if ax == axes[0]:
            patch_handles = [
                mpatches.Patch(color=METHOD_COLOR["single_tool"], label="Single tool"),
                mpatches.Patch(color=METHOD_COLOR["baseline"],    label="MajorityVote"),
                mpatches.Patch(color=METHOD_COLOR["ensemble"],    label="WeightedConsensus"),
                Line2D([0],[0], color=C_CEILING, lw=1.5, ls="--", label="Tool ceiling"),
            ]
            ax.legend(handles=patch_handles, loc="lower right",
                      fontsize=FONT_SMALL, framealpha=0.9)

    fig.suptitle("Fig A2 — Tool-limitation ceiling: all consensus strategies vs best single tool\n"
                 "95% Wilson confidence intervals shown. Dashed line = best single-tool ceiling.",
                 fontsize=FONT_TITLE, y=1.02)
    fig.tight_layout()
    save_figure(fig, outdir, "figure_a2_tool_ceiling")
    plt.close(fig)


# ─── Figure A3: Weight decomposition stacked bars ─────────────────────────────

def fig_a3_weight_decomposition(outdir):
    """
    Stacked horizontal bars per (tool, modality) showing reliability vs confidence components.
    Only tools with non-zero final_weight are shown.
    """
    all_rows = []
    for mod in MODALITY_ORDER:
        path = BENCH[MODALITY_BENCH[mod]] / "tool_confidence_weights.tsv"
        if not path.exists():
            continue
        for r in read_tsv(path):
            all_rows.append(r)

    fig, axes = plt.subplots(1, 3, figsize=(14, 5), sharey=False)
    ALPHA_BASE = 0.7
    BETA_CONF  = 0.3

    for ax, mod in zip(axes, MODALITY_ORDER):
        rows = [r for r in all_rows if r["modality"] == mod]
        rows = [r for r in rows if f(r.get("final_weight")) and f(r.get("final_weight")) > 0]
        rows_sorted = sorted(rows, key=lambda r: f(r.get("final_weight")) or 0, reverse=False)

        labels, base_parts, conf_parts, bar_colors = [], [], [], []
        for r in rows_sorted:
            base = f(r.get("base_reliability")) or 0
            fw   = f(r.get("final_weight")) or 0
            status = r.get("guardrail_status", "")
            labels.append(r["tool"])

            base_component = ALPHA_BASE * base
            if status == "applied":
                eff = f(r.get("effective_confidence")) or base
                conf_component = BETA_CONF * eff
            else:
                conf_component = 0.0  # blocked or no_confidence

            base_parts.append(base_component)
            conf_parts.append(conf_component)
            bar_colors.append(MODALITY_COLOR[mod])

        y_pos = np.arange(len(labels))

        ax.barh(y_pos, base_parts, color=MODALITY_COLOR[mod],
                height=0.62, alpha=0.85, edgecolor="white",
                label=f"0.7 × base reliability", zorder=3)
        ax.barh(y_pos, conf_parts, left=base_parts,
                color="#F0A500", height=0.62, alpha=0.9, edgecolor="white",
                label=f"0.3 × effective confidence", zorder=3)

        # Hatching for blocked tools
        for i, r in enumerate(rows_sorted):
            if r.get("guardrail_status") in ("poor_calibration", "no_confidence", "coverage_too_low"):
                ax.barh(i, base_parts[i] + conf_parts[i],
                        color="none", height=0.62,
                        edgecolor=C_CEILING, linewidth=1.2,
                        hatch="///", zorder=4, alpha=0.6)

        ax.set_yticks(y_pos)
        ax.set_yticklabels(labels, fontsize=FONT_SMALL)
        ax.set_xlabel("Weight component")
        ax.set_title(f"{MODALITY_LABEL[mod]}", fontsize=FONT_TITLE)
        ax.set_xlim(0, 1.05)
        ax.axvline(0, color=C_GREY, lw=0.5)

        if ax == axes[0]:
            blocked_patch = mpatches.Patch(
                facecolor="none", edgecolor=C_CEILING, hatch="///", linewidth=1.2,
                label="Blocked by guardrail (///)  ")
            ax.legend(
                handles=[
                    mpatches.Patch(color=MODALITY_COLOR[mod], alpha=0.85,
                                   label="0.7 × base reliability"),
                    mpatches.Patch(color="#F0A500", alpha=0.9,
                                   label="0.3 × eff. confidence"),
                    blocked_patch,
                ],
                loc="lower right", fontsize=FONT_SMALL, framealpha=0.9)

    fig.suptitle("Fig A3 — Ensemble weight decomposition: reliability component vs confidence boost\n"
                 "Hatched border = tool blocked by calibration guardrail",
                 fontsize=FONT_TITLE, y=1.02)
    fig.tight_layout()
    save_figure(fig, outdir, "figure_a3_weight_decomposition")
    plt.close(fig)


# ─── Figure A8: Reliability diagram with CI bands ─────────────────────────────

def fig_a8_calibration_ci(outdir):
    """
    Reliability diagrams (mean confidence vs observed accuracy per bin) with
    Wilson 95% CI bands on observed accuracy. One panel per modality.
    """
    fig, axes = plt.subplots(1, 3, figsize=(13, 4.5))

    TOOL_COLORS = {
        "OptiType":   "#4CAF50",
        "HLA-HD":     "#2196F3",
        "ArcasHLA":   "#FF9800",
        "T1K":        "#E53935",
        "Kourami":    "#8E24AA",
        "SpecHLA":    "#00897B",
        "POLYSOLVER": "#F57F17",
        "Seq2HLA":    "#37474F",
    }

    for ax, mod in zip(axes, MODALITY_ORDER):
        path = BENCH[MODALITY_BENCH[mod]] / "confidence_bin_summary.tsv"
        if not path.exists():
            ax.set_visible(False)
            continue
        rows = read_tsv(path)

        # Perfect calibration line
        ax.plot([0, 1], [0, 1], color=C_GREY, lw=0.8, ls=":", zorder=1, label="_nolegend_")
        ax.fill_between([0, 1], [0, 0], [1, 1],
                        alpha=0.04, color=C_GREY, zorder=0)

        from collections import defaultdict
        bins_by_tool = defaultdict(list)
        for r in rows:
            tool = r["tool"]
            mc   = f(r.get("mean_confidence"))
            oa   = f(r.get("observed_accuracy"))
            n    = int(r.get("n_rows", 0) or 0)
            if mc is None or oa is None or n == 0:
                continue
            k   = round(oa * n)
            lo, hi = wilson_ci(k, n)
            bins_by_tool[tool].append((mc, oa, lo, hi, n))

        for tool, bin_data in sorted(bins_by_tool.items()):
            if len(bin_data) < 2:
                continue
            bin_data.sort(key=lambda x: x[0])
            xs  = [b[0] for b in bin_data]
            ys  = [b[1] for b in bin_data]
            los = [b[2] for b in bin_data]
            his = [b[3] for b in bin_data]
            ns  = [b[4] for b in bin_data]
            color = TOOL_COLORS.get(tool, C_GREY)

            ax.fill_between(xs, los, his, alpha=0.18, color=color, zorder=2)
            ax.plot(xs, ys, color=color, lw=1.3, marker="o",
                    ms=max(3, min(7, min(ns) / 30)), zorder=3, label=tool)

        ax.set_xlabel("Mean predicted confidence")
        ax.set_ylabel("Observed accuracy" if ax == axes[0] else "")
        ax.set_title(f"{MODALITY_LABEL[mod]}", fontsize=FONT_TITLE)
        ax.set_xlim(-0.02, 1.02)
        ax.set_ylim(-0.02, 1.02)
        ax.legend(loc="upper left", fontsize=FONT_SMALL - 0.5, framealpha=0.9)

    fig.suptitle("Fig A8 — Reliability diagrams with 95% Wilson CI bands per calibration bin\n"
                 "Shaded band = binomial uncertainty. Dotted diagonal = perfect calibration.",
                 fontsize=FONT_TITLE, y=1.02)
    fig.tight_layout()
    save_figure(fig, outdir, "figure_a8_calibration_ci")
    plt.close(fig)


# ─── Figure A9: Weight sensitivity ────────────────────────────────────────────

def fig_a9_weight_sensitivity(outdir):
    """
    Grouped bar chart: WeightedConsensus accuracy across α values per modality.
    MajorityVote reference line per modality from all_samples benchmark.
    Directly answers W4: does pure-reliability (α=1.0) also underperform MajorityVote?
    """
    sens_path = WEIGHT_SENS / "sensitivity_comparison.tsv"
    if not sens_path.exists():
        print("  SKIP figure_a9: sensitivity_comparison.tsv not found")
        return

    sens_rows = read_tsv(sens_path)

    # Majority vote baselines from all_samples benchmarks
    mv_baselines = {}
    for mod, bench_dir in [("WES", BENCH["wes"]), ("WGS", BENCH["wgs"]), ("RNA", BENCH["rna"])]:
        mod_filter = {"WES": "wes", "WGS": "wgs", "RNA": "rnaseq"}[mod]
        rows = read_tsv(bench_dir / "method_comparison.tsv")
        mv_rows = [r for r in rows if r.get("method") == "MajorityVote"
                   and (r.get("modality", mod_filter) == mod_filter)]
        if mv_rows:
            mv_baselines[mod] = f(mv_rows[0].get("overall_correct_call_rate"))

    # Build data structure: {modality: {label: rate}}
    data = {}
    for r in sens_rows:
        mod   = r.get("modality", "")
        label = r.get("label", "")
        rate  = f(r.get("overall_correct_call_rate"))
        if mod and label and rate is not None:
            data.setdefault(mod, {})[label] = rate

    modalities_present = [m for m in ["WGS", "WES", "RNA"] if m in data]
    label_order = ["confidence_only", "equal_weight", "default", "reliability_only"]
    label_display = {
        "confidence_only":  "α=0.0 (conf only)",
        "equal_weight":     "α=0.5",
        "default":          "α=0.7 (default)",
        "reliability_only": "α=1.0 (rel only)",
    }
    label_colors = {
        "confidence_only":  "#BBBBBB",
        "equal_weight":     "#888888",
        "default":          C_WES,
        "reliability_only": C_WGS,
    }

    fig, axes = plt.subplots(1, len(modalities_present),
                             figsize=(4.5 * len(modalities_present), 5), sharey=False)
    if len(modalities_present) == 1:
        axes = [axes]

    for ax, mod in zip(axes, modalities_present):
        labels_here = [l for l in label_order if l in data.get(mod, {})]
        vals   = [data[mod][l] for l in labels_here]
        colors = [label_colors[l] for l in labels_here]
        x_pos  = np.arange(len(labels_here))

        bars = ax.bar(x_pos, vals, color=colors, width=0.6,
                      edgecolor="white", alpha=0.9, zorder=3)
        for bar, val in zip(bars, vals):
            ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.008,
                    f"{val:.3f}", ha="center", va="bottom", fontsize=FONT_SMALL)

        if mod in mv_baselines and mv_baselines[mod]:
            mv_val = mv_baselines[mod]
            ax.axhline(mv_val, color=C_CEILING, lw=1.5, ls="--", zorder=4)
            ax.text(len(labels_here) - 0.5, mv_val + 0.008,
                    f"MajorityVote {mv_val:.3f}", fontsize=FONT_SMALL, color=C_CEILING, ha="right")

        ax.set_xticks(x_pos)
        ax.set_xticklabels([label_display.get(l, l) for l in labels_here],
                           rotation=25, ha="right", fontsize=FONT_SMALL)
        ax.set_ylabel("Overall correct-call rate" if ax == axes[0] else "")
        ax.set_title(f"{mod}", fontsize=FONT_TITLE)
        ax.set_ylim(0, 1.0)
        mod_color = MODALITY_COLOR.get({"WGS": "wgs", "WES": "wes", "RNA": "rnaseq"}[mod], C_GREY)
        for bar in bars:
            if bar.get_height() < (mv_baselines.get(mod) or 0):
                bar.set_edgecolor(C_CEILING)
                bar.set_linewidth(1.2)

    fig.suptitle(
        "Fig A9 — WeightedConsensus accuracy across α values vs MajorityVote baseline (dashed)\n"
        "Key finding: WES regression persists at α=1.0 (pure reliability) — weight-learning failure, not confidence issue.",
        fontsize=FONT_TITLE, y=1.03)
    fig.tight_layout()
    save_figure(fig, outdir, "figure_a9_weight_sensitivity")
    plt.close(fig)


# ─── Verify ───────────────────────────────────────────────────────────────────

def verify():
    checks = {
        "benchmark_wgs_all_samples/method_comparison.tsv":       BENCH["wgs"] / "method_comparison.tsv",
        "benchmark_wes_all_samples/method_comparison.tsv":       BENCH["wes"] / "method_comparison.tsv",
        "benchmark_rna_all_samples/method_comparison.tsv":       BENCH["rna"] / "method_comparison.tsv",
        "benchmark_wgs_all_samples/tool_confidence_weights.tsv": BENCH["wgs"] / "tool_confidence_weights.tsv",
        "benchmark_wes_all_samples/tool_confidence_weights.tsv": BENCH["wes"] / "tool_confidence_weights.tsv",
        "benchmark_rna_all_samples/tool_confidence_weights.tsv": BENCH["rna"] / "tool_confidence_weights.tsv",
        "benchmark_wgs_all_samples/confidence_bin_summary.tsv":  BENCH["wgs"] / "confidence_bin_summary.tsv",
        "weight_sensitivity/sensitivity_comparison.tsv":         WEIGHT_SENS / "sensitivity_comparison.tsv",
    }
    ok = miss = 0
    for label, path in checks.items():
        exists = path.exists()
        print(f"  {'OK  ' if exists else 'MISS'} {label}")
        if exists:
            ok += 1
        else:
            miss += 1
    print(f"\n{ok}/{ok+miss} input files found")


# ─── Orchestrator ─────────────────────────────────────────────────────────────

def build(outdir):
    print("Generating analysis figures...")
    fig_a1_guardrail_scatter(outdir)
    fig_a2_tool_ceiling(outdir)
    fig_a3_weight_decomposition(outdir)
    fig_a8_calibration_ci(outdir)
    fig_a9_weight_sensitivity(outdir)
    print(f"\nDone. Figures written to {outdir}")


# ─── CLI ──────────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    ap = argparse.ArgumentParser(description="Generate MVHLA analysis figures")
    ap.add_argument("--outdir", default=str(DEFAULT_OUTDIR))
    ap.add_argument("--verify", action="store_true")
    args = ap.parse_args()

    if args.verify:
        verify()
    else:
        build(args.outdir)
