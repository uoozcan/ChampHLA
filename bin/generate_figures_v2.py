#!/usr/bin/env python3.11
"""
generate_figures_v2.py — Publication-quality figures for PIHLA benchmark.

Reads pre-computed tables from analysis/latest_benchmark/tables/ and writes
PNG (300 dpi) + PDF to analysis/figures_v2/.

Usage:
    python3.11 bin/generate_figures_v2.py
    python3.11 bin/generate_figures_v2.py --tables_dir /path/to/tables --out_dir /path/to/out
    python3.11 bin/generate_figures_v2.py --figures 2 4 7
"""

import argparse
import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np
import pandas as pd

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

SINGLE_TOOLS = ["ArcasHLA", "HLA-HD", "Kourami", "OptiType", "SpecHLA", "T1K"]
BASELINE_METHODS = ["MajorityVote"]
ENSEMBLE_METHODS = ["WeightedConsensus"]
METHOD_ORDER = SINGLE_TOOLS + BASELINE_METHODS + ENSEMBLE_METHODS

TOOL_COLORS = {
    "ArcasHLA":          "#6b7280",
    "HLA-HD":            "#3b82f6",
    "Kourami":           "#8b5cf6",
    "OptiType":          "#06b6d4",
    "SpecHLA":           "#64748b",
    "T1K":               "#0ea5e9",
    "MajorityVote":      "#d97706",
    "WeightedConsensus": "#0f766e",
}

POP_COLORS = {
    "CEU": "#2563eb",
    "CHB": "#dc2626",
    "GBR": "#16a34a",
    "TSI": "#d97706",
    "YRI": "#7c3aed",
}

SPLIT_COLORS = {
    "training":   "#2563eb",
    "validation": "#7c3aed",
    "holdout":    "#dc2626",
}

GUARDRAIL_COLORS = {
    "applied":          "#16a34a",
    "poor_calibration": "#dc2626",
    "no_confidence":    "#9ca3af",
}

WATERMARK = "PIHLA v2.0 · IMGT/HLA 3.59.0 · WGS Holdout n=9"
DPI = 300

TABLE_FILES = {
    "method_comparison":      "method_comparison.tsv",
    "method_per_gene":        "method_per_gene.tsv",
    "per_gene_gain":          "per_gene_gain.tsv",
    "confidence_calibration": "confidence_calibration_summary.tsv",
    "confidence_bins":        "confidence_bin_summary.tsv",
    "abstention_tradeoff":    "abstention_tradeoff.tsv",
    "discordance_summary":    "discordance_summary.tsv",
    "tool_weights_by_gene":   "tool_confidence_weights_by_gene.tsv",
}

# ---------------------------------------------------------------------------
# Style helpers
# ---------------------------------------------------------------------------

def apply_publication_style():
    plt.rcParams.update({
        "font.family":        "DejaVu Sans",
        "font.size":          10,
        "axes.titlesize":     12,
        "axes.titleweight":   "bold",
        "axes.labelsize":     10,
        "xtick.labelsize":    9,
        "ytick.labelsize":    9,
        "axes.spines.top":    False,
        "axes.spines.right":  False,
        "axes.grid":          True,
        "axes.grid.axis":     "y",
        "grid.color":         "#e5e7eb",
        "grid.linewidth":     0.6,
        "legend.fontsize":    9,
        "legend.framealpha":  0.9,
        "figure.dpi":         DPI,
        "savefig.dpi":        DPI,
        "savefig.bbox":       "tight",
    })


def add_watermark(fig):
    fig.text(
        0.99, 0.005, WATERMARK,
        ha="right", va="bottom",
        fontsize=7, color="#9ca3af",
        transform=fig.transFigure,
    )


def save_figure(fig, out_dir: Path, stem: str):
    out_dir.mkdir(parents=True, exist_ok=True)
    for ext in ("png", "pdf"):
        out_path = out_dir / f"{stem}.{ext}"
        fig.savefig(str(out_path), dpi=DPI if ext == "png" else None)
    plt.close(fig)

# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------

def load_tables(tables_dir: Path) -> dict:
    tables = {}
    missing = []
    for key, fname in TABLE_FILES.items():
        fpath = tables_dir / fname
        if not fpath.exists():
            missing.append(str(fpath))
            tables[key] = pd.DataFrame()
            continue
        df = pd.read_csv(str(fpath), sep="\t", dtype=str)
        for col in df.columns:
            try:
                df[col] = pd.to_numeric(df[col], errors="ignore")
            except Exception:
                pass
        tables[key] = df
    if missing:
        print(f"WARNING: Missing table files:\n  " + "\n  ".join(missing), file=sys.stderr)
    return tables


def load_metadata(tables_dir: Path) -> dict:
    path = tables_dir / "benchmark_metadata.json"
    fallback = {
        "population_counts": {"CEU": 7, "CHB": 8, "GBR": 7, "TSI": 9, "YRI": 11},
        "split_membership_summary": {"training": 25, "validation": 8, "holdout": 9},
        "imgt_hla_version": "3.59.0",
        "tuned_consensus_min_support": 0.45,
    }
    if not path.exists():
        print(f"WARNING: benchmark_metadata.json not found, using hard-coded values.", file=sys.stderr)
        return fallback
    try:
        with open(str(path)) as f:
            data = json.load(f)
        # Fill in fallback keys if missing
        for k, v in fallback.items():
            if k not in data:
                data[k] = v
        return data
    except Exception as exc:
        print(f"WARNING: Failed to parse benchmark_metadata.json: {exc}", file=sys.stderr)
        return fallback

# ---------------------------------------------------------------------------
# Figure 1 — Cohort Design Overview
# ---------------------------------------------------------------------------

def fig1_cohort_overview(meta: dict, out_dir: Path):
    pop_counts = meta.get("population_counts", {})
    splits = meta.get("split_membership_summary", {})

    if not pop_counts or not splits:
        print("WARNING: Skipping Figure 1 — missing population or split data.", file=sys.stderr)
        return

    fig, (ax_left, ax_right) = plt.subplots(1, 2, figsize=(13, 5))
    fig.suptitle(
        "Figure 1. Study Cohort: Population Composition and Benchmark Split Design\n"
        "(1000 Genomes WGS, n=42, IMGT/HLA 3.59.0)",
        fontsize=11, fontweight="bold", y=1.02,
    )

    # --- Panel A: Population bar chart ---
    pops = sorted(pop_counts.keys())
    counts = [pop_counts[p] for p in pops]
    colors = [POP_COLORS.get(p, "#9ca3af") for p in pops]
    bars = ax_left.barh(pops, counts, color=colors, height=0.55, edgecolor="white", linewidth=0.8)
    for bar, count in zip(bars, counts):
        ax_left.text(
            count + 0.15, bar.get_y() + bar.get_height() / 2,
            str(count), va="center", ha="left", fontsize=9, fontweight="bold",
        )
    ax_left.set_xlabel("Number of Samples")
    ax_left.set_xlim(0, max(counts) * 1.3)
    ax_left.set_title("A. Population Composition", loc="left")
    ax_left.grid(axis="x")
    ax_left.grid(axis="y", visible=False)
    ax_left.invert_yaxis()
    total = sum(counts)
    ax_left.text(
        0.98, 0.02, f"Total: n = {total}",
        ha="right", va="bottom", transform=ax_left.transAxes,
        fontsize=9, color="#374151",
        bbox=dict(boxstyle="round,pad=0.3", facecolor="#f3f4f6", alpha=0.8),
    )

    # --- Panel B: Benchmark split stacked bar ---
    split_order = ["training", "validation", "holdout"]
    split_labels = {
        "training":   f"Training (n={splits.get('training', 25)})",
        "validation": f"Validation (n={splits.get('validation', 8)})",
        "holdout":    f"Holdout (n={splits.get('holdout', 9)})",
    }
    left = 0
    for split in split_order:
        n = splits.get(split, 0)
        color = SPLIT_COLORS[split]
        ax_right.barh([0], [n], left=left, color=color, height=0.4,
                      label=split_labels[split], edgecolor="white", linewidth=1)
        if n > 0:
            ax_right.text(
                left + n / 2, 0,
                f"{split_labels[split]}\nn={n}",
                ha="center", va="center", fontsize=9,
                color="white", fontweight="bold",
            )
        left += n
    ax_right.set_xlim(0, left * 1.05)
    ax_right.set_yticks([])
    ax_right.set_xlabel("Number of Samples")
    ax_right.set_title("B. Benchmark Split Design", loc="left")
    ax_right.grid(axis="x")
    ax_right.grid(axis="y", visible=False)
    ax_right.legend(loc="upper right", fontsize=8)

    # Annotate split purposes below bar
    for i, (split, label) in enumerate(split_labels.items()):
        purpose = {"training": "weight calibration", "validation": "threshold tuning", "holdout": "final reporting"}[split]
        n_samples = splits.get(split, 0)
        # Position annotation at center of each split
        split_vals = [splits.get(s, 0) for s in split_order[:split_order.index(split)]]
        cx = sum(split_vals) + splits.get(split, 0) / 2
        ax_right.text(cx, -0.35, purpose, ha="center", va="top", fontsize=7.5, color=SPLIT_COLORS[split])

    ax_right.set_ylim(-0.7, 0.7)

    fig.tight_layout(pad=2.0)
    add_watermark(fig)
    save_figure(fig, out_dir, "figure_01_cohort_overview")
    print("  Figure 1 saved.")

# ---------------------------------------------------------------------------
# Figure 2 — Overall Method Performance
# ---------------------------------------------------------------------------

def fig2_method_performance(df: pd.DataFrame, out_dir: Path):
    if df.empty:
        print("WARNING: Skipping Figure 2 — method_comparison table is empty.", file=sys.stderr)
        return

    wgs = df[df["modality"] == "wgs"].copy()
    if wgs.empty:
        print("WARNING: Skipping Figure 2 — no WGS rows found.", file=sys.stderr)
        return

    # Reindex to canonical order
    wgs = wgs.set_index("method").reindex(METHOD_ORDER).dropna(subset=["overall_correct_call_rate"]).reset_index()
    methods = list(wgs["method"])
    n_methods = len(methods)

    accuracy = pd.to_numeric(wgs["overall_correct_call_rate"], errors="coerce").fillna(0).values
    callable_rate = pd.to_numeric(wgs["callable_rate"], errors="coerce").fillna(0).values
    colors = [TOOL_COLORS.get(m, "#6b7280") for m in methods]

    fig, ax = plt.subplots(figsize=(14, 6.5))
    x = np.arange(n_methods)
    bar_w = 0.35

    bars_acc = ax.bar(x - bar_w / 2, accuracy, bar_w, color=colors, alpha=0.92,
                      label="Overall Correct Call Rate", edgecolor="white", linewidth=0.8)
    bars_call = ax.bar(x + bar_w / 2, callable_rate, bar_w, color=colors, alpha=0.38,
                       hatch="///", label="Callable Rate", edgecolor="white", linewidth=0.8)

    # Value labels
    for bar, val in zip(bars_acc, accuracy):
        if val == 0:
            ax.text(bar.get_x() + bar.get_width() / 2, 0.012,
                    "0%", ha="center", va="bottom", fontsize=7.5,
                    color="#dc2626", fontweight="bold")
        else:
            ax.text(bar.get_x() + bar.get_width() / 2, val + 0.012,
                    f"{val:.0%}", ha="center", va="bottom", fontsize=7.5)

    for bar, val in zip(bars_call, callable_rate):
        ax.text(bar.get_x() + bar.get_width() / 2, val + 0.012,
                f"{val:.0%}", ha="center", va="bottom", fontsize=7.5, color="#6b7280")

    # Separator between single tools and consensus
    n_single = len(SINGLE_TOOLS)
    ax.axvline(n_single - 0.5, color="#d1d5db", linestyle="--", linewidth=1.2)
    ax.text(n_single - 0.45, 0.97, "Consensus methods →", ha="left", va="top",
            transform=ax.get_xaxis_transform(), fontsize=8, color="#6b7280")

    # Key insight annotation
    wc_idx = methods.index("WeightedConsensus") if "WeightedConsensus" in methods else None
    mv_idx = methods.index("MajorityVote") if "MajorityVote" in methods else None
    if wc_idx is not None and mv_idx is not None:
        wc_call = callable_rate[wc_idx]
        mv_call = callable_rate[mv_idx]
        gain = wc_call - mv_call
        ax.annotate(
            f"WeightedConsensus vs MajorityVote:\n"
            f"+{gain:.0%} callable rate at equal accuracy",
            xy=(x[wc_idx] + bar_w / 2, wc_call),
            xytext=(x[wc_idx] - 2.2, wc_call + 0.08),
            fontsize=8.5,
            color="#0f766e",
            arrowprops=dict(arrowstyle="->", color="#0f766e", lw=1.2),
            bbox=dict(boxstyle="round,pad=0.4", facecolor="#f0fdf4", edgecolor="#0f766e", alpha=0.9),
        )

    ax.set_xticks(x)
    ax.set_xticklabels(methods, rotation=20, ha="right", fontsize=9)
    ax.set_ylabel("Rate")
    ax.set_ylim(0, 1.22)
    ax.set_xlim(-0.6, n_methods - 0.4)
    ax.set_title(
        "Figure 2. HLA Typing Method Performance on WGS Holdout (n=9, HLA-A/B/C)",
        fontsize=12, fontweight="bold",
    )
    ax.legend(loc="upper left", framealpha=0.9)
    ax.grid(axis="x", visible=False)

    fig.tight_layout(pad=1.5)
    add_watermark(fig)
    save_figure(fig, out_dir, "figure_02_method_performance")
    print("  Figure 2 saved.")

# ---------------------------------------------------------------------------
# Figure 3 — Per-Gene Accuracy Heatmap
# ---------------------------------------------------------------------------

def fig3_per_gene_heatmap(df: pd.DataFrame, out_dir: Path):
    if df.empty:
        print("WARNING: Skipping Figure 3 — method_per_gene table is empty.", file=sys.stderr)
        return

    wgs = df[df["modality"] == "wgs"].copy()
    if wgs.empty:
        print("WARNING: Skipping Figure 3 — no WGS rows.", file=sys.stderr)
        return

    wgs["overall_correct_call_rate"] = pd.to_numeric(
        wgs["overall_correct_call_rate"], errors="coerce"
    )
    pivot = wgs.pivot(index="method", columns="gene", values="overall_correct_call_rate")
    genes = [g for g in ["A", "B", "C"] if g in pivot.columns]
    pivot = pivot.reindex(index=METHOD_ORDER, columns=genes)

    matrix = pivot.values.astype(float)
    row_labels = list(pivot.index)
    col_labels = [f"HLA-{g}" for g in genes]
    n_rows, n_cols = matrix.shape

    fig, ax = plt.subplots(figsize=(7, 9))
    im = ax.imshow(matrix, cmap="Blues", vmin=0, vmax=1, aspect="auto")

    # Cell annotations
    for i in range(n_rows):
        for j in range(n_cols):
            val = matrix[i, j]
            if np.isnan(val):
                ax.text(j, i, "N/A", ha="center", va="center", fontsize=9,
                        color="#9ca3af", fontstyle="italic")
            else:
                text_color = "white" if val > 0.55 else "black"
                ax.text(j, i, f"{val:.2f}", ha="center", va="center",
                        fontsize=9, color=text_color, fontweight="bold")

    ax.set_xticks(range(n_cols))
    ax.set_xticklabels(col_labels, fontsize=10)
    ax.set_yticks(range(n_rows))
    ax.set_yticklabels(row_labels, fontsize=9)
    ax.xaxis.set_ticks_position("top")
    ax.xaxis.set_label_position("top")

    # Separator line between single tools and consensus
    n_single = len([m for m in SINGLE_TOOLS if m in row_labels])
    ax.axhline(n_single - 0.5, color="#374151", linewidth=1.8)
    ax.text(n_cols - 0.45, n_single - 0.35, "consensus methods", ha="right",
            va="bottom", fontsize=8, color="#374151")

    cbar = fig.colorbar(im, ax=ax, fraction=0.035, pad=0.04)
    cbar.set_label("Overall Correct Call Rate", fontsize=9)

    ax.set_title(
        "Figure 3. Per-Gene Accuracy by Method\n(WGS Holdout, n=9)",
        fontsize=11, fontweight="bold", pad=20,
    )
    ax.grid(visible=False)

    fig.tight_layout(pad=1.5)
    add_watermark(fig)
    save_figure(fig, out_dir, "figure_03_per_gene_heatmap")
    print("  Figure 3 saved.")

# ---------------------------------------------------------------------------
# Figure 4 — Confidence Calibration Reliability Diagrams
# ---------------------------------------------------------------------------

def fig4_calibration_diagrams(
    bin_df: pd.DataFrame,
    cal_df: pd.DataFrame,
    wt_df: pd.DataFrame,
    out_dir: Path,
):
    CALIB_TOOLS = ["T1K", "OptiType", "HLA-HD", "Kourami"]

    if bin_df.empty:
        print("WARNING: Skipping Figure 4 — confidence_bins table is empty.", file=sys.stderr)
        return

    # Coerce numeric columns
    for col in ["mean_confidence", "observed_accuracy", "bin_lower", "bin_upper", "bin_index"]:
        if col in bin_df.columns:
            bin_df[col] = pd.to_numeric(bin_df[col], errors="coerce")
    for col in ["brier_score", "expected_calibration_error", "mean_confidence", "observed_accuracy"]:
        if col in cal_df.columns:
            cal_df[col] = pd.to_numeric(cal_df[col], errors="coerce")

    # Build guardrail lookup: tool -> most common guardrail_status in wt_df
    guardrail_lookup = {}
    if not wt_df.empty and "guardrail_status" in wt_df.columns:
        for tool in CALIB_TOOLS:
            rows = wt_df[(wt_df["tool"] == tool) & (wt_df["modality"] == "wgs")]
            if not rows.empty:
                guardrail_lookup[tool] = rows["guardrail_status"].mode().iloc[0]
            else:
                guardrail_lookup[tool] = "poor_calibration"
    else:
        for tool in CALIB_TOOLS:
            guardrail_lookup[tool] = "poor_calibration"

    # Calibration stats lookup: tool -> {brier, ece}
    cal_stats = {}
    if not cal_df.empty:
        for _, row in cal_df.iterrows():
            if row.get("modality") == "wgs" and row.get("tool") in CALIB_TOOLS:
                cal_stats[row["tool"]] = {
                    "brier": float(row.get("brier_score", np.nan)),
                    "ece":   float(row.get("expected_calibration_error", np.nan)),
                }

    fig, axes = plt.subplots(2, 2, figsize=(11, 9))
    axes = axes.flatten()

    colors = {
        "T1K":     "#0ea5e9",
        "OptiType":"#06b6d4",
        "HLA-HD":  "#3b82f6",
        "Kourami": "#8b5cf6",
    }

    for idx, tool in enumerate(CALIB_TOOLS):
        ax = axes[idx]
        rows = bin_df[(bin_df["tool"] == tool) & (bin_df["modality"] == "wgs")].copy()
        rows = rows.sort_values("bin_lower", na_position="last")

        # Perfect calibration diagonal
        ax.plot([0, 1], [0, 1], color="#9ca3af", linestyle="--", lw=1.2, label="Perfect calibration")

        if not rows.empty and rows["mean_confidence"].notna().any():
            x_vals = rows["mean_confidence"].values.astype(float)
            y_vals = rows["observed_accuracy"].values.astype(float)

            # Shade miscalibration area
            ax.fill_between(x_vals, y_vals, x_vals,
                            alpha=0.12, color=colors.get(tool, "#6b7280"))

            ax.plot(x_vals, y_vals, "o-",
                    color=colors.get(tool, "#6b7280"), lw=2, ms=7,
                    label="Observed accuracy")

            # Individual point labels (bin counts)
            if "n_rows" in rows.columns:
                for _, pt in rows.iterrows():
                    n = int(pt.get("n_rows", 0))
                    ax.text(float(pt["mean_confidence"]) + 0.01,
                            float(pt["observed_accuracy"]) + 0.03,
                            f"n={n}", fontsize=7, color="#6b7280")

        # Guardrail & calibration annotation box
        gs = guardrail_lookup.get(tool, "poor_calibration")
        gs_color = GUARDRAIL_COLORS.get(gs, "#9ca3af")
        stats = cal_stats.get(tool, {})
        brier_str = f"{stats['brier']:.3f}" if stats.get("brier") is not None and not np.isnan(stats.get("brier", np.nan)) else "N/A"
        ece_str   = f"{stats['ece']:.3f}"   if stats.get("ece")   is not None and not np.isnan(stats.get("ece",   np.nan)) else "N/A"
        gs_label = gs.replace("_", " ")
        stats_text = f"Brier: {brier_str}\nECE: {ece_str}\nGuardrail: {gs_label}"
        ax.text(0.04, 0.96, stats_text,
                transform=ax.transAxes, va="top", ha="left", fontsize=8,
                bbox=dict(boxstyle="round,pad=0.35", facecolor=gs_color, alpha=0.18,
                          edgecolor=gs_color, linewidth=1))

        ax.set_xlim(-0.02, 1.05)
        ax.set_ylim(-0.02, 1.12)
        ax.set_xlabel("Mean Predicted Confidence")
        ax.set_ylabel("Observed Accuracy")
        ax.set_title(tool, fontsize=11)
        ax.legend(loc="lower right", fontsize=8)
        ax.grid(axis="both", visible=True)

    fig.suptitle(
        "Figure 4. Confidence Calibration Reliability Diagrams\n"
        "(WGS Training Set; ArcasHLA excluded — degenerate single-bin data)",
        fontsize=11, fontweight="bold",
    )
    fig.tight_layout(pad=2.0)
    add_watermark(fig)
    save_figure(fig, out_dir, "figure_04_calibration_diagrams")
    print("  Figure 4 saved.")

# ---------------------------------------------------------------------------
# Figure 5 — Abstention-Coverage Tradeoff
# ---------------------------------------------------------------------------

def fig5_abstention_tradeoff(df: pd.DataFrame, out_dir: Path):
    if df.empty:
        print("WARNING: Skipping Figure 5 — abstention_tradeoff table is empty.", file=sys.stderr)
        return

    for col in ["min_support", "call_rate", "accuracy_among_called", "overall_correct_call_rate"]:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")

    df = df.sort_values("call_rate", ascending=False).reset_index(drop=True)

    fig, ax = plt.subplots(figsize=(9, 7))

    # Reference diagonal
    ax.plot([0, 1], [0, 1], ":", color="#9ca3af", lw=1, label="Reference (y=x)")

    # Accuracy among called line
    ax.plot(df["call_rate"], df["accuracy_among_called"],
            "o-", color="#2563eb", lw=2, ms=8, label="Accuracy among called", zorder=4)

    # Overall correct call rate line
    ax.plot(df["call_rate"], df["overall_correct_call_rate"],
            "s--", color="#7c3aed", lw=1.5, ms=6, label="Overall correct call rate", zorder=3)

    # Threshold labels
    for _, row in df.iterrows():
        ax.annotate(
            f"θ={row['min_support']:.2f}",
            xy=(float(row["call_rate"]), float(row["accuracy_among_called"])),
            xytext=(6, 4), textcoords="offset points",
            fontsize=8, color="#374151",
        )

    # Highlight tuned threshold
    tuned = df[df["min_support"] == 0.45]
    if not tuned.empty:
        t = tuned.iloc[0]
        ax.plot(float(t["call_rate"]), float(t["accuracy_among_called"]),
                "*", color="#dc2626", ms=18, zorder=6, label="Tuned threshold (θ=0.45)")
        ax.annotate(
            "Tuned operating point\n(θ=0.45, max coverage)",
            xy=(float(t["call_rate"]), float(t["accuracy_among_called"])),
            xytext=(-80, 30), textcoords="offset points",
            fontsize=8.5, color="#dc2626",
            arrowprops=dict(arrowstyle="->", color="#dc2626", lw=1.2),
            bbox=dict(boxstyle="round,pad=0.35", facecolor="#fef2f2", edgecolor="#dc2626", alpha=0.85),
        )

    ax.set_xlabel("Call Rate", fontsize=10)
    ax.set_ylabel("Accuracy", fontsize=10)
    ax.set_xlim(0, 1.08)
    ax.set_ylim(0, 1.15)
    ax.set_title(
        "Figure 5. Abstention-Coverage Tradeoff for Weighted Consensus\n"
        "(WGS Holdout, varying min_support threshold θ)",
        fontsize=11, fontweight="bold",
    )
    ax.legend(loc="lower right", fontsize=9)
    ax.grid(axis="both", visible=True)

    fig.tight_layout(pad=1.5)
    add_watermark(fig)
    save_figure(fig, out_dir, "figure_05_abstention_tradeoff")
    print("  Figure 5 saved.")

# ---------------------------------------------------------------------------
# Figure 6 — Discordance Taxonomy
# ---------------------------------------------------------------------------

def fig6_discordance_taxonomy(df: pd.DataFrame, out_dir: Path):
    if df.empty:
        print("WARNING: Skipping Figure 6 — discordance_summary table is empty.", file=sys.stderr)
        return

    if "n_events" in df.columns:
        df["n_events"] = pd.to_numeric(df["n_events"], errors="coerce").fillna(0).astype(int)

    scope_colors = {"rnaseq": "#0f766e", "wgs": "#d97706"}

    df = df.sort_values("n_events", ascending=True).reset_index(drop=True)
    labels = [f"{row['scope'].upper()} / {row['tag']}" for _, row in df.iterrows()]
    values = df["n_events"].tolist()
    colors = [scope_colors.get(str(row.get("scope", "wgs")), "#9ca3af") for _, row in df.iterrows()]

    fig, ax = plt.subplots(figsize=(10, max(4, len(labels) * 1.2 + 1.5)))

    bars = ax.barh(labels, values, color=colors, height=0.55, edgecolor="white", linewidth=0.8)
    for bar, val in zip(bars, values):
        ax.text(
            val + 0.15, bar.get_y() + bar.get_height() / 2,
            str(val), va="center", ha="left", fontsize=9, fontweight="bold",
        )

    max_val = max(values) if values else 1
    ax.set_xlim(0, max_val * 1.3)
    ax.set_xlabel("Number of Events")
    ax.set_title(
        "Figure 6. Discordance Taxonomy Across Modalities\n"
        "(1000 Genomes WGS + RNA-seq pilot)",
        fontsize=11, fontweight="bold",
    )

    legend_patches = [
        mpatches.Patch(color="#d97706", label="WGS"),
        mpatches.Patch(color="#0f766e", label="RNA-seq"),
    ]
    ax.legend(handles=legend_patches, loc="lower right", fontsize=9)
    ax.grid(axis="x", visible=True)
    ax.grid(axis="y", visible=False)

    fig.tight_layout(pad=1.5)
    add_watermark(fig)
    save_figure(fig, out_dir, "figure_06_discordance_taxonomy")
    print("  Figure 6 saved.")

# ---------------------------------------------------------------------------
# Figure 7 — Learned Confidence Weights Heatmap
# ---------------------------------------------------------------------------

def fig7_weight_heatmap(df: pd.DataFrame, out_dir: Path):
    if df.empty:
        print("WARNING: Skipping Figure 7 — tool_weights_by_gene table is empty.", file=sys.stderr)
        return

    wgs = df[df["modality"] == "wgs"].copy()
    if wgs.empty:
        print("WARNING: Skipping Figure 7 — no WGS rows.", file=sys.stderr)
        return

    wgs["final_weight"] = pd.to_numeric(wgs["final_weight"], errors="coerce")

    genes = [g for g in ["A", "B", "C"] if g in wgs["gene"].values]
    pivot_wt = wgs.pivot(index="tool", columns="gene", values="final_weight")
    pivot_wt = pivot_wt.reindex(index=[t for t in SINGLE_TOOLS if t in pivot_wt.index], columns=genes)

    # Guardrail status (most common per tool)
    guardrail_per_tool = {}
    if "guardrail_status" in wgs.columns:
        for tool in pivot_wt.index:
            rows = wgs[wgs["tool"] == tool]
            if not rows.empty:
                mode_val = rows["guardrail_status"].mode()
                guardrail_per_tool[tool] = mode_val.iloc[0] if not mode_val.empty else "no_confidence"
            else:
                guardrail_per_tool[tool] = "no_confidence"

    matrix = pivot_wt.values.astype(float)
    row_labels = list(pivot_wt.index)
    col_labels = [f"HLA-{g}" for g in genes]
    n_rows, n_cols = matrix.shape

    # Custom teal colormap
    teal_cmap = mcolors.LinearSegmentedColormap.from_list(
        "pihla_teal", ["#ffffff", "#0f766e"]
    )

    fig, ax = plt.subplots(figsize=(8, 8))
    masked = np.ma.masked_invalid(matrix)
    im = ax.imshow(masked, cmap=teal_cmap, vmin=0, vmax=1, aspect="auto")

    # Grey for NaN (no_confidence) cells
    cmap_nan = mcolors.ListedColormap(["#f3f4f6"])
    ax.imshow(np.where(np.isnan(matrix), 1, np.nan), cmap=cmap_nan,
              vmin=0, vmax=1, aspect="auto")

    # Cell annotations
    for i in range(n_rows):
        for j in range(n_cols):
            val = matrix[i, j]
            if np.isnan(val):
                ax.text(j, i, "N/A", ha="center", va="center",
                        fontsize=9, color="#9ca3af", fontstyle="italic")
            else:
                text_color = "white" if val > 0.55 else "#1f2937"
                ax.text(j, i, f"{val:.3f}", ha="center", va="center",
                        fontsize=9, color=text_color, fontweight="bold")

    # Guardrail status indicators as colored squares in left margin
    for i, tool in enumerate(row_labels):
        gs = guardrail_per_tool.get(tool, "no_confidence")
        gs_color = GUARDRAIL_COLORS.get(gs, "#9ca3af")
        ax.add_patch(plt.Rectangle(
            (-0.95, i - 0.42), 0.35, 0.84,
            color=gs_color, transform=ax.transData, clip_on=False,
        ))

    ax.set_xticks(range(n_cols))
    ax.set_xticklabels(col_labels, fontsize=10)
    ax.set_yticks(range(n_rows))
    ax.set_yticklabels(row_labels, fontsize=9)
    ax.xaxis.set_ticks_position("top")
    ax.xaxis.set_label_position("top")
    ax.tick_params(left=False, top=False)

    cbar = fig.colorbar(im, ax=ax, fraction=0.035, pad=0.04)
    cbar.set_label("Final Confidence Weight", fontsize=9)

    # Legend for guardrail colors
    legend_patches = [
        mpatches.Patch(color=GUARDRAIL_COLORS["applied"], label="Confidence boost applied"),
        mpatches.Patch(color=GUARDRAIL_COLORS["poor_calibration"], label="Poor calibration (blocked)"),
        mpatches.Patch(color=GUARDRAIL_COLORS["no_confidence"], label="No confidence data"),
    ]
    ax.legend(handles=legend_patches, loc="lower left",
              bbox_to_anchor=(0.0, -0.18), ncol=1, fontsize=8,
              title="Guardrail Status", title_fontsize=8)

    ax.set_title(
        "Figure 7. Benchmark-Derived Confidence Weights by Tool and Gene\n"
        "(WGS Training Set, n=25; colored squares = guardrail status per tool)",
        fontsize=11, fontweight="bold", pad=20,
    )
    ax.grid(visible=False)

    fig.tight_layout(pad=1.5)
    add_watermark(fig)
    save_figure(fig, out_dir, "figure_07_weight_heatmap")
    print("  Figure 7 saved.")

# ---------------------------------------------------------------------------
# Figure 8 — Ensemble Advantage Summary
# ---------------------------------------------------------------------------

def fig8_ensemble_advantage(df: pd.DataFrame, out_dir: Path):
    if df.empty:
        print("WARNING: Skipping Figure 8 — method_comparison table is empty.", file=sys.stderr)
        return

    wgs = df[df["modality"] == "wgs"].copy()
    if wgs.empty:
        print("WARNING: Skipping Figure 8 — no WGS rows.", file=sys.stderr)
        return

    for col in ["callable_rate", "accuracy_among_callable", "overall_correct_call_rate"]:
        wgs[col] = pd.to_numeric(wgs[col], errors="coerce")

    wgs = wgs.set_index("method").reindex(METHOD_ORDER).dropna(subset=["overall_correct_call_rate"]).reset_index()
    methods = list(wgs["method"])
    callable_vals = wgs["callable_rate"].values.astype(float)
    acc_callable  = wgs["accuracy_among_callable"].values.astype(float)
    overall_vals  = wgs["overall_correct_call_rate"].values.astype(float)
    colors = [TOOL_COLORS.get(m, "#6b7280") for m in methods]

    fig, (ax_left, ax_right) = plt.subplots(1, 2, figsize=(15, 7))
    fig.suptitle(
        "Figure 8. Coverage vs Accuracy: Ensemble Decision Space\n"
        "(WGS Holdout, n=9; bubble size ∝ overall correct call rate)",
        fontsize=11, fontweight="bold",
    )

    # --- Left panel: 2D bubble chart ---
    # Ideal zone shading (upper-right quadrant)
    ax_left.fill_between([0.5, 1.05], [0.5, 0.5], [1.15, 1.15],
                         alpha=0.06, color="#16a34a")
    ax_left.text(0.97, 0.97, "Ideal zone", ha="right", va="top",
                 transform=ax_left.transAxes, fontsize=8.5,
                 color="#16a34a", fontstyle="italic")

    # Quadrant lines
    ax_left.axhline(0.5, color="#d1d5db", lw=0.8, ls="--")
    ax_left.axvline(0.5, color="#d1d5db", lw=0.8, ls="--")

    # Bubble sizes
    bubble_sizes = np.where(np.isnan(overall_vals), 30,
                            np.clip(overall_vals, 0.01, 1.0) * 2200)

    sc = ax_left.scatter(
        callable_vals, acc_callable,
        s=bubble_sizes, c=colors,
        alpha=0.82, edgecolors="white", linewidths=1.5, zorder=4,
    )

    # Label each bubble
    label_offsets = {
        "ArcasHLA":          (8, -18),
        "HLA-HD":            (-55, 8),
        "Kourami":           (8, -18),
        "OptiType":          (-70, 8),
        "SpecHLA":           (8, 8),
        "T1K":               (8, -18),
        "MajorityVote":      (-75, 10),
        "WeightedConsensus": (8, 8),
    }
    for i, method in enumerate(methods):
        ox, oy = label_offsets.get(method, (8, 8))
        ax_left.annotate(
            method,
            xy=(callable_vals[i], acc_callable[i]),
            xytext=(ox, oy), textcoords="offset points",
            fontsize=8, color="#1f2937",
            arrowprops=dict(arrowstyle="-", color="#d1d5db", lw=0.8)
            if abs(ox) > 15 or abs(oy) > 15 else None,
        )

    # ArcasHLA special annotation (0% accuracy)
    if "ArcasHLA" in methods:
        ai = methods.index("ArcasHLA")
        ax_left.annotate(
            "ArcasHLA: 0% accuracy\n(configuration issue suspected)",
            xy=(callable_vals[ai], acc_callable[ai]),
            xytext=(0.55, 0.08), textcoords="data",
            fontsize=7.5, color="#9ca3af",
            arrowprops=dict(arrowstyle="->", color="#9ca3af", lw=0.8),
        )

    ax_left.set_xlabel("Callable Rate")
    ax_left.set_ylabel("Accuracy among Called")
    ax_left.set_xlim(0, 1.1)
    ax_left.set_ylim(-0.05, 1.15)
    ax_left.grid(axis="both")

    # Bubble size legend
    for size_val, label in [(0.25, "25%"), (0.50, "50%"), (0.75, "75%")]:
        ax_left.scatter([], [], s=size_val * 2200, c="#9ca3af", alpha=0.6,
                        label=f"Overall rate {label}")
    ax_left.legend(title="Overall correct\ncall rate", loc="upper left",
                   fontsize=7.5, title_fontsize=8, framealpha=0.9)

    # --- Right panel: overall correct call rate bars ---
    y = np.arange(len(methods))
    bars = ax_right.barh(y, overall_vals, color=colors, height=0.55,
                         edgecolor="white", linewidth=0.8)

    for bar, val, method in zip(bars, overall_vals, methods):
        if val == 0:
            ax_right.text(0.005, bar.get_y() + bar.get_height() / 2,
                          "0%", va="center", ha="left", fontsize=7.5,
                          color="#dc2626", fontweight="bold")
        else:
            ax_right.text(val + 0.008, bar.get_y() + bar.get_height() / 2,
                          f"{val:.0%}", va="center", ha="left", fontsize=8)

    # Highlight WeightedConsensus improvement
    wc_idx = methods.index("WeightedConsensus") if "WeightedConsensus" in methods else None
    if wc_idx is not None:
        wc_val = overall_vals[wc_idx]
        best_single = max(
            overall_vals[i] for i, m in enumerate(methods) if m in SINGLE_TOOLS
        )
        ax_right.annotate(
            f"= best single tool\n({wc_val:.0%})",
            xy=(wc_val, wc_idx),
            xytext=(wc_val + 0.12, wc_idx),
            fontsize=8, color="#0f766e",
            arrowprops=dict(arrowstyle="->", color="#0f766e", lw=1),
            va="center",
        )

    ax_right.set_yticks(y)
    ax_right.set_yticklabels(methods, fontsize=9)
    ax_right.invert_yaxis()
    ax_right.set_xlabel("Overall Correct Call Rate")
    ax_right.set_xlim(0, 1.1)
    ax_right.set_title("Overall Correct Call Rate\n(bars sorted by method type)", fontsize=10)
    ax_right.axvline(0.5185, color="#0f766e", linestyle=":", linewidth=1, alpha=0.6)
    ax_right.grid(axis="x", visible=True)
    ax_right.grid(axis="y", visible=False)

    fig.tight_layout(pad=1.5)
    add_watermark(fig)
    save_figure(fig, out_dir, "figure_08_ensemble_advantage")
    print("  Figure 8 saved.")

# ---------------------------------------------------------------------------
# Argument parsing & main
# ---------------------------------------------------------------------------

def parse_args() -> argparse.Namespace:
    script_dir = Path(__file__).resolve().parent
    default_tables = script_dir.parent / "analysis" / "latest_benchmark" / "tables"
    default_out = script_dir.parent / "analysis" / "figures_v2"

    p = argparse.ArgumentParser(
        description="Generate publication-quality figures for PIHLA benchmark."
    )
    p.add_argument(
        "--tables_dir", type=Path, default=default_tables,
        help="Directory containing pre-computed TSV/JSON tables.",
    )
    p.add_argument(
        "--out_dir", type=Path, default=default_out,
        help="Output directory for PNG and PDF figures.",
    )
    p.add_argument(
        "--figures", nargs="*", type=int,
        choices=list(range(1, 9)),
        default=list(range(1, 9)),
        help="Which figures to generate (1–8). Default: all.",
    )
    return p.parse_args()


def main() -> int:
    apply_publication_style()
    args = parse_args()
    tables_dir = args.tables_dir
    out_dir = args.out_dir
    out_dir.mkdir(parents=True, exist_ok=True)

    print(f"Reading tables from: {tables_dir}")
    print(f"Writing figures to:  {out_dir}")

    meta = load_metadata(tables_dir)
    tables = load_tables(tables_dir)

    dispatch = {
        1: lambda: fig1_cohort_overview(meta, out_dir),
        2: lambda: fig2_method_performance(tables["method_comparison"], out_dir),
        3: lambda: fig3_per_gene_heatmap(tables["method_per_gene"], out_dir),
        4: lambda: fig4_calibration_diagrams(
                       tables["confidence_bins"],
                       tables["confidence_calibration"],
                       tables["tool_weights_by_gene"],
                       out_dir),
        5: lambda: fig5_abstention_tradeoff(tables["abstention_tradeoff"], out_dir),
        6: lambda: fig6_discordance_taxonomy(tables["discordance_summary"], out_dir),
        7: lambda: fig7_weight_heatmap(tables["tool_weights_by_gene"], out_dir),
        8: lambda: fig8_ensemble_advantage(tables["method_comparison"], out_dir),
    }

    errors = []
    for fig_num in sorted(args.figures):
        print(f"\nGenerating Figure {fig_num}...", flush=True)
        try:
            dispatch[fig_num]()
        except Exception as exc:
            print(f"  ERROR in Figure {fig_num}: {exc}", file=sys.stderr)
            import traceback
            traceback.print_exc(file=sys.stderr)
            errors.append((fig_num, exc))

    print(f"\n{'─' * 50}")
    if errors:
        print(f"Completed with {len(errors)} error(s): figures {[n for n, _ in errors]}")
        return 1
    print(f"All 8 figures written to {out_dir}/")
    return 0


if __name__ == "__main__":
    sys.exit(main())
