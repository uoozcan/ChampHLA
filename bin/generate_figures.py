#!/usr/bin/env python3.11
"""
generate_figures.py — Publication-quality figure generator for PIHLA.

Reads benchmark TSV tables produced by hla_benchmark.py and writes refreshed
benchmark figures plus optional exploratory cohort/coverage figures as
PDF + SVG + PNG (300 dpi).

Usage:
    python3.11 bin/generate_figures.py \\
        [--tables-dir PATH]   # default: analysis/1000g_realdata/tables
        [--figures-dir PATH]  # default: analysis/1000g_realdata/figures
        [--figures 1,2,3,...] # default: all

Requirements: matplotlib>=3.5, seaborn>=0.12, numpy, pandas
    python3.11 -m pip install --user matplotlib numpy pandas seaborn
"""

import argparse
import json
import sys
import warnings
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
from typing import Optional, List  # non-interactive backend

import matplotlib.patches as mpatches
import matplotlib.patheffects as pe
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
import numpy as np
import pandas as pd
import seaborn as sns

# ---------------------------------------------------------------------------
# Publication style
# ---------------------------------------------------------------------------

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Avenir Next", "Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 9,
    "axes.titlesize": 10.5,
    "axes.titleweight": "bold",
    "axes.labelsize": 9,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
    "legend.fontsize": 8,
    "legend.framealpha": 0.95,
    "legend.edgecolor": "#d8dee9",
    "figure.dpi": 150,
    "figure.facecolor": "#fcfcf8",
    "axes.facecolor": "#fcfcf8",
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.grid": True,
    "axes.grid.axis": "y",
    "grid.alpha": 0.28,
    "grid.color": "#d7dde5",
    "savefig.bbox": "tight",
    "savefig.dpi": 300,
    "savefig.pad_inches": 0.05,
})

# Colour palette
C_WES     = "#275dad"   # cobalt
C_WGS     = "#d97706"   # amber
C_RNASEQ  = "#0f766e"   # teal
C_SINGLE  = "#94a3b8"   # grey-blue (single-tool baseline)
C_MAJORITY = "#b45309"  # darker amber
C_ENSEMBLE = "#0f766e"  # teal
C_GAIN    = "#0f766e"
C_LOSS    = "#dc2626"
C_PIHLA   = "#275dad"   # brand blue
C_NEUTRAL = "#475569"
C_ACCENT  = "#e11d48"

MODALITY_COLORS = {"wes": C_WES, "wgs": C_WGS, "rnaseq": C_RNASEQ}
MODALITY_LABELS = {"wes": "WES", "wgs": "WGS", "rnaseq": "RNA-seq"}

DISCORD_COLORS = {
    "technical_conflict":    "#E05252",
    "low_evidence_conflict": "#F97316",
    "dna_rna_discordance":   "#9333EA",
    "possible_expression_bias": "#0EA5E9",
    "no_evidence":           "#94A3B8",
}

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _load(tables_dir: Path, fname: str) -> Optional[pd.DataFrame]:
    path = tables_dir / fname
    if not path.exists():
        warnings.warn(f"[SKIP] {fname} not found — skipping dependent figure(s).")
        return None
    df = pd.read_csv(path, sep="\t", dtype=str)
    if df.empty:
        warnings.warn(f"[SKIP] {fname} is empty — skipping dependent figure(s).")
        return None
    return df


def _save(fig: plt.Figure, figures_dir: Path, stem: str) -> None:
    for fmt in ("pdf", "svg", "png"):
        out = figures_dir / f"{stem}.{fmt}"
        fig.savefig(out, format=fmt)
    plt.close(fig)
    print(f"  ✓  {stem}  (pdf / svg / png)")


def _load_json(tables_dir: Path, fname: str) -> dict:
    path = tables_dir / fname
    if not path.exists():
        return {}
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def _pct(val: float) -> str:
    return f"{val * 100:.1f}"


def _infer_source_dir(tables_dir: Path, source_dir: Optional[Path]) -> Path:
    if source_dir is not None:
        return source_dir
    if tables_dir.name == "tables":
        return tables_dir.parent
    return tables_dir


def _source_table(source_dir: Path, tables_dir: Path, fname: str) -> Optional[pd.DataFrame]:
    for path in (
        source_dir / fname,
        source_dir / "manifests_v1" / fname,
        tables_dir / fname,
    ):
        if path.exists():
            df = pd.read_csv(path, sep="	", dtype=str)
            if not df.empty:
                return df
    return None


def _empty_state_figure(figures_dir: Path, stem: str, title: str, message: str) -> None:
    fig, ax = plt.subplots(figsize=(7.2, 3.8))
    ax.axis("off")
    ax.set_title(title, fontsize=10.5, fontweight="bold", pad=14)
    ax.text(0.5, 0.56, message, ha="center", va="center", fontsize=10, color=C_NEUTRAL, wrap=True)
    ax.text(0.5, 0.22, "The corresponding benchmark table was present but contained no calibration rows.", ha="center", va="center", fontsize=7.5, color="#64748b", wrap=True)
    fig.tight_layout()
    _save(fig, figures_dir, stem)


def _write_captions(figures_dir: Path, tables_dir: Path, include_exploratory: bool) -> None:
    metadata = _load_json(tables_dir, "benchmark_metadata.json")
    source_dir = _infer_source_dir(tables_dir, None)
    source_manifest = _load_json(source_dir / "metadata", "source_manifest.json")
    truth_source = metadata.get("truth_source") or metadata.get("reference", {}).get("truth_source") or "unspecified"
    supported_loci = metadata.get("supported_loci") or []
    if isinstance(supported_loci, list):
        supported_loci_text = ", ".join(supported_loci) if supported_loci else "unspecified"
    else:
        supported_loci_text = str(supported_loci)
    tuned_support = metadata.get("tuned_consensus_min_support", "unspecified")
    lines = [
        "# Figure Captions",
        "",
        "## Benchmark Figure Set",
        "",
        f"Truth source: {truth_source}",
        f"Supported loci: {supported_loci_text}",
        f"Tuned minimum support: {tuned_support}",
    ]
    if source_manifest:
        lines.extend([
            "",
            "## Authoritative Benchmark Roots",
            "",
            f"WGS authoritative source: {source_manifest.get('wgs_authoritative_root', source_manifest.get('wgs_root', 'unspecified'))}",
            f"WES authoritative source: {source_manifest.get('wes_authoritative_root', source_manifest.get('wes_root', 'unspecified'))}",
            f"RNA authoritative source: {source_manifest.get('rna_authoritative_root', source_manifest.get('rna_root', 'unspecified'))}",
            f"Trimodal secondary source: {source_manifest.get('trimodal_root', 'unspecified')}",
            "",
            "Note: figures 1-7 are generated from a staged multi-root publication input built from the authoritative modality-specific benchmarks above.",
        ])
    lines.extend([
        "",
        "### Figure 2. Method comparison",
        "Compares single-tool baselines, majority vote, weighted consensus, and routed baselines across sequencing modalities.",
        "",
        "### Figure 3. Per-gene gains",
        "Shows weighted-consensus gains over majority vote and the best single-tool baseline by gene and modality.",
        "",
        "### Figure 4. Confidence calibration",
        "Plots observed correctness against mean predicted confidence for each tool and modality.",
        "",
        "### Figure 5. Abstention tradeoff",
        "Shows how higher support thresholds trade overall call rate against accuracy among emitted calls.",
        "",
        "### Figure 6. Discordance taxonomy",
        "Summarizes the discordance categories found in the benchmark cohort.",
        "",
        "### Figure 7. Confidence-calibrated weights",
        "Visualizes learned per-tool benchmark weights across modalities.",
    ])
    if include_exploratory:
        lines.extend([
            "",
            "## Exploratory Figure Set",
            "",
            "### Figure X1. Modality coverage overview",
            "Summarizes how many samples have truth support and completed WGS, WES, and RNA-seq coverage.",
            "",
            "### Figure X2. Cohort exclusion reasons",
            "Shows the exclusion reasons preventing candidate samples from entering the final tri-modal benchmark cohort.",
            "",
            "### Figure X3. Truth locus coverage",
            "Summarizes truth-supported locus counts across the cohort and the per-sample completeness of truth annotations.",
        ])
    (figures_dir / "captions.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("  ✓  captions.md")


# ---------------------------------------------------------------------------
# Figure 1 — Workflow Architecture
# ---------------------------------------------------------------------------

def figure_1(figures_dir: Path) -> None:
    """Programmatic workflow diagram using matplotlib patches."""
    fig, ax = plt.subplots(figsize=(14, 5))
    ax.set_xlim(0, 14)
    ax.set_ylim(0, 5)
    ax.axis("off")
    ax.set_title(
        "Figure 1. PIHLA workflow and ensemble architecture",
        fontsize=12, fontweight="bold", pad=10,
    )

    def _box(ax, x, y, w, h, label, sublabel="", color="#f0f4f8", textcolor="#1e293b",
             fontsize=8.5, bold=False):
        box = FancyBboxPatch(
            (x, y), w, h,
            boxstyle="round,pad=0.06",
            facecolor=color, edgecolor="#94A3B8", linewidth=0.8,
            zorder=2,
        )
        ax.add_patch(box)
        cy = y + h / 2
        weight = "bold" if bold else "normal"
        ax.text(x + w / 2, cy + (0.13 if sublabel else 0), label,
                ha="center", va="center", fontsize=fontsize,
                color=textcolor, fontweight=weight, zorder=3)
        if sublabel:
            ax.text(x + w / 2, cy - 0.18, sublabel,
                    ha="center", va="center", fontsize=7,
                    color="#64748b", zorder=3)

    def _arrow(ax, x1, y1, x2, y2, color="#64748b"):
        ax.annotate(
            "", xy=(x2, y2), xytext=(x1, y1),
            arrowprops=dict(
                arrowstyle="-|>",
                color=color,
                lw=1.2,
                connectionstyle="arc3,rad=0.0",
            ),
            zorder=4,
        )

    # ── Column x positions ──────────────────────────────────────────────
    X_INPUT  = 0.25
    X_TOOLS  = 2.8
    X_HARM   = 5.5
    X_CAL    = 7.9
    X_ENS    = 9.9
    X_OUT    = 12.2

    # ── 1. Input block ───────────────────────────────────────────────────
    ax.text(X_INPUT + 0.85, 4.72, "Sequencing\ninputs", ha="center",
            fontsize=8, fontweight="bold", color="#1e293b")
    input_defs = [
        ("WGS", C_WGS,    3.5),
        ("WES", C_WES,    2.55),
        ("RNA-seq", C_RNASEQ, 1.6),
    ]
    for label, col, y0 in input_defs:
        _box(ax, X_INPUT, y0, 1.7, 0.75, label, color=col, textcolor="white",
             fontsize=8.5, bold=True)

    # ── 2. Tool layer ────────────────────────────────────────────────────
    ax.text(X_TOOLS + 0.85, 4.72, "HLA callers", ha="center",
            fontsize=8, fontweight="bold", color="#1e293b")
    tools = [
        ("OptiType", "WES / WGS / RNA-seq", 3.5),
        ("ArcasHLA",  "WES / WGS / RNA-seq", 2.55),
        ("SpecHLA",   "WES / WGS",           1.6),
    ]
    for name, sub, y0 in tools:
        _box(ax, X_TOOLS, y0, 1.7, 0.75, name, sublabel=sub,
             color="#dbeafe", fontsize=8)

    # ── arrows: inputs → tools ───────────────────────────────────────────
    for _, _, iy in input_defs:
        for _, _, ty in tools:
            _arrow(ax, X_INPUT + 1.7, iy + 0.375, X_TOOLS, ty + 0.375,
                   color="#cbd5e1")

    # ── 3. Harmonization ─────────────────────────────────────────────────
    _box(ax, X_HARM, 1.4, 2.1, 3.0,
         "Call harmonization &\nconfidence extraction",
         color="#fef9c3", fontsize=8.5, bold=True)
    ax.text(X_HARM + 1.05, 2.35, "• Unified HLA call schema\n"
            "• Tool-native confidence\n"
            "• Read-support sidecars\n"
            "• Allele normalization",
            fontsize=7, va="center", color="#374151")
    for _, _, ty in tools:
        _arrow(ax, X_TOOLS + 1.7, ty + 0.375, X_HARM, 2.9, color="#94a3b8")

    # ── 4. Calibration / weight learning ─────────────────────────────────
    _box(ax, X_CAL, 2.85, 1.7, 1.55,
         "Benchmark-derived\nweight learning",
         sublabel="(training split)",
         color="#dcfce7", fontsize=8, bold=True)
    ax.text(X_CAL + 0.85, 2.95,
            "per-tool / gene /\nmodality weights",
            fontsize=7, ha="center", va="top", color="#374151")
    _arrow(ax, X_HARM + 2.1, 2.9, X_CAL, 3.425)

    # ── 5. Consensus engine ───────────────────────────────────────────────
    _box(ax, X_ENS, 1.7, 2.0, 2.4,
         "Weighted ensemble\nconsensus + abstention",
         color="#ede9fe", fontsize=8.5, bold=True)
    states = [("called", C_RNASEQ), ("low_confidence", C_WGS), ("no_call", C_LOSS)]
    for i, (lbl, col) in enumerate(states):
        bx = X_ENS + 0.15
        by = 2.65 - i * 0.52
        sb = FancyBboxPatch((bx, by), 1.7, 0.38,
                             boxstyle="round,pad=0.03",
                             facecolor=col, edgecolor="none", alpha=0.25, zorder=3)
        ax.add_patch(sb)
        ax.text(bx + 0.85, by + 0.19, lbl, ha="center", va="center",
                fontsize=7.5, color=col, fontweight="bold", zorder=4)
    _arrow(ax, X_CAL + 1.7, 3.425, X_ENS, 3.425)
    _arrow(ax, X_HARM + 2.1, 2.2, X_ENS, 2.9)

    # ── 6. Output block ──────────────────────────────────────────────────
    ax.text(X_OUT + 0.85, 4.72, "Outputs", ha="center",
            fontsize=8, fontweight="bold", color="#1e293b")
    outputs = [
        ("Consensus calls\n(JSON / TSV)", "#e0f2fe", 3.5),
        ("Discordance tags\n& taxonomy",   "#fce7f3", 2.55),
        ("Calibration\nmetrics",           "#f0fdf4", 1.6),
    ]
    for label, col, y0 in outputs:
        _box(ax, X_OUT, y0, 1.7, 0.75, label, color=col, fontsize=7.5)
        _arrow(ax, X_ENS + 2.0, y0 + 0.375, X_OUT, y0 + 0.375)

    # ── Footer note ──────────────────────────────────────────────────────
    ax.text(7, 0.18,
            "PIHLA · Nextflow-based HLA ensemble typing  ·  IMGT/HLA v3.59.0",
            ha="center", va="bottom", fontsize=7, color="#94a3b8", style="italic")

    fig.patch.set_facecolor("white")
    _save(fig, figures_dir, "figure_1_workflow_architecture")


# ---------------------------------------------------------------------------
# Figure 2 — Method Comparison (grouped bar chart)
# ---------------------------------------------------------------------------

def figure_2(tables_dir: Path, figures_dir: Path) -> None:
    df = _load(tables_dir, "method_comparison.tsv")
    if df is None:
        return
    df["overall_correct_call_rate"] = df["overall_correct_call_rate"].astype(float)

    modality_order = ["wes", "wgs", "rnaseq"]
    # Only keep modalities present in data
    modality_order = [m for m in modality_order if m in df["modality"].values]

    # Method display order and colors
    method_order = []
    for mt in df["method_type"].unique():
        for m in df.loc[df["method_type"] == mt, "method"].unique():
            if m not in method_order:
                method_order.append(m)

    def _method_color(method: str, mtype: str) -> str:
        if mtype == "ensemble":
            return C_ENSEMBLE
        if mtype == "baseline":
            return C_MAJORITY
        # single-tool: vary shade of steelblue
        tools = df.loc[df["method_type"] == "single_tool", "method"].unique().tolist()
        idx = tools.index(method) if method in tools else 0
        blues = ["#60A5FA", "#3B82F6", "#1D4ED8", "#1e40af", "#93c5fd"]
        return blues[idx % len(blues)]

    fig, ax = plt.subplots(figsize=(10, 5))

    n_methods_per_group = df.groupby("modality")["method"].nunique().max()
    group_width = 0.8
    bar_width = group_width / (n_methods_per_group + 0.5)
    gap = 0.35  # gap between modality groups

    xtick_positions = []
    xtick_labels = []
    legend_handles = {}

    for g_idx, modality in enumerate(modality_order):
        sub = df[df["modality"] == modality].copy()
        methods_here = sub["method"].tolist()
        n = len(methods_here)
        group_center = g_idx * (1.0 + gap)
        offsets = np.linspace(-(n - 1) * bar_width / 2, (n - 1) * bar_width / 2, n)

        xtick_positions.append(group_center)
        xtick_labels.append(MODALITY_LABELS.get(modality, modality.upper()))

        for offset, (_, row) in zip(offsets, sub.iterrows()):
            method = row["method"]
            mtype  = row["method_type"]
            val    = row["overall_correct_call_rate"] * 100
            color  = _method_color(method, mtype)
            x      = group_center + offset
            bar    = ax.bar(x, val, width=bar_width * 0.92, color=color,
                            edgecolor="white", linewidth=0.5, zorder=3)
            # value label
            ax.text(x, val + 0.8, f"{val:.1f}", ha="center", va="bottom",
                    fontsize=6.5, color="#374151")
            if method not in legend_handles:
                legend_handles[method] = mpatches.Patch(color=color, label=method)

        # Modality group label above bars
        ax.text(group_center, 103, MODALITY_LABELS.get(modality, modality.upper()),
                ha="center", fontsize=9, fontweight="bold",
                color=MODALITY_COLORS.get(modality, "#1e293b"))

    ax.axhline(100, color="#9CA3AF", linestyle="--", linewidth=0.8, zorder=1)
    ax.set_xticks(xtick_positions)
    ax.set_xticklabels(xtick_labels, fontsize=9)
    ax.set_ylabel("Overall correct call rate (%)", fontsize=9)
    ax.set_ylim(0, 110)
    ax.yaxis.set_major_locator(mticker.MultipleLocator(25))
    ax.set_title("Figure 2. Method comparison across modalities", fontsize=10, fontweight="bold")
    ax.legend(handles=list(legend_handles.values()), title="Method",
              loc="lower right", framealpha=0.9, fontsize=7.5)
    ax.set_axisbelow(True)
    fig.tight_layout()
    _save(fig, figures_dir, "figure_2_accuracy_comparison")


# ---------------------------------------------------------------------------
# Figure 3 — Per-Gene Gains (diverging horizontal bars)
# ---------------------------------------------------------------------------

def figure_3(tables_dir: Path, figures_dir: Path) -> None:
    df = _load(tables_dir, "per_gene_gain.tsv")
    if df is None:
        return
    df["gain_vs_majority"] = df["gain_vs_majority"].astype(float)

    modality_order = ["wes", "wgs", "rnaseq"]
    modality_order = [m for m in modality_order if m in df["modality"].values]
    gene_order = ["A", "B", "C", "DRB1", "DQB1"]
    gene_order = [g for g in gene_order if g in df["gene"].values]

    ncols = len(modality_order)
    fig, axes = plt.subplots(1, ncols, figsize=(4 * ncols, max(3, 0.6 * len(gene_order) + 1.5)),
                              sharey=True)
    if ncols == 1:
        axes = [axes]

    for ax, modality in zip(axes, modality_order):
        sub = df[df["modality"] == modality].set_index("gene").reindex(gene_order).dropna()
        gains = sub["gain_vs_majority"].values * 100
        genes = sub.index.tolist()
        colors = [C_GAIN if g >= 0 else C_LOSS for g in gains]

        bars = ax.barh(genes, gains, color=colors, edgecolor="white",
                       linewidth=0.5, height=0.6, zorder=3)
        ax.axvline(0, color="#374151", linewidth=0.8, zorder=4)

        for bar, val in zip(bars, gains):
            ha = "left" if val >= 0 else "right"
            xpos = val + (0.3 if val >= 0 else -0.3)
            ax.text(xpos, bar.get_y() + bar.get_height() / 2,
                    f"{val:+.1f}%", ha=ha, va="center",
                    fontsize=7, color="#374151")

        ax.set_title(MODALITY_LABELS.get(modality, modality.upper()),
                     fontsize=9, fontweight="bold",
                     color=MODALITY_COLORS.get(modality, "#1e293b"))
        ax.set_xlabel("Accuracy gain (pp)", fontsize=8.5)
        ax.set_xlim(-25, 65)
        ax.xaxis.set_major_locator(mticker.MultipleLocator(20))
        ax.set_axisbelow(True)
        ax.grid(axis="x", alpha=0.35)
        ax.grid(axis="y", alpha=0)

    axes[0].set_ylabel("HLA gene locus", fontsize=9)
    fig.suptitle(
        "Figure 3. Per-gene accuracy gain of weighted consensus over majority vote",
        fontsize=10, fontweight="bold", y=1.01,
    )
    gain_patch = mpatches.Patch(color=C_GAIN, label="Gain (weighted > majority)")
    loss_patch = mpatches.Patch(color=C_LOSS, label="Loss (weighted < majority)")
    axes[-1].legend(handles=[gain_patch, loss_patch], loc="lower right",
                    fontsize=7.5, framealpha=0.9)
    fig.tight_layout()
    _save(fig, figures_dir, "figure_3_per_gene_gains")


# ---------------------------------------------------------------------------
# Figure 4 — Confidence Calibration (scatter + diagonal)
# ---------------------------------------------------------------------------

def figure_4(tables_dir: Path, figures_dir: Path) -> None:
    bins_path = tables_dir / "confidence_bin_summary.tsv"
    bins_df = _load(tables_dir, "confidence_bin_summary.tsv")
    summ_df = _load(tables_dir, "confidence_calibration_summary.tsv")
    if bins_df is None:
        if bins_path.exists():
            _empty_state_figure(
                figures_dir,
                "figure_4_confidence_calibration",
                "Figure 4. Confidence calibration curves",
                "No non-empty confidence calibration bins were available for this benchmark run.",
            )
        return

    bins_df["mean_confidence"]   = bins_df["mean_confidence"].astype(float)
    bins_df["observed_accuracy"] = bins_df["observed_accuracy"].astype(float)
    bins_df["n_rows"]            = bins_df["n_rows"].astype(float)

    combos = bins_df.groupby(["tool", "modality"]).size().reset_index()[["tool", "modality"]]
    n = len(combos)
    ncols = min(3, n)
    nrows = int(np.ceil(n / ncols))

    fig, axes = plt.subplots(nrows, ncols, figsize=(4.5 * ncols, 4 * nrows),
                              squeeze=False)

    for idx, (_, row) in enumerate(combos.iterrows()):
        tool, mod = row["tool"], row["modality"]
        ax = axes[idx // ncols][idx % ncols]
        sub = bins_df[(bins_df["tool"] == tool) & (bins_df["modality"] == mod)]

        # Perfect calibration diagonal
        ax.plot([0, 1], [0, 1], linestyle="--", color="#9CA3AF", linewidth=1.0,
                label="Perfect calibration", zorder=1)
        # ±0.1 calibration band
        ax.fill_between([0, 1], [-0.1, 0.9], [0.1, 1.1], color="#E5E7EB",
                        alpha=0.4, zorder=0)

        sizes = (sub["n_rows"] / sub["n_rows"].max()) * 250 + 30
        scatter = ax.scatter(sub["mean_confidence"], sub["observed_accuracy"],
                             s=sizes, c=MODALITY_COLORS.get(mod, C_PIHLA),
                             edgecolors="white", linewidths=0.5, zorder=3,
                             alpha=0.85)

        # Annotate Brier + ECE from summary
        if summ_df is not None:
            s = summ_df[(summ_df["tool"] == tool) & (summ_df["modality"] == mod)]
            if not s.empty:
                brier = float(s.iloc[0]["brier_score"])
                ece   = float(s.iloc[0]["expected_calibration_error"])
                ax.text(0.04, 0.93,
                        f"Brier = {brier:.4f}\nECE = {ece:.3f}",
                        transform=ax.transAxes,
                        fontsize=7, va="top", color="#374151",
                        bbox=dict(boxstyle="round,pad=0.25", facecolor="white",
                                  edgecolor="#E5E7EB", linewidth=0.6))

        ax.set_xlim(-0.05, 1.05)
        ax.set_ylim(-0.05, 1.15)
        ax.set_xlabel("Mean predicted confidence", fontsize=8)
        ax.set_ylabel("Observed accuracy", fontsize=8)
        ax.set_title(f"{tool}  ·  {MODALITY_LABELS.get(mod, mod.upper())}",
                     fontsize=8.5, fontweight="bold",
                     color=MODALITY_COLORS.get(mod, "#1e293b"))
        ax.set_aspect("equal", adjustable="box")
        ax.grid(alpha=0.35)
        ax.legend(fontsize=7, loc="upper left")

    # Hide unused panels
    for i in range(n, nrows * ncols):
        axes[i // ncols][i % ncols].set_visible(False)

    fig.suptitle(
        "Figure 4. Confidence calibration curves",
        fontsize=11, fontweight="bold", y=1.01,
    )
    # Size legend
    sizes_legend = [20, 50, 100]
    handles = [plt.scatter([], [], s=(s / bins_df["n_rows"].max()) * 250 + 30,
                           c="grey", alpha=0.6, label=str(s))
               for s in sizes_legend if s <= bins_df["n_rows"].max()]
    if handles:
        fig.legend(handles=handles, title="n rows in bin",
                   loc="lower right", fontsize=7, scatterpoints=1)
    fig.tight_layout()
    _save(fig, figures_dir, "figure_4_confidence_calibration")


# ---------------------------------------------------------------------------
# Figure 5 — Abstention Tradeoff (dual-Y line chart)
# ---------------------------------------------------------------------------

def figure_5(tables_dir: Path, figures_dir: Path) -> None:
    df = _load(tables_dir, "abstention_tradeoff.tsv")
    if df is None:
        return
    df["min_support"]          = df["min_support"].astype(float)
    df["call_rate"]            = df["call_rate"].astype(float)
    df["accuracy_among_called"] = df["accuracy_among_called"].astype(float)

    fig, ax1 = plt.subplots(figsize=(7, 4.5))
    ax2 = ax1.twinx()

    # Check if there is a modality column for per-modality lines
    has_modality = "modality" in df.columns and df["modality"].nunique() > 1

    if has_modality:
        for modality, sub in df.groupby("modality"):
            col = MODALITY_COLORS.get(modality, "#64748b")
            lbl = MODALITY_LABELS.get(modality, modality.upper())
            ax1.plot(sub["min_support"], sub["accuracy_among_called"],
                     color=col, linewidth=2.0, marker="o", ms=5,
                     label=f"{lbl} accuracy")
            ax2.plot(sub["min_support"], sub["call_rate"],
                     color=col, linewidth=1.5, linestyle="--", marker="s", ms=4,
                     label=f"{lbl} call rate", alpha=0.7)
    else:
        ax1.plot(df["min_support"], df["accuracy_among_called"],
                 color=C_PIHLA, linewidth=2.0, marker="o", ms=5,
                 label="Accuracy among called")
        ax2.plot(df["min_support"], df["call_rate"],
                 color="#F97316", linewidth=1.5, linestyle="--", marker="s", ms=4,
                 label="Call rate", alpha=0.85)
        # Shaded region between curves
        ax1.fill_between(df["min_support"],
                         df["accuracy_among_called"], df["call_rate"],
                         alpha=0.08, color=C_PIHLA)

    # Tuned threshold marker
    tuned = 0.55
    ax1.axvline(tuned, color="#64748B", linestyle=":", linewidth=1.2)
    ax1.text(tuned + 0.005, ax1.get_ylim()[0] + 0.02,
             f"tuned\nthreshold\n{tuned}", fontsize=7, color="#64748B", va="bottom")

    ax1.set_xlabel("Minimum support threshold", fontsize=9)
    ax1.set_ylabel("Accuracy among called", fontsize=9, color=C_PIHLA)
    ax2.set_ylabel("Call rate", fontsize=9, color="#F97316")
    ax1.set_ylim(0, 1.1)
    ax2.set_ylim(0, 1.1)
    ax1.tick_params(axis="y", labelcolor=C_PIHLA)
    ax2.tick_params(axis="y", labelcolor="#F97316")
    ax2.spines["right"].set_visible(True)
    ax2.spines["top"].set_visible(False)

    # Combined legend
    lines1, labels1 = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax1.legend(lines1 + lines2, labels1 + labels2,
               loc="lower left", fontsize=7.5, framealpha=0.9)

    ax1.set_title("Figure 5. Abstention tradeoff: call rate vs. accuracy",
                  fontsize=10, fontweight="bold")
    ax1.grid(axis="y", alpha=0.35)
    ax1.grid(axis="x", alpha=0.2)
    fig.tight_layout()
    _save(fig, figures_dir, "figure_5_abstention_tradeoff")


# ---------------------------------------------------------------------------
# Figure 6 — Discordance Taxonomy (stacked bar)
# ---------------------------------------------------------------------------

def figure_6(tables_dir: Path, figures_dir: Path) -> None:
    df = _load(tables_dir, "discordance_summary.tsv")
    if df is None:
        return
    df["n_events"] = df["n_events"].astype(int)

    scopes = df["scope"].unique().tolist()
    tags   = df["tag"].unique().tolist()

    fig, ax = plt.subplots(figsize=(max(5, 1.8 * len(scopes)), 4.5))

    bottoms = np.zeros(len(scopes))
    scope_idx = {s: i for i, s in enumerate(scopes)}

    for tag in tags:
        vals = []
        for scope in scopes:
            sub = df[(df["scope"] == scope) & (df["tag"] == tag)]
            vals.append(int(sub["n_events"].iloc[0]) if not sub.empty else 0)
        color = DISCORD_COLORS.get(tag, "#94A3B8")
        bars = ax.bar(range(len(scopes)), vals, bottom=bottoms,
                      color=color, edgecolor="white", linewidth=0.5,
                      label=tag.replace("_", " ").title(), zorder=3)
        # Value labels
        for i, (val, bot) in enumerate(zip(vals, bottoms)):
            if val > 0:
                ax.text(i, bot + val / 2, str(val),
                        ha="center", va="center", fontsize=8,
                        color="white" if val > 0 else "#374151", fontweight="bold")
        bottoms += np.array(vals)

    ax.set_xticks(range(len(scopes)))
    scope_labels = [s.replace("_", "-") for s in scopes]
    ax.set_xticklabels(scope_labels, fontsize=8.5)
    ax.set_ylabel("Number of discordance events", fontsize=9)
    ax.set_title("Figure 6. Discordance taxonomy across modalities",
                 fontsize=10, fontweight="bold")
    ax.legend(loc="upper right", fontsize=7.5, framealpha=0.9)
    ax.set_axisbelow(True)
    total = df["n_events"].sum()
    ax.text(0.98, 0.97, f"Total events: {total}",
            transform=ax.transAxes, fontsize=8, ha="right", va="top", color="#64748B")
    ax.text(0.01, -0.14,
            "Note: counts aggregate the current authoritative modality-specific benchmark roots.",
            transform=ax.transAxes, fontsize=6.5, color="#94a3b8", style="italic")
    fig.tight_layout()
    _save(fig, figures_dir, "figure_6_discordance_taxonomy")


# ---------------------------------------------------------------------------
# Figure 7 — Learned Confidence Weights (annotated heatmap)
# ---------------------------------------------------------------------------

def figure_7(tables_dir: Path, figures_dir: Path) -> None:
    df = _load(tables_dir, "tool_confidence_weights.tsv")
    if df is None:
        return

    for col in ["base_reliability", "calibrated_confidence", "final_weight"]:
        df[col] = pd.to_numeric(df[col], errors="coerce")

    tools      = df["tool"].unique().tolist()
    modalities = ["wes", "wgs", "rnaseq"]
    modalities = [m for m in modalities if m in df["modality"].values]

    def _pivot(col: str) -> pd.DataFrame:
        return df.pivot_table(index="tool", columns="modality",
                              values=col, aggfunc="first").reindex(
            index=tools, columns=modalities)

    fw_pivot  = _pivot("final_weight")
    br_pivot  = _pivot("base_reliability")
    cc_pivot  = _pivot("calibrated_confidence")

    fig = plt.figure(figsize=(10, max(3, 1.1 * len(tools) + 1.5)))
    gs  = fig.add_gridspec(1, 4, width_ratios=[5, 1.5, 1.5, 0.4], wspace=0.05)

    ax_main  = fig.add_subplot(gs[0])
    ax_br    = fig.add_subplot(gs[1], sharey=ax_main)
    ax_cc    = fig.add_subplot(gs[2], sharey=ax_main)
    ax_cbar  = fig.add_subplot(gs[3])

    # Main heatmap
    mask = fw_pivot.isna()
    im = ax_main.imshow(fw_pivot.values.astype(float),
                        cmap="YlOrRd", vmin=0, vmax=1, aspect="auto")
    # Annotate cells
    for i, tool in enumerate(tools):
        for j, mod in enumerate(modalities):
            val = fw_pivot.iloc[i, j]
            if pd.isna(val):
                ax_main.add_patch(plt.Rectangle((j - 0.5, i - 0.5), 1, 1,
                                                fill=True, color="#e2e8f0",
                                                hatch="////", edgecolor="#94a3b8",
                                                linewidth=0.3))
                ax_main.text(j, i, "n/a", ha="center", va="center",
                             fontsize=8, color="#94a3b8")
            else:
                txt_col = "white" if val > 0.75 else "#1e293b"
                ax_main.text(j, i, f"{val:.3f}", ha="center", va="center",
                             fontsize=9, color=txt_col, fontweight="bold")

    ax_main.set_xticks(range(len(modalities)))
    ax_main.set_xticklabels([MODALITY_LABELS.get(m, m.upper()) for m in modalities],
                             fontsize=9)
    ax_main.set_yticks(range(len(tools)))
    ax_main.set_yticklabels(tools, fontsize=9)
    ax_main.set_title("Final weight\n(0.7 × reliability + 0.3 × confidence)",
                      fontsize=8.5, pad=6)
    for spine in ax_main.spines.values():
        spine.set_visible(False)

    # Side bars: base_reliability
    for i, tool in enumerate(tools):
        vals = br_pivot.iloc[i].dropna()
        mean_val = vals.mean() if not vals.empty else np.nan
        color = C_ENSEMBLE if not np.isnan(mean_val) else "#94a3b8"
        if not np.isnan(mean_val):
            ax_br.barh(i, mean_val, color=color, alpha=0.85, height=0.6,
                       edgecolor="white")
            ax_br.text(mean_val + 0.01, i, f"{mean_val:.2f}",
                       va="center", fontsize=7)
    ax_br.set_xlim(0, 1.25)
    ax_br.set_title("Base\nreliability", fontsize=8)
    ax_br.set_xlabel("")
    ax_br.tick_params(left=False, labelleft=False)
    ax_br.grid(axis="x", alpha=0.3)
    ax_br.spines["left"].set_visible(False)

    # Side bars: calibrated_confidence
    for i, tool in enumerate(tools):
        vals = cc_pivot.iloc[i].dropna()
        mean_val = vals.mean() if not vals.empty else np.nan
        if not np.isnan(mean_val):
            ax_cc.barh(i, mean_val, color=C_PIHLA, alpha=0.75, height=0.6,
                       edgecolor="white")
            ax_cc.text(mean_val + 0.01, i, f"{mean_val:.2f}",
                       va="center", fontsize=7)
        else:
            ax_cc.text(0.05, i, "—", va="center", fontsize=8, color="#94a3b8")
    ax_cc.set_xlim(0, 1.25)
    ax_cc.set_title("Calibrated\nconfidence", fontsize=8)
    ax_cc.tick_params(left=False, labelleft=False)
    ax_cc.grid(axis="x", alpha=0.3)
    ax_cc.spines["left"].set_visible(False)

    plt.colorbar(im, cax=ax_cbar, label="Final weight")
    ax_cbar.set_ylabel("Final weight", fontsize=8)

    fig.suptitle("Figure 7. Learned benchmark-derived confidence weights",
                 fontsize=10, fontweight="bold", y=1.02)
    fig.tight_layout()
    _save(fig, figures_dir, "figure_7_confidence_weights")


# ---------------------------------------------------------------------------
# Figure S1 — Per-Gene Accuracy by Tool and Modality
# ---------------------------------------------------------------------------

def figure_s1(tables_dir: Path, figures_dir: Path) -> None:
    df = _load(tables_dir, "summary_per_gene.tsv")
    if df is None:
        return
    df["overall_correct_call_rate"] = df["overall_correct_call_rate"].astype(float)

    gene_order = ["A", "B", "C", "DRB1", "DQB1"]
    gene_order = [g for g in gene_order if g in df["gene"].values]

    ncols = min(3, len(gene_order))
    nrows = int(np.ceil(len(gene_order) / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(5 * ncols, 4 * nrows),
                              squeeze=False)

    modality_order = ["wes", "wgs", "rnaseq"]

    for idx, gene in enumerate(gene_order):
        ax = axes[idx // ncols][idx % ncols]
        sub = df[df["gene"] == gene]
        tools = sub["tool"].unique().tolist()
        modalities_here = [m for m in modality_order if m in sub["modality"].values]

        n_tools = len(tools)
        n_mod   = len(modalities_here)
        bar_width = 0.8 / n_mod

        for m_idx, mod in enumerate(modalities_here):
            vals = []
            for t in tools:
                row = sub[(sub["tool"] == t) & (sub["modality"] == mod)]
                vals.append(float(row["overall_correct_call_rate"].iloc[0]) * 100
                            if not row.empty else 0.0)
            x = np.arange(n_tools) + m_idx * bar_width
            ax.bar(x, vals, width=bar_width * 0.9,
                   color=MODALITY_COLORS.get(mod, "#94A3B8"),
                   edgecolor="white", linewidth=0.5,
                   label=MODALITY_LABELS.get(mod, mod.upper()),
                   zorder=3)

        ax.set_xticks(np.arange(n_tools) + (n_mod - 1) * bar_width / 2)
        ax.set_xticklabels(tools, fontsize=8, rotation=20, ha="right")
        ax.set_ylabel("Correct call rate (%)", fontsize=8)
        ax.set_ylim(0, 115)
        ax.yaxis.set_major_locator(mticker.MultipleLocator(25))
        ax.set_title(f"HLA-{gene}", fontsize=9, fontweight="bold")
        ax.legend(fontsize=7, loc="lower right")
        ax.set_axisbelow(True)

    for i in range(len(gene_order), nrows * ncols):
        axes[i // ncols][i % ncols].set_visible(False)

    fig.suptitle("Figure S1. Per-gene accuracy by tool and modality",
                 fontsize=11, fontweight="bold", y=1.01)
    fig.tight_layout()
    _save(fig, figures_dir, "figure_s1_per_gene_accuracy")


# ---------------------------------------------------------------------------
# Figure S2 — Calibration Brier Score / ECE Heatmap
# ---------------------------------------------------------------------------

def figure_s2(tables_dir: Path, figures_dir: Path) -> None:
    summary_path = tables_dir / "confidence_calibration_summary.tsv"
    df = _load(tables_dir, "confidence_calibration_summary.tsv")
    if df is None:
        if summary_path.exists():
            _empty_state_figure(
                figures_dir,
                "figure_s2_calibration_heatmap",
                "Figure S2. Confidence calibration quality metrics (Brier score & ECE)",
                "No calibration summary rows were available for the supplementary heatmap.",
            )
        return
    df["brier_score"] = pd.to_numeric(df["brier_score"], errors="coerce")
    df["expected_calibration_error"] = pd.to_numeric(
        df["expected_calibration_error"], errors="coerce")

    tools      = df["tool"].unique().tolist()
    modalities = ["wes", "wgs", "rnaseq"]
    modalities = [m for m in modalities if m in df["modality"].values]

    def _pivot_metric(col: str) -> pd.DataFrame:
        return df.pivot_table(index="tool", columns="modality",
                              values=col, aggfunc="first").reindex(
            index=tools, columns=modalities)

    brier_p = _pivot_metric("brier_score")
    ece_p   = _pivot_metric("expected_calibration_error")

    fig, axes = plt.subplots(1, 2, figsize=(10, max(3, 0.9 * len(tools) + 1.8)))
    col_labels = [MODALITY_LABELS.get(m, m.upper()) for m in modalities]

    for ax, pivot, title, cmap in [
        (axes[0], brier_p, "Brier score\n(lower = better)", "YlOrRd"),
        (axes[1], ece_p,   "Expected calibration error\n(lower = better)", "YlOrRd"),
    ]:
        vmax = pivot.values.astype(float)
        vmax = float(np.nanmax(vmax)) if not np.all(np.isnan(vmax)) else 1.0
        im = ax.imshow(pivot.values.astype(float), cmap=cmap,
                       vmin=0, vmax=vmax, aspect="auto")
        for i, t in enumerate(tools):
            for j, m in enumerate(modalities):
                val = pivot.iloc[i, j]
                if pd.isna(val):
                    ax.add_patch(plt.Rectangle((j - 0.5, i - 0.5), 1, 1,
                                               fill=True, color="#e2e8f0",
                                               hatch="////", edgecolor="#94a3b8"))
                    ax.text(j, i, "n/a", ha="center", va="center",
                            fontsize=8, color="#94a3b8")
                else:
                    txt_col = "white" if val > vmax * 0.6 else "#1e293b"
                    ax.text(j, i, f"{val:.4f}", ha="center", va="center",
                            fontsize=8.5, color=txt_col)
        ax.set_xticks(range(len(modalities)))
        ax.set_xticklabels(col_labels, fontsize=9)
        ax.set_yticks(range(len(tools)))
        ax.set_yticklabels(tools, fontsize=9)
        ax.set_title(title, fontsize=9, pad=5)
        for sp in ax.spines.values():
            sp.set_visible(False)
        plt.colorbar(im, ax=ax, fraction=0.04, pad=0.04)

    fig.suptitle("Figure S2. Confidence calibration quality metrics (Brier score & ECE)",
                 fontsize=10, fontweight="bold", y=1.02)
    fig.tight_layout()
    _save(fig, figures_dir, "figure_s2_calibration_heatmap")


# ---------------------------------------------------------------------------
# Figure S3 — 2-Field vs 3-Field Accuracy (paired bars)
# ---------------------------------------------------------------------------

def figure_s3(tables_dir: Path, figures_dir: Path) -> None:
    df = _load(tables_dir, "ambiguity_summary.tsv")
    if df is None:
        return
    df["exact_2field_rate"] = pd.to_numeric(df["exact_2field_rate"], errors="coerce")
    df["exact_3field_rate"] = pd.to_numeric(df["exact_3field_rate"], errors="coerce")

    # Label: method + modality
    df["label"] = df["tool"] + "\n" + df["modality"].map(
        lambda m: MODALITY_LABELS.get(m, m.upper()))

    fig, ax = plt.subplots(figsize=(max(8, 1.4 * len(df)), 5))
    x    = np.arange(len(df))
    w    = 0.35
    bars1 = ax.bar(x - w / 2, df["exact_2field_rate"] * 100, width=w,
                   color=C_PIHLA, label="2-field exact match", zorder=3)
    bars2 = ax.bar(x + w / 2, df["exact_3field_rate"] * 100, width=w,
                   color=C_WGS, label="3-field exact match",
                   hatch="///", edgecolor="white", zorder=3)

    # Connecting lines between pairs
    for xi, (v2, v3) in enumerate(zip(df["exact_2field_rate"], df["exact_3field_rate"])):
        if not (pd.isna(v2) or pd.isna(v3)):
            ax.plot([xi - w / 2, xi + w / 2], [v2 * 100, v3 * 100],
                    color="#94A3B8", linewidth=0.8, zorder=4)

    ax.set_xticks(x)
    ax.set_xticklabels(df["label"], fontsize=8)
    ax.set_ylabel("Accuracy (%)", fontsize=9)
    ax.set_ylim(0, 115)
    ax.yaxis.set_major_locator(mticker.MultipleLocator(25))
    ax.legend(fontsize=8, framealpha=0.9)
    ax.set_title(
        "Figure S3. Two-field vs. three-field accuracy comparison",
        fontsize=10, fontweight="bold",
    )
    ax.text(0.01, -0.18,
            "Note: 3-field rate reported as 0 when truth/calls do not carry comparable"
            " 3-field alleles (SpecHLA WGS/WES).",
            transform=ax.transAxes, fontsize=6.5, color="#94a3b8", style="italic")
    ax.set_axisbelow(True)
    fig.tight_layout()
    _save(fig, figures_dir, "figure_s3_resolution_comparison")



# ---------------------------------------------------------------------------
# Exploratory Figures — 1000G cohort composition and coverage
# ---------------------------------------------------------------------------

def figure_x1(source_dir: Path, tables_dir: Path, figures_dir: Path) -> None:
    candidate = _source_table(source_dir, tables_dir, "candidate_seed_cohort.tsv")
    sequencing = _source_table(source_dir, tables_dir, "sequencing_manifest.tsv")
    truth = _source_table(source_dir, tables_dir, "truth_manifest.tsv")

    labels = ["Truth", "WGS", "WES", "RNA-seq"]
    values = []
    total = None
    if candidate is not None and {"truth_available", "wgs_done", "wes_done", "rnaseq_done"}.issubset(candidate.columns):
        total = len(candidate.index)
        values = [
            int(candidate["truth_available"].eq("1").sum()),
            int(candidate["wgs_done"].eq("1").sum()),
            int(candidate["wes_done"].eq("1").sum()),
            int(candidate["rnaseq_done"].eq("1").sum()),
        ]
    elif sequencing is not None and truth is not None:
        truth_samples = set(truth["sample"].dropna())
        total = len(truth_samples)
        values = [
            total,
            int(sequencing.loc[sequencing["modality"].eq("wgs") & sequencing["available"].eq("1"), "sample"].nunique()),
            int(sequencing.loc[sequencing["modality"].eq("wes") & sequencing["available"].eq("1"), "sample"].nunique()),
            int(sequencing.loc[sequencing["modality"].eq("rnaseq") & sequencing["available"].eq("1"), "sample"].nunique()),
        ]
    else:
        warnings.warn("[SKIP] exploratory modality coverage inputs not found.")
        return

    fig, ax = plt.subplots(figsize=(6.8, 4.2))
    bars = ax.bar(labels, values, color=[C_NEUTRAL, C_WGS, C_WES, C_RNASEQ], width=0.62, edgecolor="white", linewidth=1.4)
    ax.set_title("Figure X1. 1000G modality coverage overview", fontsize=10.5, fontweight="bold")
    ax.set_ylabel("Samples")
    ax.set_axisbelow(True)
    ymax = max(values + [1])
    ax.set_ylim(0, ymax * 1.22)
    for bar, value in zip(bars, values):
        suffix = f" ({value / total:.0%})" if total else ""
        ax.text(bar.get_x() + bar.get_width() / 2, value + ymax * 0.03, f"{value}{suffix}", ha="center", va="bottom", fontsize=8, color=C_NEUTRAL)
    ax.text(0.01, -0.16, "Counts summarize candidate cohort truth support and per-modality availability/completion.", transform=ax.transAxes, fontsize=7, color="#64748b")
    fig.tight_layout()
    _save(fig, figures_dir, "figure_x1_modality_coverage")


def figure_x2(source_dir: Path, tables_dir: Path, figures_dir: Path) -> None:
    cohort = _source_table(source_dir, tables_dir, "cohort_manifest.tsv")
    if cohort is None or "excluded_reason" not in cohort.columns:
        warnings.warn("[SKIP] cohort_manifest.tsv not found — skipping figure_x2.")
        return
    include_col = cohort["include"] if "include" in cohort.columns else pd.Series(["0"] * len(cohort))
    excluded = cohort.loc[include_col.fillna("0") != "1", "excluded_reason"].fillna("unspecified")
    if excluded.empty:
        summary = pd.Series([0], index=["no exclusions"], dtype=int)
    else:
        summary = excluded.value_counts().sort_values(ascending=True)

    fig, ax = plt.subplots(figsize=(8.4, max(3.4, 0.6 * len(summary.index) + 1.7)))
    ax.barh(summary.index.tolist(), summary.values.tolist(), color=C_ACCENT, alpha=0.9)
    ax.set_title("Figure X2. Cohort exclusion reasons", fontsize=10.5, fontweight="bold")
    ax.set_xlabel("Samples")
    ax.set_axisbelow(True)
    xmax = max(summary.values.tolist() + [1])
    ax.set_xlim(0, xmax * 1.18)
    for idx, value in enumerate(summary.values.tolist()):
        ax.text(value + xmax * 0.02, idx, str(value), va="center", fontsize=8, color=C_NEUTRAL)
    ax.grid(axis="x")
    fig.tight_layout()
    _save(fig, figures_dir, "figure_x2_exclusion_reasons")


def figure_x3(source_dir: Path, tables_dir: Path, figures_dir: Path) -> None:
    truth_long = _source_table(source_dir, tables_dir, "truth_long.tsv")
    truth_manifest = _source_table(source_dir, tables_dir, "truth_manifest.tsv")
    if truth_long is None or not {"sample", "gene"}.issubset(truth_long.columns):
        warnings.warn("[SKIP] truth_long.tsv not found — skipping figure_x3.")
        return

    gene_counts = truth_long.groupby("gene")["sample"].nunique().sort_values(ascending=False)
    sample_gene_counts = truth_long.groupby("sample")["gene"].nunique().sort_values()
    supported_loci = 0
    if truth_manifest is not None and "truth_supported_loci" in truth_manifest.columns and not truth_manifest.empty:
        supported_loci = max(len(str(value).split(",")) for value in truth_manifest["truth_supported_loci"].dropna())

    fig, axes = plt.subplots(1, 2, figsize=(11.4, 4.2), gridspec_kw={"width_ratios": [1.2, 1.0]})
    axes[0].bar(gene_counts.index.tolist(), gene_counts.values.tolist(), color=C_PIHLA, width=0.62)
    axes[0].set_title("Per-locus truth support", fontsize=10, fontweight="bold")
    axes[0].set_ylabel("Samples with truth")
    axes[0].tick_params(axis="x", rotation=0)
    top = max(gene_counts.values.tolist() + [1])
    for idx, value in enumerate(gene_counts.values.tolist()):
        axes[0].text(idx, value + top * 0.03, str(value), ha="center", va="bottom", fontsize=8)

    bins = np.arange(sample_gene_counts.min(), sample_gene_counts.max() + 2) - 0.5
    axes[1].hist(sample_gene_counts.values.astype(int), bins=bins, color=C_RNASEQ, edgecolor="white", linewidth=1.2)
    axes[1].set_title("Per-sample truth completeness", fontsize=10, fontweight="bold")
    axes[1].set_xlabel("Truth-supported loci per sample")
    axes[1].set_ylabel("Samples")
    if supported_loci:
        axes[1].axvline(supported_loci, color=C_ACCENT, linestyle="--", linewidth=1.4)
        axes[1].text(supported_loci + 0.08, axes[1].get_ylim()[1] * 0.9, f"Expected loci: {supported_loci}", color=C_ACCENT, fontsize=7)

    fig.suptitle("Figure X3. 1000G truth locus coverage", fontsize=10.5, fontweight="bold", y=1.02)
    fig.tight_layout()
    _save(fig, figures_dir, "figure_x3_truth_locus_coverage")

# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

FIGURE_MAP = {
    "1":  lambda t, f: figure_1(f),
    "2":  lambda t, f: figure_2(t, f),
    "3":  lambda t, f: figure_3(t, f),
    "4":  lambda t, f: figure_4(t, f),
    "5":  lambda t, f: figure_5(t, f),
    "6":  lambda t, f: figure_6(t, f),
    "7":  lambda t, f: figure_7(t, f),
    "s1": lambda t, f: figure_s1(t, f),
    "s2": lambda t, f: figure_s2(t, f),
    "s3": lambda t, f: figure_s3(t, f),
}

EXPLORATORY_MAP = {
    "x1": figure_x1,
    "x2": figure_x2,
    "x3": figure_x3,
}


def main() -> None:
    repo_root = Path(__file__).resolve().parent.parent
    default_tables = repo_root / "analysis/1000g_realdata/tables"
    default_figures = repo_root / "analysis/1000g_realdata/figures"

    parser = argparse.ArgumentParser(
        description="Generate PIHLA publication figures from benchmark TSV tables.")
    parser.add_argument("--tables-dir", type=Path, default=default_tables,
                        help="Directory containing benchmark TSV output tables.")
    parser.add_argument("--figures-dir", type=Path, default=default_figures,
                        help="Directory to write figure files into.")
    parser.add_argument("--source-dir", type=Path, default=None,
                        help="Directory containing optional source TSVs such as truth_long.tsv and candidate_seed_cohort.tsv.")
    parser.add_argument("--figures", type=str, default=None,
                        help="Comma-separated list of benchmark figures to generate, e.g. 2,3,s1. Default: all benchmark figures.")
    parser.add_argument("--exploratory-figures", type=str, default="x1,x2,x3",
                        help="Comma-separated list of exploratory figures to generate when enabled.")
    parser.add_argument("--include-exploratory", choices=("auto", "always", "never"), default="auto",
                        help="Generate exploratory 1000G cohort figures when source TSVs are available.")
    parser.add_argument("--skip-captions", action="store_true",
                        help="Do not write captions.md.")
    args = parser.parse_args()

    tables_dir = args.tables_dir
    figures_dir = args.figures_dir
    source_dir = _infer_source_dir(tables_dir, args.source_dir)
    figures_dir.mkdir(parents=True, exist_ok=True)

    if args.figures:
        keys = [k.strip().lower() for k in args.figures.split(",") if k.strip()]
    else:
        keys = list(FIGURE_MAP.keys())

    unknown = [k for k in keys if k not in FIGURE_MAP]
    if unknown:
        print(f"Unknown figure keys: {unknown}. Valid: {list(FIGURE_MAP.keys())}")
        sys.exit(1)

    exploratory_requested = args.include_exploratory != "never"
    exploratory_inputs_present = any((source_dir / name).exists() for name in (
        "candidate_seed_cohort.tsv",
        "truth_long.tsv",
        "sequencing_source_from_results.tsv",
        "cohort_manifest.tsv",
    )) or (tables_dir / "cohort_manifest.tsv").exists()
    run_exploratory = exploratory_requested and (args.include_exploratory == "always" or exploratory_inputs_present)
    xkeys = [k.strip().lower() for k in args.exploratory_figures.split(",") if k.strip()]

    print(f"Tables dir : {tables_dir}")
    print(f"Figures dir: {figures_dir}")
    print(f"Source dir : {source_dir}")
    print(f"Generating benchmark figures: {', '.join(keys)}")
    if run_exploratory:
        print(f"Generating exploratory figures: {', '.join(xkeys)}")
    print()

    for key in keys:
        print(f"→ Figure {key.upper()}")
        try:
            FIGURE_MAP[key](tables_dir, figures_dir)
        except Exception as exc:
            print(f"  ✗  ERROR generating figure {key}: {exc}")
            import traceback
            traceback.print_exc()

    if run_exploratory:
        for key in xkeys:
            if key not in EXPLORATORY_MAP:
                print(f"  ✗  Unknown exploratory figure key: {key}")
                continue
            print(f"→ Figure {key.upper()}")
            try:
                EXPLORATORY_MAP[key](source_dir, tables_dir, figures_dir)
            except Exception as exc:
                print(f"  ✗  ERROR generating figure {key}: {exc}")
                import traceback
                traceback.print_exc()

    if not args.skip_captions:
        _write_captions(figures_dir, tables_dir, include_exploratory=run_exploratory)

    print("\nDone.")


if __name__ == "__main__":
    main()
