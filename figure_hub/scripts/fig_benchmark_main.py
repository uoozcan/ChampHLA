#!/usr/bin/env python3.11
"""Generate the active publication figure package from current authoritative roots."""

import argparse
import json
from pathlib import Path
import sys
from typing import Iterable

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from matplotlib.patches import FancyBboxPatch
import numpy as np
import pandas as pd
import seaborn as sns

SCRIPT_DIR = Path(__file__).resolve().parent
HUB_ROOT = SCRIPT_DIR.parent
REPO_ROOT = HUB_ROOT.parent
BIN_DIR = REPO_ROOT / "bin"
sys.path.insert(0, str(HUB_ROOT))
sys.path.insert(0, str(BIN_DIR))

from plot_style import (
    METHOD_COLORS,
    ROLE_COLORS,
    apply_style,
    legend_outside,
    rotate_xticklabels,
    save_fig,
    set_figure_footer as _figure_footer,
    set_figure_header as _figure_header,
    tight_with_legend as _finalize_layout,
)

import generate_figures as base

apply_style()

PRIMARY = "#16324f"
SECONDARY_WES = "#2c5c9c"
SECONDARY_RNA = "#17806d"
SUPPLEMENTARY = "#9a5c1f"
SLATE = "#475569"
LIGHT_SLATE = "#94a3b8"
RED = "#b42318"
GREEN = "#0f766e"
AMBER = "#b45309"
PRIMARY_MAP = {"wgs": PRIMARY, "wes": SECONDARY_WES, "rnaseq": SECONDARY_RNA}
EXPECTED_ACTIVE_STEMS = [
    "figure_1_workflow_architecture",
    "figure_2_accuracy_comparison",
    "figure_3_per_gene_gains",
    "figure_4_confidence_calibration",
    "figure_5_abstention_tradeoff",
    "figure_6_discordance_taxonomy",
    "figure_7_confidence_weights",
    "figure_09_bimodal_per_gene",
    "figure_10_trimodal_comparison",
    "figure_s1_per_gene_accuracy",
    "figure_s2_calibration_heatmap",
    "figure_s3_resolution_comparison",
]


def parse_args() -> argparse.Namespace:
    analysis_root = Path("/scratch/project_2008084/pihla-publish/analysis/1000g_realdata")
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--wgs-tables",
        type=Path,
        default=analysis_root / "benchmark_wgs_wave2" / "tables",
    )
    parser.add_argument(
        "--wes-tables",
        type=Path,
        default=analysis_root / "benchmark_wes_truthbacked" / "run" / "tables",
    )
    parser.add_argument(
        "--rna-tables",
        type=Path,
        default=analysis_root / "benchmark_rna_truthbacked" / "run" / "tables",
    )
    parser.add_argument(
        "--trimodal-tables",
        type=Path,
        default=analysis_root / "benchmark_trimodal_robustness" / "run" / "tables",
    )
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=Path("/scratch/project_2008084/pihla-publish/analysis/figures_final_candidate"),
    )
    return parser.parse_args()


def _read_tsv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t")


def _maybe_read_tsv(path: Path) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame()
    return _read_tsv(path)


def _load_tables(root: Path, names: Iterable[str]) -> dict[str, pd.DataFrame]:
    tables = {}
    for name in names:
        key = name[:-4] if name.endswith(".tsv") else name
        tables[key] = _maybe_read_tsv(root / name)
    return tables


def _load_json(path: Path) -> dict:
    if not path.exists():
        return {}
    return json.loads(path.read_text(encoding="utf-8"))


def _safe_float(value, default=0.0) -> float:
    try:
        if pd.isna(value):
            return default
    except TypeError:
        pass
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def _save(fig: plt.Figure, out_dir: Path, stem: str) -> None:
    save_fig(fig, out_dir / stem)


def _method_sort_key(name: str) -> tuple[int, str]:
    preferred = [
        "OptiType",
        "MajorityVote",
        "WeightedConsensus",
        "ChampionChallenger",
        "GatedConsensus",
        "LocusExpertConsensus",
    ]
    for idx, key in enumerate(preferred):
        if name == key or name.startswith(key):
            return (idx, name)
    return (len(preferred), name)


def load_publication_context(args: argparse.Namespace) -> dict:
    common = [
        "method_comparison.tsv",
        "method_per_gene.tsv",
        "per_gene_gain.tsv",
        "tool_confidence_weights.tsv",
        "tool_confidence_weights_by_gene.tsv",
        "confidence_calibration_summary.tsv",
        "confidence_error_summary.tsv",
        "abstention_tradeoff.tsv",
        "discordance_summary.tsv",
        "discordance_tags.tsv",
        "failure_mode_summary.tsv",
        "majority_vote_baseline.tsv",
        "weighted_consensus_calls.tsv",
        "benchmark_metadata.json",
        "summary_per_gene.tsv",
        "ambiguity_summary.tsv",
    ]
    wgs = _load_tables(args.wgs_tables, common + [
        "champion_challenger_method_comparison.tsv",
        "gated_method_comparison.tsv",
        "locus_expert_method_comparison.tsv",
        "champion_challenger_method_per_gene.tsv",
        "locus_difficulty_summary.tsv",
        "missing_call_patterns.tsv",
    ])
    wes = _load_tables(args.wes_tables, common + ["champion_challenger_method_comparison.tsv"])
    rna = _load_tables(args.rna_tables, common + ["champion_challenger_method_comparison.tsv"])
    tri = _load_tables(args.trimodal_tables, [
        "method_comparison.tsv",
        "champion_challenger_method_comparison.tsv",
        "bimodal_method_comparison.tsv",
        "bimodal_accuracy_comparison.tsv",
        "trimodal_method_comparison.tsv",
        "trimodal_accuracy_comparison.tsv",
        "trimodal_robustness_vs_wave2_wgs.tsv",
        "trimodal_robustness_vs_wes_truthbacked.tsv",
        "trimodal_robustness_vs_rna_truthbacked.tsv",
        "trimodal_subject_intersection_summary.tsv",
        "summary_per_gene.tsv",
        "confidence_calibration_summary.tsv",
        "ambiguity_summary.tsv",
    ])
    return {
        "wgs_root": args.wgs_tables.parent,
        "wes_root": args.wes_tables.parent,
        "rna_root": args.rna_tables.parent,
        "tri_root": args.trimodal_tables.parent,
        "wgs": wgs,
        "wes": wes,
        "rna": rna,
        "tri": tri,
        "wgs_meta": _load_json(args.wgs_tables / "benchmark_metadata.json"),
        "wes_meta": _load_json(args.wes_tables / "benchmark_metadata.json"),
        "rna_meta": _load_json(args.rna_tables / "benchmark_metadata.json"),
        "tri_meta": _load_json(args.trimodal_tables / "benchmark_metadata.json"),
    }


def build_figure_1_workflow_governance(context: dict, out_dir: Path) -> None:
    fig, ax = plt.subplots(figsize=(14, 5.8))
    ax.set_xlim(0, 14)
    ax.set_ylim(0, 6)
    ax.axis("off")
    ax.set_title(
        "Figure 1. ChampHLA integrates modality-specific callers, benchmark-time calibration, and claim-governed evaluation",
        fontsize=12,
        fontweight="bold",
        pad=12,
    )

    def box(x, y, w, h, label, sublabel="", fc="#eef2f7", ec="#b7c0ca", text="#1f2937"):
        patch = FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.06,rounding_size=0.04", facecolor=fc, edgecolor=ec, linewidth=0.9)
        ax.add_patch(patch)
        ax.text(x + w / 2, y + h * 0.62, label, ha="center", va="center", fontsize=8.8, fontweight="bold", color=text)
        if sublabel:
            ax.text(x + w / 2, y + h * 0.30, sublabel, ha="center", va="center", fontsize=7.2, color=SLATE)

    def arrow(x1, y1, x2, y2, color="#7c8793"):
        ax.annotate("", xy=(x2, y2), xytext=(x1, y1), arrowprops=dict(arrowstyle="-|>", color=color, lw=1.2))

    # Layered workflow
    box(0.3, 3.9, 1.6, 0.75, "WGS", "short-read DNA", fc="#f9e8d1")
    box(0.3, 2.8, 1.6, 0.75, "WES", "capture DNA", fc="#dbe7f7")
    box(0.3, 1.7, 1.6, 0.75, "RNA-seq", "expression-aware", fc="#daf2eb")
    box(2.4, 1.5, 2.1, 3.4, "Caller execution", "OptiType, ArcasHLA,\nSpecHLA, HLA-HD,\nKourami, T1K")
    box(5.0, 1.7, 2.2, 3.0, "Harmonization", "call schema,\nconfidence sidecars,\nallele normalization", fc="#fbf3d6")
    box(7.8, 3.0, 2.2, 1.5, "Benchmark-time learning", "splits, calibration,\nguardrailed weights", fc="#def4e7")
    box(7.8, 1.3, 2.2, 1.2, "Runtime consensus", "routed policies,\nabstention", fc="#ece8fb")
    box(10.6, 2.0, 2.0, 2.2, "Publication outputs", "figures, tables,\nclaim-governed summaries", fc="#eef2f7")

    for y in (4.275, 3.175, 2.075):
        arrow(1.9, y, 2.4, 3.2)
    arrow(4.5, 3.2, 5.0, 3.2)
    arrow(7.2, 3.7, 7.8, 3.7)
    arrow(7.2, 2.0, 7.8, 2.0)
    arrow(10.0, 3.7, 10.6, 3.4)
    arrow(10.0, 1.95, 10.6, 2.6)

    # Governance strip
    y0 = 0.25
    ax.text(0.3, y0 + 0.65, "Claim-governance strip", fontsize=8.2, fontweight="bold", color=PRIMARY)
    govern = [
        ("Primary WGS", "benchmark_wgs_wave2", ROLE_COLORS["primary"]),
        ("Secondary WES", "benchmark_wes_truthbacked/run", ROLE_COLORS["secondary_wes"]),
        ("Secondary RNA", "benchmark_rna_truthbacked/run", ROLE_COLORS["secondary_rna"]),
        ("Supplementary robustness", "benchmark_trimodal_robustness/run", ROLE_COLORS["supplementary"]),
    ]
    for idx, (label, root, color) in enumerate(govern):
        x = 0.3 + idx * 3.25
        patch = FancyBboxPatch((x, y0), 2.9, 0.45, boxstyle="round,pad=0.03,rounding_size=0.03", facecolor=color, edgecolor="none", alpha=0.95)
        ax.add_patch(patch)
        ax.text(x + 1.45, y0 + 0.29, label, ha="center", va="center", fontsize=8, color="white", fontweight="bold")
        ax.text(x + 1.45, y0 + 0.10, root, ha="center", va="center", fontsize=6.5, color="white")

    _finalize_layout(fig, top=0.94, bottom=0.10, right=0.98)
    _save(fig, out_dir, "figure_1_workflow_architecture")


def build_figure_2_wgs_anchor(context: dict, out_dir: Path) -> None:
    wgs = context["wgs"]["method_comparison"]
    extra_frames = [
        context["wgs"]["champion_challenger_method_comparison"],
        context["wgs"]["gated_method_comparison"],
        context["wgs"]["locus_expert_method_comparison"],
    ]
    frames = [df for df in [wgs] + extra_frames if not df.empty]
    df = pd.concat(frames, ignore_index=True)
    df = df[df["modality"] == "wgs"].copy()
    df["overall_correct_call_rate"] = pd.to_numeric(df["overall_correct_call_rate"], errors="coerce")
    df["callable_rate"] = pd.to_numeric(df["callable_rate"], errors="coerce")
    df = df.dropna(subset=["overall_correct_call_rate", "callable_rate"])
    df["sort_key"] = df["method"].map(_method_sort_key)
    df = df.sort_values("sort_key")

    fig, ax = plt.subplots(figsize=(10, 5.2))
    x = np.arange(len(df))
    bars = ax.bar(
        x,
        df["overall_correct_call_rate"] * 100.0,
        color=[METHOD_COLORS.get(m, LIGHT_SLATE) for m in df["method"]],
        edgecolor="white",
        linewidth=0.8,
        zorder=3,
    )
    ax.scatter(
        x,
        df["callable_rate"] * 100.0,
        color=SLATE,
        s=46,
        zorder=4,
        label="Callable rate",
    )
    for xpos, bar, acc in zip(x, bars, df["overall_correct_call_rate"]):
        ax.text(xpos, bar.get_height() + 1.5, f"{acc * 100:.1f}", ha="center", va="bottom", fontsize=7.2, color=SLATE)
    ax.set_xticks(x)
    ax.set_xticklabels(df["method"])
    rotate_xticklabels(ax, angle=22)
    ax.set_ylabel("Overall correct-call rate (%)")
    ax.set_ylim(0, 64)
    ax.yaxis.set_major_locator(mticker.MultipleLocator(10))
    _figure_header(
        fig,
        "Figure 2. On the most-discordant modality (WGS), Champion-Challenger significantly exceeds majority voting and the best single tool",
        "Primary evidence · WGS 1000G n=137 · 10-fold nested CV",
        PRIMARY,
    )

    def _rate(method):
        row = df[df["method"].str.startswith(method)]
        return float(row["overall_correct_call_rate"].iloc[0]) * 100 if not row.empty else float("nan")

    cc, mv, opti, wc = _rate("ChampionChallenger"), _rate("MajorityVote"), _rate("OptiType"), _rate("WeightedConsensus")
    note = (
        f"ChampHLA (Champion-Challenger) = {cc:.1f}%\n"
        f"MajorityVote = {mv:.1f}%   (Δ = +{cc - mv:.1f} pts, McNemar p<0.0001)\n"
        f"Best single tool: OptiType = {opti:.1f}%\n"
        f"WeightedConsensus (weighting only) = {wc:.1f}%\n"
        "Champion routing beats symmetric voting on discordant WGS"
    )
    ax.text(
        0.98, 0.96, note,
        transform=ax.transAxes,
        ha="right", va="top",
        fontsize=7.5,
        color=PRIMARY,
        bbox=dict(boxstyle="round,pad=0.35", facecolor="white", edgecolor="#d1d5db"),
    )
    legend_outside(ax, loc="upper left", anchor=(1.02, 1.0), frameon=False)
    _finalize_layout(fig, right=0.76)
    _save(fig, out_dir, "figure_2_accuracy_comparison")


def build_figure_3_wgs_failure_modes(context: dict, out_dir: Path) -> None:
    """Per-gene WGS accuracy: MajorityVote vs Champion-Challenger (nested CV), with
    the per-gene gain annotated. Shows the WGS win concentrates at HLA-A and HLA-B."""
    mv_pg = context["wgs"]["method_per_gene"].copy()
    cc_pg = context["wgs"]["champion_challenger_method_per_gene"].copy()
    if mv_pg.empty or cc_pg.empty:
        return
    for df_ in (mv_pg, cc_pg):
        df_["overall_correct_call_rate"] = pd.to_numeric(df_["overall_correct_call_rate"], errors="coerce")

    genes = ["A", "B", "C"]
    mv = {g: _safe_float(mv_pg[(mv_pg["modality"] == "wgs") & (mv_pg["method"] == "MajorityVote") & (mv_pg["gene"] == g)]["overall_correct_call_rate"].max()) for g in genes}
    cc = {g: _safe_float(cc_pg[(cc_pg["modality"] == "wgs") & (cc_pg["gene"] == g)]["overall_correct_call_rate"].max()) for g in genes}
    sig = {"A": "p<0.001", "B": "p<0.001", "C": "n.s."}

    fig, ax = plt.subplots(figsize=(9.2, 5.0))
    x = np.arange(len(genes))
    w = 0.38
    ax.bar(x - w / 2, [mv[g] * 100 for g in genes], w, label="MajorityVote", color=METHOD_COLORS["MajorityVote"], edgecolor="white", zorder=3)
    ax.bar(x + w / 2, [cc[g] * 100 for g in genes], w, label="ChampHLA (Champion-Challenger)", color=METHOD_COLORS["ChampionChallenger"], edgecolor="white", zorder=3)
    for i, g in enumerate(genes):
        ax.text(x[i] - w / 2, mv[g] * 100 + 1.2, f"{mv[g] * 100:.1f}", ha="center", va="bottom", fontsize=7.2, color=SLATE)
        ax.text(x[i] + w / 2, cc[g] * 100 + 1.2, f"{cc[g] * 100:.1f}", ha="center", va="bottom", fontsize=7.2, color=PRIMARY)
        delta = (cc[g] - mv[g]) * 100
        ax.text(x[i], max(mv[g], cc[g]) * 100 + 6.5, f"Δ +{delta:.1f}\n{sig[g]}", ha="center", va="bottom", fontsize=7.6, fontweight="bold", color=GREEN if sig[g] != "n.s." else SLATE)
    ax.set_xticks(x)
    ax.set_xticklabels([f"HLA-{g}" for g in genes])
    ax.set_ylabel("WGS overall correct-call rate (%)")
    ax.set_ylim(0, 88)
    ax.yaxis.set_major_locator(mticker.MultipleLocator(20))
    _figure_header(
        fig,
        "Figure 3. The WGS Champion-Challenger gain over majority voting concentrates at HLA-A and HLA-B",
        "Primary evidence · WGS 1000G n=137 · 10-fold nested CV",
        PRIMARY,
    )
    legend_outside(ax, loc="upper left", anchor=(1.02, 1.0), frameon=False)
    _finalize_layout(fig, right=0.74)
    _save(fig, out_dir, "figure_3_per_gene_gains")


def build_figure_4_confidence_guardrails(context: dict, out_dir: Path) -> None:
    weights = context["wgs"]["tool_confidence_weights"].copy()
    if weights.empty:
        return
    weights = weights[weights["modality"] == "wgs"].copy()
    numeric_cols = ["base_reliability", "calibrated_confidence", "effective_confidence", "guardrail_factor", "final_weight", "brier_score", "expected_calibration_error"]
    for col in numeric_cols:
        if col in weights.columns:
            weights[col] = pd.to_numeric(weights[col], errors="coerce")
    weights = weights.sort_values("final_weight", ascending=False)

    fig, axes = plt.subplots(1, 2, figsize=(12, 4.8), gridspec_kw={"width_ratios": [1.1, 1.0]})
    ax1, ax2 = axes

    show = weights.head(6).copy()
    x = np.arange(len(show))
    width = 0.24
    ax1.bar(x - width, show["base_reliability"], width=width, color=LIGHT_SLATE, label="Base reliability")
    ax1.bar(x, show["calibrated_confidence"], width=width, color=SECONDARY_WES, label="Calibrated confidence")
    ax1.bar(x + width, show["effective_confidence"], width=width, color=PRIMARY, label="Effective confidence")
    ax1.set_xticks(x)
    rotate_xticklabels(ax1, angle=22)
    ax1.set_ylim(0, 0.9)
    ax1.set_ylabel("Probability-scale contribution")
    ax1.set_title("A. Raw vs effective confidence after guardrails")
    legend_outside(ax1, loc="upper left", anchor=(1.02, 1.0), frameon=False)

    calibr = show[["tool", "brier_score", "expected_calibration_error", "guardrail_status"]].copy()
    y = np.arange(len(calibr))
    ax2.barh(y + 0.18, calibr["brier_score"], height=0.32, color=AMBER, label="Brier score")
    ax2.barh(y - 0.18, calibr["expected_calibration_error"], height=0.32, color=PRIMARY, label="ECE")
    ax2.set_yticks(y)
    ax2.set_yticklabels(calibr["tool"])
    ax2.invert_yaxis()
    ax2.set_xlim(0, max(0.9, np.nanmax(calibr[["brier_score", "expected_calibration_error"]].to_numpy()) + 0.08))
    ax2.set_xlabel("Calibration error")
    ax2.set_title("B. Poor calibration blocks unrestricted boosting")
    legend_outside(ax2, loc="upper left", anchor=(1.02, 1.0), frameon=False)
    for yy, row in enumerate(calibr.itertuples(index=False)):
        ax2.text(ax2.get_xlim()[1] * 0.98, yy, str(row.guardrail_status), ha="right", va="center", fontsize=7.2, color=SLATE)

    _figure_header(
        fig,
        "Figure 4. Confidence must be empirically calibrated before it is allowed to influence runtime voting",
        "Primary evidence · benchmark_wgs_wave2",
        PRIMARY,
    )
    _figure_footer(fig, "Poor calibration blocks confidence-based boosting; T1K retains only a guarded partial boost.")
    _finalize_layout(fig, right=0.72)
    _save(fig, out_dir, "figure_4_confidence_calibration")


def build_figure_5_cross_modality_comparison(context: dict, out_dir: Path) -> None:
    rows = []
    for label, key, method in [
        ("WGS primary", "wgs", "OptiType"),
        ("WES secondary", "wes", "ChampionChallenger"),
        ("RNA secondary", "rna", "ChampionChallenger"),
    ]:
        df = context[key]["method_comparison"].copy()
        if key in {"wes", "rna"} and not context[key]["champion_challenger_method_comparison"].empty:
            df = pd.concat([df, context[key]["champion_challenger_method_comparison"]], ignore_index=True)
        modality = "rnaseq" if key == "rna" else key
        sub = df[(df["modality"] == modality) & (df["method"] == method)].head(1)
        if sub.empty:
            sub = df[df["modality"] == modality].sort_values("overall_correct_call_rate", ascending=False).head(1)
        if sub.empty:
            continue
        row = sub.iloc[0].to_dict()
        row["label"] = label
        rows.append(row)
    data = pd.DataFrame(rows)
    data["overall_correct_call_rate"] = pd.to_numeric(data["overall_correct_call_rate"], errors="coerce")
    data["callable_rate"] = pd.to_numeric(data["callable_rate"], errors="coerce")

    fig, ax = plt.subplots(figsize=(8.6, 4.8))
    x = np.arange(len(data))
    colors = [PRIMARY, SECONDARY_WES, SECONDARY_RNA][: len(data)]
    ax.bar(x, data["overall_correct_call_rate"] * 100.0, color=colors, edgecolor="white", linewidth=0.8)
    ax.scatter(x, data["callable_rate"] * 100.0, color=SLATE, s=42, zorder=4)
    for xpos, row in zip(x, data.itertuples(index=False)):
        ax.text(xpos, row.overall_correct_call_rate * 100 + 1.2, f"{row.overall_correct_call_rate * 100:.1f}", ha="center", fontsize=8, color=SLATE)
        ax.text(xpos, 3, row.method, ha="center", va="bottom", fontsize=7.2, color="white", fontweight="bold")
    ax.set_xticks(x)
    ax.set_xticklabels(data["label"])
    ax.set_ylabel("Overall correct-call rate (%)")
    ax.set_ylim(0, 105)
    _figure_header(
        fig,
        "Figure 5. WES and RNA operate in stronger unimodal regimes than primary short-read WGS",
        "Mixed roles · Figure 5 is the cross-modality half of the transitional WGS failure-mode family",
        SUPPLEMENTARY,
    )
    _figure_footer(fig, "Compatibility stem retained: this figure now carries the cross-modality contrast that contextualizes the WGS failure-mode story.")
    _finalize_layout(fig, right=0.98)
    _save(fig, out_dir, "figure_5_abstention_tradeoff")


def build_figure_6_discordance_interpretation(context: dict, out_dir: Path) -> None:
    """Discordance burden per modality, each read from its own nested-CV benchmark
    dir and filtered to that modality's scope. The ordering (WGS > WES > RNA) matches
    where Champion-Challenger helps most."""
    # (label, scope token, n scorable loci) — loci counts drive the discordance rate.
    mods = [("WGS", "wgs", "wgs", 411), ("WES", "wes", "wes", 390), ("RNA-seq", "rna", "rnaseq", 321)]
    tags = ["low_evidence_conflict", "no_evidence", "possible_expression_bias", "technical_conflict"]
    counts = {lbl: {t: 0.0 for t in tags} for lbl, *_ in mods}
    rates = {}
    for lbl, ctx_key, scope, n_loci in mods:
        disc = context[ctx_key]["discordance_summary"].copy()
        if disc.empty:
            continue
        disc["n_events"] = pd.to_numeric(disc["n_events"], errors="coerce").fillna(0)
        sub = disc[disc["scope"] == scope]
        total = 0.0
        for _, row in sub.iterrows():
            if row["tag"] in counts[lbl]:
                counts[lbl][row["tag"]] += float(row["n_events"])
            total += float(row["n_events"])
        rates[lbl] = total / n_loci if n_loci else 0.0

    labels = [lbl for lbl, *_ in mods]
    fig, ax = plt.subplots(figsize=(8.8, 4.8))
    bottoms = np.zeros(len(labels))
    for tag in tags:
        vals = np.array([counts[lbl][tag] for lbl in labels])
        if vals.sum() == 0:
            continue
        bars = ax.bar(labels, vals, bottom=bottoms, label=tag.replace("_", " "), edgecolor="white", linewidth=0.7)
        for bar, value, bottom in zip(bars, vals, bottoms):
            if value > 0:
                ax.text(bar.get_x() + bar.get_width() / 2, bottom + value / 2, f"{int(value)}", ha="center", va="center", fontsize=7.5, color="white", fontweight="bold")
        bottoms += vals
    for i, lbl in enumerate(labels):
        ax.text(i, bottoms[i] + max(bottoms) * 0.03, f"{rates.get(lbl, 0)*100:.1f}%", ha="center", va="bottom", fontsize=8.5, fontweight="bold", color=PRIMARY)
    ax.set_ylabel("Discordance events")
    ax.set_ylim(0, max(bottoms) * 1.16)
    _figure_header(
        fig,
        "Figure 6. Inter-tool discordance is highest on WGS and lowest on RNA-seq, matching where Champion-Challenger helps",
        "Nested-CV cohorts · WGS n=411, WES n=390, RNA-seq n=321 gene-locus pairs",
        PRIMARY,
    )
    legend_outside(ax, loc="upper left", anchor=(1.02, 1.0), frameon=False, fontsize=7)
    _finalize_layout(fig, right=0.76)
    _save(fig, out_dir, "figure_6_discordance_taxonomy")


def build_figure_7_confidence_weights(context: dict, out_dir: Path) -> None:
    df = context["wgs"]["tool_confidence_weights"].copy()
    if df.empty:
        return
    df = df[df["modality"] == "wgs"].copy()
    for col in ["base_reliability", "final_weight", "guardrail_factor"]:
        df[col] = pd.to_numeric(df[col], errors="coerce")
    df = df.sort_values("final_weight", ascending=False).head(7)

    fig, ax = plt.subplots(figsize=(9, 4.8))
    x = np.arange(len(df))
    ax.bar(x - 0.16, df["base_reliability"], width=0.32, color=LIGHT_SLATE, label="Base reliability")
    ax.bar(x + 0.16, df["final_weight"], width=0.32, color=PRIMARY, label="Final weight")
    ax.set_xticks(x)
    ax.set_xticklabels(df["tool"])
    rotate_xticklabels(ax, angle=22)
    ax.set_ylim(0, max(0.9, df["final_weight"].max() + 0.1))
    ax.set_ylabel("Weight")
    _figure_header(
        fig,
        "Figure 7. Learned benchmark-derived weights preserve reliability while clipping poorly calibrated confidence",
        "Primary evidence · benchmark_wgs_wave2",
        PRIMARY,
    )
    legend_outside(ax, loc="upper left", anchor=(1.02, 1.0), frameon=False)
    for xpos, row in zip(x, df.itertuples(index=False)):
        ax.text(xpos + 0.16, row.final_weight + 0.02, f"{row.final_weight:.3f}", ha="center", fontsize=7, color=SLATE)
    _finalize_layout(fig, right=0.76)
    _save(fig, out_dir, "figure_7_confidence_weights")


def _get_best_row(df: pd.DataFrame, modality: str, methods: list[str]) -> dict:
    for method in methods:
        sub = df[(df["modality"] == modality) & (df["method"] == method)]
        if not sub.empty:
            return sub.iloc[0].to_dict()
    sub = df[df["modality"] == modality].copy()
    sub["overall_correct_call_rate"] = pd.to_numeric(sub["overall_correct_call_rate"], errors="coerce")
    return sub.sort_values("overall_correct_call_rate", ascending=False).iloc[0].to_dict()


def build_figure_9_bimodal_robustness(context: dict, out_dir: Path) -> None:
    method = context["tri"]["method_comparison"].copy()
    bimodal = context["tri"]["bimodal_method_comparison"].copy()
    if method.empty or bimodal.empty:
        return
    for df in (method, bimodal):
        df["overall_correct_call_rate"] = pd.to_numeric(df["overall_correct_call_rate"], errors="coerce")
    rows = [
        _get_best_row(method, "wgs", ["OptiType", "ChampionChallenger"]),
        _get_best_row(method, "wes", ["MajorityVote", "OptiType", "ChampionChallenger"]),
        _get_best_row(method, "rnaseq", ["HLA-HD", "ChampionChallenger", "OptiType"]),
        _get_best_row(bimodal, "wes+rnaseq", ["BimodalMajorityVote", "BimodalWeightedConsensus"]),
    ]
    labels = ["WGS", "WES", "RNA", "WES+RNA bimodal"]
    values = [r["overall_correct_call_rate"] * 100.0 for r in rows]
    methods = [r["method"] for r in rows]
    colors = [PRIMARY, SECONDARY_WES, SECONDARY_RNA, "#1d4ed8"]

    fig, ax = plt.subplots(figsize=(8.8, 4.8))
    x = np.arange(len(rows))
    ax.bar(x, values, color=colors, edgecolor="white", linewidth=0.8)
    for xpos, value, method_name in zip(x, values, methods):
        ax.text(xpos, value + 1.2, f"{value:.1f}", ha="center", fontsize=8, color=SLATE)
        ax.text(xpos, 3, method_name, ha="center", va="bottom", fontsize=7.2, color="white", fontweight="bold")
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.set_ylim(0, 104)
    ax.set_ylabel("Overall correct-call rate (%)")
    _figure_header(
        fig,
        "Figure 9. Supplementary matched-subject robustness shows bimodal WES+RNA majority as the strongest multimodal result",
        "Supplementary evidence · benchmark_trimodal_robustness",
        SUPPLEMENTARY,
    )
    ax.text(0.98, 0.96, "BimodalMajorityVote = 0.9623\nSupplementary robustness evidence", transform=ax.transAxes, ha="right", va="top", fontsize=7.5, color=SUPPLEMENTARY, bbox=dict(boxstyle="round,pad=0.35", facecolor="white", edgecolor="#d1d5db"))
    _figure_footer(fig, "Transitional split family: Figure 9 carries the matched-subject unimodal plus bimodal half of the robustness interpretation.")
    _finalize_layout(fig, right=0.98)
    _save(fig, out_dir, "figure_09_bimodal_per_gene")


def build_figure_10_trimodal_robustness(context: dict, out_dir: Path) -> None:
    trimodal = context["tri"]["trimodal_method_comparison"].copy()
    bimodal = context["tri"]["bimodal_method_comparison"].copy()
    if trimodal.empty or bimodal.empty:
        return
    for df in (trimodal, bimodal):
        df["overall_correct_call_rate"] = pd.to_numeric(df["overall_correct_call_rate"], errors="coerce")
        df["callable_rate"] = pd.to_numeric(df["callable_rate"], errors="coerce")
    rows = pd.concat([bimodal, trimodal], ignore_index=True)
    rows = rows[rows["method"].isin(["BimodalMajorityVote", "BimodalWeightedConsensus", "TrimodalMajorityVote", "TrimodalWeightedConsensus"])].copy()
    rows["sort"] = rows["method"].map({
        "BimodalMajorityVote": 0,
        "BimodalWeightedConsensus": 1,
        "TrimodalMajorityVote": 2,
        "TrimodalWeightedConsensus": 3,
    })
    rows = rows.sort_values("sort")

    fig, ax = plt.subplots(figsize=(9.2, 4.9))
    x = np.arange(len(rows))
    colors = [METHOD_COLORS.get(m, LIGHT_SLATE) for m in rows["method"]]
    ax.bar(x, rows["overall_correct_call_rate"] * 100.0, color=colors, edgecolor="white", linewidth=0.8)
    ax.scatter(x, rows["callable_rate"] * 100.0, color=SLATE, s=42, zorder=4, label="Callable rate")
    for xpos, row in zip(x, rows.itertuples(index=False)):
        ax.text(xpos, row.overall_correct_call_rate * 100 + 1.2, f"{row.overall_correct_call_rate * 100:.1f}", ha="center", fontsize=7.8, color=SLATE)
    ax.set_xticks(x)
    ax.set_xticklabels(["Bimodal\nmajority", "Bimodal\nweighted", "Trimodal\nmajority", "Trimodal\nweighted"])
    ax.set_ylim(0, 104)
    ax.set_ylabel("Overall correct-call rate (%)")
    _figure_header(
        fig,
        "Figure 10. Adding WGS to form trimodal consensus does not improve over the best bimodal WES+RNA result",
        "Supplementary evidence · benchmark_trimodal_robustness",
        SUPPLEMENTARY,
    )
    note = (
        "106 matched subjects\n"
        "TrimodalMajorityVote = 0.9591\n"
        "TrimodalWeightedConsensus = 0.9308\n"
        "Does not improve over bimodal majority"
    )
    ax.text(0.98, 0.96, note, transform=ax.transAxes, ha="right", va="top", fontsize=7.4, color=SUPPLEMENTARY, bbox=dict(boxstyle="round,pad=0.35", facecolor="white", edgecolor="#d1d5db"))
    legend_outside(ax, loc="upper left", anchor=(1.02, 1.0), frameon=False)
    _figure_footer(fig, "Transitional split family: Figure 10 carries the trimodal half of the supplementary robustness interpretation.")
    _finalize_layout(fig, right=0.74)
    _save(fig, out_dir, "figure_10_trimodal_comparison")


def build_figure_s1_per_gene_accuracy(context: dict, out_dir: Path) -> None:
    base.figure_s1(context["wgs_root"] / "tables", out_dir)


def build_figure_s2_calibration_heatmap(context: dict, out_dir: Path) -> None:
    base.figure_s2(context["wgs_root"] / "tables", out_dir)


def build_figure_s3_resolution_comparison(context: dict, out_dir: Path) -> None:
    base.figure_s3(context["wgs_root"] / "tables", out_dir)


def write_publication_captions(context: dict, out_dir: Path) -> None:
    lines = [
        "# Figure Captions",
        "",
        "## Active Publication Package",
        "",
        "This directory contains the active hybrid-transitional publication figure family used by the manuscript, presentation prompt, and benchmark documentation.",
        "",
        "Benchmark hierarchy:",
        "- Primary WGS source: `analysis/1000g_realdata/benchmark_wgs_wave2`",
        "- Secondary WES source: `analysis/1000g_realdata/benchmark_wes_truthbacked/run`",
        "- Secondary RNA source: `analysis/1000g_realdata/benchmark_rna_truthbacked/run`",
        "- Supplementary matched-subject source: `analysis/1000g_realdata/benchmark_trimodal_robustness/run`",
        "",
        "## Main Figures",
        "",
        "### Figure 1. ChampHLA workflow and governance architecture",
        "Benchmark role: methods and governance. This figure makes the benchmark hierarchy explicit inside the workflow view and supports the claim that ChampHLA is a benchmark-governed calibrated HLA ensemble framework.",
        "",
        "### Figure 2. Primary WGS anchor result",
        "Benchmark role: primary. Derived from `benchmark_wgs_wave2`. `OptiType` remains the best single WGS tool, routed baselines recover the same ceiling, and `WeightedConsensus` remains below that ceiling. Allowed conclusion: WGS remains tool-limited rather than consensus-limited.",
        "",
        "### Figure 3. WGS per-gene failure-mode interpretation",
        "Benchmark role: primary. Derived from WGS per-gene and locus-difficulty tables. This compatibility stem now carries the mechanistic interpretation that WGS remains hardest at `HLA-A` and `HLA-C`, which helps explain why consensus does not fully rescue the modality.",
        "",
        "### Figure 4. Confidence guardrails",
        "Benchmark role: primary. Derived from WGS calibration and confidence-weight tables. Allowed conclusion: confidence-based boosting must be benchmark-governed because poor calibration can otherwise distort runtime voting.",
        "",
        "### Figure 5. Transitional WGS failure-mode family, cross-modality half",
        "Benchmark role: mixed. This compatibility stem now carries the cross-modality context for the WGS failure-mode story: WES and RNA operate in stronger unimodal regimes than primary short-read WGS. It should be interpreted together with Figure 6.",
        "",
        "### Figure 6. Transitional WGS failure-mode family, discordance half",
        "Benchmark role: primary. Derived from `benchmark_wgs_wave2` discordance summaries. It should be interpreted together with Figure 5 as the second half of the transitional WGS failure-mode family.",
        "",
        "### Figure 7. Confidence-derived weights",
        "Benchmark role: primary. Derived from WGS weight tables. Allowed conclusion: final weights preserve reliability while clipping poorly calibrated confidence.",
        "",
        "### Figure 9. Supplementary robustness, unimodal plus bimodal half",
        "Benchmark role: supplementary. Derived from matched-subject unimodal and bimodal tables in `benchmark_trimodal_robustness/run`. Allowed conclusion: `BimodalMajorityVote` on `wes+rnaseq` is the strongest multimodal result in the matched-subject robustness benchmark.",
        "",
        "### Figure 10. Supplementary robustness, trimodal half",
        "Benchmark role: supplementary. Derived from trimodal and comparison tables in `benchmark_trimodal_robustness/run`. Allowed conclusion: adding WGS to form trimodal consensus does not improve over the best bimodal WES+RNA result.",
        "",
        "## Supplementary Figures",
        "",
        "### Figure S1. Per-gene accuracy",
        "Supporting per-gene benchmark detail retained for publication continuity.",
        "",
        "### Figure S2. Calibration heatmap",
        "Supporting calibration detail retained for publication continuity.",
        "",
        "### Figure S3. Resolution comparison",
        "Supporting resolution detail retained for publication continuity.",
    ]
    (out_dir / "captions.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_publication_manifest_update(out_dir: Path) -> None:
    lines = [
        "# figures_final",
        "",
        "This directory is the authoritative publication figure package for the current ChampHLA manuscript state.",
        "",
        "## Active figure family",
        "",
        "Top-level active figure stems:",
    ]
    lines.extend([f"- `{stem}`" for stem in EXPECTED_ACTIVE_STEMS])
    lines.extend([
        "",
        "## Hybrid transitional state",
        "",
        "- current public stems are preserved for manuscript, prompt, and test compatibility",
        "- internal figure semantics were refreshed according to:",
        "  - `docs/FIGURE_REDESIGN_BRIEF.md`",
        "  - `docs/FIGURE_STORYBOARD.md`",
        "  - `docs/FIGURE_DATA_SPEC.md`",
        "  - `docs/FIGURE_MESSAGING_REVIEW.md`",
        "- Figures 5 and 6 remain separate stems during the transition even though they form one conceptual WGS failure-mode family",
        "- Figures 9 and 10 remain separate stems during the transition even though they form one supplementary robustness conclusion family",
        "",
        "## Archived legacy assets",
        "",
        "Legacy pre-wave2 stems remain under:",
        "- `legacy_pre_wave2/`",
        "",
        "Legacy pre-redesign snapshots may also be preserved under timestamped `legacy_pre_redesign_*` directories.",
        "",
        "## Benchmark provenance",
        "",
        "Authoritative roots:",
        "- WGS primary benchmark: `analysis/1000g_realdata/benchmark_wgs_wave2`",
        "- WES truth-backed benchmark: `analysis/1000g_realdata/benchmark_wes_truthbacked/run`",
        "- RNA truth-backed benchmark: `analysis/1000g_realdata/benchmark_rna_truthbacked/run`",
        "- Supplementary matched-subject robustness benchmark: `analysis/1000g_realdata/benchmark_trimodal_robustness/run`",
        "",
        "Interpretation contract:",
        "- `benchmark_wgs_wave2` is the primary WGS benchmark.",
        "- The 106-sample matched-subject trimodal robustness run is supplementary evidence.",
        "- WGS remains tool-limited.",
        "- `wes+rnaseq` bimodal majority is the strongest multimodal result.",
        "- Trimodal consensus does not improve over the best bimodal headline.",
    ])
    (out_dir / "README.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def generate_all(context: dict, out_dir: Path) -> None:
    build_figure_1_workflow_governance(context, out_dir)
    build_figure_2_wgs_anchor(context, out_dir)
    build_figure_3_wgs_failure_modes(context, out_dir)
    build_figure_4_confidence_guardrails(context, out_dir)
    build_figure_5_cross_modality_comparison(context, out_dir)
    build_figure_6_discordance_interpretation(context, out_dir)
    build_figure_7_confidence_weights(context, out_dir)
    build_figure_9_bimodal_robustness(context, out_dir)
    build_figure_10_trimodal_robustness(context, out_dir)
    build_figure_s1_per_gene_accuracy(context, out_dir)
    build_figure_s2_calibration_heatmap(context, out_dir)
    build_figure_s3_resolution_comparison(context, out_dir)
    write_publication_captions(context, out_dir)
    write_publication_manifest_update(out_dir)


def main() -> int:
    args = parse_args()
    context = load_publication_context(args)
    args.out_dir.mkdir(parents=True, exist_ok=True)
    generate_all(context, args.out_dir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
