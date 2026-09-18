#!/usr/bin/env python3.11
"""Champion-Challenger vs Majority-Voting explainer figures for the ChampHLA (MVHLA) manuscript.

Generates three layout variants from authoritative WES truth-backed benchmark tables so the
user can pick one for the manuscript:

  figure_11_champion_challenger_combined -- one multi-panel figure (mechanism + MV-vs-CC contrast)
  figure_11a_cc_mechanism                -- the Champion-Challenger mechanism schematic only
  figure_11b_mv_vs_cc                    -- the MV-vs-CC contrast + real-number accuracy/audit

Real numbers are read from:
  <wes-tables>/method_comparison.tsv                 (OptiType / MajorityVote / WeightedConsensus)
  <wes-tables>/../sweeps/wes_champion_override_sweep.tsv  (tuned Champion-Challenger optimum + audit)

If a table is missing, documented constants from CHAMPHLA_MANUSCRIPT_V1.md are used as fallback.
Style (Arial, Wong palette, 300 DPI, legends outside axes) comes from figure_hub/plot_style.py.
"""

import argparse
from pathlib import Path
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
import numpy as np
import pandas as pd

SCRIPT_DIR = Path(__file__).resolve().parent
HUB_ROOT = SCRIPT_DIR.parent
sys.path.insert(0, str(HUB_ROOT))

from plot_style import (  # noqa: E402
    METHOD_COLORS,
    apply_style,
    save_fig,
    set_figure_footer as _figure_footer,
    set_figure_header as _figure_header,
)

apply_style()

# Palette (mirrors fig_benchmark_main.py local constants)
PRIMARY = "#16324f"
SLATE = "#475569"
LIGHT_SLATE = "#94a3b8"
INK = "#1f2937"
GREEN = "#0f766e"
RED = "#b42318"
AMBER = "#b45309"
MV_COLOR = METHOD_COLORS["MajorityVote"]            # #b45309
CC_COLOR = METHOD_COLORS["ChampionChallenger"]      # #2d728f
OPTI_COLOR = LIGHT_SLATE

# Nested-CV held-out numbers (10-fold, out-of-fold operating points). WGS is the
# headline (significant), WES/RNA match majority voting. Fallbacks only; live values
# are loaded from the nested-CV tables in load_numbers().
ANALYSIS = Path("/scratch/project_2008084/pihla-publish/analysis")
FALLBACK_BY_MODALITY = {
    "wgs": dict(mv=0.3844, cc=0.5012, cc_lo=0.4531, cc_hi=0.5493, best=0.4975, best_name="OptiType",
                delta=0.1168, mcnemar_p=1e-4, b=60, c=12, n=411),
    "wes": dict(mv=0.9359, cc=0.9436, cc_lo=0.9161, cc_hi=0.9625, best=0.9231, best_name="OptiType",
                delta=0.0077, mcnemar_p=0.4531, b=5, c=2, n=390),
    "rna": dict(mv=0.9502, cc=0.9408, cc_lo=0.9094, cc_hi=0.9618, best=0.9429, best_name="HLA-HD",
                delta=-0.0093, mcnemar_p=0.4531, b=2, c=5, n=321),
}
# WGS per-gene gains (for the per-gene panel) and champions (for the mechanism schematic).
FALLBACK_WGS_PER_GENE = {"A": (0.3577, 0.4818, "p<0.001"), "B": (0.4964, 0.6934, "p<0.001"), "C": (0.2993, 0.3285, "n.s.")}
WGS_CHAMPIONS = {"A": "OptiType", "B": "OptiType", "C": "T1K"}


# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------
def parse_args() -> argparse.Namespace:
    analysis_root = Path(
        "/scratch/project_2008084/pihla-publish/analysis/1000g_realdata"
    )
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--wes-tables",
        type=Path,
        default=analysis_root / "benchmark_wes_truthbacked" / "run" / "tables",
        help="dir containing method_comparison.tsv (and ../sweeps/wes_champion_override_sweep.tsv)",
    )
    parser.add_argument(
        "--sweep",
        type=Path,
        default=None,
        help="explicit path to wes_champion_override_sweep.tsv (default: <wes-tables>/../sweeps/...)",
    )
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=Path("/scratch/project_2008084/pihla-publish/analysis/figures_final_candidate"),
    )
    return parser.parse_args()


def _maybe_read_tsv(path: Path) -> pd.DataFrame:
    if path and path.exists():
        return pd.read_csv(path, sep="\t")
    return pd.DataFrame()


def _tool_best(mc: pd.DataFrame, modality: str):
    tools = mc[(mc["modality"] == modality) & (mc["method_type"] == "single_tool")]
    if tools.empty:
        return None, None
    tools = tools.sort_values("overall_correct_call_rate", ascending=False)
    top = tools.iloc[0]
    return float(top["overall_correct_call_rate"]), str(top["method"])


def load_numbers(args: argparse.Namespace) -> dict:
    """Load per-modality nested-CV held-out numbers (MV, Champion-Challenger, best tool,
    McNemar) plus the WGS per-gene gains, preferring the on-disk tables over fallbacks."""
    data = {"by_modality": {m: dict(v) for m, v in FALLBACK_BY_MODALITY.items()},
            "wgs_per_gene": dict(FALLBACK_WGS_PER_GENE),
            "champion": dict(WGS_CHAMPIONS), "_source": "fallback constants"}

    for lbl in ("wgs", "wes", "rna"):
        mod = "rnaseq" if lbl == "rna" else lbl
        cv = ANALYSIS / f"benchmark_{lbl}_cv_recalibrated" / "tables"
        ncv = ANALYSIS / "nested_cv_champion_challenger" / lbl
        mc = _maybe_read_tsv(cv / "method_comparison.tsv")
        cc = _maybe_read_tsv(cv / "champion_challenger_method_comparison.tsv")
        mci = _maybe_read_tsv(ncv / "nested_cv_mcnemar.tsv")
        d = data["by_modality"][lbl]
        if not mc.empty:
            mvrow = mc[(mc["modality"] == mod) & (mc["method"] == "MajorityVote")]
            if not mvrow.empty:
                d["mv"] = round(float(mvrow.iloc[0]["overall_correct_call_rate"]), 4)
                d["n"] = int(mvrow.iloc[0]["gene_rows"])
            best, best_name = _tool_best(mc, mod)
            if best is not None:
                d["best"], d["best_name"] = round(best, 4), best_name
        if not cc.empty:
            r = cc.iloc[0]
            d["cc"] = round(float(r["overall_correct_call_rate"]), 4)
            d["cc_lo"] = float(r["overall_correct_call_rate_ci_lo"])
            d["cc_hi"] = float(r["overall_correct_call_rate_ci_hi"])
        if not mci.empty:
            r = mci.iloc[0]
            d["mcnemar_p"] = float(r["p_value"]); d["b"] = int(r["b"]); d["c"] = int(r["c"])
        d["delta"] = round(d["cc"] - d["mv"], 4)

    # WGS per-gene (MV vs CC)
    cv = ANALYSIS / "benchmark_wgs_cv_recalibrated" / "tables"
    mvpg = _maybe_read_tsv(cv / "method_per_gene.tsv")
    ccpg = _maybe_read_tsv(cv / "champion_challenger_method_per_gene.tsv")
    if not mvpg.empty and not ccpg.empty:
        sig = {"A": "p<0.001", "B": "p<0.001", "C": "n.s."}
        for g in ("A", "B", "C"):
            mv_r = mvpg[(mvpg["modality"] == "wgs") & (mvpg["method"] == "MajorityVote") & (mvpg["gene"] == g)]
            cc_r = ccpg[(ccpg["modality"] == "wgs") & (ccpg["gene"] == g)]
            if not mv_r.empty and not cc_r.empty:
                data["wgs_per_gene"][g] = (round(float(mv_r.iloc[0]["overall_correct_call_rate"]), 4),
                                           round(float(cc_r.iloc[0]["overall_correct_call_rate"]), 4), sig[g])
        data["_source"] = str(cv)
    return data


# ---------------------------------------------------------------------------
# Schematic primitives (axes use normalised 0..1 coordinates)
# ---------------------------------------------------------------------------
def _box(ax, x, y, w, h, label, sub="", fc="#eef2f7", ec="#b7c0ca", text=INK,
         lw=0.9, fs=8.8, sub_fs=7.0, sub_color=SLATE):
    ax.add_patch(FancyBboxPatch(
        (x, y), w, h, boxstyle="round,pad=0.012,rounding_size=0.02",
        facecolor=fc, edgecolor=ec, linewidth=lw, transform=ax.transData,
        clip_on=False, mutation_aspect=0.55))
    ly = y + h * (0.60 if sub else 0.5)
    ax.text(x + w / 2, ly, label, ha="center", va="center",
            fontsize=fs, fontweight="bold", color=text, zorder=5)
    if sub:
        ax.text(x + w / 2, y + h * 0.27, sub, ha="center", va="center",
                fontsize=sub_fs, color=sub_color, zorder=5)


def _arrow(ax, x1, y1, x2, y2, color="#7c8793", lw=1.4, style="-|>"):
    ax.annotate("", xy=(x2, y2), xytext=(x1, y1),
                arrowprops=dict(arrowstyle=style, color=color, lw=lw),
                annotation_clip=False)


# ---------------------------------------------------------------------------
# Panel: Champion-Challenger mechanism
# ---------------------------------------------------------------------------
def draw_mechanism(ax, d: dict) -> None:
    champ = d.get("champion", {"A": "OptiType", "B": "OptiType", "C": "T1K"})
    if len(set(champ.values())) == 1:
        champ_sub = f"{next(iter(champ.values()))} for HLA-A/-B/-C\nhighest training accuracy → default call"
    else:
        champ_sub = "  ".join(f"{g}:{t}" for g, t in champ.items()) + "\n(per gene; WGS shown) → default call"

    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")

    # --- top pipeline ---
    _box(ax, 0.020, 0.70, 0.205, 0.20, "Per-gene champion", champ_sub,
         fc="#dbe7f7", ec="#9bb4d6", text=PRIMARY)
    _arrow(ax, 0.225, 0.80, 0.275, 0.80)
    _box(ax, 0.275, 0.70, 0.235, 0.20, "Weighted ensemble vote",
         "w = 0.7·reliability + 0.3·confidence\nguardrail: Brier & ECE < 0.35",
         fc="#fbf3d6", ec="#d9c98a", text=INK)
    _arrow(ax, 0.510, 0.80, 0.560, 0.80)
    _box(ax, 0.560, 0.70, 0.205, 0.20, "Challenger pair\n≠ champion?",
         "ensemble's top-weighted\nallele pair", fc="#ece8fb", ec="#c3b8ec", text="#4c1d95")

    # challenger differs -> down into the gate test
    _arrow(ax, 0.620, 0.70, 0.580, 0.615, color="#7c8793")
    ax.text(0.628, 0.660, "yes, differs", fontsize=7, color=SLATE, ha="left", va="center")
    # challenger == champion -> champion retained (clean right-side drop, clears gate box)
    _arrow(ax, 0.748, 0.70, 0.727, 0.325, color=PRIMARY, lw=1.2)
    ax.text(0.752, 0.52, "no", fontsize=7, color=SLATE, ha="left", va="center")
    # champion uncallable -> fallback (clean left-side drop, clears gate box)
    _arrow(ax, 0.120, 0.70, 0.120, 0.325, color=AMBER, lw=1.2)
    ax.text(0.135, 0.52, "champion\nuncallable", fontsize=7, color=AMBER, ha="left", va="center")

    # --- four-gate override test (centred container) ---
    ax.text(0.550, 0.640, "Override test — all 4 gates must pass",
            fontsize=8.4, fontweight="bold", color=PRIMARY, ha="center")
    gates = [
        "1.  Support fraction  ≥ min",
        "2.  Weight margin  ≥ min",
        "3.  Supporting tools  ≥ min",
        "4.  Challenger allele unambiguous",
    ]
    gy = 0.560
    for g in gates:
        _box(ax, 0.395, gy, 0.310, 0.048, g, fc="#f3f5f8", ec="#c4ccd6",
             text=INK, fs=7.8)
        gy -= 0.060
    # gate cluster -> outcomes
    _arrow(ax, 0.500, 0.375, 0.455, 0.325, color=GREEN, lw=1.6)
    ax.text(0.440, 0.350, "all pass", fontsize=7, color=GREEN, ha="right", va="center")
    _arrow(ax, 0.640, 0.375, 0.690, 0.325, color=RED, lw=1.4)
    ax.text(0.655, 0.350, "any fails", fontsize=7, color=RED, ha="left", va="center")

    # --- outcomes (decision-trace categories), no crossing arrows ---
    _box(ax, 0.020, 0.150, 0.230, 0.175, "champion_missing",
         "champion uncallable →\nweighted-consensus fallback",
         fc="#fdebd8", ec="#e6bd86", text=AMBER, fs=8.4)
    _box(ax, 0.300, 0.150, 0.250, 0.175, "challenger_override",
         "adopt consensus allele pair", fc="#dff3e8", ec="#86c7a6", text=GREEN, fs=8.6)
    _box(ax, 0.610, 0.150, 0.250, 0.175, "champion_retained",
         "keep champion's call", fc="#dbe7f7", ec="#9bb4d6", text=PRIMARY, fs=8.6)

    # --- majority-voting contrast strip ---
    ax.add_patch(FancyBboxPatch(
        (0.010, 0.010, ), 0.980, 0.085, boxstyle="round,pad=0.006,rounding_size=0.01",
        facecolor="#fbeede", edgecolor=MV_COLOR, linewidth=1.0, alpha=0.65,
        transform=ax.transData, clip_on=False, mutation_aspect=0.25))
    ax.text(0.5, 0.052,
            "Majority Voting (baseline):  every callable tool casts 1 equal vote  →  most-voted "
            "allele pair wins  →  tie = no_call.   No champion, no reliability weights, no calibration.",
            ha="center", va="center", fontsize=8.0, color="#7a3d00", fontweight="bold")


# ---------------------------------------------------------------------------
# Panel: MV vs CC mini-flows
# ---------------------------------------------------------------------------
def draw_mv_cc_flows(ax, d: dict) -> None:
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    ax.text(0.5, 0.97, "Decision logic", fontsize=9.5, fontweight="bold",
            color=PRIMARY, ha="center")

    # Majority Voting column (left)
    ax.text(0.25, 0.90, "Majority Voting", fontsize=9, fontweight="bold",
            color=MV_COLOR, ha="center")
    mv_steps = [
        ("All callable tools", "8 tools, harmonised calls"),
        ("Equal vote each", "1 tool = 1 vote"),
        ("Most-voted pair", "tie → no_call"),
    ]
    y = 0.74
    for i, (lab, sub) in enumerate(mv_steps):
        _box(ax, 0.05, y, 0.40, 0.13, lab, sub, fc="#fbeede", ec="#e6bd86",
             text="#7a3d00", fs=8.4)
        if i < len(mv_steps) - 1:
            _arrow(ax, 0.25, y, 0.25, y - 0.06, color=MV_COLOR)
        y -= 0.235

    # Champion-Challenger column (right)
    ax.text(0.75, 0.90, "Champion-Challenger", fontsize=9, fontweight="bold",
            color=CC_COLOR, ha="center")
    cc_steps = [
        ("Reliability-weighted vote", "per-gene champion = default"),
        ("Gated override test", "4 simultaneous criteria"),
        ("Audited consensus", "retain / override / fallback"),
    ]
    y = 0.74
    for i, (lab, sub) in enumerate(cc_steps):
        _box(ax, 0.55, y, 0.40, 0.13, lab, sub, fc="#dceef3", ec="#9cc6d4",
             text=CC_COLOR, fs=8.4)
        if i < len(cc_steps) - 1:
            _arrow(ax, 0.75, y, 0.75, y - 0.06, color=CC_COLOR)
        y -= 0.235

    ax.text(0.5, 0.045,
            "Symmetric & unweighted  vs.  routed, weight-gated & auditable",
            ha="center", va="center", fontsize=7.6, color=SLATE, style="italic")


# ---------------------------------------------------------------------------
# Panel: accuracy bars
# ---------------------------------------------------------------------------
def draw_accuracy_bars(ax, d: dict) -> None:
    """Champion-Challenger vs Majority Voting across the three modalities (nested CV).
    WGS is the significant win; WES/RNA match MV."""
    bm = d["by_modality"]
    order = [("wgs", "WGS\n(n=411)"), ("wes", "WES\n(n=390)"), ("rna", "RNA-seq\n(n=321)")]
    x = np.arange(len(order))
    w = 0.38
    mv = [bm[k]["mv"] * 100 for k, _ in order]
    cc = [bm[k]["cc"] * 100 for k, _ in order]
    ax.bar(x - w / 2, mv, w, label="Majority Voting", color=MV_COLOR, edgecolor="white", zorder=3)
    ax.bar(x + w / 2, cc, w, label="Champion-Challenger", color=CC_COLOR, edgecolor="white", zorder=3)
    for i, (k, _) in enumerate(order):
        ax.text(x[i] - w / 2, mv[i] + 1.5, f"{mv[i]:.1f}", ha="center", va="bottom", fontsize=7.2, color="#7a3d00")
        ax.text(x[i] + w / 2, cc[i] + 1.5, f"{cc[i]:.1f}", ha="center", va="bottom", fontsize=7.2, color=CC_COLOR)
        p = bm[k]["mcnemar_p"]
        tag = "***" if p < 0.001 else ("*" if p < 0.05 else "n.s.")
        ax.text(x[i], max(mv[i], cc[i]) + 7.5, f"Δ{(cc[i]-mv[i]):+.1f}\n{tag}", ha="center", va="bottom",
                fontsize=7.6, fontweight="bold", color=GREEN if p < 0.05 else SLATE)
    ax.set_xticks(x)
    ax.set_xticklabels([lbl for _, lbl in order], fontsize=8)
    ax.set_ylim(0, 108)
    ax.set_ylabel("Overall correct-call rate (%)", fontsize=8.5)
    ax.grid(axis="y", alpha=0.28)
    ax.legend(loc="upper left", fontsize=7.6, frameon=False)
    ax.set_title("Accuracy vs majority voting (10-fold nested CV)", fontsize=9.5, color=PRIMARY)


# ---------------------------------------------------------------------------
# Panel: override audit
# ---------------------------------------------------------------------------
def draw_override_audit(ax, d: dict) -> None:
    """Per-gene WGS accuracy (MV vs CC): the modality win concentrates at HLA-A and -B,
    where champion routing recovers the per-gene single-tool ceiling."""
    pg = d["wgs_per_gene"]
    genes = ["A", "B", "C"]
    x = np.arange(len(genes))
    w = 0.38
    mv = [pg[g][0] * 100 for g in genes]
    cc = [pg[g][1] * 100 for g in genes]
    ax.bar(x - w / 2, mv, w, color=MV_COLOR, edgecolor="white", zorder=3)
    ax.bar(x + w / 2, cc, w, color=CC_COLOR, edgecolor="white", zorder=3)
    for i, g in enumerate(genes):
        ax.text(x[i], max(mv[i], cc[i]) + 3, f"+{cc[i]-mv[i]:.0f}\n{pg[g][2]}", ha="center", va="bottom",
                fontsize=7.0, fontweight="bold", color=GREEN if pg[g][2] != "n.s." else SLATE)
    ax.set_xticks(x)
    ax.set_xticklabels([f"HLA-{g}" for g in genes], fontsize=8)
    ax.set_ylim(0, 92)
    ax.set_ylabel("WGS correct-call (%)", fontsize=8.5)
    ax.grid(axis="y", alpha=0.28)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    ax.set_title("WGS per-gene: routing wins at HLA-A & -B", fontsize=9, color=PRIMARY)


# ---------------------------------------------------------------------------
# Figure builders
# ---------------------------------------------------------------------------
def _save(fig, out_dir: Path, stem: str) -> None:
    save_fig(fig, out_dir / stem)


def build_figure_11a(d: dict, out_dir: Path) -> None:
    fig, ax = plt.subplots(figsize=(13.0, 6.4))
    _figure_header(
        fig,
        "Figure 11. ChampHLA Champion-Challenger consensus mechanism",
        "A benchmark champion is the default; the weighted ensemble overrides it only when all four "
        "evidence gates pass",
        PRIMARY, title_y=0.975, role_y=0.93)
    draw_mechanism(ax, d)
    _figure_footer(
        fig,
        "Per-gene champions and override thresholds selected out of fold under nested cross-validation "
        "(WGS: HLA-A/-B = OptiType, HLA-C = T1K).",
        y=0.02)
    fig.subplots_adjust(left=0.01, right=0.99, top=0.90, bottom=0.06)
    _save(fig, out_dir, "figure_11a_cc_mechanism")


def build_figure_11b(d: dict, out_dir: Path) -> None:
    fig = plt.figure(figsize=(13.0, 6.0))
    gs = fig.add_gridspec(
        2, 2, width_ratios=[1.0, 1.15], height_ratios=[1.25, 1.0],
        wspace=0.28, hspace=0.55, left=0.06, right=0.97, top=0.86, bottom=0.12)
    ax_flow = fig.add_subplot(gs[:, 0])
    ax_bars = fig.add_subplot(gs[0, 1])
    ax_audit = fig.add_subplot(gs[1, 1])

    _figure_header(
        fig,
        "Figure 11. Majority Voting vs Champion-Challenger consensus",
        "Champion-Challenger significantly beats majority voting on the discordant WGS modality and "
        "matches it on WES/RNA (10-fold nested CV)",
        PRIMARY, title_y=0.975, role_y=0.925)
    draw_mv_cc_flows(ax_flow, d)
    draw_accuracy_bars(ax_bars, d)
    draw_override_audit(ax_audit, d)
    _figure_footer(
        fig,
        "10-fold nested cross-validation; operating points selected out of fold, scores pooled over "
        "held-out folds. *** exact McNemar p<0.001.",
        y=0.02)
    _save(fig, out_dir, "figure_11b_mv_vs_cc")


def build_figure_11_combined(d: dict, out_dir: Path) -> None:
    fig = plt.figure(figsize=(14.5, 10.2))
    outer = fig.add_gridspec(
        2, 1, height_ratios=[1.05, 1.0], hspace=0.20,
        left=0.04, right=0.97, top=0.92, bottom=0.06)
    ax_mech = fig.add_subplot(outer[0])
    lower = outer[1].subgridspec(
        2, 2, width_ratios=[1.0, 1.15], height_ratios=[1.25, 1.0],
        wspace=0.28, hspace=0.6)
    ax_flow = fig.add_subplot(lower[:, 0])
    ax_bars = fig.add_subplot(lower[0, 1])
    ax_audit = fig.add_subplot(lower[1, 1])

    _figure_header(
        fig,
        "Figure 11. How ChampHLA's Champion-Challenger consensus works and how it differs from majority voting",
        "(A) per-gene champion with a four-gate ensemble override   ·   (B) contrast with majority "
        "voting and the nested-CV accuracy outcome across modalities",
        PRIMARY, title_y=0.982, role_y=0.952)
    draw_mechanism(ax_mech, d)
    ax_mech.text(-0.005, 1.02, "A", transform=ax_mech.transAxes,
                 fontsize=14, fontweight="bold", color=PRIMARY, va="bottom")
    draw_mv_cc_flows(ax_flow, d)
    ax_flow.text(-0.02, 1.04, "B", transform=ax_flow.transAxes,
                 fontsize=14, fontweight="bold", color=PRIMARY, va="bottom")
    draw_accuracy_bars(ax_bars, d)
    draw_override_audit(ax_audit, d)
    _figure_footer(
        fig,
        "10-fold nested cross-validation; per-gene champions and override thresholds selected out of "
        "fold (WGS: HLA-A/-B = OptiType, HLA-C = T1K). *** exact McNemar p<0.001.",
        y=0.015)
    _save(fig, out_dir, "figure_11_champion_challenger_combined")


def write_captions(d: dict, out_dir: Path) -> None:
    bm = d["by_modality"]
    w, e, r = bm["wgs"], bm["wes"], bm["rna"]
    lines = [
        "# Figure 11 captions — Champion-Challenger vs Majority Voting",
        "",
        f"_Data source: {d.get('_source')} (10-fold nested cross-validation)_",
        "",
        "**Figure 11 (combined).** How ChampHLA's Champion-Challenger consensus works and how it "
        "differs from majority voting. (A) For each HLA gene a benchmark-designated champion "
        "(selected per gene and modality out of fold; WGS: HLA-A/-B = OptiType, HLA-C = T1K) provides "
        "the default call; a reliability- and confidence-weighted ensemble (w = 0.7·reliability + "
        "0.3·confidence, gated by a Brier/ECE < 0.35 calibration guardrail) may override it only when "
        "all four gates pass: minimum support fraction, weight margin, supporting tools (thresholds "
        "selected out of fold), and an unambiguous challenger allele. Every locus is recorded as "
        "champion_retained, challenger_override, or champion_missing (weighted-consensus fallback). "
        "Majority voting, by contrast, treats every callable tool as one equal vote with ties returning "
        "no_call. (B) Under 10-fold nested cross-validation, Champion-Challenger significantly beats "
        f"majority voting on the discordant WGS modality ({w['cc']:.4f} vs {w['mv']:.4f}; +{w['delta']*100:.1f} "
        f"points; exact McNemar p<0.001) and the best single tool ({w['best_name']} {w['best']:.4f}), while "
        f"matching majority voting on the high-accuracy WES ({e['cc']:.4f} vs {e['mv']:.4f}) and RNA-seq "
        f"({r['cc']:.4f} vs {r['mv']:.4f}) modalities (both p=0.45). The WGS gain concentrates at HLA-A "
        "and HLA-B and arises from champion routing (the override gate does not fire on WGS).",
        "",
        "**Figure 11a.** Champion-Challenger mechanism schematic (panel A above, standalone).",
        "",
        "**Figure 11b.** Majority Voting vs Champion-Challenger contrast: nested-CV accuracy across "
        "modalities and the WGS per-gene breakdown (panel B above, standalone).",
        "",
    ]
    (out_dir / "captions_champion_challenger.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    args = parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    d = load_numbers(args)
    print(f"[fig11] data source: {d.get('_source')}")
    for m in ("wgs", "wes", "rna"):
        b = d["by_modality"][m]
        print(f"[fig11] {m}: MV={b['mv']} CC={b['cc']} best={b.get('best_name')} {b.get('best')} "
              f"delta={b['delta']} p={b['mcnemar_p']}")
    build_figure_11a(d, args.out_dir)
    build_figure_11b(d, args.out_dir)
    build_figure_11_combined(d, args.out_dir)
    write_captions(d, args.out_dir)
    print(f"[fig11] wrote figures + captions to {args.out_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
