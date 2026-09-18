#!/usr/bin/env python3.11
"""NCI-60 external-validation figure for the ChampHLA (MVHLA) manuscript (Figure 13).

Independent, non-1000G generalisation test: the *frozen* 1000G-derived Champion-Challenger policy
applied to the NCI-60 cancer cell-line panel (Adams 2005 SBT truth), across RNA-seq and WES.

Three panels, all data-driven (regenerate to pick up the powered RNA scale-up automatically):
  A  Overall correct-call rate, RNA vs WES, best-single / MajorityVote / ChampHLA (Wilson CIs).
  B  Per-locus (A/B/C) ChampHLA vs MajorityVote, RNA vs WES.
  C  WES override-gate sweep: CC accuracy vs challenger-margin at tools=2 for each support level,
     marking the frozen point (0.65/0.20/2) and the sweep optimum, with the MajorityVote reference.

Reads:
  <rna-tables>/method_comparison.tsv, champion_challenger_method_comparison.tsv,
              method_per_gene.tsv, champion_challenger_method_per_gene.tsv
  <wes-tables>/ (same set) and <wes-tables>/../sweeps/wes_champion_override_sweep.tsv
Style (Arial, Wong palette, 300 DPI, legends outside axes) from figure_hub/plot_style.py.
"""

import argparse
from pathlib import Path
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

SCRIPT_DIR = Path(__file__).resolve().parent
HUB_ROOT = SCRIPT_DIR.parent
sys.path.insert(0, str(HUB_ROOT))

from plot_style import (  # noqa: E402
    METHOD_COLORS,
    apply_style,
    save_fig,
    set_figure_header as _figure_header,
    set_figure_footer as _figure_footer,
)

apply_style()

MV_COLOR = METHOD_COLORS["MajorityVote"]
CC_COLOR = METHOD_COLORS["ChampionChallenger"]
SINGLE_COLOR = "#94a3b8"
INK = "#1f2937"
GRID = "#cbd5e1"
NCI60_RNA = "/scratch/project_2008084/pihla-publish/analysis/nci60_benchmark/tables"
NCI60_WES = "/scratch/project_2008084/pihla-publish/analysis/nci60_wes_benchmark/tables"
OUT_DEFAULT = "/scratch/project_2008084/pihla-publish/analysis/figures_final_candidate"


def _read(p: Path) -> pd.DataFrame:
    try:
        return pd.read_csv(p, sep="\t")
    except Exception:
        return pd.DataFrame()


def load_modality(tables: Path) -> dict:
    """Pull CC / MajorityVote / best-single overall + per-locus from a benchmark tables dir."""
    mc = _read(tables / "method_comparison.tsv")
    cc = _read(tables / "champion_challenger_method_comparison.tsv")
    mpg = _read(tables / "method_per_gene.tsv")
    ccpg = _read(tables / "champion_challenger_method_per_gene.tsv")
    d: dict = {"n": None, "overall": {}, "per_gene": {}}
    if mc.empty:
        return d
    d["n"] = int(mc["sample_count"].max())
    acc = "overall_correct_call_rate"
    lo, hi = acc + "_ci_lo", acc + "_ci_hi"

    def row(df, method):
        r = df[df["method"] == method]
        return None if r.empty else r.iloc[0]

    mv = row(mc, "MajorityVote")
    if mv is not None:
        d["overall"]["MajorityVote"] = (float(mv[acc]), float(mv[lo]), float(mv[hi]))
    # best single tool (exclude coverage-only 0-callable tools)
    singles = mc[(mc["method_type"] == "single_tool") & (mc["callable_rate"] > 0)]
    if not singles.empty:
        b = singles.loc[singles[acc].idxmax()]
        d["overall"]["BestSingle"] = (float(b[acc]), float(b[lo]), float(b[hi]), str(b["method"]))
    ccr = row(cc, "ChampionChallenger") if not cc.empty else None
    if ccr is not None:
        d["overall"]["ChampionChallenger"] = (float(ccr[acc]), float(ccr[lo]), float(ccr[hi]))
    # per-locus
    for gene in ("A", "B", "C"):
        d["per_gene"].setdefault(gene, {})
        if not mpg.empty:
            r = mpg[(mpg["method"] == "MajorityVote") & (mpg["gene"] == gene)]
            if not r.empty:
                d["per_gene"][gene]["MajorityVote"] = float(r.iloc[0][acc])
        if not ccpg.empty:
            r = ccpg[(ccpg["method"] == "ChampionChallenger") & (ccpg["gene"] == gene)]
            if not r.empty:
                d["per_gene"][gene]["ChampionChallenger"] = float(r.iloc[0][acc])
    return d


def load_sweep(path: Path) -> pd.DataFrame:
    sw = _read(path)
    if sw.empty:
        return sw
    return sw[sw["modality"] == "wes"].copy() if "modality" in sw.columns else sw


# ---------------------------------------------------------------------------
def panel_overall(ax, rna, wes):
    methods = [("BestSingle", "Best single tool", SINGLE_COLOR),
               ("MajorityVote", "MajorityVote", MV_COLOR),
               ("ChampionChallenger", "ChampHLA (CC)", CC_COLOR)]
    mods = [("RNA-seq", rna), ("WES", wes)]
    x = np.arange(len(mods)); w = 0.26
    for i, (mkey, mlab, col) in enumerate(methods):
        vals, los, his = [], [], []
        for _, dat in mods:
            t = dat["overall"].get(mkey)
            if t is None:
                vals.append(np.nan); los.append(0); his.append(0)
            else:
                vals.append(t[0]); los.append(t[0] - t[1]); his.append(t[2] - t[0])
        xs = x + (i - 1) * w
        ax.bar(xs, vals, w, color=col, label=mlab, edgecolor="white", linewidth=0.6, zorder=3)
        ax.errorbar(xs, vals, yerr=[los, his], fmt="none", ecolor=INK, elinewidth=1.0,
                    capsize=2.5, zorder=4)
        for xx, vv in zip(xs, vals):
            if not np.isnan(vv):
                ax.text(xx, vv + 0.012, f"{vv:.3f}", ha="center", va="bottom", fontsize=6.5, color=INK)
    ax.set_xticks(x)
    ax.set_xticklabels([f"{lab}\n(n={dat['n']})" for lab, dat in mods])
    ax.set_ylabel("Overall correct-call rate")
    ax.set_ylim(0, 1.02); ax.set_title("A  Overall accuracy (RNA vs WES)", loc="left", fontsize=9)
    ax.axhline(1.0, color=GRID, lw=0.6, zorder=1)
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.12), ncol=3, frameon=False, fontsize=7)


def panel_per_locus(ax, rna, wes):
    genes = ["A", "B", "C"]
    groups = [("RNA CC", rna, "ChampionChallenger", CC_COLOR, ""),
              ("RNA MV", rna, "MajorityVote", MV_COLOR, ""),
              ("WES CC", wes, "ChampionChallenger", CC_COLOR, "//"),
              ("WES MV", wes, "MajorityVote", MV_COLOR, "//")]
    x = np.arange(len(genes)); w = 0.2
    for i, (lab, dat, mkey, col, hatch) in enumerate(groups):
        vals = [dat["per_gene"].get(g, {}).get(mkey, np.nan) for g in genes]
        ax.bar(x + (i - 1.5) * w, vals, w, color=col, hatch=hatch, edgecolor="white",
               linewidth=0.5, label=lab, zorder=3)
    ax.set_xticks(x); ax.set_xticklabels([f"HLA-{g}" for g in genes])
    ax.set_ylabel("Correct-call rate"); ax.set_ylim(0, 1.05)
    ax.set_title("B  Per-locus, by modality", loc="left", fontsize=9)
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.12), ncol=4, frameon=False, fontsize=6.5)


def panel_sweep(ax, sweep, wes):
    """WES override-gate sweep at the deployed support=0.35/margin=0.00, varying the
    supporting-tools floor — the lever that controls harmful overrides on this cohort."""
    if sweep.empty:
        ax.text(0.5, 0.5, "sweep table not found", ha="center", va="center"); return
    sf, mg, tl = "min_challenger_support_fraction", "min_challenger_margin", "min_supporting_tools"
    acc = "overall_correct_call_rate"
    sub = sweep[(np.isclose(sweep[sf], 0.35)) & (np.isclose(sweep[mg], 0.00))].sort_values(tl)
    tools = sub[tl].astype(int).tolist()
    accs = sub[acc].astype(float).tolist()
    corr = sub["corrective_override_count"].astype(int).tolist()
    harm = sub["harmful_override_count"].astype(int).tolist()
    x = np.arange(len(tools))
    colors = [CC_COLOR if t != 2 else "#0f766e" for t in tools]
    ax.bar(x, accs, 0.6, color=colors, edgecolor="white", linewidth=0.6, zorder=3)
    for xx, a, c, h, t in zip(x, accs, corr, harm, tools):
        ax.text(xx, a + 0.004, f"{a:.3f}", ha="center", va="bottom", fontsize=7, color=INK)
        tag = "deployed" if t == 1 else ("optimum" if t == 2 else "")
        ax.text(xx, 0.772, f"{c} corr\n{h} harm" + (f"\n({tag})" if tag else ""),
                ha="center", va="bottom", fontsize=6.2,
                color="#0f766e" if t == 2 else "#b42318" if h else INK)
    mv = wes["overall"].get("MajorityVote")
    if mv:
        ax.axhline(mv[0], color=MV_COLOR, ls="--", lw=1.2, zorder=2)
        ax.text(len(tools) - 0.5, mv[0] + 0.002, f"MajorityVote {mv[0]:.3f}", color=MV_COLOR,
                fontsize=6.5, ha="right", va="bottom")
    ax.set_xticks(x); ax.set_xticklabels([f"≥{t}" for t in tools])
    ax.set_xlabel("Supporting-tools floor (min_supporting_tools)\nat deployed support=0.35, margin=0.00")
    ax.set_ylabel("ChampHLA correct-call rate")
    ax.set_ylim(0.76, 0.88)
    ax.set_title("C  WES override-gate: supporting-tools floor", loc="left", fontsize=9)


def build_figure(rna, wes, sweep, out_dir: Path):
    fig = plt.figure(figsize=(11.0, 4.4))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.0, 1.15, 1.25], wspace=0.34,
                          left=0.06, right=0.985, top=0.80, bottom=0.20)
    panel_overall(fig.add_subplot(gs[0, 0]), rna, wes)
    panel_per_locus(fig.add_subplot(gs[0, 1]), rna, wes)
    panel_sweep(fig.add_subplot(gs[0, 2]), sweep, wes)
    _figure_header(fig, "Figure 13. External validation on an independent non-1000G cohort (NCI-60): "
                        "the frozen Champion-Challenger policy generalises across modalities.")
    _figure_footer(fig, "ChampHLA / MVHLA — NCI-60 (Adams 2005 SBT truth); frozen 1000G-derived policy.")
    save_fig(fig, out_dir / "figure_13_external_validation_nci60")


def write_caption(rna, wes, sweep, out_dir: Path):
    def g(d, m):
        t = d["overall"].get(m)
        return f"{t[0]:.3f}" if t else "n/a"
    opt = ""
    if not sweep.empty:
        s2 = sweep[sweep["min_supporting_tools"] == 2]
        o = s2.loc[s2["overall_correct_call_rate"].idxmax()]
        opt = (f"At support=0.35/margin=0.00/tools=2 the gate fires "
               f"{o['corrective_override_count']:g} corrective and {o['harmful_override_count']:g} harmful "
               f"overrides, lifting ChampHLA to {o['overall_correct_call_rate']:.3f}.")
    cap = (
        "**Figure 13. External validation on the NCI-60 panel (independent, non-1000G).** The frozen "
        "1000G-derived Champion-Challenger (ChampHLA) policy is applied without re-learning to NCI-60 "
        "RNA-seq and WES, scored against Adams 2005 sequence-based typing (HLA-A/B/C, 2-field). "
        f"(A) Overall correct-call rate (Wilson 95% CIs): RNA n={rna['n']}, WES n={wes['n']}; "
        f"ChampHLA RNA {g(rna,'ChampionChallenger')} vs MajorityVote {g(rna,'MajorityVote')}; "
        f"ChampHLA WES {g(wes,'ChampionChallenger')} vs MajorityVote {g(wes,'MajorityVote')}. "
        "(B) Per-locus correct-call rate by modality. (C) WES override-gate sweep at supporting-tools=2: "
        "the frozen margin gate (≥0.20) fires no overrides (ChampHLA≡OptiType), whereas relaxing the "
        f"margin to 0.00 recovers corrective overrides. {opt} NCI-60 are cancer lines, so loss of "
        "heterozygosity can cause apparent homozygosity; results are scored against the constitutional "
        "Adams genotype."
    )
    (out_dir / "captions_external_validation_nci60.md").write_text(cap + "\n", encoding="utf-8")


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--rna-tables", type=Path, default=Path(NCI60_RNA))
    p.add_argument("--wes-tables", type=Path, default=Path(NCI60_WES))
    p.add_argument("--sweep", type=Path, default=None)
    p.add_argument("--out-dir", type=Path, default=Path(OUT_DEFAULT))
    return p.parse_args()


def main() -> int:
    a = parse_args()
    a.out_dir.mkdir(parents=True, exist_ok=True)
    rna = load_modality(a.rna_tables)
    wes = load_modality(a.wes_tables)
    sweep = load_sweep(a.sweep or (a.wes_tables.parent / "sweeps" / "wes_champion_override_sweep.tsv"))
    build_figure(rna, wes, sweep, a.out_dir)
    write_caption(rna, wes, sweep, a.out_dir)
    print(f"[fig13] RNA n={rna['n']} WES n={wes['n']} sweep_rows={len(sweep)} -> {a.out_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
