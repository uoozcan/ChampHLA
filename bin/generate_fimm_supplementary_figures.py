#!/usr/bin/env python3
"""
generate_fimm_supplementary_figures.py

Generates publication-quality supplementary figures S4–S7 for the ChampHLA manuscript:

  S4 — FIMM cross-modality HLA concordance (scRNA vs BulkRNA vs WES)
  S5 — LOH/allele dropout QC heatmap (tools × modalities)
  S6 — FIMM WES LOH candidate classification (stacked bar per gene)
  S7 — Overall survival Kaplan-Meier curves (FIMM AML/MDS cohort)

Usage:
  module load gcc/11.3.0 biopythontools/11.3.0_3.10.6
  python3 bin/generate_fimm_supplementary_figures.py [--outdir PATH]

All figures use the hub figure standards (Arial, 300 DPI, Wong palette).
"""

import argparse
import math
import sys
import warnings
from collections import Counter
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "figure_hub"))
from plot_style import (
    apply_style, legend_outside, rotate_xticklabels, save_fig,
    tight_with_legend, WONG, WONG_CYCLE, TOOL_COLORS, MOD_COLORS,
)

apply_style()

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns

# ── Paths ─────────────────────────────────────────────────────────────────────
_LOCAL = Path("/scratch/project_2008084/pihla_local")
_PUBLISH = Path("/scratch/project_2008084/pihla-publish")

CALLS_DEIDENT = _LOCAL / "fimm_results" / "fimm_hla_calls_deident.tsv"
CALLS_IDENT   = _LOCAL / "fimm_results" / "fimm_hla_calls.tsv"

LOH_DEIDENT = _LOCAL / "analysis" / "loh_analysis" / "fimm_loh_candidates_wes_deident.tsv"
LOH_IDENT   = _LOCAL / "analysis" / "loh_analysis" / "fimm_loh_candidates_wes.tsv"

HOMOZYGOSITY = _PUBLISH / "analysis" / "loh_analysis" / "benchmark_homozygosity_error_summary.tsv"

PATIENT_DATA = Path("/scratch/project_2008084/patient_data_v3.xlsx")
BOOK2        = Path("/scratch/project_2008084/Book2.xlsx")
HLA_CALLS    = _LOCAL / "fimm_results" / "fimm_hla_calls.tsv"

DEFAULT_OUTDIR = _PUBLISH / "analysis" / "figures_final"

GENES = ["A", "B", "C"]
CUTOFF_DATE = datetime(2026, 5, 22)

# ── S4: Cross-Modality Concordance ────────────────────────────────────────────

MODALITIES = {
    "scrna":   ["scrna_arcashla", "scrna_optitype"],
    "bulkrna": ["bulkrna_OptiType", "bulkrna_arcasHLA", "bulkrna_SpecHLA"],
    "wes":     ["wes_OptiType", "wes_arcasHLA", "wes_SpecHLA"],
}


def _pairwise_concordance(df, col_a, col_b):
    if col_a == col_b:
        valid = df[[col_a]].dropna()
        return float("nan") if valid.empty else 1.0
    valid = df[[col_a, col_b]].dropna()
    if valid.empty:
        return float("nan")
    return (valid[col_a] == valid[col_b]).mean()


def figure_s4_concordance(df, outdir):
    """Cross-modality HLA concordance bar chart."""
    mod_pairs = [
        ("scrna", "wes", "scRNA vs WES"),
        ("scrna", "bulkrna", "scRNA vs Bulk RNA"),
        ("bulkrna", "wes", "Bulk RNA vs WES"),
    ]
    fig, axes = plt.subplots(1, len(GENES), figsize=(12, 4.5), sharey=True)
    cmap_list = [WONG["blue"], WONG["orange"], WONG["green"]]

    for gi, gene in enumerate(GENES):
        ax = axes[gi]
        ys, labels, colors_list = [], [], []

        for pi, (ma, mb, label) in enumerate(mod_pairs):
            col_a = next((c for c in df.columns if c.startswith(ma) and c.endswith(f"_{gene}")), None)
            col_b = next((c for c in df.columns if c.startswith(mb) and c.endswith(f"_{gene}")), None)
            if col_a and col_b:
                conc = _pairwise_concordance(df, col_a, col_b)
                ys.append(conc)
                labels.append(label)
                colors_list.append(cmap_list[pi])

        bars = ax.bar(range(len(ys)), ys, color=colors_list, edgecolor="white")
        for bar, y in zip(bars, ys):
            if not np.isnan(y):
                ax.text(bar.get_x() + bar.get_width() / 2, y + 0.02,
                        f"{y:.2f}", ha="center", va="bottom", fontsize=7)
        ax.set_xticks(range(len(labels)))
        ax.set_xticklabels(labels, rotation=30, ha="right", fontsize=7)
        ax.set_ylim(0, 1.15)
        ax.set_title(f"HLA-{gene}", fontweight="bold")
        if gi == 0:
            ax.set_ylabel("Concordance rate")

    fig.suptitle("Figure S4. FIMM Cross-Modality HLA Concordance",
                 fontsize=11, fontweight="bold", y=1.02)
    tight_with_legend(fig, right=0.98)
    save_fig(fig, outdir / "figure_s4_fimm_concordance")
    print("  Saved figure_s4_fimm_concordance.{pdf,svg,png}")


# ── S5: LOH / Allele Dropout Heatmap ─────────────────────────────────────────

TOOL_ORDER = ["ArcasHLA", "HLA-HD", "Kourami", "OptiType", "POLYSOLVER",
              "Seq2HLA", "SpecHLA", "T1K"]
MOD_LABELS = {"wgs": "WGS", "wes": "WES", "rnaseq": "RNA"}
MODALITY_ORDER = ["wgs", "wes", "rnaseq"]


def figure_s5_dropout(outdir):
    """Allele dropout heatmap and stacked bar."""
    if not HOMOZYGOSITY.exists():
        print(f"  SKIP figure_s5 — input not found: {HOMOZYGOSITY}")
        return

    df = pd.read_csv(HOMOZYGOSITY, sep="\t")
    df["false_duplicate_rate"] = pd.to_numeric(df["false_duplicate_rate"], errors="coerce")
    df["possible_allele_dropout"] = pd.to_numeric(df["possible_allele_dropout"], errors="coerce").fillna(0)

    col_keys = [(m, g) for m in MODALITY_ORDER for g in GENES]
    col_labels = [f"{MOD_LABELS.get(m, m)}\n{g}" for m, g in col_keys]
    tools_present = [t for t in TOOL_ORDER if t in df["tool"].values]

    heatmap = np.full((len(tools_present), len(col_keys)), np.nan)
    for i, tool in enumerate(tools_present):
        for j, (mod, gene) in enumerate(col_keys):
            row = df[(df["tool"] == tool) & (df["modality"] == mod) & (df["gene"] == gene)]
            if not row.empty:
                heatmap[i, j] = row["false_duplicate_rate"].iloc[0]

    bar_data = {}
    for tool in tools_present:
        counts = []
        for mod in MODALITY_ORDER:
            sub = df[(df["tool"] == tool) & (df["modality"] == mod)]
            counts.append(sub["possible_allele_dropout"].sum())
        bar_data[tool] = counts

    fig = plt.figure(figsize=(15, 5.5))
    gs = fig.add_gridspec(1, 2, width_ratios=[2.6, 1], wspace=0.38)
    ax_heat = fig.add_subplot(gs[0])
    ax_bar = fig.add_subplot(gs[1])

    cmap = plt.cm.RdYlGn_r
    cmap.set_bad(color="#e8e8e8")
    im = ax_heat.imshow(heatmap, cmap=cmap, vmin=0.0, vmax=1.0,
                        aspect="auto", interpolation="nearest")
    for x in np.arange(-0.5, len(col_keys), 1):
        ax_heat.axvline(x, color="white", lw=0.5)
    for y in np.arange(-0.5, len(tools_present), 1):
        ax_heat.axhline(y, color="white", lw=0.5)
    for sep in [2.5, 5.5]:
        ax_heat.axvline(sep, color="#555555", lw=1.2)
    for i in range(len(tools_present)):
        for j in range(len(col_keys)):
            val = heatmap[i, j]
            if np.isnan(val):
                continue
            txt = f"{val:.2f}" if val > 0 else "0"
            color = "white" if val > 0.6 else "black"
            ax_heat.text(j, i, txt, ha="center", va="center", fontsize=6.5, color=color)
    ax_heat.set_xticks(range(len(col_keys)))
    ax_heat.set_xticklabels(col_labels, fontsize=7.5)
    ax_heat.set_yticks(range(len(tools_present)))
    ax_heat.set_yticklabels(tools_present, fontsize=8.5)
    ax_heat.set_title("A  False Duplicate Rate", fontsize=9.5, loc="left", fontweight="bold")
    cbar = fig.colorbar(im, ax=ax_heat, fraction=0.025, pad=0.02)
    cbar.set_label("False duplicate rate", fontsize=8)
    cbar.ax.tick_params(labelsize=7)

    x = np.arange(len(MODALITY_ORDER))
    bottom = np.zeros(len(MODALITY_ORDER))
    bar_handles = []
    for tool in tools_present:
        vals = np.array(bar_data[tool])
        ax_bar.bar(x, vals, bottom=bottom,
                   color=TOOL_COLORS.get(tool, "#aaaaaa"),
                   label=tool, edgecolor="white", linewidth=0.4)
        bar_handles.append(mpatches.Patch(color=TOOL_COLORS.get(tool, "#aaaaaa"), label=tool))
        bottom += vals
    ax_bar.set_xticks(x)
    ax_bar.set_xticklabels([MOD_LABELS.get(m, m) for m in MODALITY_ORDER])
    ax_bar.set_ylabel("Possible allele dropout calls")
    ax_bar.set_title("B  Allele Dropout by Modality", fontsize=9.5, loc="left", fontweight="bold")
    legend_outside(ax_bar, loc="upper left", anchor=(1.02, 1.0),
                   handles=bar_handles, fontsize=7)

    fig.suptitle("Figure S5. HLA Allele Dropout Characterisation Across Tools and Modalities",
                 fontsize=10, y=1.03, fontweight="bold")
    rotate_xticklabels(ax_heat, angle=0, ha="center")
    tight_with_legend(fig, right=0.78, top=0.90, bottom=0.16, left=0.06)
    save_fig(fig, outdir / "figure_s5_allele_dropout_heatmap")
    print("  Saved figure_s5_allele_dropout_heatmap.{pdf,svg,png}")


# ── S6: FIMM WES LOH Status ──────────────────────────────────────────────────

STATUS_ORDER = ["candidate_loh", "allelic_imbalance_review",
                "balanced_heterozygous", "insufficient_evidence"]
STATUS_COLORS = {
    "candidate_loh":            WONG["vermil"],
    "allelic_imbalance_review": WONG["orange"],
    "balanced_heterozygous":    WONG["green"],
    "insufficient_evidence":    "#aaaaaa",
}
STATUS_LABELS = {
    "candidate_loh":            "Candidate LOH",
    "allelic_imbalance_review": "Allelic Imbalance (Review)",
    "balanced_heterozygous":    "Balanced Heterozygous",
    "insufficient_evidence":    "Insufficient Evidence",
}


def figure_s6_loh(outdir):
    """FIMM WES LOH status distribution per gene."""
    candidates_path = LOH_DEIDENT if LOH_DEIDENT.exists() else LOH_IDENT
    if not candidates_path.exists():
        print(f"  SKIP figure_s6 — input not found: {candidates_path}")
        return

    df = pd.read_csv(candidates_path, sep="\t", dtype=str)
    print(f"  LOH data: {len(df)} rows from {candidates_path.name}")

    counts = {g: {s: 0 for s in STATUS_ORDER} for g in GENES}
    for _, row in df.iterrows():
        gene = str(row.get("gene", ""))
        status = str(row.get("loh_status_wes", "insufficient_evidence"))
        if gene in counts and status in counts[gene]:
            counts[gene][status] += 1

    fig, ax = plt.subplots(figsize=(6, 4.5))
    x = np.arange(len(GENES))
    bottom = np.zeros(len(GENES))

    for status in STATUS_ORDER:
        vals = np.array([counts[g][status] for g in GENES], dtype=float)
        ax.bar(x, vals, bottom=bottom,
               color=STATUS_COLORS[status],
               label=STATUS_LABELS[status],
               edgecolor="white", linewidth=0.5)
        for xi, (v, b) in enumerate(zip(vals, bottom)):
            if v > 0:
                ax.text(xi, b + v / 2, str(int(v)),
                        ha="center", va="center", fontsize=8, color="white", fontweight="bold")
        bottom += vals

    ax.set_xticks(x)
    ax.set_xticklabels([f"HLA-{g}" for g in GENES])
    ax.set_ylabel("Number of samples")
    ax.set_title("Figure S6. FIMM WES LOH Status Distribution per Gene",
                 loc="left", fontweight="bold", fontsize=9.5)
    legend_outside(ax, loc="upper left", anchor=(1.02, 1.0))

    ax.text(0.01, -0.18,
            "WES allele frequencies less reliable for LOH than WGS.\n"
            "WGS data required for confirmation of candidate LOH events.",
            transform=ax.transAxes, fontsize=6.5, color=WONG["vermil"], style="italic")

    tight_with_legend(fig, right=0.68)
    save_fig(fig, outdir / "figure_s6_fimm_loh_wes")
    print("  Saved figure_s6_fimm_loh_wes.{pdf,svg,png}")


# ── S7: Survival KM Curves ───────────────────────────────────────────────────

_AA_PROPS = {
    "A": (0,   8.1,  31),  "R": (0.65, 10.5, 124), "N": (0.23, 11.6,  56),
    "D": (0.23, 13.0,  54), "C": (0.77,  5.5,  55), "Q": (0.50, 10.5,  85),
    "E": (0.42, 12.3,  83), "G": (0,     9.0,   3), "H": (0.58, 10.4,  96),
    "I": (0,    5.2,  111), "L": (0,     4.9, 111), "K": (0.33, 11.3, 119),
    "M": (0,    5.7,  105), "F": (0,     5.2, 132), "P": (0.39,  8.0,  32.5),
    "S": (0.70,  9.2,  32), "T": (0.71,  8.6,  61), "W": (0.13,  5.4, 170),
    "Y": (0.20,  6.2, 136), "V": (0,     5.9,  84),
}

HLA_A_ABG = {
    "A*01:01": "YTLFYERKTQ", "A*02:01": "YTLFYDRKTQ", "A*02:05": "YTLFYDRKTQ",
    "A*03:01": "YTLFYERKTH", "A*11:01": "YTLFYERKTH", "A*24:02": "YTHFYDRKTQ",
    "A*25:01": "YTLFYDRKTH", "A*29:01": "YTLFYERKAQ", "A*30:01": "YTLFYERKAQ",
    "A*31:01": "YTLFYDRKTH", "A*32:01": "YTLFYDRKTQ", "A*68:01": "YTLFYDRKTQ",
}
HLA_B_ABG = {
    "B*07:02": "YNLFYNRTQH", "B*08:01": "YTHFYDRTQR", "B*13:02": "YTLFYDRTQH",
    "B*15:01": "YTLFYNRTQH", "B*18:01": "YTLFYNRTQH", "B*18:03": "YTLFYNRTQH",
    "B*27:05": "YTHFYNRTQH", "B*35:01": "YTLFYNRTQH", "B*35:08": "YTLFYNRTQH",
    "B*37:01": "YTLFYDRTQH", "B*38:01": "YTLFYNRTQH", "B*39:01": "YTLFYNRTQH",
    "B*40:01": "YTLFYDRTQH", "B*41:01": "YTLFYNRTQH", "B*44:02": "YTLFYDRTQH",
    "B*44:08": "YTLFYDRTQH", "B*51:01": "YTLFYNRTQH", "B*56:01": "YTLFYNRTQH",
    "B*57:01": "YTLFYDRTQH", "B*58:01": "YTLFYDRTQH",
}
HLA_C_ABG = {
    "C*01:02": "YTLFYRRTH", "C*01:06": "YTLFYRRTH", "C*01:08": "YTLFYRRTH",
    "C*02:02": "YTLFYRRTQ", "C*03:02": "YTLFYRRTH", "C*03:03": "YTLFYRRTH",
    "C*03:04": "YTLFYRRTH", "C*04:01": "YTLFYRRTQ", "C*05:01": "YTLFYRRTQ",
    "C*06:02": "YTLFYRRTQ", "C*07:01": "YTLFYRRTQ", "C*07:02": "YTLFYRRTQ",
    "C*12:03": "YTLFYRRTH", "C*14:02": "YTLFYRRTH", "C*17:01": "YTLFYRRTQ",
}
_LOCUS_TABLES = {"A": HLA_A_ABG, "B": HLA_B_ABG, "C": HLA_C_ABG}


def grantham(aa1, aa2):
    if aa1 == aa2:
        return 0.0
    if aa1 not in _AA_PROPS or aa2 not in _AA_PROPS:
        return 50.0
    c1, p1, v1 = _AA_PROPS[aa1]
    c2, p2, v2 = _AA_PROPS[aa2]
    return math.sqrt(0.1018 * (c1-c2)**2 + 0.000399 * (p1-p2)**2 + 0.000791 * (v1-v2)**2)


def hed_for_pair(allele_str, locus):
    if pd.isna(allele_str):
        return float("nan")
    parts = str(allele_str).split("|")
    if len(parts) != 2:
        return float("nan")
    a1, a2 = parts[0].strip(), parts[1].strip()
    if a1 == a2:
        return 0.0
    table = _LOCUS_TABLES.get(locus, {})
    seq1, seq2 = table.get(a1), table.get(a2)
    if seq1 is None:
        grp = a1.split(":")[0]
        candidates = [v for k, v in table.items() if k.startswith(grp)]
        seq1 = candidates[0] if candidates else None
    if seq2 is None:
        grp = a2.split(":")[0]
        candidates = [v for k, v in table.items() if k.startswith(grp)]
        seq2 = candidates[0] if candidates else None
    if seq1 is None or seq2 is None:
        return float("nan")
    n = min(len(seq1), len(seq2))
    dists = [grantham(seq1[i], seq2[i]) for i in range(n)]
    return sum(dists) / n if dists else 0.0


def _load_survival_data():
    """Load and merge FIMM survival data (mirrors generate_fimm_survival_hed_report.py)."""
    from lifelines import KaplanMeierFitter
    from lifelines.statistics import logrank_test

    hla = pd.read_csv(HLA_CALLS, sep="\t")
    hla["patient_id"] = hla["sample_id"].apply(lambda s: ".".join(s.split("_")[:2]))

    b2 = pd.read_excel(BOOK2, header=0)
    b2["donor_str"] = b2["donor"].astype(str)
    fimm_pids = set(hla["patient_id"].unique())
    b2_fimm = b2[b2["donor_str"].isin(fimm_pids)].copy()
    b2_fimm["date"] = pd.to_datetime(b2_fimm["date"], errors="coerce")
    b2_earliest = (
        b2_fimm.sort_values("date")
        .groupby("donor_str", as_index=False)
        .first()
        [["donor_str", "diagnosis (text)", "disease stage", "date", "months since diagnosis"]]
    )
    b2_earliest.rename(columns={
        "donor_str": "patient_id", "diagnosis (text)": "diagnosis",
        "disease stage": "stage", "date": "sampling_date",
        "months since diagnosis": "months_since_dx",
    }, inplace=True)

    pd3 = pd.read_excel(PATIENT_DATA, sheet_name="scrna")
    pd3["study_id"] = pd3["study id"].astype(str)
    pd3_surv = pd3[pd3["study_id"].isin(fimm_pids)][
        ["study_id", "date of death", "cause of death", "gender"]
    ].copy()
    pd3_surv.rename(columns={
        "study_id": "patient_id", "date of death": "date_of_death",
        "cause of death": "cause_of_death",
    }, inplace=True)
    pd3_surv["date_of_death"] = pd.to_datetime(pd3_surv["date_of_death"], errors="coerce")

    hla_pat = hla.sort_values("sample_id").drop_duplicates("patient_id", keep="first").copy()

    def consensus_allele_pair(row, locus):
        pairs = []
        for col in [f"wes_OptiType_{locus}", f"wes_arcasHLA_{locus}", f"wes_SpecHLA_{locus}"]:
            val = row.get(col, "")
            if pd.notna(val) and str(val).strip() and "|" in str(val):
                a1, a2 = str(val).split("|", 1)
                pairs.append(tuple(sorted([a1.strip(), a2.strip()])))
        if not pairs:
            return None, False
        most_common, count = Counter(pairs).most_common(1)[0]
        concordant = (count == len(pairs) and len(pairs) == 3)
        return f"{most_common[0]}|{most_common[1]}", concordant

    for locus in GENES:
        consensus_vals, concordant_vals = zip(
            *[consensus_allele_pair(row, locus) for _, row in hla_pat.iterrows()]
        )
        hla_pat[f"consensus_{locus}"] = list(consensus_vals)
        hla_pat[f"HED_{locus}"] = hla_pat[f"consensus_{locus}"].apply(
            lambda x: hed_for_pair(x, locus)
        )

    hla_pat["HED_total"] = hla_pat[["HED_A", "HED_B", "HED_C"]].sum(axis=1, min_count=1)

    df = pd3_surv.merge(b2_earliest, on="patient_id", how="outer")
    df = df.merge(
        hla_pat[["patient_id", "consensus_A", "consensus_B", "consensus_C",
                 "HED_A", "HED_B", "HED_C", "HED_total"]],
        on="patient_id", how="outer"
    )

    def _cohort(pid):
        prefix = str(pid).split(".")[0]
        return "FIMM" if prefix in ("FH", "FHRB") else prefix

    df["cohort"] = df["patient_id"].apply(_cohort)
    df["event"] = df["date_of_death"].notna().astype(int)
    df["os_from_dx"] = np.where(
        df["event"] == 1,
        df["months_since_dx"] + (df["date_of_death"] - df["sampling_date"]).dt.days / 30.44,
        df["months_since_dx"] + (CUTOFF_DATE - df["sampling_date"]).dt.days / 30.44,
    )

    df_surv = df[df["os_from_dx"].notna() & (df["os_from_dx"] > 0)].copy()
    return df, df_surv


def figure_s7_survival(outdir):
    """Two-panel KM figure: by cohort and by HED group."""
    if not PATIENT_DATA.exists() or not BOOK2.exists() or not HLA_CALLS.exists():
        print("  SKIP figure_s7 — survival data files not found")
        return

    from lifelines import KaplanMeierFitter
    from lifelines.statistics import logrank_test

    print("  Loading survival data...")
    df, df_surv = _load_survival_data()

    HED_COLORS = {"High HED": WONG["vermil"], "Low HED": WONG["blue"]}

    def simplify_dx(d):
        d = str(d)
        if "MDS" in d and "AML" in d:
            return "MDS->AML"
        if "AML" in d:
            return "AML"
        if "MDS" in d:
            return "MDS"
        return "Other"

    fig, axes = plt.subplots(1, 2, figsize=(14, 5.5))

    # Panel A: KM by diagnosis group
    ax = axes[0]
    df2 = df_surv.copy()
    df2["dx_grp"] = df2["diagnosis"].apply(simplify_dx)
    DX_COLORS = {"AML": WONG["blue"], "MDS->AML": WONG["orange"],
                 "MDS": WONG["green"], "Other": "#aaaaaa"}
    kmf = KaplanMeierFitter()
    for grp in ["AML", "MDS->AML", "MDS"]:
        mask = df2["dx_grp"] == grp
        if mask.sum() < 2:
            continue
        kmf.fit(df2.loc[mask, "os_from_dx"], df2.loc[mask, "event"],
                label=f"{grp} (n={mask.sum()})")
        kmf.plot_survival_function(ax=ax, ci_show=True,
                                   color=DX_COLORS.get(grp, WONG["blue"]))
    ax.set_xlabel("OS from diagnosis (months)")
    ax.set_ylabel("Survival probability")
    ax.set_title("A  Overall Survival by Diagnosis", fontweight="bold", fontsize=9.5)
    ax.set_ylim(-0.05, 1.05)
    ax.legend(frameon=False, fontsize=8)

    # Panel B: KM by HED group
    ax = axes[1]
    df3 = df_surv.dropna(subset=["HED_total", "os_from_dx"]).copy()
    if len(df3) >= 4:
        median_hed = df3["HED_total"].median()
        df3["HED_group"] = np.where(df3["HED_total"] >= median_hed, "High HED", "Low HED")
        kmf2 = KaplanMeierFitter()
        for grp in ["High HED", "Low HED"]:
            mask = df3["HED_group"] == grp
            if mask.sum() < 2:
                continue
            kmf2.fit(df3.loc[mask, "os_from_dx"], df3.loc[mask, "event"],
                     label=f"{grp} (n={mask.sum()})")
            kmf2.plot_survival_function(ax=ax, ci_show=True,
                                        color=HED_COLORS.get(grp, WONG["blue"]))

        hi = df3[df3["HED_group"] == "High HED"]
        lo = df3[df3["HED_group"] == "Low HED"]
        if len(hi) >= 2 and len(lo) >= 2:
            result = logrank_test(lo["os_from_dx"], hi["os_from_dx"],
                                  lo["event"], hi["event"])
            ax.text(0.05, 0.15, f"Log-rank p = {result.p_value:.3f}",
                    transform=ax.transAxes, fontsize=9, color=WONG["blue"])

        ax.text(0.05, 0.07,
                f"Split at median HED = {median_hed:.1f} — exploratory (n={len(df3)})",
                transform=ax.transAxes, fontsize=7, color="#888888", style="italic")

    ax.set_xlabel("OS from diagnosis (months)")
    ax.set_ylabel("Survival probability")
    ax.set_title("B  OS by HED Group (median split)", fontweight="bold", fontsize=9.5)
    ax.set_ylim(-0.05, 1.05)
    ax.legend(frameon=False, fontsize=8)

    fig.suptitle("Figure S7. FIMM Cohort Overall Survival and HLA Evolutionary Divergence",
                 fontsize=11, fontweight="bold", y=1.02)
    fig.text(0.01, -0.04,
             "Exploratory analysis — small cohort (n~28). Interpret with caution. "
             "HED computed from WES consensus alleles at HLA-A, -B, -C using Grantham distances.",
             fontsize=7, color="#888888")
    tight_with_legend(fig, right=0.95, top=0.90, bottom=0.12)
    save_fig(fig, outdir / "figure_s7_fimm_survival_hed")
    print("  Saved figure_s7_fimm_survival_hed.{pdf,svg,png}")


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    args = ap.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    # S4: Concordance
    calls_path = CALLS_DEIDENT if CALLS_DEIDENT.exists() else CALLS_IDENT
    if calls_path.exists():
        df_calls = pd.read_csv(calls_path, sep="\t", dtype=str)
        print(f"Loaded {len(df_calls)} rows from {calls_path.name}")
        figure_s4_concordance(df_calls, args.outdir)
    else:
        print(f"SKIP S4 — input not found: {calls_path}")

    # S5: Allele dropout
    figure_s5_dropout(args.outdir)

    # S6: LOH
    figure_s6_loh(args.outdir)

    # S7: Survival
    figure_s7_survival(args.outdir)

    print(f"\nAll supplementary figures saved to {args.outdir}")


if __name__ == "__main__":
    main()
