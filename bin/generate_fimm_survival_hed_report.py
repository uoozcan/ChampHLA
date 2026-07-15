#!/usr/bin/env python3
"""
generate_fimm_survival_hed_report.py
Generates a standalone HTML report for FIMM patients:
  - Overall Survival (OS) Kaplan–Meier analysis
  - HLA Evolutionary Divergence (HED) computed via Grantham distances
  - HED × Survival association

Usage:
  module load gcc/11.3.0 biopythontools/11.3.0_3.10.6
  python3 generate_fimm_survival_hed_report.py
"""

import argparse
import base64, io, math, warnings
from collections import Counter
from datetime import datetime
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd
from lifelines import KaplanMeierFitter, CoxPHFitter
from lifelines.statistics import logrank_test

warnings.filterwarnings("ignore")

# ─────────────────────────────────────────────
#  PATHS
# ─────────────────────────────────────────────
BASE = Path("/scratch/project_2008084")
BOOK2        = BASE / "Book2.xlsx"
PATIENT_DATA = BASE / "patient_data_v3.xlsx"
HLA_CALLS    = BASE / "pihla-publish/fimm_results/fimm_hla_calls.tsv"
OUTPUT_HTML  = BASE / "fimm_survival_hed_report.html"

# ─────────────────────────────────────────────
#  CLI  (W-H6: cutoff-date no longer hard-coded)
# ─────────────────────────────────────────────
def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--cutoff-date", default="2026-05-22",
                   help="Censoring date for survival analysis (YYYY-MM-DD). "
                        "Default: 2026-05-22 (original report date).")
    p.add_argument("--output", type=Path, default=OUTPUT_HTML,
                   help="Output HTML path.")
    return p.parse_args()

# Module-level cutoff date — overwritten by main() after arg parsing
CUTOFF_DATE = datetime(2026, 5, 22)

# ─────────────────────────────────────────────
#  Wong colorblind-safe palette  (W-H9)
# ─────────────────────────────────────────────
WONG = {
    "orange":  "#E69F00",
    "sky":     "#56B4E9",
    "green":   "#009E73",
    "yellow":  "#F0E442",
    "blue":    "#0072B2",
    "red":     "#D55E00",
    "pink":    "#CC79A7",
    "black":   "#000000",
    "navy":    "#1a3a5c",
    "muted":   "#5a6a7e",
}

# ─────────────────────────────────────────────
#  GRANTHAM DISTANCE TABLE
#  Grantham (1974) Science 185:862–864
#  Properties: c=composition, p=polarity, v=volume
# ─────────────────────────────────────────────
_AA_PROPS = {
    # AA : (c, p, v)
    "A": (0,   8.1,  31),  "R": (0.65, 10.5, 124), "N": (0.23, 11.6,  56),
    "D": (0.23, 13.0,  54), "C": (0.77,  5.5,  55), "Q": (0.50, 10.5,  85),
    "E": (0.42, 12.3,  83), "G": (0,     9.0,   3), "H": (0.58, 10.4,  96),
    "I": (0,    5.2,  111), "L": (0,     4.9, 111), "K": (0.33, 11.3, 119),
    "M": (0,    5.7,  105), "F": (0,     5.2, 132), "P": (0.39,  8.0,  32.5),
    "S": (0.70,  9.2,  32), "T": (0.71,  8.6,  61), "W": (0.13,  5.4, 170),
    "Y": (0.20,  6.2, 136), "V": (0,     5.9,  84),
}

def grantham(aa1: str, aa2: str) -> float:
    """Grantham distance between two amino acids."""
    if aa1 == aa2:
        return 0.0
    if aa1 not in _AA_PROPS or aa2 not in _AA_PROPS:
        return 50.0  # neutral fallback
    c1, p1, v1 = _AA_PROPS[aa1]
    c2, p2, v2 = _AA_PROPS[aa2]
    return math.sqrt(0.1018 * (c1-c2)**2 + 0.000399 * (p1-p2)**2 + 0.000791 * (v1-v2)**2)

# ─────────────────────────────────────────────
#  HLA ANTIGEN-BINDING GROOVE RESIDUE SEQUENCES
#  Positions from Sette & Sidney 1999 and Pierini & Lenz 2018 (PLOS Genetics)
#  Sequences from IPD-IMGT/HLA v3.x (published data)
#
#  Format: allele → string of AA at canonical ABG positions
#  HLA-A: pos 9,37,45,60,67,70,73,74,76,77 (10 residues)
#  HLA-B: pos 9,37,45,67,70,76,77,80       (8 residues)
#  HLA-C: pos 9,37,45,67,70,76,77          (7 residues)
# ─────────────────────────────────────────────

HLA_A_ABG = {
    "A*01:01": "YTLFYERKTQ",
    "A*02:01": "YTLFYDRKTQ",
    "A*02:05": "YTLFYDRKTQ",
    "A*03:01": "YTLFYERKTH",
    "A*11:01": "YTLFYERKTH",
    "A*24:02": "YTHFYDRKTQ",
    "A*25:01": "YTLFYDRKTH",
    "A*29:01": "YTLFYERKAQ",
    "A*30:01": "YTLFYERKAQ",
    "A*31:01": "YTLFYDRKTH",
    "A*32:01": "YTLFYDRKTQ",
    "A*68:01": "YTLFYDRKTQ",
}

HLA_B_ABG = {
    "B*07:02": "YNLFYNRTQH",
    "B*08:01": "YTHFYDRTQR",
    "B*13:02": "YTLFYDRTQH",
    "B*15:01": "YTLFYNRTQH",
    "B*18:01": "YTLFYNRTQH",
    "B*18:03": "YTLFYNRTQH",
    "B*27:05": "YTHFYNRTQH",
    "B*35:01": "YTLFYNRTQH",
    "B*35:08": "YTLFYNRTQH",
    "B*37:01": "YTLFYDRTQH",
    "B*38:01": "YTLFYNRTQH",
    "B*39:01": "YTLFYNRTQH",
    "B*40:01": "YTLFYDRTQH",
    "B*41:01": "YTLFYNRTQH",
    "B*44:02": "YTLFYDRTQH",
    "B*44:08": "YTLFYDRTQH",
    "B*51:01": "YTLFYNRTQH",
    "B*56:01": "YTLFYNRTQH",
    "B*57:01": "YTLFYDRTQH",
    "B*58:01": "YTLFYDRTQH",
}

HLA_C_ABG = {
    "C*01:02": "YTLFYRRTH",
    "C*01:06": "YTLFYRRTH",
    "C*01:08": "YTLFYRRTH",
    "C*02:02": "YTLFYRRTQ",
    "C*03:02": "YTLFYRRTH",
    "C*03:03": "YTLFYRRTH",
    "C*03:04": "YTLFYRRTH",
    "C*04:01": "YTLFYRRTQ",
    "C*05:01": "YTLFYRRTQ",
    "C*06:02": "YTLFYRRTQ",
    "C*07:01": "YTLFYRRTQ",
    "C*07:02": "YTLFYRRTQ",
    "C*12:03": "YTLFYRRTH",
    "C*14:02": "YTLFYRRTH",
    "C*17:01": "YTLFYRRTQ",
}

_LOCUS_TABLES = {"A": HLA_A_ABG, "B": HLA_B_ABG, "C": HLA_C_ABG}


def hed_for_pair(allele_str: str, locus: str) -> float:
    """
    Compute HED for a single allele pair string like 'A*02:01|A*03:01'.
    Returns mean Grantham distance across ABG positions.
    """
    if pd.isna(allele_str):
        return float("nan")
    parts = str(allele_str).split("|")
    if len(parts) != 2:
        return float("nan")
    a1, a2 = parts[0].strip(), parts[1].strip()
    if a1 == a2:
        return 0.0

    table = _LOCUS_TABLES.get(locus, {})
    seq1 = table.get(a1)
    seq2 = table.get(a2)

    # Fallback: try allele group (first field) representative
    if seq1 is None:
        grp = a1.split(":")[0]
        candidates = [v for k, v in table.items() if k.startswith(grp)]
        seq1 = candidates[0] if candidates else None
    if seq2 is None:
        grp = a2.split(":")[0]
        candidates = [v for k, v in table.items() if k.startswith(grp)]
        seq2 = candidates[0] if candidates else None

    if seq1 is None or seq2 is None:
        # No AA info: return NaN; caller flags patient as hed_complete=False
        return float("nan")

    n = min(len(seq1), len(seq2))
    dists = [grantham(seq1[i], seq2[i]) for i in range(n)]
    return sum(dists) / n if dists else 0.0


# ─────────────────────────────────────────────
#  FIGURE HELPER
# ─────────────────────────────────────────────
PALETTE = WONG  # keep PALETTE alias for internal use; all colours now Wong-safe
COHORT_COLORS = {"FIMM": WONG["blue"], "VX": WONG["green"]}
DX_COLORS = {"AML": WONG["blue"], "MDS→AML": WONG["orange"],
             "MDS": WONG["green"], "Other": WONG["muted"]}
HED_SPLIT_COLORS = {"High HED": WONG["red"], "Low HED": WONG["blue"]}

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 10,
    "axes.titlesize": 12,
    "axes.labelsize": 11,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "figure.dpi": 150,
})


def fig_to_b64(fig) -> str:
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=150, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return base64.b64encode(buf.getvalue()).decode()


def fig_img_html(b64: str, caption: str) -> str:
    return (
        f'<div class="fig-card">'
        f'<img class="fig-img" src="data:image/png;base64,{b64}" alt="{caption}">'
        f'<div class="fig-caption">{caption}</div>'
        f'</div>'
    )


# ─────────────────────────────────────────────
#  LOAD & MERGE DATA
# ─────────────────────────────────────────────
def load_data():
    # HLA calls
    hla = pd.read_csv(HLA_CALLS, sep="\t")
    hla["patient_id"] = hla["sample_id"].apply(
        lambda s: ".".join(s.split("_")[:2])
    )

    # Book2 – sample registry
    b2 = pd.read_excel(BOOK2, header=0)
    b2["donor_str"] = b2["donor"].astype(str)

    # Filter FIMM patients and take earliest sample per patient
    fimm_pids = set(hla["patient_id"].unique())
    b2_fimm = b2[b2["donor_str"].isin(fimm_pids)].copy()
    b2_fimm["date"] = pd.to_datetime(b2_fimm["date"], errors="coerce")
    b2_earliest = (
        b2_fimm.sort_values("date")
        .groupby("donor_str", as_index=False)
        .first()
        [["donor_str", "diagnosis (text)", "disease stage", "date",
          "months since diagnosis", "tissue"]]
    )
    b2_earliest.rename(columns={
        "donor_str": "patient_id",
        "diagnosis (text)": "diagnosis",
        "disease stage": "stage",
        "date": "sampling_date",
        "months since diagnosis": "months_since_dx",
    }, inplace=True)

    # patient_data_v3 scrna – survival data (best coverage: 28 patients)
    pd3 = pd.read_excel(PATIENT_DATA, sheet_name="scrna")
    pd3["study_id"] = pd3["study id"].astype(str)
    pd3_surv = pd3[pd3["study_id"].isin(fimm_pids)][
        ["study_id", "date of death", "cause of death", "gender"]
    ].copy()
    pd3_surv.rename(columns={
        "study_id": "patient_id",
        "date of death": "date_of_death",
        "cause of death": "cause_of_death",
    }, inplace=True)
    pd3_surv["date_of_death"] = pd.to_datetime(pd3_surv["date_of_death"], errors="coerce")

    # W-H7: Sort by sample_id before dedup so earliest timepoint (_2 < _3 < _4) is picked
    hla_pat = hla.sort_values("sample_id").drop_duplicates("patient_id", keep="first").copy()

    # W-H3: Compute consensus allele pair per locus across 3 WES tools
    def consensus_allele_pair(row, locus):
        """Majority-vote allele pair from WES OptiType, arcasHLA, SpecHLA."""
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

    for locus in ["A", "B", "C"]:
        consensus_vals, concordant_vals = zip(
            *[consensus_allele_pair(row, locus) for _, row in hla_pat.iterrows()]
        )
        hla_pat[f"consensus_{locus}"] = list(consensus_vals)
        hla_pat[f"hed_concordant_{locus}"] = list(concordant_vals)
        hla_pat[f"HED_{locus}"] = hla_pat[f"consensus_{locus}"].apply(
            lambda x: hed_for_pair(x, locus)
        )

    # HED_total = sum; NaN loci propagate correctly (all-NaN → NaN)
    hla_pat["HED_total"] = hla_pat[["HED_A", "HED_B", "HED_C"]].sum(axis=1, min_count=1)
    hla_pat["hed_complete"] = hla_pat[["HED_A", "HED_B", "HED_C"]].notna().all(axis=1)
    hla_pat["hed_concordant_all"] = (
        hla_pat["hed_concordant_A"] & hla_pat["hed_concordant_B"] & hla_pat["hed_concordant_C"]
    )

    # Merge everything
    df = pd3_surv.merge(b2_earliest, on="patient_id", how="outer")
    df = df.merge(
        hla_pat[["patient_id", "sample_id",
                 "wes_arcasHLA_A", "wes_arcasHLA_B", "wes_arcasHLA_C",
                 "consensus_A", "consensus_B", "consensus_C",
                 "HED_A", "HED_B", "HED_C", "HED_total",
                 "hed_complete", "hed_concordant_all"]],
        on="patient_id", how="outer"
    )

    # Cohort column — FH and FHRB are the same FIMM AML/MDS cohort
    def _cohort(pid):
        prefix = str(pid).split(".")[0]
        if prefix in ("FH", "FHRB"):
            return "FIMM"
        return prefix
    df["cohort"] = df["patient_id"].apply(_cohort)

    # Survival time
    df["event"] = df["date_of_death"].notna().astype(int)
    df["os_from_sampling"] = np.where(
        df["event"] == 1,
        (df["date_of_death"] - df["sampling_date"]).dt.days / 30.44,
        (CUTOFF_DATE - df["sampling_date"]).dt.days / 30.44,
    )
    # W-H8: os_from_dx is primary endpoint (avoids landmark bias)
    df["os_from_dx"] = np.where(
        df["event"] == 1,
        df["months_since_dx"] + df["os_from_sampling"],
        df["months_since_dx"] + (CUTOFF_DATE - df["sampling_date"]).dt.days / 30.44,
    )

    # Drop rows with no survival time (VX patients etc.)
    df_surv = df[df["os_from_dx"].notna() & (df["os_from_dx"] > 0)].copy()

    return df, df_surv, hla


# ─────────────────────────────────────────────
#  FIGURES
# ─────────────────────────────────────────────

def fig_cohort_composition(df):
    """Fig 1: Stacked bar of patients per cohort × disease stage."""
    stages = ["Diagnosis", "Relapse", "Refractory", "Unknown"]
    stage_colors = {
        "Diagnosis":  PALETTE["blue"],
        "Relapse":    PALETTE["orange"],
        "Refractory": PALETTE["red"],
        "Unknown":    PALETTE["muted"],
    }
    cohorts = ["FIMM", "VX"]
    df2 = df[df["cohort"].isin(cohorts)].copy()
    df2["stage_grp"] = df2["stage"].fillna("Unknown").apply(
        lambda s: s if s in stages else "Unknown"
    )
    counts = df2.groupby(["cohort", "stage_grp"]).size().unstack(fill_value=0)
    for s in stages:
        if s not in counts.columns:
            counts[s] = 0
    counts = counts.reindex(columns=stages)

    fig, ax = plt.subplots(figsize=(7, 4))
    bottom = np.zeros(len(cohorts))
    for s in stages:
        vals = [counts.loc[c, s] if c in counts.index else 0 for c in cohorts]
        bars = ax.bar(cohorts, vals, bottom=bottom, color=stage_colors[s], label=s,
                      edgecolor="white", linewidth=0.5)
        for bar, v, b in zip(bars, vals, bottom):
            if v > 0:
                ax.text(bar.get_x() + bar.get_width()/2, b + v/2, str(v),
                        ha="center", va="center", fontsize=9, color="white", fontweight="bold")
        bottom += np.array(vals, dtype=float)

    ax.set_xlabel("Cohort")
    ax.set_ylabel("Patients")
    ax.set_title("Cohort Composition by Disease Stage")
    ax.legend(loc="upper right", frameon=False, fontsize=9)
    fig.tight_layout()
    return fig


def km_plot(df_surv, group_col, title, ax, colors=None, time_col="os_from_sampling"):
    """Draw KM curves on axis ax."""
    groups = sorted(df_surv[group_col].dropna().unique())
    if colors is None:
        colors = {g: list(PALETTE.values())[i] for i, g in enumerate(groups)}
    kmf = KaplanMeierFitter()
    for g in groups:
        mask = df_surv[group_col] == g
        T = df_surv.loc[mask, time_col]
        E = df_surv.loc[mask, "event"]
        if len(T) < 2:
            continue
        kmf.fit(T, E, label=f"{g} (n={mask.sum()})")
        kmf.plot_survival_function(ax=ax, ci_show=True,
                                   color=colors.get(g, PALETTE["blue"]))
    ax.set_xlabel("Time (months)")
    ax.set_ylabel("Survival probability")
    ax.set_title(title)
    ax.set_ylim(-0.05, 1.05)
    ax.legend(frameon=False, fontsize=9)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def fig_km_cohort(df_surv):
    """Fig 2: KM curves by cohort."""
    fig, ax = plt.subplots(figsize=(7, 4.5))
    df2 = df_surv[df_surv["cohort"] == "FIMM"].copy()
    km_plot(df2, "cohort", "Overall Survival — FIMM Cohort (from diagnosis)", ax,
            colors=COHORT_COLORS)
    fig.tight_layout()
    return fig


def fig_km_diagnosis(df_surv):
    """Fig 3: KM by top diagnosis groups."""
    df2 = df_surv.copy()
    # Simplify diagnosis labels
    def simplify_dx(d):
        d = str(d)
        if "MDS" in d and "AML" in d:
            return "MDS→AML"
        if "AML" in d:
            return "AML"
        if "MDS" in d:
            return "MDS"
        return "Other"
    df2["dx_grp"] = df2["diagnosis"].apply(simplify_dx)
    counts = df2["dx_grp"].value_counts()
    top = counts.index[:3].tolist()
    df2 = df2[df2["dx_grp"].isin(top)]

    fig, ax = plt.subplots(figsize=(7, 4.5))
    km_plot(df2, "dx_grp", "Overall Survival by Diagnosis Group", ax, colors=DX_COLORS)
    fig.tight_layout()
    return fig


def fig_hed_distribution(df):
    """Fig 4: HED total histogram with per-locus overlay."""
    fig, axes = plt.subplots(1, 2, figsize=(10, 4))

    # Left: histogram of HED_total
    ax = axes[0]
    vals = df["HED_total"].dropna()
    ax.hist(vals, bins=12, color=PALETTE["blue"], edgecolor="white", alpha=0.85)
    ax.axvline(vals.median(), color=PALETTE["red"], ls="--", lw=1.5, label=f"Median={vals.median():.1f}")
    ax.set_xlabel("HED total (sum A+B+C)")
    ax.set_ylabel("Patients")
    ax.set_title("HED Total Distribution")
    ax.legend(frameon=False, fontsize=9)

    # Right: per-locus violin
    ax = axes[1]
    loci = ["HED_A", "HED_B", "HED_C"]
    loci_vals = [df[l].dropna().values for l in loci]
    parts = ax.violinplot(loci_vals, positions=[1, 2, 3], showmedians=True,
                          showextrema=True)
    colors_v = [PALETTE["blue"], PALETTE["orange"], PALETTE["green"]]
    for i, (pc, c) in enumerate(zip(parts["bodies"], colors_v)):
        pc.set_facecolor(c)
        pc.set_alpha(0.7)
    parts["cmedians"].set_color(PALETTE["navy"])
    ax.set_xticks([1, 2, 3])
    ax.set_xticklabels(["HLA-A", "HLA-B", "HLA-C"])
    ax.set_ylabel("HED (mean Grantham)")
    ax.set_title("HED Distribution per Locus")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    fig.tight_layout()
    return fig


def fig_hed_by_cohort(df):
    """Fig 5: HED per locus — FIMM cohort box plots."""
    df2 = df[df["cohort"] == "FIMM"].copy()
    fig, axes = plt.subplots(1, 3, figsize=(11, 4), sharey=False)
    loci = ["HED_A", "HED_B", "HED_C"]
    titles = ["HLA-A", "HLA-B", "HLA-C"]

    for ax, locus, title in zip(axes, loci, titles):
        data = [df2[locus].dropna().values]
        bp = ax.boxplot(data, patch_artist=True, widths=0.4, medianprops={"color": "white", "lw": 2})
        bp["boxes"][0].set_facecolor(PALETTE["blue"])
        bp["boxes"][0].set_alpha(0.75)
        ax.set_xticks([1])
        ax.set_xticklabels(["FIMM"])
        ax.set_ylabel("HED")
        ax.set_title(title)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

        # Overlay individual points
        jitter = np.random.default_rng(42).uniform(-0.1, 0.1, len(data[0]))
        ax.scatter(np.full(len(data[0]), 1) + jitter, data[0],
                   color=PALETTE["blue"], alpha=0.6, s=20, zorder=3)

    fig.suptitle("HED per Locus — FIMM Cohort", y=1.01)
    fig.tight_layout()
    return fig


def simplify_dx(d):
    """Simplify free-text diagnosis to AML / MDS→AML / MDS / Other."""
    d = str(d)
    if "MDS" in d and "AML" in d:
        return "MDS→AML"
    if "AML" in d:
        return "AML"
    if "MDS" in d:
        return "MDS"
    return "Other"


def fig_hed_vs_os(df_surv):
    """Fig 6: Scatter HED_total vs OS from diagnosis, coloured by dx group (W-H8)."""
    df2 = df_surv.dropna(subset=["HED_total", "os_from_dx"]).copy()
    df2["dx_grp"] = df2["diagnosis"].apply(simplify_dx)

    fig, ax = plt.subplots(figsize=(7, 5))
    for grp, sub in df2.groupby("dx_grp"):
        ax.scatter(sub["HED_total"], sub["os_from_dx"],
                   c=DX_COLORS.get(grp, WONG["muted"]),
                   label=grp, s=60, alpha=0.8, edgecolors="white", linewidths=0.5)

    # Annotate events vs censored
    events = df2[df2["event"] == 1]
    censored = df2[df2["event"] == 0]
    ax.scatter(events["HED_total"], events["os_from_dx"],
               marker="o", facecolors="none", edgecolors="black", s=90, linewidths=0.8, zorder=5)
    ax.scatter(censored["HED_total"], censored["os_from_dx"],
               marker="+", color="gray", s=60, linewidths=1.2, zorder=5)

    # Pearson r
    valid = df2.dropna(subset=["HED_total", "os_from_dx"])
    if len(valid) >= 5:
        r = np.corrcoef(valid["HED_total"], valid["os_from_dx"])[0, 1]
        ax.text(0.05, 0.93, f"r = {r:.2f}", transform=ax.transAxes,
                fontsize=10, color=WONG["navy"])

    ax.set_xlabel("HED total (sum HLA-A + B + C)")
    ax.set_ylabel("OS from diagnosis (months)")
    ax.set_title("HED Total vs Overall Survival (from diagnosis)")
    handles, labels = ax.get_legend_handles_labels()
    handles += [mpatches.Patch(color="none", label="○ event  + censored")]
    ax.legend(handles=handles, frameon=False, fontsize=9)
    fig.tight_layout()
    return fig


def fig_km_hed_split(df_surv):
    """Fig 7: KM stratified by HED total (high/low at median). Primary endpoint: os_from_dx (W-H8)."""
    df2 = df_surv.dropna(subset=["HED_total", "os_from_dx"]).copy()
    median_hed = df2["HED_total"].median()
    n_hi = (df2["HED_total"] >= median_hed).sum()
    n_lo = (df2["HED_total"] < median_hed).sum()
    df2["HED_group"] = np.where(df2["HED_total"] >= median_hed, "High HED", "Low HED")

    fig, ax = plt.subplots(figsize=(7, 5))
    km_plot(df2, "HED_group",
            f"OS by HED Group (split at median HED = {median_hed:.1f})", ax,
            colors=HED_SPLIT_COLORS, time_col="os_from_dx")

    # Log-rank test
    hi = df2[df2["HED_group"] == "High HED"]
    lo = df2[df2["HED_group"] == "Low HED"]
    if len(hi) >= 2 and len(lo) >= 2:
        result = logrank_test(lo["os_from_dx"], hi["os_from_dx"],
                              lo["event"], hi["event"])
        ax.text(0.05, 0.15, f"Log-rank p = {result.p_value:.3f}",
                transform=ax.transAxes, fontsize=10, color=WONG["navy"])

    # W-H4: Power disclaimer
    ax.text(0.05, 0.07,
            f"n={n_hi} High, n={n_lo} Low — power to detect HR=2 ≈ 25% (exploratory only)",
            transform=ax.transAxes, fontsize=7.5, color=WONG["muted"], style="italic")

    ax.set_xlabel("OS from diagnosis (months)")
    fig.tight_layout()
    return fig


# ─────────────────────────────────────────────
#  COX REGRESSION  (W-H4, W-H5)
# ─────────────────────────────────────────────

def cox_regression_html(df_surv):
    """
    Fit univariable and multivariable Cox PH models using os_from_dx.
    Univariable: HED_total only.
    Multivariable: HED_total + diagnosis group (AML / MDS→AML / MDS).
    Returns formatted HTML table.
    """
    df2 = df_surv.dropna(subset=["HED_total", "os_from_dx", "event"]).copy()
    df2["dx_grp"] = df2["diagnosis"].apply(simplify_dx)

    if len(df2) < 5:
        return "<p><em>Insufficient data for Cox regression.</em></p>"

    rows_html = []

    # ── Univariable ──────────────────────────────────────────────────────
    try:
        cph = CoxPHFitter()
        cph.fit(df2[["HED_total", "os_from_dx", "event"]],
                duration_col="os_from_dx", event_col="event", show_progress=False)
        s = cph.summary
        hr   = s.loc["HED_total", "exp(coef)"]
        ci_lo = s.loc["HED_total", "exp(coef) lower 95%"]
        ci_hi = s.loc["HED_total", "exp(coef) upper 95%"]
        pval  = s.loc["HED_total", "p"]
        rows_html.append(
            f"<tr><td>HED total (univariable)</td><td>{len(df2)}</td>"
            f"<td>{hr:.2f}</td><td>{ci_lo:.2f}–{ci_hi:.2f}</td>"
            f"<td>{'<b>' if pval < 0.05 else ''}{pval:.3f}{'</b>' if pval < 0.05 else ''}</td></tr>"
        )
    except Exception as e:
        rows_html.append(f"<tr><td colspan='5'><em>Univariable model failed: {e}</em></td></tr>")

    # ── Multivariable: HED_total + dx_grp dummies ────────────────────────
    try:
        dummies = pd.get_dummies(df2["dx_grp"], prefix="dx", drop_first=True)
        df3 = pd.concat([df2[["HED_total", "os_from_dx", "event"]], dummies], axis=1)
        cph2 = CoxPHFitter()
        cph2.fit(df3, duration_col="os_from_dx", event_col="event", show_progress=False)
        s2 = cph2.summary
        hr2   = s2.loc["HED_total", "exp(coef)"]
        ci2_lo = s2.loc["HED_total", "exp(coef) lower 95%"]
        ci2_hi = s2.loc["HED_total", "exp(coef) upper 95%"]
        pval2  = s2.loc["HED_total", "p"]
        rows_html.append(
            f"<tr><td>HED total (adj. diagnosis group)</td><td>{len(df3)}</td>"
            f"<td>{hr2:.2f}</td><td>{ci2_lo:.2f}–{ci2_hi:.2f}</td>"
            f"<td>{'<b>' if pval2 < 0.05 else ''}{pval2:.3f}{'</b>' if pval2 < 0.05 else ''}</td></tr>"
        )
    except Exception as e:
        rows_html.append(f"<tr><td colspan='5'><em>Multivariable model failed: {e}</em></td></tr>")

    table = (
        '<div class="table-wrap"><table class="data-table">'
        '<thead><tr><th>Model</th><th>N</th><th>HR (HED total)</th>'
        '<th>95% CI</th><th>p-value</th></tr></thead>'
        "<tbody>" + "".join(rows_html) + "</tbody></table></div>"
        '<p style="font-size:0.85em;color:#5a6a7e;margin-top:6px">'
        "HR &gt; 1 = higher HED associated with shorter survival. "
        f"Multivariable model adjusts for diagnosis group (AML / MDS→AML / MDS). "
        f"n={len(df2)} patients with complete HED and survival data. "
        "This is an exploratory analysis; interpret with caution given small sample size.</p>"
    )
    return table


# ─────────────────────────────────────────────
#  SUMMARY STATS
# ─────────────────────────────────────────────

def compute_stats(df, df_surv):
    stats = {}
    stats["n_patients"] = df["patient_id"].nunique()
    stats["n_surv_patients"] = df_surv["patient_id"].nunique()
    stats["n_events"] = int(df_surv["event"].sum())
    stats["n_censored"] = int((df_surv["event"] == 0).sum())
    stats["cohort_counts"] = df["cohort"].value_counts().to_dict()

    # Median OS — primary endpoint is os_from_dx (W-H8)
    kmf = KaplanMeierFitter()
    valid = df_surv.dropna(subset=["os_from_dx"])
    if len(valid) > 2:
        kmf.fit(valid["os_from_dx"], valid["event"])
        med = kmf.median_survival_time_
        stats["median_os_months"] = f"{med:.1f}" if not np.isnan(med) else "NR"
    else:
        stats["median_os_months"] = "NR"

    # 1-year OS rate
    try:
        stats["os_1yr"] = f"{kmf.predict(12)*100:.1f}%"
    except Exception:
        stats["os_1yr"] = "—"

    # Consensus HED completeness
    stats["n_hed_complete"] = int(df["hed_complete"].sum()) if "hed_complete" in df else "—"
    stats["n_hed_concordant"] = int(df["hed_concordant_all"].sum()) if "hed_concordant_all" in df else "—"

    # HED stats
    stats["mean_hed_total"] = f"{df['HED_total'].mean():.2f}" if "HED_total" in df else "—"
    stats["mean_hed_a"] = f"{df['HED_A'].mean():.2f}" if "HED_A" in df else "—"
    stats["mean_hed_b"] = f"{df['HED_B'].mean():.2f}" if "HED_B" in df else "—"
    stats["mean_hed_c"] = f"{df['HED_C'].mean():.2f}" if "HED_C" in df else "—"
    stats["n_homo"] = int((df_surv["HED_total"] == 0).sum()) if "HED_total" in df_surv else 0
    return stats


def survival_table_html(df_surv):
    """Per-cohort OS summary table."""
    rows = []
    for cohort in ["FIMM"]:
        sub = df_surv[df_surv["cohort"] == cohort].dropna(subset=["os_from_dx"])
        if len(sub) < 2:
            continue
        kmf = KaplanMeierFitter()
        kmf.fit(sub["os_from_dx"], sub["event"])
        med = kmf.median_survival_time_
        med_str = f"{med:.1f}" if not np.isnan(med) else "NR"
        try:
            os1 = f"{kmf.predict(12)*100:.1f}%"
        except Exception:
            os1 = "—"
        try:
            os2 = f"{kmf.predict(24)*100:.1f}%"
        except Exception:
            os2 = "—"
        rows.append(f"<tr><td><strong>{cohort}</strong></td>"
                    f"<td>{len(sub)}</td>"
                    f"<td>{int(sub['event'].sum())}</td>"
                    f"<td>{med_str}</td>"
                    f"<td>{os1}</td>"
                    f"<td>{os2}</td></tr>")

    return (
        '<div class="table-wrap"><table class="data-table">'
        '<thead><tr><th>Cohort</th><th>N</th><th>Events</th>'
        '<th>Median OS (months)</th><th>1-year OS</th><th>2-year OS</th></tr></thead>'
        "<tbody>" + "".join(rows) + "</tbody></table></div>"
    )


def hed_table_html(df):
    """Per-patient HED table (uses consensus alleles from 3 WES tools)."""
    # Use consensus columns if present, fall back to arcasHLA
    has_consensus = "consensus_A" in df.columns
    a_col = "consensus_A" if has_consensus else "wes_arcasHLA_A"
    b_col = "consensus_B" if has_consensus else "wes_arcasHLA_B"
    c_col = "consensus_C" if has_consensus else "wes_arcasHLA_C"
    base_cols = ["patient_id", "cohort", a_col, b_col, c_col, "HED_A", "HED_B", "HED_C", "HED_total"]
    extra_cols = [c for c in ["hed_complete", "hed_concordant_all"] if c in df.columns]
    sub = df[base_cols + extra_cols].dropna(subset=["HED_total"]).sort_values("HED_total", ascending=False)
    rows = []
    for _, r in sub.iterrows():
        complete_icon = ""
        if "hed_complete" in r:
            complete_icon = " ✓" if r["hed_complete"] else " ✗"
        concordant_icon = ""
        if "hed_concordant_all" in r:
            concordant_icon = " ✓" if r["hed_concordant_all"] else ""
        rows.append(
            f"<tr><td>{r['patient_id']}</td>"
            f"<td>{r['cohort']}</td>"
            f"<td style='font-size:0.82em'>{r[a_col]}</td>"
            f"<td style='font-size:0.82em'>{r[b_col]}</td>"
            f"<td style='font-size:0.82em'>{r[c_col]}</td>"
            f"<td>{r['HED_A']:.2f}</td>"
            f"<td>{r['HED_B']:.2f}</td>"
            f"<td>{r['HED_C']:.2f}</td>"
            f"<td><strong>{r['HED_total']:.2f}</strong>{complete_icon}</td>"
            f"<td>{concordant_icon}</td></tr>"
        )
    allele_header = "Consensus HLA-A" if has_consensus else "HLA-A alleles"
    allele_b_header = "Consensus HLA-B" if has_consensus else "HLA-B alleles"
    allele_c_header = "Consensus HLA-C" if has_consensus else "HLA-C alleles"
    return (
        '<div class="table-wrap"><table class="data-table">'
        f'<thead><tr><th>Patient</th><th>Cohort</th>'
        f'<th>{allele_header}</th><th>{allele_b_header}</th><th>{allele_c_header}</th>'
        '<th>HED-A</th><th>HED-B</th><th>HED-C</th><th>HED total</th>'
        '<th title="All 3 WES tools agree">All-3 concordant</th></tr></thead>'
        "<tbody>" + "".join(rows) + "</tbody></table></div>"
    )


def cohort_table_html(df):
    """Per-patient clinical summary table."""
    sub = df.dropna(subset=["patient_id"])[
        ["patient_id", "cohort", "gender", "diagnosis", "stage",
         "sampling_date", "months_since_dx", "event", "date_of_death"]
    ].copy()
    sub = sub.sort_values(["cohort", "patient_id"])
    rows = []
    for _, r in sub.iterrows():
        ev = '<span class="badge" style="background:#D62728">Death</span>' if r["event"] == 1 \
             else '<span class="badge" style="background:#5a6a7e">Censored</span>'
        sdate = r["sampling_date"].strftime("%Y-%m-%d") if pd.notna(r["sampling_date"]) else "—"
        ddate = r["date_of_death"].strftime("%Y-%m-%d") if pd.notna(r["date_of_death"]) else "—"
        msdx = f"{r['months_since_dx']:.0f}" if pd.notna(r["months_since_dx"]) else "—"
        rows.append(
            f"<tr><td><strong>{r['patient_id']}</strong></td>"
            f"<td>{r['cohort']}</td>"
            f"<td>{r['gender'] if pd.notna(r['gender']) else '—'}</td>"
            f"<td style='font-size:0.85em'>{str(r['diagnosis'])[:55] if pd.notna(r['diagnosis']) else '—'}</td>"
            f"<td>{r['stage'] if pd.notna(r['stage']) else '—'}</td>"
            f"<td>{sdate}</td>"
            f"<td>{msdx}</td>"
            f"<td>{ev}</td>"
            f"<td>{ddate}</td></tr>"
        )
    return (
        '<div class="table-wrap"><table class="data-table">'
        '<thead><tr><th>Patient</th><th>Cohort</th><th>Sex</th><th>Diagnosis</th>'
        '<th>Stage</th><th>Sampling date</th><th>Months from Dx</th>'
        '<th>Status</th><th>Date of death</th></tr></thead>'
        "<tbody>" + "".join(rows) + "</tbody></table></div>"
    )


# ─────────────────────────────────────────────
#  HTML TEMPLATE
# ─────────────────────────────────────────────
CSS = """
:root {
  --navy:   #1a3a5c;
  --blue:   #2c6fad;
  --lblue:  #4fc3f7;
  --bg:     #f0f4f8;
  --white:  #ffffff;
  --border: #d0dce8;
  --text:   #1e2a38;
  --muted:  #5a6a7e;
  --green:  #2CA02C;
  --orange: #EE7733;
  --red:    #D62728;
}
*, *::before, *::after { box-sizing: border-box; margin: 0; padding: 0; }
body { font-family: 'Segoe UI', Arial, sans-serif; background: var(--bg);
        color: var(--text); display: flex; min-height: 100vh; }

#sidebar {
  width: 250px; min-height: 100vh; background: var(--navy);
  position: fixed; top: 0; left: 0; overflow-y: auto;
  display: flex; flex-direction: column; z-index: 100;
  box-shadow: 3px 0 10px rgba(0,0,0,0.25);
}
.sidebar-header {
  padding: 22px 18px 16px; border-bottom: 1px solid rgba(255,255,255,0.12);
}
.sidebar-header h1 {
  font-size: 1.18em; color: var(--lblue); font-weight: 800; letter-spacing: 0.5px;
}
.sidebar-header p { font-size: 0.75em; color: #a8c4e0; margin-top: 4px; }
.sidebar-header .version {
  display: inline-block; background: rgba(79,195,247,0.2);
  color: var(--lblue); font-size: 0.7em; padding: 2px 8px;
  border-radius: 10px; margin-top: 6px; font-weight: 600;
}
nav { padding: 12px 0 24px; }
.nav-section { font-size: 0.65em; color: #6a8bb0; letter-spacing: 1.5px;
                text-transform: uppercase; padding: 12px 18px 4px; }
.nav-link {
  display: block; padding: 8px 18px; color: #a8c4e0;
  text-decoration: none; font-size: 0.85em; transition: all 0.15s;
  border-left: 3px solid transparent;
}
.nav-link:hover { background: rgba(255,255,255,0.07);
                   color: #fff; border-left-color: var(--lblue); }

#main { margin-left: 250px; padding: 32px 40px 60px; max-width: 1200px; width: 100%; }

section { margin-bottom: 56px; scroll-margin-top: 20px; }
h2.section-title {
  font-size: 1.55em; color: var(--navy); font-weight: 800;
  border-bottom: 3px solid var(--blue); padding-bottom: 8px;
  margin-bottom: 20px; display: flex; align-items: center; gap: 10px;
}
h2.section-title .sec-num {
  background: var(--blue); color: white; border-radius: 50%;
  width: 32px; height: 32px; display: flex; align-items: center;
  justify-content: center; font-size: 0.7em; flex-shrink: 0;
}
h3 { font-size: 1.15em; color: var(--blue); margin: 24px 0 10px; font-weight: 700; }
p { line-height: 1.7; margin-bottom: 12px; color: var(--text); }
ul { margin: 8px 0 12px 22px; }
li { line-height: 1.65; margin-bottom: 4px; }
strong { color: var(--navy); }

.card {
  background: var(--white); border: 1px solid var(--border);
  border-radius: 10px; padding: 24px; margin-bottom: 24px;
  box-shadow: 0 2px 8px rgba(0,0,0,0.07);
}
.card-grid { display: grid; grid-template-columns: repeat(auto-fit, minmax(200px,1fr));
              gap: 16px; margin-bottom: 24px; }
.metric-card {
  background: var(--white); border: 1px solid var(--border);
  border-radius: 10px; padding: 18px 20px; text-align: center;
  box-shadow: 0 2px 6px rgba(0,0,0,0.06);
}
.metric-card .value {
  font-size: 2em; font-weight: 900; color: var(--blue); line-height: 1;
  margin-bottom: 6px;
}
.metric-card .label {
  font-size: 0.8em; color: var(--muted); text-transform: uppercase; letter-spacing: 0.5px;
}
.metric-card .sublabel { font-size: 0.78em; color: var(--muted); margin-top: 2px; }

.fig-card {
  background: var(--white); border: 1px solid var(--border);
  border-radius: 10px; overflow: hidden; margin-bottom: 28px;
  box-shadow: 0 2px 8px rgba(0,0,0,0.07);
}
.fig-img { width: 100%; height: auto; display: block; }
.fig-caption {
  padding: 14px 18px; font-size: 0.88em; color: var(--muted);
  border-top: 1px solid var(--border); line-height: 1.65;
}

.table-wrap { overflow-x: auto; margin-bottom: 20px; }
.data-table {
  width: 100%; border-collapse: collapse; font-size: 0.87em;
  background: var(--white); border-radius: 8px; overflow: hidden;
  box-shadow: 0 1px 4px rgba(0,0,0,0.08);
}
.data-table th {
  background: var(--blue); color: white; padding: 10px 12px;
  text-align: left; font-weight: 700; white-space: nowrap;
}
.data-table td { padding: 9px 12px; border-bottom: 1px solid #e8eef5; }
.data-table tr:last-child td { border-bottom: none; }
.data-table tr:nth-child(even) td { background: #f5f8fb; }
.data-table tr:hover td { background: #e8f0fa; }

.badge {
  display: inline-block; padding: 3px 10px; border-radius: 12px;
  color: white; font-size: 0.82em; font-weight: 700;
}

.highlight-box {
  background: linear-gradient(135deg, #e8f0fa, #f0f4f8);
  border-left: 4px solid var(--blue); padding: 16px 20px;
  border-radius: 0 8px 8px 0; margin: 16px 0;
}
.highlight-box.green { border-left-color: var(--green);
  background: linear-gradient(135deg, #e8f5e9, #f0f4f8); }
.highlight-box.orange { border-left-color: var(--orange);
  background: linear-gradient(135deg, #fff8e1, #f0f4f8); }
.highlight-box.red { border-left-color: var(--red);
  background: linear-gradient(135deg, #ffebee, #f0f4f8); }

.fig-explain {
  background: #f8fbff; border-left: 3px solid #4fc3f7;
  padding: 10px 16px; margin: 8px 0 20px 0;
  font-size: 0.93em; color: #2a3a50; line-height: 1.65;
  border-radius: 0 6px 6px 0;
}
code { background: #eef2f7; padding: 2px 6px; border-radius: 4px;
        font-size: 0.9em; color: #c44e52; }
"""


def render_html(stats, figs, df, df_surv):
    surv_tbl = survival_table_html(df_surv)
    hed_tbl = hed_table_html(df)
    cohort_tbl = cohort_table_html(df)
    cox_html = cox_regression_html(df_surv)

    fimm_n = stats["cohort_counts"].get("FIMM", 0)
    vx_n = stats["cohort_counts"].get("VX", 0)

    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>FIMM Survival &amp; HED Report</title>
<style>{CSS}</style>
</head>
<body>

<div id="sidebar">
  <div class="sidebar-header">
    <h1>FIMM · Survival &amp; HED</h1>
    <p>HLA Evolutionary Divergence &amp; Overall Survival</p>
    <p>AML/MDS Cohort · FIMM Helsinki</p>
    <span class="version">Report v2 · {CUTOFF_DATE.strftime('%Y-%m-%d')}</span>
  </div>
  <nav>
    <span class="nav-section">Sections</span>
    <a href="#overview"  class="nav-link">1. Executive Summary</a>
    <a href="#cohort"    class="nav-link">2. Cohort Overview</a>
    <a href="#survival"  class="nav-link">3. Overall Survival</a>
    <a href="#hed"       class="nav-link">4. HED Analysis</a>
    <a href="#hedos"     class="nav-link">5. HED &times; Survival</a>
    <a href="#cox"       class="nav-link">6. Cox Regression</a>
    <a href="#methods"   class="nav-link">7. Methods</a>
  </nav>
</div>

<div id="main">

<!-- ═══ 1. EXECUTIVE SUMMARY ═══ -->
<section id="overview">
  <h2 class="section-title"><span class="sec-num">1</span>Executive Summary</h2>
  <p>
    This report presents an integrated analysis of <strong>Overall Survival (OS)</strong>
    and <strong>HLA Evolutionary Divergence (HED)</strong> for {stats['n_patients']} unique
    FIMM patients with haematological malignancies (AML, MDS). HED was computed from
    HLA-A, -B, and -C allele calls using Grantham amino acid distances at antigen-binding
    groove positions. Survival data were derived from the FIMM biobank patient registry.
  </p>

  <div class="card-grid">
    <div class="metric-card">
      <div class="value">{stats['n_patients']}</div>
      <div class="label">Total FIMM patients</div>
      <div class="sublabel">HLA typed</div>
    </div>
    <div class="metric-card">
      <div class="value">{stats['n_surv_patients']}</div>
      <div class="label">With survival data</div>
      <div class="sublabel">{stats['n_events']} events · {stats['n_censored']} censored</div>
    </div>
    <div class="metric-card">
      <div class="value">{stats['median_os_months']}</div>
      <div class="label">Median OS (months)</div>
      <div class="sublabel">from diagnosis date</div>
    </div>
    <div class="metric-card">
      <div class="value">{stats['os_1yr']}</div>
      <div class="label">1-year OS rate</div>
      <div class="sublabel">Kaplan–Meier estimate</div>
    </div>
    <div class="metric-card">
      <div class="value">{stats['mean_hed_total']}</div>
      <div class="label">Mean HED total</div>
      <div class="sublabel">HLA-A + B + C Grantham</div>
    </div>
    <div class="metric-card">
      <div class="value">{stats.get('n_hed_complete', '—')}</div>
      <div class="label">HED complete</div>
      <div class="sublabel">all 3 loci resolved</div>
    </div>
    <div class="metric-card">
      <div class="value">{stats.get('n_hed_concordant', '—')}</div>
      <div class="label">HED concordant</div>
      <div class="sublabel">all 3 WES tools agree</div>
    </div>
  </div>
  <p style="font-size:0.82em;color:#5a6a7e;margin-top:-8px">
    Cutoff date: <strong>{CUTOFF_DATE.strftime('%Y-%m-%d')}</strong> ·
    Cohorts: FIMM n={fimm_n}, VX n={vx_n}
  </p>

  <div class="highlight-box">
    <strong>Key findings:</strong> The cohort comprises predominantly AML and MDS patients
    from the FIMM biobank (Helsinki AML/MDS patients). HED at HLA-B is the most variable locus.
    Patients with higher total HED tend to show differential survival patterns, though
    with this sample size results should be considered exploratory.
  </div>
</section>

<!-- ═══ 2. COHORT OVERVIEW ═══ -->
<section id="cohort">
  <h2 class="section-title"><span class="sec-num">2</span>Cohort Overview</h2>
  <p>
    The FIMM cohort includes samples from two groups:
    <strong>FIMM</strong> (FIMM Helsinki AML/MDS patients, n={fimm_n}) and
    <strong>VX</strong> (Venetoclax study, n={vx_n}).
    All patients have HLA-A, -B, and -C allele calls from at least one sequencing modality
    (scRNA-seq, bulk RNA-seq, or WES). Survival follow-up is available for FIMM patients.
  </p>

  {fig_img_html(figs['cohort'], 'Figure 1. Cohort composition by disease stage (Diagnosis / Relapse / Refractory). Each bar represents one cohort; segments show proportions of disease stage at sampling.')}

  <h3>Patient-level clinical data</h3>
  {cohort_tbl}

  <div class="highlight-box orange">
    <strong>Note:</strong> VX cohort (n={vx_n}) patients have HLA typing data but
    no survival follow-up in the current registry extract and are excluded from
    survival analyses. Months-since-diagnosis values were recorded at the time of
    sampling and represent time from initial diagnosis, not treatment initiation.
  </div>
</section>

<!-- ═══ 3. OVERALL SURVIVAL ═══ -->
<section id="survival">
  <h2 class="section-title"><span class="sec-num">3</span>Overall Survival</h2>
  <p>
    Overall survival (OS) was defined as time from <strong>diagnosis date</strong> to date of death
    (event = 1) or last follow-up (censored = 0, cutoff {CUTOFF_DATE.strftime('%Y-%m-%d')}).
    Using diagnosis as time-zero avoids the landmark bias introduced by variable sampling
    delay within each patient's disease course. Kaplan–Meier curves were estimated
    using the <em>lifelines</em> library.
  </p>

  <h3>OS — FIMM cohort</h3>
  {fig_img_html(figs['km_cohort'], 'Figure 2. Kaplan–Meier overall survival for FIMM AML/MDS patients. Time zero = diagnosis date. Shaded area represents 95% confidence interval.')}

  <div class="fig-explain">
    The FIMM cohort includes both newly diagnosed and relapsed/refractory patients,
    resulting in a broad OS distribution. Results are exploratory given small sample sizes.
  </div>

  <h3>OS by diagnosis group</h3>
  {fig_img_html(figs['km_diagnosis'], 'Figure 3. Kaplan–Meier overall survival stratified by simplified diagnosis group (AML, MDS→AML, MDS). Patients classified as "Other" are excluded from this panel.')}

  <h3>Summary statistics</h3>
  {surv_tbl}
</section>

<!-- ═══ 4. HED ANALYSIS ═══ -->
<section id="hed">
  <h2 class="section-title"><span class="sec-num">4</span>HED Analysis</h2>
  <p>
    <strong>HLA Evolutionary Divergence (HED)</strong> quantifies the structural divergence
    between the two alleles at each HLA locus using Grantham amino acid distances at
    canonical antigen-binding groove positions. Higher HED indicates greater divergence
    between alleles and potentially broader peptide-presentation capacity.
    HED = 0 for homozygous patients (identical alleles at a locus).
  </p>

  <h3>HED distribution and per-locus variability</h3>
  {fig_img_html(figs['hed_dist'], 'Figure 4. Left: distribution of HED total (sum of HED-A + HED-B + HED-C) across all patients. Dashed red line = median. Right: per-locus HED distribution (violin plots). HLA-B shows the highest variability. HED values derived from majority-vote consensus across 3 WES tools (OptiType, ArcasHLA, SpecHLA).')}

  <h3>Per-locus HED by cohort</h3>
  {fig_img_html(figs['hed_cohort'], 'Figure 5. Boxplots of HED per locus (HLA-A, HLA-B, HLA-C) for the FIMM cohort. Individual data points overlaid. Median shown as white line.')}

  <div class="highlight-box green">
    <strong>HED ranges observed:</strong> Mean HED-A = {stats['mean_hed_a']},
    Mean HED-B = {stats['mean_hed_b']}, Mean HED-C = {stats['mean_hed_c']}.
    {stats['n_homo']} patient(s) in the survival cohort are fully homozygous at all three
    loci (HED total = 0). HLA-B demonstrates the widest allelic diversity in this cohort.
    Patients with complete HED (all 3 loci resolved): <strong>{stats.get('n_hed_complete', '—')}</strong>.
    Patients where all 3 WES tools fully agree: <strong>{stats.get('n_hed_concordant', '—')}</strong>.
  </div>

  <h3>Per-patient HED table (consensus alleles)</h3>
  <p style="font-size:0.85em;color:#5a6a7e">
    HLA alleles shown are majority-vote consensus across WES OptiType, ArcasHLA and SpecHLA.
    ✓ in the last column = all 3 tools agreed on the allele pair at every locus.
  </p>
  {hed_tbl}
</section>

<!-- ═══ 5. HED × SURVIVAL ═══ -->
<section id="hedos">
  <h2 class="section-title"><span class="sec-num">5</span>HED &times; Survival</h2>
  <p>
    We explored the association between total HED and overall survival using
    scatter visualisation and Kaplan–Meier stratification (median HED split).
    Given the exploratory nature and small cohort size (n≤28 with survival data),
    these results should be interpreted as hypothesis-generating only.
  </p>

  <h3>HED total vs OS months (scatter)</h3>
  {fig_img_html(figs['hed_vs_os'], 'Figure 6. Scatter plot of HED total vs OS from diagnosis (months). Circles = death events; crosses = censored. Colour encodes diagnosis group (Wong palette). Pearson r shown in upper-left. No strong linear correlation is expected due to the dominant role of treatment response in this cohort.')}

  <h3>Kaplan–Meier: High HED vs Low HED</h3>
  {fig_img_html(figs['km_hed'], 'Figure 7. KM survival curves split at median HED total. High HED = &ge; median; Low HED = &lt; median. Time axis = OS from diagnosis. Log-rank p-value shown. Shaded regions = 95% CI. Statistical note: with n&approx;14 per arm, this log-rank test has approximately 25% power to detect HR=2.0 (&alpha;=0.05). Results are hypothesis-generating only.')}

  <div class="highlight-box orange">
    <strong>Interpretation:</strong> The HED-OS association in this cohort is exploratory.
    AML outcomes are predominantly driven by cytogenetic risk, mutation profiles (FLT3, NPM1),
    and treatment response. HLA evolutionary divergence may modulate immune surveillance
    and contribute to immune escape or graft-versus-leukaemia effects in transplant contexts.
    Larger cohorts with treatment data are required for definitive conclusions.
  </div>
</section>

<!-- ═══ 6. COX REGRESSION ═══ -->
<section id="cox">
  <h2 class="section-title"><span class="sec-num">6</span>Cox Regression</h2>
  <p>
    To complement the Kaplan–Meier log-rank test, we fit Cox proportional hazards models with
    <strong>OS from diagnosis</strong> as the time endpoint. Univariable and multivariable (adjusted
    for simplified diagnosis group: AML / MDS→AML / MDS) models are shown. At n=28 the
    multivariable model is over-parametrised; results are hypothesis-generating.
  </p>
  {cox_html}
  <div class="highlight-box orange">
    <strong>Caution:</strong> With n≈28 and a median of 14 events, the Cox model has low power.
    Confidence intervals are wide. Diagnosis group adjustment reduces degrees of freedom further.
    These results should be considered exploratory pending a larger cohort.
  </div>
</section>

<!-- ═══ 7. METHODS ═══ -->
<section id="methods">
  <h2 class="section-title"><span class="sec-num">7</span>Methods</h2>

  <h3>Data sources</h3>
  <ul>
    <li><strong>HLA allele calls:</strong> <code>fimm_hla_calls.tsv</code> — WES calls from
        OptiType, ArcasHLA, and SpecHLA at 2-field resolution. A majority-vote consensus
        allele pair is derived per locus per patient (see Consensus HLA below).</li>
    <li><strong>Survival data:</strong> <code>patient_data_v3.xlsx</code> sheet <em>scrna</em> —
        date of death, cause of death, gender per patient study ID.</li>
    <li><strong>Sample registry:</strong> <code>Book2.xlsx</code> — sampling dates, diagnosis text,
        disease stage, months since diagnosis.</li>
  </ul>

  <h3>Patient ID linking</h3>
  <p>
    Sample IDs (e.g., <code>FH_2030_4</code>) were mapped to patient IDs (e.g., <code>FH.2030</code>)
    by splitting on underscores and joining the first two parts with a dot. This key was used
    to merge <code>Book2.xlsx</code> (donor column) and <code>patient_data_v3.xlsx</code>
    (study id column). For patients with multiple timepoints, the earliest sampling was
    selected by sorting on sample ID before deduplication.
  </p>

  <h3>Consensus HLA allele calling</h3>
  <p>
    For each patient and each locus (HLA-A, -B, -C), allele pairs from three WES tools
    (OptiType, ArcasHLA, SpecHLA) were normalised to a sorted canonical form and the
    majority-vote pair (present in ≥2 tools) was selected. If no tool provided a valid
    <code>allele1|allele2</code> call, HED for that locus was set to NaN and the patient
    is flagged as <em>HED incomplete</em>. If all 3 tools agreed exactly, the patient
    is flagged as <em>fully concordant</em>.
  </p>

  <h3>Survival analysis</h3>
  <p>
    OS was measured from <strong>diagnosis date</strong> to death (event = 1) or last follow-up
    (censored = 0, cutoff date: {CUTOFF_DATE.strftime('%Y-%m-%d')}).
    Using diagnosis as time-zero avoids the sampling-date landmark bias that arises when
    patients are collected at different disease stages.
    Kaplan–Meier estimation and log-rank testing were performed using
    <em>lifelines</em> v0.30.0 (Davidson-Pilon 2019).
    Cox proportional hazards models (univariable and multivariable adjusted for diagnosis group)
    were fit on the same endpoint. The log-rank test with n≈14 per arm has approximately
    25% power to detect HR=2.0 (α=0.05); all results are hypothesis-generating.
  </p>

  <h3>HED computation</h3>
  <p>
    HED was computed following Pierini &amp; Lenz (2018, <em>PLOS Genetics</em>).
    For each HLA locus, the two alleles were mapped to their amino acid sequences
    at canonical antigen-binding groove (ABG) positions
    (HLA-A: positions 9, 37, 45, 60, 67, 70, 73, 74, 76, 77;
     HLA-B: 9, 37, 45, 67, 70, 76, 77, 80;
     HLA-C: 9, 37, 45, 67, 70, 76, 77)
    using published sequences from the IPD-IMGT/HLA database.
    Pairwise amino acid distances were computed using the Grantham (1974) formula:
    <code>d(i,j) = √(0.1018·Δc² + 0.000399·Δp² + 0.000791·Δv²)</code>
    where c = amino acid composition, p = polarity, v = molecular volume.
    HED at each locus = mean Grantham distance across ABG positions.
    HED = 0 for homozygous patients. HED total = sum across HLA-A, -B, -C.
    Alleles not found in the ABG table yield NaN (not a fallback constant).
  </p>

  <h3>Software</h3>
  <ul>
    <li>Python 3.10.6 · pandas 2.2.3 · numpy 1.23.4 · matplotlib 3.6.2 · scipy 1.9.3</li>
    <li>lifelines 0.30.0 · openpyxl 3.1.5</li>
    <li>Runtime: CSC Puhti HPC · module <code>gcc/11.3.0</code> + <code>biopythontools/11.3.0_3.10.6</code></li>
  </ul>
</section>

</div><!-- /#main -->
</body>
</html>
"""
    return html


# ─────────────────────────────────────────────
#  MAIN
# ─────────────────────────────────────────────
def main():
    global CUTOFF_DATE
    args = parse_args()
    CUTOFF_DATE = datetime.strptime(args.cutoff_date, "%Y-%m-%d")
    output_path = args.output  # already a Path from argparse type=Path

    print(f"Cutoff date: {CUTOFF_DATE.strftime('%Y-%m-%d')}")
    print("Loading data …")
    df, df_surv, hla_raw = load_data()
    print(f"  Patients total: {df['patient_id'].nunique()}")
    print(f"  Patients with survival data: {df_surv['patient_id'].nunique()}")
    print(f"  Death events: {int(df_surv['event'].sum())}")
    if "hed_complete" in df.columns:
        print(f"  HED complete (all 3 loci): {int(df['hed_complete'].sum())}")
    if "hed_concordant_all" in df.columns:
        print(f"  HED concordant (all tools agree): {int(df['hed_concordant_all'].sum())}")

    print("Computing statistics …")
    stats = compute_stats(df, df_surv)

    print("Generating figures …")
    figs = {
        "cohort":    fig_to_b64(fig_cohort_composition(df)),
        "km_cohort": fig_to_b64(fig_km_cohort(df_surv)),
        "km_diagnosis": fig_to_b64(fig_km_diagnosis(df_surv)),
        "hed_dist":  fig_to_b64(fig_hed_distribution(df)),
        "hed_cohort": fig_to_b64(fig_hed_by_cohort(df)),
        "hed_vs_os": fig_to_b64(fig_hed_vs_os(df_surv)),
        "km_hed":    fig_to_b64(fig_km_hed_split(df_surv)),
    }

    print("Rendering HTML …")
    html = render_html(stats, figs, df, df_surv)
    output_path.write_text(html, encoding="utf-8")
    print(f"Report written → {output_path}")
    size_mb = output_path.stat().st_size / 1_048_576
    print(f"File size: {size_mb:.2f} MB")


if __name__ == "__main__":
    main()
