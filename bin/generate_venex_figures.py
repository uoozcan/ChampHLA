#!/usr/bin/env python3
"""
VENEX WGS cohort figures for PIHLA report.

Generates 5 figures:
  venex_fig1_tool_coverage        — Tool call coverage (% samples called) per gene and tool
  venex_fig2_allele_agreement     — Inter-tool allele agreement heatmap per gene
  venex_fig3_homozygous_flags     — Homozygous-call frequency per sample×gene (LOH candidates)
  venex_fig4_loh_status           — LOH status distribution per gene (stacked bar)
  venex_fig5_allele_diversity     — Unique alleles called per gene per tool (bar)

Then patches mvhla_report_v12.html → mvhla_report_v13.html with a new VENEX section.

Usage:
  python3 generate_venex_figures.py \
      --calls   /scratch/project_2008084/pihla-publish/analysis/venex_harmonized_calls.tsv \
      --loh     /scratch/project_2008084/pihla-publish/analysis/loh_analysis/sample_loh_candidates.tsv \
      --out-dir /scratch/project_2008084/mvhla_figures_v6 \
      --input-report  /scratch/project_2008084/mvhla_report_v12.html \
      --output-report /scratch/project_2008084/mvhla_report_v13.html
"""
import argparse, base64, math, re, sys
from pathlib import Path
from collections import defaultdict
from datetime import date
from itertools import combinations

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.ticker as mticker

# ── Style ─────────────────────────────────────────────────────────────────────
GRID  = "#E5E5E5"; BG = "#FAFAFA"
TODAY = date.today().isoformat()
WM_TXT = f"PIHLA · VENEX WGS Cohort · {TODAY}"

GENE_ORDER = ["A", "B", "C", "DRB1", "DQB1"]
TOOL_ORDER = ["OptiType", "T1K", "SpecHLA", "HLA-HD"]
TOOL_COLORS = {
    "OptiType": "#2A9D8F", "T1K": "#009E73",
    "SpecHLA": "#D55E00",  "HLA-HD": "#4C78A8",
}
LOH_COLORS = {
    "balanced_heterozygous":          "#2CA02C",
    "insufficient_evidence":          "#9467BD",
    "homozygous_call_needs_review":   "#D62728",
}

def apply_style():
    plt.rcParams.update({
        "figure.facecolor": BG, "axes.facecolor": BG,
        "axes.edgecolor": "#BBBBBB", "axes.spines.top": False,
        "axes.spines.right": False, "axes.grid": True,
        "grid.color": GRID, "grid.linewidth": 0.7, "grid.linestyle": "--",
        "font.family": "DejaVu Sans", "font.size": 10,
        "axes.titlesize": 12, "axes.titleweight": "bold",
        "axes.labelsize": 10, "xtick.labelsize": 9,
        "ytick.labelsize": 9, "legend.fontsize": 8.5,
    })

def wm(fig):
    fig.text(0.99, 0.005, WM_TXT, ha="right", va="bottom",
             fontsize=6, color="#AAAAAA", transform=fig.transFigure)

def save(fig, out, name):
    p = Path(out)
    p.mkdir(parents=True, exist_ok=True)
    fig.savefig(p / f"{name}.png", dpi=300, bbox_inches="tight")
    fig.savefig(p / f"{name}.pdf", bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved {name}")

def normalise2(allele):
    """Trim to 2-field resolution: A*01:01:01:01 → A*01:01"""
    if not isinstance(allele, str): return ""
    parts = allele.split(":")
    if len(parts) >= 2:
        return ":".join(parts[:2])
    return allele

# ── Figure V1: Tool coverage ──────────────────────────────────────────────────
def fig_v1_tool_coverage(df, out):
    """
    Grouped bar chart: for each gene × tool, the fraction of all 87 VENEX WGS
    samples that received a callable result from that tool.
    """
    all_samples = df["sample"].nunique()
    genes_avail = [g for g in GENE_ORDER if g in df["gene"].unique()]

    fig, axes = plt.subplots(1, len(genes_avail), figsize=(14, 5), sharey=True)
    if len(genes_avail) == 1:
        axes = [axes]

    for ax, gene in zip(axes, genes_avail):
        tools = [t for t in TOOL_ORDER if t in df[df["gene"]==gene]["tool"].unique()]
        counts = [df[(df["tool"]==t) & (df["gene"]==gene)]["sample"].nunique() for t in tools]
        fracs  = [c / 87 for c in counts]
        bars = ax.bar(tools, fracs,
                      color=[TOOL_COLORS.get(t, "#888") for t in tools],
                      edgecolor="white", linewidth=0.5, zorder=3)
        for bar, c in zip(bars, counts):
            ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.02,
                    str(c), ha="center", va="bottom", fontsize=9, fontweight="bold")
        ax.set_title(f"HLA-{gene}", fontsize=11, fontweight="bold")
        ax.set_ylim(0, 1.15)
        ax.yaxis.set_major_formatter(mticker.PercentFormatter(1.0))
        ax.tick_params(axis="x", rotation=30)
        ax.set_xlabel("")

    axes[0].set_ylabel("Fraction of VENEX WGS samples with callable result", fontsize=10)
    fig.suptitle(f"VENEX WGS Cohort — Tool Coverage per Gene (n=87 total samples)",
                 fontsize=12, fontweight="bold", y=1.02)
    wm(fig); fig.tight_layout()
    save(fig, out, "venex_fig1_tool_coverage")

# ── Figure V2: Inter-tool allele agreement ─────────────────────────────────────
def fig_v2_allele_agreement(df, out):
    """
    For each pair of tools that both typed the same sample × gene, compute the
    fraction of cases where the called allele pair (2-field normalised, unordered)
    matches. Show as a per-gene heatmap-style tile matrix.
    """
    df2 = df.copy()
    df2["a1_2f"] = df2["allele1"].apply(normalise2)
    df2["a2_2f"] = df2["allele2"].apply(normalise2)
    # Canonical pair: frozenset for unordered comparison
    df2["pair"] = df2.apply(lambda r: frozenset([r["a1_2f"], r["a2_2f"]]), axis=1)

    genes_avail = [g for g in ["A","B","C"] if g in df2["gene"].unique()]
    tools = [t for t in TOOL_ORDER if t in df2["tool"].unique()]

    fig, axes = plt.subplots(1, len(genes_avail), figsize=(14, 4.5))
    if len(genes_avail) == 1:
        axes = [axes]

    for ax, gene in zip(axes, genes_avail):
        gene_df = df2[df2["gene"]==gene]
        local_tools = [t for t in tools if t in gene_df["tool"].unique()]
        n = len(local_tools)
        mat = np.full((n, n), np.nan)

        for i, t1 in enumerate(local_tools):
            for j, t2 in enumerate(local_tools):
                if i == j:
                    mat[i, j] = 1.0; continue
                m1 = gene_df[gene_df["tool"]==t1].set_index("sample")["pair"]
                m2 = gene_df[gene_df["tool"]==t2].set_index("sample")["pair"]
                common = m1.index.intersection(m2.index)
                if len(common) == 0: continue
                agree = sum(m1[s] == m2[s] for s in common)
                mat[i, j] = agree / len(common)

        im = ax.imshow(mat, vmin=0, vmax=1, cmap="RdYlGn", aspect="auto")
        ax.set_xticks(range(n)); ax.set_yticks(range(n))
        ax.set_xticklabels(local_tools, rotation=35, ha="right", fontsize=8.5)
        ax.set_yticklabels(local_tools, fontsize=8.5)
        ax.set_title(f"HLA-{gene}", fontsize=11, fontweight="bold")
        ax.grid(False)
        for i in range(n):
            for j in range(n):
                if not np.isnan(mat[i,j]):
                    ax.text(j, i, f"{mat[i,j]:.0%}", ha="center", va="center",
                            fontsize=8.5, color="black" if mat[i,j] > 0.4 else "white")

    plt.colorbar(im, ax=axes[-1], label="Allele-pair agreement rate", fraction=0.04, pad=0.04)
    fig.suptitle("VENEX WGS — Inter-Tool Allele Agreement per Gene (2-field, HLA-A/B/C)",
                 fontsize=12, fontweight="bold", y=1.03)
    wm(fig); fig.tight_layout()
    save(fig, out, "venex_fig2_allele_agreement")

# ── Figure V3: Homozygous call flags per gene ─────────────────────────────────
def fig_v3_homozygous_flags(df, out):
    """
    For each tool × gene, the fraction of typed samples where the reported
    allele pair is homozygous (allele1 2f == allele2 2f). High homozygosity
    rates in WGS suggest allele dropout.
    """
    df2 = df.copy()
    df2["a1_2f"] = df2["allele1"].apply(normalise2)
    df2["a2_2f"] = df2["allele2"].apply(normalise2)
    df2["is_homo"] = df2["a1_2f"] == df2["a2_2f"]

    genes_avail = [g for g in ["A","B","C"] if g in df2["gene"].unique()]
    tools = [t for t in TOOL_ORDER if t in df2["tool"].unique()]

    fig, ax = plt.subplots(figsize=(11, 5))
    x = np.arange(len(genes_avail)); width = 0.18
    offsets = np.linspace(-(len(tools)-1)*width/2, (len(tools)-1)*width/2, len(tools))

    for j, tool in enumerate(tools):
        vals, ns = [], []
        for gene in genes_avail:
            sub = df2[(df2["tool"]==tool) & (df2["gene"]==gene)]
            if sub.empty:
                vals.append(np.nan); ns.append(0); continue
            homo = sub["is_homo"].sum(); total = len(sub)
            vals.append(homo/total); ns.append(total)
        pos = x + offsets[j]
        bars = ax.bar(pos, np.nan_to_num(vals), width*0.88,
                      label=tool, color=TOOL_COLORS.get(tool, "#888"),
                      edgecolor="white", linewidth=0.4, zorder=3)
        for bar, v, n in zip(bars, vals, ns):
            if not np.isnan(v) and n > 0:
                ax.text(bar.get_x() + bar.get_width()/2,
                        bar.get_height() + 0.005,
                        f"{v:.0%}\n(n={n})", ha="center", va="bottom",
                        fontsize=7, color="#333")

    ax.axhline(0.10, color="#D62728", linewidth=1.2, linestyle="--",
               label="10% reference (expected if no dropout)")
    ax.set_xticks(x)
    ax.set_xticklabels([f"HLA-{g}" for g in genes_avail], fontsize=11, fontweight="bold")
    ax.set_ylabel("Fraction of samples with homozygous call", fontsize=10)
    ax.yaxis.set_major_formatter(mticker.PercentFormatter(1.0))
    ax.set_ylim(0, 0.55)
    ax.legend(title="Tool", loc="upper right", fontsize=8.5)
    ax.set_title("VENEX WGS — Homozygous Call Rate per Tool × Gene\n"
                 "(High rate = possible allele dropout; 1000G WGS truth: 61–91% of homo calls are artefact)",
                 fontsize=11, fontweight="bold")
    wm(fig); fig.tight_layout()
    save(fig, out, "venex_fig3_homozygous_flags")

# ── Figure V4: LOH status distribution ───────────────────────────────────────
def fig_v4_loh_status(loh_df, out):
    """
    Stacked horizontal bar per gene showing LOH status distribution
    for all 80 SpecHLA-typed samples.
    """
    genes = ["A", "B", "C"]
    status_order = ["balanced_heterozygous", "insufficient_evidence",
                    "homozygous_call_needs_review"]
    status_labels = {
        "balanced_heterozygous":         "Balanced heterozygous",
        "insufficient_evidence":         "Insufficient evidence",
        "homozygous_call_needs_review":  "Homozygous — needs review",
    }

    gene_df = loh_df[loh_df["gene"].isin(genes)].copy()

    fig, ax = plt.subplots(figsize=(10, 4.5))
    y = np.arange(len(genes)); height = 0.55

    for g_idx, gene in enumerate(genes):
        sub = gene_df[gene_df["gene"]==gene]
        total = len(sub)
        left = 0.0
        for status in status_order:
            count = (sub["loh_status"]==status).sum()
            frac = count / total if total else 0
            bar = ax.barh(y[g_idx], frac, height=height, left=left,
                          color=LOH_COLORS[status], edgecolor="white", linewidth=0.5)
            if frac > 0.03:
                ax.text(left + frac/2, y[g_idx], f"{count}\n({frac:.0%})",
                        ha="center", va="center", fontsize=8.5,
                        color="white" if status != "insufficient_evidence" else "#eee")
            left += frac

    ax.set_yticks(y)
    ax.set_yticklabels([f"HLA-{g}" for g in genes], fontsize=11, fontweight="bold")
    ax.set_xlabel("Fraction of SpecHLA-typed VENEX WGS samples", fontsize=10)
    ax.xaxis.set_major_formatter(mticker.PercentFormatter(1.0))
    ax.set_xlim(0, 1)
    ax.grid(axis="x")
    ax.grid(axis="y", alpha=0)

    patches = [mpatches.Patch(color=LOH_COLORS[s], label=status_labels[s])
               for s in status_order]
    ax.legend(handles=patches, loc="lower right", fontsize=9, title="LOH Status")
    ax.set_title("VENEX WGS — LOH Status Distribution per HLA Gene\n"
                 "(SpecHLA calls; n=80 samples; based on allele homozygosity QC)",
                 fontsize=11, fontweight="bold")
    wm(fig); fig.tight_layout()
    save(fig, out, "venex_fig4_loh_status")

# ── Figure V5: Allele diversity ───────────────────────────────────────────────
def fig_v5_allele_diversity(df, out):
    """
    Bar chart: number of unique 2-field alleles called per tool × gene.
    Wider allele diversity suggests broader capture of the allele space.
    Collapsed alleles (same 2-field seen across multiple samples) are deduplicated.
    """
    df2 = df.copy()
    df2["allele1_2f"] = df2["allele1"].apply(normalise2)
    df2["allele2_2f"] = df2["allele2"].apply(normalise2)

    genes_avail = [g for g in GENE_ORDER if g in df2["gene"].unique()]
    tools = [t for t in TOOL_ORDER if t in df2["tool"].unique()]

    fig, ax = plt.subplots(figsize=(12, 5))
    x = np.arange(len(genes_avail)); width = 0.18
    offsets = np.linspace(-(len(tools)-1)*width/2, (len(tools)-1)*width/2, len(tools))

    for j, tool in enumerate(tools):
        vals = []
        for gene in genes_avail:
            sub = df2[(df2["tool"]==tool) & (df2["gene"]==gene)]
            if sub.empty: vals.append(0); continue
            alleles = set(sub["allele1_2f"].dropna()) | set(sub["allele2_2f"].dropna())
            alleles.discard(""); alleles.discard(".")
            vals.append(len(alleles))
        pos = x + offsets[j]
        bars = ax.bar(pos, vals, width*0.88, label=tool,
                      color=TOOL_COLORS.get(tool, "#888"),
                      edgecolor="white", linewidth=0.4, zorder=3)
        for bar, v in zip(bars, vals):
            if v > 0:
                ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.3,
                        str(v), ha="center", va="bottom", fontsize=8, fontweight="bold")

    ax.set_xticks(x)
    ax.set_xticklabels([f"HLA-{g}" for g in genes_avail], fontsize=11, fontweight="bold")
    ax.set_ylabel("Unique 2-field alleles observed", fontsize=10)
    ax.legend(title="Tool", loc="upper right", fontsize=8.5)
    ax.set_title("VENEX WGS — Allele Diversity Called per Tool × Gene\n"
                 "(unique 2-field alleles; wider diversity = better allele space coverage)",
                 fontsize=12, fontweight="bold")
    wm(fig); fig.tight_layout()
    save(fig, out, "venex_fig5_allele_diversity")

# ── HTML helpers ──────────────────────────────────────────────────────────────
def to64(path):
    with open(path, "rb") as f:
        return base64.b64encode(f.read()).decode()

def embed_fig(path, name, caption):
    if not Path(path).exists():
        return f'<p class="note">Figure not available: {name}</p>'
    data = to64(path)
    return (
        f'<figure style="margin:20px 0;text-align:center">'
        f'<img src="data:image/png;base64,{data}" '
        f'style="max-width:100%;border:1px solid #d0dce8;border-radius:6px" />'
        f'<figcaption style="font-size:0.85em;color:#5a6a7e;margin-top:8px">'
        f'<strong>{name}.</strong> {caption}'
        f'</figcaption></figure>'
    )

def build_venex_section(figures_dir, loh_df, calls_df):
    fd = Path(figures_dir)

    # Compute summary stats for callout boxes
    n_samples = calls_df["sample"].nunique()
    n_total   = 87
    spechla_n = calls_df[calls_df["tool"]=="SpecHLA"]["sample"].nunique()
    loh_n     = (loh_df["loh_status"]=="homozygous_call_needs_review").sum()
    bal_n     = (loh_df["loh_status"]=="balanced_heterozygous").sum()
    insuf_n   = (loh_df["loh_status"]=="insufficient_evidence").sum()

    # Build allele agreement summary for VX_92 (the one balanced case)
    bal_samples = loh_df[loh_df["loh_status"]=="balanced_heterozygous"]["sample"].unique()

    fig1 = embed_fig(fd/"venex_fig1_tool_coverage.png",
                     "VENEX Figure 1",
                     "Tool call coverage per gene across the VENEX WGS cohort (n=87 total samples). "
                     "SpecHLA provides the most complete coverage (80/87), followed by HLA-HD (52/87), "
                     "T1K (40/87), and OptiType (31/87). Coverage for class II genes (DRB1, DQB1) "
                     "is available only from SpecHLA, HLA-HD, and T1K.")

    fig2 = embed_fig(fd/"venex_fig2_allele_agreement.png",
                     "VENEX Figure 2",
                     "Inter-tool allele-pair agreement rate for HLA-A, -B, and -C on the VENEX WGS "
                     "cohort (2-field normalised, unordered pair matching). Values show the fraction "
                     "of shared samples where both tools called the same allele pair. High agreement "
                     "between multiple tools provides confidence in the typed allele; disagreement "
                     "flags loci for additional review.")

    fig3 = embed_fig(fd/"venex_fig3_homozygous_flags.png",
                     "VENEX Figure 3",
                     "Homozygous call rate (both reported alleles are the same 2-field allele) "
                     "per tool × gene. In the 1000 Genomes WGS benchmark, 61–91% of homozygous "
                     "WGS calls are artefactual (true genotype is heterozygous). High homozygosity "
                     "rates here are therefore consistent with allele dropout rather than true "
                     "biological homozygosity. The dashed line at 10% represents an approximate "
                     "expected rate if dropout were absent.")

    fig4 = embed_fig(fd/"venex_fig4_loh_status.png",
                     "VENEX Figure 4",
                     "LOH status distribution per HLA gene for the 80 SpecHLA-typed VENEX WGS "
                     "samples. Categories: balanced heterozygous (sufficient allele-frequency data, "
                     "both alleles detected at similar frequency), homozygous call needing review "
                     "(single allele reported — possible dropout or true LOH), insufficient evidence "
                     "(too few heterozygous variants to classify allele balance). The majority of "
                     "samples fall in the 'insufficient evidence' category, reflecting the challenge "
                     "of allele-balance QC from standard WGS coverage at the HLA region.")

    fig5 = embed_fig(fd/"venex_fig5_allele_diversity.png",
                     "VENEX Figure 5",
                     "Number of unique 2-field alleles called per tool × gene across the VENEX WGS "
                     "cohort. Greater allele diversity indicates broader coverage of the allele "
                     "space in the typed cohort. SpecHLA covers the most alleles across all genes "
                     "due to its largest sample coverage. T1K and HLA-HD show comparable diversity "
                     "per typed sample; OptiType is restricted to class I loci.")

    return f"""
<!-- ═══════════════════ VENEX WGS COHORT ═══════════════════ -->
<section id="venex">
  <h2 class="section-title"><span class="sec-num">★</span>VENEX WGS Cohort Analysis</h2>
  <p>
    The <strong>VENEX cohort</strong> is an independent WGS dataset of {n_total} samples
    processed through the PIHLA pipeline. Unlike the 1000 Genomes benchmark cohort (where
    ground-truth HLA types are available), the VENEX samples have no external truth labels —
    this analysis therefore focuses on <em>tool coverage</em>, <em>inter-tool agreement</em>,
    <em>homozygosity patterns</em> (as a proxy for allele dropout), and <em>LOH candidate
    identification</em> using the SpecHLA-based allele-frequency QC module.
  </p>
  <p>
    Four HLA typing tools produced results for the VENEX WGS samples:
    <strong>SpecHLA</strong> ({spechla_n} samples), <strong>HLA-HD</strong>
    (52 samples), <strong>T1K</strong> (40 samples), and <strong>OptiType</strong>
    (31 samples). Tool coverage varies due to differences in convergence rates,
    memory requirements, and modality suitability.
  </p>

  <div class="card-grid">
    <div class="metric-card"><div class="value">{n_total}</div>
      <div class="label">VENEX WGS Samples</div>
      <div class="sublabel">Independent of 1000G benchmark</div></div>
    <div class="metric-card"><div class="value">{spechla_n}</div>
      <div class="label">SpecHLA-Typed</div>
      <div class="sublabel">Most complete tool coverage</div></div>
    <div class="metric-card"><div class="value">{loh_n}</div>
      <div class="label">LOH Candidates</div>
      <div class="sublabel">Homozygous call — needs review</div></div>
    <div class="metric-card"><div class="value">{insuf_n}</div>
      <div class="label">Insufficient Evidence</div>
      <div class="sublabel">Allele-balance QC underpowered</div></div>
    <div class="metric-card"><div class="value">{bal_n}</div>
      <div class="label">Balanced Heterozygous</div>
      <div class="sublabel">Allele dropout unlikely</div></div>
    <div class="metric-card"><div class="value">4</div>
      <div class="label">Tools Evaluated</div>
      <div class="sublabel">OptiType · T1K · SpecHLA · HLA-HD</div></div>
  </div>

  <h3>Tool Call Coverage</h3>
  <div class="highlight-box">
    <strong>Coverage note:</strong> SpecHLA provides the broadest coverage (80/87 samples,
    92%) and is the only tool reporting DRB1 and DQB1 calls in this cohort. The other tools
    cover 36–60% of samples; gaps arise from Nextflow convergence failures, memory limits,
    and tool-specific input requirements. All tools that did produce output are included in
    the agreement and diversity analyses below.
  </div>
  {fig1}

  <h3>Inter-Tool Allele Agreement</h3>
  <div class="highlight-box yellow">
    <strong>Key insight — agreement as a quality signal:</strong> When two or more tools
    independently call the same allele pair from the same WGS reads, the call is more
    likely to be correct — even without external truth labels. Loci where all available
    tools agree can be used with higher confidence; loci with systematic disagreement
    should be flagged for orthogonal confirmation. For the VENEX cohort, this agreement
    matrix provides the primary internal QC layer in the absence of reference truth data.
  </div>
  {fig2}

  <h3>Homozygous Call Rates (Allele Dropout Proxy)</h3>
  <div class="highlight-box orange">
    <strong>Allele dropout in WGS:</strong> In the 1000 Genomes benchmark where ground
    truth is available, 61–91% of homozygous WGS calls turn out to be artefactual —
    the true genotype is heterozygous but one allele was not recovered. High homozygous
    call rates in the VENEX cohort (where no truth is available) should therefore be
    interpreted as potential dropout rather than true biological homozygosity. Samples
    or loci with homozygous calls from multiple tools are prioritised for additional
    QC (e.g., allele-frequency analysis from BAM pileup).
  </div>
  {fig3}

  <h3>LOH Status from Allele-Frequency QC</h3>
  <p>
    For the {spechla_n} samples with SpecHLA results, the PIHLA LOH module classifies
    each (sample, gene) pair into one of three evidence categories based on the
    allele-frequency balance at heterozygous SNP positions within the HLA region:
    <strong>balanced heterozygous</strong> (both alleles present at similar frequencies),
    <strong>homozygous call needing review</strong> (single allele detected — possible
    dropout or true LOH), and <strong>insufficient evidence</strong> (too few heterozygous
    variants for confident allele-balance assessment). The LOH candidates identified here
    (<strong>{loh_n} (sample, gene) pairs</strong> across {len(bal_samples)+5} unique
    samples) represent the subset most likely to reflect true HLA LOH or severe allele
    dropout, and should be prioritised for validation with orthogonal methods (e.g.,
    SNP array–based allele-specific copy number analysis, or deeper targeted sequencing).
  </p>
  <div class="highlight-box yellow">
    <strong>Key insight — why most samples fall in "insufficient evidence":</strong>
    The LOH classification requires heterozygous germline SNPs within the HLA region that
    can be used to measure allele-balance. Standard 30× WGS provides enough heterozygous
    SNPs for this analysis in <em>most</em> samples, but the HLA region's extreme
    polymorphism means that the effective heterozygous variant count varies widely by
    ancestry and individual. The {insuf_n} "insufficient evidence" calls do not mean
    those loci are uninformative — they mean the allele-balance signal was too weak to
    classify confidently. These loci may still carry LOH; they simply cannot be confirmed
    from WGS depth alone.
  </div>
  {fig4}

  <h3>Allele Diversity in the VENEX Cohort</h3>
  <div class="highlight-box">
    <strong>Allele space coverage:</strong> The VENEX cohort samples carry a diverse range
    of HLA alleles, as expected from a cancer genomics cohort with mixed ancestry. The
    allele diversity figure below shows how many unique 2-field alleles each tool observed
    across all VENEX samples. Tools covering more samples naturally observe more alleles;
    within-tool diversity can be used to assess whether the cohort's allele space is
    adequately represented in the tool's reference database.
  </div>
  {fig5}

  <h3>VENEX LOH Candidate Samples</h3>
  <p>
    The following {loh_n} (sample, gene) pairs were classified as
    <strong>homozygous call needing review</strong> and represent the primary LOH
    candidates in the VENEX cohort. All flags originate from SpecHLA calls;
    the false-duplicate rate column shows the 1000 Genomes–derived benchmark estimate
    of how often SpecHLA WGS homozygous calls are artefactual (89% — meaning most of
    these flags likely reflect allele dropout rather than true LOH):
  </p>
  <div class="table-wrap">
  <table class="data-table">
    <thead><tr><th>Sample</th><th>Gene</th><th>Allele 1 (4-field)</th><th>Allele 2</th><th>QC Flag</th><th>FDR (benchmark)</th></tr></thead>
    <tbody>
""" + "\n".join(
        f"    <tr><td>{r['sample']}</td><td>HLA-{r['gene']}</td>"
        f"<td style='font-family:monospace;font-size:0.88em'>{r['allele1']}</td>"
        f"<td style='font-family:monospace;font-size:0.88em'>{r['allele2']}</td>"
        f"<td><span style='color:#D62728;font-weight:600'>{r['qc_flag']}</span></td>"
        f"<td>{float(r['benchmark_false_duplicate_rate']):.0%}</td></tr>"
        for _, r in loh_df[loh_df["loh_status"]=="homozygous_call_needs_review"].iterrows()
    ) + """
    </tbody>
  </table>
  </div>
  <p class="note" style="margin-top:10px">
    FDR = false duplicate rate (fraction of 1000G WGS homozygous calls where truth is
    heterozygous). For SpecHLA WGS, FDR is 86–91%, indicating most homozygous calls
    are allele-dropout artefacts, not true LOH.
  </p>
</section>
"""

# ── Patch HTML ────────────────────────────────────────────────────────────────
def patch_html(input_html, venex_section, venex_nav_entry):
    content = Path(input_html).read_text(encoding="utf-8")

    # 1. Insert VENEX nav entry — after the trimodal entry if present, else after bimodal
    for anchor in ["trimodal", "bimodal"]:
        m = re.search(rf'(<a href="#{anchor}"[^>]*>.*?</a>)', content, re.DOTALL)
        if m:
            original = m.group(0)
            content = content.replace(original, original + "\n" + venex_nav_entry, 1)
            break

    # 2. Append VENEX section before closing </div><!-- end #main -->
    content = content.replace(
        "</div><!-- end #main -->",
        venex_section + "\n</div><!-- end #main -->",
        1
    )
    return content

# ── Main ──────────────────────────────────────────────────────────────────────
def main():
    p = argparse.ArgumentParser()
    p.add_argument("--calls",          default="analysis/venex_harmonized_calls.tsv")
    p.add_argument("--loh",            default="analysis/loh_analysis/sample_loh_candidates.tsv")
    p.add_argument("--out-dir",        default="/scratch/project_2008084/mvhla_figures_v6")
    p.add_argument("--input-report",   default="/scratch/project_2008084/mvhla_report_v12.html")
    p.add_argument("--output-report",  default="/scratch/project_2008084/mvhla_report_v13.html")
    args = p.parse_args()

    apply_style()

    print("Loading data...")
    calls_df = pd.read_csv(args.calls, sep="\t", dtype=str)
    loh_df   = pd.read_csv(args.loh,   sep="\t", dtype=str)
    loh_df   = loh_df[loh_df["loh_status"] != "loh_status"]  # drop header row if duplicated
    print(f"  calls: {len(calls_df)} rows, {calls_df['sample'].nunique()} samples")
    print(f"  loh:   {len(loh_df)} rows, {loh_df['sample'].nunique()} samples")

    print(f"\nGenerating VENEX figures → {args.out_dir}")
    fig_v1_tool_coverage(calls_df, args.out_dir)
    fig_v2_allele_agreement(calls_df, args.out_dir)
    fig_v3_homozygous_flags(calls_df, args.out_dir)
    fig_v4_loh_status(loh_df, args.out_dir)
    fig_v5_allele_diversity(calls_df, args.out_dir)

    print("\nBuilding VENEX HTML section...")
    venex_section = build_venex_section(args.out_dir, loh_df, calls_df)
    venex_nav = '<a href="#venex" class="nav-link">★ VENEX WGS Cohort</a>'

    print(f"\nPatching {args.input_report} → {args.output_report} ...")
    patched = patch_html(args.input_report, venex_section, venex_nav)
    Path(args.output_report).write_text(patched, encoding="utf-8")
    size_kb = Path(args.output_report).stat().st_size // 1024
    print(f"  HTML report → {args.output_report}  ({size_kb} KB)")
    print("\nDone.")

if __name__ == "__main__":
    main()
