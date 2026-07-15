#!/usr/bin/env python3
"""Generate ChampHLA two-week progress report (Jun 15-29, 2026).

Self-contained HTML with embedded base64 PNG figures, sidebar navigation,
and CSS-variable styling following the project's existing report pattern.

Usage:
    module load gcc/11.3.0 biopythontools/11.3.0_3.10.6
    python3 bin/generate_progress_report_jun2026.py
"""

import base64
import html
from datetime import datetime
from pathlib import Path

import pandas as pd

BASE = Path("/scratch/project_2008084/pihla-publish")
OUT = BASE / "analysis" / "champhla_progress_report_jun2026.html"
TIMING_TSV = BASE / "analysis" / "performance_benchmark_n30" / "tables" / "timing_summary_n30.tsv"
FIGS_FINAL = BASE / "analysis" / "figures_final"
FIGS_CANDIDATE = BASE / "analysis" / "figures_final_candidate"

FIGURES = [
    {
        "path": FIGS_FINAL / "figure_8_computational_performance.png",
        "label": "Figure 8",
        "title": "Computational resource requirements",
        "caption": "Three-panel figure showing peak RAM (GB), CPU utilisation (%), and wall-clock time (hours) per tool across WGS, WES, and RNA-seq. n=30 per modality; RNA tool-level n varies (T1K/Seq2HLA n=30, ArcasHLA/HLA-HD/SpecHLA n=29, OptiType n=26).",
        "date": "Jun 19, 2026",
        "tag": "Updated",
    },
    {
        "path": FIGS_FINAL / "figure_11_champion_challenger_combined.png",
        "label": "Figure 11",
        "title": "Champion-Challenger ensemble mechanism",
        "caption": "Three-panel composite: (A) per-gene champion assignment and override gate logic, (B) MV vs CC accuracy comparison, (C) confidence-derived weights.",
        "date": "Jun 24, 2026",
        "tag": "New",
    },
    {
        "path": FIGS_FINAL / "figure_11a_cc_mechanism.png",
        "label": "Figure 11a",
        "title": "Champion-Challenger mechanism detail",
        "caption": "Per-gene champion assignment and override gate logic showing how the frozen 1000G policy routes calls.",
        "date": "Jun 24, 2026",
        "tag": "New",
    },
    {
        "path": FIGS_FINAL / "figure_11b_mv_vs_cc.png",
        "label": "Figure 11b",
        "title": "Majority Voting vs Champion-Challenger",
        "caption": "Direct comparison of majority voting vs Champion-Challenger accuracy across modalities.",
        "date": "Jun 24, 2026",
        "tag": "New",
    },
    {
        "path": FIGS_FINAL / "figure_12_silver_truth_hprc.png",
        "label": "Figure 12",
        "title": "Silver-standard truth validation (HPRC)",
        "caption": "Locityper + HPRC pangenome + Immuannot silver-truth concordance vs 1000G gold truth. 94.4% allele-level concordance (two-field) at per-sample depth; 81.1% two-field / 82.2% G-group across 30 samples.",
        "date": "Jun 26, 2026",
        "tag": "New",
    },
    {
        "path": FIGS_CANDIDATE / "figure_13_external_validation_nci60.png",
        "label": "Figure 13",
        "title": "External validation on NCI-60",
        "caption": "Champion-Challenger performance on the independent NCI-60 cancer cell-line panel (RNA-seq, n=4 pilot). Frozen 1000G policy generalises with 2 corrective overrides, 0 harmful.",
        "date": "Jun 28, 2026",
        "tag": "Candidate",
    },
    {
        "path": FIGS_FINAL / "figure_s4_fimm_concordance.png",
        "label": "Figure S4",
        "title": "FIMM cross-modality HLA concordance",
        "caption": "Concordance rates for HLA-A, HLA-B, HLA-C across scRNA-seq, bulk RNA-seq, and WES in the FIMM AML/MDS cohort.",
        "date": "Jun 18, 2026",
        "tag": "New",
    },
    {
        "path": FIGS_FINAL / "figure_s5_allele_dropout_heatmap.png",
        "label": "Figure S5",
        "title": "HLA allele dropout characterisation",
        "caption": "Two-panel: (A) false duplicate rate heatmap across tools/modalities/genes; (B) allele dropout counts by modality from 1000G truth-backed benchmarks.",
        "date": "Jun 18, 2026",
        "tag": "New",
    },
    {
        "path": FIGS_FINAL / "figure_s6_fimm_loh_wes.png",
        "label": "Figure S6",
        "title": "FIMM WES LOH status distribution",
        "caption": "Stacked bar of LOH classification per HLA gene from SpecHLA WES allele frequencies. Categories: candidate LOH, allelic imbalance, balanced heterozygous, insufficient evidence.",
        "date": "Jun 18, 2026",
        "tag": "New",
    },
    {
        "path": FIGS_FINAL / "figure_s7_fimm_survival_hed.png",
        "label": "Figure S7",
        "title": "FIMM overall survival and HLA evolutionary divergence",
        "caption": "Two-panel: (A) KM survival by diagnosis group (AML, MDS, MDS→AML); (B) KM by HED group (median split). n≈28, exploratory.",
        "date": "Jun 18, 2026",
        "tag": "New",
    },
]


def safe(text):
    return html.escape(str(text), quote=True)


def img_file_b64(path):
    with open(path, "rb") as f:
        return base64.b64encode(f.read()).decode()


def render_table(headers, rows, caption=None):
    head = "".join(f"<th>{safe(h)}</th>" for h in headers)
    body = ""
    for row in rows:
        body += "<tr>" + "".join(f"<td>{cell}</td>" for cell in row) + "</tr>"
    cap = f"<caption>{safe(caption)}</caption>" if caption else ""
    return (
        f"<div class='table-wrap'><table>{cap}"
        f"<thead><tr>{head}</tr></thead>"
        f"<tbody>{body}</tbody></table></div>"
    )


def build_timing_table(tsv_path):
    df = pd.read_csv(tsv_path, sep="\t")
    headers = ["Tool", "Modality", "n", "Wall-clock (h)", "Peak RAM (GB)", "CPU (%)"]
    rows = []
    for _, r in df.iterrows():
        mod = r["modality"]
        mod_lower = safe(mod.lower())
        mod_display = safe(mod)
        rows.append([
            safe(r["tool"]),
            f"<span class='modality-badge {mod_lower}'>{mod_display}</span>",
            safe(str(int(r["n_samples"]))),
            f"{r['median_runtime_hours']:.3f}",
            f"{r['median_peak_ram_gb']:.2f}",
            f"{r['median_cpu_pct']:.0f}",
        ])
    return render_table(headers, rows, "Median computational resource usage per tool and modality (n=30 benchmark)")


def build_figure_gallery(figures):
    cards = []
    for fig in figures:
        if not fig["path"].exists():
            continue
        b64 = img_file_b64(fig["path"])
        tag_class = fig["tag"].lower()
        cards.append(f"""
        <article class="figure-card panel">
          <div class="figure-header">
            <div>
              <p class="eyebrow">{safe(fig['label'])}</p>
              <h3>{safe(fig['title'])}</h3>
            </div>
            <div style="display:flex;gap:8px;align-items:center;">
              <span class="fig-tag {tag_class}">{safe(fig['tag'])}</span>
              <span class="figure-date">{safe(fig['date'])}</span>
            </div>
          </div>
          <div class="figure-box">
            <img src="data:image/png;base64,{b64}" alt="{safe(fig['title'])}" />
          </div>
          <div class="figure-copy">
            <p>{safe(fig['caption'])}</p>
          </div>
        </article>
        """)
    return "\n".join(cards)


SECTIONS = [
    "Executive Summary",
    "New Figures",
    "External Validation",
    "Pipeline Changes",
    "Computational Performance",
    "Documentation",
    "Manuscript Updates",
    "Next Steps",
]


def section_id(title):
    return title.lower().replace(" ", "-").replace("&", "and")


def build_html():
    generated = datetime.now().strftime("%Y-%m-%d %H:%M")

    nav_html = "".join(
        f"<a href='#{section_id(s)}'>{safe(s)}</a>" for s in SECTIONS
    )

    figure_gallery = build_figure_gallery(FIGURES)
    timing_table = build_timing_table(TIMING_TSV)
    n_figures = sum(1 for f in FIGURES if f["path"].exists())

    return f"""<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8" />
  <meta name="viewport" content="width=device-width, initial-scale=1.0" />
  <title>ChampHLA Progress Report — June 2026</title>
  <style>
    :root {{
      --bg: #f3f6fb;
      --surface: #ffffff;
      --ink: #102033;
      --muted: #5b6472;
      --line: #d7e1ee;
      --navy: #0f2744;
      --blue: #275dad;
      --teal: #0f766e;
      --amber: #b45309;
      --green: #065f46;
      --rose: #b91c1c;
      --shadow: 0 20px 44px rgba(15, 39, 68, 0.10);
      --radius: 22px;
      --sidebar-width: 300px;
    }}
    * {{ box-sizing: border-box; }}
    html {{ scroll-behavior: smooth; }}
    body {{
      margin: 0;
      font-family: Arial, Helvetica, sans-serif;
      color: var(--ink);
      background: linear-gradient(180deg, #eef3f9 0%, #f8fbfd 100%);
    }}
    .layout {{
      display: grid;
      grid-template-columns: var(--sidebar-width) minmax(0, 1fr);
      min-height: 100vh;
    }}
    aside {{
      position: sticky; top: 0; align-self: start; height: 100vh;
      padding: 28px 22px;
      background: linear-gradient(180deg, #102844 0%, #174576 55%, #1f5ea8 100%);
      color: #fff;
      border-right: 1px solid rgba(255,255,255,0.10);
      overflow-y: auto;
    }}
    .brand {{ margin-bottom: 22px; padding-bottom: 18px; border-bottom: 1px solid rgba(255,255,255,0.14); }}
    .brand .kicker {{
      margin: 0 0 10px; color: #c8dbf2; font-size: 12px;
      letter-spacing: 0.08em; text-transform: uppercase; font-weight: 700;
    }}
    .brand h1 {{ margin: 0 0 10px; font-size: 24px; line-height: 1.15; }}
    .brand p {{ margin: 0; color: #d7e6f7; font-size: 13px; line-height: 1.5; }}
    .sidebar-meta {{ display: grid; gap: 10px; margin: 18px 0 24px; }}
    .sidebar-chip {{
      padding: 10px 12px; border: 1px solid rgba(255,255,255,0.14);
      border-radius: 16px; background: rgba(255,255,255,0.07);
      font-size: 13px; color: #e7f0fa;
    }}
    .sidebar-chip strong {{ display: block; margin-bottom: 4px; color: #fff; font-size: 11px; text-transform: uppercase; letter-spacing: 0.06em; }}
    .nav {{ display: grid; gap: 6px; }}
    .nav a {{
      color: #d9e8f7; text-decoration: none; padding: 9px 12px;
      border-radius: 14px; border: 1px solid transparent;
      font-size: 13px; line-height: 1.35; transition: all 0.15s;
    }}
    .nav a:hover, .nav a.active {{
      background: rgba(255,255,255,0.10); border-color: rgba(255,255,255,0.16); color: #fff;
    }}
    main {{ padding: 30px 32px 48px; }}
    .hero {{
      background: linear-gradient(135deg, #102844 0%, #174576 55%, #275dad 100%);
      color: #fff; border-radius: 28px; box-shadow: var(--shadow);
      padding: 32px 34px; margin-bottom: 24px;
    }}
    .hero .eyebrow, .eyebrow {{
      margin: 0 0 12px; font-size: 12px; letter-spacing: 0.08em;
      text-transform: uppercase; font-weight: 700; color: #c9dcf4;
    }}
    .hero h2 {{ margin: 0 0 10px; font-size: 34px; line-height: 1.08; }}
    .hero p {{ margin: 0; max-width: 1040px; font-size: 16px; line-height: 1.6; color: #e3eef9; }}
    section {{ margin-top: 28px; }}
    .section-header {{ margin-bottom: 16px; }}
    .section-header h2 {{ margin: 0 0 8px; font-size: 26px; color: var(--navy); }}
    .section-header p {{ margin: 0; color: var(--muted); font-size: 15px; line-height: 1.6; max-width: 960px; }}
    .grid {{ display: grid; gap: 18px; }}
    .metrics-grid {{ grid-template-columns: repeat(auto-fit, minmax(180px, 1fr)); }}
    .two-col {{ grid-template-columns: 1fr 1fr; align-items: start; }}
    .panel, .metric {{
      background: linear-gradient(180deg, #fbfdff 0%, #ffffff 100%);
      border: 1px solid var(--line); border-radius: var(--radius);
      box-shadow: 0 8px 24px rgba(15, 39, 68, 0.05);
    }}
    .panel {{ padding: 22px; }}
    .panel h3 {{ margin: 0 0 12px; font-size: 18px; color: var(--navy); }}
    .panel p, .panel li {{ margin: 0; font-size: 14px; line-height: 1.65; color: var(--ink); }}
    .panel ul {{ margin: 8px 0 0; padding-left: 18px; display: grid; gap: 8px; }}
    .metric {{ padding: 18px; text-align: center; }}
    .metric .label {{ font-size: 11px; text-transform: uppercase; letter-spacing: 0.06em; color: var(--muted); font-weight: 700; margin-bottom: 6px; }}
    .metric .value {{ font-size: 28px; font-weight: 700; color: var(--navy); line-height: 1.1; margin-bottom: 6px; }}
    .metric .note {{ font-size: 12px; line-height: 1.4; color: var(--muted); }}
    .callout {{
      padding: 18px 20px; border-left: 5px solid var(--blue);
      background: #f5f9ff; border-radius: 16px; color: var(--ink);
    }}
    .callout.success {{ border-left-color: var(--teal); background: #f0fdf4; }}
    .callout.warn {{ border-left-color: var(--amber); background: #fffaf2; }}
    .table-wrap {{ overflow-x: auto; border: 1px solid var(--line); border-radius: 18px; background: #fff; }}
    table {{ width: 100%; border-collapse: collapse; min-width: 600px; }}
    caption {{ text-align: left; padding: 16px 18px 0; font-size: 13px; color: var(--muted); }}
    thead th {{ position: sticky; top: 0; background: #eef5fd; color: var(--navy); font-size: 11px; text-transform: uppercase; letter-spacing: 0.05em; }}
    th, td {{ padding: 10px 14px; border-bottom: 1px solid #e7eef6; text-align: left; font-size: 13px; vertical-align: top; }}
    tbody tr:nth-child(even) {{ background: #fbfdff; }}
    tbody tr:hover {{ background: #f0f5fc; }}
    .modality-badge {{
      display: inline-block; padding: 3px 8px; border-radius: 999px;
      font-size: 11px; font-weight: 700; text-transform: uppercase; letter-spacing: 0.04em;
    }}
    .modality-badge.wgs {{ color: #1e40af; background: #dbeafe; }}
    .modality-badge.wes {{ color: #065f46; background: #d1fae5; }}
    .modality-badge.rna {{ color: #9a3412; background: #ffedd5; }}
    .figure-card {{ padding: 22px; margin-bottom: 4px; }}
    .figure-header {{ display: flex; justify-content: space-between; gap: 16px; align-items: flex-start; margin-bottom: 14px; flex-wrap: wrap; }}
    .figure-header h3 {{ margin: 0; font-size: 18px; }}
    .figure-date {{ font-size: 12px; color: var(--muted); background: #f2f6fb; border: 1px solid var(--line); padding: 6px 10px; border-radius: 999px; white-space: nowrap; }}
    .figure-box {{
      border: 1px solid var(--line); border-radius: 20px; background: #f9fbfe;
      padding: 12px; display: flex; justify-content: center; align-items: center;
      min-height: 200px; margin-bottom: 14px;
    }}
    .figure-box img {{ display: block; max-width: 100%; width: 100%; height: auto; border-radius: 14px; }}
    .figure-copy p {{ font-size: 13px; color: var(--muted); line-height: 1.55; }}
    .fig-tag {{
      display: inline-block; padding: 4px 10px; border-radius: 999px;
      font-size: 11px; font-weight: 700; text-transform: uppercase; letter-spacing: 0.04em;
    }}
    .fig-tag.new {{ color: #065f46; background: #dff7ef; border: 1px solid #b6e7d7; }}
    .fig-tag.updated {{ color: #1e40af; background: #dbeafe; border: 1px solid #bfdbfe; }}
    .fig-tag.candidate {{ color: #9a3412; background: #fff1e8; border: 1px solid #f7d4c1; }}
    .changelog-card {{ margin-bottom: 12px; }}
    .changelog-card .ch-date {{ font-size: 11px; color: var(--muted); text-transform: uppercase; letter-spacing: 0.04em; }}
    .changelog-card .ch-title {{ font-size: 16px; font-weight: 700; color: var(--navy); margin: 4px 0 6px; }}
    .changelog-card .ch-body {{ font-size: 13px; color: var(--ink); line-height: 1.6; }}
    .changelog-card code {{ background: #eef3fa; padding: 2px 6px; border-radius: 6px; font-size: 12px; }}
    .result-table {{ width: 100%; border-collapse: collapse; min-width: 500px; }}
    .result-table th {{ background: #eef5fd; color: var(--navy); font-size: 12px; text-transform: uppercase; padding: 10px 12px; text-align: left; }}
    .result-table td {{ padding: 10px 12px; border-bottom: 1px solid #e7eef6; font-size: 13px; }}
    .result-table tr:hover {{ background: #f5f9ff; }}
    .highlight {{ font-weight: 700; color: var(--navy); }}
    footer {{
      margin-top: 40px; padding: 24px; text-align: center;
      font-size: 12px; color: var(--muted); border-top: 1px solid var(--line);
    }}
    @media (max-width: 1180px) {{
      .layout {{ grid-template-columns: 1fr; }}
      aside {{ position: relative; height: auto; }}
      .metrics-grid {{ grid-template-columns: repeat(2, 1fr); }}
      .two-col {{ grid-template-columns: 1fr; }}
    }}
    @media (max-width: 720px) {{
      main {{ padding: 16px; }}
      .hero h2 {{ font-size: 26px; }}
      .metrics-grid {{ grid-template-columns: 1fr; }}
    }}
  </style>
</head>
<body>
<div class="layout">
  <aside>
    <div class="brand">
      <p class="kicker">Progress Report</p>
      <h1>ChampHLA</h1>
      <p>Two-week development summary<br/>June 15 &ndash; 29, 2026</p>
    </div>
    <div class="sidebar-meta">
      <div class="sidebar-chip"><strong>Period</strong>Jun 15 &ndash; 29, 2026</div>
      <div class="sidebar-chip"><strong>Figures</strong>{n_figures} new/updated</div>
      <div class="sidebar-chip"><strong>Generated</strong>{safe(generated)}</div>
    </div>
    <nav class="nav" id="sidebar-nav">{nav_html}</nav>
  </aside>

  <main>
    <!-- HERO -->
    <section class="hero">
      <p class="eyebrow">ChampHLA &middot; Majority Voting HLA Typing Pipeline</p>
      <h2>Two-Week Progress Report</h2>
      <p>Summary of all new results, figures, pipeline changes, and documentation produced between June 15 and June 29, 2026. This was one of the most productive periods in the project &mdash; spanning external validation on two independent cohorts, ten new publication figures, deployment documentation, and manuscript finalisation.</p>
    </section>

    <!-- EXECUTIVE SUMMARY -->
    <section id="{section_id('Executive Summary')}">
      <div class="section-header">
        <h2>Executive Summary</h2>
        <p>Key deliverables completed in the last two weeks.</p>
      </div>
      <div class="grid metrics-grid">
        <div class="metric"><div class="label">New Figures</div><div class="value">{n_figures}</div><div class="note">Publication-quality (Figs 8, 11, 12, 13, S4&ndash;S7)</div></div>
        <div class="metric"><div class="label">External Cohorts</div><div class="value">2</div><div class="note">NCI-60 (RNA) + GIAB HG002 (WGS)</div></div>
        <div class="metric"><div class="label">Bug Fixes</div><div class="value">2</div><div class="note">Seq2HLA routing + Fig 8 layout</div></div>
        <div class="metric"><div class="label">New Docs</div><div class="value">4</div><div class="note">INSTALL, USAGE, README, Scouting</div></div>
        <div class="metric"><div class="label">Manuscript</div><div class="value">V3</div><div class="note">110.7 KB &mdash; all sections updated</div></div>
        <div class="metric"><div class="label">New Modules</div><div class="value">3</div><div class="note">Locityper, Immuannot, Silver-truth</div></div>
      </div>
    </section>

    <!-- NEW FIGURES -->
    <section id="{section_id('New Figures')}">
      <div class="section-header">
        <h2>New Figures</h2>
        <p>All figures created or updated during this period. Each figure is embedded directly in this report for offline viewing.</p>
      </div>
      <div class="grid">
        {figure_gallery}
      </div>
    </section>

    <!-- EXTERNAL VALIDATION -->
    <section id="{section_id('External Validation')}">
      <div class="section-header">
        <h2>External Validation Results</h2>
        <p>Two independent, non-1000G cohorts were used to test whether the frozen Champion-Challenger policy generalises beyond the training population.</p>
      </div>

      <div class="grid two-col">
        <!-- NCI-60 -->
        <div class="panel">
          <h3>NCI-60 Cancer Cell Lines (RNA-seq)</h3>
          <div class="callout success" style="margin-bottom:14px;">
            <strong>Key finding:</strong> The frozen 1000G-derived RNA Champion-Challenger policy fired 2 overrides on this external cohort &mdash; both corrective, 0 harmful &mdash; lifting champion-only accuracy from 0.75 to 0.917.
          </div>
          <div class="table-wrap">
            <table class="result-table">
              <thead><tr><th>Method</th><th>Accuracy</th><th>Score</th></tr></thead>
              <tbody>
                <tr><td class="highlight">ChampHLA (CC)</td><td>11/12</td><td class="highlight">0.917</td></tr>
                <tr><td>MajorityVote</td><td>11/12</td><td>0.917</td></tr>
                <tr><td>OptiType</td><td>11/12</td><td>0.917</td></tr>
                <tr><td>HLA-HD</td><td>11/12</td><td>0.917</td></tr>
                <tr><td>T1K</td><td>8/12</td><td>0.667</td></tr>
                <tr><td>ArcasHLA</td><td>7/12</td><td>0.583</td></tr>
              </tbody>
            </table>
          </div>
          <p style="margin-top:10px;font-size:12px;color:var(--muted);">
            n=4 pilot (SK-MEL-28, RPMI-8226, NCI-H23, OVCAR-8). Truth: Adams et al. 2005 SBT class I. 95% CI: 0.65&ndash;0.99.
          </p>
        </div>

        <!-- HPRC -->
        <div class="panel">
          <h3>GIAB/HPRC HG002 (WGS)</h3>
          <div class="callout success" style="margin-bottom:14px;">
            <strong>Key finding:</strong> Perfect three-field recovery on the field-standard reference genome, scored against an independent clinical gold standard (Chin et al. 2020).
          </div>
          <div class="table-wrap">
            <table class="result-table">
              <thead><tr><th>Method</th><th>Accuracy</th><th>Score</th></tr></thead>
              <tbody>
                <tr><td class="highlight">ChampHLA (CC)</td><td>3/3</td><td class="highlight">1.000</td></tr>
                <tr><td>MajorityVote</td><td>3/3</td><td>1.000</td></tr>
                <tr><td>OptiType</td><td>3/3</td><td>1.000</td></tr>
                <tr><td>T1K</td><td>3/3</td><td>1.000</td></tr>
                <tr><td>SpecHLA</td><td>3/3</td><td>1.000</td></tr>
                <tr><td>HLA-HD</td><td>2/3</td><td>0.667</td></tr>
              </tbody>
            </table>
          </div>
          <p style="margin-top:10px;font-size:12px;color:var(--muted);">
            n=1 (HG002). Gold truth: A*01:01/26:01 &middot; B*35:08/38:01 &middot; C*04:01/12:03. 95% CI: 0.44&ndash;1.0.
          </p>
        </div>
      </div>
    </section>

    <!-- PIPELINE CHANGES -->
    <section id="{section_id('Pipeline Changes')}">
      <div class="section-header">
        <h2>Pipeline &amp; Infrastructure Changes</h2>
        <p>Code changes, bug fixes, and new modules added to the ChampHLA Nextflow pipeline.</p>
      </div>
      <div class="grid">
        <div class="panel changelog-card">
          <p class="ch-date">Jun 19, 2026</p>
          <p class="ch-title">Seq2HLA BAM-path routing fix</p>
          <p class="ch-body">
            Fixed a bug in <code>main.nf</code> where Seq2HLA was unconditionally skipped for BAM inputs.
            The BAM&rarr;FASTQ conversion path did not include Seq2HLA, causing it to silently produce no results in the RNA benchmark.
            Fix: added Seq2HLA to the <code>need_fastq</code> condition and created a dedicated routing block for RNA BAM inputs.
            After the fix, Seq2HLA was re-run on all 30 RNA benchmark samples (SLURM job 35195038, all 30 completed).
          </p>
        </div>
        <div class="panel changelog-card">
          <p class="ch-date">Jun 19, 2026</p>
          <p class="ch-title">Figure 8 redesign: 3&times;3 &rarr; 1&times;3</p>
          <p class="ch-body">
            Replaced the misleading 3&times;3 layout (rows=HLA genes with fabricated per-gene multiplier) with an honest 1&times;3 layout
            showing three panels: Peak RAM (GB), CPU utilisation (%), and wall-clock time (hours) per tool. Timing summary regenerated
            from fresh Nextflow execution traces including the Seq2HLA re-run data.
          </p>
        </div>
        <div class="panel changelog-card">
          <p class="ch-date">Jun 18, 2026</p>
          <p class="ch-title">SpecHLA containerisation</p>
          <p class="ch-body">
            Created <code>containers/spechla/Dockerfile</code> &mdash; a purpose-built image for SpecHLA v1.1 using Miniconda base
            with compiled dependencies and internal reference databases. Resolves the historical SpecHLA containerisation blocker.
            Container reference added to <code>nextflow.config</code>.
          </p>
        </div>
        <div class="panel changelog-card">
          <p class="ch-date">Jun 24, 2026</p>
          <p class="ch-title">Silver-truth infrastructure</p>
          <p class="ch-body">
            Three new Nextflow modules integrated: <code>modules/locityper.nf</code> (pangenome-aware HLA genotyping),
            <code>modules/immuannot.nf</code> (assembly&rarr;HLA annotation), and <code>modules/silver_truth.nf</code>
            (orchestration workflow). Supporting scripts: <code>bin/generate_silver_truth.py</code>,
            <code>bin/parse_locityper_results.py</code>, <code>bin/parse_immuannot_gtf.py</code>.
            Addresses the manuscript's 2014 two-field resolution ceiling by enabling silver-standard truth derivation
            for cohorts lacking gold-standard typing.
          </p>
        </div>
        <div class="panel changelog-card">
          <p class="ch-date">Jun 18&ndash;19, 2026</p>
          <p class="ch-title">RNA benchmark completion</p>
          <p class="ch-body">
            Completed performance benchmark across all 6 RNA-capable tools (ArcasHLA, HLA-HD, OptiType, Seq2HLA, SpecHLA, T1K)
            on 30 samples. OptiType timed out on 4 samples (n=26); all others completed n=29&ndash;30.
            Timing summary regenerated with correct RNA trace sources.
          </p>
        </div>
      </div>
    </section>

    <!-- COMPUTATIONAL PERFORMANCE -->
    <section id="{section_id('Computational Performance')}">
      <div class="section-header">
        <h2>Computational Performance</h2>
        <p>Resource profiling from the n=30 benchmark across WGS, WES, and RNA-seq modalities. Data from Nextflow execution traces.</p>
      </div>
      <div class="grid metrics-grid" style="margin-bottom:18px;">
        <div class="metric"><div class="label">WES Fastest</div><div class="value">0.008 h</div><div class="note">ArcasHLA &mdash; 0.55 GB RAM</div></div>
        <div class="metric"><div class="label">RNA Lightest</div><div class="value">0.14 GB</div><div class="note">Seq2HLA &mdash; 0.049 h</div></div>
        <div class="metric"><div class="label">WGS Heaviest</div><div class="value">5.9 GB</div><div class="note">Kourami &mdash; 1065% CPU</div></div>
        <div class="metric"><div class="label">RNA Slowest</div><div class="value">0.452 h</div><div class="note">HLA-HD &mdash; 8.7 GB RAM</div></div>
      </div>
      <div class="panel">
        {timing_table}
      </div>
    </section>

    <!-- DOCUMENTATION -->
    <section id="{section_id('Documentation')}">
      <div class="section-header">
        <h2>Documentation &amp; Deployment</h2>
        <p>New user-facing documentation created to support ChampHLA deployment and reproducibility.</p>
      </div>
      <div class="grid two-col">
        <div class="panel">
          <h3>INSTALLATION.md</h3>
          <p>Cross-platform installation guide (9.8 KB) covering HPC, Linux, macOS, and Windows. Includes Docker, Singularity, and local-SpecHLA deployment modes. Created Jun 18.</p>
        </div>
        <div class="panel">
          <h3>USAGE.md</h3>
          <p>Usage guide (7.3 KB) with four deployment scenarios: RNA-seq FASTQ, WES BAM with SpecHLA, DNA+consensus, and SLURM/Puhti batch. Created Jun 18.</p>
        </div>
        <div class="panel">
          <h3>README.md</h3>
          <p>Updated with "Get Running in 5 Minutes" quick-start section (15.4 KB total). Includes tool overview, parameter reference table, and CSC Puhti container paths. Updated Jun 18.</p>
        </div>
        <div class="panel">
          <h3>NEW_HLA_BENCHMARK_DATASETS.md</h3>
          <p>Comprehensive scouting report (13.6 KB) documenting the search for non-1000G external validation cohorts. Rejected options with reasons (consHLA, GeT-RM, 829 WES, 652 RNA). Accepted: NCI-60, GIAB/HPRC, SweHLA. Created Jun 29.</p>
        </div>
      </div>
    </section>

    <!-- MANUSCRIPT UPDATES -->
    <section id="{section_id('Manuscript Updates')}">
      <div class="section-header">
        <h2>Manuscript Updates</h2>
        <p>CHAMPHLA_MANUSCRIPT_V3.md is now the authoritative manuscript version (110.7 KB).</p>
      </div>
      <div class="panel">
        <h3>Key sections added or updated</h3>
        <ul>
          <li><strong>Results &sect;2 &mdash; Computational performance:</strong> All tool timing numbers updated with new TSV medians from the completed RNA benchmark. Seq2HLA data now included.</li>
          <li><strong>Methods &sect;7 &mdash; Trace collection:</strong> Updated to reflect dedicated RNA benchmark, tool-level n-count variation, and OptiType timeout documentation.</li>
          <li><strong>Figure 8 caption:</strong> Changed from "3&times;3 compound figure (rows: HLA-A, -B, -C)" to "Three-panel figure showing peak RAM, CPU utilisation, and wall-clock time" with correct per-tool n-counts.</li>
          <li><strong>Supplementary &sect;S7 (S7.1&ndash;S7.4):</strong> FIMM clinical application sections added &mdash; cross-modality concordance, allele dropout, LOH analysis, and survival/HED.</li>
          <li><strong>Champion-Challenger method:</strong> Override policy, auditable traces, and frozen-policy generalisation results integrated.</li>
          <li><strong>Silver-standard truth methodology:</strong> Locityper + HPRC pangenome + Immuannot approach documented.</li>
          <li><strong>External validation:</strong> NCI-60 and GIAB/HPRC results sections added.</li>
          <li><strong>Deployment planning:</strong> 16 GB / resource-constrained environment guidance added.</li>
        </ul>
      </div>
    </section>

    <!-- NEXT STEPS -->
    <section id="{section_id('Next Steps')}">
      <div class="section-header">
        <h2>Next Steps</h2>
        <p>Ongoing and planned work following this two-week period.</p>
      </div>
      <div class="grid two-col">
        <div class="callout">
          <strong>In progress</strong>
          <ul style="margin-top:8px;padding-left:18px;">
            <li>NCI-60 WES full expansion (60 samples) &mdash; SLURM jobs running with frozen WES Champion-Challenger policy</li>
            <li>NCI-60 RNA full expansion (60 samples) &mdash; manifests built, batched runner ready</li>
          </ul>
        </div>
        <div class="callout warn">
          <strong>On hold (pending decision)</strong>
          <ul style="margin-top:8px;padding-left:18px;">
            <li>HPRC cohort expansion (HG003/HG004 AJ parents + HG005&ndash;HG007 Han trio) &mdash; documented and ready to resume</li>
            <li>SweHLA WGS &mdash; controlled access application required (via NBIS)</li>
            <li>NCI-60 Interim34 results aggregation &mdash; awaiting full-scale job completion</li>
          </ul>
        </div>
      </div>
    </section>

    <footer>
      ChampHLA Progress Report &middot; Generated {safe(generated)} &middot; Figures embedded as base64 &middot; Self-contained offline HTML
    </footer>
  </main>
</div>

<script>
(function() {{
  const nav = document.getElementById('sidebar-nav');
  if (!nav) return;
  const links = nav.querySelectorAll('a[href^="#"]');
  const sections = [];
  links.forEach(link => {{
    const id = link.getAttribute('href').slice(1);
    const el = document.getElementById(id);
    if (el) sections.push({{ link, el }});
  }});
  if (!sections.length) return;
  const observer = new IntersectionObserver(entries => {{
    entries.forEach(entry => {{
      const match = sections.find(s => s.el === entry.target);
      if (match) match.link.classList.toggle('active', entry.isIntersecting);
    }});
  }}, {{ rootMargin: '-20% 0px -70% 0px' }});
  sections.forEach(s => observer.observe(s.el));
}})();
</script>
</body>
</html>"""


def main():
    html_content = build_html()
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(html_content, encoding="utf-8")
    size_mb = OUT.stat().st_size / (1024 * 1024)
    print(f"Report written to {OUT} ({size_mb:.1f} MB)")


if __name__ == "__main__":
    main()
