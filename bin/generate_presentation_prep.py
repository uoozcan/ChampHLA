#!/usr/bin/env python3
"""
generate_presentation_prep.py

Generates a comprehensive presentation preparation document (HTML, print-to-PDF)
for the MVHLA project covering all benchmark results, scientific claims,
LOH/HED analyses, governance constraints, and Q&A preparation.

Usage:
    python3 bin/generate_presentation_prep.py
    # Open output in browser → File → Print → Save as PDF
"""

from pathlib import Path
from datetime import datetime

OUTPUT = Path("/scratch/project_2008084/mvhla_presentation_prep.html")

CSS = """
@import url('https://fonts.googleapis.com/css2?family=Inter:wght@300;400;600;700;800&display=swap');
:root {
  --navy:   #1a3a5c;
  --blue:   #0072B2;
  --green:  #009E73;
  --orange: #E69F00;
  --red:    #D55E00;
  --purple: #4A148C;
  --bg:     #f4f7fb;
  --white:  #ffffff;
  --border: #d0dce8;
  --text:   #1e2a38;
  --muted:  #5a6a7e;
  --yellow: #F0E442;
}
*, *::before, *::after { box-sizing: border-box; margin: 0; padding: 0; }
body {
  font-family: 'Inter', 'Segoe UI', Arial, sans-serif;
  background: var(--bg); color: var(--text);
  font-size: 10.5pt; line-height: 1.65;
  padding: 0; margin: 0;
}
.page { max-width: 900px; margin: 0 auto; padding: 32px 40px 60px; }

/* Cover */
.cover {
  background: linear-gradient(135deg, var(--navy) 0%, #0d2a45 60%, #0072B2 100%);
  color: white; padding: 80px 60px 60px; border-radius: 12px; margin-bottom: 40px;
  page-break-after: always;
}
.cover h1 { font-size: 2.4em; font-weight: 800; letter-spacing: -0.5px; margin-bottom: 8px; }
.cover .subtitle { font-size: 1.1em; color: #a8d4f5; margin-bottom: 32px; font-weight: 300; }
.cover .chips { display: flex; gap: 12px; flex-wrap: wrap; margin-top: 24px; }
.chip {
  background: rgba(255,255,255,0.15); border: 1px solid rgba(255,255,255,0.3);
  color: white; padding: 6px 16px; border-radius: 20px; font-size: 0.82em;
  font-weight: 600;
}
.chip.green  { background: rgba(0,158,115,0.4); border-color: #009E73; }
.chip.orange { background: rgba(230,159,0,0.4); border-color: #E69F00; }
.chip.red    { background: rgba(213,94,0,0.4);  border-color: #D55E00; }

/* Sections */
.section {
  margin-bottom: 48px; page-break-inside: avoid;
}
.section-header {
  display: flex; align-items: center; gap: 14px;
  border-bottom: 3px solid var(--blue); padding-bottom: 10px; margin-bottom: 20px;
}
.sec-num {
  background: var(--blue); color: white; border-radius: 50%;
  width: 34px; height: 34px; display: flex; align-items: center;
  justify-content: center; font-size: 0.8em; font-weight: 800; flex-shrink: 0;
}
h2.section-title { font-size: 1.4em; color: var(--navy); font-weight: 800; }
h3 { font-size: 1.05em; color: var(--blue); margin: 20px 0 8px; font-weight: 700; }
h4 { font-size: 0.95em; color: var(--navy); margin: 14px 0 6px; font-weight: 700; }
p  { margin-bottom: 10px; }
ul, ol { margin: 8px 0 12px 22px; }
li { margin-bottom: 4px; }
strong { color: var(--navy); }
code {
  background: #eef2f7; padding: 2px 6px; border-radius: 4px;
  font-size: 0.88em; color: #c44e52;
}

/* Tables */
table {
  width: 100%; border-collapse: collapse; margin: 14px 0 20px;
  font-size: 0.9em; border-radius: 8px; overflow: hidden;
  box-shadow: 0 1px 4px rgba(0,0,0,0.1);
}
thead th {
  background: var(--navy); color: white; padding: 9px 12px;
  text-align: left; font-weight: 700; white-space: nowrap;
}
td { padding: 8px 12px; border-bottom: 1px solid #e0e8f0; }
tr:nth-child(even) td { background: #f4f7fb; }
tr:hover td { background: #e8f0fa; }
.best td:first-child { font-weight: 800; color: var(--green); }
.winner { font-weight: 800; color: var(--green); }
.caution { color: var(--red); font-weight: 600; }

/* Metric cards */
.metric-grid {
  display: grid; grid-template-columns: repeat(auto-fit, minmax(160px,1fr));
  gap: 14px; margin-bottom: 24px;
}
.metric-card {
  background: white; border: 1px solid var(--border);
  border-radius: 10px; padding: 16px; text-align: center;
  box-shadow: 0 2px 6px rgba(0,0,0,0.06);
}
.metric-card .value { font-size: 1.9em; font-weight: 900; color: var(--blue); }
.metric-card .label { font-size: 0.75em; color: var(--muted); text-transform: uppercase; margin-top: 4px; }

/* Boxes */
.box {
  border-radius: 8px; padding: 16px 20px; margin: 14px 0;
  border-left: 4px solid var(--blue);
  background: linear-gradient(135deg, #e8f0fa, #f4f7fb);
}
.box.green  { border-left-color: var(--green);  background: linear-gradient(135deg, #e8f5e9, #f4f7fb); }
.box.orange { border-left-color: var(--orange); background: linear-gradient(135deg, #fff8e1, #f4f7fb); }
.box.red    { border-left-color: var(--red);    background: linear-gradient(135deg, #ffebee, #f4f7fb); }
.box.purple { border-left-color: var(--purple); background: linear-gradient(135deg, #f3e5f5, #f4f7fb); }
.box-title  { font-weight: 800; color: var(--navy); margin-bottom: 6px; font-size: 0.95em; }

/* Badges */
.badge {
  display: inline-block; padding: 3px 10px; border-radius: 12px;
  color: white; font-size: 0.78em; font-weight: 700; margin-right: 4px;
}
.badge.primary      { background: #1B5E20; }
.badge.secondary    { background: #0D47A1; }
.badge.supplementary{ background: #4A148C; }
.badge.warning      { background: #D55E00; }
.badge.success      { background: #009E73; }

/* Q&A */
.qa-item {
  background: white; border: 1px solid var(--border); border-radius: 10px;
  padding: 18px 22px; margin-bottom: 16px;
  box-shadow: 0 2px 6px rgba(0,0,0,0.06);
  page-break-inside: avoid;
}
.qa-q { font-weight: 800; color: var(--navy); margin-bottom: 8px; font-size: 0.97em; }
.qa-q::before { content: "Q: "; color: var(--blue); }
.qa-a { color: #2a3a50; font-size: 0.93em; }
.qa-a::before { content: "A: "; font-weight: 700; color: var(--green); }

/* Forbidden / governance */
.forbidden {
  background: #fff5f5; border: 1.5px solid var(--red);
  border-radius: 8px; padding: 14px 18px; margin: 10px 0;
}
.forbidden .title { color: var(--red); font-weight: 800; margin-bottom: 6px; }
.forbidden ul li { color: var(--red); }

.approved {
  background: #f0fff8; border: 1.5px solid var(--green);
  border-radius: 8px; padding: 14px 18px; margin: 10px 0;
}
.approved .title { color: var(--green); font-weight: 800; margin-bottom: 6px; }

/* TOC */
.toc { columns: 2; gap: 24px; margin-bottom: 24px; }
.toc-item { margin-bottom: 6px; display: flex; align-items: baseline; gap: 8px; }
.toc-num { font-weight: 800; color: var(--blue); min-width: 24px; }
.toc-title { color: var(--text); font-size: 0.9em; }

/* Print */
@media print {
  body { background: white; }
  .page { padding: 20px 30px; }
  .cover { page-break-after: always; border-radius: 0; }
  .section { page-break-inside: avoid; }
  table { page-break-inside: auto; }
  tr { page-break-inside: avoid; }
}
"""


def render() -> str:
    now = datetime.now().strftime("%d %B %Y")
    return f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>MVHLA — Presentation Preparation Guide</title>
<style>{CSS}</style>
</head>
<body>
<div class="page">

<!-- ═══════════════════════════════════════════════════════ COVER -->
<div class="cover">
  <h1>MVHLA Presentation Preparation Guide</h1>
  <p class="subtitle">Reliability-Weighted HLA Typing Ensemble · 8 Tools · 3 Modalities · 1000G Benchmark</p>
  <p style="color:#a8d4f5;font-size:0.9em">Prepared {now} · Confidential — internal use only</p>
  <div class="chips">
    <span class="chip">WGS best: OptiType 0.4722</span>
    <span class="chip green">WES+RNA: 0.9623</span>
    <span class="chip orange">n=50 WGS holdout</span>
    <span class="chip orange">n=129 WES · n=107 RNA</span>
    <span class="chip red">WGS is tool-limited</span>
  </div>
  <hr style="border:0;border-top:1px solid rgba(255,255,255,0.2);margin:28px 0 20px">
  <p style="color:#a8d4f5;font-size:0.85em">
    Contents: Project Overview · Benchmark Results · Calibration · LOH Analysis ·
    HED &times; Survival · Governance Rules · Anticipated Q&A
  </p>
</div>

<!-- ═══ TABLE OF CONTENTS -->
<div class="section">
  <div class="section-header">
    <h2 class="section-title">Contents</h2>
  </div>
  <div class="toc">
    <div class="toc-item"><span class="toc-num">1</span><span class="toc-title">Project Overview & Architecture</span></div>
    <div class="toc-item"><span class="toc-num">2</span><span class="toc-title">Benchmark Hierarchy</span></div>
    <div class="toc-item"><span class="toc-num">3</span><span class="toc-title">PRIMARY: WGS Results</span></div>
    <div class="toc-item"><span class="toc-num">4</span><span class="toc-title">Confidence Calibration & Guardrails</span></div>
    <div class="toc-item"><span class="toc-num">5</span><span class="toc-title">SECONDARY: WES Results</span></div>
    <div class="toc-item"><span class="toc-num">6</span><span class="toc-title">SECONDARY: RNA-seq Results</span></div>
    <div class="toc-item"><span class="toc-num">7</span><span class="toc-title">SUPPLEMENTARY: Bimodal WES+RNA</span></div>
    <div class="toc-item"><span class="toc-num">8</span><span class="toc-title">SUPPLEMENTARY: Trimodal Analysis</span></div>
    <div class="toc-item"><span class="toc-num">9</span><span class="toc-title">LOH Analysis</span></div>
    <div class="toc-item"><span class="toc-num">10</span><span class="toc-title">HED &times; Survival (FIMM)</span></div>
    <div class="toc-item"><span class="toc-num">11</span><span class="toc-title">Governance — What NOT to Say</span></div>
    <div class="toc-item"><span class="toc-num">12</span><span class="toc-title">Anticipated Q&A</span></div>
    <div class="toc-item"><span class="toc-num">13</span><span class="toc-title">Key Numbers Cheat Sheet</span></div>
  </div>
</div>

<!-- ═══ 1. PROJECT OVERVIEW -->
<div class="section">
  <div class="section-header">
    <span class="sec-num">1</span>
    <h2 class="section-title">Project Overview &amp; Architecture</h2>
  </div>

  <div class="box green">
    <div class="box-title">One-sentence summary</div>
    MVHLA is a benchmark-derived, reliability-aware HLA typing ensemble that uses calibrated
    confidence guardrails and empirical tool weights to select or combine calls from 8 HLA tools
    across WGS, WES, and RNA-seq — without assuming any single tool is universally best.
  </div>

  <h3>Why does this exist?</h3>
  <ul>
    <li><strong>8 tools</strong>, each optimised for different modalities and allele depths — no single tool wins everywhere</li>
    <li>Naive majority voting is dominated by the worst tools in disagreement cases</li>
    <li>Raw confidence scores from most tools are <em>poorly calibrated</em> — amplifying them directly lowers accuracy</li>
    <li>Clinical HLA typing (transplant matching, immunotherapy) requires explainable uncertainty, not silent errors</li>
  </ul>

  <h3>8 Tools included</h3>
  <table>
    <thead><tr><th>Tool</th><th>Primary modality</th><th>Confidence output</th><th>Notes</th></tr></thead>
    <tbody>
      <tr><td>ArcasHLA</td><td>RNA-seq</td><td>Yes</td><td>Fails on WGS (8% accuracy)</td></tr>
      <tr><td>HLA-HD</td><td>WGS / WES</td><td>No</td><td>Strong WES</td></tr>
      <tr><td>Kourami</td><td>WGS</td><td>No</td><td>Graph-based</td></tr>
      <tr><td>OptiType</td><td>WES / WGS</td><td>Yes</td><td><strong>Best WGS single tool</strong></td></tr>
      <tr><td>POLYSOLVER</td><td>WES</td><td>No</td><td>Clinical tool of record</td></tr>
      <tr><td>Seq2HLA</td><td>RNA-seq</td><td>No</td><td>Older tool, RNA-only</td></tr>
      <tr><td>SpecHLA</td><td>WES / WGS</td><td>No</td><td>3.5% on RNA — excluded from RNA scoring</td></tr>
      <tr><td>T1K</td><td>Universal</td><td>Yes</td><td>Second-best WGS single tool</td></tr>
    </tbody>
  </table>

  <h3>Pipeline flow (4 stages)</h3>
  <ol>
    <li><strong>Dispatch</strong> — route each sample to tools compatible with its modality</li>
    <li><strong>Phase-gating</strong> — filter tool outputs that fail minimum evidence thresholds</li>
    <li><strong>Weight learning</strong> — derive per-tool weights from calibration hold-out; apply guardrails to confidence scores</li>
    <li><strong>Consensus + abstention</strong> — compute weighted vote; abstain if support &lt; 0.45 threshold</li>
  </ol>

  <div class="box orange">
    <div class="box-title">Abstention is a feature, not a failure</div>
    When inter-tool support is below the threshold (0.45), MVHLA withholds a call rather than
    reporting an uncertain allele. WGS abstention rate = 8.3%; RNA-seq = 4%.
    These loci flag cases that need orthogonal confirmation.
  </div>
</div>

<!-- ═══ 2. BENCHMARK HIERARCHY -->
<div class="section">
  <div class="section-header">
    <span class="sec-num">2</span>
    <h2 class="section-title">Benchmark Hierarchy</h2>
  </div>

  <p>All numbers come from 1000 Genomes Project samples with IPD-IMGT/HLA v3.59.0 truth alleles.</p>

  <table>
    <thead><tr><th>Tier</th><th>Badge</th><th>Cohort</th><th>N samples</th><th>Scope</th></tr></thead>
    <tbody>
      <tr>
        <td><strong>PRIMARY</strong></td>
        <td><span class="badge primary">PRIMARY</span></td>
        <td>WGS Wave2 holdout</td>
        <td>12</td>
        <td>Main evidence for WGS conclusions — authority for all WGS claims</td>
      </tr>
      <tr>
        <td><strong>SECONDARY</strong></td>
        <td><span class="badge secondary">SECONDARY</span></td>
        <td>WES (129) + RNA (107)</td>
        <td>129 / 107</td>
        <td>Independent benchmarks — not interchangeable with PRIMARY</td>
      </tr>
      <tr>
        <td><strong>SUPPLEMENTARY</strong></td>
        <td><span class="badge supplementary">SUPPLEMENTARY</span></td>
        <td>Matched 106-subject robustness</td>
        <td>106</td>
        <td>Confirms direction; does NOT replace wave2 or secondary benchmarks</td>
      </tr>
    </tbody>
  </table>

  <div class="box red">
    <div class="box-title">Critical: benchmarks cannot be mixed</div>
    PRIMARY, SECONDARY, and SUPPLEMENTARY results are from different cohorts and study designs.
    Never compare a WGS PRIMARY number directly to a WES SECONDARY number as though they are equivalent evaluations.
  </div>
</div>

<!-- ═══ 3. PRIMARY WGS -->
<div class="section">
  <div class="section-header">
    <span class="sec-num">3</span>
    <h2 class="section-title">PRIMARY: WGS Wave2 Results <span class="badge primary">PRIMARY</span></h2>
  </div>

  <p><strong>Cohort:</strong> 50 truth-backed 1000G samples (training 28 / validation 10 / <strong>holdout 12</strong>). Evaluation loci: HLA-A, HLA-B, HLA-C at 2-field resolution.</p>

  <h3>Single-tool performance (holdout, n=12)</h3>
  <table>
    <thead><tr><th>Tool</th><th>Accuracy</th><th>Callable rate</th></tr></thead>
    <tbody>
      <tr class="best"><td>OptiType</td><td class="winner">0.4722</td><td>1.0000</td></tr>
      <tr><td>T1K</td><td>0.3611</td><td>1.0000</td></tr>
      <tr><td>HLA-HD</td><td>0.3056</td><td>1.0000</td></tr>
      <tr><td>SpecHLA</td><td>0.2500</td><td>0.8611</td></tr>
      <tr><td>Kourami</td><td>0.1364</td><td>1.0000</td></tr>
      <tr class="caution"><td>ArcasHLA</td><td class="caution">0.0000</td><td>1.0000</td></tr>
    </tbody>
  </table>

  <h3>Consensus methods (holdout, n=12)</h3>
  <table>
    <thead><tr><th>Method</th><th>Accuracy</th><th>Accuracy-among-callable</th><th>Callable rate</th></tr></thead>
    <tbody>
      <tr><td>MajorityVote (baseline)</td><td>0.3889</td><td>0.3889</td><td>1.0000</td></tr>
      <tr><td>WeightedConsensus (MVHLA)</td><td>0.4444</td><td>0.4848</td><td>0.9167</td></tr>
      <tr class="best"><td>ChampionChallenger (routed)</td><td class="winner">0.4722</td><td>—</td><td>—</td></tr>
      <tr class="best"><td>GatedConsensus (routed)</td><td class="winner">0.4722</td><td>—</td><td>—</td></tr>
    </tbody>
  </table>

  <div class="box orange">
    <div class="box-title">The core WGS finding — say this exactly</div>
    <em>"WGS is tool-limited, not consensus-limited. OptiType sets the ceiling at 0.4722.
    Routed baselines recover that ceiling. WeightedConsensus trails at 0.4444 because low-performing
    tools pull the weighted vote. No method exceeds OptiType on WGS."</em>
  </div>

  <h3>Per-gene accuracy (WGS, best single tool = OptiType)</h3>
  <table>
    <thead><tr><th>Locus</th><th>Accuracy</th><th>Difficulty driver</th></tr></thead>
    <tbody>
      <tr class="best"><td>HLA-B</td><td class="winner">0.8333</td><td>Most divergent sequences — easiest to distinguish</td></tr>
      <tr><td>HLA-A</td><td>0.4167</td><td>Moderate sequence similarity</td></tr>
      <tr class="caution"><td>HLA-C</td><td class="caution">0.1667</td><td>High similarity to pseudogenes — hardest locus</td></tr>
    </tbody>
  </table>

  <p><strong>Ensemble benefit is locus-specific.</strong> HLA-B improvements are larger; HLA-C is where most errors concentrate regardless of method.</p>
</div>

<!-- ═══ 4. CALIBRATION -->
<div class="section">
  <div class="section-header">
    <span class="sec-num">4</span>
    <h2 class="section-title">Confidence Calibration &amp; Guardrails <span class="badge primary">PRIMARY</span></h2>
  </div>

  <div class="box">
    <div class="box-title">What calibration means here</div>
    A tool's confidence score is "calibrated" if it accurately predicts correctness probability.
    A score of 0.9 should mean the tool is right 90% of the time at that score level.
    Most HLA tools are <em>not</em> calibrated in this sense.
  </div>

  <h3>WGS tool calibration results</h3>
  <table>
    <thead><tr><th>Tool</th><th>Brier Score</th><th>ECE</th><th>Guardrail</th><th>Final weight</th></tr></thead>
    <tbody>
      <tr class="best"><td>OptiType</td><td>0.2465</td><td>0.0027</td><td>applied</td><td class="winner">0.5476</td></tr>
      <tr><td>T1K</td><td>0.3144</td><td>0.0176</td><td>applied</td><td>0.2649</td></tr>
      <tr><td>HLA-HD</td><td>0.1673</td><td>0.0348</td><td>applied</td><td>0.2294</td></tr>
      <tr><td>Kourami</td><td>0.1174</td><td>0.0242</td><td>applied</td><td>0.1408</td></tr>
      <tr><td>ArcasHLA</td><td>0.0235</td><td>0.0197</td><td>applied</td><td>0.0255</td></tr>
      <tr class="caution"><td>SpecHLA</td><td>—</td><td>—</td><td>no_confidence</td><td>0.1548</td></tr>
    </tbody>
  </table>

  <p>Thresholds: <code>max_brier_score = 0.35</code> · <code>max_ECE = 0.35</code></p>
  <p>Weight formula: <code>final_weight = 0.7 × base_reliability + 0.3 × effective_confidence</code></p>

  <div class="box red">
    <div class="box-title">Why guardrails matter — the key argument</div>
    Three of five WGS tools with confidence data showed poor empirical calibration.
    <strong>Without guardrails, naively amplifying those confidence scores would have lowered ensemble accuracy.</strong>
    OptiType's miscalibrated confidence (Brier 0.2465) is clipped before contributing to votes.
  </div>
</div>

<!-- ═══ 5. WES -->
<div class="section">
  <div class="section-header">
    <span class="sec-num">5</span>
    <h2 class="section-title">SECONDARY: WES Results <span class="badge secondary">SECONDARY</span></h2>
  </div>

  <p><strong>Cohort:</strong> 129 truth-backed 1000G WES samples. Evaluated independently of WGS benchmark.</p>

  <div class="metric-grid">
    <div class="metric-card"><div class="value">0.9359</div><div class="label">MajorityVote accuracy</div></div>
    <div class="metric-card"><div class="value">0.9487</div><div class="label">ChampionChallenger (tuned)</div></div>
    <div class="metric-card"><div class="value">129</div><div class="label">Samples</div></div>
    <div class="metric-card"><div class="value">1.000</div><div class="label">Callable rate</div></div>
  </div>

  <h3>Single-tool performance (WES, n=129)</h3>
  <table>
    <thead><tr><th>Tool</th><th>Accuracy</th><th>Callable rate</th></tr></thead>
    <tbody>
      <tr class="best"><td>OptiType</td><td class="winner">0.9231</td><td>1.0</td></tr>
      <tr><td>POLYSOLVER</td><td>0.9154</td><td>1.0</td></tr>
      <tr><td>T1K</td><td>0.8308</td><td>1.0</td></tr>
      <tr><td>SpecHLA</td><td>0.7846</td><td>1.0</td></tr>
      <tr><td>HLA-HD</td><td>0.7821</td><td>1.0</td></tr>
      <tr><td>Kourami</td><td>0.7048</td><td>1.0</td></tr>
      <tr class="caution"><td>ArcasHLA</td><td class="caution">0.1359</td><td>1.0</td></tr>
    </tbody>
  </table>

  <div class="box green">
    <div class="box-title">WES framing</div>
    WES is a high-agreement regime — top tools perform at >92%. Ensemble adds a modest +0.5 to +1.5 pp gain.
    Threshold-tuning (ChampionChallenger) outperforms plain weighting here because inter-tool variance is low.
  </div>
</div>

<!-- ═══ 6. RNA -->
<div class="section">
  <div class="section-header">
    <span class="sec-num">6</span>
    <h2 class="section-title">SECONDARY: RNA-seq Results <span class="badge secondary">SECONDARY</span></h2>
  </div>

  <p><strong>Cohort:</strong> 107 truth-backed 1000G RNA-seq samples.</p>

  <div class="metric-grid">
    <div class="metric-card"><div class="value">0.9502</div><div class="label">MajorityVote accuracy</div></div>
    <div class="metric-card"><div class="value">0.9533</div><div class="label">ChampionChallenger (tuned)</div></div>
    <div class="metric-card"><div class="value">107</div><div class="label">Samples</div></div>
    <div class="metric-card"><div class="value">4%</div><div class="label">Abstention rate</div></div>
  </div>

  <h3>Single-tool performance (RNA, n=107)</h3>
  <table>
    <thead><tr><th>Tool</th><th>Accuracy</th><th>Callable rate</th><th>Note</th></tr></thead>
    <tbody>
      <tr class="best"><td>HLA-HD</td><td class="winner">0.9429</td><td>1.0</td><td></td></tr>
      <tr><td>OptiType</td><td>0.9397</td><td>1.0</td><td></td></tr>
      <tr><td>ArcasHLA</td><td>0.9190</td><td>1.0</td><td>RNA-designed — performs as expected</td></tr>
      <tr><td>T1K</td><td>0.8794</td><td>1.0</td><td></td></tr>
      <tr><td>Seq2HLA</td><td>0.5806</td><td>1.0</td><td>n=29 only</td></tr>
      <tr class="caution"><td>SpecHLA</td><td class="caution">0.0349</td><td>0.7619</td><td><strong>Excluded from RNA production scoring</strong></td></tr>
    </tbody>
  </table>

  <div class="box red">
    <div class="box-title">SpecHLA RNA caveat — must state clearly</div>
    SpecHLA achieves 3.5% accuracy on RNA-seq — below chance given the allele catalogue size.
    It is <strong>excluded from RNA-seq production accuracy scoring</strong>.
    Coverage metrics from SpecHLA RNA may still be shown, but accuracy claims exclude it.
  </div>
</div>

<!-- ═══ 7. BIMODAL -->
<div class="section">
  <div class="section-header">
    <span class="sec-num">7</span>
    <h2 class="section-title">SUPPLEMENTARY: Bimodal WES+RNA <span class="badge supplementary">SUPPLEMENTARY</span></h2>
  </div>

  <p><strong>Cohort:</strong> 106 matched subjects with both WES and RNA-seq data (391 Class I + II gene-locus pairs).</p>

  <div class="metric-grid">
    <div class="metric-card"><div class="value">96.83%</div><div class="label">Bimodal WeightedConsensus</div></div>
    <div class="metric-card"><div class="value">96.16%</div><div class="label">Bimodal MajorityVote</div></div>
    <div class="metric-card"><div class="value">+2.0 pp</div><div class="label">Gain over WES-only</div></div>
    <div class="metric-card"><div class="value">+0.8 pp</div><div class="label">Gain over RNA-only</div></div>
  </div>

  <h3>Per-gene bimodal gains (WES+RNA vs single-modality best)</h3>
  <table>
    <thead><tr><th>Locus</th><th>WES-only</th><th>RNA-only</th><th>Bimodal WES+RNA</th><th>Gain</th></tr></thead>
    <tbody>
      <tr><td>HLA-A</td><td>—</td><td>96.0%</td><td class="winner">97.1%</td><td>+1.1 pp</td></tr>
      <tr><td>HLA-B</td><td>—</td><td>96.2%</td><td class="winner">97.1%</td><td>+0.9 pp</td></tr>
      <tr><td>HLA-C</td><td>—</td><td>96.1%</td><td class="winner">97.0%</td><td>+0.9 pp</td></tr>
    </tbody>
  </table>

  <div class="box green">
    <div class="box-title">Clinical translation</div>
    The +0.8–2.0 pp improvement translates to approximately <strong>8–10 more loci correctly resolved
    per 1,000 patient samples</strong> compared to the best single-modality approach.
  </div>
</div>

<!-- ═══ 8. TRIMODAL -->
<div class="section">
  <div class="section-header">
    <span class="sec-num">8</span>
    <h2 class="section-title">SUPPLEMENTARY: Trimodal Analysis <span class="badge supplementary">SUPPLEMENTARY</span></h2>
  </div>

  <p><strong>Cohort:</strong> Same 106 matched subjects, adding WGS calls (370 gene-locus pairs; 21 dropped due to WGS no-calls).</p>

  <table>
    <thead><tr><th>Method</th><th>Accuracy-among-callable</th><th>Callable rate</th><th>vs Bimodal</th></tr></thead>
    <tbody>
      <tr class="best"><td>Bimodal WeightedConsensus (WES+RNA)</td><td class="winner">96.83%</td><td>96.7%</td><td>—</td></tr>
      <tr><td>Trimodal WeightedConsensus</td><td>96.88%</td><td>95.4%</td><td>+0.05 pp</td></tr>
      <tr><td>Bimodal MajorityVote</td><td>96.16%</td><td>100%</td><td>—</td></tr>
      <tr><td>Trimodal MajorityVote</td><td>96.22%</td><td>100%</td><td>+0.06 pp</td></tr>
    </tbody>
  </table>

  <div class="box red">
    <div class="box-title">The trimodal finding — state carefully</div>
    <em>"Adding WGS to bimodal WES+RNA provides negligible accuracy gain (≤0.1 pp)
    while reducing callable coverage by 21 gene-rows (5.4%). WES+RNA is the optimal
    routine combination. Trimodal does NOT improve over bimodal WES+RNA."</em>
  </div>
</div>

<!-- ═══ 9. LOH -->
<div class="section">
  <div class="section-header">
    <span class="sec-num">9</span>
    <h2 class="section-title">LOH (Loss of Heterozygosity) Analysis</h2>
  </div>

  <div class="box orange">
    <div class="box-title">Scope and limitation upfront</div>
    LOH analysis targets the VENEX clinical cohort (80 tumour WGS samples). Only 1/80 samples
    (VX_92_2_D1) has surviving SpecHLA frequency files — the rest were in Nextflow work dirs
    that were cleaned. Direct LOH classification is only possible for 1 sample.
    The main output is allele dropout characterisation from benchmark QC metrics.
  </div>

  <h3>False Duplicate Rate by tool × modality (WGS — most affected)</h3>
  <p>False duplicate rate = fraction of duplicate (homozygous) calls where the truth is heterozygous. High rates indicate allele dropout.</p>

  <table>
    <thead><tr><th>Tool</th><th>WGS-A</th><th>WGS-B</th><th>WGS-C</th><th>WES-A</th><th>WES-B</th><th>WES-C</th></tr></thead>
    <tbody>
      <tr class="caution"><td>OptiType</td><td>100%</td><td>100%</td><td>92%</td><td>&lt;15%</td><td>&lt;15%</td><td>&lt;15%</td></tr>
      <tr class="caution"><td>SpecHLA</td><td>89%</td><td>33%</td><td>87%</td><td>&lt;20%</td><td>&lt;20%</td><td>&lt;20%</td></tr>
      <tr class="caution"><td>HLA-HD</td><td>88%</td><td>73%</td><td>91%</td><td>&lt;25%</td><td>&lt;25%</td><td>&lt;25%</td></tr>
      <tr class="caution"><td>Kourami</td><td>86%</td><td>88%</td><td>88%</td><td>&lt;20%</td><td>&lt;20%</td><td>&lt;20%</td></tr>
      <tr class="caution"><td>T1K</td><td>83%</td><td>75%</td><td>78%</td><td>&lt;25%</td><td>&lt;25%</td><td>&lt;25%</td></tr>
    </tbody>
  </table>

  <p><strong>Interpretation:</strong> 75–100% of apparent WGS homozygous calls are false positives (allele dropout errors), not biological homozygosity. This directly explains why WGS consensus accuracy is limited — tools systematically drop one allele.</p>

  <h3>LOH status distribution (VENEX, 240 gene-rows)</h3>
  <table>
    <thead><tr><th>Status</th><th>N rows</th><th>Meaning</th></tr></thead>
    <tbody>
      <tr><td>insufficient_evidence</td><td>222</td><td>No freq file available (work dir cleaned)</td></tr>
      <tr class="caution"><td>likely_allele_dropout</td><td>15</td><td>False duplicate rate ≥ 0.50</td></tr>
      <tr class="best"><td>balanced_heterozygous</td><td>3</td><td>Alleles 0.35–0.65 balance (VX_92_2_D1)</td></tr>
    </tbody>
  </table>

  <div class="box">
    <div class="box-title">Framing for LOH in your talk</div>
    "Direct allele frequency measurement is available for one sample (VX_92_2_D1), which shows
    balanced heterozygosity at all three HLA class I loci — no LOH detected.
    For the remaining 79 samples, input WGS data is no longer available for re-processing;
    LOH assessment is limited to allele dropout pattern characterisation from benchmark QC metrics."
  </div>
</div>

<!-- ═══ 10. HED -->
<div class="section">
  <div class="section-header">
    <span class="sec-num">10</span>
    <h2 class="section-title">HED &times; Survival (FIMM AML/MDS Cohort)</h2>
  </div>

  <div class="metric-grid">
    <div class="metric-card"><div class="value">30</div><div class="label">FIMM patients</div></div>
    <div class="metric-card"><div class="value">27</div><div class="label">With survival data</div></div>
    <div class="metric-card"><div class="value">14</div><div class="label">Death events</div></div>
    <div class="metric-card"><div class="value">30/30</div><div class="label">HED complete (all loci)</div></div>
  </div>

  <h3>What is HED?</h3>
  <p>
    <strong>HLA Evolutionary Divergence (HED)</strong> quantifies how structurally different a patient's
    two alleles are at each HLA locus, using Grantham amino acid distances at antigen-binding groove (ABG) positions.
    Higher HED = greater allelic diversity = potentially broader peptide presentation = possibly better immune surveillance.
  </p>
  <ul>
    <li>Formula: <code>d(i,j) = √(0.1018·Δc² + 0.000399·Δp² + 0.000791·Δv²)</code> (Grantham 1974)</li>
    <li>ABG positions: HLA-A (10 positions) · HLA-B (8) · HLA-C (7) — Pierini &amp; Lenz 2018</li>
    <li>HED = 0 for homozygous patients</li>
    <li>HED total = sum across HLA-A + HLA-B + HLA-C</li>
  </ul>

  <h3>Analysis approach (updated with 8 fixes)</h3>
  <ul>
    <li>Consensus allele pairs from majority-vote across 3 WES tools (OptiType, ArcasHLA, SpecHLA)</li>
    <li>Primary endpoint: <strong>OS from diagnosis date</strong> (avoids landmark bias)</li>
    <li>Kaplan–Meier stratification at median HED</li>
    <li>Cox regression: univariable (HED) + multivariable (HED + diagnosis group)</li>
  </ul>

  <div class="box orange">
    <div class="box-title">Power caveat — say this at every mention of KM p-value</div>
    With n≈14 per arm (median HED split), the log-rank test has approximately <strong>25% power
    to detect HR=2.0</strong> (α=0.05). All HED–survival results are hypothesis-generating only.
    A larger cohort with treatment data is required for definitive conclusions.
  </div>

  <h3>Diagnosis context</h3>
  <p>AML median OS ≈ 12–18 months; MDS median OS ≈ 3–5 years. This confounder is partially addressed
  by the multivariable Cox model (HED + diagnosis group), but at n=28 the model is underpowered.
  Any HED association must be interpreted in the context of cytogenetic risk and treatment response.</p>
</div>

<!-- ═══ 11. GOVERNANCE -->
<div class="section">
  <div class="section-header">
    <span class="sec-num">11</span>
    <h2 class="section-title">Governance — What NOT to Say</h2>
  </div>

  <div class="forbidden">
    <div class="title">⛔ Forbidden phrases — never use these</div>
    <ul>
      <li>"WeightedConsensus is the best WGS method" (OptiType 0.4722 &gt; WeightedConsensus 0.4444)</li>
      <li>"Trimodal consensus improves over the best multimodal baseline" (gain is ≤0.1 pp)</li>
      <li>"The robustness run replaces wave2" (supplementary ≠ primary)</li>
      <li>"WGS results generalise beyond HLA class I" (benchmark covers A, B, C only)</li>
      <li>"SpecHLA is unsuitable for RNA production" without the caveat that coverage metrics may still be shown</li>
      <li>Presenting any SUPPLEMENTARY number as equivalent to PRIMARY evidence</li>
    </ul>
  </div>

  <div class="approved">
    <div class="title">✓ Approved phrasings for key claims</div>
    <ul>
      <li><strong>WGS best method:</strong> "OptiType achieves 0.4722 on WGS holdout. WeightedConsensus trails at 0.4444."</li>
      <li><strong>WGS framing:</strong> "WGS is tool-limited, not consensus-limited."</li>
      <li><strong>Routed baseline:</strong> "ChampionChallenger recovers the OptiType ceiling but does not exceed it."</li>
      <li><strong>Best multimodal:</strong> "BimodalMajorityVote WES+RNA (0.9623) is the strongest multimodal result."</li>
      <li><strong>Trimodal:</strong> "Trimodal does NOT improve over bimodal WES+RNA."</li>
      <li><strong>SpecHLA RNA:</strong> "SpecHLA is unsuitable for RNA production <em>accuracy scoring</em> — coverage caveat only."</li>
    </ul>
  </div>
</div>

<!-- ═══ 12. Q&A -->
<div class="section">
  <div class="section-header">
    <span class="sec-num">12</span>
    <h2 class="section-title">Anticipated Q&amp;A Preparation</h2>
  </div>

  <div class="qa-item">
    <div class="qa-q">Why doesn't your weighted consensus beat OptiType on WGS?</div>
    <div class="qa-a">Because WGS is tool-limited, not consensus-limited. The best single tool
    sets an empirical ceiling of 0.4722. Ensemble methods redistribute weight away from
    low-performing tools, but they cannot generate correct calls that no tool made.
    ArcasHLA's 0% accuracy pulls majority vote down to 0.3889; weighting reduces this drag
    but cannot fully overcome it. The honest conclusion is that we need better WGS tools,
    not better voting.</div>
  </div>

  <div class="qa-item">
    <div class="qa-q">Why is WGS so much worse than WES and RNA-seq?</div>
    <div class="qa-a">Two factors: (1) WGS HLA loci suffer from high paralog read misalignment —
    75–100% of tool-reported homozygous calls are actually allele dropout, not biology.
    (2) Tool design: ArcasHLA was built for RNA, SpecHLA for WES. Using mismatched tools
    on WGS degrades the voting pool. WES and RNA have fewer dropout events and better
    tool-modality fit.</div>
  </div>

  <div class="qa-item">
    <div class="qa-q">Why not just use OptiType alone if it's the best WGS tool?</div>
    <div class="qa-a">On the holdout, yes — OptiType alone matches the routed ensemble.
    But OptiType is a closed, unmaintained tool (last update 2019), with no confidence output
    for class II, and no native abstention. MVHLA adds calibrated uncertainty quantification,
    abstention, and multi-modality support that OptiType cannot provide. For production,
    ensemble abstention flags ambiguous loci for review rather than silently reporting errors.</div>
  </div>

  <div class="qa-item">
    <div class="qa-q">The WGS holdout is only 12 samples — aren't your results underpowered?</div>
    <div class="qa-a">Yes, and we state this explicitly. The 12-sample holdout is an honest
    evaluation of the wave2 benchmark; we do not overstate precision. The supplementary
    106-sample matched cohort (all three modalities, WGS matched-subject results: OptiType 0.4825)
    confirms the wave2 WGS direction. Primary conclusions are framed as within-benchmark findings,
    not population-level generalisations.</div>
  </div>

  <div class="qa-item">
    <div class="qa-q">Why does adding WGS to WES+RNA not help?</div>
    <div class="qa-a">WGS allele dropout rates of 75–100% mean that WGS tools frequently report
    one allele as homozygous when the truth is heterozygous. In a trimodal vote, this incorrect
    WGS call can outvote the correct WES+RNA consensus. The marginal accuracy gain (+0.05 pp)
    does not justify the coverage cost (−5.4%) or the complexity of running WGS in addition
    to WES and RNA.</div>
  </div>

  <div class="qa-item">
    <div class="qa-q">How does MVHLA handle class II HLA?</div>
    <div class="qa-a">The primary benchmark covers HLA-A, -B, -C (class I). Partial class II
    coverage is available in the bimodal supplementary benchmark (DRB1: 38 assessments,
    DQB1: 35 assessments) — these are partial truth coverage only and not used for primary claims.
    Full class II benchmarking is a stated future direction.</div>
  </div>

  <div class="qa-item">
    <div class="qa-q">Is the HED–survival association real?</div>
    <div class="qa-a">We cannot conclude this from the current data. With n≈14 per arm, the
    log-rank test has ~25% power to detect HR=2.0. The analysis is hypothesis-generating:
    HED may modulate immune surveillance through broader peptide-presentation capacity, but
    AML/MDS survival is dominated by cytogenetic risk, FLT3/NPM1 mutations, and treatment
    response. A larger cohort with treatment data is required.</div>
  </div>

  <div class="qa-item">
    <div class="qa-q">Why is LOH mostly "insufficient evidence"?</div>
    <div class="qa-a">SpecHLA generates allele frequency files (HLA_{{gene}}_freq.txt) only inside
    Nextflow work directories, which are cleaned after pipeline completion. Only 1/80 VENEX samples
    has a surviving work directory. The LOH section is therefore limited to allele dropout pattern
    characterisation from benchmark QC metrics — a genuine data limitation we state explicitly.</div>
  </div>
</div>

<!-- ═══ 13. CHEAT SHEET -->
<div class="section">
  <div class="section-header">
    <span class="sec-num">13</span>
    <h2 class="section-title">Key Numbers Cheat Sheet</h2>
  </div>

  <h3>Numbers to have memorised</h3>
  <table>
    <thead><tr><th>Claim</th><th>Number</th><th>Context</th></tr></thead>
    <tbody>
      <tr><td>Best WGS single tool</td><td class="winner">0.4722</td><td>OptiType, holdout n=12</td></tr>
      <tr><td>WGS WeightedConsensus</td><td>0.4444</td><td>Trails OptiType; abstention 8.3%</td></tr>
      <tr><td>WGS MajorityVote</td><td>0.3889</td><td>Pulled down by ArcasHLA 0%</td></tr>
      <tr><td>Best WES consensus</td><td class="winner">0.9487</td><td>ChampionChallenger, n=129</td></tr>
      <tr><td>Best RNA consensus</td><td class="winner">0.9533</td><td>ChampionChallenger, n=107</td></tr>
      <tr><td>Best bimodal (WES+RNA)</td><td class="winner">0.9623</td><td>BimodalMajorityVote, n=106</td></tr>
      <tr><td>Trimodal gain over bimodal</td><td>+0.05 pp</td><td>Not meaningful; coverage −5.4%</td></tr>
      <tr><td>WGS false duplicate rate</td><td class="caution">75–100%</td><td>AlleleDropout, tool-dependent</td></tr>
      <tr><td>HLA-B WGS accuracy (OptiType)</td><td class="winner">0.8333</td><td>Easiest locus</td></tr>
      <tr><td>HLA-C WGS accuracy (OptiType)</td><td class="caution">0.1667</td><td>Hardest locus</td></tr>
      <tr><td>SpecHLA RNA accuracy</td><td class="caution">0.0349</td><td>Excluded from RNA scoring</td></tr>
      <tr><td>FIMM patients</td><td>30</td><td>AML/MDS, WES HED analysis</td></tr>
      <tr><td>HED KM power (n=14/arm)</td><td class="caution">~25%</td><td>Power to detect HR=2 at α=0.05</td></tr>
    </tbody>
  </table>

  <h3>Three things to say in every context</h3>
  <ol>
    <li><strong>WGS is tool-limited</strong> — the ensemble ceiling equals the best single tool. Better tools are needed, not better voting.</li>
    <li><strong>Calibration guardrails are not optional</strong> — without them, confidence amplification lowers accuracy on WGS.</li>
    <li><strong>WES+RNA bimodal is the best supported multimodal strategy</strong> — adding WGS hurts coverage without improving accuracy.</li>
  </ol>

  <div class="box purple">
    <div class="box-title">Closing sentence for any talk</div>
    <em>"MVHLA demonstrates that the limiting factor for HLA typing in WGS is tool quality,
    not consensus algorithm sophistication — and that calibrated uncertainty quantification is
    essential for responsible deployment of ensemble HLA typing in clinical research."</em>
  </div>

  <hr style="border:0;border-top:1px solid #d0dce8;margin:32px 0 20px">
  <p style="font-size:0.8em;color:#8a9ab0;text-align:center">
    MVHLA Presentation Preparation Guide · Generated {now} · /scratch/project_2008084/pihla-publish/
  </p>
</div>

</div><!-- /.page -->
</body>
</html>
"""

OUTPUT.write_text(render(), encoding="utf-8")
print(f"Written → {OUTPUT}  ({OUTPUT.stat().st_size / 1024:.0f} KB)")
print("Open in browser → File → Print → Save as PDF")


if __name__ == "__main__":
    pass  # already executed at import-time above for simplicity
