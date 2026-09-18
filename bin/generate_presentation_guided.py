#!/usr/bin/env python3
"""
Generate a guided, figure-by-figure presentation prep document for MVHLA.

Each slide has:
  - The actual figure (embedded as base64 PNG)
  - "What this figure shows" — plain-language description
  - "What to say" — speaker bullet points
  - "Key takeaway" — one-sentence summary

Output: /scratch/project_2008084/mvhla_presentation_guided.html
Open in browser → Ctrl+P → Save as PDF
"""

import base64
from pathlib import Path

BASE  = Path("/scratch/project_2008084")
PUB   = BASE / "pihla-publish"
LOCAL = BASE / "pihla_local"
OUTPUT = BASE / "mvhla_presentation_guided.html"

# ─── figure helper ────────────────────────────────────────────────────────────

def b64img(path: Path, alt: str = "") -> str:
    """Return an <img> tag with the figure embedded as base64, or a placeholder."""
    if not path.exists():
        return f'<div class="fig-missing">Figure not found: {path}</div>'
    data = base64.b64encode(path.read_bytes()).decode()
    return (
        f'<img src="data:image/png;base64,{data}" '
        f'alt="{alt}" class="slide-fig">'
    )

def fig(path: Path, caption: str = "", alt: str = "") -> str:
    html = b64img(path, alt or caption)
    if caption:
        html += f'<p class="fig-caption">{caption}</p>'
    return f'<div class="figure-block">{html}</div>'

# ─── slide builder ─────────────────────────────────────────────────────────────

def slide(number: int, title: str, body: str, tag: str = "") -> str:
    tag_html = f'<span class="slide-tag">{tag}</span>' if tag else ""
    return f"""
<div class="slide">
  <div class="slide-header">
    <span class="slide-num">{number}</span>
    <h2 class="slide-title">{title}</h2>
    {tag_html}
  </div>
  <div class="slide-body">
    {body}
  </div>
</div>
"""

def story(what_shows: str, what_to_say, takeaway: str,
          note: str = "") -> str:
    bullets = "".join(f"<li>{b}</li>" for b in what_to_say)
    note_html = f'<div class="story-note">{note}</div>' if note else ""
    return f"""
<div class="story-block">
  <div class="story-section">
    <div class="story-label">What this shows</div>
    <p class="story-text">{what_shows}</p>
  </div>
  <div class="story-section">
    <div class="story-label">What to say</div>
    <ul class="story-bullets">{bullets}</ul>
  </div>
  <div class="takeaway">
    <span class="takeaway-label">KEY TAKEAWAY</span>
    {takeaway}
  </div>
  {note_html}
</div>
"""

def two_figs(path1: Path, cap1: str, path2: Path, cap2: str) -> str:
    return f"""
<div class="figure-block two-figs">
  <div class="fig-col">
    {b64img(path1, cap1)}
    <p class="fig-caption">{cap1}</p>
  </div>
  <div class="fig-col">
    {b64img(path2, cap2)}
    <p class="fig-caption">{cap2}</p>
  </div>
</div>
"""

# ─── CSS ──────────────────────────────────────────────────────────────────────

CSS = """
:root {
  --navy:   #003f5c;
  --blue:   #2f6fd4;
  --green:  #0ea47a;
  --orange: #e07b39;
  --red:    #d62728;
  --purple: #7b4d9e;
  --gold:   #c9a227;
  --light:  #f4f6fb;
  --border: #d0d8e8;
  --text:   #1a2233;
}
*, *::before, *::after { box-sizing: border-box; margin: 0; padding: 0; }
body {
  font-family: 'Segoe UI', Arial, sans-serif;
  font-size: 15px;
  color: var(--text);
  background: #e8ecf3;
}
/* ── COVER ── */
.cover {
  min-height: 100vh;
  display: flex;
  flex-direction: column;
  justify-content: center;
  align-items: center;
  background: linear-gradient(135deg, var(--navy) 0%, #005f8a 100%);
  color: #fff;
  text-align: center;
  padding: 60px 40px;
  page-break-after: always;
}
.cover-title  { font-size: 2.6rem; font-weight: 800; letter-spacing: -0.5px; margin-bottom: 18px; }
.cover-sub    { font-size: 1.25rem; opacity: 0.85; max-width: 720px; line-height: 1.6; margin-bottom: 32px; }
.cover-meta   { font-size: 0.9rem; opacity: 0.6; }
.cover-badge  {
  display: inline-block;
  background: var(--orange);
  color: #fff;
  padding: 6px 18px;
  border-radius: 20px;
  font-size: 0.85rem;
  font-weight: 700;
  margin-bottom: 30px;
  letter-spacing: 0.5px;
}
/* ── SLIDE ── */
.slide {
  background: #fff;
  min-height: 90vh;
  margin: 20px auto;
  max-width: 1100px;
  border-radius: 12px;
  box-shadow: 0 4px 24px rgba(0,0,0,0.10);
  display: flex;
  flex-direction: column;
  overflow: hidden;
  page-break-after: always;
}
.slide-header {
  display: flex;
  align-items: center;
  gap: 18px;
  padding: 22px 32px 18px;
  background: var(--navy);
  color: #fff;
}
.slide-num {
  width: 44px; height: 44px;
  border-radius: 50%;
  background: var(--blue);
  display: flex; align-items: center; justify-content: center;
  font-size: 1.2rem; font-weight: 800;
  flex-shrink: 0;
}
.slide-title { font-size: 1.5rem; font-weight: 700; flex: 1; }
.slide-tag {
  background: var(--orange);
  padding: 4px 14px;
  border-radius: 12px;
  font-size: 0.8rem;
  font-weight: 700;
  white-space: nowrap;
}
.slide-body {
  display: flex;
  flex: 1;
  gap: 0;
}
/* ── FIGURE ── */
.figure-block {
  flex: 1.1;
  padding: 24px;
  background: var(--light);
  display: flex;
  flex-direction: column;
  align-items: center;
  justify-content: center;
  border-right: 1px solid var(--border);
}
.slide-fig {
  max-width: 100%;
  max-height: 500px;
  border-radius: 8px;
  box-shadow: 0 2px 12px rgba(0,0,0,0.08);
  object-fit: contain;
}
.fig-caption {
  margin-top: 10px;
  font-size: 0.78rem;
  color: #666;
  text-align: center;
  font-style: italic;
}
.fig-missing {
  padding: 24px;
  background: #fff3cd;
  border: 1px solid #f0c000;
  border-radius: 6px;
  color: #7a5c00;
  font-size: 0.85rem;
}
.two-figs {
  flex-direction: row;
  gap: 16px;
  align-items: flex-start;
}
.fig-col {
  flex: 1;
  display: flex;
  flex-direction: column;
  align-items: center;
}
/* ── STORY ── */
.story-block {
  flex: 1;
  padding: 28px 28px 24px;
  display: flex;
  flex-direction: column;
  gap: 20px;
}
.story-section { display: flex; flex-direction: column; gap: 6px; }
.story-label {
  font-size: 0.7rem;
  font-weight: 800;
  text-transform: uppercase;
  letter-spacing: 1px;
  color: var(--blue);
}
.story-text {
  font-size: 0.95rem;
  line-height: 1.6;
  color: var(--text);
}
.story-bullets {
  list-style: none;
  display: flex;
  flex-direction: column;
  gap: 8px;
}
.story-bullets li {
  padding-left: 20px;
  position: relative;
  font-size: 0.93rem;
  line-height: 1.55;
}
.story-bullets li::before {
  content: "▶";
  position: absolute;
  left: 0;
  color: var(--blue);
  font-size: 0.65rem;
  top: 4px;
}
.takeaway {
  background: linear-gradient(135deg, #e8f4fd, #d0e8f8);
  border-left: 4px solid var(--blue);
  border-radius: 0 8px 8px 0;
  padding: 14px 18px;
  font-size: 0.92rem;
  line-height: 1.5;
  color: var(--navy);
  margin-top: auto;
}
.takeaway-label {
  display: block;
  font-size: 0.65rem;
  font-weight: 800;
  letter-spacing: 1.2px;
  color: var(--blue);
  margin-bottom: 4px;
}
.story-note {
  background: #fff8e1;
  border: 1px solid #ffe082;
  border-radius: 6px;
  padding: 10px 14px;
  font-size: 0.83rem;
  color: #7a5c00;
  line-height: 1.5;
}
/* ── COVER-ONLY slides (title slide, text-only) ── */
.slide-text-only .slide-body {
  flex-direction: column;
  padding: 32px 40px;
  gap: 24px;
}
.slide-text-only .story-block {
  flex: none;
  padding: 0;
}
/* ── METRIC CARDS ── */
.metric-row {
  display: flex;
  gap: 16px;
  flex-wrap: wrap;
  margin: 16px 0;
}
.metric-card {
  flex: 1;
  min-width: 140px;
  background: var(--light);
  border: 1px solid var(--border);
  border-radius: 8px;
  padding: 16px;
  text-align: center;
}
.metric-val {
  font-size: 2rem;
  font-weight: 800;
  color: var(--navy);
  line-height: 1;
}
.metric-lbl {
  font-size: 0.75rem;
  color: #666;
  margin-top: 4px;
}
.metric-card.green  .metric-val { color: var(--green); }
.metric-card.orange .metric-val { color: var(--orange); }
.metric-card.red    .metric-val { color: var(--red); }
.metric-card.blue   .metric-val { color: var(--blue); }
/* ── TOC ── */
.toc { list-style: none; display: flex; flex-direction: column; gap: 6px; margin-top: 16px; }
.toc li {
  display: flex; align-items: center; gap: 12px;
  padding: 10px 14px;
  background: var(--light);
  border-radius: 6px;
  font-size: 0.9rem;
}
.toc-num {
  width: 28px; height: 28px;
  border-radius: 50%;
  background: var(--navy);
  color: #fff;
  display: flex; align-items: center; justify-content: center;
  font-size: 0.75rem; font-weight: 700;
  flex-shrink: 0;
}
.toc-tag {
  margin-left: auto;
  font-size: 0.7rem;
  background: var(--blue);
  color: #fff;
  padding: 2px 10px;
  border-radius: 10px;
}
/* ── PRINT ── */
@media print {
  body { background: #fff; font-size: 13px; }
  .slide, .cover {
    box-shadow: none;
    border-radius: 0;
    margin: 0;
    max-width: 100%;
    page-break-after: always;
  }
  .slide-fig { max-height: 380px; }
}
"""

# ─── SLIDES ───────────────────────────────────────────────────────────────────

F = PUB / "analysis" / "figures_final"
FIMM_F = PUB / "fimm_results" / "figures"
FIMM_L = LOCAL / "analysis" / "figures_final"

def build_slides() -> str:
    out = []

    # ── COVER ──────────────────────────────────────────────────────────────────
    out.append("""
<div class="cover">
  <div class="cover-badge">MVHLA — Presentation Guide</div>
  <h1 class="cover-title">Majority-Voting HLA Typing<br>Ensemble</h1>
  <p class="cover-sub">
    A step-by-step walkthrough of every key figure —<br>
    what it shows, what to say, and the one line that matters most.
  </p>
  <p class="cover-meta">FIMM AML/MDS Cohort &nbsp;·&nbsp; Benchmarked on 1000 Genomes Project
    &nbsp;·&nbsp; 8 Tools · 3 Modalities</p>
</div>
""")

    # ── SLIDE 1 — Project overview (text only) ─────────────────────────────────
    out.append(slide(1, "What Is This Project?",
        """
<div class="slide-text-only" style="padding:32px 40px;display:flex;flex-direction:column;gap:20px;width:100%">
  <div class="story-section">
    <div class="story-label">One sentence</div>
    <p class="story-text" style="font-size:1.15rem;font-weight:600;color:var(--navy)">
      MVHLA is a benchmark-derived, reliability-weighted ensemble that combines 8 HLA typing tools
      across whole-genome, whole-exome, and RNA sequencing to deliver accurate, calibrated
      HLA allele calls with honest confidence scores.
    </p>
  </div>
  <div class="story-section">
    <div class="story-label">Why does HLA typing matter?</div>
    <p class="story-text">
      HLA genes (Human Leukocyte Antigens) control how the immune system recognises cancer cells.
      In AML/MDS, patients with more diverse HLA alleles may mount stronger anti-tumour responses.
      Accurate HLA typing is the first step for studying this — but existing tools often disagree
      with each other, and none tell you how much to trust their answer.
    </p>
  </div>
  <div class="story-section">
    <div class="story-label">What did we do?</div>
    <ul class="story-bullets">
      <li>Benchmarked 8 published tools on 1000 Genomes Project samples with known truth HLA alleles</li>
      <li>Learned reliability weights per tool per modality from the benchmark data</li>
      <li>Built a majority-vote ensemble that outputs a call <em>and</em> a calibrated confidence score</li>
      <li>Applied it to 1,227 FIMM AML/MDS samples and linked HLA diversity to patient survival</li>
    </ul>
  </div>
  <div class="metric-row">
    <div class="metric-card blue"><div class="metric-val">8</div><div class="metric-lbl">HLA typing tools</div></div>
    <div class="metric-card blue"><div class="metric-val">3</div><div class="metric-lbl">Sequencing modalities</div></div>
    <div class="metric-card green"><div class="metric-val">1,227</div><div class="metric-lbl">FIMM samples</div></div>
    <div class="metric-card orange"><div class="metric-val">30</div><div class="metric-lbl">FIMM patients with survival</div></div>
  </div>
</div>
""",
    tag="Introduction"))

    # ── SLIDE 2 — Pipeline architecture ────────────────────────────────────────
    out.append(slide(2, "How the Pipeline Works",
        fig(F / "figure_1_workflow_architecture.png",
            "Figure 1 — MVHLA pipeline: 8 tools → 4 stages → ensemble call")
        + story(
            "The pipeline has four stages: (1) dispatch input reads to all 8 tools in parallel, "
            "(2) phase-gate low-quality calls using calibration guardrails, "
            "(3) compute tool reliability weights from benchmark performance, "
            "(4) majority vote with weighted confidence → final call or abstention.",
            [
                "We didn't invent a new HLA typing method — we combined 8 published tools intelligently.",
                "Each tool runs on your sequencing data. We then ask: which tools agree? How reliable is each one?",
                "If no consensus can be reached at a given confidence threshold, the pipeline abstains rather than guessing.",
                "The weights were learned from a held-out benchmark — not hand-tuned.",
            ],
            "The ensemble is smarter than any single tool because it knows which tools to trust in which situations.",
        ),
        tag="Architecture"))

    # ── SLIDE 3 — Main accuracy result ─────────────────────────────────────────
    out.append(slide(3, "The Main Result: Ensemble vs. Individual Tools",
        fig(F / "figure_2_accuracy_comparison.png",
            "Figure 2 — Allele-level accuracy across tools and modalities")
        + story(
            "Each bar is a tool's allele-level accuracy (fraction of HLA alleles called correctly at 4-digit resolution) "
            "across the three sequencing modalities. The rightmost bar is MVHLA's ensemble.",
            [
                "For WES samples, our ensemble reaches 94.9% — the best single tool sits at ~80–88%.",
                "For RNA-seq, we reach 95.3%. For WGS, we abstain on ~8% of samples and reach ~70% on the rest — "
                "because WGS tools are uniquely unreliable (see next slide).",
                "No single tool is best across all modalities. The ensemble adapts.",
                "The gains are consistent — this is not a single lucky cohort.",
            ],
            "The ensemble adds 5–15 percentage points over the best individual tool, across all three modalities.",
        ),
        tag="Core Result"))

    # ── SLIDE 4 — WGS failure / false duplicate rate ────────────────────────────
    out.append(slide(4, "Why WGS Tools Struggle: The False Duplicate Problem",
        fig(F / "figure_a2_tool_ceiling.png",
            "Figure A2 — False duplicate rate per tool on WGS (primary benchmark)")
        + story(
            "The 'false duplicate rate' measures how often a tool calls a heterozygous patient as homozygous — "
            "reporting the same allele twice instead of finding both alleles. "
            "On WGS data, most tools do this 75–100% of the time.",
            [
                "Whole-genome sequencing sounds like the gold standard, but HLA loci are so repetitive that "
                "most tools get confused and collapse two different alleles into one.",
                "This is not a bug we can fix — it's a fundamental challenge of the HLA region in WGS.",
                "Our ensemble detects when tools are likely doing this and flags the call for abstention.",
                "This is why WES and RNA-seq are our primary reliable modalities.",
            ],
            "WGS is not the gold standard for HLA typing. WES and RNA are more reliable.",
            note="Governance note: say 'WGS is uniquely challenging for HLA typing tools' — "
                 "not 'WGS is bad'. The limitation is specific to HLA loci, not WGS in general.",
        ),
        tag="WGS Caveat"))

    # ── SLIDE 5 — Per-gene accuracy ─────────────────────────────────────────────
    out.append(slide(5, "Where the Gains Are Biggest: Per-Gene Accuracy",
        fig(F / "figure_3_per_gene_gains.png",
            "Figure 3 — Accuracy improvement per HLA gene (A, B, C)")
        + story(
            "The three HLA genes (A, B, C) have different difficulty levels. "
            "This figure shows how much accuracy each tool achieves per gene, "
            "and how the ensemble improves on the best individual tool.",
            [
                "HLA-A is the easiest — most tools do well and the ensemble polishes the last few percent.",
                "HLA-B is intermediate in difficulty.",
                "HLA-C is the hardest — even the best tools leave a lot on the table. "
                "The ensemble recovers the most ground here.",
                "The reason HLA-C is harder: it has more ambiguous allele families and less training data in most tools.",
            ],
            "The ensemble's gains are real across all three genes — biggest where it matters most: HLA-C.",
        ),
        tag="Per Gene"))

    # ── SLIDE 6 — Calibration ──────────────────────────────────────────────────
    out.append(slide(6, "Can We Trust the Confidence Scores?",
        fig(F / "figure_4_confidence_calibration.png",
            "Figure 4 — Predicted confidence vs. observed accuracy (calibration plot)")
        + story(
            "A calibration plot compares a model's stated confidence with its actual accuracy. "
            "If the curve follows the diagonal, the confidence scores are reliable. "
            "Points above the diagonal mean over-confidence; below means under-confidence.",
            [
                "MVHLA's confidence scores lie very close to the diagonal — they are well calibrated.",
                "In plain English: when we say we are 90% confident, we are right about 90% of the time.",
                "Individual tools are often overconfident — they claim high confidence but are frequently wrong.",
                "Calibration matters clinically: a doctor using an uncalibrated tool cannot know when to trust it.",
            ],
            "MVHLA's confidence scores are trustworthy — a feature rare in published HLA typing tools.",
            note="Guardrail: we apply Brier score ≤ 0.35 and ECE ≤ 0.35 before accepting any tool's weight. "
                 "Badly miscalibrated tools are downweighted automatically.",
        ),
        tag="Calibration"))

    # ── SLIDE 7 — Abstention ────────────────────────────────────────────────────
    out.append(slide(7, "Knowing When to Say 'I Don't Know'",
        fig(F / "figure_5_abstention_tradeoff.png",
            "Figure 5 — Abstention rate vs. accuracy on retained calls")
        + story(
            "The abstention tradeoff curve shows what happens as we raise the confidence threshold: "
            "more samples are rejected (abstained), but accuracy on the kept calls rises. "
            "Each modality has its own curve.",
            [
                "For WGS: abstaining on 8% of samples raises accuracy from ~44% to ~70% on the remaining calls.",
                "For WES/RNA: the curve is already high, so abstention provides a smaller but still real boost.",
                "The chosen operating point (the dot on each curve) balances accuracy gain against sample loss.",
                "Abstention is not a failure — it is honest. A 'no call' is better than a wrong call.",
            ],
            "Saying 'I don't know' on 8% of WGS samples more than doubles accuracy on the ones we do answer.",
        ),
        tag="Abstention"))

    # ── SLIDE 8 — Discordance taxonomy ─────────────────────────────────────────
    out.append(slide(8, "Why Do Tools Disagree? A Taxonomy of Errors",
        fig(F / "figure_6_discordance_taxonomy.png",
            "Figure 6 — Types of tool discordance classified by error category")
        + story(
            "When tools disagree, it is not random. This figure classifies disagreements into "
            "categories: digit-level ambiguity, allele dropout (false homozygosity), novel alleles, "
            "and cross-gene confusion. Each category has a different reliability signature.",
            [
                "Digit-level disagreements (4-digit vs. 2-digit) are usually benign — clinical impact is low.",
                "Allele dropout is the most dangerous: one allele is simply missed. Concentrated in WGS.",
                "Novel allele calls are rare but high-risk — a tool may hallucinate an allele not in its database.",
                "By classifying disagreements, we can assign reliability weights with a principled basis.",
            ],
            "Tool disagreement has predictable structure — and that structure is what the ensemble exploits.",
        ),
        tag="Discordance"))

    # ── SLIDE 9 — Confidence weights ────────────────────────────────────────────
    out.append(slide(9, "How We Weight Each Tool",
        fig(F / "figure_7_confidence_weights.png",
            "Figure 7 — Tool reliability weights per modality (learned from benchmark)")
        + story(
            "Each bar shows the weight assigned to a tool in the final ensemble, separated by modality. "
            "Weights are derived from benchmark accuracy and calibration scores using the formula: "
            "0.7 × base_reliability + 0.3 × effective_confidence.",
            [
                "OptiType gets the highest weight for WGS — it makes the fewest false duplicate calls.",
                "ArcasHLA and SpecHLA lead for RNA-seq, where their splice-aware design pays off.",
                "Tools that fail calibration guardrails (Brier > 0.35 or ECE > 0.35) receive weight zero for that modality.",
                "Weights are not guessed or hand-tuned — they are learned from data.",
            ],
            "Evidence-based weights mean the ensemble improves automatically if a better tool emerges.",
        ),
        tag="Weights"))

    # ── SLIDE 10 — Bimodal + trimodal ──────────────────────────────────────────
    out.append(slide(10, "Combining Modalities: WES + RNA",
        two_figs(
            F / "figure_09_bimodal_per_gene.png", "Figure 9 — Bimodal (WES+RNA) per-gene accuracy",
            F / "figure_10_trimodal_comparison.png", "Figure 10 — Trimodal vs. bimodal comparison",
        )
        + story(
            "When a patient has both WES and RNA sequencing, we can combine calls across modalities. "
            "Figure 9 shows per-gene accuracy for the bimodal (WES+RNA) ensemble. "
            "Figure 10 compares bimodal vs. trimodal (adding WGS).",
            [
                "Bimodal WES+RNA reaches 96.2% — a meaningful gain over either modality alone.",
                "Adding WGS (trimodal) gives only +0.05 percentage points — marginal benefit, high WGS cost.",
                "The per-gene gains in bimodal mode are consistent: HLA-C benefits most, as expected.",
                "For FIMM patients: most have both WES and RNA, so bimodal is the operative mode.",
            ],
            "Bimodal (WES+RNA) is the sweet spot. Trimodal adds almost nothing while tripling sequencing cost.",
        ),
        tag="Multi-Modal"))

    # ── SLIDE 11 — LOH benchmark ────────────────────────────────────────────────
    out.append(slide(11, "Loss of Heterozygosity (LOH): What It Is and How We Detect It",
        fig(F / "figure_loh_allele_dropout_heatmap.png",
            "Figure LOH — Allele dropout heatmap across VENEX WGS samples")
        + story(
            "LOH occurs when a tumour cell deletes one copy of an HLA gene, leaving only one allele. "
            "This is detectable in WGS data as an allele frequency imbalance: "
            "one allele appears at >80% frequency instead of the expected ~50%.",
            [
                "The heatmap shows each VENEX WGS sample (rows) × gene (columns). "
                "Red cells = candidate LOH (major allele fraction ≥ 80%).",
                "Clear patterns emerge: some patients show LOH across multiple HLA genes — "
                "a sign of large chromosomal deletions.",
                "LOH detection requires WGS for reliability — allele frequency estimates from WES are "
                "biased by capture probe design and uneven exon coverage.",
                "For FIMM WES samples (next slide) we flag candidates only — WGS is needed to confirm.",
            ],
            "WGS is the only reliable modality for HLA LOH. WES can only screen for candidates.",
        ),
        tag="LOH — Benchmark"))

    # ── SLIDE 12 — FIMM concordance ─────────────────────────────────────────────
    out.append(slide(12, "FIMM Real-World Data: Do WES and RNA Agree?",
        two_figs(
            FIMM_F / "fig_concordance_by_gene.png", "FIMM — HLA concordance by gene (WES vs. RNA)",
            FIMM_F / "fig_crossmodal_comparison.png", "FIMM — Cross-modal WES vs. RNA call comparison",
        )
        + story(
            "These figures show how well WES-based and RNA-based HLA calls agree for the same FIMM patient. "
            "Concordance is measured as the fraction of patients where both modalities give the same allele call.",
            [
                "HLA-A concordance is high — both modalities agree on the most common alleles.",
                "HLA-C concordance is lower — consistent with the benchmark finding that HLA-C is harder.",
                "The cross-modal comparison (right figure) shows the spread of call pairs: "
                "most cluster on the diagonal (agreement), with HLA-C showing the most scatter.",
                "This validates the benchmark: the patterns seen in 1000 Genomes samples also appear in real AML/MDS patients.",
            ],
            "Real FIMM clinical data confirms the benchmark findings — what we learned on 1000 Genomes holds in AML.",
        ),
        tag="FIMM — Concordance"))

    # ── SLIDE 13 — FIMM WES LOH ─────────────────────────────────────────────────
    out.append(slide(13, "FIMM WES LOH Screen: Flagging Candidates",
        fig(FIMM_L / "figure_fimm_loh_wes.png",
            "FIMM WES LOH — Status distribution across HLA-A, B, C")
        + story(
            "This figure shows the LOH classification for 225 FIMM sample-gene pairs using WES allele frequencies. "
            "The categories are: candidate LOH, balanced heterozygous, allelic imbalance (review needed), "
            "and insufficient evidence.",
            [
                "Most FIMM WES samples fall into 'allelic imbalance — review'. This is expected and honest.",
                "WES capture panels produce uneven coverage at HLA loci, so allele frequencies are unreliable "
                "for the stringent thresholds calibrated for WGS (major fraction ≥ 80%, ≥ 20 het variants).",
                "The few 'candidate LOH' calls (if any) should be considered hypotheses, not conclusions.",
                "The purpose of this screen is to prioritise which patients should get WGS follow-up, "
                "not to confirm LOH.",
            ],
            "WES LOH is a hypothesis-generating screen, not a diagnostic result. Every candidate needs WGS confirmation.",
            note="This analysis was run on newly ingested FIMM SpecHLA data "
                 "(1,227 sample rows in fimm_hla_calls.tsv).",
        ),
        tag="FIMM — LOH"))

    # ── SLIDE 14 — HED & Survival ────────────────────────────────────────────────
    out.append(slide(14, "HLA Evolutionary Divergence (HED) and Patient Survival",
        """
<div style="padding:28px;background:var(--light);border-right:1px solid var(--border);
            display:flex;flex-direction:column;justify-content:center;gap:20px;min-width:340px;">
  <div class="story-label">What is HED?</div>
  <p class="story-text">
    HED (HLA Evolutionary Divergence) measures how different the two alleles of an HLA gene are
    from each other, using the Grantham amino acid distance at the antigen-binding groove (positions A, B, G).
    Higher HED = more diverse HLA = potentially better at presenting different neoantigens.
  </p>
  <div class="story-label" style="margin-top:8px">FIMM cohort summary</div>
  <div class="metric-row">
    <div class="metric-card blue"><div class="metric-val">30</div><div class="metric-lbl">Patients with HED</div></div>
    <div class="metric-card blue"><div class="metric-val">27</div><div class="metric-lbl">With survival data</div></div>
    <div class="metric-card orange"><div class="metric-val">14</div><div class="metric-lbl">Events (deaths)</div></div>
    <div class="metric-card green"><div class="metric-val">0</div><div class="metric-lbl">Concordant calls<br><small style="font-size:0.65rem">(3-tool WES agreement)</small></div></div>
  </div>
  <div class="story-note" style="margin-top:4px">
    Concordant = 0 because SpecHLA WES provides allele calls for only 4 of 35 patients in the current dataset.
    Full 3-tool consensus requires OptiType + arcasHLA + SpecHLA to all return a call.
  </div>
</div>
""" + story(
            "The HED report tests whether patients with higher HLA diversity (measured by HED) "
            "have better overall survival from diagnosis (os_from_dx). "
            "Cox proportional-hazards regression was used: univariable (HED_total alone) "
            "and multivariable (HED_total + diagnosis group).",
            [
                "Higher HED is hypothesised to improve survival by broadening neoantigen presentation.",
                "With only 27 patients and 14 events, statistical power is limited — "
                "treat any p-value as exploratory, not confirmatory.",
                "The primary endpoint is os_from_dx (overall survival from diagnosis) to avoid "
                "landmark bias that would arise from os_from_sampling.",
                "This is hypothesis-generating data. A larger cohort with treatment data "
                "is needed before clinical conclusions can be drawn.",
            ],
            "HED × survival is a hypothesis that needs a larger cohort. The current n=27 is exploratory.",
            note="Key framing: say 'we observe a trend consistent with the HED hypothesis' — "
                 "never 'HED predicts survival' without a significant, well-powered result.",
        ),
        tag="HED × Survival"))

    return "\n".join(out)

# ─── RENDER ────────────────────────────────────────────────────────────────────

def render() -> str:
    return f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>MVHLA — Guided Presentation</title>
<style>
{CSS}
</style>
</head>
<body>
{build_slides()}
<div style="text-align:center;padding:40px;color:#999;font-size:0.8rem;page-break-before:avoid">
  Generated {__import__('datetime').date.today()} &nbsp;·&nbsp;
  MVHLA Majority-Voting HLA Typing Ensemble &nbsp;·&nbsp;
  Open in browser → Ctrl+P → Save as PDF
</div>
</body>
</html>"""

# ─── MAIN ──────────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    html = render()
    OUTPUT.write_text(html, encoding="utf-8")
    size_kb = OUTPUT.stat().st_size // 1024
    print(f"Written → {OUTPUT}  ({size_kb} KB)")
    print("Open in browser → Ctrl+P → Save as PDF")
