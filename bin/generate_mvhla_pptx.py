#!/usr/bin/env python3.11
"""Generate MVHLA presentation PPTX with all current publication figures."""

from pptx import Presentation
from pptx.util import Inches, Pt, Emu
from pptx.dml.color import RGBColor
from pptx.enum.text import PP_ALIGN
from pptx.util import Inches, Pt
import os
from datetime import date

FIGURES_DIR = "/scratch/project_2008084/pihla-publish/analysis/figures_final_candidate"
OUT_PATH = "/scratch/project_2008084/pihla-publish/docs/MVHLA_Presentation_2026-06-04.pptx"

# Color palette (Wong 2011 colorblind-safe)
C_BLUE   = RGBColor(0, 114, 178)
C_ORANGE = RGBColor(230, 159, 0)
C_GREEN  = RGBColor(0, 158, 115)
C_RED    = RGBColor(213, 94, 0)
C_SKY    = RGBColor(86, 180, 233)
C_WHITE  = RGBColor(255, 255, 255)
C_DARK   = RGBColor(30, 30, 40)
C_LGRAY  = RGBColor(245, 246, 250)
C_MGRAY  = RGBColor(180, 185, 200)

SLIDE_W = Inches(13.33)
SLIDE_H = Inches(7.5)

prs = Presentation()
prs.slide_width  = SLIDE_W
prs.slide_height = SLIDE_H

BLANK = prs.slide_layouts[6]   # completely blank

def rgb(r, g, b):
    return RGBColor(r, g, b)

def add_rect(slide, left, top, width, height, fill_color, line_color=None):
    shape = slide.shapes.add_shape(1, left, top, width, height)
    shape.fill.solid()
    shape.fill.fore_color.rgb = fill_color
    if line_color:
        shape.line.color.rgb = line_color
        shape.line.width = Pt(0.75)
    else:
        shape.line.fill.background()
    return shape

def add_text(slide, text, left, top, width, height,
             font_size=18, bold=False, color=C_DARK,
             align=PP_ALIGN.LEFT, wrap=True, italic=False):
    txBox = slide.shapes.add_textbox(left, top, width, height)
    tf = txBox.text_frame
    tf.word_wrap = wrap
    p = tf.paragraphs[0]
    p.alignment = align
    run = p.add_run()
    run.text = text
    run.font.size = Pt(font_size)
    run.font.bold = bold
    run.font.italic = italic
    run.font.color.rgb = color
    return txBox

def header_bar(slide, title, subtitle=None):
    """Blue header bar at top."""
    add_rect(slide, 0, 0, SLIDE_W, Inches(1.1), C_BLUE)
    add_text(slide, title,
             Inches(0.35), Inches(0.08), Inches(12.5), Inches(0.6),
             font_size=28, bold=True, color=C_WHITE)
    if subtitle:
        add_text(slide, subtitle,
                 Inches(0.35), Inches(0.68), Inches(12.5), Inches(0.4),
                 font_size=14, bold=False, color=C_SKY)

def footer_bar(slide, text="mvHLA  |  Majority-Voting HLA Ensemble  |  2026"):
    add_rect(slide, 0, Inches(7.2), SLIDE_W, Inches(0.3), C_DARK)
    add_text(slide, text, Inches(0.3), Inches(7.21), Inches(12.5), Inches(0.25),
             font_size=9, color=C_MGRAY, align=PP_ALIGN.LEFT)

def add_figure_slide(fig_png, title, subtitle, caption_lines, caption_note=None):
    slide = prs.slides.add_slide(BLANK)
    add_rect(slide, 0, 0, SLIDE_W, SLIDE_H, C_LGRAY)
    header_bar(slide, title, subtitle)
    footer_bar(slide)

    # Figure image — left panel
    fig_left   = Inches(0.25)
    fig_top    = Inches(1.2)
    fig_width  = Inches(8.5)
    fig_height = Inches(5.7)
    if os.path.exists(fig_png):
        slide.shapes.add_picture(fig_png, fig_left, fig_top, fig_width, fig_height)
    else:
        add_rect(slide, fig_left, fig_top, fig_width, fig_height, C_MGRAY)
        add_text(slide, f"[figure not found]\n{fig_png}",
                 fig_left + Inches(0.1), fig_top + Inches(0.1),
                 fig_width - Inches(0.2), fig_height - Inches(0.2),
                 font_size=11, color=C_DARK)

    # Caption panel — right side
    cap_left  = Inches(9.0)
    cap_top   = Inches(1.2)
    cap_width = Inches(4.1)

    add_rect(slide, cap_left - Inches(0.15), cap_top,
             cap_width + Inches(0.3), Inches(5.7),
             C_WHITE, line_color=C_MGRAY)

    y = cap_top + Inches(0.15)
    for line in caption_lines:
        bold_flag = line.startswith("**") and line.endswith("**")
        txt = line.strip("*") if bold_flag else line
        add_text(slide, txt,
                 cap_left, y, cap_width, Inches(0.5),
                 font_size=11, bold=bold_flag, color=C_DARK, wrap=True)
        y += Inches(0.42) if bold_flag else Inches(0.36)

    if caption_note:
        add_rect(slide, cap_left - Inches(0.1), y + Inches(0.05),
                 cap_width + Inches(0.2), Inches(0.6), C_LGRAY)
        add_text(slide, caption_note,
                 cap_left, y + Inches(0.1), cap_width, Inches(0.6),
                 font_size=9, italic=True, color=C_MGRAY, wrap=True)


# ── SLIDE 1: Title ───────────────────────────────────────────────────────────
s1 = prs.slides.add_slide(BLANK)
add_rect(s1, 0, 0, SLIDE_W, SLIDE_H, C_DARK)
add_rect(s1, 0, Inches(2.4), SLIDE_W, Inches(2.8), C_BLUE)

add_text(s1, "mvHLA",
         Inches(0.6), Inches(0.5), Inches(12), Inches(1.0),
         font_size=54, bold=True, color=C_WHITE, align=PP_ALIGN.LEFT)

add_text(s1, "Majority-Voting HLA Ensemble",
         Inches(0.6), Inches(1.45), Inches(12), Inches(0.6),
         font_size=24, bold=False, color=C_SKY, align=PP_ALIGN.LEFT)

add_text(s1,
         "A benchmark-derived reliability-weighted multi-tool ensemble for HLA typing\n"
         "with calibration-gated confidence integration across WGS, WES, and RNA-seq",
         Inches(0.6), Inches(2.55), Inches(12), Inches(1.0),
         font_size=16, bold=False, color=C_WHITE, align=PP_ALIGN.LEFT, wrap=True)

add_text(s1,
         "HLA-A  ·  HLA-B  ·  HLA-C  |  8 tools  |  WGS · WES · RNA-seq  |  1000 Genomes truth-backed cohort",
         Inches(0.6), Inches(3.65), Inches(12), Inches(0.5),
         font_size=13, bold=False, color=C_ORANGE, align=PP_ALIGN.LEFT)

add_text(s1,
         "Target venue: Bioinformatics\n"
         f"Presentation date: {date.today().strftime('%B %d, %Y')}",
         Inches(0.6), Inches(6.5), Inches(8), Inches(0.7),
         font_size=12, color=C_MGRAY, align=PP_ALIGN.LEFT)


# ── SLIDE 2: Scope & Cohort ───────────────────────────────────────────────────
s2 = prs.slides.add_slide(BLANK)
add_rect(s2, 0, 0, SLIDE_W, SLIDE_H, C_LGRAY)
header_bar(s2, "Study Scope & Benchmark Cohort",
           "Eight HLA typing tools · Three sequencing modalities · Truth-backed 1000 Genomes samples")
footer_bar(s2)

bullets = [
    ("Platform", "Nextflow DSL2, Docker/Singularity, CSC Puhti HPC"),
    ("Tools integrated", "OptiType, ArcasHLA, SpecHLA, HLA-HD, POLYSOLVER, Kourami, T1K, Seq2HLA"),
    ("Loci reported", "HLA-A, HLA-B, HLA-C (two-field resolution; IMGT/HLA 3.59.0)"),
    ("Primary cohort", "WGS  n=50  |  splits: 28 train / 10 val / 12 holdout  |  5 population groups"),
    ("Secondary cohorts", "WES  n=129  ·  RNA-seq  n=107  ·  Cross-modal (WES+RNA)  n=106"),
    ("Truth source", "1000 Genomes HLA dataset (Gourraud et al. 2014)  —  public reference tier"),
    ("Weight formula", "final_weight = 0.7 × base_reliability + 0.3 × effective_confidence"),
    ("Guardrail threshold", "ECE or Brier > 0.35 → poor_calibration → confidence boost blocked"),
]

col_x = [Inches(0.4), Inches(3.1)]
y0 = Inches(1.25)
for i, (label, val) in enumerate(bullets):
    y = y0 + i * Inches(0.68)
    add_rect(s2, col_x[0], y, Inches(2.6), Inches(0.52), C_BLUE)
    add_text(s2, label, col_x[0] + Inches(0.1), y + Inches(0.08),
             Inches(2.4), Inches(0.4), font_size=12, bold=True, color=C_WHITE)
    add_rect(s2, col_x[1], y, Inches(9.8), Inches(0.52), C_WHITE, line_color=C_MGRAY)
    add_text(s2, val, col_x[1] + Inches(0.1), y + Inches(0.08),
             Inches(9.6), Inches(0.4), font_size=12, color=C_DARK)


# ── SLIDE 3: Key Results Summary ─────────────────────────────────────────────
s3 = prs.slides.add_slide(BLANK)
add_rect(s3, 0, 0, SLIDE_W, SLIDE_H, C_LGRAY)
header_bar(s3, "Key Quantitative Results",
           "WGS is tool-limited; WES and RNA-seq approach ensemble ceiling via routing")
footer_bar(s3)

panels = [
    ("WGS  (n=50, holdout n=12)", C_BLUE, [
        "OptiType (best single tool):  0.4722  callable 1.0",
        "WeightedConsensus:            0.4444  callable 0.917",
        "  accuracy-among-callable:    0.4848",
        "MajorityVote:                 0.3889",
        "ChampionChallenger:           0.4722  (recovers ceiling)",
        "",
        "Bottleneck: HLA-C (0.17 for all methods)",
        "HLA-B easiest (0.83); HLA-A intermediate (0.42)",
    ]),
    ("WES  (n=129)", C_GREEN, [
        "MajorityVote:                 0.9359  callable 1.0",
        "WeightedConsensus (tuned):    0.9359  callable 0.985",
        "ChampionChallenger (tuned):   0.9487  callable 1.0",
        "",
        "Inter-tool agreement high → ensemble near ceiling",
    ]),
    ("RNA-seq  (n=107)", C_ORANGE, [
        "MajorityVote:                 0.9502  callable 1.0",
        "WeightedConsensus (tuned):    0.9502  callable 0.994",
        "ChampionChallenger (tuned):   0.9533  callable 1.0",
        "",
        "Abstention rate: 4% RNA · 7% WES · 8% WGS",
        "ArcasHLA: 91.9% RNA vs 8.0% WGS — modality matters",
    ]),
]

for i, (title, color, lines) in enumerate(panels):
    px = Inches(0.3) + i * Inches(4.35)
    py = Inches(1.2)
    pw = Inches(4.15)
    ph = Inches(5.7)
    add_rect(s3, px, py, pw, ph, C_WHITE, line_color=C_MGRAY)
    add_rect(s3, px, py, pw, Inches(0.42), color)
    add_text(s3, title, px + Inches(0.1), py + Inches(0.06),
             pw - Inches(0.2), Inches(0.35),
             font_size=12, bold=True, color=C_WHITE)
    ty = py + Inches(0.52)
    for line in lines:
        mono = line.startswith("  ")
        add_text(s3, line.strip(), px + Inches(0.12), ty,
                 pw - Inches(0.24), Inches(0.38),
                 font_size=10, color=C_MGRAY if line == "" else C_DARK,
                 bold=False, italic=mono)
        ty += Inches(0.38)


# ── FIGURE SLIDES ────────────────────────────────────────────────────────────
FIGURES = [
    (
        "figure_1_workflow_architecture.png",
        "Figure 1 — MVHLA Workflow Architecture",
        "Methods and governance: benchmark-governed calibrated HLA ensemble framework",
        [
            "**Main message**",
            "MVHLA is a Nextflow DSL2 platform that harmonises 8 HLA typing tools into",
            "a unified benchmark-governed calibrated ensemble.",
            "",
            "**Benchmark role:** Methods / governance",
            "",
            "Inputs: BAM · CRAM · FASTQ",
            "Outputs: consensus calls, abstention flags,",
            "discordance tags, weight artifacts",
            "",
            "Supports Docker, Singularity, SLURM,",
            "and CSC Puhti profiles.",
        ],
        "Benchmark hierarchy is explicit in the architecture view."
    ),
    (
        "figure_2_accuracy_comparison.png",
        "Figure 2 — Primary WGS Accuracy Comparison",
        "WGS holdout (n=12): OptiType ceiling; ensemble recovers but does not exceed",
        [
            "**Main message**",
            "WeightedConsensus (0.4444 overall) does not exceed",
            "OptiType (0.4722) on the WGS holdout.",
            "ChampionChallenger recovers the ceiling (0.4722).",
            "",
            "**Benchmark role:** Primary",
            "Source: benchmark_wgs_wave2",
            "",
            "Tool spread is extreme (ArcasHLA 0.0 → OptiType 0.4722).",
            "MajorityVote falls to 0.3889 — low-performing tools",
            "outvote the strongest caller on conflicted loci.",
            "",
            "WGS bottleneck = tool landscape, not ensemble formula.",
        ],
        "Routed baselines (GatedConsensus, locus-expert panel) also recover 0.4722."
    ),
    (
        "figure_3_per_gene_gains.png",
        "Figure 3 — Per-Gene WGS Failure-Mode Interpretation",
        "HLA-B easiest · HLA-A intermediate · HLA-C hardest (gene-dense MHC region)",
        [
            "**Main message**",
            "Per-gene holdout accuracy varies markedly across class I loci.",
            "",
            "HLA-B:  0.8333  (OptiType & ChampionChallenger)",
            "HLA-A:  0.4167  (OptiType, T1K, ChampionChallenger)",
            "HLA-C:  0.1667  (all leading WGS methods)",
            "",
            "**Benchmark role:** Primary (WGS per-gene tables)",
            "",
            "HLA-C sits in the most gene-dense MHC segment",
            "with highest pseudogene sequence similarity.",
            "",
            "Also shows: WES/RNA operate in stronger",
            "unimodal regimes than WGS (Fig 5 context).",
        ],
        "Cross-modality context for WGS failure-mode story."
    ),
    (
        "figure_4_confidence_calibration.png",
        "Figure 4 — Confidence Calibration & Guardrails",
        "Raw tool confidence cannot be used naively — calibration gating is required",
        [
            "**Main message**",
            "Several tools are strongly overconfident relative to empirical truth.",
            "Calibration checks block naive confidence amplification.",
            "",
            "**Benchmark role:** Primary",
            "",
            "Poor calibration (ECE/Brier > 0.35) → blocked:",
            "  OptiType:  Brier 0.2465 / ECE 0.0027  applied",
            "  HLA-HD:    Brier 0.1673 / ECE 0.0348  blocked",
            "  Kourami:   Brier 0.1174 / ECE 0.0242  blocked",
            "  T1K:       Brier 0.309  / ECE 0.296   partial boost",
            "  ArcasHLA:  Brier 0.079  / ECE 0.076   min WGS weight",
            "",
            "SpecHLA: no_confidence → base-reliability only.",
        ],
        "Confidence expansion alone was not enough; guardrail layer was required."
    ),
    (
        "figure_5_abstention_tradeoff.png",
        "Figure 5 — Abstention–Accuracy Tradeoff (Cross-Modality)",
        "Abstention is an adaptive diagnostic signal, not a failure",
        [
            "**Main message**",
            "WeightedConsensus abstains selectively when the ensemble",
            "is genuinely conflicted. Abstention rates scale with",
            "ensemble difficulty.",
            "",
            "**Benchmark role:** Mixed (WGS primary; WES/RNA secondary)",
            "",
            "WGS holdout callable rate:  0.9167",
            "WES callable rate (tuned):  0.9846",
            "RNA callable rate (tuned):  0.9938",
            "",
            "Abstention rates:  8% WGS · 7% WES · 4% RNA",
            "",
            "Clinical interpretation: abstained locus → orthogonal",
            "confirmation (Sanger / PCR-SSO); committed call",
            "carries 48% among-callable accuracy on WGS.",
        ],
        "Abstention mechanism behaves adaptively — activates on genuine conflict."
    ),
    (
        "figure_6_discordance_taxonomy.png",
        "Figure 6 — Discordance Taxonomy",
        "Five-category structured disagreement labelling across modalities",
        [
            "**Main message**",
            "Disagreement patterns are structurally distinct per modality.",
            "",
            "**Benchmark role:** Primary (WGS discordance summaries)",
            "",
            "WGS:  low_evidence_conflict dominant (n=157)",
            "      → wide tool-accuracy spread",
            "WES:  low-evidence conflicts rare (n=26)",
            "RNA:  possible_expression_bias prominent (n=292)",
            "      → allele-specific expression effects",
            "",
            "Five categories:",
            "  technical_conflict",
            "  low_evidence_conflict",
            "  dna_rna_discordance",
            "  possible_expression_bias",
            "  no_evidence",
        ],
        "WGS false-duplicate rates 61–91% vs WES <31% → allele dropout artefact."
    ),
    (
        "figure_7_confidence_weights.png",
        "Figure 7 — Benchmark-Derived Confidence Weights",
        "Final runtime weights: OptiType dominant; T1K partial boost; HLA-HD/Kourami clipped",
        [
            "**Main message**",
            "Final weights preserve reliability while clipping poorly",
            "calibrated confidence.",
            "",
            "**Benchmark role:** Primary (WGS weight tables)",
            "",
            "OptiType:  0.5476  (dominant reliability)",
            "T1K:       0.2649  (partial confidence boost retained)",
            "HLA-HD:    0.2294  (clipped to base-reliability)",
            "Kourami:   clipped",
            "SpecHLA:   base-reliability only (no_confidence)",
            "ArcasHLA:  minimal WGS weight",
            "",
            "Without guardrails: broader confidence ingestion",
            "worsened WeightedConsensus performance.",
        ],
        "Guardrail is not merely theoretical — it materially improves ensemble behaviour."
    ),
]

for fig_file, title, subtitle, captions, note in FIGURES:
    fig_path = os.path.join(FIGURES_DIR, fig_file)
    add_figure_slide(fig_path, title, subtitle, captions, note)


# ── Supplementary Section Divider ────────────────────────────────────────────
sdiv = prs.slides.add_slide(BLANK)
add_rect(sdiv, 0, 0, SLIDE_W, SLIDE_H, C_DARK)
add_rect(sdiv, 0, Inches(2.8), SLIDE_W, Inches(0.08), C_ORANGE)
add_text(sdiv, "Supplementary Results",
         Inches(0.8), Inches(3.0), Inches(12), Inches(1.0),
         font_size=48, bold=True, color=C_WHITE, align=PP_ALIGN.CENTER)
add_text(sdiv, "Bimodal & Trimodal Cross-Modal Robustness  ·  Per-Gene Detail  ·  Calibration Heatmaps  ·  Resolution Comparison",
         Inches(0.8), Inches(4.1), Inches(12), Inches(0.5),
         font_size=15, bold=False, color=C_SKY, align=PP_ALIGN.CENTER)


# ── Supplementary Figure Slides ──────────────────────────────────────────────
SUPP_FIGURES = [
    (
        "figure_09_bimodal_per_gene.png",
        "Figure 9 — Bimodal Robustness: Unimodal + Bimodal Half",
        "Matched-subject WES+RNA: BimodalMajorityVote is strongest multimodal result",
        [
            "**Main message**",
            "BimodalMajorityVote on WES+RNA is the strongest",
            "multimodal result in the matched-subject robustness benchmark.",
            "",
            "**Benchmark role:** Supplementary",
            "Source: benchmark_trimodal_robustness/run",
            "Matched subjects n=106 (WES+RNA both available)",
            "",
            "Headline loci: HLA-A, -B, -C",
            "DRB1 / DQB1 present with partial truth",
            "(n=38/35) — excluded from headline claims.",
            "",
            "VENEX cohort: internal QC only,",
            "no truth labels — excluded from all figures.",
        ],
        "Interpret together with Figure 10 (trimodal half)."
    ),
    (
        "figure_10_trimodal_comparison.png",
        "Figure 10 — Trimodal Robustness Comparison",
        "Adding WGS to WES+RNA does not improve over best bimodal result",
        [
            "**Main message**",
            "TrimodalMajorityVote does not improve over",
            "BimodalMajorityVote (WES+RNA) in matched-subject analysis.",
            "",
            "**Benchmark role:** Supplementary",
            "Source: benchmark_trimodal_robustness/run",
            "",
            "WGS tool landscape is too heterogeneous to",
            "contribute net benefit when pooled with",
            "high-agreement WES+RNA tools.",
            "",
            "Conclusion: cross-modal gain is already captured",
            "by WES+RNA; WGS adds noise not signal at",
            "current cohort size and tool quality.",
        ],
        "NOTE: benchmark_trimodal_all_samples/ WGS n=136 (older) — do NOT use for single-modality claims."
    ),
    (
        "figure_s1_per_gene_accuracy.png",
        "Figure S1 — Per-Gene Accuracy (All Tools, All Modalities)",
        "Supplementary detail retained for publication continuity",
        [
            "**Main message**",
            "Per-gene accuracy across all integrated tools",
            "and modalities at two-field resolution.",
            "",
            "**Benchmark role:** Supplementary supporting detail",
            "",
            "Supports modality-appropriate tool selection argument.",
            "",
            "Key observations:",
            "  ArcasHLA: 91.9% RNA vs 8.0% WGS",
            "  SpecHLA:  high WGS accuracy, 3.5% RNA",
            "",
            "Tool selection must be modality-driven,",
            "not availability-driven.",
        ],
        None
    ),
    (
        "figure_s2_calibration_heatmap.png",
        "Figure S2 — Calibration Heatmap (ECE & Brier Score)",
        "Per-tool, per-modality calibration quality across the full benchmark cohort",
        [
            "**Main message**",
            "Calibration quality varies widely across tools",
            "and modalities — supporting the need for",
            "per-tool guardrail decisions.",
            "",
            "**Benchmark role:** Supplementary supporting detail",
            "",
            "ECE / Brier thresholds:  0.35 cutoff",
            "",
            "Tools passing WGS calibration guard:",
            "  T1K: Brier 0.309 / ECE 0.296 ✓",
            "  ArcasHLA: Brier 0.079 / ECE 0.076 ✓",
            "     (but near-zero WGS base-reliability)",
            "",
            "IMGT/HLA 3.59.0 pinned across all runs.",
        ],
        None
    ),
    (
        "figure_s3_resolution_comparison.png",
        "Figure S3 — Resolution Comparison (2-Field vs 3-Field)",
        "Primary endpoint: 2-field exact pair match; secondary: 3-field and G/P group",
        [
            "**Main message**",
            "Primary accuracy endpoint is exact two-field allele-pair",
            "match rate. Three-field and G/P-group rates are secondary.",
            "",
            "**Benchmark role:** Supplementary supporting detail",
            "",
            "A call is counted correct ONLY when both alleles",
            "of the returned pair match ground truth exactly",
            "under IMGT/HLA 3.59.0.",
            "",
            "No partial credit for heterozygous loci.",
            "Secondary metrics blank when group-level",
            "encodings absent from truth or call.",
        ],
        None
    ),
]

for fig_file, title, subtitle, captions, note in SUPP_FIGURES:
    fig_path = os.path.join(FIGURES_DIR, fig_file)
    add_figure_slide(fig_path, title, subtitle, captions, note)


# ── FINAL SLIDE: Conclusions ─────────────────────────────────────────────────
sc = prs.slides.add_slide(BLANK)
add_rect(sc, 0, 0, SLIDE_W, SLIDE_H, C_LGRAY)
header_bar(sc, "Conclusions & Next Steps",
           "Solid: WGS wave2 + WES/RNA results · Remaining: supplementary assembly before submission")
footer_bar(sc)

conclusions = [
    ("WGS is tool-limited, not ensemble-limited",
     "WeightedConsensus recovers to OptiType ceiling (0.4722) via routing; cannot exceed it. "
     "HLA-C is the hardest locus (0.167) due to MHC gene-density and pseudogene similarity."),
    ("Calibration guardrail is essential",
     "Three of five tools with confidence scores fail WGS calibration checks. "
     "Without guardrails, broader confidence ingestion worsens ensemble performance."),
    ("WES & RNA-seq operate near the ensemble ceiling",
     "ChampionChallenger reaches 0.9487 WES / 0.9533 RNA. "
     "Abstention rates are low after threshold tuning (1.5% RNA, 1.5% WES)."),
    ("Bimodal (WES+RNA) is the best multimodal strategy",
     "TrimodalMajorityVote does not improve over BimodalMajorityVote. "
     "WGS adds noise not signal at current cohort size and tool quality."),
    ("Remaining before submission",
     "Supplementary assembly (Tables S1–S3, Figures S1–S3). "
     "Figure 1 (architecture diagram — manual, outside AI scope)."),
]

y0 = Inches(1.25)
for i, (heading, body) in enumerate(conclusions):
    y = y0 + i * Inches(1.12)
    add_rect(sc, Inches(0.3), y, Inches(0.3), Inches(0.9), C_BLUE)
    add_text(sc, heading,
             Inches(0.75), y + Inches(0.04), Inches(11.8), Inches(0.38),
             font_size=13, bold=True, color=C_BLUE)
    add_text(sc, body,
             Inches(0.75), y + Inches(0.42), Inches(11.8), Inches(0.56),
             font_size=11, color=C_DARK, wrap=True)


# ── Save ─────────────────────────────────────────────────────────────────────
prs.save(OUT_PATH)
print(f"Saved: {OUT_PATH}")
