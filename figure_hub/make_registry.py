#!/usr/bin/env python3.11
"""
make_registry.py — Generate the MVHLA master figure registry.

Outputs (both in figure_hub/):
  figure_registry.tsv   — tab-separated, importable into Excel / Numbers
  figure_registry.md    — markdown table for quick scanning

Update this script whenever a figure is added, removed, or changed.
Re-run to refresh both outputs:
    python3.11 figure_hub/make_registry.py
"""

from pathlib import Path
import csv
import io

HUB = Path(__file__).resolve().parent
TSV_OUT = HUB / "figure_registry.tsv"
MD_OUT  = HUB / "figure_registry.md"

# ---------------------------------------------------------------------------
# REGISTRY — edit this list to update the registry
# ---------------------------------------------------------------------------
# Column keys (in display order):
#   figure_id, filename_stem, manuscript_role, manuscript_section,
#   active, source_analysis, hub_script, original_script,
#   arial_font, legend_outside, formats, last_updated, notes

REGISTRY = [
    # ── Active main figures ─────────────────────────────────────────────────
    dict(
        figure_id        = "Fig 1",
        filename_stem    = "figure_1_workflow_architecture",
        manuscript_role  = "main",
        manuscript_section = "Methods",
        active           = "yes",
        source_analysis  = "manual / architecture diagram",
        hub_script       = "—",
        original_script  = "manual (outside AI scope)",
        arial_font       = "unknown",
        legend_outside   = "N/A",
        formats          = "pdf+svg+png",
        last_updated     = "2026-06-01",
        notes            = "Architecture schematic; manually prepared. Benchmark hierarchy explicit in view.",
    ),
    dict(
        figure_id        = "Fig 2",
        filename_stem    = "figure_2_accuracy_comparison",
        manuscript_role  = "main",
        manuscript_section = "Results",
        active           = "yes",
        source_analysis  = "analysis/1000g_realdata/benchmark_wgs_wave2",
        hub_script       = "scripts/fig_benchmark_main.py",
        original_script  = "bin/generate_figures.py",
        arial_font       = "yes",
        legend_outside   = "yes",
        formats          = "pdf+svg+png",
        last_updated     = "2026-06-03",
        notes            = "Primary WGS anchor result (n=50, holdout n=12). Regenerated 2026-06-03 with fig_benchmark_main.py; Arial enforced via apply_style().",
    ),
    dict(
        figure_id        = "Fig 3",
        filename_stem    = "figure_3_per_gene_gains",
        manuscript_role  = "main",
        manuscript_section = "Results",
        active           = "yes",
        source_analysis  = "analysis/1000g_realdata/benchmark_wgs_wave2 (per-gene tables)",
        hub_script       = "scripts/fig_benchmark_main.py",
        original_script  = "bin/generate_figures.py",
        arial_font       = "yes",
        legend_outside   = "yes",
        formats          = "pdf+svg+png",
        last_updated     = "2026-06-03",
        notes            = "Per-gene WGS failure-mode. HLA-B easiest (0.83), HLA-C hardest (0.17). Regenerated 2026-06-03; Arial enforced.",
    ),
    dict(
        figure_id        = "Fig 4",
        filename_stem    = "figure_4_confidence_calibration",
        manuscript_role  = "main",
        manuscript_section = "Results",
        active           = "yes",
        source_analysis  = "analysis/1000g_realdata/benchmark_wgs_wave2 (calibration tables)",
        hub_script       = "scripts/fig_benchmark_main.py",
        original_script  = "bin/generate_figures.py",
        arial_font       = "yes",
        legend_outside   = "yes",
        formats          = "pdf+svg+png",
        last_updated     = "2026-06-03",
        notes            = "Confidence guardrails: OptiType/HLA-HD/Kourami blocked; T1K partial boost. Regenerated 2026-06-03; Arial enforced.",
    ),
    dict(
        figure_id        = "Fig 5",
        filename_stem    = "figure_5_abstention_tradeoff",
        manuscript_role  = "main",
        manuscript_section = "Results",
        active           = "yes",
        source_analysis  = "analysis/1000g_realdata/benchmark_wgs_wave2 + WES + RNA (calibration/confidence)",
        hub_script       = "scripts/fig_benchmark_main.py",
        original_script  = "bin/generate_figures.py",
        arial_font       = "yes",
        legend_outside   = "yes",
        formats          = "pdf+svg+png",
        last_updated     = "2026-06-03",
        notes            = "Abstention–accuracy tradeoff. Rates: 8% WGS, 7% WES, 4% RNA. Regenerated 2026-06-03; Arial enforced.",
    ),
    dict(
        figure_id        = "Fig 6",
        filename_stem    = "figure_6_discordance_taxonomy",
        manuscript_role  = "main",
        manuscript_section = "Results",
        active           = "yes",
        source_analysis  = "analysis/1000g_realdata/benchmark_wgs_wave2 (discordance summaries)",
        hub_script       = "scripts/fig_benchmark_main.py",
        original_script  = "bin/generate_figures.py",
        arial_font       = "yes",
        legend_outside   = "yes",
        formats          = "pdf+svg+png",
        last_updated     = "2026-06-03",
        notes            = "5-category discordance taxonomy. WGS low_evidence n=157; RNA possible_expression_bias n=292. Regenerated 2026-06-03; Arial enforced.",
    ),
    dict(
        figure_id        = "Fig 7",
        filename_stem    = "figure_7_confidence_weights",
        manuscript_role  = "main",
        manuscript_section = "Results",
        active           = "yes",
        source_analysis  = "analysis/1000g_realdata/benchmark_wgs_wave2 (weight tables)",
        hub_script       = "scripts/fig_benchmark_main.py",
        original_script  = "bin/generate_figures.py",
        arial_font       = "yes",
        legend_outside   = "yes",
        formats          = "pdf+svg+png",
        last_updated     = "2026-06-03",
        notes            = "Final runtime weights: OptiType 0.5476, T1K 0.2649, HLA-HD 0.2294. Regenerated 2026-06-03; Arial enforced.",
    ),
    dict(
        figure_id        = "Fig 8 (file: 09)",
        filename_stem    = "figure_09_bimodal_per_gene",
        manuscript_role  = "main",
        manuscript_section = "Results (multimodal section)",
        active           = "yes",
        source_analysis  = "analysis/1000g_realdata/benchmark_trimodal_robustness/run (WES+RNA matched n=106)",
        hub_script       = "scripts/fig_benchmark_main.py",
        original_script  = "bin/generate_figures_pub.py",
        arial_font       = "yes",
        legend_outside   = "yes",
        formats          = "pdf+svg+png",
        last_updated     = "2026-06-03",
        notes            = "Bimodal WES+RNA: strongest multimodal result. Regenerated 2026-06-03 with fig_benchmark_main.py; Arial enforced via apply_style().",
    ),
    dict(
        figure_id        = "Fig 9 (file: 10)",
        filename_stem    = "figure_10_trimodal_comparison",
        manuscript_role  = "main",
        manuscript_section = "Results (multimodal section)",
        active           = "yes",
        source_analysis  = "analysis/1000g_realdata/benchmark_trimodal_robustness/run",
        hub_script       = "scripts/fig_benchmark_main.py",
        original_script  = "bin/generate_figures_pub.py",
        arial_font       = "yes",
        legend_outside   = "yes",
        formats          = "pdf+svg+png",
        last_updated     = "2026-06-03",
        notes            = "Trimodal vs bimodal: adding WGS does not improve over WES+RNA. Regenerated 2026-06-03 with fig_benchmark_main.py; Arial enforced.",
    ),
    # ── Active supplementary ────────────────────────────────────────────────
    dict(
        figure_id        = "Fig S1",
        filename_stem    = "figure_s1_per_gene_accuracy",
        manuscript_role  = "supplementary",
        manuscript_section = "Supplementary",
        active           = "yes",
        source_analysis  = "per-modality benchmark runs (WGS/WES/RNA all samples)",
        hub_script       = "scripts/fig_benchmark_main.py",
        original_script  = "bin/generate_figures_pub.py",
        arial_font       = "yes",
        legend_outside   = "yes",
        formats          = "pdf+svg+png",
        last_updated     = "2026-06-03",
        notes            = "Per-gene accuracy all tools × modalities. Regenerated 2026-06-03 with fig_benchmark_main.py; Arial enforced.",
    ),
    dict(
        figure_id        = "Fig S2",
        filename_stem    = "figure_s2_calibration_heatmap",
        manuscript_role  = "supplementary",
        manuscript_section = "Supplementary",
        active           = "yes",
        source_analysis  = "per-modality calibration tables (ECE, Brier scores)",
        hub_script       = "scripts/fig_benchmark_main.py",
        original_script  = "bin/generate_figures_pub.py",
        arial_font       = "yes",
        legend_outside   = "yes",
        formats          = "pdf+svg+png",
        last_updated     = "2026-06-03",
        notes            = "Calibration heatmap per tool × modality. ECE/Brier threshold 0.35. Regenerated 2026-06-03; Arial enforced.",
    ),
    dict(
        figure_id        = "Fig S3",
        filename_stem    = "figure_s3_resolution_comparison",
        manuscript_role  = "supplementary",
        manuscript_section = "Supplementary",
        active           = "yes",
        source_analysis  = "per-modality benchmark tables (2-field vs 3-field, G/P group)",
        hub_script       = "scripts/fig_benchmark_main.py",
        original_script  = "bin/generate_figures_pub.py",
        arial_font       = "yes",
        legend_outside   = "yes",
        formats          = "pdf+svg+png",
        last_updated     = "2026-06-03",
        notes            = "Resolution comparison (2-field primary; 3-field/G-group secondary). Regenerated 2026-06-03; Arial enforced.",
    ),
    # ── FIMM-specific figures ───────────────────────────────────────────────
    dict(
        figure_id        = "FIMM-CONC-A",
        filename_stem    = "fig_concordance_by_gene",
        manuscript_role  = "fimm_only",
        manuscript_section = "N/A (internal / supplementary candidate)",
        active           = "yes",
        source_analysis  = "pihla_local/fimm_results/fimm_hla_calls_deident.tsv",
        hub_script       = "scripts/fig_fimm_concordance.py",
        original_script  = "bin/analyze_fimm_hla.py",
        arial_font       = "yes",
        legend_outside   = "yes",
        formats          = "pdf+svg+png",
        last_updated     = "2026-06-03",
        notes            = "Per-gene tool-pair concordance heatmap (genes A, B, C). Regenerated 2026-06-03 with fig_fimm_concordance.py using fimm_hla_calls_deident.tsv; Arial enforced.",
    ),
    dict(
        figure_id        = "FIMM-CROSS",
        filename_stem    = "fig_crossmodal_comparison",
        manuscript_role  = "fimm_only",
        manuscript_section = "N/A",
        active           = "yes",
        source_analysis  = "pihla_local/fimm_results/fimm_hla_calls_deident.tsv",
        hub_script       = "scripts/fig_fimm_concordance.py",
        original_script  = "bin/analyze_fimm_hla.py",
        arial_font       = "yes",
        legend_outside   = "yes",
        formats          = "pdf+svg+png",
        last_updated     = "2026-06-03",
        notes            = "Cross-modality concordance bars (scRNA vs WES vs BulkRNA). Regenerated 2026-06-03 with fig_fimm_concordance.py; Arial enforced.",
    ),
    dict(
        figure_id        = "FIMM-SCRNA",
        filename_stem    = "fig_scrna_concordance_heatmap",
        manuscript_role  = "fimm_only",
        manuscript_section = "N/A",
        active           = "yes",
        source_analysis  = "pihla_local/fimm_results/fimm_hla_calls_deident.tsv",
        hub_script       = "scripts/fig_fimm_concordance.py",
        original_script  = "bin/analyze_fimm_hla.py",
        arial_font       = "yes",
        legend_outside   = "yes",
        formats          = "pdf+svg+png",
        last_updated     = "2026-06-03",
        notes            = "Per-sample scRNA vs WES concordance heatmap. Regenerated 2026-06-03 with fig_fimm_concordance.py; Arial enforced.",
    ),
    dict(
        figure_id        = "FIMM-LOH",
        filename_stem    = "figure_fimm_loh_wes",
        manuscript_role  = "fimm_only",
        manuscript_section = "N/A (exploratory / supplementary candidate)",
        active           = "yes",
        source_analysis  = "pihla_local/fimm_results/fimm_loh_candidates_wes_deident.tsv",
        hub_script       = "scripts/fig_fimm_loh_wes.py",
        original_script  = "bin/analysis/analyze_fimm_loh_wes.py",
        arial_font       = "yes",
        legend_outside   = "yes",
        formats          = "pdf+svg+png",
        last_updated     = "2026-06-03",
        notes            = "LOH status stacked bar per HLA gene. WES exploratory (WGS confirmation needed). Regenerated 2026-06-03 with fig_fimm_loh_wes.py; Arial enforced.",
    ),
    # ── Removed / Legacy ────────────────────────────────────────────────────
    dict(
        figure_id        = "Removed",
        filename_stem    = "figure_a1_guardrail_scatter",
        manuscript_role  = "removed",
        manuscript_section = "N/A",
        active           = "no",
        source_analysis  = "WGS calibration tables",
        hub_script       = "—",
        original_script  = "bin/generate_figures_analysis.py",
        arial_font       = "partial",
        legend_outside   = "unknown",
        formats          = "pdf+png",
        last_updated     = "2026-06-01",
        notes            = "Guardrail scatter — archived in legacy_pre_redesign_2026-06-01/. Superseded by Fig 4.",
    ),
    dict(
        figure_id        = "Removed",
        filename_stem    = "figure_a2_tool_ceiling",
        manuscript_role  = "removed",
        manuscript_section = "N/A",
        active           = "no",
        source_analysis  = "WGS benchmark tables",
        hub_script       = "—",
        original_script  = "bin/generate_figures_analysis.py",
        arial_font       = "partial",
        legend_outside   = "unknown",
        formats          = "pdf+png",
        last_updated     = "2026-06-01",
        notes            = "Tool ceiling plot — archived in legacy. Content merged into Fig 2.",
    ),
    dict(
        figure_id        = "Removed",
        filename_stem    = "figure_a3_weight_decomposition",
        manuscript_role  = "removed",
        manuscript_section = "N/A",
        active           = "no",
        source_analysis  = "WGS weight tables",
        hub_script       = "—",
        original_script  = "bin/generate_figures_analysis.py",
        arial_font       = "partial",
        legend_outside   = "unknown",
        formats          = "pdf+png",
        last_updated     = "2026-06-01",
        notes            = "Weight decomposition — archived. Content merged into Fig 7.",
    ),
    dict(
        figure_id        = "Removed",
        filename_stem    = "figure_a8_calibration_ci",
        manuscript_role  = "removed",
        manuscript_section = "N/A",
        active           = "no",
        source_analysis  = "WGS calibration CI tables",
        hub_script       = "—",
        original_script  = "bin/generate_figures_analysis.py",
        arial_font       = "partial",
        legend_outside   = "unknown",
        formats          = "pdf+png",
        last_updated     = "2026-06-01",
        notes            = "Calibration confidence intervals — archived. Content in Fig S2.",
    ),
    dict(
        figure_id        = "Removed",
        filename_stem    = "figure_a9_weight_sensitivity",
        manuscript_role  = "removed",
        manuscript_section = "N/A",
        active           = "no",
        source_analysis  = "analysis/weight_sensitivity/",
        hub_script       = "—",
        original_script  = "bin/generate_figures_analysis.py",
        arial_font       = "partial",
        legend_outside   = "unknown",
        formats          = "pdf+png",
        last_updated     = "2026-06-01",
        notes            = "Weight sensitivity analysis — archived. Future-work note in manuscript.",
    ),
    dict(
        figure_id        = "Removed",
        filename_stem    = "figure_f4_abstention_phasespace",
        manuscript_role  = "removed",
        manuscript_section = "N/A",
        active           = "no",
        source_analysis  = "WGS/WES/RNA benchmark runs",
        hub_script       = "—",
        original_script  = "bin/analysis/generate_supplementary_figures.py",
        arial_font       = "partial",
        legend_outside   = "unknown",
        formats          = "pdf+png",
        last_updated     = "2026-06-01",
        notes            = "Abstention phase-space — archived. Concept retained in Fig 5.",
    ),
    dict(
        figure_id        = "Removed",
        filename_stem    = "figure_f6_discordance_tags",
        manuscript_role  = "removed",
        manuscript_section = "N/A",
        active           = "no",
        source_analysis  = "WGS discordance summaries",
        hub_script       = "—",
        original_script  = "bin/analysis/generate_supplementary_figures.py",
        arial_font       = "partial",
        legend_outside   = "unknown",
        formats          = "pdf+png",
        last_updated     = "2026-06-01",
        notes            = "Discordance tag breakdown — archived. Merged into Fig 6.",
    ),
    dict(
        figure_id        = "Removed",
        filename_stem    = "figure_f7_sample_heatmap",
        manuscript_role  = "removed",
        manuscript_section = "N/A",
        active           = "no",
        source_analysis  = "WGS per-sample results",
        hub_script       = "—",
        original_script  = "bin/analysis/generate_supplementary_figures.py",
        arial_font       = "partial",
        legend_outside   = "unknown",
        formats          = "pdf+png",
        last_updated     = "2026-06-01",
        notes            = "Per-sample heatmap — archived (too granular for manuscript).",
    ),
    dict(
        figure_id        = "Removed",
        filename_stem    = "figure_loh_allele_dropout_heatmap",
        manuscript_role  = "removed",
        manuscript_section = "N/A",
        active           = "no",
        source_analysis  = "analysis/loh_analysis/benchmark_homozygosity_error_summary.tsv",
        hub_script       = "scripts/fig_loh_dropout_heatmap.py",
        original_script  = "bin/analysis/figure_loh_allele_dropout_heatmap.py",
        arial_font       = "no",
        legend_outside   = "partial",
        formats          = "pdf+png",
        last_updated     = "2026-06-01",
        notes            = "Benchmark LOH/allele-dropout heatmap (2-panel). Hub script exists. Removed from active set; may be reinstated as supplementary.",
    ),
    dict(
        figure_id        = "Removed",
        filename_stem    = "figure_per_gene_barplot",
        manuscript_role  = "removed",
        manuscript_section = "N/A",
        active           = "no",
        source_analysis  = "pre-wave2 per-gene benchmark tables",
        hub_script       = "—",
        original_script  = "bin/figure_per_gene_barplot.py",
        arial_font       = "partial",
        legend_outside   = "unknown",
        formats          = "pdf+png",
        last_updated     = "2026-04-27",
        notes            = "Pre-wave2 per-gene barplot — archived in legacy_pre_wave2/. Superseded by Fig 3.",
    ),
]

COLUMNS = [
    "figure_id", "filename_stem", "manuscript_role", "manuscript_section",
    "active", "source_analysis", "hub_script", "original_script",
    "arial_font", "legend_outside", "formats", "last_updated", "notes",
]

# ---------------------------------------------------------------------------
# Write TSV
# ---------------------------------------------------------------------------
with TSV_OUT.open("w", newline="", encoding="utf-8") as fh:
    writer = csv.DictWriter(fh, fieldnames=COLUMNS, delimiter="\t",
                            extrasaction="ignore")
    writer.writeheader()
    writer.writerows(REGISTRY)

print(f"TSV  -> {TSV_OUT}  ({len(REGISTRY)} rows)")

# ---------------------------------------------------------------------------
# Write Markdown
# ---------------------------------------------------------------------------

def md_row(rec: dict) -> str:
    cells = [rec.get(c, "") for c in COLUMNS]
    return "| " + " | ".join(cells) + " |"

header = "| " + " | ".join(COLUMNS) + " |"
sep    = "| " + " | ".join(["---"] * len(COLUMNS)) + " |"

lines = [
    "# MVHLA Figure Registry",
    "",
    f"**{len(REGISTRY)} figures total** — generated by `figure_hub/make_registry.py`  ",
    f"Last regenerated: 2026-06-02",
    "",
    "## Active — Main manuscript",
    "",
    header, sep,
]
for r in REGISTRY:
    if r["manuscript_role"] == "main":
        lines.append(md_row(r))

lines += ["", "## Active — Supplementary", "", header, sep]
for r in REGISTRY:
    if r["manuscript_role"] == "supplementary":
        lines.append(md_row(r))

lines += ["", "## Active — FIMM only", "", header, sep]
for r in REGISTRY:
    if r["manuscript_role"] == "fimm_only":
        lines.append(md_row(r))

lines += ["", "## Removed / Legacy", "", header, sep]
for r in REGISTRY:
    if r["manuscript_role"] == "removed":
        lines.append(md_row(r))

lines += [
    "",
    "---",
    "",
    "## Arial font status legend",
    "",
    "| Value | Meaning |",
    "|---|---|",
    "| `yes` | Script explicitly sets Arial as primary font via hub `plot_style.py` |",
    "| `partial` | Script uses `Avenir Next → Arial → Helvetica` fallback; Arial renders on Linux HPC if installed |",
    "| `no` | Script uses DejaVu Sans or system sans-serif; Arial not in font list |",
    "| `unknown` | Manually produced figure; font not controlled by script |",
    "",
    "Figures marked `partial` or `no` should be **regenerated** using the corresponding `figure_hub/scripts/fig_*.py`.",
]

MD_OUT.write_text("\n".join(lines) + "\n", encoding="utf-8")
print(f"MD   -> {MD_OUT}")
