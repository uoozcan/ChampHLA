# PIHLA v3.0.0 — Release Summary

**Date:** 2026-04-18
**Status:** Publication-ready benchmark complete

---

## What's New in v3.0.0

### 8-Tool HLA Typing Ensemble

| Tool | Modalities | Notes |
|------|-----------|-------|
| OptiType | WGS, WES, RNA-seq | Best WGS single-tool (43.3%) |
| ArcasHLA | WGS, WES, RNA-seq | RNA-seq optimised; WGS weight 0.042 (94% shrinkage) |
| SpecHLA | WGS, WES | Variant-aware typing |
| HLA-HD | WGS, WES | High-resolution; WGS poor_calibration |
| Kourami | WGS | Graph-based; WGS poor_calibration |
| T1K | WGS, WES, RNA-seq | RNA-seq best at 90.5% |
| POLYSOLVER | WGS, WES | no_confidence on most modalities |
| Seq2HLA | RNA-seq | Transcript-based |

### read_confidence_v2 Weighting

```
final_weight = 0.7 × base_reliability + 0.3 × effective_confidence_score
```

Guardrail statuses per tool × modality:
- `applied` — confidence integration active with shrinkage
- `poor_calibration` — Brier/ECE > 0.35; confidence score discounted
- `no_confidence` — <50% coverage; weight falls back to base reliability

### Tri-modal Benchmark

| Modality | Typed | Holdout | WeightedConsensus |
|----------|-------|---------|-------------------|
| WGS | 99 | 20 | 26.7% (callable 56.7%) |
| WES | 51 | 12 | 91.7% |
| RNA-seq | 50 | 7 | 100.0% |

Split: 60/20/20 (training/validation/holdout). Truth: Gourraud 2014, loci A/B/C.
IMGT/HLA: v3.59.0. Benchmark scope: `partial_realdata_benchmark`.

---

## Core Pipeline Files

| File | Description |
|------|-------------|
| `main.nf` | Main Nextflow pipeline (DSL2) |
| `nextflow.config` | Pipeline configuration with CSC Puhti profile |
| `conf/benchmark_1000g_full_cohort.yaml` | Tri-modal benchmark config |
| `conf/benchmark_1000g_wes.yaml` | WES-only benchmark config |
| `conf/benchmark_1000g_rna.yaml` | RNA-seq-only benchmark config |
| `bin/run_1000g_benchmark.py` | Benchmark driver script |
| `bin/generate_figures_v2.py` | Figure generation (PNG 300dpi + PDF) |
| `analysis/generate_report_v4.py` | HTML report generator |

## Key Output Files

| File | Description |
|------|-------------|
| `analysis/benchmark_trimodal_50samples/tables/` | All benchmark TSVs |
| `analysis/figures_final/` | 8 publication figures (PNG + PDF) |
| `analysis/pihla_benchmark_report_v4.html` | Standalone benchmark report (3.0 MB) |
| `docs/MANUSCRIPT_DRAFT_V1.md` | Manuscript draft with final numbers |

## Publication Figures

| Figure | Description |
|--------|-------------|
| figure_01_cohort_overview | Cohort design and sample counts |
| figure_02_method_performance | Method comparison across modalities |
| figure_03_per_gene_heatmap | Per-gene accuracy heatmap |
| figure_04_calibration_diagrams | Confidence calibration curves |
| figure_05_abstention_tradeoff | Callable rate vs accuracy |
| figure_06_discordance_taxonomy | DNA/RNA discordance taxonomy |
| figure_07_weight_heatmap | Confidence weight heatmap |
| figure_08_ensemble_advantage | Ensemble advantage over single-tool baselines |

---

## SLURM Scripts

| Script | Purpose |
|--------|---------|
| `slurm_benchmark_trimodal.sh` | Run full tri-modal benchmark |
| `slurm_benchmark_wes_50.sh` | Run WES benchmark |
| `slurm_benchmark_rna_50.sh` | Run RNA-seq benchmark |
| `slurm_generate_figures.sh` | Regenerate all figures |
| `slurm_wgs_wave2_typing.sh` | Type 49 wave-2 WGS samples (array, %10 throttle) |
| `slurm_download_wgs_wave2.sh` | Download NYGC 30x CRAMs from EBI (array, %10 throttle) |

---

## Known Limitations

- Benchmark loci limited to HLA-A, -B, -C (Gourraud 2014 truth source)
- DRB1/DQB1 excluded: no truth labels available in current truth source
- WGS accuracy low by design (whole-genome depth; callable rate 56.7% at threshold 0.45)
- Holdout sizes: WGS n=20, WES n=12, RNA-seq n=7

---

## HPC Environment

- Platform: CSC Puhti (SLURM)
- Account: project_2008084
- Container runtime: Singularity
- Python environment: `module load python-data`
- Nextflow: ≥23.04.0
