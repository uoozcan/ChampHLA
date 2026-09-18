# Benchmark Figures

The scientific benchmark workflow now targets only real 1000 Genomes samples with public HLA ground truth and matched WGS, WES, and RNA-seq availability. The synthetic fixture remains in the repo only for CI and parser validation.

## Real-Data Workflow

1. Build canonical manifests:
```bash
python3 bin/build_1000g_benchmark_manifests.py \
  --truth /path/to/1000g_hla_truth.tsv \
  --sequencing /path/to/1000g_sequencing_source.tsv \
  --output-dir results/1000g_manifests \
  --acquisition-date 2014-07-25 \
  --supported-loci A,B,C,DRB1,DQB1
```

2. Run the split-aware benchmark:
```bash
python3 bin/run_1000g_benchmark.py \
  --config conf/benchmark_1000g_config.example.yaml \
  --output-dir results/1000g_benchmark
```

The runner:
- assembles `truth_manifest.tsv`, `sequencing_manifest.tsv`, and `cohort_manifest.tsv`
- restricts analysis to strict tri-modal truth-backed samples
- learns confidence weights on the training split
- optionally tunes abstention support on validation
- reports headline benchmark outputs from holdout only

## Truth and Splits

Truth source priority for this workflow is fixed to the public 1000 Genomes HLA dataset from Gourraud et al. The real-data benchmark uses:
- `training`: confidence calibration and weight learning
- `validation`: abstention threshold tuning when available
- `holdout`: final reported performance only

Primary benchmark loci:
- `A`
- `B`
- `C`
- `DRB1`
- `DQB1`

Primary resolution:
- exact two-field allele-pair match

Secondary resolution:
- exact three-field allele-pair match where derivable and comparable

## Supported Benchmark Parsers

Generic parser families:
- `spec_hla_result`
- `optitype_tsv`
- `arcashla_json`
- `wide_hla_table`
- `long_hla_table`
- `auto`

Tool aliases accepted by the benchmark layer:
- `hlahd_table`
- `polysolver_table`
- `kourami_table`
- `t1k_table`
- `seq2hla_table`

These aliases support benchmark-time ingestion of precomputed result tables. The current PIHLA execution layer now includes native workflow routes for HLA-HD, POLYSOLVER, Kourami, T1K, and Seq2HLA as well, but phase-gated reporting still applies whenever a tool fails to complete successfully in a given modality.

## Current Primary Results Set

The current primary quantitative benchmark is the expanded 50-sample WGS-only wave in:
- `/scratch/project_2008084/pihla-publish/analysis/1000g_realdata/benchmark_wgs_wave2`

This wave should now be treated as the authoritative WGS benchmark because it provides:
- stable sample-level training/validation/holdout splits (`28/10/12`)
- a larger truth-backed six-tool WGS cohort than the original 42-sample wave
- non-empty calibration outputs for `OptiType`, `T1K`, `ArcasHLA`, `HLA-HD`, and `Kourami`
- finalized baseline and routed-method comparison tables
- a direct comparison artifact against the earlier 42-sample wave

The 42-sample WGS wave remains useful as a development and calibration reference, but it should no longer be the main manuscript-facing quantitative benchmark. The earlier tri-modal three-sample run remains useful for workflow integration and multi-modality demonstrations, but not as the primary statistical comparison set.

## WGS Wave 2 Figure Mapping

### Figure 2. WGS cohort assembly and split design
Data source:
- `analysis/1000g_realdata/wgs_wave2_inputs/truth_manifest.tsv`
- `analysis/1000g_realdata/wgs_wave2_inputs/cohort_manifest.tsv`
- `analysis/1000g_realdata/benchmark_wgs_wave2/tables/benchmark_metadata.json`
Primary message:
- the primary WGS benchmark now uses 50 truth-backed samples with deterministic, population-aware splits.

### Figure 3. WGS holdout method comparison
Data source:
- `analysis/1000g_realdata/benchmark_wgs_wave2/tables/method_comparison.tsv`
- `analysis/1000g_realdata/benchmark_wgs_wave2/tables/champion_challenger_method_comparison.tsv`
- `analysis/1000g_realdata/benchmark_wgs_wave2/tables/locus_expert_method_comparison.tsv`
Primary message:
- the earlier 42-sample baseline tie does not survive the expanded cohort: `WeightedConsensus` falls behind `OptiType`, while the strongest routed baselines recover but do not exceed the `OptiType` ceiling.

### Figure 4. WGS per-gene tool behavior
Data source:
- `analysis/1000g_realdata/benchmark_wgs_wave2/tables/method_per_gene.tsv`
- `analysis/1000g_realdata/benchmark_wgs_wave2/tables/champion_challenger_method_per_gene.tsv`
- `analysis/1000g_realdata/benchmark_wgs_wave2/tables/locus_difficulty_summary.tsv`
Primary message:
- WGS performance remains strongly gene-dependent, with `HLA-B` strongest, `HLA-A` intermediate, and `HLA-C` the dominant unresolved failure locus.

### Figure 5. WGS calibration and confidence reliability
Data source:
- `analysis/1000g_realdata/benchmark_wgs_wave2/tables/confidence_calibration_summary.tsv`
- `analysis/1000g_realdata/benchmark_wgs_wave2/tables/tool_confidence_weights.tsv`
- `analysis/1000g_realdata/benchmark_wgs_wave2/tables/wgs_wave2_vs_wave1_diff.tsv`
Primary message:
- the calibration story is stable after cohort expansion: recalibration remains methodologically necessary, but the larger wave does not rescue weighted consensus above the best single tool.

### Figure 6. WGS disagreement and routed-baseline interpretation
Data source:
- `analysis/1000g_realdata/benchmark_wgs_wave2/tables/consensus_decision_trace.tsv`
- `analysis/1000g_realdata/benchmark_wgs_wave2/tables/champion_challenger_trace.tsv`
- `analysis/1000g_realdata/benchmark_wgs_wave2/tables/tool_vs_truth_error_taxonomy.tsv`
Primary message:
- routed baselines can recover the strongest single-tool decisions on difficult loci, but disagreement structure still points to a tool-limited WGS regime rather than a remaining consensus-architecture gap.

### Figure 7. WGS learned confidence weights
Data source:
- `analysis/1000g_realdata/benchmark_wgs_wave2/tables/tool_confidence_weights.tsv`
- `analysis/1000g_realdata/benchmark_wgs_wave2/tables/tool_confidence_weights_by_gene.tsv`
Primary message:
- runtime weights remain stable after expansion, with `OptiType` still dominant and other tools contributing limited or guarded signal.

## Main Outputs

Tables emitted by the real-data benchmark:
- `truth_manifest.tsv`
- `sequencing_manifest.tsv`
- `cohort_manifest.tsv`
- `harmonized_benchmark_rows.tsv`
- `cohort_overview.tsv`
- `summary_full_cohort.tsv`
- `summary_per_gene.tsv`
- `ambiguity_summary.tsv`
- `ambiguity_summary_by_gene.tsv`
- `tool_confidence_weights.tsv`
- `tool_confidence_weights_by_gene.tsv`
- `majority_vote_baseline.tsv`
- `weighted_consensus_calls.tsv`
- `method_comparison.tsv`
- `method_ambiguity_summary.tsv`
- `method_ambiguity_summary_by_gene.tsv`
- `confidence_calibration_summary.tsv`
- `confidence_error_summary.tsv`
- `confidence_error_summary_by_gene.tsv`
- `abstention_tradeoff.tsv`
- `discordance_summary.tsv`
- `reference_metadata.tsv`
- `benchmark_metadata.json`
- `consensus_runtime_weights.json`

Figures:
- `figure_2_accuracy_comparison.svg`
- `figure_3_per_gene_gains.svg`
- `figure_4_confidence_calibration.svg`
- `figure_5_abstention_tradeoff.svg`
- `figure_6_discordance_taxonomy.svg`
- `figure_7_confidence_weights.svg`
- `captions.md`

## Metadata Requirements

`benchmark_metadata.json` records:
- truth source and acquisition date
- final tri-modal cohort size
- per-population counts
- per-modality sample counts before and after filtering
- split membership summary
- excluded-sample counts by reason
- supported loci
- tuned consensus support threshold

## Figure Shells

### Figure 1. Workflow and ensemble architecture
Data source:
- workflow DAG and schema docs
Panels:
- tool execution layer by modality
- harmonized call/confidence schema
- benchmark calibration and weight learning
- runtime consensus and abstention outputs
Primary message:
- PIHLA is a calibrated workflow platform, not only a tool wrapper.

### Figure 2. Cohort assembly and benchmark design
Data source:
- `truth_manifest.tsv`
- `sequencing_manifest.tsv`
- `cohort_manifest.tsv`
- `benchmark_metadata.json`
Panels:
- truth ingestion
- tri-modal intersection
- exclusions by reason
- final split assignment
Primary message:
- the scientific benchmark is a strict real-data matched cohort with explicit filtering.

### Figure 3. Method comparison on holdout
Data source:
- `method_comparison.tsv`
- `summary_full_cohort.tsv`
Panels:
- per-modality A/B/C accuracy
- callable-rate overlay
- weighted consensus vs majority vote vs best single tool
Primary message:
- calibrated consensus should be compared against both single-tool and equal-weight baselines.

### Figure 4. Per-gene gains
Data source:
- `summary_per_gene.tsv`
- `method_ambiguity_summary_by_gene.tsv`
Panels:
- A, B, and C gains by modality
- weighted consensus gain over majority vote
Primary message:
- ensemble benefit is gene-specific and should not be summarized only as a cohort-wide mean.

### Figure 5. Calibration and confidence-stratified error
Data source:
- `confidence_calibration_summary.tsv`
- `confidence_error_summary.tsv`
- `confidence_error_summary_by_gene.tsv`
Panels:
- calibration curve by modality
- Brier/ECE bars
- confidence-bin error rates
Primary message:
- confidence is meaningful only when tied to observed correctness and passed through explicit calibration-aware guardrails before runtime use.

### Figure 6. Abstention and disagreement interpretation
Data source:
- `abstention_tradeoff.tsv`
- `discordance_summary.tsv`
Panels:
- abstention tradeoff curves
- discordance tag counts by modality
Primary message:
- uncertainty handling should separate weak evidence from hard disagreement.

### Figure 7. Learned confidence weights
Data source:
- `tool_confidence_weights.tsv`
- `tool_confidence_weights_by_gene.tsv`
- `consensus_runtime_weights.json`
Panels:
- tool-level weights by modality
- gene-level weight heatmap
Primary message:
- benchmark-derived weighting changes runtime voting behavior in a traceable way, and the guardrail layer makes those changes scientifically interpretable rather than purely score-driven.

## Table Shells

- Table 1. Supported tools, modalities, confidence evidence, and native runtime integration.
- Table 2. Cohort assembly, public truth source, supported loci, and split design.
- Table 3. Holdout method comparison across single tools, majority vote, and weighted consensus.
- Table 4. Per-gene A/B/C performance and gain over majority vote.
- Table 5. Calibration, confidence-coverage, and abstention summary metrics.
- Table 6. Discordance taxonomy and phase-gated tool availability notes.
