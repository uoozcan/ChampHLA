# Benchmark Figures

This workflow builds manuscript-style benchmark tables and figures for single-tool baselines, majority vote, and weighted HLA consensus across WES, WGS, and RNA-seq.

## Truth Hierarchy

Use the strongest truth source available for each cohort in this order:
1. orthogonal clinical typing
2. targeted HLA NGS
3. public reference truth

The active hierarchy is recorded in `benchmark_metadata.json`.

## Benchmark Splits

Recommended split roles:
- `training`: fit confidence calibration and benchmark-derived weights
- `validation`: tune abstention thresholds
- `holdout`: final reporting only

The active split labels are recorded in `benchmark_metadata.json`.

## Frozen Schemas

See [HLA_ENSEMBLE_SCHEMAS.md](/users/ozcanumu/scratch/project_2008084/pihla-publish/docs/HLA_ENSEMBLE_SCHEMAS.md) for the frozen TSV/JSON interfaces used by the benchmark layer.

## Supported Result Parsers

- `spec_hla_result`
- `optitype_tsv`
- `arcashla_json`
- `wide_hla_table`
- `long_hla_table`
- `auto`

## Supported Confidence Parsers

- `optitype_result_objective`
- `long_confidence_table`
- `wide_confidence_table`
- `json_confidence_table`

## Main Benchmark Outputs

### Tables
- `harmonized_benchmark_rows.tsv`
- `summary_full_cohort.tsv`
- `summary_per_gene.tsv`
- `tool_confidence_weights.tsv`
- `tool_confidence_weights_by_gene.tsv`
- `ambiguity_summary.tsv`
- `ambiguity_summary_by_gene.tsv`
- `reference_metadata.tsv`
- `majority_vote_baseline.tsv`
- `weighted_consensus_calls.tsv`
- `method_comparison.tsv`
- `method_per_gene.tsv`
- `method_ambiguity_summary.tsv`
- `method_ambiguity_summary_by_gene.tsv`
- `per_gene_gain.tsv`
- `confidence_bin_summary.tsv`
- `confidence_calibration_summary.tsv`
- `abstention_tradeoff.tsv`
- `discordance_tags.tsv`
- `discordance_summary.tsv`
- `benchmark_metadata.json`
- `consensus_runtime_weights.json`

### Figures
- `figure_2_accuracy_comparison.svg`
- `figure_3_per_gene_gains.svg`
- `figure_4_confidence_calibration.svg`
- `figure_5_abstention_tradeoff.svg`
- `figure_6_discordance_taxonomy.svg`
- `figure_7_confidence_weights.svg`
- `captions.md`

## Metric Definitions

- Primary resolution: two-field allele resolution
- Secondary ambiguity-aware layers: three-field, G-group, and P-group summaries when comparable data are available
- Primary success metric: exact correct allele pair at the configured resolution
- IMGT/HLA version and secondary resolution targets are recorded in `reference_metadata.tsv` and `benchmark_metadata.json`
- Majority-vote baseline: strict plurality across tools within each modality
- Weighted consensus: benchmark-derived per-tool/per-gene/per-modality weights with abstention
- Abstention thresholds:
  - `min_support`: minimum support fraction for the winning weighted allele pair
  - `min_margin`: minimum margin over the runner-up weighted allele pair
- Calibration summary:
  - `brier_score`
  - `expected_calibration_error`

## Discordance Taxonomy

The current benchmark layer emits the following high-level tags:
- `technical_conflict`
- `low_evidence_conflict`
- `dna_rna_discordance`
- `possible_expression_bias`
- `no_evidence`
