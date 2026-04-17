# HLA Ensemble Schemas

## Headline vs Detailed Outputs

Headline benchmark and consensus outputs are currently scoped to `HLA-A`, `HLA-B`, and `HLA-C`. This applies to benchmark-facing tables, manuscript-ready summaries, and default consensus interpretation.

Detailed per-tool output directories may still contain additional loci such as `DQB1`, `DRB1`, and other tool-native HLA genes. These richer outputs are preserved for inspection, troubleshooting, and future scope expansion, but they are not part of the current headline benchmark claim set unless explicitly stated.

Modality is resolved once at workflow runtime and then propagated explicitly into aggregation, weighting, and consensus layers so samplesheet and directory-based runs behave consistently. Runtime weight consumption supports both the native runtime-weight schema and imported `raw_accuracy` calibration JSON files.

## Unified HLA Call TSV Schema

One row per `sample x modality x tool x gene`.

Required columns:
- `sample`
- `modality`
- `tool`
- `gene`
- `truth_allele1_raw`
- `truth_allele2_raw`
- `allele1_raw`
- `allele2_raw`
- `truth_allele1`
- `truth_allele2`
- `allele1`
- `allele2`
- `truth_allele1_3field`
- `truth_allele2_3field`
- `allele1_3field`
- `allele2_3field`
- `call_status`
- `correct_status`
- `is_callable`
- `is_correct`
- `is_correct_2field`
- `is_correct_3field`
- `is_correct_g_group`
- `is_correct_p_group`
- `match_grade`
- `imgt_hla_version`
- `runtime_hours`
- `max_ram_gb`
- `source_file`

Semantics:
- `call_status`: `callable` or `missing`
- `correct_status`: `correct` or `incorrect`
- `is_callable`: `1` only when both alleles are emitted after normalization at the configured primary resolution
- `is_correct`: `1` only when the normalized allele pair matches truth exactly at the configured primary resolution
- `is_correct_2field`: exact pair match at two-field resolution
- `is_correct_3field`: exact pair match at three-field resolution
- `is_correct_g_group`: exact pair match at G-group level when both truth and call provide comparable G-group alleles; blank otherwise
- `is_correct_p_group`: exact pair match at P-group level when both truth and call provide comparable P-group alleles; blank otherwise
- `match_grade`: one of `missing`, `exact_3field`, `exact_2field`, `g_group`, `p_group`, or `mismatch`
- `imgt_hla_version`: IMGT/HLA release pinned for the benchmark run

## Unified Confidence TSV Schema

One row per `sample x modality x tool x gene` after confidence extraction/calibration.

Required columns:
- `sample`
- `modality`
- `tool`
- `gene`
- `confidence_score`
- `confidence_source`
- `raw_confidence`
- `read_support`

Semantics:
- `confidence_score`: normalized to `[0,1]`
- `confidence_source`: one of `optitype_result_objective`, `long_confidence_table`, `wide_confidence_table`, `json_confidence_table`, or `missing`
- `raw_confidence`: tool-native confidence before normalization, when present
- `read_support`: sidecar read-support count used for saturation-based confidence mapping when direct confidence is unavailable

## Consensus Benchmark Schema

### Majority-vote baseline
One row per `sample x modality x gene`.

Columns:
- `sample`
- `modality`
- `gene`
- `method`
- `truth_allele1`
- `truth_allele2`
- `allele1`
- `allele2`
- `truth_allele1_3field`
- `truth_allele2_3field`
- `allele1_3field`
- `allele2_3field`
- `call_status`
- `is_callable`
- `is_correct`
- `is_correct_2field`
- `is_correct_3field`
- `is_correct_g_group`
- `is_correct_p_group`
- `match_grade`
- `imgt_hla_version`
- `agreeing_tools`
- `contributing_tools`
- `support_fraction`
- `support_margin`
- `discordance_tag`

### Weighted consensus benchmark
Same grouping and core columns as majority vote, plus:
- `total_weight`

Runtime weight inputs currently support two compatible machine-readable forms:
- native runtime schema with `tool_weights` and `gene_weights` containing `final_weight`
- imported calibration schema with `raw_accuracy`, converted at runtime into tool-level and gene-level weights

Method-level ambiguity outputs:
- `method_ambiguity_summary.tsv`
- `method_ambiguity_summary_by_gene.tsv`

Confidence-stratified error outputs:
- `confidence_error_summary.tsv`
- `confidence_error_summary_by_gene.tsv`

Consensus `call_status` values:
- `called`
- `low_confidence`
- `no_call`

Current default consensus gene scope:
- `A`
- `B`
- `C`

Discordance tags currently emitted by the benchmark layer:
- `consensus_call`
- `technical_conflict`
- `low_evidence_conflict`
- `no_evidence`
- `dna_rna_discordance`
- `possible_expression_bias`


## Real-Data Cohort Manifest Schemas

### Truth manifest
One row per truth-backed sample.

Required columns:
- `sample`
- `population`
- `truth_source`
- `acquisition_date`
- `truth_supported_loci`
- `truth_gene_count`

### Sequencing manifest
One row per `sample x modality` availability record.

Required columns:
- `sample`
- `population`
- `modality`
- `data_locator`
- `available`

Semantics:
- `modality`: one of `wgs`, `wes`, `rnaseq`
- `available`: `1` only when the modality is present and benchmarkable for that sample

### Cohort manifest
One row per sample after intersecting truth and sequencing availability.

Required columns:
- `sample`
- `population`
- `include`
- `split`
- `truth_supported_loci`
- `wgs_available`
- `wes_available`
- `rnaseq_available`
- `excluded_reason`

Semantics:
- `include`: `1` only for samples retained in the strict tri-modal benchmark cohort
- `split`: `training`, `validation`, or `holdout` for included samples
- `excluded_reason`: blank for included samples; otherwise records the exclusion rule that removed the sample
