# HLA Ensemble Schemas

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

Discordance tags currently emitted by the benchmark layer:
- `consensus_call`
- `technical_conflict`
- `low_evidence_conflict`
- `no_evidence`
- `dna_rna_discordance`
- `possible_expression_bias`
