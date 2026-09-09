# Supplementary material: plurality-centered ChampHLA benchmark

## S1. Governance and evidence map

The active manuscript is the plurality-centered benchmark. The conditional
Guarded CC manuscript is archived as superseded, while all historical method
artifacts remain immutable. The retrospective primary amendment records that
method framing and development denominators were visible at selection time; it
cannot authorize prospective confirmation until signed. Every result has an
explicit evidence role, validity state, canonical row key, source artifact, and
checksum in `result_registry.tsv`.

## S2. Complete caller and integration head-to-head

The machine-generated development comparison is stored in
`artifacts/generated/development_evaluation/head_to_head.tsv`. It contains every
deployable method in the modality-specific comparator family, paired
comparator-minus-plurality differences, donor-clustered simultaneous intervals,
unadjusted and Holm-adjusted tests, and missing-row validity flags. The compact
view is `artifacts/generated/consensus_head_to_head.md`.

The caller chosen as best after observing outcomes is explicitly descriptive.
Inference is retained for every frozen caller; it is not rerun only against the
selected winner. A comparator with an absent locus row is marked incomplete,
its comparison is not qualifying, and the prospective gates fail closed.

## S3. Development plurality results

Pair-level plurality called 365 of 390 eligible WES loci correctly
[RESULT:DEV_WES_CONSENSUS]. The two-thirds abstaining baseline called 328 of 390
correctly [RESULT:DEV_TT_WES_BASE], while the archived Guarded CC analysis
recorded 369 of 390 [RESULT:DEV_TT_WES_GCC]. These are development artifacts
created before the partial-call correction and corrected end-to-end rerun.

For bulk RNA-seq, plurality called 306 of 321 eligible loci correctly
[RESULT:DEV_RNA_CONSENSUS]. The two-thirds baseline called 296 of 321 correctly
[RESULT:DEV_TT_RNA_BASE], while the archived Guarded CC artifact also recorded
306 of 321 [RESULT:DEV_TT_RNA_GCC]. These values are provisional development
evidence and do not establish prospective equivalence.

## S4. Baseline ladder and historical integration ablations

The complete comparator manifest is `configs/comparator_manifest.json`.
Integration and historical-ablation counts are generated from its `kind` field.
It retains raw Champion–Challenger, Guarded CC, MV-floor, MetaConsensus,
WeightedConsensus, RefFormer, and EvidenceGatedCC regardless of outcome. Simple
two-thirds consensus remains a coverage-accuracy sensitivity baseline rather
than the primary reference, and the candidate-set oracle remains outside the
deployable family.

RefFormer did not pass its frozen retention gate [RESULT:REFFORMER_GATE].
EvidenceGatedCC did not pass its development gate [RESULT:EGCC_GATE]. NCI-60
MetaConsensus transfer analyses remain exploratory and are not presented as
plurality confirmation [RESULT:META_NCI60_WES] [RESULT:META_NCI60_RNA].

## S5. Tie, partial-call, and missingness sensitivities

`artifacts/generated/development_evaluation/caller_call_status.tsv` reports
complete, partial, and missing records for each caller, assay, and locus.
`secondary_endpoints.tsv` reports plurality call rate, called-only accuracy,
allele accuracy, tie count, post-selection best-caller summary, and the
candidate-generation oracle. The corrected rerun regenerates these tables from
native outputs; legacy zero partial-call counts are not interpreted as evidence
that native partial calls were absent because the historical parser collapsed
the missing-second-allele state.

Sensitivity analyses retain the fixed eligible-truth denominator. Partial and
missing calls remain incorrect for genotype accuracy; allele accuracy may
credit the one retained allele through multiset matching. Tie resolution is
tested for caller-order invariance and reported separately from support.

## S6. Per-gene, population, and external-stratum outputs

`per_gene_results.tsv` and `stratified_results.tsv` are generated directly by
the evaluator. `external_strata_results.tsv` reports donor-independent and
new-library overlapping-donor cohorts separately. The evaluator refuses an
input containing multiple independence strata unless the prospective design
names one primary stratum, and it marks all other strata as not pooled into the
primary estimate. Repeated tissues and libraries use the donor identifier as
the resampling cluster.

## S7. Invalid WGS diagnostics

The legacy WGS lane is retained only as an invalid diagnostic. Pair-level
plurality recorded 164 of 411 loci [RESULT:DEV_WGS_CONSENSUS], the two-thirds
baseline recorded 79 of 411 [RESULT:DEV_WGS_LEGACY_BASE], and Guarded CC
recorded 205 of 411 [RESULT:DEV_WGS_LEGACY_GCC]. MetaConsensus from the same
defective source is likewise invalid [RESULT:META_WGS_INTERNAL]. None of these
records enters a gate, pooled estimate, primary figure, abstract, or WGS
performance conclusion.

The defect was an input slice that omitted HLA alternate-contig mappings and
mates outside the chromosome-restricted interval. The corrected production lane
starts from each full CRAM, performs mate-aware extraction, retains the chr6 MHC
and HLA-A/B/C alternate contigs, runs the frozen five-caller WGS panel, and
requires complete native-output, parser, reference, environment, and named
human-review audits before truth joining.

## S8. Dataset inventory and evidence boundaries

The version-controlled discovery package contains 11 dataset/assay candidates,
10 truth sources, 391 crosswalk rows, and 111 truth-free pilot rows. These are
inventory counts, not performance estimates. Its audit records per-lane states:
inventory, metadata verified, truth verified, executable, prediction frozen,
and evaluated.

HPRC Release 2 has a feasible truth-blind WGS roster but still requires
independent phased-assembly truth generation and audit. NCI-60 RNA is an
exploratory strict pilot; CCLE WES has verified metadata but no complete
accuracy-eligible A/B/C truth. AFGR, FNLCR, DICE, transplant WES, and other
controlled candidates remain blocked until authorization, sample mapping,
truth independence, and access terms are verified. Rejected or circular truth
sources remain visible and cannot enter denominators.

## S9. Runtime and resource reporting

The corrected execution manifest records caller wall time, CPU time, peak
memory, scheduler status, input bytes, output bytes, container digest, command,
and reference checksums. Empty output or an absent scheduler process is an
error, not a missing call. Runtime tables will be generated only from the
corrected frozen runs so historical pilot resources are not mistaken for the
production benchmark.

## S10. Count reconciliation and reproducibility

`decisions/20260908_count_reconciliation.json` records the known historical
discrepancies and the neutral rule for resolving them from corrected canonical
rows. `artifacts/registered_results_long.tsv` is the canonical current result
artifact. Registry validation checks identity and values, while an independent
evaluator must reproduce each headline count before release.

Local verification runs both pytest entry forms, regenerates the registry and
head-to-head tables, validates the discovery inventory, and audits this
supplement together with the main text. Ubuntu and Windows CI exercise the same
suite, including cross-platform POSIX serialization of freeze paths.
