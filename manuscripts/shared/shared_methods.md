# Shared verified methods

## Scope and caller panels

ChampHLA is a research-only workflow for HLA-A, HLA-B, and HLA-C typing from
the evaluated short-read WGS, WES, and bulk RNA-seq designs. Each assay uses a
frozen compatible caller panel. Tool versions, executable artifacts, reference
build, IMGT/HLA release, configuration, and neutral disposition are recorded in
`configs/comparator_manifest.json`. A tool is not removed because its observed
performance is unfavorable. The candidate-set oracle is reported separately
because it is not deployable.

## Harmonisation and partial calls

Native outputs are retained and checksummed. Alleles are normalized to
canonical two-field nomenclature without changing the source text. Genotypes
are represented as unordered allele pairs. An explicit repeated allele or a
native homozygosity declaration produces a homozygous pair. A single allele
with a missing second allele, including the native token `-`, produces
`call_status=partial`: allele 1 is retained and allele 2 remains empty. The
workflow never converts a partial call into a homozygous genotype. Partial and
missing predictions are incorrect for fixed-denominator genotype accuracy;
partial-call frequency is reported by caller, assay, and locus.

## Primary algorithm

`SimplePluralityLex` is the stable machine identifier; scientific text calls it
“pair-level plurality consensus.” `MajorityVote` is accepted only as a legacy
alias. Each intended-use caller that returns a complete canonical unordered
allele pair contributes one equal vote. Partial and missing calls do not vote.
The pair with the largest support is selected, and lexicographic pair ordering
resolves an exact top-support tie. If no caller supplies a complete pair, the
result is `no_evidence`.

The audit record includes the complete, partial, and missing caller counts; top
support and its fraction among complete calls; tie indicator and tied pairs;
supporting callers; method version; and native source hashes. The rule is a
transparent operating default, not a novel voting algorithm.

## Endpoints and inference

The primary endpoint is exact unordered two-field genotype correctness at the
fixed eligible-truth denominator. Only truth loci resolving uniquely to a
two-field pair enter this denominator. Missing, partial, abstained, and
`no_evidence` predictions are incorrect. Secondary endpoints are callability,
called-only genotype accuracy, allele accuracy, tie frequency, partial calls,
missingness, per-gene and population results, and runtime/resources.

Every deployable comparator is paired with plurality on the same eligible
loci. Repeated libraries and tissues are clustered by donor. Holm correction is
applied separately within each modality's frozen deployable-comparator family.
Simultaneous subject-clustered bootstrap intervals assess the prospective
noninferiority margin. The caller selected as best after outcomes are observed
is descriptive; all prespecified caller comparisons remain in the corrected
family. Donor-independent cohorts and new-library overlapping-donor cohorts are
reported as separate strata and are never pooled into one headline estimate.

The three machine-readable gates are `no_comparator_holm_superior`, a
superiority screen; `consensus_noninferior_2pp`, the prospective simultaneous
margin assessment; and `three_modality_claim_ready`, which additionally
requires valid WGS, WES, and RNA-seq, frozen prospective design, sample-size
minimums, complete comparator records, and the WGS production audit. Failure to
reject superiority is not interpreted as equivalence.

## Truth firewall and rerun integrity

Predictions are created without truth-bearing columns. Code, configurations,
panels, native-source hashes, and predictions are frozen and revalidated before
truth is joined once. The join refuses overwrite and requires exact prediction
and truth locus sets. Registry entries then resolve to unique canonical result
rows and validate identity fields, counts, accuracy, validity, evidence role,
and source hashes.

The corrected rerun starts from full, mate-aware WGS CRAM extraction retaining
the chr6 MHC and HLA-A/B/C alternate contigs, plus end-to-end WES and RNA-seq
runs through the same parser. All lanes fail closed on a process error, empty
native output, resource mismatch, or incomplete caller-by-locus record set.
Before truth joining, WGS additionally requires a complete production audit,
parser round trips, environment/reference manifests, and named stratified human
review. The existing WES and RNA-seq eligible-truth denominators can change only
through a separately documented truth-eligibility defect.
