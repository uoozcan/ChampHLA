# ChampHLA: auditable plurality consensus for cross-assay short-read HLA typing

## Abstract

HLA genotyping from short-read sequencing is sensitive to assay design,
reference representation, and caller-specific failure modes. We present
ChampHLA, an auditable workflow that executes compatible caller panels,
preserves native evidence, harmonises complete, partial, and missing calls, and
applies pair-level plurality consensus as a selection-free default. The
plurality rule itself is intentionally simple and is not claimed as a novel
voting algorithm. The contributions are its reproducible implementation, a
truth-blind evidence architecture, cross-assay benchmarking, and explicit tests
of where simple consensus succeeds or fails. Development evidence supports the
use of plurality for WES and bulk RNA-seq, while historical WGS results were
invalidated by an input-extraction defect. A corrected three-modality rerun is
therefore a prerequisite for the final performance claim. ChampHLA is intended
for research use in HLA-A/B/C typing, not clinical decision-making.

## Introduction

Class-I HLA typing is a demanding short-read inference problem. Extensive
polymorphism, sequence similarity among loci, alternative reference contigs,
multimapping, and uneven assay coverage can produce errors that differ across
callers. A tool that performs well in exome data may not lead in RNA-seq or
whole-genome data. Selecting one caller globally therefore risks hiding an
assay-specific tradeoff, while a complex integration model can add instability
without improving an already strong consensus.

ChampHLA is organized around a simpler operational question: can a transparent
pair-level vote provide a robust default without choosing a winner in advance?
The workflow runs only modality-compatible tools, preserves their native
outputs and provenance, and reduces complete calls to canonical unordered
two-field HLA-A/B/C pairs. Every complete caller receives equal weight. The
highest-support pair is emitted and exact ties are resolved deterministically.
This rule is deliberately unsurprising; scientific value comes from making the
entire comparison traceable and from locating the assay boundary at which the
default stops being reliable.

The study addresses three linked questions. First, how much does individual
caller performance vary by assay? Second, how closely does selection-free
plurality track the strongest observed caller without selecting that caller?
Third, do more elaborate integration policies reliably improve on plurality?
The comparison is guarded by a truth firewall, fixed denominators, donor-level
inference, modality-specific multiplicity families, and explicit invalid-data
states.

## Results

### An auditable workflow separates evidence generation from evaluation

ChampHLA emits one record for each caller, subject, assay, and locus, including
the native source hash and normalized call state. A native missing second allele
remains partial instead of being converted to homozygosity. Partial and missing
outputs cannot vote, and each plurality result reports contributing callers,
support, ties, missingness, parser version, and source hashes. Prediction files
contain no truth and must be frozen before the one-time truth join. The result
registry resolves every reported result to one canonical artifact row and its
source artifact rather than trusting a checksum alone.

### Individual-caller performance is assay dependent

The development comparison showed that the caller with the highest observed
accuracy differed by modality. Because “best observed caller” is selected after
viewing outcomes, it is reported descriptively; confirmatory inference retains
every frozen deployable caller in the modality-specific family. Exhaustive
caller results and the candidate-generation oracle are assigned to the
supplement so the main comparison does not turn a post-selection summary into a
prospective claim.

### Pair-level plurality is the recommended selection-free default

In the existing development artifacts, pair-level plurality produced 365
correct WES genotypes among 390 eligible loci [RESULT:DEV_WES_CONSENSUS] and 306
correct RNA-seq genotypes among 321 eligible loci
[RESULT:DEV_RNA_CONSENSUS]. These denominators and the choice to foreground
plurality were visible before the primary amendment; the estimates are
therefore provisional development evidence rather than prospective
confirmation.

| Modality | Pair-level plurality | Best observed caller | Strongest integration comparator | Interpretation |
|---|---:|---|---|---|
| WES | 365/390 [RESULT:DEV_WES_CONSENSUS] | generated in the corrected head-to-head | generated in the corrected head-to-head | provisional development evidence |
| bulk RNA-seq | 306/321 [RESULT:DEV_RNA_CONSENSUS] | generated in the corrected head-to-head | generated in the corrected head-to-head | provisional development evidence |
| WGS | excluded | awaiting corrected rerun | awaiting corrected rerun | no valid performance result |

The final version of this table is generated from the registry and prospective
evaluation output. It will include fixed-denominator accuracy, call rate,
paired comparator-minus-plurality difference, simultaneous confidence interval,
and Holm-adjusted inference. No WGS estimate enters the main table until the
full-CRAM production audit passes.

### More complex integration did not establish a reliable advantage

The historical ablations include raw Champion–Challenger, Guarded CC,
MV-floor, MetaConsensus, RefFormer, EvidenceGatedCC, and weighted variants. The
number of retained integration and historical-ablation methods is generated
from the comparator manifest rather than typed into the manuscript. RefFormer
failed its frozen retention decision [RESULT:REFFORMER_GATE], and
EvidenceGatedCC failed its development decision [RESULT:EGCC_GATE]. These
results support preserving the methods as informative ablations, not promoting
them to a headline algorithm. Complete head-to-head and baseline-ladder results
remain visible in the supplement.

### The historical WGS lane identified a technical boundary, not a method result

The historical WGS input was constructed from a chromosome-restricted slice
that omitted HLA alternate-contig alignments and mates outside the selected
interval. All performance records derived from that lane are marked invalid.
They remain discoverable for diagnosis but are excluded from gates, pooled
estimates, primary figures, and performance language. The corrected WGS lane
starts from the full CRAM, keeps relevant alternate contigs, preserves mates,
and runs the frozen five-caller panel before any truth access.

### Three explicit gates replace an ambiguous retention decision

The evaluation reports three outcomes. `no_comparator_holm_superior` asks
whether any frozen comparator has a Holm-significant positive advantage over
plurality; this is a superiority screen and does not establish equivalence.
`consensus_noninferior_2pp` evaluates the predeclared prospective margin using
simultaneous donor-clustered upper confidence bounds. Finally,
`three_modality_claim_ready` requires all modalities to be valid, prospective,
adequately sized, complete across comparator records, and technically audited.
The current development analysis does not satisfy the final gate.

## Methods

### Study design and evidence roles

Every dataset is assigned a role before predictions are evaluated:
development, same-resource confirmation, independent validation, exploratory,
or invalid. Donor-independent cohorts and new-library cohorts that reuse donors
are analyzed separately. Repeated tissues and libraries share a donor cluster.
Computational, imputed, circular, and RNA-informed truth sources cannot enter an
accuracy denominator for the affected assay.

The imported discovery program tracks dataset lanes through inventory,
metadata verification, truth verification, executability, prediction freeze,
and evaluation. HPRC Release 2 supplies a feasible truth-blind WGS roster, but
evaluation awaits independently audited phased-assembly truth. The strict
NCI-60 RNA panel is exploratory and under the confirmation target. CCLE WES
overlap is metadata-verified but lacks complete accuracy-eligible A/B/C truth.
Controlled candidates remain blocked until authorization, mapping, truth
independence, and data-use terms are documented.

### Harmonisation and plurality

The full algorithm, call-state contract, endpoints, and truth firewall are
defined in `manuscripts/shared/shared_methods.md`. Briefly, each complete
canonical pair receives one equal vote from its caller; partial and missing
calls receive none. Maximum support wins, lexicographic ordering resolves an
exact top tie, and absence of a complete pair returns `no_evidence`.
`SimplePluralityLex` remains the machine identifier. “Pair-level plurality
consensus” is the scientific name, and `MajorityVote` is a legacy input alias.

### Outcomes and statistical analysis

The primary outcome is exact unordered two-field genotype accuracy at the
fixed eligible-truth denominator. Secondary outcomes include callability,
called-only accuracy, allele accuracy, ties, partial calls, missingness,
per-gene and population results, and resource use. Each method is compared with
plurality on paired loci. Tests and simultaneous intervals cluster by donor;
Holm adjustment is performed within each frozen modality-specific deployable
family. The candidate oracle is descriptive and separate.

### Corrected rerun and release control

WGS, WES, and RNA-seq are rerun end to end after the partial-call correction.
Execution fails closed for missing processes, empty outputs, resource
mismatches, or incomplete caller/locus records. WGS additionally requires a
complete native-output audit and named stratified human review before truth can
be joined. Code, configurations, panels, predictions, and provenance are then
frozen; truth is joined once; registry rows, tables, prose, and release hashes
are regenerated from the evaluated artifacts.

## Discussion

The current evidence motivates pair-level plurality as ChampHLA's recommended
default for WES and bulk RNA-seq development data. That recommendation is about
transparency and robustness to caller choice, not novelty of the voting rule.
The data also show why a cross-assay conclusion cannot be inferred from WES and
RNA-seq alone: input construction and caller compatibility can dominate the WGS
comparison.

The revised architecture turns that limitation into a testable boundary. If
the corrected prospective WGS, WES, and RNA-seq lanes all satisfy the margin,
the conclusion will be that plurality was not materially exceeded across the
evaluated HLA short-read assays. If valid WGS fails, plurality will be described
as a strong WES/RNA-seq default with a demonstrated WGS limitation. If WGS
integrity fails again, the paper will remain a software and benchmark report
without a finalized cross-assay performance claim.

The scope is restricted to research-only HLA-A/B/C typing from the evaluated
short-read designs. Class II loci, long reads, single-cell pseudobulk,
transplant outcomes, loss of heterozygosity, HLA evolutionary divergence, and
clinical safety are outside the primary claim.

## Data and software availability

The release bundle will include source code, compact fixtures, comparator and
environment manifests, dataset roles, native-output checksums, frozen
predictions, one-time join evidence, registry-backed tables, statistical
outputs, and claim-audit records. A DOI-ready bundle is released only after the
author amendment, corrected rerun, WGS technical and human-review gates,
cross-platform tests, manuscript audits, and clean repository freeze all pass.
