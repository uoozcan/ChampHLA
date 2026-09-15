# ChampHLA: auditable plurality consensus for cross-assay short-read HLA typing

## Abstract

HLA genotyping from short-read sequencing is sensitive to assay design,
reference representation, and caller-specific failure modes. We present
ChampHLA, an auditable workflow that executes compatible caller panels,
preserves native evidence, harmonises complete, partial, and missing calls, and
applies ChampHLA-Consensus, a pair-level plurality consensus, as a
selection-free default. The
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
callers [@Robinson2020; @Szolek2014]. A tool that performs well in exome data may not lead in RNA-seq or
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
source artifact rather than trusting a checksum alone (Figure 1).

![Figure 1. Workflow and truth firewall. Truth-free WGS, WES, and RNA-seq inputs flow through assay-compatible callers, native-output preservation, call-state harmonisation, pair-level plurality, and an immutable prediction freeze. Segregated truth can enter only the one-time non-overwriting join; the dashed path marks forbidden truth access during prediction.](../figures/main/Figure1_workflow_firewall.png)

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

The corrected same-resource version of this table is generated from the
registry and same-resource evaluation output; later donor-independent results
are registered and reported as a separate prospective evidence layer. The
table will include fixed-denominator accuracy, call rate,
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

The same-resource benchmark uses the 1000 Genomes resource [@1000Genomes2015].
The independent WGS lane is based on the Human Pangenome Reference Consortium
resource [@Liao2023], and the small independent RNA-seq feasibility lane uses
NCI-60 material [@Adams2005]. These source citations identify the cohorts; they
do not convert an unfinished lane into validation evidence.

The imported discovery program tracks dataset lanes through inventory,
metadata verification, truth verification, executability, prediction freeze,
and evaluation. HPRC Release 2 supplies a feasible truth-blind WGS roster, but
evaluation awaits independently audited phased-assembly truth. The strict
NCI-60 RNA panel is exploratory and under the confirmation target. CCLE WES
overlap is metadata-verified but lacks complete accuracy-eligible A/B/C truth.
Controlled candidates remain blocked until authorization, mapping, truth
independence, and data-use terms are documented.

Figure 2 maps these evidence roles and current lane states without presenting
an aspirational cohort as completed or pooling across independence strata.

![Figure 2. Cohort and evidence map. Dataset roles, independence strata, available or target scale, truth source, and current lane state are shown separately. Color encodes evidence role rather than performance; same-resource, donor-independent, exploratory, and invalid evidence are not pooled.](../figures/main/Figure2_cohort_evidence_map.png)

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

The corrected same-resource WGS, WES, and RNA-seq lanes are rerun end to end
after the partial-call correction. Remote WGS CRAM extraction permits two
five-hour attempts only for allowlisted transport failures. Each attempt has an
immutable staging-ledger row and separate stderr; unsuccessful attempts and
their partial BAMs are quarantined and excluded. Because full CRAMs are
streamed rather than retained, their manifest checksums are recorded as
upstream-declared identities, not as locally verified values; extracted HLA
BAMs and indexes are hashed locally. Execution fails closed for missing
processes, empty outputs, resource
mismatches, or incomplete caller/locus records. WGS additionally requires a
complete native-output audit and named stratified human review before truth can
be joined. Code, configurations, panels, predictions, and provenance are then
frozen; truth is joined once; registry rows, tables, prose, and release hashes
are regenerated from the evaluated artifacts.

Caller databases are recorded component by component rather than assigned a
single repository-wide IPD-IMGT/HLA release. The frozen panel intentionally
contains different caller-specific releases, including a heterogeneous
SpecHLA bundle; the versioned caller/reference attestation is authoritative.
The frozen assay-compatible panels comprise HLA-HD [@Kawaguchi2017], Kourami
[@Lee2018], OptiType [@Szolek2014], POLYSOLVER [@Shukla2015], SpecHLA
[@Wang2023], T1K [@Song2023], and arcasHLA [@Orenbuch2020]. Each caller is used
only in the modalities declared in the frozen panel manifest.

## Discussion

The current evidence motivates pair-level plurality as ChampHLA's recommended
default for WES and bulk RNA-seq development data. That recommendation is about
transparency and robustness to caller choice, not novelty of the voting rule.
The data also show why a cross-assay conclusion cannot be inferred from WES and
RNA-seq alone: input construction and caller compatibility can dominate the WGS
comparison.

The revised architecture turns that limitation into a testable boundary. The
corrected same-resource WGS, WES, and RNA-seq lanes can establish benchmark
readiness, but cannot establish prospective noninferiority. The predeclared
two-percentage-point noninferiority conclusion is reserved for genuinely
prospective, donor-independent cohorts. If valid corrected WGS underperforms,
plurality will be described as a strong WES/RNA-seq default with a demonstrated
WGS limitation. If WGS integrity fails again, the paper will remain a software
and benchmark report without a finalized cross-assay performance claim.

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

## References

Adams S, Robbins F-M, Chen D, et al. HLA class I and II genotype of the NCI-60
cell lines. *Journal of Translational Medicine*. 2005;3:11.
doi:10.1186/1479-5876-3-11.

The 1000 Genomes Project Consortium. A global reference for human genetic
variation. *Nature*. 2015;526:68–74. doi:10.1038/nature15393.

Kawaguchi S, Higasa K, Shimizu M, Yamada R, Matsuda F. HLA-HD: an accurate HLA
typing algorithm for next-generation sequencing data. *Human Mutation*.
2017;38:788–797. doi:10.1002/humu.23230.

Lee H, Kingsford C. Kourami: graph-guided assembly for novel human leukocyte
antigen allele discovery. *Genome Biology*. 2018;19:16.
doi:10.1186/s13059-018-1388-2.

Liao W-W, Asri M, Ebler J, et al. A draft human pangenome reference. *Nature*.
2023;617:312–324. doi:10.1038/s41586-023-05896-x.

Orenbuch R, Filip I, Comito D, Shaman J, Pe'er I, Rabadan R. arcasHLA: high
resolution HLA typing from RNAseq. *Bioinformatics*. 2020;36:33–40.
doi:10.1093/bioinformatics/btz474.

Robinson J, Barker DJ, Georgiou X, Cooper MA, Flicek P, Marsh SGE. IPD-IMGT/HLA
Database. *Nucleic Acids Research*. 2020;48:D948–D955.
doi:10.1093/nar/gkz950.

Shukla SA, Rooney MS, Rajasagi M, et al. Comprehensive analysis of cancer-associated
somatic mutations in class I HLA genes. *Nature Biotechnology*.
2015;33:1152–1158. doi:10.1038/nbt.3344.

Song L, Bai G, Liu XS, Li B, Li H. T1K: efficient and accurate KIR and HLA
genotyping with RNA-seq. *Genome Research*. 2023;33:923–931.
doi:10.1101/gr.277585.122.

Szolek A, Schubert B, Mohr C, Sturm M, Feldhahn M, Kohlbacher O. OptiType:
precision HLA typing from next-generation sequencing data. *Bioinformatics*.
2014;30:3310–3316. doi:10.1093/bioinformatics/btu548.

Wang S, Wang M, Chen L, Pan G, Wang Y, Li S. SpecHLA enables full-resolution HLA
typing from sequencing data. *Cell Reports Methods*. 2023;3:100589.
doi:10.1016/j.crmeth.2023.100589.
