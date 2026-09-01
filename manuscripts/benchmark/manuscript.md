# ChampHLA: an auditable multi-caller workflow and multimodal benchmark for short-read HLA typing

## Abstract

Short-read HLA typing varies across sequencing assays and computational
callers. ChampHLA provides a reproducible workflow for running, harmonising,
and auditing modality-appropriate HLA callers. We evaluate exact unordered
two-field HLA-A, HLA-B, and HLA-C calls while preserving missingness and fixed
denominators. The study shows that caller-output consensus is close to
saturation in WES and RNA-seq, whereas WGS is more sensitive to caller panel and
input preparation. A frozen Guarded Champion–Challenger policy preserves
high-agreement consensus and can resolve lower-agreement loci, but its
development performance is not presented as external validation. Learned and
read-evidence extensions did not pass their prespecified retention gates.
ChampHLA is therefore positioned as a research-only, auditable benchmarking and
integration framework pending completion of frozen unseen-subject confirmation.

## Introduction

HLA typing from short reads is difficult because of extreme polymorphism,
reference bias, multimapping, assay-specific coverage, and heterogeneous caller
outputs. Individual tools use different reference and inference strategies, so
their errors and missingness differ. Combining them may improve robustness, but
the value of integration must be assessed against both strong individual tools
and transparent consensus rules.

ChampHLA integrates intended-use callers for WGS, WES, and bulk RNA-seq and
records the evidence behind every consensus decision. This paper asks three
questions: whether a simple consensus protects against weak callers, where
caller disagreement leaves room for selection, and whether more complex
selection methods provide reproducible benefit beyond transparent baselines.

## Methods

The shared harmonisation, method, endpoint, comparator, and truth-firewall
definitions are maintained in `manuscripts/shared/shared_methods.md`.

The current 1000 Genomes analysis is development evidence. Subject-unseen 2014
laboratory-typed 1000G samples are reserved for same-resource confirmation;
they are not described as an independent cohort. HPRC phased assemblies define
the independent WGS validation arm. NCI-60 and E-MTAB-197 are secondary
transfer or independent-library analyses and are not pooled into the primary
three-modality family.

## Results

### Reproducible internal analysis

The frozen two-thirds baseline was correct at 328 of 390 WES loci
[RESULT:DEV_TT_WES_BASE]. Guarded CC was correct at 369 loci
[RESULT:DEV_TT_WES_GCC]. In RNA-seq the corresponding counts were 296 of 321
[RESULT:DEV_TT_RNA_BASE] and 306 of 321 [RESULT:DEV_TT_RNA_GCC]. These are
development results. The interpretation is resolution of loci on which the
abstaining baseline has insufficient agreement, not improved conditional
accuracy on protected high-agreement calls.

### WGS input-integrity finding

Historical WGS accuracy estimates are excluded because the source slice omitted
HLA alternate-contig mappings and mates outside the chromosome-6 interval. WGS
performance will be reported only after the full-CRAM, mate-aware five-caller
audit and frozen confirmation are complete.

### Negative method ablations

RefFormer failed its prespecified candidate-reranker gate
[RESULT:REFFORMER_GATE]. EvidenceGatedCC also failed its development gate
[RESULT:EGCC_GATE]. These negative results constrain the method claim: caller
selection or added read evidence did not reliably improve already saturated
WES/RNA consensus.

## Discussion

The present evidence supports an auditable software and benchmark contribution.
It does not yet support universal Champion–Challenger superiority or completed
three-modality external validation. The abstaining two-thirds baseline is useful
for measuring resolution of consensus failures, but always-call plurality and
stronger learned comparators must remain visible. The study is restricted to
research-only class-I, exact two-field typing; class II, long reads, clinical
safety, FIMM LOH, HED, and survival claims are outside the main manuscript.

## Data and software availability

The submission bundle will contain code, compact fixtures, manifests, fold and
cohort assignments, frozen predictions, source-table checksums, statistical
outputs, and figure-source tables. A stable DOI will be assigned only after the
confirmation and release freezes validate.

