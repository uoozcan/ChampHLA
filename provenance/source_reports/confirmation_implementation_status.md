# Implementation status

The isolated confirmation workflow is implemented and locally verified. No
existing ChampHLA, RefFormer, MetaConsensus, EGCC, or manuscript file was
modified.

## Completed gates

- Frozen intended-use caller panels and two-thirds baseline.
- Truth-blind `TwoThirdsGuardedCC` without the truth-derived homozygosity rule.
- Truth-bearing runtime-column rejection and checksum-protected join-once truth
  firewall.
- Exact subject-cluster sign-flip inference, subject bootstrap, Holm correction,
  per-gene harm, plurality non-regression, secondary comparator, candidate-oracle,
  capacity, and power-projection outputs.
- Alias/relative overlap crosswalk, external cohort schemas, source provenance,
  caller-native WGS audit, and 50-record manual-review sheet.
- Independent evaluator and discovery reproduction gate.
- Thirty-two passing unit/integration tests, including caller-order
  invariance, parser round trips, truth-column rejection, immutable freezes,
  code-only freeze scope, HPRC selection, truth-free IHWG registry sanitation,
  and strict five-caller pilot collection.

## Discovery results only

| Modality | Two-thirds baseline | Guarded CC | Difference |
|---|---:|---:|---:|
| WES | 328/390 | 365/390 | +9.49 points |
| RNA-seq | 296/321 | 304/321 | +2.49 points |
| WGS | 79/411 | 205/411 | +30.66 points |

The original 30-subject WGS native outputs have now passed automated
round-trip parsing (450/450 expected caller/locus records) and a 50-record
stratified manual review. That audit also detected and fixed an isolated audit
parser bug that could mistake HLA-DMA/DOB tokens for class-I alleles. However,
the old input BAMs were built by a chromosome-6 interval stream that omitted
HLA alternate-contig reads and discarded mates outside the slice. Therefore,
their WGS accuracy remains invalid for confirmation.

A truth-blind full-CRAM repair pilot is defined in
`configs/full_cram_wgs_pilot.tsv`. It uses mate-aware retrieval of the full MHC
interval plus every HLA-A/B/C GRCh38DH contig and then runs exactly HLA-HD,
Kourami, OptiType, SpecHLA, and T1K. This provides a bounded capacity and input
integrity gate without materializing whole-genome BAM/FASTQ files. The software
still cannot label any result `three_modality_confirmed` unless the repaired
WGS audit passes, Guarded CC emits every top call, the truth-blind capacity gate
passes, external minimum sample sizes are reached, and the immutable freeze
validates.

Roihu jobs explicitly initialize CSC's non-interactive Lmod environment and pin
`bio-apps/v202603`, `samtools/1.21`, and `nextflow/25.10.2-standalone`; the first
submission exposed that plain batch shells do not provide `module`. That failed
job exited before reading or writing sequencing data and is retained in the job
history as an environment-gate failure.

The first full-CRAM subject, HG00096, has now passed the extraction gate. Its
mate-aware BAM and index pass SHA-256 verification and `samtools quickcheck`;
the BAM contains 1,640,886 reads, including mappings on chr6 and HLA-A/B/C
alternate contigs. The initial 68-minute stream was preserved after two
samtools long-option compatibility failures, then resumed without downloading
again. The frozen five-caller truth-blind job is submitted as Roihu job 977640.
HG00097 remains gated until the first caller run validates.

The official HPRC Release 2 metadata at commit
`5a939042026331a823a6307fe36a3d7e0188a6e0` yields 166 truth-free,
development-nonoverlapping subjects with two assembly haplotypes and paired
Illumina WGS, so the WGS target of 120 is feasible. The verified local IHWG
inventory contains only two WES and two study-provenance-valid RNA subjects;
the locked minima of 89 and 130 are not met. The three-modality confirmation
claim is therefore unavailable unless genuinely new, orthogonally typed public
WES/RNA cohorts are identified.

A systematic truth-free search has therefore replaced the earlier ad hoc name
search. The source IHWG map was reduced to 366 unique records containing only
IHW number, primary line name, and ancestry; its sanitized SHA-256 is
`1ddbfb9be35fed17e490167850a30b1d875bcb63723bd6ca780eaa05cefe86f0`.
Roihu job 977650 completed 366 official ENA read-run queries with zero failures.
The fail-closed technical/name screen found only four possible WES subjects and
22 possible RNA-seq subjects across 81 runs. Those are absolute upper bounds
before study-level provenance and overlap exclusions, and are already far below
the locked minima of 89 and 130. Therefore the current public inventory cannot
support the three-modality confirmation. Every matching run remains
`requires_study_level_review`; short names, sample titles, and library metadata
alone cannot establish that the sequenced material is the reference line.

The discovery projection at 120 WGS, 89 WES, and 150 RNA subjects is promising,
but is deliberately marked invalid until WGS repair passes. External truth has
not been opened or used for policy changes.

## Final truth-free raw CC policy

The historical evaluator mixed inference code with truth-bearing result rows.
The confirmation copy now freezes the same development reliability/calibrated-
confidence weighting formula into `models/frozen_raw_cc_policy.json` and runs a
standalone inference scorer that rejects truth columns. Against historical
outer-fold raw-CC outputs, the final frozen policy agrees at 383/390 WES,
318/321 RNA, and 409/411 WGS loci; exact identity is not expected because the
historical calls used fold-specific fitted weights.

As a development-only diagnostic, the final bundle produces Guarded-CC accuracy
of 369/390 WES, 306/321 RNA, and 205/411 WGS versus two-thirds baselines of
328/390, 296/321, and 79/411. It has no losses relative to the baseline by
construction and passes the plurality non-regression diagnostic. These are
resubstitution/discovery results, not confirmation evidence, and the evaluator
therefore reports `discovery_only` and `headline_retained=false`.

Revalidation with the stricter external-freeze rule correctly rejects the old
development bundle: four loci in the invalid legacy WGS input have neither an
intended-use caller call nor a raw CC call, so Guarded CC cannot emit the
mandatory top call. The legacy freeze is forbidden for external use. The freeze
implementation was also narrowed to code/config/model/template paths plus
explicit named inputs; it no longer sweeps development truth tables into the
code checksum manifest.
