# ChampHLA HLA-ground-truth dataset assessment

Audit date: 2026-09-07. This is an inventory and pilot-preparation result. No raw reads or controlled files were downloaded, no access application was submitted, and no HLA predictions were run.

## Decision summary

The registry is structurally complete and fail-closed, but the full plan's data-readiness acceptance gate is not yet met.

| Candidate or component | Status | Evidence-backed disposition |
|---|---|---|
| NCI-60 RNA, PRJNA433861 | **verified** | 60 public paired RNA runs; 58 map to SBT-table lines. Eleven canonical lines have complete strict two-field A/B/C truth. These can form a partial RNA pilot after the normal prediction freeze, but do not meet the 12-line target. |
| Original NCI-60 WES | **rejected** | The raw endpoint cited by the WES publication no longer resolves. Processed CellMiner variant calls are not raw WES input. |
| CCLE WES overlap, PRJNA523380 | **verified metadata; conditional pilot** | Three NCI-60 lines have public paired WES: K-562, SF-295, and HS 578T. All three have incomplete or ambiguous strict A/B/C SBT truth, so none can enter the primary accuracy denominator yet. |
| AFGR MKK RNA | **conditional** | ENCODE resolves 166 MKK RNA experiments and their file accessions. The controlled EGA HLA file has not been opened, so exact truth mapping remains gated. |
| AFGR MKK WGS | **conditional** | The AFGR publication reports 164 genomes, but the audited ENCODE AFGR collection contains RNA-seq and ATAC-seq, not WGS. A raw WGS repository accession is unresolved. |
| FNLCR PBMC RNA, phs003177 | **conditional** | Controlled 96-donor orthogonal Sanger candidate. The original calls must be distinguished from the ten RNA-informed reconciliations; the latter are circular for RNA evaluation. |
| DICE RNA, phs001703/PRJNA494278 | **conditional** | Controlled multi-cell-type candidate. The exact DNA HLA file, donor mapping, relationships, layout, and access terms need verification. |
| Transplant WES, phs003394/PRJNA1006802 | **conditional data; rejected truth claim** | The study has 576 germline WES samples, but “HLA-matched” is not sample-level allele truth. It does not close the WES minimum unless the controlled manifest proves independent two-field A/B/C genotypes and unique donor counts. |
| E-MTAB-197 | **conditional** | Public LCL RNA plus reported PCR-SSOP A/B/C truth. Exact FASTQ mapping and CEU/HapMap overlap with ChampHLA development donors remain unresolved. |
| PRJEB6763/IHWG | **conditional** | Public targeted long-range-PCR MHC data from about 95 reference lines. It is a targeted-MHC stress set, not WES, and cannot contribute to the WES capacity requirement. |
| GSE120221 | **conditional** | Public 10x single-cell data with reported PCR/NGS HLA truth. It belongs in an exploratory pseudobulk lane until a 10x adapter is implemented and validated. |

The definitive WES determination is therefore **no verified cohort currently closes the 89-subject WES gap**. The controlled transplant study remains a lead, not evidence that the gap is closed.

## Prepared pilots

`pilot_manifest.tsv` is truth-free and organized as one row per subject/modality/locus, matching the existing `cohort/subject/modality/gene/source` contract.

- **NCI60_PUBLIC_PILOT:** eleven non-derivative lines with public RNA and complete strict A/B/C SBT truth are marked `ready`. SF-295 is retained as the twelfth, tissue-diversifying dual-modality stress line; its RNA and WES rows are `blocked_truth_ambiguity` and cannot be scored in the primary denominator. The planned 12-line matched public WES/RNA pilot is thus prepared but not currently executable as a strict accuracy pilot.
- **AFGR_MKK_GATED_PILOT:** twelve donors have RNA and WGS rows. RNA is `blocked_truth_crosswalk`; WGS is `blocked_raw_accession`. It may execute only when a sample-level EGA truth crosswalk, access terms, and raw WGS accessions are verified.

Selection does not inspect HLA allele values. NCI-60 prioritizes complete-call eligibility, removes derivative duplicates, and uses tissue diversity for the blocked twelfth stress line. Allele-diversity optimization is intentionally deferred because using truth alleles to select the pilot would weaken the truth firewall.

## Acceptance audit

The machine-readable `audit_summary.json` reports:

| Gate | Result |
|---|---:|
| Public RNA pilot ready | FAIL (11/12 strict subjects) |
| Public WES pilot ready | FAIL |
| Public matched WES/RNA pilot ready | FAIL |
| Matched WGS/RNA candidate accession-resolved | FAIL |
| Verified cohort closes WES minimum | FAIL |
| Overlap classification complete | PASS |
| Relatedness classification complete | FAIL (AFGR Coriell pedigree audit pending) |
| Truth circularity classification complete | PASS |

An overall `passed: false` is the intended fail-closed outcome. It prevents inventory breadth from being confused with confirmation readiness.

## Truth and overlap rules

- Primary truth is an unordered exact two-field HLA-A/B/C genotype derived from direct Sanger, SSOP, targeted NGS, or a multi-laboratory reference consensus.
- Original higher-resolution strings, suffixes, ambiguity codes, null alleles, and disputed calls must be preserved. The exact IMGT/HLA release is not yet present in this workspace and must be pinned before normalization.
- Computational or imputed calls are inventory-only and have `accuracy_eligible=no`. RNA-informed reconciled labels are excluded from RNA accuracy denominators.
- Donor-independent cohorts and new-library/overlapping-donor cohorts receive equal prominence but separate estimates. They are never pooled into a single accuracy estimate.
- Repeated tissues and libraries remain visible but are clustered by donor. Related or derivative cell lines are flagged; a single canonical representative is used in the primary panel.
- Class II, partial-locus truth, targeted-MHC data, and pseudobulk single-cell data remain separate exploratory lanes.

## Comparison design

`comparison_design.json` freezes the intended comparison without executing it. Pair-level plurality (`SimplePluralityLex`) is primary. It retains Guarded CC, raw CC, MV-floor, nested-CV best-tool selection, MetaConsensus, simple two-thirds consensus, the separately reported candidate-set oracle, and every compatible component caller: HLA-HD/Kourami/OptiType/SpecHLA/T1K for WGS; HLA-HD/OptiType/POLYSOLVER/SpecHLA/T1K for WES; and arcasHLA/HLA-HD/OptiType/T1K for RNA.

Required outputs are fixed-denominator genotype accuracy, allele accuracy, call/abstention rate, partial-call frequency, tie frequency, caller agreement, generation-versus-selection error, per-locus results, and runtime/resources. Inference uses paired method differences, simultaneous donor-clustered confidence intervals, and Holm correction within each modality's frozen deployable-comparator family. The candidate oracle and the caller selected as best after observing outcomes are descriptive, not members of a post-selection confirmatory test.

## Access instructions and gates

Public accessions may be inspected through ENA/SRA, ENCODE, GEO, and BioStudies. This inventory stores official metadata and direct public FASTQ URLs only for the NCI/CCLE crosswalk; it does not fetch those files.

Controlled candidates require an authorized institutional environment:

1. Request AFGR truth through EGA datasets `EGAD00001011379` and `EGAD00010002577`, then verify permitted processing and the donor-level file manifest.
2. Review dbGaP study manifests for `phs003177`, `phs001703`, and `phs003394` under the approved data-use agreement. Do not move files into an unapproved cloud or local environment.
3. Record only accession/mapping evidence in this registry until authorization exists. Keep truth files outside prediction input paths and join only after predictions are checksummed.

## Unresolved and rejection log

- **NCI-60 WES:** rejected as a public raw-data route because the historical endpoint is unavailable. Reopen only if NCI restores raw BAM/FASTQ access or supplies a replacement accession.
- **CCLE overlap:** exact raw access is resolved, but all three overlapping lines fail strict complete A/B/C truth. They remain method stress cases, not accuracy cases.
- **AFGR WGS:** resolve an official study/file accession and confirm the same donor identifiers used by ENCODE RNA and EGA truth.
- **AFGR relatedness:** audit Coriell pedigrees before selecting twelve independent donors.
- **FNLCR:** obtain original Sanger calls and isolate the ten disputed/reconciled subjects.
- **DICE:** identify the exact DNA HLA truth file and map each sorted-cell library to the typed donor.
- **Transplant WES:** prove unique donor counts and sample-level two-field A/B/C truth; otherwise reject it as a ground-truth cohort.
- **E-MTAB-197:** extract accession-level FASTQ and SSOP identifiers and run the existing ChampHLA alias/relative overlap logic.
- **PRJEB6763:** correct any downstream record calling the assay “WES”; validate reference-line aliases, null alleles, and run mapping.
- **GSE120221:** implement and test a 10x-to-pseudobulk adapter before execution.

## Recommended next discovery targets

1. Resolve the AFGR WGS accession through the publication authors/data-coordination record; this is the most direct route to a large ancestry-diverse matched WGS/RNA panel.
2. Inspect the authorized `phs003394` phenotype manifest specifically for allele-level HLA-A/B/C fields. This gives the fastest definitive controlled-data answer for WES capacity.
3. Crosswalk public GeT-RM and IHIW/Coriell reference panels against SRA/ENA WES projects, requiring study-level provenance rather than cell-line name matches.
4. Audit E-MTAB-197 as a new-library overlap cohort and FNLCR as the highest-value controlled independent RNA cohort.

## Reproduction

The source builder consumes already-downloaded official metadata snapshots and never downloads reads:

```bash
PYTHONPATH=src python3 scripts/build_dataset_discovery_sources.py \
  --nci-hla-html /path/to/PMC555742.html \
  --nci-rna-ena /path/to/PRJNA433861.read_run.tsv \
  --ccle-wes-ena /path/to/PRJNA523380.wxs.read_run.tsv \
  --afgr-supplement-zip /path/to/afgr-media-2.zip \
  --afgr-encode-json /path/to/encode-afgr.json \
  --development discovery/inputs/harmonized_calls_truth_blind.tsv \
  --output-dir external/dataset_discovery
```

Validate all four interfaces and regenerate the acceptance summary:

```bash
audit_dataset_discovery_registry \
  --datasets external/dataset_discovery/dataset_registry.tsv \
  --truth external/dataset_discovery/truth_registry.tsv \
  --crosswalk external/dataset_discovery/sample_crosswalk.tsv \
  --pilot external/dataset_discovery/pilot_manifest.tsv \
  --output external/dataset_discovery/audit_summary.json
```

The audit command exits `2` while any readiness gate is false. That is expected for this inventory snapshot.

## Evidence sources

- AFGR study and supplement: <https://pmc.ncbi.nlm.nih.gov/articles/PMC10659267/>
- MKK HLA study: <https://www.nature.com/articles/s41591-024-02944-5>
- AFGR HLA data: <https://ega-archive.org/datasets/EGAD00001011379>
- NCI-60 SBT truth: <https://pmc.ncbi.nlm.nih.gov/articles/PMC555742/>
- NCI-60 WES publication: <https://pmc.ncbi.nlm.nih.gov/articles/PMC4102467/>
- NCI-60 RNA: <https://www.ncbi.nlm.nih.gov/bioproject/PRJNA433861>
- CCLE WES: <https://www.ncbi.nlm.nih.gov/bioproject/PRJNA523380>
- FNLCR: <https://pmc.ncbi.nlm.nih.gov/articles/PMC9883133/>
- DICE: <https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/study.cgi?study_id=phs001703.v6.p1>
- Transplant WES: <https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/study.cgi?study_id=phs003394.v1.p1>
- E-MTAB-197: <https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-197>
- PRJEB6763: <https://www.ebi.ac.uk/ena/browser/view/PRJEB6763>
- GSE120221: <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE120221>
