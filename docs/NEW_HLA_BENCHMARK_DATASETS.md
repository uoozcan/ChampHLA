# New HLA Ground-Truth Datasets — Scouting Shortlist

**Purpose.** The current ChampHLA / MVHLA benchmark relies entirely on 1000 Genomes truth
(Gourraud et al. 2014, IMGT/HLA 3.59.0). For the next analysis we want **independent, non-1000G**
cohorts that pair raw **WES / WGS / RNA-seq** reads with **gold-standard HLA typing** so we can test
generalisation of the Champion-Challenger consensus beyond 1000G.

**Status.** This is a shortlist + acquisition runbook only — **no data has been downloaded**.
Accessions marked _(verify)_ need confirmation against the live portal before use. Building Nextflow
input manifests is a separate follow-up task.

**Execution status (2026-06-25):**
- **NCI-60 RNA — DONE / VALIDATED (wave 1).** 4-line frozen-RNA-policy pilot complete (job 35255953):
  ChampHLA=MajorityVote=0.917; the frozen CC override gate fired twice, both corrective, 0 harmful
  (champion-only 0.75→0.917). See `analysis/nci60_benchmark/README.md`.
- **NCI-60 WES — RUNNING (wave 1b).** Same 5 lines, second modality; reads = open SRA **SRP150855**
  (NCI watson server was dead, but SRA hosts the exomes). Frozen WES CC policy. `conf/benchmark_nci60_wes.yaml`.
- **GIAB Ashkenazi-TRIO WES+RNA — DONE (2026-07-01).** HG002/HG003/HG004, non-1000G GERMLINE,
  clinical SBT gold truth (Stanford; F1000 8:1751 Table 1; Mendelian-consistent). **WES CC=MV=1.000
  (9/9 perfect); RNA CC=MV=0.889 (8/9).** ChampHLA's strong modalities, no LOH. Reads RNA
  SRR15909917-19 / WES SRR2962669/2692/2694. `analysis/giab_trio_benchmark/README.md`.
- **SweHLA / SweGen — AWAITING ACCESS.** User requested EGA (study **EGAS50000000906**, CRAM/GRCh38)
  + NBIS SweHLA truth; build ready (`analysis/swehla_benchmark/PHASE0_STATUS.md`).
- **829 WES — REJECTED** (cohort *is* 1000 Genomes — not independent).
- **652 RNA — REJECTED** for class-I bulk (Geuvadis/Montgomery = 1000G/HapMap; non-1000G subsets are
  class-II-only / mono-allelic / 10x scRNA). See `analysis/pcrsbt_benchmark/PHASE0_STATUS.md`.
- **IHWG B-LCL — PARKED** (truth public, but independent class-I reads sparse — HapMap wall).
- **HPRC — DEFERRED** (needs assembly→HLA truth derivation; HG002 clinical truth secured as a control).
- **consHLA / ZERO Childhood Cancer — REJECTED (integrity, verified 2026-06-28).** The open 76-patient
  Suppl. CSV is consHLA's **in-silico** consensus (HLA-HD-based), not experimental truth; the only
  clinical SSO truth (10 haem patients, Table 4) is **aggregate-only (no per-sample genotypes)** and
  those patients have **no deposited reads**. No (experimental truth × accessible reads) pairing exists.
  See `analysis/conshla_benchmark/PHASE0_STATUS.md`.
- **CCLE/DepMap — READ SOURCE ONLY** (open reads, but its HLA "types" are in-silico → pair with
  experimental truth, e.g. NCI-60 overlap). **BeatAML / Arab-trios — REJECTED** (no experimental truth).

**Fresh scout (2026-06-28) — best new lead is GeT-RM, but the reads wall holds:**
- **GeT-RM 108-panel — BEST EXPERIMENTAL TRUTH, but no matched short-read reads for the independent
  subset.** Bettinotti et al., *J Mol Diagn* 2018 (PMC**6939753**): 108 Coriell NIGMS DNA reference
  materials characterised for all 11 classical HLA loci by **multi-lab PCR-SSO + high-resolution SBT**,
  reported at **three-field** resolution (e.g. `A*02:01:01`) — higher than our 2-field Gourraud-2014
  truth, and openly published in the paper's Tables 1–3. **Independence:** spot-check of 12 listed
  `NA#####` IDs → 10 are **not** in our 1000G truth (only NA12878, NA12273 overlap), so a large
  non-1000G subset exists. **Reads BLOCKER (verified via ENA):** the non-1000G NIGMS lines lack public
  short-read WGS/WES/RNA — e.g. NA07439 has only Oxford-Nanopore **HLA-amplicon** runs (PRJEB61855,
  wrong modality), NA10005/NA16688 have none. The lines that *do* have deep reads (NA12878, NA12273…)
  are 1000G/HapMap. **CONFIRMED 2026-06-28 (truth built, reads-blocked):** the authoritative CDC
  consolidated table yields 108 NA-line 3-field genotypes (`bin/build_getrm_truth.py` →
  `analysis/getrm_benchmark/truth_long{,_3field}.tsv`), but **0/108 overlap our runs, 1/108 our 1000G
  truth**, and an ENA sweep finds only **2/108 with usable short reads** (NA12273 RNA, NA17221 WGS; 33
  amplicon-only). → excellent open experimental 3-field truth but **no matched short reads → not a viable
  benchmark cohort**; truth retained as a resource. See `analysis/getrm_benchmark/PHASE0_STATUS.md`.
- **A549 multi-modality (PRJNA667475, *HLA* 2025, Jiang et al. tan.70049) — single line, already covered.**
  Open short-read RNA-seq + long-read Iso-Seq for A549, which we already hold via NCI-60 (truth A549-ATCC,
  Adams SBT). Marginal value (one line; Iso-Seq could confirm its truth) — not a new cohort.
- **Confirmed dead ends (re-verified):** the 682-RNA-seq benchmark (PMC10827116, Mangul lab) class-I sets
  = Geuvadis/Montgomery 1000G/HapMap; the Feb-2026 medRxiv 29-biospecimen multi-platform WGS benchmark is
  **dbGaP controlled (phs004346)**.

**Second scout round (2026-06-28) — all hit the same four walls:**
- **ABraOM / SABE-WGS-1171** (Brazil, admixed non-1000G, 1171 WGS; Naslavsky 2022): HLA calls are
  **WGS-derived in-silico**, not experimental SBT; raw-read availability in SRA/ENA unconfirmed (portal
  serves variants). → in-silico truth → REJECT.
- **TCLP / Boegel cancer cell-line HLA catalog** + the 2025 pediatric cancer cell-line HLA compendium
  (iScience): RNA-seq **in-silico** (seq2HLA) types → REJECT (same as CCLE/DepMap).
- **precisionFDA Truth Challenge (HG002/3/4, matched Illumina+HiFi+ONT, MHC-focused)** = the **GIAB**
  Ashkenazi trio — i.e. the already-cleared GIAB option, not a new cohort.
- **PacBio reference-grade HLA, 46 Japanese DNAs** (PMC6180199): experimental reference-grade alleles,
  but reads are **HLA-amplicon** (ION PGM + SMRT), not WGS/WES/RNA → wrong modality; likely JPT overlap.
- **UK Biobank WES HLA-calling** (Nature Commun Biol 2023): huge, but HLA labels are tool-called
  (in-silico) and reads are application-access controlled → REJECT.

**Structural conclusion (after multi-session, multi-angle scouting).** Open *experimental* HLA truth +
open *short-read* WGS/WES/RNA + *non-1000G* is genuinely scarce: deeply-sequenced clinically-typed
germline samples are almost all 1000G/HapMap/Coriell; every fresh candidate falls into one of four buckets
— (a) in-silico "truth" (CCLE/TCLP/ABraOM/TCGA/UKB), (b) amplicon-only reads (donor registries, Japanese
reference, GeT-RM NA07439), (c) 1000G/HapMap overlap, or (d) controlled access (dbGaP phs004346, EGA
SweHLA/consHLA). The only clean independent open options remain **NCI-60** (DONE, RNA+WES) and **GIAB/HPRC**
(HG002 clinical SBT gold + open NIST WGS — cleared). The highest-ROI *new* contribution is GeT-RM's
**3-field orthogonal truth** re-scored on existing 1000G-overlap samples (a truth-quality upgrade
addressing the manuscript's 2014 two-field ceiling, not a new cohort). Controlled-access experimental-truth
cohorts (SweHLA WGS, the Feb-2026 dbGaP multi-platform set) require the user to obtain access first.

**Selection criteria.** (1) public or controlled-but-obtainable raw reads; (2) high-resolution HLA
truth (ideally Sanger SBT / NGS consensus, ≥2-field); (3) class I (A/B/C) at minimum, class II a
bonus; (4) **not** drawn from 1000 Genomes (to avoid truth/cohort circularity).

**CWD / CIWD catalogues — NOT a cohort; adopted as a commonness-stratification + plausibility layer
(2026-07-09).** Prompted by the Table 1 resource list in a 2024 *Best Pract Res Clin Haematol* HLA-typing
review (ScienceDirect S1521692624000252), we assessed the **CWD 2.0.0** (Mack 2013) and **CIWD 3.0.0**
(Hurley 2020, HLA 95:516; PMC7317522) catalogues. **Verdict: they are allele-classification catalogues,
not ground-truth datasets** — each row is an HLA allele tagged common / intermediate / well-documented /
not-CIWD across seven population groups, compiled from ~8M donor-registry typings. They carry **no
per-sample genotypes and no reads**, so they cannot pair with sequence data as a benchmark cohort (they
fall into the same donor-registry "amplicon-only reads" bucket already rejected above). → **Not added to
the cohort list.** Instead vendored as an *auxiliary analytical layer* over the existing ground-truth
benchmarks: `assets/ciwd_3.0.0.tsv` (+ `bin/ciwd.py`) drives (i) concordance **stratified by allele
commonness** (`summary_ciwd_stratified.tsv`, Fig. `fig_ciwd_stratified`) — showing whether the consensus
holds on rare/well-documented alleles, not just common ones — and (ii) a **biological-plausibility QC
flag** for not-CIWD/novel calls (`summary_ciwd_plausibility.tsv`). All additive; overall benchmark
numbers unchanged. See `assets/README_CIWD.md`.

---

## 1. SweHLA — high-confidence HLA from 1000 Swedish genomes  ·  **WGS**

- **What:** Population HLA bio-resource derived from the SweGen 1000 Swedish whole genomes. High-
  confidence consensus genotypes (class I **and** class II) determined by an *n−1 concordance rule*
  across multiple callers.
- **Why for us:** Large, independent **WGS** cohort with class I+II truth — directly stresses our
  weakest modality (WGS) on a non-1000G population.
- **Truth:** Multi-tool consensus, high-confidence subset published as a table; 2-field+.
- **Access:**
  - Paper: *SweHLA*, Eur. J. Hum. Genet. 2020 — https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7170882/
    (preprint: https://www.biorxiv.org/content/10.1101/660241).
  - HLA truth table: paper supplementary.
  - Raw WGS: **SweGen** — controlled access via NBIS/SciLifeLab (Swedish Bioinformatics);
    application required _(verify current data-access route)_.
- **Caveats:** WGS reads are controlled-access (data-access agreement). Population is European-only.

---

## 2. IHWG / 17th IHIW reference B-lymphoblastoid cell lines  ·  **cross-modality (RNA-seq + WGS/WES)**

- **What:** 382–406 International Histocompatibility Working Group reference B-LCLs typed by NGS in a
  17-lab single-blind study (17th International HLA & Immunogenetics Workshop). The de-facto
  **gold-standard reference panel** for evaluating HLA reagents/platforms; ≥88 unique haplotypes,
  60–95% of European allelic/haplotype diversity.
- **Why for us:** The cleanest possible HLA truth (multi-lab consensus, full-locus, high-resolution),
  and many of these immortalised lines have **public sequencing** (RNA-seq and/or WGS/WES) in
  ENA/SRA, enabling a true cross-modality test.
- **Truth:** 17th IHIW consensus NGS typing, class I + II, high resolution.
  - Paper: Hum. Immunol. 2019, "Next-generation HLA typing of 382 IHWG reference B-LCLs" —
    https://www.sciencedirect.com/science/article/abs/pii/S0198885919302277
  - Panel registry: https://bioregistry.io/registry/ihw ; IHWG @ Fred Hutch
    https://www.fredhutch.org/en/research/institutes-networks-ircs/international-histocompatibility-working-group.html
  - Ongoing: 18th IHIW reference cell-line project https://www.ihiw18.org/component-immunogenetics/reference-cell-lines/
- **Access:** Cell lines via Coriell / ECACC / IHWG repository. Sequencing reads: search ENA/SRA by
  IHW#### / cell-line ID _(verify per-line availability)_.
- **Caveats:** ⚠️ **Overlap check required** — some B-LCLs coincide with 1000G/HapMap lines (e.g.
  GM/NA IDs). Exclude any line already in our 1000G truth set to keep the test independent.

---

## 3. Human Pangenome Reference Consortium (HPRC)  ·  **WGS (+ long-read truth)**

- **What:** Release 1 = 47 diploid samples; **Release 2 (May 2025) = 200+** genetically diverse
  individuals with phased, >99% accurate long-read assemblies (PacBio HiFi + ONT) plus high-coverage
  Illumina short reads.
- **Why for us:** HLA truth can be derived directly from the **phased assemblies** (assembly-based
  HLA calling, e.g. via IPD-IMGT/HLA alignment / immuannot / locityper), giving an independent,
  modern, multi-ancestry **WGS** benchmark with matched Illumina reads to feed our short-read tools.
- **Truth:** Not pre-tabulated as HLA calls — **derive** from assemblies (extract MHC haplotypes →
  map to IMGT/HLA). Long-read T1K support in our pipeline is also directly testable here.
- **Access:**
  - Portal: https://humanpangenome.org/ ; Ensembl https://projects.ensembl.org/hprc/
  - Data: public AWS S3 (`s3://human-pangenome-reference`) / GitHub HPRC index _(verify bucket/path)_.
  - Draft pangenome paper: Nature 2023 — https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10172123/
- **Caveats:** Requires an assembly→HLA truth-extraction step (extra work, but fully public). HG002
  is shared with GIAB and well-characterised — good positive control.

---

## 4. Published PCR-SBT cohorts (cancer / cell-line panels)  ·  **WES + RNA-seq**

Smaller, paired-modality cohorts with orthogonal (PCR-SBT) truth — fastest path to a matched
**WES+RNA-seq** test.

- **NCI-60 cell-line panel (~58 lines):** high HLA-region coverage; used to validate arcasHLA &
  OptiType. RNA-seq (and some WES) public in GEO/SRA; HLA truth from prior typing.
  - arcasHLA validation context: https://pmc.ncbi.nlm.nih.gov/articles/PMC6956775/
- **Clinical cancer cohort (~28 patients), WES + RNA-seq, PCR-SBT gold standard** — surfaced in
  HLA-tool benchmarking literature _(identify exact paper + accession before use)_.
- **652 RNA-seq benchmark set** (gold-standard alleles): "A rigorous benchmarking of alignment-based
  HLA typing algorithms for RNA-seq data" — https://www.biorxiv.org/content/10.1101/2023.05.22.541750
- **829 WES cohort** (HLA-A,-B,-C,-DRB1,-DQB1): "Benchmarking freely available HLA typing
  algorithms…" — https://pubmed.ncbi.nlm.nih.gov/36426357/ ,
  https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9679531/
- **MHC genotyping benchmark (class I + II):** https://pmc.ncbi.nlm.nih.gov/articles/PMC10170851/
- **Caveats:** Truth resolution/method varies per study; some cohorts are controlled-access
  (dbGaP/EGA). Confirm truth method (PCR-SBT vs in-silico) per dataset — in-silico "truth" is not
  acceptable as ground truth for us.

---

## Recommended acquisition order

1. **IHWG 17th IHIW B-LCLs** — best truth quality; start by pulling the published consensus typing
   table and screening ENA/SRA for lines with public RNA-seq/WGS, **after removing 1000G/HapMap
   overlaps**. Fastest route to a clean cross-modality test.
2. **NCI-60 / small PCR-SBT cohorts** — quickest matched **WES+RNA-seq** validation; modest size,
   mostly open-access.
3. **SweHLA** — strongest independent **WGS** test (our weak modality); begin the SweGen
   data-access application early (controlled access has lead time).
4. **HPRC** — highest long-term value (large, modern, multi-ancestry) but needs an assembly→HLA
   truth-derivation pipeline; schedule once 1–3 are moving.

## Integrity notes (per project policy)
- Use only **experimentally determined** truth (Sanger SBT / NGS multi-lab consensus); do **not**
  treat in-silico calls as ground truth.
- Record IMGT/HLA version for each truth source and normalise to the pipeline-pinned version before
  scoring.
- Flag and exclude any sample overlapping the existing 1000G truth set to preserve independence.

_Sources gathered via web search, June 2026; accessions marked (verify) pending confirmation._
