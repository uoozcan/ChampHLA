# ChampHLA Manuscript Draft V3

## Title
**ChampHLA: a reproducible Nextflow ensemble for champion–challenger HLA typing across WGS, WES, and RNA sequencing**

**Running title:** ChampHLA ensemble HLA typing

**Authors:** Umut Onur Özcan¹, Nemo Ikonen¹, Philipp Sergeev¹, Ella Sinervuori¹, Khalid Saeed², Ziaurrehman Tanoli¹˒³*, Markus Vähä-Koskela¹˒³*

**Affiliations:**
¹Institute for Molecular Medicine Finland (FIMM), HiLIFE, University of Helsinki, Helsinki, Finland
²Joint AstraZeneca–Cancer Research Horizons Functional Genomics Centre, Cambridge, UK
³iCAN Digital Precision Cancer Medicine Flagship, University of Helsinki and Helsinki University Hospital, Helsinki, Finland

*Corresponding authors

**Availability:** https://github.com/uoozcan/ChampHLA | IMGT/HLA 3.59.0 | Pipeline version 2.0.0

---

## Abstract

Human leukocyte antigen (HLA) genes encode the polymorphic cell-surface molecules that present peptide antigens to T lymphocytes and are therefore central to adaptive immunity, allograft survival, autoimmune disease susceptibility, pharmacogenomic safety, and cancer neoantigen presentation. Accurate HLA typing at the population scale—increasingly at high, allele-level resolution—is now a prerequisite for transplantation matching, immunotherapy design, and pharmacogenomics screening, yet the field currently relies on a fragmented ecosystem of independently developed computational typing tools, each with distinct input format requirements, resolution capabilities, and accuracy profiles across sequencing modalities. No tool performs uniformly well across whole-genome sequencing (WGS), whole-exome sequencing (WES), and RNA-seq, and practitioners who run multiple tools routinely encounter discordant allele calls for the same sample, with no principled mechanism for resolving the disagreement.

Here we present **ChampHLA** (formerly PIHLA/mvHLA), a reproducible Nextflow DSL2 pipeline that integrates eight state-of-the-art HLA typing tools—OptiType, HLA-HD, ArcasHLA, T1K, SpecHLA, Kourami, POLYSOLVER, and Seq2HLA—into a single unified framework with containerised execution, harmonised output schemas, benchmark-derived reliability weights, and a novel **Champion-Challenger** ensemble consensus method. ChampHLA eliminates per-tool installation burden, routes each tool to the appropriate input format automatically, and resolves cross-tool disagreement through a principled override policy: each HLA gene is assigned a benchmark-designated champion, and the ensemble may override the champion's call only when the aggregate evidence meets simultaneous criteria for support breadth, weight margin, tool diversity, and allele unambiguity. The champions and override thresholds are selected under nested cross-validation—learned on training folds and evaluated strictly out of fold—so that no reported accuracy is inflated by tuning on the data it is measured against, and every called locus carries a fully auditable decision trace.

In a benchmarking study spanning three sequencing modalities across cohorts derived from the 1000 Genomes Project, evaluated under **10-fold population-stratified nested cross-validation** with operating points selected strictly out of fold, ChampHLA's benefit scaled with inter-tool discordance. In the difficult WGS modality it achieved an overall correct-call rate of 0.5012 (95% CI: 0.453–0.549; n=137), **significantly exceeding majority voting** (0.3844; +11.7 percentage points; exact McNemar p<0.0001) and matching **the best single tool** (OptiType, 0.4975; McNemar p=0.74, not significant), by routing each locus to its benchmark-designated champion—recovering the single-tool accuracy ceiling that symmetric voting dilutes. In the high-accuracy WES and RNA-seq modalities, where strong tools already largely agree, it **matched majority voting** within overlapping confidence intervals (WES 0.9436 vs 0.9359; RNA-seq 0.9408 vs 0.9502; McNemar p=0.45, not significant). An ablation confirmed that this gain arises from champion routing rather than the weighting scheme: pure reliability-weighted voting underperformed majority voting in every modality. In a supplementary matched-subject cross-modal analysis of 106 subjects, bimodal WES+RNA majority voting achieved the best multimodal result (0.9623 overall correct-call rate)—which Champion-Challenger matches exactly while exceeding the best single tool—so Champion-Challenger can serve as the default consensus in the recommended bimodal configuration at no accuracy cost; adding WGS to form a trimodal consensus did not improve performance (0.9591), identifying short-read WGS single-tool quality as the current accuracy bottleneck. To extend benchmarking to cohorts that lack gold-standard typing, ChampHLA also generates an orthogonal, agreement-gated **silver-standard truth**—Locityper genotyping against a Human Pangenome Reference Consortium (HPRC) pangenome database, annotated by Immuannot—which reproduced 1000 Genomes gold calls at 81.1% (two-field) / 82.2% (G-group) allele-level concordance across 30 samples—rising to 94.4% on a smaller per-sample-depth pilot—directly addressing the two-field resolution ceiling of the 2014 truth panel. Finally, to test generalisation beyond the cohort on which the policy was calibrated, we applied the frozen 1000G-derived weights and Champion-Challenger policy—without any re-learning—to the independent, non-1000-Genomes NCI-60 cancer cell-line panel scored against Adams et al. (2005) sequence-based typing: ChampHLA matched majority voting on whole-exome data (0.842) and on RNA-seq (0.793 versus 0.817, overlapping confidence intervals) with an active, predominantly corrective override gate (15 of 20 RNA overrides corrective), indicating that the Champion-Challenger framework transfers out of distribution.

ChampHLA is freely available at https://github.com/uoozcan/ChampHLA and requires only Nextflow (≥23.04.0) and Docker or Singularity; no per-tool installation is needed. Pre-configured profiles are provided for Docker-based local execution and SLURM-based HPC deployment, including a dedicated configuration for the CSC Puhti supercomputer.

---

## Abstract (Short, ≤200 words)

Human leukocyte antigen (HLA) typing is essential for transplantation, pharmacogenomic safety, autoimmune genetics, and cancer immunotherapy, yet NGS-based typing remains fragmented across tools with incompatible inputs and no principled way to reconcile discordant predictions. We present **ChampHLA**, a reproducible Nextflow DSL2 pipeline integrating eight HLA typing tools (OptiType, HLA-HD, ArcasHLA, T1K, SpecHLA, Kourami, POLYSOLVER, Seq2HLA) in a containerised framework with a novel **Champion-Challenger** ensemble: a per-gene benchmark champion whose call is overridden only when weighted ensemble evidence is simultaneously broad, marginal, tool-diverse, and unambiguous, with a fully auditable decision trace. Under 10-fold nested cross-validation with out-of-fold operating points, ChampHLA's advantage scaled with inter-tool discordance—on WGS it significantly beat majority voting (0.5012 vs 0.3844; exact McNemar p<0.0001) and matched the best single tool via champion routing, while on the high-accuracy WES and RNA-seq modalities it matched majority voting (0.9436, 0.9408; not significant); an ablation showed the gain comes from routing, not weighting. ChampHLA also generates agreement-gated silver-standard truth (Locityper + HPRC pangenome + Immuannot) reproducing 1000 Genomes gold calls at 81.1% across 30 samples for cohorts lacking gold typing, and under frozen weights transferred to a non-1000G panel (NCI-60), matching majority voting. ChampHLA is freely available at https://github.com/uoozcan/ChampHLA.

---

## Introduction

### 1. HLA Typing: Biological Role and Clinical Importance

The major histocompatibility complex (MHC), encoded in humans on chromosome 6p21.3 and commonly termed the human leukocyte antigen (HLA) system, constitutes the most polymorphic locus in the entire human genome. The classical HLA class I genes—HLA-A, HLA-B, and HLA-C—encode transmembrane glycoproteins that present intracellular peptide fragments to CD8⁺ cytotoxic T lymphocytes, whereas the classical HLA class II genes—HLA-DR, HLA-DQ, and HLA-DP—present extracellular peptides to CD4⁺ helper T lymphocytes. Together, these molecules serve as the primary mechanism through which the adaptive immune system surveys cellular proteomes for evidence of infection, malignant transformation, or graft non-self. At the latest count, more than 40,000 distinct alleles have been catalogued across all classical HLA loci in the IMGT/HLA database (Robinson et al., 2020), a count that continues to grow with the expansion of large-scale sequencing programs.

The clinical stakes of accurate HLA typing are considerable. In haematopoietic stem cell transplantation (HSCT), recipients and donors are matched for at least HLA-A, -B, -C, and -DRB1 (8/8 match) or preferably also for HLA-DQB1 (10/10 match); each additional mismatch at these loci increases the risk of acute graft-versus-host disease and graft failure, while paradoxically also influencing graft-versus-tumour effect strength (Petersdorf et al., 2007). Beyond transplantation, HLA alleles underlie some of the strongest pharmacogenomic safety signals known: HLA-B*57:01 is a near-perfect predictor of abacavir hypersensitivity reaction (Mallal et al., 2008), and HLA-B*15:02 predicts severe cutaneous adverse reactions to carbamazepine in Southeast Asian populations (Chung et al., 2004). In cancer immunotherapy, HLA genotype determines which tumour-specific neoantigens can be processed and presented to cytotoxic lymphocytes; errors in HLA typing therefore propagate directly into personalised vaccine and T-cell receptor therapy design pipelines, potentially selecting ineffective or immunologically invisible neoantigen candidates (Schumacher and Schreiber, 2015). In genetic association studies, the HLA locus harbours hundreds of the strongest known disease associations for autoimmune and inflammatory conditions, including type 1 diabetes, rheumatoid arthritis, ankylosing spondylitis, and coeliac disease (Matzaraki et al., 2017). These applications collectively demand HLA typing that is simultaneously accurate, high-resolution, reproducible, and scalable to the cohort sizes now routinely encountered in clinical genomics and biobank-linked studies.

### 2. Computational Challenges of HLA Locus Typing

Although next-generation sequencing (NGS) has largely displaced serological and sequence-specific oligonucleotide typing for population-scale HLA studies, the HLA locus poses unique computational challenges that distinguish it from the rest of the human genome. The extreme polymorphism of the locus means that short reads generated from standard whole-genome or whole-exome protocols must be aligned to an enormously diverse reference space: the IMGT/HLA database describes allele-level variation for HLA-A alone across sequences that differ from one another by tens to hundreds of base pairs at polymorphic sites. Standard whole-genome alignment strategies, which map reads to a single linear reference genome, systematically mismap or lose reads that originate from alleles poorly represented in the reference, introducing coverage and allele-frequency biases that distort downstream typing calls.

The structural organisation of the HLA region compounds this difficulty. The six-megabase MHC on chromosome 6p21.3 contains not only the classical HLA genes but also pseudogenes, class I chain-related genes (MICA, MICB), complement component genes, and numerous non-classical HLA-like loci. Short reads originating from polymorphic sites in classical loci frequently align with comparable or higher confidence to paralogous non-classical sequences, producing spurious allele calls that are difficult to distinguish from true mappings without allele-specific databases and scoring models. At the RNA level, allele-specific expression adds a further complication: heterozygous individuals may express one allele at a substantially lower level than the other, causing RNA-seq-based callers to under-represent or entirely miss the minor allele in low-expression conditions. Collectively, these challenges mean that no single alignment strategy, reference space, or scoring model is optimal across all HLA genes, all sequencing modalities, and all population backgrounds.

### 3. State-of-the-Art HLA Typing Tools

The past decade has produced a substantial ecosystem of computational HLA typing tools, each targeting a distinct subset of the typing problem. **OptiType** applies integer linear programming over reads aligned to an HLA allele reference, maximising the number of reads explained jointly by an allele pair; it was originally designed for class I typing from WES and WGS data and delivers high precision when sufficient HLA-region coverage is available (Szolek et al., 2014). **ArcasHLA** reconstructs allele pairs from RNA-seq-aligned BAM files using a diploid reconstruction algorithm that leverages the sensitivity of transcriptomic coverage at expressed alleles; it performs strongly on RNA-seq but its accuracy degrades substantially on DNA inputs where allele expression is not informative (Orenbuch et al., 2019). **HLA-HD** performs bidirectional alignment of FASTQ reads to a bowtie2-indexed reference derived from IMGT/HLA and supports multi-locus high-resolution class I and class II typing; it is well suited to WGS and WES inputs and achieves competitive accuracy at four-field resolution (Mishima et al., 2019). **T1K** uses a k-mer-based typing strategy that scales to large cohorts, supports both short-read and long-read inputs (PacBio HiFi and Oxford Nanopore), and extends coverage to classical class II loci (Song et al., 2023). **SpecHLA** performs HLA typing through variant-aware alignment with an exon-only mode specifically designed for WES data; it provides higher resolution in exon-targeted regions but requires native installation because no public Docker image is maintained. **Kourami** assembles HLA allele sequences from graph-guided read alignment against a curated allele-graph database and is optimised for WGS BAM inputs (Lee and Kingsford, 2018). **POLYSOLVER** uses a Bayesian statistical scoring model over WES read alignments against an allele reference and has been applied extensively in cancer genomics pipelines (Shukla et al., 2015). **Seq2HLA** types HLA alleles from RNA-seq reads by aligning against an HLA allele index and applying a statistical assignment model; it is rapid but restricted to RNA-seq inputs and class I loci (Boegel et al., 2012). Each of these tools has demonstrated utility in its target context, but none performs uniformly well across all sequencing modalities, population backgrounds, and resolution requirements.

### 4. Gap and Motivation

Despite the strength of individual HLA typing tools, the field lacks a unified, reproducible pipeline that integrates multiple tools into a single harmonised analysis and resolves the disagreements between them in a principled, auditable manner. Current practice requires practitioners to install each tool independently—a process often complicated by conflicting dependencies, tool-specific reference database paths, and environment-specific compilation requirements—then run each tool separately, parse heterogeneous output formats, and manually adjudicate disagreements with no systematic framework for deciding which tool to trust on any given gene or sample. Input format requirements are particularly heterogeneous: some tools accept only FASTQ reads, others require coordinate-sorted BAM alignments, and still others have been validated only against specific reference genome builds or HLA database versions, meaning that applying multiple tools to the same sample frequently requires custom preprocessing and format conversion scripts. Resolution also varies: some tools report alleles at two-field resolution, others at four-field or G-group notation, and the absence of a common normalisation layer makes cross-tool accuracy comparisons unreliable. When tools disagree—which they do non-trivially, particularly at HLA-A and HLA-C in WGS data, and at heterozygous RNA-seq loci affected by allele-specific expression—the practitioner has no basis for choosing between conflicting predictions without a reliability reference grounded in empirical benchmarking. ChampHLA was developed to address these gaps: it provides a single Nextflow DSL2 pipeline that encapsulates all eight tools in containerised modules, routes inputs to each tool's required format automatically, normalises outputs to a common allele schema, and applies a benchmark-calibrated Champion-Challenger ensemble method to resolve disagreements and generate reliable consensus calls without requiring any per-tool installation or manual output curation. A further gap concerns evaluation rather than typing: outside a small number of public panels such as the 1000 Genomes Project, gold-standard HLA truth is scarce, which prevents practitioners from benchmarking ensemble accuracy on their own cohorts—a gap ChampHLA addresses by generating an orthogonal, agreement-gated silver-standard truth from pangenome- and assembly-based callers.

### 5. Study Overview

Here we describe ChampHLA version 2.0.0, a reproducible Nextflow ensemble pipeline for HLA typing that integrates eight tools across WGS, WES, RNA-seq, and long-read sequencing. ChampHLA makes three distinct methodological contributions beyond tool integration. First, an empirical calibration guardrail system evaluates each tool's confidence signal against a benchmark truth set and blocks miscalibrated confidence from inflating consensus weights; only tools whose confidence proxies correlate reliably with call correctness (Brier score and expected calibration error both below 0.35) receive a confidence-based weight boost, while miscalibrated tools fall back to base-reliability weighting. Second, the Champion-Challenger ensemble method designates a per-gene benchmark champion and applies a parameterised, nested-cross-validated override policy (thresholds selected on training folds and applied out of fold) that only replaces the champion's call when the weighted ensemble evidence for an alternative is simultaneously sufficient in breadth, magnitude, tool diversity, and allele unambiguity. Third, to support benchmarking on cohorts that lack gold-standard typing, ChampHLA generates an orthogonal, agreement-gated silver-standard truth from two pangenome- and assembly-based callers (Locityper and Immuannot) that are held separate from the eight-tool consensus panel. The result is a framework in which practitioners can obtain reliable, auditable HLA consensus calls from a single command across any supported sequencing modality, using only Nextflow and a container runtime, with every decision traceable from raw input through harmonised calls to final consensus output.

---

## Materials and Methods

### 1. ChampHLA Nextflow Pipeline Architecture

ChampHLA is implemented in Nextflow DSL2 (Di Tommaso et al., 2017) as a modular pipeline in which each of the eight HLA typing tools occupies a self-contained process module (**Figure 1**). Each module encapsulates the complete execution environment for one tool—including the container image, input preprocessing, tool invocation, and output normalisation—so that adding or updating a tool requires only a change to its isolated module without modifying the core workflow logic. Containerisation is supported via Docker (for local workstations) and Singularity (for HPC environments where Docker is not available due to privilege requirements), with pre-built container images pulled automatically at runtime or from a user-specified Singularity SIF cache directory. The pipeline is version-controlled alongside all configuration files, and pipeline runs are anchored to a versioned manifest of input samples, ensuring that results are independently reproducible from the same inputs.

The pipeline accepts three input formats: paired-end FASTQ files, sorted and indexed BAM files, and CRAM files with a supplied reference FASTA. Input format discovery proceeds via either a directory scan (all files matching the specified format in `--input`) or an explicit samplesheet CSV that pairs sample identifiers with file paths and optional metadata. CRAM files are first decoded to BAM using the supplied reference; BAM files are then routed to tools based on their requirements. For WES BAM inputs, the pipeline optionally extracts only reads mapping to the HLA region (chromosomal coordinates chr6:28,477,797–34,448,354) before format conversion, substantially reducing intermediate disk space and processing time without loss of typing sensitivity. For tools requiring FASTQ input (OptiType, T1K, HLA-HD, Seq2HLA), BAM files are converted via `samtools fastq` with duplicate-read removal. For tools requiring coordinate-sorted BAM (ArcasHLA, SpecHLA, POLYSOLVER, Kourami), BAM files are passed directly or re-sorted as needed.

The pipeline supports four sequencing type modes: `dna` (WGS or WES), `rna` (bulk RNA-seq), `longreads_hifi` (PacBio HiFi), and `longreads_ont` (Oxford Nanopore Technology). Sequencing type is declared at runtime via `--seq_type` and determines which subset of tools is activated: POLYSOLVER and Kourami are skipped for RNA-seq runs because they accept only DNA BAM inputs; Seq2HLA and ArcasHLA receive the `--seq_type rna` flag activating their RNA-specific alignment modes; long-read runs activate only T1K, which is the sole tool with validated long-read support in the current version. The active tool set can be further restricted by the user via the `--tools` parameter, which accepts a comma-separated list of tool names (e.g., `--tools optitype,arcashla,hlahd`).

**HPC Execution.** For SLURM-based clusters, ChampHLA provides a generic `slurm` profile and a pre-tuned `puhti` profile for the CSC Puhti supercomputer. The SLURM profile submits each tool module as a separate job using `executor = 'slurm'`, with resource labels (`process_low`, `process_medium`, `process_high`) that map to 2–16 CPUs, 8–64 GB RAM, and 2–24 hour time limits. Retry logic with linear resource scaling is built in: if a job fails due to memory overflow or wall-time limit, Nextflow automatically resubmits with multiplied resources up to the pipeline maximum (controlled by `--max_cpus`, `--max_memory`, `--max_time`). Singularity image resolution is handled via a user-specified cache directory (`--singularity_cache_dir`), avoiding repeated download of large container images across runs. A typical HPC invocation is:

```bash
nextflow run champhla/main.nf \
  --input_samplesheet samples.csv \
  --input_type bam \
  --seq_type dna \
  --tools optitype,hlahd,t1k,arcashla,spechla,kourami,polysolver \
  --enable_majority_voting \
  --weighting calibrated \
  --mv_genes A,B,C \
  --max_cpus 16 \
  --max_memory 64.GB \
  -profile puhti,singularity \
  -resume
```

The `-resume` flag enables Nextflow's content-hash-based caching, so that completed tool steps are not rerun when adding new samples or updating a single module in a subsequent run.

**Local Execution.** For local workstations with Docker installed, the `docker` profile is used. Minimum recommended hardware is 8 CPUs and 32 GB RAM; running the full eight-tool panel on a single WES sample requires approximately 32–64 GB peak RAM and 4–12 hours wall time depending on which tools are enabled. Memory-constrained environments can use `--extract_hla_region` to reduce the BAM input to the HLA region before tool execution. A typical local invocation for an RNA-seq experiment is:

```bash
nextflow run champhla/main.nf \
  --input fastq_files/ \
  --input_type fastq \
  --seq_type rna \
  --tools optitype,arcashla,hlahd,t1k,seq2hla \
  --optitype_seq_type rna \
  --enable_majority_voting \
  -profile docker
```

**Output Schema.** Each run emits a standardised output directory containing: `harmonized_calls.tsv` (one row per tool × sample × gene, with allele pair, call status, and correctness flags when truth is available); `tool_availability_summary.tsv` (callable status per tool × sample); `consensus_calls.tsv` (per-sample consensus allele pair per gene); `runtime_weights.json` (the modality-specific tool weights applied at runtime); and `benchmark_metadata.json` (pipeline version, IMGT/HLA version, parameters, and cohort metadata). When truth data are provided via a benchmark YAML configuration, additional artifacts are emitted: `method_comparison.tsv`, `confidence_calibration_summary.tsv`, `tool_confidence_weights.tsv`, `discordance_summary.tsv`, and `champion_challenger_method_comparison.tsv`. Per-tool computational resource profiles (peak RAM, CPU utilisation, wall-clock time by modality) are shown in **Supplementary Figure S14**.

### 2. Champion-Challenger Ensemble Method

The central methodological contribution of ChampHLA is the **Champion-Challenger** (CC) consensus algorithm, which frames ensemble HLA calling as a structured decision problem rather than a symmetric vote. The algorithm operates per gene, per sample, after all tool outputs have been collected and harmonised. **Figure 6** contrasts the routed, evidence-gated Champion-Challenger decision flow with the symmetric majority-voting baseline and summarises its effect across the three modality benchmarks under nested cross-validation.

For each HLA gene (HLA-A, -B, -C; extensible to other loci), a single tool is designated as the **champion** based on benchmarking evidence—specifically, the tool with the highest overall correct-call rate for that gene and modality on the training folds under nested cross-validation. The champion's call is the default output. The ensemble may produce a **challenger override** only if all four of the following conditions are satisfied simultaneously:

1. **Support fraction:** the aggregate benchmark-derived weight supporting the challenger allele pair constitutes at least a minimum fraction of the total weight across all candidate allele pairs (`min_challenger_support_fraction`);
2. **Weight margin:** the challenger's aggregate weight exceeds that of the next-best allele pair candidate by at least a minimum margin (`min_challenger_margin`);
3. **Tool diversity:** at least a minimum number of distinct callable tools independently support the challenger allele pair (`min_supporting_tools`);
4. **Ambiguity guard:** if enabled, the challenger allele pair must be non-ambiguous under the current IMGT/HLA nomenclature, blocking overrides driven by nomenclature-ambiguous calls.

When the champion tool is not callable for a given sample–gene pair (e.g., due to a tool failure or low read coverage abstention), the algorithm falls back to weighted consensus rather than leaving the locus uncalled, preserving callable rate. Every call is recorded with a full decision trace indicating whether the outcome was `champion_retained`, `champion_missing` (fallback used), or `challenger_override`. All triggered overrides are post-hoc classified as **corrective** (override correct, champion wrong), **harmful** (override wrong, champion correct), or **neutral** (both calls wrong) using available truth labels, enabling transparent auditing of ensemble decision quality across parameter settings.

The four threshold parameters are not set by hand. To avoid the optimistic bias of selecting an operating point on the same data it is evaluated against, they are chosen under **nested cross-validation**. Within each of 10 population-stratified folds, a grid of values—`min_challenger_support_fraction` ∈ {0.20, 0.35, 0.50, 0.65}; `min_challenger_margin` ∈ {0.0, 0.05, 0.10, 0.20}; `min_supporting_tools` ∈ {1, 2, 3}—is evaluated on the training folds only; the combination maximising training-fold correct-call rate, together with the per-gene champions (likewise designated from the training folds as the highest-accuracy single tool per gene), is frozen, and the held-out fold is scored with that frozen policy. Held-out predictions are pooled across folds for reporting, and the paired Champion-Challenger-versus-majority-voting comparison is assessed with an exact McNemar test on discordant loci (`bin/nested_cv_champion_challenger.py`). For end users, the pipeline additionally writes the full-cohort sweep surface to machine-readable TSV files (`wes_champion_override_sweep.tsv`, `rna_sweep_summary.json`); this surface is an in-sample tuning aid for selecting an institution-specific operating point and is not itself a performance estimate. The accuracy figures reported in this study are the held-out nested-cross-validation estimates produced by this re-analysis (`bin/nested_cv_champion_challenger.py`, applied to the harmonised call tables). The released Nextflow pipeline applies a single frozen operating point for production—the bundled per-modality weights (`tool_weights_{wgs,wes,rna}.json`) and default override policy, which are the full-cohort analogue of the per-fold policy—and is not intended to re-estimate held-out accuracy; users benchmarking on their own truth-backed cohorts should use the nested-CV re-analysis to obtain unbiased estimates.

### 3. Benchmark-Derived Confidence Weights and Calibration Guardrail

ChampHLA learns per-tool reliability weights from benchmark data rather than treating all tools as equally trustworthy. For each tool, modality, and gene, a **base reliability** weight is computed as the fraction of calls matching the ground-truth allele pair at two-field resolution on the training folds of the nested cross-validation, so that weights applied to a held-out fold are never learned from it. Tools with high base reliability receive proportionally higher influence on the weighted consensus vote.

In addition to base reliability, each tool's **confidence signal** is evaluated empirically. Tool-native confidence proxies are extracted and normalised to [0, 1]: OptiType's integer-linear-programming objective score is divided by the maximum observed training-split value; ArcasHLA's posterior allele-pair probability is used directly; HLA-HD, T1K, and Kourami use per-allele read-support counts divided by a 50-read normalisation target; SpecHLA, POLYSOLVER, and Seq2HLA produce no per-call confidence output and are unconditionally assigned `no_confidence` status. Calibrated probabilities are derived from proxy scores using Platt scaling (logistic regression on training-split proxy scores and correctness labels with leave-one-out cross-validation). Calibration quality is assessed by Brier score and expected calibration error (ECE) computed over 10 equal-frequency confidence bins.

A **calibration guardrail** blocks miscalibrated confidence from distorting ensemble weights: only tools with both Brier score and ECE below the guardrail threshold (0.35 in the current study) are permitted to use confidence for weight boosting. For tools that pass the guardrail, **effective_confidence** is the Platt-scaled calibrated probability output from the leave-one-out logistic regression applied to the tool's native confidence proxy on training-split calls. For blocked tools, effective_confidence is set equal to base_reliability, collapsing the blended weight to base_reliability entirely. Note that the guardrail threshold (0.35) and the blend coefficients (0.7/0.3) are independent parameters: the threshold determines whether confidence enters the blending formula at all, while the coefficients determine how heavily confidence is weighted versus base_reliability when it does. The final blended weight is:

```
final_weight = 0.7 × base_reliability + 0.3 × effective_confidence
```

where the 0.7/0.3 coefficients represent a conservative benchmark-derived default that keeps reliability the dominant signal. This formula is applied at runtime when `--weighting calibrated` is specified; uniform equal weights are applied under `--weighting equal`. Pre-computed modality-specific weight files (`tool_weights_wgs.json`, `tool_weights_wes.json`, `tool_weights_rna.json`) are bundled with the pipeline and applied by default when no in-run benchmark data are available.

### 4. Tool Harmonisation and Allele Normalisation

Outputs from all tools are parsed by tool-specific parsers into a unified allele call schema. Each record encodes: `sample`, `tool`, `modality`, `gene`, `allele1`, `allele2` (the two-field normalised pair), `raw_allele1`, `raw_allele2` (the original caller output for audit), `call_status` (called / no_call / low_confidence), `is_callable` (0/1), `is_correct` (0/1 when truth is available), and `imgt_hla_version`. Allele strings are normalised by stripping colon-delimited fields beyond the second (e.g., `HLA-A*02:01:01` → `HLA-A*02:01`), applying IMGT/HLA nomenclature rules for null-allele and expression-variant suffix notations, and resolving ambiguity-code expansions. The IMGT/HLA database version is pinned in the benchmark configuration (currently 3.59.0) and propagated through `benchmark_metadata.json`, ensuring that truth normalisation and call normalisation operate under identical nomenclature assumptions.

### 5. Discordance Taxonomy

Post-consensus disagreement is characterised using a five-category taxonomy applied as an interpretation layer over harmonised calls; category assignments do not alter accuracy statistics. `technical_conflict` is assigned when at least two callable tools return divergent allele pairs without a clear supermajority. `low_evidence_conflict` is applied when total benchmark-derived weight across callable tools falls below the minimum support threshold. `dna_rna_discordance` labels disagreements between DNA-modality (WGS/WES) and RNA-modality calls for the same sample and gene in a multimodal context, consistent with somatic loss of heterozygosity or allele-specific expression. `possible_expression_bias` labels intra-RNA disagreements consistent with allele-specific expression causing one allele to be under-represented in the RNA-seq read pool. `no_evidence` is applied when no callable output is available from any tool for a given sample–gene pair. Discordance labels are written to `discordance_summary.tsv` and `discordance_tags.tsv` without modifying consensus output.

### 6. 1000 Genomes Project Benchmark Cohort

The primary WGS benchmark cohort comprised 137 truth-backed 1000 Genomes Project (1000G) samples spanning multiple superpopulation groups, with six WGS-capable tools (OptiType, ArcasHLA, HLA-HD, T1K, SpecHLA, Kourami) contributing benchmark-ingestable outputs. Rather than a single fixed holdout, WGS accuracy was evaluated by 10-fold population-stratified nested cross-validation: within each fold the per-gene champions, confidence weights, and Champion-Challenger override policy were learned on the training folds and applied to the held-out fold, with predictions pooled across folds for reporting (Methods §2). A 50-sample "wave2" subset with a frozen 28/10/12 train/validation/holdout split is additionally retained for the confidence-calibration-guardrail characterisation (Results §4, Table 5), which is a per-tool property estimated on a dedicated training split.

The secondary WES benchmark comprised 130 truth-backed 1000G samples from populations including CEU (n=12), FIN (n=7), GBR (n=6), TSI (n=6), YRI (n=7), and 92 additional samples. Seven tools contributed to the WES benchmark: ArcasHLA, HLA-HD, Kourami, OptiType, POLYSOLVER, SpecHLA, and T1K; Kourami produced callable outputs for 120 of 130 samples. The secondary RNA-seq benchmark comprised 107 truth-backed 1000G samples. Six tools ran on RNA: ArcasHLA, HLA-HD, OptiType, Seq2HLA, SpecHLA, and T1K. Seq2HLA was callable on 29 of 107 samples due to a container configuration failure; SpecHLA achieved a callable rate of 76.2% (80/105 samples with HLA-region reads) and was excluded from default production accuracy scoring due to its near-zero RNA-seq accuracy (0.035); the RNA-seq Champion-Challenger accuracy reported in this study uses a SpecHLA-inclusive scoring run to benchmark the ensemble against the largest possible tool panel. A supplementary matched-subject trimodal robustness cohort of 106 samples was defined as the intersection of samples with WGS, WES, and RNA-seq outputs and ground-truth HLA calls.

HLA ground-truth labels were derived from the 2014 1000 Genomes HLA typing study (Gourraud et al., 2014), normalised to IMGT/HLA 3.59.0 at two-field resolution for HLA-A, -B, and -C. The IMGT/HLA version was pinned throughout all benchmark runs and propagated through `benchmark_metadata.json`. Primary accuracy was defined as exact two-field allele-pair correct-call rate: a locus was counted correct only when both alleles of the returned pair matched the ground-truth pair exactly at two-field resolution. Wilson score 95% confidence intervals were computed for all proportional accuracy estimates. Because every method is scored on the same gene-locus pairs, the pre-specified primary comparison—Champion-Challenger versus majority voting—was assessed per modality (three tests) with an exact McNemar test on discordant loci. The per-gene WGS McNemar tests and all secondary stratifications (allele commonness, ancestry, resolution) are reported as descriptive robustness checks rather than confirmatory tests; we note that the two significant per-gene WGS results (HLA-A and HLA-B) remain significant after Holm–Bonferroni correction across the three per-gene tests (both adjusted p<0.005). All primary and stratified estimates derive from a single, frozen analysis provenance: the harmonised nested-cross-validation call tables under `analysis/nested_cv_champion_challenger/` and `analysis/benchmark_{wgs,wes,rna}_cv_recalibrated/`, regenerable from the deposited inputs by `bin/nested_cv_champion_challenger.py` (Data Availability), avoiding the "garden of forking paths" that an unpinned, in-sample operating-point search would introduce.

### 7. Computational Performance Benchmarking

Per-tool wall-clock time, peak resident memory (RAM), and CPU utilisation were measured from Nextflow execution traces (`execution_trace.txt`) produced automatically by the ChampHLA pipeline during 1000 Genomes Project benchmark runs. The trace captures three fields per completed process: `realtime` (wall-clock duration), `peak_rss` (peak resident set size in MB), and `%cpu` (average CPU utilisation). Timing was collected across n=30 independent samples per modality. WGS timing data were collected from a dedicated n=30 performance benchmark run using population-stratified samples (6 per superpopulation: CEU, CHB, GBR, TSI, YRI) from the WGS cohort. RNA-seq timing data were collected from a dedicated n=30 performance benchmark run using HLA-region-extracted BAMs from the truth-backed RNA-seq cohort; tool-level sample counts vary (T1K and Seq2HLA n=30; ArcasHLA, HLA-HD, and SpecHLA n=29; OptiType n=26 due to four samples exceeding the 8 h wall-clock limit). WES timing was measured from 30 samples selected from the truth-backed WES cohort (n=130 total), choosing those with the most tools completing successfully to maximise multi-tool trace coverage. All applicable tools were run together in a single Nextflow invocation per sample. Median and interquartile range (P25–P75) are reported for each tool × modality combination. Timing data were aggregated using `bin/parse_nextflow_trace_timing.py` and stored in `analysis/performance_benchmark_n30/tables/timing_summary_n30.tsv`.

### 8. Orthogonal Silver-Standard Truth Generation and Pangenome Validation

A practical barrier to benchmarking HLA ensembles on new cohorts is that gold-standard HLA typing is rarely available outside a handful of public resources. To enable benchmarking on truth-less cohorts, ChampHLA additionally integrates two methodologically orthogonal, assembly- and pangenome-based tools as optional containerised modules (`modules/locityper.nf`, `modules/immuannot.nf`)—Locityper (Prodanov et al., 2025) and Immuannot (Zhou et al., 2024)—and generates an **agreement-gated silver-standard truth** that is consumed by the same benchmark harness as gold truth. These two tools are deliberately kept outside the eight-tool consensus panel: their role is to construct an independent reference, not to vote. An anti-circularity guardrail (`truth.generated_from` in `bin/hla_benchmark.py`) automatically excludes any tool that contributed to a silver-truth label from being scored against that label.

Locityper genotypes complex loci by recruiting short reads to a locus and selecting the haplotype pair from a pangenome that jointly maximises the read-alignment likelihood, insert-size distribution, and read-depth profile. We constructed a locus database with `locityper target -v` directly from the Human Pangenome Reference Consortium v1.1 Minigraph-Cactus pangenome (Liao et al., 2023; Hickey et al., 2024): the GRCh38 major histocompatibility complex region (chr6:28–34 Mb; 72,315 variants across 45 assembled samples) was extracted with `bcftools`, yielding 47, 70, and 48 real population haplotypes for HLA-A, -B, and -C respectively. Per-sample background read depth was estimated with `locityper preproc` over an N-free autosomal region (chr14:50.0–53.6 Mb). Because Locityper reports genotypes as pangenome haplotype identifiers, each database haplotype was annotated with Immuannot against the bundled IPD-IMGT/HLA reference to obtain a per-locus haplotype→allele crosswalk (for example, the GRCh38 and CHM13 reference haplotypes map to HLA-A\*03:01 and HLA-A\*01:01, respectively). For contrast, a parallel database was built from the full-length IPD-IMGT/HLA 3.64.0 genomic allele set (A: 1,966; B: 892; C: 1,458 alleles after removing partial sequences). The silver-truth reference databases use the latest IPD-IMGT/HLA release available at database-build time (3.64.0); this version affects only the silver-truth pipeline's internal allele set, whereas all gold-truth normalisation and primary accuracy scoring are pinned to IMGT/HLA 3.59.0 (Methods §6).

Silver-standard truth was generated by `bin/generate_silver_truth.py`, which emits a two-field call for a sample and locus only when the Locityper and Immuannot predictions agree and the Locityper genotype quality exceeds a configurable threshold, abstaining otherwise; the output `truth_long.tsv` is schema-compatible with the gold-truth path and required no changes to the benchmark ingestion code. To quantify silver-truth reliability against gold truth, Locityper was run on 99 HLA-region 1000 Genomes WGS samples carrying Gourraud et al. (2014) gold labels; concordance was evaluated at both two-field and G-group resolution for 30 samples (sharing a single depth profile) and for a three-sample pilot with per-sample depth profiles. Reads for the pilot full-WGS samples were obtained by streaming only the MHC and background regions from the remote 30× CRAMs, and all genotyping was executed in containerised SLURM jobs on the CSC Puhti supercomputer.

### 9. External Validation Cohorts

To test generalisation beyond the 1000 Genomes Project—on which the consensus weights and the Champion-Challenger operating points were derived—we assembled an independent, non-1000-Genomes external-validation cohort from the NCI-60 cancer cell-line panel. Constitutional gold-standard HLA genotypes for the panel were taken from the sequence-based typing of Adams et al. (2005, PMC555742) and normalised to the on-disk truth schema (`sample`, `gene`, `allele1`, `allele2`) at two-field resolution for HLA-A, -B and -C; because the Adams typing is mixed-resolution, 42 lines carry two-field-scorable class-I truth (HLA-A *n*=38–39, -B *n*=26, -C *n*=18). Two sequencing modalities were obtained from public archives: RNA-seq from SRA study `SRP133178` (BioProject PRJNA433861; Reinhold et al., 2019) and whole-exome capture from SRA study `SRP150855`, with cell-line names harmonised to the Adams truth identifiers. Deep total-RNA libraries were uniformly downsampled to 50 million read-pairs (streaming `seqtk` subsampling) before typing, to normalise depth across lines and bound per-tool temporary-storage use. Crucially, NCI-60 was treated as a **held-out generalisation test**: we applied the *frozen* 1000G-derived per-modality confidence weights (`tool_weights_rna.json`, `tool_weights_wes.json`) and the already-tuned Champion-Challenger policy (per-modality champions and override gate) without any re-learning on NCI-60, scoring ChampHLA, MajorityVote and the individual tools against the Adams genotype with the same harness used for the 1000G benchmark. As NCI-60 lines are tumour-derived, loss of heterozygosity and aneuploidy can produce apparent HLA homozygosity at some loci; calls were scored against the constitutional Adams genotype and LOH-affected loci flagged through the existing discordance taxonomy rather than removed. To probe the calibration of the frozen WES override gate at scale, we additionally swept the Champion-Challenger operating point on the WES cohort over support-fraction, challenger-margin and supporting-tool thresholds.

To broaden the external test beyond a single cancer-cell-line panel—covering all three modalities and germline as well as tumour-derived material—we assembled three further independent, non-1000-Genomes cohorts with experimental truth, in each case applying the *frozen* 1000G-derived per-modality weights and Champion-Challenger policy without re-learning. First, the **GIAB Ashkenazi trio** (HG002, HG003 and HG004; germline, so free of the loss-of-heterozygosity confound that affects NCI-60) was scored in ChampHLA's two strong modalities against clinical sequence-based typing from the Stanford Blood Center (Table 1 of Llamas et al., 2019, *F1000Research* 8:1751; HG002 cross-confirmed by Chin et al., 2020, *Nat Commun*, PMC7508831, Supplementary Table 4), which we hard-coded into the on-disk truth schema (`bin/build_giab_trio_truth.py`) and validated for Mendelian consistency across the trio. Open ENA reads were used for both arms—RNA-seq runs `SRR15909917`/`SRR15909918`/`SRR15909919` and whole-exome runs `SRR2962669`/`SRR2962692`/`SRR2962694`—applying the frozen RNA policy (champions A=OptiType, B/-C=ArcasHLA; support 0.35, margin 0.0, ≥2 tools) and WES policy (champions A/-B/-C=OptiType; support 0.35, margin 0.00, ≥1 tool). Second, as a whole-genome gold-truth anchor, we scored **GIAB HG002 WGS** from the open `HG002.GRCh38.2x250` alignment (chr6 MHC region sliced remotely without full download) against the same clinical typing, using the frozen WGS policy (champions A=T1K, B/-C=OptiType; support 0.65, margin 0.20, ≥2 tools). Third, we scored the **International HLA and Immunogenetics Workshop (IHWG) MHC-reference cells**—study-verified Illumina WGS from BioProject `PRJNA764575` (SSTO, DBB, APD and QBL, all homozygous MHC typing cells)—against the IPD-IMGT/HLA IHIW multi-laboratory reference-cell consensus, with the same frozen WGS policy. For the IHWG cohort we adopted a read-provenance guardrail: candidate runs were accepted only when the submitting study's context confirmed cell identity, because matching the free-text ENA `cell_line` field alone yielded collisions with unrelated material (embryos, organoids and sarcoma lines) and was discarded.

### 10. Flow-Cytometry HLA-A2 Serotyping (FIMM AML cohort)

As an orthogonal, protein-level ground truth independent of sequence-based typing, surface HLA-A2 expression was measured by flow cytometry with the anti–HLA-A2 monoclonal antibody clone BB7.2 (Parham and Brodsky, 1981) on bone-marrow samples from the FIMM AML cohort (Supplementary S7.1), each scored HLA-A2 positive or negative. Because HLA-A2 is a germline trait, flow samples were matched to ChampHLA calls at the donor level (collapsing sequencing timepoints and replicates), and flow status was verified single-valued within each donor. An HLA-A2–positive genotype was defined as either HLA-A allele in the `A*02` allele group (strict definition); a cross-reactive definition additionally counted `A*68` and `A*69`, reflecting reported BB7.2 recognition of these A2-related antigens. Serotype calls were derived from the multi-tool consensus HLA-A genotype (majority across available modality × tool calls) and, separately, per modality × tool. Concordance with flow cytometry was summarised as sensitivity, specificity and accuracy with Wilson score 95% confidence intervals (`bin/validate_hla_a2_flow.py`; per-donor genotype detail via `bin/champhla_detail_flow_cohort.py`).

---

## Results and Discussion

### 1. HLA Typing Accuracy: ChampHLA versus Individual Tools

All accuracy comparisons in this section were computed under **10-fold, population-stratified nested cross-validation**: the per-gene champions and the override-policy operating point are selected on the training folds only, and every sample is scored strictly out of fold, so no accuracy estimate is inflated by tuning on the data it is measured against (Methods §2). Under this design, ChampHLA's Champion-Challenger consensus delivered its largest — and only statistically significant — accuracy gain in the most difficult modality. In the WGS benchmark (n=137 samples, 411 gene-locus pairs at HLA-A, -B, and -C), Champion-Challenger achieved an overall correct-call rate of **0.5012** (95% CI: 0.453–0.549), exceeding MajorityVote (0.3844) by **11.7 percentage points** and matching the best single tool (OptiType, 0.4975; Table 3). The paired improvement over MajorityVote was highly significant (exact McNemar test: b=60 discordant loci resolved correctly by Champion-Challenger versus c=12 resolved correctly by MajorityVote; **p<0.0001**). The paired comparison against the best single tool, by contrast, was **not** significant (0.5012 vs 0.4975; exact McNemar p=0.74; b=19, c=16): on WGS, champion routing does not so much beat the strongest individual tool as recover the single-tool accuracy ceiling that symmetric majority voting dilutes across the many low-accuracy WGS callers. The mechanism is per-gene champion **routing** rather than the override gate: with the out-of-fold-selected champions (HLA-A and -B: OptiType; HLA-C: T1K) and the conservative override policy, no override fired on any held-out WGS locus, so the entire gain over MajorityVote arises from routing each locus to its benchmark-designated champion instead of a symmetric vote that is diluted by the many low-accuracy WGS tools (ArcasHLA 0.080, Kourami 0.159, SpecHLA 0.165).

On the two high-accuracy modalities, where the strong tools already largely agree, Champion-Challenger **matched** MajorityVote within overlapping confidence intervals, with differences that were not statistically significant. In WES (n=130 samples, 390 gene-locus pairs) it reached 0.9436 (95% CI: 0.916–0.963) versus 0.9359 for MajorityVote (McNemar p=0.45; b=5, c=2) and 0.9231 for the best single tool (OptiType); the override gate fired on a small number of loci and was net corrective. In RNA-seq (n=107 samples, 321 gene-locus pairs) it reached 0.9408 (95% CI: 0.909–0.962) versus 0.9502 for MajorityVote (McNemar p=0.45; b=2, c=5) and 0.9429 for the best single tool, HLA-HD. The consistent picture across modalities is that the benefit of Champion-Challenger **scales with the degree of inter-tool discordance**: it is decisive where tools disagree most and majority voting fails (WGS) and neutral where they concur (WES/RNA). (The RNA-seq scoring run forcibly includes SpecHLA to benchmark against the largest possible tool panel; SpecHLA is excluded from the default production configuration because of its 76.2% RNA callable rate and near-zero RNA accuracy, 0.035.)

This modality pattern is explained by an ablation that isolates the two ingredients of the ensemble. Pure reliability-weighted voting **without** champion routing (WeightedConsensus, Table 3) underperforms MajorityVote in every modality — 0.3017 versus 0.3844 in WGS, 0.8949 versus 0.9359 in WES, and 0.9221 versus 0.9502 in RNA-seq — because it abstains on low-support loci (callable rate 0.58 / 0.93 / 0.96) and, in WGS, spreads weight across discordant tools. Champion routing recovers and exceeds the majority-voting level precisely where weighting alone fails, confirming that the accuracy contribution comes from **routing to a benchmark-validated champion, not from the weighting scheme per se**. The weighting scheme's distinct value lies elsewhere — a calibration guardrail that prevents miscalibrated tool confidence from corrupting the weights (Results §4), and a principled abstention option that raises accuracy among called loci (Supplementary Section S8).

These results are presented in full in **Table 3**, which reports callable rate, accuracy among callable loci, overall correct-call rate, 95% confidence intervals, and the paired McNemar test for all methods across the WGS, WES, and RNA-seq nested-CV benchmarks. **Figure 2** displays the method-level accuracy comparison across modalities, and **Figure 3** shows per-gene (HLA-A, HLA-B, HLA-C) accuracy changes relative to the MajorityVote baseline; extended per-gene accuracy detail for all tools and modalities is provided in **Supplementary Figure S1**. Per-gene accuracy breakdowns and the bimodal/trimodal analyses are presented in **Table 4**. The abstention–accuracy tradeoff between callable rate and accuracy among callable loci for the weighting-only WeightedConsensus is a secondary property, reported in **Supplementary Figure S8**.

---

**Table 1. Supported callers, modalities, evidence types, and benchmark-derived base reliability**

*Base reliability = overall correct-call rate on the respective modality benchmark: the full 137-sample cohort (WGS) or the full truth-backed cohort (WES/RNA-seq).*

| Tool | WGS | WES | RNA-seq | Confidence type | Base reliability (WGS) | Base reliability (WES) | Base reliability (RNA) |
|---|---|---|---|---|---|---|---|
| OptiType | ✓ | ✓ | ✓ | Objective function score | 0.498 | 0.923 | 0.940 |
| HLA-HD | ✓ | ✓ | ✓ | Read-support count | 0.260 | 0.782 | 0.943 |
| T1K | ✓ | ✓ | ✓ | Read-support count | 0.328 | 0.831 | 0.879 |
| ArcasHLA | ✓ | ✓ | ✓ | Native posterior probability | 0.080 | 0.136 | 0.919 |
| SpecHLA | ✓ | ✓ | ✓* | No confidence | 0.165 | 0.785 | 0.035* |
| Kourami | ✓ | ✓ | — | Read-coverage metric | 0.159 | 0.705 | — |
| POLYSOLVER | — | ✓ | — | No confidence | — | 0.915 | — |
| Seq2HLA | — | — | ✓† | No confidence | — | — | 0.581† |

*RNA-seq SpecHLA: 76.2% callable rate; excluded from default production accuracy scoring. †Seq2HLA RNA-seq: n=29 callable samples due to container failure.

---

**Table 2. Cohort design, benchmark role, and split structure**

| Benchmark | Modality | Role | n (cohort) | n (evaluation) | Truth source | IMGT/HLA version | Populations |
|---|---|---|---|---|---|---|---|
| WGS truth-backed | WGS | Primary | 137 | 137 (10-fold nested CV) | Gourraud 2014 | 3.59.0 | Mixed 1000G populations |
| WES truth-backed | WES | Primary | 130 | 130 (10-fold nested CV) | Gourraud 2014 | 3.59.0 | CEU=12, FIN=7, GBR=6, TSI=6, YRI=7, other=92 |
| RNA truth-backed | RNA-seq | Primary | 107 | 107 (10-fold nested CV) | Gourraud 2014 | 3.59.0 | Mixed 1000G populations |
| Trimodal robustness | WGS+WES+RNA | Supplementary | 106 | 106 (matched) | Gourraud 2014 | 3.59.0 | Matched intersection |

*Primary accuracy (Table 3), the per-gene, discordance, and resolution analyses (Tables 4, 6, 7) all use the 137-sample WGS cohort under 10-fold population-stratified nested cross-validation with strictly out-of-fold operating-point selection. A 50-sample WGS "wave2" subset with a frozen 28/10/12 train/validation/holdout split is retained only for the confidence-calibration-guardrail characterisation (Table 5).*

---

**Table 3. Accuracy comparison under 10-fold nested cross-validation — ChampHLA (Champion-Challenger) versus baselines and individual tools across modalities**

**WGS (n=137 samples, 411 gene-locus pairs)** — *Champion-Challenger vs MajorityVote: +11.7 pp, exact McNemar p<0.0001 (b=60, c=12)*

| Method | Type | Callable rate | Acc. among callable | Overall acc. | 95% CI |
|---|---|---|---|---|---|
| OptiType | Single tool | 0.990 | 0.5025 | 0.4975 | [0.449, 0.546] |
| T1K | Single tool | 1.000 | 0.3280 | 0.3280 | [0.283, 0.377] |
| HLA-HD | Single tool | 1.000 | 0.2598 | 0.2598 | [0.218, 0.306] |
| MajorityVote | Baseline | 0.998 | 0.3854 | 0.3844 | [0.339, 0.432] |
| WeightedConsensus (weighting-only ablation) | Ablation | 0.577 | 0.5232 | 0.3017 | [0.259, 0.348] |
| **ChampHLA (Champion-Challenger, nested-CV)** | **Routed** | **0.993** | **0.5049** | **0.5012** | **[0.453, 0.549]** |

**WES (n=130 samples, 390 gene-locus pairs)** — *Champion-Challenger vs MajorityVote: +0.8 pp, McNemar p=0.45 (b=5, c=2), n.s.*

| Method | Type | Callable rate | Acc. among callable | Overall acc. | 95% CI |
|---|---|---|---|---|---|
| OptiType | Single tool | 1.000 | 0.9231 | 0.9231 | [0.892, 0.946] |
| POLYSOLVER | Single tool | 1.000 | 0.9154 | 0.9154 | [0.884, 0.939] |
| MajorityVote | Baseline | 1.000 | 0.9359 | 0.9359 | [0.907, 0.956] |
| WeightedConsensus (weighting-only ablation) | Ablation | 0.931 | 0.9614 | 0.8949 | [0.861, 0.922] |
| **ChampHLA (Champion-Challenger, nested-CV)** | **Routed** | **1.000** | **0.9436** | **0.9436** | **[0.916, 0.963]** |

**RNA-seq (n=107 samples, 321 gene-locus pairs)** — *Champion-Challenger vs MajorityVote: −0.9 pp, McNemar p=0.45 (b=2, c=5), n.s.*

| Method | Type | Callable rate | Acc. among callable | Overall acc. | 95% CI |
|---|---|---|---|---|---|
| HLA-HD | Single tool | 1.000 | 0.9429 | 0.9429 | [0.912, 0.964] |
| OptiType | Single tool | 1.000 | 0.9397 | 0.9397 | [0.908, 0.961] |
| ArcasHLA | Single tool | 1.000 | 0.9190 | 0.9190 | [0.884, 0.944] |
| MajorityVote | Baseline | 1.000 | 0.9502 | 0.9502 | [0.921, 0.969] |
| WeightedConsensus (weighting-only ablation) | Ablation | 0.960 | 0.9610 | 0.9221 | [0.888, 0.947] |
| **ChampHLA (Champion-Challenger, nested-CV)** | **Routed** | **1.000** | **0.9408** | **0.9408** | **[0.909, 0.962]** |

*Champion-Challenger operating points (per-gene champions and override policy) were selected strictly out of fold under 10-fold population-stratified nested cross-validation and pooled over held-out folds (Methods §2; `bin/nested_cv_champion_challenger.py`). MajorityVote has no tunable parameters, so its nested-CV score equals its full-cohort score. WeightedConsensus is shown as the weighting-only ablation (fixed default policy) and abstains on low-support loci, which lowers its overall correct-call rate despite the highest accuracy among callable loci. Paired Champion-Challenger-vs-MajorityVote comparison by exact McNemar test on discordant loci (b = loci ChampHLA calls correctly and MajorityVote does not; c = the reverse). The paired Champion-Challenger-vs-best-single-tool exact McNemar test (`nested_cv_mcnemar_vs_best_tool.tsv`) is not significant on WGS (0.5012 vs OptiType 0.4975; p=0.74; b=19, c=16) or RNA-seq (0.9408 vs HLA-HD 0.9429; p=1.0) and only nominally favours Champion-Challenger on WES (0.9436 vs OptiType 0.9231; p=0.039; b=10, c=2), which does not survive Holm correction across the three modality tests; Champion-Challenger's sole robust paired advantage is therefore over majority voting on WGS. RNA-seq scoring includes SpecHLA to benchmark against the full panel; the default production configuration excludes it (76.2% callable, 0.035 accuracy).*

---

**Figure 1.** *ChampHLA workflow and benchmark governance architecture.* The pipeline accepts BAM, CRAM, or FASTQ inputs and routes each sample through eight containerised HLA typing modules, with automatic input format conversion and tool selection based on sequencing type. Harmonised outputs are aggregated into a unified call schema and passed through the Champion-Challenger ensemble module, which applies benchmark-derived per-tool reliability weights, confidence calibration guardrails, and a nested-cross-validated override policy to produce a final auditable consensus call per gene per sample. The benchmark governance strip at the bottom of the figure indicates the role hierarchy of benchmark cohorts: out-of-fold weight learning and operating-point selection under 10-fold nested cross-validation for primary accuracy reporting, with a retained training split for confidence-calibration characterisation. *(Source: `figures_final/figure_1_workflow_architecture`)*

**Figure 2.** *Accuracy comparison across methods and sequencing modalities (10-fold nested cross-validation).* Bar plots show overall correct-call rate (y-axis) for individual tools (grey), the majority-vote baseline (blue), and ChampHLA Champion-Challenger (red) across the WGS (n=137; left), WES (n=130; centre), and RNA-seq (n=107; right) benchmarks; all Champion-Challenger operating points are selected strictly out of fold. Error bars represent Wilson score 95% confidence intervals. ChampHLA **significantly exceeds majority voting in WGS** — the most discordant modality — through per-gene champion routing (0.5012 vs 0.3844; exact McNemar p<0.0001), thereby recovering the best-single-tool ceiling (OptiType 0.4975; matched, p=0.74), and matches majority voting in the high-accuracy WES and RNA-seq modalities where the strong tools already agree. The contrast illustrates that the ensemble's benefit scales with inter-tool discordance. *(Source: `figures_final/figure_2_accuracy_comparison`)*

**Figure 3.** *Per-gene accuracy gains and losses relative to MajorityVote baseline.* Signed bar plots show the change in overall correct-call rate for each method relative to MajorityVote, separately for HLA-A, HLA-B, and HLA-C, across all three modalities. Positive bars indicate methods that outperform MajorityVote at that locus; negative bars indicate regressions. ChampHLA Champion-Challenger shows its largest positive per-gene gains over MajorityVote in WGS, where champion routing recovers the single-tool ceiling that symmetric voting dilutes; in WES and RNA-seq the per-gene changes are small, consistent with the modality-level parity. *(Source: `figures_final/figure_3_per_gene_gains`)*

### 2. Computational Resource Requirements

The eight HLA typing tools span a wide range in peak RAM and wall-clock time across sequencing modalities (**Supplementary Figure S14**; measured from n=30 1000G samples per modality via Nextflow execution trace, see Methods §7; RNA-seq tool-level n varies from 26 to 30 due to individual tool failures). All values are medians; error whiskers show the P25–P75 interquartile range.

In the WES benchmark, **ArcasHLA** was the fastest tool (median 0.008 h wall-clock, 0.55 GB peak RAM, ~207% CPU), followed by **OptiType** (0.013 h, 0.59 GB, ~101% CPU). **POLYSOLVER** required the most wall-clock time in WES (0.105 h, 1.5 GB), while **Kourami** (0.016 h, 3.9 GB) and **T1K** (0.034 h, 3.9 GB) had the highest peak RAM demands in WES. **SpecHLA** required 0.089 h and 0.96 GB in WES mode. **HLA-HD** ran in 0.015 h with 2.2 GB peak RAM, benefiting from its multi-threaded alignment strategy (~789% CPU utilisation). These WES measurements reflect HLA-region-extracted BAM inputs rather than full exome BAMs, which accounts for the relatively short wall-clock times compared to whole-genome inputs.

In the RNA-seq benchmark, wall-clock times are longer than WES, reflecting the read volume and alignment complexity of transcriptomic data processed from HLA-region-extracted BAMs. **ArcasHLA** required a median of 0.034 h and 1.5 GB peak RAM; **Seq2HLA** was the most memory-efficient RNA tool (0.049 h, 0.14 GB, ~469% CPU). **HLA-HD** was the most time-consuming RNA-seq tool (0.452 h, 8.7 GB), followed by **SpecHLA** (0.393 h, 1.8 GB, n=29). **T1K** (0.150 h, 10.0 GB) had the highest peak RAM among RNA tools. **OptiType** required 0.098 h and 5.3 GB peak RAM (n=26; four samples exceeded the 8 h wall-clock limit). Interquartile ranges indicate moderate sample-to-sample variability in RNA-seq timing, consistent with variable read depth and expressed allele coverage.

In the WGS benchmark, wall-clock times were comparable to WES due to HLA-region BAM extraction prior to tool execution. **OptiType** remained the fastest tool (median 0.010 h, 0.45 GB, ~98% CPU), followed by **ArcasHLA** (0.021 h, 0.94 GB, ~213% CPU). **HLA-HD** (0.033 h, 1.9 GB, ~411% CPU), **POLYSOLVER** (0.034 h, 0.92 GB, ~101% CPU), and **T1K** (0.031 h, 2.6 GB, ~550% CPU) clustered in the 0.03–0.05 h range. **Kourami** required the highest peak RAM in WGS mode (5.9 GB, 0.051 h, ~1065% CPU), while **SpecHLA** was the slowest WGS tool (0.243 h, 1.5 GB). Interquartile ranges were narrow for WGS, reflecting the consistency of HLA-region-extracted inputs.

These profiles inform deployment planning: environments with fewer than 16 GB RAM can comfortably run OptiType, ArcasHLA, HLA-HD, POLYSOLVER, and SpecHLA across all applicable modalities. Kourami (5.9 GB WGS), T1K (10.0 GB RNA), and HLA-HD in RNA mode (8.7 GB) have the highest peak RAM demands. For time-constrained workflows, the OptiType+ArcasHLA+HLA-HD panel delivers ensemble consensus at reduced wall-clock overhead. Full per-tool resource profiles across all modalities and metrics are shown in **Supplementary Figure S14**.

### 3. Result Concordance and Heterogeneity Across Tools

A central motivation for ChampHLA is that different HLA typing tools frequently predict different allele pairs for the same sample. To quantify this heterogeneity, ChampHLA computes structured discordance labels for every gene-locus pair and writes them to `discordance_summary.tsv` (**Supplementary Figure S12**; **Table 6**). Across the three modalities the discordance burden tracks exactly the ordering in which Champion-Challenger helps: it is highest in WGS, where the ensemble delivers a large accuracy gain, and lowest in RNA-seq, where the ensemble is neutral.

WGS is by a wide margin the most discordant modality: **174 discordance events across 411 gene-locus pairs (42.3%)**—173 `low_evidence_conflict` and one `no_evidence`—reflecting the wide spread of per-tool WGS accuracy (0.08–0.50), so that many loci carry no allele pair with strong weighted support. This is precisely the regime in which routing each locus to its benchmark-designated champion, rather than taking a symmetric vote across mostly weak tools, produces the +11.7-percentage-point accuracy gain reported in Results §1.

RNA-seq, in contrast, is the **least** discordant modality: only 13 `low_evidence_conflict` events across 321 gene-locus pairs (4.0%), because the strong RNA tools (HLA-HD, OptiType, ArcasHLA at 0.92–0.95) largely agree. Where intra-RNA disagreement does occur it is consistent with allele-specific expression at heterozygous loci—one allele expressed at lower read depth is resolved differently by tools with different sensitivity thresholds—but at this low rate there is little for an ensemble to correct, consistent with Champion-Challenger matching majority voting on RNA-seq.

WES sits between the two: 27 `low_evidence_conflict` events across 390 gene-locus pairs (6.9%). Under nested cross-validation the Champion-Challenger override gate acted on a small number of these discordant loci and was net corrective (Results §1); every override is preserved in the audit trail (`champion_challenger_method_comparison.tsv`) for post-hoc review without altering the primary accuracy statistics.

**Table 6.** *Discordance summary by modality and category (10-fold nested-CV cohorts)*

| Modality | low_evidence_conflict | no_evidence | possible_expression_bias | Total events | Discordance rate |
|---|---|---|---|---|---|
| WGS (n=411 loci) | 173 | 1 | 0 | 174 | 42.3% |
| WES (n=390 loci) | 27 | 0 | 0 | 27 | 6.9% |
| RNA-seq (n=321 loci) | 13 | 0 | 0 | 13 | 4.0% |

*RNA-seq rate exceeds 100% because multiple discordance categories can be tagged per locus. Category counts extracted from `discordance_summary.tsv` per benchmark root (IMGT/HLA 3.59.0).*

### 4. Confidence Calibration and Benchmark-Derived Weights

The confidence-calibration characterisation in this section was computed on the full 137-sample WGS cohort, consistent with the primary benchmark. Because the ablation (Results §1) shows that weighting does not drive the reported accuracy, the guardrail functions as a soundness safeguard on the ensemble weights rather than as a contributor to the correct-call rates. The calibration guardrail system revealed that raw tool confidence signals are not directly comparable across HLA callers and cannot safely be used to amplify ensemble weights without empirical validation. In the WGS cohort, five of eight tools contributed confidence signals; three of five were blocked by the guardrail. OptiType's integer-linear-programming objective score saturated at 1.0 across virtually all calls regardless of whether the call was correct (Brier score=0.503, ECE=0.503), indicating that OptiType's confidence proxy is uninformative for weight boosting despite the tool's strong base accuracy (0.498). HLA-HD similarly showed saturation-driven miscalibration (Brier=0.550, ECE=0.583), and Kourami showed the most extreme miscalibration of all five tools (Brier=0.817, ECE=0.827). Under the guardrail, all three are assigned final weights equal to their base reliabilities (OptiType=0.498, HLA-HD=0.260, Kourami=0.159).

The practical stakes of this miscalibration are best illustrated by a counterfactual. Without the guardrail, OptiType's saturated mean confidence of approximately 1.0 would push its blended final weight to approximately 0.7 × 0.4975 + 0.3 × 1.0 = 0.648—a ~30% relative increase over its base reliability, driven entirely by an uninformative proxy that is equally high for correct and incorrect calls. HLA-HD, with mean confidence 0.876 but base reliability 0.260, would be similarly over-weighted relative to its actual accuracy. The calibration guardrail prevents both inflations, ensuring that the weight ordering correctly reflects the tool accuracy ranking rather than confidence-proxy magnitude (**Supplementary Figure S13**).

Two tools passed the guardrail: T1K (Brier=0.312, ECE=0.297), which retained a partial confidence boost resulting in a final weight of 0.339 (versus base reliability 0.328), and ArcasHLA (Brier=0.079, ECE=0.076), which showed excellent calibration but contributed negligible weight due to near-zero WGS accuracy (base reliability=0.080; ArcasHLA is an RNA-seq-optimised tool applied out of distribution on WGS data). These calibration results are summarised in **Table 5** and visualised in **Figure 4**. A per-tool calibration heatmap showing bin-level predicted-versus-observed accuracy across all confidence-producing tools is provided in **Supplementary Figure S2**.

**Table 4.** *Per-gene HLA-A, HLA-B, HLA-C accuracy — WGS (137-sample nested CV) and bimodal/trimodal robustness*

**WGS (n=137 samples, 10-fold nested CV; per-gene overall correct-call rate)**

| Method | HLA-A | HLA-B | HLA-C |
|---|---|---|---|
| MajorityVote | 0.3577 | 0.4964 | 0.2993 |
| Best single tool | 0.4853 (OptiType) | 0.6985 (OptiType) | 0.3413 (T1K) |
| **ChampHLA (Champion-Challenger)** | **0.4818** | **0.6934** | **0.3285** |
| CC vs MV, exact McNemar *p* | **<0.001** | **<0.001** | 0.48 (n.s.) |

*The WGS accuracy gain concentrates at HLA-A (+12.4 percentage points) and HLA-B (+19.7 pp), both highly significant; HLA-C — the hardest WGS locus for every tool — shows a small, non-significant gain. In each gene Champion-Challenger tracks the per-gene single-tool ceiling (routing), confirming that the modality-level gain in Table 3 is the sum of per-gene champion routing rather than the override gate.*

**Bimodal WES+RNA (n=106 matched samples)**

| Method | HLA-A (overall) | HLA-B (overall) | HLA-C (overall) |
|---|---|---|---|
| BimodalMajorityVote | 0.9717 | 0.9623 | 0.9528 |
| **BimodalChampionChallenger (nested-CV)** | 0.9623 | 0.9528 | 0.9717 |
| BimodalWeightedConsensus | 0.9434 | 0.9528 | 0.9245 |

*Bimodal Champion-Challenger (10-fold nested CV over the 106 matched subjects) matches BimodalMajorityVote overall (0.9623 vs 0.9623; exact McNemar p=1.0, two discordant loci each way) and exceeds the best single tool (0.9423); per gene it trails majority voting at HLA-A/-B and leads at HLA-C, a net wash (Results §7, Figure 5).*

**Table 5.** *Confidence calibration, guardrail status, and final weights (137-sample WGS cohort)*

| Tool | Confidence type | Brier score | ECE | Guardrail | Final weight (WGS) |
|---|---|---|---|---|---|
| ArcasHLA | Posterior probability | 0.079 | 0.076 | Pass | 0.081 |
| T1K | Read-support proxy | 0.312 | 0.297 | Pass | 0.339 |
| OptiType | Objective score | 0.503 | 0.503 | Fail (blocked) | 0.498 |
| HLA-HD | Read-support proxy | 0.550 | 0.583 | Fail (blocked) | 0.260 |
| Kourami | Read-coverage proxy | 0.817 | 0.827 | Fail (blocked) | 0.159 |
| SpecHLA | No confidence | — | — | No confidence | 0.165 |

*Guardrail threshold: Brier and ECE both < 0.35 required for confidence boost. Blocked tools: effective_confidence = base_reliability, giving final_weight = base_reliability. Computed on the full 137-sample WGS cohort; the guardrail decisions are identical to those obtained on the wave2 training split, confirming that the calibration ranking is cohort-stable.*

**Figure 4.** *Confidence calibration curves and guardrail outcomes by tool.* Calibration plots (reliability diagrams) show mean predicted confidence (x-axis) versus observed accuracy (y-axis) in 10 equal-frequency bins for each tool with a non-null confidence signal. The diagonal represents perfect calibration. OptiType, HLA-HD, and Kourami cluster near the top-right corner (mean confidence ≈ 1.0 regardless of accuracy), indicating severe proxy saturation. T1K and ArcasHLA fall near the diagonal, indicating usable calibration. Blocked tools are highlighted with a red border; passing tools with a green border. *(Source: `figures_final/figure_4_confidence_calibration`)*

### 5. Resolution Analysis

Secondary resolution metrics—three-field exact match rate (exact_3field), G-group compatibility rate, P-group compatibility rate, and ambiguity-compatible rate—were evaluated for all tools across the WGS (n=137) and WES truth-backed benchmarks using `method_comparison_multiresolution.tsv` and `summary_full_cohort_multiresolution.tsv`. Four key findings emerge.

**Finding 1: Only OptiType returns verified ≥3-field allele calls.** For OptiType in WGS, exact_3field = exact_2field = 0.4975, confirming that OptiType's native output format encodes alleles at three-field resolution and that the 2014 truth set supports verification at this level. All other tools return essentially two-field outputs: exact_3field ≈ 0 for ArcasHLA, T1K, and SpecHLA, with only marginal non-zero values (HLA-HD 0.0105 and Kourami 0.0082 in WGS; POLYSOLVER 0.0026 in WES), attributable to incidental allele-level coincidence rather than true three-field resolution. For three-field or four-field clinical applications, OptiType is currently the only tool in the ChampHLA panel whose output contributes at extended resolution.

**Finding 2: Ambiguity-compatible analysis reveals that apparent WGS accuracy gaps are substantially an artefact of the 2014 truth set.** When calls are evaluated using ambiguity-compatible matching—where a call is scored correct if it falls within the IMGT/HLA ambiguity equivalence class of the truth allele—apparent WGS accuracy increases substantially for all methods. MajorityVote: exact_2field = 0.3844 versus ambiguity_compatible = 0.7932, a 40.9 percentage-point gap. OptiType: 0.4975 exact versus 0.8775 ambiguity_compatible. Kourami: 0.1592 exact versus 0.5429 ambiguity_compatible in WGS. These gaps indicate that a large proportion of apparent WGS typing errors are cases where the tool's output is a valid allele within the ambiguity group of the truth allele, but does not match the specific representative allele recorded in the 2014 truth set. The ambiguity-compatible rate therefore provides a more accurate estimate of clinical utility for WGS typing workflows.

**Finding 3: G-group and P-group match rates are uniformly zero across all tools and modalities.** g_group_match_rate = 0.0 and p_group_match_rate = 0.0 for all tools in all modalities. This is the expected truth-ceiling result: the 2014 1000G truth set carries no G-group or P-group suffix annotations. Evaluation of G-group and P-group accuracy requires a more recent truth set derived from high-resolution Sanger sequencing or long-read validation.

**Finding 4: POLYSOLVER output is capped at two-field resolution despite WES optimisation.** In the WES benchmark, POLYSOLVER achieves exact_2field = 0.9154 but exact_3field = 0.0026. This confirms that POLYSOLVER's output format is two-field-truncated by design, consistent with published documentation. For applications requiring three-field or four-field WES typing, POLYSOLVER cannot contribute beyond the two-field ensemble.

These resolution findings are summarised in **Table 7** and visualised in **Supplementary Figure S3**.

**Table 7.** *Resolution metric summary — WGS (137-sample cohort) and WES truth-backed benchmark*

| Tool | Modality | exact_2field | exact_3field | ambiguity_compatible | g_group | p_group |
|---|---|---|---|---|---|---|
| OptiType | WGS | 0.4975 | 0.4975 | 0.8775 | 0.0 | 0.0 |
| MajorityVote | WGS | 0.3844 | — | 0.7932 | 0.0 | 0.0 |
| WeightedConsensus | WGS | 0.3017 | — | 0.5450 | 0.0 | 0.0 |
| Kourami | WGS | 0.1592 | 0.0082 | 0.5429 | 0.0 | 0.0 |
| HLA-HD | WGS | 0.2598 | 0.0105 | 0.1942 | 0.0 | 0.0 |
| T1K | WGS | 0.3280 | 0.0 | 0.2725 | 0.0 | 0.0 |
| OptiType | WES | 0.9231 | 0.9231 | — | 0.0 | 0.0 |
| POLYSOLVER | WES | 0.9154 | 0.0026 | — | 0.0 | 0.0 |
| MajorityVote | WES | 0.9359 | 0.0 | — | 0.0 | 0.0 |

*exact_3field = exact_2field for OptiType only (sole tool with three-field native output). All others: exact_3field ≈ 0. g_group and p_group = 0.0 uniformly across all tools and modalities (truth-ceiling: 2014 1000G truth carries no G-group/P-group annotations). ambiguity_compatible shown for WGS where the gap from exact_2field is largest and most clinically informative. Sources: `method_comparison_multiresolution.tsv`, `summary_full_cohort_multiresolution.tsv` per benchmark root (IMGT/HLA 3.59.0).*

Primary accuracy in all ChampHLA benchmarks was evaluated using exact two-field allele-pair match rate against the 2014 1000 Genomes truth set, which provides ground-truth calls at two-field resolution for HLA-A, -B, and -C for most samples. Secondary accuracy metrics—three-field exact match, G-group compatibility, P-group compatibility, and ambiguity-compatible rate—were computed for all tools but are sparsely populated for the majority of WGS loci because the 2014 truth set does not carry G-group or P-group suffix notations for most alleles, and most loci do not have three-field truth encodings. As a result, secondary metric columns are blank for the majority of WGS rows, and all cross-method comparisons in this study use only two-field exact-match accuracy. Secondary resolution comparison results are presented in **Supplementary Figure S3** for the subset of loci where truth encodings permit.

This limitation is important to contextualise. The two-field resolution at which all primary benchmarks are evaluated is clinically relevant for most transplantation matching decisions at the allele-group level; however, high-resolution clinical typing for HSCT matching increasingly requires four-field or at minimum three-field resolution at HLA-B and HLA-C loci where allele-level mismatches have stronger clinical associations than group-level mismatches. The extent to which different tools—and ChampHLA—perform at three-field and four-field resolution in WES and RNA-seq modalities with truth sets supporting those resolutions is a high-value question for future benchmark cycles.

**Robustness across allele commonness and ancestry.** A concern for any consensus HLA caller is that a headline accuracy gain may be dominated by common alleles, which are over-represented in every cohort and easiest to type, or may not generalise across ancestries. We tested both on the held-out nested-CV calls. For commonness, we stratified each genotype by the rarer of its two truth alleles using the **Common, Intermediate and Well-Documented (CIWD) version 3.0.0** catalogue (Hurley et al., 2020, *HLA* 95:516–531; an allele-frequency classification compiled from ~8 million donor-registry typings—not a ground-truth genotype resource, and used here purely as an annotation layer). The WGS Champion-Challenger gain is **not** an artefact of rare alleles: within the common-allele stratum (n=397 loci) Champion-Challenger improves over majority voting from 0.390 to 0.504—essentially the full +11-point margin seen overall—while on the common-dominated WES and RNA-seq benchmarks the two methods stay within one point (WES 0.955 vs 0.952; RNA-seq 0.946 vs 0.955). For ancestry, we stratified held-out accuracy by 1000 Genomes superpopulation mapped to continental groups (EUR/AFR): Champion-Challenger was greater than or equal to majority voting in every ancestry stratum and modality tested—most visibly for the WGS African-ancestry (0.571 vs 0.476) and unlabelled (0.515 vs 0.364) subsets—so routing to a benchmark-designated champion does not disadvantage any ancestry group. The small African-ancestry sample (n=21 loci) and incomplete ancestry labelling of the WGS cohort mean per-population weight learning (Limitations) remains a valuable extension rather than a resolved question. Full breakdowns are in **Supplementary Table S9** (commonness) and **Supplementary Table S10** (ancestry) (`analysis/nested_cv_champion_challenger/*/nested_cv_by_{ciwd,ancestry}.tsv`). As a complementary quality-control view, any callable call whose allele is not-CIWD or absent from the catalogue is flagged as biologically implausible—a candidate typing error or genuinely novel allele (`summary_ciwd_plausibility.tsv`). These stratifications are additive: the overall two-field figures reported elsewhere are unchanged. *(Catalogue vendored at `assets/ciwd_3.0.0.tsv`; classifier `bin/ciwd.py`; see `assets/README_CIWD.md`.)*

### 6. Orthogonal Silver-Standard Validation with the HPRC Pangenome

The resolution ceiling described in §5—the absence of G-group and high-resolution truth in the 2014 1000 Genomes panel—motivates an independent, assembly-derived reference. ChampHLA's orthogonal silver-truth pipeline (Locityper genotyping against an HPRC pangenome database, with Immuannot supplying the haplotype→allele crosswalk) reproduces 1000 Genomes gold calls with high concordance and supplies exactly the kind of G-group-resolved, long-read-derived validation that §5 identifies as missing. Against gold truth, the HPRC-pangenome database achieved **94.4%** allele-level concordance on the three-sample pilot with per-sample depth profiles (HLA-A 83.3%, HLA-B and HLA-C both 100%) and **81.1% two-field / 82.2% G-group** concordance across 30 samples sharing a single depth profile (per-locus two-field: HLA-A 73.3%, HLA-B 76.7%, HLA-C 93.3%; **Table 8**, **Supplementary Figure S16A**). The strong agreement at G-group resolution—where the 1000 Genomes gold panel itself carries no annotations—confirms that the silver-truth pipeline recovers clinically meaningful allele groups rather than merely two-field strings.

Database construction, rather than the genotyping algorithm, was the dominant determinant of accuracy. A control database built from the full IPD-IMGT/HLA allele set reached only **46.1% two-field / 50.6% G-group** across the same 30 samples (55.6% on the pilot), because the thousands of near-identical and partial genomic alleles drive Locityper to over-resolve calls to rare neighbouring alleles. Replacing the allele set with the 47–70 real assembled HPRC haplotypes per locus nearly doubled concordance under identical reads, depth profiles, and software. This result is a practical caution for assembly-based truth derivation generally: the choice of reference panel matters more than the choice of genotyper.

The two databases also differed sharply in confidence calibration (**Supplementary Figure S16B**), with a direct bearing on the silver-truth design. Under the HPRC database, Locityper genotype quality tracked correctness—concordance rose monotonically from 0.811 over all calls to 0.905 for the highest-confidence subset (genotype quality ≥ 20)—whereas under the IMGT-allele database, quality was uninformative or anti-correlated (0.461 over all calls versus 0.432 for genotype quality ≥ 10): the genotyper was *confidently wrong*. This miscalibration under a weak reference is precisely why ChampHLA gates silver truth on agreement between two methodologically orthogonal tools rather than on any single tool's self-reported confidence. Orthogonal agreement remains a valid safeguard even when an individual tool's quality score does not, which is the central design rationale for the agreement-gated truth generator.

Two factors explain the gap between the 94.4% pilot and the 81.1% 30-sample concordance, and both are addressable. First, the 30-sample run shared a single depth profile, whereas the pilot used per-sample preprocessing and recovered roughly 13 percentage points; region-restricted input BAMs that cannot be preprocessed individually therefore understate achievable accuracy. Second, the 45-sample HPRC v1.1 panel does not carry every rare allele—for example, a gold HLA-A\*29:02 allele absent from the panel was assigned the nearest available haplotype (A\*32:01)—and the forthcoming HPRC v2 release (200+ samples) is expected to close most of the residual HLA-A gap. Even at current scale, the agreement-gated pipeline is accurate enough to serve as silver-standard truth for benchmarking ChampHLA on cohorts lacking gold-standard HLA typing, extending the framework beyond the public 1000 Genomes panel.

**Table 8. HPRC pangenome silver-standard truth concordance versus 1000 Genomes gold truth.** *Per-locus values are exact two-field allele concordance; overall is reported at both two-field and G-group resolution. Locityper genotyping with the HPRC v1.1 pangenome database versus a control database built from IPD-IMGT/HLA 3.64.0 alleles, evaluated on 30 1000 Genomes WGS samples (shared depth profile) and a three-sample per-sample-depth pilot.*

| Database | Samples | Depth profile | HLA-A | HLA-B | HLA-C | Overall (2-field) | Overall (G-group) |
|---|---|---|---|---|---|---|---|
| HPRC pangenome | 30 | shared | 0.733 | 0.767 | 0.933 | **0.811** | **0.822** |
| HPRC pangenome | 3 (pilot) | per-sample | 0.833 | 1.000 | 1.000 | **0.944** | **0.944** |
| IPD-IMGT/HLA alleles | 30 | shared | 0.333 | 0.550 | 0.500 | 0.461 | 0.506 |
| IPD-IMGT/HLA alleles | 3 (pilot) | per-sample | 0.333 | 0.833 | 0.500 | 0.556 | 0.556 |

### 7. Cross-Modality Robustness: Bimodal and Trimodal Analysis

The supplementary 106-sample matched-subject trimodal robustness analysis tested whether combining HLA evidence across WGS, WES, and RNA-seq improves ensemble accuracy over any single modality on the same subjects. The results were unambiguous: bimodal WES+RNA majority voting achieved the best multimodal result at **0.9623** overall correct-call rate (callable rate=0.9937; 95% CI: 0.935–0.978; **Figure 5**), outperforming all unimodal methods. At the per-gene level, bimodal MajorityVote achieved 0.9717 (HLA-A), 0.9623 (HLA-B), and 0.9528 (HLA-C), with near-complete callable coverage at all three loci. These results confirm that WES and RNA-seq are mutually reinforcing modalities for ensemble HLA typing. Applying the Champion-Challenger consensus to the same bimodal call set under 10-fold nested cross-validation reproduced the majority-voting result exactly (0.9623 vs 0.9623; exact McNemar p=1.0, two discordant loci each way) while exceeding the best single tool (HLA-HD RNA-seq, 0.9423; Table 4): in the recommended bimodal configuration Champion-Challenger therefore matches majority voting without regression and beats any individual tool, so it can serve as the default consensus—providing calibration-guarded, auditable calls at no accuracy cost—rather than being reserved for the discordant WGS modality where it produces a significant gain.

Adding WGS as a third modality did not improve the bimodal headline: TrimodalMajorityVote achieved **0.9591** overall correct-call rate (callable rate=0.9906; **Supplementary Figure S15**), 3.2 percentage points below the bimodal WES+RNA result. Per-gene trimodal accuracy (A=0.9717, B=0.9623, C=0.9434) is largely unchanged from the bimodal result at HLA-A and HLA-B but shows a marginal improvement at HLA-C (+0.9 percentage points) that does not offset the overall dilution. This negative result for WGS addition is itself informative: it demonstrates that short-read WGS HLA typing in its current form adds more noise than signal when combined with an already-strong WES+RNA bimodal ensemble. The most accurate multimodal strategy for HLA typing with current short-read technology is therefore WES+RNA bimodal consensus, with WGS reserved for cases where WES and RNA-seq data are unavailable.

**Figure 5.** *Bimodal WES+RNA per-gene comparison (n=106 matched subjects; 10-fold nested cross-validation).* Grouped bars show overall correct-call rate at HLA-A, HLA-B, and HLA-C for bimodal MajorityVote, bimodal Champion-Challenger, and the best single tool (per gene, across WES and RNA-seq callers). Overall, bimodal Champion-Challenger matches bimodal majority voting exactly (0.9623 vs 0.9623; exact McNemar p=1.0) and both exceed the best single tool (0.9423); per gene the two consensus methods trade the lead—majority voting at HLA-A and -B, Champion-Challenger at HLA-C—a net wash. In the recommended bimodal configuration Champion-Challenger can therefore serve as the default consensus at no accuracy cost while adding calibration-guarded, auditable calls. *(Source: `figures_final/figure_09_bimodal_per_gene`)*

**Figure 6.** *How the Champion-Challenger consensus works and how it differs from majority voting.* **(A)** Mechanism schematic: for each HLA gene a benchmark-designated champion (selected per gene and modality — e.g. OptiType for HLA-A/-B and T1K for HLA-C in WGS) provides the default call, while a reliability- and confidence-weighted ensemble (`final_weight = 0.7 × base_reliability + 0.3 × effective_confidence`, gated by a Brier/ECE < 0.35 calibration guardrail) may override the champion only when all four gates pass simultaneously — a minimum challenger support fraction, weight margin, number of supporting tools (thresholds selected strictly out of fold), and an unambiguous challenger allele — producing one of three audited decision-trace outcomes (`challenger_override`, `champion_retained`, `champion_missing` fallback). Majority voting, by contrast, treats every callable tool as one equal vote and returns `no_call` on ties, with no champion, reliability weights, or calibration. **(B)** Under 10-fold nested cross-validation, Champion-Challenger's advantage over majority voting scales with inter-tool discordance: it is decisive in WGS (0.5012 vs 0.3844; +11.7 percentage points; exact McNemar p<0.0001), where per-gene champion routing recovers the single-tool ceiling that symmetric voting dilutes, and neutral in the high-accuracy WES (0.9436 vs 0.9359) and RNA-seq (0.9408 vs 0.9502) modalities where the strong tools already agree. All values are pooled over held-out folds. *(Source: `figures_final/figure_11_champion_challenger_combined`; alternate layouts `figure_11a_cc_mechanism`, `figure_11b_mv_vs_cc`)*

### 8. External Validation on Independent Non-1000G Cohorts

A central question for any benchmark-tuned ensemble is whether it generalises beyond the cohort on which it was calibrated. We therefore applied the frozen 1000G-derived weights and Champion-Challenger policy, without re-learning, to the independent NCI-60 cancer cell-line panel scored against the Adams et al. (2005) sequence-based constitutional genotypes (Supplementary Figure S17).

On **RNA-seq** (powered, *n*=42; HLA-A *n*=38, -B *n*=26, -C *n*=18), the frozen policy transferred as an actively corrective mechanism but did not beat the baseline. ChampHLA reached an overall correct-call rate of 0.793 (per-locus HLA-A 0.763, -B 0.923, -C 0.667; Wilson 95% CI 0.693–0.866), tracking just below majority voting (0.817; CI 0.720–0.886) and its best single tool OptiType (0.805). Importantly, the override gate stayed strongly corrective out of distribution: 20 overrides fired, of which **15 were corrective** (concentrated at HLA-B and -C), 1 harmful and 4 neutral—93.8% corrective precision among non-neutral decisions. The reason ChampHLA lands marginally below majority voting rather than above it is structural: the frozen RNA policy assigns the HLA-B/-C champion to ArcasHLA (the strongest B/C tool on the 1000G training cohort, but comparatively weak on NCI-60, where it scores 0.494 overall), so the champion baseline for two of three loci is low and the corrective overrides recover much—but not quite all—of that gap. The honest reading, consistent with WES, is that the Champion-Challenger gate generalises as a predominantly corrective mechanism while overall accuracy matches majority voting within overlapping confidence intervals; a modality-specific re-selection of the RNA HLA-B/-C champion (rather than freezing the 1000G choice) would be the lever to recover the shortfall, mirroring the WES supporting-tools-floor result below.

On **WES** (n=42; HLA-A *n*=38, -B *n*=26, -C *n*=18), the deployed policy also generalised, with an instructive caveat. At the operating point ChampHLA actually ships (support fraction 0.35, challenger margin 0.00, ≥1 supporting tool—the WES optimum from the 1000G sweep), ChampHLA reached an overall correct-call rate of 0.842, exactly matching majority voting (0.842; Wilson 95% CI 0.747–0.905) and above its champion OptiType (0.817). Importantly, the override gate remained active out of distribution: it fired seven overrides, of which three were corrective—all at HLA-A, where it correctly resolved apparent-homozygous calls in EKVX, NCI-H522 and TK-10 that OptiType had split—while one was harmful and three were neutral, all at HLA-C and all involving spuriously high-resolution challenger alleles proposed by a single tool. The net effect was parity with majority voting. A sweep of the operating point on the WES cohort (Supplementary Figure S17C) localised the issue to the supporting-tools floor rather than the support or margin gates: holding support=0.35 and margin=0.00 and raising the floor from one to two supporting tools eliminates every harmful and neutral override (the spurious single-tool HLA-C calls) while retaining all three corrective HLA-A overrides, lifting ChampHLA to 0.854 and overtaking majority voting; requiring three supporting tools is over-strict and suppresses all overrides (reverting to OptiType, 0.817).

To test whether this behaviour holds beyond a single cancer-cell-line panel, we scored three further independent non-1000G cohorts under the same frozen policies, extending the external check to germline material (free of the loss-of-heterozygosity confound) and to all three modalities. The overall picture across all six external arms is summarised in **Table 9**.

On the **GIAB Ashkenazi trio** (HG002/HG003/HG004; germline, clinical sequence-based gold truth) in ChampHLA's two strong modalities, the frozen policies transferred cleanly. On **WES** ChampHLA matched majority voting at a perfect correct-call rate (1.000, 9/9; OptiType and HLA-HD likewise 9/9), with all champions retained. On **RNA-seq** ChampHLA again matched majority voting (0.889, 8/9), with HG002 and HG004 typed perfectly at all six alleles; the single miss was a second-field slip at HG003 HLA-C (`06:08` vs truth `06:02`) on which both the frozen HLA-C champion and OptiType err. On this cleanest, LOH-free, gold-truth cohort ChampHLA is therefore perfect on WES and near-perfect on RNA-seq, consistent with the NCI-60 reading that the framework transfers out of distribution and tracks majority voting.

As a whole-genome gold-truth anchor, **GIAB HG002 WGS** (scored against the same clinical typing) was recovered exactly by ChampHLA at all three class-I loci (1.000, 3/3—indeed at three-field resolution: A\*01:01/26:01, B\*35:08/38:01, C\*04:01/12:03), equal to majority voting; all three champions were already correct so no override was needed, and the only single-tool error (HLA-HD at HLA-C) was avoided by the OptiType champion routing. This is a clean recovery in ChampHLA's weakest modality, though n=1 makes it an anchor rather than a powered comparison.

On the **IHWG MHC-reference cells** (study-verified WGS of the homozygous typing cells SSTO, DBB, APD and QBL, scored against the IPD-IMGT/HLA multi-laboratory consensus), ChampHLA reached 0.833 (10/12) against 0.917 (11/12) for majority voting, with three of four cells typed perfectly at all six alleles. As on NCI-60 RNA, the entire shortfall traces to a single frozen, modality-specific champion: the WGS HLA-A champion T1K is strong on 1000G WGS but weak on these reference cells (HLA-A 0.50), whereas the OptiType-championed HLA-B and -C are perfect (1.00); the one miss is a second-field slip at QBL's homozygous HLA-A that the (weaker) T1K would have recovered. The result is a non-1000G, reference-grade external check that reproduces the cross-cohort pattern, with the caveats that these are homozygous typing cells (which inflate agreement), carry mild reference-allele circularity, and give a wide confidence interval at n=4.

**Table 9. Cross-cohort external validation under frozen 1000G policies (overall two-field correct-call rate).** All arms apply the frozen 1000G-derived per-modality weights and Champion-Challenger policy without re-learning to independent, non-1000G cohorts with experimental truth. NCI-60 is the powered comparison; the germline arms are small-n confirmatory checks (wide confidence intervals). *n* is samples / scorable HLA-A/-B/-C gene-rows.

| Cohort | Material | Modality | *n* (samples / gene-rows) | ChampHLA (CC) | MajorityVote | Truth source |
|---|---|---|---|---|---|---|
| NCI-60 | cancer cell lines | RNA-seq | 42 / 81 | 0.793 | 0.817 | Adams et al. (2005) SBT |
| NCI-60 | cancer cell lines | WES | 42 / 82 | 0.842 | 0.842 | Adams et al. (2005) SBT |
| GIAB Ashkenazi trio | germline | RNA-seq | 3 / 9 | 0.889 | 0.889 | Stanford clinical SBT (gold) |
| GIAB Ashkenazi trio | germline | WES | 3 / 9 | 1.000 | 1.000 | Stanford clinical SBT (gold) |
| GIAB HG002 | germline | WGS | 1 / 3 | 1.000 | 1.000 | Chin et al. (2020) clinical SBT |
| IHWG MHC-reference cells | germline LCL | WGS | 4 / 12 | 0.833 | 0.917 | IPD-IMGT/HLA IHIW multi-lab |

Two points temper interpretation. First, at these sample sizes the confidence intervals are wide and the ChampHLA–majority-voting differences are not individually significant; the value of the external test lies in the *direction and mechanism* of the effect rather than the absolute gap. Second, NCI-60 are tumour lines, so loss of heterozygosity inflates apparent homozygosity at some loci—indeed the gate's corrective HLA-A overrides recover exactly such constitutional homozygous genotypes, which OptiType had mis-split. Taken together, the external validation supports the claim that the Champion-Challenger framework generalises out of distribution—matching majority voting on both WES and RNA-seq (within overlapping confidence intervals) with an active, mostly-corrective override gate—and that where it trails majority voting the cause is a *frozen* modality-specific choice rather than a failure of the consensus design: on WES a too-strict supporting-tools floor (recoverable, below), and on RNA a frozen HLA-B/-C champion (ArcasHLA) that is strong on 1000G but weak on NCI-60. Both are addressable by transparent, modality-specific recalibration rather than by abandoning the design. This reading is reinforced across the full set of six external arms and all three modalities (Table 9): in every arm ChampHLA matches majority voting within overlapping confidence intervals, and the two arms where it trails—NCI-60 RNA-seq and the IHWG WGS reference cells—share the identical mechanism, a single frozen modality-specific *champion* (ArcasHLA and T1K respectively) that is strong on 1000G but weak out of distribution while the remaining champions perform well. Conversely, on the cleanest arm—the germline, LOH-free, clinical-gold-truth GIAB trio—ChampHLA is perfect on WES and near-perfect on RNA-seq, and it recovers the HG002 whole-genome anchor exactly. The consistent conclusion is that the Champion-Challenger consensus design transfers out of distribution, with any residual shortfall traceable to a frozen champion choice rather than to the consensus logic, and recoverable by modality-specific re-selection.

### 9. Orthogonal Flow-Cytometry HLA-A2 Validation in a Clinical AML Cohort

To validate ChampHLA against a ground truth fully orthogonal to sequencing—surface protein rather than DNA/RNA—we compared computationally derived HLA-A2 serotype calls with flow-cytometry HLA-A2 staining (clone BB7.2) in the FIMM AML cohort. Flow results were available for 84 bone-marrow samples (49 HLA-A2 positive, 35 negative), of which 16 donors also had ChampHLA typing; the remainder lacked matched sequencing and were out of scope. Though small, the set was class-balanced (8 positive, 8 negative donors).

Under the strict `A*02` definition, the multi-tool consensus HLA-A2 call agreed with flow cytometry in every informative donor: sensitivity 8/8 = 1.00 (95% CI 0.68–1.00), specificity 7/7 = 1.00 (0.65–1.00), and accuracy 1.00 (one uninformative tie excluded) (**Table 10**; **Figure 7A**). Concordance held across A2 subtypes, including a donor typed `A*02:06`—an A2 antigen outside the common `A*02:01`—correctly predicted positive. Per method, the RNA-based tools and WES ArcasHLA were essentially concordant with flow, whereas WES OptiType had the lowest specificity (0.40; 3 of 5 negative donors mis-called positive), systematically over-calling `A*02` (**Figure 7B**). Extending the definition to the cross-reactive `A*68`/`A*69` group lowered consensus specificity to 0.86 (one `A*68:01` donor re-classified positive), indicating that for BB7.2 the strict `A*02` mapping is the appropriate genotype-to-serotype rule in this cohort.

Beyond accuracy, the orthogonal assay exposed two data-integrity problems that sequence-only cross-tool comparison left unresolved. In one donor the consensus HLA-A and HLA-B genotypes were irreconcilable across sequencing samples (homozygous `A*02:01` in the WES and most replicates versus `A*03:01/A*68:01` in the bulk-RNA sample, sharing no allele group); the flow HLA-A2–positive result identified the `A*02:01` genotype as correct and flagged the discordant RNA sample as a probable swap. In a second donor a single WES-OptiType `A*02:01` call was contradicted by every other tool and modality (consensus `A*03:01/A*25:01`) and by a flow-negative result, isolating it as a tool-level artefact. These cases illustrate a practical benefit of an orthogonal protein-level assay—it adjudicates both typing accuracy and sample provenance. Given the small sample size (wide confidence intervals), we present this as a proof-of-concept external validation rather than a powered accuracy estimate.

**Table 10. Flow-cytometry HLA-A2 concordance of ChampHLA consensus and individual methods (FIMM AML cohort, strict `A*02` definition).** *Sensitivity (HLA-A2–positive donors), specificity (negative donors) and accuracy with Wilson score 95% confidence intervals; n = informative donors with a call for that method. Consensus is the per-donor majority across available modality × tool calls.*

| Method | n | Sensitivity (95% CI) | Specificity (95% CI) | Accuracy |
|---|---|---|---|---|
| **ChampHLA consensus** | 15 | 1.00 (0.68–1.00) | 1.00 (0.65–1.00) | 1.00 |
| scRNA ArcasHLA | 8 | 0.80 (0.38–0.96) | 1.00 (0.44–1.00) | 0.88 |
| scRNA OptiType | 5 | 1.00 (0.44–1.00) | 1.00 (0.34–1.00) | 1.00 |
| bulkRNA ArcasHLA | 15 | 1.00 (0.68–1.00) | 1.00 (0.65–1.00) | 1.00 |
| bulkRNA OptiType | 8 | 0.80 (0.38–0.96) | 1.00 (0.44–1.00) | 0.88 |
| bulkRNA SpecHLA | 8 | 0.50 (0.15–0.85) | 1.00 (0.51–1.00) | 0.75 |
| WES ArcasHLA | 12 | 1.00 (0.61–1.00) | 1.00 (0.61–1.00) | 1.00 |
| WES OptiType | 10 | 1.00 (0.57–1.00) | 0.40 (0.12–0.77) | 0.70 |
| WES SpecHLA | 10 | 0.75 (0.30–0.95) | 1.00 (0.61–1.00) | 0.90 |

**Figure 7.** *Orthogonal flow-cytometry HLA-A2 validation (FIMM AML cohort).* **(A)** Confusion matrix of the ChampHLA multi-tool consensus HLA-A2 call (strict `A*02`) versus flow-cytometry BB7.2 staining across 15 informative donors (8 positive, 7 negative). **(B)** Sensitivity and specificity by method (consensus and each modality × tool) with Wilson score 95% confidence intervals; WES OptiType shows markedly reduced specificity, reflecting systematic over-calling of `A*02`. *(Source: `figures_final/figure_14_flow_hla_a2`)*

---

## Recommendations

### Practical Guidance for ChampHLA End Users

The following section describes four concrete usage scenarios covering the most common deployment contexts. Fully annotated command-line examples with configuration templates are provided in the Supplementary Section and the repository README at https://github.com/uoozcan/ChampHLA.

#### Scenario 1: Running ChampHLA on an HPC Cluster (SLURM)

For users running ChampHLA on a SLURM-managed HPC cluster, the recommended approach uses Singularity for container execution. Before the first run, pre-pull container images to a shared cache directory to avoid repeated downloads:

```bash
# Set Singularity image cache path in your params file
singularity_cache_dir = "/path/to/shared/singularity_cache"
```

A complete WES run with Champion-Challenger ensemble on SLURM:

```bash
nextflow run champhla/main.nf \
  --input_samplesheet /path/to/samples.csv \
  --input_type bam \
  --seq_type dna \
  --tools optitype,hlahd,t1k,arcashla,spechla,kourami,polysolver \
  --extract_hla_region \
  --enable_majority_voting \
  --weighting calibrated \
  --mv_genes A,B,C \
  --mv_resolution 2 \
  --outdir /path/to/results \
  --max_cpus 16 \
  --max_memory 64.GB \
  --max_time 24.h \
  -profile slurm,singularity \
  -resume
```

Resource requirements per tool module are declared as `process_medium` (8 CPUs, 32 GB, 8 h) by default, with automatic retry at higher resources if a job fails. Runtime scales approximately linearly with cohort size; see **Supplementary Figure S14** for per-tool resource profiles.

#### Scenario 2: Running ChampHLA on a Local Computer (Docker)

For local execution, ensure Docker is installed and that at least 32 GB RAM and 8 CPUs are available. The `--extract_hla_region` flag is strongly recommended for WES BAM inputs to reduce intermediate file sizes and peak memory.

```bash
nextflow run champhla/main.nf \
  --input /path/to/fastq_files/ \
  --input_type fastq \
  --seq_type dna \
  --tools optitype,hlahd,t1k,arcashla \
  --enable_majority_voting \
  --weighting calibrated \
  --outdir ./results \
  -profile docker
```

If RAM is constrained (< 32 GB), consider running a subset of tools optimised for memory efficiency. OptiType, T1K, and ArcasHLA collectively provide strong ensemble performance at lower memory requirements than the full eight-tool panel; see **Supplementary Figure S14** (Results §2) for per-tool resource profiles.

#### Scenario 3: RNA-seq Input for Transcriptomic HLA Typing

For bulk RNA-seq samples, use `--seq_type rna` and `--optitype_seq_type rna` to activate the RNA-specific tool modes. POLYSOLVER and Kourami are automatically excluded for RNA-seq inputs. HLA-HD (the strongest single tool in the RNA-seq benchmark; Table 3) is strongly recommended for RNA-seq runs. ArcasHLA is recommended in addition when a native probabilistic posterior is needed for downstream interpretation.

```bash
nextflow run champhla/main.nf \
  --input_samplesheet rna_samples.csv \
  --input_type bam \
  --seq_type rna \
  --tools optitype,arcashla,hlahd,t1k \
  --optitype_seq_type rna \
  --enable_majority_voting \
  --weighting calibrated \
  -profile docker
```

Be aware that RNA-seq HLA typing is susceptible to allele-specific expression bias at heterozygous loci (see Results §2). The `discordance_summary.tsv` output can identify loci with possible_expression_bias events, which should be flagged for manual review in clinical contexts.

#### Scenario 4: Clinical WES Typing with High-Confidence Output

For clinical or research settings where reliable allele calls are required and false positives (harmful Champion-Challenger overrides) are particularly consequential, we recommend a more conservative override policy—raising `min_challenger_support_fraction` and `min_challenger_margin` above the nested-CV-selected default and requiring at least two supporting tools. The full in-sample sweep surface (`wes_champion_override_sweep.tsv`) is provided as a tuning aid for choosing an institution-specific operating point (it is a tuning aid, not a performance estimate). This trades some corrective-override recall for increased precision:

```bash
nextflow run champhla/main.nf \
  ... \
  --cc_min_support_fraction 0.50 \
  --cc_min_margin 0.05 \
  --cc_min_tools 2 \
  --enable_majority_voting \
  --weighting calibrated \
  -profile docker
```

The full sweep TSV (`wes_champion_override_sweep.tsv`) is included in the benchmark outputs and allows users to inspect the complete accuracy-override tradeoff surface for their own cohort, enabling institution-specific policy tuning without modifying the core algorithm.

### Limitations and Future Directions

The most important current limitation of ChampHLA is the low absolute accuracy of short-read WGS HLA typing: even after Champion-Challenger routing lifts the overall correct-call rate to 0.5012 (Table 3), roughly half of WGS loci remain miscalled, and HLA-C stays hard (0.3285) for every tool and consensus method. Champion routing significantly improves over majority voting on WGS (Results §1) but cannot exceed the accuracy of the best available per-gene tool; the priority for future WGS improvement therefore lies in better individual caller quality, particularly at HLA-A and HLA-C. Expanding to higher-resolution truth and incorporating long-read WGS data (where T1K already provides a pathway) are the most impactful near-term improvements. A substantial fraction of the apparent WGS error is also an artefact of the two-field 2014 truth set (Results §5): under ambiguity-compatible scoring, WGS accuracy rises to 0.79–0.88.

Seven of the eight integrated tools are distributed as pre-built Docker container images that are automatically pulled at first execution and can be converted to Singularity images for HPC environments without root access. SpecHLA (v1.1) requires a more complex dependency chain—including conda-managed bioinformatics libraries and internal reference databases—that historically precluded containerisation. The pipeline now provides a purpose-built Dockerfile (`containers/spechla/Dockerfile`) that packages the full SpecHLA installation into a single container image, as well as a local installation mode (`--use_local_spechla true --spechla_path /path/to/SpecHLA`) for environments where containerisation is impractical, such as HPC systems where SpecHLA's internal Singularity calls conflict with the host Singularity runtime. Users select which tools to include via the `--tools` parameter, and only the corresponding containers and databases need to be installed; SpecHLA can be omitted entirely without affecting the functionality of the remaining seven-tool panel. Seq2HLA and SpecHLA experienced container configuration failures in the RNA-seq benchmark, reducing their callable rates to 29/107 and 80/107 respectively; these failures reflect pipeline engineering gaps rather than inherent tool limitations and will be resolved in future benchmark cycles. A practical consequence is that the effective RNA-seq panel is smaller than eight tools: Seq2HLA (29/107 callable) and SpecHLA (76.2% callable at 0.035 accuracy) contribute little to RNA-seq typing, and the RNA-seq Champion-Challenger accuracy is reported from a SpecHLA-inclusive scoring run to benchmark against the largest possible panel—default production RNA-seq typing relies on the OptiType/HLA-HD/ArcasHLA/T1K core. The full eight-tool panel is realised on DNA (WGS/WES) inputs.

The ground-truth HLA labels used in all benchmarks are derived from the 2014 1000 Genomes typing study, which provides only two-field resolution for most loci. This limits resolution analysis and may systematically penalise high-resolution callers that correctly identify alleles post-dating the 2014 reference. Future benchmark cycles should incorporate complementary truth sources with higher-resolution labelling, including orthogonal clinical typing data and newer high-resolution reference panels. The orthogonal silver-standard truth capability introduced here (Methods §8; Results §6) is a first step in this direction—it reproduces gold calls at G-group resolution, which the 2014 panel cannot encode—but its accuracy is currently bounded by the 45-sample HPRC v1.1 panel's allele coverage and, for region-restricted input BAMs, by shared rather than per-sample depth estimation; the larger HPRC v2 panel (200+ samples) and per-sample preprocessing are expected to lift both bounds.

Benchmark-derived reliability weights are currently estimated as global averages across mixed-ancestry cohorts. Given that HLA allele frequencies differ substantially across ancestry groups, and that tools validated primarily on European reference panels may perform differently on non-European samples, per-population weight learning is a high-priority extension for clinical deployment contexts. The Champion-Challenger method currently uses a single champion per gene applied uniformly across all samples; a per-population champion designation—e.g., using HLA-HD as the WGS champion for CHB samples if it outperforms OptiType in that subgroup—could improve accuracy for ancestrally stratified cohorts.

Several of the confirmatory analyses are deliberately small and are presented as proof-of-concept rather than powered estimates: the orthogonal flow-cytometry validation covers 15 informative donors (with two data-integrity catches), the silver-standard-truth concordance is powered at 30 samples (with a three-sample per-depth pilot), and the germline external cohorts (GIAB trio, HG002, IHWG reference cells) comprise three to four samples each. Their value lies in orthogonality and mechanism—protein-level truth, assembly-based truth, and out-of-distribution transfer of the frozen policy—rather than in narrow confidence intervals; larger replications are planned.

The current headline evaluation is restricted to HLA-A, -B, and -C. DRB1 and DQB1 are implemented in the framework and are available for consensus voting via `--mv_genes A,B,C,DQA1,DQB1,DRB1`, but truth-backed benchmarking at these loci awaits complementary truth resources for the 1000 Genomes cohorts used in this study. Class II typing (clinically important for transplantation matching) is therefore left to future work as a priority for the next benchmark cycle. Extension to five-locus (A, B, C, DRB1, DQB1) truth-backed evaluation would substantially increase the clinical relevance of the benchmark. Tighter integration of long-read and pangenome-graph backends into the consensus panel (the silver-truth pipeline already employs a pangenome-graph genotyper), per-population weight stratification, allele-frequency priors in the override policy, and integration with downstream immunoinformatics tools (neoantigen prediction pipelines, HLA ligandome databases, transplant compatibility scoring) are all planned for the ChampHLA version 2.1 roadmap. Community contributions via the GitHub repository are welcome and will shape the feature priority of the next release.

A concise summary of limitations and future directions is provided above; full technical detail, extended scenario walkthroughs, and the version 2 roadmap are available in the Supplementary Section.

---

## Conclusion

ChampHLA addresses one of the most persistent operational challenges in computational HLA typing: the fragmentation of tool ecosystems, the heterogeneity of individual tool outputs, and the absence of a principled, reproducible framework for resolving disagreements between tools. By integrating eight state-of-the-art HLA typing tools into a single Nextflow DSL2 pipeline with containerised execution, ChampHLA eliminates per-tool installation burden and routes inputs automatically to each tool's required format, removing the engineering barriers that currently prevent many research groups from routinely applying ensemble HLA typing at scale.

The core methodological contribution—the Champion-Challenger ensemble algorithm—was motivated by the observation that different tools predict different alleles for the same sample, and that current ensemble methods either treat all tools equally (majority voting) or trust confidence signals that are not empirically calibrated (naive weighted voting). Champion-Challenger addresses both limitations: it designates a benchmark-validated champion for each gene and allows the ensemble to override the champion only when the aggregate evidence is simultaneously strong in weight, broad in tool support, and unambiguous in allele identity. Under nested cross-validation with strictly out-of-fold operating points, this approach produced its largest—and only statistically significant—accuracy gain in the hardest, most-discordant modality: in WGS it reached 0.5012, significantly above majority voting (0.3844; +11.7 percentage points; exact McNemar p<0.0001) and level with the best single tool (OptiType 0.4975; p=0.74), through per-gene champion routing. In the high-accuracy WES and RNA-seq modalities, where the strong tools already agree, it matched majority voting within overlapping confidence intervals (0.9436 and 0.9408; differences not significant). An ablation confirmed that the accuracy contribution comes from champion routing rather than the weighting scheme, whose distinct value is a calibration guardrail that blocks miscalibrated tool confidence from corrupting ensemble weights. The evidence-gated override itself rarely changed the outcome—it triggered no overrides on WGS and did not significantly alter overall accuracy in any modality—and is therefore best understood not as an accuracy booster but as an auditable, predominantly corrective safety mechanism that retains the champion's call unless the weighted ensemble evidence is simultaneously broad, marginal, and unambiguous. Likewise, the calibration guardrail is a soundness contribution—it prevents an uninformative confidence proxy from distorting the weights—rather than a driver of the reported accuracy.

The usability contribution is equally important: a practitioner who can run Nextflow and has Docker or Singularity installed can obtain harmonised HLA calls from the full eight-tool panel with a single command, without reading eight separate tool documentation pages, managing eight separate reference databases, or writing custom output-parsing scripts. The harmonised output schema, structured discordance labels, and calibration-guardrail audit trail provide interpretable quality metrics at every stage of the pipeline, from individual tool outputs through ensemble decisions to final consensus calls.

We anticipate that ChampHLA will be most immediately useful to research groups conducting population-scale HLA typing for pharmacogenomic biomarker discovery, transplantation cohort studies, and cancer immunogenomics projects where both accuracy and reproducibility at scale are essential. The matched-subject robustness benchmark identifies WES+RNA bimodal consensus as the strongest tested multimodal configuration under current short-read technology; WGS addition did not improve accuracy in this cohort, suggesting that WGS single-tool quality improvement is the more impactful near-term priority for WGS-equipped workflows. These findings also delineate when each part of ChampHLA matters. Where paired WES and RNA-seq are available, bimodal consensus is the recommended configuration; there Champion-Challenger matches bimodal majority voting exactly (0.9623 vs 0.9623; McNemar p=1.0) while beating the best single tool (0.9423), so it can serve as the default consensus—delivering calibration-guarded, auditable calls at no accuracy cost—rather than being reserved for difficult modalities. Where only WGS is available—as in many biobank, archival, or DNA-only cohorts—champion routing provides a significant and fully auditable accuracy gain over majority voting, recovering the best-single-tool accuracy ceiling that symmetric voting dilutes (Results §1). The contribution of ChampHLA is therefore the combination of a reproducible multi-tool substrate, a calibration-governed and auditable consensus that does no harm on concordant data and helps on discordant data, and an orthogonal silver-standard-truth capability for cohorts without gold-standard typing—rather than a uniformly higher accuracy than majority voting across all modalities.

---

## Supplementary Materials

All supplementary items are provided in a separate document, **ChampHLA — Supplementary Materials** (`docs/CHAMPHLA_MANUSCRIPT_V3_SUPPLEMENTARY.md` / `.docx`): Supplementary Figures S1–S8 and S12–S17, Supplementary Tables S9–S11, and Supplementary Notes S1–S11. Each supplementary figure is embedded there with its source stem under `analysis/figures_final/`.

## Acknowledgements

The authors thank the 1000 Genomes Project Consortium for generating and sharing the HLA truth-backed dataset used as the primary benchmark in this study (Gourraud et al., 2014). Computational analyses were performed on the CSC Puhti supercomputer (CSC — IT Center for Science, Finland). Funding information will be provided prior to submission. The authors declare no competing interests.

## Data Availability

The ChampHLA Nextflow pipeline (version 2.0.0) is openly available at https://github.com/uoozcan/ChampHLA. The nested cross-validation re-analysis and all derived tables—`analysis/nested_cv_champion_challenger/` (per-modality method comparison, per-gene accuracy, per-gene McNemar tests, and the CIWD- and ancestry-stratified TSVs), the `analysis/benchmark_{wgs,wes,rna}_cv_recalibrated/` harmonised call tables, calibration artifacts, confidence weight files, and the (in-sample) Champion-Challenger sweep surfaces—together with the regenerated main and supplementary figures, will be deposited to Zenodo at the time of submission and assigned a permanent DOI (cited here as [Zenodo DOI, to be inserted]). The primary and stratified results are regenerable from the deposited inputs with `bin/reproduce_manuscript_numbers.sh` (which runs `bin/nested_cv_champion_challenger.py` for all three modalities, `bin/populate_cc_stubs_from_nested_cv.py`, and the two figure generators); an HPC driver is provided as `slurm_nested_cv.sh`, and a caption-consistent, journal-ready figure set is assembled by `bin/build_submission_figures.py`. The bimodal WES+RNA Champion-Challenger result (Table 4) is produced by `bin/nested_cv_champion_challenger.py --bimodal`. All benchmark runs are anchored to IMGT/HLA database version 3.59.0 (pinned in `benchmark_metadata.json` for each benchmark root). The 1000 Genomes HLA truth panel used for all benchmarks is available at the 1000 Genomes Project FTP archive (Gourraud et al., 2014). The orthogonal silver-standard truth pipeline (`bin/generate_silver_truth.py`, `bin/validate_silver_truth.py`, and the `modules/locityper.nf`, `modules/immuannot.nf`, and `modules/silver_truth.nf` Nextflow modules) is distributed with the pipeline; the Human Pangenome Reference Consortium v1.1 Minigraph-Cactus pangenome used to build the Locityper database is publicly available from the HPRC open-data archive (Liao et al., 2023), and the silver-truth genotyping outputs and concordance summaries will be deposited to Zenodo upon acceptance. The independent NCI-60 external-validation cohort uses openly available data: constitutional HLA truth from Adams et al. (2005, PMC555742), NCI-60 RNA-seq from SRA study SRP133178 (BioProject PRJNA433861), and NCI-60 whole-exome capture from SRA study SRP150855; the derived `truth_long.tsv`, frozen-policy benchmark tables, and the WES Champion-Challenger override sweep will be deposited to Zenodo upon acceptance. The additional germline external-validation cohorts also use openly available data: the GIAB Ashkenazi trio (HG002/HG003/HG004) RNA-seq runs SRR15909917–SRR15909919 and whole-exome runs SRR2962669/SRR2962692/SRR2962694, scored against Stanford Blood Center clinical sequence-based typing (Llamas et al., 2019, F1000Research 8:1751; Chin et al., 2020); the GIAB HG002 whole-genome alignment (`HG002.GRCh38.2x250`); and the IHWG MHC-reference-cell whole-genome sequencing from BioProject PRJNA764575, scored against the IPD-IMGT/HLA IHIW multi-laboratory reference-cell consensus. The derived `truth_long.tsv` files and frozen-policy benchmark tables for these cohorts will likewise be deposited to Zenodo upon acceptance. Flow-cytometry HLA-A2 results and the derived concordance tables for the FIMM AML cohort are provided in `pihla-publish/fimm_results/hla_a2_flow/` (source spreadsheet `AB_HLA_A2_typing.xlsx`); the analysis scripts are `bin/validate_hla_a2_flow.py`, `bin/champhla_detail_flow_cohort.py`, and `bin/plot_hla_a2_flow.py`.

---

## References

Adams, S.D. et al. (2005). Ambiguous allele combinations in HLA Class I and II sequence-based typing: when precise nucleotide sequencing leads to imprecise allele identification. *Journal of Translational Medicine*, 3, 30.

Boegel, S. et al. (2012). HLA typing from RNA-seq sequence reads. *Genome Medicine*, 4, 102.

Chin, C.S. et al. (2020). A diploid assembly-based benchmark for variants in the major histocompatibility complex. *Nature Communications*, 11, 4794.

Chung, W.H. et al. (2004). Medical genetics: a marker for Stevens-Johnson syndrome. *Nature*, 428, 486.

Di Tommaso, P. et al. (2017). Nextflow enables reproducible computational workflows. *Nature Biotechnology*, 35, 316–319.

Gourraud, P.A. et al. (2014). HLA diversity in the 1000 genomes dataset. *PLoS ONE*, 9, e97282.

Hickey, G. et al. (2024). Pangenome graph construction from genome alignments with Minigraph-Cactus. *Nature Biotechnology*, 42, 663–673.

Lee, H. & Kingsford, C. (2018). Kourami: graph-guided assembly for novel human leukocyte antigen allele discovery. *Genome Biology*, 19, 16.

Liao, W.W. et al. (2023). A draft human pangenome reference. *Nature*, 617, 312–324.

Llamas, B. et al. (2019). A strategy for building and using a human reference pangenome. *F1000Research*, 8, 1751.

Mallal, S. et al. (2008). HLA-B*5701 screening for hypersensitivity to abacavir. *New England Journal of Medicine*, 358, 568–579.

Matzaraki, V. et al. (2017). The MHC locus and genetic susceptibility to autoimmune and infectious diseases. *Genome Biology*, 18, 76.

Mishima, H. et al. (2019). A practical guide for the evaluation of next-generation sequencing-based HLA typing. *Human Immunology*, 80, 891–898.

Orenbuch, R. et al. (2019). arcasHLA: high-resolution HLA typing from RNAseq. *Bioinformatics*, 36, 33–40.

Parham, P. & Brodsky, F.M. (1981). Partial purification and some properties of BB7.2, a cytotoxic monoclonal antibody with specificity for HLA-A2 and a variant of HLA-A28. *Human Immunology*, 3, 277–299.

Petersdorf, E.W. et al. (2007). Major histocompatibility complex variation and outcomes in hematopoietic cell transplantation. *Immunological Reviews*, 214, 229–244.

Prodanov, T. et al. (2025). Locityper enables targeted genotyping of complex polymorphic genes. *Nature Genetics*, 57, 2901–2908.

Robinson, J. et al. (2020). IPD-IMGT/HLA database. *Nucleic Acids Research*, 48, D948–D955.

Schumacher, T.N. & Schreiber, R.D. (2015). Neoantigens in cancer immunotherapy. *Science*, 348, 69–74.

Shukla, S.A. et al. (2015). Comprehensive analysis of cancer-associated somatic mutations in class I HLA genes. *Nature Biotechnology*, 33, 1152–1158.

Song, L. et al. (2023). Efficient and accurate KIR and HLA genotyping with massively parallel sequencing data. *Genome Research*, 33, 923–931.

Szolek, A. et al. (2014). OptiType: precision HLA typing from next-generation sequencing data. *Bioinformatics*, 30, 3310–3316.

Zhou, Y. et al. (2024). Full-resolution HLA and KIR gene annotations for human genome assemblies. *Genome Research*, 34, 1931–1941.
