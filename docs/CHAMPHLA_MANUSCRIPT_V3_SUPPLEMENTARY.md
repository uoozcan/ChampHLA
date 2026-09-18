# ChampHLA — Supplementary Materials

**Supplementary material for:** *ChampHLA: a reproducible Nextflow ensemble for champion–challenger HLA typing across WGS, WES, and RNA sequencing.* Özcan et al.

This document contains all supplementary figures, tables, and notes referenced by the main manuscript. Supplementary Figures are numbered S1–S14, Supplementary Tables S1–S3, and Supplementary Notes S1–S8, each series numbered independently. All figures are regenerable from the deposited benchmark tables via the scripts noted in Data Availability; each caption records its source figure stem under `analysis/figures_final/`.

---

## Supplementary Figures

**Supplementary Figure S1.** *Extended per-gene, per-tool accuracy across sequencing modalities.* Per-gene (HLA-A, HLA-B, HLA-C) overall correct-call rate for every individual tool, majority voting, and Champion-Challenger in the WGS, WES, and RNA-seq nested cross-validation benchmarks. Extends the summary per-gene view of main-text Figure 3 to the full tool panel, showing the wide spread of single-tool accuracy—particularly in WGS—that motivates per-gene champion routing. *(Source: `figures_final/figure_s1_per_gene_accuracy`)*

**Supplementary Figure S2.** *Per-tool confidence calibration heatmap.* Brier score and expected calibration error (ECE) for each tool's native confidence proxy across modalities, with the Brier/ECE < 0.35 guardrail threshold marked. Tools whose proxy fails the guardrail (OptiType, HLA-HD, Kourami) are blocked from receiving a confidence boost; only T1K and ArcasHLA pass. Complements main-text Figure 4. *(Source: `figures_final/figure_s2_calibration_heatmap`)*

**Supplementary Figure S3.** *Multi-resolution match-rate comparison by tool.* Exact two-field, exact three-field, G-group, P-group, and ambiguity-compatible match rates by tool for the WGS and WES truth-backed benchmarks, illustrating the two-field ceiling of the 2014 1000 Genomes truth panel (near-zero three-field/G-group/P-group rates for all tools except OptiType). Supports main-text Table 7 and Results §5. *(Source: `figures_final/figure_s3_resolution_comparison`)*

**Supplementary Figure S4.** *Cross-modality HLA concordance in the FIMM AML/MDS cohort.* Per-gene (HLA-A, -B, -C) allele-call concordance between all pairs of sequencing modalities (single-cell RNA-seq, bulk RNA-seq, WES) applied to a clinical AML/MDS cohort. Concordance is highest between bulk RNA and WES, consistent with the WES+RNA complementarity observed in the 1000 Genomes benchmark; scRNA-vs-WES concordance is lower, particularly at HLA-C, reflecting reduced HLA-region read depth in droplet scRNA-seq. See Supplementary Note S7.1. *(Source: `figures_final/figure_s4_fimm_concordance`)*

**Supplementary Figure S5.** *HLA allele dropout and homozygosity QC by tool and modality.* False duplicate rate (fraction of heterozygous truth loci called homozygous) per tool and modality from the 1000 Genomes truth-backed benchmarks. ArcasHLA shows the highest false-duplicate rate in RNA-seq (up to 1.0 at HLA-A/-B) under allele-specific expression, whereas WES shows near-zero rates for most tools—supporting WES-based calls as the primary input for downstream LOH and HED analyses. See Supplementary Note S7.2. *(Source: `figures_final/figure_s5_allele_dropout_heatmap`)*

**Supplementary Figure S6.** *FIMM WES loss-of-heterozygosity candidate classification.* Classification of each patient–gene pair into candidate LOH (major-allele fraction ≥ 0.80 with ≥ 20 heterozygous variant sites), allelic imbalance requiring review, balanced heterozygous, or insufficient evidence, from SpecHLA WES allele-frequency outputs. These classifications are exploratory: WES exon-capture allele frequencies may not transfer thresholds calibrated for WGS, and WGS confirmation is required before clinical interpretation. See Supplementary Note S7.3. *(Source: `figures_final/figure_s6_fimm_loh_wes`)*

**Supplementary Figure S7.** *Overall survival stratified by HLA evolutionary divergence (FIMM cohort).* Kaplan–Meier overall-survival curves for 28 FIMM patients with survival follow-up, stratified by diagnosis group (AML, MDS, MDS→AML) and by median total HLA evolutionary divergence (HED) computed from WES consensus calls at HLA-A/-B/-C. The analysis is exploratory (n ≈ 14 per arm; ~25% power for HR = 2.0) and hypothesis-generating. See Supplementary Note S7.4. *(Source: `figures_final/figure_s7_fimm_survival_hed`)*

**Supplementary Figure S8.** *Abstention–accuracy tradeoff curves across modalities.* Line plots show overall correct-call rate (y-axis) against callable rate (x-axis) as the WeightedConsensus support threshold varies from 0 (always call) to 1 (always abstain), for WGS, WES, and RNA-seq benchmarks. Each point represents one threshold setting; the operating point selected out of fold under nested cross-validation is marked. This abstention behaviour is a secondary reliability property of the weighting-only consensus and is not part of the primary Champion-Challenger accuracy comparison. The WES and RNA-seq curves show a clear tradeoff: accepting a modest callable-rate reduction (from 1.0 to 0.93 in WES, 0.96 in RNA-seq) yields accuracy-among-callable estimates above MajorityVote accuracy. *(Source: `figures_final/figure_5_abstention_tradeoff`)*

**Supplementary Figure S9.** *Discordance taxonomy counts by modality.* Bar charts show the number of discordance events per category (`low_evidence_conflict`, `no_evidence`, `possible_expression_bias`, `technical_conflict`) for the WGS (n=411 loci), WES (n=390), and RNA-seq (n=321) benchmarks. WGS carries by far the highest discordance (42.3%, dominated by `low_evidence_conflict` from the wide spread of per-tool WGS accuracy), WES is intermediate (6.9%), and RNA-seq the lowest (4.0%), where the strong tools largely agree. The discordance ordering matches the modality ordering in which Champion-Challenger improves over majority voting. *(Source: `figures_final/figure_6_discordance_taxonomy`)*

**Supplementary Figure S10.** *Benchmark-derived reliability weights by tool and modality.* Horizontal bar charts show the final blended weight (`final_weight = 0.7 × base_reliability + 0.3 × effective_confidence`) for each tool in the WGS, WES, and RNA-seq modalities. Tool-level weights are shown separately for HLA-A, HLA-B, and HLA-C where gene-level variation exists. Guardrail status is encoded by bar fill (solid = confidence boost applied; hatched = blocked). OptiType dominates WGS weighting; HLA-HD and OptiType lead in RNA-seq; POLYSOLVER and OptiType lead in WES. ArcasHLA's near-zero WGS weight and near-average RNA-seq weight illustrate why modality-specific weight learning is essential. *(Source: `figures_final/figure_7_confidence_weights`)*

**Supplementary Figure S11.** *Computational resource requirements per HLA typing tool across sequencing modalities.* Three-panel figure showing peak RAM (GB), CPU utilisation (%), and wall-clock time (hours). Within each panel, grouped bars show median resource values per tool for WGS (orange), WES (blue), and RNA-seq (green) from 1000 Genomes Project samples (n=30 per modality for WGS and WES; RNA-seq: T1K and Seq2HLA n=30, ArcasHLA/HLA-HD/SpecHLA n=29, OptiType n=26 due to four samples exceeding the 8 h wall-clock limit). Error whiskers indicate P25–P75 interquartile range. Tools not applicable to a given modality are shown with hatched bars. All values measured from Nextflow execution traces (`realtime`, `peak_rss`, `%cpu` fields) from dedicated n=30 benchmark runs. *(Source: `figures_final/figure_8_computational_performance`)*

**Supplementary Figure S12.** *Trimodal WGS+WES+RNA robustness analysis (n=106 matched subjects; supplementary robustness benchmark).* Two-panel comparison of trimodal and bimodal ensemble methods. Panel A: accuracy among callable loci for bimodal MajorityVote, bimodal WeightedConsensus, trimodal MajorityVote, and trimodal WeightedConsensus, with exact percentage annotations and Wilson score 95% confidence intervals. Panel B: callable rate (fraction of loci with a call) for the same four methods, showing the coverage cost of adding WGS. Hatched bars indicate trimodal methods. Adding WGS as a third modality does not improve over bimodal WES+RNA and reduces overall accuracy slightly, confirming that short-read WGS adds noise rather than signal to an already-strong WES+RNA ensemble. This finding is specific to the current benchmark cohort and tool set; it does not constitute a general recommendation to exclude WGS from all HLA typing workflows, but identifies WGS single-tool quality improvement as the most impactful open problem for WGS-equipped HLA typing pipelines. *(Source: `figures_final/figure_10_trimodal_comparison`)*

**Supplementary Figure S13.** *Orthogonal silver-standard truth validation against 1000 Genomes gold truth.* (A) Allele-level concordance per locus (HLA-A, -B, -C) and overall, at two-field and G-group resolution, for a Locityper genotyping database built from the HPRC v1.1 Minigraph-Cactus pangenome versus a control database built from IPD-IMGT/HLA alleles, on 30 samples (shared depth profile) and a three-sample per-sample-depth pilot. The HPRC pangenome database nearly doubles concordance (overall two-field 0.811 versus 0.461 at 30 samples; 0.944 versus 0.556 in the pilot). (B) Confidence calibration: concordance as a function of Locityper genotype-quality threshold for the 30-sample run. Concordance increases monotonically with quality under the HPRC database (0.811→0.905) but is uninformative under the IMGT-allele database (0.461→0.432), demonstrating that genotype confidence is well-calibrated only with an adequate pangenome reference and motivating the agreement-gated silver-truth design. *(Source: `figures_final/figure_12_silver_truth_hprc`)*

**Supplementary Figure S14.** *External validation on an independent non-1000G cohort (NCI-60).* The frozen 1000G-derived Champion-Challenger policy is applied without re-learning to NCI-60 RNA-seq and WES, scored against Adams (2005) sequence-based typing (HLA-A/-B/-C, two-field). **(A)** Overall correct-call rate with Wilson 95% confidence intervals for the best single tool, majority voting, and ChampHLA, by modality (RNA-seq *n*=42, WES *n*=42); on WES ChampHLA (0.842) matches majority voting and exceeds OptiType (0.817), and on RNA-seq ChampHLA (0.793) tracks just below majority voting (0.817) with a predominantly corrective gate (15 of 20 overrides corrective). **(B)** Per-locus correct-call rate (HLA-A/-B/-C) for ChampHLA and majority voting in each modality. **(C)** WES override-gate behaviour at the deployed support=0.35/margin=0.00 as the supporting-tools floor varies: the deployed floor of ≥1 tool gives 0.842 (three corrective, one harmful, three neutral overrides, equal to majority voting); raising it to ≥2 removes all harmful/neutral overrides and yields the optimum (0.854, three corrective / zero harmful, above the majority-voting reference 0.842); ≥3 suppresses all overrides (0.817, equal to OptiType). *(Source: `figures_final_candidate/figure_13_external_validation_nci60`)*

---

## Supplementary Tables

**Supplementary Table S1. Accuracy by allele commonness (CIWD 3.0.0).** Held-out nested-CV correct-call rate stratified by the commonness of the rarer truth allele (CIWD 3.0.0; `analysis/nested_cv_champion_challenger/*/nested_cv_by_ciwd.tsv`). The common stratum dominates every cohort; intermediate/well-documented strata carry too few loci (n ≤ 3) for inference and are omitted here. The WGS Champion-Challenger gain is retained within the common stratum, confirming it is not a rare-allele artefact.

| Modality | Stratum | n loci | MajorityVote | ChampHLA (CC) |
|---|---|---|---|---|
| WGS | common | 397 | 0.390 | **0.504** |
| WES | common | 377 | 0.952 | 0.955 |
| RNA-seq | common | 312 | 0.955 | 0.946 |

**Supplementary Table S2. Accuracy by continental ancestry.** Held-out nested-CV correct-call rate stratified by 1000 Genomes superpopulation mapped to continental ancestry (EUR = CEU/FIN/GBR/TSI; AFR = YRI; `analysis/nested_cv_champion_challenger/*/nested_cv_by_ancestry.tsv`). Champion-Challenger is greater than or equal to majority voting in every ancestry stratum and modality; the African-ancestry and unlabelled strata carry small n and wide intervals, so per-population weight learning (Limitations) remains a valuable extension.

| Modality | Ancestry | n loci | MajorityVote | ChampHLA (CC) |
|---|---|---|---|---|
| WGS | AFR | 21 | 0.476 | **0.571** |
| WGS | EUR | 93 | 0.430 | 0.441 |
| WGS | unlabelled | 297 | 0.364 | **0.515** |
| WES | AFR | 21 | 0.905 | **1.000** |
| WES | EUR | 93 | 0.968 | 0.968 |
| WES | unlabelled | 276 | 0.928 | 0.931 |
| RNA-seq | AFR | 21 | 0.952 | 0.905 |
| RNA-seq | EUR | 93 | 0.979 | 0.979 |
| RNA-seq | unlabelled | 207 | 0.937 | 0.928 |

**Supplementary Table S3. Attempted external-validation cohorts excluded on integrity grounds.** Beyond the external cohorts reported in Results §8 (NCI-60, GIAB trio, HG002, IHWG), we scouted five further non-1000G candidates and rejected each after Phase-0 provenance verification, before any benchmarking. We report them here because the exclusion criteria—experimental (not in-silico) per-sample truth, independence from the 1000 Genomes calibration cohort, and matched open short reads—are the same standards that make an external validation meaningful, and because published per-sample HLA tables are frequently in-silico tool output rather than orthogonal truth (verified twice here). Per-cohort Phase-0 records are in `analysis/{swehla,pcrsbt,conshla,getrm}_benchmark/PHASE0_STATUS.md`.

| Candidate cohort (reference) | Material / modality | Reason excluded | Standard upheld |
|---|---|---|---|
| SweHLA / SweGen (Nabais Sá et al., *EJHG* 2020) | germline WGS | Truth is access-gated: the open DOI record (10.17044/NBIS/G000009) is metadata-only with zero downloadable files; per-sample SweHLA genotypes require NBIS/SciLifeLab (swefreq) registration | Truth must be obtainable |
| PCR-SBT 829-WES (Yu et al. 2022, PMC9679531) | germline WES | The 829 WES samples **are** 1000 Genomes Phase 3 → circular with our primary 1000G truth; not an independent cohort | Independence from the calibration cohort |
| PCR-SBT 652-RNA (Mangul Lab, PMC10827116) | bulk / single-cell RNA | Every diploid class-I subset is 1000G/HapMap (Geuvadis, Montgomery); the genuinely independent subsets are class-II-only, artificial mono-allelic (B721.221), 10× scRNA, or unreleased | Independent diploid class-I bulk reads |
| consHLA / ZERO Childhood Cancer (*BMC Bioinformatics* 2025, PMC12363109) | WGS + RNA | Published per-sample "truth" is consHLA's own in-silico HLA-HD calls; the experimental clinical truth (10 patients) is aggregate-only with no per-sample genotypes and no deposited reads | Truth must be experimental and per-sample, not in-silico |
| GeT-RM (Bettinotti et al., *J Mol Diagn* 2018, PMC6939753) | Coriell reference DNA | Integrity-clean three-field PCR-SSO/SBT truth was built for all 108 lines, but only 2/108 have usable public short reads (NA12273 RNA-seq, NA17221 WGS) → n=2, not a benchmark; the truth table is retained for future use | Truth needs matched open short reads |

This discipline—rejecting a self-referential in-silico "truth", catching two cohorts that are covertly the 1000 Genomes resource we calibrate on, and declining an access-gated or read-less panel—is why the external-validation cohorts we do report (Results §8) are independent, experimentally-truthed, and openly reproducible.

---

## Supplementary Notes

### Supplementary Note S1. Full HPC Usage Guide

Detailed SLURM configuration with example `nextflow.config` overrides, Singularity SIF cache setup instructions, CSC Puhti-specific parameter file, resource sizing guide by cohort size and sequencing modality, and annotated example commands for WGS, WES, and RNA-seq cohort runs.

### Supplementary Note S2. Full Local Computer Usage Guide

Docker installation prerequisites, per-tool resource profiles (RAM, CPU, disk), recommended tool subsets for memory-constrained environments (< 16 GB RAM, < 32 GB RAM), Docker Compose alternative for multi-sample parallelisation, and expected runtimes on consumer hardware.

### Supplementary Note S3. Full Scenario Walkthroughs

Step-by-step walkthroughs for all four user scenarios described in the Recommendations section, including: (S3.1) SLURM WES cohort run with full eight-tool panel and Champion-Challenger output; (S3.2) local Docker RNA-seq run with HLA-HD+ArcasHLA+OptiType+T1K panel; (S3.3) samplesheet CSV preparation and validation; (S3.4) custom benchmark YAML configuration for a new cohort with user-provided truth data.

### Supplementary Note S4. Extended Resolution and Ambiguity Analysis

Extended per-locus resolution analysis for the WGS (n=137) and WES truth-backed benchmarks, including per-gene breakdowns of exact_3field, ambiguity_compatible, g_group, and p_group match rates for all eight tools. Supplements the summary findings in Results §5 (main text Table 7 and Supplementary Figure S3). Includes discussion of the 2014 truth-set ambiguity ceiling and its implications for interpreting apparent WGS accuracy gaps.

### Supplementary Note S5. Extended Limitations and ChampHLA v2 Roadmap

Full discussion of each limitation identified in the main text, including: training-split size and its effect on weight uncertainty; proxy-based confidence representation and its limitations relative to calibrated probability models; single-champion-per-gene design and its implications for population-stratified cohorts; truth-set ceiling at two-field resolution; and WGS tool-limitation as the central barrier to WGS ensemble accuracy improvement. Planned features for version 2.1 include: DRB1/DQB1 five-locus evaluation with complementary truth resources; per-population weight stratification; allele-frequency priors in override policy; isotonic regression confidence calibration to replace proxy-based Platt scaling; graph-aware and long-read backends via T1K; and integration hooks for neoantigen prediction and transplant compatibility scoring downstream tools.

### Supplementary Note S6. Champion-Challenger Sweep Surface

Description of the sweep TSV file format, how to read the accuracy-override tradeoff table, and guidance for selecting a custom operating point for institution-specific precision requirements. **The sweep surface is an in-sample tuning aid for choosing a production operating point, not a held-out performance estimate:** its accuracies are computed on the full benchmark cohort without a train/test split, so the accuracy at any swept point is optimistically biased and must not be read as the method's expected accuracy. All accuracy figures reported in the main manuscript instead come from the nested cross-validation path (Methods §2), in which the operating point is selected strictly out of fold; the sweep is provided only so users can trade correct-call rate against override aggressiveness for their own precision requirements.

### Supplementary Note S7. FIMM Clinical HLA Typing and Loss of Heterozygosity Analysis

#### S7.1 Cross-Modality HLA Concordance (FIMM AML/MDS Cohort)

To assess ChampHLA's clinical utility on real-world cancer samples, we applied the pipeline to a cohort of AML and MDS patients from the Institute for Molecular Medicine Finland (FIMM). HLA typing was performed across three sequencing modalities—single-cell RNA-seq (scRNA-seq), bulk RNA-seq, and whole-exome sequencing (WES)—using ChampHLA tools available for each modality (OptiType, ArcasHLA, and SpecHLA for WES and bulk RNA; OptiType and ArcasHLA for scRNA). Cross-modality concordance rates were computed for each HLA gene (A, B, C) by comparing allele calls between all modality pairs (**Supplementary Figure S4**).

Concordance was highest between bulk RNA and WES across all three genes, consistent with the strong complementarity of these modalities observed in the 1000 Genomes benchmark (Results §7). scRNA vs WES concordance was lower, particularly at HLA-C, reflecting the reduced read depth at HLA loci typical of droplet-based scRNA-seq protocols. These cross-modality patterns are consistent with the benchmark finding that bimodal WES+RNA consensus produces the highest ensemble accuracy.

#### S7.2 HLA Allele Dropout and Homozygosity QC

Before interpreting allele-level results, we characterised allele dropout patterns across tools and modalities using the 1000 Genomes truth-backed benchmarks (**Supplementary Figure S5**). The false duplicate rate (fraction of heterozygous truth loci called as homozygous) varied substantially by tool and modality. ArcasHLA exhibited the highest false duplicate rate in RNA-seq (1.0 at HLA-A and HLA-B), consistent with its tendency to call a single allele when the minor allele falls below its detection threshold under allele-specific expression. In contrast, WES showed near-zero false duplicate rates for most tools, supporting the use of WES-based allele calls as the primary input for downstream analyses such as LOH classification and HED computation.

#### S7.3 FIMM WES LOH Candidate Classification

Loss of heterozygosity (LOH) at HLA loci is a tumour immune escape mechanism in haematological malignancies. Using SpecHLA WES allele frequency outputs from the FIMM cohort, we classified each patient–gene pair into four categories (**Supplementary Figure S6**): candidate LOH (major allele fraction ≥ 0.80 with ≥ 20 heterozygous variant sites), allelic imbalance requiring review, balanced heterozygous, or insufficient evidence. These classifications are exploratory: WES allele frequencies reflect exon-capture read depth at HLA loci rather than genome-wide allele balance, and thresholds calibrated for WGS may not transfer reliably to WES data. WGS confirmation is required before drawing clinical conclusions from candidate LOH events identified in WES data.

#### S7.4 Overall Survival and HLA Evolutionary Divergence

In a subset of 28 FIMM patients with available survival follow-up, we computed HLA Evolutionary Divergence (HED)—the mean Grantham amino acid distance at antigen-binding groove positions—from WES consensus allele calls at HLA-A, -B, and -C. Kaplan–Meier survival curves were stratified by diagnosis group (AML, MDS, MDS→AML) and by median HED total (**Supplementary Figure S7**). The analysis is exploratory given the small cohort size (n ≈ 14 per arm for the HED split), which provides approximately 25% power to detect a hazard ratio of 2.0. Results should be interpreted as hypothesis-generating for future larger-scale studies investigating the relationship between HLA diversity and outcomes in haematological malignancies.

### Supplementary Note S8. Abstention–Accuracy Tradeoff (weighting-only consensus)

**Supplementary Figure S8** shows the abstention–accuracy tradeoff for the weighting-only WeightedConsensus across modalities: as the support threshold rises, callable rate falls and accuracy among called loci rises. This is a secondary reliability property (a "defer rather than err" option) of the weighting scheme and is not part of the primary Champion-Challenger accuracy comparison; the overall correct-call rate of WeightedConsensus remains below majority voting in every modality (main-text Table 3, ablation row).
