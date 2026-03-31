# Manuscript Draft V1

## Working Title
A calibrated uncertainty-aware multi-tool ensemble for HLA typing from WGS, WES, and RNA-seq

## Running Title
Calibrated HLA ensemble typing

## Core Claim
We present a Nextflow-based HLA ensemble platform that harmonizes multiple HLA callers, calibrates tool confidence against truth-backed benchmarks, integrates evidence across genes and sequencing modalities, supports abstention when evidence is weak, and reports structured disagreement patterns with reproducible execution metadata.

## Abstract
### Background
HLA typing from next-generation sequencing is algorithmically rich but operationally fragmented. Different callers perform unevenly across genes, sequencing modalities, and confidence regimes, while simple consensus rules often ignore uncertainty and locus-specific behavior.

### Results
We developed a reproducible Nextflow-based ensemble framework for HLA typing across WGS, WES, and RNA-seq. The framework harmonizes heterogeneous caller outputs into a unified schema, calibrates confidence from tool-native scores and read-support sidecars, learns benchmark-derived per-tool, per-gene, and per-modality weights, and generates weighted consensus calls with abstention and discordance tagging. In the current internal fixture benchmark, weighted consensus matched the best-performing callers on WES, improved WGS overall correct-call rate from 0.8333 under strict majority vote to 1.0, and matched the strongest RNA-seq baseline at 0.8333. Confidence calibration summaries also showed modality-specific differences, with OptiType performing well on WES and substantially less well calibrated on RNA-seq and WGS than on WES.

### Conclusions
This project reframes HLA typing as a calibrated evidence-integration problem rather than a simple caller-selection or majority-voting exercise. By combining ensemble inference, HLA-specific evaluation rigor, and workflow reproducibility, the platform is designed to support both methodological benchmarking and translational immunogenomics workflows.

## Introduction
### 1. Clinical and biological motivation
HLA typing underpins transplantation, disease association analysis, pharmacogenomics, and cancer immunotherapy. Because the HLA locus is highly polymorphic and clinically consequential, both accuracy and reproducibility are critical.

### 2. Computational landscape
Many HLA callers exist for short-read DNA and RNA sequencing, including OptiType, HLA*LA, HISAT-genotype, arcasHLA, Kourami, and related tools. However, performance varies by gene, input modality, resolution, and reference assumptions.

### 3. Workflow gap
Although HLA methods are algorithmically diverse, the surrounding workflow layer remains fragmented. Reproducible execution, harmonized outputs, versioned references, and standardized benchmarking are often underdeveloped.

### 4. Consensus gap
Existing consensus strategies improve practical usability, but many rely on simple agreement rules and do not explicitly calibrate confidence, quantify uncertainty, or distinguish between weak evidence and true biological discordance.

### 5. Study objective
Our goal is to build a calibrated, uncertainty-aware, multi-tool HLA ensemble platform that integrates heterogeneous caller outputs across WGS, WES, and RNA-seq and supports reproducible benchmarking, ambiguity-aware evaluation, and downstream immunoinformatics use.

## Methods
### 1. Workflow architecture
- Nextflow-based execution.
- Containerized tool integration.
- Profile-aware deployment across local, HPC, and portable environments.
- Standardized output structure for harmonized calls, benchmark summaries, consensus outputs, and manuscript figures.

### 2. Tool harmonization layer
- Parse heterogeneous tool outputs into a unified HLA call schema.
- Preserve raw allele strings and normalized comparison-ready forms.
- Normalize by gene and allele resolution.
- Record runtime and provenance metadata.

### 3. Confidence modeling
- Extract tool-native confidence where available.
- Accept sidecar read-support or confidence files when native confidence is absent.
- Convert confidence evidence to normalized scores in `[0,1]`.
- Summarize calibration using confidence-bin tables, Brier score, and expected calibration error.

### 4. Weight learning
- Learn benchmark-derived weights at tool, gene, and modality levels.
- Generate machine-readable runtime weight artifacts.
- Record coverage of confidence evidence and fallback behavior when native confidence is unavailable.
- Use benchmark-derived weights during weighted consensus over normalized allele pairs.

### 5. Consensus and abstention
- Implement strict majority-vote baseline.
- Implement weighted consensus over normalized allele pairs.
- Add abstention based on support and support-margin thresholds.
- Return `called`, `low_confidence`, or `no_call` states.

### 6. Discordance taxonomy
- Label disagreement structure using tags such as:
  - `technical_conflict`
  - `low_evidence_conflict`
  - `dna_rna_discordance`
  - `possible_expression_bias`
  - `no_evidence`
- Treat discordance interpretation as an analysis layer distinct from primary truth benchmarking.

### 7. Ambiguity-aware HLA evaluation
- 2-field evaluation remains the current primary endpoint.
- The benchmark now also emits 3-field exact-match rates plus G-group and P-group match rates when truth and calls provide comparable group-suffixed alleles.
- IMGT/HLA database version is pinned in the benchmark configuration and propagated through `benchmark_metadata.json`, `reference_metadata.tsv`, and `harmonized_benchmark_rows.tsv`.
- Group-level metrics are intentionally blank when comparable G-group or P-group encodings are absent.

### 8. Benchmark design
- Use truth hierarchy:
  1. orthogonal clinical typing
  2. targeted HLA NGS
  3. public reference truth
- Use split hierarchy:
  - training for calibration and weight learning,
  - validation for abstention threshold tuning,
  - holdout for final reporting.

### 9. Current internal benchmark cohort
The current internal fixture benchmark used during development contains two synthetic truth-backed samples across genes A, B, and C, with WES, WGS, and RNA-seq tool outputs available for the supported benchmark callers. These results are useful for validating method behavior and output structure, but they should be treated as provisional until larger holdout benchmarking is complete.

### 10. Metrics
Primary:
- 2-field exact pair accuracy
- callable rate
- overall correct-call rate

Secondary:
- 3-field accuracy
- ambiguity-aware agreement
- per-gene performance
- calibration error
- abstention tradeoff
- discordance burden

## Results
### 1. Harmonized benchmark framework
The current benchmark layer produces harmonized call tables, confidence summaries, per-tool/per-gene weight artifacts, strict majority-vote baselines, weighted consensus calls, calibration tables, abstention tradeoff summaries, and discordance tables from one unified benchmark configuration.

### 2. Benchmark-derived confidence weights
The current internal benchmark successfully differentiates confidence-weighted behavior across tools and modalities. OptiType shows strong apparent confidence and accuracy on WES, while calibration quality declines on WGS and RNA-seq. ArcasHLA uses sidecar read-support-derived confidence and remains reasonably calibrated on the current RNA-seq fixture set.

### 3. Comparison with single tools and majority vote
In the current fixture benchmark:
- WES performance is perfect for the strongest tools, strict majority vote, and weighted consensus.
- WGS weighted consensus improves overall correct-call rate to 1.0, compared with 0.8333 for OptiType alone and 0.8333 for strict majority vote.
- RNA-seq weighted consensus currently matches the strongest RNA-seq baseline at 0.8333 overall correct-call rate rather than improving beyond it.

These results suggest that the benchmark-trained weighting layer already adds value in the mixed-confidence WGS setting, while RNA-seq remains a harder regime where more evidence sources or richer caller diversity may be needed.

### 4. Confidence calibration
Current calibration summaries show clear modality dependence. OptiType on WES has mean confidence 0.985 with observed accuracy 1.0, whereas OptiType on WGS and RNA-seq shows larger calibration error. The current RNA-seq ArcasHLA sidecar-based confidence profile produces observed accuracy 0.8333 with mean confidence 0.7533. These early outputs support the motivation for calibrated rather than raw cross-tool confidence use.

### 5. Ambiguity-aware and version-aware evaluation
The benchmark now records IMGT/HLA version 3.59.0 in both machine-readable metadata and row-level benchmark outputs. The new ambiguity summary layer shows that some tools retain perfect two-field performance while losing agreement at three-field resolution. In the current fixture benchmark, SpecHLA remains at 1.0 exact two-field rate on both WES and WGS but drops to 0.0 exact three-field rate, while OptiType retains 1.0 exact three-field rate on WES and 0.8333 on WGS. These results validate the need for resolution-aware reporting rather than a single exact-match metric.

### 6. Abstention behavior
The benchmark now emits abstention tradeoff tables across support thresholds. In the current small fixture cohort, weighted consensus mostly emits calls rather than abstentions, which is expected because several loci remain unambiguous in this development set. Larger holdout cohorts will be required to characterize abstention utility more realistically.

### 7. Discordance interpretation
The current discordance summary identifies one cross-modality `dna_rna_discordance` event in the fixture benchmark. This validates the basic disagreement taxonomy wiring and provides a starting point for richer biological and technical discordance analysis in larger cohorts.

### 8. Reproducibility and workflow portability
The benchmark layer now produces machine-readable tables and figure-ready outputs in a consistent structure. This is useful not only for manuscript generation but also for portable regression testing and future workflow hardening.

## Discussion
### 1. Main contribution
The project contributes a calibrated HLA ensemble framework rather than only another wrapper around existing callers.

### 2. Why calibration matters
Raw tool confidence is not directly comparable across HLA callers. Calibration and weight learning make confidence operationally useful.

### 3. Why HLA-specific rigor matters
Nomenclature consistency, resolution-aware evaluation, and database version control are necessary for fair interpretation of HLA typing results.

### 4. Why workflow engineering matters
Reproducibility, portability, and provenance are part of the scientific contribution because they determine whether the method can be trusted and reused.

### 5. Relation to consHLA and prior work
Unlike simpler rule-based consensus frameworks, this platform is designed as a calibrated multi-tool ensemble with explicit uncertainty handling and a modular workflow substrate.

### 6. Limitations
- the current reported quantitative results come from a small internal fixture benchmark,
- G-group and P-group evaluation still depends on truth/call inputs that explicitly carry comparable group-suffixed alleles,
- long-read and graph-first backends remain future extensions,
- RNA-seq improvement over the strongest individual baseline is not yet demonstrated in the current internal benchmark.

### 7. Future directions
- expand graph-aware and long-read backends,
- add broader ancestry/population stratification,
- extend downstream immunoinformatics interfaces,
- add richer tumor-specific discordance interpretation,
- expand ambiguity-aware benchmarking onto larger holdout datasets with richer three-field and group-suffixed truth labels.

## Conclusion
This project positions HLA typing as a calibrated evidence-integration task implemented within a reproducible workflow platform. Its long-term value lies in combining method development, HLA-specific rigor, and operational reproducibility.

## Planned Main Figures
- Figure 1. Workflow and ensemble architecture
- Figure 2. Accuracy comparison across single-tool, majority-vote, and weighted-consensus methods
- Figure 3. Per-gene gains of weighted consensus over majority vote
- Figure 4. Confidence calibration curves and summaries
- Figure 5. Abstention tradeoff curves
- Figure 6. Discordance taxonomy across modalities
- Figure 7. Learned benchmark-derived confidence weights

## Planned Main Tables
- Table 1. Supported callers, modalities, and evidence types
- Table 2. Cohort design, truth hierarchy, and benchmark splits
- Table 3. Main benchmark comparison results
- Table 4. Per-gene performance and gains
- Table 5. Calibration and abstention summaries
- Table 6. Discordance taxonomy summary
