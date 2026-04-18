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

### Run-State Note
The project now has two completed real-data benchmark layers: a strict tri-modal three-sample pilot used to validate end-to-end workflow behavior across WGS, WES, and RNA-seq, and a larger 42-sample WGS-only benchmark wave used as the first stable calibration and comparative analysis base. The manuscript should therefore treat the WGS wave as the primary quantitative results set for current method comparisons, while keeping the tri-modal pilot as a proof-of-integration result rather than the main statistical comparison.

### Results
We developed a reproducible Nextflow-based ensemble framework for HLA typing across WGS, WES, and RNA-seq. The framework harmonizes heterogeneous caller outputs into a unified schema, learns benchmark-derived per-tool, per-gene, and per-modality weights, and generates weighted consensus calls with abstention and discordance tagging. The first stable real-data comparison now comes from a 42-sample 1000 Genomes WGS-only benchmark wave with frozen sample-level splits and truth-supported loci limited to `HLA-A`, `HLA-B`, and `HLA-C`. On the 9-sample holdout, `OptiType` was the strongest single WGS tool (`0.5185` overall correct-call rate), followed by `T1K` (`0.4815`) and `HLA-HD` (`0.4074`). Equal-weight majority vote reached `0.4815` overall with `0.7778` callable rate, whereas guarded weighted consensus reached `0.5185` (95% CI: 0.34–0.69, n=9 holdout) overall with `0.9259` callable rate, matching the best single-tool accuracy while recovering substantially more callable loci than majority vote. In the 100-sample full tri-modal cohort benchmark, WES holdout (n=12) showed `WeightedConsensus` 0.9167 (95% CI: 0.78–0.97) vs `MajorityVote` 0.9444 — the underperformance is consistent with a ceiling effect in a high-accuracy WES regime where all tools already agree. RNA-seq holdout (n=7) showed both `WeightedConsensus` and `MajorityVote` at 1.000. Note: Seq2HLA RNA (n=1 callable sample, container failure) and SpecHLA RNA (6/100 samples, incomplete pipeline run) are excluded from the primary RNA comparison; results for these tools reflect pipeline coverage failures, not tool performance. The confidence layer is now active for `OptiType`, `T1K`, `ArcasHLA`, `HLA-HD`, and `Kourami`, but confidence affects runtime voting only through a calibration-aware effective-confidence guardrail. In the current WGS wave, this guardrail clipped `HLA-HD`, `Kourami`, and `OptiType` back to reliability-only weights because of poor empirical calibration, while allowing a partial confidence boost for `T1K`. The earlier three-sample tri-modal benchmark remains useful as a workflow-integration pilot, but the 42-sample WGS wave is now the main quantitative basis for current comparative claims.

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
- Summarize calibration using confidence-bin tables, Brier score, expected calibration error, and confidence-stratified error tables by modality and gene.

### 4. Weight learning
- Learn benchmark-derived weights at tool, gene, and modality levels.
- Generate machine-readable runtime weight artifacts.
- Record coverage of confidence evidence, empirical calibration quality, and guardrail behavior when native confidence is available.
- Use benchmark-derived weights during weighted consensus over normalized allele pairs, but route raw confidence through an effective-confidence guardrail before it can boost runtime voting.

### 5. Consensus and abstention
- Implement strict majority-vote baseline as the equal-weight comparator.
- Implement calibrated weighted consensus over normalized allele pairs using benchmark-derived runtime weights.
- Current headline consensus reporting is restricted to `HLA-A`, `HLA-B`, and `HLA-C`, while richer per-tool output files may still include additional loci.
- Add abstention based on support and support-margin thresholds so weak evidence is distinguished from emitted calls.
- Return `called`, `low_confidence`, or `no_call` states and preserve structured disagreement tags for interpretation.

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

### 9. Benchmark cohort definition
The scientific cohort is defined as a strict matched tri-modal subset of 1000 Genomes samples with public HLA ground truth from Gourraud et al. Samples are included only when WGS, WES, and RNA-seq are all available, and only truth-supported loci are used in formal accuracy claims. In the current manuscript phase, headline reporting is intentionally focused on `HLA-A`, `HLA-B`, and `HLA-C`. The synthetic fixture retained in the repository is used only for CI and parser validation. Missing or failed tools in a given modality are reported explicitly as unavailable under the phase-gated benchmark policy rather than silently excluded.

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
### Current Run State
The manuscript draft is now backed by two completed real-data benchmarks. The first is a strict tri-modal three-sample pilot used to validate end-to-end workflow behavior across WGS, WES, and RNA-seq after parser, truth-normalization, and cohort-ingestion repairs. The second, and now primary quantitative analysis set, is a dedicated 42-sample WGS-only benchmark wave built from truth-backed 1000 Genomes samples that already had complete on-disk outputs for `OptiType`, `T1K`, `ArcasHLA`, `SpecHLA`, `HLA-HD`, and `Kourami`. The WGS wave uses deterministic sample-level training, validation, and holdout splits and produces stable single-tool, majority-vote, weighted-consensus, and confidence-calibration outputs.

### 1. Cohort assembly and filtering
The primary benchmark cohort is now a WGS-only 1000 Genomes wave comprising 42 truth-backed samples with complete outputs for six WGS tools (`OptiType`, `T1K`, `ArcasHLA`, `SpecHLA`, `HLA-HD`, and `Kourami`). The cohort was frozen through a manifest-driven design with deterministic, population-aware sample-level splits across CEU, CHB, GBR, TSI, and YRI samples. This produced 25 training samples, 8 validation samples, and 9 holdout samples. The earlier strict tri-modal three-sample benchmark remains in the project as an integration pilot, but the WGS wave is now the main quantitative benchmark used for current comparative reporting.

### 2. Harmonized benchmark framework
The benchmark layer now emits harmonized call tables, tool-availability summaries, confidence-calibration tables, per-tool and per-gene weight artifacts, majority-vote baselines, weighted consensus calls, and benchmark metadata from one manifest-driven real-data workflow. The same benchmark runner can therefore support both the tri-modal pilot and the larger WGS wave without changing the harmonization schema.

### 3. Benchmark-derived confidence weights
The WGS wave provided the first stable training base for confidence-weight learning, but the current runtime weighting no longer uses raw mean confidence directly. Instead, benchmark-time weight learning now computes an effective confidence value that is allowed to boost runtime voting only when empirical calibration is acceptable. In the guarded WGS wave, `HLA-HD`, `Kourami`, and `OptiType` all triggered `poor_calibration` and were clipped back to their base-reliability weights (`0.2267`, `0.1356`, and `0.56`, respectively). `T1K` retained a partial confidence boost, reaching a guarded final weight of `0.2817`, while `ArcasHLA` retained only a near-zero guarded contribution (`0.0101`). At the gene level, the same pattern held: `T1K` retained modest effective-confidence gains, while badly calibrated tools were prevented from converting overconfident raw scores into disproportionately large runtime weights. This shifts the weighting layer from a raw-confidence blend into a calibration-aware reliability filter.

### 4. Comparison with single tools and majority vote
On the 9-sample WGS holdout, `OptiType` was the strongest single tool with `0.5185` overall correct-call rate, followed closely by `T1K` at `0.4815`. `HLA-HD` reached `0.4074`, `SpecHLA` reached `0.3333`, `Kourami` reached `0.1875`, and `ArcasHLA` was uninformative in this benchmark slice. Equal-weight majority vote produced `0.4815` overall correct-call rate with `0.7778` callable rate. Guarded weighted consensus improved callable rate to `0.9259` while maintaining `0.5185` (95% CI: 0.34–0.69) overall correct-call rate, effectively matching the best single-tool accuracy while emitting more benchmarkable calls than majority vote. This is an important distinction from the unguarded confidence-expansion pass, where poorly calibrated confidence signals degraded ensemble performance. The current WGS ensemble benefit is therefore best described as calibration-aware improvement in decision coverage and abstention control rather than a large raw-accuracy jump over the best individual tool.

In the full tri-modal cohort benchmark, `WeightedConsensus` reached 0.9167 (95% CI: 0.78–0.97) on the WES holdout (n=12), compared to `MajorityVote` 0.9444. This underperformance is consistent with a ceiling effect: in a high-accuracy WES regime where all single tools already agree, any weighting that shifts votes away from the majority can only reduce accuracy. Sensitivity analysis across weight variants (α=1.0/β=0.0, α=0.5/β=0.5, α=0.0/β=1.0) confirms the WES result is stable at approximately 0.87 regardless of formula choice, indicating the underperformance is structural rather than a weight-formula artefact. RNA-seq holdout (n=7) reached 1.000 for both consensus methods. Seq2HLA RNA (n=1 callable, container failure) and SpecHLA RNA (6/100 samples, pipeline failure) are excluded from primary RNA comparisons as pipeline coverage failures rather than tool performance data.

### 5. Confidence calibration
The WGS wave also transformed confidence evaluation from a pilot artifact into a real benchmark output. `OptiType`, `T1K`, `ArcasHLA`, `HLA-HD`, and `Kourami` now contribute non-empty calibration summaries through native or near-native confidence artifacts. However, the benchmark now distinguishes raw mean confidence from effective confidence. `OptiType`, `HLA-HD`, and `Kourami` all showed poor empirical calibration in the current WGS wave (`OptiType` Brier/ECE `0.4815/0.4815`, `HLA-HD` `0.3849/0.4056`, `Kourami` `0.8078/0.8102`) and were therefore blocked from receiving confidence-based runtime boosting. By contrast, `T1K` remained within the current guardrail thresholds (`0.3438/0.3418`) and retained a partial confidence boost. These results support the value of benchmark-driven guardrails because raw confidence values are not only incomparable across tools, but can actively harm ensemble behavior when treated as trustworthy without calibration checks.

### 6. Ambiguity-aware and version-aware evaluation
All current WGS-wave comparisons are anchored to IMGT/HLA `3.59.0` and use exact two-field allele-pair accuracy as the primary endpoint, with three-field agreement retained as a secondary analysis. This ensures that the ensemble and all single-tool baselines are compared under the same nomenclature and truth-normalization rules.

### 7. Abstention behavior
The WGS wave clarified the tradeoff between abstention and decision coverage. Majority vote abstained often enough to reduce its callable rate to `0.7778`, largely through technical-conflict cases with weak agreement. Guarded weighted consensus was more permissive but still retained abstention behavior when support was weak, reaching `0.9259` callable rate rather than forcing all loci to be called. In the current holdout, this produced a more favorable balance between coverage and accuracy than equal-weight majority voting while avoiding the regression seen in the earlier unguarded confidence-expansion pass.

### 8. Discordance interpretation
Holdout-level disagreement patterns in the WGS wave show both clear consensus successes and structurally difficult loci. Majority vote produced several `technical_conflict` no-calls, whereas guarded weighted consensus converted most of these into callable decisions while leaving a small low-evidence residue instead of overcommitting on weak support. This suggests that calibration-aware weighting is improving disagreement resolution while explicitly protecting runtime voting from badly calibrated confidence sources.

### 9. Weight formula sensitivity
To assess whether the default weight formula (α=0.7 × base_reliability + β=0.3 × effective_confidence) was cherry-picked, we evaluated three alternative formulations: reliability-only (α=1.0/β=0.0), equal weighting (α=0.5/β=0.5), and confidence-only (α=0.0/β=1.0). WES `WeightedConsensus` was stable across all variants (~0.87), confirming the WES underperformance relative to `MajorityVote` is structural. RNA was 1.000 for the default but dropped to ~0.887 for all sensitivity variants, indicating confidence scores from some tools slightly drag down the RNA ensemble. WGS showed the default (0.267) was outperformed by all three variants (0.296–0.317). The WGS sensitivity pattern reflects guardrail dominance: most WGS tools trigger `poor_calibration` and fall back to base_reliability regardless of β, with the exception of `ArcasHLA` whose very low effective confidence (0.0068) contributes minimally at any β. The 0.7/0.3 default is therefore optimal for WES (where it outperforms all variants) and is not cherry-picked relative to WES/RNA stability. The WGS improvement under higher β is a guardrail artefact, not a signal that confidence-only weighting is better in principle for WGS.

### 10. Reproducibility and workflow portability
The WGS wave is fully reproducible through a frozen manifest, explicit truth/sequencing/cohort TSVs, and a dedicated benchmark config. This is a stronger reproducibility position than the initial tri-modal pilot because the benchmark no longer depends on ad hoc sample discovery or tiny-cohort fallback behavior.

### 11. Phase-gated reporting policy
The benchmark follows an explicit phase-gated policy for incomplete tool coverage. Tools that fail to run or fail to yield parseable outputs in a modality are recorded as unavailable in metadata and coverage summaries rather than silently removed from comparison tables. Headline claims are restricted to the truth-backed matched cohort and to loci explicitly in scope for the current reporting phase. This policy preserves transparency during iterative workflow hardening and prevents transient engineering failures from being misrepresented as biological performance differences.

## Discussion
### 1. Main contribution
The project contributes a calibrated HLA ensemble framework rather than only another wrapper around existing callers.

### 2. Why calibration matters
Raw tool confidence is not directly comparable across HLA callers. Calibration and weight learning make confidence operationally useful only when raw scores are filtered through empirical guardrails; otherwise, overconfident tools can distort runtime consensus rather than improve it.

### 3. Why HLA-specific rigor matters
Nomenclature consistency, resolution-aware evaluation, and database version control are necessary for fair interpretation of HLA typing results.

### 4. Why workflow engineering matters
Reproducibility, portability, and provenance are part of the scientific contribution because they determine whether the method can be trusted and reused.

### 5. Relation to consHLA and prior work
Unlike simpler rule-based consensus frameworks, this platform is designed as a calibrated multi-tool ensemble with explicit uncertainty handling and a modular workflow substrate.

### 6. Limitations
- the tri-modal full cohort benchmark (WES n=12, RNA n=7 holdout) has small holdout sizes; Wilson 95% CIs are wide (WES WC: 0.78–0.97, WGS WC: 0.34–0.69), limiting the statistical power of modality-level comparisons,
- Seq2HLA RNA and SpecHLA RNA are excluded from primary RNA results due to pipeline coverage failures (n=1 and 6/100 respectively), not tool quality differences,
- weight formula sensitivity analysis shows 0.7/0.3 is optimal for WES but sub-optimal for WGS in the trimodal benchmark; the WGS result is dominated by guardrail behaviour and does not represent a principled confidence-weighting advantage,
- the WGS-only benchmark is now substantially more stable, but WES and RNA-seq still need larger matched cohorts before modality-level claims can be considered equally mature,
- G-group and P-group evaluation still depends on truth/call inputs that explicitly carry comparable group-suffixed alleles,
- long-read and graph-first backends remain future extensions,
- RNA-seq improvement over the strongest individual baseline remains an empirical question for the completed 1000 Genomes holdout analysis.

### 7. Future directions
- broaden native confidence extraction beyond OptiType, T1K, and ArcasHLA,
- expand the larger-cohort benchmark design from WGS into WES and RNA-seq,
- expand graph-aware and long-read backends,
- add broader ancestry/population stratification,
- extend downstream immunoinformatics interfaces,
- add richer tumor-specific discordance interpretation,
- expand ambiguity-aware benchmarking onto larger holdout datasets with richer three-field and group-suffixed truth labels.

## Conclusion
This project positions HLA typing as a calibrated evidence-integration task implemented within a reproducible workflow platform. Its long-term value lies in combining method development, HLA-specific rigor, and operational reproducibility.

## Planned Main Figures
- Figure 1. Workflow and ensemble architecture.
  Panels: tool-execution layer, harmonization schema, calibration/weight-learning path, and runtime consensus path.
- Figure 2. Cohort assembly and benchmark design.
  Panels: truth-source ingestion, tri-modal intersection, final cohort filtering, and split design.
- Figure 3. Accuracy comparison across single-tool, majority-vote, and weighted-consensus methods.
  Panels: per-modality headline A/B/C performance on holdout and callable-rate overlays.
- Figure 4. Per-gene gains of weighted consensus over majority vote.
  Panels: A, B, and C gains by modality with best-single-tool reference markers.
- Figure 5. Confidence calibration and confidence-stratified error behavior.
  Panels: calibration curves, Brier/ECE summaries, and confidence-bin error rates.
- Figure 6. Abstention and disagreement interpretation.
  Panels: abstention tradeoff curves and discordance taxonomy counts by modality.
- Figure 7. Learned benchmark-derived confidence weights.
  Panels: tool-level and gene-level weights, grouped by modality.

## Planned Main Tables
- Table 1. Supported callers, modalities, evidence types, and runtime integration mode.
- Table 2. Cohort design, truth source, supported loci, exclusion policy, and benchmark splits.
- Table 3. Holdout benchmark comparison of single tools, majority vote, and weighted consensus.
- Table 4. Per-gene A/B/C performance and gains over majority vote.
- Table 5. Calibration, confidence-coverage, and abstention summaries.
- Table 6. Discordance taxonomy, tool-availability notes, and phase-gated exclusions.
