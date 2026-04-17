# Phase 2 Execution Plan

## Purpose

This document defines the next execution phase for PIHLA after the first real-data benchmark wave. Phase 2 is focused on turning the current project state into a stronger scientific and technical package by:

- expanding the real-data evidence base,
- improving multimodal robustness,
- refining calibration-aware consensus,
- and tightening manuscript-ready outputs.

The main principle is simple: the largest current limitation is no longer whether the workflow can run, but whether the benchmark evidence is broad and stable enough to support stronger scientific claims.

## Priority Labels

- `P0`: must-have for the next strong scientific milestone
- `P1`: important and high-value, but can follow after `P0`
- `P2`: useful later once core benchmark quality is stronger

## Phase 2 Goals

1. Expand the WGS benchmark beyond the current 42-sample wave.
2. Stabilize WES and RNA so they become true benchmark cohorts rather than pilot-only tracks.
3. Upgrade confidence weighting from guardrailed raw normalization to stronger empirical recalibration.
4. Convert the improved benchmark outputs into a cleaner manuscript and figure package.

## Workstreams

### Workstream A: WGS Benchmark Expansion

Priority: `P0`

Why this matters:
- WGS is currently the strongest evidence base in the project.
- Larger WGS cohorts will stabilize tool rankings, per-gene results, and calibration behavior.
- This gives the fastest path to stronger scientific claims.

Tasks:
1. Audit additional truth-backed 1000 Genomes WGS samples beyond the current 42-sample wave.
2. Freeze an expanded WGS cohort manifest with explicit inclusion rules.
3. Rebuild `truth_manifest.tsv`, `sequencing_manifest.tsv`, and `cohort_manifest.tsv` for the expanded wave.
4. Run the expanded WGS benchmark using the same split-aware logic.
5. Compare the expanded-wave outputs against the current 42-sample wave.

Exact deliverables:
- updated WGS cohort manifest
- expanded WGS benchmark output directory
- refreshed:
  - `method_comparison.tsv`
  - `summary_full_cohort.tsv`
  - `summary_per_gene.tsv`
  - `tool_confidence_weights.tsv`
  - `confidence_calibration_summary.tsv`
  - `benchmark_metadata.json`
- comparison note between 42-sample and expanded WGS waves

Success criteria:
- expanded WGS cohort remains truth-backed and reproducible
- holdout sample count increases or at minimum training size grows substantially
- weighted-consensus behavior remains interpretable under the larger cohort

### Workstream B: WES and RNA Stabilization

Priority: `P0`

Why this matters:
- PIHLA’s multimodal claim is currently supported mainly by a small tri-modal integration pilot.
- WES and RNA need higher tool completeness and cleaner ingestion before they can support stronger comparative results.

Tasks:
1. Build a tool-by-tool WES failure matrix from completed outputs and logs.
2. Build a tool-by-tool RNA failure matrix from completed outputs and logs.
3. Separate runtime failures from parser-ingestion failures.
4. Fix the highest-impact tool gaps first, especially tools that already emit outputs but are not benchmark-ingested cleanly.
5. Reduce phase-gated missingness in WES and RNA benchmark tables.
6. Assemble a second real multimodal benchmark wave once coverage is improved.

Exact deliverables:
- `WES_TOOL_FAILURE_MATRIX.tsv`
- `RNA_TOOL_FAILURE_MATRIX.tsv`
- parser coverage audit note
- list of fixed tool-ingestion routes
- refreshed WES/RNA benchmarkable result set
- second multimodal benchmark run directory

Success criteria:
- more WES/RNA tools contribute benchmark-ingested rows
- fewer samples are excluded for modality/tool-availability reasons
- multimodal comparison becomes stronger than a proof-of-integration pilot

### Workstream C: Confidence Recalibration

Priority: `P0`

Why this matters:
- Confidence expansion exposed a real methodological issue: raw or weakly normalized confidence can be badly miscalibrated.
- The current guardrail is strong and necessary, but the next step is a true recalibration layer.

Tasks:
1. Define a recalibration strategy for tools with enough training data.
2. Compare raw normalized confidence, effective confidence, and recalibrated confidence.
3. Add recalibrated-confidence evaluation outputs.
4. Test whether recalibrated runtime weights improve weighted consensus over the current guardrail-only version.
5. Preserve the guardrail as a safety layer even after recalibration is added.

Exact deliverables:
- recalibration design note
- updated benchmark code for recalibrated confidence support
- new outputs such as:
  - recalibrated confidence summary table
  - raw vs recalibrated confidence comparison table
  - runtime weight comparison table
- benchmark comparison:
  - majority vote
  - current guarded weighted consensus
  - recalibrated weighted consensus

Success criteria:
- recalibration improves confidence reliability for at least some tools
- recalibrated weighted consensus is at least as strong as the current guarded version
- the method remains interpretable and benchmark-derived

### Workstream D: Abstention and Discordance Maturity

Priority: `P1`

Why this matters:
- Abstention and discordance are already useful differentiators for PIHLA.
- These outputs make the ensemble more interpretable and more clinically or operationally realistic.

Tasks:
1. Refine discordance taxonomy labels where current tags are too broad.
2. Add per-gene abstention summaries.
3. Add clearer reporting of low-evidence versus tool-conflict cases.
4. Align abstention and discordance outputs with manuscript figure language.

Exact deliverables:
- refined discordance taxonomy note
- per-gene abstention summary table
- improved discordance summary tables
- updated figure captions and Results text

Success criteria:
- disagreement patterns are easier to interpret biologically or technically
- abstention outputs are understandable without reading code

### Workstream E: Manuscript and Figure Polish

Priority: `P1`

Why this matters:
- The benchmark now has a real scientific core.
- Better figures and tighter text will increase the project’s publishability without changing the underlying science.

Tasks:
1. Update the main manuscript draft after the expanded WGS wave.
2. Refresh figure text to match the latest calibration and guardrail behavior.
3. Make the WGS wave the primary Results layer.
4. Keep the tri-modal benchmark framed as integration validation and multimodal readiness.
5. Sharpen differentiation from simple multi-tool aggregation and ConshLA-like approaches.

Exact deliverables:
- updated `MANUSCRIPT_DRAFT_V1.md`
- updated `BENCHMARK_FIGURES.md`
- updated HTML presentation/report files
- final figure/table checklist

Success criteria:
- all reported numbers come from the latest real-data benchmark outputs
- the manuscript clearly explains why PIHLA is more than a workflow wrapper

### Workstream F: Workflow and User-Facing Hardening

Priority: `P1`

Why this matters:
- Strong science is more useful when the software is easier to run and easier to understand.
- Better usability will help both adoption and reproducibility.

Tasks:
1. Improve samplesheet-driven execution further.
2. Simplify benchmark configuration examples.
3. Clarify output directory structure and benchmark artifacts.
4. Improve docs for common execution paths and expected outputs.

Exact deliverables:
- improved example configs
- cleaner README benchmark section
- output interpretation guide
- user-facing run examples for WGS, WES, and RNA

Success criteria:
- a new user can understand how to run and interpret PIHLA with less manual guidance

### Workstream G: Advanced Extensions

Priority: `P2`

Why this matters:
- These are high-upside ideas, but they should not distract from strengthening the current evidence base.

Possible topics:
- richer multimodal consensus logic
- clinical or report-style exports
- GL String or HML outputs
- more advanced ensemble models
- broader patient-level integration layers

Exact deliverables:
- exploratory design notes only for now

## Recommended Execution Order

1. `P0` WGS benchmark expansion
2. `P0` WES/RNA stabilization
3. `P0` confidence recalibration
4. `P1` abstention and discordance refinement
5. `P1` manuscript and figure polish
6. `P1` workflow hardening
7. `P2` advanced extensions

## Immediate Next Package

The most practical next package to execute immediately is:

1. freeze an expanded WGS cohort manifest
2. run the expanded WGS benchmark
3. compare the expanded WGS results against the current 42-sample wave
4. produce WES and RNA tool-failure matrices
5. design the recalibration layer using the larger WGS training set

## Definition of Phase 2 Success

Phase 2 will be considered successful if the project reaches all of the following:

- a larger and more stable WGS real-data benchmark than the current 42-sample wave
- improved WES and RNA benchmark coverage beyond pilot-level integration status
- a calibration-aware consensus model that remains robust after confidence expansion
- a manuscript-ready Results and Figures package based on the stronger benchmark wave

## Guiding Principle

The core value of PIHLA is no longer just that it runs many HLA tools. Its real value is that it benchmarks them on truth-backed real data, learns from their strengths and weaknesses, and uses calibration-aware logic to decide when their evidence should or should not shape consensus. Phase 2 should strengthen that identity rather than dilute it.
