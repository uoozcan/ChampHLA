# Project Strategy V2

## Purpose
This document captures the working strategy that emerged from the planning and literature-review discussions for the HLA typing pipeline project. It is a structured project record rather than a verbatim chat transcript.

## Conversation Summary
The project began as a Nextflow-based multi-tool HLA typing pipeline with a consensus mechanism intended to choose the best allele prediction across tools. The main strategic concern was whether that concept was sufficiently novel compared with existing solutions, especially consHLA.

The strategy discussion then evolved in three steps:
1. assess whether read-based confidence weighting is scientifically useful,
2. identify what would make the project publishable and differentiated,
3. turn that into a practical roadmap and manuscript direction.

The main conclusion was that the strongest scientific positioning is not "a pipeline that runs multiple HLA tools," but "a calibrated, uncertainty-aware, multi-tool HLA ensemble platform" for WGS, WES, and RNA-seq.

## What Was Retained From Plan V1
Plan V1 had several strong ideas that remain central and should not be lost.

### Strengths from V1
- Build a benchmark-centered calibration layer rather than only adding runtime heuristics.
- Learn per-tool confidence weights from truth-backed benchmark data.
- Allow sidecar confidence/read-support inputs when tools do not emit native confidence.
- Produce runtime-consumable weight artifacts from benchmark outputs.
- Add a formal weighted consensus consumer instead of relying only on majority vote.
- Benchmark single tools against consensus approaches rather than making informal claims.
- Add abstention and deterministic tie-breaking.
- Treat disagreement structure as important output, not only final calls.

### Why V1 Was Valuable
V1 turned the project from a general pipeline into a method-oriented benchmark and consensus framework. That was the first major step toward novelty.

## How Plan V2 Differs From Plan V1
Plan V2 keeps the best methodological core of V1, but changes the framing and priorities in important ways.

### 1. V2 broadens the project from "consensus method" to "HLA ensemble platform"
Plan V1 focused mainly on confidence weighting and runtime consensus consumption.
Plan V2 keeps that, but adds two equally important pillars:
- HLA-specific rigor around nomenclature, IMGT/HLA versioning, and ambiguity-aware evaluation.
- workflow hardening, portability, provenance, and downstream immunoinformatics integration.

### 2. V2 makes HLA nomenclature and database rigor first-class requirements
Plan V1 did not emphasize IMGT/HLA database pinning, resolution-aware normalization, or G-group/P-group comparisons.
Plan V2 adds them immediately because the literature makes clear that these are essential to fair HLA evaluation.

### 3. V2 changes the evaluation philosophy
Plan V1 mainly targeted exact truth-backed accuracy plus confidence weighting.
Plan V2 adds ambiguity-aware evaluation and resolution-aware scoring so the benchmark reflects real HLA practice rather than only exact-string matching.

### 4. V2 integrates the literature trends directly into the roadmap
Plan V1 differentiated mainly from consHLA.
Plan V2 still does that, but also responds to broader literature themes:
- graph-aware and pangenome-aware HLA methods,
- RNA-seq and multimodal typing,
- long-read extensibility,
- reproducible workflow engineering,
- downstream immunoinformatics utility.

### 5. V2 turns workflow engineering into part of the contribution
Plan V1 treated workflow integration mostly as implementation context.
Plan V2 treats reproducibility, portability, versioning, containers, testing, and modular backend support as part of the scientific and software contribution.

## Final Strategic Positioning
The project should be positioned as:

> A calibrated, uncertainty-aware, multi-tool HLA ensemble platform for WGS, WES, and RNA-seq, with benchmark-trained evidence weighting, ambiguity-aware evaluation, and reproducible Nextflow execution.

This is stronger than:
- a simple multi-tool HLA pipeline,
- a majority-vote workflow,
- or a best-tool selector.

## Core Differentiation From consHLA
consHLA is best treated as an important reference point, but not the only one.

### consHLA-style features that overlap
- consensus as a central design goal,
- clinically meaningful integration across assay types,
- a final interpretive output rather than raw caller output only.

### Main differences our project should emphasize
- multi-tool ensemble instead of primarily one-backend consensus logic,
- benchmark-trained confidence calibration rather than rule-based agreement alone,
- explicit abstention and uncertainty states,
- structured discordance taxonomy,
- ambiguity-aware HLA evaluation,
- modular backend roadmap including graph-aware and future long-read methods,
- stronger workflow reproducibility and downstream interoperability.

## Final Prioritized Plan

### Pillar 1: Calibrated Ensemble Benchmark Core
Highest priority.
- unified call schema and confidence schema,
- benchmark-derived weights,
- weighted consensus,
- majority-vote baseline,
- abstention,
- discordance tagging,
- confidence calibration outputs.

### Pillar 2: HLA Rigor and Fair Evaluation
Immediate next priority.
- IMGT/HLA version pinning,
- nomenclature normalization,
- 2-field / 3-field / G-group / P-group evaluation,
- ambiguity-aware scoring,
- stratified robustness analyses.

### Pillar 3: Workflow Hardening and Downstream Hooks
Next layer after scientific core is stable.
- provenance and config/version capture,
- containers and test hardening,
- modular backend support,
- NetMHCpan / pVACtools-ready exports,
- report auditability.

## What We Should Keep Doing Next
1. Finish stabilizing the calibrated ensemble benchmark core in code and tests.
2. Add HLA database and nomenclature rigor before widening feature scope.
3. Build the manuscript around the method plus rigor plus reproducibility story.
4. Treat long-read and broader graph backends as modular expansion targets, not blockers for the first paper.

## Working Manuscript Message
The paper should argue that existing HLA workflows are either algorithmically strong but operationally fragmented, or operationally useful but methodologically simplistic in how they combine evidence. This project addresses that gap by combining:
- calibrated ensemble inference,
- HLA-aware evaluation rigor,
- and reproducible workflow engineering.
