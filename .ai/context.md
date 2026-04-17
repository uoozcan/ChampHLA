**Objective:** Operate PIHLA as a calibrated, uncertainty-aware HLA ensemble workflow covering benchmark execution, manuscript support, and figure generation.

**Authority model:** Repo-local `.ai/` files are authoritative for AI coordination. Claude memory is a synchronized copy.

**Workspace path:** `/users/ozcanumu/scratch/project_2008084/pihla-publish/`

**Runtime path:** `/scratch/project_2008084/pihla-publish/` is the execution truth for real runs, benchmark outputs, and scientific reporting.

## Architecture
- Workflow engine: Nextflow DSL2 pipeline with module-per-tool execution.
- Execution profiles: Docker, Singularity, Puhti, and generic SLURM.
- Major scripts: `bin/run_1000g_benchmark.py`, `bin/run_1000g_figure_workflow.py`, `bin/build_1000g_benchmark_manifests.py`, `bin/hla_benchmark.py`.
- Scientific artifacts: benchmark tables, figures, manuscript drafts, and analysis docs under `analysis/` and `docs/`.

## Key commands
- `nextflow run main.nf -params-file params.yaml -profile puhti,singularity -resume`
- `python3 bin/run_1000g_benchmark.py --config conf/benchmark_1000g_realdata.yaml --output-dir results/1000g_benchmark`
- `python3 bin/run_1000g_figure_workflow.py --config conf/benchmark_1000g_realdata.yaml --output-dir analysis/1000g_realdata`
- `pytest tests/ -v`

## Coordination defaults
- Claude Code is the default implementer for pipeline, parser, script, figure, and test changes.
- Claude Code also owns direct manuscript edits.
- Codex is used mainly for brainstorming, scientific critique, evaluation, tradeoff analysis, and `/review`-style verification.
- Handoffs must be written to `.ai/handoffs/` before ownership changes.

## Active risks
- The workspace copy may drift from `/scratch/project_2008084/pihla-publish`; use the runtime path for execution truth.
- Fixture-backed outputs must not be used as scientific evidence.
- Benchmark and manuscript claims must match `benchmark_metadata.json` scope and reporting policy.

## Highest-value AI workflows
- Claude-first benchmark/debug fixes
- Codex option analysis before major scientific or workflow choices
- Claude manuscript or figure edits informed by Codex critique
- Codex verification before scientific claims when a second opinion is useful
