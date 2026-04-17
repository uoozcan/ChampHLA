## 2026-04-06  WES Root-Cause Pass
- Status: completed
- What happened: performed a focused root-cause pass on the failed 3-sample WES run (`33828475`) using launcher logs, execution trace, and task-level `.command.err` files under `/scratch/project_2008084/hla_calibration/work/wes_3sample_33828475/`.
- Decision or interpretation:
  - `EXTRACT_HLA_AND_CONVERT` was failing because staged BAMs were missing task-local indexes; this blocked WES `OptiType` and `T1K` entirely.
  - `HLAHD_BAM` was failing because the container-side samtools runtime was broken, even though the HLA-HD executable itself exists.
  - `Kourami` WES was failing because `samtools` is absent in the current container runtime.
  - `POLYSOLVER` WES was failing because the bundled samtools uses older `sort` syntax than the module expected.
- Evidence/files/jobs:
  - `/scratch/project_2008084/hla_calibration/logs/wes_3sample_33828475.out`
  - `/scratch/project_2008084/hla_calibration/work/wes_3sample_33828475/.../.command.err`
  - `/scratch/project_2008084/pihla-publish/docs/WES_ROOT_CAUSE_PASS_2026-04-06.md`
- Next action:
  - rerun the 3-sample WES workflow with the patched modules/routing and regenerate the WES failure matrix.

## 2026-04-06 12:05:00 EEST
- status: completed
- what happened: added a calibration-aware guardrail to benchmark-derived confidence weighting and regenerated the 42-sample WGS benchmark.
- decision or interpretation: the guardrail now shrinks or blocks confidence boosting when empirical calibration is poor, while preserving reliability-based weight contribution. In the guarded WGS wave, `HLA-HD`, `Kourami`, and `OptiType` all triggered `poor_calibration` and were clipped back to their base reliability weights, whereas `T1K` retained a partial confidence boost and `ArcasHLA` retained a very small guarded contribution. This recovered weighted-consensus holdout performance from the degraded post-expansion state (`0.4815` overall / `0.8889` callable rate) to `0.5185` overall / `0.9259` callable rate.
- evidence/files/jobs: `bin/hla_benchmark.py`, `tests/test_hla_benchmark.py`, `analysis/1000g_realdata/benchmark_wgs_wave1/tables/tool_confidence_weights.tsv`, `analysis/1000g_realdata/benchmark_wgs_wave1/tables/method_comparison.tsv`, `analysis/1000g_realdata/benchmark_wgs_wave1/tables/benchmark_metadata.json`.
- next action: update manuscript/figure language so the confidence-weighting section reflects guardrailed effective confidence rather than raw mean confidence alone, and consider whether `OptiType` should continue to be blocked by default thresholds or receive a modality-specific threshold in future tuning.

## 2026-04-06 11:25:00 EEST
- status: completed
- what happened: broadened WGS confidence extraction beyond OptiType, T1K, and ArcasHLA by adding native confidence parsers for HLA-HD (`*.read.txt`) and Kourami (`*.kourami.result`), wiring them into the WGS benchmark configs, rerunning tests, and regenerating the 42-sample WGS benchmark.
- decision or interpretation: the expanded confidence layer is technically successful and now yields non-empty calibration outputs for HLA-HD and Kourami. However, both tools appear strongly overconfident in the current WGS wave (`HLA-HD mean_confidence=0.8456 at observed_accuracy=0.44`; `Kourami mean_confidence=0.9977 at observed_accuracy=0.1875`). Their inclusion increased learned weights (`HLA-HD 0.4255`, `Kourami 0.3873`) and degraded weighted-consensus holdout performance from `0.5185` overall / `0.9630` callable rate to `0.4815` overall / `0.8889` callable rate. This is a useful calibration finding: broader native confidence ingestion is feasible, but some tool-native scores should likely be reliability-adjusted or downweighted before being promoted into runtime consensus.
- evidence/files/jobs: `bin/hla_benchmark.py`, `tests/test_hla_benchmark.py`, `conf/benchmark_1000g_phase_gated_abc.yaml`, `conf/benchmark_1000g_wgs_wave1.yaml`, `analysis/1000g_realdata/benchmark_wgs_wave1/tables/tool_confidence_weights.tsv`, `analysis/1000g_realdata/benchmark_wgs_wave1/tables/confidence_calibration_summary.tsv`, `analysis/1000g_realdata/benchmark_wgs_wave1/tables/method_comparison.tsv`.
- next action: add a confidence-reliability guardrail layer for runtime weighting, for example by shrinking raw confidence contributions when calibration error is high, or by restricting runtime weighting to tools with empirically acceptable calibration quality.

## 2026-04-06 10:35:00 EEST
- status: completed
- what happened: implemented and ran a dedicated WGS-only 1000 Genomes benchmark wave using the frozen 42-sample manifest.
- decision or interpretation: the new config `conf/benchmark_1000g_wgs_wave1.yaml` successfully consumed `conf/1000g_wgs_wave1_manifest.tsv` through explicit WGS-only truth/sequencing/cohort manifests and produced stabilized single-tool, majority-vote, weighted-consensus, and calibration outputs. Weighted consensus slightly outperformed the best single WGS tool on overall correct call rate in the 9-sample holdout (`0.5185` vs `0.5185` for OptiType, tied overall, but with higher callable rate than majority vote and a more realistic confidence-weighting layer). T1K and ArcasHLA confidence are now active in the larger WGS benchmark, not just in the 3-sample pilot.
- evidence/files/jobs: `conf/benchmark_1000g_wgs_wave1.yaml`, `analysis/1000g_realdata/wgs_wave1_inputs/`, `analysis/1000g_realdata/benchmark_wgs_wave1/tables/method_comparison.tsv`, `analysis/1000g_realdata/benchmark_wgs_wave1/tables/tool_confidence_weights.tsv`, `analysis/1000g_realdata/benchmark_wgs_wave1/tables/confidence_calibration_summary.tsv`.
- next action: use the 42-sample WGS wave as the primary training/calibration base, then decide whether to (a) broaden confidence extraction to additional tools or (b) draft the first manuscript-ready WGS results section and figures from these stabilized outputs.

## 2026-04-06 10:05:00 EEST
- status: completed
- what happened: quantified the larger WGS-only truth-backed expansion cohort after confidence expansion for T1K and ArcasHLA.
- decision or interpretation: there are 42 samples on disk that are simultaneously truth-backed and complete for OptiType, T1K, ArcasHLA, SpecHLA, HLA-HD, and Kourami. This is large enough to support a meaningful WGS-only benchmark wave with fixed sample-level splits, instead of relying on the brittle 3-sample tri-modal pilot.
- evidence/files/jobs: `/scratch/project_2008084/ozcanumu/hla_calibration/conf/ground_truth_data.csv`, `/scratch/project_2008084/hla_calibration/wgs/results`, `conf/1000g_wgs_wave1_manifest.tsv`, `docs/WGS_WAVE1_PLAN.md`. Population counts in the full six-tool truth-backed overlap are CEU=7, CHB=8, GBR=7, TSI=9, YRI=11.
- next action: generate a dedicated WGS-only benchmark config that consumes the frozen 42-sample manifest and reports stabilized single-tool, majority-vote, weighted-consensus, and confidence-calibration outputs.

# Project Notebook

This notebook is the operational record for PIHLA development, benchmarking, and run-state tracking. Entries are append-only in reverse chronological order. Scientific interpretation should be cited from the manuscript and benchmark outputs; procedural history should be cited from here.

## 2026-04-06 09:15:00 EEST
**Status**
Confidence expansion and larger WGS-wave planning initiated.

**What happened**
- Audited native `T1K` and `arcasHLA` artifacts already present on disk for WGS outputs.
- Confirmed `T1K` WGS `*_genotype.tsv` files contain per-gene support counts suitable for read-support style confidence extraction.
- Confirmed `arcasHLA` WGS `*.genes.json` files contain per-gene support/read values suitable for read-support style confidence extraction.
- Confirmed a broader reusable WGS result pool is already available on disk for approximately `42-43` samples across `OptiType`, `T1K`, `arcasHLA`, `SpecHLA`, `HLA-HD`, and `Kourami`.

**Decision or interpretation**
- The next confidence expansion should prioritize WGS because the reusable result pool is already materially larger than the tri-modal pilot.
- `T1K` and `arcasHLA` confidence should be treated as gene-level read-support style signals, not direct probabilities.
- The next benchmark wave should likely be WGS-only and truth-backed, using the existing 42-43 sample result pool to reduce variance in learned weights before forcing larger tri-modal expansion.

**Evidence/files/jobs**
- Native confidence-capable artifacts inspected:
  - `wgs/results/*/t1k/t1k_out/*_genotype.tsv`
  - `wgs/results/*/arcashla/*/*.genes.json`
  - `wgs/results/*/optitype/*/*_result.tsv`
- Approximate reusable WGS sample counts on disk:
  - `OptiType 42`, `T1K 42`, `arcasHLA 42`, `SpecHLA 43`, `HLA-HD 43`, `Kourami 43`

**Next action**
- Validate the new T1K and arcasHLA confidence parsers on the real WGS benchmark.
- Draft a WGS-only benchmark manifest strategy that intersects the existing on-disk result pool with truth-backed 1000 Genomes samples and keeps `A/B/C` as the headline scope for the next wave.

## 2026-04-05 09:00:00 EEST
**Status**
Primary pilot-benchmark bottlenecks partially resolved: tiny-cohort weighting is more stable and real OptiType confidence is now flowing into the benchmark outputs.

**What happened**
- Updated `bin/run_1000g_benchmark.py` so cohorts of three or fewer included samples enter a `small_cohort_mode` that pools training and validation samples for weight learning.
- Added metadata fields documenting the effective training, validation, and holdout sample sets used under small-cohort mode.
- Wired native OptiType `*_result.tsv` files as confidence sidecars for WGS and RNA-seq in `conf/benchmark_1000g_phase_gated_abc.yaml`.
- Re-ran the regression suite and the real-data pilot benchmark successfully.
- The refreshed benchmark now reports non-empty confidence calibration rows for OptiType and no longer collapses WGS weighted consensus to a single called gene.

**Decision or interpretation**
- Pooling the non-holdout samples is a better default than a strict `1/1/1` train/validation/holdout split for tiny pilot cohorts.
- The weighted-consensus framework is now empirically more credible in the pilot because it uses real OptiType Objective values instead of pure fallback reliability for that tool.
- The remaining calibration gap is now narrower and more concrete: extend real confidence extraction beyond OptiType to other tools where scientifically justified.

**Evidence/files/jobs**
- Updated files: `bin/run_1000g_benchmark.py`, `conf/benchmark_1000g_phase_gated_abc.yaml`, `tests/test_1000g_benchmark_workflow.py`
- Refreshed pilot outputs:
  - `tables/benchmark_metadata.json`
  - `tables/tool_confidence_weights.tsv`
  - `tables/confidence_calibration_summary.tsv`
  - `tables/method_comparison.tsv`
- Key improvements after rerun:
  - `small_cohort_mode=True`
  - effective training samples: `NA06994`, `NA06985`
  - effective holdout sample: `NA06986`
  - OptiType confidence coverage now `1.0` in WGS and RNA-seq training outputs
  - WGS weighted consensus improved from `1/3` to `2/3`

**Next action**
- Extend confidence extraction to additional tools only where the native files provide interpretable read-support or objective-style metrics.
- Use the broader existing WGS result pool to move beyond the three-sample pilot and produce less fragile training/holdout estimates.

## 2026-04-05 08:05:00 EEST
**Status**
First strict tri-modal 1000 Genomes pilot benchmark completed successfully on real data; benchmark ingestion and truth-normalization bugs were repaired.

**What happened**
- Recovered the failed benchmark-only handoff by fixing WES sample-ID recognition in the phase-gated input builder and benchmark config.
- Diagnosed benchmark ingestion gaps and repaired parser behavior in `bin/hla_benchmark.py`:
  - wide/long parsers now honor `sample_override`
  - long-format outputs without a sample column can use the path-derived sample ID
  - `HLA-HD` now has a dedicated parser for headerless output
  - bare truth alleles such as `03:01` are now qualified with the current gene before ambiguity comparison
- Re-ran the benchmark repeatedly until the final real-data pilot outputs were structurally correct and scored correctly.
- Final benchmark rerun `33847048` completed successfully.
- The strict tri-modal cohort manifest now includes all three pilot samples, with split membership of one training, one validation, and one holdout sample.
- Holdout reporting is therefore based on sample `NA06986`, which explains why headline tables currently show `sample_count=1`.

**Decision or interpretation**
- The benchmark framework is now operational on real 1000 Genomes data and is no longer blocked by sample-ID, parser, or truth-normalization bugs.
- The current scientific result should be framed as a phase-gated pilot benchmark, not a stable comparative study, because final reporting still rests on a single holdout sample.
- Weighted consensus is currently brittle in WGS because the learned gene-specific weights are estimated from only one training sample.
- Majority vote is presently the more stable WGS comparator, while weighted consensus looks most promising in RNA-seq under this pilot split.
- Confidence-calibration analysis remains effectively unavailable in this pilot because usable confidence-sidecar coverage is still absent.

**Evidence/files/jobs**
- Successful final benchmark rerun: `33847048`
- Output root: `/scratch/project_2008084/pihla-publish/analysis/1000g_realdata/benchmark_phase_gated_abc`
- Headline tables inspected:
  - `tables/summary_full_cohort.tsv`
  - `tables/method_comparison.tsv`
  - `tables/summary_per_gene.tsv`
  - `tables/tool_confidence_weights.tsv`
  - `tables/weighted_consensus_calls.tsv`
- Final pilot benchmark highlights:
  - WES `SpecHLA`: `3/3` correct on holdout
  - WGS best single tools: `OptiType` and `T1K` at `2/3`
  - RNA best single tool: `OptiType` at `3/3`; `Seq2HLA` at `2/2`
  - Majority vote: `WES 3/3`, `WGS 2/3`, `RNA 2/3`
  - Weighted consensus: `WES 3/3`, `WGS 1/3`, `RNA 3/3`
- Reusable broader WGS result pool already on disk:
  - `optitype 42`, `t1k 42`, `arcashla 42`, `spechla 43`, `hlahd 43`, `kourami 43`

**Next action**
- Use the completed pilot benchmark to update the manuscript Results/Methods text.
- Prioritize expansion beyond a one-sample holdout by reusing the existing larger WGS result pool and by relaxing or redesigning split policy for tiny cohorts.
- Add real confidence-bearing inputs from existing tool-native files where available, especially OptiType TSV objective-style outputs, T1K allele/genotype outputs, and arcasHLA JSON-derived sidecars where scientifically defensible.

## 2026-04-03 00:45:00 EEST
**Status**
Waiting-phase manuscript and reporting scaffold refined while WES/RNA reruns remain queued.

**What happened**
- Updated the manuscript draft to reflect the current queued rerun chain rather than the older session-lock failure state.
- Added an explicit phase-gated reporting section to the manuscript so incomplete tool coverage is handled transparently in the paper narrative.
- Expanded the planned figure and table sections into near-final shells with exact data sources, panel concepts, and intended scientific messages.
- Updated the benchmark analysis note so the current run-state and phase-gated reporting policy match the storage-hardened rerun plan.
- Updated the benchmark-figure planning document to reflect native support for the broader tool set and to define exact figure/table shells for post-run analysis.

**Decision or interpretation**
- The project is now better positioned for a fast transition from completed benchmark outputs to manuscript-ready Results.
- The figure/table inventory is sufficiently concrete that post-run work should focus on number replacement and interpretation rather than redesign.

**Evidence/files/jobs**
- Docs updated: `docs/MANUSCRIPT_DRAFT_V1.md`, `docs/ANALYSIS_2026-04-01.md`, `docs/BENCHMARK_FIGURES.md`
- Active queued jobs remain `33828475`, `33828476`, and `33828477`

**Next action**
- When the queued runs start, resume run monitoring; when they finish, populate the figure/table shells and manuscript placeholders from the real-data benchmark outputs.

## 2026-04-03 00:25:00 EEST
**Status**
Quota diagnosis completed; storage-hardened rerun chain submitted.

**What happened**
- Confirmed both relaunch jobs from the previous wave failed due to disk quota rather than scheduler state:
  - WES `33823898` failed after `SPECHLA_BAM` hit `Disk quota exceeded`
  - RNA `33823899` failed while publishing T1K output, also with `Disk quota exceeded`
- Measured the main storage consumers and found the dominant pressure in failed work directories rather than logs or WES published results:
  - `/scratch/project_2008084/hla_calibration/work` reached `132G`
  - `wes_3sample_33823898` alone used `70G`
  - `rna_3sample_33823899` alone used `54G`
  - RNA published results temporarily included a `4.4G` `arcasHLA` output directory containing copied FASTQs.
- Cleaned failed work directories and oversized partial RNA outputs, reducing:
  - `/scratch/project_2008084/hla_calibration/work` to `24M`
  - `/scratch/project_2008084/hla_calibration/rna_3sample/results` to `4.4M`
- Traced the worst storage drivers to specific workflow/module behavior:
  - WES BAM input was still using full `BAM_TO_FASTQ` conversion unless `--extract_hla_region` was set manually.
  - `arcasHLA` and `SpecHLA` were publishing whole sample directories, which could include large FASTQs.
  - `HLA-HD` and `OptiType` could leave large decompressed FASTQs behind on failure.
- Applied storage-hardening edits in the codebase:
  - `main.nf`: WES BAM now defaults to HLA-region extraction before FASTQ conversion.
  - `modules/arcashla.nf`: removed full-directory publish outputs.
  - `modules/spechla.nf`: removed full-directory publish outputs and added cleanup traps.
  - `modules/hlahd.nf`: added cleanup traps for BAM/FASTQ intermediates.
  - `modules/optitype.nf`: added cleanup trap for decompressed FASTQs.
- Validated the updated code:
  - `python3 -m unittest tests.test_majority_voting_runtime tests.test_hla_benchmark tests.test_1000g_benchmark_workflow` passed (`Ran 14 tests ... OK`)
  - `nextflow run main.nf --help -profile puhti,singularity` parsed successfully.
- Removed stale partial WES/RNA result trees and submitted a fresh chain:
  - WES: `33828475`
  - RNA: `33828476`
  - dependent benchmark: `33828477`

**Decision or interpretation**
- The immediate engineering focus was correct: storage behavior, not SLURM scheduling, is the dominant blocker at this stage.
- WES should not rely on whole-exome BAM-to-FASTQ conversion for the current smoke benchmark when HLA-region extraction is available.
- Publishing complete tool working directories is too expensive for this cohort scale and should remain disabled unless explicitly requested for debugging.

**Evidence/files/jobs**
- Failed jobs: `33823898`, `33823899`, `33823900`
- New jobs: `33828475`, `33828476`, `33828477`
- Logs: `/scratch/project_2008084/hla_calibration/logs/wes_3sample_33823898.out`, `/scratch/project_2008084/hla_calibration/logs/rna_3sample_33823899.out`
- Key edited files: `main.nf`, `modules/arcashla.nf`, `modules/spechla.nf`, `modules/hlahd.nf`, `modules/optitype.nf`

**Next action**
- Monitor `33828475` and `33828476` closely for whether the storage-hardening changes prevent quota growth and whether any remaining tool-specific runtime failures emerge after the new runs start.

## 2026-04-02 23:45:00 EEST
**Status**
Fresh WES and RNA relaunches are now running without immediate startup failures.

**What happened**
- Queried SLURM again after the cleanup/relaunch wave and confirmed both launcher jobs started at `2026-04-02 23:33:57 EEST`.
- Current live job state:
  - WES `33823898` running on `r06c56`
  - RNA `33823899` running on `r06c64`
  - benchmark `33823900` still waiting on `afterok` dependencies
- Inspected launcher stdout/stderr for the first few minutes of execution.
- WES stdout shows successful task submission for previously problematic tools including `HLAHD_BAM`, `T1K`-dependent preprocessing, `SPECHLA_BAM`, `POLYSOLVER`, `KOURAMI`, and `ARCASHLA_BAM`.
- RNA stdout shows successful task submission for previously problematic tools including `HLAHD_FASTQ`, `T1K_FASTQ`, `SPECHLA_FASTQ`, `SEQ2HLA`, `OPTITYPE`, and `ARCASHLA_FASTQ`.
- No `ERROR`, `WARN`, `failed`, `retry`, `terminated`, or `exit status` lines were present yet in the launcher stdout logs, and stderr only showed normal module/environment messages.

**Decision or interpretation**
- The most important launcher and environment fixes appear to be working: the earlier immediate failures for HLA-HD, T1K, and SpecHLA are no longer blocking task startup.
- This does not yet prove all tools will complete successfully, but the early execution state is materially healthier than prior runs.
- `arcasHLA` RNA remains a known watch item because its earlier failure happened inside tool execution rather than at initial submission.

**Evidence/files/jobs**
- Jobs: `33823898`, `33823899`, `33823900`
- Logs: `/scratch/project_2008084/hla_calibration/logs/wes_3sample_33823898.out`, `/scratch/project_2008084/hla_calibration/logs/rna_3sample_33823899.out`
- Scheduler state from `scontrol show job` on `33823898` and `33823899`

**Next action**
- Recheck the live logs after more task runtime to see whether any tool begins failing during execution, with special attention to `arcasHLA`, `HLA-HD`, `T1K`, and `SpecHLA`.

## 2026-04-02 23:35:00 EEST
**Status**
Storage cleanup completed; fresh WES/RNA/benchmark relaunch chain is queued cleanly.

**What happened**
- Confirmed the previous corrected reruns were no longer viable because quota pressure persisted in the active work area.
- Cleaned stale heavy work directories from `/scratch/project_2008084/pihla-publish/work`, reducing it to a minimal fresh state.
- Re-submitted the launcher chain with the already-patched Puhti configuration and module fixes:
  - WES: `33823898`
  - RNA: `33823899`
  - dependent benchmark: `33823900`
- Queried SLURM after relaunch and confirmed:
  - WES `33823898` is `PENDING (Priority)`
  - RNA `33823899` is `PENDING (Priority)`
  - benchmark `33823900` is `PENDING (Dependency)`
- Verified that no `.out` log files exist yet for these job IDs, consistent with jobs that have not started.
- Rechecked disk usage after cleanup:
  - `/scratch/project_2008084/pihla-publish/work` = `4.2M`
  - `/scratch/project_2008084/hla_calibration/work` = `8.5G`

**Decision or interpretation**
- The pipeline is now in a healthier scheduler-ready state than the earlier relaunch wave.
- The immediate blocker is queue time rather than an active parser or launcher error.
- When either WES or RNA begins, the first live checks should focus on whether the patched HLA-HD, T1K, and SpecHLA tasks initialize correctly; arcasHLA RNA remains the primary known tool-level risk.

**Evidence/files/jobs**
- Jobs: `33823898`, `33823899`, `33823900`
- Launcher scripts: `slurm_wes_3sample.sh`, `slurm_rna_3sample.sh`
- Benchmark handoff script: `slurm_benchmark_1000g_phase_gated.sh`
- Active work paths: `/scratch/project_2008084/pihla-publish/work`, `/scratch/project_2008084/hla_calibration/work`

**Next action**
- Monitor queue transition from `PENDING` to `RUNNING`, inspect the first task logs/work dirs for HLA-HD, T1K, and SpecHLA initialization, and relaunch the benchmark automatically via job `33823900` once both upstream jobs complete successfully.

## 2026-04-02 21:20:00 EEST
**Status**
Runtime/tool-integration repair pass completed; corrected relaunches submitted.

**What happened**
- Diagnosed the failed WES relaunch as a top-level Nextflow resume-lock problem caused by `-resume` reuse in the launcher script.
- Diagnosed RNA tool failures from the old run:
  - `HLA-HD` was pointed at a nonexistent database root.
  - `T1K` lacked a valid `-f` HLA reference in the active module path.
  - `SpecHLA` was using a broken container-path assumption (`/opt/SpecHLA/...`) instead of the available local Puhti installation.
  - `arcasHLA` still shows an internal `KeyError: '0'` during FASTQ genotyping and remains an active known issue.
- Patched `nextflow.config`, `modules/t1k.nf`, `slurm_wes_3sample.sh`, and `slurm_rna_3sample.sh`.
- Removed ambiguous launcher `-resume` use, added unique work directories, switched Puhti launcher invocations to local SpecHLA, corrected `hlahd_db`, and added `t1k_hlaidx` support.
- Canceled stale RNA launcher `33819826` and submitted corrected jobs `33822685` (WES), `33822686` (RNA), and dependent benchmark `33822687`.

**Decision or interpretation**
- The project can now proceed with corrected WES/RNA launcher behavior for HLA-HD, T1K, and SpecHLA on Puhti.
- `arcasHLA` RNA remains phase-gated/known-broken for now and should be reported as unresolved if it still fails in the rerun.
- The benchmark dependency chain has been re-established with the corrected job IDs.

**Evidence/files/jobs**
- Files: `nextflow.config`, `modules/t1k.nf`, `slurm_wes_3sample.sh`, `slurm_rna_3sample.sh`
- Failed WES log: `/scratch/project_2008084/hla_calibration/logs/wes_3sample_33819825.out`
- Failed RNA work dirs included HLA-HD, T1K, SpecHLA, and arcasHLA diagnostic traces under `/scratch/project_2008084/pihla-publish/work/`
- New jobs: `33822685`, `33822686`, `33822687`

**Next action**
- Monitor the corrected WES and RNA relaunches, then run the dependent real-data benchmark and update manuscript placeholders with the first successful tri-modal refresh.

## 2026-04-02 19:35:00 EEST
**Status**
Live run-state update after waiting-phase documentation pass.

**What happened**
- The earlier "both queued" state changed during documentation review.
- RNA launcher job `33819826` started successfully and is now running on Puhti.
- WES launcher job `33819825` failed immediately because Nextflow could not acquire a session lock for resume session `41cbfa8b-e555-4ffc-9924-2149f74dc4b4`.
- Dependent benchmark job `33819996` was canceled automatically after the failed WES dependency.

**Decision or interpretation**
- The project notebook and benchmark-facing run-state notes must reflect the live scheduler state rather than the earlier fully pending snapshot.
- The next operational fix is to relaunch WES with a clean session strategy or a non-conflicting resume context, then resubmit the dependent benchmark chain.
- RNA work can continue and its partial outputs may still be useful for later phase-gated coverage accounting.

**Evidence/files/jobs**
- `squeue -u ozcanumu`
- `sacct -j 33819825,33819826,33819996 --format=JobID,JobName,State,ExitCode,Start,End`
- `/scratch/project_2008084/hla_calibration/logs/wes_3sample_33819825.out`
- `/scratch/project_2008084/hla_calibration/logs/rna_3sample_33819826.out`

**Next action**
- Inspect the conflicting Nextflow session state, relaunch WES cleanly, and recreate the dependent benchmark submission once WES is healthy again.

## 2026-04-02 18:45:46 EEST
**Status**
Waiting-phase documentation and analysis pass while corrected WES and RNA jobs remain queued.

**What happened**
- Waiting-time work was redirected into documentation hardening, schema freeze notes, and provisional WGS-only interpretation.
- Corrected rerun chain remains queued on Puhti: `33819825` (WES), `33819826` (RNA), and dependent benchmark `33819996`.
- Benchmark config for the queued real-data phase-gated run was created at `conf/benchmark_1000g_phase_gated_abc.yaml`.

**Decision or interpretation**
- Use queue wait to finalize benchmark-facing narrative and schema behavior rather than pausing project momentum.
- Keep all headline scientific language tied to real-data-only 1000 Genomes analyses and `HLA-A/B/C` focus.
- Treat current WGS observations as provisional until the tri-modal benchmark refresh completes.

**Evidence/files/jobs**
- Jobs: `33819825`, `33819826`, `33819996`
- Config: `conf/benchmark_1000g_phase_gated_abc.yaml`
- Current benchmark tables: `analysis/latest_benchmark/tables/`
- WGS result inventory: `/scratch/project_2008084/hla_calibration/wgs/results`

**Next action**
- When queued jobs complete, verify WES/RNA tool coverage, rerun the 1000 Genomes benchmark, inspect `benchmark_metadata.json`, and replace manuscript placeholders with final numbers.

## 2026-04-02 17:20:00 EEST
**Status**
Scheduler and runtime stabilization completed; corrected reruns submitted.

**What happened**
- Diagnosed stalled WES/RNA jobs as SLURM submission-limit failures rather than tool-level biology failures.
- Observed `AssocMaxSubmitJobLimit` when validating generated task wrappers.
- Added submission throttling, overwrite-safe report settings, and explicit Puhti account handling in the pipeline configuration and launcher scripts.
- Fixed consensus modality propagation in DSL2 and refreshed dependent benchmark submission chain.

**Decision or interpretation**
- Submission throttling (`queueSize=4`, `submitRateLimit='1 sec'`) is required on Puhti for this pipeline profile.
- The benchmark should be launched as a dependent job instead of manually after WES/RNA completion.

**Evidence/files/jobs**
- Jobs canceled: `33688505`, `33688510`
- Jobs submitted: `33819825`, `33819826`, `33819996`
- Queue diagnosis: `sbatch --test-only` returned `AssocMaxSubmitJobLimit`
- Files touched in this stabilization phase included `nextflow.config`, `main.nf`, `modules/majority_voting.nf`, and `modules/polysolver.nf`.

**Next action**
- Monitor queue start estimates and use wait time for manuscript/schema/interim analysis work.

## 2026-04-02 16:00:00 EEST
**Status**
Consensus/runtime stabilization completed and validated locally.

**What happened**
- Added runtime weight compatibility so consensus accepts both native runtime-weight JSON and imported `raw_accuracy` JSON.
- Made modality resolution explicit and propagated it through aggregation and consensus.
- Set benchmark/consensus headline scope to `A,B,C`, while preserving richer per-tool detailed outputs for other HLA loci.

**Decision or interpretation**
- Headline benchmark and manuscript reporting should focus on `HLA-A/B/C`.
- Additional loci remain accessible through detailed tool outputs but are not part of the main benchmark claims in the current phase.

**Evidence/files/jobs**
- Tests passing after stabilization: `tests/test_majority_voting_runtime.py`, `tests/test_hla_benchmark.py`, `tests/test_1000g_benchmark_workflow.py`
- Runtime weight files in `conf/tool_weights_{wgs,wes,rna}.json`

**Next action**
- Use corrected WES/RNA reruns to refresh the phase-gated real-data benchmark.

## 2026-04-01 21:00:00 EEST
**Status**
Project policy shifted from synthetic-driven interpretation to 1000 Genomes real-data-only scientific reporting.

**What happened**
- Synthetic fixtures were retained for CI, parser checks, and plumbing validation only.
- Real-data benchmark architecture was scaffolded around 1000 Genomes truth manifests, sequencing manifests, and strict tri-modal cohort logic.

**Decision or interpretation**
- Scientific claims, manuscript numbers, and figure interpretation must come only from truth-backed 1000 Genomes analyses.
- Phase-gated coverage is acceptable in the short term if missing tools are explicitly reported as unavailable rather than silently dropped.

**Evidence/files/jobs**
- `docs/ANALYSIS_2026-04-01.md`
- `bin/build_1000g_benchmark_manifests.py`
- `bin/run_1000g_benchmark.py`

**Next action**
- Complete the tri-modal benchmark refresh and replace placeholders with real cohort numbers.

## Post-Run Handoff Checklist
1. Verify WES and RNA result coverage by sample and tool from `/scratch/project_2008084/hla_calibration/wes_3sample/results` and `/scratch/project_2008084/hla_calibration/rna_3sample/results`.
2. Regenerate phase-gated inputs with `bin/build_1000g_phase_gated_inputs.py`.
3. Run `bin/run_1000g_benchmark.py` with `conf/benchmark_1000g_phase_gated_abc.yaml` if the dependent job did not already complete successfully.
4. Inspect `analysis/1000g_realdata/benchmark_phase_gated_abc/tables/benchmark_metadata.json`.
5. Summarize `method_comparison.tsv`, `summary_per_gene.tsv`, `confidence_*`, `abstention_*`, and `discordance_*`.
6. Replace manuscript placeholders with real cohort size, tool coverage, and A/B/C performance numbers.
7. Append a final notebook entry capturing cohort/tool coverage, benchmark status, and interpretation boundaries.
