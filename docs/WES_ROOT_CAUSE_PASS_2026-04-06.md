# WES Root-Cause Pass (2026-04-06)

## Scope
Focused debugging pass for the failed 3-sample WES run under:
- `/scratch/project_2008084/hla_calibration/logs/wes_3sample_33828475.out`
- `/scratch/project_2008084/hla_calibration/work/wes_3sample_33828475/`

Goal:
- identify true runtime causes of WES tool failures
- separate quick code fixes from deeper environment/container problems
- apply the safest fixes first to improve the next WES rerun

## Confirmed Root Causes

### 1. `EXTRACT_HLA_AND_CONVERT` failed before `OptiType` and `T1K` could launch
Evidence from task stderr:
- `Could not retrieve index file`
- `Random alignment retrieval only works for indexed ... BAM or CRAM`

Interpretation:
- staged BAM files inside the task did not always carry a staged `.bai`
- region extraction therefore failed before the FASTQ-producing branch completed
- because `OptiType` and `T1K` depend on `ch_fastq`, they were never launched in WES

Fix applied:
- `modules/bam_to_fastq.nf`
- `EXTRACT_HLA_REGION` and `EXTRACT_HLA_AND_CONVERT` now create a BAM index inside the task if the staged `.bai` is absent

Expected impact:
- WES `OptiType` and `T1K` should now be able to launch
- WES HLA-region FASTQ generation should become stable for indexed BAM inputs

### 2. `HLAHD_BAM` failed due container-side samtools runtime issues, not missing HLA-HD itself
Evidence:
- container probe showed HLA-HD executable exists at `/app/hlahd.1.4.0/bin/hlahd.sh`
- task stderr showed:
  - `samtools: error while loading shared libraries: libcurl.so.4: cannot open shared object file`

Interpretation:
- the HLA-HD image does not ship a working `samtools`
- the current Puhti strategy of binding a host `samtools` binary into the container is not reliable because the binary's host-linked shared libraries are not resolved safely in the container runtime

Fixes applied:
- `modules/hlahd.nf`
  - both BAM and FASTQ modules now resolve `hlahd.sh` robustly with fallback to `/app/hlahd.1.4.0/bin/hlahd.sh`
- `main.nf`
  - for BAM/WES runs, `HLA-HD` is now routed through `HLAHD_FASTQ` using the same extracted FASTQ channel as `OptiType`/`T1K`
  - `HLAHD_BAM` remains available for non-WES BAM paths

Expected impact:
- WES HLA-HD should no longer depend on containerized BAM-mode `samtools`
- WES HLA-HD should become recoverable through the FASTQ path once extraction succeeds

### 3. `Kourami` failed because `samtools` is absent in the container runtime
Evidence:
- container probe: `bwa` and `java` exist, `samtools` does not
- task stderr: `samtools: command not found`

Interpretation:
- current WES Kourami failures are real runtime-environment failures, not parser issues
- naive host-binary bind attempts were not stable enough to treat as a safe production fix in this pass

Fix applied:
- `main.nf`
- WES `Kourami` is now explicitly phase-gated with a warning on Puhti instead of failing at runtime

Expected impact:
- cleaner WES runs with honest reporting
- avoids repeated known-bad Kourami failures in WES until a stable samtools source is available in-container or via a validated wrapper

### 4. `POLYSOLVER` failed because the bundled samtools uses older sort syntax
Evidence from task stderr:
- `sort: invalid option -- '@'`
- `fail to open file 8`
- subsequent segfault in `fixmate`

Interpretation:
- the POLYSOLVER container bundles an older samtools where `sort` does not support the newer `-@` / `-o` syntax used in the module
- the broken `sort` command produced invalid intermediate input, and the later `fixmate` segfault was a downstream symptom

Fix applied:
- `modules/polysolver.nf`
- preprocessing now uses the older samtools sort form:
  - `samtools sort -n <in.bam> <out.prefix>`
  - `samtools sort <in.bam> <out.prefix>`

Expected impact:
- POLYSOLVER preprocessing should now be compatible with the container's bundled samtools
- the prior exit-139 failure mode should be reduced or eliminated if no second tool-specific issue remains

## Validation Performed
- `python3 -m unittest tests.test_majority_voting_runtime tests.test_hla_benchmark tests.test_1000g_benchmark_workflow`
  - passed (`Ran 16 tests ... OK`)
- `module load nextflow && nextflow run main.nf --help -profile puhti,singularity`
  - passed

## Current WES Status After This Pass
Fixed now:
- WES FASTQ extraction index handling
- WES HLA-HD routing strategy
- POLYSOLVER samtools syntax mismatch
- explicit Kourami WES phase-gating on Puhti

Still remaining for later WES follow-up:
- verify next WES rerun produces real `OptiType`, `T1K`, and `HLA-HD` outputs
- inspect whether POLYSOLVER now completes end-to-end after the sort syntax fix
- revisit Kourami-on-WES only after a stable samtools runtime strategy exists
- separately debug `ArcasHLA` WES exit `1`

## Recommended Next Action
- rerun the 3-sample WES workflow with the patched code
- rebuild the WES failure matrix from the new run
- then decide whether WES is ready for the next multimodal benchmark wave or needs one more focused tool pass
