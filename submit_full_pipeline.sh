#!/bin/bash
# submit_full_pipeline.sh
#
# Submits the complete pipeline:
#   1. Immediate cleanup of already-processed RNA intermediates
#   2. WES batches (10 samples at a time, 5 batches)
#   3. RNA batches (10 samples at a time, 5 batches)
#   4. Trimodal benchmark (after all batches complete)
#   5. Figure generation (after benchmark)
#
# Usage: bash submit_full_pipeline.sh

set -euo pipefail

REPO=/scratch/project_2008084/pihla-publish
cd "$REPO"

echo "[$(date)] Submitting full PIHLA pipeline"
echo ""

# ── Cleanup (immediate, no dependency) ───────────────────────────────────
CLEANUP=$(sbatch --parsable slurm_cleanup_rna_intermediates.sh)
echo "Cleanup:        job $CLEANUP"

# ── WES missing samples (49 total, 5 batches of ~10) ─────────────────────
WES_MISSING=(
    HG00096 HG00097 HG00099 HG00100 HG00101 HG00102 HG00103 HG00104 HG00105 HG00106
    HG00107 HG00108 HG00109 NA10847 NA10851 NA11830 NA11840 NA11920 NA12004 NA12044
    NA12154 NA18561 NA18562 NA18563 NA18564 NA19093 NA19099 NA19116 NA19119 NA19129
    NA19130 NA19131 NA19152 NA19153 NA19209 NA19210 NA19238 NA19239 NA20502 NA20504
    NA20505 NA20506 NA20508 NA20509 NA20510 NA20512 NA20519 NA20521 NA20528
)

echo ""
echo "Submitting WES batches (${#WES_MISSING[@]} samples):"
PREV_WES=""
for (( i=0; i<${#WES_MISSING[@]}; i+=10 )); do
    batch=("${WES_MISSING[@]:$i:10}")
    if [[ -z "$PREV_WES" ]]; then
        JOB=$(sbatch --parsable "$REPO/slurm_wes_batch10.sh" "${batch[@]}")
    else
        JOB=$(sbatch --parsable --dependency=afterok:${PREV_WES} "$REPO/slurm_wes_batch10.sh" "${batch[@]}")
    fi
    PREV_WES="$JOB"
    echo "  WES batch $((i/10+1)): job $JOB → samples ${batch[*]}"
done
LAST_WES="$PREV_WES"

# ── RNA missing samples (50 total, 5 batches of 10) ──────────────────────
RNA_MISSING=(
    HG00107 HG00551 HG00553 HG00554 HG00637 HG00638 HG00640 HG00641 HG00734 NA07037
    NA07048 NA07051 NA07056 NA07347 NA07357 NA12044 NA12154 NA12751 NA18501 NA18504
    NA18507 NA18516 NA18517 NA18522 NA18523 NA18526 NA18530 NA18532 NA18534 NA18536
    NA18537 NA18542 NA18543 NA18544 NA18545 NA18561 NA18562 NA18563 NA18564 NA18858
    NA18861 NA19131 NA19153 NA19209 NA19210 NA19238 NA19239 NA20519 NA20521 NA20528
)

echo ""
echo "Submitting RNA batches (${#RNA_MISSING[@]} samples, after cleanup):"
PREV_RNA="$CLEANUP"
for (( i=0; i<${#RNA_MISSING[@]}; i+=10 )); do
    batch=("${RNA_MISSING[@]:$i:10}")
    JOB=$(sbatch --parsable --dependency=afterok:${PREV_RNA} "$REPO/slurm_rna_batch10.sh" "${batch[@]}")
    PREV_RNA="$JOB"
    echo "  RNA batch $((i/10+1)): job $JOB → samples ${batch[*]}"
done
LAST_RNA="$PREV_RNA"

# ── Trimodal benchmark (after last WES batch AND last RNA batch) ──────────
echo ""
BENCH=$(sbatch --parsable --dependency=afterok:${LAST_WES}:${LAST_RNA} \
    "$REPO/slurm_benchmark_trimodal.sh")
echo "Benchmark:      job $BENCH (after WES job $LAST_WES and RNA job $LAST_RNA)"

# ── Figures (after benchmark) ─────────────────────────────────────────────
FIG1=$(sbatch --parsable --dependency=afterok:${BENCH} "$REPO/slurm_generate_figures.sh")
FIG2=$(sbatch --parsable --dependency=afterok:${BENCH} "$REPO/slurm_generate_new_figures.sh")
echo "Figures (v2):   job $FIG1"
echo "Figures (new):  job $FIG2"

echo ""
echo "[$(date)] All jobs submitted. Monitor with:"
echo "  squeue -u ozcanumu -o '%.10i %.22j %.8T %.10M %R'"
