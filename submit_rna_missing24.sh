#!/bin/bash
# submit_rna_missing24.sh — RNA typing for 24 trimodal-eligible samples missing RNA.
#
# Submits a SLURM array job (3 tasks, one task at a time) where each task
# processes one pre-generated batch file of samples.  After the array
# finishes the trimodal benchmark is submitted automatically.
#
# Usage: bash submit_rna_missing24.sh

set -euo pipefail
REPO=/scratch/project_2008084/pihla-publish
BATCH_DIR="$REPO/analysis/rna_missing24_batches"
mkdir -p "$BATCH_DIR"

# ── 24 samples with WGS+WES done, RNA missing ────────────────────────────
RNA_MISSING24=(
    HG00551 HG00553 HG00554 HG00637 HG00638 HG00640 HG00641 HG00734
    NA18501 NA18504 NA18507 NA18516 NA18522 NA18523 NA18526 NA18530
    NA18532 NA18534 NA18536 NA18537 NA18542 NA18543 NA18544 NA18545
)

# ── Write batch files (10 / 10 / 4) ─────────────────────────────────────
echo "[$(date)] Writing batch files to $BATCH_DIR"
BATCH_SIZE=10
BATCH_IDX=0
for (( i=0; i<${#RNA_MISSING24[@]}; i+=BATCH_SIZE )); do
    slice=("${RNA_MISSING24[@]:$i:$BATCH_SIZE}")
    printf '%s\n' "${slice[@]}" > "$BATCH_DIR/batch_${BATCH_IDX}.txt"
    echo "  batch_${BATCH_IDX}.txt: ${#slice[@]} samples"
    BATCH_IDX=$(( BATCH_IDX + 1 ))
done
N_BATCHES=$BATCH_IDX
LAST_IDX=$(( N_BATCHES - 1 ))

# ── Submit array job (one task at a time to avoid disk pressure) ──────────
echo ""
echo "[$(date)] Submitting RNA array job (0-${LAST_IDX}%1)..."
ARRAY_JOB=$(sbatch --parsable \
    --array=0-${LAST_IDX}%1 \
    "$REPO/slurm_rna_array.sh" "$BATCH_DIR")
echo "  RNA array job: $ARRAY_JOB"

# ── Trimodal benchmark after all array tasks finish ───────────────────────
BENCH=$(sbatch --parsable \
    --dependency=afterok:${ARRAY_JOB} \
    "$REPO/slurm_benchmark_trimodal.sh")
echo "  Benchmark job: $BENCH"

echo ""
echo "[$(date)] Done. Monitor:"
echo "  squeue -u ozcanumu -o '%.10i %.22j %.8T %.10M %R'"
