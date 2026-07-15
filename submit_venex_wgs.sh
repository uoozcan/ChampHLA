#!/bin/bash
# submit_venex_wgs.sh — VENEX WGS array submission for both Allas batches.
#
# Submits two independent SLURM array jobs (one per Allas bucket), each
# running one task at a time to control disk usage.  Each array task
# processes BATCH_SIZE=5 samples: download DLO segments via S3 (no Swift
# auth needed) → decrypt crypt4gh → type (PIHLA) → upload results → clean.
#
# Prerequisite:
#   export C4GH_PASSPHRASE="<passphrase for all_data_csc_key.sec>"
#
# Usage: bash submit_venex_wgs.sh

set -euo pipefail

[[ -n "${C4GH_PASSPHRASE:-}" ]] || {
    echo "ERROR: C4GH_PASSPHRASE not set." >&2
    echo "       Run: export C4GH_PASSPHRASE='<passphrase for all_data_csc_key.sec>'" >&2
    exit 1
}

REPO=/scratch/project_2008084/pihla-publish
BATCH_SIZE=5

B1_LIST="$REPO/venex_batch1_samples.txt"
B12_LIST="$REPO/venex_batch1_2_samples.txt"

for f in "$B1_LIST" "$B12_LIST"; do
    [[ -f "$f" ]] || { echo "ERROR: $f not found" >&2; exit 1; }
done

# ── Compute array bounds ──────────────────────────────────────────────────
n_b1=$(wc -l < "$B1_LIST")
n_b12=$(wc -l < "$B12_LIST")
last_b1=$(( (n_b1 + BATCH_SIZE - 1) / BATCH_SIZE - 1 ))
last_b12=$(( (n_b12 + BATCH_SIZE - 1) / BATCH_SIZE - 1 ))

echo "[$(date)] VENEX WGS array submission"
echo "  Batch1:   $n_b1 samples → array 0-${last_b1} ($(( last_b1 + 1 )) tasks, %1)"
echo "  Batch1_2: $n_b12 samples → array 0-${last_b12} ($(( last_b12 + 1 )) tasks, %1)"
echo ""

# ── Batch1 (VenEx_DNA_NONHUS_Batch1) ─────────────────────────────────────
JOB_B1=$(sbatch --parsable \
    --export=ALL \
    --array=0-${last_b1}%1 \
    "$REPO/slurm_venex_array.sh" \
    --s3bucket "s3allas:2014061-DNA_seq/DNA_seq/DRAGEN/VenEx_DNA_NONHUS_Batch1" \
    --sample-list "$B1_LIST" \
    --modality wgs \
    --batch-size $BATCH_SIZE)
echo "  Batch1 array job:   $JOB_B1"

# ── Batch1_2 (VenEx_DNA_NONHUS_Batch1_2) — runs in parallel ──────────────
JOB_B12=$(sbatch --parsable \
    --export=ALL \
    --array=0-${last_b12}%1 \
    "$REPO/slurm_venex_array.sh" \
    --s3bucket "s3allas:2014061-DNA_seq/DNA_seq/DRAGEN/VenEx_DNA_NONHUS_Batch1_2" \
    --sample-list "$B12_LIST" \
    --modality wgs \
    --batch-size $BATCH_SIZE)
echo "  Batch1_2 array job: $JOB_B12"

echo ""
echo "[$(date)] Done. Monitor:"
echo "  squeue -u ozcanumu -o '%.10i %.22j %.8T %.10M %R'"
