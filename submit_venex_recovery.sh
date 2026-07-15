#!/bin/bash
# Recovery submission for VENEX samples with incomplete results (<30 result files).
# Uses --batch-size 2 to limit peak concurrent disk usage (~120 GB per task).
#
# Prerequisite:
#   export C4GH_PASSPHRASE="<passphrase for all_data_csc_key.sec>"
#
# Usage: bash submit_venex_recovery.sh
set -euo pipefail

[[ -n "${C4GH_PASSPHRASE:-}" ]] || {
    echo "ERROR: C4GH_PASSPHRASE not set." >&2
    echo "       Run: export C4GH_PASSPHRASE='<passphrase>'" >&2
    exit 1
}

REPO=/scratch/project_2008084/pihla-publish
BATCH_SIZE=2

echo "[$(date)] Generating recovery sample lists..."
bash "${REPO}/create_venex_recovery_lists.sh"

B1_LIST="${REPO}/venex_batch1_recovery.txt"
B12_LIST="${REPO}/venex_batch1_2_recovery.txt"

n_b1=$(wc -l < "$B1_LIST")
n_b12=$(wc -l < "$B12_LIST")

echo ""
echo "[$(date)] VENEX recovery submission"
echo "  Batch1:   $n_b1 samples, batch-size $BATCH_SIZE"
echo "  Batch1_2: $n_b12 samples, batch-size $BATCH_SIZE"
echo ""

if [[ "$n_b1" -gt 0 ]]; then
    last_b1=$(( (n_b1 + BATCH_SIZE - 1) / BATCH_SIZE - 1 ))
    JOB_B1=$(sbatch --parsable --export=ALL \
        --array=0-${last_b1}%1 \
        "${REPO}/slurm_venex_array.sh" \
        --s3bucket "s3allas:2014061-DNA_seq/DNA_seq/DRAGEN/VenEx_DNA_NONHUS_Batch1" \
        --sample-list "$B1_LIST" \
        --modality wgs \
        --batch-size $BATCH_SIZE)
    echo "  Batch1 recovery array job: $JOB_B1"
else
    echo "  Batch1: no samples need recovery"
fi

if [[ "$n_b12" -gt 0 ]]; then
    last_b12=$(( (n_b12 + BATCH_SIZE - 1) / BATCH_SIZE - 1 ))
    JOB_B12=$(sbatch --parsable --export=ALL \
        --array=0-${last_b12}%1 \
        "${REPO}/slurm_venex_array.sh" \
        --s3bucket "s3allas:2014061-DNA_seq/DNA_seq/DRAGEN/VenEx_DNA_NONHUS_Batch1_2" \
        --sample-list "$B12_LIST" \
        --modality wgs \
        --batch-size $BATCH_SIZE)
    echo "  Batch1_2 recovery array job: $JOB_B12"
else
    echo "  Batch1_2: no samples need recovery"
fi

echo ""
echo "[$(date)] Done. Monitor:"
echo "  squeue -u ozcanumu -o '%.14i %.22j %.8T %.10M %R'"
