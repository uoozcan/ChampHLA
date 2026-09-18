#!/bin/bash
# Submit VENEX FASTQ recovery jobs for samples that lack FASTQ-pipeline results.
# Downloads R1+R2 FASTQ pairs from Allas (~14 GB/sample) and runs:
#   optitype, t1k, hlahd, spechla  (polysolver is BAM-only, skipped)
#
# Results land in:
#   /scratch/project_2008084/hla_calibration/venex/results/wgs_fastq/
# and are uploaded to:
#   s3allas:2014061-DNA_seq_results/fastq/{sample}/
#
# Prerequisite:
#   export C4GH_PASSPHRASE="<passphrase for all_data_csc_key.sec>"
#
# Usage: bash submit_venex_fastq_recovery.sh
set -euo pipefail

[[ -n "${C4GH_PASSPHRASE:-}" ]] || {
    echo "ERROR: C4GH_PASSPHRASE not set." >&2
    echo "       Run: export C4GH_PASSPHRASE='<passphrase>'" >&2
    exit 1
}

REPO=/scratch/project_2008084/pihla-publish
BATCH_SIZE=2

# Regenerate the BAM-pipeline recovery lists (samples with <30 BAM results)
# so we target the same set of incomplete samples.
echo "[$(date)] Generating recovery sample lists..."
bash "${REPO}/create_venex_recovery_lists.sh"

B1_LIST="${REPO}/venex_batch1_recovery.txt"
B12_LIST="${REPO}/venex_batch1_2_recovery.txt"

n_b1=$(grep -c . "$B1_LIST" 2>/dev/null || echo 0)
n_b12=$(grep -c . "$B12_LIST" 2>/dev/null || echo 0)

echo ""
echo "[$(date)] VENEX FASTQ recovery submission"
echo "  Batch1:   $n_b1 samples, batch-size $BATCH_SIZE"
echo "  Batch1_2: $n_b12 samples, batch-size $BATCH_SIZE"
echo ""

if [[ "$n_b1" -gt 0 ]]; then
    last_b1=$(( (n_b1 + BATCH_SIZE - 1) / BATCH_SIZE - 1 ))
    JOB_B1=$(sbatch --parsable --export=ALL \
        --array=0-${last_b1}%1 \
        "${REPO}/slurm_venex_fastq_array.sh" \
        --s3bucket "s3allas:2014061-DNA_seq/DNA_seq/DRAGEN/VenEx_DNA_NONHUS_Batch1" \
        --sample-list "$B1_LIST" \
        --batch-size $BATCH_SIZE)
    echo "  Batch1 FASTQ recovery array job: $JOB_B1"
else
    echo "  Batch1: no samples need recovery"
fi

if [[ "$n_b12" -gt 0 ]]; then
    last_b12=$(( (n_b12 + BATCH_SIZE - 1) / BATCH_SIZE - 1 ))
    JOB_B12=$(sbatch --parsable --export=ALL \
        --array=0-${last_b12}%1 \
        "${REPO}/slurm_venex_fastq_array.sh" \
        --s3bucket "s3allas:2014061-DNA_seq/DNA_seq/DRAGEN/VenEx_DNA_NONHUS_Batch1_2" \
        --sample-list "$B12_LIST" \
        --batch-size $BATCH_SIZE)
    echo "  Batch1_2 FASTQ recovery array job: $JOB_B12"
else
    echo "  Batch1_2: no samples need recovery"
fi

echo ""
echo "[$(date)] Done. Monitor:"
echo "  squeue -u ozcanumu -o '%.14i %.22j %.8T %.10M %R'"
echo ""
echo "Results will appear in:"
echo "  /scratch/project_2008084/hla_calibration/venex/results/wgs_fastq/"
