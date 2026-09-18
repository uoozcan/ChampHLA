#!/bin/bash
# submit_fimm_scrna.sh — FIMM scRNA array submission.
#
# !! Before running this script:
#
#   1. Find the scRNA files in Allas:
#      RCLONE=/appl/opt/csc-cli-utils/allas-cli-utils/rclone
#      $RCLONE ls s3allas:psergeev-2012380-MISC/ | grep -i "scrna\|single\|10x\|cellranger"
#
#   2. Update S3_BUCKET in this script and in slurm_fimm_scrna_array.sh.
#
#   3. Check read length of R1 to choose mode:
#      - Standard paired FASTQs (R1 = full-length cDNA): MODE=standard
#      - 10X Chromium (R1 = 28bp barcode+UMI):           MODE=tenx
#
#   4. Confirm KEY_REMOTE_NAME after: rclone ls s3allas:2012380-keys/
#      Update KEY_REMOTE_NAME in slurm_fimm_scrna_array.sh if needed.
#
# Prerequisite:
#   export C4GH_PASSPHRASE="<passphrase for FIMM decryption key>"
#
# Usage: bash submit_fimm_scrna.sh [standard|tenx]

set -euo pipefail

[[ -n "${C4GH_PASSPHRASE:-}" ]] || {
    echo "ERROR: C4GH_PASSPHRASE not set." >&2
    echo "       Run: export C4GH_PASSPHRASE='<passphrase for FIMM key>'" >&2
    exit 1
}

REPO=/scratch/project_2008084/pihla-publish
BATCH_SIZE=3
SAMPLE_LIST="$REPO/fimm_samples.txt"
MODE="${1:-standard}"

# !! UPDATE to the correct scRNA path after discovery
S3_BUCKET="s3allas:psergeev-2012380-MISC/fimm_ga2_heckman"

[[ -f "$SAMPLE_LIST" ]] || { echo "ERROR: $SAMPLE_LIST not found" >&2; exit 1; }
[[ "$MODE" == "standard" || "$MODE" == "tenx" ]] || {
    echo "ERROR: mode must be 'standard' or 'tenx'" >&2
    echo "Usage: $0 [standard|tenx]" >&2
    exit 1
}

mkdir -p "${REPO}/analysis/logs/fimm"

n=$(wc -l < "$SAMPLE_LIST")
last=$(( (n + BATCH_SIZE - 1) / BATCH_SIZE - 1 ))

echo "[$(date)] FIMM scRNA array submission (mode: $MODE)"
echo "  Sample list: $SAMPLE_LIST ($n samples)"
echo "  Batch size:  $BATCH_SIZE"
echo "  S3 bucket:   $S3_BUCKET"
echo "  Array:       0-${last} ($(( last + 1 )) tasks, %1)"
echo ""

JOB=$(sbatch --parsable \
    --export=ALL \
    --array=0-${last}%1 \
    "$REPO/slurm_fimm_scrna_array.sh" \
    --sample-list "$SAMPLE_LIST" \
    --batch-size $BATCH_SIZE \
    --mode "$MODE" \
    --s3bucket "$S3_BUCKET")
echo "  scRNA array job: $JOB (mode: $MODE)"
echo ""
echo "[$(date)] Submitted. Monitor:"
echo "  squeue -u ozcanumu -o '%.10i %.28j %.8T %.10M %R'"
echo "  tail -f ${REPO}/analysis/logs/fimm/scrna_arr_${JOB}_*.out"
