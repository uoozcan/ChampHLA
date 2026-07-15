#!/bin/bash
# submit_fimm_bulkrna.sh — FIMM bulkRNA array submission.
#
# Downloads encrypted FASTQ pairs from Allas (s3allas:psergeev-2012380-MISC),
# decrypts with crypt4gh, runs PIHLA RNA pipeline (arcashla, optitype,
# seq2hla, t1k, spechla, hlahd), uploads results, and cleans local files.
#
# !! Prerequisite: list the key bucket first to confirm KEY_REMOTE_NAME
#   /appl/opt/csc-cli-utils/allas-cli-utils/rclone ls s3allas:2012380-keys/
#   Then update KEY_REMOTE_NAME in slurm_fimm_bulkrna_array.sh if needed.
#
# Prerequisite:
#   export C4GH_PASSPHRASE="<passphrase for FIMM decryption key>"
#
# Usage: bash submit_fimm_bulkrna.sh

set -euo pipefail

[[ -n "${C4GH_PASSPHRASE:-}" ]] || {
    echo "ERROR: C4GH_PASSPHRASE not set." >&2
    echo "       Run: export C4GH_PASSPHRASE='<passphrase for FIMM key>'" >&2
    exit 1
}

REPO=/scratch/project_2008084/pihla-publish
BATCH_SIZE=3
SAMPLE_LIST="$REPO/fimm_samples.txt"

[[ -f "$SAMPLE_LIST" ]] || { echo "ERROR: $SAMPLE_LIST not found" >&2; exit 1; }

mkdir -p "${REPO}/analysis/logs/fimm"

n=$(wc -l < "$SAMPLE_LIST")
last=$(( (n + BATCH_SIZE - 1) / BATCH_SIZE - 1 ))

echo "[$(date)] FIMM bulkRNA array submission"
echo "  Sample list: $SAMPLE_LIST ($n samples)"
echo "  Batch size:  $BATCH_SIZE"
echo "  Array:       0-${last} ($(( last + 1 )) tasks, %1)"
echo ""

JOB=$(sbatch --parsable \
    --export=ALL \
    --array=0-${last}%1 \
    "$REPO/slurm_fimm_bulkrna_array.sh" \
    --sample-list "$SAMPLE_LIST" \
    --batch-size $BATCH_SIZE)
echo "  bulkRNA array job: $JOB"
echo ""
echo "[$(date)] Submitted. Monitor:"
echo "  squeue -u ozcanumu -o '%.10i %.28j %.8T %.10M %R'"
echo "  tail -f ${REPO}/analysis/logs/fimm/bulkrna_arr_${JOB}_*.out"
