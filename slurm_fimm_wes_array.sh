#!/bin/bash
#SBATCH --job-name=pihla_fimm_wes
#SBATCH --account=project_2008084
#SBATCH --partition=small
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --output=/scratch/project_2008084/pihla-publish/analysis/logs/fimm/wes_arr_%A_%a.out
#SBATCH --error=/scratch/project_2008084/pihla-publish/analysis/logs/fimm/wes_arr_%A_%a.err

# FIMM WES array task — each task processes one batch of samples from Allas.
# SLURM_ARRAY_TASK_ID * BATCH_SIZE gives the offset into the sample list.
#
# Allas paths:
#   WES files: s3allas:psergeev-2012380-MISC/fimm_ga2_heckman/vcp/fimm_ga2_heckman/
#   Keys:      s3allas:2012380-keys/
#
# !! Before running: list the keys bucket to find the key filename:
#   /appl/opt/csc-cli-utils/allas-cli-utils/rclone ls s3allas:2012380-keys/
#   Then update KEY_REMOTE_NAME below if it differs from the default.
#
# Files are crypt4gh-encrypted BAMs (.bam.c4gh). The script handles both
# Swift DLO objects and regular S3 objects.
#
# Prerequisites (set before submitting via submit_fimm_wes.sh):
#   export C4GH_PASSPHRASE="<passphrase for FIMM key>"
#
# Usage (do not call directly — use submit_fimm_wes.sh):
#   sbatch --array=0-N%1 --export=ALL slurm_fimm_wes_array.sh \
#     --sample-list <file> [--batch-size 3]

set -euo pipefail
trap 'echo "ERR at line $LINENO: $BASH_COMMAND" >&2' ERR

# ── Parse arguments ───────────────────────────────────────────────────────
SAMPLE_LIST=""
BATCH_SIZE=3

while [[ $# -gt 0 ]]; do
    case "$1" in
        --sample-list) SAMPLE_LIST="$2"; shift 2 ;;
        --batch-size)  BATCH_SIZE="$2";  shift 2 ;;
        *) echo "Unknown arg: $1" >&2; exit 1 ;;
    esac
done

[[ -z "$SAMPLE_LIST" ]] && { echo "ERROR: --sample-list required" >&2; exit 1; }
[[ -z "${C4GH_PASSPHRASE:-}" ]] && { echo "ERROR: C4GH_PASSPHRASE not set" >&2; exit 1; }

OFFSET=$(( SLURM_ARRAY_TASK_ID * BATCH_SIZE ))

# ── Paths ─────────────────────────────────────────────────────────────────
REPO=/scratch/project_2008084/pihla-publish
SCRATCH=/scratch/project_2008084/hla_calibration
LOCAL_INPUT=${SCRATCH}/fimm/inputs/wes
LOCAL_RESULTS=${SCRATCH}/fimm/results/wes
WORK_DIR=${SCRATCH}/fimm/work/wes
LOG_DIR=${REPO}/analysis/logs/fimm

S3_BUCKET="s3allas:psergeev-2012380-MISC/fimm_ga2_heckman/vcp/fimm_ga2_heckman"
# Segments bucket derived from main bucket (for DLO fallback)
S3_SEGS_BUCKET="${S3_BUCKET/psergeev-2012380-MISC/psergeev-2012380-MISC_segments}"
S3_KEYS="s3allas:2012380-keys"

# !! UPDATE THIS after running: rclone ls s3allas:2012380-keys/
KEY_REMOTE_NAME="all_data_csc_key.sec"

RCLONE=/appl/opt/csc-cli-utils/allas-cli-utils/rclone
CRYPT4GH=/appl/soft/bio/biopython/gcc_11.3.0/3.10.6_sqlite/bin/python3
KEY_FILE="${TMPDIR:-/tmp}/fimm_key.sec"

mkdir -p "$LOCAL_INPUT" "$LOCAL_RESULTS" "$WORK_DIR" "$LOG_DIR"

module load nextflow
module load samtools

# ── Fetch decryption key from S3 ─────────────────────────────────────────
echo "[$(date)] Fetching decryption key (${KEY_REMOTE_NAME})..."
"$RCLONE" copyto "${S3_KEYS}/${KEY_REMOTE_NAME}" "$KEY_FILE"
if [[ ! -f "$KEY_FILE" ]]; then
    "$RCLONE" copy "${S3_KEYS}/${KEY_REMOTE_NAME}" "$(dirname "$KEY_FILE")/"
    downloaded_key="$(dirname "$KEY_FILE")/${KEY_REMOTE_NAME}"
    [[ -f "$downloaded_key" ]] && mv "$downloaded_key" "$KEY_FILE"
fi
[[ -f "$KEY_FILE" ]] || { echo "ERROR: failed to download decryption key" >&2; exit 1; }
chmod 600 "$KEY_FILE"

# ── Select this task's batch ──────────────────────────────────────────────
mapfile -t ALL_SAMPLES < "$SAMPLE_LIST"
TOTAL=${#ALL_SAMPLES[@]}
batch=("${ALL_SAMPLES[@]:$OFFSET:$BATCH_SIZE}")

echo "[$(date)] FIMM WES array task ${SLURM_ARRAY_TASK_ID}"
echo "  S3 bucket: $S3_BUCKET"
echo "  Offset:    $OFFSET / $TOTAL"
echo "  Samples:   ${batch[*]}"
echo ""

[[ ${#batch[@]} -eq 0 ]] && { echo "No samples for this task — exiting"; exit 0; }

# ── Helper: download a DLO file via S3 segment concatenation ─────────────
download_dlo() {
    local seg_prefix="$1"
    local out_file="$2"
    local tmpdir="${out_file}.segtmp"
    mkdir -p "$tmpdir"

    local seg_id
    seg_id=$("$RCLONE" lsf "${seg_prefix}/" 2>/dev/null | head -1 | tr -d '/') || seg_id=""
    if [[ -z "$seg_id" ]]; then
        rm -rf "$tmpdir"
        echo "  ERROR: no segments found under ${seg_prefix}" >&2
        return 1
    fi

    echo "  Downloading DLO segments (id: $seg_id)..."
    "$RCLONE" copy "${seg_prefix}/${seg_id}/" "$tmpdir/" --transfers 8
    echo "  Concatenating segments..."
    find "$tmpdir" -maxdepth 1 -type f | sort -V | xargs cat > "$out_file"
    rm -rf "$tmpdir"
}

# ── Helper: download a file (regular S3 or DLO fallback) ─────────────────
download_file() {
    local s3_path="$1"
    local out_file="$2"
    local size

    "$RCLONE" copyto "$s3_path" "$out_file" 2>/dev/null && \
        size=$(stat -c%s "$out_file" 2>/dev/null || echo 0) && \
        [[ "$size" -gt 0 ]] && return 0

    echo "  Direct copy produced empty file — trying DLO segment download..."
    local seg_prefix="${s3_path/$S3_BUCKET/$S3_SEGS_BUCKET}"
    download_dlo "$seg_prefix" "$out_file"
}

# ── Step 1: Download + decrypt each sample's BAM ─────────────────────────
echo "[$(date)] Step 1: downloading and decrypting BAM files from Allas..."
downloaded=()

for sample in "${batch[@]}"; do
    n_done=$(find "${LOCAL_RESULTS}/${sample}/results/${sample}" \
        \( -name "*.txt" -o -name "*.tsv" \) 2>/dev/null | wc -l) || n_done=0
    if [[ "$n_done" -ge 10 ]]; then
        echo "  SKIP $sample (already typed: $n_done result files)"
        continue
    fi

    # Discover BAM filename for this sample
    echo "  Listing Allas for $sample..."
    sample_files=$("$RCLONE" lsf "${S3_BUCKET}/${sample}/" 2>/dev/null) || sample_files=""
    if [[ -z "$sample_files" ]]; then
        # Try flat listing (no per-sample subdirectory)
        sample_files=$("$RCLONE" lsf "${S3_BUCKET}/" 2>/dev/null | grep "^${sample}" || true)
    fi
    if [[ -z "$sample_files" ]]; then
        echo "  WARNING: $sample — no files found at ${S3_BUCKET}/${sample}/" >&2
        continue
    fi

    bam_enc_name=$(echo "$sample_files" | grep "\.bam\.c4gh$" | head -1 || true)
    if [[ -z "$bam_enc_name" ]]; then
        echo "  WARNING: $sample — no .bam.c4gh file found in listing" >&2
        echo "  Found: $sample_files" >&2
        continue
    fi

    # Strip leading path components to get just the filename
    bam_basename=$(basename "$bam_enc_name")
    bam_stem="${bam_basename%.c4gh}"
    bam_dec="${LOCAL_INPUT}/${bam_stem}"
    bam_enc="${LOCAL_INPUT}/${bam_basename}"

    if [[ -f "$bam_dec" ]]; then
        echo "  SKIP $sample (decrypted BAM already local)"
        downloaded+=("$sample"); continue
    fi

    # S3 path (with or without per-sample subdirectory)
    if [[ "$bam_enc_name" == */* ]]; then
        s3_path="${S3_BUCKET}/${bam_enc_name}"
    else
        s3_path="${S3_BUCKET}/${sample}/${bam_enc_name}"
    fi

    echo "  Downloading $sample BAM (${bam_enc_name})..."
    if ! download_file "$s3_path" "$bam_enc"; then
        echo "  ERROR: $sample BAM download failed" >&2
        rm -f "$bam_enc"; continue
    fi
    if [[ ! -s "$bam_enc" ]]; then
        echo "  ERROR: $sample — empty file after download" >&2
        rm -f "$bam_enc"; continue
    fi
    echo "  Downloaded $(du -sh "$bam_enc" | cut -f1) encrypted"

    echo "  Decrypting $sample..."
    C4GH_PASSPHRASE="$C4GH_PASSPHRASE" \
        "$CRYPT4GH" -m crypt4gh decrypt --sk "$KEY_FILE" \
        < "$bam_enc" > "$bam_dec"
    rm -f "$bam_enc"
    if [[ ! -s "$bam_dec" ]]; then
        echo "  ERROR: $sample — empty file after decryption (wrong passphrase?)" >&2
        rm -f "$bam_dec"; continue
    fi
    echo "  Decrypted $(du -sh "$bam_dec" | cut -f1)"

    echo "  Indexing $sample BAM..."
    samtools index -@ 4 "$bam_dec"

    downloaded+=("$sample")
    echo "  $sample ready"
done

rm -f "$KEY_FILE"

if [[ ${#downloaded[@]} -eq 0 ]]; then
    echo "WARNING: no samples ready for this task — check Allas paths and C4GH_PASSPHRASE"
    exit 1
fi

# ── Step 2: Submit PIHLA WES pipeline jobs ───────────────────────────────
echo "[$(date)] Step 2: submitting PIHLA WES jobs..."
SAMPLE_SCRIPT="$REPO/slurm_wes_sample.sh"
job_id_list=""
for sample in "${downloaded[@]}"; do
    job_id=$(sbatch --parsable "$SAMPLE_SCRIPT" \
        "$sample" "$LOCAL_INPUT" "$LOCAL_RESULTS" "$WORK_DIR")
    job_id_list="${job_id_list} ${job_id}"
    echo "  Submitted $sample → job $job_id"
done

# ── Step 3: Wait for pipeline jobs ───────────────────────────────────────
if [[ -n "$job_id_list" ]]; then
    ids=$(echo "$job_id_list" | tr ' ' ',' | sed 's/^,//')
    echo "[$(date)] Step 3: waiting for jobs ($ids)"
    while squeue --jobs="$ids" -h 2>/dev/null | grep -q .; do
        sleep 60
    done
    echo "[$(date)] Pipeline jobs finished"
fi

# ── Step 4: Upload results to Allas ──────────────────────────────────────
ALLAS_RESULTS="s3allas:psergeev-2012380-MISC-results/fimm/wes"
echo "[$(date)] Step 4: uploading results..."
for sample in "${downloaded[@]}"; do
    result_dir="${LOCAL_RESULTS}/${sample}/results/${sample}"
    [[ ! -d "$result_dir" ]] && continue
    n_results=$(find "$result_dir" \( -name "*.txt" -o -name "*.tsv" \) 2>/dev/null | wc -l)
    if [[ "$n_results" -gt 0 ]]; then
        echo "  Uploading $sample ($n_results files)..."
        "$RCLONE" copy "$result_dir" "${ALLAS_RESULTS}/${sample}/" \
            --include "*.txt" --include "*.tsv"
        echo "  Uploaded $sample"
    else
        echo "  WARNING: $sample has no results — skipping upload"
    fi
done

# ── Step 5: Clean local BAM files ────────────────────────────────────────
echo "[$(date)] Step 5: cleaning local BAM files..."
for sample in "${downloaded[@]}"; do
    result_dir="${LOCAL_RESULTS}/${sample}/results/${sample}"
    n_results=$(find "$result_dir" \( -name "*.txt" -o -name "*.tsv" \) 2>/dev/null | wc -l)
    if [[ "$n_results" -gt 0 ]]; then
        find "$LOCAL_INPUT" -maxdepth 1 \
            \( -name "${sample}*.bam" -o -name "${sample}*.bai" \) \
            -delete 2>/dev/null || true
        find "${LOCAL_RESULTS}/${sample}" \
            \( -name "*.bam" -o -name "*.fq.gz" -o -name "*.fastq.gz" \
               -o -name "*.alignment.p" \) -delete 2>/dev/null || true
        rm -rf "${WORK_DIR}/wes_${sample}" 2>/dev/null || true
        echo "  Cleaned $sample"
    else
        echo "  KEPT $sample inputs (no results — investigate)"
    fi
done

echo "[$(date)] Task ${SLURM_ARRAY_TASK_ID} complete. Disk:"
df -h /scratch/project_2008084/ | tail -1
