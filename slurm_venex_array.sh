#!/bin/bash
#SBATCH --job-name=pihla_venex_arr
#SBATCH --account=project_2008084
#SBATCH --partition=small
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --output=/scratch/project_2008084/pihla-publish/analysis/logs/venex_arr_%A_%a.out
#SBATCH --error=/scratch/project_2008084/pihla-publish/analysis/logs/venex_arr_%A_%a.err

# VENEX WGS array task — each task processes one batch of samples from Allas.
# SLURM_ARRAY_TASK_ID * BATCH_SIZE gives the offset into the sample list.
#
# Data format in Allas:
#   Each sample lives in a per-sample subdirectory; all files are crypt4gh-
#   encrypted (.c4gh). Large files (BAM, FASTQ) are stored as Swift DLOs,
#   meaning S3 sees them as 0-byte manifest objects. This script downloads
#   the DLO segment files directly from the _segments S3 bucket and cats
#   them together — no OS_AUTH_TOKEN or Swift auth required.
#
# Prerequisites (set before submitting via submit_venex_wgs.sh):
#   export C4GH_PASSPHRASE="<passphrase for all_data_csc_key.sec>"
#
# Usage (do not call directly — use submit_venex_wgs.sh):
#   sbatch --array=0-N%1 --export=ALL slurm_venex_array.sh \
#     --s3bucket s3allas:2014061-DNA_seq/DNA_seq/DRAGEN/VenEx_DNA_NONHUS_Batch1 \
#     --sample-list <file>    \
#     [--modality wgs]        \
#     [--batch-size 5]

set -euo pipefail
trap 'echo "ERR at line $LINENO: $BASH_COMMAND" >&2' ERR

# ── Parse arguments ───────────────────────────────────────────────────────
S3_BUCKET=""
SAMPLE_LIST=""
MODALITY="wgs"
BATCH_SIZE=5

while [[ $# -gt 0 ]]; do
    case "$1" in
        --s3bucket)    S3_BUCKET="$2";  shift 2 ;;
        --sample-list) SAMPLE_LIST="$2"; shift 2 ;;
        --modality)    MODALITY="$2";    shift 2 ;;
        --batch-size)  BATCH_SIZE="$2";  shift 2 ;;
        *) echo "Unknown arg: $1" >&2; exit 1 ;;
    esac
done

[[ -z "$S3_BUCKET" ]]    && { echo "ERROR: --s3bucket required" >&2; exit 1; }
[[ -z "$SAMPLE_LIST" ]]  && { echo "ERROR: --sample-list required" >&2; exit 1; }
[[ -z "${C4GH_PASSPHRASE:-}" ]] && { echo "ERROR: C4GH_PASSPHRASE not set" >&2; exit 1; }

OFFSET=$(( SLURM_ARRAY_TASK_ID * BATCH_SIZE ))

# ── Paths ─────────────────────────────────────────────────────────────────
REPO=/scratch/project_2008084/pihla-publish
SCRATCH=/scratch/project_2008084/hla_calibration
LOCAL_INPUT=${SCRATCH}/venex/inputs/${MODALITY}
LOCAL_RESULTS=${SCRATCH}/venex/results/${MODALITY}
WORK_DIR=${SCRATCH}/venex/work/${MODALITY}
LOG_DIR=${REPO}/analysis/logs/venex

# Derive the _segments bucket from the main bucket name
# s3allas:2014061-DNA_seq/path  →  s3allas:2014061-DNA_seq_segments/path
S3_SEGS_BUCKET="${S3_BUCKET/DNA_seq\//DNA_seq_segments/}"

mkdir -p "$LOCAL_INPUT" "$LOCAL_RESULTS" "$WORK_DIR" "$LOG_DIR"

RCLONE=/appl/opt/csc-cli-utils/allas-cli-utils/rclone
CRYPT4GH=/appl/soft/bio/biopython/gcc_11.3.0/3.10.6_sqlite/bin/python3
KEY_FILE="${TMPDIR:-/tmp}/all_data_csc_key.sec"

module load nextflow
module load samtools

# ── Fetch decryption key from S3 (no Swift auth needed) ──────────────────
echo "[$(date)] Fetching decryption key..."
"$RCLONE" copy s3allas:2014061-keys/all_data_csc_key.sec "$(dirname "$KEY_FILE")/"
[[ -f "$KEY_FILE" ]] || { echo "ERROR: failed to download decryption key" >&2; exit 1; }
chmod 600 "$KEY_FILE"

# ── Select this task's batch ──────────────────────────────────────────────
mapfile -t ALL_SAMPLES < "$SAMPLE_LIST"
TOTAL=${#ALL_SAMPLES[@]}
batch=("${ALL_SAMPLES[@]:$OFFSET:$BATCH_SIZE}")

echo "[$(date)] VENEX array task ${SLURM_ARRAY_TASK_ID}"
echo "  S3 bucket: $S3_BUCKET"
echo "  Modality:  $MODALITY"
echo "  Offset:    $OFFSET / $TOTAL"
echo "  Samples:   ${batch[*]}"
echo ""

[[ ${#batch[@]} -eq 0 ]] && { echo "No samples for this task — exiting"; exit 0; }

# ── Helper: download a DLO file via S3 segment concatenation ─────────────
# Usage: download_dlo <s3_segs_bucket/prefix> <local_output_file>
download_dlo() {
    local seg_prefix="$1"
    local out_file="$2"
    local tmpdir
    tmpdir="${out_file}.segtmp"
    mkdir -p "$tmpdir"

    # Find the segment batch ID (one subdirectory under the DLO prefix)
    local seg_id
    seg_id=$("$RCLONE" lsf "${seg_prefix}/" 2>/dev/null | head -1 | tr -d '/') || seg_id=""
    if [[ -z "$seg_id" ]]; then
        rm -rf "$tmpdir"
        echo "  ERROR: no segments found under ${seg_prefix}" >&2
        return 1
    fi

    echo "  Downloading segments (id: $seg_id)..."
    "$RCLONE" copy "${seg_prefix}/${seg_id}/" "$tmpdir/" --transfers 8

    echo "  Concatenating segments..."
    # Sort numerically and concatenate
    find "$tmpdir" -maxdepth 1 -type f | sort -V | xargs cat > "$out_file"
    rm -rf "$tmpdir"
}

# ── Step 1: Download + decrypt each sample's BAM ─────────────────────────
echo "[$(date)] Step 1: downloading and decrypting from Allas..."
downloaded=()

for sample in "${batch[@]}"; do
    bam_dec="${LOCAL_INPUT}/${sample}_tumor.bam"
    bam_enc="${LOCAL_INPUT}/${sample}_tumor.bam.c4gh"

    if [[ -f "$bam_dec" ]]; then
        echo "  SKIP $sample (decrypted BAM already local)"
        downloaded+=("$sample"); continue
    fi

    sdir="${LOCAL_RESULTS}/${sample}/results/${sample}"
    [[ -d "${sdir}/hlahd" ]]    && _hh=1 || _hh=0
    [[ -d "${sdir}/optitype" ]] && _op=1 || _op=0
    [[ -d "${sdir}/t1k" ]]      && _t1=1 || _t1=0
    if [[ "$_hh" -eq 1 && "$_op" -eq 1 && "$_t1" -eq 1 ]]; then
        echo "  SKIP $sample (already typed: hlahd+optitype+t1k present)"
        continue
    fi

    # Check segments exist for this sample's BAM
    seg_prefix="${S3_SEGS_BUCKET}/${sample}/${sample}_tumor.bam.c4gh"
    if ! "$RCLONE" lsf "${seg_prefix}/" >/dev/null 2>&1; then
        echo "  WARNING: $sample — no BAM segments at ${seg_prefix}, skipping" >&2
        continue
    fi

    echo "  Downloading $sample encrypted BAM (~15 GB, patience)..."
    if ! download_dlo "$seg_prefix" "$bam_enc"; then
        echo "  ERROR: $sample download failed, skipping" >&2
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
    echo "  All samples in this task already have results — skipping."
    exit 0
fi

# ── Step 2: Submit PIHLA pipeline jobs ────────────────────────────────────
echo "[$(date)] Step 2: submitting PIHLA jobs..."
SAMPLE_SCRIPT="$REPO/slurm_wgs_sample.sh"
job_id_list=""
for sample in "${downloaded[@]}"; do
    job_id=$(sbatch --parsable "$SAMPLE_SCRIPT" \
        "$sample" "$LOCAL_INPUT" "$LOCAL_RESULTS" "$WORK_DIR")
    job_id_list="${job_id_list} ${job_id}"
    echo "  Submitted $sample → job $job_id"
done

# ── Step 3: Wait for pipeline jobs ────────────────────────────────────────
if [[ -n "$job_id_list" ]]; then
    ids=$(echo "$job_id_list" | tr ' ' ',' | sed 's/^,//')
    echo "[$(date)] Step 3: waiting for jobs ($ids)"
    while squeue --jobs="$ids" -h 2>/dev/null | grep -q .; do
        sleep 60
    done
    echo "[$(date)] Pipeline jobs finished"
fi

# ── Step 4: Upload results via S3 ────────────────────────────────────────
ALLAS_RESULTS="s3allas:2014061-DNA_seq_results"
echo "[$(date)] Step 4: uploading results..."
for sample in "${downloaded[@]}"; do
    result_dir="$LOCAL_RESULTS/$sample/results/$sample"
    [[ ! -d "$result_dir" ]] && continue
    n_results=$(find "$result_dir" -name "*.txt" 2>/dev/null | wc -l)
    if [[ "$n_results" -gt 0 ]]; then
        echo "  Uploading $sample ($n_results files)..."
        "$RCLONE" copy "$result_dir" "${ALLAS_RESULTS}/${sample}/" \
            --include "*.txt" --include "*.tsv"
        echo "  Uploaded $sample"
    else
        echo "  WARNING: $sample has no results — skipping upload"
    fi
done

# ── Step 5: Clean local files ─────────────────────────────────────────────
echo "[$(date)] Step 5: cleaning local files..."
for sample in "${downloaded[@]}"; do
    result_dir="$LOCAL_RESULTS/$sample/results/$sample"
    n_results=$(find "$result_dir" -name "*.txt" 2>/dev/null | wc -l)
    if [[ "$n_results" -gt 0 ]]; then
        find "$LOCAL_INPUT" -maxdepth 1 \
            \( -name "${sample}_tumor.bam" -o -name "${sample}_tumor.bam.bai" \) \
            -delete 2>/dev/null || true
        rm -rf "$WORK_DIR/wgs_${sample}" 2>/dev/null || true
        find "$LOCAL_RESULTS/$sample" \
            \( -name "*.bam" -o -name "*.fq.gz" -o -name "*.fastq.gz" \
               -o -name "*.alignment.p" \) -delete 2>/dev/null || true
        echo "  Cleaned $sample"
    else
        echo "  KEPT $sample inputs (no results — investigate)"
    fi
done

echo "[$(date)] Task ${SLURM_ARRAY_TASK_ID} complete. Disk:"
df -h /scratch/project_2008084/ | tail -1
