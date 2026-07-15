#!/bin/bash
#SBATCH --job-name=pihla_venex_fq_arr
#SBATCH --account=project_2008084
#SBATCH --partition=small
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --output=/scratch/project_2008084/pihla-publish/analysis/logs/venex_fq_arr_%A_%a.out
#SBATCH --error=/scratch/project_2008084/pihla-publish/analysis/logs/venex_fq_arr_%A_%a.err

# VENEX WGS array task — downloads FASTQ pairs from Allas and runs the
# FASTQ-input HLA typing pipeline (optitype, t1k, hlahd, spechla).
#
# FASTQ files in Allas follow the DRAGEN naming convention:
#   {sample_id}_S{nn}[_L004]_R{1,2}_001.fastq.gz.c4gh
# The _S{nn} index varies per sample, so filenames are discovered dynamically.
# Both R1 and R2 are Swift DLOs stored in the _segments bucket.
#
# Disk footprint: ~14 GB per sample pair (vs 30–80 GB for full BAM → FASTQ),
# so batch-size 2 peaks at ~28 GB — well within scratch quota.
#
# Prerequisites (set before submitting via submit_venex_fastq_recovery.sh):
#   export C4GH_PASSPHRASE="<passphrase for all_data_csc_key.sec>"
#
# Usage (do not call directly — use submit_venex_fastq_recovery.sh):
#   sbatch --array=0-N%1 --export=ALL slurm_venex_fastq_array.sh \
#     --s3bucket s3allas:2014061-DNA_seq/DNA_seq/DRAGEN/VenEx_DNA_NONHUS_Batch1 \
#     --sample-list <file>    \
#     [--batch-size 2]

set -euo pipefail
trap 'echo "ERR at line $LINENO: $BASH_COMMAND" >&2' ERR

# ── Parse arguments ───────────────────────────────────────────────────────
S3_BUCKET=""
SAMPLE_LIST=""
BATCH_SIZE=2

while [[ $# -gt 0 ]]; do
    case "$1" in
        --s3bucket)    S3_BUCKET="$2";  shift 2 ;;
        --sample-list) SAMPLE_LIST="$2"; shift 2 ;;
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
LOCAL_INPUT=${SCRATCH}/venex/inputs/fastq
LOCAL_RESULTS=${SCRATCH}/venex/results/wgs_fastq
WORK_DIR=${SCRATCH}/venex/work/wgs_fastq
LOG_DIR=${REPO}/analysis/logs/venex

# Derive the _segments bucket from the main bucket name
# s3allas:2014061-DNA_seq/path  →  s3allas:2014061-DNA_seq_segments/path
S3_SEGS_BUCKET="${S3_BUCKET/DNA_seq\//DNA_seq_segments/}"

mkdir -p "$LOCAL_INPUT" "$LOCAL_RESULTS" "$WORK_DIR" "$LOG_DIR"

RCLONE=/appl/opt/csc-cli-utils/allas-cli-utils/rclone
CRYPT4GH=/appl/soft/bio/biopython/gcc_11.3.0/3.10.6_sqlite/bin/python3
KEY_FILE="${TMPDIR:-/tmp}/all_data_csc_key.sec"

module load nextflow

# ── Fetch decryption key from S3 ─────────────────────────────────────────
echo "[$(date)] Fetching decryption key..."
"$RCLONE" copy s3allas:2014061-keys/all_data_csc_key.sec "$(dirname "$KEY_FILE")/"
[[ -f "$KEY_FILE" ]] || { echo "ERROR: failed to download decryption key" >&2; exit 1; }
chmod 600 "$KEY_FILE"

# ── Select this task's batch ──────────────────────────────────────────────
mapfile -t ALL_SAMPLES < "$SAMPLE_LIST"
TOTAL=${#ALL_SAMPLES[@]}
batch=("${ALL_SAMPLES[@]:$OFFSET:$BATCH_SIZE}")

echo "[$(date)] VENEX FASTQ array task ${SLURM_ARRAY_TASK_ID}"
echo "  S3 bucket: $S3_BUCKET"
echo "  Offset:    $OFFSET / $TOTAL"
echo "  Samples:   ${batch[*]}"
echo ""

[[ ${#batch[@]} -eq 0 ]] && { echo "No samples for this task — exiting"; exit 0; }

# ── Helper: download a DLO file via S3 segment concatenation ─────────────
download_dlo() {
    local seg_prefix="$1"
    local out_file="$2"
    local tmpdir
    tmpdir="${out_file}.segtmp"
    mkdir -p "$tmpdir"

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
    find "$tmpdir" -maxdepth 1 -type f | sort -V | xargs cat > "$out_file"
    rm -rf "$tmpdir"
}

# ── Step 1: Download + decrypt each sample's FASTQ pair ──────────────────
echo "[$(date)] Step 1: downloading and decrypting FASTQs from Allas..."
declare -A SAMPLE_FQ1
declare -A SAMPLE_FQ2

for sample in "${batch[@]}"; do
    fq_r1="${LOCAL_INPUT}/${sample}_R1.fastq.gz"
    fq_r2="${LOCAL_INPUT}/${sample}_R2.fastq.gz"

    # Skip if already decrypted locally
    if [[ -f "$fq_r1" && -f "$fq_r2" ]]; then
        echo "  SKIP $sample (decrypted FASTQs already local)"
        SAMPLE_FQ1[$sample]="$fq_r1"
        SAMPLE_FQ2[$sample]="$fq_r2"
        continue
    fi

    # Skip if FASTQ results already complete.
    # Threshold = 10: T1K+SpecHLA alone produce 4 files, which is insufficient —
    # HLA-HD and/or OptiType must also have run. Samples with HLA-HD produce 50+
    # files; OptiType adds 1 more. A count ≥10 means at least one heavy tool ran.
    n_done=$(find "${LOCAL_RESULTS}/${sample}/results/${sample}" \
        \( -name "*.txt" -o -name "*.tsv" \) 2>/dev/null | wc -l) || n_done=0
    if [[ "$n_done" -ge 60 ]]; then
        echo "  SKIP $sample (FASTQ results already present: $n_done files)"
        continue
    fi

    # Discover FASTQ filenames dynamically (DRAGEN naming: {sample}_S{nn}[_L004]_R{1,2}_001.fastq.gz.c4gh)
    r1_enc_name=$("$RCLONE" lsf "${S3_BUCKET}/${sample}/" 2>/dev/null | grep "_R1_001\.fastq\.gz\.c4gh$" | head -1) || r1_enc_name=""
    r2_enc_name=$("$RCLONE" lsf "${S3_BUCKET}/${sample}/" 2>/dev/null | grep "_R2_001\.fastq\.gz\.c4gh$" | head -1) || r2_enc_name=""

    if [[ -z "$r1_enc_name" || -z "$r2_enc_name" ]]; then
        echo "  WARNING: $sample — FASTQ files not found at ${S3_BUCKET}/${sample}/, skipping" >&2
        continue
    fi
    echo "  Found: $r1_enc_name / $r2_enc_name"

    # Download R1 — try DLO segment concatenation first; fall back to direct copy
    # (some samples store FASTQs as regular S3 objects instead of Swift DLOs)
    seg_prefix_r1="${S3_SEGS_BUCKET}/${sample}/${r1_enc_name}"
    fq_enc_r1="${LOCAL_INPUT}/${sample}_R1.fastq.gz.c4gh"
    echo "  Downloading $sample R1 (~7 GB, patience)..."
    if ! download_dlo "$seg_prefix_r1" "$fq_enc_r1"; then
        echo "  DLO not found — trying direct copy for $sample R1..."
        if "$RCLONE" copy "${S3_BUCKET}/${sample}/${r1_enc_name}" "$LOCAL_INPUT/"; then
            mv "${LOCAL_INPUT}/${r1_enc_name}" "$fq_enc_r1" 2>/dev/null || true
        else
            echo "  ERROR: $sample R1 download failed (DLO and direct both failed), skipping" >&2
            rm -f "$fq_enc_r1" "${LOCAL_INPUT}/${r1_enc_name}" 2>/dev/null || true; continue
        fi
    fi
    if [[ ! -s "$fq_enc_r1" ]]; then
        echo "  ERROR: $sample R1 — empty after download, skipping" >&2
        continue
    fi

    # Download R2 — same DLO-first, direct-copy fallback
    seg_prefix_r2="${S3_SEGS_BUCKET}/${sample}/${r2_enc_name}"
    fq_enc_r2="${LOCAL_INPUT}/${sample}_R2.fastq.gz.c4gh"
    echo "  Downloading $sample R2 (~7 GB, patience)..."
    if ! download_dlo "$seg_prefix_r2" "$fq_enc_r2"; then
        echo "  DLO not found — trying direct copy for $sample R2..."
        if "$RCLONE" copy "${S3_BUCKET}/${sample}/${r2_enc_name}" "$LOCAL_INPUT/"; then
            mv "${LOCAL_INPUT}/${r2_enc_name}" "$fq_enc_r2" 2>/dev/null || true
        else
            echo "  ERROR: $sample R2 download failed (DLO and direct both failed), skipping" >&2
            rm -f "$fq_enc_r1" "$fq_enc_r2" "${LOCAL_INPUT}/${r2_enc_name}" 2>/dev/null || true; continue
        fi
    fi
    if [[ ! -s "$fq_enc_r2" ]]; then
        echo "  ERROR: $sample R2 — empty after download, skipping" >&2
        rm -f "$fq_enc_r1"; continue
    fi

    # Decrypt R1
    echo "  Decrypting $sample R1..."
    C4GH_PASSPHRASE="$C4GH_PASSPHRASE" \
        "$CRYPT4GH" -m crypt4gh decrypt --sk "$KEY_FILE" \
        < "$fq_enc_r1" > "$fq_r1"
    rm -f "$fq_enc_r1"
    if [[ ! -s "$fq_r1" ]]; then
        echo "  ERROR: $sample R1 — empty after decryption (wrong passphrase?)" >&2
        rm -f "$fq_r1"; continue
    fi

    # Decrypt R2
    echo "  Decrypting $sample R2..."
    C4GH_PASSPHRASE="$C4GH_PASSPHRASE" \
        "$CRYPT4GH" -m crypt4gh decrypt --sk "$KEY_FILE" \
        < "$fq_enc_r2" > "$fq_r2"
    rm -f "$fq_enc_r2"
    if [[ ! -s "$fq_r2" ]]; then
        echo "  ERROR: $sample R2 — empty after decryption (wrong passphrase?)" >&2
        rm -f "$fq_r1" "$fq_r2"; continue
    fi

    echo "  Decrypted $(du -sh "$fq_r1" | cut -f1) R1 + $(du -sh "$fq_r2" | cut -f1) R2"
    SAMPLE_FQ1[$sample]="$fq_r1"
    SAMPLE_FQ2[$sample]="$fq_r2"
    echo "  $sample ready"
done

rm -f "$KEY_FILE"

n_ready=0
for _k in "${!SAMPLE_FQ1[@]}"; do (( n_ready++ )) || true; done
if [[ $n_ready -eq 0 ]]; then
    echo "  All samples in this task already have results — skipping."
    exit 0
fi

# ── Step 2: Submit PIHLA FASTQ pipeline jobs ──────────────────────────────
echo "[$(date)] Step 2: submitting PIHLA FASTQ jobs..."
SAMPLE_SCRIPT="$REPO/slurm_venex_fastq_sample.sh"
job_id_list=""
for sample in "${!SAMPLE_FQ1[@]}"; do
    job_id=$(sbatch --parsable "$SAMPLE_SCRIPT" \
        "$sample" "${SAMPLE_FQ1[$sample]}" "${SAMPLE_FQ2[$sample]}" \
        "$LOCAL_RESULTS" "$WORK_DIR")
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

# ── Step 4: Upload results to Allas ──────────────────────────────────────
ALLAS_RESULTS="s3allas:2014061-DNA_seq_results/fastq"
echo "[$(date)] Step 4: uploading results..."
for sample in "${!SAMPLE_FQ1[@]}"; do
    result_dir="$LOCAL_RESULTS/$sample/results/$sample"
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

# ── Step 5: Clean local FASTQ files ──────────────────────────────────────
echo "[$(date)] Step 5: cleaning local files..."
for sample in "${!SAMPLE_FQ1[@]}"; do
    result_dir="$LOCAL_RESULTS/$sample/results/$sample"
    n_results=$(find "$result_dir" \( -name "*.txt" -o -name "*.tsv" \) 2>/dev/null | wc -l)
    if [[ "$n_results" -ge 60 ]]; then
        # Sufficient results (≥10 files means at least HLA-HD or OptiType ran) —
        # safe to clean FASTQs and work dir.
        rm -f "${SAMPLE_FQ1[$sample]}" "${SAMPLE_FQ2[$sample]}" 2>/dev/null || true
        rm -rf "$WORK_DIR/${sample}" 2>/dev/null || true
        find "$LOCAL_RESULTS/$sample" \
            \( -name "*.bam" -o -name "*.fq.gz" -o -name "*.fastq.gz" \
               -o -name "*.alignment.p" \) -delete 2>/dev/null || true
        echo "  Cleaned $sample"
    elif [[ "$n_results" -gt 0 ]]; then
        # Partial results (T1K/SpecHLA only) — keep FASTQs and work dir so the
        # next recovery pass can resume from cache and only re-run the failed tools.
        find "$LOCAL_RESULTS/$sample" \
            \( -name "*.bam" -o -name "*.alignment.p" \) -delete 2>/dev/null || true
        echo "  PARTIAL $sample ($n_results files) — kept FASTQs and work dir for next recovery"
    else
        echo "  KEPT $sample inputs (no results — investigate)"
    fi
done

echo "[$(date)] Task ${SLURM_ARRAY_TASK_ID} complete. Disk:"
df -h /scratch/project_2008084/ | tail -1
