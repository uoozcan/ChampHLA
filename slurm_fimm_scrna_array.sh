#!/bin/bash
#SBATCH --job-name=pihla_fimm_scrna
#SBATCH --account=project_2008084
#SBATCH --partition=small
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --output=/scratch/project_2008084/pihla-publish/analysis/logs/fimm/scrna_arr_%A_%a.out
#SBATCH --error=/scratch/project_2008084/pihla-publish/analysis/logs/fimm/scrna_arr_%A_%a.err

# FIMM scRNA array task — each task processes one batch of samples from Allas.
#
# Two pipeline modes (set via --mode):
#
#   standard (default):
#     Files are standard paired FASTQs (R1/R2 same read length).
#     Runs PIHLA RNA pipeline: arcashla, optitype, seq2hla, t1k, spechla, hlahd.
#     Use when scRNA data is stored as cDNA R1/R2 pairs of equal length
#     (e.g., Smart-seq2 or pseudo-bulk FASTQs).
#
#   tenx:
#     Files are 10X Chromium format (R1 = 28bp barcode+UMI, R2 = cDNA).
#     Runs T1K standalone in scRNA barcode mode — skips tools that require
#     standard paired reads. T1K called directly (not via Nextflow).
#     Use when R1 is short (~28 bp).
#
# !! Before running:
#   1. Discover scRNA file location:
#      rclone ls s3allas:psergeev-2012380-MISC/ | grep -i "scrna\|single\|10x\|cellranger"
#   2. Update S3_BUCKET below with the correct path.
#   3. Check R1 read length to confirm the mode (standard vs tenx).
#   4. Confirm KEY_REMOTE_NAME after listing s3allas:2012380-keys/
#
# Prerequisites (set before submitting via submit_fimm_scrna.sh):
#   export C4GH_PASSPHRASE="<passphrase for FIMM key>"
#
# Usage (do not call directly — use submit_fimm_scrna.sh):
#   sbatch --array=0-N%1 --export=ALL slurm_fimm_scrna_array.sh \
#     --sample-list <file> [--batch-size 3] [--mode standard|tenx] \
#     [--s3bucket s3allas:psergeev-2012380-MISC/path/to/scrna]

set -euo pipefail
trap 'echo "ERR at line $LINENO: $BASH_COMMAND" >&2' ERR

# ── Parse arguments ───────────────────────────────────────────────────────
SAMPLE_LIST=""
BATCH_SIZE=3
MODE="standard"
# !! UPDATE S3_BUCKET to the correct scRNA path after discovery
S3_BUCKET="s3allas:psergeev-2012380-MISC/fimm_ga2_heckman"

while [[ $# -gt 0 ]]; do
    case "$1" in
        --sample-list) SAMPLE_LIST="$2"; shift 2 ;;
        --batch-size)  BATCH_SIZE="$2";  shift 2 ;;
        --mode)        MODE="$2";        shift 2 ;;
        --s3bucket)    S3_BUCKET="$2";   shift 2 ;;
        *) echo "Unknown arg: $1" >&2; exit 1 ;;
    esac
done

[[ -z "$SAMPLE_LIST" ]] && { echo "ERROR: --sample-list required" >&2; exit 1; }
[[ -z "${C4GH_PASSPHRASE:-}" ]] && { echo "ERROR: C4GH_PASSPHRASE not set" >&2; exit 1; }
[[ "$MODE" == "standard" || "$MODE" == "tenx" ]] || \
    { echo "ERROR: --mode must be 'standard' or 'tenx'" >&2; exit 1; }

OFFSET=$(( SLURM_ARRAY_TASK_ID * BATCH_SIZE ))

# ── Paths ─────────────────────────────────────────────────────────────────
REPO=/scratch/project_2008084/pihla-publish
SCRATCH=/scratch/project_2008084/hla_calibration
LOCAL_INPUT=${SCRATCH}/fimm/inputs/scrna
LOCAL_RESULTS=${SCRATCH}/fimm/results/scrna
WORK_DIR=${SCRATCH}/fimm/work/scrna
LOG_DIR=${REPO}/analysis/logs/fimm

S3_SEGS_BUCKET="${S3_BUCKET/psergeev-2012380-MISC/psergeev-2012380-MISC_segments}"
S3_KEYS="s3allas:2012380-keys"

# !! UPDATE THIS after running: rclone ls s3allas:2012380-keys/
KEY_REMOTE_NAME="all_data_csc_key.sec"

RCLONE=/appl/opt/csc-cli-utils/allas-cli-utils/rclone
CRYPT4GH=/appl/soft/bio/biopython/gcc_11.3.0/3.10.6_sqlite/bin/python3
KEY_FILE="${TMPDIR:-/tmp}/fimm_key.sec"

# T1K binary location (for tenx mode)
T1K_BIN=/appl/soft/bio/t1k/bin/run-t1k
T1K_HLAIDX=/scratch/project_2008084/hla_references/t1k_hlaidx

mkdir -p "$LOCAL_INPUT" "$LOCAL_RESULTS" "$WORK_DIR" "$LOG_DIR"
module load nextflow

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

echo "[$(date)] FIMM scRNA array task ${SLURM_ARRAY_TASK_ID} (mode: $MODE)"
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

# ── Step 1: Download + decrypt each sample's FASTQ pair ──────────────────
echo "[$(date)] Step 1: downloading and decrypting scRNA FASTQ files from Allas..."
downloaded=()

for sample in "${batch[@]}"; do
    fastq1_dec="${LOCAL_INPUT}/${sample}_R1.fastq.gz"
    fastq2_dec="${LOCAL_INPUT}/${sample}_R2.fastq.gz"

    if [[ -f "$fastq1_dec" && -f "$fastq2_dec" ]]; then
        echo "  SKIP $sample (decrypted FASTQs already local)"
        downloaded+=("$sample"); continue
    fi

    n_done=$(find "${LOCAL_RESULTS}/${sample}" \
        \( -name "*.txt" -o -name "*.tsv" \) 2>/dev/null | wc -l) || n_done=0
    if [[ "$n_done" -ge 5 ]]; then
        echo "  SKIP $sample (already processed: $n_done result files)"
        continue
    fi

    echo "  Listing Allas for $sample..."
    sample_files=$("$RCLONE" lsf "${S3_BUCKET}/${sample}/" 2>/dev/null) || sample_files=""
    if [[ -z "$sample_files" ]]; then
        sample_files=$("$RCLONE" lsf "${S3_BUCKET}/" 2>/dev/null | grep "^${sample}" || true)
    fi
    if [[ -z "$sample_files" ]]; then
        echo "  WARNING: $sample — no files found at ${S3_BUCKET}/${sample}/" >&2
        continue
    fi

    r1_enc_name=$(echo "$sample_files" | grep -i "_R1_\|_R1\." | grep "\.c4gh$" | head -1 || true)
    r2_enc_name=$(echo "$sample_files" | grep -i "_R2_\|_R2\." | grep "\.c4gh$" | head -1 || true)

    if [[ -z "$r1_enc_name" || -z "$r2_enc_name" ]]; then
        echo "  WARNING: $sample — could not identify R1/R2 encrypted FASTQ files" >&2
        echo "  Found: $sample_files" >&2
        continue
    fi

    for read_num in 1 2; do
        [[ "$read_num" -eq 1 ]] && enc_name="$r1_enc_name" && dec_file="$fastq1_dec"
        [[ "$read_num" -eq 2 ]] && enc_name="$r2_enc_name" && dec_file="$fastq2_dec"

        enc_file="${LOCAL_INPUT}/${sample}_R${read_num}.fastq.gz.c4gh"
        s3_path="${S3_BUCKET}/${enc_name}"

        echo "  Downloading $sample R${read_num} ($enc_name)..."
        if ! download_file "$s3_path" "$enc_file"; then
            echo "  ERROR: $sample R${read_num} download failed" >&2
            rm -f "$enc_file"; continue 2
        fi
        if [[ ! -s "$enc_file" ]]; then
            echo "  ERROR: $sample R${read_num} — empty after download" >&2
            rm -f "$enc_file"; continue 2
        fi

        echo "  Decrypting $sample R${read_num}..."
        C4GH_PASSPHRASE="$C4GH_PASSPHRASE" \
            "$CRYPT4GH" -m crypt4gh decrypt --sk "$KEY_FILE" \
            < "$enc_file" > "$dec_file"
        rm -f "$enc_file"
        if [[ ! -s "$dec_file" ]]; then
            echo "  ERROR: $sample R${read_num} — empty after decryption" >&2
            rm -f "$dec_file"; continue 2
        fi
        echo "  Decrypted R${read_num}: $(du -sh "$dec_file" | cut -f1)"
    done

    [[ -f "$fastq1_dec" && -f "$fastq2_dec" ]] && downloaded+=("$sample") && \
        echo "  $sample ready"
done

rm -f "$KEY_FILE"

if [[ ${#downloaded[@]} -eq 0 ]]; then
    echo "WARNING: no samples ready for this task — check Allas paths and C4GH_PASSPHRASE"
    exit 1
fi

# ── Step 2a: Standard mode — submit PIHLA RNA pipeline ───────────────────
if [[ "$MODE" == "standard" ]]; then
    echo "[$(date)] Step 2 (standard): submitting PIHLA RNA jobs..."
    SAMPLE_SCRIPT="$REPO/slurm_rna_sample.sh"
    job_id_list=""
    for sample in "${downloaded[@]}"; do
        job_id=$(sbatch --parsable "$SAMPLE_SCRIPT" \
            "$sample" "$LOCAL_INPUT" "$LOCAL_RESULTS" "$WORK_DIR")
        job_id_list="${job_id_list} ${job_id}"
        echo "  Submitted $sample → job $job_id"
    done

    if [[ -n "$job_id_list" ]]; then
        ids=$(echo "$job_id_list" | tr ' ' ',' | sed 's/^,//')
        echo "[$(date)] Waiting for jobs ($ids)"
        while squeue --jobs="$ids" -h 2>/dev/null | grep -q .; do
            sleep 60
        done
        echo "[$(date)] Pipeline jobs finished"
    fi
fi

# ── Step 2b: 10X mode — run T1K in scRNA barcode mode ────────────────────
if [[ "$MODE" == "tenx" ]]; then
    echo "[$(date)] Step 2 (tenx): running T1K in scRNA barcode mode..."
    for sample in "${downloaded[@]}"; do
        out_dir="${LOCAL_RESULTS}/${sample}/t1k_scrna"
        mkdir -p "$out_dir"

        echo "  Running T1K scRNA for $sample..."
        # R1 = barcode+UMI (28bp for 10X 3' v3: 16bp barcode + 12bp UMI)
        # R2 = cDNA read
        # --barcode: 10xv3 chemistry (16+12bp)
        # Adjust --barcode to match the library prep chemistry used
        "$T1K_BIN" \
            -1 "${LOCAL_INPUT}/${sample}_R1.fastq.gz" \
            -2 "${LOCAL_INPUT}/${sample}_R2.fastq.gz" \
            --ref "${T1K_HLAIDX}/hlaidx_rna_seq.fa" \
            --od "$out_dir" \
            -t 4 \
            --barcode 10xv3 \
            --preset hla \
            2>&1 | tee "${out_dir}/${sample}_t1k_scrna.log"

        echo "  T1K scRNA complete for $sample"
    done
fi

# ── Step 3: Upload results to Allas ──────────────────────────────────────
ALLAS_RESULTS="s3allas:psergeev-2012380-MISC-results/fimm/scrna"
echo "[$(date)] Step 3: uploading results..."
for sample in "${downloaded[@]}"; do
    if [[ "$MODE" == "standard" ]]; then
        result_dir="${LOCAL_RESULTS}/${sample}/results/${sample}"
    else
        result_dir="${LOCAL_RESULTS}/${sample}/t1k_scrna"
    fi
    [[ ! -d "$result_dir" ]] && continue
    n_results=$(find "$result_dir" \( -name "*.txt" -o -name "*.tsv" -o -name "*.log" \) \
        2>/dev/null | wc -l)
    if [[ "$n_results" -gt 0 ]]; then
        echo "  Uploading $sample ($n_results files)..."
        "$RCLONE" copy "$result_dir" "${ALLAS_RESULTS}/${sample}/" \
            --include "*.txt" --include "*.tsv" --include "*.log"
        echo "  Uploaded $sample"
    else
        echo "  WARNING: $sample has no results — skipping upload"
    fi
done

# ── Step 4: Clean local FASTQ files ──────────────────────────────────────
echo "[$(date)] Step 4: cleaning local FASTQ files..."
for sample in "${downloaded[@]}"; do
    if [[ "$MODE" == "standard" ]]; then
        result_dir="${LOCAL_RESULTS}/${sample}/results/${sample}"
        min_files=5
    else
        result_dir="${LOCAL_RESULTS}/${sample}/t1k_scrna"
        min_files=1
    fi
    n_results=$(find "$result_dir" \( -name "*.txt" -o -name "*.tsv" -o -name "*.log" \) \
        2>/dev/null | wc -l)
    if [[ "$n_results" -ge "$min_files" ]]; then
        rm -f "${LOCAL_INPUT}/${sample}_R1.fastq.gz" \
               "${LOCAL_INPUT}/${sample}_R2.fastq.gz"
        find "${LOCAL_RESULTS}/${sample}" \
            \( -name "*.fq.gz" -o -name "*.fastq.gz" -o -name "*.alignment.p" \) \
            -delete 2>/dev/null || true
        echo "  Cleaned $sample"
    else
        echo "  KEPT $sample inputs (insufficient results — investigate)"
    fi
done

echo "[$(date)] Task ${SLURM_ARRAY_TASK_ID} complete. Disk:"
df -h /scratch/project_2008084/ | tail -1
