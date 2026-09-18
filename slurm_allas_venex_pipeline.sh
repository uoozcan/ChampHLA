#!/bin/bash
#SBATCH --job-name=pihla_venex
#SBATCH --account=project_2008084
#SBATCH --partition=small
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --output=/scratch/project_2008084/pihla-publish/analysis/logs/venex_%j.out
#SBATCH --error=/scratch/project_2008084/pihla-publish/analysis/logs/venex_%j.err

# ────────────────────────────────────────────────────────────────────────────
# PIHLA batch pipeline for VENEX samples stored in CSC Allas.
#
# Workflow per batch:
#   1. Authenticate with Allas (uses pre-saved token or prompts)
#   2. Download CRAM/BAM/FASTQ files for BATCH_SIZE samples
#   3. Run PIHLA HLA typing pipeline (Nextflow)
#   4. Upload result TXT files back to Allas results bucket
#   5. Delete local inputs and large intermediates
#   6. Repeat for next batch
#
# Prerequisites:
#   - Run 'allas-conf' interactively once to save authentication token
#   - ALLAS_BUCKET: Allas bucket containing VENEX input files
#   - SAMPLE_LIST: text file listing one sample ID per line
#   - File naming convention in Allas: {sample_id}.{cram|bam|fastq.gz}
#
# Usage:
#   sbatch slurm_allas_venex_pipeline.sh \
#     --bucket allas:venex-samples \
#     --sample-list /path/to/venex_samples.txt \
#     --modality wgs \
#     --batch-size 5 \
#     [--offset 0]   # start from this sample index (for resuming)
# ────────────────────────────────────────────────────────────────────────────

set -euo pipefail

# ── Parse arguments ───────────────────────────────────────────────────────────
ALLAS_BUCKET=""
SAMPLE_LIST=""
MODALITY="wgs"
BATCH_SIZE=5
OFFSET=0

while [[ $# -gt 0 ]]; do
    case "$1" in
        --bucket)      ALLAS_BUCKET="$2"; shift 2 ;;
        --sample-list) SAMPLE_LIST="$2";  shift 2 ;;
        --modality)    MODALITY="$2";     shift 2 ;;
        --batch-size)  BATCH_SIZE="$2";   shift 2 ;;
        --offset)      OFFSET="$2";       shift 2 ;;
        *) echo "Unknown arg: $1" >&2; exit 1 ;;
    esac
done

[[ -z "$ALLAS_BUCKET" ]] && { echo "ERROR: --bucket required" >&2; exit 1; }
[[ -z "$SAMPLE_LIST" ]]  && { echo "ERROR: --sample-list required" >&2; exit 1; }

# ── Paths ─────────────────────────────────────────────────────────────────────
REPO=/scratch/project_2008084/pihla-publish
SCRATCH=/scratch/project_2008084/hla_calibration
LOCAL_INPUT=${SCRATCH}/venex/inputs/${MODALITY}
LOCAL_RESULTS=${SCRATCH}/venex/results/${MODALITY}
WORK_DIR=${SCRATCH}/venex/work/${MODALITY}
ALLAS_RESULTS="${ALLAS_BUCKET}_results"
LOG_DIR=${REPO}/analysis/logs/venex

mkdir -p "$LOCAL_INPUT" "$LOCAL_RESULTS" "$WORK_DIR" "$LOG_DIR"
RCLONE=/appl/opt/csc-cli-utils/allas-cli-utils/rclone

# ── Load modules ──────────────────────────────────────────────────────────────
module load allas
module load nextflow
module load samtools

echo "[$(date)] VENEX PIHLA pipeline"
echo "  Bucket:     $ALLAS_BUCKET"
echo "  Modality:   $MODALITY"
echo "  Batch size: $BATCH_SIZE"
echo "  Offset:     $OFFSET"
echo ""

# ── Read sample list ──────────────────────────────────────────────────────────
mapfile -t ALL_SAMPLES < "$SAMPLE_LIST"
TOTAL=${#ALL_SAMPLES[@]}
echo "[$(date)] Total samples in list: $TOTAL"

# ── Process batches ───────────────────────────────────────────────────────────
for (( i=OFFSET; i<TOTAL; i+=BATCH_SIZE )); do
    batch=("${ALL_SAMPLES[@]:$i:$BATCH_SIZE}")
    batch_num=$(( i/BATCH_SIZE + 1 ))
    echo ""
    echo "══════════════════════════════════════════════════"
    echo "[$(date)] BATCH ${batch_num}: samples ${batch[*]}"
    echo "══════════════════════════════════════════════════"

    # ── Step 1: Download from Allas ───────────────────────────────────────────
    echo "[$(date)] Step 1: downloading from Allas..."
    downloaded=()
    for sample in "${batch[@]}"; do
        # Try common file extensions
        for ext in cram bam fastq.gz fq.gz; do
            remote="${ALLAS_BUCKET}/${sample}.${ext}"
            local_file="${LOCAL_INPUT}/${sample}.${ext}"
            if [[ -f "$local_file" ]]; then
                echo "  SKIP $sample.$ext (already local)"
                downloaded+=("$sample")
                break
            fi
            if "$RCLONE" ls "$remote" >/dev/null 2>&1; then
                echo "  Downloading $sample.$ext..."
                "$RCLONE" copy "$remote" "$LOCAL_INPUT/" \
                    --progress --transfers 4 2>/dev/null
                # For CRAM: also download index if available
                "$RCLONE" copy "${remote%.cram}.cram.crai" "$LOCAL_INPUT/" 2>/dev/null || true
                "$RCLONE" copy "${remote%.bam}.bam.bai"    "$LOCAL_INPUT/" 2>/dev/null || true
                echo "  Downloaded $sample.$ext ($(du -sh "$local_file" | cut -f1))"
                downloaded+=("$sample")
                break
            fi
        done
        # Check paired FASTQs (R1/R2 convention)
        r1="${LOCAL_INPUT}/${sample}_R1.fastq.gz"
        if [[ ! -f "$r1" ]]; then
            for r1_remote in "${ALLAS_BUCKET}/${sample}_R1.fastq.gz" \
                             "${ALLAS_BUCKET}/${sample}_1.fastq.gz"; do
                if "$RCLONE" ls "$r1_remote" >/dev/null 2>&1; then
                    echo "  Downloading $sample FASTQ pair..."
                    r2_remote="${r1_remote/_R1./_R2.}"; r2_remote="${r1_remote/_1./_2.}"
                    "$RCLONE" copy "$r1_remote" "$LOCAL_INPUT/" 2>/dev/null
                    "$RCLONE" copy "$r2_remote" "$LOCAL_INPUT/" 2>/dev/null || true
                    downloaded+=("$sample")
                    break
                fi
            done
        fi
    done

    [[ ${#downloaded[@]} -eq 0 ]] && echo "  WARNING: no files downloaded for this batch" && continue

    # ── Step 2: Submit PIHLA per-sample jobs ──────────────────────────────────
    echo "[$(date)] Step 2: submitting PIHLA pipeline jobs..."
    job_id_list=""
    for sample in "${downloaded[@]}"; do
        # Find the local input file
        input_file=$(find "$LOCAL_INPUT" -maxdepth 1 \
            \( -name "${sample}.cram" -o -name "${sample}.bam" \
               -o -name "${sample}_R1.fastq.gz" \) 2>/dev/null | head -1)
        [[ -z "$input_file" ]] && echo "  SKIP $sample: no local file found" && continue

        result_dir="$LOCAL_RESULTS/$sample/results"
        n_done=$(find "$result_dir" -name "*_optitype.txt" 2>/dev/null | wc -l)
        [[ "$n_done" -gt 0 ]] && echo "  SKIP $sample: already has results" && continue

        # Determine sample script based on modality and input type
        ext="${input_file##*.}"
        if [[ "$ext" == "cram" || "$ext" == "bam" ]]; then
            SAMPLE_SCRIPT="$REPO/slurm_${MODALITY}_sample.sh"
        else
            SAMPLE_SCRIPT="$REPO/slurm_rna_sample.sh"
        fi

        if [[ ! -f "$SAMPLE_SCRIPT" ]]; then
            echo "  WARNING: $SAMPLE_SCRIPT not found for $sample"
            continue
        fi

        job_id=$(sbatch --parsable "$SAMPLE_SCRIPT" \
            "$sample" "$LOCAL_INPUT" "$LOCAL_RESULTS" "$WORK_DIR")
        job_id_list="${job_id_list} ${job_id}"
        echo "  Submitted $sample → job $job_id"
    done

    # ── Step 3: Wait for pipeline jobs ────────────────────────────────────────
    if [[ -n "$job_id_list" ]]; then
        ids=$(echo "$job_id_list" | tr ' ' ',' | sed 's/^,//')
        echo "[$(date)] Step 3: waiting for pipeline jobs ($ids)..."
        while squeue --jobs="$ids" -h 2>/dev/null | grep -q .; do
            sleep 60
        done
        echo "[$(date)] Pipeline jobs finished"
    fi

    # ── Step 4: Upload results to Allas ──────────────────────────────────────
    echo "[$(date)] Step 4: uploading results to Allas..."
    for sample in "${downloaded[@]}"; do
        result_dir="$LOCAL_RESULTS/$sample/results/$sample"
        [[ ! -d "$result_dir" ]] && continue

        n_results=$(find "$result_dir" -name "*.txt" 2>/dev/null | wc -l)
        if [[ "$n_results" -gt 0 ]]; then
            echo "  Uploading $sample results ($n_results files)..."
            "$RCLONE" copy "$result_dir" \
                "${ALLAS_RESULTS}/${sample}/" \
                --include "*.txt" --include "*.tsv" \
                --progress 2>/dev/null
            echo "  Uploaded $sample ✓"
        else
            echo "  WARNING: $sample has no result files — skipping upload"
        fi
    done

    # ── Step 5: Clean up local files ─────────────────────────────────────────
    echo "[$(date)] Step 5: cleaning local files..."
    for sample in "${downloaded[@]}"; do
        result_dir="$LOCAL_RESULTS/$sample/results/$sample"
        n_results=$(find "$result_dir" -name "*.txt" 2>/dev/null | wc -l)

        if [[ "$n_results" -gt 0 ]]; then
            # Remove input files
            find "$LOCAL_INPUT" -maxdepth 1 \
                \( -name "${sample}.*" -o -name "${sample}_R*.fastq.gz" \) \
                -delete 2>/dev/null || true
            # Remove large Nextflow intermediates
            rm -rf "$WORK_DIR/${MODALITY}_${sample}" 2>/dev/null || true
            rm -rf "$WORK_DIR/rna_${sample}"         2>/dev/null || true
            # Remove large pipeline intermediate outputs (keep result TXTs)
            find "$LOCAL_RESULTS/$sample" -name "*.bam" -delete 2>/dev/null || true
            find "$LOCAL_RESULTS/$sample" -name "*.fq.gz" -delete 2>/dev/null || true
            find "$LOCAL_RESULTS/$sample" -name "*.fastq.gz" -delete 2>/dev/null || true
            find "$LOCAL_RESULTS/$sample" -name "*.alignment.p" -delete 2>/dev/null || true
            echo "  Cleaned $sample"
        else
            echo "  KEPT input for $sample (no results — investigate)"
        fi
    done

    echo "[$(date)] Batch $batch_num complete"
    df -h /scratch/project_2008084/ | tail -1
done

echo ""
echo "[$(date)] All batches complete"
echo "Results uploaded to: $ALLAS_RESULTS"
echo "Local results dir:   $LOCAL_RESULTS"
