#!/bin/bash
#SBATCH --job-name=pihla_wes_b10
#SBATCH --account=project_2008084
#SBATCH --partition=small
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --output=/scratch/project_2008084/pihla-publish/analysis/logs/wes_batch10_%j.out
#SBATCH --error=/scratch/project_2008084/pihla-publish/analysis/logs/wes_batch10_%j.err

# Downloads (if needed), types, and cleans 10 WES samples.
# Usage: sbatch slurm_wes_batch10.sh <sample1> <sample2> ... <sample10>
#
# WES BAMs are small (HLA region only, ~100MB each) so this is mainly
# for consistency with the RNA batching approach.

set -euo pipefail

WES_INDEX_DIR=/scratch/project_2008084/hla_calibration/wes/index
WES_INDICES=(
    "$WES_INDEX_DIR/sample_bam_urls_batch5.tsv"
    "$WES_INDEX_DIR/sample_bam_urls_batch2.tsv"
    "$WES_INDEX_DIR/sample_bam_urls_batch3.tsv"
    "$WES_INDEX_DIR/sample_bam_urls_batch4.tsv"
    "$WES_INDEX_DIR/sample_bam_urls_recovery_wave001.tsv"
    "$WES_INDEX_DIR/sample_bam_urls.tsv"
)
BAM_ROOT=/scratch/project_2008084/hla_calibration/wes/bams
BATCH_ROOT=/scratch/project_2008084/hla_calibration/wes_batches
WORK_ROOT=/scratch/project_2008084/hla_calibration/work
PIPELINE_DIR=/scratch/project_2008084/pihla-publish
SAMTOOLS=/projappl/project_2008084/bin/samtools

mkdir -p "$BAM_ROOT" /scratch/project_2008084/hla_calibration/logs

samples=("$@")
if [[ ${#samples[@]} -eq 0 ]]; then
    echo "Usage: sbatch slurm_wes_batch10.sh <sample1> [sample2 ...]" >&2
    exit 1
fi

echo "[$(date)] WES batch: ${#samples[@]} samples: ${samples[*]}"
module load samtools 2>/dev/null || true
if "${SAMTOOLS}" --version >/dev/null 2>&1; then SAM="${SAMTOOLS}"; else SAM=samtools; fi

# ── Step 1: Download BAMs (HLA region only) ───────────────────────────────
echo "[$(date)] Step 1: downloading HLA-region BAMs"
for sample in "${samples[@]}"; do
    bam_count=$(find "$BAM_ROOT" -maxdepth 1 -name "${sample}*.bam" 2>/dev/null | wc -l)
    [[ "$bam_count" -gt 0 ]] && echo "  SKIP $sample (BAM exists)" && continue

    line=""
    for idx_file in "${WES_INDICES[@]}"; do
        [[ -f "$idx_file" ]] || continue
        line=$(grep "^${sample}"$'\t' "$idx_file" 2>/dev/null | head -1 || true)
        [[ -n "$line" ]] && break
    done
    [[ -z "$line" ]] && echo "  SKIP $sample: not in any WES index" && continue

    primary_url=$(echo "$line" | cut -f2)
    alt_urls=$(echo "$line" | cut -f3 | tr '|' ' ')
    bam_filename=$(basename "$primary_url")
    out_bam="$BAM_ROOT/$bam_filename"

    echo "  Downloading HLA region for $sample..."
    tmp_bam="${out_bam}.tmp"
    downloaded=0
    for try_url in "$primary_url" $alt_urls; do
        [[ -z "$try_url" ]] && continue
        for region in "6:28000000-34000000" "chr6:28,000,000-34,000,000" "chr6:28000000-34000000"; do
            rm -f "$tmp_bam"
            "$SAM" view -b -o "$tmp_bam" "$try_url" "$region" 2>/dev/null || continue
            read_count=$("$SAM" view -c "$tmp_bam" 2>/dev/null || echo "0")
            if [[ "$read_count" -gt 0 ]]; then
                mv "$tmp_bam" "$out_bam"
                "$SAM" index "$out_bam"
                echo "  Downloaded $sample: $read_count reads ($(du -sh "$out_bam" | cut -f1)) [url: $(basename $try_url), region: $region]"
                downloaded=1
                break 2
            fi
        done
    done
    rm -f "$tmp_bam"
    [[ "$downloaded" -eq 0 ]] && echo "  WARNING: all regions produced 0 reads for $sample — skipping"
done
echo "[$(date)] BAM downloads complete"

# ── Step 2: Submit pipeline jobs ──────────────────────────────────────────
echo "[$(date)] Step 2: submitting pipeline jobs"
job_id_list=""
for sample in "${samples[@]}"; do
    bam_count=$(find "$BAM_ROOT" -maxdepth 1 -name "${sample}*.bam" 2>/dev/null | wc -l)
    [[ "$bam_count" -eq 0 ]] && echo "  SKIP $sample: no BAM" && continue

    n_done=$(find "$BATCH_ROOT" -path "*/$sample/results/**/*_optitype.txt" 2>/dev/null | wc -l)
    [[ "$n_done" -gt 0 ]] && echo "  SKIP $sample: already has results" && continue

    job_id=$(sbatch --parsable \
        "$PIPELINE_DIR/slurm_wes_sample.sh" \
        "$sample" "$BAM_ROOT" "$BATCH_ROOT" "$WORK_ROOT")
    job_id_list="${job_id_list} ${job_id}"
    echo "  Submitted $sample → job $job_id"
done

# ── Step 3: Wait ──────────────────────────────────────────────────────────
if [[ -n "${job_id_list}" ]]; then
    echo "[$(date)] Step 3: waiting for jobs"
    while squeue --jobs="$(echo "$job_id_list" | tr ' ' ',' | sed 's/^,//')" -h 2>/dev/null | grep -q .; do
        sleep 60
    done
    echo "[$(date)] All pipeline jobs finished"
fi

# ── Step 4: Clean BAMs (WES BAMs are small, but keep disk tidy) ──────────
echo "[$(date)] Step 4: cleaning downloaded BAMs"
for sample in "${samples[@]}"; do
    n_done=$(find "$BATCH_ROOT" -path "*/$sample/results/**/*_optitype.txt" 2>/dev/null | wc -l)
    if [[ "$n_done" -gt 0 ]]; then
        find "$BAM_ROOT" -maxdepth 1 -name "${sample}*.bam" -delete 2>/dev/null || true
        find "$BAM_ROOT" -maxdepth 1 -name "${sample}*.bai" -delete 2>/dev/null || true
        echo "  Cleaned BAMs for $sample"
    else
        echo "  WARNING: $sample has no results — keeping BAM for inspection"
    fi
done

echo "[$(date)] Batch complete. Disk:"
df -h /scratch/project_2008084/ | tail -1
