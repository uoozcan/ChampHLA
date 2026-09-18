#!/bin/bash
#SBATCH --job-name=pihla_wgs_b11
#SBATCH --account=project_2008084
#SBATCH --partition=small
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --output=/scratch/project_2008084/pihla-publish/analysis/logs/wgs_batch11_%j.out
#SBATCH --error=/scratch/project_2008084/pihla-publish/analysis/logs/wgs_batch11_%j.err

# Downloads HLA-region BAM from NYGC 30x CRAM, types HLA, then cleans the BAM.
# Usage: sbatch slurm_wgs_batch11.sh <sample1> [sample2 ...]

set -euo pipefail

INDEX=/scratch/project_2008084/hla_calibration/wgs/index/sample_cram_urls_batch3.tsv
BAM_ROOT=/scratch/project_2008084/hla_calibration/wgs/bams_new
BATCH_ROOT=/scratch/project_2008084/hla_calibration/wgs_batches
WORK_ROOT=/scratch/project_2008084/hla_calibration/work
PIPELINE_DIR=/scratch/project_2008084/pihla-publish
REF=/scratch/project_2008084/references/GRCh38.fa
REGION_CHR="chr6:28000000-34000000"
REGION_NUM="6:28000000-34000000"

mkdir -p "$BAM_ROOT" /scratch/project_2008084/hla_calibration/logs

samples=("$@")
[[ ${#samples[@]} -eq 0 ]] && { echo "Usage: $0 <sample1> [...]" >&2; exit 1; }

echo "[$(date)] WGS batch: ${#samples[@]} samples: ${samples[*]}"
module load samtools 2>/dev/null || true
SAM=samtools

# ── Step 1: Download HLA region from CRAM ────────────────────────────────
echo "[$(date)] Step 1: downloading HLA-region BAMs from NYGC CRAMs"
for sample in "${samples[@]}"; do
    bam_count=$(find "$BAM_ROOT" -maxdepth 1 -name "${sample}*.bam" 2>/dev/null | wc -l)
    [[ "$bam_count" -gt 0 ]] && echo "  SKIP $sample (BAM exists)" && continue

    line=$(grep "^${sample}"$'\t' "$INDEX" 2>/dev/null | head -1 || true)
    [[ -z "$line" ]] && echo "  SKIP $sample: not in index" && continue

    cram_url=$(echo "$line" | cut -f3)
    out_bam="$BAM_ROOT/${sample}.hla.bam"
    tmp_bam="${out_bam}.tmp"

    echo "  Downloading HLA region for $sample..."
    downloaded=0
    for region in "$REGION_CHR" "$REGION_NUM"; do
        if "$SAM" view -b -T "$REF" -o "$tmp_bam" "$cram_url" "$region" 2>/dev/null && \
           [[ -s "$tmp_bam" ]]; then
            mv "$tmp_bam" "$out_bam"
            "$SAM" index "$out_bam"
            nreads=$("$SAM" view -c "$out_bam" 2>/dev/null || echo "?")
            echo "  Downloaded $sample: $nreads reads ($(du -sh "$out_bam" | cut -f1))"
            downloaded=1
            break
        fi
    done
    rm -f "$tmp_bam"
    [[ "$downloaded" -eq 0 ]] && echo "  WARNING: download failed for $sample"
done
echo "[$(date)] Downloads complete"

# ── Step 2: Submit pipeline jobs ──────────────────────────────────────────
echo "[$(date)] Step 2: submitting pipeline jobs"
job_id_list=""
for sample in "${samples[@]}"; do
    bam_count=$(find "$BAM_ROOT" -maxdepth 1 -name "${sample}*.bam" 2>/dev/null | wc -l)
    [[ "$bam_count" -eq 0 ]] && echo "  SKIP $sample: no BAM" && continue

    n_done=$(find "$BATCH_ROOT" -path "*/results/${sample}/optitype/*" 2>/dev/null | wc -l)
    [[ "$n_done" -gt 0 ]] && echo "  SKIP $sample: already typed" && continue

    job_id=$(sbatch --parsable \
        "$PIPELINE_DIR/slurm_wgs_sample.sh" \
        "$sample" "$BAM_ROOT" "$BATCH_ROOT" "$WORK_ROOT")
    job_id_list="${job_id_list} ${job_id}"
    echo "  Submitted $sample → job $job_id"
done

# ── Step 3: Wait ──────────────────────────────────────────────────────────
if [[ -n "${job_id_list}" ]]; then
    ids=$(echo "$job_id_list" | tr ' ' ',' | sed 's/^,//')
    echo "[$(date)] Step 3: waiting for jobs: $ids"
    while squeue --jobs="$ids" -h 2>/dev/null | grep -q .; do
        sleep 60
    done
    echo "[$(date)] Pipeline jobs finished"
fi

# ── Step 4: Clean HLA BAMs ────────────────────────────────────────────────
echo "[$(date)] Step 4: cleaning HLA BAMs"
for sample in "${samples[@]}"; do
    n_done=$(find "$BATCH_ROOT" -path "*/results/${sample}/optitype/*" 2>/dev/null | wc -l)
    if [[ "$n_done" -gt 0 ]]; then
        find "$BAM_ROOT" -maxdepth 1 -name "${sample}*.bam" -delete 2>/dev/null || true
        find "$BAM_ROOT" -maxdepth 1 -name "${sample}*.bai" -delete 2>/dev/null || true
        # Remove Nextflow work directory for this sample
        rm -rf "${WORK_ROOT}/wgs_${sample}" 2>/dev/null || true
        rm -rf "${WORK_ROOT}/wgs_arcashla_${sample}" 2>/dev/null || true
        echo "  Cleaned $sample (BAM + work dir)"
    else
        echo "  WARNING: $sample has no results — keeping BAM"
    fi
done

echo "[$(date)] Batch complete. Disk:"
df -h /scratch/project_2008084/ | tail -1
