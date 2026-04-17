#!/bin/bash
#SBATCH --job-name=pihla_wgs_cohort
#SBATCH --account=project_2008084
#SBATCH --partition=small
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --output=/scratch/project_2008084/hla_calibration/logs/wgs_cohort_%j.out
#SBATCH --error=/scratch/project_2008084/hla_calibration/logs/wgs_cohort_%j.err

# Submits per-sample WGS pipeline jobs for all samples in the CRAM index.
# Run AFTER slurm_download_wgs_cram.sh has completed.

set -euo pipefail

INDEX=/scratch/project_2008084/hla_calibration/wgs/index/sample_cram_urls.tsv
INPUT_ROOT=/scratch/project_2008084/hla_calibration/wgs/bams_new
OUTPUT_ROOT=/scratch/project_2008084/hla_calibration/wgs_batches
WORK_ROOT=/scratch/project_2008084/hla_calibration/work
PIPELINE_DIR=/scratch/project_2008084/pihla-publish

mkdir -p /scratch/project_2008084/hla_calibration/logs

echo "[$(date)] Submitting WGS pipeline jobs for all samples"
submitted=0
skipped=0

while IFS=$'\t' read -r sample_id _rest; do
  [[ -z "${sample_id}" ]] && continue

  # Skip if results already exist from old or new pipeline
  if [[ -d "/scratch/project_2008084/hla_calibration/wgs/results/${sample_id}" ]] || \
     [[ -d "${OUTPUT_ROOT}/${sample_id}/results" ]]; then
    echo "  SKIP ${sample_id}: results already exist"
    (( skipped++ )) || true
    continue
  fi

  # Check that BAM exists
  bam_count=$(find "${INPUT_ROOT}" -maxdepth 1 -name "${sample_id}*.bam" 2>/dev/null | wc -l)
  if [[ "${bam_count}" -eq 0 ]]; then
    echo "  SKIP ${sample_id}: no BAM found under ${INPUT_ROOT}"
    (( skipped++ )) || true
    continue
  fi

  # Rate-limit: wait if >= 15 WGS sample jobs already running/pending
  while [[ $(squeue -u ozcanumu -h --name=pihla_wgs_sample 2>/dev/null | wc -l) -ge 15 ]]; do
    sleep 60
  done

  job_id=$(sbatch --parsable \
    "${PIPELINE_DIR}/slurm_wgs_sample.sh" \
    "${sample_id}" \
    "${INPUT_ROOT}" \
    "${OUTPUT_ROOT}" \
    "${WORK_ROOT}")

  echo "  Submitted ${sample_id} → job ${job_id}"
  (( submitted++ )) || true
done < "${INDEX}"

echo "[$(date)] Done: ${submitted} jobs submitted, ${skipped} samples skipped"
