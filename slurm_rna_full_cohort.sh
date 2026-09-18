#!/bin/bash
#SBATCH --job-name=pihla_rna_cohort
#SBATCH --account=project_2008084
#SBATCH --partition=small
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --output=/scratch/project_2008084/hla_calibration/logs/rna_cohort_%j.out
#SBATCH --error=/scratch/project_2008084/hla_calibration/logs/rna_cohort_%j.err

# Submits per-sample RNA pipeline jobs for all samples in the download index.
# Run AFTER slurm_download_rna.sh has completed successfully.
# Usage:
#   sbatch slurm_rna_full_cohort.sh                  # submit immediately
#   sbatch --dependency=afterok:<DL_JOB_ID> slurm_rna_full_cohort.sh  # after download job

set -euo pipefail

INDEX=/scratch/project_2008084/hla_calibration/rna/index/sample_fastq_urls.tsv
INPUT_ROOT=/scratch/project_2008084/hla_calibration/rna/fastqs
OUTPUT_ROOT=/scratch/project_2008084/hla_calibration/rna_batches
WORK_ROOT=/scratch/project_2008084/hla_calibration/work
PIPELINE_DIR=/scratch/project_2008084/pihla-publish

mkdir -p /scratch/project_2008084/hla_calibration/logs

echo "[$(date)] Submitting RNA pipeline jobs for all samples"
submitted=0
skipped=0

while IFS=$'\t' read -r sample_id _rest; do
  [[ -z "${sample_id}" ]] && continue

  r1="${INPUT_ROOT}/${sample_id}_R1.fastq.gz"
  r2="${INPUT_ROOT}/${sample_id}_R2.fastq.gz"

  if [[ ! -f "${r1}" || ! -f "${r2}" ]]; then
    echo "  SKIP ${sample_id}: FASTQs not found (${r1})"
    (( skipped++ )) || true
    continue
  fi

  # Rate-limit: wait if >= 15 RNA sample jobs already running/pending
  while [[ $(squeue -u ozcanumu -h --name=pihla_rna_sample 2>/dev/null | wc -l) -ge 15 ]]; do
    sleep 60
  done

  job_id=$(sbatch --parsable \
    "${PIPELINE_DIR}/slurm_rna_sample.sh" \
    "${sample_id}" \
    "${INPUT_ROOT}" \
    "${OUTPUT_ROOT}" \
    "${WORK_ROOT}")

  echo "  Submitted ${sample_id} → job ${job_id}"
  (( submitted++ )) || true
done < "${INDEX}"

echo "[$(date)] Done: ${submitted} jobs submitted, ${skipped} samples skipped (missing FASTQs)"
