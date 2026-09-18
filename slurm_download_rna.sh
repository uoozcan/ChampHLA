#!/bin/bash
#SBATCH --job-name=pihla_dl_rna
#SBATCH --account=project_2008084
#SBATCH --partition=small
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --array=1-50
#SBATCH --output=/scratch/project_2008084/hla_calibration/logs/dl_rna_%a_%j.out
#SBATCH --error=/scratch/project_2008084/hla_calibration/logs/dl_rna_%a_%j.err

# Downloads GEUVADIS RNA-seq FASTQ pairs from EBI for all samples in the index file.
# Uses SLURM array: task N downloads the Nth sample.

set -euo pipefail

INDEX=/scratch/project_2008084/hla_calibration/rna/index/sample_fastq_urls.tsv
OUTDIR=/scratch/project_2008084/hla_calibration/rna/fastqs

mkdir -p "${OUTDIR}"
mkdir -p /scratch/project_2008084/hla_calibration/logs

# Read the Nth line (1-based, no header in this file)
line=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${INDEX}")
if [[ -z "${line}" ]]; then
  echo "No data for array task ${SLURM_ARRAY_TASK_ID}" >&2
  exit 0
fi

sample_id=$(echo "${line}" | cut -f1)
# Primary R1/R2 are columns 2 and 3 (may contain pipe-separated alternatives; use first)
r1_raw=$(echo "${line}" | cut -f2 | cut -d'|' -f1)
r2_raw=$(echo "${line}" | cut -f3 | cut -d'|' -f1)

# Prepend https:// if the URL has no protocol
[[ "${r1_raw}" == http* || "${r1_raw}" == ftp* ]] && r1_url="${r1_raw}" || r1_url="https://${r1_raw}"
[[ "${r2_raw}" == http* || "${r2_raw}" == ftp* ]] && r2_url="${r2_raw}" || r2_url="https://${r2_raw}"

out_r1="${OUTDIR}/${sample_id}_R1.fastq.gz"
out_r2="${OUTDIR}/${sample_id}_R2.fastq.gz"

echo "[$(date)] Task ${SLURM_ARRAY_TASK_ID}: ${sample_id}"
echo "  R1: ${r1_url}"
echo "  R2: ${r2_url}"

if [[ -f "${out_r1}" ]]; then
  echo "  R1 already exists, skipping download"
else
  wget -q --no-check-certificate -c -O "${out_r1}" "${r1_url}"
  echo "  R1 done ($(du -sh "${out_r1}" | cut -f1))"
fi

if [[ -f "${out_r2}" ]]; then
  echo "  R2 already exists, skipping download"
else
  wget -q --no-check-certificate -c -O "${out_r2}" "${r2_url}"
  echo "  R2 done ($(du -sh "${out_r2}" | cut -f1))"
fi

echo "[$(date)] ${sample_id} complete"
