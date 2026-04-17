#!/bin/bash
#SBATCH --job-name=pihla_wgs_dl
#SBATCH --account=project_2008084
#SBATCH --partition=small
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --array=1-8
#SBATCH --output=/scratch/project_2008084/hla_calibration/logs/wgs_dl_%A_%a.out
#SBATCH --error=/scratch/project_2008084/hla_calibration/logs/wgs_dl_%A_%a.err

# Downloads HLA-region BAMs for the 8 WGS samples missing from the previous run.
# Streams chr6:28000000-34000000 from NYGC 30x CRAM files at EBI.

set -euo pipefail

module load samtools 2>/dev/null || true

SAMPLES_FILE=/scratch/project_2008084/hla_calibration/wgs/wgs_samples_missing.txt
OUT_DIR=/scratch/project_2008084/hla_calibration/wgs/bams_new
URL_MAP=/scratch/project_2008084/hla_calibration/wgs/index/sample_cram_urls.tsv
REF_CACHE=/scratch/project_2008084/hla_calibration/wgs/cram_ref_cache

export REF_PATH="https://www.ebi.ac.uk/ena/cram/md5/%s"
export REF_CACHE="${REF_CACHE}/%2s/%2s/%s"
mkdir -p "${OUT_DIR}" "${REF_CACHE}"

sample_id=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${SAMPLES_FILE}")
[[ -z "${sample_id}" ]] && { echo "[ERROR] Empty sample for task ${SLURM_ARRAY_TASK_ID}"; exit 1; }

out_bam="${OUT_DIR}/${sample_id}_hla.bam"

if [[ -f "${out_bam}" ]] && [[ -f "${out_bam}.bai" ]]; then
  echo "[SKIP] Already done: ${sample_id}"
  exit 0
fi

CRAM_URL=$(awk -v s="${sample_id}" '$1==s{print $3}' "${URL_MAP}")
[[ -z "${CRAM_URL}" ]] && { echo "[ERROR] No CRAM URL for ${sample_id}" >&2; exit 1; }

echo "[INFO] Sample: ${sample_id}"
echo "[INFO] CRAM:   ${CRAM_URL}"
echo "[INFO] Date:   $(date)"

TMP="${out_bam}.tmp"
SUCCESS=0
for REGION in "chr6:28000000-34000000" "6:28000000-34000000"; do
  if samtools view -b -h -@ "${SLURM_CPUS_PER_TASK}" -o "${TMP}" "${CRAM_URL}" "${REGION}" 2>/dev/null; then
    NREADS=$(samtools view -c "${TMP}" 2>/dev/null || echo 0)
    if [[ "${NREADS}" -gt 0 ]]; then
      mv "${TMP}" "${out_bam}"
      samtools index "${out_bam}"
      echo "[OK] ${sample_id}: ${NREADS} HLA reads"
      SUCCESS=1
      break
    fi
  fi
  rm -f "${TMP}"
done

[[ "${SUCCESS}" -eq 0 ]] && { echo "[ERROR] Failed to extract reads for ${sample_id}" >&2; exit 1; }
