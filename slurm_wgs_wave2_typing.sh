#!/bin/bash
#SBATCH --job-name=pihla_wgs_type_wave2
#SBATCH --account=project_2008084
#SBATCH --partition=small
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --array=1-49%10
#SBATCH --output=/scratch/project_2008084/hla_calibration/logs/wgs_type_wave2_%A_%a.out
#SBATCH --error=/scratch/project_2008084/hla_calibration/logs/wgs_type_wave2_%A_%a.err

# Runs PIHLA Nextflow WGS HLA typing for the 49 wave-2 samples (those present in
# WES/RNA batches but not yet in wgs_batches). At most 10 tasks run concurrently
# to avoid disk space exhaustion.
# Run AFTER slurm_download_wgs_wave2.sh completes.

set -euo pipefail

SAMPLES_FILE=/scratch/project_2008084/hla_calibration/wgs/wgs_wave2_samples.txt
INPUT_DIR=/scratch/project_2008084/hla_calibration/wgs/bams_new
OUTPUT_ROOT=/scratch/project_2008084/hla_calibration/wgs_batches
WORK_ROOT=/scratch/project_2008084/hla_calibration/work

module load nextflow
module load samtools 2>/dev/null || true

mkdir -p /scratch/project_2008084/hla_calibration/logs "${WORK_ROOT}"

sample_id=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${SAMPLES_FILE}")
[[ -z "${sample_id}" ]] && { echo "[ERROR] Empty sample for task ${SLURM_ARRAY_TASK_ID}"; exit 1; }

bam_path="${INPUT_DIR}/${sample_id}_hla.bam"
if [[ ! -f "${bam_path}" ]]; then
  echo "[ERROR] BAM not found: ${bam_path}" >&2
  exit 1
fi

# Skip if already typed (check for at least 6 tool result dirs)
result_dir="${OUTPUT_ROOT}/${sample_id}/results/${sample_id}"
if [[ -d "${result_dir}" ]]; then
  n_tools=$(ls -d "${result_dir}"/*/2>/dev/null | wc -l)
  if [[ "${n_tools}" -ge 6 ]]; then
    echo "[SKIP] ${sample_id} already has ${n_tools} tool results"
    exit 0
  fi
fi

[[ -f "${bam_path}.bai" ]] || samtools index -@ 8 "${bam_path}"

samplesheet="${OUTPUT_ROOT}/${sample_id}/wgs_${sample_id}_samplesheet.csv"
mkdir -p "${OUTPUT_ROOT}/${sample_id}"
cat > "${samplesheet}" <<CSV
sample_id,bam_path
${sample_id},${bam_path}
CSV

LAUNCH_DIR="${OUTPUT_ROOT}/${sample_id}/nxf_launch"
mkdir -p "${LAUNCH_DIR}"
cd "${LAUNCH_DIR}"

nextflow run /scratch/project_2008084/pihla-publish/main.nf \
  -resume \
  -name "pihla_wgs_${sample_id}_${SLURM_JOB_ID}" \
  -w "${WORK_ROOT}/wgs_${sample_id}" \
  --input_samplesheet "${samplesheet}" \
  --input_type bam \
  --tools spechla,hlahd,optitype,polysolver,kourami,t1k,arcashla \
  --spechla_exon_only 0 \
  --seq_type dna \
  --run_modality wgs \
  --outdir "${OUTPUT_ROOT}/${sample_id}/results" \
  --slurm_account project_2008084 \
  --use_local_spechla true \
  --spechla_path /projappl/project_2008084/SpecHLAx \
  -profile puhti,singularity

echo "[DONE] ${sample_id} at $(date)"
