#!/bin/bash
#SBATCH --job-name=seq2hla_rna_bench
#SBATCH --account=project_2008084
#SBATCH --partition=small
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --array=1-30%10
#SBATCH --output=/scratch/project_2008084/hla_calibration/logs/seq2hla_rna_%A_%a.out
#SBATCH --error=/scratch/project_2008084/hla_calibration/logs/seq2hla_rna_%A_%a.err

# Seq2HLA-only RNA benchmark — runs Seq2HLA on the same 30 RNA samples
# used in the main RNA performance benchmark.
# Output goes to a SEPARATE directory so existing traces are preserved.

set -euo pipefail

SAMPLE_LIST=/scratch/project_2008084/pihla-publish/analysis/performance_benchmark_n30/rna_n30_sample_list.txt
SAMPLESHEET_DIR=/scratch/project_2008084/pihla-publish/analysis/performance_benchmark_n30/samplesheets/rna
OUTPUT_ROOT=/scratch/project_2008084/pihla-publish/analysis/performance_benchmark_n30/runs/rna_hla_seq2hla
WORK_ROOT=/scratch/project_2008084/hla_calibration/work/seq2hla_rna_bench

module load nextflow
module load samtools 2>/dev/null || true

mkdir -p /scratch/project_2008084/hla_calibration/logs "${WORK_ROOT}"

sample_id=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${SAMPLE_LIST}")
[[ -z "${sample_id}" ]] && { echo "[ERROR] Empty sample for task ${SLURM_ARRAY_TASK_ID}"; exit 1; }

samplesheet="${SAMPLESHEET_DIR}/${sample_id}.csv"
if [[ ! -f "${samplesheet}" ]]; then
  echo "[ERROR] Samplesheet not found: ${samplesheet}" >&2
  exit 1
fi

bam_file="/scratch/project_2008084/pihla-publish/analysis/performance_benchmark_n30/rna_hla_bams/${sample_id}_rna_hla.bam"
if [[ ! -f "${bam_file}" ]]; then
  echo "[ERROR] HLA-region BAM not found: ${bam_file}" >&2
  exit 1
fi

outdir="${OUTPUT_ROOT}/${sample_id}"
mkdir -p "${outdir}"

LAUNCH_DIR="${WORK_ROOT}/nxf_${sample_id}"
mkdir -p "${LAUNCH_DIR}"
cd "${LAUNCH_DIR}"

echo "[START] ${sample_id} at $(date)"

nextflow run /scratch/project_2008084/pihla-publish/main.nf \
  -name "seq2hla_rna_${sample_id}_${SLURM_JOB_ID}" \
  -w "${WORK_ROOT}/work_${sample_id}" \
  -params-file /scratch/project_2008084/pihla-publish/conf/puhti_params.yaml \
  --input_samplesheet "${samplesheet}" \
  --input_type bam \
  --tools seq2hla \
  --seq_type rna \
  --run_modality rnaseq \
  --outdir "${outdir}" \
  -profile puhti,singularity

echo "[DONE] ${sample_id} at $(date)"
echo "[TRACE] $(ls ${outdir}/pipeline_info/execution_trace.txt 2>/dev/null || echo 'no trace found')"
