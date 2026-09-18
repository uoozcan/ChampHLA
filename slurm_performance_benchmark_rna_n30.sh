#!/bin/bash
#SBATCH --job-name=perf_bench_rna_n30
#SBATCH --account=project_2008084
#SBATCH --partition=small
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --array=1-30%10
#SBATCH --output=/scratch/project_2008084/hla_calibration/logs/perf_bench_rna_%A_%a.out
#SBATCH --error=/scratch/project_2008084/hla_calibration/logs/perf_bench_rna_%A_%a.err

# RNA-seq computational performance benchmark — n=30 samples, HLA-region-extracted BAMs.
# One Nextflow job per sample so each produces one execution_trace.txt.
#
# Prerequisites:
#   - Run slurm_download_rna_hla_bams.sh first to download HLA-region BAMs
#   - BAMs must exist at analysis/performance_benchmark_n30/rna_hla_bams/
#
# After all 30 array jobs complete, run:
#   python3.11 bin/parse_nextflow_trace_timing.py \
#     --rna-traces analysis/performance_benchmark_n30/runs/rna_hla/*/pipeline_info/execution_trace.txt \
#     --wgs-traces analysis/performance_benchmark_n30/runs/wgs/*/pipeline_info/execution_trace.txt \
#     --n-per-modality 30 --out-dir analysis/performance_benchmark_n30/tables/

set -euo pipefail

SAMPLE_LIST=/scratch/project_2008084/pihla-publish/analysis/performance_benchmark_n30/rna_n30_sample_list.txt
SAMPLESHEET_DIR=/scratch/project_2008084/pihla-publish/analysis/performance_benchmark_n30/samplesheets/rna
OUTPUT_ROOT=/scratch/project_2008084/pihla-publish/analysis/performance_benchmark_n30/runs/rna_hla
WORK_ROOT=/scratch/project_2008084/hla_calibration/work/perf_bench_rna

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
  -resume \
  -name "perf_bench_rna_${sample_id}_${SLURM_JOB_ID}" \
  -w "${WORK_ROOT}/work_${sample_id}" \
  -params-file /scratch/project_2008084/pihla-publish/conf/puhti_params.yaml \
  --input_samplesheet "${samplesheet}" \
  --input_type bam \
  --tools spechla,hlahd,optitype,t1k,arcashla,seq2hla \
  --spechla_exon_only 0 \
  --seq_type rna \
  --run_modality rnaseq \
  --optitype_seq_type rna \
  --outdir "${outdir}" \
  --use_local_spechla true \
  --spechla_path /projappl/project_2008084/SpecHLAx \
  -profile puhti,singularity

echo "[DONE] ${sample_id} at $(date)"
echo "[TRACE] $(ls ${outdir}/pipeline_info/execution_trace.txt 2>/dev/null || echo 'no trace found')"
