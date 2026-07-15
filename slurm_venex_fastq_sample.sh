#!/bin/bash
#SBATCH --job-name=pihla_venex_fq
#SBATCH --account=project_2008084
#SBATCH --partition=small
#SBATCH --time=36:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --output=/scratch/project_2008084/hla_calibration/logs/venex_fq_%x_%j.out
#SBATCH --error=/scratch/project_2008084/hla_calibration/logs/venex_fq_%x_%j.err

set -euo pipefail

if [[ $# -lt 3 ]]; then
  echo "Usage: $0 <sample_id> <fastq1> <fastq2> [output_root] [work_root]" >&2
  exit 1
fi

sample_id="$1"
fastq1="$2"
fastq2="$3"
output_root="${4:-/scratch/project_2008084/hla_calibration/venex/results/wgs_fastq}"
work_root="${5:-/scratch/project_2008084/hla_calibration/venex/work/wgs_fastq}"

module load nextflow

mkdir -p /scratch/project_2008084/hla_calibration/logs
mkdir -p "${output_root}/${sample_id}"
mkdir -p "${work_root}"

if [[ ! -f "${fastq1}" || ! -f "${fastq2}" ]]; then
  echo "Missing FASTQ pair for ${sample_id}: ${fastq1} / ${fastq2}" >&2
  exit 1
fi

samplesheet="${output_root}/${sample_id}/venex_${sample_id}_fastq_samplesheet.csv"
cat > "${samplesheet}" <<CSV
sample_id,fastq_1,fastq_2
${sample_id},${fastq1},${fastq2}
CSV

LAUNCH_DIR="${output_root}/${sample_id}/nxf_launch"
mkdir -p "${LAUNCH_DIR}"
cd "${LAUNCH_DIR}"

nextflow run /scratch/project_2008084/pihla-publish/main.nf \
  -params-file /scratch/project_2008084/pihla-publish/conf/puhti_params.yaml \
  -c /scratch/project_2008084/pihla-publish/conf/wgs_disk_fix.config \
  -resume \
  -name pihla_venex_fq_${sample_id}_${SLURM_JOB_ID} \
  -w "${work_root}/${sample_id}" \
  --input_samplesheet "${samplesheet}" \
  --input_type fastq \
  --tools optitype,t1k,hlahd \
  --optitype_seq_type dna \
  --seq_type dna \
  --run_modality wgs \
  --outdir "${output_root}/${sample_id}/results" \
  --slurm_account project_2008084 \
  --singularity_cache_dir /scratch/project_2008084/hla_references/singularity_cache/containers \
  -profile puhti,singularity
