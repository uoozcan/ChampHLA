#!/bin/bash
#SBATCH --job-name=pihla_rna_sample
#SBATCH --account=project_2008084
#SBATCH --partition=small
#SBATCH --time=16:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --output=/scratch/project_2008084/hla_calibration/logs/rna_sample_%x_%j.out
#SBATCH --error=/scratch/project_2008084/hla_calibration/logs/rna_sample_%x_%j.err

set -euo pipefail

if [[ $# -lt 1 ]]; then
  echo "Usage: $0 <sample_id> [input_root] [output_root] [work_root]" >&2
  exit 1
fi

sample_id="$1"
input_root="${2:-/scratch/project_2008084/hla_calibration/rna_3sample_input_named}"
output_root="${3:-/scratch/project_2008084/hla_calibration/rna_batches}"
work_root="${4:-/scratch/project_2008084/hla_calibration/work}"

module load nextflow

mkdir -p /scratch/project_2008084/hla_calibration/logs
mkdir -p "${output_root}/${sample_id}"
mkdir -p "${work_root}"
mkdir -p "${input_root}"

fastq1="${input_root}/${sample_id}_R1.fastq.gz"
fastq2="${input_root}/${sample_id}_R2.fastq.gz"
if [[ ! -f "${fastq1}" || ! -f "${fastq2}" ]]; then
  echo "Missing FASTQ pair for ${sample_id} under ${input_root}" >&2
  exit 1
fi

samplesheet="${output_root}/${sample_id}/rna_${sample_id}_samplesheet.csv"
cat > "${samplesheet}" <<CSV
sample_id,fastq_1,fastq_2
${sample_id},${fastq1},${fastq2}
CSV

# Use a per-sample launch directory so concurrent jobs don't share .nextflow/cache locks
LAUNCH_DIR="${output_root}/${sample_id}/nxf_launch"
mkdir -p "${LAUNCH_DIR}"
cd "${LAUNCH_DIR}"

nextflow run /scratch/project_2008084/pihla-publish/main.nf \
  -resume \
  -name pihla_rna_${sample_id}_${SLURM_JOB_ID} \
  -w "${work_root}/rna_${sample_id}" \
  --input_samplesheet "${samplesheet}" \
  --input_type fastq \
  --tools arcashla,optitype,seq2hla,t1k,spechla,hlahd \
  --optitype_seq_type rna \
  --seq_type rna \
  --run_modality rnaseq \
  --outdir "${output_root}/${sample_id}/results" \
  --slurm_account project_2008084 \
  --use_local_spechla true \
  --spechla_path /projappl/project_2008084/SpecHLAx \
  -profile puhti,singularity
