#!/bin/bash
#SBATCH --job-name=pihla_wes_sample
#SBATCH --account=project_2008084
#SBATCH --partition=small
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --output=/scratch/project_2008084/hla_calibration/logs/wes_sample_%x_%j.out
#SBATCH --error=/scratch/project_2008084/hla_calibration/logs/wes_sample_%x_%j.err

set -euo pipefail

if [[ $# -lt 1 ]]; then
  echo "Usage: $0 <sample_id> [input_root] [output_root] [work_root]" >&2
  exit 1
fi

sample_id="$1"
input_root="${2:-/scratch/project_2008084/hla_calibration/wes_3sample_input}"
output_root="${3:-/scratch/project_2008084/hla_calibration/wes_batches}"
work_root="${4:-/scratch/project_2008084/hla_calibration/work}"

module load nextflow
module load samtools

mkdir -p /scratch/project_2008084/hla_calibration/logs
mkdir -p "${output_root}/${sample_id}"
mkdir -p "${work_root}"

bam_path="$(find "${input_root}" -maxdepth 1 -type f -name "${sample_id}*.bam" | head -n 1)"
if [[ -z "${bam_path}" ]]; then
  echo "No BAM found for ${sample_id} under ${input_root}" >&2
  exit 1
fi

[[ -f "${bam_path}.bai" ]] || samtools index -@ 8 "${bam_path}"

samplesheet="${output_root}/${sample_id}/wes_${sample_id}_samplesheet.csv"
cat > "${samplesheet}" <<CSV
sample_id,bam_path
${sample_id},${bam_path}
CSV

# Use a per-sample launch directory so concurrent jobs don't share .nextflow/cache locks
LAUNCH_DIR="${output_root}/${sample_id}/nxf_launch"
mkdir -p "${LAUNCH_DIR}"
cd "${LAUNCH_DIR}"

nextflow run /scratch/project_2008084/pihla-publish/main.nf \
  -params-file /scratch/project_2008084/pihla-publish/conf/puhti_params.yaml \
  -resume \
  -name pihla_wes_${sample_id}_${SLURM_JOB_ID} \
  -w "${work_root}/wes_${sample_id}" \
  --input_samplesheet "${samplesheet}" \
  --input_type bam \
  --tools spechla,hlahd,optitype,polysolver,kourami,t1k,arcashla \
  --spechla_exon_only 1 \
  --seq_type dna \
  --run_modality wes \
  --outdir "${output_root}/${sample_id}/results" \
  --slurm_account project_2008084 \
  --use_local_spechla true \
  --spechla_path /projappl/project_2008084/SpecHLAx \
  --singularity_cache_dir /scratch/project_2008084/hla_references/singularity_cache/containers \
  -profile puhti,singularity
