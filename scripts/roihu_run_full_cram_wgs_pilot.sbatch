#!/bin/bash
#SBATCH --job-name=cc_fullcram_callers
#SBATCH --account=project_2008084
#SBATCH --partition=small
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --array=1-2%1
#SBATCH --output=/scratch/project_2008084/champhla_cc_confirmation/runs/full_cram_wgs_pilot/logs/callers_%A_%a.out
#SBATCH --error=/scratch/project_2008084/champhla_cc_confirmation/runs/full_cram_wgs_pilot/logs/callers_%A_%a.err

# Runs only the frozen five-caller WGS panel. This stage has no truth input and
# writes exclusively inside the isolated confirmation project.
set -eo pipefail

ROOT=/scratch/project_2008084/champhla_cc_confirmation
PIPELINE=/scratch/project_2008084/pihla-publish
MANIFEST=${ROOT}/configs/full_cram_wgs_pilot.tsv
RUN_ROOT=${ROOT}/runs/full_cram_wgs_pilot
INPUT_ROOT=${RUN_ROOT}/inputs
OUTPUT_ROOT=${RUN_ROOT}/caller_outputs
WORK_ROOT=${RUN_ROOT}/work
LOG_ROOT=${RUN_ROOT}/logs

export CSC_ENV_INIT_NON_INTERACTIVE=yes
export NXF_OPTS='-Xms256m -Xmx4g'
source /etc/profile.d/zz-csc-env.sh
module load bio-apps/v202603
module load samtools/1.21
module load nextflow/25.10.2-standalone
set -u
mkdir -p "${OUTPUT_ROOT}" "${WORK_ROOT}" "${LOG_ROOT}"

line_number=$((SLURM_ARRAY_TASK_ID + 1))
IFS=$'\t' read -r sample_id cram_url expected_md5 source_record < <(sed -n "${line_number}p" "${MANIFEST}")
input_bam=${INPUT_ROOT}/${sample_id}/${sample_id}.mate_aware_hla.bam
if [[ ! -f "${INPUT_ROOT}/${sample_id}/EXTRACTION_COMPLETE" || ! -s "${input_bam}" || ! -s "${input_bam}.bai" ]]; then
  echo "Validated extracted input is absent for ${sample_id}" >&2
  exit 2
fi

sample_root=${OUTPUT_ROOT}/${sample_id}
mkdir -p "${sample_root}"
completion_marker=${sample_root}/CALLERS_COMPLETE
if [[ -e "${completion_marker}" ]]; then
  mv "${completion_marker}" "${completion_marker}.superseded.${SLURM_JOB_ID}"
fi
{
  nextflow -version
  samtools --version
  apptainer --version
  module list
} > "${sample_root}/caller_environment.txt" 2>&1
samplesheet=${sample_root}/${sample_id}.samplesheet.csv
printf 'sample_id,bam_path\n%s,%s\n' "${sample_id}" "${input_bam}" > "${samplesheet}"
sha256sum "${input_bam}" "${input_bam}.bai" "${samplesheet}" > "${sample_root}/pre_run_sha256.txt"
sha256sum \
  /scratch/project_2008084/hla_references/singularity_cache/containers/hlahd.sif \
  /scratch/project_2008084/hla_references/singularity_cache/containers/kourami.sif \
  /scratch/project_2008084/hla_references/singularity_cache/containers/optitype.sif \
  /scratch/project_2008084/hla_references/singularity_cache/containers/t1k.sif \
  > "${sample_root}/container_sha256.txt"
sha256sum \
  "${PIPELINE}/main.nf" \
  "${PIPELINE}/nextflow.config" \
  "${PIPELINE}/conf/roihu_params.yaml" \
  "${PIPELINE}/conf/wgs_disk_fix.config" \
  "${ROOT}/configs/wgs_pilot_resource_fix.config" \
  "${ROOT}/models/frozen_raw_cc_policy.json" \
  /projappl/project_2008084/SpecHLAx/script/whole/SpecHLA.sh \
  > "${sample_root}/software_bundle_sha256.txt"
git -C /projappl/project_2008084/SpecHLAx rev-parse HEAD \
  > "${sample_root}/spechla_git_commit.txt" 2>/dev/null || printf 'not_a_git_checkout\n' \
  > "${sample_root}/spechla_git_commit.txt"

launch_root=${sample_root}/nxf_launch
mkdir -p "${launch_root}"
cd "${launch_root}"

nextflow run "${PIPELINE}/main.nf" \
  -resume \
  -name "cc_fullcram_${sample_id}_${SLURM_JOB_ID}" \
  -w "${WORK_ROOT}/${sample_id}" \
  -params-file "${PIPELINE}/conf/roihu_params.yaml" \
  -c "${PIPELINE}/conf/wgs_disk_fix.config" \
  -c "${ROOT}/configs/wgs_pilot_resource_fix.config" \
  --input_samplesheet "${samplesheet}" \
  --input_type bam \
  --tools spechla,hlahd,optitype,kourami,t1k \
  --hlahd_db /app/hlahd.1.4.0 \
  --extract_hla_region false \
  --spechla_exon_only 0 \
  --seq_type dna \
  --run_modality wgs \
  --outdir "${sample_root}/results" \
  --slurm_account project_2008084 \
  --use_local_spechla true \
  --spechla_path /projappl/project_2008084/SpecHLAx \
  --singularity_cache_dir /scratch/project_2008084/hla_references/singularity_cache/containers \
  -profile roihu,singularity

validation=${sample_root}/caller_output_validation.tsv
printf 'caller\tnative_summary\tsha256\n' > "${validation}"
for caller_suffix in spechla hlahd optitype kourami t1k; do
  mapfile -t matches < <(find "${sample_root}/results" -type f \
    -name "${sample_id}_${caller_suffix}.txt" -size +0c -print)
  if [[ ${#matches[@]} -ne 1 ]]; then
    echo "Expected exactly one nonempty ${caller_suffix} summary, found ${#matches[@]}" >&2
    exit 3
  fi
  if grep -q '^# No results generated' "${matches[0]}"; then
    echo "${caller_suffix} emitted a no-results placeholder" >&2
    exit 4
  fi
  if ! grep -Eq '(^|[^[:alnum:]])(HLA-)?[ABC]\*[0-9A-Za-z]+:[0-9A-Za-z]+' "${matches[0]}"; then
    echo "${caller_suffix} summary contains no parseable class-I two-field allele" >&2
    exit 5
  fi
  printf '%s\t%s\t%s\n' "${caller_suffix}" "${matches[0]}" \
    "$(sha256sum "${matches[0]}" | cut -d ' ' -f 1)" >> "${validation}"
done
find "${sample_root}/results" -type f -print0 | sort -z | xargs -0 sha256sum > "${sample_root}/prediction_files_sha256.txt"
touch "${completion_marker}"
