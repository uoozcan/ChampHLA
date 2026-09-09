#!/bin/bash
#SBATCH --job-name=champhla_stage
#SBATCH --account=project_2008084
#SBATCH --partition=small
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G

# Stage exactly one truth-free manifest row. For WGS, create the frozen
# mate-aware HLA BAM; for WES/RNA, retain the checksummed source input.
set -euo pipefail

CODE_ROOT=${CHAMPHLA_CODE_ROOT:-/scratch/project_2008084/champhla_plurality}
source "${CHAMPHLA_PATH_CONFIG:-${CODE_ROOT}/configs/roihu_paths.env}"
manifest=${1:?usage: roihu_stage_inputs.sbatch MANIFEST}
line_number=$((SLURM_ARRAY_TASK_ID + 1))
IFS=$'\t' read -r cohort sample_id donor_id modality input_type input_uri index_uri index_checksum source_checksum reference_build read_layout independence_stratum evidence_role \
  < <(sed -n "${line_number}p" "${manifest}")

if [[ -z "${sample_id:-}" || -z "${modality:-}" ]]; then
  echo "Missing manifest row for array task ${SLURM_ARRAY_TASK_ID}" >&2
  exit 2
fi
sample_root=${CHAMPHLA_INPUT_ROOT}/${modality}/${sample_id}
mkdir -p "${sample_root}"
test ! -e "${sample_root}/STAGE_COMPLETE"

export CSC_ENV_INIT_NON_INTERACTIVE=yes
source /etc/profile.d/zz-csc-env.sh
module load bio-apps/v202603
module load samtools/1.21

fetch_one() {
  local uri=$1
  local target=$2
  if [[ "${uri}" == /* ]]; then
    ln -s "${uri}" "${target}"
  elif [[ "${uri}" == s3://* ]]; then
    aws s3 --no-sign-request cp "${uri}" "${target}"
  else
    curl --fail --location --retry 4 --output "${target}" "${uri}"
  fi
}

verify_one() {
  local expected=$1
  local target=$2
  expected=${expected#sha256:}
  if [[ ${#expected} -eq 64 ]]; then
    printf '%s  %s\n' "${expected}" "${target}" | sha256sum --check --status
  else
    expected=${expected#md5:}
    printf '%s  %s\n' "${expected}" "${target}" | md5sum --check --status
  fi
}

if [[ "${input_type}" == fastq_pair ]]; then
  IFS=';' read -r uri1 uri2 <<< "${input_uri}"
  IFS=';' read -r sha1 sha2 <<< "${source_checksum}"
  fetch_one "${uri1}" "${sample_root}/${sample_id}_R1.fastq.gz"
  fetch_one "${uri2}" "${sample_root}/${sample_id}_R2.fastq.gz"
  verify_one "${sha1}" "${sample_root}/${sample_id}_R1.fastq.gz"
  verify_one "${sha2}" "${sample_root}/${sample_id}_R2.fastq.gz"
elif [[ "${input_type}" == cram || "${input_type}" == bam ]]; then
  raw=${sample_root}/${sample_id}.${input_type}
  idx=${raw}.$([[ "${input_type}" == cram ]] && printf 'crai' || printf 'bai')
  fetch_one "${input_uri}" "${raw}"
  fetch_one "${index_uri}" "${idx}"
  verify_one "${source_checksum}" "${raw}"
  verify_one "${index_checksum}" "${idx}"
  samtools quickcheck -v "${raw}"
else
  echo "Unsupported input type: ${input_type}" >&2
  exit 3
fi

if [[ "${modality}" == wgs ]]; then
  [[ "${input_type}" == cram ]]
  reference=${CHAMPHLA_REFERENCE_ROOT}/genomes/GRCh38_full_analysis_set_plus_decoy_hla.fa
  regions=${sample_root}/grch38dh_hla_regions.bed
  awk -F '\t' 'BEGIN {OFS="\t"} $1 == "chr6" {print $1,27999999,34000000} $1 ~ /^HLA-[ABC]\*/ {print $1,0,$2}' \
    "${reference}.fai" > "${regions}"
  [[ $(wc -l < "${regions}") -ge 4 ]]
  tmp_bam=${sample_root}/${sample_id}.mate_aware_hla.unsorted.bam
  final_bam=${sample_root}/${sample_id}.mate_aware_hla.bam
  samtools view --threads "${SLURM_CPUS_PER_TASK}" --reference "${reference}" \
    --fetch-pairs --target-file "${regions}" --bam --output "${tmp_bam}" \
    "${sample_root}/${sample_id}.cram"
  samtools quickcheck -v "${tmp_bam}"
  samtools sort --threads "${SLURM_CPUS_PER_TASK}" -o "${final_bam}" "${tmp_bam}"
  samtools index --threads "${SLURM_CPUS_PER_TASK}" "${final_bam}"
  samtools quickcheck -v "${final_bam}"
  rm -f "${tmp_bam}"
fi

find "${sample_root}" -maxdepth 1 -type f -print0 | sort -z | xargs -0 sha256sum \
  > "${sample_root}/staged_files_sha256.txt"
touch "${sample_root}/STAGE_COMPLETE"
