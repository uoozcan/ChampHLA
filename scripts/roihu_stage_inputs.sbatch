#!/bin/bash
#SBATCH --job-name=champhla_stage
#SBATCH --account=project_2008084
#SBATCH --partition=small
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G

# Stage exactly one truth-free manifest row.
#   WGS: stream the source CRAM over the network and keep only the frozen
#        mate-aware HLA BAM. The 18 GB source is never written to disk.
#   WES/RNA: fetch and verify the source, which the callers consume directly.
#        It is released after validation by roihu_prune_sample_work.sh.
# nounset is enabled after the CSC environment is sourced: zz-csc-env.sh reads
# PS1, which is unbound in a batch shell, and would abort the job immediately.
set -eo pipefail

CODE_ROOT=${CHAMPHLA_CODE_ROOT:-/scratch/project_2008084/champhla_plurality}
source "${CHAMPHLA_PATH_CONFIG:-${CODE_ROOT}/configs/roihu_paths.env}"
manifest=${1:?usage: roihu_stage_inputs.sbatch MANIFEST}
line_number=$((SLURM_ARRAY_TASK_ID + 1))
# Tab is IFS whitespace, so `IFS=$'\t' read` collapses runs of tabs and silently drops
# empty columns, shifting every later field. RNA rows have no index, so index_uri and
# index_checksum are empty and source_checksum was receiving read_layout. mapfile -d does
# not collapse, so empty columns are preserved.
mapfile -t -d $'\t' _row < <(printf '%s' "$(sed -n "${line_number}p" "${manifest}")")
if [[ ${#_row[@]} -ne 13 ]]; then
  echo "Manifest row ${line_number} has ${#_row[@]} fields, expected 13" >&2
  exit 2
fi
cohort=${_row[0]}
sample_id=${_row[1]}
donor_id=${_row[2]}
modality=${_row[3]}
input_type=${_row[4]}
input_uri=${_row[5]}
index_uri=${_row[6]}
index_checksum=${_row[7]}
source_checksum=${_row[8]}
reference_build=${_row[9]}
read_layout=${_row[10]}
independence_stratum=${_row[11]}
evidence_role=${_row[12]}

if [[ -z "${sample_id:-}" || -z "${modality:-}" ]]; then
  echo "Missing manifest row for array task ${SLURM_ARRAY_TASK_ID}" >&2
  exit 2
fi
sample_root=${CHAMPHLA_INPUT_ROOT}/${modality}/${sample_id}
stage_root=${sample_root}.staging.${SLURM_JOB_ID}
quarantine=${CHAMPHLA_RUN_ROOT}/quarantine/staging/${modality}/${sample_id}.${SLURM_JOB_ID}
test ! -e "${sample_root}"
test ! -e "${stage_root}"
test ! -e "${quarantine}"
mkdir -p "${stage_root}"
on_stage_failure() {
  local code=$?
  mkdir -p "$(dirname "${quarantine}")"
  if [[ -d "${stage_root}" ]]; then
    mv "${stage_root}" "${quarantine}"
  fi
  exit "${code}"
}
trap on_stage_failure ERR

export CSC_ENV_INIT_NON_INTERACTIVE=yes
source /etc/profile.d/zz-csc-env.sh
set -u
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

build_hla_regions() {
  local reference=$1
  local regions=$2
  awk -F '\t' 'BEGIN {OFS="\t"} $1 == "chr6" {print $1,27999999,34000000} $1 ~ /^HLA-[ABC]\*/ {print $1,0,$2}' \
    "${reference}.fai" > "${regions}"
  [[ $(wc -l < "${regions}") -ge 4 ]]
}

if [[ "${input_type}" == fastq_pair ]]; then
  IFS=';' read -r uri1 uri2 <<< "${input_uri}"
  IFS=';' read -r sha1 sha2 <<< "${source_checksum}"
  fetch_one "${uri1}" "${stage_root}/${sample_id}_R1.fastq.gz"
  fetch_one "${uri2}" "${stage_root}/${sample_id}_R2.fastq.gz"
  verify_one "${sha1}" "${stage_root}/${sample_id}_R1.fastq.gz"
  verify_one "${sha2}" "${stage_root}/${sample_id}_R2.fastq.gz"
elif [[ "${modality}" == wgs ]]; then
  # Stream the extraction: a 30x CRAM is roughly 18 GB and the callers never see
  # it, only the mate-aware HLA BAM derived from it. Retaining 137 of them does
  # not fit the project allocation, so the source is read over the network and
  # never written to disk. samtools 1.21 is built with libcurl=yes S3=yes GCS=yes
  # (plugins=no, so the handlers are compiled in and HTS_PATH is not needed).
  [[ "${input_type}" == cram ]]
  [[ "${input_uri}" == http://* || "${input_uri}" == https://* || "${input_uri}" == s3://* ]]
  reference=${CHAMPHLA_REFERENCE_ROOT}/genomes/GRCh38_full_analysis_set_plus_decoy_hla.fa
  regions=${stage_root}/grch38dh_hla_regions.bed
  build_hla_regions "${reference}" "${regions}"

  # The index is small, so it is fetched and verified rather than streamed. Passing
  # it explicitly also stops samtools caching a remote index into the working
  # directory, which previously cluttered $HOME.
  idx=${stage_root}/${sample_id}.cram.crai
  fetch_one "${index_uri}" "${idx}"
  verify_one "${index_checksum}" "${idx}"

  tmp_bam=${stage_root}/${sample_id}.mate_aware_hla.unsorted.bam
  final_bam=${stage_root}/${sample_id}.mate_aware_hla.bam
  samtools view --threads "${SLURM_CPUS_PER_TASK}" --reference "${reference}" \
    --fetch-pairs --target-file "${regions}" --bam --output "${tmp_bam}" \
    -X "${input_uri}" "${idx}"
  samtools quickcheck -v "${tmp_bam}"
  samtools sort --threads "${SLURM_CPUS_PER_TASK}" -o "${final_bam}" "${tmp_bam}"
  samtools index --threads "${SLURM_CPUS_PER_TASK}" "${final_bam}"
  samtools quickcheck -v "${final_bam}"
  rm -f "${tmp_bam}"

  # The source checksum cannot be reverified because the CRAM is never on disk.
  # Record it as the declared upstream identity, and say so, rather than implying
  # a verification that did not happen.
  cat > "${stage_root}/streamed_source.json" <<JSON
{
  "sample_id": "${sample_id}",
  "modality": "${modality}",
  "source_uri": "${input_uri}",
  "source_checksum_declared": "${source_checksum}",
  "source_checksum_verified": false,
  "source_retained": false,
  "index_uri": "${index_uri}",
  "index_checksum_verified": true,
  "extracted_bam_sha256": "$(sha256sum "${final_bam}" | cut -d ' ' -f 1)",
  "extraction": "samtools view --fetch-pairs --target-file <chr6 MHC + HLA-A/B/C alt contigs>",
  "note": "The source CRAM was read over the network and never written to disk, so its checksum is the value declared in the frozen manifest and was not independently reverified. The extracted BAM is checksummed here and by staged_files_sha256.txt."
}
JSON
elif [[ "${input_type}" == cram || "${input_type}" == bam ]]; then
  # WES keeps the download path: its callers consume the CRAM itself, so swapping
  # in an extracted BAM would change what every WES caller sees. The staged input
  # is released after validation instead, by roihu_prune_sample_work.sh.
  raw=${stage_root}/${sample_id}.${input_type}
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

find "${stage_root}" -maxdepth 1 -type f ! -name staged_files_sha256.txt -print0 \
  | sort -z | xargs -0 sha256sum > "${stage_root}/staged_files_sha256.txt"
test -s "${stage_root}/staged_files_sha256.txt"
touch "${stage_root}/STAGE_COMPLETE"
mv "${stage_root}" "${sample_root}"
trap - ERR
