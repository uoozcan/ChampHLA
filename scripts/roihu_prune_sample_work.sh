#!/bin/bash
# Release a validated sample's reclaimable storage in low-storage mode.
#
# Two things are reclaimable once a sample is validated:
#   the Nextflow work directory, which is reproducible from the frozen inputs; and
#   for WES and RNA-seq, the staged source input, which the frozen manifest pins by
#   URI and checksum and which can therefore be re-fetched and re-verified.
#
# The WGS mate-aware BAM is never released. Under the streaming staging path it is
# the only surviving artifact of the extraction, and it is small.
#
# Native caller outputs, prediction evidence, logs and every checksum record are
# always retained.
set -euo pipefail

run_root=${1:?usage: roihu_prune_sample_work.sh RUN_ROOT MODALITY SAMPLE}
modality=${2:?usage: roihu_prune_sample_work.sh RUN_ROOT MODALITY SAMPLE}
sample_id=${3:?usage: roihu_prune_sample_work.sh RUN_ROOT MODALITY SAMPLE}
[[ "${run_root}" == /scratch/project_2008084/champhla_plurality_runs* ]]
[[ "${modality}" =~ ^(wgs|wes|rnaseq)$ ]]
[[ "${sample_id}" =~ ^[A-Za-z0-9._-]+$ ]]

sample_root=${run_root}/caller_outputs/${modality}/${sample_id}
work_root=${run_root}/work/${modality}/${sample_id}
test -f "${sample_root}/CALLERS_COMPLETE"
test -s "${sample_root}/caller_output_validation.tsv"
test -s "${sample_root}/prediction_files_sha256.txt"
expected=$([[ "${modality}" == rnaseq ]] && printf '4' || printf '5')
observed=$(( $(wc -l < "${sample_root}/caller_output_validation.tsv") - 1 ))
[[ "${observed}" -eq "${expected}" ]]
test -d "${work_root}"

find "${work_root}" -depth -type f -delete
find "${work_root}" -depth -type l -delete
find "${work_root}" -depth -type d -empty -delete
test ! -e "${work_root}"
printf 'work_pruned_after_validation\n' > "${sample_root}/LOW_STORAGE_WORK_PRUNED"

# Release the staged source input for the modalities whose callers consumed it
# directly. WGS is excluded: its staged artifact is the extracted BAM, not a source.
if [[ "${modality}" == wes || "${modality}" == rnaseq ]]; then
  input_root=${CHAMPHLA_INPUT_ROOT:?CHAMPHLA_INPUT_ROOT must be set}/${modality}/${sample_id}
  [[ "${input_root}" == /scratch/project_2008084/champhla_plurality_inputs* ]]
  test -d "${input_root}"
  # Refuse to release anything unless the record of what was staged survives.
  test -s "${input_root}/staged_files_sha256.txt"
  test -f "${input_root}/STAGE_COMPLETE"

  released=$(find "${input_root}" -maxdepth 1 -type f \
    ! -name staged_files_sha256.txt ! -name STAGE_COMPLETE \
    ! -name INPUT_RELEASED_AFTER_VALIDATION -printf '%s\n' | awk '{sum+=$1} END{print sum+0}')
  find "${input_root}" -maxdepth 1 -type f \
    ! -name staged_files_sha256.txt ! -name STAGE_COMPLETE \
    ! -name INPUT_RELEASED_AFTER_VALIDATION -delete
  find "${input_root}" -maxdepth 1 -type l -delete

  {
    printf 'released_bytes=%s\n' "${released}"
    printf 'released_after=CALLERS_COMPLETE_and_validation\n'
    printf 'recoverable_from=frozen manifest input_uri and source_checksum\n'
    printf 'retained_record=staged_files_sha256.txt\n'
  } > "${input_root}/INPUT_RELEASED_AFTER_VALIDATION"
  printf 'input_released_after_validation bytes=%s\n' "${released}" \
    > "${sample_root}/LOW_STORAGE_INPUT_RELEASED"
fi
