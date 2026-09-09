#!/bin/bash
# Remove only a validated sample's reproducible Nextflow work in low-storage mode.
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
