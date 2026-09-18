#!/bin/bash

set -euo pipefail

sample_list="${1:-/scratch/project_2008084/pihla-publish/conf/1000g_smoke_samples.txt}"
analysis_root="${2:-/scratch/project_2008084/pihla-publish/analysis/1000g_realdata}"
wes_input_root="${3:-/scratch/project_2008084/hla_calibration/wes_3sample_input}"
rna_input_root="${4:-/scratch/project_2008084/hla_calibration/rna_3sample_input_named}"
wes_output_root="/scratch/project_2008084/hla_calibration/wes_batches"
rna_output_root="/scratch/project_2008084/hla_calibration/rna_batches"
work_root="/scratch/project_2008084/hla_calibration/work"
job_manifest="${analysis_root}/batch_job_manifest.tsv"

mkdir -p "${analysis_root}" /scratch/project_2008084/hla_calibration/logs "${wes_output_root}" "${rna_output_root}" "${work_root}"

printf 'sample_id\tmodality\tjob_id\tbatch_output_root\twork_root\n' > "${job_manifest}"

deps=()
while IFS= read -r sample_id; do
  [[ -n "${sample_id}" ]] || continue

  wes_submit="$(sbatch /scratch/project_2008084/pihla-publish/slurm_wes_sample.sh "${sample_id}" "${wes_input_root}" "${wes_output_root}" "${work_root}")"
  wes_job_id="${wes_submit##* }"
  printf '%s\t%s\t%s\t%s\t%s\n' "${sample_id}" "wes" "${wes_job_id}" "${wes_output_root}" "${work_root}" >> "${job_manifest}"
  deps+=("${wes_job_id}")

  rna_submit="$(sbatch /scratch/project_2008084/pihla-publish/slurm_rna_sample.sh "${sample_id}" "${rna_input_root}" "${rna_output_root}" "${work_root}")"
  rna_job_id="${rna_submit##* }"
  printf '%s\t%s\t%s\t%s\t%s\n' "${sample_id}" "rnaseq" "${rna_job_id}" "${rna_output_root}" "${work_root}" >> "${job_manifest}"
  deps+=("${rna_job_id}")
done < "${sample_list}"

dep_expr="$(IFS=:; printf '%s' "${deps[*]}")"
merge_submit="$(sbatch --dependency="afterany:${dep_expr}" /scratch/project_2008084/pihla-publish/slurm_merge_and_benchmark_1000g.sh "${job_manifest}" "${sample_list}" "${analysis_root}")"
merge_job_id="${merge_submit##* }"

echo "Job manifest: ${job_manifest}"
echo "Merge/benchmark job: ${merge_job_id}"
