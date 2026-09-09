#!/bin/bash
# Retry exactly one quarantined failed sample without overwriting either attempt.
set -euo pipefail

CODE_ROOT=${CHAMPHLA_CODE_ROOT:-/scratch/project_2008084/champhla_plurality}
source "${CHAMPHLA_PATH_CONFIG:-${CODE_ROOT}/configs/roihu_paths.env}"
manifest=${1:?usage: roihu_retry_sample.sh MANIFEST MODALITY SAMPLE FAILED_JOB_ID}
modality=${2:?usage: roihu_retry_sample.sh MANIFEST MODALITY SAMPLE FAILED_JOB_ID}
sample_id=${3:?usage: roihu_retry_sample.sh MANIFEST MODALITY SAMPLE FAILED_JOB_ID}
failed_job=${4:?usage: roihu_retry_sample.sh MANIFEST MODALITY SAMPLE FAILED_JOB_ID}
export PYTHONPATH=${CHAMPHLA_CODE_ROOT}/src

wave=$(basename "${manifest}" .tsv)
ledger=${CHAMPHLA_RUN_ROOT}/manifests/${wave}.ledger.tsv
test -s "${ledger}"
test -f "${CHAMPHLA_INPUT_ROOT}/${modality}/${sample_id}/STAGE_COMPLETE"
test -d "${CHAMPHLA_RUN_ROOT}/quarantine/${modality}/${sample_id}.${failed_job}"
test ! -e "${CHAMPHLA_RUN_ROOT}/caller_outputs/${modality}/${sample_id}"
task_id=$(awk -F '\t' -v m="${modality}" -v s="${sample_id}" 'NR>1 && $4==m && $2==s {print NR-1}' "${manifest}")
[[ "${task_id}" =~ ^[0-9]+$ ]]
attempt=$(python3 -c 'import csv,sys; rows=[r for r in csv.DictReader(open(sys.argv[1]),delimiter="\t") if r["modality"]==sys.argv[2] and r["sample_id"]==sys.argv[3]]; assert rows and all(r["state"]=="failed" and r["job_id"]==sys.argv[4] for r in rows); print(max(int(r["attempt"]) for r in rows)+1)' \
  "${ledger}" "${modality}" "${sample_id}" "${failed_job}")
retry_job=$(sbatch --parsable --hold --array="${task_id}" \
  --export=ALL,CHAMPHLA_EXECUTION_MODE="${CHAMPHLA_EXECUTION_MODE:-full_scale}" \
  --output="${CHAMPHLA_RUN_ROOT}/logs/retry_${modality}_${sample_id}_%A_%a.out" \
  --error="${CHAMPHLA_RUN_ROOT}/logs/retry_${modality}_${sample_id}_%A_%a.err" \
  "${CHAMPHLA_CODE_ROOT}/scripts/roihu_run_sample.sbatch" "${manifest}")
python3 -c 'from champhla_confirmation.cli import transition_run_ledger_main; raise SystemExit(transition_run_ledger_main())' \
  --ledger "${ledger}" --cohort "$(awk -F '\t' -v n=$((task_id+1)) 'NR==n{print $1}' "${manifest}")" \
  --sample "${sample_id}" --modality "${modality}" --state resubmitted \
  --job-id "${retry_job}" --attempt "${attempt}" --supersedes-job-id "${failed_job}"
scontrol release "${retry_job}"
printf 'retry_job=%s attempt=%s supersedes=%s\n' "${retry_job}" "${attempt}" "${failed_job}"
