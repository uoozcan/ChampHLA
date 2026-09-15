#!/bin/bash
#SBATCH --job-name=champhla_finalize
#SBATCH --account=project_2008084
#SBATCH --partition=small
#SBATCH --time=00:10:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G

set -eo pipefail
CODE_ROOT=${CHAMPHLA_CODE_ROOT:-/scratch/project_2008084/champhla_plurality}
source "${CHAMPHLA_PATH_CONFIG:-${CODE_ROOT}/configs/roihu_paths.env}"
export PYTHONPATH=${CHAMPHLA_CODE_ROOT}/src
manifest=${1:?usage: finalize MANIFEST STAGE_LEDGER CALLER_LEDGER STAGE_JOB CALLER_JOB}
stage_ledger=${2:?}
caller_ledger=${3:?}
stage_job=${4:?}
caller_job=${5:?}
line_number=$((SLURM_ARRAY_TASK_ID + 1))
mapfile -t -d $'\t' _row < <(printf '%s' "$(sed -n "${line_number}p" "${manifest}")")
[[ ${#_row[@]} -eq 13 ]]
cohort=${_row[0]}; sample_id=${_row[1]}; modality=${_row[3]}

scheduler_state() {
  local job=$1 observed= poll
  # afterany may release before accounting is queryable. Bound the wait to one
  # minute and never leave an active ledger row merely because sacct lagged.
  for poll in $(seq 1 12); do
    observed=$(sacct -j "${job}_${SLURM_ARRAY_TASK_ID}" --noheader --parsable2 --format=State \
      | head -1 | cut -d'|' -f1 | cut -d+ -f1)
    if [[ -n "${observed}" ]]; then printf '%s\n' "${observed}"; return 0; fi
    sleep 5
  done
  printf 'terminal_unknown\n'
}
stage_state=$(scheduler_state "${stage_job}")
caller_state=$(scheduler_state "${caller_job}")
completion=${CHAMPHLA_RUN_ROOT}/caller_outputs/${CHAMPHLA_RUN_ID}/${CHAMPHLA_RUN_ROLE}/${modality}/${sample_id}/CALLERS_COMPLETE
python3 - "${stage_ledger}" "${caller_ledger}" "${cohort}" "${sample_id}" "${modality}" \
  "${stage_state}" "${caller_state}" "${completion}" <<'PY'
import sys
from champhla_confirmation.staging import finalize_sample
result = finalize_sample(*sys.argv[1:])
if not result["passed"]:
    raise SystemExit("sample not fully validated: " + "; ".join(result["failures"]))
PY
