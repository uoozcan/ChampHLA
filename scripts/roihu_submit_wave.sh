#!/bin/bash
set -euo pipefail

CODE_ROOT=${CHAMPHLA_CODE_ROOT:-/scratch/project_2008084/champhla_plurality}
source "${CHAMPHLA_PATH_CONFIG:-${CODE_ROOT}/configs/roihu_paths.env}"
manifest=${1:?usage: roihu_submit_wave.sh MANIFEST MODALITY pilot|capacity|production}
modality=${2:?usage: roihu_submit_wave.sh MANIFEST MODALITY pilot|capacity|production}
wave=${3:?usage: roihu_submit_wave.sh MANIFEST MODALITY pilot|capacity|production}
export PYTHONPATH=${CHAMPHLA_CODE_ROOT}/src

audit=${CHAMPHLA_RUN_ROOT}/manifests/$(basename "${manifest}").audit.json
mkdir -p "${CHAMPHLA_RUN_ROOT}/manifests" "${CHAMPHLA_RUN_ROOT}/logs"
test -s "${CHAMPHLA_RUN_ROOT}/environment_inventory.json"
python3 -c 'import json,sys; assert json.load(open(sys.argv[1]))["passed"] is True' \
  "${CHAMPHLA_RUN_ROOT}/environment_inventory.json"
python3 -c 'from champhla_confirmation.cli import audit_run_manifest_main; raise SystemExit(audit_run_manifest_main())' \
  --manifest "${manifest}" --output "${audit}"

case "${wave}" in
  pilot) limit=2; concurrency=1 ;;
  capacity) limit=10; concurrency=2 ;;
  production)
    python3 -c 'import json,sys; a=json.load(open(sys.argv[1])); c=json.load(open(sys.argv[2])); assert a["status"] == "SIGNED_BY_AUTHOR" and a["signed_by"]; assert c["status"] == "FROZEN"' \
      "${CHAMPHLA_CODE_ROOT}/decisions/20260908_consensus_primary_amendment.json" \
      "${CHAMPHLA_CODE_ROOT}/configs/comparator_manifest.json"
    test -s "${CHAMPHLA_RUN_ROOT}/storage_gate.json"
    python3 -c 'import json,sys; assert json.load(open(sys.argv[1]))["passed"] is True' \
      "${CHAMPHLA_RUN_ROOT}/storage_gate.json"
    limit=1000000; concurrency=4
    ;;
  *) echo "Unknown wave: ${wave}" >&2; exit 2 ;;
esac

subset=${CHAMPHLA_RUN_ROOT}/manifests/${modality}_${wave}.tsv
awk -F '\t' -v modality="${modality}" -v limit="${limit}" 'BEGIN{OFS="\t"} NR==1{print;next} $4==modality && count<limit{print;count++}' "${manifest}" > "${subset}"
records=$(( $(wc -l < "${subset}") - 1 ))
[[ ${records} -gt 0 ]]
if [[ "${wave}" == production && "${manifest}" == *same_resource* ]]; then
  expected=$(python3 -c 'import json,sys; print(json.load(open(sys.argv[1]))[sys.argv[2]])' \
    "${CHAMPHLA_CODE_ROOT}/configs/roihu_storage_targets.json" "${modality}")
  [[ "${records}" -eq "${expected}" ]]
fi
stage_job=$(sbatch --parsable --array="1-${records}%${concurrency}" \
  --output="${CHAMPHLA_RUN_ROOT}/logs/stage_${modality}_${wave}_%A_%a.out" \
  --error="${CHAMPHLA_RUN_ROOT}/logs/stage_${modality}_${wave}_%A_%a.err" \
  "${CHAMPHLA_CODE_ROOT}/scripts/roihu_stage_inputs.sbatch" "${subset}")
caller_job=$(sbatch --parsable --dependency="afterok:${stage_job}" --array="1-${records}%${concurrency}" \
  --output="${CHAMPHLA_RUN_ROOT}/logs/callers_${modality}_${wave}_%A_%a.out" \
  --error="${CHAMPHLA_RUN_ROOT}/logs/callers_${modality}_${wave}_%A_%a.err" \
  "${CHAMPHLA_CODE_ROOT}/scripts/roihu_run_sample.sbatch" "${subset}")
ledger=${CHAMPHLA_RUN_ROOT}/manifests/${modality}_${wave}.ledger.tsv
python3 -c 'from champhla_confirmation.cli import initialize_run_ledger_main; raise SystemExit(initialize_run_ledger_main())' \
  --manifest "${subset}" --output "${ledger}" --job-id "${caller_job}"
printf 'stage_job=%s caller_job=%s records=%s subset=%s ledger=%s\n' \
  "${stage_job}" "${caller_job}" "${records}" "${subset}" "${ledger}"
