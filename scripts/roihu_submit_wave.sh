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
python3 -c 'import json,subprocess,sys; data=json.load(open(sys.argv[1])); current=subprocess.check_output(["git","-C",sys.argv[2],"rev-parse","HEAD"],text=True).strip(); assert data["repository"]["clean"] is True and data["repository"]["commit"] == current' \
  "${CHAMPHLA_RUN_ROOT}/environment_inventory.json" "${CHAMPHLA_CODE_ROOT}"
python3 -c 'from champhla_confirmation.roihu import validate_workflow_lock; import sys; result=validate_workflow_lock(sys.argv[1],sys.argv[2]); assert result["passed"], result["failures"]' \
  "${CHAMPHLA_CODE_ROOT}/configs/roihu_workflow_lock.json" "${CHAMPHLA_CODE_ROOT}"
python3 -c 'from champhla_confirmation.manifests import validate_comparator_manifest; import sys; failures=validate_comparator_manifest(sys.argv[1], True); assert not failures, failures' \
  "${CHAMPHLA_CODE_ROOT}/configs/comparator_manifest.json"
if [[ "${manifest}" == *hprc* ]]; then
  python3 -c 'from champhla_confirmation.manifests import validate_hprc_truth_protocol; import sys; failures=validate_hprc_truth_protocol(sys.argv[1], True); assert not failures, failures' \
    "${CHAMPHLA_CODE_ROOT}/configs/hprc_truth_protocol.json"
fi
python3 -c 'from champhla_confirmation.cli import audit_run_manifest_main; raise SystemExit(audit_run_manifest_main())' \
  --manifest "${manifest}" --output "${audit}"

execution_mode=full_scale
case "${wave}" in
  pilot)
    # Pilots run in the same mode production will use, so the retained figure the
    # storage gate consumes describes a state that actually persists. Peak disk is
    # still measured before the release, so both projections stay honest.
    limit=2
    concurrency=1
    execution_mode=sequential_low_storage
    ;;
  capacity)
    limit=10
    test -s "${CHAMPHLA_RUN_ROOT}/storage_gate.json"
    execution_mode=$(python3 -c 'import json,sys; g=json.load(open(sys.argv[1])); assert g["passed"] is True; print(g["execution_mode"])' \
      "${CHAMPHLA_RUN_ROOT}/storage_gate.json")
    concurrency=$([[ "${execution_mode}" == full_scale ]] && printf '2' || printf '1')
    ;;
  production)
    python3 -c 'import json,sys; a=json.load(open(sys.argv[1])); assert a["status"] == "SIGNED_BY_AUTHOR" and a["signed_by"] and a["signed_at_utc"]' \
      "${CHAMPHLA_CODE_ROOT}/decisions/20260908_consensus_primary_amendment.json"
    test -s "${CHAMPHLA_RUN_ROOT}/storage_gate.json"
    python3 -c 'import json,sys; assert json.load(open(sys.argv[1]))["passed"] is True' \
      "${CHAMPHLA_RUN_ROOT}/storage_gate.json"
    python3 -c 'import hashlib,json,sys; gate=json.load(open(sys.argv[1])); observed=hashlib.sha256(open(sys.argv[2],"rb").read()).hexdigest(); assert gate["availability_source"] == "project_allocation" and gate["environment_inventory_sha256"] == observed' \
      "${CHAMPHLA_RUN_ROOT}/storage_gate.json" "${CHAMPHLA_RUN_ROOT}/environment_inventory.json"
    execution_mode=$(python3 -c 'import json,sys; g=json.load(open(sys.argv[1])); assert g["execution_mode"] in {"full_scale","sequential_low_storage"}; print(g["execution_mode"])' \
      "${CHAMPHLA_RUN_ROOT}/storage_gate.json")
    limit=1000000
    concurrency=$([[ "${execution_mode}" == full_scale ]] && printf '4' || printf '1')
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
ledger=${CHAMPHLA_RUN_ROOT}/manifests/${modality}_${wave}.ledger.tsv
test ! -e "${ledger}"
python3 -c 'from champhla_confirmation.cli import initialize_run_ledger_main; raise SystemExit(initialize_run_ledger_main())' \
  --manifest "${subset}" --output "${ledger}" --job-id "submission_pending"

if [[ "${execution_mode}" == sequential_low_storage ]]; then
  previous_job=
  caller_jobs=()
  for task_id in $(seq 1 "${records}"); do
    dependency=()
    if [[ -n "${previous_job}" ]]; then
      dependency=(--dependency="afterok:${previous_job}")
    fi
    stage_job=$(sbatch --parsable "${dependency[@]}" --array="${task_id}" \
      --output="${CHAMPHLA_RUN_ROOT}/logs/stage_${modality}_${wave}_%A_%a.out" \
      --error="${CHAMPHLA_RUN_ROOT}/logs/stage_${modality}_${wave}_%A_%a.err" \
      "${CHAMPHLA_CODE_ROOT}/scripts/roihu_stage_inputs.sbatch" "${subset}")
    caller_job=$(sbatch --parsable --dependency="afterok:${stage_job}" --array="${task_id}" \
      --export=ALL,CHAMPHLA_EXECUTION_MODE=sequential_low_storage \
      --output="${CHAMPHLA_RUN_ROOT}/logs/callers_${modality}_${wave}_%A_%a.out" \
      --error="${CHAMPHLA_RUN_ROOT}/logs/callers_${modality}_${wave}_%A_%a.err" \
      "${CHAMPHLA_CODE_ROOT}/scripts/roihu_run_sample.sbatch" "${subset}")
    caller_jobs+=("${caller_job}")
    previous_job=${caller_job}
  done
  printf 'execution_mode=%s caller_jobs=%s records=%s subset=%s ledger=%s\n' \
    "${execution_mode}" "${caller_jobs[*]}" "${records}" "${subset}" "${ledger}"
else
  stage_job=$(sbatch --parsable --array="1-${records}%${concurrency}" \
    --output="${CHAMPHLA_RUN_ROOT}/logs/stage_${modality}_${wave}_%A_%a.out" \
    --error="${CHAMPHLA_RUN_ROOT}/logs/stage_${modality}_${wave}_%A_%a.err" \
    "${CHAMPHLA_CODE_ROOT}/scripts/roihu_stage_inputs.sbatch" "${subset}")
  caller_job=$(sbatch --parsable --dependency="afterok:${stage_job}" --array="1-${records}%${concurrency}" \
    --export=ALL,CHAMPHLA_EXECUTION_MODE=full_scale \
    --output="${CHAMPHLA_RUN_ROOT}/logs/callers_${modality}_${wave}_%A_%a.out" \
    --error="${CHAMPHLA_RUN_ROOT}/logs/callers_${modality}_${wave}_%A_%a.err" \
    "${CHAMPHLA_CODE_ROOT}/scripts/roihu_run_sample.sbatch" "${subset}")
  printf 'execution_mode=%s stage_job=%s caller_job=%s records=%s subset=%s ledger=%s\n' \
    "${execution_mode}" "${stage_job}" "${caller_job}" "${records}" "${subset}" "${ledger}"
fi
