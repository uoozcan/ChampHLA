#!/bin/bash
set -euo pipefail

CODE_ROOT=${CHAMPHLA_CODE_ROOT:-/scratch/project_2008084/champhla_plurality}
source "${CHAMPHLA_PATH_CONFIG:-${CODE_ROOT}/configs/roihu_paths.env}"
manifest=${1:?usage: roihu_submit_wave.sh MANIFEST MODALITY WAVE RUN_ID [SIZE INDEX]}
modality=${2:?usage: roihu_submit_wave.sh MANIFEST MODALITY WAVE RUN_ID [SIZE INDEX]}
wave=${3:?usage: roihu_submit_wave.sh MANIFEST MODALITY WAVE RUN_ID [SIZE INDEX]}
run_id=${4:?usage: roihu_submit_wave.sh MANIFEST MODALITY WAVE RUN_ID [SIZE INDEX]}
batch_size=${5:-0}
batch_index=${6:-0}
[[ "${run_id}" =~ ^[a-z0-9][a-z0-9._-]{2,63}$ ]] || {
  echo "unsafe run_id: ${run_id}" >&2
  exit 2
}
[[ "${modality}" =~ ^(wgs|wes|rnaseq)$ ]] || { echo "unsupported modality" >&2; exit 2; }
case "${wave}" in
  pilot) run_role=technical_pilot ;;
  capacity) run_role=capacity_validation ;;
  production|batch) run_role=production ;;
  *) echo "Unknown wave: ${wave}" >&2; exit 2 ;;
esac
export PYTHONPATH=${CHAMPHLA_CODE_ROOT}/src

audit_label=${wave}
if [[ "${wave}" == batch ]]; then audit_label=batch${batch_index}; fi
audit=${CHAMPHLA_RUN_ROOT}/manifests/${run_id}_${modality}_${audit_label}.source_manifest.audit.json
mkdir -p "${CHAMPHLA_RUN_ROOT}/manifests" "${CHAMPHLA_RUN_ROOT}/logs"
test ! -e "${audit}"
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
chain_afterok=
# Only the two-sample technical pilots may run before the author signs. Capacity
# and every production form, including deterministic batches, share this gate.
if [[ "${run_role}" != technical_pilot ]]; then
  python3 -c 'import json,sys; a=json.load(open(sys.argv[1])); assert a["status"] == "SIGNED_BY_AUTHOR" and a["signed_by"] and a["signed_at_utc"]' \
    "${CHAMPHLA_CODE_ROOT}/decisions/20260908_consensus_primary_amendment.json"
fi
case "${wave}" in
  pilot)
    # Pilots run in the same mode production will use, so the retained figure the
    # storage gate consumes describes a state that actually persists. Peak disk is
    # still measured before the release, so both projections stay honest.
    limit=2
    concurrency=1
    execution_mode=sequential_low_storage
    ;;
  batch)
    # A batch is a deterministic slice of this modality's rows in frozen-manifest order.
    [[ "${batch_size}" =~ ^[1-9][0-9]*$ ]] || { echo "batch requires SIZE >= 1" >&2; exit 2; }
    [[ "${batch_index}" =~ ^[0-9]+$ ]] || { echo "batch requires INDEX >= 0" >&2; exit 2; }
    limit=${batch_size}
    concurrency=1
    execution_mode=sequential_low_storage
    # Only one sample is in flight at a time, so the gate is asked about a batch-sized
    # target rather than the whole cohort.
    test -s "${CHAMPHLA_RUN_ROOT}/storage_gate.json"
    python3 -c 'import json,sys; g=json.load(open(sys.argv[1])); assert g["passed"] is True, g["execution_mode"]' \
      "${CHAMPHLA_RUN_ROOT}/storage_gate.json"
    # A batch may not start while the previous one still has unvalidated rows.
    previous=${CHAMPHLA_RUN_ROOT}/manifests/${run_id}_${modality}_batch$((batch_index - 1)).ledger.tsv
    if [[ "${batch_index}" -gt 0 ]]; then
      test -s "${previous}"
      python3 -c 'import csv,sys; rows=list(csv.DictReader(open(sys.argv[1]),delimiter="\t")); bad=[r["sample_id"] for r in rows if r["state"] not in {"validated","frozen"}]; assert not bad, "previous batch has unvalidated samples: " + ",".join(sorted(set(bad)))' \
        "${previous}"
      previous_submission=${CHAMPHLA_RUN_ROOT}/manifests/${run_id}_${modality}_batch$((batch_index - 1)).submission.txt
      test -s "${previous_submission}"
      chain_afterok=$(awk -F= '$1=="last_finalizer_job"{print $2}' "${previous_submission}")
      [[ "${chain_afterok}" =~ ^[0-9]+$ ]]
    fi
    ;;
  capacity)
    limit=10
    test -s "${CHAMPHLA_RUN_ROOT}/storage_gate.json"
    execution_mode=$(python3 -c 'import json,sys; g=json.load(open(sys.argv[1])); assert g["passed"] is True; print(g["execution_mode"])' \
      "${CHAMPHLA_RUN_ROOT}/storage_gate.json")
    concurrency=$([[ "${execution_mode}" == full_scale ]] && printf '2' || printf '1')
    ;;
  production)
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
esac

label=${wave}
offset=0
if [[ "${wave}" == batch ]]; then
  label=batch${batch_index}
  offset=$((batch_index * batch_size))
fi
artifact_prefix=${run_id}_${modality}_${label}
subset=${CHAMPHLA_RUN_ROOT}/manifests/${artifact_prefix}.tsv
test ! -e "${subset}"
awk -F '\t' -v modality="${modality}" -v limit="${limit}" -v offset="${offset}" \
  'BEGIN{OFS="\t"} NR==1{print;next} $4==modality{ if (seen++ < offset) next; if (count<limit){print;count++} }' \
  "${manifest}" > "${subset}"
records=$(( $(wc -l < "${subset}") - 1 ))
[[ ${records} -gt 0 ]]
if [[ "${wave}" == production && "${manifest}" == *same_resource* ]]; then
  expected=$(python3 -c 'import json,sys; print(json.load(open(sys.argv[1]))[sys.argv[2]])' \
    "${CHAMPHLA_CODE_ROOT}/configs/roihu_storage_targets.json" "${modality}")
  [[ "${records}" -eq "${expected}" ]]
fi
ledger=${CHAMPHLA_RUN_ROOT}/manifests/${artifact_prefix}.ledger.tsv
stage_ledger=${CHAMPHLA_RUN_ROOT}/manifests/${artifact_prefix}.stage_ledger.tsv
submission_record=${CHAMPHLA_RUN_ROOT}/manifests/${artifact_prefix}.submission.txt
test ! -e "${ledger}"
test ! -e "${stage_ledger}"
test ! -e "${submission_record}"
python3 -c 'from champhla_confirmation.cli import initialize_run_ledger_main; raise SystemExit(initialize_run_ledger_main())' \
  --manifest "${subset}" --output "${ledger}" --run-id "${run_id}" \
  --run-role "${run_role}"
python3 -c 'from champhla_confirmation.cli import initialize_stage_ledger_main; raise SystemExit(initialize_stage_ledger_main())' \
  --manifest "${subset}" --output "${stage_ledger}" --run-id "${run_id}" \
  --run-role "${run_role}"

if [[ "${execution_mode}" == sequential_low_storage ]]; then
  previous_job=${chain_afterok}
  caller_jobs=()
  stage_jobs=()
  finalizer_jobs=()
  for task_id in $(seq 1 "${records}"); do
    dependency=()
    if [[ -n "${previous_job}" ]]; then
      dependency=(--dependency="afterok:${previous_job}")
    fi
    stage_job=$(sbatch --parsable --hold --kill-on-invalid-dep=yes "${dependency[@]}" --array="${task_id}" \
      --export=ALL,CHAMPHLA_RUN_ID="${run_id}",CHAMPHLA_RUN_ROLE="${run_role}" \
      --output="${CHAMPHLA_RUN_ROOT}/logs/stage_${artifact_prefix}_%A_%a.out" \
      --error="${CHAMPHLA_RUN_ROOT}/logs/stage_${artifact_prefix}_%A_%a.err" \
      "${CHAMPHLA_CODE_ROOT}/scripts/roihu_stage_inputs.sbatch" "${subset}" "${stage_ledger}")
    caller_job=$(sbatch --parsable --hold --kill-on-invalid-dep=yes --dependency="afterok:${stage_job}" --array="${task_id}" \
      --export=ALL,CHAMPHLA_EXECUTION_MODE=sequential_low_storage,CHAMPHLA_RUN_ID="${run_id}",CHAMPHLA_RUN_ROLE="${run_role}" \
      --output="${CHAMPHLA_RUN_ROOT}/logs/callers_${artifact_prefix}_%A_%a.out" \
      --error="${CHAMPHLA_RUN_ROOT}/logs/callers_${artifact_prefix}_%A_%a.err" \
      "${CHAMPHLA_CODE_ROOT}/scripts/roihu_run_sample.sbatch" "${subset}" "${ledger}" "${stage_ledger}")
    finalizer_job=$(sbatch --parsable --hold --dependency="afterany:${stage_job}:${caller_job}" --array="${task_id}" \
      --export=ALL,CHAMPHLA_RUN_ID="${run_id}",CHAMPHLA_RUN_ROLE="${run_role}" \
      --output="${CHAMPHLA_RUN_ROOT}/logs/finalize_${artifact_prefix}_%A_%a.out" \
      --error="${CHAMPHLA_RUN_ROOT}/logs/finalize_${artifact_prefix}_%A_%a.err" \
      "${CHAMPHLA_CODE_ROOT}/scripts/roihu_finalize_sample.sbatch" "${subset}" "${stage_ledger}" "${ledger}" \
      "${stage_job}" "${caller_job}")
    mapfile -t -d $'\t' submission_row < <(printf '%s' "$(sed -n "$((task_id+1))p" "${subset}")")
    cohort_id=${submission_row[0]}
    sample_id=${submission_row[1]}
    python3 -c 'from champhla_confirmation.cli import transition_stage_ledger_main; raise SystemExit(transition_stage_ledger_main())' \
      --ledger "${stage_ledger}" --cohort "${cohort_id}" --sample "${sample_id}" --modality "${modality}" \
      --attempt 1 --state submitted --scheduler-job-id "${stage_job}_${task_id}"
    python3 -c 'from champhla_confirmation.cli import transition_run_ledger_main; raise SystemExit(transition_run_ledger_main())' \
      --ledger "${ledger}" --cohort "${cohort_id}" --sample "${sample_id}" --modality "${modality}" \
      --state submitted --job-id "${caller_job}_${task_id}"
    scontrol release "${stage_job}" "${caller_job}" "${finalizer_job}"
    stage_jobs+=("${stage_job}")
    caller_jobs+=("${caller_job}")
    finalizer_jobs+=("${finalizer_job}")
    previous_job=${finalizer_job}
  done
  {
    printf 'run_id=%s\nrun_role=%s\nwave=%s\n' "${run_id}" "${run_role}" "${wave}"
    printf 'last_finalizer_job=%s\nstage_jobs=%s\ncaller_jobs=%s\nfinalizer_jobs=%s\n' "${previous_job}" "${stage_jobs[*]}" "${caller_jobs[*]}" "${finalizer_jobs[*]}"
    printf 'subset_sha256=%s\nledger=%s\nstage_ledger=%s\n' "$(sha256sum "${subset}" | cut -d ' ' -f 1)" "${ledger}" "${stage_ledger}"
  } > "${submission_record}"
  printf 'execution_mode=%s caller_jobs=%s records=%s subset=%s ledger=%s\n' \
    "${execution_mode}" "${caller_jobs[*]}" "${records}" "${subset}" "${ledger}"
else
  stage_job=$(sbatch --parsable --hold --array="1-${records}%${concurrency}" \
    --export=ALL,CHAMPHLA_RUN_ID="${run_id}",CHAMPHLA_RUN_ROLE="${run_role}" \
    --output="${CHAMPHLA_RUN_ROOT}/logs/stage_${artifact_prefix}_%A_%a.out" \
    --error="${CHAMPHLA_RUN_ROOT}/logs/stage_${artifact_prefix}_%A_%a.err" \
    "${CHAMPHLA_CODE_ROOT}/scripts/roihu_stage_inputs.sbatch" "${subset}" "${stage_ledger}")
  caller_job=$(sbatch --parsable --hold --kill-on-invalid-dep=yes --dependency="afterok:${stage_job}" --array="1-${records}%${concurrency}" \
    --export=ALL,CHAMPHLA_EXECUTION_MODE=full_scale,CHAMPHLA_RUN_ID="${run_id}",CHAMPHLA_RUN_ROLE="${run_role}" \
    --output="${CHAMPHLA_RUN_ROOT}/logs/callers_${artifact_prefix}_%A_%a.out" \
    --error="${CHAMPHLA_RUN_ROOT}/logs/callers_${artifact_prefix}_%A_%a.err" \
    "${CHAMPHLA_CODE_ROOT}/scripts/roihu_run_sample.sbatch" "${subset}" "${ledger}" "${stage_ledger}")
  finalizer_job=$(sbatch --parsable --hold --dependency="afterany:${stage_job}:${caller_job}" --array="1-${records}%${concurrency}" \
    --export=ALL,CHAMPHLA_RUN_ID="${run_id}",CHAMPHLA_RUN_ROLE="${run_role}" \
    --output="${CHAMPHLA_RUN_ROOT}/logs/finalize_${artifact_prefix}_%A_%a.out" \
    --error="${CHAMPHLA_RUN_ROOT}/logs/finalize_${artifact_prefix}_%A_%a.err" \
    "${CHAMPHLA_CODE_ROOT}/scripts/roihu_finalize_sample.sbatch" "${subset}" "${stage_ledger}" "${ledger}" \
    "${stage_job}" "${caller_job}")
  python3 - "${subset}" "${stage_ledger}" "${ledger}" "${modality}" "${stage_job}" "${caller_job}" <<'PY'
import csv, sys
from champhla_confirmation.staging import transition_stage_attempt
from champhla_confirmation.roihu import transition_run_sample
manifest, stages, callers, modality, stage_job, caller_job = sys.argv[1:]
for index, row in enumerate(csv.DictReader(open(manifest), delimiter="\t"), 1):
    transition_stage_attempt(stages, row["cohort"], row["sample_id"], modality, 1,
                             "submitted", scheduler_job_id=f"{stage_job}_{index}")
    transition_run_sample(callers, row["cohort"], row["sample_id"], modality,
                          "submitted", job_id=f"{caller_job}_{index}")
PY
  scontrol release "${stage_job}" "${caller_job}" "${finalizer_job}"
  {
    printf 'run_id=%s\nrun_role=%s\nwave=%s\n' "${run_id}" "${run_role}" "${wave}"
    printf 'stage_job=%s\nlast_caller_job=%s\nfinalizer_job=%s\n' "${stage_job}" "${caller_job}" "${finalizer_job}"
    printf 'subset_sha256=%s\nledger=%s\nstage_ledger=%s\n' "$(sha256sum "${subset}" | cut -d ' ' -f 1)" "${ledger}" "${stage_ledger}"
  } > "${submission_record}"
  printf 'execution_mode=%s stage_job=%s caller_job=%s records=%s subset=%s ledger=%s\n' \
    "${execution_mode}" "${stage_job}" "${caller_job}" "${records}" "${subset}" "${ledger}"
fi
