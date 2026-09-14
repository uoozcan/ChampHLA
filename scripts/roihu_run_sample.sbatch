#!/bin/bash
#SBATCH --job-name=champhla_callers
#SBATCH --account=project_2008084
#SBATCH --partition=small
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=180G

# Run the frozen assay panel for one staged sample. Nextflow placeholders and
# masked process failures are rejected by the native-output checks below.
# nounset is enabled after the CSC environment is sourced: zz-csc-env.sh reads
# PS1, which is unbound in a batch shell, and would abort the job immediately.
set -eo pipefail

CODE_ROOT=${CHAMPHLA_CODE_ROOT:-/scratch/project_2008084/champhla_plurality}
source "${CHAMPHLA_PATH_CONFIG:-${CODE_ROOT}/configs/roihu_paths.env}"
export PYTHONPATH=${CHAMPHLA_CODE_ROOT}/src
manifest=${1:?usage: roihu_run_sample.sbatch MANIFEST LEDGER STAGE_LEDGER}
ledger=${2:?usage: roihu_run_sample.sbatch MANIFEST LEDGER STAGE_LEDGER}
stage_ledger=${3:?usage: roihu_run_sample.sbatch MANIFEST LEDGER STAGE_LEDGER}
run_id=${CHAMPHLA_RUN_ID:?CHAMPHLA_RUN_ID is required}
run_role=${CHAMPHLA_RUN_ROLE:?CHAMPHLA_RUN_ROLE is required}
[[ "${run_id}" =~ ^[a-z0-9][a-z0-9._-]{2,63}$ ]]
[[ "${run_role}" =~ ^(technical_pilot|capacity_validation|production)$ ]]
if [[ "${run_role}" != technical_pilot ]]; then
  python3 -c 'import json,sys; a=json.load(open(sys.argv[1])); assert a["status"] == "SIGNED_BY_AUTHOR" and a["signed_by"] and a["signed_at_utc"]' \
    "${CHAMPHLA_CODE_ROOT}/decisions/20260908_consensus_primary_amendment.json"
fi
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
input_root=${CHAMPHLA_INPUT_ROOT}/${modality}/${sample_id}
sample_root=${CHAMPHLA_RUN_ROOT}/caller_outputs/${run_id}/${run_role}/${modality}/${sample_id}
work_root=${CHAMPHLA_RUN_ROOT}/work/${run_id}/${run_role}/${modality}/${sample_id}
test -f "${ledger}"
test -f "${stage_ledger}"
test -f "${input_root}/STAGE_COMPLETE"
python3 - "${stage_ledger}" "${cohort}" "${sample_id}" "${modality}" <<'PY'
import csv, sys
ledger, cohort, sample, modality = sys.argv[1:]
rows = [r for r in csv.DictReader(open(ledger), delimiter="\t")
        if (r["cohort"], r["sample_id"], r["modality"]) == (cohort, sample, modality)]
assert rows, "missing staging ledger row"
latest = max(rows, key=lambda r: int(r["attempt"]))
assert latest["state"] == "validated", "staging ledger is not validated"
assert latest["extracted_sha256"], "staging ledger lacks extracted hash"
PY
test ! -e "${sample_root}/CALLERS_COMPLETE"
mkdir -p "${sample_root}/results" "${sample_root}/launch" "${work_root}"
test -s "${CHAMPHLA_RUN_ROOT}/environment_inventory.json"
python3 -c 'import json,sys; assert json.load(open(sys.argv[1]))["passed"] is True' \
  "${CHAMPHLA_RUN_ROOT}/environment_inventory.json"

export CSC_ENV_INIT_NON_INTERACTIVE=yes
export NXF_OPTS='-Xms256m -Xmx4g'
source /etc/profile.d/zz-csc-env.sh
set -u
module load bio-apps/v202603
module load samtools/1.21
module load nextflow/25.10.2-standalone

ledger_transition() {
  python3 -c 'from champhla_confirmation.cli import transition_run_ledger_main; raise SystemExit(transition_run_ledger_main())' \
    --ledger "${ledger}" --cohort "${cohort}" --sample "${sample_id}" \
    --modality "${modality}" "$@"
}
scheduler_task_id=${SLURM_ARRAY_JOB_ID:-${SLURM_JOB_ID}}_${SLURM_ARRAY_TASK_ID}
ledger_transition --state running --job-id "${scheduler_task_id}"
preserve_task_diagnostics() {
  # Nextflow keeps each task's stderr, stdout and script in its work directory, and the
  # parent log records only a truncated excerpt. Releasing work without copying these
  # first leaves a failure undiagnosable, so they are saved into the quarantine. Each is
  # capped, so an enormous log cannot defeat the point of reclaiming the space.
  local work_root=$1 dest=$2 task rel name
  [[ -d "${work_root}" ]] || return 0
  mkdir -p "${dest}"
  while IFS= read -r -d '' task; do
    rel=$(basename "$(dirname "${task}")")$(basename "${task}")
    for name in .command.sh .command.err .command.log .command.out .exitcode; do
      if [[ -f "${task}/${name}" ]]; then
        tail -c 1000000 "${task}/${name}" > "${dest}/${rel}${name}" 2>/dev/null || true
      fi
    done
  done < <(find "${work_root}" -mindepth 2 -maxdepth 2 -type d -print0 2>/dev/null)
}

reap_nextflow_children() {
  # Nextflow submits each process as its own Slurm job named nf-<PROCESS>_(<tag>), and
  # they survive the parent. Cancel only those carrying this sample's tag, so a
  # concurrently running sample is never affected.
  local orphans owner
  # USER can be unset in a batch shell, and this runs under `set -u` inside the trap.
  owner=${USER:-$(whoami)}
  orphans=$(squeue -h -u "${owner}" -o '%i %j' 2>/dev/null \
    | awk -v tag="(${sample_id}" '$2 ~ /^nf-/ && index($2, tag) {print $1}')
  if [[ -n "${orphans}" ]]; then
    echo "reaping nextflow child jobs for ${sample_id}: ${orphans}" >&2
    # shellcheck disable=SC2086
    scancel ${orphans} 2>/dev/null || true
  fi
}

on_failure() {
  local code=$?
  reap_nextflow_children

  # A ledger write can fail for reasons that have nothing to do with this job: a full
  # filesystem already left one row reading "running" while the job had failed. Keep going
  # so the quarantine still happens, but make the desync visible instead of inferred.
  local ledger_failed=
  ledger_transition --state failed --exit-code "${code}" || ledger_failed=1

  quarantine=${CHAMPHLA_RUN_ROOT}/quarantine/${run_id}/${run_role}/${modality}/${sample_id}.${SLURM_JOB_ID}
  mkdir -p "$(dirname "${quarantine}")"
  test ! -e "${quarantine}"
  if [[ -d "${sample_root}" ]]; then
    mv "${sample_root}" "${quarantine}"
  fi
  if [[ -n "${ledger_failed}" ]]; then
    echo "LEDGER WRITE FAILED for ${modality}/${sample_id}; ledger is out of sync" >&2
    mkdir -p "${quarantine}"
    printf 'state=failed exit_code=%s job_id=%s\n' "${code}" "${SLURM_JOB_ID}" \
      > "${quarantine}/LEDGER_WRITE_FAILED"
  fi

  # Record the work size before releasing it, so the storage measurement survives the
  # cleanup. Work is reproducible from the frozen inputs; leaving it behind on every
  # failure is what filled the allocation.
  if [[ -d "${work_root}" ]]; then
    mkdir -p "${quarantine}"
    du -sb "${work_root}" 2>/dev/null | awk '{print "failed_work_bytes=" $1}' \
      > "${quarantine}/WORK_SIZE_AT_FAILURE" || true
    preserve_task_diagnostics "${work_root}" "${quarantine}/work_diagnostics"
    if [[ "${CHAMPHLA_EXECUTION_MODE:-full_scale}" == sequential_low_storage ]]; then
      rm -rf "${work_root}"
    fi
  fi
  exit "${code}"
}
trap on_failure ERR

case "${modality}" in
  wgs)
    tools=spechla,hlahd,optitype,kourami,t1k
    seq_type=dna
    input_type_effective=bam
    printf 'sample_id,bam_path\n%s,%s\n' "${sample_id}" "${input_root}/${sample_id}.mate_aware_hla.bam" > "${sample_root}/samplesheet.csv"
    extra=(--spechla_exon_only 0 --extract_hla_region false)
    callers=(SpecHLA HLA-HD OptiType Kourami T1K)
    suffixes=(spechla hlahd optitype kourami t1k)
    ;;
  wes)
    tools=hlahd,optitype,polysolver,spechla,t1k
    seq_type=wes
    [[ "${input_type}" == bam || "${input_type}" == cram ]]
    input_type_effective=${input_type}
    printf 'sample_id,%s_path\n%s,%s\n' "${input_type}" "${sample_id}" \
      "${input_root}/${sample_id}.${input_type}" > "${sample_root}/samplesheet.csv"
    extra=(--spechla_exon_only 1 --extract_hla_region false)
    callers=(HLA-HD OptiType POLYSOLVER SpecHLA T1K)
    suffixes=(hlahd optitype polysolver spechla t1k)
    ;;
  rnaseq)
    tools=arcashla,hlahd,optitype,t1k
    seq_type=rna
    input_type_effective=fastq
    printf 'sample_id,fastq_1,fastq_2\n%s,%s,%s\n' "${sample_id}" "${input_root}/${sample_id}_R1.fastq.gz" "${input_root}/${sample_id}_R2.fastq.gz" > "${sample_root}/samplesheet.csv"
    extra=(--optitype_seq_type rna --extract_hla_region false)
    callers=(ArcasHLA HLA-HD OptiType T1K)
    suffixes=(arcashla hlahd optitype t1k)
    ;;
  *) echo "Unsupported modality: ${modality}" >&2; exit 2 ;;
esac

printf 'run_id=%s\nrun_role=%s\n' "${run_id}" "${run_role}" > "${sample_root}/run_identity.txt"
sha256sum "${sample_root}/samplesheet.csv" "${sample_root}/run_identity.txt" "${manifest}" \
  "${CHAMPHLA_RUN_ROOT}/environment_inventory.json" \
  "${CHAMPHLA_CODE_ROOT}/configs/comparator_manifest.json" \
  "${CHAMPHLA_CODE_ROOT}/configs/rerun_protocol.json" \
  "${CHAMPHLA_CODE_ROOT}/configs/roihu_workflow_lock.json" \
  "${CHAMPHLA_PIPELINE_ROOT}/main.nf" > "${sample_root}/pre_run_sha256.txt"
cd "${sample_root}/launch"
# GNU time is not installed on Roihu compute nodes and no cgroup peak-memory file is
# readable there, so elapsed time comes from the shell and peak RSS from Slurm accounting.
# Both sources are recorded so a reader can tell how each number was obtained.
run_started=${SECONDS}
nextflow run "${CHAMPHLA_PIPELINE_ROOT}/main.nf" \
  -name "plurality_${run_id}_${modality}_${sample_id}_${SLURM_JOB_ID}" \
  -w "${work_root}" \
  -params-file "${CHAMPHLA_PIPELINE_ROOT}/conf/roihu_params.yaml" \
  --input_samplesheet "${sample_root}/samplesheet.csv" --input_type "${input_type_effective}" \
  --tools "${tools}" --seq_type "${seq_type}" --run_modality "${modality}" \
  --outdir "${sample_root}/results" --slurm_account "${CHAMPHLA_ACCOUNT}" \
  --singularity_cache_dir "${CHAMPHLA_CONTAINER_ROOT}" "${extra[@]}" \
  -profile roihu,singularity
run_elapsed=$(( SECONDS - run_started ))

peak_rss_raw=$(sstat --noheader --parsable2 --format=MaxRSS -j "${SLURM_JOB_ID}.batch" 2>/dev/null | head -1 || true)
peak_rss_kb=$(printf '%s' "${peak_rss_raw}" | awk '
  {
    v=$0; u=substr(v, length(v), 1);
    if (u ~ /^[KMGT]$/) { n=substr(v, 1, length(v)-1)+0 } else { n=v+0; u="K" }
    m=(u=="K")?1:(u=="M")?1024:(u=="G")?1048576:(u=="T")?1073741824:1;
    printf "%d", n*m
  }')
[[ "${peak_rss_kb}" =~ ^[0-9]+$ ]] || peak_rss_kb=0
{
  printf 'elapsed_seconds=%s\n' "${run_elapsed}"
  printf 'peak_memory_kbytes=%s\n' "${peak_rss_kb}"
  printf 'elapsed_source=bash SECONDS\n'
  if [[ "${peak_rss_kb}" -gt 0 ]]; then
    printf 'peak_memory_source=slurm sstat MaxRSS\n'
  else
    printf 'peak_memory_source=unavailable\n'
  fi
} > "${sample_root}/resource_usage.txt"

validation=${sample_root}/caller_output_validation.tsv
printf 'caller\tnative_summary\tsha256\n' > "${validation}"
for index in "${!suffixes[@]}"; do
  suffix=${suffixes[$index]}
  caller=${callers[$index]}
  mapfile -t matches < <(find "${sample_root}/results" -type f -name "${sample_id}_${suffix}.txt" -size +0c -print)
  [[ ${#matches[@]} -eq 1 ]]
  ! grep -Eqi 'placeholder|stub|no results generated|warning:.*no output' "${matches[0]}"
  python3 - "${caller}" "${matches[0]}" <<'PY'
import sys
from champhla_confirmation.parsers import parse_caller_call
caller, path = sys.argv[1:]
native = open(path, encoding="utf-8", errors="strict").read()
for gene in ("A", "B", "C"):
    call = parse_caller_call(caller, native, gene)
    assert call["call_status"] in {"callable", "partial", "missing"}
PY
  printf '%s\t%s\t%s\n' "${suffix}" "${matches[0]}" "$(sha256sum "${matches[0]}" | cut -d ' ' -f 1)" >> "${validation}"
done
find "${sample_root}/results" -type f -print0 | sort -z | xargs -0 sha256sum > "${sample_root}/prediction_files_sha256.txt"
touch "${sample_root}/CALLERS_COMPLETE"
runtime=$(awk -F= '$1=="elapsed_seconds"{print $2}' "${sample_root}/resource_usage.txt")
peak_kb=$(awk -F= '$1=="peak_memory_kbytes"{print $2}' "${sample_root}/resource_usage.txt")
peak_memory=$((peak_kb * 1024))
# Peak is measured while everything is still on disk, so it reflects the true
# high-water mark for this sample.
peak_disk=$(du -sb "${input_root}" "${sample_root}" "${work_root}" | awk '{sum+=$1} END{print sum}')

# Release before measuring what is retained. The storage gate multiplies retained
# bytes by the full cohort size, so a figure taken before the release would
# describe a state that never persists and would make the projection meaningless.
if [[ "${CHAMPHLA_EXECUTION_MODE:-full_scale}" == sequential_low_storage ]]; then
  bash "${CHAMPHLA_CODE_ROOT}/scripts/roihu_prune_sample_work.sh" \
    "${CHAMPHLA_RUN_ROOT}" "${run_id}" "${run_role}" "${modality}" "${sample_id}"
fi
retained_disk=$(du -sb "${input_root}" "${sample_root}" | awk '{sum+=$1} END{print sum}')

output_hash=$(sha256sum "${sample_root}/prediction_files_sha256.txt" | cut -d ' ' -f 1)
ledger_transition --state succeeded --exit-code 0 --runtime-seconds "${runtime}" \
  --peak-memory-bytes "${peak_memory}" --peak-disk-bytes "${peak_disk}" \
  --retained-disk-bytes "${retained_disk}" \
  --output-sha256 "${output_hash}"
ledger_transition --state validated
trap - ERR
