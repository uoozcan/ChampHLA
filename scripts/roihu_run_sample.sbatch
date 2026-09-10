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
manifest=${1:?usage: roihu_run_sample.sbatch MANIFEST}
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
sample_root=${CHAMPHLA_RUN_ROOT}/caller_outputs/${modality}/${sample_id}
wave=$(basename "${manifest}" .tsv)
ledger=${CHAMPHLA_RUN_ROOT}/manifests/${wave}.ledger.tsv
test -f "${input_root}/STAGE_COMPLETE"
test ! -e "${sample_root}/CALLERS_COMPLETE"
mkdir -p "${sample_root}/results" "${sample_root}/launch" "${CHAMPHLA_RUN_ROOT}/work/${modality}/${sample_id}"
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
ledger_transition --state running --job-id "${SLURM_JOB_ID}"
on_failure() {
  local code=$?
  ledger_transition --state failed --exit-code "${code}" || true
  quarantine=${CHAMPHLA_RUN_ROOT}/quarantine/${modality}/${sample_id}.${SLURM_JOB_ID}
  mkdir -p "$(dirname "${quarantine}")"
  test ! -e "${quarantine}"
  if [[ -d "${sample_root}" ]]; then
    mv "${sample_root}" "${quarantine}"
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

sha256sum "${sample_root}/samplesheet.csv" "${manifest}" \
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
  -name "plurality_${modality}_${sample_id}_${SLURM_JOB_ID}" \
  -w "${CHAMPHLA_RUN_ROOT}/work/${modality}/${sample_id}" \
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
peak_disk=$(du -sb "${input_root}" "${sample_root}" "${CHAMPHLA_RUN_ROOT}/work/${modality}/${sample_id}" | awk '{sum+=$1} END{print sum}')

# Release before measuring what is retained. The storage gate multiplies retained
# bytes by the full cohort size, so a figure taken before the release would
# describe a state that never persists and would make the projection meaningless.
if [[ "${CHAMPHLA_EXECUTION_MODE:-full_scale}" == sequential_low_storage ]]; then
  bash "${CHAMPHLA_CODE_ROOT}/scripts/roihu_prune_sample_work.sh" \
    "${CHAMPHLA_RUN_ROOT}" "${modality}" "${sample_id}"
fi
retained_disk=$(du -sb "${input_root}" "${sample_root}" | awk '{sum+=$1} END{print sum}')

output_hash=$(sha256sum "${sample_root}/prediction_files_sha256.txt" | cut -d ' ' -f 1)
ledger_transition --state succeeded --exit-code 0 --runtime-seconds "${runtime}" \
  --peak-memory-bytes "${peak_memory}" --peak-disk-bytes "${peak_disk}" \
  --retained-disk-bytes "${retained_disk}" \
  --output-sha256 "${output_hash}"
ledger_transition --state validated
trap - ERR
