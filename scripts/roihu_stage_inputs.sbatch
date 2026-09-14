#!/bin/bash
#SBATCH --job-name=champhla_stage
#SBATCH --account=project_2008084
#SBATCH --partition=small
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G

# Stage exactly one truth-free manifest row.
#   WGS: stream the source CRAM over the network and keep only the frozen
#        mate-aware HLA BAM. The 18 GB source is never written to disk.
#   WES/RNA: fetch and verify the source, which the callers consume directly.
#        It is released after validation by roihu_prune_sample_work.sh.
# nounset is enabled after the CSC environment is sourced: zz-csc-env.sh reads
# PS1, which is unbound in a batch shell, and would abort the job immediately.
set -eo pipefail

CODE_ROOT=${CHAMPHLA_CODE_ROOT:-/scratch/project_2008084/champhla_plurality}
source "${CHAMPHLA_PATH_CONFIG:-${CODE_ROOT}/configs/roihu_paths.env}"
manifest=${1:?usage: roihu_stage_inputs.sbatch MANIFEST STAGE_LEDGER}
stage_ledger=${2:?usage: roihu_stage_inputs.sbatch MANIFEST STAGE_LEDGER}
run_id=${CHAMPHLA_RUN_ID:?CHAMPHLA_RUN_ID is required}
run_role=${CHAMPHLA_RUN_ROLE:?CHAMPHLA_RUN_ROLE is required}
[[ "${run_id}" =~ ^[a-z0-9][a-z0-9._-]{2,63}$ ]]
[[ "${run_role}" =~ ^(technical_pilot|capacity_validation|production)$ ]]
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

if [[ -z "${sample_id:-}" || -z "${modality:-}" ]]; then
  echo "Missing manifest row for array task ${SLURM_ARRAY_TASK_ID}" >&2
  exit 2
fi
sample_root=${CHAMPHLA_INPUT_ROOT}/${modality}/${sample_id}
stage_root=${sample_root}.staging.${SLURM_JOB_ID}
quarantine=${CHAMPHLA_RUN_ROOT}/quarantine/staging/${run_id}/${run_role}/${modality}/${sample_id}.${SLURM_JOB_ID}
attempt=1
stage_started=${SECONDS}
scheduler_task_id=${SLURM_ARRAY_JOB_ID:-${SLURM_JOB_ID}}_${SLURM_ARRAY_TASK_ID}
export PYTHONPATH=${CHAMPHLA_CODE_ROOT}/src

stage_transition() {
  python3 -c 'from champhla_confirmation.cli import transition_stage_ledger_main; raise SystemExit(transition_stage_ledger_main())' \
    --ledger "${stage_ledger}" --cohort "${cohort}" --sample "${sample_id}" \
    --modality "${modality}" --attempt "${attempt}" "$@"
}
stage_transition --state running --scheduler-job-id "${scheduler_task_id}"

# Staging is idempotent. A batch that is re-run after a partial failure will find its
# earlier inputs still in place, and re-downloading a verified input wastes the transfer
# without changing anything. Accept an existing stage only if it is complete and every
# file still matches the checksum recorded when it was staged.
if [[ -e "${sample_root}" ]]; then
  # Strip any leading directory so both the relative format written now and the legacy
  # absolute-path format verify against the files where they actually are.
  if [[ "${modality}" != wgs && -f "${sample_root}/STAGE_COMPLETE" && -s "${sample_root}/staged_files_sha256.txt" ]] \
     && (cd "${sample_root}" && sed 's|  .*/|  |' staged_files_sha256.txt \
         | sha256sum --check --status); then
    echo "already staged and verified: ${modality}/${sample_id}"
    staged_hash=$(sha256sum "${sample_root}/staged_files_sha256.txt" | cut -d ' ' -f 1)
    stage_transition --state validated --exit-code 0 --extracted-sha256 "${staged_hash}" \
      --runtime-seconds "$((SECONDS-stage_started))" --peak-storage-bytes "$(du -sb "${sample_root}" | cut -f1)" \
      --source-checksum-verified true --index-checksum-verified true
    exit 0
  fi
  # Present but not verifiable: a broken stage, not a resumable one.
  echo "existing stage for ${modality}/${sample_id} failed verification; quarantining" >&2
  mkdir -p "$(dirname "${quarantine}")"
  mv "${sample_root}" "${quarantine}.unverified"
fi
test ! -e "${sample_root}"
test ! -e "${stage_root}"
test ! -e "${quarantine}"
mkdir -p "${stage_root}"
on_stage_failure() {
  local code=$?
  set +e
  mkdir -p "$(dirname "${quarantine}")"
  if [[ -d "${stage_root}" ]]; then
    mv "${stage_root}" "${quarantine}"
  fi
  stage_transition --state failed --exit-code "${code}" --failure-class deterministic \
    --runtime-seconds "$((SECONDS-stage_started))" --peak-storage-bytes 0
  exit "${code}"
}
trap on_stage_failure ERR

export CSC_ENV_INIT_NON_INTERACTIVE=yes
source /etc/profile.d/zz-csc-env.sh
set -u
module load bio-apps/v202603
module load samtools/1.21

fetch_one() {
  local uri=$1
  local target=$2
  if [[ "${uri}" == /* ]]; then
    ln -s "${uri}" "${target}"
  elif [[ "${uri}" == s3://* ]]; then
    aws s3 --no-sign-request cp "${uri}" "${target}"
  else
    curl --fail --location --retry 4 --output "${target}" "${uri}"
  fi
}

verify_one() {
  local expected=$1
  local target=$2
  expected=${expected#sha256:}
  if [[ ${#expected} -eq 64 ]]; then
    printf '%s  %s\n' "${expected}" "${target}" | sha256sum --check --status
  else
    expected=${expected#md5:}
    printf '%s  %s\n' "${expected}" "${target}" | md5sum --check --status
  fi
}

build_hla_regions() {
  local reference=$1
  local regions=$2
  awk -F '\t' 'BEGIN {OFS="\t"} $1 == "chr6" {print $1,27999999,34000000} $1 ~ /^HLA-[ABC]\*/ {print $1,0,$2}' \
    "${reference}.fai" > "${regions}"
  [[ $(wc -l < "${regions}") -ge 4 ]]
}

if [[ "${input_type}" == fastq_pair ]]; then
  IFS=';' read -r uri1 uri2 <<< "${input_uri}"
  IFS=';' read -r sha1 sha2 <<< "${source_checksum}"
  fetch_one "${uri1}" "${stage_root}/${sample_id}_R1.fastq.gz"
  fetch_one "${uri2}" "${stage_root}/${sample_id}_R2.fastq.gz"
  verify_one "${sha1}" "${stage_root}/${sample_id}_R1.fastq.gz"
  verify_one "${sha2}" "${stage_root}/${sample_id}_R2.fastq.gz"
elif [[ "${modality}" == wgs ]]; then
  # Stream the extraction: a 30x CRAM is roughly 18 GB and the callers never see
  # it, only the mate-aware HLA BAM derived from it. Retaining 137 of them does
  # not fit the project allocation, so the source is read over the network and
  # never written to disk. samtools 1.21 is built with libcurl=yes S3=yes GCS=yes
  # (plugins=no, so the handlers are compiled in and HTS_PATH is not needed).
  [[ "${input_type}" == cram ]]
  [[ "${input_uri}" == http://* || "${input_uri}" == https://* || "${input_uri}" == s3://* ]]
  reference=${CHAMPHLA_REFERENCE_ROOT}/genomes/GRCh38_full_analysis_set_plus_decoy_hla.fa
  regions=${stage_root}/grch38dh_hla_regions.bed
  build_hla_regions "${reference}" "${regions}"

  # The index is small, so it is fetched and verified rather than streamed. Passing
  # it explicitly also stops samtools caching a remote index into the working
  # directory, which previously cluttered $HOME.
  idx=${stage_root}/${sample_id}.cram.crai
  fetch_one "${index_uri}" "${idx}"
  verify_one "${index_checksum}" "${idx}"

  tmp_bam=${stage_root}/${sample_id}.mate_aware_hla.unsorted.bam
  final_bam=${stage_root}/${sample_id}.mate_aware_hla.bam
  max_attempts=${CHAMPHLA_STAGE_MAX_ATTEMPTS:-2}
  attempt_timeout=${CHAMPHLA_STAGE_ATTEMPT_TIMEOUT_SECONDS:-18000}
  retry_backoff=${CHAMPHLA_STAGE_RETRY_BACKOFF_SECONDS:-60}
  [[ "${max_attempts}" -eq 2 && "${attempt_timeout}" -eq 18000 ]]
  extraction_ok=false
  for attempt in $(seq 1 "${max_attempts}"); do
    attempt_root=${stage_root}/attempt_${attempt}
    test ! -e "${attempt_root}"
    mkdir -p "${attempt_root}"
    attempt_bam=${attempt_root}/${sample_id}.mate_aware_hla.unsorted.bam
    attempt_stderr=${attempt_root}/samtools.stderr
    attempt_started=${SECONDS}
    if [[ "${attempt}" -gt 1 ]]; then
      stage_transition --state submitted --append-attempt --scheduler-job-id "${scheduler_task_id}" \
        --superseded-attempt "$((attempt-1))"
      stage_transition --state running --scheduler-job-id "${scheduler_task_id}"
    fi
    set +e
    timeout --signal=TERM "${attempt_timeout}" samtools view \
      --threads "${SLURM_CPUS_PER_TASK}" --reference "${reference}" \
      --fetch-pairs --target-file "${regions}" --bam --output "${attempt_bam}" \
      -X "${input_uri}" "${idx}" 2>"${attempt_stderr}"
    extraction_code=$?
    set -e
    if [[ "${extraction_code}" -eq 0 ]]; then
      mv "${attempt_bam}" "${tmp_bam}"
      mv "${attempt_stderr}" "${stage_root}/extraction_attempt_${attempt}.stderr"
      rmdir "${attempt_root}"
      python3 - "${stage_root}/extraction_attempt_${attempt}.json" "${scheduler_task_id}" \
        "${attempt}" "${input_uri}" "${index_uri}" "$((SECONDS-attempt_started))" <<'PY'
import json, sys
path, job, attempt, source, index, runtime = sys.argv[1:]
with open(path, "w", encoding="utf-8", newline="\n") as handle:
    json.dump({"schema_version": "champhla-stage-attempt-1", "scheduler_job_id": job,
               "attempt": int(attempt), "state": "succeeded", "exit_code": 0,
               "source_uri": source, "index_uri": index,
               "runtime_seconds": int(runtime)}, handle, indent=2, sort_keys=True)
    handle.write("\n")
PY
      extraction_ok=true
      break
    fi
    classification=${attempt_root}/classification.json
    python3 -c 'from champhla_confirmation.cli import classify_stage_failure_main; raise SystemExit(classify_stage_failure_main())' \
      --exit-code "${extraction_code}" --stderr "${attempt_stderr}" \
      --source-uri "${input_uri}" --output "${classification}"
    failure_class=$(python3 -c 'import json,sys; print(json.load(open(sys.argv[1]))["failure_class"])' "${classification}")
    retryable=$(python3 -c 'import json,sys; print(str(json.load(open(sys.argv[1]))["retryable"]).lower())' "${classification}")
    partial_size=$(du -sb "${attempt_root}" | cut -f1)
    partial_bam_hash=""
    if [[ -f "${attempt_bam}" ]]; then partial_bam_hash=$(sha256sum "${attempt_bam}" | cut -d ' ' -f 1); fi
    python3 - "${attempt_root}/attempt_audit.json" "${scheduler_task_id}" "${attempt}" \
      "${input_uri}" "${index_uri}" "${extraction_code}" "${failure_class}" \
      "$((SECONDS-attempt_started))" "${partial_size}" "${partial_bam_hash}" <<'PY'
import json, sys
path, job, attempt, source, index, code, kind, runtime, size, partial_hash = sys.argv[1:]
with open(path, "w", encoding="utf-8", newline="\n") as handle:
    json.dump({"schema_version": "champhla-stage-attempt-1", "scheduler_job_id": job,
               "attempt": int(attempt), "state": "failed", "exit_code": int(code),
               "failure_class": kind, "source_uri": source, "index_uri": index,
               "runtime_seconds": int(runtime), "partial_bytes": int(size),
               "partial_bam_sha256": partial_hash, "cleanup_eligible": False},
              handle, indent=2, sort_keys=True)
    handle.write("\n")
PY
    partial_hash=$(find "${attempt_root}" -type f -print0 | sort -z | xargs -0 sha256sum | sha256sum | cut -d ' ' -f 1)
    attempt_quarantine=${quarantine}.attempt${attempt}
    test ! -e "${attempt_quarantine}"
    mv "${attempt_root}" "${attempt_quarantine}"
    stage_transition --state failed --exit-code "${extraction_code}" \
      --failure-class "${failure_class}" --extracted-sha256 "${partial_hash}" \
      --runtime-seconds "$((SECONDS-stage_started))" --peak-storage-bytes "${partial_size}"
    if [[ "${retryable}" != true || "${attempt}" -eq "${max_attempts}" ]]; then
      context_quarantine=${quarantine}.context
      test ! -e "${context_quarantine}"
      mv "${stage_root}" "${context_quarantine}"
      trap - ERR
      exit "${extraction_code}"
    fi
    sleep "${retry_backoff}"
  done
  [[ "${extraction_ok}" == true ]]
  samtools quickcheck -v "${tmp_bam}"
  samtools sort --threads "${SLURM_CPUS_PER_TASK}" -o "${final_bam}" "${tmp_bam}"
  samtools index --threads "${SLURM_CPUS_PER_TASK}" "${final_bam}"
  samtools quickcheck -v "${final_bam}"
  rm -f "${tmp_bam}"

  # The source checksum cannot be reverified because the CRAM is never on disk.
  # Record it as the declared upstream identity, and say so, rather than implying
  # a verification that did not happen.
  cat > "${stage_root}/streamed_source.json" <<JSON
{
  "sample_id": "${sample_id}",
  "modality": "${modality}",
  "source_uri": "${input_uri}",
  "source_checksum_declared": "${source_checksum}",
  "source_checksum_verified": false,
  "source_retained": false,
  "index_uri": "${index_uri}",
  "index_checksum_verified": true,
  "extracted_bam_sha256": "$(sha256sum "${final_bam}" | cut -d ' ' -f 1)",
  "extraction": "samtools view --fetch-pairs --target-file <chr6 MHC + HLA-A/B/C alt contigs>",
  "note": "The source CRAM was read over the network and never written to disk, so its checksum is the value declared in the frozen manifest and was not independently reverified. The extracted BAM is checksummed here and by staged_files_sha256.txt."
}
JSON
elif [[ "${input_type}" == cram || "${input_type}" == bam ]]; then
  # WES keeps the download path: its callers consume the CRAM itself, so swapping
  # in an extracted BAM would change what every WES caller sees. The staged input
  # is released after validation instead, by roihu_prune_sample_work.sh.
  raw=${stage_root}/${sample_id}.${input_type}
  idx=${raw}.$([[ "${input_type}" == cram ]] && printf 'crai' || printf 'bai')
  fetch_one "${input_uri}" "${raw}"
  fetch_one "${index_uri}" "${idx}"
  verify_one "${source_checksum}" "${raw}"
  verify_one "${index_checksum}" "${idx}"
  samtools quickcheck -v "${raw}"
else
  echo "Unsupported input type: ${input_type}" >&2
  exit 3
fi

# Relative names, so the record still verifies after the directory is renamed below.
( cd "${stage_root}" && find . -maxdepth 1 -type f ! -name staged_files_sha256.txt -printf '%P\0' \
  | sort -z | xargs -0 sha256sum > staged_files_sha256.txt )
test -s "${stage_root}/staged_files_sha256.txt"
touch "${stage_root}/STAGE_COMPLETE"
mv "${stage_root}" "${sample_root}"
staged_hash=$(sha256sum "${sample_root}/staged_files_sha256.txt" | cut -d ' ' -f 1)
if [[ "${modality}" == wgs ]]; then
  staged_hash=$(sha256sum "${sample_root}/${sample_id}.mate_aware_hla.bam" | cut -d ' ' -f 1)
  source_verified=false
else
  source_verified=true
fi
stage_transition --state validated --exit-code 0 --extracted-sha256 "${staged_hash}" \
  --runtime-seconds "$((SECONDS-stage_started))" --peak-storage-bytes "$(du -sb "${sample_root}" | cut -f1)" \
  --source-checksum-verified "${source_verified}" \
  --index-checksum-verified "$([[ -n "${index_uri}" ]] && printf true || printf not_applicable)"
trap - ERR
