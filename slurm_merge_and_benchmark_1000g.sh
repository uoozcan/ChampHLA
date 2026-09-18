#!/bin/bash
#SBATCH --job-name=pihla_merge_bench
#SBATCH --account=project_2008084
#SBATCH --partition=small
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --output=/scratch/project_2008084/hla_calibration/logs/merge_benchmark_%j.out
#SBATCH --error=/scratch/project_2008084/hla_calibration/logs/merge_benchmark_%j.err

set -euo pipefail

if [[ $# -lt 1 ]]; then
  echo "Usage: $0 <job_manifest.tsv> [sample_list] [analysis_root]" >&2
  exit 1
fi

job_manifest="$1"
sample_list="${2:-/scratch/project_2008084/pihla-publish/conf/1000g_smoke_samples.txt}"
analysis_root="${3:-/scratch/project_2008084/pihla-publish/analysis/1000g_realdata}"

wes_batch_root="/scratch/project_2008084/hla_calibration/wes_batches"
rna_batch_root="/scratch/project_2008084/hla_calibration/rna_batches"
wes_canonical="/scratch/project_2008084/hla_calibration/wes_3sample/results"
rna_canonical="/scratch/project_2008084/hla_calibration/rna_3sample/results"
status_tsv="${analysis_root}/batch_status.tsv"
phase_inputs="${analysis_root}/phase_gated_inputs"
benchmark_out="${analysis_root}/benchmark_phase_gated_abc"

mkdir -p /scratch/project_2008084/hla_calibration/logs
mkdir -p "${analysis_root}"

if [[ ! -f "${job_manifest}" ]]; then
  echo "Missing job manifest: ${job_manifest}" >&2
  exit 1
fi

mapfile -t samples < "${sample_list}"

readarray -t manifest_rows < <(tail -n +2 "${job_manifest}")
declare -A job_id_by_key
for row in "${manifest_rows[@]}"; do
  IFS=$'\t' read -r sample_id modality job_id batch_output_root _ <<< "${row}"
  job_id_by_key["${sample_id}|${modality}"]="${job_id}|${batch_output_root}"
done

get_expected_tools() {
  local modality="$1"
  if [[ "${modality}" == "wes" ]]; then
    printf '%s' "SpecHLA,HLA-HD,ArcasHLA,OptiType,POLYSOLVER,Kourami,T1K"
  else
    printf '%s' "ArcasHLA,OptiType,Seq2HLA,T1K,SpecHLA,HLA-HD"
  fi
}

collect_present_tools() {
  local results_dir="$1"
  python3 - "$results_dir" <<'PY'
import sys
from pathlib import Path

root = Path(sys.argv[1])
mapping = {
    "_spechla.txt": "SpecHLA",
    "_hlahd.txt": "HLA-HD",
    "_arcashla.txt": "ArcasHLA",
    "_optitype.txt": "OptiType",
    "_polysolver.txt": "POLYSOLVER",
    "_kourami.txt": "Kourami",
    "_t1k.txt": "T1K",
    "_seq2hla.txt": "Seq2HLA",
}
present = []
if root.exists():
    for suffix, label in mapping.items():
        if any(root.rglob(f"*{suffix}")):
            present.append(label)
print(",".join(sorted(set(present))))
PY
}

slurm_state() {
  local job_id="$1"
  local state
  state="$(sacct -X -n -P -j "${job_id}" -o State | head -n 1 | tr -d '[:space:]')"
  printf '%s' "${state:-UNKNOWN}"
}

rm -rf "${wes_canonical}" "${rna_canonical}"
mkdir -p "${wes_canonical}" "${rna_canonical}"

printf 'sample_id\tmodality\tjob_id\tslurm_state\tbatch_status\ttools_expected\ttools_present\tbatch_results_dir\n' > "${status_tsv}"

all_ok=1
for sample_id in "${samples[@]}"; do
  for modality in wes rnaseq; do
    key="${sample_id}|${modality}"
    if [[ -z "${job_id_by_key[${key}]:-}" ]]; then
      printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        "${sample_id}" "${modality}" "NA" "MISSING" "failed" \
        "$(get_expected_tools "${modality}")" "" "" >> "${status_tsv}"
      all_ok=0
      continue
    fi

    IFS='|' read -r job_id batch_output_root <<< "${job_id_by_key[${key}]}"
    batch_results_dir="${batch_output_root}/${sample_id}/results"
    state="$(slurm_state "${job_id}")"
    present_tools="$(collect_present_tools "${batch_results_dir}")"
    batch_status="failed"

    if [[ -d "${batch_results_dir}/${sample_id}" && -n "${present_tools}" ]]; then
      batch_status="success"
      target_root="${rna_canonical}"
      [[ "${modality}" == "wes" ]] && target_root="${wes_canonical}"
      mkdir -p "${target_root}"
      rm -rf "${target_root:?}/${sample_id}"
      cp -a "${batch_results_dir}/${sample_id}" "${target_root}/${sample_id}"
    else
      all_ok=0
    fi

    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
      "${sample_id}" "${modality}" "${job_id}" "${state}" "${batch_status}" \
      "$(get_expected_tools "${modality}")" "${present_tools}" "${batch_results_dir}" >> "${status_tsv}"
  done
done

if [[ "${all_ok}" -ne 1 ]]; then
  echo "At least one required per-sample batch failed; benchmark will not run." >&2
  exit 1
fi

cd /scratch/project_2008084/pihla-publish

python3 bin/build_1000g_phase_gated_inputs.py \
  --truth-csv /scratch/project_2008084/ozcanumu/hla_calibration/conf/ground_truth_data.csv \
  --wgs-results /scratch/project_2008084/hla_calibration/wgs/results \
  --wes-results "${wes_canonical}" \
  --rnaseq-results "${rna_canonical}" \
  --samples "$(paste -sd, "${sample_list}")" \
  --supported-loci A,B,C \
  --output-dir "${phase_inputs}"

python3 bin/run_1000g_benchmark.py \
  --config /scratch/project_2008084/pihla-publish/conf/benchmark_1000g_phase_gated_abc.yaml \
  --output-dir "${benchmark_out}"
