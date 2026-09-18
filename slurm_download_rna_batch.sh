#!/bin/bash
#SBATCH --job-name=pihla_dl_rna_batch
#SBATCH --account=project_2008084
#SBATCH --partition=small
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --output=/scratch/project_2008084/hla_calibration/logs/dl_rna_batch_%j.out
#SBATCH --error=/scratch/project_2008084/hla_calibration/logs/dl_rna_batch_%j.err

set -euo pipefail

INDEX=${INDEX:-/scratch/project_2008084/hla_calibration/rna/index/sample_fastq_urls.tsv}
OUTDIR=${OUTDIR:-/scratch/project_2008084/hla_calibration/rna/fastqs}
STATUS_DIR=${STATUS_DIR:-/scratch/project_2008084/hla_calibration/rna/batch_status}
START_IDX=${START_IDX:?START_IDX is required}
END_IDX=${END_IDX:?END_IDX is required}
BATCH_LABEL=${BATCH_LABEL:-batch_${START_IDX}_${END_IDX}}

mkdir -p "${OUTDIR}"
mkdir -p "${STATUS_DIR}"
mkdir -p /scratch/project_2008084/hla_calibration/logs

status_file="${STATUS_DIR}/${BATCH_LABEL}_${SLURM_JOB_ID}.tsv"
printf 'sample\tindex\tstatus\tr1_bytes\tr2_bytes\tnote\n' > "${status_file}"

fetch_one() {
  local url="$1"
  local target="$2"
  local tmp="${target}.part"

  rm -f "${tmp}"
  wget -q --no-check-certificate -c -O "${tmp}" "${url}"
  if [[ ! -s "${tmp}" ]]; then
    rm -f "${tmp}"
    return 1
  fi
  mv -f "${tmp}" "${target}"
}

any_failed=0

echo "[$(date)] RNA batch download start: ${BATCH_LABEL} (indices ${START_IDX}-${END_IDX})"

for idx in $(seq "${START_IDX}" "${END_IDX}"); do
  line=$(sed -n "${idx}p" "${INDEX}" || true)
  if [[ -z "${line}" ]]; then
    printf 'NA\t%s\tskipped\t0\t0\tno_data_for_index\n' "${idx}" >> "${status_file}"
    continue
  fi

  sample_id=$(echo "${line}" | cut -f1)
  r1_raw=$(echo "${line}" | cut -f2 | cut -d'|' -f1)
  r2_raw=$(echo "${line}" | cut -f3 | cut -d'|' -f1)

  [[ "${r1_raw}" == http* || "${r1_raw}" == ftp* ]] && r1_url="${r1_raw}" || r1_url="http://${r1_raw}"
  [[ "${r2_raw}" == http* || "${r2_raw}" == ftp* ]] && r2_url="${r2_raw}" || r2_url="http://${r2_raw}"

  out_r1="${OUTDIR}/${sample_id}_R1.fastq.gz"
  out_r2="${OUTDIR}/${sample_id}_R2.fastq.gz"

  echo "[$(date)] ${sample_id} (index ${idx})"
  echo "  R1: ${r1_url}"
  echo "  R2: ${r2_url}"

  if [[ -s "${out_r1}" && -s "${out_r2}" ]]; then
    printf '%s\t%s\talready_complete\t%s\t%s\tcomplete_pair_present\n' \
      "${sample_id}" "${idx}" "$(stat -c%s "${out_r1}")" "$(stat -c%s "${out_r2}")" >> "${status_file}"
    echo "  complete pair already present, skipping"
    continue
  fi

  rm -f "${out_r1}" "${out_r2}" "${out_r1}.part" "${out_r2}.part"

  if ! fetch_one "${r1_url}" "${out_r1}"; then
    any_failed=1
    printf '%s\t%s\tfailed\t0\t0\tr1_download_failed\n' "${sample_id}" "${idx}" >> "${status_file}"
    echo "  R1 failed"
    continue
  fi

  if ! fetch_one "${r2_url}" "${out_r2}"; then
    any_failed=1
    rm -f "${out_r1}" "${out_r2}" "${out_r1}.part" "${out_r2}.part"
    printf '%s\t%s\tfailed\t0\t0\tr2_download_failed\n' "${sample_id}" "${idx}" >> "${status_file}"
    echo "  R2 failed"
    continue
  fi

  printf '%s\t%s\tcompleted\t%s\t%s\tdownload_ok\n' \
    "${sample_id}" "${idx}" "$(stat -c%s "${out_r1}")" "$(stat -c%s "${out_r2}")" >> "${status_file}"
  echo "  completed: R1=$(du -sh "${out_r1}" | cut -f1) R2=$(du -sh "${out_r2}" | cut -f1)"
done

echo "[$(date)] RNA batch download end: ${BATCH_LABEL}"
echo "Status file: ${status_file}"

exit ${any_failed}
