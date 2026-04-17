#!/bin/bash
set -euo pipefail

INDEX=${INDEX:-/scratch/project_2008084/hla_calibration/rna/index/sample_fastq_urls.tsv}
BATCH_SIZE=${BATCH_SIZE:-4}
SCRIPT=${SCRIPT:-/scratch/project_2008084/pihla-publish/slurm_download_rna_batch.sh}

if [[ ! -f "${INDEX}" ]]; then
  echo "Missing index: ${INDEX}" >&2
  exit 1
fi
if [[ ! -x "${SCRIPT}" ]]; then
  echo "Missing batch script: ${SCRIPT}" >&2
  exit 1
fi

n=$(wc -l < "${INDEX}")
echo "Submitting RNA download batches for ${n} samples with batch size ${BATCH_SIZE}"

start=1
batch_num=1
while [[ ${start} -le ${n} ]]; do
  end=$((start + BATCH_SIZE - 1))
  if [[ ${end} -gt ${n} ]]; then
    end=${n}
  fi
  label=$(printf 'rna_batch_%02d_%03d_%03d' "${batch_num}" "${start}" "${end}")
  cmd=(sbatch --export=ALL,START_IDX=${start},END_IDX=${end},BATCH_LABEL=${label} "${SCRIPT}")
  echo "${cmd[*]}"
  "${cmd[@]}"
  start=$((end + 1))
  batch_num=$((batch_num + 1))
done
