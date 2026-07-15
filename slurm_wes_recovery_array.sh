#!/bin/bash
#SBATCH --job-name=pihla_wes_recovery
#SBATCH --account=project_2008084
#SBATCH --partition=small
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --output=/scratch/project_2008084/hla_calibration/logs/wes_recovery_%A_%a.out
#SBATCH --error=/scratch/project_2008084/hla_calibration/logs/wes_recovery_%A_%a.err

set -euo pipefail

if [[ $# -lt 1 ]]; then
  echo "Usage: $0 <sample_manifest> [input_root] [output_root] [work_root]" >&2
  exit 1
fi

if [[ -z "${SLURM_ARRAY_TASK_ID:-}" ]]; then
  echo "SLURM_ARRAY_TASK_ID is required for array execution" >&2
  exit 1
fi

PIPELINE_DIR="${PIHLA_PIPELINE_DIR:-/scratch/project_2008084/pihla-publish}"
SAMPLE_SCRIPT="${PIPELINE_DIR}/slurm_wes_sample.sh"
sample_manifest="$1"
input_root="${2:-/scratch/project_2008084/hla_calibration/wes/bams}"
output_root="${3:-/scratch/project_2008084/hla_calibration/wes_batches}"
work_root="${4:-/scratch/project_2008084/hla_calibration/work}"

if [[ ! -f "${sample_manifest}" ]]; then
  echo "Missing sample manifest: ${sample_manifest}" >&2
  exit 1
fi

if [[ ! -x "${SAMPLE_SCRIPT}" ]]; then
  echo "Missing or non-executable sample launcher: ${SAMPLE_SCRIPT}" >&2
  exit 1
fi

sample_id="$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "${sample_manifest}" | tr -d '\r')"
if [[ -z "${sample_id}" ]]; then
  echo "No sample found for task index ${SLURM_ARRAY_TASK_ID} in ${sample_manifest}" >&2
  exit 1
fi

exec "${SAMPLE_SCRIPT}" "${sample_id}" "${input_root}" "${output_root}" "${work_root}"
