#!/bin/bash
#SBATCH --job-name=mvhla_wes_cv
#SBATCH --account=project_2008084
#SBATCH --partition=small
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --output=/scratch/project_2008084/pihla-publish/analysis/logs/wes_cv_recalibrated_%j.out
#SBATCH --error=/scratch/project_2008084/pihla-publish/analysis/logs/wes_cv_recalibrated_%j.err

# Purpose: Re-run WES benchmark with probabilistic_recalibrated mode.
# Generates cross_validation_weight_summary.tsv with LOO-CV calibration metrics for WES.
# On WES, most tools pass calibration (applied status) — CV estimates will validate
# that the retained tools are stably calibrated, not just lucky on training split.
# Addresses W2 for WES modality.

set -euo pipefail

mkdir -p /scratch/project_2008084/pihla-publish/analysis/logs

cd /scratch/project_2008084/pihla-publish

OUTDIR=/scratch/project_2008084/pihla-publish/analysis/benchmark_wes_cv_recalibrated

echo "[$(date)] Starting WES CV-recalibrated benchmark"
echo "Config:     conf/benchmark_1000g_wes.yaml"
echo "Mode:       probabilistic_recalibrated (Platt scaling + LOO CV)"
echo "Output dir: ${OUTDIR}"

python3 bin/run_1000g_benchmark.py \
    --config conf/benchmark_1000g_wes.yaml \
    --output-dir "${OUTDIR}" \
    --benchmark-mode probabilistic_recalibrated

echo "[$(date)] Done. Key output: ${OUTDIR}/tables/cross_validation_weight_summary.tsv"
