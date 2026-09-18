#!/bin/bash
#SBATCH --job-name=mvhla_wgs_cv
#SBATCH --account=project_2008084
#SBATCH --partition=small
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --output=/scratch/project_2008084/pihla-publish/analysis/logs/wgs_cv_recalibrated_%j.out
#SBATCH --error=/scratch/project_2008084/pihla-publish/analysis/logs/wgs_cv_recalibrated_%j.err

# Purpose: Re-run WGS benchmark with probabilistic_recalibrated mode.
# Generates cross_validation_weight_summary.tsv with LOO-CV ECE/Brier metrics.
# Addresses W2: replaces single-split calibration estimates with cross-validated ones.
# Output compared against benchmark_wgs_all_samples/ to verify guardrail decision stability.

set -euo pipefail

mkdir -p /scratch/project_2008084/pihla-publish/analysis/logs

cd /scratch/project_2008084/pihla-publish

OUTDIR=/scratch/project_2008084/pihla-publish/analysis/benchmark_wgs_cv_recalibrated

echo "[$(date)] Starting WGS CV-recalibrated benchmark"
echo "Config:     conf/benchmark_1000g_wgs.yaml"
echo "Mode:       probabilistic_recalibrated (Platt scaling + LOO CV)"
echo "Output dir: ${OUTDIR}"

python3 bin/run_1000g_benchmark.py \
    --config conf/benchmark_1000g_wgs.yaml \
    --output-dir "${OUTDIR}" \
    --benchmark-mode probabilistic_recalibrated

echo "[$(date)] Done. Key output: ${OUTDIR}/tables/cross_validation_weight_summary.tsv"
echo ""
echo "To compare CV vs legacy calibration metrics:"
echo "  diff ${OUTDIR}/tables/tool_confidence_weights.tsv \\"
echo "       analysis/benchmark_wgs_all_samples/tables/tool_confidence_weights.tsv"
