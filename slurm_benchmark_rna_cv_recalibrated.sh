#!/bin/bash
#SBATCH --job-name=mvhla_rna_cv
#SBATCH --account=project_2008084
#SBATCH --partition=small
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --output=/scratch/project_2008084/pihla-publish/analysis/logs/rna_cv_recalibrated_%j.out
#SBATCH --error=/scratch/project_2008084/pihla-publish/analysis/logs/rna_cv_recalibrated_%j.err

# Purpose: Re-run RNA-seq benchmark with probabilistic_recalibrated mode.
# Generates cross_validation_weight_summary.tsv with LOO-CV calibration metrics for RNA.
# ArcasHLA (guardrail_status=poor_calibration on RNA) is the critical case —
# CV will confirm whether this is a stable finding or a training-split artifact.
# Addresses W2 for RNA modality.

set -euo pipefail

mkdir -p /scratch/project_2008084/pihla-publish/analysis/logs

cd /scratch/project_2008084/pihla-publish

OUTDIR=/scratch/project_2008084/pihla-publish/analysis/benchmark_rna_cv_recalibrated

echo "[$(date)] Starting RNA-seq CV-recalibrated benchmark"
echo "Config:     conf/benchmark_1000g_rna.yaml"
echo "Mode:       probabilistic_recalibrated (Platt scaling + LOO CV)"
echo "Output dir: ${OUTDIR}"

python3 bin/run_1000g_benchmark.py \
    --config conf/benchmark_1000g_rna.yaml \
    --output-dir "${OUTDIR}" \
    --benchmark-mode probabilistic_recalibrated

echo "[$(date)] Done. Key output: ${OUTDIR}/tables/cross_validation_weight_summary.tsv"
