#!/bin/bash
#SBATCH --job-name=pihla_bench_trimodal
#SBATCH --account=project_2008084
#SBATCH --partition=small
#SBATCH --time=6:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --output=/scratch/project_2008084/pihla-publish/analysis/logs/benchmark_trimodal_%j.out
#SBATCH --error=/scratch/project_2008084/pihla-publish/analysis/logs/benchmark_trimodal_%j.err

set -euo pipefail

mkdir -p /scratch/project_2008084/pihla-publish/analysis/logs

cd /scratch/project_2008084/pihla-publish

# Tri-modal cohort: WGS+WES+RNA-seq, all included samples evaluated together.
# Weights are computed from and applied to the full cohort (no holdout split).

echo "[$(date)] Starting tri-modal benchmark (WGS+WES+RNA-seq)"
echo "Config:     conf/benchmark_1000g_full_cohort.yaml"
echo "Output dir: /scratch/project_2008084/pihla-publish/analysis/benchmark_trimodal_all_samples"

python3 bin/run_1000g_benchmark.py \
    --config conf/benchmark_1000g_full_cohort.yaml \
    --output-dir /scratch/project_2008084/pihla-publish/analysis/benchmark_trimodal_all_samples

echo "[$(date)] Tri-modal benchmark complete"
