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

# NOTE: This benchmark has already been run (2026-04-18, benchmark_trimodal_50samples).
# Re-submit only if the input data changes (e.g. additional samples typed).
# Tri-modal cohort: WGS n=99, WES n=51, RNA-seq n=50; 60/20/20 split.

echo "[$(date)] Starting tri-modal benchmark (WGS+WES+RNA-seq)"
echo "Config:     conf/benchmark_1000g_full_cohort.yaml"
echo "Output dir: /scratch/project_2008084/pihla-publish/analysis/benchmark_trimodal_50samples"

python3 bin/run_1000g_benchmark.py \
    --config conf/benchmark_1000g_full_cohort.yaml \
    --output-dir /scratch/project_2008084/pihla-publish/analysis/benchmark_trimodal_50samples

echo "[$(date)] Tri-modal benchmark complete"
