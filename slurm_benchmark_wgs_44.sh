#!/bin/bash
#SBATCH --job-name=pihla_bench_wgs
#SBATCH --account=project_2008084
#SBATCH --partition=small
#SBATCH --time=4:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --output=/scratch/project_2008084/pihla-publish/analysis/logs/benchmark_wgs44_%j.out
#SBATCH --error=/scratch/project_2008084/pihla-publish/analysis/logs/benchmark_wgs44_%j.err

set -euo pipefail

mkdir -p /scratch/project_2008084/pihla-publish/analysis/logs

cd /scratch/project_2008084/pihla-publish

echo "[$(date)] Starting WGS benchmark (all included samples)"
echo "Config:     conf/benchmark_1000g_wgs.yaml"
echo "Output dir: /scratch/project_2008084/pihla-publish/analysis/benchmark_wgs_all_samples"

python3 bin/run_1000g_benchmark.py \
    --config conf/benchmark_1000g_wgs.yaml \
    --output-dir /scratch/project_2008084/pihla-publish/analysis/benchmark_wgs_all_samples

echo "[$(date)] WGS benchmark complete"
