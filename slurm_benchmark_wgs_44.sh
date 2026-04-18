#!/bin/bash
#SBATCH --job-name=pihla_bench_wgs44
#SBATCH --account=project_2008084
#SBATCH --partition=small
#SBATCH --time=4:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --output=/scratch/project_2008084/pihla-publish/analysis/logs/benchmark_wgs44_%j.out
#SBATCH --error=/scratch/project_2008084/pihla-publish/analysis/logs/benchmark_wgs44_%j.err

set -euo pipefail

mkdir -p /scratch/project_2008084/pihla-publish/analysis/logs

cd /users/ozcanumu/scratch/project_2008084/pihla-publish

echo "[$(date)] Starting WGS 44-sample (interim) benchmark"
echo "Config:     conf/benchmark_wgs_44only.yaml"
echo "Output dir: /scratch/project_2008084/pihla-publish/analysis/benchmark_wgs_44samples"

python3 bin/run_1000g_benchmark.py \
    --config conf/benchmark_wgs_44only.yaml \
    --output-dir /scratch/project_2008084/pihla-publish/analysis/benchmark_wgs_44samples

echo "[$(date)] WGS (44-sample interim) benchmark complete"
