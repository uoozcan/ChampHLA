#!/bin/bash
#SBATCH --job-name=pihla_bench_wgs_wave1
#SBATCH --account=project_2008084
#SBATCH --partition=small
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --output=/scratch/project_2008084/hla_calibration/logs/benchmark_wgs_wave1_%j.out
#SBATCH --error=/scratch/project_2008084/hla_calibration/logs/benchmark_wgs_wave1_%j.err

set -euo pipefail

module load nextflow >/dev/null 2>&1 || true
mkdir -p /scratch/project_2008084/hla_calibration/logs

cd /scratch/project_2008084/pihla-publish

python3 bin/run_1000g_benchmark.py \
  --config /scratch/project_2008084/pihla-publish/conf/benchmark_1000g_wgs_wave1.yaml \
  --output-dir /scratch/project_2008084/pihla-publish/analysis/1000g_realdata/benchmark_wgs_wave1
