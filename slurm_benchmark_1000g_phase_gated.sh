#!/bin/bash
#SBATCH --job-name=pihla_bench_1000g
#SBATCH --account=project_2008084
#SBATCH --partition=small
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --output=/scratch/project_2008084/hla_calibration/logs/benchmark_1000g_%j.out
#SBATCH --error=/scratch/project_2008084/hla_calibration/logs/benchmark_1000g_%j.err

set -euo pipefail

module load nextflow >/dev/null 2>&1 || true
mkdir -p /scratch/project_2008084/hla_calibration/logs

cd /scratch/project_2008084/pihla-publish

python3 bin/build_1000g_phase_gated_inputs.py \
  --truth-csv /scratch/project_2008084/ozcanumu/hla_calibration/conf/ground_truth_data.csv \
  --wgs-results /scratch/project_2008084/hla_calibration/wgs/results \
  --wes-results /scratch/project_2008084/hla_calibration/wes_3sample/results \
  --rnaseq-results /scratch/project_2008084/hla_calibration/rna_3sample/results \
  --samples NA06985,NA06986,NA06994 \
  --supported-loci A,B,C \
  --output-dir /scratch/project_2008084/pihla-publish/analysis/1000g_realdata/phase_gated_inputs

python3 bin/run_1000g_benchmark.py \
  --config /scratch/project_2008084/pihla-publish/conf/benchmark_1000g_phase_gated_abc.yaml \
  --output-dir /scratch/project_2008084/pihla-publish/analysis/1000g_realdata/benchmark_phase_gated_abc
