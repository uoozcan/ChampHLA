#!/bin/bash
#SBATCH --job-name=pihla_figures
#SBATCH --account=project_2008084
#SBATCH --partition=small
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --output=logs/slurm_figures_%j.out
#SBATCH --error=logs/slurm_figures_%j.err

# Regenerate publication-quality figures (PNG + PDF) for key benchmark runs.
# Requires generate_figures_v2.py and a complete Python env with matplotlib + seaborn.
#
# Usage:
#   sbatch slurm_generate_figures.sh
#
# Outputs: <benchmark_dir>/figures_v2/{figure_01.png, figure_01.pdf, ...}

set -euo pipefail

REPO=/scratch/project_2008084/pihla-publish
SCRIPT=$REPO/bin/generate_figures_v2.py

mkdir -p $REPO/logs

# Load Python environment (adjust module if needed)
module load python-data/3.9

echo "=== Generating figures for benchmark_wgs_wave1 ==="
python3 $SCRIPT \
  --tables_dir $REPO/analysis/1000g_realdata/benchmark_wgs_wave1/tables \
  --out_dir    $REPO/analysis/1000g_realdata/benchmark_wgs_wave1/figures_v2

echo "=== Generating figures for full_cohort_benchmark ==="
python3 $SCRIPT \
  --tables_dir $REPO/analysis/full_cohort_benchmark/tables \
  --out_dir    $REPO/analysis/full_cohort_benchmark/figures_v2

echo "=== Generating figures for benchmark_wes_50samples ==="
python3 $SCRIPT \
  --tables_dir $REPO/analysis/benchmark_wes_50samples/tables \
  --out_dir    $REPO/analysis/benchmark_wes_50samples/figures_v2

echo "=== Generating figures for benchmark_rna_50samples ==="
python3 $SCRIPT \
  --tables_dir $REPO/analysis/benchmark_rna_50samples/tables \
  --out_dir    $REPO/analysis/benchmark_rna_50samples/figures_v2

echo "=== All figure runs complete ==="
