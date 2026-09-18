#!/bin/bash
#SBATCH --job-name=pihla_new_figs
#SBATCH --account=project_2008084
#SBATCH --partition=small
#SBATCH --time=1:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --output=/scratch/project_2008084/pihla-publish/analysis/logs/new_figures_%j.out
#SBATCH --error=/scratch/project_2008084/pihla-publish/analysis/logs/new_figures_%j.err

set -euo pipefail

REPO=/scratch/project_2008084/pihla-publish
mkdir -p $REPO/analysis/logs

module load python-data

cd $REPO

echo "[$(date)] Generating new figures (A–E) and tables (A–B)"

python3 bin/generate_new_figures.py \
  --tables-dir $REPO/analysis/benchmark_trimodal_all_samples/tables \
  --figures-dir $REPO/analysis/figures_new

echo "[$(date)] Done"
ls $REPO/analysis/figures_new/
