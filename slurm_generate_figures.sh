#!/bin/bash
#SBATCH --job-name=pihla_figures
#SBATCH --account=project_2008084
#SBATCH --partition=small
#SBATCH --time=1:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --output=/scratch/project_2008084/pihla-publish/analysis/logs/figures_%j.out
#SBATCH --error=/scratch/project_2008084/pihla-publish/analysis/logs/figures_%j.err

# Regenerate publication-quality figures (PNG 300dpi + PDF) from the final
# tri-modal benchmark (benchmark_trimodal_50samples).
#
# Usage:
#   sbatch slurm_generate_figures.sh

set -euo pipefail

REPO=/scratch/project_2008084/pihla-publish
SCRIPT=$REPO/bin/generate_figures_v2.py

mkdir -p $REPO/analysis/logs

# Load Python environment (CSC Puhti)
module load python-data

cd $REPO

echo "[$(date)] Generating figures from tri-modal benchmark tables"

python3 $SCRIPT \
  --tables_dir $REPO/analysis/benchmark_trimodal_50samples/tables \
  --out_dir    $REPO/analysis/figures_final

echo "[$(date)] Figure generation complete"
ls $REPO/analysis/figures_final/
