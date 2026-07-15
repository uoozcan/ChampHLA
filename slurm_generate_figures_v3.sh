#!/bin/bash
#SBATCH --job-name=pihla_figs_v3
#SBATCH --account=project_2008084
#SBATCH --partition=small
#SBATCH --time=1:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --output=/scratch/project_2008084/pihla-publish/analysis/logs/figures_v3_%j.out
#SBATCH --error=/scratch/project_2008084/pihla-publish/analysis/logs/figures_v3_%j.err

set -euo pipefail
REPO=/scratch/project_2008084/pihla-publish
mkdir -p $REPO/analysis/logs
module load python-data
cd $REPO
echo "[$(date)] Generating majority-voting figures (v3)"
python3 bin/generate_figures_v3.py \
  --tables-dir $REPO/analysis/benchmark_trimodal_all_samples/tables \
  --out-dir    $REPO/analysis/figures_v3
echo "[$(date)] Done"
ls $REPO/analysis/figures_v3/
