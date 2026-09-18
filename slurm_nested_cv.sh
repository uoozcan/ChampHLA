#!/bin/bash
#SBATCH --job-name=champhla_nested_cv
#SBATCH --account=project_2008084
#SBATCH --partition=small
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --output=/scratch/project_2008084/pihla-publish/analysis/logs/nested_cv_%j.out
#SBATCH --error=/scratch/project_2008084/pihla-publish/analysis/logs/nested_cv_%j.err

# Official held-out benchmark entry point for ChampHLA.
#
# Reproduces the manuscript's primary (held-out, nested cross-validation) accuracy
# numbers and figures from the harmonised benchmark tables — no HLA tools are re-run.
# This is the authoritative path for the reported accuracy; the released Nextflow
# pipeline applies a single frozen operating point for production calls (see Methods §2).
#
# Usage:  sbatch slurm_nested_cv.sh      (or: bash slurm_nested_cv.sh on a login node)

set -euo pipefail
mkdir -p /scratch/project_2008084/pihla-publish/analysis/logs
cd /scratch/project_2008084/pihla-publish
module load python-data/3.12 2>/dev/null || true

echo "[$(date)] ChampHLA nested-CV re-analysis"
bash bin/reproduce_manuscript_numbers.sh

# Bimodal WES+RNA Champion-Challenger (recommended-configuration check)
echo "[$(date)] Bimodal WES+RNA nested-CV"
python3 bin/nested_cv_champion_challenger.py \
  --harmonized analysis/benchmark_trimodal_all_samples/tables/harmonized_benchmark_rows.tsv \
  --bimodal --genes A,B,C --folds 10 --seed 42 \
  --out analysis/nested_cv_champion_challenger/bimodal

echo "[$(date)] Done. Held-out method comparisons under analysis/nested_cv_champion_challenger/."
