#!/bin/bash
#SBATCH --job-name=pihla_cleanup_work
#SBATCH --account=project_2008084
#SBATCH --partition=small
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=4G
#SBATCH --output=/scratch/project_2008084/pihla-publish/analysis/logs/cleanup_work_%j.out
#SBATCH --error=/scratch/project_2008084/pihla-publish/analysis/logs/cleanup_work_%j.err

# Removes Nextflow pipeline work directories (~201G).
# Safe to run once pipeline outputs are in wgs_batches/, wes_batches/, rna_batches/.

set -euo pipefail
WORK=/scratch/project_2008084/hla_calibration/work

echo "[$(date)] Disk before:"
df -h /scratch/project_2008084/ | tail -1

echo "[$(date)] Removing Nextflow work dir: $WORK"
rm -rf "${WORK:?}"/*
echo "[$(date)] Done"

echo "[$(date)] Disk after:"
df -h /scratch/project_2008084/ | tail -1
