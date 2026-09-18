#!/bin/bash
#SBATCH --job-name=pihla_weight_sens
#SBATCH --account=project_2008084
#SBATCH --partition=small
#SBATCH --time=4:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --output=/scratch/project_2008084/pihla-publish/analysis/logs/weight_sensitivity_%j.out
#SBATCH --error=/scratch/project_2008084/pihla-publish/analysis/logs/weight_sensitivity_%j.err

# Weight formula sensitivity analysis for PIHLA benchmark.
# Runs 3 variants of the trimodal benchmark with different alpha/beta weights
# to demonstrate that the default 0.7/0.3 split is not cherry-picked.
#
# Variants:
#   alpha=1.0, beta=0.0  -> reliability-only (no confidence term)
#   alpha=0.5, beta=0.5  -> equal weighting
#   alpha=0.0, beta=1.0  -> confidence-only (no reliability term)
#
# Output: analysis/weight_sensitivity/{alpha_1.0_beta_0.0, alpha_0.5_beta_0.5, alpha_0.0_beta_1.0}/

set -euo pipefail

module load python-data

REPO=/scratch/project_2008084/pihla-publish
CONFIG=$REPO/conf/benchmark_1000g_full_cohort.yaml
BASE_OUT=$REPO/analysis/weight_sensitivity

mkdir -p $REPO/analysis/logs

cd $REPO

echo "[$(date)] Starting weight formula sensitivity analysis"

for ALPHA in 1.0 0.5 0.0; do
    BETA=$(python3 -c "print(round(1.0 - $ALPHA, 1))")
    OUT_DIR=$BASE_OUT/alpha_${ALPHA}_beta_${BETA}
    echo "[$(date)] Running alpha=$ALPHA beta=$BETA -> $OUT_DIR"
    python3 bin/hla_benchmark.py \
        --config $CONFIG \
        --output-dir $OUT_DIR \
        --weight-alpha $ALPHA \
        --weight-beta $BETA
    echo "[$(date)] Done: alpha=$ALPHA beta=$BETA"
done

echo "[$(date)] All 3 weight variants complete"
echo "Results:"
for DIR in $BASE_OUT/alpha_*; do
    echo "  $DIR"
    cut -f1-3,6-9 $DIR/tables/method_comparison.tsv 2>/dev/null | head -5
    echo "---"
done
