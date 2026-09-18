#!/bin/bash
#SBATCH --job-name=mvhla_sweep_wes
#SBATCH --account=project_2008084
#SBATCH --partition=small
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --output=/scratch/project_2008084/pihla-publish/analysis/logs/weight_sweep_wes_%j.out
#SBATCH --error=/scratch/project_2008084/pihla-publish/analysis/logs/weight_sweep_wes_%j.err

# Purpose: WES alpha sweep — 6 additional alpha values filling in the response surface.
# Existing weight_sensitivity already has alpha 0.0, 0.5, 1.0 (from trimodal_50samples).
# This run adds 0.2, 0.4, 0.6, 0.8 on the full WES all_samples benchmark (N=129).
# Critical finding to confirm: WeightedConsensus underperforms MajorityVote at ALL alpha
# values, proving the WES regression is a weight-learning failure, not an alpha-choice artifact.
# Addresses W4 with a full response surface. Feeds Figure A9 update.

set -euo pipefail

mkdir -p /scratch/project_2008084/pihla-publish/analysis/logs
mkdir -p /scratch/project_2008084/pihla-publish/analysis/weight_sensitivity_wes

cd /scratch/project_2008084/pihla-publish

BASEDIR=/scratch/project_2008084/pihla-publish/analysis/weight_sensitivity_wes

for ALPHA in 0.0 0.2 0.4 0.6 0.8 1.0; do
    BETA=$(python3 -c "print(round(1.0 - ${ALPHA}, 1))")
    LABEL="alpha_${ALPHA}_beta_${BETA}"
    OUTDIR="${BASEDIR}/${LABEL}"

    echo "[$(date)] Running WES sweep: alpha=${ALPHA} beta=${BETA}"
    echo "  Output: ${OUTDIR}"

    python3 bin/run_1000g_benchmark.py \
        --config conf/benchmark_1000g_wes.yaml \
        --output-dir "${OUTDIR}" \
        --weight-alpha "${ALPHA}" \
        --weight-beta  "${BETA}"

    echo "[$(date)] Done: ${LABEL}"
done

# Collect results into a single comparison table
echo "[$(date)] Collecting results..."
python3 - <<'PYEOF'
import csv, json
from pathlib import Path

BASEDIR = Path("/scratch/project_2008084/pihla-publish/analysis/weight_sensitivity_wes")
# Also include existing sensitivity points for comparison
EXTRA_DIRS = [
    ("0.0", "1.0", "confidence_only",  Path("/scratch/project_2008084/pihla-publish/analysis/weight_sensitivity/alpha_0.0_beta_1.0")),
    ("0.5", "0.5", "equal_weight",     Path("/scratch/project_2008084/pihla-publish/analysis/weight_sensitivity/alpha_0.5_beta_0.5")),
    ("1.0", "0.0", "reliability_only", Path("/scratch/project_2008084/pihla-publish/analysis/weight_sensitivity/alpha_1.0_beta_0.0")),
]

def read_tsv(p):
    with open(p) as f:
        return list(csv.DictReader(f, delimiter="\t"))

rows = []
for alpha_dir in sorted(BASEDIR.iterdir()):
    mc_path = alpha_dir / "tables" / "method_comparison.tsv"
    if not mc_path.exists():
        continue
    parts = alpha_dir.name.split("_")
    alpha = parts[1]
    beta  = parts[3] if len(parts) > 3 else str(round(1 - float(alpha), 1))
    for r in read_tsv(mc_path):
        mod_map = {"wes": "WES", "wgs": "WGS", "rnaseq": "RNA"}
        mod = mod_map.get(r.get("modality", ""), r.get("modality", ""))
        rows.append({
            "weight_alpha": alpha,
            "weight_beta":  beta,
            "label":        f"alpha_{alpha}",
            "modality":     mod,
            "method":       r.get("method", ""),
            "method_type":  r.get("method_type", ""),
            "overall_correct_call_rate": r.get("overall_correct_call_rate", ""),
            "callable_rate": r.get("callable_rate", ""),
            "ci_lo": r.get("overall_correct_call_rate_ci_lo", ""),
            "ci_hi": r.get("overall_correct_call_rate_ci_hi", ""),
            "source": str(mc_path),
        })

outpath = BASEDIR / "wes_weight_sweep_comparison.tsv"
fieldnames = ["weight_alpha","weight_beta","label","modality","method","method_type",
              "overall_correct_call_rate","callable_rate","ci_lo","ci_hi","source"]
with open(outpath, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
    w.writeheader()
    w.writerows(rows)
print(f"Wrote {len(rows)} rows to {outpath}")
PYEOF

echo "[$(date)] WES weight sweep complete."
echo "Results: ${BASEDIR}/wes_weight_sweep_comparison.tsv"
