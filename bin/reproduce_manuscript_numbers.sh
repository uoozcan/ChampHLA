#!/bin/bash
# Reproduce the ChampHLA manuscript's primary and stratified numbers + figures from the
# deposited harmonised benchmark tables (no HLA tools are re-run). See Data Availability.
#
# Regenerates, for WGS/WES/RNA:
#   - nested cross-validation held-out Champion-Challenger vs MajorityVote (Table 3, Fig 2/11)
#   - per-gene accuracy + per-gene McNemar (Table 4, Fig 3)
#   - CIWD-commonness and ancestry stratifications (Supplementary S9/S10)
#   - populates the champion_challenger_* stubs consumed by the figure generators
#   - the main + supplementary figures (into analysis/figures_repro/ to avoid clobbering figures_final)
#
# Usage:  bash bin/reproduce_manuscript_numbers.sh
set -euo pipefail

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO"
module load python-data/3.12 2>/dev/null || true   # Roihu; no-op elsewhere if python3 has pandas/matplotlib

ANALYSIS="$REPO/analysis"
NCV="$ANALYSIS/nested_cv_champion_challenger"
FIGREPRO="$ANALYSIS/figures_repro"
mkdir -p "$FIGREPRO"

echo "[1/4] Nested cross-validation (per modality) ..."
for pair in wgs:wgs wes:wes rna:rnaseq; do
  lbl="${pair%%:*}"; mod="${pair##*:}"
  D="$ANALYSIS/benchmark_${lbl}_cv_recalibrated/tables"
  python3 bin/nested_cv_champion_challenger.py \
    --harmonized "$D/harmonized_benchmark_rows.tsv" \
    --weights    "$D/consensus_runtime_weights.json" \
    --modality "$mod" --genes A,B,C --folds 10 --seed 42 \
    --out "$NCV/$lbl"
  # Weighting sensitivity control (Table 3, sensitivity + override-audit panels):
  # equal-weights and champion-routing-only should leave the significant WGS result
  # unchanged (0 overrides) and move WES/RNA only within their confidence intervals.
  python3 bin/nested_cv_champion_challenger.py \
    --harmonized "$D/harmonized_benchmark_rows.tsv" \
    --weights    "$D/consensus_runtime_weights.json" \
    --modality "$mod" --genes A,B,C --folds 10 --seed 42 \
    --equal-weights --out "$NCV/${lbl}_equalw"
  python3 bin/nested_cv_champion_challenger.py \
    --harmonized "$D/harmonized_benchmark_rows.tsv" \
    --weights    "$D/consensus_runtime_weights.json" \
    --modality "$mod" --genes A,B,C --folds 10 --seed 42 \
    --no-override --out "$NCV/${lbl}_nooverride"
done

echo "[2/4] Populate champion_challenger_* stubs from the nested-CV outputs ..."
python3 bin/populate_cc_stubs_from_nested_cv.py

echo "[3/4] Regenerate main + supplementary benchmark figures ..."
python3 figure_hub/scripts/fig_benchmark_main.py \
  --wgs-tables "$ANALYSIS/benchmark_wgs_cv_recalibrated/tables" \
  --wes-tables "$ANALYSIS/benchmark_wes_cv_recalibrated/tables" \
  --rna-tables "$ANALYSIS/benchmark_rna_cv_recalibrated/tables" \
  --out-dir "$FIGREPRO"

echo "[4/4] Regenerate the Champion-Challenger figure (Fig 11) ..."
python3 figure_hub/scripts/fig_champion_challenger.py --out-dir "$FIGREPRO"

echo
echo "Done. Held-out method comparison (should match Table 3):"
for lbl in wgs wes rna; do
  echo "  --- $lbl ---"
  column -t -s$'\t' "$NCV/$lbl/nested_cv_method_comparison.tsv" 2>/dev/null || cat "$NCV/$lbl/nested_cv_method_comparison.tsv"
done
echo
echo "Figures written to: $FIGREPRO  (promote to analysis/figures_final/ only after visual check)."
