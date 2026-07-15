#!/bin/bash
# submit_trimodal_expansion.sh — Expand trimodal sample set by submitting
# WES and RNA-seq SLURM jobs for all currently-missing samples, then
# re-running the trimodal benchmark as a dependency.
#
# Usage: bash submit_trimodal_expansion.sh [--dry-run]
#
# What it does:
#   1. Runs find_trimodal_candidates.py to refresh gap manifests
#   2. Submits slurm_wes_batch10.sh for WES-missing samples
#   3. Submits slurm_rna_array.sh (array) for RNA-missing samples
#   4. Chains slurm_benchmark_trimodal.sh after both complete

set -euo pipefail

REPO=/scratch/project_2008084/pihla-publish
CALIB_ROOT=/scratch/project_2008084/hla_calibration
DRY_RUN=0
[[ "${1:-}" == "--dry-run" ]] && DRY_RUN=1 && echo "[dry-run] No jobs will be submitted"

RNA_BATCH5="$CALIB_ROOT/rna/index/sample_fastq_urls_batch5.tsv"
WES_MANIFEST="$REPO/analysis/wes_priority_manifest.txt"
BATCH_DIR="$REPO/analysis/rna_expansion_batches"

# ── Step 1: Refresh gap manifests ─────────────────────────────────────────
echo "[$(date)] Step 1: auditing current trimodal gaps"
if [[ "$DRY_RUN" -eq 1 ]]; then
    python3 "$REPO/bin/find_trimodal_candidates.py" --dry-run
else
    python3 "$REPO/bin/find_trimodal_candidates.py"
fi

# ── Step 2: WES jobs ──────────────────────────────────────────────────────
WES_JOB=""
if [[ -f "$WES_MANIFEST" ]] && [[ -s "$WES_MANIFEST" ]]; then
    mapfile -t wes_samples < "$WES_MANIFEST"
    echo ""
    echo "[$(date)] Step 2: submitting WES jobs for ${#wes_samples[@]} samples"
    for s in "${wes_samples[@]}"; do echo "  $s"; done

    if [[ "$DRY_RUN" -eq 0 ]]; then
        WES_JOB=$(sbatch --parsable \
            "$REPO/slurm_wes_batch10.sh" "${wes_samples[@]}")
        echo "  WES batch job: $WES_JOB"
    else
        echo "  [dry-run] would run: sbatch slurm_wes_batch10.sh ${wes_samples[*]}"
    fi
else
    echo ""
    echo "[$(date)] Step 2: no WES samples missing — skipping WES submission"
fi

# ── Step 3: RNA jobs ──────────────────────────────────────────────────────
ARRAY_JOB=""
if [[ -f "$RNA_BATCH5" ]] && [[ -s "$RNA_BATCH5" ]]; then
    # Extract sample IDs (column 1) from the batch5 URL manifest
    mapfile -t rna_samples < <(cut -f1 "$RNA_BATCH5")
    echo ""
    echo "[$(date)] Step 3: preparing RNA batch files for ${#rna_samples[@]} samples"

    if [[ "$DRY_RUN" -eq 0 ]]; then
        rm -rf "$BATCH_DIR"
        mkdir -p "$BATCH_DIR"
    fi

    BATCH_SIZE=10
    BATCH_IDX=0
    for (( i=0; i<${#rna_samples[@]}; i+=BATCH_SIZE )); do
        slice=("${rna_samples[@]:$i:$BATCH_SIZE}")
        batch_file="$BATCH_DIR/batch_${BATCH_IDX}.txt"
        echo "  batch_${BATCH_IDX}.txt: ${#slice[@]} samples (${slice[*]})"
        if [[ "$DRY_RUN" -eq 0 ]]; then
            printf '%s\n' "${slice[@]}" > "$batch_file"
        fi
        BATCH_IDX=$(( BATCH_IDX + 1 ))
    done
    LAST_IDX=$(( BATCH_IDX - 1 ))

    if [[ "$DRY_RUN" -eq 0 ]]; then
        ARRAY_JOB=$(sbatch --parsable \
            --array=0-${LAST_IDX}%1 \
            "$REPO/slurm_rna_array.sh" "$BATCH_DIR")
        echo "  RNA array job: $ARRAY_JOB"
    else
        echo "  [dry-run] would run: sbatch --array=0-${LAST_IDX}%1 slurm_rna_array.sh $BATCH_DIR"
    fi
else
    echo ""
    echo "[$(date)] Step 3: no RNA samples missing — skipping RNA submission"
fi

# ── Step 4: Chain benchmark ───────────────────────────────────────────────
echo ""
if [[ -z "$WES_JOB" && -z "$ARRAY_JOB" ]]; then
    echo "[$(date)] Nothing to do — already at trimodal ceiling."
    echo "  Run:  python3 bin/find_trimodal_candidates.py --dry-run  to confirm."
    exit 0
fi

# Build dependency string from whichever job IDs are non-empty
DEP_IDS=""
[[ -n "$WES_JOB"   ]] && DEP_IDS="${DEP_IDS}:${WES_JOB}"
[[ -n "$ARRAY_JOB" ]] && DEP_IDS="${DEP_IDS}:${ARRAY_JOB}"
DEP_IDS="${DEP_IDS#:}"   # strip leading colon

echo "[$(date)] Step 4: chaining benchmark re-run (dependency: afterok:${DEP_IDS})"
if [[ "$DRY_RUN" -eq 0 ]]; then
    BENCH=$(sbatch --parsable \
        --dependency=afterok:${DEP_IDS} \
        "$REPO/slurm_benchmark_trimodal.sh")
    echo "  Benchmark job: $BENCH"
else
    echo "  [dry-run] would run: sbatch --dependency=afterok:${DEP_IDS} slurm_benchmark_trimodal.sh"
fi

# ── Summary ───────────────────────────────────────────────────────────────
echo ""
echo "[$(date)] Submitted. Monitor with:"
echo "  squeue -u ozcanumu -o '%.10i %.22j %.8T %.10M %R'"
echo ""
echo "After completion, regenerate figures:"
echo "  cd $REPO"
echo "  python3 bin/generate_figures_v3.py \\"
echo "      --tables-dir analysis/benchmark_trimodal_all_samples/tables \\"
echo "      --out-dir /scratch/project_2008084/mvhla_figures_v5"
echo "  python3 bin/generate_html_report_v7.py \\"
echo "      --tables-dir analysis/benchmark_trimodal_all_samples/tables \\"
echo "      --figures-dir /scratch/project_2008084/mvhla_figures_v5 \\"
echo "      --venex-dir /scratch/project_2008084/hla_calibration/venex/results/wgs \\"
echo "      --output /scratch/project_2008084/mvhla_report_v9.html"
