#!/bin/bash
# submit_wgs_wave3.sh — Submit WGS typing for 6 samples not yet in the CRAM index
# (NA11894, NA12234, NA12489, NA12763, NA12775, NA12812), then chain the
# benchmark re-run as a SLURM dependency.
#
# These samples have WES + RNA already done; adding WGS brings the trimodal
# cohort from 100 to 106 (the hard ceiling given current NYGC 30x availability).
#
# Usage: bash submit_wgs_wave3.sh [--dry-run]

set -euo pipefail

REPO=/scratch/project_2008084/pihla-publish
DRY_RUN=0
[[ "${1:-}" == "--dry-run" ]] && DRY_RUN=1 && echo "[dry-run] No jobs will be submitted"

SAMPLES=(NA11894 NA12234 NA12489 NA12763 NA12775 NA12812)

echo "[$(date)] WGS wave3: ${#SAMPLES[@]} samples — ${SAMPLES[*]}"

if [[ "$DRY_RUN" -eq 0 ]]; then
    WGS_JOB=$(sbatch --parsable "$REPO/slurm_wgs_batch11.sh" "${SAMPLES[@]}")
    echo "  WGS wave3 job: $WGS_JOB"

    BENCH=$(sbatch --parsable \
        --dependency=afterok:${WGS_JOB} \
        "$REPO/slurm_benchmark_trimodal.sh")
    echo "  Benchmark job: $BENCH"
else
    echo "  [dry-run] would run: sbatch slurm_wgs_batch11.sh ${SAMPLES[*]}"
    echo "  [dry-run] would chain: sbatch --dependency=afterok:<WGS_JOB> slurm_benchmark_trimodal.sh"
fi

echo ""
echo "Monitor with:  squeue -u ozcanumu -o '%.10i %.22j %.8T %.10M %R'"
echo ""
echo "After benchmark completes, regenerate figures:"
echo "  cd $REPO"
echo "  python3 bin/generate_figures_v3.py \\"
echo "      --tables-dir analysis/benchmark_trimodal_all_samples/tables \\"
echo "      --out-dir /scratch/project_2008084/mvhla_figures_v5"
echo "  python3 bin/generate_html_report_v7.py \\"
echo "      --tables-dir analysis/benchmark_trimodal_all_samples/tables \\"
echo "      --figures-dir /scratch/project_2008084/mvhla_figures_v5 \\"
echo "      --venex-dir /scratch/project_2008084/hla_calibration/venex/results/wgs \\"
echo "      --output /scratch/project_2008084/mvhla_report_v9.html"
