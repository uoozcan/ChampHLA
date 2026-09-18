#!/bin/bash
#SBATCH --job-name=pihla_cleanup_rna
#SBATCH --account=project_2008084
#SBATCH --partition=small
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=4G
#SBATCH --output=/scratch/project_2008084/pihla-publish/analysis/logs/cleanup_rna_%j.out
#SBATCH --error=/scratch/project_2008084/pihla-publish/analysis/logs/cleanup_rna_%j.err

# Removes large intermediate files from already-processed RNA samples.
# Keeps all HLA result TXT files and confidence data.
# Safe to run once results have been verified in rna_batches/*/results/.
#
# Deletes per sample (~4G freed each):
#   - Extracted FASTQ copies used by ArcasHLA   (*.fq.gz, ~1.9G x2)
#   - ArcasHLA alignment index                   (*.alignment.p, ~266M)
#   - T1K candidate reads                        (*_candidate_*.fq, ~79M x2)
#   - T1K aligned reads                          (*_aligned_*.fa, ~9M x2)
#   - Original downloaded FASTQs                 (rna/fastqs/*.fastq.gz)

set -euo pipefail

RNA_BATCHES=/scratch/project_2008084/hla_calibration/rna_batches
RNA_FASTQS=/scratch/project_2008084/hla_calibration/rna/fastqs

mkdir -p /scratch/project_2008084/pihla-publish/analysis/logs

echo "[$(date)] Starting RNA intermediate cleanup"
echo "  Scanning: $RNA_BATCHES"

freed=0

for sample_dir in "$RNA_BATCHES"/*/; do
    sample=$(basename "$sample_dir")
    results_dir="$sample_dir/results/$sample"
    [[ -d "$results_dir" ]] || continue

    # Only clean if result files exist (safety check)
    n_results=$(find "$results_dir" -name "*.txt" -path "*/arcashla/*" 2>/dev/null | wc -l)
    if [[ "$n_results" -eq 0 ]]; then
        echo "  SKIP $sample: no arcashla result txt found"
        continue
    fi

    # ArcasHLA extracted FASTQs
    find "$results_dir/arcashla" -name "*.fq.gz" -delete 2>/dev/null || true
    echo "    $sample: removed arcashla FASTQs"

    # ArcasHLA alignment index
    find "$results_dir/arcashla" -name "*.alignment.p" -delete 2>/dev/null || true
    echo "    $sample: removed arcashla alignment"

    # T1K candidate and aligned reads (directory may not exist for all samples)
    [[ -d "$results_dir/t1k/t1k_out" ]] && \
        find "$results_dir/t1k/t1k_out" -name "*_candidate_*.fq" -delete 2>/dev/null || true
    [[ -d "$results_dir/t1k/t1k_out" ]] && \
        find "$results_dir/t1k/t1k_out" -name "*_aligned_*.fa"   -delete 2>/dev/null || true
    echo "    $sample: removed t1k intermediates (if present)"

    ((freed++)) || true
done

echo "[$(date)] Cleaned intermediates for $freed RNA samples"
echo ""

# Remove original downloaded FASTQs (results safely stored in rna_batches)
echo "[$(date)] Removing original RNA FASTQs from $RNA_FASTQS"
fastq_count=$(find "$RNA_FASTQS" -name "*.fastq.gz" 2>/dev/null | wc -l)
echo "  Found $fastq_count FASTQ files"
find "$RNA_FASTQS" -name "*.fastq.gz" -delete 2>/dev/null
echo "[$(date)] Original FASTQs removed"

echo ""
echo "[$(date)] Cleanup complete. Checking disk usage:"
df -h /scratch/project_2008084/
