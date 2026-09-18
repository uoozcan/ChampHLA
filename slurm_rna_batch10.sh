#!/bin/bash
#SBATCH --job-name=pihla_rna_b10
#SBATCH --account=project_2008084
#SBATCH --partition=small
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --output=/scratch/project_2008084/pihla-publish/analysis/logs/rna_batch10_%j.out
#SBATCH --error=/scratch/project_2008084/pihla-publish/analysis/logs/rna_batch10_%j.err

# Downloads, types, and cleans 10 RNA samples.
# Usage: sbatch slurm_rna_batch10.sh <sample1> <sample2> ... <sample10>
#
# Each sample:
#   1. Downloads FASTQ pair from EBI
#   2. Submits HLA typing pipeline job
#   3. Waits for all 10 pipeline jobs to finish
#   4. Deletes FASTQs and large intermediates (keeps result TXT files)

set -euo pipefail

INDEX=/scratch/project_2008084/hla_calibration/rna/index/sample_fastq_urls_batch2.tsv
FASTQ_ROOT=/scratch/project_2008084/hla_calibration/rna/fastqs
BATCH_ROOT=/scratch/project_2008084/hla_calibration/rna_batches
WORK_ROOT=/scratch/project_2008084/hla_calibration/work
PIPELINE_DIR=/scratch/project_2008084/pihla-publish
SAMTOOLS=/projappl/project_2008084/bin/samtools

mkdir -p "$FASTQ_ROOT" /scratch/project_2008084/hla_calibration/logs

samples=("$@")
if [[ ${#samples[@]} -eq 0 ]]; then
    echo "Usage: sbatch slurm_rna_batch10.sh <sample1> [sample2 ...]" >&2
    exit 1
fi

echo "[$(date)] RNA batch: ${#samples[@]} samples: ${samples[*]}"
module load python-data 2>/dev/null || true

# ── Step 1: Download FASTQs ───────────────────────────────────────────────
echo "[$(date)] Step 1: downloading FASTQs"
for sample in "${samples[@]}"; do
    r1="$FASTQ_ROOT/${sample}_R1.fastq.gz"
    r2="$FASTQ_ROOT/${sample}_R2.fastq.gz"
    [[ -f "$r1" && -f "$r2" ]] && echo "  SKIP $sample (FASTQs exist)" && continue

    # Find URLs in index
    line=$(grep "^${sample}"$'\t' "$INDEX" 2>/dev/null | head -1 || true)
    [[ -z "$line" ]] && echo "  SKIP $sample: not in index" && continue

    url_r1=$(echo "$line" | cut -f2)
    url_r2=$(echo "$line" | cut -f3)

    echo "  Downloading $sample R1..."
    wget -q -O "$r1" "$url_r1" || curl -s -o "$r1" "$url_r1"
    echo "  Downloading $sample R2..."
    wget -q -O "$r2" "$url_r2" || curl -s -o "$r2" "$url_r2"
    echo "  Downloaded $sample"
done
echo "[$(date)] Downloads complete"

# ── Step 2: Submit pipeline jobs ──────────────────────────────────────────
echo "[$(date)] Step 2: submitting pipeline jobs"
job_id_list=""
for sample in "${samples[@]}"; do
    r1="$FASTQ_ROOT/${sample}_R1.fastq.gz"
    [[ ! -f "$r1" ]] && echo "  SKIP $sample: no FASTQ" && continue

    # Check if already has results
    if [[ -d "$BATCH_ROOT/$sample/results" ]]; then
        n_done=$(find "$BATCH_ROOT/$sample/results" -name "*_arcashla.txt" 2>/dev/null | wc -l)
    else
        n_done=0
    fi

    job_id=$(sbatch --parsable \
        "$PIPELINE_DIR/slurm_rna_sample.sh" \
        "$sample" "$FASTQ_ROOT" "$BATCH_ROOT" "$WORK_ROOT")
    job_id_list="${job_id_list} ${job_id}"
    echo "  Submitted $sample → job $job_id"
done

# ── Step 3: Wait for all pipeline jobs ────────────────────────────────────
if [[ -n "${job_id_list}" ]]; then
    echo "[$(date)] Step 3: waiting for jobs"
    # Poll until all jobs leave the queue
    while squeue --jobs="$(echo "$job_id_list" | tr ' ' ',' | sed 's/^,//')" -h 2>/dev/null | grep -q .; do
        sleep 60
    done
    echo "[$(date)] All pipeline jobs finished"
fi

# ── Step 4: Clean intermediates and FASTQs ────────────────────────────────
echo "[$(date)] Step 4: cleaning intermediates"
for sample in "${samples[@]}"; do
    results_dir="$BATCH_ROOT/$sample/results/$sample"
    [[ -d "$results_dir" ]] || continue

    n_results=$(find "$results_dir" -name "*.txt" -path "*/arcashla/*" 2>/dev/null | wc -l)
    if [[ "$n_results" -eq 0 ]]; then
        echo "  WARNING: $sample has no arcashla result — keeping files for inspection"
        continue
    fi

    find "$results_dir/arcashla" -name "*.fq.gz"          -delete 2>/dev/null || true
    find "$results_dir/arcashla" -name "*.alignment.p"    -delete 2>/dev/null || true
    [[ -d "$results_dir/t1k/t1k_out" ]] && \
        find "$results_dir/t1k/t1k_out" -name "*_candidate_*.fq" -delete 2>/dev/null || true
    [[ -d "$results_dir/t1k/t1k_out" ]] && \
        find "$results_dir/t1k/t1k_out" -name "*_aligned_*.fa"   -delete 2>/dev/null || true
    rm -f "$FASTQ_ROOT/${sample}_R1.fastq.gz" "$FASTQ_ROOT/${sample}_R2.fastq.gz"
    echo "  Cleaned $sample"
done

echo "[$(date)] Batch complete. Disk:"
df -h /scratch/project_2008084/ | tail -1
