#!/bin/bash
#SBATCH --job-name=nci60_wes_full
#SBATCH --account=project_2008084
#SBATCH --partition=small
#SBATCH --time=36:00:00
#SBATCH --cpus-per-task=20
#SBATCH --mem=180G
#SBATCH --output=/scratch/project_2008084/hla_calibration/logs/nci60_wes_full_%j.out
#SBATCH --error=/scratch/project_2008084/hla_calibration/logs/nci60_wes_full_%j.err

# NCI-60 WES SCALE-UP (42 lines with 2-field A/B/C truth). Reads: SRA SRP150855 (exome capture).
# Truth: Adams 2005 (PMC555742). Frozen WES Champion-Challenger policy. Powers HLA-A (n=39).
# After this completes:
#   python3 bin/run_1000g_benchmark.py --config conf/benchmark_nci60_wes.yaml \
#           --output-dir analysis/nci60_wes_benchmark
set -euo pipefail

module load nextflow

BASE=/scratch/project_2008084/hla_calibration/nci60/wes
FASTQ=$BASE/fastqs
INDEX=$BASE/index/full_fastq_urls.tsv
mkdir -p "$FASTQ" /scratch/project_2008084/hla_calibration/logs /scratch/project_2008084/hla_calibration/work

# ── Download all WES FASTQs, named <SAMPLE>_R1/_R2 so main.nf derives the truth sample IDs ──
while IFS=$'\t' read -r SAMPLE R1 R2; do
  [ -z "$SAMPLE" ] && continue
  echo "[download] $SAMPLE"
  wget -c -q -O "$FASTQ/${SAMPLE}_R1.fastq.gz" "$R1"
  wget -c -q -O "$FASTQ/${SAMPLE}_R2.fastq.gz" "$R2"
done < "$INDEX"
echo "[download] complete:"; du -sh "$FASTQ"

# ── Run the ChampHLA pipeline (DNA/WES tools) on all lines ──
cd /scratch/project_2008084/pihla-publish
nextflow run main.nf \
  -params-file /scratch/project_2008084/pihla-publish/conf/puhti_params.yaml \
  -resume \
  -name nci60_wes_full_${SLURM_JOB_ID} \
  -w /scratch/project_2008084/hla_calibration/work/nci60_wes_full \
  --input "$FASTQ" \
  --input_type fastq \
  --tools optitype,t1k,spechla,hlahd \
  --optitype_seq_type dna \
  --seq_type dna \
  --run_modality wes \
  --outdir "$BASE/results" \
  --slurm_account project_2008084 \
  --use_local_spechla true \
  --spechla_path /projappl/project_2008084/SpecHLAx \
  --singularity_cache_dir /scratch/project_2008084/hla_references/singularity_cache/containers \
  -profile puhti,singularity

echo "[scale] pipeline finished. Results in $BASE/results"

# ── Disk-safe: drop raw FASTQs once per-tool outputs exist (results are tiny txt files) ──
if [ -d "$BASE/results" ] && [ "$(find "$BASE/results" -name '*_optitype.txt' | wc -l)" -ge 30 ]; then
  echo "[cleanup] removing raw FASTQs to free disk"; rm -f "$FASTQ"/*.fastq.gz
fi
echo "[done]"
