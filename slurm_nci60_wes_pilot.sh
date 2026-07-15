#!/bin/bash
#SBATCH --job-name=nci60_wes_pilot
#SBATCH --account=project_2008084
#SBATCH --partition=small
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=20
#SBATCH --mem=180G
#SBATCH --output=/scratch/project_2008084/hla_calibration/logs/nci60_wes_pilot_%j.out
#SBATCH --error=/scratch/project_2008084/hla_calibration/logs/nci60_wes_pilot_%j.err

# NCI-60 WES external-validation PILOT (5 lines with full A/B/C 2-field truth).
# Reads: SRA SRP150855 (NCI-60 exome capture). Truth: Adams 2005 (PMC555742). Pipeline: ChampHLA main.nf (DNA/WES).
# Second modality on the same lines as the validated RNA pilot (NCI-H23, OVCAR-8, RPMI-8226, SK-MEL-28) + EKVX.
# After this completes, build sequencing_source.tsv and run:
#   python3 bin/run_1000g_benchmark.py --config conf/benchmark_nci60_wes.yaml \
#           --output-dir analysis/nci60_wes_benchmark
set -euo pipefail

module load nextflow

BASE=/scratch/project_2008084/hla_calibration/nci60/wes
FASTQ=$BASE/fastqs
INDEX=$BASE/index/pilot_fastq_urls.tsv
mkdir -p "$FASTQ" /scratch/project_2008084/hla_calibration/logs /scratch/project_2008084/hla_calibration/work

# ── Download pilot FASTQs, named <SAMPLE>_R1/_R2 so main.nf derives the truth sample IDs ──
while IFS=$'\t' read -r SAMPLE R1 R2; do
  [ -z "$SAMPLE" ] && continue
  echo "[download] $SAMPLE"
  wget -c -q -O "$FASTQ/${SAMPLE}_R1.fastq.gz" "$R1"
  wget -c -q -O "$FASTQ/${SAMPLE}_R2.fastq.gz" "$R2"
done < "$INDEX"
echo "[download] complete:"; ls -lh "$FASTQ"

# ── Run the ChampHLA pipeline (DNA/WES tools) ──
cd /scratch/project_2008084/pihla-publish
nextflow run main.nf \
  -params-file /scratch/project_2008084/pihla-publish/conf/puhti_params.yaml \
  -resume \
  -name nci60_wes_pilot_${SLURM_JOB_ID} \
  -w /scratch/project_2008084/hla_calibration/work/nci60_wes_pilot \
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

echo "[pilot] pipeline finished. Results in $BASE/results"
