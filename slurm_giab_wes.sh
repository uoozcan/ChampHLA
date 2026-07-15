#!/bin/bash
#SBATCH --job-name=giab_wes
#SBATCH --account=project_2008084
#SBATCH --partition=small
#SBATCH --time=36:00:00
# HEAD/ORCHESTRATOR ONLY (puhti profile = slurm executor; tools run as sub-jobs).
#SBATCH --cpus-per-task=6
#SBATCH --mem=24G
#SBATCH --output=/scratch/project_2008084/hla_calibration/logs/giab_wes_%j.out
#SBATCH --error=/scratch/project_2008084/hla_calibration/logs/giab_wes_%j.err

# GIAB Ashkenazi-trio WES (HG002/HG003/HG004) vs Stanford clinical SBT gold truth. Deep exome
# (~10 GB/sample). Frozen/deployed WES Champion-Challenger policy. Post-type FASTQ cleanup to bound disk.
# After: python3 bin/run_1000g_benchmark.py --config conf/benchmark_giab_wes.yaml \
#          --output-dir analysis/giab_trio_benchmark/wes
set -uo pipefail
module load nextflow

BASE=/scratch/project_2008084/hla_calibration/giab_trio/wes
STAGE=$BASE/fastqs; RESULTS=$BASE/results
INDEX=/scratch/project_2008084/pihla-publish/analysis/giab_trio_benchmark/index/wes_fastq_urls.tsv
WORK=/scratch/project_2008084/hla_calibration/work/giab_wes
rm -rf "$WORK"; mkdir -p "$STAGE" "$RESULTS" /scratch/project_2008084/hla_calibration/logs

dl(){  # url out
  for t in 1 2 3 4; do
    wget -q --tries=3 --timeout=300 -O "$2" "$1" && [ -s "$2" ] && gzip -t "$2" 2>/dev/null && return 0
    echo "   [retry $t] $2"; rm -f "$2"; sleep 10
  done
  return 1
}
while IFS=$'\t' read -r S R1 R2; do
  [ -z "$S" ] && continue
  echo "[download] $S"
  if ! dl "$R1" "$STAGE/${S}_R1.fastq.gz"; then echo "[SKIP] $S R1 failed"; rm -f "$STAGE/${S}"_R*; continue; fi
  if ! dl "$R2" "$STAGE/${S}_R2.fastq.gz"; then echo "[SKIP] $S R2 failed"; rm -f "$STAGE/${S}"_R*; continue; fi
done < "$INDEX"
echo "[staged] $(ls "$STAGE"/*_R1.fastq.gz 2>/dev/null | wc -l) samples, $(du -sh "$STAGE" 2>/dev/null | cut -f1)"

cd /scratch/project_2008084/pihla-publish
nextflow run main.nf \
  -params-file /scratch/project_2008084/pihla-publish/conf/puhti_params.yaml \
  -resume -name giab_wes_${SLURM_JOB_ID} \
  -w "$WORK" \
  --input "$STAGE" \
  --input_type fastq \
  --tools optitype,t1k,spechla,hlahd \
  --optitype_seq_type dna \
  --seq_type dna \
  --run_modality wes \
  --outdir "$RESULTS" \
  --slurm_account project_2008084 \
  --use_local_spechla true \
  --spechla_path /projappl/project_2008084/SpecHLAx \
  --singularity_cache_dir /scratch/project_2008084/hla_references/singularity_cache/containers \
  -profile puhti,singularity

echo "[giab-wes] pipeline finished. Results in $RESULTS"
# disk-safe: drop raw FASTQs once tool outputs exist
if [ "$(find "$RESULTS" -name '*_optitype.txt' | wc -l)" -ge 2 ]; then rm -f "$STAGE"/*.fastq.gz; fi
