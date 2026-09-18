#!/bin/bash
#SBATCH --job-name=ihwg_wgs
#SBATCH --account=project_2008084
#SBATCH --partition=small
#SBATCH --time=24:00:00
# HEAD/ORCHESTRATOR ONLY (puhti profile = slurm executor; tools run as sub-jobs).
#SBATCH --cpus-per-task=6
#SBATCH --mem=24G
#SBATCH --output=/scratch/project_2008084/hla_calibration/logs/ihwg_wgs_%j.out
#SBATCH --error=/scratch/project_2008084/hla_calibration/logs/ihwg_wgs_%j.err

# IHWG/IHIW WGS pilot — n=4 study-verified MHC reference cells (SSTO, DBB, APD, QBL; PRJNA764575,
# cell_line-annotated Illumina paired short-read). Truth: IPD-IMGT IHIW multi-lab 4-field consensus.
# Small MHC-enriched libraries (~50 MB/run) → no subsampling needed.
# After: python3 bin/run_1000g_benchmark.py --config conf/benchmark_ihwg_wgs.yaml \
#          --output-dir analysis/ihwg_benchmark/wgs
set -uo pipefail   # NOT -e: a bad download skips its line, not the whole job
module load nextflow

BASE=/scratch/project_2008084/hla_calibration/ihwg/wgs
STAGE=$BASE/stage; RESULTS=$BASE/results
INDEX=/scratch/project_2008084/pihla-publish/analysis/ihwg_benchmark/index/wgs_fastq_urls.tsv
WORK=/scratch/project_2008084/hla_calibration/work/ihwg_wgs
rm -rf "$STAGE" "$WORK"; mkdir -p "$STAGE" "$RESULTS" /scratch/project_2008084/hla_calibration/logs

dl(){  # url out
  for t in 1 2 3 4; do
    wget -q --tries=3 --timeout=180 -O "$2" "$1" && [ -s "$2" ] && gzip -t "$2" 2>/dev/null && return 0
    echo "   [retry $t] $2"; rm -f "$2"; sleep 8
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
  -resume -name ihwg_wgs_${SLURM_JOB_ID} \
  -w "$WORK" \
  --input "$STAGE" \
  --input_type fastq \
  --tools optitype,t1k,spechla,hlahd \
  --optitype_seq_type dna \
  --seq_type dna \
  --run_modality wgs \
  --outdir "$RESULTS" \
  --slurm_account project_2008084 \
  --use_local_spechla true \
  --spechla_path /projappl/project_2008084/SpecHLAx \
  --singularity_cache_dir /scratch/project_2008084/hla_references/singularity_cache/containers \
  -profile puhti,singularity

echo "[ihwg-wgs] pipeline finished. Results in $RESULTS"
