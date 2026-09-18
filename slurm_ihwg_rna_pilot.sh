#!/bin/bash
#SBATCH --job-name=ihwg_rna_pilot
#SBATCH --account=project_2008084
#SBATCH --partition=small
#SBATCH --time=24:00:00
# HEAD/ORCHESTRATOR ONLY (puhti profile = slurm executor; tools run as sub-jobs). Modest resources.
#SBATCH --cpus-per-task=6
#SBATCH --mem=32G
#SBATCH --output=/scratch/project_2008084/hla_calibration/logs/ihwg_rna_pilot_%j.out
#SBATCH --error=/scratch/project_2008084/hla_calibration/logs/ihwg_rna_pilot_%j.err

# IHWG/IHIW RNA pilot — 5 clean non-1000G reference B-LCLs (JY, MOLT-4, BJ, PAR, BER).
# Reads: cell_line-verified public RNA-seq (ENA). Truth: IPD-IMGT/HLA IHIW multi-lab consensus.
# Subsample deep libraries to 50M read-pairs (NCI-60 lesson: deep RNA overruns OptiType /tmp).
# After: python3 bin/run_1000g_benchmark.py --config conf/benchmark_ihwg_rna.yaml \
#          --output-dir analysis/ihwg_benchmark/rna
set -uo pipefail   # NOT -e: a single bad download must skip its line, not abort the job
module load nextflow
module load seqtk

BASE=/scratch/project_2008084/hla_calibration/ihwg/rna
STAGE=$BASE/stage; RESULTS=$BASE/results
INDEX=/scratch/project_2008084/pihla-publish/analysis/ihwg_benchmark/index/rna_pilot_fastq_urls.tsv
WORK=/scratch/project_2008084/hla_calibration/work/ihwg_rna_pilot
SUBN=50000000
rm -rf "$STAGE" "$WORK"; mkdir -p "$STAGE" "$RESULTS" /scratch/project_2008084/hla_calibration/logs

# robust download: retry, require non-empty + gzip-valid
dl(){  # url out
  for t in 1 2 3 4; do
    wget -q --tries=3 --timeout=180 -O "$2" "$1" && [ -s "$2" ] && gzip -t "$2" 2>/dev/null && return 0
    echo "   [retry $t] $2"; rm -f "$2"; sleep 8
  done
  return 1
}
while IFS=$'\t' read -r S R1 R2 RC; do
  [ -z "$S" ] && continue
  echo "[download] $S"
  if ! dl "$R1" "$STAGE/${S}_R1.full.fastq.gz"; then echo "[SKIP] $S R1 download failed"; rm -f "$STAGE/${S}"_R*; continue; fi
  if ! dl "$R2" "$STAGE/${S}_R2.full.fastq.gz"; then echo "[SKIP] $S R2 download failed"; rm -f "$STAGE/${S}"_R*; continue; fi
  # cap depth at ~50M pairs using FRACTION-based seqtk (streaming, O(1) memory; a count-based
  # subsample reservoir-holds 50M reads in RAM and OOM-kills the head job on deep libraries).
  FRAC=$(awk -v rc="${RC:-0}" -v n="$SUBN" 'BEGIN{ if(rc<=0||rc<=n) print "1.0"; else printf "%.4f", n/rc }')
  echo "[subsample frac=$FRAC (rc=${RC:-NA})] $S"
  if [ "$FRAC" = "1.0" ]; then
    mv "$STAGE/${S}_R1.full.fastq.gz" "$STAGE/${S}_R1.fastq.gz"
    mv "$STAGE/${S}_R2.full.fastq.gz" "$STAGE/${S}_R2.fastq.gz"
  elif seqtk sample -s100 "$STAGE/${S}_R1.full.fastq.gz" "$FRAC" | gzip > "$STAGE/${S}_R1.fastq.gz" \
       && seqtk sample -s100 "$STAGE/${S}_R2.full.fastq.gz" "$FRAC" | gzip > "$STAGE/${S}_R2.fastq.gz"; then
    rm -f "$STAGE/${S}_R1.full.fastq.gz" "$STAGE/${S}_R2.full.fastq.gz"
  else
    echo "[SKIP] $S subsample failed"; rm -f "$STAGE/${S}"_R*; continue
  fi
done < "$INDEX"
echo "[staged] $(ls "$STAGE"/*_R1.fastq.gz 2>/dev/null | wc -l) samples, $(du -sh "$STAGE" 2>/dev/null | cut -f1)"

cd /scratch/project_2008084/pihla-publish
nextflow run main.nf \
  -params-file /scratch/project_2008084/pihla-publish/conf/puhti_params.yaml \
  -resume -name ihwg_rna_pilot_${SLURM_JOB_ID} \
  -w "$WORK" \
  --input "$STAGE" \
  --input_type fastq \
  --tools arcashla,optitype,seq2hla,t1k,hlahd \
  --optitype_seq_type rna \
  --seq_type rna \
  --run_modality rnaseq \
  --outdir "$RESULTS" \
  --slurm_account project_2008084 \
  --use_local_spechla true \
  --spechla_path /projappl/project_2008084/SpecHLAx \
  --singularity_cache_dir /scratch/project_2008084/hla_references/singularity_cache/containers \
  -profile puhti,singularity

echo "[ihwg-rna] pipeline finished. Results in $RESULTS"
rm -rf "$STAGE" "$WORK"
