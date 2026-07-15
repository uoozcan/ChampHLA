#!/bin/bash
#SBATCH --job-name=nci60_rna_full
#SBATCH --account=project_2008084
#SBATCH --partition=small
# HEAD/ORCHESTRATOR JOB ONLY: the puhti profile uses executor='slurm', so Nextflow submits each
# tool as its own SLURM sub-job. This parent only runs the Nextflow head + wget + seqtk|gzip, so it
# needs modest resources. Over-provisioning (was 20 cpu/180G) buried it behind a 2-day backfill wait.
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=6
#SBATCH --mem=24G
#SBATCH --output=/scratch/project_2008084/hla_calibration/logs/nci60_rna_full_%j.out
#SBATCH --error=/scratch/project_2008084/hla_calibration/logs/nci60_rna_full_%j.err

# NCI-60 RNA SCALE-UP (39 truth lines not in the pilot; +4 pilot = 43). Reads: ENA PRJNA433861.
# DISK-SAFE BATCHED: process CHUNK samples at a time (download -> type -> delete chunk fastqs+work),
# so peak disk stays bounded despite deep total-RNA libraries. Results accumulate in $BASE/results.
# After this completes:
#   python3 bin/run_1000g_benchmark.py --config conf/benchmark_nci60.yaml --output-dir analysis/nci60_benchmark
set -uo pipefail   # NOT -e: a transient bad download must skip its line, not abort the whole job
module load nextflow
module load seqtk

# robust download: retry, require non-empty + gzip-valid (a 0-byte/corrupt mate previously segfaulted seqtk)
dl(){  # url out
  for t in 1 2 3 4; do
    wget -q --tries=3 --timeout=180 -O "$2" "$1" && [ -s "$2" ] && gzip -t "$2" 2>/dev/null && return 0
    echo "   [retry $t] $2"; rm -f "$2"; sleep 8
  done
  return 1
}

BASE=/scratch/project_2008084/hla_calibration/nci60/rna
INDEX=$BASE/index/remaining6_fastq_urls.tsv
RESULTS=$BASE/results
CHUNK=3   # deep total-RNA: ~27 GB/sample (fastq+work); 3/chunk ≈ 80 GB peak, safe under the 1 TB quota
mkdir -p "$RESULTS" /scratch/project_2008084/hla_calibration/logs

mapfile -t LINES < "$INDEX"
N=${#LINES[@]}
echo "[rna-scale] $N samples, chunk size $CHUNK"

i=0; chunk=0
while [ $i -lt $N ]; do
  chunk=$((chunk+1))
  STAGE=$BASE/stage_c${chunk}
  WORK=/scratch/project_2008084/hla_calibration/work/nci60_rna_full_c${chunk}
  rm -rf "$STAGE" "$WORK"; mkdir -p "$STAGE"
  # ── stage this chunk's FASTQs, SUBSAMPLED to 50M read-pairs (pilot depth) ──
  # Deep total-RNA (>100M reads) overruns OptiType/razers3 container /tmp ("quota exceeded").
  # seqtk subsampling to a uniform 50M pairs matches the pilot depth and keeps the footprint small.
  SUBN=50000000
  for j in $(seq 0 $((CHUNK-1))); do
    idx=$((i+j)); [ $idx -ge $N ] && break
    IFS=$'\t' read -r S R1 R2 <<< "${LINES[$idx]}"
    [ -z "$S" ] && continue
    echo "[chunk $chunk][download] $S"
    if ! dl "$R1" "$STAGE/${S}_R1.full.fastq.gz"; then echo "[SKIP] $S R1 failed"; rm -f "$STAGE/${S}"_R*; continue; fi
    if ! dl "$R2" "$STAGE/${S}_R2.full.fastq.gz"; then echo "[SKIP] $S R2 failed"; rm -f "$STAGE/${S}"_R*; continue; fi
    echo "[chunk $chunk][subsample->${SUBN}] $S"
    if seqtk sample -s100 "$STAGE/${S}_R1.full.fastq.gz" $SUBN | gzip > "$STAGE/${S}_R1.fastq.gz" \
       && seqtk sample -s100 "$STAGE/${S}_R2.full.fastq.gz" $SUBN | gzip > "$STAGE/${S}_R2.fastq.gz"; then
      rm -f "$STAGE/${S}_R1.full.fastq.gz" "$STAGE/${S}_R2.full.fastq.gz"   # free the deep originals
    else
      echo "[SKIP] $S subsample failed"; rm -f "$STAGE/${S}"_R*; continue
    fi
  done
  echo "[chunk $chunk] staged: $(du -sh "$STAGE" | cut -f1)"
  # ── type this chunk ──
  cd /scratch/project_2008084/pihla-publish
  nextflow run main.nf \
    -params-file /scratch/project_2008084/pihla-publish/conf/puhti_params.yaml \
    -name nci60_rna_full_c${chunk}_${SLURM_JOB_ID} \
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
  # ── free disk before next chunk ──
  echo "[chunk $chunk] done; cleaning stage + work"
  rm -rf "$STAGE" "$WORK"
  i=$((i+CHUNK))
done
echo "[rna-scale] all chunks done. Results in $RESULTS"
