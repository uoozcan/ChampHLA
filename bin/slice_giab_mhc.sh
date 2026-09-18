#!/bin/bash
# Slice the chr6 MHC region from GIAB GRCh38 Illumina WGS BAMs (remote, no full download) and
# emit per-sample paired FASTQ for the ChampHLA WGS pipeline. Disk-safe: only MHC reads are fetched.
#
# GRCh38 MHC window (with flank): chr6:28,477,797-33,448,354. The GIAB analysis-set BAMs use 'chr6'.
set -euo pipefail
module load samtools 2>/dev/null || true

REGION="chr6:28477797-33448354"
OUT=/scratch/project_2008084/hla_calibration/hprc/wgs/fastqs
TMP=/scratch/project_2008084/hla_calibration/hprc/wgs/slices
mkdir -p "$OUT" "$TMP"

GIAB=https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data
# sample -> remote BAM (GRCh38 2x250 novoalign analysis-set)
declare -A BAM=(
  [HG002]="$GIAB/AshkenazimTrio/HG002_NA24385_son/NIST_Illumina_2x250bps/novoalign_bams/HG002.GRCh38.2x250.bam"
  [HG003]="$GIAB/AshkenazimTrio/HG003_NA24149_father/NIST_Illumina_2x250bps/novoalign_bams/HG003.GRCh38.2x250.bam"
  [HG004]="$GIAB/AshkenazimTrio/HG004_NA24143_mother/NIST_Illumina_2x250bps/novoalign_bams/HG004.GRCh38.2x250.bam"
  [HG005]="$GIAB/ChineseTrio/HG005_NA24631_son/NIST_Illumina_2x250bps/novoalign_bams/HG005.GRCh38.2x250.bam"
  [HG006]="$GIAB/ChineseTrio/HG006_NA24694-huref_father/NIST_Illumina_2x250bps/novoalign_bams/HG006.GRCh38.2x250.bam"
  [HG007]="$GIAB/ChineseTrio/HG007_NA24695-hu38168_mother/NIST_Illumina_2x250bps/novoalign_bams/HG007.GRCh38.2x250.bam"
)

SAMPLES=("${@:-HG002}")   # default HG002; pass sample IDs to slice more
for S in "${SAMPLES[@]}"; do
  URL="${BAM[$S]:-}"
  [ -z "$URL" ] && { echo "[slice] no BAM URL for $S; skip"; continue; }
  echo "[slice] $S  $REGION"
  # name-sort the MHC slice so samtools fastq pairs correctly
  samtools view -b "$URL" "$REGION" 2>/dev/null \
    | samtools sort -n -@4 -o "$TMP/${S}_mhc.namesorted.bam" - 2>/dev/null
  samtools fastq -@4 \
    -1 "$OUT/${S}_R1.fastq.gz" -2 "$OUT/${S}_R2.fastq.gz" \
    -0 /dev/null -s /dev/null -n "$TMP/${S}_mhc.namesorted.bam" 2>/dev/null
  rm -f "$TMP/${S}_mhc.namesorted.bam"
  echo "[slice] $S -> $(du -sh "$OUT/${S}_R1.fastq.gz" "$OUT/${S}_R2.fastq.gz" | cut -f1 | tr '\n' ' ')"
done
echo "[slice] done -> $OUT"
