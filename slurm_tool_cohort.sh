#!/bin/bash
#SBATCH --job-name=pihla_tool_cohort
#SBATCH --account=project_2008084
#SBATCH --partition=small
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --output=/scratch/project_2008084/hla_calibration/logs/tool_cohort_%x_%j.out
#SBATCH --error=/scratch/project_2008084/hla_calibration/logs/tool_cohort_%x_%j.err

# Submit per-sample, single-tool Nextflow jobs for an entire cohort in batches.
#
# Usage (as a SLURM job or directly on login node):
#   sbatch slurm_tool_cohort.sh --tool hlahd --modality wes [--batch-size 10]
#   bash   slurm_tool_cohort.sh --tool arcashla --modality wgs
#
# Arguments:
#   --tool         Tool name: optitype|hlahd|t1k|arcashla|polysolver|kourami|spechla|seq2hla
#   --modality     rna|wes|wgs
#   --batch-size   Max concurrent jobs (default: 10)
#   --input-dir    Override default input directory
#   --output-dir   Override default output root
#   --work-dir     Override default Nextflow work root

set -euo pipefail

TOOL=""
MODALITY=""
BATCH_SIZE=10
INPUT_DIR=""
OUTPUT_DIR=""
WORK_ROOT=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --tool)        TOOL="$2";       shift 2 ;;
    --modality)    MODALITY="$2";   shift 2 ;;
    --batch-size)  BATCH_SIZE="$2"; shift 2 ;;
    --input-dir)   INPUT_DIR="$2";  shift 2 ;;
    --output-dir)  OUTPUT_DIR="$2"; shift 2 ;;
    --work-dir)    WORK_ROOT="$2";  shift 2 ;;
    *) echo "Unknown arg: $1" >&2; exit 1 ;;
  esac
done

[[ -z "$TOOL" || -z "$MODALITY" ]] && {
  echo "Usage: $0 --tool TOOL --modality rna|wes|wgs [--batch-size N]" >&2
  exit 1
}

BASE=/scratch/project_2008084/hla_calibration
PIPELINE=/scratch/project_2008084/pihla-publish/main.nf

# Defaults by modality
case "$MODALITY" in
  rna)
    INPUT_DIR="${INPUT_DIR:-${BASE}/rna/fastqs}"
    OUTPUT_DIR="${OUTPUT_DIR:-${BASE}/rna_batches}"
    SEQ_TYPE=rna
    RUN_MODALITY=rnaseq
    SPECHLA_EXON=0
    INPUT_TYPE=fastq
    ;;
  wes)
    INPUT_DIR="${INPUT_DIR:-${BASE}/wes/bams}"
    OUTPUT_DIR="${OUTPUT_DIR:-${BASE}/wes_batches}"
    SEQ_TYPE=dna
    RUN_MODALITY=wes
    SPECHLA_EXON=1
    INPUT_TYPE=bam
    ;;
  wgs)
    INPUT_DIR="${INPUT_DIR:-${BASE}/wgs/bams_new}"
    OUTPUT_DIR="${OUTPUT_DIR:-${BASE}/wgs_batches}"
    SEQ_TYPE=dna
    RUN_MODALITY=wgs
    SPECHLA_EXON=0
    INPUT_TYPE=bam
    ;;
  *)
    echo "Unknown modality: $MODALITY" >&2; exit 1 ;;
esac

WORK_ROOT="${WORK_ROOT:-${BASE}/work}"
LOGS_DIR="${BASE}/logs"
mkdir -p "$LOGS_DIR"

JOB_NAME="pihla_${MODALITY}_${TOOL}"

# Build sample list from input dir
declare -a SAMPLES=()
if [[ "$INPUT_TYPE" == "fastq" ]]; then
  while IFS= read -r r1; do
    s=$(basename "$r1" _R1.fastq.gz)
    SAMPLES+=("$s")
  done < <(find "$INPUT_DIR" -maxdepth 1 -name "*_R1.fastq.gz" | sort)
else
  while IFS= read -r bam; do
    s=$(basename "$bam" .bam | grep -oP '^[A-Z0-9]+')
    [[ -n "$s" ]] && SAMPLES+=("$s")
  done < <(find "$INPUT_DIR" -maxdepth 1 -name "*.bam" | sort)
fi

echo "[$(date)] tool=$TOOL modality=$MODALITY batch_size=$BATCH_SIZE"
echo "[$(date)] Found ${#SAMPLES[@]} samples in $INPUT_DIR"

submitted=0; skipped=0

for sample in "${SAMPLES[@]}"; do
  # Skip if result already published
  sentinel="${OUTPUT_DIR}/${sample}/results/${sample}/${TOOL}/${sample}_${TOOL}.txt"
  if [[ -f "$sentinel" ]]; then
    (( skipped++ )) || true
    continue
  fi

  # Find input file(s)
  if [[ "$INPUT_TYPE" == "fastq" ]]; then
    r1="${INPUT_DIR}/${sample}_R1.fastq.gz"
    r2="${INPUT_DIR}/${sample}_R2.fastq.gz"
    [[ -f "$r1" && -f "$r2" ]] || { echo "  SKIP $sample: FASTQs missing"; (( skipped++ )) || true; continue; }
    SS_CONTENT="sample_id,fastq_1,fastq_2
${sample},${r1},${r2}"
  else
    bam=$(find "$INPUT_DIR" -maxdepth 1 -name "${sample}*.bam" | head -1)
    [[ -n "$bam" ]] || { echo "  SKIP $sample: BAM missing"; (( skipped++ )) || true; continue; }
    SS_CONTENT="sample_id,bam_path
${sample},${bam}"
  fi

  # Rate-limit: wait until fewer than BATCH_SIZE jobs with this name are running/pending
  while [[ $(squeue -u ozcanumu -h --name="$JOB_NAME" 2>/dev/null | wc -l) -ge $BATCH_SIZE ]]; do
    sleep 30
  done

  # Create per-tool launch dir and samplesheet
  LAUNCH_DIR="${OUTPUT_DIR}/${sample}/nxf_${TOOL}_launch"
  mkdir -p "$LAUNCH_DIR"
  SS="${OUTPUT_DIR}/${sample}/${MODALITY}_${sample}_${TOOL}_samplesheet.csv"
  echo "$SS_CONTENT" > "$SS"

  WORK_DIR="${WORK_ROOT}/${MODALITY}_${TOOL}_${sample}"

  job_id=$(sbatch --parsable \
    --job-name="$JOB_NAME" \
    --account=project_2008084 \
    --partition=small \
    --time=12:00:00 \
    --cpus-per-task=2 \
    --mem=8G \
    --output="${LOGS_DIR}/${MODALITY}_${TOOL}_${sample}_%j.out" \
    --error="${LOGS_DIR}/${MODALITY}_${TOOL}_${sample}_%j.err" \
    --wrap="
      cd '${LAUNCH_DIR}'
      module load nextflow
      module load samtools
      nextflow run '${PIPELINE}' \
        -name pihla_${MODALITY}_${TOOL}_${sample}_\${SLURM_JOB_ID} \
        -w '${WORK_DIR}' \
        --input_samplesheet '${SS}' \
        --input_type '${INPUT_TYPE}' \
        --tools '${TOOL}' \
        --seq_type '${SEQ_TYPE}' \
        --run_modality '${RUN_MODALITY}' \
        --spechla_exon_only '${SPECHLA_EXON}' \
        --outdir '${OUTPUT_DIR}/${sample}/results' \
        --slurm_account project_2008084 \
        --use_local_spechla true \
        --spechla_path /projappl/project_2008084/SpecHLAx \
        -profile puhti,singularity
    ")

  echo "  Submitted $sample → job $job_id"
  (( submitted++ )) || true
done

echo "[$(date)] Done: $submitted submitted, $skipped skipped (already have results or missing input)"
