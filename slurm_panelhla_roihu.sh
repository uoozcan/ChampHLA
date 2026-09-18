#!/bin/bash
#SBATCH --job-name=panelhla
#SBATCH --account=project_2008084
#SBATCH --partition=small
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#
# Run PanelHLA on CSC Roihu.
#
#   sbatch slurm_panelhla_roihu.sh <modality> <input_dir> <input_type> <tools> <outdir> <workdir>
#
# e.g.
#   sbatch slurm_panelhla_roihu.sh wgs  /path/to/bam_dir  bam   optitype,hlahd,t1k  /path/out  /path/work
#   sbatch slurm_panelhla_roihu.sh rnaseq /path/to/fastqs fastq arcashla,seq2hla    /path/out  /path/work
#
# This is the Roihu launcher. The older slurm_*.sh scripts are Puhti-era: they say
# `module load nextflow`, which cannot work here -- see bin/roihu_env.sh for why.
set -euo pipefail

MODALITY="${1:?modality: wgs, wes or rnaseq}"
INPUT="${2:?input directory}"
INPUT_TYPE="${3:?input type: bam, cram or fastq}"
TOOLS="${4:?comma-separated tool list}"
OUTDIR="${5:?output directory}"
WORKDIR="${6:?work directory}"

PIPE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck disable=SC1091
source "${PIPE}/bin/roihu_env.sh"

SEQ_TYPE=dna
EXON_ONLY=0
case "${MODALITY}" in
  rnaseq) SEQ_TYPE=rna ;;
  wes)    EXON_ONLY=1 ;;   # SpecHLA: 1 selects exon-only, per its code
esac

mkdir -p "${OUTDIR}" "${WORKDIR}"

echo "[panelhla] modality=${MODALITY} tools=${TOOLS}"
echo "[panelhla] nextflow=$(command -v nextflow) samtools=$(command -v samtools)"

nextflow run "${PIPE}/main.nf" \
  -params-file "${PIPE}/conf/roihu_params.yaml" \
  -name "panelhla_${MODALITY}_${SLURM_JOB_ID}" \
  -w "${WORKDIR}" \
  --input "${INPUT}" \
  --input_type "${INPUT_TYPE}" \
  --tools "${TOOLS}" \
  --seq_type "${SEQ_TYPE}" \
  --optitype_seq_type "${SEQ_TYPE}" \
  --spechla_exon_only "${EXON_ONLY}" \
  --run_modality "${MODALITY}" \
  --outdir "${OUTDIR}" \
  --slurm_account project_2008084 \
  -profile roihu,singularity

echo "[panelhla] done. results in ${OUTDIR}"
