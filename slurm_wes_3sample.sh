#!/bin/bash
#SBATCH --job-name=pihla_wes_3
#SBATCH --account=project_2008084
#SBATCH --partition=small
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=20
#SBATCH --mem=180G
#SBATCH --output=/scratch/project_2008084/hla_calibration/logs/wes_3sample_%j.out
#SBATCH --error=/scratch/project_2008084/hla_calibration/logs/wes_3sample_%j.err

set -euo pipefail

module load nextflow
module load samtools
mkdir -p /scratch/project_2008084/hla_calibration/logs
mkdir -p /scratch/project_2008084/hla_calibration/work

for bam in /scratch/project_2008084/hla_calibration/wes_3sample_input/*.bam; do
  [ -f "${bam}.bai" ] || samtools index -@ 8 "$bam"
done

cd /scratch/project_2008084/pihla-publish

nextflow run main.nf \
  -name pihla_wes_3_${SLURM_JOB_ID} \
  -w /scratch/project_2008084/hla_calibration/work/wes_3sample_${SLURM_JOB_ID} \
  --input /scratch/project_2008084/hla_calibration/wes_3sample_input \
  --input_type bam \
  --tools spechla,hlahd,optitype,polysolver,kourami,t1k,arcashla \
  --spechla_exon_only 1 \
  --seq_type dna \
  --run_modality wes \
  --outdir /scratch/project_2008084/hla_calibration/wes_3sample/results \
  --slurm_account project_2008084 \
  --use_local_spechla true \
  --spechla_path /projappl/project_2008084/SpecHLAx \
  -profile puhti,singularity
