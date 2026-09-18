#!/bin/bash
#SBATCH --job-name=pihla_rna_3
#SBATCH --account=project_2008084
#SBATCH --partition=small
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=20
#SBATCH --mem=180G
#SBATCH --output=/scratch/project_2008084/hla_calibration/logs/rna_3sample_%j.out
#SBATCH --error=/scratch/project_2008084/hla_calibration/logs/rna_3sample_%j.err

set -euo pipefail

module load nextflow
mkdir -p /scratch/project_2008084/hla_calibration/logs
mkdir -p /scratch/project_2008084/hla_calibration/work
mkdir -p /scratch/project_2008084/hla_calibration/rna_3sample_input_named
cd /scratch/project_2008084/hla_calibration/rna_3sample_input_named

ln -sf /scratch/project_2008084/hla_calibration/rna_3sample_input/ERR188327_1.fastq.gz NA06985_R1.fastq.gz
ln -sf /scratch/project_2008084/hla_calibration/rna_3sample_input/ERR188327_2.fastq.gz NA06985_R2.fastq.gz
ln -sf /scratch/project_2008084/hla_calibration/rna_3sample_input/ERR188213_1.fastq.gz NA06986_R1.fastq.gz
ln -sf /scratch/project_2008084/hla_calibration/rna_3sample_input/ERR188213_2.fastq.gz NA06986_R2.fastq.gz
ln -sf /scratch/project_2008084/hla_calibration/rna_3sample_input/ERR188047_1.fastq.gz NA06994_R1.fastq.gz
ln -sf /scratch/project_2008084/hla_calibration/rna_3sample_input/ERR188047_2.fastq.gz NA06994_R2.fastq.gz

cd /scratch/project_2008084/pihla-publish

nextflow run main.nf \
  -name pihla_rna_3_${SLURM_JOB_ID} \
  -w /scratch/project_2008084/hla_calibration/work/rna_3sample_${SLURM_JOB_ID} \
  --input /scratch/project_2008084/hla_calibration/rna_3sample_input_named \
  --input_type fastq \
  --tools arcashla,optitype,seq2hla,t1k,spechla,hlahd \
  --optitype_seq_type rna \
  --seq_type rna \
  --run_modality rnaseq \
  --outdir /scratch/project_2008084/hla_calibration/rna_3sample/results \
  --slurm_account project_2008084 \
  --use_local_spechla true \
  --spechla_path /projappl/project_2008084/SpecHLAx \
  -profile puhti,singularity
