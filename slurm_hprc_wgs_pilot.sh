#!/bin/bash
#SBATCH --job-name=hprc_wgs
#SBATCH --account=project_2008084
#SBATCH --partition=small
#SBATCH --time=36:00:00
# HEAD/ORCHESTRATOR ONLY: puhti profile uses executor='slurm' (tools run as their own sub-jobs),
# so the parent just runs the Nextflow head — modest resources (the RNA-scale queue lesson).
#SBATCH --cpus-per-task=6
#SBATCH --mem=24G
#SBATCH --output=/scratch/project_2008084/hla_calibration/logs/hprc_wgs_%j.out
#SBATCH --error=/scratch/project_2008084/hla_calibration/logs/hprc_wgs_%j.err

# GIAB/HPRC external validation (WGS). Reads = chr6 MHC slices of open GIAB GRCh38 Illumina WGS,
# pre-staged by bin/slice_giab_mhc.sh (no download here). Frozen WGS Champion-Challenger policy.
# After this completes:
#   python3 bin/run_1000g_benchmark.py --config conf/benchmark_hprc.yaml --output-dir analysis/hprc_benchmark
set -euo pipefail
module load nextflow

BASE=/scratch/project_2008084/hla_calibration/hprc/wgs
FASTQ=$BASE/fastqs
mkdir -p /scratch/project_2008084/hla_calibration/logs /scratch/project_2008084/hla_calibration/work
echo "[hprc] staged samples:"; ls "$FASTQ"/*_R1.fastq.gz 2>/dev/null | sed 's#.*/##;s/_R1.*//'

cd /scratch/project_2008084/pihla-publish
nextflow run main.nf \
  -params-file /scratch/project_2008084/pihla-publish/conf/puhti_params.yaml \
  -resume \
  -name hprc_wgs_${SLURM_JOB_ID} \
  -w /scratch/project_2008084/hla_calibration/work/hprc_wgs \
  --input "$FASTQ" \
  --input_type fastq \
  --tools optitype,t1k,spechla,hlahd \
  --optitype_seq_type dna \
  --seq_type dna \
  --run_modality wgs \
  --outdir "$BASE/results" \
  --slurm_account project_2008084 \
  --use_local_spechla true \
  --spechla_path /projappl/project_2008084/SpecHLAx \
  --singularity_cache_dir /scratch/project_2008084/hla_references/singularity_cache/containers \
  -profile puhti,singularity

echo "[hprc] pipeline finished. Results in $BASE/results"
