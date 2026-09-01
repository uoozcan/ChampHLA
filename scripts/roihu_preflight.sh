#!/bin/bash
set -euo pipefail

ROOT=/scratch/project_2008084/champhla_publication_candidate
REFERENCE=/scratch/project_2008084/hla_references/genomes/GRCh38_full_analysis_set_plus_decoy_hla.fa

export CSC_ENV_INIT_NON_INTERACTIVE=yes
source /etc/profile.d/zz-csc-env.sh
module load bio-apps/v202603
module load samtools/1.21
module load nextflow/25.10.2-standalone

test -d "${ROOT}"
test -s "${REFERENCE}"
test -s "${REFERENCE}.fai"
python3 --version
samtools --version | head -n 1
nextflow -version
apptainer --version
df -h /scratch

for image in hlahd.sif kourami.sif optitype.sif t1k.sif; do
  test -s "/scratch/project_2008084/hla_references/singularity_cache/containers/${image}"
done

sha256sum "${ROOT}/configs/confirmation_protocol.json" \
  "${ROOT}/configs/runtime_versions.json"

