#!/bin/bash
set -eo pipefail

ROOT=${CHAMPHLA_CODE_ROOT:-/scratch/project_2008084/champhla_plurality}
PATH_CONFIG=${CHAMPHLA_PATH_CONFIG:-${ROOT}/configs/roihu_paths.env}

export CSC_ENV_INIT_NON_INTERACTIVE=yes
source /etc/profile.d/zz-csc-env.sh
set -u
module load bio-apps/v202603
module load samtools/1.21
module load nextflow/25.10.2-standalone

test -s "${PATH_CONFIG}"
source "${PATH_CONFIG}"
test -d "${CHAMPHLA_CODE_ROOT}"
test -d "${CHAMPHLA_RUN_ROOT}"
test -d "${CHAMPHLA_INPUT_ROOT}"
test -d "${CHAMPHLA_TRUTH_ROOT}"
test -s "${CHAMPHLA_REFERENCE_ROOT}/genomes/GRCh38_full_analysis_set_plus_decoy_hla.fa"
test -s "${CHAMPHLA_REFERENCE_ROOT}/genomes/GRCh38_full_analysis_set_plus_decoy_hla.fa.fai"
python3 --version
samtools --version | head -n 1
nextflow -version
apptainer --version
df -h /scratch
du -sh "${CHAMPHLA_CODE_ROOT}" "${CHAMPHLA_RUN_ROOT}" "${CHAMPHLA_INPUT_ROOT}" \
  "${CHAMPHLA_TRUTH_ROOT}" "${CHAMPHLA_REFERENCE_ROOT}"
module list
sinfo --format='%P %a %l %D %c %m %G'

export PYTHONPATH="${CHAMPHLA_CODE_ROOT}/src"
python3 -c 'from champhla_confirmation.cli import audit_roihu_environment_main; raise SystemExit(audit_roihu_environment_main())' \
  --site-config "${CHAMPHLA_CODE_ROOT}/configs/roihu_site.json" \
  --output "${CHAMPHLA_RUN_ROOT}/environment_inventory.json"

sha256sum "${CHAMPHLA_CODE_ROOT}/configs/confirmation_protocol.json" \
  "${CHAMPHLA_CODE_ROOT}/configs/runtime_versions.json" \
  "${CHAMPHLA_CODE_ROOT}/configs/roihu_site.json"
