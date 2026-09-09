#!/bin/bash
set -euo pipefail

CODE_ROOT=${CHAMPHLA_CODE_ROOT:-/scratch/project_2008084/champhla_plurality}
source "${CHAMPHLA_PATH_CONFIG:-${CODE_ROOT}/configs/roihu_paths.env}"
run_id=${1:?usage: roihu_export_compact.sh RUN_ID FREEZE_VALIDATION_JSON}
freeze_validation=${2:?usage: roihu_export_compact.sh RUN_ID FREEZE_VALIDATION_JSON}
run_root=${CHAMPHLA_RUN_ROOT}/${run_id}
export_root=${CHAMPHLA_RUN_ROOT}/exports
test -d "${run_root}"
test -s "${freeze_validation}"
python3 -c 'import json,sys; assert json.load(open(sys.argv[1]))["valid"] is True' "${freeze_validation}"
mkdir -p "${export_root}"
file_list=${export_root}/${run_id}.files.txt
archive=${export_root}/${run_id}.compact.tar.gz
test ! -e "${archive}"

find "${run_root}" -type f \
  \( -name '*.md' -o -name '*.tsv' -o -name '*.json' -o -name '*.png' \
     -o -name '*.svg' -o -name '*.pdf' -o -name '*.log' -o -name '*.txt' \) \
  -size -100M -printf '%P\n' | LC_ALL=C sort > "${file_list}"
[[ -s "${file_list}" ]]
! grep -Ei '\.(cram|crai|bam|bai|fastq|fq)(\.gz)?$|\.sif$|\.agc$' "${file_list}"
tar --directory "${run_root}" --files-from "${file_list}" --create --gzip --file "${archive}"
sha256sum "${archive}" "${file_list}" > "${archive}.sha256"
