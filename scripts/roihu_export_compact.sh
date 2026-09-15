#!/bin/bash
set -euo pipefail

CODE_ROOT=${CHAMPHLA_CODE_ROOT:-/scratch/project_2008084/champhla_plurality}
source "${CHAMPHLA_PATH_CONFIG:-${CODE_ROOT}/configs/roihu_paths.env}"
run_id=${1:?usage: roihu_export_compact.sh RUN_ID FREEZE_VALIDATION_JSON}
freeze_validation=${2:?usage: roihu_export_compact.sh RUN_ID FREEZE_VALIDATION_JSON}
export_root=${CHAMPHLA_RUN_ROOT}/exports
policy=${CHAMPHLA_RELEASE_ALLOWLIST:-${CODE_ROOT}/configs/release_allowlist.json}
test -d "${CHAMPHLA_RUN_ROOT}"
test -s "${freeze_validation}"
python3 -c 'import json,sys; assert json.load(open(sys.argv[1]))["valid"] is True' "${freeze_validation}"
mkdir -p "${export_root}"
file_list=${export_root}/${run_id}.files.txt
export_manifest=${export_root}/${run_id}.export_manifest.json
archive=${export_root}/${run_id}.compact.tar.gz
test ! -e "${archive}" && test ! -e "${file_list}" && test ! -e "${export_manifest}"

python3 -c 'from champhla_recovery.cli import build_compact_export_manifest_main; raise SystemExit(build_compact_export_manifest_main())' \
  --root "${CHAMPHLA_RUN_ROOT}" --run-id "${run_id}" --policy "${policy}" \
  --output "${export_manifest}" --list-output "${file_list}"
python3 -c 'from champhla_recovery.cli import verify_compact_export_manifest_main; raise SystemExit(verify_compact_export_manifest_main())' \
  --root "${CHAMPHLA_RUN_ROOT}" --manifest "${export_manifest}"
tar --directory "${CHAMPHLA_RUN_ROOT}" --verbatim-files-from --files-from "${file_list}" \
  --create --gzip --file "${archive}"
sha256sum "${archive}" "${file_list}" "${export_manifest}" > "${archive}.sha256"
