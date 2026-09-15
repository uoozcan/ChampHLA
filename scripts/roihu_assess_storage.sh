#!/bin/bash
set -euo pipefail

CODE_ROOT=${CHAMPHLA_CODE_ROOT:-/scratch/project_2008084/champhla_plurality}
source "${CHAMPHLA_PATH_CONFIG:-${CODE_ROOT}/configs/roihu_paths.env}"
gate_id=${1:?usage: roihu_assess_storage.sh GATE_ID WGS_LEDGER WES_LEDGER RNA_LEDGER}
wgs_ledger=${2:?}
wes_ledger=${3:?}
rna_ledger=${4:?}
[[ "${gate_id}" =~ ^[a-z0-9][a-z0-9._-]{2,63}$ ]] || {
  echo "unsafe gate_id: ${gate_id}" >&2
  exit 2
}
export PYTHONPATH=${CHAMPHLA_CODE_ROOT}/src

preflight_pointer=${CHAMPHLA_RUN_ROOT}/preflight/CURRENT
test -s "${preflight_pointer}"
test ! -L "${preflight_pointer}"
[[ "$(wc -l < "${preflight_pointer}")" -eq 1 ]]
preflight_id=$(sed -n '1p' "${preflight_pointer}")
[[ "${preflight_id}" =~ ^[a-z0-9][a-z0-9._-]{2,63}$ ]]
inventory=${CHAMPHLA_RUN_ROOT}/preflight/${preflight_id}/environment_inventory.json
test -s "${inventory}"
test ! -L "${inventory}"

gate_parent=${CHAMPHLA_RUN_ROOT}/storage_gates
gate_root=${gate_parent}/${gate_id}
mkdir -p "${gate_parent}"
test ! -e "${gate_root}"
mkdir "${gate_root}"
python3 -c 'from champhla_confirmation.cli import assess_roihu_storage_main; raise SystemExit(assess_roihu_storage_main())' \
  --pilot-ledger "${wgs_ledger}" --pilot-ledger "${wes_ledger}" \
  --pilot-ledger "${rna_ledger}" \
  --targets "${CHAMPHLA_CODE_ROOT}/configs/roihu_storage_targets.json" \
  --environment-inventory "${inventory}" --evidence-id "${gate_id}" \
  --output "${gate_root}/storage_gate.json"
(cd "${gate_root}" && sha256sum storage_gate.json > evidence.sha256)
printf 'storage gate created but not promoted: %s\n' "${gate_root}"
