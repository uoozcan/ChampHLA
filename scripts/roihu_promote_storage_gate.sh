#!/bin/bash
set -euo pipefail

CODE_ROOT=${CHAMPHLA_CODE_ROOT:-/scratch/project_2008084/champhla_plurality}
source "${CHAMPHLA_PATH_CONFIG:-${CODE_ROOT}/configs/roihu_paths.env}"
gate_id=${1:?usage: roihu_promote_storage_gate.sh GATE_ID}
[[ "${gate_id}" =~ ^[a-z0-9][a-z0-9._-]{2,63}$ ]] || {
  echo "unsafe gate_id: ${gate_id}" >&2
  exit 2
}

gate_parent=${CHAMPHLA_RUN_ROOT}/storage_gates
gate_root=${gate_parent}/${gate_id}
gate=${gate_root}/storage_gate.json
test -s "${gate}"
test ! -L "${gate_root}"
test ! -L "${gate}"
(cd "${gate_root}" && sha256sum --check --status evidence.sha256)

preflight_pointer=${CHAMPHLA_RUN_ROOT}/preflight/CURRENT
test -s "${preflight_pointer}"
test ! -L "${preflight_pointer}"
[[ "$(wc -l < "${preflight_pointer}")" -eq 1 ]]
preflight_id=$(sed -n '1p' "${preflight_pointer}")
[[ "${preflight_id}" =~ ^[a-z0-9][a-z0-9._-]{2,63}$ ]]
inventory=${CHAMPHLA_RUN_ROOT}/preflight/${preflight_id}/environment_inventory.json
test -s "${inventory}"
python3 -c 'import hashlib,json,sys; g=json.load(open(sys.argv[1])); observed=hashlib.sha256(open(sys.argv[2],"rb").read()).hexdigest(); assert g["schema_version"] == "champhla-roihu-storage-gate-4" and g["evidence_id"] == sys.argv[3] and g["passed"] is True and g["environment_inventory_sha256"] == observed' \
  "${gate}" "${inventory}" "${gate_id}"

pointer=${gate_parent}/CURRENT
pointer_tmp=$(mktemp "${gate_parent}/.CURRENT.XXXXXX")
printf '%s\n' "${gate_id}" > "${pointer_tmp}"
mv -f "${pointer_tmp}" "${pointer}"
printf 'promoted storage gate: %s\n' "${gate_id}"
