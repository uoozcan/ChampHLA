#!/bin/bash
set -euo pipefail

CODE_ROOT=${CHAMPHLA_CODE_ROOT:-/scratch/project_2008084/champhla_plurality}
source "${CHAMPHLA_PATH_CONFIG:-${CODE_ROOT}/configs/roihu_paths.env}"
evidence_id=${1:?usage: roihu_promote_preflight.sh EVIDENCE_ID}
[[ "${evidence_id}" =~ ^[a-z0-9][a-z0-9._-]{2,63}$ ]] || {
  echo "unsafe evidence_id: ${evidence_id}" >&2
  exit 2
}
export PYTHONPATH=${CHAMPHLA_CODE_ROOT}/src

preflight_parent=${CHAMPHLA_RUN_ROOT}/preflight
evidence_root=${preflight_parent}/${evidence_id}
inventory=${evidence_root}/environment_inventory.json
workflow_audit=${evidence_root}/workflow_lock_audit.json
checksum_manifest=${evidence_root}/evidence.sha256
test -s "${inventory}"
test -s "${workflow_audit}"
test -s "${checksum_manifest}"
(cd "${evidence_root}" && sha256sum --check --status evidence.sha256)
python3 -c 'from champhla_confirmation.roihu import validate_environment_inventory; import sys; f=validate_environment_inventory(sys.argv[1], project_root=sys.argv[2], max_age_seconds=3600); assert not f, f' \
  "${inventory}" "${CHAMPHLA_CODE_ROOT}"
python3 -c 'import json,sys; a=json.load(open(sys.argv[1])); assert a["passed"] is True' \
  "${workflow_audit}"
python3 -c 'import json,sys; assert json.load(open(sys.argv[1]))["evidence_id"] == sys.argv[2]' \
  "${inventory}" "${evidence_id}"

pointer=${preflight_parent}/CURRENT
pointer_tmp=$(mktemp "${preflight_parent}/.CURRENT.XXXXXX")
printf '%s\n' "${evidence_id}" > "${pointer_tmp}"
mv -f "${pointer_tmp}" "${pointer}"
printf 'promoted preflight evidence: %s\n' "${evidence_id}"
