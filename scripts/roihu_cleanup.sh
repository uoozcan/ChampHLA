#!/bin/bash
set -euo pipefail

plan=${1:?usage: roihu_cleanup.sh CLEANUP_PLAN [--execute]}
mode=${2:---dry-run}
python3 -c 'import json,sys; p=json.load(open(sys.argv[1])); assert p["schema_version"] == "champhla-roihu-cleanup-plan-1"; assert p["executable"]; assert all(x["path"] in {"work","tmp"} for x in p["eligible_paths"])' "${plan}"
run_root=$(python3 -c 'import json,sys; print(json.load(open(sys.argv[1]))["run_root"])' "${plan}")
[[ "${run_root}" == /scratch/project_2008084/champhla_plurality_runs/* ]]
if [[ "${mode}" != --execute ]]; then
  echo "Dry run only; eligible paths:"
  python3 -c 'import json,sys; [print(x["absolute_path"]) for x in json.load(open(sys.argv[1]))["eligible_paths"]]' "${plan}"
  exit 0
fi
while IFS= read -r target; do
  [[ "${target}" == "${run_root}/work" || "${target}" == "${run_root}/tmp" ]] || {
    echo "Refusing unsafe cleanup target: ${target}" >&2
    exit 3
  }
  rm -rf -- "${target}"
done < <(python3 -c 'import json,sys; [print(x["absolute_path"]) for x in json.load(open(sys.argv[1]))["eligible_paths"]]' "${plan}")
