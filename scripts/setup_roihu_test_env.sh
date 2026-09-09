#!/bin/bash
set -eo pipefail

ROOT=${CHAMPHLA_CODE_ROOT:-/scratch/project_2008084/champhla_plurality}
PATH_CONFIG=${CHAMPHLA_PATH_CONFIG:-${ROOT}/configs/roihu_paths.env}
EXPECTED_PYTHON=3.11.15

export CSC_ENV_INIT_NON_INTERACTIVE=yes
source /etc/profile.d/zz-csc-env.sh
set -u
module load bio-apps/v202603
module load samtools/1.21
module load nextflow/25.10.2-standalone
test -s "${PATH_CONFIG}"
source "${PATH_CONFIG}"
ROOT=${CHAMPHLA_CODE_ROOT}

actual_python="$(python3 -c 'import platform; print(platform.python_version())')"
if [[ "${actual_python}" != "${EXPECTED_PYTHON}" ]]; then
  echo "Expected Python ${EXPECTED_PYTHON}; found ${actual_python}" >&2
  exit 2
fi

if [[ -x "${ROOT}/.venv/bin/python" ]]; then
  venv_python="$(${ROOT}/.venv/bin/python -c 'import platform; print(platform.python_version())')"
  if [[ "${venv_python}" != "${EXPECTED_PYTHON}" ]]; then
    echo "Existing .venv uses Python ${venv_python}; preserve or relocate it before retrying" >&2
    exit 3
  fi
else
  python3 -m venv "${ROOT}/.venv"
fi

export PIP_CACHE_DIR="${ROOT}/.cache/pip"
"${ROOT}/.venv/bin/python" -m pip install --disable-pip-version-check -e "${ROOT}[test]"
"${ROOT}/.venv/bin/pytest" -q "${ROOT}/tests"
"${ROOT}/.venv/bin/python" -m pytest -q "${ROOT}/tests"
