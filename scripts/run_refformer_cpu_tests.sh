#!/bin/bash
set -euo pipefail

ROOT=$(cd "$(dirname "$0")/.." && pwd)
cd "${ROOT}/ablations/refformer"
python3 -c 'import torch; assert torch.cuda.is_available() is False or isinstance(torch.__version__, str)'
PYTHONPATH=src python3 -m pytest -q

