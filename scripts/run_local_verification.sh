#!/bin/bash
set -euo pipefail

ROOT=$(cd "$(dirname "$0")/.." && pwd)
cd "${ROOT}"
export PYTHONPATH="${ROOT}/src"

python3 -m pytest -q
python3 -m champhla_recovery.cli render-manuscript-tables \
  --registry result_registry.tsv \
  --root "${ROOT}" \
  --output manuscripts/shared/generated_results.md
python3 -m champhla_recovery.cli audit-manuscript-claims \
  --manuscript manuscripts/benchmark/manuscript.md \
  --registry result_registry.tsv \
  --claims manuscripts/claim_audit.tsv \
  --root "${ROOT}" \
  --output artifacts/benchmark_manuscript_audit.json || true
python3 -m champhla_recovery.cli audit-manuscript-claims \
  --manuscript manuscripts/method_conditional/manuscript.md \
  --registry result_registry.tsv \
  --claims manuscripts/claim_audit.tsv \
  --root "${ROOT}" \
  --output artifacts/method_manuscript_audit.json || true

