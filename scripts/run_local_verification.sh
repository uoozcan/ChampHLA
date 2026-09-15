#!/bin/bash
set -euo pipefail

ROOT=$(cd "$(dirname "$0")/.." && pwd)
cd "${ROOT}"
export PYTHONPATH="${ROOT}/src"

python3 -c 'import sys; print("WARNING: local Python is below the supported >=3.10 boundary; results are advisory") if sys.version_info < (3, 10) else None'

expect_blocked() {
  set +e
  "$@"
  status=$?
  set -e
  if [ "${status}" -ne 2 ]; then
    echo "Expected fail-closed exit 2, observed ${status}: $*" >&2
    return 1
  fi
}

pytest -q
python3 -m pytest -q
python3 scripts/verify_discovery_reproduction.py

python3 -m champhla_recovery.cli render-manuscript-tables \
  --registry result_registry.tsv \
  --root "${ROOT}" \
  --output manuscripts/shared/generated_results.md

RENDER_CHECK=$(mktemp -d)
trap 'rm -rf -- "${RENDER_CHECK}"' EXIT
python3 scripts/consensus_head_to_head.py \
  --registry result_registry.tsv \
  --comparators configs/comparator_manifest.json \
  --evaluation-dir artifacts/generated/development_evaluation \
  --output "${RENDER_CHECK}/consensus_head_to_head.md" \
  --summary "${RENDER_CHECK}/consensus_head_to_head.json"

expect_blocked python3 -c \
  'from champhla_confirmation.cli import audit_dataset_discovery_registry_main; raise SystemExit(audit_dataset_discovery_registry_main())' \
  --datasets external/dataset_discovery/dataset_registry.tsv \
  --truth external/dataset_discovery/truth_registry.tsv \
  --crosswalk external/dataset_discovery/sample_crosswalk.tsv \
  --pilot external/dataset_discovery/pilot_manifest.tsv \
  --output external/dataset_discovery/audit_summary.json

expect_blocked python3 -m champhla_recovery.cli audit-manuscript-claims \
  --manuscript manuscripts/benchmark/manuscript.md \
  --registry result_registry.tsv \
  --claims manuscripts/claim_audit.tsv \
  --root "${ROOT}" \
  --output artifacts/benchmark_manuscript_audit.json

expect_blocked python3 -m champhla_recovery.cli audit-manuscript-claims \
  --manuscript manuscripts/benchmark/manuscript.md \
  --supplement manuscripts/benchmark/supplementary.md \
  --registry result_registry.tsv \
  --claims manuscripts/claim_audit.tsv \
  --root "${ROOT}" \
  --output artifacts/benchmark_combined_audit.json

python3 -m champhla_recovery.cli audit-manuscript-claims \
  --manuscript manuscripts/method_conditional/manuscript.md \
  --registry result_registry.tsv \
  --claims manuscripts/claim_audit.tsv \
  --root "${ROOT}" \
  --output artifacts/method_manuscript_audit.json

expect_blocked python3 -m champhla_recovery.cli audit-release-readiness \
  --project-root "${ROOT}" \
  --config configs/release_requirements.json \
  --output artifacts/release_readiness.json

expect_blocked python3 -m champhla_recovery.cli validate-figure-manifest \
  --project-root "${ROOT}" \
  --manifest manuscripts/figures/figure_manifest.tsv \
  --output artifacts/figure_audit.json

expect_blocked python3 -m champhla_recovery.cli build-submission-package \
  --project-root "${ROOT}" \
  --main manuscripts/benchmark/manuscript.md \
  --supplement manuscripts/benchmark/supplementary.md \
  --bibliography manuscripts/submission/references.bib \
  --declarations manuscripts/submission/declarations_template.md \
  --output-dir manuscripts/submission \
  --manifest manuscripts/submission/submission_manifest.tsv \
  --parity-audit artifacts/submission_docx_parity.json \
  --bibliography-audit artifacts/submission_bibliography_audit.json \
  --declarations-audit artifacts/submission_declarations_audit.json

echo "Local verification passed; prospective manuscript/release gates remain correctly blocked."
