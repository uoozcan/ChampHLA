#!/bin/bash
set -euo pipefail

CODE_ROOT=${CHAMPHLA_CODE_ROOT:-/scratch/project_2008084/champhla_plurality}
source "${CHAMPHLA_PATH_CONFIG:-${CODE_ROOT}/configs/roihu_paths.env}"
export PYTHONPATH=${CHAMPHLA_CODE_ROOT}/src
output_root=${CHAMPHLA_RUN_ROOT}/manifests/same_resource
mkdir -p "${output_root}"

# Only missing public CRAI files are streamed to calculate their checksums.
# No CRAM, BAM, or FASTQ data are downloaded by this metadata-only command.
python3 -c 'from champhla_confirmation.cli import resolve_same_resource_index_checksums_main; raise SystemExit(resolve_same_resource_index_checksums_main())' \
  --roster "${CHAMPHLA_CODE_ROOT}/cohorts/same_resource_truth_free_roster.tsv" \
  --assay-manifest "${CHAMPHLA_CODE_ROOT}/cohorts/official_1000g_assay_manifest.tsv" \
  --output "${output_root}/public_index_checksums.tsv"

python3 -c 'from champhla_confirmation.cli import build_same_resource_run_manifest_main; raise SystemExit(build_same_resource_run_manifest_main())' \
  --roster "${CHAMPHLA_CODE_ROOT}/cohorts/same_resource_truth_free_roster.tsv" \
  --assay-manifest "${CHAMPHLA_CODE_ROOT}/cohorts/official_1000g_assay_manifest.tsv" \
  --ena-report "${CHAMPHLA_CODE_ROOT}/cohorts/sources/geuvadis_ena_fastq_report.tsv" \
  --index-checksums "${output_root}/public_index_checksums.tsv" \
  --output "${output_root}/same_resource_truth_free.tsv"

python3 -c 'from champhla_confirmation.cli import freeze_run_manifest_main; raise SystemExit(freeze_run_manifest_main())' \
  --manifest "${output_root}/same_resource_truth_free.tsv" \
  --expected-counts "${CHAMPHLA_CODE_ROOT}/configs/roihu_storage_targets.json" \
  --output "${output_root}/same_resource_truth_free.freeze.json"
