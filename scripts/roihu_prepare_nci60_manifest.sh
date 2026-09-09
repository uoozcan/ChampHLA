#!/bin/bash
set -euo pipefail

CODE_ROOT=${CHAMPHLA_CODE_ROOT:-/scratch/project_2008084/champhla_plurality}
source "${CHAMPHLA_PATH_CONFIG:-${CODE_ROOT}/configs/roihu_paths.env}"
export PYTHONPATH=${CHAMPHLA_CODE_ROOT}/src
output_root=${CHAMPHLA_RUN_ROOT}/manifests/nci60
mkdir -p "${output_root}"

# Fetch metadata only (not reads) from the official ENA Portal API. The MD5s
# become part of the truth-free run manifest and are checked after staging.
curl --fail --location --retry 4 --get \
  --data-urlencode 'accession=PRJNA433861' \
  --data-urlencode 'result=read_run' \
  --data-urlencode 'fields=run_accession,fastq_ftp,fastq_md5' \
  --data-urlencode 'format=tsv' \
  --output "${output_root}/ena_file_report.tsv" \
  'https://www.ebi.ac.uk/ena/portal/api/filereport'

python3 -c 'from champhla_confirmation.cli import build_nci60_run_manifest_main; raise SystemExit(build_nci60_run_manifest_main())' \
  --pilot "${CHAMPHLA_CODE_ROOT}/external/dataset_discovery/pilot_manifest.tsv" \
  --ena-report "${output_root}/ena_file_report.tsv" \
  --output "${output_root}/nci60_rna_truth_free.tsv"
python3 -c 'from champhla_confirmation.cli import freeze_run_manifest_main; raise SystemExit(freeze_run_manifest_main())' \
  --manifest "${output_root}/nci60_rna_truth_free.tsv" \
  --expected-counts "${CHAMPHLA_CODE_ROOT}/configs/nci60_expected_counts.json" \
  --output "${output_root}/nci60_rna_truth_free.freeze.json"
