#!/bin/bash
set -eo pipefail

ROOT=${CHAMPHLA_CODE_ROOT:-/scratch/project_2008084/champhla_plurality}
PATH_CONFIG=${CHAMPHLA_PATH_CONFIG:-${ROOT}/configs/roihu_paths.env}

export CSC_ENV_INIT_NON_INTERACTIVE=yes
source /etc/profile.d/zz-csc-env.sh
set -u
module load bio-apps/v202603
module load samtools/1.21
module load nextflow/25.10.2-standalone

test -s "${PATH_CONFIG}"
source "${PATH_CONFIG}"
test -d "${CHAMPHLA_CODE_ROOT}"
test -d "${CHAMPHLA_RUN_ROOT}"
test -d "${CHAMPHLA_INPUT_ROOT}"
test -d "${CHAMPHLA_TRUTH_ROOT}"
test -s "${CHAMPHLA_REFERENCE_ROOT}/genomes/GRCh38_full_analysis_set_plus_decoy_hla.fa"
test -s "${CHAMPHLA_REFERENCE_ROOT}/genomes/GRCh38_full_analysis_set_plus_decoy_hla.fa.fai"
test -s "${CHAMPHLA_REFERENCE_ROOT}/t1k_hlaidx/_dna_seq.fa"
test -s "${CHAMPHLA_REFERENCE_ROOT}/t1k_hlaidx/_rna_seq.fa"
test -s "${CHAMPHLA_CODE_ROOT}/configs/roihu_workflow_lock.json"
test -s "${CHAMPHLA_PIPELINE_ROOT}/main.nf"
test -s "${CHAMPHLA_PIPELINE_ROOT}/nextflow.config"
apptainer exec "${CHAMPHLA_CONTAINER_ROOT}/hlahd.sif" \
  test -s /app/hlahd.1.4.0/HLA_gene.split.txt
apptainer exec "${CHAMPHLA_CONTAINER_ROOT}/hlahd.sif" \
  test -d /app/hlahd.1.4.0/dictionary
apptainer exec "${CHAMPHLA_CONTAINER_ROOT}/hlahd.sif" \
  test -d /app/hlahd.1.4.0/freq_data

# Reconstruct both permitted POLYSOLVER wrapper variants inside the pinned container.
# This validates the source hash, exact three-token chr6 correction, independent TMP
# substitution, and the derived-wrapper audit without writing into the Git checkout.
patch_tmp=$(mktemp -d "${CHAMPHLA_RUN_ROOT}/polysolver-preflight.XXXXXX")
cleanup_patch_tmp() {
  rm -rf "${patch_tmp}"
}
trap cleanup_patch_tmp EXIT
for contig in 6 chr6; do
  output=${patch_tmp}/wrapper.${contig}
  audit=${patch_tmp}/audit.${contig}.json
  apptainer exec \
    --bind "${CHAMPHLA_CODE_ROOT}:${CHAMPHLA_CODE_ROOT}:ro" \
    --bind "${patch_tmp}:${patch_tmp}" \
    "${CHAMPHLA_CONTAINER_ROOT}/polysolver.sif" \
    python3 "${CHAMPHLA_CODE_ROOT}/workflow/bin/patch_polysolver_wrapper.py" \
      --source /home/polysolver/scripts/shell_call_hla_type \
      --output "${output}" --audit "${audit}" \
      --spec "${CHAMPHLA_CODE_ROOT}/workflow/conf/polysolver_wrapper_patch.json" \
      --contig "${contig}" --tmp-dir /frozen/picard_tmp
  apptainer exec \
    --bind "${CHAMPHLA_CODE_ROOT}:${CHAMPHLA_CODE_ROOT}:ro" \
    --bind "${patch_tmp}:${patch_tmp}" \
    "${CHAMPHLA_CONTAINER_ROOT}/polysolver.sif" \
    python3 "${CHAMPHLA_CODE_ROOT}/workflow/bin/patch_polysolver_wrapper.py" --verify-only \
      --source /home/polysolver/scripts/shell_call_hla_type \
      --output "${output}" --audit "${audit}" \
      --spec "${CHAMPHLA_CODE_ROOT}/workflow/conf/polysolver_wrapper_patch.json" \
      --contig "${contig}" --tmp-dir /frozen/picard_tmp
done
python3 --version
samtools --version | head -n 1
nextflow -version
apptainer --version
df -h /scratch
du -sh "${CHAMPHLA_CODE_ROOT}" "${CHAMPHLA_RUN_ROOT}" "${CHAMPHLA_INPUT_ROOT}" \
  "${CHAMPHLA_TRUTH_ROOT}" "${CHAMPHLA_REFERENCE_ROOT}"
module list
sinfo --format='%P %a %l %D %c %m %G'

export PYTHONPATH="${CHAMPHLA_CODE_ROOT}/src"
python3 -c 'from champhla_confirmation.cli import validate_workflow_lock_main; raise SystemExit(validate_workflow_lock_main())' \
  --project-root "${CHAMPHLA_CODE_ROOT}" \
  --lock "${CHAMPHLA_CODE_ROOT}/configs/roihu_workflow_lock.json" \
  --output "${CHAMPHLA_RUN_ROOT}/workflow_lock_audit.json"
python3 -c 'from champhla_confirmation.manifests import validate_caller_reference_attestation; import sys; failures=validate_caller_reference_attestation(sys.argv[1], True); assert not failures, failures' \
  "${CHAMPHLA_CODE_ROOT}/configs/roihu_caller_reference_attestation.json"
python3 -c 'from champhla_confirmation.cli import audit_roihu_environment_main; raise SystemExit(audit_roihu_environment_main())' \
  --site-config "${CHAMPHLA_CODE_ROOT}/configs/roihu_site.json" \
  --output "${CHAMPHLA_RUN_ROOT}/environment_inventory.json"

sha256sum "${CHAMPHLA_CODE_ROOT}/configs/confirmation_protocol.json" \
  "${CHAMPHLA_CODE_ROOT}/configs/runtime_versions.json" \
  "${CHAMPHLA_CODE_ROOT}/configs/roihu_site.json" \
  "${CHAMPHLA_CODE_ROOT}/configs/roihu_workflow_lock.json"
