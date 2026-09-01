#!/bin/bash
#SBATCH --job-name=cc_kourami_repair
#SBATCH --account=project_2008084
#SBATCH --partition=small
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=180G
#SBATCH --output=/scratch/project_2008084/champhla_cc_confirmation/runs/full_cram_wgs_pilot/logs/kourami_repair_%j.out
#SBATCH --error=/scratch/project_2008084/champhla_cc_confirmation/runs/full_cram_wgs_pilot/logs/kourami_repair_%j.err

# Fail-closed repair for the production Kourami wrapper, which can mask a
# killed Java process with `|| true`. This script consumes the complete,
# immutable mate-aware HLA BAM and writes only inside the isolated pilot.
set -eo pipefail

ROOT=/scratch/project_2008084/champhla_cc_confirmation
RUN_ROOT=${ROOT}/runs/full_cram_wgs_pilot
CONTAINER=/scratch/project_2008084/hla_references/singularity_cache/containers/kourami.sif
PARSER=/scratch/project_2008084/pihla-publish/bin/parse_kourami_results.py
sample_id=${1:?usage: sbatch roihu_repair_kourami_full_bam.sbatch SAMPLE_ID}
input_bam=${RUN_ROOT}/inputs/${sample_id}/${sample_id}.mate_aware_hla.bam
repair_root=${RUN_ROOT}/kourami_repair/${sample_id}
work=${repair_root}/attempt_${SLURM_JOB_ID}
accepted=${repair_root}/accepted

export CSC_ENV_INIT_NON_INTERACTIVE=yes
source /etc/profile.d/zz-csc-env.sh
module load bio-apps/v202603
module load samtools/1.21
set -u

if [[ ! -s "${input_bam}" || ! -s "${input_bam}.bai" || ! -s "${CONTAINER}" ]]; then
  echo "Validated input BAM, index, or pinned Kourami container is absent" >&2
  exit 2
fi
samtools quickcheck -v "${input_bam}"
mkdir -p "${work}"
cd "${work}"

{
  samtools --version
  apptainer --version
  module list
} > environment.txt 2>&1
sha256sum "${input_bam}" "${input_bam}.bai" "${CONTAINER}" "${PARSER}" \
  "${ROOT}/scripts/roihu_repair_kourami_full_bam.sbatch" > pre_run_sha256.txt

r1=${sample_id}.all_hla_1.fq.gz
r2=${sample_id}.all_hla_2.fq.gz
singletons=${sample_id}.all_hla_singletons.fq.gz
panel_bam=${sample_id}.panel.bam
prefix=${sample_id}.kourami
summary=${sample_id}_kourami.txt

# Use every record in the mate-aware HLA BAM, including HLA alternate-contig
# mappings. Name collation preserves pairs without modifying the input.
if [[ ! -s "${r1}" || ! -s "${r2}" ]]; then
  samtools collate -O -u -@ "${SLURM_CPUS_PER_TASK}" "${input_bam}" | \
    samtools fastq -@ "${SLURM_CPUS_PER_TASK}" -n \
      -1 "${r1}" -2 "${r2}" -0 /dev/null -s "${singletons}" -
fi

# Stream SAM into BAM to avoid the large text intermediate that was killed in
# the production wrapper. Two BWA threads bound its per-task memory footprint.
if [[ ! -s "${panel_bam}" ]]; then
  apptainer exec "${CONTAINER}" bwa mem -t 2 \
    /usr/local/bin/kourami-0.9.6/db/All_FINAL_with_Decoy.fa.gz \
    "${r1}" "${r2}" | samtools view -b -o "${panel_bam}" -
fi
samtools quickcheck -v "${panel_bam}"

apptainer exec "${CONTAINER}" bash -lc \
  'jar=/usr/local/bin/kourami-0.9.6/build/Kourami.jar; [[ -f "$jar" ]] || jar=/usr/local/bin/kourami-0.9.6/target/Kourami.jar; [[ -f "$jar" ]] || jar=/usr/local/bin/Kourami.jar; exec java -Xmx24g -jar "$jar" -d /usr/local/bin/kourami-0.9.6/db "$1" -o "$2"' \
  bash "${panel_bam}" "${prefix}"

if [[ ! -s "${prefix}.result" ]]; then
  echo "Kourami produced no nonempty result" >&2
  exit 3
fi
python3 "${PARSER}" --input "${prefix}.result" --sample "${sample_id}" --output "${summary}"
if ! grep -Eq '(^|[^[:alnum:]])(HLA-)?[ABC]\*[0-9A-Za-z]+:[0-9A-Za-z]+' "${summary}"; then
  echo "Kourami repair contains no parseable class-I two-field allele" >&2
  exit 4
fi

sha256sum "${prefix}.result" "${summary}" "${panel_bam}" > output_sha256.txt
mkdir -p "${accepted}"
cp -p "${prefix}.result" "${summary}" environment.txt pre_run_sha256.txt \
  output_sha256.txt "${accepted}/"
sha256sum "${accepted}/${prefix}.result" "${accepted}/${summary}" \
  > "${accepted}/accepted_sha256.txt"
touch "${accepted}/REPAIR_COMPLETE"
