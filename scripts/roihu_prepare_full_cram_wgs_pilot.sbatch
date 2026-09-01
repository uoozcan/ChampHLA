#!/bin/bash
#SBATCH --job-name=cc_fullcram_extract
#SBATCH --account=project_2008084
#SBATCH --partition=small
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --array=1-2%1
#SBATCH --output=/scratch/project_2008084/champhla_cc_confirmation/runs/full_cram_wgs_pilot/logs/extract_%A_%a.out
#SBATCH --error=/scratch/project_2008084/champhla_cc_confirmation/runs/full_cram_wgs_pilot/logs/extract_%A_%a.err

# Truth-blind, mate-aware extraction from the original public WGS CRAM.  The
# interval list contains the primary-assembly MHC plus every HLA-A/B/C contig in
# the exact GRCh38DH reference.  No truth file is read by this script.
set -eo pipefail

ROOT=/scratch/project_2008084/champhla_cc_confirmation
MANIFEST=${ROOT}/configs/full_cram_wgs_pilot.tsv
REFERENCE=/scratch/project_2008084/hla_references/genomes/GRCh38_full_analysis_set_plus_decoy_hla.fa
RUN_ROOT=${ROOT}/runs/full_cram_wgs_pilot
INPUT_ROOT=${RUN_ROOT}/inputs
LOG_ROOT=${RUN_ROOT}/logs

export CSC_ENV_INIT_NON_INTERACTIVE=yes
source /etc/profile.d/zz-csc-env.sh
module load bio-apps/v202603
module load samtools/1.21
set -u
mkdir -p "${INPUT_ROOT}" "${LOG_ROOT}"

line_number=$((SLURM_ARRAY_TASK_ID + 1))
IFS=$'\t' read -r sample_id cram_url expected_md5 source_record < <(sed -n "${line_number}p" "${MANIFEST}")
if [[ -z "${sample_id:-}" || -z "${cram_url:-}" ]]; then
  echo "Missing pilot-manifest row for task ${SLURM_ARRAY_TASK_ID}" >&2
  exit 2
fi

sample_root=${INPUT_ROOT}/${sample_id}
mkdir -p "${sample_root}"
{
  samtools --version
  module list
} > "${sample_root}/extraction_environment.txt" 2>&1
regions=${sample_root}/grch38dh_hla_regions.bed
tmp_bam=${sample_root}/${sample_id}.mate_aware_hla.unsorted.bam
final_bam=${sample_root}/${sample_id}.mate_aware_hla.bam

# The first interval captures classical chromosome-6 mappings. Every A/B/C HLA
# alternate contig is included in full; D-class and pseudogene contigs are not.
awk -F '\t' 'BEGIN {OFS="\t"}
  $1 == "chr6" {print $1,27999999,34000000}
  $1 ~ /^HLA-[ABC]\*/ {print $1,0,$2}
' "${REFERENCE}.fai" > "${regions}"

if [[ $(wc -l < "${regions}") -lt 4 ]]; then
  echo "Reference interval construction failed: ${regions}" >&2
  exit 3
fi

if [[ ! -s "${tmp_bam}" ]]; then
  samtools view \
    --threads "${SLURM_CPUS_PER_TASK}" \
    --reference "${REFERENCE}" \
    --fetch-pairs \
    --target-file "${regions}" \
    --bam \
    --output "${tmp_bam}" \
    "${cram_url}"
fi

samtools quickcheck -v "${tmp_bam}"
samtools sort --threads "${SLURM_CPUS_PER_TASK}" -o "${final_bam}" "${tmp_bam}"
samtools index --threads "${SLURM_CPUS_PER_TASK}" "${final_bam}"
samtools quickcheck -v "${final_bam}"
samtools flagstat --threads "${SLURM_CPUS_PER_TASK}" "${final_bam}" > "${sample_root}/flagstat.txt"
samtools idxstats "${final_bam}" > "${sample_root}/idxstats.tsv"
sha256sum "${final_bam}" "${final_bam}.bai" "${regions}" > "${sample_root}/input_sha256.txt"

{
  printf 'sample_id\tcram_url\tdeclared_cram_md5\tsource_record\treference\treference_sha256\textraction_policy\n'
  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "${sample_id}" "${cram_url}" "${expected_md5}" "${source_record}" \
    "${REFERENCE}" "$(sha256sum "${REFERENCE}" | cut -d ' ' -f 1)" \
    'chr6:28000000-34000000 plus complete HLA-A/B/C contigs; samtools --fetch-pairs'
} > "${sample_root}/source_provenance.tsv"

# The unsorted stream is retained until all validation succeeds, then removed to
# preserve capacity. The validated coordinate-sorted BAM remains immutable input.
rm -f "${tmp_bam}"
touch "${sample_root}/EXTRACTION_COMPLETE"
