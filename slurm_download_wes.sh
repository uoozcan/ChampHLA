#!/bin/bash
#SBATCH --job-name=pihla_dl_wes
#SBATCH --account=project_2008084
#SBATCH --partition=small
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --array=1-50
#SBATCH --output=/scratch/project_2008084/hla_calibration/logs/dl_wes_%a_%j.out
#SBATCH --error=/scratch/project_2008084/hla_calibration/logs/dl_wes_%a_%j.err

# Downloads the HLA region (chr6:28-34 Mb) from 1000G exome BAMs via samtools streaming.
# Produces small HLA-region BAMs (~50-150 MB each) rather than full exome BAMs.
# Tries primary URL first, then alt-date fallbacks (col 3 of TSV, pipe-separated).
# Also tries both chr6 and 6 chromosome notation since 1000G BAMs vary.

set -euo pipefail

INDEX=/scratch/project_2008084/hla_calibration/wes/index/sample_bam_urls.tsv
OUTDIR=/scratch/project_2008084/hla_calibration/wes/bams
SAMTOOLS=/projappl/project_2008084/bin/samtools

mkdir -p "${OUTDIR}"
mkdir -p /scratch/project_2008084/hla_calibration/logs

module load samtools 2>/dev/null || true

# Prefer project samtools if it can actually run (binary may have missing shared libs on some nodes)
if "${SAMTOOLS}" --version >/dev/null 2>&1; then
  SAM="${SAMTOOLS}"
else
  SAM=samtools
fi

# Read the Nth line
line=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${INDEX}")
if [[ -z "${line}" ]]; then
  echo "No data for array task ${SLURM_ARRAY_TASK_ID}" >&2
  exit 0
fi

sample_id=$(echo "${line}" | cut -f1)
# Primary URL (col 2) and alt-date fallbacks (col 3, pipe-separated)
PRIMARY_URL=$(echo "${line}" | cut -f2)
ALT_URLS=$(echo "${line}" | cut -f3 | tr '|' ' ')

# Derive output filename from primary URL (preserves 1000G naming convention)
bam_filename=$(basename "${PRIMARY_URL}")
out_bam="${OUTDIR}/${bam_filename}"
out_bai="${out_bam}.bai"

echo "[$(date)] Task ${SLURM_ARRAY_TASK_ID}: ${sample_id}"
echo "  Primary URL: ${PRIMARY_URL}"
echo "  Output: ${out_bam}"

if [[ -f "${out_bam}" && -f "${out_bai}" ]]; then
  echo "  Already exists and indexed, skipping"
  exit 0
fi

# Try each URL with both chr notations; accept first combination that yields reads
SUCCESS=0
for URL in ${PRIMARY_URL} ${ALT_URLS}; do
  for REGION in "chr6:28000000-34000000" "6:28000000-34000000"; do
    echo "  Trying: ${URL} region ${REGION}"
    TMP_BAM="${out_bam}.tmp"
    if "${SAM}" view -b -@ "${SLURM_CPUS_PER_TASK}" -o "${TMP_BAM}" "${URL}" "${REGION}" 2>/dev/null; then
      NREADS=$("${SAM}" view -c "${TMP_BAM}" 2>/dev/null || echo 0)
      if [[ "${NREADS}" -gt 0 ]]; then
        mv "${TMP_BAM}" "${out_bam}"
        "${SAM}" index "${out_bam}"
        echo "  [OK] ${NREADS} reads — BAM: $(du -sh "${out_bam}" | cut -f1)"
        SUCCESS=1
        break 2
      else
        echo "  [WARN] 0 reads from ${URL} region ${REGION}"
      fi
    else
      echo "  [WARN] samtools failed for ${URL} region ${REGION}"
    fi
    rm -f "${TMP_BAM}"
  done
done

if [[ "${SUCCESS}" -eq 0 ]]; then
  echo "  [FAIL] No reads found from any URL/region combination for ${sample_id}" >&2
  exit 1
fi

echo "[$(date)] ${sample_id} complete"
