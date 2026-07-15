#!/bin/bash
#SBATCH --job-name=dl_rna_hla
#SBATCH --account=project_2008084
#SBATCH --partition=small
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --array=1-30%10
#SBATCH --output=/scratch/project_2008084/hla_calibration/logs/dl_rna_hla_%A_%a.out
#SBATCH --error=/scratch/project_2008084/hla_calibration/logs/dl_rna_hla_%A_%a.err

# Stream-extract HLA region (chr6:28-34 Mb) from remote Geuvadis RNA-seq BAMs.
# Downloads only ~50-200 MB per sample instead of 5-10 GB full BAMs.

set -euo pipefail

MANIFEST=/scratch/project_2008084/pihla-publish/analysis/performance_benchmark_n30/rna_n30_url_manifest.tsv
OUTDIR=/scratch/project_2008084/pihla-publish/analysis/performance_benchmark_n30/rna_hla_bams
HLA_REGION="chr6:28000000-34000000"

module load samtools 2>/dev/null || true

mkdir -p "${OUTDIR}" /scratch/project_2008084/hla_calibration/logs

line=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${MANIFEST}")
sample_id=$(echo "${line}" | cut -f1)
bam_url=$(echo "${line}" | cut -f2)

[[ -z "${sample_id}" ]] && { echo "[ERROR] Empty sample for task ${SLURM_ARRAY_TASK_ID}"; exit 1; }
[[ -z "${bam_url}" ]] && { echo "[ERROR] No URL for ${sample_id}"; exit 1; }

outbam="${OUTDIR}/${sample_id}_rna_hla.bam"

if [[ -f "${outbam}" && -f "${outbam}.bai" ]]; then
    reads=$(samtools view -c "${outbam}" 2>/dev/null || echo 0)
    if [[ "${reads}" -gt 1000 ]]; then
        echo "[SKIP] ${sample_id} already downloaded (${reads} reads)"
        exit 0
    fi
fi

echo "[START] ${sample_id} at $(date)"
echo "[URL] ${bam_url}"

samtools view -b -h "${bam_url}" ${HLA_REGION} \
    | samtools sort -@ 1 -o "${outbam}"

samtools index "${outbam}"

reads=$(samtools view -c "${outbam}")
echo "[DONE] ${sample_id} — ${reads} reads extracted, $(du -h "${outbam}" | cut -f1)"

if [[ "${reads}" -lt 100 ]]; then
    echo "[WARN] Very few reads (${reads}) — may indicate wrong region name or empty BAM"
fi
