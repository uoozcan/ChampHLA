#!/bin/bash
# submit_path_b_pipeline.sh — Path B: 100 trimodal samples, 10-sample batches throughout.
#
# All WGS, WES, and RNA batches are strictly 10 samples to control disk usage.
# Each batch downloads, types, and cleans before the next batch starts.
#
# Usage: bash submit_path_b_pipeline.sh

set -euo pipefail
REPO=/scratch/project_2008084/pihla-publish
cd "$REPO"

echo "[$(date)] Path B: 10-sample batches for 100 trimodal samples"
echo ""

# ── Sample lists ──────────────────────────────────────────────────────────

# 32 new GEUVADIS+truth samples (WGS + WES + RNA all needed)
NEW_32=(
    NA06984 NA06989 NA11829 NA11831 NA11832 NA11843
    HG00183 HG00185 HG00186 HG00187 HG00188 HG00189 HG00266
    HG00110 HG00111 HG00112 HG00114 HG00115 HG00116
    NA20503 NA20507 NA20513 NA20514 NA20515 NA20516
    NA18486 NA18488 NA18489 NA18498 NA18499 NA18510 NA18511
)

# Existing cohort missing WES (49 samples)
WES_MISSING=(
    HG00096 HG00097 HG00099 HG00100 HG00101 HG00102 HG00103 HG00104 HG00105 HG00106
    HG00107 HG00108 HG00109 NA10847 NA10851 NA11830 NA11840 NA11920 NA12004 NA12044
    NA12154 NA18561 NA18562 NA18563 NA18564 NA19093 NA19099 NA19116 NA19119 NA19129
    NA19130 NA19131 NA19152 NA19153 NA19209 NA19210 NA19238 NA19239 NA20502 NA20504
    NA20505 NA20506 NA20508 NA20509 NA20510 NA20512 NA20519 NA20521 NA20528
)

# Existing cohort missing RNA (19 samples with GEUVADIS URLs)
RNA_MISSING=(
    NA07037 NA07048 NA07051 NA07056 NA07347 NA07357 NA12044 NA12154 NA12751
    NA18517 NA18522 NA18523 NA18526 NA18530 NA18532 NA18536 NA18537 NA18542 NA18543
)

# ── Helper: submit batches of exactly 10 ─────────────────────────────────
# Prints progress to stderr; echoes only the final job ID to stdout.
submit_batches() {
    local script="$1"; local dep="$2"; shift 2
    local samples=("$@")
    local prev="$dep"
    local batch_num=1
    for (( i=0; i<${#samples[@]}; i+=10 )); do
        local batch=("${samples[@]:$i:10}")
        if [[ -z "$prev" ]]; then
            JOB=$(sbatch --parsable "$script" "${batch[@]}")
        else
            JOB=$(sbatch --parsable --dependency=afterok:${prev} "$script" "${batch[@]}")
        fi
        prev="$JOB"
        echo "  batch $((batch_num++)): job $JOB [${#batch[@]} samples]" >&2
    done
    echo "$prev"
}

# ── WGS batches: 32 new samples, 10 per batch ────────────────────────────
echo "WGS batches (32 new samples, batch3 index):"
LAST_WGS=$(submit_batches "$REPO/slurm_wgs_batch10.sh" "" "${NEW_32[@]}")
echo "  → last job: $LAST_WGS"

# ── WES batches: 49 existing missing, 10 per batch (batch2 index) ─────────
echo ""
echo "WES batches (49 existing missing, batch2 index):"
LAST_WES_EX=$(submit_batches "$REPO/slurm_wes_batch10.sh" "" "${WES_MISSING[@]}")
echo "  → last job: $LAST_WES_EX"

# ── WES batches: 32 new samples, 10 per batch (batch3 index) ─────────────
echo ""
echo "WES batches (32 new samples, batch3 index):"
TMP_WES=$(mktemp /tmp/wes_b3_XXXXXX.sh)
sed 's|sample_bam_urls_batch2\.tsv|sample_bam_urls_batch3.tsv|' \
    "$REPO/slurm_wes_batch10.sh" > "$TMP_WES"
chmod +x "$TMP_WES"
LAST_WES_NEW=$(submit_batches "$TMP_WES" "" "${NEW_32[@]}")
echo "  → last job: $LAST_WES_NEW"

# ── RNA batches: 19 existing GEUVADIS, 10 per batch (batch2 index) ────────
echo ""
echo "RNA batches (19 existing GEUVADIS, batch2 index):"
LAST_RNA_EX=$(submit_batches "$REPO/slurm_rna_batch10.sh" "" "${RNA_MISSING[@]}")
echo "  → last job: $LAST_RNA_EX"

# ── RNA batches: 32 new samples, 10 per batch (batch3 index) ─────────────
echo ""
echo "RNA batches (32 new samples, batch3 index):"
TMP_RNA=$(mktemp /tmp/rna_b3_XXXXXX.sh)
sed 's|sample_fastq_urls_batch2\.tsv|sample_fastq_urls_batch3.tsv|' \
    "$REPO/slurm_rna_batch10.sh" > "$TMP_RNA"
chmod +x "$TMP_RNA"
LAST_RNA_NEW=$(submit_batches "$TMP_RNA" "" "${NEW_32[@]}")
echo "  → last job: $LAST_RNA_NEW"

# ── Benchmark: after ALL typing complete ─────────────────────────────────
echo ""
ALL_DEPS="${LAST_WGS}:${LAST_WES_EX}:${LAST_WES_NEW}:${LAST_RNA_EX}:${LAST_RNA_NEW}"
BENCH=$(sbatch --parsable --dependency=afterok:${ALL_DEPS} \
    "$REPO/slurm_benchmark_trimodal.sh")
echo "Benchmark:     job $BENCH"

FIG1=$(sbatch --parsable --dependency=afterok:${BENCH} "$REPO/slurm_generate_figures.sh")
FIG2=$(sbatch --parsable --dependency=afterok:${BENCH} "$REPO/slurm_generate_new_figures.sh")
echo "Figures (v2):  job $FIG1"
echo "Figures (new): job $FIG2"

echo ""
echo "[$(date)] All jobs submitted."
echo "Monitor: squeue -u ozcanumu -o '%.10i %.22j %.8T %.10M %R'"
