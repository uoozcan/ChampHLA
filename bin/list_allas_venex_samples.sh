#!/bin/bash
# List all sample IDs available in an Allas bucket for VENEX
# Usage: bash list_allas_venex_samples.sh allas:bucket-name [output_sample_list.txt]

BUCKET="${1:-allas:venex-samples}"
OUTPUT="${2:-venex_samples.txt}"
RCLONE=/appl/opt/csc-cli-utils/allas-cli-utils/rclone

module load allas 2>/dev/null || true

echo "Listing samples in $BUCKET..."
"$RCLONE" ls "$BUCKET" 2>/dev/null | \
    awk '{print $2}' | \
    grep -E '\.(cram|bam|fastq\.gz|fq\.gz)$' | \
    sed 's/\.\(cram\|bam\|fastq\.gz\|fq\.gz\|_R[12]\.fastq\.gz\)$//' | \
    sed 's/_R[12]$//' | \
    sort -u > "$OUTPUT"

echo "Found $(wc -l < "$OUTPUT") unique samples → $OUTPUT"
head -5 "$OUTPUT"
