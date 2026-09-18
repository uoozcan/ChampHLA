#!/bin/bash
# Generates venex_batch1_recovery.txt and venex_batch1_2_recovery.txt
# containing samples missing any of hlahd, optitype, or t1k tool results.
set -euo pipefail

REPO=/scratch/project_2008084/pihla-publish
RESULTS=/scratch/project_2008084/hla_calibration/venex/results/wgs

for BATCH in batch1 batch1_2; do
    LIST="${REPO}/venex_${BATCH}_samples.txt"
    OUT="${REPO}/venex_${BATCH}_recovery.txt"
    > "$OUT"
    while IFS= read -r sample; do
        sdir="${RESULTS}/${sample}/results/${sample}"
        [[ -d "${sdir}/hlahd" ]]    && has_hlahd=1    || has_hlahd=0
        [[ -d "${sdir}/optitype" ]] && has_optitype=1 || has_optitype=0
        [[ -d "${sdir}/t1k" ]]      && has_t1k=1      || has_t1k=0
        if [[ "$has_hlahd" -eq 0 || "$has_optitype" -eq 0 || "$has_t1k" -eq 0 ]]; then
            echo "$sample" >> "$OUT"
        fi
    done < "$LIST"
    echo "$(wc -l < "$OUT") samples need recovery in $BATCH → $OUT"
done
