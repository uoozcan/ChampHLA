"""
make_performance_benchmark_samplesheets.py

Generates per-sample samplesheet CSV files for the WGS n=30 computational
performance benchmark. Selects 30 samples stratified by population from the
WGS wave2 manifest, verifies BAM files exist, and writes one samplesheet per sample.

Usage:
    python3.11 bin/make_performance_benchmark_samplesheets.py

Outputs:
    analysis/performance_benchmark_n30/samplesheets/wgs/{SAMPLE}.csv
    analysis/performance_benchmark_n30/wgs_n30_sample_list.txt
"""

import csv
import sys
from collections import defaultdict
from pathlib import Path

MANIFEST = Path(
    "/scratch/project_2008084/pihla-publish/analysis/1000g_realdata"
    "/wgs_wave2_inputs/sequencing_manifest.tsv"
)
BAM_DIR = Path("/scratch/project_2008084/hla_calibration/wgs/bams_new")
OUT_DIR = Path(
    "/scratch/project_2008084/pihla-publish/analysis"
    "/performance_benchmark_n30/samplesheets/wgs"
)
SAMPLE_LIST = Path(
    "/scratch/project_2008084/pihla-publish/analysis"
    "/performance_benchmark_n30/wgs_n30_sample_list.txt"
)
N_TOTAL = 30

# ---------------------------------------------------------------------------
# Load manifest
# ---------------------------------------------------------------------------
samples_by_pop: dict[str, list[str]] = defaultdict(list)

with open(MANIFEST) as fh:
    reader = csv.DictReader(fh, delimiter="\t")
    for row in reader:
        if row.get("available", "0") != "1":
            continue
        sample = row["sample"]
        pop = row["population"]
        bam = BAM_DIR / f"{sample}_hla.bam"
        if not bam.exists():
            print(f"  SKIP {sample}: BAM not found at {bam}", file=sys.stderr)
            continue
        samples_by_pop[pop].append(sample)

print(f"Populations: {dict({k: len(v) for k, v in samples_by_pop.items()})}", file=sys.stderr)

# ---------------------------------------------------------------------------
# Stratified selection: 6 per population (5 pops × 6 = 30)
# ---------------------------------------------------------------------------
pops = sorted(samples_by_pop)
n_per_pop = N_TOTAL // len(pops)
remainder = N_TOTAL % len(pops)

selected: list[str] = []
for i, pop in enumerate(pops):
    n = n_per_pop + (1 if i < remainder else 0)
    chosen = samples_by_pop[pop][:n]
    selected.extend(chosen)
    print(f"  {pop}: selected {len(chosen)}/{len(samples_by_pop[pop])}: {chosen}", file=sys.stderr)

print(f"\nTotal selected: {len(selected)} samples", file=sys.stderr)

# ---------------------------------------------------------------------------
# Write one samplesheet per sample
# ---------------------------------------------------------------------------
OUT_DIR.mkdir(parents=True, exist_ok=True)

for sample in selected:
    bam = BAM_DIR / f"{sample}_hla.bam"
    bai = BAM_DIR / f"{sample}_hla.bam.bai"
    out_csv = OUT_DIR / f"{sample}.csv"
    with open(out_csv, "w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(["sample_id", "bam_path"])
        writer.writerow([sample, str(bam)])
    print(f"  Wrote {out_csv}", file=sys.stderr)

# Write sample list
SAMPLE_LIST.parent.mkdir(parents=True, exist_ok=True)
with open(SAMPLE_LIST, "w") as fh:
    fh.write("\n".join(selected) + "\n")

print(f"\nSample list written to {SAMPLE_LIST}")
print(f"Samplesheets written to {OUT_DIR}/")
print(f"Run slurm_performance_benchmark_wgs_n30.sh to submit the array job.")
