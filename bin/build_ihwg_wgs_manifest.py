#!/usr/bin/env python3
"""Build the IHWG WGS sequencing_source.tsv from typed results. Run after the pilot; then score with:
  python3 bin/run_1000g_benchmark.py --config conf/benchmark_ihwg_wgs.yaml --output-dir analysis/ihwg_benchmark/wgs
NOTE: clear analysis/ihwg_benchmark/wgs/{tables,manifests} before re-scoring a grown cohort (the runner
reuses a cached cohort_manifest — see the NCI-60 gotcha).
"""
import os
BASE = "/scratch/project_2008084/pihla-publish/analysis/ihwg_benchmark"
RESULTS = "/scratch/project_2008084/hla_calibration/ihwg/wgs/results"
TRUTH = os.path.join(BASE, "truth_long.tsv")
OUT = os.path.join(BASE, "wgs", "sequencing_source.tsv")

truth = set()
with open(TRUTH) as fh:
    next(fh)
    for line in fh:
        s = line.split("\t", 1)[0].strip()
        if s:
            truth.add(s)

rows = []
for s in sorted(truth):
    if os.path.isdir(os.path.join(RESULTS, s, "optitype")):
        rows.append((s, "", "wgs", os.path.join(RESULTS, s), "1"))

os.makedirs(os.path.dirname(OUT), exist_ok=True)
with open(OUT, "w") as fh:
    fh.write("sample\tpopulation\tmodality\tdata_locator\tavailable\n")
    for r in rows:
        fh.write("\t".join(r) + "\n")
print(f"[ihwg-wgs-manifest] wrote {len(rows)} typed samples -> {OUT}")
