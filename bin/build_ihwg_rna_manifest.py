#!/usr/bin/env python3
"""Build the IHWG RNA sequencing_source.tsv from whichever lines have been typed.

A sample is included only if (a) it has IHIW truth and (b) its pipeline results dir has the champion
tool output (optitype/). Run after the IHWG RNA pilot/scale completes, then score with:
  python3 bin/run_1000g_benchmark.py --config conf/benchmark_ihwg_rna.yaml --output-dir analysis/ihwg_benchmark/rna
"""
import os

BASE = "/scratch/project_2008084/pihla-publish/analysis/ihwg_benchmark"
RESULTS = "/scratch/project_2008084/hla_calibration/ihwg/rna/results"
TRUTH = os.path.join(BASE, "truth_long.tsv")
OUT = os.path.join(BASE, "rna", "sequencing_source.tsv")

truth_samples = set()
with open(TRUTH) as fh:
    next(fh)
    for line in fh:
        s = line.split("\t", 1)[0].strip()
        if s:
            truth_samples.add(s)

rows, skipped = [], []
for s in sorted(truth_samples):
    rdir = os.path.join(RESULTS, s)
    if os.path.isdir(os.path.join(rdir, "optitype")):
        rows.append((s, "", "rnaseq", rdir, "1"))
    elif os.path.isdir(rdir):
        skipped.append(s)

os.makedirs(os.path.dirname(OUT), exist_ok=True)
with open(OUT, "w") as fh:
    fh.write("sample\tpopulation\tmodality\tdata_locator\tavailable\n")
    for r in rows:
        fh.write("\t".join(r) + "\n")

print(f"[ihwg-manifest] wrote {len(rows)} typed samples -> {OUT}")
if skipped:
    print(f"[ihwg-manifest] {len(skipped)} have a results dir but no optitype/ yet: {', '.join(skipped)}")
