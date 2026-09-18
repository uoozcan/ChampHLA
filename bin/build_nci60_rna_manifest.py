#!/usr/bin/env python3
"""Build the NCI-60 RNA sequencing_source.tsv from whatever results dirs exist.

A sample is included only if (a) it has 2-field truth in truth_long.tsv and
(b) its pipeline results dir is present AND has calling-tool output (optitype/).
Run after the RNA scale job completes to produce the powered (~43-line) manifest.
"""
import os, sys

BASE = "/scratch/project_2008084/pihla-publish/analysis/nci60_benchmark"
RESULTS = "/scratch/project_2008084/hla_calibration/nci60/rna/results"
TRUTH = os.path.join(BASE, "truth_long.tsv")
OUT = os.path.join(BASE, "sequencing_source.tsv")

# truth sample names (skip header)
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
    has_call = os.path.isdir(os.path.join(rdir, "optitype")) or os.path.isdir(rdir)
    # require the champion tool output so we never benchmark a half-typed line
    if os.path.isdir(os.path.join(rdir, "optitype")):
        rows.append((s, "", "rnaseq", rdir, "1"))
    else:
        skipped.append(s)

with open(OUT, "w") as fh:
    fh.write("sample\tpopulation\tmodality\tdata_locator\tavailable\n")
    for r in rows:
        fh.write("\t".join(r) + "\n")

print(f"[manifest] wrote {len(rows)} samples -> {OUT}")
if skipped:
    print(f"[manifest] {len(skipped)} truth lines NOT yet typed (no optitype/): {', '.join(skipped)}")
