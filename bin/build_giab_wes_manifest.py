#!/usr/bin/env python3
"""GIAB trio WES sequencing_source.tsv. Clear tables/manifests before re-scoring."""
import os
BASE="/scratch/project_2008084/pihla-publish/analysis/giab_trio_benchmark"
RES="/scratch/project_2008084/hla_calibration/giab_trio/wes/results"
truth={l.split("\t",1)[0] for i,l in enumerate(open(BASE+"/truth_long.tsv")) if i}
rows=[(s,"","wes",os.path.join(RES,s),"1") for s in sorted(truth) if os.path.isdir(os.path.join(RES,s,"optitype"))]
os.makedirs(BASE+"/wes",exist_ok=True)
with open(BASE+"/wes/sequencing_source.tsv","w") as fh:
    fh.write("sample\tpopulation\tmodality\tdata_locator\tavailable\n")
    for r in rows: fh.write("\t".join(r)+"\n")
print(f"[giab-wes-manifest] wrote {len(rows)} samples")
