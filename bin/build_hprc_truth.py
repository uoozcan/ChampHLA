#!/usr/bin/env python3
"""Build GIAB/HPRC HLA truth (analysis/hprc_benchmark/truth_long.tsv).

Truth sources (both experimentally grounded, not in-silico short-read calls):
  - HG002: clinical sequence-based typing (SBT) GOLD, Chin et al. 2020 Nat Commun (PMC7508831,
    Suppl. Table 4) — hard-coded below. Doubles as the Immuannot validation control.
  - HG003/HG004/HG005/HG006/HG007: derived from long-read phased assemblies via Immuannot
    (YingZhou001/Immuannot, IPD-IMGT/HLA). Parsed from Immuannot GTF outputs when present under
    --immuannot-dir; samples without an Immuannot file are simply skipped (report realised n).

Schema: sample<TAB>gene<TAB>allele1<TAB>allele2, two-field, HLA-A/-B/-C.
"""
import argparse
import glob
import os
import re

# HG002 clinical SBT gold (two-field), Chin 2020 PMC7508831 Suppl. Table 4
HG002_GOLD = {
    "A": ("01:01", "26:01"),
    "B": ("35:08", "38:01"),
    "C": ("04:01", "12:03"),
}

LOCI = ("A", "B", "C")
ALLELE_RE = re.compile(r"([ABC])\*(\d+:\d+)")  # capture gene + two-field from e.g. A*01:01:01


def two_field(allele):
    parts = allele.split(":")
    return ":".join(parts[:2]) if len(parts) >= 2 else allele


def parse_immuannot_gtf(path):
    """Return {gene: [alleles...]} of HLA-A/-B/-C two-field calls from an Immuannot GTF.

    Immuannot annotates each assembled haplotype contig with an `allele`/`consensus` attribute like
    `HLA-A*01:01:01:01`. We collect the two-field call per locus (one per haplotype).
    """
    calls = {g: [] for g in LOCI}
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            # find HLA-[ABC]*NN:NN in the attribute column
            for m in re.finditer(r"HLA-([ABC])\*(\d+:\d+)", line):
                gene, tf = m.group(1), m.group(2)
                if tf not in calls[gene]:
                    calls[gene].append(tf)
    return calls


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    here = os.path.dirname(os.path.abspath(__file__))
    repo = os.path.dirname(here)
    ap.add_argument("--out", default=os.path.join(repo, "analysis/hprc_benchmark/truth_long.tsv"))
    ap.add_argument("--immuannot-dir", default=None,
                    help="dir with per-sample Immuannot GTFs named <SAMPLE>*.gtf (HG005, parents)")
    args = ap.parse_args()

    rows = []
    samples = []

    # HG002 gold
    for g in LOCI:
        a1, a2 = HG002_GOLD[g]
        rows.append(("HG002", g, a1, a2))
    samples.append("HG002")

    # Immuannot-derived samples (HG003/4/5/6/7) when GTFs are present
    if args.immuannot_dir and os.path.isdir(args.immuannot_dir):
        for gtf in sorted(glob.glob(os.path.join(args.immuannot_dir, "*.gtf"))):
            sample = os.path.basename(gtf).split(".")[0].split("_")[0]
            if sample == "HG002":
                continue  # gold already used; HG002 GTF is only for the control comparison
            calls = parse_immuannot_gtf(gtf)
            wrote = False
            for g in LOCI:
                alleles = [two_field(a) for a in calls[g]]
                if not alleles:
                    continue
                a1 = alleles[0]
                a2 = alleles[1] if len(alleles) > 1 else alleles[0]  # homozygous
                rows.append((sample, g, a1, a2))
                wrote = True
            if wrote:
                samples.append(sample)

    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    with open(args.out, "w") as fh:
        fh.write("sample\tgene\tallele1\tallele2\n")
        for r in rows:
            fh.write("\t".join(r) + "\n")
    print(f"[hprc-truth] wrote {len(rows)} rows for {len(samples)} samples -> {args.out}")
    print(f"[hprc-truth] samples: {', '.join(samples)}")
    if len(samples) == 1:
        print("[hprc-truth] NOTE: only HG002 gold present; run Immuannot (HG005/parents) and re-run "
              "with --immuannot-dir to expand the cohort.")


if __name__ == "__main__":
    main()
