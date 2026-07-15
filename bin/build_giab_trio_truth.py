#!/usr/bin/env python3
"""Build GIAB Ashkenazi-trio HLA truth (analysis/giab_trio_benchmark/truth_long.tsv).

Truth: clinical sequence-based HLA typing (Stanford Blood Center, 2016-12-16) of the GIAB Ashkenazi
trio HG002/HG003/HG004, published in Table 1 of "A strategy for building and using a human reference
pangenome" (F1000Research 2019, 8:1751) — full diploid, unphased. HG002 cross-confirmed by Chin et al.
2020 (Nat Commun MHC benchmark, Suppl. Table 4). Mendelian-consistent across the trio (verified).
Two-field HLA-A/-B/-C. Germline (no LOH), non-1000G.
"""
import os

TRIO = {
    "HG002": {"A": ("01:01", "26:01"), "B": ("35:08", "38:01"), "C": ("04:01", "12:03")},
    "HG003": {"A": ("26:01", "30:01"), "B": ("13:02", "38:01"), "C": ("06:02", "12:03")},
    "HG004": {"A": ("01:01", "33:01"), "B": ("14:02", "35:08"), "C": ("04:01", "08:02")},
}
OUT = "/scratch/project_2008084/pihla-publish/analysis/giab_trio_benchmark/truth_long.tsv"


def main():
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    n = 0
    with open(OUT, "w") as fh:
        fh.write("sample\tgene\tallele1\tallele2\n")
        for s in ("HG002", "HG003", "HG004"):
            for g in ("A", "B", "C"):
                a1, a2 = TRIO[s][g]
                fh.write(f"{s}\t{g}\t{a1}\t{a2}\n")
                n += 1
    print(f"[giab-trio-truth] wrote {n} rows (3 samples x A/B/C) -> {OUT}")


if __name__ == "__main__":
    main()
