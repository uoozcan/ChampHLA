#!/usr/bin/env python3
"""Build GeT-RM HLA truth tables from the CDC consolidated PGx-HLA Excel.

Source: CDC GeT-RM "Consolidated PGx-HLA table" (sheet 'PGx HLA Genotypes'), the open
publication of Bettinotti et al., J Mol Diagn 2018 (PMC6939753) — 108 Coriell NIGMS DNA
reference materials characterised at all 11 classical HLA loci by multi-lab PCR-SSO + SBT,
reported at three-field resolution (e.g. *02:01:01).

Emits two long-format truth files (schema: sample<TAB>gene<TAB>allele1<TAB>allele2) for
HLA-A/-B/-C:
  - truth_long_3field.tsv : full sub-allele resolution (e.g. 02:01:01)
  - truth_long.tsv        : collapsed to two-field (e.g. 02:01)

NOTE (read availability): these 108 lines have essentially no public short-read WGS/WES/RNA
(ENA: only NA12273 RNA-Seq + NA17221 WGS usable; 33 amplicon-only). The truth is retained as
an integrity-clean experimental resource; it is not a viable standalone benchmark cohort.
See analysis/getrm_benchmark/PHASE0_STATUS.md.
"""
import os
import openpyxl

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
XLSX = os.path.join(REPO, "analysis/getrm_benchmark/getrm_consolidated_PGx-HLA_table.xlsx")
OUT_DIR = os.path.join(REPO, "analysis/getrm_benchmark")
SHEET = "PGx HLA Genotypes"

# 0-based column indices in the 'PGx HLA Genotypes' sheet (header on row 2)
COL_ID = 1
LOCUS_COLS = {"A": (74, 75), "B": (77, 78), "C": (80, 81)}  # (allele1, allele2)


def norm(raw):
    """'*02:01:01' / '02:01:01N' -> '02:01:01' (strip leading * and surrounding space)."""
    if raw is None:
        return None
    s = str(raw).strip().lstrip("*").strip()
    if not s or s.lower() in ("not typed", "nt", "na", "none", "-"):
        return None
    return s


def to_two_field(allele):
    if allele is None:
        return None
    parts = allele.split(":")
    return ":".join(parts[:2]) if len(parts) >= 2 else allele


def main():
    wb = openpyxl.load_workbook(XLSX, read_only=True, data_only=True)
    ws = wb[SHEET]
    rows3, rows2 = [], []
    n_samples = set()
    for r in ws.iter_rows(min_row=3, values_only=True):
        sid = r[COL_ID]
        if not sid:
            continue
        sid = str(sid).strip()
        wrote = False
        for gene, (i1, i2) in LOCUS_COLS.items():
            a1, a2 = norm(r[i1]), norm(r[i2])
            if a1 is None and a2 is None:
                continue
            # homozygous: a single typed allele fills both
            a1 = a1 or a2
            a2 = a2 or a1
            rows3.append((sid, gene, a1, a2))
            rows2.append((sid, gene, to_two_field(a1), to_two_field(a2)))
            wrote = True
        if wrote:
            n_samples.add(sid)

    os.makedirs(OUT_DIR, exist_ok=True)
    for fname, rows in (("truth_long_3field.tsv", rows3), ("truth_long.tsv", rows2)):
        path = os.path.join(OUT_DIR, fname)
        with open(path, "w") as fh:
            fh.write("sample\tgene\tallele1\tallele2\n")
            for row in rows:
                fh.write("\t".join(row) + "\n")
        print(f"[getrm] wrote {len(rows)} rows ({len(n_samples)} samples) -> {path}")


if __name__ == "__main__":
    main()
