#!/usr/bin/env python3
"""Build IHWG/IHIW reference-panel HLA truth from the IPD-IMGT/HLA IHIW cell typing table.

Source: IPD-IMGT/HLA "IHIW Cell Typing Table" (open xlsx,
https://www.ebi.ac.uk/ipd/imgt/hla/cells/ihiw/), sheet 'Expressed Typings' — 366 IHWG reference
B-LCLs, multi-lab consensus, FULL (4-field) resolution. Encoding: comma = the two genotype alleles;
a single value = homozygous (allele1 == allele2).

Emits (schema sample<TAB>gene<TAB>allele1<TAB>allele2, sample = IHW number), HLA-A/-B/-C (+ DRB1/DQB1):
  - truth_long_4field.tsv : full resolution (e.g. 24:02:01:01)
  - truth_long.tsv        : two-field collapse (e.g. 24:02)
Plus ihiw_line_map.tsv (IHW_number, primary_name, ancestry, homozygous) for matching ENA reads -> truth.
"""
import os
import openpyxl

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
XLSX = os.path.join(REPO, "analysis/ihwg_benchmark/ihiw_cell_typing_table.xlsx")
OUT = os.path.join(REPO, "analysis/ihwg_benchmark")
SHEET = "Expressed Typings"
# 0-based columns in 'Expressed Typings' (header row 5; data from row 6)
COL = {"ihw": 0, "name": 1, "ancestry": 2, "homozygous": 4,
       "A": 6, "B": 7, "C": 8, "DRB1": 13, "DQB1": 19}
LOCI = ("A", "B", "C")


def parse_geno(raw):
    """'a, b' -> (a,b); single 'a' -> (a,a) homozygous; >2 -> first two; '' -> None."""
    if raw is None:
        return None
    s = str(raw).strip().replace("/", ",")  # treat the rare slash like a separator
    if not s:
        return None
    parts = [p.strip() for p in s.split(",") if p.strip()]
    if not parts:
        return None
    if len(parts) == 1:
        return (parts[0], parts[0])
    return (parts[0], parts[1])


def two_field(a):
    p = a.split(":")
    return ":".join(p[:2]) if len(p) >= 2 else a


def main():
    wb = openpyxl.load_workbook(XLSX, read_only=True, data_only=True)
    ws = wb[SHEET]
    rows4, rows2, lmap = [], [], []
    samples = set()
    for r in ws.iter_rows(min_row=6, values_only=True):
        ihw = r[COL["ihw"]]
        if not ihw or not str(ihw).strip().upper().startswith("IHW"):
            continue
        ihw = str(ihw).strip()
        lmap.append((ihw, str(r[COL["name"]] or "").strip(),
                     str(r[COL["ancestry"]] or "").strip(), str(r[COL["homozygous"]] or "").strip()))
        wrote = False
        for g in LOCI:
            geno = parse_geno(r[COL[g]])
            if geno is None:
                continue
            a1, a2 = geno
            rows4.append((ihw, g, a1, a2))
            rows2.append((ihw, g, two_field(a1), two_field(a2)))
            wrote = True
        if wrote:
            samples.add(ihw)

    os.makedirs(OUT, exist_ok=True)
    for fname, rows in (("truth_long_4field.tsv", rows4), ("truth_long.tsv", rows2)):
        with open(os.path.join(OUT, fname), "w") as fh:
            fh.write("sample\tgene\tallele1\tallele2\n")
            for x in rows:
                fh.write("\t".join(x) + "\n")
        print(f"[ihwg] wrote {len(rows)} rows ({len(samples)} samples) -> {fname}")
    with open(os.path.join(OUT, "ihiw_line_map.tsv"), "w") as fh:
        fh.write("ihw_number\tprimary_name\tancestry\thomozygous\n")
        for x in lmap:
            fh.write("\t".join(x) + "\n")
    print(f"[ihwg] wrote line map ({len(lmap)} lines)")


if __name__ == "__main__":
    main()
