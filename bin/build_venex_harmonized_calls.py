#!/usr/bin/env python3
"""
Build venex_harmonized_calls.tsv from PIHLA pipeline output files.

Scans the VENEX results directory and parses all available tool outputs
(OptiType, T1K, HLA-HD, SpecHLA) into a long-format TSV.

Output columns: sample  tool  gene  allele1  allele2  depth1  depth2  reads

Usage:
    python3 bin/build_venex_harmonized_calls.py \
        --results /scratch/project_2008084/hla_calibration/venex/results/wgs \
        --out analysis/venex_harmonized_calls.tsv
"""

import argparse
import csv
import re
from pathlib import Path


CLASSICAL_GENES = {
    "A", "B", "C", "DRB1", "DRB3", "DRB4", "DRB5",
    "DQA1", "DQB1", "DPA1", "DPB1",
}


def strip_hla_prefix(s: str) -> str:
    """Remove leading 'HLA-' from gene or allele strings."""
    if isinstance(s, str) and s.upper().startswith("HLA-"):
        return s[4:]
    return s


def normalise_gene(raw: str) -> str:
    g = strip_hla_prefix(raw.strip())
    return g if g in CLASSICAL_GENES else ""


def parse_optitype(path: Path, sid: str):
    """Parse PIHLA OptiType output: Gene/Allele1/Allele2 TSV with # comments."""
    rows = []
    with path.open(encoding="utf-8") as fh:
        reader = None
        for line in fh:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            if reader is None:
                # First non-comment line is the header
                header = line.split("\t")
                # Expect Gene, Allele1, Allele2
                reader = True
                continue
            parts = line.split("\t")
            if len(parts) < 3:
                continue
            gene = normalise_gene(parts[0])
            if not gene:
                continue
            rows.append({
                "sample": sid, "tool": "OptiType", "gene": gene,
                "allele1": strip_hla_prefix(parts[1].strip()),
                "allele2": strip_hla_prefix(parts[2].strip()),
                "depth1": "", "depth2": "", "reads": "",
            })
    return rows


def parse_t1k(path: Path, sid: str):
    """Parse PIHLA T1K output: same Gene/Allele1/Allele2 format as OptiType."""
    rows = []
    with path.open(encoding="utf-8") as fh:
        header_seen = False
        for line in fh:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            if not header_seen:
                header_seen = True  # skip header row
                continue
            parts = line.split("\t")
            if len(parts) < 3:
                continue
            gene = normalise_gene(parts[0])
            if not gene:
                continue
            rows.append({
                "sample": sid, "tool": "T1K", "gene": gene,
                "allele1": strip_hla_prefix(parts[1].strip()),
                "allele2": strip_hla_prefix(parts[2].strip()),
                "depth1": "", "depth2": "", "reads": "",
            })
    return rows


def parse_hlahd(path: Path, sid: str):
    """Parse PIHLA HLA-HD output: gene/allele1/allele2/depth1/depth2, no header."""
    rows = []
    with path.open(encoding="utf-8") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 3:
                continue
            gene = normalise_gene(parts[0])
            if not gene:
                continue
            depth1 = parts[3].strip() if len(parts) > 3 else ""
            depth2 = parts[4].strip() if len(parts) > 4 else ""
            rows.append({
                "sample": sid, "tool": "HLA-HD", "gene": gene,
                "allele1": strip_hla_prefix(parts[1].strip()),
                "allele2": strip_hla_prefix(parts[2].strip()),
                "depth1": depth1, "depth2": depth2, "reads": "",
            })
    return rows


def parse_spechla(path: Path, sid: str):
    """Parse PIHLA SpecHLA output: wide TSV, columns like HLA_A_1 / HLA_A_2."""
    rows = []
    header = None
    with path.open(encoding="utf-8") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if header is None:
                header = parts
                continue
            # Data row
            record = dict(zip(header, parts))
            # Collect per-gene alleles from column names like HLA_A_1, HLA_DRB1_2
            gene_alleles: dict = {}
            for col, val in record.items():
                # Match HLA_<gene>_<1|2>
                m = re.match(r"^HLA_(.+)_([12])$", col)
                if not m:
                    continue
                gene = m.group(1)
                slot = int(m.group(2))
                if gene not in CLASSICAL_GENES:
                    continue
                gene_alleles.setdefault(gene, ["", ""])
                gene_alleles[gene][slot - 1] = val.strip()
            for gene, (a1, a2) in gene_alleles.items():
                rows.append({
                    "sample": sid, "tool": "SpecHLA", "gene": gene,
                    "allele1": strip_hla_prefix(a1),
                    "allele2": strip_hla_prefix(a2),
                    "depth1": "", "depth2": "", "reads": "",
                })
    return rows


def parse_sample(sample_dir: Path) -> list:
    sid = sample_dir.name
    result_dir = sample_dir / "results" / sid
    all_rows = []

    tool_specs = [
        ("optitype", f"{sid}_optitype.txt",  parse_optitype),
        ("t1k",      f"{sid}_t1k.txt",        parse_t1k),
        ("hlahd",    f"{sid}_hlahd.txt",       parse_hlahd),
        ("spechla",  f"{sid}_spechla.txt",     parse_spechla),
    ]

    for tool_dir, fname, parser in tool_specs:
        fpath = result_dir / tool_dir / fname
        if fpath.exists():
            try:
                rows = parser(fpath, sid)
                all_rows.extend(rows)
            except Exception as exc:
                print(f"  WARNING: {sid}/{tool_dir}: {exc}")

    return all_rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--results",
        type=Path,
        default=Path("/scratch/project_2008084/hla_calibration/venex/results/wgs"),
    )
    parser.add_argument(
        "--out",
        type=Path,
        default=Path("analysis/venex_harmonized_calls.tsv"),
    )
    args = parser.parse_args()

    results_dir = args.results
    if not results_dir.is_dir():
        print(f"ERROR: results directory not found: {results_dir}")
        raise SystemExit(1)

    all_rows = []
    samples_seen = 0
    for sample_dir in sorted(results_dir.iterdir()):
        if not sample_dir.is_dir():
            continue
        rows = parse_sample(sample_dir)
        if rows:
            all_rows.extend(rows)
            samples_seen += 1

    args.out.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = ["sample", "tool", "gene", "allele1", "allele2", "depth1", "depth2", "reads"]
    with args.out.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(all_rows)

    print(f"Wrote {len(all_rows)} rows from {samples_seen} samples → {args.out}")

    # Summary by tool
    from collections import Counter
    tool_counts = Counter(r["tool"] for r in all_rows)
    for tool, n in sorted(tool_counts.items()):
        n_samples = len({r["sample"] for r in all_rows if r["tool"] == tool})
        print(f"  {tool}: {n} rows, {n_samples} samples")


if __name__ == "__main__":
    main()
