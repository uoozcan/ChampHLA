#!/usr/bin/env python3
"""Build NCI-60 HLA class I truth (truth_long.tsv) from the Adams et al. 2005 SBT table.

Truth source: Adams S et al., "HLA class I and II genotype of the NCI-60 cell lines",
J Transl Med 2005;3:11 (PMC555742), Table 2 (sequence-based typing of class I loci).
The Adams table is mixed-resolution: many alleles are reported only to 1 field (e.g. `03a`,
group-level) or `N.R.` (not resolved). For the 2-field exact-pair benchmark we keep a locus only
when BOTH alleles are resolvable to 2 fields; lower-resolution loci are dropped (not scored) rather
than guessed. Cell-line names are harmonised to the NCI-60 RNA-seq sample aliases (ENA PRJNA433861).

Output: <out>/truth_long.tsv  with columns  sample<TAB>gene<TAB>allele1<TAB>allele2  (2-field).
"""
import argparse
import re
from pathlib import Path


def two_field(tok: str):
    """Return 'NN:NN' if the Adams token is 2-field resolvable, else None."""
    tok = tok.strip()
    if tok in ("", "N.R.", "NR", "-", "n.r."):
        return None
    if tok.endswith("new") or tok.endswith("a"):   # novel or 1-field group only
        return None
    if re.fullmatch(r"\d{4,6}", tok):              # e.g. 020101 -> 02:01 ; 5501 -> 55:01
        return f"{tok[:2]}:{tok[2:4]}"
    return None


def locus_pair(cell: str):
    toks = [t for t in cell.split(",") if t.strip()]
    tf = [two_field(t) for t in toks]
    if not tf or any(v is None for v in tf):
        return None
    if len(tf) == 1:
        return (tf[0], tf[0])                       # homozygous
    return tuple(sorted(tf[:2]))


def canon(name: str) -> str:
    """Canonicalise a cell-line name for matching: upper, strip non-alnum, drop ATCC suffix.

    NCI-60 naming quirks handled: trailing 'ATCC' tag (RNA 'A549-ATCC' = Adams 'A549'); the
    '786-O' (letter O) vs '786-0' (zero) discrepancy between Adams and the ENA aliases.
    """
    c = re.sub(r"[^A-Z0-9]", "", name.upper())
    c = re.sub(r"ATCC$", "", c)
    c = c.replace("786O", "7860")          # Adams '786-O' (O) -> ENA '786-0' (zero)
    return c


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--adams-raw", type=Path,
                    default=Path("/scratch/project_2008084/hla_calibration/nci60/nci60_classI_raw.tsv"))
    ap.add_argument("--rna-index", type=Path,
                    default=Path("/scratch/project_2008084/hla_calibration/nci60/rna/index/sample_fastq_urls.tsv"))
    ap.add_argument("--out", type=Path,
                    default=Path("/scratch/project_2008084/pihla-publish/analysis/nci60_benchmark"))
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)

    # RNA sample names (canonical sequencing aliases) -> map canon -> rna name
    rna_names = [l.split("\t")[0] for l in args.rna_index.read_text().splitlines() if l.strip()]
    rna_by_canon = {canon(n): n for n in rna_names}

    rows = [l.split("\t") for l in args.adams_raw.read_text().splitlines()]
    data = [r for r in rows[1:] if len(r) >= 6]

    truth_lines = ["sample\tgene\tallele1\tallele2"]
    sample_map = ["adams_name\trna_name\tmatched"]
    matched, unmatched_adams = set(), []
    per_locus = {"A": 0, "B": 0, "C": 0}
    for r in data:
        adams_name = r[0].strip()
        rna = rna_by_canon.get(canon(adams_name))
        sample_map.append(f"{adams_name}\t{rna or ''}\t{'1' if rna else '0'}")
        if not rna:
            unmatched_adams.append(adams_name)
            continue
        matched.add(rna)
        for gene, col in (("A", 3), ("B", 4), ("C", 5)):
            pair = locus_pair(r[col])
            if pair:
                per_locus[gene] += 1
                truth_lines.append(f"{rna}\t{gene}\t{pair[0]}\t{pair[1]}")

    (args.out / "truth_long.tsv").write_text("\n".join(truth_lines) + "\n")
    (args.out / "nci60_sample_map.tsv").write_text("\n".join(sample_map) + "\n")

    unmatched_rna = sorted(set(rna_names) - matched)
    print(f"[nci60-truth] Adams data rows: {len(data)}")
    print(f"[nci60-truth] matched to RNA samples: {len(matched)}")
    print(f"[nci60-truth] 2-field scorable loci: A={per_locus['A']} B={per_locus['B']} C={per_locus['C']}"
          f"  (total truth rows={len(truth_lines)-1})")
    if unmatched_adams:
        print(f"[nci60-truth] UNMATCHED Adams names ({len(unmatched_adams)}): {unmatched_adams}")
    if unmatched_rna:
        print(f"[nci60-truth] RNA samples with no Adams match ({len(unmatched_rna)}): {unmatched_rna}")
    print(f"[nci60-truth] wrote {args.out/'truth_long.tsv'} and nci60_sample_map.tsv")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
