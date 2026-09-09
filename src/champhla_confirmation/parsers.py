from __future__ import annotations

import csv
import io
import re

from .io import canonical_allele, canonical_pair, normalize_gene

ALLELE_PATTERN = re.compile(
    r"(?<![A-Za-z0-9])(?:HLA-)?[ABC]\*[0-9A-Za-z]+:[0-9A-Za-z]+(?::[0-9A-Za-z]+)*(?:[NLSQCA])?",
    re.I,
)


def _alleles(text: str, gene: str) -> list[str]:
    target = normalize_gene(gene)
    values = []
    for match in ALLELE_PATTERN.findall(text):
        try:
            allele = canonical_allele(match, target)
        except ValueError:
            continue
        values.append(allele)
    return values


def _line_call(text: str, gene: str):
    for line in text.splitlines():
        values = _alleles(line, gene)
        if len(values) == 2:
            pair = canonical_pair(values[0], values[1], gene)
            return {"allele1": pair[0], "allele2": pair[1], "call_status": "callable"}
        if len(values) == 1 and re.search(r"(?:homo|hom|same)", line, re.I):
            pair = canonical_pair(values[0], values[0], gene)
            return {"allele1": pair[0], "allele2": pair[1], "call_status": "callable"}
        if len(values) == 1 and re.search(r"(?:^|[\t ,;])-(?:$|[\t ,;])", line):
            return {"allele1": values[0], "allele2": "", "call_status": "partial"}
    return None


def _optitype_matrix(text: str, gene: str):
    rows = [row for row in csv.reader(io.StringIO(text), delimiter="\t") if row]
    if len(rows) < 2:
        return None
    header = [value.strip().upper() for value in rows[0]]
    target = normalize_gene(gene)
    indices = [i for i, value in enumerate(header) if value in {f"{target}1", f"{target}2"}]
    if len(indices) != 2:
        return None
    values = rows[-1]
    if max(indices) >= len(values):
        return None
    alleles = [values[index].strip() for index in indices]
    present = [value for value in alleles if value not in {"", "-", "."}]
    if len(present) == 2:
        pair = canonical_pair(present[0], present[1], target)
        return {"allele1": pair[0], "allele2": pair[1], "call_status": "callable"}
    if len(present) == 1:
        return {"allele1": canonical_allele(present[0], target), "allele2": "",
                "call_status": "partial"}
    return {"allele1": "", "allele2": "", "call_status": "missing"}


def parse_caller_call(caller: str, text: str, gene: str) -> dict[str, str]:
    """Parse a complete, partial, or missing call from a native summary.

    Wrappers may use a normalized line-oriented summary. OptiType's standard
    matrix form is handled explicitly; other pinned wrappers are validated by
    gene-specific allele lines. A lone allele is retained as partial and is
    never promoted to homozygous without explicit native evidence.
    """
    if caller == "OptiType":
        call = _optitype_matrix(text, gene)
        if call:
            return call
    call = _line_call(text, gene)
    if call:
        return call
    values = _alleles(text, gene)
    if len(values) == 1:
        return {"allele1": values[0], "allele2": "", "call_status": "partial"}
    if len(values) == 2:
        pair = canonical_pair(values[0], values[1], gene)
        return {"allele1": pair[0], "allele2": pair[1], "call_status": "callable"}
    return {"allele1": "", "allele2": "", "call_status": "missing"}


def parse_caller_pair(caller: str, text: str, gene: str):
    """Compatibility adapter returning only complete diploid pairs."""
    call = parse_caller_call(caller, text, gene)
    if call["call_status"] != "callable":
        return None
    return (call["allele1"], call["allele2"])
