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


def _line_pair(text: str, gene: str):
    for line in text.splitlines():
        values = _alleles(line, gene)
        if len(values) == 2:
            return canonical_pair(values[0], values[1], gene)
        if len(values) == 1 and re.search(r"(?:homo|hom|same)", line, re.I):
            return canonical_pair(values[0], values[0], gene)
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
    return canonical_pair(values[indices[0]], values[indices[1]], target)


def parse_caller_pair(caller: str, text: str, gene: str):
    """Parse a diploid pair from the pinned caller-native result summaries.

    Wrappers may use a normalized line-oriented summary. OptiType's standard
    matrix form is handled explicitly; other pinned wrappers are validated by
    gene-specific allele lines. Unknown or ambiguous content fails closed.
    """
    if caller == "OptiType":
        pair = _optitype_matrix(text, gene)
        if pair:
            return pair
    pair = _line_pair(text, gene)
    if pair:
        return pair
    values = _alleles(text, gene)
    if len(values) == 1:
        return canonical_pair(values[0], values[0], gene)
    if len(values) == 2:
        return canonical_pair(values[0], values[1], gene)
    return None
