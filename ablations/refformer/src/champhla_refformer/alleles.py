"""HLA normalization utilities shared by indexing, training, and inference."""

from __future__ import annotations

import re
from typing import Iterable, Tuple


SUFFIX_RE = re.compile(r"([NLSCAQ])$")


def clean(value) -> str:
    return "" if value is None else str(value).strip()


def normalize_allele(value, gene: str | None = None, fields: int = 2) -> str:
    text = clean(value).upper().replace("HLA-", "")
    if not text or text in {"-", "NA", "N/A", "NONE", "NULL", "."}:
        return ""
    text = text.split("/", 1)[0].split(";", 1)[0].strip()
    suffix = ""
    match = SUFFIX_RE.search(text)
    if match:
        suffix = match.group(1)
        text = text[:-1]
    if "*" not in text:
        if not gene:
            return ""
        text = f"{gene.upper()}*{text}"
    allele_gene, number = text.split("*", 1)
    allele_gene = allele_gene.strip()
    if gene and allele_gene != gene.upper():
        return ""
    number = number.replace("_", ":")
    if ":" not in number and number.isdigit() and len(number) >= 4:
        number = ":".join(number[index:index + 2] for index in range(0, len(number), 2))
    parts = [part for part in number.split(":") if part]
    if len(parts) < 2:
        return ""
    normalized = f"{allele_gene}*{':'.join(parts[:fields])}"
    if fields > 2 and suffix:
        normalized += suffix
    return normalized


def canonical_pair(allele1, allele2, gene: str | None = None) -> Tuple[str, str]:
    pair = [normalize_allele(allele1, gene), normalize_allele(allele2, gene)]
    if not all(pair):
        return ("", "")
    return tuple(sorted(pair))


def pair_key(pair: Iterable[str]) -> str:
    values = tuple(sorted(clean(item) for item in pair))
    return "+".join(values)


def allele_gene(allele: str) -> str:
    return normalize_allele(allele).split("*", 1)[0] if normalize_allele(allele) else ""

