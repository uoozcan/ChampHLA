from __future__ import annotations

import csv
import hashlib
import json
import re
from pathlib import Path
from typing import Iterable

TRUTH_TOKENS = {
    "truth", "truth_allele1", "truth_allele2", "label", "labels",
    "is_correct", "is_correct_2field", "correct_status", "candidate_set_oracle",
    "match_grade", "compatibility_grade",
}
UNCALLABLE = {"", "missing", "uncallable", "no_call", "no_consensus", "abstained", "failed"}
ALLELE_RE = re.compile(r"^(?:HLA-)?([A-Za-z0-9]+)\*([0-9A-Za-z]+)(?::([0-9A-Za-z]+))?(?::.*)?$")


def sha256(path: str | Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def sha256_text(value: str) -> str:
    return hashlib.sha256(value.encode("utf-8")).hexdigest()


def read_tsv(path: str | Path) -> list[dict[str, str]]:
    with Path(path).open(newline="", encoding="utf-8", errors="replace") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: str | Path, rows: Iterable[dict], fields: list[str] | None = None) -> None:
    rows = list(rows)
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    if fields is None:
        fields = list(rows[0]) if rows else []
    with target.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields,
                                extrasaction="ignore", lineterminator="\n")
        if fields:
            writer.writeheader()
            writer.writerows(rows)


def read_json(path: str | Path):
    return json.loads(Path(path).read_text(encoding="utf-8"))


def write_json(path: str | Path, value: object) -> None:
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def reject_truth_columns(rows: list[dict[str, str]], context: str) -> None:
    if not rows:
        return
    columns = {column.lower() for column in rows[0]}
    forbidden = sorted(columns & TRUTH_TOKENS)
    if forbidden:
        raise ValueError(f"{context} contains forbidden truth-bearing columns: {forbidden}")


def normalize_modality(value: str) -> str:
    token = value.strip().lower().replace("-", "")
    aliases = {"rna": "rnaseq", "rnaseq": "rnaseq", "wgs": "wgs", "wes": "wes"}
    if token not in aliases:
        raise ValueError(f"unsupported modality: {value!r}")
    return aliases[token]


def normalize_gene(value: str) -> str:
    token = value.strip().upper()
    if token.startswith("HLA-"):
        token = token[4:]
    if token not in {"A", "B", "C"}:
        raise ValueError(f"unsupported primary gene: {value!r}")
    return token


def canonical_allele(value: str, gene: str | None = None) -> str:
    token = value.strip().replace("HLA-", "")
    if not token:
        return ""
    if "*" not in token and gene:
        token = f"{normalize_gene(gene)}*{token}"
    match = ALLELE_RE.match(token)
    if not match or match.group(3) is None:
        raise ValueError(f"allele is not resolvable to two fields: {value!r}")
    allele_gene = normalize_gene(match.group(1))
    if gene and allele_gene != normalize_gene(gene):
        raise ValueError(f"allele/gene mismatch: {value!r} vs {gene!r}")
    return f"{allele_gene}*{match.group(2)}:{match.group(3)}"


def canonical_pair(a: str, b: str, gene: str | None = None) -> tuple[str, str]:
    first = canonical_allele(a, gene)
    second = canonical_allele(b, gene)
    if not first or not second:
        raise ValueError("a callable diploid pair requires two alleles")
    return tuple(sorted((first, second)))


def is_callable(row: dict[str, str]) -> bool:
    status = row.get("call_status", "").strip().lower()
    return status not in UNCALLABLE and bool(row.get("allele1", "").strip() and row.get("allele2", "").strip())


def row_key(row: dict[str, str]) -> tuple[str, str, str, str]:
    return (row.get("cohort", ""), row["subject"], normalize_modality(row["modality"]),
            normalize_gene(row["gene"]))


def require_columns(rows: list[dict[str, str]], required: set[str], context: str) -> None:
    if not rows:
        raise ValueError(f"{context} is empty")
    missing = sorted(required - set(rows[0]))
    if missing:
        raise ValueError(f"{context} missing columns: {missing}")
