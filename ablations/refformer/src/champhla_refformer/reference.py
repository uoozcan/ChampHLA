"""Pinned IPD-IMGT/HLA reference indexing and alias resolution."""

from __future__ import annotations

import csv
import hashlib
import json
import re
from collections import defaultdict
from pathlib import Path

from . import GENES
from .alleles import normalize_allele


RENAME_RE = re.compile(r"(?:renamed|identical to)\s+(?:HLA-)?([A-Z0-9]+\*[0-9:]+[NLSCAQ]?)", re.I)
ALIGNMENT_ALLELE_RE = re.compile(r"^\s+([A-Z0-9]+\*\S+)\s+(.+?)\s*$")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def parse_fasta(path: Path):
    header = None
    sequence = []
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(sequence)
                header, sequence = line[1:], []
            else:
                sequence.append(line.upper())
    if header is not None:
        yield header, "".join(sequence)


def fasta_allele(header: str) -> str:
    fields = header.split()
    return fields[1] if len(fields) > 1 else ""


def parse_nucleotide_alignment(path: Path, gene: str) -> tuple[dict[str, dict], dict]:
    """Expand IMGT dash notation into a shared, gap-preserving cDNA coordinate."""
    blocks: list[dict[str, str]] = []
    current: dict[str, str] = {}
    with path.open(encoding="utf-8") as handle:
        for raw_line in handle:
            if raw_line.startswith(" cDNA"):
                if current:
                    blocks.append(current)
                    current = {}
                continue
            match = ALIGNMENT_ALLELE_RE.match(raw_line.rstrip("\n"))
            if not match or not match.group(1).startswith(gene + "*"):
                continue
            allele, rendered = match.groups()
            current[allele] = "".join(rendered.split()).replace("|", "")
    if current:
        blocks.append(current)

    expanded: dict[str, list[str]] = defaultdict(list)
    usable_blocks, ignored_blocks, total_columns = 0, 0, 0
    for block in blocks:
        reference_name = next(iter(block), None)
        if reference_name is None:
            continue
        reference = block[reference_name]
        # A terminal allele-specific extension has no shared reference and must
        # not be projected onto the common coordinate.
        if len(block) < 10 or not any(base in reference for base in "ACGT"):
            ignored_blocks += 1
            continue
        if any(len(chunk) != len(reference) for chunk in block.values()):
            raise ValueError(f"inconsistent IMGT alignment block width: {path}")
        for allele, chunk in block.items():
            resolved = []
            for reference_base, base in zip(reference, chunk):
                if base == "-":
                    base = reference_base
                if base == ".":
                    base = "-"
                elif base == "*" or base not in "ACGTN-":
                    base = "N"
                resolved.append(base)
            expanded[allele].append("".join(resolved))
        total_columns += len(reference)
        usable_blocks += 1

    records = {}
    for allele, chunks in expanded.items():
        aligned = "".join(chunks)
        if len(aligned) != total_columns:
            continue
        records[allele] = {
            "aligned_sequence": aligned,
            "alignment_mask": "".join("1" if base in "ACGT" else "0" for base in aligned),
        }
    metadata = {
        "shared_columns": total_columns,
        "usable_blocks": usable_blocks,
        "ignored_allele_specific_blocks": ignored_blocks,
        "aligned_alleles": len(records),
        "gap_character": "-",
        "unknown_character": "N",
        "dash_notation_expanded": True,
    }
    return records, metadata


def load_aliases(history_path: Path) -> dict[str, str]:
    aliases = {}
    with history_path.open(newline="", encoding="utf-8") as handle:
        rows = (line for line in handle if line and not line.startswith("#"))
        reader = csv.reader(rows)
        header = next(reader)
        current_index = header.index("3590")
        for row in reader:
            if len(row) <= current_index:
                continue
            current = normalize_allele(row[current_index])
            current_full = normalize_allele(row[current_index], fields=4)
            if not current:
                continue
            for historical in row[1:]:
                old = normalize_allele(historical)
                old_full = normalize_allele(historical, fields=4)
                # Preserve higher-field history even when both names collapse to
                # the same two-field group. Runtime output remains two-field.
                if old_full and old_full != current_full:
                    aliases.setdefault(old_full, current)
                if old and old != current:
                    aliases.setdefault(old, current)
    return aliases


def load_deleted(path: Path) -> tuple[set[str], dict[str, str]]:
    deleted, renamed = set(), {}
    with path.open(newline="", encoding="utf-8") as handle:
        rows = (line for line in handle if line and not line.startswith("#"))
        for row in csv.DictReader(rows):
            old = normalize_allele(row.get("Allele"))
            if not old:
                continue
            deleted.add(old)
            match = RENAME_RE.search(row.get("Description", ""))
            if match:
                new = normalize_allele(match.group(1))
                if new:
                    renamed[old] = new
    return deleted, renamed


def build_index(source_root: Path, output: Path) -> dict:
    source_root = source_root.resolve()
    groups = defaultdict(lambda: {"sequences": [], "source_alleles": [], "sequence_masks": [],
                                  "aligned_sequences": [], "alignment_masks": [], "partial": False})
    files = []
    alignment_metadata = {}
    for gene in GENES:
        fasta = source_root / "fasta" / f"{gene}_nuc.fasta"
        alignment = source_root / "alignments" / f"{gene}_nuc.txt"
        if not fasta.exists() or not alignment.exists():
            raise FileNotFoundError(f"missing pinned IMGT asset for {gene}")
        files.extend([fasta, alignment])
        aligned_records, alignment_metadata[gene] = parse_nucleotide_alignment(alignment, gene)
        for header, sequence in parse_fasta(fasta):
            full = fasta_allele(header)
            two_field = normalize_allele(full, gene)
            if not two_field or not sequence:
                continue
            record = groups[two_field]
            if sequence not in record["sequences"]:
                record["sequences"].append(sequence)
                record["sequence_masks"].append("1" * len(sequence))
            aligned = aligned_records.get(full)
            if aligned and aligned["aligned_sequence"] not in record["aligned_sequences"]:
                record["aligned_sequences"].append(aligned["aligned_sequence"])
                record["alignment_masks"].append(aligned["alignment_mask"])
            record["source_alleles"].append(full)
            record["partial"] = record["partial"] or "partial" in header.lower() or "N" in sequence

    history = source_root / "Allelelist_history.txt"
    deleted_path = source_root / "Deleted_alleles.txt"
    release = source_root / "release_version.txt"
    licence = source_root / "LICENCE.md"
    files.extend([history, deleted_path, release, licence, source_root / "Allelelist.txt"])
    aliases = load_aliases(history)
    deleted, renamed = load_deleted(deleted_path)
    aliases.update(renamed)
    payload = {
        "schema_version": "imgt-reference-index-2",
        "release": "IPD-IMGT/HLA 3.59.0",
        "git_commit": "f03bbe951d3106ff00f1801407a7a13278c8866b",
        "genes": list(GENES),
        "groups": {key: value for key, value in sorted(groups.items())},
        "aliases": dict(sorted(aliases.items())),
        "deleted_two_field": sorted(deleted),
        "alignment_metadata": alignment_metadata,
        "files": [{"path": str(path.relative_to(source_root)), "bytes": path.stat().st_size,
                   "sha256": sha256(path)} for path in sorted(files)],
        "alignment_policy": "IMGT dash notation is expanded on a shared per-gene cDNA coordinate; gaps and unsequenced positions are explicitly masked",
    }
    encoded = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    payload["index_sha256"] = hashlib.sha256(encoded).hexdigest()
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return payload


class ReferenceIndex:
    def __init__(self, path: str | Path):
        self.path = Path(path)
        self.payload = json.loads(self.path.read_text(encoding="utf-8"))
        checksum = self.payload.get("index_sha256")
        unsigned = dict(self.payload)
        unsigned.pop("index_sha256", None)
        observed = hashlib.sha256(json.dumps(unsigned, sort_keys=True, separators=(",", ":")).encode()).hexdigest()
        if checksum != observed:
            raise ValueError("IMGT index checksum mismatch")

    @property
    def checksum(self) -> str:
        return self.payload["index_sha256"]

    def resolve(self, allele: str) -> tuple[str, dict]:
        full = normalize_allele(allele, fields=4)
        normalized = normalize_allele(allele)
        original = normalized
        flags = {"alias_resolved": False, "deleted": False, "unresolved": False, "partial": False}
        alias_key = full if full in self.payload["aliases"] else normalized
        if alias_key in self.payload["aliases"]:
            normalized = self.payload["aliases"][alias_key]
            flags["alias_resolved"] = True
        flags["deleted"] = original in self.payload["deleted_two_field"] or normalized in self.payload["deleted_two_field"]
        record = self.payload["groups"].get(normalized)
        if record is None:
            flags["unresolved"] = True
        else:
            flags["partial"] = bool(record.get("partial"))
        return normalized, flags

    def sequences(self, allele: str) -> tuple[list[str], dict]:
        resolved, flags = self.resolve(allele)
        record = self.payload["groups"].get(resolved)
        return (list(record.get("sequences", [])) if record else []), flags
