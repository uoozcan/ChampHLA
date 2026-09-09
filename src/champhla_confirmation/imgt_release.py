"""Identify an IPD-IMGT/HLA release from database content rather than from a filename.

A caller ships an HLA database without always declaring which IPD-IMGT/HLA release it was
built from. A filename such as ``a_complete.3100.new.eb.fasta`` is suggestive but it is not
provenance, and in this project two such filenames turned out to be wrong. This module
identifies the release from what a database actually contains.

Comparison runs at **two-field resolution**. Comparing full allele designations across
sources is unsound: IPD-IMGT/HLA renames a two-field allele to three fields once subtypes
are discovered, so ``A*32:55`` in one source is ``A*32:55:01`` in another while both denote
the same two-field type. Two fields is also the resolution this project scores at.

The identification is the earliest release that contains every two-field type in the
database: a database cannot contain types its source release had not yet published.
Coverage -- the share of that release's types the database carries -- distinguishes an
unfiltered database (coverage near 1.0, falling sharply at the next release) from a
deliberately filtered one such as a complete-sequences-only subset.

Comparison runs against a pinned local cache of official per-release allele lists so that
identification is reproducible offline and every input is checksummed. When no release
contains the set the result is ``NO_MATCH``; the module never guesses and never lets a
filename decide.
"""
from __future__ import annotations

import re
from pathlib import Path

from .io import sha256

RELEASE_FILE_RE = re.compile(r"Allelelist\.(\d+)\.(\d+)\.(\d+)\.txt$")
RELEASE_RE = re.compile(r"\b(\d+)\.(\d+)\.(\d+)\b")
EXPRESSION_SUFFIX = re.compile(r"[A-Z]$")
POLYSOLVER_FIELD = re.compile(r"(\d+)([a-zA-Z])")

SOURCE_KINDS = ("allele_list", "fasta_headers", "polysolver_headers", "declared_header")
EVIDENCE_KINDS = ("in_file_header", "content_match", "version_file", "upstream_commit")

HETEROGENEOUS = "HETEROGENEOUS"
NO_MATCH = "NO_MATCH"
NOT_RELEASE_BEARING = "NOT_RELEASE_BEARING"

# A database carrying at least this share of its release's two-field types is treated as
# unfiltered, so its release is identified sharply rather than as a lower bound.
UNFILTERED_COVERAGE = 0.99


def release_key(release: str) -> tuple[int, int, int]:
    """Order releases numerically so that 3.9.0 precedes 3.10.0."""
    match = RELEASE_RE.search(release)
    if not match:
        raise ValueError(f"not a release designation: {release!r}")
    return (int(match.group(1)), int(match.group(2)), int(match.group(3)))


def two_field(name: str) -> str | None:
    """Reduce an allele designation to its gene and first two numeric fields."""
    cleaned = name.strip().upper()
    if cleaned.startswith("HLA-"):
        cleaned = cleaned[4:]
    if "*" not in cleaned:
        return None
    gene, _, rest = cleaned.partition("*")
    if not gene:
        return None
    fields = [EXPRESSION_SUFFIX.sub("", field) for field in rest.split(":")[:2]]
    if not fields or not all(field.isdigit() for field in fields):
        return None
    return f"{gene}*{':'.join(fields)}"


def extract_two_field_types(path: str | Path, kind: str) -> set[str]:
    """Extract the two-field types a database source actually contains."""
    if kind not in SOURCE_KINDS or kind == "declared_header":
        raise ValueError(f"unsupported allele source kind: {kind}")
    types: set[str] = set()
    for line in Path(path).read_text(encoding="utf-8", errors="replace").splitlines():
        stripped = line.strip()
        if not stripped or stripped.startswith("#"):
            continue
        if kind == "polysolver_headers":
            parts = stripped.lstrip(">").split("_")
            if len(parts) < 4 or parts[0].lower() != "hla":
                continue
            fields = parts[2:]
            match = POLYSOLVER_FIELD.fullmatch(fields[-1])
            if match:
                fields[-1] = match.group(1)
            candidate = two_field(f"{parts[1]}*{':'.join(fields)}")
            if candidate:
                types.add(candidate)
            continue
        if kind == "fasta_headers" and not stripped.startswith(">"):
            continue
        for token in re.split(r"[,\s]+", stripped.lstrip(">")):
            if "*" in token:
                candidate = two_field(token)
                if candidate:
                    types.add(candidate)
                break
    return types


def declared_release(path: str | Path) -> str | None:
    """Read a release from a source that declares its own version in its header block.

    Only the leading commented header is inspected, so a release-like number in the body of
    a file cannot be mistaken for a declaration.
    """
    for line in Path(path).read_text(encoding="utf-8", errors="replace").splitlines():
        stripped = line.strip()
        if not stripped:
            continue
        if not stripped.startswith("#"):
            break
        match = RELEASE_RE.search(stripped)
        if match:
            return ".".join(match.groups())
    return None


def load_release_index(cache_dir: str | Path) -> dict[str, dict]:
    """Load the pinned per-release allele lists, recording each file's checksum."""
    root = Path(cache_dir)
    if not root.is_dir():
        raise ValueError(f"allele-list cache is not a directory: {cache_dir}")
    index: dict[str, dict] = {}
    for path in sorted(root.iterdir()):
        match = RELEASE_FILE_RE.search(path.name)
        if not path.is_file() or not match:
            continue
        release = ".".join(match.groups())
        if release in index:
            raise ValueError(f"duplicate allele list for release {release}")
        index[release] = {
            "types": extract_two_field_types(path, "allele_list"),
            "source": path.name,
            "sha256": sha256(path),
        }
    if not index:
        raise ValueError(f"allele-list cache contains no release files: {cache_dir}")
    return index


def identify_release(types: set[str], index: dict[str, dict]) -> dict:
    """Identify the release a database was built from, at two-field resolution."""
    if not types:
        raise ValueError("cannot identify a release from an empty type set")
    genes = {value.split("*", 1)[0] for value in types}
    ordered = sorted(index, key=release_key)

    measured = []
    for release in ordered:
        scoped = {value for value in index[release]["types"] if value.split("*", 1)[0] in genes}
        if not scoped:
            continue
        measured.append({
            "release": release,
            "missing": len(types - scoped),
            "coverage": len(types & scoped) / len(scoped),
            "release_types": len(scoped),
        })

    contained = [row for row in measured if row["missing"] == 0]
    if not contained:
        closest = min(measured, key=lambda row: row["missing"])
        return {
            "release": NO_MATCH,
            "evidence": "content_match",
            "two_field_types": len(types),
            "candidate_window": [measured[0]["release"], measured[-1]["release"]],
            "closest_release": closest["release"],
            "missing_from_closest": closest["missing"],
        }

    match = contained[0]
    position = measured.index(match)
    successor = measured[position + 1] if position + 1 < len(measured) else None
    filtered = match["coverage"] < UNFILTERED_COVERAGE

    record = {
        "release": match["release"],
        "evidence": "content_match",
        "two_field_types": len(types),
        "candidate_window": [measured[0]["release"], measured[-1]["release"]],
        "coverage_at_release": round(match["coverage"], 6),
        "release_types_in_scope": match["release_types"],
        "filtered_subset": filtered,
        "allele_list_source": index[match["release"]]["source"],
        "allele_list_sha256": index[match["release"]]["sha256"],
    }
    if position:
        previous = measured[position - 1]
        record["predecessor"] = previous["release"]
        record["missing_from_predecessor"] = previous["missing"]
    if successor:
        record["successor"] = successor["release"]
        record["coverage_at_successor"] = round(successor["coverage"], 6)
    if filtered:
        record["interpretation"] = (
            "database is a filtered subset of its release; the release is the earliest that "
            "contains every type it carries"
        )
    return record


def caller_release_summary(component_releases: list[str]) -> str:
    """Summarise a caller from its components without hiding disagreement.

    A caller whose components come from different releases is heterogeneous. That is a
    property of the deployed database, not a failure to determine it, so it is reported as
    such rather than collapsed into one release.
    """
    informative = {
        value for value in component_releases
        if value and value != NOT_RELEASE_BEARING
    }
    if not informative:
        return NO_MATCH
    if NO_MATCH in informative:
        return NO_MATCH
    if len(informative) == 1:
        return next(iter(informative))
    return HETEROGENEOUS


def resolve_source(path: str | Path, kind: str, cache_dir: str | Path | None = None) -> dict:
    """Resolve one database source, preferring a declared release over an inferred one."""
    if kind == "declared_header":
        release = declared_release(path)
        return {
            "release": release or NO_MATCH,
            "evidence": "in_file_header",
            "component_sha256": sha256(path),
        }
    if cache_dir is None:
        raise ValueError("content matching requires a pinned allele-list cache")
    record = identify_release(extract_two_field_types(path, kind), load_release_index(cache_dir))
    record["component_sha256"] = sha256(path)
    return record
