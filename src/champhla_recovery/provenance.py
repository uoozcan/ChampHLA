from __future__ import annotations

from pathlib import Path

from .io import read_tsv, sha256, write_json


def build_originals_manifest(source_map: str, project_root: str, output: str) -> dict:
    root = Path(project_root).resolve()
    records = []
    failures = []
    for mapping in read_tsv(source_map):
        original = Path(mapping["original_path"])
        copied = root / mapping["copied_path"]
        if not original.is_file():
            failures.append(f"missing original: {original}")
            continue
        if not copied.is_file():
            failures.append(f"missing copy: {copied}")
            continue
        original_hash = sha256(original)
        copied_hash = sha256(copied)
        if original_hash != copied_hash:
            failures.append(f"copy differs from original: {copied}")
        records.append({
            "source_role": mapping["source_role"],
            "original_path": str(original.resolve()),
            "copied_path": mapping["copied_path"],
            "original_bytes": original.stat().st_size,
            "original_mtime_ns": original.stat().st_mtime_ns,
            "original_sha256": original_hash,
            "copied_sha256": copied_hash,
            "hashes_match": original_hash == copied_hash,
            "source_git_state": mapping.get("source_git_state", ""),
        })
    result = {
        "schema_version": "champhla-publication-originals-1",
        "policy": "read_only_sources_and_manuscript",
        "records": records,
        "failures": failures,
        "passed": not failures,
    }
    write_json(output, result)
    return result

