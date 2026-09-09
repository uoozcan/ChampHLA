from __future__ import annotations

from datetime import datetime, timezone
import os
from pathlib import Path

from .io import read_json, sha256, write_json


CODE_DIRS = {"src", "scripts", "configs", "models", "templates"}
CODE_FILES = {"pyproject.toml", "README.md", "IMPLEMENTATION_STATUS.md", "ORIGINALS_MANIFEST.json"}


def _tree_hashes(root: Path) -> list[dict]:
    records = []
    for path in sorted(root.rglob("*")):
        if not path.is_file() or any(part in {"__pycache__", ".pytest_cache", ".git"} for part in path.parts):
            continue
        relative = path.relative_to(root)
        if not (relative.name in CODE_FILES and len(relative.parts) == 1
                or relative.parts[0] in CODE_DIRS):
            continue
        records.append({"path": relative.as_posix(), "bytes": path.stat().st_size,
                        "sha256": sha256(path)})
    return records


def freeze_bundle(project_root: str | Path, named_inputs: dict[str, str | Path],
                  output_manifest: str | Path) -> dict:
    root = Path(project_root).resolve()
    target = Path(output_manifest).resolve()
    if not root.is_dir():
        raise ValueError(f"project root does not exist: {root}")
    inputs = {}
    for name, raw_path in sorted(named_inputs.items()):
        path = Path(raw_path).resolve()
        if not path.is_file():
            raise ValueError(f"freeze input {name} is not a file: {path}")
        try:
            relative = path.relative_to(root).as_posix()
        except ValueError as error:
            raise ValueError(f"freeze input {name} is outside the project root: {path}") from error
        inputs[name] = {"path": relative, "path_scope": "repository_relative",
                        "bytes": path.stat().st_size, "sha256": sha256(path)}
    payload = {
        "schema_version": "confirmation-freeze-2",
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "truth_joined": False,
        "policy": "immutable predictions and policy before external truth",
        "project_hash_scope": {"directories": sorted(CODE_DIRS), "files": sorted(CODE_FILES)},
        "project_root": ".",
        "project_root_locator": Path(os.path.relpath(root, target.parent)).as_posix(),
        "project_files": _tree_hashes(root),
        "inputs": inputs,
    }
    write_json(target, payload)
    return payload


def validate_freeze(manifest_path: str | Path) -> dict:
    manifest = read_json(manifest_path)
    failures = []
    manifest_file = Path(manifest_path).resolve()
    if manifest.get("schema_version") == "confirmation-freeze-2":
        root = (manifest_file.parent / manifest["project_root_locator"]).resolve()
    else:
        root = Path(manifest.get("project_root", ""))
    for name, record in manifest.get("inputs", {}).items():
        path = (root / record["path"] if record.get("path_scope") == "repository_relative"
                else Path(record["path"]))
        if not path.is_file():
            failures.append(f"{name}:missing")
        elif sha256(path) != record["sha256"]:
            failures.append(f"{name}:checksum_mismatch")
    for record in manifest.get("project_files", []):
        path = root / record["path"]
        if not path.is_file():
            failures.append(f"project:{record['path']}:missing")
        elif sha256(path) != record["sha256"]:
            failures.append(f"project:{record['path']}:checksum_mismatch")
    return {"valid": not failures, "failures": failures, "checked_inputs": len(manifest.get("inputs", {})),
            "checked_project_files": len(manifest.get("project_files", []))}
