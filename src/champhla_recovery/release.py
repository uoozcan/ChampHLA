from __future__ import annotations

import subprocess
from pathlib import Path

from .io import sha256, write_json


EXCLUDED_PARTS = {".git", "__pycache__", ".pytest_cache", ".venv", "runs", "runtime_data", "truth_sealed"}


def freeze_release(root: str, output: str) -> dict:
    project = Path(root).resolve()
    target = Path(output).resolve()
    files = []
    for path in sorted(project.rglob("*")):
        if not path.is_file() or path.resolve() == target:
            continue
        relative = path.relative_to(project)
        if any(part in EXCLUDED_PARTS for part in relative.parts):
            continue
        files.append({"path": str(relative), "bytes": path.stat().st_size, "sha256": sha256(path)})
    try:
        commit = subprocess.run(
            ["git", "-C", str(project), "rev-parse", "HEAD"], check=True,
            capture_output=True, text=True,
        ).stdout.strip()
        dirty = bool(subprocess.run(
            ["git", "-C", str(project), "status", "--porcelain"], check=True,
            capture_output=True, text=True,
        ).stdout.strip())
    except (subprocess.CalledProcessError, FileNotFoundError):
        commit, dirty = "not_a_git_repository", True
    result = {
        "schema_version": "champhla-release-freeze-1",
        "project_root": str(project), "git_commit": commit, "git_dirty": dirty,
        "files": files, "file_count": len(files),
        "truth_data_included": False,
    }
    write_json(target, result)
    return result

