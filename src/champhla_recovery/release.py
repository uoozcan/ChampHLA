from __future__ import annotations

import fnmatch
import gzip
import io
import json
import os
import re
import shutil
import subprocess
import tarfile
import tempfile
from pathlib import Path

from .io import read_json, sha256, write_json


HEX64 = re.compile(r"^[0-9a-f]{64}$")
SAFE_RUN_ID = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]{0,127}$")
HOUSEKEEPING_PARTS = {".git", "__pycache__", ".pytest_cache", ".venv"}
SECRET_CONTENT = (
    re.compile(r"-----BEGIN (?:OPENSSH |RSA |EC )?PRIVATE KEY-----"),
    re.compile(r"(?i)aws_secret_access_key\s*[:=]"),
)


def _git(project: Path, *args: str) -> str:
    return subprocess.run(
        ["git", "-C", str(project), *args], check=True,
        capture_output=True, text=True,
    ).stdout


def _load_profile(policy_path: str | Path, profile: str) -> tuple[dict, str]:
    policy_file = Path(policy_path).resolve()
    policy = read_json(policy_file)
    if policy.get("schema_version") != "champhla-release-allowlist-1":
        raise ValueError("unsupported release allow-list schema")
    profiles = policy.get("profiles", {})
    if profile not in profiles:
        raise ValueError(f"unknown release allow-list profile: {profile}")
    return profiles[profile], sha256(policy_file)


def _matches(path: str, patterns: list[str], run_id: str | None = None) -> bool:
    for pattern in patterns:
        if run_id is not None:
            pattern = pattern.replace("{run_id}", run_id)
        variants = {pattern, pattern.replace("/**/", "/")}
        if any(fnmatch.fnmatchcase(path, candidate) for candidate in variants):
            return True
    return False


def _deny_reason(relative: Path, profile: dict) -> str | None:
    posix = relative.as_posix()
    lower = posix.lower()
    parts = {part.lower() for part in relative.parts}
    denied_parts = {part.lower() for part in profile.get("denied_parts", [])}
    if parts & denied_parts:
        return "prohibited path component"
    if any(lower.endswith(suffix.lower()) for suffix in profile.get("denied_suffixes", [])):
        return "prohibited payload suffix"
    if any(re.search(pattern, posix, flags=re.IGNORECASE)
           for pattern in profile.get("denied_name_patterns", [])):
        return "prohibited or secret-bearing name"
    return None


def _content_secret(path: Path, profile: dict) -> bool:
    if path.suffix.lower() not in set(profile.get("content_scan_suffixes", [])):
        return False
    if path.stat().st_size > int(profile.get("content_scan_max_bytes", 5_000_000)):
        return False
    try:
        text = path.read_text(encoding="utf-8")
    except (UnicodeDecodeError, OSError):
        return False
    return any(pattern.search(text) for pattern in SECRET_CONTENT)


def _audit_file(path: Path, relative: Path, profile: dict,
                run_id: str | None = None) -> list[str]:
    failures: list[str] = []
    posix = relative.as_posix()
    if path.is_symlink():
        failures.append(f"{posix}: symlinks are prohibited")
        return failures
    reason = _deny_reason(relative, profile)
    if reason:
        failures.append(f"{posix}: {reason}")
    if path.stat().st_size > int(profile["max_file_bytes"]):
        failures.append(f"{posix}: exceeds maximum compact-file size")
    if not _matches(posix, profile["allowed_patterns"], run_id):
        failures.append(f"{posix}: path is not explicitly allow-listed")
    if _content_secret(path, profile):
        failures.append(f"{posix}: private-key or credential material detected")
    return failures


def _tracked_entries(project: Path) -> list[tuple[Path, int]]:
    raw = subprocess.run(
        ["git", "-C", str(project), "ls-files", "--stage", "-z"], check=True,
        capture_output=True,
    ).stdout
    entries: list[tuple[Path, int]] = []
    for item in raw.split(b"\0"):
        if not item:
            continue
        metadata, name = item.split(b"\t", 1)
        git_mode = metadata.split(b" ", 1)[0]
        entries.append((Path(name.decode("utf-8")), 0o755 if git_mode == b"100755" else 0o644))
    return entries


def _scan_local_prohibited(project: Path, profile: dict, target: Path) -> list[str]:
    """Detect prohibited payloads even when Git ignores or does not track them."""
    failures: list[str] = []
    for directory, names, filenames in os.walk(project, topdown=True, followlinks=False):
        names[:] = sorted(name for name in names if name not in HOUSEKEEPING_PARTS)
        for name in sorted(filenames):
            path = Path(directory) / name
            relative = path.relative_to(project)
            if path.resolve() == target:
                continue
            reason = _deny_reason(relative, profile)
            if path.is_symlink():
                failures.append(f"{relative.as_posix()}: symlink present in release tree")
            elif reason:
                failures.append(f"{relative.as_posix()}: {reason} present in release tree")
            elif _content_secret(path, profile):
                failures.append(
                    f"{relative.as_posix()}: private-key or credential material present in release tree"
                )
    return failures


def freeze_release(root: str, output: str, policy: str | None = None) -> dict:
    """Freeze only clean, tracked, explicitly allow-listed repository artifacts."""
    project = Path(root).resolve()
    target = Path(output).resolve()
    policy_path = policy or str(project / "configs" / "release_allowlist.json")
    profile, policy_sha256 = _load_profile(policy_path, "repository")
    failures: list[str] = []
    files: list[dict] = []
    try:
        commit = _git(project, "rev-parse", "HEAD").strip()
        dirty = bool(_git(project, "status", "--porcelain").strip())
        tracked = _tracked_entries(project)
    except (subprocess.CalledProcessError, FileNotFoundError):
        commit, dirty, tracked = "not_a_git_repository", True, []
        failures.append("project root is not a readable Git repository")
    for relative, mode in tracked:
        path = project / relative
        if path.resolve() == target:
            continue
        if not path.exists():
            failures.append(f"{relative.as_posix()}: tracked file is missing")
            continue
        file_failures = _audit_file(path, relative, profile)
        failures.extend(file_failures)
        if not file_failures:
            files.append({
                "path": relative.as_posix(),
                "bytes": path.stat().st_size,
                "sha256": sha256(path),
                "mode": mode,
            })
    failures.extend(_scan_local_prohibited(project, profile, target))
    if dirty:
        failures.append("Git worktree is not clean")
    result = {
        "schema_version": "champhla-release-freeze-3",
        "project_root": ".",
        "project_root_locator": Path(os.path.relpath(project, target.parent)).as_posix(),
        "git_commit": commit,
        "git_dirty": dirty,
        "allowlist_sha256": policy_sha256,
        "files": files,
        "file_count": len(files),
        "truth_data_included": False,
        "violations": sorted(set(failures)),
        "valid": not failures,
    }
    write_json(target, result)
    return result


def build_compact_export_manifest(root: str, run_id: str, policy: str,
                                  output: str, list_output: str) -> dict:
    """Select a Roihu run's compact evidence using an explicit path policy."""
    if not SAFE_RUN_ID.fullmatch(run_id):
        raise ValueError("unsafe run_id")
    run_root = Path(root).resolve()
    target = Path(output).resolve()
    listing = Path(list_output).resolve()
    if target.exists() or listing.exists():
        raise FileExistsError("compact export evidence is immutable")
    profile, policy_sha256 = _load_profile(policy, "roihu_compact_export")
    selected: list[dict] = []
    failures: list[str] = []
    for path in sorted(run_root.rglob("*")):
        if not (path.is_file() or path.is_symlink()):
            continue
        relative = path.relative_to(run_root)
        posix = relative.as_posix()
        if not _matches(posix, profile["candidate_patterns"], run_id):
            continue
        file_failures = _audit_file(path, relative, profile, run_id)
        failures.extend(file_failures)
        if not file_failures:
            selected.append({"path": posix, "bytes": path.stat().st_size, "sha256": sha256(path)})
    if not selected:
        failures.append("no compact run evidence matched the explicit allow-list")
    result = {
        "schema_version": "champhla-compact-export-manifest-1",
        "run_id": run_id,
        "root": ".",
        "allowlist_sha256": policy_sha256,
        "files": selected,
        "file_count": len(selected),
        "violations": sorted(set(failures)),
        "valid": not failures,
    }
    write_json(target, result)
    if failures:
        return result
    listing.parent.mkdir(parents=True, exist_ok=True)
    with listing.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write("".join(f"{row['path']}\n" for row in selected))
    return result


def verify_compact_export_manifest(root: str, manifest_path: str) -> list[str]:
    root_path = Path(root).resolve()
    manifest = read_json(manifest_path)
    failures: list[str] = []
    if manifest.get("schema_version") != "champhla-compact-export-manifest-1":
        failures.append("unsupported compact export manifest schema")
    if manifest.get("valid") is not True:
        failures.append("compact export manifest is not valid")
    for row in manifest.get("files", []):
        relative = Path(str(row.get("path", "")))
        if relative.is_absolute() or ".." in relative.parts:
            failures.append(f"unsafe archive path: {relative}")
            continue
        path = root_path / relative
        if not path.is_file() or path.is_symlink():
            failures.append(f"missing or unsafe compact artifact: {relative.as_posix()}")
            continue
        if path.stat().st_size != row.get("bytes"):
            failures.append(f"size drift: {relative.as_posix()}")
        expected = str(row.get("sha256", "")).lower()
        if not HEX64.fullmatch(expected) or sha256(path) != expected:
            failures.append(f"checksum drift: {relative.as_posix()}")
    return sorted(set(failures))


def _safe_archive_name(name: str) -> bool:
    path = Path(name)
    return bool(name) and not path.is_absolute() and ".." not in path.parts and "\\" not in name


def create_release_archive(root: str, freeze_manifest: str, archive: str) -> dict:
    """Create a deterministic archive from a valid release freeze."""
    project = Path(root).resolve()
    manifest_path = Path(freeze_manifest).resolve()
    archive_path = Path(archive).resolve()
    if archive_path.exists():
        raise FileExistsError("release archive is immutable")
    manifest = read_json(manifest_path)
    if manifest.get("schema_version") != "champhla-release-freeze-3" or manifest.get("valid") is not True:
        raise ValueError("release freeze is not valid")
    rows = manifest.get("files", [])
    seen: set[str] = set()
    for row in rows:
        name = str(row.get("path", ""))
        if not _safe_archive_name(name) or name in seen:
            raise ValueError(f"unsafe or duplicate release path: {name!r}")
        seen.add(name)
        path = project / name
        if not path.is_file() or path.is_symlink():
            raise ValueError(f"release file missing or unsafe: {name}")
        if path.stat().st_size != row.get("bytes") or sha256(path) != row.get("sha256"):
            raise ValueError(f"release file drift: {name}")
    archive_path.parent.mkdir(parents=True, exist_ok=True)
    with archive_path.open("xb") as raw:
        with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as compressed:
            with tarfile.open(fileobj=compressed, mode="w", format=tarfile.PAX_FORMAT) as bundle:
                for row in sorted(rows, key=lambda item: item["path"]):
                    name = row["path"]
                    path = project / name
                    info = tarfile.TarInfo(name=name)
                    info.size = path.stat().st_size
                    info.mode = int(row.get("mode", 0o644))
                    info.mtime = 0
                    info.uid = info.gid = 0
                    info.uname = info.gname = ""
                    with path.open("rb") as handle:
                        bundle.addfile(info, handle)
                manifest_bytes = json.dumps(
                    manifest, indent=2, sort_keys=True, ensure_ascii=False,
                ).encode("utf-8") + b"\n"
                info = tarfile.TarInfo(name="RELEASE_MANIFEST.json")
                info.size = len(manifest_bytes)
                info.mode = 0o644
                info.mtime = 0
                info.uid = info.gid = 0
                info.uname = info.gname = ""
                bundle.addfile(info, io.BytesIO(manifest_bytes))
    return {
        "schema_version": "champhla-release-archive-1",
        "archive": archive_path.name,
        "archive_sha256": sha256(archive_path),
        "archive_bytes": archive_path.stat().st_size,
        "release_manifest_sha256": sha256(manifest_path),
        "file_count": len(rows),
    }


def verify_release_archive(archive: str, output: str | None = None) -> dict:
    """Extract into a temporary directory and validate the exact frozen contents."""
    archive_path = Path(archive).resolve()
    failures: list[str] = []
    with tempfile.TemporaryDirectory(prefix="champhla-release-verify-") as temporary:
        extraction = Path(temporary)
        with tarfile.open(archive_path, mode="r:gz") as bundle:
            members = bundle.getmembers()
            names = [member.name for member in members]
            if len(names) != len(set(names)):
                failures.append("archive contains duplicate paths")
            for member in members:
                if not _safe_archive_name(member.name):
                    failures.append(f"unsafe archive path: {member.name!r}")
                if not member.isfile():
                    failures.append(f"non-regular archive member: {member.name!r}")
            if not failures:
                for member in members:
                    source = bundle.extractfile(member)
                    if source is None:
                        failures.append(f"unreadable archive member: {member.name}")
                        continue
                    destination = extraction / member.name
                    destination.parent.mkdir(parents=True, exist_ok=True)
                    with destination.open("wb") as handle:
                        shutil.copyfileobj(source, handle)
        manifest_path = extraction / "RELEASE_MANIFEST.json"
        manifest = read_json(manifest_path) if manifest_path.is_file() else {}
        if manifest.get("schema_version") != "champhla-release-freeze-3":
            failures.append("archive release manifest is missing or unsupported")
        expected = {str(row.get("path", "")): row for row in manifest.get("files", [])}
        observed = {
            path.relative_to(extraction).as_posix()
            for path in extraction.rglob("*") if path.is_file()
        } - {"RELEASE_MANIFEST.json"}
        if observed != set(expected):
            failures.append("archive member set differs from release freeze")
        for name, row in expected.items():
            path = extraction / name
            if not path.is_file():
                failures.append(f"archive member missing: {name}")
                continue
            if path.stat().st_size != row.get("bytes"):
                failures.append(f"archive member size drift: {name}")
            if sha256(path) != row.get("sha256"):
                failures.append(f"archive member checksum drift: {name}")
    result = {
        "schema_version": "champhla-release-archive-verification-1",
        "archive": archive_path.name,
        "archive_sha256": sha256(archive_path),
        "file_count": len(expected) if "expected" in locals() else 0,
        "failures": sorted(set(failures)),
        "passed": not failures,
        "temporary_extraction_removed": True,
    }
    if output:
        output_path = Path(output)
        if output_path.exists():
            raise FileExistsError("release archive verification is immutable")
        write_json(output_path, result)
    return result
