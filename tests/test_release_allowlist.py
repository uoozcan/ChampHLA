import json
import io
import subprocess
import tarfile
from pathlib import Path

import pytest

from champhla_recovery.release import (
    build_compact_export_manifest,
    create_release_archive,
    freeze_release,
    verify_compact_export_manifest,
    verify_release_archive,
)


ROOT = Path(__file__).resolve().parents[1]
POLICY = ROOT / "configs" / "release_allowlist.json"


def _git(root: Path, *args: str) -> None:
    subprocess.run(["git", "-C", str(root), *args], check=True, capture_output=True)


def _repository(tmp_path: Path) -> Path:
    root = tmp_path / "repository"
    (root / "configs").mkdir(parents=True)
    (root / "src").mkdir()
    (root / "configs" / "release_allowlist.json").write_bytes(POLICY.read_bytes())
    (root / "src" / "method.py").write_text("VALUE = 1\n", encoding="utf-8")
    (root / ".gitignore").write_text("*.bam\n.venv/\n", encoding="utf-8")
    _git(root, "init", "-q")
    _git(root, "config", "user.name", "Test")
    _git(root, "config", "user.email", "test@example.invalid")
    _git(root, "add", ".")
    _git(root, "commit", "-qm", "fixture")
    return root


def test_repository_freeze_uses_clean_tracked_allowlist(tmp_path: Path):
    root = _repository(tmp_path)
    result = freeze_release(str(root), str(tmp_path / "release.json"))
    assert result["valid"] is True
    assert result["git_dirty"] is False
    assert {row["path"] for row in result["files"]} == {
        ".gitignore", "configs/release_allowlist.json", "src/method.py",
    }
    assert {row["mode"] for row in result["files"]} == {0o644}


def test_release_archive_is_deterministic_and_verifies_after_extraction(tmp_path: Path):
    root = _repository(tmp_path)
    freeze_path = tmp_path / "release.json"
    result = freeze_release(str(root), str(freeze_path))
    assert result["valid"] is True
    first = create_release_archive(str(root), str(freeze_path), str(tmp_path / "first.tar.gz"))
    second = create_release_archive(str(root), str(freeze_path), str(tmp_path / "second.tar.gz"))
    assert first["archive_sha256"] == second["archive_sha256"]
    verification = verify_release_archive(str(tmp_path / "first.tar.gz"))
    assert verification["passed"] is True
    assert verification["file_count"] == result["file_count"]


def test_release_archive_verifier_rejects_path_traversal(tmp_path: Path):
    archive = tmp_path / "malicious.tar.gz"
    with tarfile.open(archive, "w:gz") as bundle:
        payload = b"escape"
        info = tarfile.TarInfo("../escape.txt")
        info.size = len(payload)
        bundle.addfile(info, io.BytesIO(payload))
    result = verify_release_archive(str(archive))
    assert result["passed"] is False
    assert any("unsafe archive path" in failure for failure in result["failures"])


@pytest.mark.parametrize(
    "relative",
    [
        "ignored.bam",
        "runtime_data/sample.fastq.gz",
        "roihu_containers/caller.sif",
        "assemblies/cohort.agc",
        "references/GRCh38.fa",
        "secrets/id_ed25519",
    ],
)
def test_repository_freeze_detects_untracked_or_ignored_prohibited_payloads(
        tmp_path: Path, relative: str):
    root = _repository(tmp_path)
    target = root / relative
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text("prohibited\n", encoding="utf-8")
    result = freeze_release(str(root), str(tmp_path / "release.json"))
    assert result["valid"] is False
    assert any(relative in failure for failure in result["violations"])


def test_repository_freeze_detects_private_key_content(tmp_path: Path):
    root = _repository(tmp_path)
    secret = root / "artifacts" / "oops.txt"
    secret.parent.mkdir()
    secret.write_text("-----BEGIN OPENSSH PRIVATE KEY-----\n", encoding="utf-8")
    result = freeze_release(str(root), str(tmp_path / "release.json"))
    assert result["valid"] is False
    assert any("credential material" in failure for failure in result["violations"])


def test_compact_export_selects_only_run_scoped_reviewable_evidence(tmp_path: Path):
    (tmp_path / "manifests").mkdir()
    (tmp_path / "logs").mkdir()
    marker = tmp_path / "caller_outputs" / "run-1" / "technical_pilot" / "wgs" / "S1"
    marker.mkdir(parents=True)
    (tmp_path / "manifests" / "run-1_wgs.ledger.tsv").write_text("state\nvalidated\n")
    (tmp_path / "logs" / "stage_run-1_1.out").write_text("ok\n")
    (marker / "CALLERS_COMPLETE").write_text("")
    output = tmp_path / "exports" / "manifest.json"
    listing = tmp_path / "exports" / "files.txt"
    result = build_compact_export_manifest(
        str(tmp_path), "run-1", str(POLICY), str(output), str(listing),
    )
    assert result["valid"] is True
    assert result["file_count"] == 3
    assert verify_compact_export_manifest(str(tmp_path), str(output)) == []
    assert listing.read_text().splitlines() == sorted(listing.read_text().splitlines())


@pytest.mark.parametrize(
    "relative",
    [
        "predictions/run-1/sample.bam",
        "logs/run-1_secret.txt",
        "evaluation/run-1/truth_join.json",
        "review/run-1/unreviewed.txt",
    ],
)
def test_compact_export_fails_closed_on_prohibited_or_unreviewed_candidate(
        tmp_path: Path, relative: str):
    allowed = tmp_path / "manifests" / "run-1_wgs.ledger.tsv"
    allowed.parent.mkdir(parents=True)
    allowed.write_text("state\nvalidated\n")
    candidate = tmp_path / relative
    candidate.parent.mkdir(parents=True, exist_ok=True)
    candidate.write_text("x\n")
    output = tmp_path / "exports" / "manifest.json"
    listing = tmp_path / "exports" / "files.txt"
    result = build_compact_export_manifest(
        str(tmp_path), "run-1", str(POLICY), str(output), str(listing),
    )
    assert result["valid"] is False
    assert not listing.exists()
    assert any(relative in failure for failure in result["violations"])


def test_compact_manifest_detects_checksum_drift_and_is_immutable(tmp_path: Path):
    allowed = tmp_path / "manifests" / "run-1_wgs.ledger.tsv"
    allowed.parent.mkdir(parents=True)
    allowed.write_text("state\nvalidated\n")
    output = tmp_path / "exports" / "manifest.json"
    listing = tmp_path / "exports" / "files.txt"
    build_compact_export_manifest(
        str(tmp_path), "run-1", str(POLICY), str(output), str(listing),
    )
    allowed.write_text("state\nchanged\n")
    assert any("drift" in failure for failure in verify_compact_export_manifest(tmp_path, output))
    with pytest.raises(FileExistsError):
        build_compact_export_manifest(
            str(tmp_path), "run-1", str(POLICY), str(output), str(listing),
        )


def test_compact_export_rejects_unsafe_run_id(tmp_path: Path):
    with pytest.raises(ValueError, match="unsafe run_id"):
        build_compact_export_manifest(
            str(tmp_path), "../escape", str(POLICY),
            str(tmp_path / "manifest.json"), str(tmp_path / "files.txt"),
        )


def test_roihu_export_script_uses_shared_manifest_not_extension_only_find():
    script = (ROOT / "scripts" / "roihu_export_compact.sh").read_text(encoding="utf-8")
    assert "build_compact_export_manifest_main" in script
    assert "verify_compact_export_manifest_main" in script
    assert "find \"${run_root}\"" not in script
