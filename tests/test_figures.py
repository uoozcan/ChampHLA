import csv
import hashlib
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

from champhla_recovery.figures import REQUIRED_MAIN, validate_figure_manifest


ROOT = Path(__file__).resolve().parents[1]


def test_current_figure_manifest_validates_generated_rows_and_exposes_blockers(tmp_path: Path):
    result = validate_figure_manifest(
        ROOT, ROOT / "manuscripts/figures/figure_manifest.tsv",
        tmp_path / "audit.json",
    )
    assert result["generated_rows_valid"] is True
    assert result["generated_figures"] == 2
    assert set(result["blocked_figures"]) == {"MAIN-3", "MAIN-4", "MAIN-5"}
    assert result["passed"] is False


def test_figure_manifest_has_all_required_main_roles_and_no_fake_blocked_outputs():
    path = ROOT / "manuscripts/figures/figure_manifest.tsv"
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert {row["figure_id"] for row in rows} == REQUIRED_MAIN
    for row in rows:
        if row["status"] == "BLOCKED_EVIDENCE":
            assert row["output_paths"] == ""
            assert row["output_sha256"] == ""
            assert row["blocker"]


def test_figure_validator_detects_source_hash_drift(tmp_path: Path):
    source = ROOT / "manuscripts/figures/figure_manifest.tsv"
    with source.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
        fields = list(rows[0])
    rows[0]["source_data_sha256"] = "0" * 64
    forged = tmp_path / "figure_manifest.tsv"
    with forged.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    result = validate_figure_manifest(ROOT, forged)
    assert result["generated_rows_valid"] is False
    assert any("source_data checksum drift" in failure for failure in result["failures"])


def test_figure_pngs_are_declared_at_publication_resolution():
    result = validate_figure_manifest(ROOT, ROOT / "manuscripts/figures/figure_manifest.tsv")
    assert result["outputs_valid"] is True


def test_static_figure_generation_is_byte_deterministic(tmp_path: Path):
    pytest.importorskip("matplotlib", reason="requires the optional publication extra")
    for relative in ("configs/confirmation_protocol.json", "configs/dataset_roles.json",
                     "manuscripts/figures/scripts/generate_static_figures.py"):
        target = tmp_path / relative
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(ROOT / relative, target)
    command = [sys.executable, str(tmp_path / "manuscripts/figures/scripts/generate_static_figures.py"),
               "--root", str(tmp_path), "--generation-commit", "fixture-commit"]
    subprocess.run(command, check=True)
    outputs = sorted((tmp_path / "manuscripts/figures/main").iterdir())
    first = {path.name: hashlib.sha256(path.read_bytes()).hexdigest() for path in outputs}
    subprocess.run(command, check=True)
    second = {path.name: hashlib.sha256(path.read_bytes()).hexdigest() for path in outputs}
    assert first == second
