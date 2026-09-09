import json
from pathlib import Path

from champhla_confirmation.manifests import validate_comparator_manifest
from champhla_recovery.readiness import audit_release_readiness
from champhla_recovery.io import sha256


ROOT = Path(__file__).resolve().parents[1]


def test_comparator_manifest_is_complete_but_not_falsely_frozen():
    path = ROOT / "configs" / "comparator_manifest.json"
    assert validate_comparator_manifest(str(path)) == []
    failures = validate_comparator_manifest(str(path), require_frozen=True)
    assert failures
    assert "comparator manifest status is not FROZEN" in failures


def test_release_audit_fails_closed_on_unfinished_external_work(tmp_path: Path):
    result = audit_release_readiness(
        str(ROOT), str(ROOT / "configs" / "release_requirements.json"),
        str(tmp_path / "readiness.json"),
    )
    assert result["release_ready"] is False
    assert "signed_amendment" in result["blockers"]
    assert "corrected_three_modality_evaluation" in result["blockers"]
    assert result["gates"]["registry_row_validation"]["passed"] is True


def test_generated_development_outputs_match_provenance_manifest():
    manifest = json.loads((
        ROOT / "provenance" / "20260908_development_plurality_evaluation.json"
    ).read_text(encoding="utf-8"))
    for relative, expected in manifest["inputs"].items():
        assert sha256(ROOT / relative) == expected
    for relative, expected in manifest["outputs"].items():
        assert sha256(ROOT / relative) == expected
