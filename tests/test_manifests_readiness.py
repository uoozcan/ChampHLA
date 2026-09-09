import json
from pathlib import Path

from champhla_confirmation.manifests import (
    validate_caller_reference_attestation,
    validate_comparator_manifest,
    validate_hprc_truth_protocol,
)
from champhla_recovery.readiness import audit_release_readiness
from champhla_recovery.io import sha256


ROOT = Path(__file__).resolve().parents[1]


def test_comparator_manifest_is_complete_but_not_falsely_frozen():
    path = ROOT / "configs" / "comparator_manifest.json"
    assert validate_comparator_manifest(str(path)) == []
    failures = validate_comparator_manifest(str(path), require_frozen=True)
    assert failures
    assert "comparator manifest status is not FROZEN" in failures
    assert "Caller:HLA-HD.imgt_hla_version is unresolved" in failures


def test_caller_reference_attestation_exposes_release_blockers():
    path = ROOT / "configs" / "roihu_caller_reference_attestation.json"
    failures = validate_caller_reference_attestation(str(path), require_ready=True)
    assert "HLA-HD: database release unresolved" in failures
    assert "OptiType: database release unresolved" in failures
    assert "POLYSOLVER: database release unresolved" in failures
    assert "SpecHLA: database release unresolved" in failures
    assert validate_caller_reference_attestation(str(path), require_ready=False) == []


def test_hprc_execution_remains_blocked_until_tools_and_archive_are_pinned():
    path = ROOT / "configs" / "hprc_truth_protocol.json"
    assert validate_hprc_truth_protocol(str(path), require_frozen=False) == []
    failures = validate_hprc_truth_protocol(str(path), require_frozen=True)
    assert "HPRC AGC archive SHA-256 is unresolved" in failures
    assert "HLA-ASM.artifact_sha256 is not an exact SHA-256" in failures
    assert "Immuannot.wrapper_sha256 is not an exact SHA-256" in failures
    assert "HPRC truth protocol is not FROZEN" in failures


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
