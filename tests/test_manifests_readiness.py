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


def test_comparator_manifest_is_frozen_and_agrees_with_the_attestation():
    path = ROOT / "configs" / "comparator_manifest.json"
    assert validate_comparator_manifest(str(path)) == []
    assert validate_comparator_manifest(str(path), require_frozen=True) == []


def test_frozen_manifest_rejects_a_release_the_attestation_does_not_support(tmp_path: Path):
    """The manifest and the attestation are written by different steps; they must not drift."""
    path = ROOT / "configs" / "comparator_manifest.json"
    manifest = json.loads(path.read_text(encoding="utf-8"))
    for row in manifest["comparators"]:
        if row["method_id"] == "Caller:HLA-HD":
            row["imgt_hla_version"] = "3.32.0"
    drifted = tmp_path / "configs" / "comparator_manifest.json"
    drifted.parent.mkdir(parents=True, exist_ok=True)
    drifted.write_text(json.dumps(manifest), encoding="utf-8")
    manifest["provenance_attestation"] = str(
        ROOT / "configs" / "roihu_caller_reference_attestation.json"
    )
    drifted.write_text(json.dumps(manifest), encoding="utf-8")
    failures = validate_comparator_manifest(str(drifted), require_frozen=True)
    assert any("disagrees with attestation evidence" in failure for failure in failures)


def test_caller_reference_attestation_resolves_every_caller():
    path = ROOT / "configs" / "roihu_caller_reference_attestation.json"
    assert validate_caller_reference_attestation(str(path), require_ready=True) == []
    assert validate_caller_reference_attestation(str(path), require_ready=False) == []


def test_heterogeneous_database_is_resolved_not_blocked():
    """SpecHLA mixes 3.38.0 and 3.51.0 components. That is a property, not a failure."""
    path = ROOT / "configs" / "roihu_caller_reference_attestation.json"
    data = json.loads(path.read_text(encoding="utf-8"))
    spechla = data["callers"]["SpecHLA"]
    assert spechla["database_release"] == "HETEROGENEOUS"
    assert spechla["release_resolved"] is True
    releases = {
        component["release"] for component in spechla["reference_components"]
    }
    assert {"3.38.0", "3.51.0"} <= releases


def test_release_resolved_cannot_be_asserted_without_component_evidence(tmp_path: Path):
    """Resolution is derived from evidence, so flipping the boolean must not pass."""
    path = ROOT / "configs" / "roihu_caller_reference_attestation.json"
    data = json.loads(path.read_text(encoding="utf-8"))
    for component in data["callers"]["OptiType"]["reference_components"]:
        component["release"] = "NO_MATCH"
    data["callers"]["OptiType"]["release_resolved"] = True
    forged = tmp_path / "attestation.json"
    forged.write_text(json.dumps(data), encoding="utf-8")
    failures = validate_caller_reference_attestation(str(forged), require_ready=True)
    assert any("release_resolved does not match" in failure for failure in failures)


def test_non_release_bearing_components_do_not_decide_a_caller_release():
    """An exon-split table carries no alleles and must not vote on the release."""
    path = ROOT / "configs" / "roihu_caller_reference_attestation.json"
    data = json.loads(path.read_text(encoding="utf-8"))
    hlahd = data["callers"]["HLA-HD"]
    split = [
        component for component in hlahd["reference_components"]
        if component["artifact_id"].startswith("HLA_gene.split")
    ]
    assert split and all(c["release"] == "NOT_RELEASE_BEARING" for c in split)
    assert hlahd["database_release"] == "3.15.0"


def test_claim_audit_gates_point_at_distinct_artifacts():
    """readiness treats these as two gates; aliasing them checks one thing twice."""
    requirements = json.loads(
        (ROOT / "configs" / "release_requirements.json").read_text(encoding="utf-8")
    )
    paths = requirements["paths"]
    assert paths["main_claim_audit"] != paths["supplement_claim_audit"]
    for left, right in requirements["distinct_path_pairs"]:
        assert paths[left] != paths[right]


def test_supplement_may_cite_an_invalid_result_only_when_labelled_diagnostic(tmp_path: Path):
    """The supplement rules exist so the WGS diagnostic can be cited without being quotable."""
    from champhla_recovery.manuscript import audit_claims

    registry = ROOT / "result_registry.tsv"
    claims = ROOT / "manuscripts" / "claim_audit.tsv"
    main = ROOT / "manuscripts" / "benchmark" / "manuscript.md"

    labelled = tmp_path / "labelled.md"
    labelled.write_text(
        "## S\n\nThe legacy lane is an invalid diagnostic and enters no gate "
        "[RESULT:DEV_WGS_CONSENSUS].\n",
        encoding="utf-8",
    )
    result = audit_claims(str(main), str(registry), str(claims),
                          str(tmp_path / "ok.json"), str(ROOT), str(labelled))
    assert result["failures"] == []

    unlabelled = tmp_path / "unlabelled.md"
    unlabelled.write_text(
        "## S\n\nWhole-genome plurality reached 164 of 411 loci "
        "[RESULT:DEV_WGS_CONSENSUS].\n",
        encoding="utf-8",
    )
    result = audit_claims(str(main), str(registry), str(claims),
                          str(tmp_path / "bad.json"), str(ROOT), str(unlabelled))
    assert any("without diagnostic labeling" in failure for failure in result["failures"])


def test_abstaining_rule_is_excluded_from_the_always_call_family():
    """Holm runs within the always-call family; the abstaining rule is reported separately."""
    manifest = json.loads(
        (ROOT / "configs" / "comparator_manifest.json").read_text(encoding="utf-8")
    )
    entry = next(
        row for row in manifest["comparators"]
        if row["method_id"] == "SimpleTwoThirdsConsensus"
    )
    assert entry["primary_family"] is False
    design = json.loads(
        (ROOT / "configs" / "consensus_evaluation.json").read_text(encoding="utf-8")
    )
    assert "SimpleTwoThirdsConsensus" in design["abstaining_comparators"]
    assert design["family_policy"]


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
