import json
from pathlib import Path

from champhla_confirmation.manifests import (
    validate_caller_reference_attestation,
    validate_comparator_manifest,
    validate_hprc_truth_protocol,
)
from champhla_recovery.readiness import audit_goal_completion, audit_release_readiness
from champhla_recovery.io import sha256


ROOT = Path(__file__).resolve().parents[1]


def test_comparator_manifest_is_frozen_and_agrees_with_the_attestation():
    path = ROOT / "configs" / "comparator_manifest.json"
    assert validate_comparator_manifest(str(path)) == []
    assert validate_comparator_manifest(str(path), require_frozen=True) == []


def test_frozen_historical_comparators_have_immutable_surviving_provenance():
    """Archived methods stay visible, but an absent executable must not mean absent evidence."""
    path = ROOT / "configs" / "comparator_manifest.json"
    manifest = json.loads(path.read_text(encoding="utf-8"))
    historical = {
        row["method_id"]: row for row in manifest["comparators"]
        if row["kind"].startswith("historical")
    }
    assert {"WeightedConsensus", "RefFormer", "EvidenceGatedCC"} <= set(historical)
    for row in historical.values():
        assert row["deployable"] is False
        assert len(row["artifact_sha256"]) == 64
        assert len(row["source_sha256"]) == 64
        assert len(row["workflow_hash"]) == 64
        assert row["reference_artifacts"]
        for reference in row["reference_artifacts"]:
            candidate = ROOT / reference["artifact_id"]
            if candidate.is_file():
                assert sha256(candidate) == reference["sha256"]


def test_frozen_manifest_rejects_missing_historical_provenance(tmp_path: Path):
    path = ROOT / "configs" / "comparator_manifest.json"
    manifest = json.loads(path.read_text(encoding="utf-8"))
    row = next(r for r in manifest["comparators"] if r["method_id"] == "RefFormer")
    row["artifact_sha256"] = ""
    row["reference_artifacts"] = []
    forged = tmp_path / "configs" / "comparator_manifest.json"
    forged.parent.mkdir(parents=True)
    forged.write_text(json.dumps(manifest), encoding="utf-8")
    manifest["provenance_attestation"] = str(
        ROOT / "configs" / "roihu_caller_reference_attestation.json"
    )
    forged.write_text(json.dumps(manifest), encoding="utf-8")
    failures = validate_comparator_manifest(str(forged), require_frozen=True)
    assert "RefFormer.artifact_sha256 is not an exact SHA-256" in failures
    assert "RefFormer.reference_artifacts is empty" in failures


def test_frozen_comparator_rows_share_the_exact_current_workflow_lock():
    path = ROOT / "configs" / "comparator_manifest.json"
    manifest = json.loads(path.read_text(encoding="utf-8"))
    reference = manifest["workflow_lock"]
    assert sha256(ROOT / reference["path"]) == reference["sha256"]
    assert {row["workflow_hash"] for row in manifest["comparators"]} == {
        reference["sha256"]
    }


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


def test_protocols_defer_to_caller_specific_imgt_release_evidence():
    confirmation = json.loads((ROOT / "configs" / "confirmation_protocol.json").read_text())
    rerun = json.loads((ROOT / "configs" / "rerun_protocol.json").read_text())
    assert "imgt_hla_version" not in confirmation
    assert confirmation["imgt_hla_release_policy"]["global_release"] is None
    assert rerun["environment"]["imgt_hla_release_policy"]["global_release"] is None
    assert "caller_reference_attestation.json" in confirmation["imgt_hla_release_policy"]["attestation"]


def test_lifecycle_records_reflect_frozen_comparators_and_validated_polysolver_transform():
    prospective = json.loads((ROOT / "configs" / "consensus_evaluation_prospective.json").read_text())
    assert prospective["status"] == "PREDECLARED_PENDING_AUTHOR_SIGNATURE"
    assert prospective["comparator_manifest_status"] == "FROZEN"
    attestation = json.loads((ROOT / "configs" / "roihu_caller_reference_attestation.json").read_text())
    transform = attestation["callers"]["POLYSOLVER"]["wrapper_compatibility_transform"]
    assert transform["status"] == "FROZEN_TECHNICALLY_VALIDATED"
    assert len(transform["technical_pilot_evidence"]) == 2
    assert all(row["hg38_contig_substitutions"] == 3 for row in transform["technical_pilot_evidence"])
    assert all(row["temporary_directory_substitutions"] == 1 for row in transform["technical_pilot_evidence"])
    decision = json.loads((ROOT / "decisions" / "20260914_polysolver_contig_compatibility.json").read_text())
    assert decision["status"] == "TECHNICALLY_VALIDATED"
    assert "not directly comparable" in decision["method_boundary"]["historical_caveat"]


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


def test_bibliographic_dois_are_not_misclassified_as_unregistered_results(tmp_path: Path):
    from champhla_recovery.manuscript import audit_claims

    manuscript = tmp_path / "manuscript.md"
    manuscript.write_text(
        "# Title\n\n## References\n\nExample. 2020;48:D948–D955. "
        "doi:10.1093/nar/gkz950.\n",
        encoding="utf-8",
    )
    result = audit_claims(
        str(manuscript), str(ROOT / "result_registry.tsv"),
        str(ROOT / "manuscripts/claim_audit.tsv"), str(tmp_path / "audit.json"), str(ROOT),
    )
    assert "main numerical result lacks a registry reference" not in result["failures"]


def test_public_validation_gap_repeat_search_preserves_locked_minima():
    decision = json.loads(
        (ROOT / "decisions/20260909_public_validation_gap_assessment.json").read_text(
            encoding="utf-8"
        )
    )
    assert decision["prospective_minima"] == {"wes_donors": 89, "rnaseq_donors": 130}
    repeat = decision["repeat_search"]
    assert repeat["criteria_unchanged"] is True
    assert len(repeat["candidate_dispositions"]) >= 6
    assert repeat["outcome"] == decision["outcome"]
    assert all(row["disposition"] for row in repeat["candidate_dispositions"])

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


def test_full_goal_audit_is_distinct_and_fails_every_unfinished_contract_area(tmp_path: Path):
    project = tmp_path / "empty-project"
    project.mkdir()
    result = audit_goal_completion(
        str(project), str(ROOT / "configs" / "goal_completion_requirements.json"),
        str(tmp_path / "goal-audit.json"),
    )
    expected = {
        "same_resource_benchmark",
        "donor_independent_three_modality",
        "wgs_native_and_named_review",
        "truth_firewall_all_lanes",
        "statistics_recount_registry",
        "identical_commit_verification",
        "issue_ledger_closed",
        "manuscript_claims_declarations_references",
        "figures_reproducible_and_reviewed",
        "markdown_docx_submission_parity",
        "release_archive_safe_and_verified",
        "clean_frozen_git_state",
    }
    assert set(result["gates"]) == expected
    assert set(result["blockers"]) == expected
    assert result["completion_ready"] is False


def test_full_goal_requirements_keep_all_independent_minima():
    config = json.loads((
        ROOT / "configs" / "goal_completion_requirements.json"
    ).read_text(encoding="utf-8"))
    assert {
        modality: row["minimum_donors"]
        for modality, row in config["independent_lanes"].items()
    } == {"wgs": 120, "wes": 89, "rnaseq": 130}
    assert set(config["firewall_lanes"]) == {
        "same_resource", "independent_wgs", "independent_wes", "independent_rnaseq",
    }


def test_generated_development_outputs_match_provenance_manifest():
    manifest = json.loads((
        ROOT / "provenance" / "20260908_development_plurality_evaluation.json"
    ).read_text(encoding="utf-8"))
    for relative, expected in manifest["inputs"].items():
        assert sha256(ROOT / relative) == expected
    for relative, expected in manifest["outputs"].items():
        assert sha256(ROOT / relative) == expected
