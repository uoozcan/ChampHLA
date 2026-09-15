from __future__ import annotations

import subprocess
import csv
from pathlib import Path

from champhla_confirmation.freeze import validate_freeze
from champhla_confirmation.io import read_json
from champhla_confirmation.manifests import (
    validate_caller_reference_attestation,
    validate_comparator_manifest,
)
from champhla_confirmation.roihu import validate_workflow_lock

from .registry import validate_registry
from .io import write_json


def _json(root: Path, relative: str) -> dict:
    path = root / relative
    return read_json(path) if path.is_file() else {}


def audit_release_readiness(root: str, config_path: str, output: str) -> dict:
    """Evaluate benchmark-release and independent-claim gates separately."""
    project = Path(root).resolve()
    config = read_json(config_path)
    paths = config["paths"]
    gates: dict[str, dict] = {}

    amendment = _json(project, paths["amendment"])
    gates["signed_amendment"] = {
        "passed": amendment.get("status") == "SIGNED_BY_AUTHOR"
        and bool(amendment.get("signed_by")) and bool(amendment.get("signed_at_utc")),
        "evidence": paths["amendment"],
    }

    comparator_path = project / paths["comparator_manifest"]
    comparator_failures = (validate_comparator_manifest(str(comparator_path), True)
                           if comparator_path.is_file() else ["manifest missing"])
    gates["frozen_comparators"] = {
        "passed": not comparator_failures,
        "evidence": paths["comparator_manifest"],
        "failures": comparator_failures,
    }

    caller_reference_path = project / paths["caller_reference_attestation"]
    caller_reference_failures = (
        validate_caller_reference_attestation(str(caller_reference_path), True)
        if caller_reference_path.is_file() else ["attestation missing"]
    )
    gates["caller_reference_provenance"] = {
        "passed": not caller_reference_failures,
        "evidence": paths["caller_reference_attestation"],
        "failures": caller_reference_failures,
    }

    workflow_path = project / paths["workflow_lock"]
    workflow = (validate_workflow_lock(workflow_path, project)
                if workflow_path.is_file() else {"passed": False, "failures": ["lock missing"]})
    gates["frozen_workflow"] = {
        "passed": bool(workflow["passed"]), "evidence": paths["workflow_lock"],
        "failures": workflow["failures"],
    }

    run_manifest = _json(project, paths["same_resource_manifest_freeze"])
    gates["frozen_same_resource_manifest"] = {
        "passed": bool(run_manifest.get("frozen")) and bool(run_manifest.get("passed"))
        and run_manifest.get("expected_samples_by_modality") == {"wgs": 137, "wes": 130, "rnaseq": 107}
        and run_manifest.get("expected_caller_locus_records") == 5289
        and run_manifest.get("expected_plurality_rows") == 1122,
        "evidence": paths["same_resource_manifest_freeze"],
    }

    nci60 = _json(project, paths["nci60_manifest_freeze"])
    gates["frozen_nci60_exploratory_manifest"] = {
        "passed": bool(nci60.get("frozen")) and nci60.get("expected_samples_by_modality") == {"rnaseq": 11},
        "evidence": paths["nci60_manifest_freeze"],
        "required_for_benchmark_release": False,
    }

    evaluation = _json(project, paths["corrected_evaluation"])
    gates["corrected_three_modality_evaluation"] = {
        "passed": bool(evaluation.get("same_resource_benchmark_ready")),
        "evidence": paths["corrected_evaluation"],
    }

    wgs = _json(project, paths["wgs_audit"])
    gates["wgs_technical_and_human_review"] = {
        "passed": bool(wgs.get("passed")) and bool(wgs.get("human_review_passed"))
        and bool(wgs.get("named_reviewer")),
        "evidence": paths["wgs_audit"],
    }

    freeze_path = project / paths["prediction_freeze"]
    freeze_result = validate_freeze(freeze_path) if freeze_path.is_file() else {
        "valid": False, "failures": ["freeze missing"]
    }
    gates["truth_blind_prediction_freeze"] = {
        "passed": bool(freeze_result["valid"]),
        "evidence": paths["prediction_freeze"],
        "failures": freeze_result["failures"],
    }

    join = _json(project, paths["truth_join"])
    gates["one_time_truth_join"] = {
        "passed": bool(join.get("join_once")) and bool(join.get("prediction_sha256"))
        and bool(join.get("truth_sha256")) and bool(join.get("joined_sha256")),
        "evidence": paths["truth_join"],
    }

    registry_path = project / paths["registry"]
    registry_failures = (validate_registry(str(registry_path), str(project))
                         if registry_path.is_file() else ["registry missing"])
    gates["registry_row_validation"] = {
        "passed": not registry_failures,
        "evidence": paths["registry"],
        "failures": registry_failures,
    }

    reconciliation = _json(project, paths["count_reconciliation"])
    gates["count_reconciliation"] = {
        "passed": reconciliation.get("status") == "RESOLVED"
        and not reconciliation.get("unexplained_discrepancies", ["not reported"]),
        "evidence": paths["count_reconciliation"],
    }

    for name in ("main_claim_audit", "supplement_claim_audit"):
        audit = _json(project, paths[name])
        gates[name] = {
            "passed": bool(audit.get("submission_ready")),
            "evidence": paths[name],
        }

    required_tests = set(config.get(
        "benchmark_required_test_attestations", config.get("test_attestations", {})
    ))
    for environment, relative in sorted(config.get("test_attestations", {}).items()):
        attestation = _json(project, relative)
        gates[f"tests_{environment}"] = {
            "passed": bool(attestation.get("passed")), "evidence": relative,
            "required_for_benchmark_release": environment in required_tests,
        }

    try:
        dirty = bool(subprocess.run(
            ["git", "-C", str(project), "status", "--porcelain"], check=True,
            capture_output=True, text=True,
        ).stdout.strip())
    except (FileNotFoundError, subprocess.CalledProcessError):
        dirty = True
    gates["clean_git_state"] = {"passed": not dirty, "evidence": "git status --porcelain"}

    benchmark_gate_names = {
        "signed_amendment", "frozen_comparators", "caller_reference_provenance",
        "frozen_workflow", "frozen_same_resource_manifest", "corrected_three_modality_evaluation",
        "wgs_technical_and_human_review", "truth_blind_prediction_freeze",
        "one_time_truth_join", "registry_row_validation", "count_reconciliation",
        "main_claim_audit", "supplement_claim_audit", "clean_git_state",
        *{f"tests_{name}" for name in required_tests},
    }
    blockers = sorted(
        name for name in benchmark_gate_names if not gates.get(name, {}).get("passed", False)
    )
    independent_path = paths.get("independent_evaluation", "")
    independent = _json(project, independent_path) if independent_path else {}
    independent_ready = bool(independent.get("three_modality_claim_ready"))
    independent_blockers = [] if independent_ready else [
        "donor-independent prospective WGS/WES/RNA minima and noninferiority gate not complete"
    ]
    result = {
        "schema_version": "champhla-release-readiness-2",
        "release_ready": not blockers,
        "same_resource_benchmark_ready": not blockers,
        "independent_three_modality_claim_ready": independent_ready,
        "gates": gates,
        "blockers": blockers,
        "independent_claim_evidence": independent_path,
        "independent_claim_blockers": independent_blockers,
        "policy": "independent validation is reported separately and cannot block an honest same-resource benchmark release",
    }
    write_json(output, result)
    return result


def _passed_json(project: Path, relative: str, field: str = "passed") -> tuple[bool, dict]:
    data = _json(project, relative)
    return bool(data.get(field)), data


def _read_issue_rows(path: Path) -> list[dict[str, str]]:
    if not path.is_file():
        return []
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def audit_goal_completion(root: str, config_path: str, output: str) -> dict:
    """Audit the full research goal without weakening benchmark-release semantics."""
    project = Path(root).resolve()
    config = read_json(config_path)
    if config.get("schema_version") != "champhla-goal-completion-requirements-1":
        raise ValueError("unsupported goal-completion requirements schema")
    paths = config["paths"]
    gates: dict[str, dict] = {}

    benchmark_ok, benchmark = _passed_json(
        project, paths["benchmark_readiness"], "release_ready",
    )
    gates["same_resource_benchmark"] = {
        "passed": benchmark_ok and bool(benchmark.get("same_resource_benchmark_ready")),
        "evidence": paths["benchmark_readiness"],
    }

    independent_ok, independent = _passed_json(
        project, paths["independent_evaluation"], "independent_three_modality_claim_ready",
    )
    lane_failures: list[str] = []
    for modality, spec in config["independent_lanes"].items():
        lane = _json(project, spec["path"])
        if not lane.get("valid"):
            lane_failures.append(f"{modality}: validation is not valid")
        if lane.get("evidence_role") != "independent_validation":
            lane_failures.append(f"{modality}: evidence role is not independent_validation")
        if int(lane.get("eligible_donors", 0)) < int(spec["minimum_donors"]):
            lane_failures.append(f"{modality}: donor minimum not met")
        if not lane.get("truth_frozen") or not lane.get("prediction_frozen"):
            lane_failures.append(f"{modality}: truth/prediction freeze is incomplete")
    gates["donor_independent_three_modality"] = {
        "passed": independent_ok and not lane_failures,
        "evidence": paths["independent_evaluation"],
        "failures": lane_failures,
    }

    wgs = _json(project, paths["wgs_audit"])
    gates["wgs_native_and_named_review"] = {
        "passed": bool(wgs.get("passed")) and bool(wgs.get("human_review_passed"))
        and bool(wgs.get("named_reviewer")) and bool(wgs.get("reviewed_at_utc")),
        "evidence": paths["wgs_audit"],
    }

    firewall_failures: list[str] = []
    for lane, spec in config["firewall_lanes"].items():
        freeze_path = project / spec["freeze"]
        if not freeze_path.is_file():
            firewall_failures.append(f"{lane}: prediction freeze missing")
        else:
            result = validate_freeze(freeze_path)
            firewall_failures.extend(f"{lane}: {failure}" for failure in result["failures"])
        join = _json(project, spec["join"])
        if not (join.get("join_once") and join.get("prediction_sha256")
                and join.get("truth_sha256") and join.get("joined_sha256")
                and join.get("non_overwriting") is True):
            firewall_failures.append(f"{lane}: one-time non-overwriting join is invalid")
    gates["truth_firewall_all_lanes"] = {
        "passed": not firewall_failures,
        "evidence": config["firewall_lanes"],
        "failures": firewall_failures,
    }

    statistics = _json(project, paths["statistical_audit"])
    gates["statistics_recount_registry"] = {
        "passed": bool(statistics.get("passed"))
        and bool(statistics.get("independent_recount_agrees"))
        and bool(statistics.get("registry_reconciled"))
        and bool(statistics.get("multiplicity_valid"))
        and bool(statistics.get("donor_clustering_valid")),
        "evidence": paths["statistical_audit"],
    }

    try:
        head = subprocess.run(
            ["git", "-C", str(project), "rev-parse", "HEAD"], check=True,
            capture_output=True, text=True,
        ).stdout.strip()
        tree = subprocess.run(
            ["git", "-C", str(project), "rev-parse", "HEAD^{tree}"], check=True,
            capture_output=True, text=True,
        ).stdout.strip()
    except (FileNotFoundError, subprocess.CalledProcessError):
        head, tree = "", ""
    test_failures: list[str] = []
    for environment, relative in config["test_attestations"].items():
        attestation = _json(project, relative)
        if not attestation.get("passed"):
            test_failures.append(f"{environment}: tests not passed")
        if attestation.get("commit") != head or attestation.get("tree") != tree:
            test_failures.append(f"{environment}: commit/tree differs from release checkout")
        if not attestation.get("pytest_entrypoint", {}).get("passed"):
            test_failures.append(f"{environment}: pytest entrypoint did not pass")
        if not attestation.get("python_module_entrypoint", {}).get("passed"):
            test_failures.append(f"{environment}: python -m pytest entrypoint did not pass")
    gates["identical_commit_verification"] = {
        "passed": not test_failures,
        "evidence": config["test_attestations"],
        "failures": test_failures,
    }

    issues = _read_issue_rows(project / paths["issue_ledger"])
    unresolved = [row.get("issue_id", "") for row in issues
                  if row.get("status") not in {"FIXED_VERIFIED", "SUPERSEDED"}]
    gates["issue_ledger_closed"] = {
        "passed": bool(issues) and not unresolved,
        "evidence": paths["issue_ledger"],
        "unresolved": unresolved,
    }

    main = _json(project, paths["main_claim_audit"])
    combined = _json(project, paths["combined_claim_audit"])
    declarations = _json(project, paths["declarations_audit"])
    bibliography = _json(project, paths["bibliography_audit"])
    gates["manuscript_claims_declarations_references"] = {
        "passed": bool(main.get("submission_ready"))
        and bool(combined.get("submission_ready"))
        and bool(declarations.get("passed")) and bool(bibliography.get("passed")),
        "evidence": [paths["main_claim_audit"], paths["combined_claim_audit"],
                     paths["declarations_audit"], paths["bibliography_audit"]],
    }

    figure = _json(project, paths["figure_audit"])
    figure_review = _json(project, paths["figure_review"])
    gates["figures_reproducible_and_reviewed"] = {
        "passed": bool(figure.get("passed")) and bool(figure.get("registry_agreement"))
        and bool(figure.get("source_hashes_valid")) and bool(figure.get("outputs_valid"))
        and bool(figure_review.get("passed")) and bool(figure_review.get("named_reviewer"))
        and bool(figure_review.get("reviewed_at_utc")),
        "evidence": [paths["figure_audit"], paths["figure_review"]],
    }

    submission = _json(project, paths["submission_audit"])
    submission_review = _json(project, paths["submission_review"])
    gates["markdown_docx_submission_parity"] = {
        "passed": bool(submission.get("passed")) and bool(submission.get("text_parity"))
        and bool(submission.get("structure_valid")) and bool(submission.get("figure_order_valid"))
        and bool(submission_review.get("passed")) and bool(submission_review.get("named_reviewer"))
        and bool(submission_review.get("reviewed_at_utc")),
        "evidence": [paths["submission_audit"], paths["submission_review"]],
    }

    release_freeze = _json(project, paths["release_freeze"])
    archive = _json(project, paths["archive_verification"])
    gates["release_archive_safe_and_verified"] = {
        "passed": benchmark_ok and bool(release_freeze.get("valid"))
        and not release_freeze.get("violations") and bool(archive.get("passed"))
        and bool(archive.get("temporary_extraction_removed")),
        "evidence": [paths["release_freeze"], paths["archive_verification"]],
    }

    try:
        dirty = bool(subprocess.run(
            ["git", "-C", str(project), "status", "--porcelain"], check=True,
            capture_output=True, text=True,
        ).stdout.strip())
    except (FileNotFoundError, subprocess.CalledProcessError):
        dirty = True
    gates["clean_frozen_git_state"] = {
        "passed": bool(head) and bool(tree) and not dirty,
        "evidence": "git HEAD, HEAD^{tree}, status --porcelain",
        "commit": head,
        "tree": tree,
    }

    blockers = sorted(name for name, gate in gates.items() if not gate["passed"])
    result = {
        "schema_version": "champhla-goal-completion-audit-1",
        "completion_ready": not blockers,
        "gates": gates,
        "blockers": blockers,
        "policy": "Full completion requires all 12 contract clauses; benchmark readiness alone is insufficient.",
    }
    write_json(output, result)
    return result
