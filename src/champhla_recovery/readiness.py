from __future__ import annotations

import subprocess
from pathlib import Path

from champhla_confirmation.freeze import validate_freeze
from champhla_confirmation.io import read_json
from champhla_confirmation.manifests import validate_comparator_manifest

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
        "signed_amendment", "frozen_comparators", "corrected_three_modality_evaluation",
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
