from __future__ import annotations

import json
from pathlib import Path

import pytest

from champhla_confirmation.evaluation import evaluate
from champhla_confirmation.io import read_tsv, write_tsv


FIELDS = [
    "cohort", "subject", "donor", "independence_stratum", "modality", "gene",
    "method", "allele1", "allele2", "call_status", "truth_allele1",
    "truth_allele2", "truth_status", "truth_source", "truth_source_sha256",
]


def row(subject: str, modality: str, method: str, predicted: tuple[str, str],
        truth: tuple[str, str], stratum: str = "donor_independent") -> dict[str, str]:
    return {
        "cohort": "C", "subject": subject, "donor": subject,
        "independence_stratum": stratum, "modality": modality, "gene": "A",
        "method": method, "allele1": predicted[0], "allele2": predicted[1],
        "call_status": "callable", "truth_allele1": truth[0],
        "truth_allele2": truth[1], "truth_status": "resolved",
        "truth_source": "orthogonal", "truth_source_sha256": "abc",
    }


def write_design(path: Path, valid: list[str], primary_stratum: str = "") -> None:
    path.write_text(json.dumps({
        "reference_method": "SimplePluralityLex", "valid_modalities": valid,
        "invalid_modalities": {"wgs": "fixture"} if "wgs" not in valid else {},
        "prospective": False, "noninferiority_margin_points": 2,
        "minimum_subjects": {"wgs": 1, "wes": 1, "rnaseq": 1},
        "primary_independence_stratum": primary_stratum,
    }))


def test_invalid_modality_cannot_defeat_valid_modality_gate(tmp_path: Path):
    truth = ("A*01:01", "A*02:01")
    wrong = ("A*03:01", "A*24:02")
    rows = []
    for modality in ("wes", "wgs"):
        for index in range(7):
            plurality = truth if modality == "wes" else wrong
            candidate = truth
            rows.extend([
                row(f"{modality}{index}", modality, "SimplePluralityLex", plurality, truth),
                row(f"{modality}{index}", modality, "Candidate", candidate, truth),
            ])
    joined = tmp_path / "joined.tsv"
    write_tsv(joined, rows, FIELDS)
    design = tmp_path / "design.json"
    write_design(design, ["wes"], "donor_independent")
    result = evaluate(str(joined), str(tmp_path / "eval"), bootstrap=200,
                      allow_partial=True, evaluation_design=str(design))
    assert result["no_comparator_holm_superior"] is True
    assert result["consensus_noninferior_2pp"] is True
    assert result["three_modality_claim_ready"] is False
    comparisons = read_tsv(tmp_path / "eval" / "head_to_head.tsv")
    assert next(r for r in comparisons if r["modality"] == "wgs")[
        "holm_significant_superiority"] == "1"


def test_multiple_independence_strata_are_never_implicitly_pooled(tmp_path: Path):
    truth = ("A*01:01", "A*02:01")
    rows = []
    for index, stratum in enumerate(("donor_independent", "new_library_overlap")):
        rows.extend([
            row(f"S{index}", "wes", "SimplePluralityLex", truth, truth, stratum),
            row(f"S{index}", "wes", "Candidate", truth, truth, stratum),
        ])
    joined = tmp_path / "joined.tsv"
    write_tsv(joined, rows, FIELDS)
    with pytest.raises(ValueError, match="cannot be pooled"):
        evaluate(str(joined), str(tmp_path / "bad"), bootstrap=100, allow_partial=True)
    design = tmp_path / "design.json"
    write_design(design, ["wes"], "donor_independent")
    evaluate(str(joined), str(tmp_path / "ok"), bootstrap=100,
             allow_partial=True, evaluation_design=str(design))
    strata = read_tsv(tmp_path / "ok" / "external_strata_results.tsv")
    assert {r["independence_stratum"] for r in strata} == {
        "donor_independent", "new_library_overlap",
    }
    assert sum(int(r["pooled_into_primary"]) for r in strata) == 1


def test_identical_legacy_alias_is_accepted_but_conflict_fails(tmp_path: Path):
    truth = ("A*01:01", "A*02:01")
    wrong = ("A*03:01", "A*24:02")
    rows = [
        row("S", "wes", "SimplePluralityLex", truth, truth),
        row("S", "wes", "MajorityVote", truth, truth),
        row("S", "wes", "Candidate", truth, truth),
    ]
    joined = tmp_path / "joined.tsv"
    write_tsv(joined, rows, FIELDS)
    evaluate(str(joined), str(tmp_path / "ok"), bootstrap=50, allow_partial=True)
    rows[1] = row("S", "wes", "MajorityVote", wrong, truth)
    write_tsv(joined, rows, FIELDS)
    with pytest.raises(ValueError, match="duplicate prediction"):
        evaluate(str(joined), str(tmp_path / "bad"), bootstrap=50, allow_partial=True)


def test_same_resource_readiness_is_separate_from_independent_claim(tmp_path: Path):
    truth = ("A*01:01", "A*02:01")
    rows = []
    for modality in ("wgs", "wes", "rnaseq"):
        rows.extend([
            row(f"{modality}1", modality, "SimplePluralityLex", truth, truth,
                "new_library_overlap"),
            row(f"{modality}1", modality, "Candidate", truth, truth,
                "new_library_overlap"),
        ])
    joined = tmp_path / "joined.tsv"
    write_tsv(joined, rows, FIELDS)
    design = tmp_path / "same.json"
    design.write_text(json.dumps({
        "reference_method": "SimplePluralityLex",
        "evidence_tier": "same_resource_confirmation",
        "valid_modalities": ["wgs", "wes", "rnaseq"], "invalid_modalities": {},
        "prospective": False, "noninferiority_margin_points": 2,
        "minimum_subjects": {"wgs": 1, "wes": 1, "rnaseq": 1},
        "primary_independence_stratum": "new_library_overlap",
    }))
    audit = tmp_path / "wgs.json"
    audit.write_text('{"passed": true}\n')
    result = evaluate(
        str(joined), str(tmp_path / "eval"), bootstrap=50, mode="external",
        wgs_audit_summary=str(audit), evaluation_design=str(design),
    )
    assert result["same_resource_benchmark_ready"] is True
    assert result["three_modality_claim_ready"] is False
    assert result["status"] == "same_resource_benchmark_complete"
