from __future__ import annotations

from pathlib import Path

from champhla_confirmation.consensus import build_consensus
from champhla_confirmation.cli import run_truth_blind_predictions_main
from champhla_confirmation.evaluation import evaluate, join_truth
from champhla_confirmation.freeze import freeze_bundle
from champhla_confirmation.io import read_tsv, sha256, write_tsv
from champhla_confirmation.schema import normalize_caller_row
from champhla_recovery.manuscript import audit_claims
from champhla_recovery.registry import validate_registry


def test_prediction_cli_runs_without_cc_inputs(tmp_path: Path, monkeypatch):
    pairs = {
        "A": ("A*01:01", "A*02:01"),
        "B": ("B*07:02", "B*08:01"),
        "C": ("C*07:01", "C*07:02"),
    }
    rows = []
    for gene, pair in pairs.items():
        for caller in ("HLA-HD", "OptiType", "POLYSOLVER", "SpecHLA", "T1K"):
            rows.append({
                "cohort": "E", "subject": "S1", "modality": "wes", "gene": gene,
                "caller": caller, "allele1": pair[0], "allele2": pair[1],
                "call_status": "callable",
            })
    calls = tmp_path / "calls.tsv"
    predictions = tmp_path / "predictions.tsv"
    manifest = tmp_path / "manifest.json"
    write_tsv(calls, rows)
    monkeypatch.setattr("sys.argv", [
        "run_truth_blind_predictions", "--calls", str(calls),
        "--predictions", str(predictions), "--manifest", str(manifest),
    ])
    assert run_truth_blind_predictions_main() == 0
    methods = {row["method"] for row in read_tsv(predictions)}
    assert "SimplePluralityLex" in methods
    assert "TwoThirdsGuardedCC" not in methods


def test_caller_to_manuscript_evidence_chain(tmp_path: Path):
    truth_pair = ("A*01:01", "A*02:01")
    calls = [
        normalize_caller_row({
            "cohort": "E", "subject": "S1", "modality": "wes", "gene": "A",
            "caller": caller, "allele1": truth_pair[0], "allele2": truth_pair[1],
            "call_status": "callable", "source_sha256": f"hash-{caller}",
        })
        for caller in ("HLA-HD", "OptiType")
    ]
    predictions = build_consensus(calls)
    prediction_path = tmp_path / "predictions.tsv"
    write_tsv(prediction_path, predictions)
    protocol = tmp_path / "protocol.json"
    protocol.write_text("{}\n", encoding="utf-8")
    freeze = tmp_path / "freeze.json"
    freeze_bundle(tmp_path, {"predictions": prediction_path, "protocol": protocol}, freeze)

    truth = tmp_path / "truth.tsv"
    write_tsv(truth, [{
        "cohort": "E", "subject": "S1", "modality": "wes", "gene": "A",
        "truth_allele1": truth_pair[0], "truth_allele2": truth_pair[1],
        "truth_status": "resolved", "truth_source": "orthogonal",
        "source_sha256": "truth-hash",
    }])
    joined = tmp_path / "joined.tsv"
    join_manifest = tmp_path / "join.json"
    join_truth(str(prediction_path), str(truth), str(freeze), str(joined), str(join_manifest))
    evaluation_dir = tmp_path / "evaluation"
    evaluate(str(joined), str(evaluation_dir), bootstrap=50, allow_partial=True)
    result = read_tsv(evaluation_dir / "primary_modality_results.tsv")[0]

    canonical = tmp_path / "results.tsv"
    write_tsv(canonical, [{
        "result_id": "PIPELINE", "cohort": "E", "modality": result["modality"],
        "method": result["reference_method"], "evidence_role": "development",
        "analysis_status": "discovery", "validity": "valid", "subjects": result["subjects"],
        "loci": result["loci"], "correct": result["correct"], "accuracy": result["accuracy"],
    }])
    registry = tmp_path / "registry.tsv"
    write_tsv(registry, [{
        "result_id": "PIPELINE", "cohort": "E", "modality": result["modality"],
        "method": result["reference_method"], "evidence_role": "development",
        "analysis_status": "discovery", "validity": "valid", "subjects": result["subjects"],
        "loci": result["loci"], "correct": result["correct"], "accuracy": result["accuracy"],
        "artifact_path": canonical.name, "artifact_sha256": sha256(canonical),
        "artifact_record_key": "result_id=PIPELINE",
        "source_artifact_path": str(evaluation_dir.relative_to(tmp_path) / "primary_modality_results.tsv"),
        "source_artifact_sha256": sha256(evaluation_dir / "primary_modality_results.tsv"),
        "abstract_allowed": "0", "claim_boundary": "integration fixture",
    }])
    assert validate_registry(str(registry), str(tmp_path)) == []

    manuscript_dir = tmp_path / "benchmark"
    manuscript_dir.mkdir()
    manuscript = manuscript_dir / "manuscript.md"
    manuscript.write_text("# Results\nPipeline result [RESULT:PIPELINE].\n", encoding="utf-8")
    claims = tmp_path / "claims.tsv"
    write_tsv(claims, [{
        "claim_id": "P1", "draft": "benchmark", "section": "Results",
        "claim": "pipeline", "status": "supported", "evidence_result_id": "PIPELINE",
        "action": "retain",
    }])
    audit = audit_claims(
        str(manuscript), str(registry), str(claims), str(tmp_path / "audit.json"), str(tmp_path),
    )
    assert audit["submission_ready"] is True
