from pathlib import Path

from champhla_recovery.io import sha256
from champhla_recovery.manuscript import audit_claims
from champhla_recovery.registry import validate_registry


def _registry(tmp_path: Path) -> Path:
    artifact = tmp_path / "artifact.tsv"
    artifact.write_text("x\n1\n", encoding="utf-8")
    registry = tmp_path / "registry.tsv"
    registry.write_text(
        "result_id\tcohort\tmodality\tmethod\tanalysis_status\tvalidity\tsubjects\tloci\tcorrect\taccuracy\tartifact_path\tartifact_sha256\tabstract_allowed\tclaim_boundary\n"
        f"GOOD\tdev\twes\tmethod\tdiscovery\tvalid\t1\t3\t3\t1.0\tartifact.tsv\t{sha256(artifact)}\t0\ttest\n"
        f"BAD\tdev\twgs\tmethod\tdiscovery\tinvalid\t1\t3\t1\t0.333\tartifact.tsv\t{sha256(artifact)}\t0\ttest\n",
        encoding="utf-8",
    )
    return registry


def test_registry_rejects_invalid_abstract_claim(tmp_path: Path):
    registry = _registry(tmp_path)
    text = registry.read_text(encoding="utf-8").replace("BAD\tdev\twgs\tmethod\tdiscovery\tinvalid", "BAD\tdev\twgs\tmethod\tdiscovery\tinvalid").replace("\t0\ttest\n", "\t1\ttest\n", 2)
    registry.write_text(text, encoding="utf-8")
    assert validate_registry(str(registry), str(tmp_path))


def test_manuscript_accepts_registered_valid_result(tmp_path: Path):
    registry = _registry(tmp_path)
    manuscript_dir = tmp_path / "benchmark"
    manuscript_dir.mkdir()
    manuscript = manuscript_dir / "manuscript.md"
    manuscript.write_text("# Results\nVerified result [RESULT:GOOD].\n", encoding="utf-8")
    claims = tmp_path / "claims.tsv"
    claims.write_text(
        "claim_id\tdraft\tsection\tclaim\tstatus\tevidence_result_id\taction\n"
        "C1\tbenchmark\tResults\tVerified\tsupported\tGOOD\tretain\n",
        encoding="utf-8",
    )
    result = audit_claims(str(manuscript), str(registry), str(claims), str(tmp_path / "audit.json"), str(tmp_path))
    assert result["submission_ready"] is True


def test_manuscript_rejects_placeholder_and_invalid_wgs_number(tmp_path: Path):
    registry = _registry(tmp_path)
    manuscript = tmp_path / "source.md"
    manuscript.write_text("[AUTHOR: confirm] WGS was 0.399 and 0.494. [RESULT:BAD]\n", encoding="utf-8")
    claims = tmp_path / "claims.tsv"
    claims.write_text(
        "claim_id\tdraft\tsection\tclaim\tstatus\tevidence_result_id\taction\n"
        "C1\tsource_original\tResults\tbad\tremove\tBAD\tremove\n",
        encoding="utf-8",
    )
    result = audit_claims(str(manuscript), str(registry), str(claims), str(tmp_path / "audit.json"), str(tmp_path))
    assert result["submission_ready"] is False
    assert any("invalid-WGS" in failure for failure in result["failures"])

