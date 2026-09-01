from pathlib import Path

from champhla_recovery.exposure import build_exposure_ledger
from champhla_recovery.io import read_json, read_tsv


def test_exposure_ledger_deduplicates_and_records_roles(tmp_path: Path):
    folds = tmp_path / "training_folds.tsv"
    folds.write_text("subject\tfold\nHG00096\t1\nHG00096\t2\nNA06985\t3\n", encoding="utf-8")
    pilot = tmp_path / "pilot.tsv"
    pilot.write_text("sample_id\nHG00097\n", encoding="utf-8")
    output = tmp_path / "ledger.tsv"
    summary = tmp_path / "summary.json"
    result = build_exposure_ledger([str(folds), str(pilot)], str(output), str(summary))
    assert result["unique_subjects"] == 3
    rows = read_tsv(output)
    assert {row["subject"] for row in rows} == {"HG00096", "HG00097", "NA06985"}
    assert next(row for row in rows if row["subject"] == "HG00096")["exposure_role"] == "development_exposed"
    assert read_json(summary)["output_sha256"]

