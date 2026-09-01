from pathlib import Path

from champhla_recovery.io import read_json, read_tsv
from champhla_recovery.truth import prepare_locked_truth


def test_2014_truth_resolves_same_two_field_ambiguity_and_rejects_conflict(tmp_path: Path):
    source = tmp_path / "truth.txt"
    source.write_text(
        '"id" "sbgroup" "A" "A.1" "B" "B.1" "C" "C.1"\n'
        '"HG00001" "GBR" "01:01:01/01:01:02" "02:01" "07:02" "08:01" "07:01" "06:02"\n'
        '"HG00002" "FIN" "01:01/02:01" "03:01" "07:02" "08:01" "07:01" "06:02"\n',
        encoding="utf-8",
    )
    truth = tmp_path / "sealed" / "truth.tsv"
    registry = tmp_path / "registry.tsv"
    manifest = tmp_path / "manifest.json"
    result = prepare_locked_truth(str(source), str(truth), str(registry), str(manifest))
    assert result["eligible_complete_subjects"] == 1
    assert result["eligible_any_exact_locus_subjects"] == 2
    assert len(read_tsv(truth)) == 15
    rows = {row["subject"]: row for row in read_tsv(registry)}
    assert rows["HG00001"]["complete_truth_available"] == "1"
    assert rows["HG00002"]["primary_eligible"] == "1"
    assert rows["HG00002"]["exact_truth_loci"] == "2"
    assert read_json(manifest)["inferred_2018_truth_allowed"] is False


def test_conflicting_duplicate_subject_is_ineligible(tmp_path: Path):
    source = tmp_path / "truth.txt"
    header = '"id" "sbgroup" "A" "A.1" "B" "B.1" "C" "C.1"\n'
    source.write_text(
        header
        + '"HG00001" "GBR" "01:01" "02:01" "07:02" "08:01" "07:01" "06:02"\n'
        + '"HG00001" "GBR" "03:01" "02:01" "07:02" "08:01" "07:01" "06:02"\n',
        encoding="utf-8",
    )
    result = prepare_locked_truth(
        str(source), str(tmp_path / "truth.tsv"), str(tmp_path / "registry.tsv"),
        str(tmp_path / "manifest.json"),
    )
    assert result["eligible_complete_subjects"] == 0
    assert result["eligible_any_exact_locus_subjects"] == 1
    assert result["duplicate_conflicts"] == 1
