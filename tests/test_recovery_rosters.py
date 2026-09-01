import json
from pathlib import Path

from champhla_recovery.io import read_tsv
from champhla_recovery.rosters import freeze_rosters


def _write(path: Path, text: str):
    path.write_text(text, encoding="utf-8")


def test_roster_is_truth_blind_deterministic_and_excludes_exposed(tmp_path: Path):
    registry = tmp_path / "registry.tsv"
    _write(registry, "subject\tpopulation\tprimary_eligible\tcomplete_truth_available\n" + "".join(
        f"HG{i:05d}\t{'EUR' if i % 2 else 'AFR'}\t1\t1\n" for i in range(1, 10)
    ))
    assay = tmp_path / "assay.tsv"
    _write(assay, "sample_id\tmodality\tpopulation\tinput_url\tindex_url\tinput_md5\tgenome_build\n" + "".join(
        f"HG{i:05d}\t{mod}\t{'EUR' if i % 2 else 'AFR'}\thttps://x/{i}\t\t\tGRCh38DH\n"
        for mod in ("wgs", "wes", "rnaseq") for i in range(1, 10)
    ))
    exposure = tmp_path / "exposure.tsv"
    _write(exposure, "subject\nHG00001\n")
    config = tmp_path / "config.json"
    config.write_text(json.dumps({
        "selection_seed": "fixed", "initial_targets": {"wgs": 3, "wes": 3, "rnaseq": 3}
    }), encoding="utf-8")
    first = tmp_path / "first.tsv"
    second = tmp_path / "second.tsv"
    freeze_rosters(str(registry), str(assay), str(exposure), str(config), str(first), str(tmp_path / "a.json"))
    freeze_rosters(str(registry), str(assay), str(exposure), str(config), str(second), str(tmp_path / "b.json"))
    assert first.read_bytes() == second.read_bytes()
    rows = read_tsv(first)
    assert len(rows) == 9
    assert "HG00001" not in {row["subject"] for row in rows}
    assert all("truth" not in field.lower() for field in rows[0])


def test_roster_records_infeasible_modality_without_weakening_truth(tmp_path: Path):
    registry = tmp_path / "registry.tsv"
    _write(registry, "subject\tpopulation\tprimary_eligible\nHG00001\tEUR\t1\n")
    assay = tmp_path / "assay.tsv"
    _write(assay, "sample_id\tmodality\tpopulation\n" + "".join(
        f"HG00001\t{modality}\tEUR\n" for modality in ("wgs", "wes", "rnaseq")
    ))
    exposure = tmp_path / "exposure.tsv"
    _write(exposure, "subject\n")
    config = tmp_path / "config.json"
    config.write_text(json.dumps({
        "selection_seed": "fixed", "initial_targets": {"wgs": 1, "wes": 1, "rnaseq": 2}
    }), encoding="utf-8")
    result = freeze_rosters(
        str(registry), str(assay), str(exposure), str(config),
        str(tmp_path / "roster.tsv"), str(tmp_path / "summary.json"),
    )
    assert result["all_minimums_met"] is False
    assert result["modalities"]["rnaseq"]["status"] == "infeasible_minimum_not_met"
    assert result["modalities"]["wgs"]["status"] == "roster_frozen"
