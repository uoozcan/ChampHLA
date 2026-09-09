"""Release identification must be decided by content, and must refuse to guess."""
from pathlib import Path

import pytest

from champhla_confirmation.imgt_release import (
    NOT_RELEASE_BEARING,
    NO_MATCH,
    caller_release_summary,
    declared_release,
    extract_two_field_types,
    identify_release,
    load_release_index,
    two_field,
)


def write_cache(root: Path, releases: dict[str, list[str]]) -> Path:
    cache = root / "cache"
    cache.mkdir()
    for release, alleles in releases.items():
        lines = ["# file: Allelelist.txt", f"# version: IPD-IMGT/HLA {release}"]
        lines += [f"HLA{index:05d},{allele}" for index, allele in enumerate(alleles, start=1)]
        (cache / f"Allelelist.{release}.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")
    return cache


def test_two_field_truncation_collapses_field_depth_differences():
    # A*32:55 and A*32:55:01 are the same two-field type; comparing full designations
    # across sources with different depths is what makes naive matching fail.
    assert two_field("A*32:55") == "A*32:55"
    assert two_field("HLA-A*32:55:01") == "A*32:55"
    assert two_field("B*15:01:01:02N") == "B*15:01"
    assert two_field("not-an-allele") is None


def test_component_resolves_to_the_release_that_first_contains_it(tmp_path: Path):
    cache = write_cache(tmp_path, {
        "3.10.0": ["A*01:01", "A*01:02"],
        "3.11.0": ["A*01:01", "A*01:02", "A*01:03"],
        "3.12.0": ["A*01:01", "A*01:02", "A*01:03", "A*01:04"],
    })
    source = tmp_path / "db.txt"
    source.write_text("A*01:01\nA*01:02\nA*01:03\n", encoding="utf-8")
    types = extract_two_field_types(source, "allele_list")
    record = identify_release(types, load_release_index(cache))
    assert record["release"] == "3.11.0"
    assert record["coverage_at_release"] == 1.0
    assert record["filtered_subset"] is False
    assert record["missing_from_predecessor"] == 1


def test_a_filtered_subset_is_flagged_rather_than_silently_treated_as_complete(tmp_path: Path):
    cache = write_cache(tmp_path, {
        "3.10.0": [f"A*01:{i:02d}" for i in range(1, 21)],
        "3.11.0": [f"A*01:{i:02d}" for i in range(1, 41)],
    })
    source = tmp_path / "db.txt"
    source.write_text("\n".join(f"A*01:{i:02d}" for i in range(1, 6)) + "\n", encoding="utf-8")
    record = identify_release(extract_two_field_types(source, "allele_list"), load_release_index(cache))
    assert record["release"] == "3.10.0"
    assert record["filtered_subset"] is True
    assert "filtered subset" in record["interpretation"]


def test_a_set_no_release_contains_is_reported_as_no_match(tmp_path: Path):
    """An unmatched database must never be resolved to its closest neighbour."""
    cache = write_cache(tmp_path, {"3.10.0": ["A*01:01"], "3.11.0": ["A*01:01", "A*01:02"]})
    source = tmp_path / "db.txt"
    source.write_text("A*01:01\nA*99:99\n", encoding="utf-8")
    record = identify_release(extract_two_field_types(source, "allele_list"), load_release_index(cache))
    assert record["release"] == NO_MATCH
    assert record["closest_release"] in {"3.10.0", "3.11.0"}
    assert record["missing_from_closest"] == 1


def test_polysolver_header_convention_is_parsed(tmp_path: Path):
    source = tmp_path / "poly.headers"
    source.write_text(">hla_a_01_01_01_01\n>hla_a_01_01_01_02n\n>hla_b_08_59\n", encoding="utf-8")
    assert extract_two_field_types(source, "polysolver_headers") == {"A*01:01", "B*08:59"}


def test_declared_release_reads_only_the_header_block(tmp_path: Path):
    source = tmp_path / "declared.txt"
    source.write_text(
        "# version: IPD-IMGT/HLA 3.38.0\nA*01:01\n# version: IPD-IMGT/HLA 9.99.9\n",
        encoding="utf-8",
    )
    assert declared_release(source) == "3.38.0"


def test_a_caller_mixing_releases_is_heterogeneous_not_unresolved():
    assert caller_release_summary(["3.38.0", "3.51.0"]) == "HETEROGENEOUS"
    assert caller_release_summary(["3.38.0", "3.38.0"]) == "3.38.0"
    assert caller_release_summary([NOT_RELEASE_BEARING, "3.15.0"]) == "3.15.0"
    assert caller_release_summary([NOT_RELEASE_BEARING]) == NO_MATCH
    assert caller_release_summary(["3.15.0", NO_MATCH]) == NO_MATCH


def test_identification_requires_a_pinned_cache_and_a_nonempty_set(tmp_path: Path):
    with pytest.raises(ValueError):
        load_release_index(tmp_path / "missing")
    empty = write_cache(tmp_path, {"3.10.0": ["A*01:01"]})
    with pytest.raises(ValueError):
        identify_release(set(), load_release_index(empty))
