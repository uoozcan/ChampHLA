from __future__ import annotations

import hashlib
import importlib.util
import json
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "workflow" / "bin" / "patch_polysolver_wrapper.py"
SPEC = importlib.util.spec_from_file_location("polysolver_wrapper_patch", SCRIPT)
assert SPEC and SPEC.loader
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


INTERVALS = [
    "6:29941260-29945884",
    "6:31353872-31357187",
    "6:31268749-31272105",
]


def make_fixture(tmp_path: Path):
    source = tmp_path / "source.sh"
    source.write_text(
        "TMP_DIR=/home/polysolver\n"
        "if [ $build == \"hg38\" ]; then\n"
        + "\n".join(f"samtools view $bam {token}" for token in INTERVALS)
        + "\nfi\n",
        encoding="utf-8",
    )
    spec = tmp_path / "spec.json"
    spec.write_text(json.dumps({
        "schema_version": "champhla-polysolver-wrapper-patch-1",
        "transformation_version": "test-1",
        "source_sha256": hashlib.sha256(source.read_bytes()).hexdigest(),
        "temporary_directory_token": "TMP_DIR=/home/polysolver",
        "hg38_interval_tokens": INTERVALS,
        "expected_temporary_directory_substitutions": 1,
        "expected_hg38_contig_substitutions_for_chr6": 3,
    }), encoding="utf-8")
    return source, spec


@pytest.mark.parametrize("contig,expected", [("6", 0), ("chr6", 3)])
def test_polysolver_transform_is_exact_and_audited(tmp_path: Path, contig: str, expected: int):
    source, spec = make_fixture(tmp_path)
    output, audit = tmp_path / "derived.sh", tmp_path / "audit.json"
    record = MODULE.transform(source, output, audit, spec, contig, "/work/picard_tmp")
    text = output.read_text(encoding="utf-8")
    assert record["hg38_contig_substitutions"] == expected
    assert text.count("TMP_DIR=/work/picard_tmp") == 1
    assert sum(text.count("chr" + token) for token in INTERVALS) == expected
    assert MODULE.verify(source, output, audit, spec, contig, "/work/picard_tmp") == record


def test_polysolver_transform_rejects_contig_source_and_count_drift(tmp_path: Path):
    source, spec = make_fixture(tmp_path)
    output, audit = tmp_path / "derived.sh", tmp_path / "audit.json"
    with pytest.raises(ValueError, match="contig"):
        MODULE.transform(source, output, audit, spec, "chromosome6", "/tmp/p")
    source.write_text(source.read_text(encoding="utf-8") + "# altered\n", encoding="utf-8")
    with pytest.raises(ValueError, match="source SHA-256"):
        MODULE.transform(source, output, audit, spec, "chr6", "/tmp/p")
    source, spec = make_fixture(tmp_path)
    content = source.read_text(encoding="utf-8").replace(INTERVALS[-1], "6:1-2")
    source.write_text(content, encoding="utf-8")
    data = json.loads(spec.read_text(encoding="utf-8"))
    data["source_sha256"] = hashlib.sha256(source.read_bytes()).hexdigest()
    spec.write_text(json.dumps(data), encoding="utf-8")
    with pytest.raises(ValueError, match="occurred 0 times"):
        MODULE.transform(source, output, audit, spec, "chr6", "/tmp/p")


def test_polysolver_derived_wrapper_hash_drift_is_rejected(tmp_path: Path):
    source, spec = make_fixture(tmp_path)
    output, audit = tmp_path / "derived.sh", tmp_path / "audit.json"
    MODULE.transform(source, output, audit, spec, "chr6", "/tmp/p")
    output.write_text(output.read_text(encoding="utf-8") + "# drift\n", encoding="utf-8")
    with pytest.raises(ValueError, match="derived_wrapper_sha256"):
        MODULE.verify(source, output, audit, spec, "chr6", "/tmp/p")


def test_frozen_polysolver_patch_source_matches_attestation():
    spec = json.loads((ROOT / "workflow/conf/polysolver_wrapper_patch.json").read_text())
    attestation = json.loads((ROOT / "configs/roihu_caller_reference_attestation.json").read_text())
    wrapper = next(
        item for item in attestation["callers"]["POLYSOLVER"]["reference_components"]
        if item["artifact_id"] == "scripts/shell_call_hla_type"
    )
    assert spec["source_sha256"] == wrapper["sha256"]
