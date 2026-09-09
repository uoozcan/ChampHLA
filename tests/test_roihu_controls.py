from __future__ import annotations

import json
from pathlib import Path

import pytest

from champhla_confirmation.hprc_truth import build_hprc_assembly_truth
from champhla_confirmation.io import read_tsv, write_tsv
from champhla_confirmation.roihu import (
    RUN_LEDGER_FIELDS,
    audit_run_manifest,
    assess_storage,
    build_cleanup_plan,
    build_nci60_run_manifest,
    initialize_run_ledger,
    freeze_run_manifest,
    transition_run_record,
    transition_run_sample,
)
from champhla_recovery.recount import independent_recount


MANIFEST_FIELDS = [
    "cohort", "sample_id", "donor_id", "modality", "input_type", "input_uri",
    "index_uri", "index_checksum", "source_checksum", "reference_build", "read_layout",
    "independence_stratum", "evidence_role",
]


def manifest_row(modality="wes", sample="S1"):
    return {
        "cohort": "C", "sample_id": sample, "donor_id": sample, "modality": modality,
        "input_type": "cram", "input_uri": f"https://example.org/{sample}.cram",
        "index_uri": f"https://example.org/{sample}.cram.crai",
        "index_checksum": "b" * 64, "source_checksum": "a" * 64,
        "reference_build": "GRCh38DH" if modality != "rnaseq" else "GRCh37",
        "read_layout": "paired", "independence_stratum": "new_library_overlap",
        "evidence_role": "same_resource_confirmation",
    }


def test_truth_free_manifest_and_ledger(tmp_path: Path):
    manifest = tmp_path / "manifest.tsv"
    write_tsv(manifest, [manifest_row()], MANIFEST_FIELDS)
    assert audit_run_manifest(manifest)["passed"] is True
    ledger = tmp_path / "ledger.tsv"
    rows = initialize_run_ledger(manifest, ledger, "deadbeef")
    assert len(rows) == 5
    running = transition_run_record(
        transition_run_record(rows[0], "submitted", job_id="1"), "running"
    )
    with pytest.raises(ValueError, match="requires exit_code=0"):
        transition_run_record(running, "succeeded")
    succeeded = transition_run_record(
        running, "succeeded", exit_code="0", output_sha256="c" * 64,
    )
    assert transition_run_record(succeeded, "validated")["state"] == "validated"
    frozen = freeze_run_manifest(manifest, tmp_path / "freeze.json", {"wes": 1})
    assert frozen["expected_caller_locus_records"] == 15
    assert frozen["expected_plurality_rows"] == 3
    with pytest.raises(ValueError, match="expected 2 unique samples"):
        freeze_run_manifest(manifest, tmp_path / "bad-freeze.json", {"wes": 2})
    submitted = initialize_run_ledger(manifest, ledger, "deadbeef", "123")
    assert all(row["state"] == "submitted" for row in submitted)
    transition_run_sample(ledger, "C", "S1", "wes", "running")
    transition_run_sample(
        ledger, "C", "S1", "wes", "succeeded", exit_code="0", output_sha256="d" * 64,
    )
    transition_run_sample(ledger, "C", "S1", "wes", "validated")
    assert all(row["state"] == "validated" for row in read_tsv(ledger))


def test_manifest_rejects_truth_duplicates_reference_and_missing_pair(tmp_path: Path):
    rows = [manifest_row(), manifest_row()]
    rows[0]["truth_allele1"] = "A*01:01"
    rows[0]["reference_build"] = "hg19"
    rows[0]["input_type"] = "fastq_pair"
    rows[0]["input_uri"] = "https://example.org/r1.fq.gz"
    rows[0]["index_uri"] = ""
    rows[0]["index_checksum"] = ""
    manifest = tmp_path / "bad.tsv"
    write_tsv(manifest, rows, MANIFEST_FIELDS + ["truth_allele1"])
    failures = audit_run_manifest(manifest)["failures"]
    assert any("truth-bearing" in failure for failure in failures)
    assert any("duplicate" in failure for failure in failures)
    assert any("requires 2 URI/checksum" in failure for failure in failures)
    assert any("expected GRCh38DH" in failure for failure in failures)


def test_storage_and_cleanup_are_fail_closed(tmp_path: Path):
    ledger = tmp_path / "ledger.tsv"
    rows = []
    for modality in ("wgs", "wes", "rnaseq"):
        row = {field: "" for field in RUN_LEDGER_FIELDS}
        row.update({"modality": modality, "state": "validated", "peak_disk_bytes": "1000"})
        rows.append(row)
    write_tsv(ledger, rows, list(RUN_LEDGER_FIELDS))
    result = assess_storage(ledger, {"wgs": 1, "wes": 1, "rnaseq": 1}, 200 * 1024 ** 3)
    assert result["passed"] is True
    run_root = tmp_path / "run"
    (run_root / "work").mkdir(parents=True)
    (run_root / "work" / "temporary").write_text("x")
    invalid = tmp_path / "invalid.json"
    invalid.write_text('{"valid": false}\n')
    plan = build_cleanup_plan(run_root, invalid, tmp_path / "cleanup.json")
    assert plan["executable"] is False
    valid = tmp_path / "valid.json"
    valid.write_text('{"valid": true}\n')
    plan = build_cleanup_plan(run_root, valid, tmp_path / "cleanup2.json")
    assert plan["executable"] is True
    assert (run_root / "work" / "temporary").exists()


def test_hprc_truth_requires_dual_method_concordance(tmp_path: Path):
    rows = []
    for gene in ("A", "B", "C"):
        for method in ("HLA-ASM", "Immuannot"):
            for haplotype, allele in (("1", f"{gene}*01:01"), ("2", f"{gene}*02:01")):
                rows.append({
                    "sample_id": "HG1", "haplotype": haplotype, "gene": gene,
                    "method": method, "allele": allele, "exon2_complete": "yes",
                    "exon3_complete": "yes", "equally_supported_conflict": "no",
                    "source_sha256": ("a" if method == "HLA-ASM" else "b") * 64,
                })
    calls = tmp_path / "assembly.tsv"
    write_tsv(calls, rows)
    truth = tmp_path / "truth.tsv"
    summary = build_hprc_assembly_truth(calls, truth, tmp_path / "audit.json")
    assert summary["resolved_loci"] == 3
    rows[-1]["allele"] = "C*03:01"
    write_tsv(calls, rows)
    summary = build_hprc_assembly_truth(calls, truth, tmp_path / "audit.json")
    assert summary["resolved_loci"] == 2
    assert next(row for row in read_tsv(truth) if row["gene"] == "C")["truth_status"] == "unresolved"


def test_independent_recount_does_not_pool_strata(tmp_path: Path):
    joined = tmp_path / "joined.tsv"
    base = {
        "cohort": "C", "modality": "wes", "gene": "A", "method": "SimplePluralityLex",
        "allele1": "A*01:01", "allele2": "A*02:01", "call_status": "callable",
        "truth_allele1": "A*01:01", "truth_allele2": "A*02:01", "truth_status": "resolved",
    }
    write_tsv(joined, [
        {**base, "subject": "S1", "donor": "D1", "independence_stratum": "donor_independent"},
        {**base, "subject": "S2", "donor": "D2", "independence_stratum": "new_library_overlap"},
    ])
    result = independent_recount(joined, tmp_path / "counts.tsv", tmp_path / "counts.json")
    assert len(result["rows"]) == 2
    assert all(row["correct"] == 1 for row in result["rows"])


def test_nci60_builder_requires_exact_ready_roster_and_checksums(tmp_path: Path):
    pilot = tmp_path / "pilot.tsv"
    pilot_rows = []
    report_rows = []
    for index in range(11):
        accession = f"SRR{index:08d}"
        for gene in ("A", "B", "C"):
            pilot_rows.append({
                "cohort": "NCI60_PUBLIC_PILOT", "subject": f"CELL{index}",
                "modality": "rnaseq", "gene": gene, "source_accession": accession,
                "gate_status": "ready",
            })
        report_rows.append({
            "run_accession": accession,
            "fastq_ftp": f"ftp.example/{accession}_1.fastq.gz;ftp.example/{accession}_2.fastq.gz",
            "fastq_md5": f"{'a' * 32};{'b' * 32}",
        })
    write_tsv(pilot, pilot_rows)
    report = tmp_path / "ena.tsv"
    write_tsv(report, report_rows)
    output = tmp_path / "nci60.tsv"
    result = build_nci60_run_manifest(pilot, report, output)
    assert result["passed"] is True
    assert result["samples_by_modality"] == {"rnaseq": 11}
    assert all(row["evidence_role"] == "exploratory" for row in read_tsv(output))
