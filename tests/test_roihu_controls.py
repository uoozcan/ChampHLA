from __future__ import annotations

import json
import os
import subprocess
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
    build_same_resource_run_manifest,
    collect_run_outputs,
    directory_tree_identity,
    freeze_workflow_lock,
    initialize_run_ledger,
    reconcile_terminal_run_sample,
    freeze_run_manifest,
    transition_run_record,
    transition_run_sample,
    validate_environment_inventory,
    validate_workflow_lock,
)
from champhla_confirmation.staging import (
    build_run_disposition,
    classify_stage_failure,
    initialize_stage_ledger,
    transition_stage_attempt,
    validate_certificate_pair,
    finalize_sample,
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


@pytest.mark.parametrize("message", [
    "Container header CRC32 failure", "connection reset by peer",
    "operation timed out", "temporary failure in name resolution",
    "HTTP response code 429", "HTTP 503 Service Unavailable",
    "Seek at offset 12003535533 failed",
    'samtools view: error closing "https://example.org/a.cram": -1',
    'samtools view: error reading file "https://example.org/a.cram"',
    "EOF marker is absent. The input is probably truncated",
])
def test_remote_stage_failure_transient_allowlist(message):
    result = classify_stage_failure(1, message, "https://example.org/a.cram")
    assert result["failure_class"] == "transient_transport"
    assert result["retryable"] is True


@pytest.mark.parametrize("message", [
    "HTTP 403 Forbidden", "HTTP 404 Not Found", "reference mismatch",
    "unknown contig chr6", "checksum mismatch", "no such file",
    "some generic failure",
])
def test_stage_failure_deterministic(message):
    result = classify_stage_failure(1, message, "https://example.org/a.cram")
    assert result["failure_class"] == "deterministic"
    assert result["retryable"] is False


def test_remote_stage_timeout_is_retryable_but_local_timeout_is_not():
    assert classify_stage_failure(124, "", "https://example.org/a.cram")["retryable"]
    assert not classify_stage_failure(124, "", "/data/a.cram")["retryable"]


def test_generic_local_read_error_remains_deterministic():
    result = classify_stage_failure(
        1, 'samtools view: error reading file "/data/a.cram"',
        "https://example.org/a.cram",
    )
    assert result["failure_class"] == "deterministic"
    assert result["retryable"] is False


def test_stage_attempt_ledger_is_append_only_and_validated(tmp_path: Path):
    manifest = tmp_path / "manifest.tsv"
    ledger = tmp_path / "stage.tsv"
    write_tsv(manifest, [manifest_row()], MANIFEST_FIELDS)
    initialize_stage_ledger(manifest, ledger, "pilot-wes-stagev2", "technical_pilot")
    transition_stage_attempt(ledger, "C", "S1", "wes", 1, "submitted",
                             scheduler_job_id="101")
    transition_stage_attempt(ledger, "C", "S1", "wes", 1, "running")
    transition_stage_attempt(ledger, "C", "S1", "wes", 1, "failed",
                             exit_code="1", failure_class="transient_transport")
    transition_stage_attempt(ledger, "C", "S1", "wes", 2, "submitted", True,
                             scheduler_job_id="101", superseded_attempt="1")
    transition_stage_attempt(ledger, "C", "S1", "wes", 2, "running")
    with pytest.raises(ValueError, match="requires exit_code=0"):
        transition_stage_attempt(ledger, "C", "S1", "wes", 2, "validated")
    transition_stage_attempt(ledger, "C", "S1", "wes", 2, "validated",
                             exit_code="0", extracted_sha256="f" * 64)
    rows = read_tsv(ledger)
    assert [row["state"] for row in rows] == ["failed", "validated"]
    with pytest.raises(ValueError, match="unique"):
        transition_stage_attempt(ledger, "C", "S1", "wes", 2, "submitted", True)


def test_certificate_pair_detects_match_mismatch_and_expiry(tmp_path: Path, monkeypatch):
    identity, certificate = tmp_path / "id", tmp_path / "id-cert.pub"
    identity.write_text("private-not-reported", encoding="utf-8")
    certificate.write_text("certificate", encoding="utf-8")
    fingerprint = "SHA256:AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
    cert_fingerprint = fingerprint

    def fake_run(command, **kwargs):
        nonlocal cert_fingerprint
        if command[1] == "-y":
            return subprocess.CompletedProcess(command, 0, "ssh-ed25519 AAAA\n", "")
        if command[1] == "-lf":
            return subprocess.CompletedProcess(command, 0, f"256 {fingerprint} key (ED25519)\n", "")
        return subprocess.CompletedProcess(
            command, 0,
            f"Public key: ED25519-CERT {cert_fingerprint}\n"
            "Valid: from 2026-09-14T00:00:00 to 2026-09-16T00:00:00\n", "")

    monkeypatch.setattr(subprocess, "run", fake_run)
    active = 1789344000  # 2026-09-15T00:00:00 UTC
    assert validate_certificate_pair(identity, certificate, active)["passed"] is True
    cert_fingerprint = "SHA256:BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB"
    mismatch = validate_certificate_pair(identity, certificate, active)
    assert mismatch["passed"] is False and mismatch["fingerprint_match"] is False
    cert_fingerprint = fingerprint
    expired = validate_certificate_pair(identity, certificate, 1789516801)
    assert expired["passed"] is False and expired["valid_now"] is False


def test_run_disposition_contains_operations_but_no_accuracy(tmp_path: Path):
    manifest = tmp_path / "manifest.tsv"
    stages, callers, output = tmp_path / "stages.tsv", tmp_path / "callers.tsv", tmp_path / "out.tsv"
    write_tsv(manifest, [manifest_row()], MANIFEST_FIELDS)
    initialize_stage_ledger(manifest, stages, "pilot-wes-stagev2", "technical_pilot")
    initialize_run_ledger(manifest, callers, "pilot-wes-stagev2", "technical_pilot", "deadbeef")
    rows = build_run_disposition([stages], [callers], output)
    assert len(rows) == 6
    assert {row["evidence_layer"] for row in rows} == {"staging", "caller"}
    assert "accuracy" not in read_tsv(output)[0]


def test_finalizer_passes_only_with_validated_stage_callers_and_marker(tmp_path: Path):
    manifest = tmp_path / "manifest.tsv"
    stages, callers = tmp_path / "stages.tsv", tmp_path / "callers.tsv"
    marker = tmp_path / "CALLERS_COMPLETE"
    write_tsv(manifest, [manifest_row()], MANIFEST_FIELDS)
    initialize_stage_ledger(manifest, stages, "pilot-wes-stagev2", "technical_pilot")
    initialize_run_ledger(manifest, callers, "pilot-wes-stagev2", "technical_pilot", "deadbeef")
    transition_stage_attempt(stages, "C", "S1", "wes", 1, "submitted", scheduler_job_id="1_1")
    transition_stage_attempt(stages, "C", "S1", "wes", 1, "running")
    transition_stage_attempt(stages, "C", "S1", "wes", 1, "validated",
                             exit_code="0", extracted_sha256="e" * 64)
    transition_run_sample(callers, "C", "S1", "wes", "submitted", job_id="2_1")
    transition_run_sample(callers, "C", "S1", "wes", "running")
    transition_run_sample(callers, "C", "S1", "wes", "succeeded",
                          exit_code="0", output_sha256="f" * 64)
    transition_run_sample(callers, "C", "S1", "wes", "validated")
    assert finalize_sample(stages, callers, "C", "S1", "wes", "COMPLETED",
                           "COMPLETED", marker)["passed"] is False
    marker.touch()
    assert finalize_sample(stages, callers, "C", "S1", "wes", "COMPLETED",
                           "COMPLETED", marker)["passed"] is True


def test_truth_free_manifest_and_ledger(tmp_path: Path):
    manifest = tmp_path / "manifest.tsv"
    write_tsv(manifest, [manifest_row()], MANIFEST_FIELDS)
    assert audit_run_manifest(manifest)["passed"] is True
    ledger = tmp_path / "ledger.tsv"
    rows = initialize_run_ledger(
        manifest, ledger, "pilot-wes-001", "technical_pilot", "deadbeef",
    )
    assert len(rows) == 5
    assert {row["run_id"] for row in rows} == {"pilot-wes-001"}
    assert {row["run_role"] for row in rows} == {"technical_pilot"}
    running = transition_run_record(
        transition_run_record(rows[0], "submitted", job_id="1"), "running"
    )
    with pytest.raises(ValueError, match="requires exit_code=0"):
        transition_run_record(running, "succeeded")
    succeeded = transition_run_record(
        running, "succeeded", exit_code="0", output_sha256="c" * 64,
    )
    assert transition_run_record(succeeded, "validated")["state"] == "validated"

    cancelled = transition_run_record(rows[1], "submitted", job_id="1302428")
    cancelled = transition_run_record(
        cancelled, "cancelled", exit_code="dependency_cancelled",
    )
    assert cancelled["state"] == "cancelled"
    with pytest.raises(ValueError, match="invalid run-ledger transition"):
        transition_run_record(cancelled, "running")

    failed = transition_run_record(running, "failed", exit_code="1")
    with pytest.raises(ValueError, match="supersedes_job_id"):
        transition_run_record(failed, "resubmitted", job_id="2", attempt="2")
    retried = transition_run_record(
        failed, "resubmitted", job_id="2", attempt="2", supersedes_job_id="1",
    )
    assert retried["attempt"] == "2"
    assert retried["supersedes_job_id"] == "1"
    frozen = freeze_run_manifest(manifest, tmp_path / "freeze.json", {"wes": 1})
    assert frozen["expected_caller_locus_records"] == 15
    assert frozen["expected_plurality_rows"] == 3
    with pytest.raises(ValueError, match="expected 2 unique samples"):
        freeze_run_manifest(manifest, tmp_path / "bad-freeze.json", {"wes": 2})
    submitted = initialize_run_ledger(
        manifest, ledger, "pilot-wes-001", "technical_pilot", "deadbeef", "123",
    )
    assert all(row["state"] == "submitted" for row in submitted)
    transition_run_sample(ledger, "C", "S1", "wes", "running")
    transition_run_sample(
        ledger, "C", "S1", "wes", "succeeded", exit_code="0", output_sha256="d" * 64,
    )
    transition_run_sample(ledger, "C", "S1", "wes", "validated")
    assert all(row["state"] == "validated" for row in read_tsv(ledger))


def test_terminal_reconciliation_closes_mixed_active_caller_rows(tmp_path: Path):
    manifest, ledger = tmp_path / "manifest.tsv", tmp_path / "ledger.tsv"
    write_tsv(manifest, [manifest_row()], MANIFEST_FIELDS)
    rows = initialize_run_ledger(manifest, ledger, "pilot-wes-stagev2",
                                 "technical_pilot", "deadbeef")
    rows[0] = transition_run_record(rows[0], "submitted", job_id="100_1")
    rows[0] = transition_run_record(rows[0], "running")
    rows[1] = transition_run_record(rows[1], "submitted", job_id="100_1")
    write_tsv(ledger, rows, list(RUN_LEDGER_FIELDS))
    reconcile_terminal_run_sample(ledger, "C", "S1", "wes", "CANCELLED by dependency")
    states = {row["state"] for row in read_tsv(ledger)}
    assert states == {"failed", "cancelled"}
    assert not states.intersection({"pending", "submitted", "running", "succeeded"})


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
    ledgers = []
    for modality in ("wgs", "wes", "rnaseq"):
        ledger = tmp_path / f"{modality}.tsv"
        rows = []
        for sample in ("S1", "S2"):
            for caller in __import__("champhla_confirmation.panels", fromlist=["PANELS"]).PANELS[modality]:
                row = {field: "" for field in RUN_LEDGER_FIELDS}
                row.update({
                    "run_id": f"pilot-{modality}-001", "run_role": "technical_pilot",
                    "cohort": "C", "sample_id": sample, "caller": caller,
                    "modality": modality, "state": "validated", "runtime_seconds": "10",
                    "peak_memory_bytes": "200", "peak_disk_bytes": "1000",
                    "retained_disk_bytes": "500", "git_commit": "deadbeef",
                })
                rows.append(row)
        write_tsv(ledger, rows, list(RUN_LEDGER_FIELDS))
        ledgers.append(ledger)
    result = assess_storage(
        ledgers, {"wgs": 1, "wes": 1, "rnaseq": 1}, 200 * 1024 ** 3,
    )
    assert result["passed"] is True
    assert result["execution_mode"] == "full_scale"
    assert result["availability_source"] == "project_allocation"
    assert result["validated_pilot_samples_by_modality"] == {
        "rnaseq": 2, "wes": 2, "wgs": 2,
    }
    mixed_rows = read_tsv(ledgers[0])
    mixed_rows[0]["git_commit"] = "different"
    write_tsv(ledgers[0], mixed_rows, list(RUN_LEDGER_FIELDS))
    mixed = assess_storage(ledgers, {"wgs": 1, "wes": 1, "rnaseq": 1}, 200 * 1024 ** 3)
    assert mixed["passed"] is False
    assert any("one nonempty Git commit" in failure for failure in mixed["failures"])
    mixed_rows[0]["git_commit"] = "deadbeef"
    write_tsv(ledgers[0], mixed_rows, list(RUN_LEDGER_FIELDS))
    rejected = assess_storage(
        ledgers, {"wgs": 1, "wes": 1, "rnaseq": 1}, 200 * 1024 ** 3,
        availability_source="filesystem_global",
    )
    assert rejected["passed"] is False
    low_space = assess_storage(
        ledgers, {"wgs": 1_000_000, "wes": 1, "rnaseq": 1},
        100 * 1024 ** 3 + 600_000_000,
    )
    assert low_space["full_scale_passed"] is False
    assert low_space["sequential_low_storage_passed"] is True
    assert low_space["execution_mode"] == "sequential_low_storage"
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


def test_run_identity_and_storage_measurements_are_fail_closed(tmp_path: Path):
    manifest = tmp_path / "manifest.tsv"
    write_tsv(manifest, [manifest_row()], MANIFEST_FIELDS)
    with pytest.raises(ValueError, match="run_id"):
        initialize_run_ledger(manifest, tmp_path / "ledger.tsv", "BAD ID", "technical_pilot")
    with pytest.raises(ValueError, match="run_role"):
        initialize_run_ledger(manifest, tmp_path / "ledger.tsv", "valid-run", "pilot")

    ledgers = []
    from champhla_confirmation.panels import PANELS
    for modality in ("wgs", "wes", "rnaseq"):
        rows = []
        for sample in ("S1", "S2"):
            for caller in PANELS[modality]:
                row = {field: "" for field in RUN_LEDGER_FIELDS}
                row.update({
                    "run_id": f"pilot-{modality}", "run_role": "technical_pilot",
                    "cohort": "C", "sample_id": sample, "modality": modality,
                    "caller": caller, "state": "validated", "runtime_seconds": "1",
                    "peak_memory_bytes": "2", "peak_disk_bytes": "100",
                    "retained_disk_bytes": "10",
                })
                rows.append(row)
        path = tmp_path / f"{modality}.tsv"
        write_tsv(path, rows, list(RUN_LEDGER_FIELDS))
        ledgers.append(path)
    rows = read_tsv(ledgers[1])
    rows[0]["peak_disk_bytes"] = "101"
    write_tsv(ledgers[1], rows, list(RUN_LEDGER_FIELDS))
    result = assess_storage(ledgers, {"wgs": 1, "wes": 1, "rnaseq": 1}, 200 * 1024 ** 3)
    assert result["passed"] is False
    assert any("inconsistent peak_disk_bytes" in failure for failure in result["failures"])
    with pytest.raises(ValueError, match="exactly three"):
        assess_storage(ledgers[:2], {"wgs": 1, "wes": 1}, 200 * 1024 ** 3)


def test_sequential_storage_uses_p95_of_paired_sample_transients(tmp_path: Path):
    """Crossed marginal maxima must not erase the larger sample transient."""
    from champhla_confirmation.panels import PANELS

    ledgers = []
    for modality in ("wgs", "wes", "rnaseq"):
        rows = []
        measurements = ((1000, 900), (900, 100))
        for sample, (peak, retained) in zip(("S1", "S2"), measurements):
            for caller in PANELS[modality]:
                row = {field: "" for field in RUN_LEDGER_FIELDS}
                row.update({
                    "run_id": f"pilot-{modality}",
                    "run_role": "technical_pilot", "cohort": "C",
                    "sample_id": sample, "modality": modality, "caller": caller,
                    "state": "validated", "runtime_seconds": "1",
                    "peak_memory_bytes": "2", "peak_disk_bytes": str(peak),
                    "retained_disk_bytes": str(retained), "git_commit": "deadbeef",
                })
                rows.append(row)
        ledger = tmp_path / f"{modality}.tsv"
        write_tsv(ledger, rows, list(RUN_LEDGER_FIELDS))
        ledgers.append(ledger)

    result = assess_storage(
        ledgers, {"wgs": 1, "wes": 1, "rnaseq": 1}, 200 * 1024 ** 3,
    )
    assert result["pilot_transient_p95_bytes_by_modality"] == {
        "wgs": 800, "wes": 800, "rnaseq": 800,
    }
    expected = 2700 + 1000 + 100 * 1024 ** 3
    assert result["sequential_required_available_bytes"] == expected


def test_canonical_collection_rejects_nonproduction_run_role(tmp_path: Path):
    manifest = tmp_path / "manifest.tsv"
    write_tsv(manifest, [manifest_row()], MANIFEST_FIELDS)
    ledger = tmp_path / "ledger.tsv"
    initialize_run_ledger(
        manifest, ledger, "pilot-wes-001", "technical_pilot", "deadbeef", "1",
    )
    summary = collect_run_outputs(
        tmp_path / "callers", manifest, [ledger], tmp_path / "calls.tsv",
        tmp_path / "summary.json",
    )
    assert summary["passed"] is False
    assert "canonical collection requires run_role=production" in summary["failures"]


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


def test_same_resource_builder_has_exact_truth_free_production_matrix(tmp_path: Path):
    root = Path(__file__).resolve().parents[1]
    output = tmp_path / "same_resource.tsv"
    roster_keys = {
        (row["subject"], row["modality"]) for row in read_tsv(
            root / "cohorts" / "same_resource_truth_free_roster.tsv"
        )
    }
    index_rows = [
        {"sample_id": row["sample_id"], "modality": row["modality"],
         "index_uri": row["index_url"], "bytes": "1", "md5": "c" * 32,
         "sha256": "d" * 64}
        for row in read_tsv(root / "cohorts" / "official_1000g_assay_manifest.tsv")
        if (row["sample_id"], row["modality"]) in roster_keys
        and row["modality"] != "rnaseq" and not row["index_md5"]
    ]
    index_registry = tmp_path / "index_checksums.tsv"
    write_tsv(index_registry, index_rows)
    result = build_same_resource_run_manifest(
        root / "cohorts" / "same_resource_truth_free_roster.tsv",
        root / "cohorts" / "official_1000g_assay_manifest.tsv",
        root / "cohorts" / "sources" / "geuvadis_ena_fastq_report.tsv",
        index_registry,
        output,
    )
    assert result["passed"] is True
    assert result["samples_by_modality"] == {"wgs": 137, "wes": 130, "rnaseq": 107}
    frozen = freeze_run_manifest(
        output, tmp_path / "same_resource.freeze.json",
        {"wgs": 137, "wes": 130, "rnaseq": 107},
    )
    assert frozen["expected_caller_locus_records"] == 5289
    assert frozen["expected_plurality_rows"] == 1122
    assert all(row["evidence_role"] == "same_resource_confirmation" for row in read_tsv(output))


def test_directory_tree_hash_uses_posix_relative_paths(tmp_path: Path):
    root = tmp_path / "reference"
    (root / "nested").mkdir(parents=True)
    (root / "nested" / "a.txt").write_text("A\n", encoding="utf-8")
    first = directory_tree_identity(root)
    assert first["present"] is True
    assert first["files"] == 1
    assert len(first["tree_sha256"]) == 64
    (root / "nested" / "a.txt").write_text("B\n", encoding="utf-8")
    assert directory_tree_identity(root)["tree_sha256"] != first["tree_sha256"]


def test_environment_inventory_requires_version_identity_freshness_and_post_pilot_time(tmp_path: Path):
    path = tmp_path / "environment.json"
    payload = {
        "schema_version": "champhla-roihu-environment-inventory-2",
        "evidence_id": "post-pilot-20260915t170000z",
        "created_at_utc": "2026-09-15T17:00:00Z",
        "passed": True,
        "project_storage": {
            "free_bytes": 123, "filesystem_global_space_ignored": True,
        },
        "repository": {
            "clean": True, "commit": "deadbeef", "tracked_tree_sha256": "a" * 64,
        },
        "workflow_lock": {"passed": True},
    }
    path.write_text(json.dumps(payload), encoding="utf-8")
    assert validate_environment_inventory(
        path, max_age_seconds=3600,
        not_before_utc="2026-09-15T16:59:59Z",
        now_utc="2026-09-15T17:30:00Z",
    ) == []

    payload["created_at_utc"] = "2026-09-15T16:00:00Z"
    path.write_text(json.dumps(payload), encoding="utf-8")
    failures = validate_environment_inventory(
        path, max_age_seconds=3600,
        not_before_utc="2026-09-15T16:30:00Z",
        now_utc="2026-09-15T17:30:01Z",
    )
    assert "environment inventory predates the pilot evidence" in failures
    assert "environment inventory is stale" in failures

    payload.pop("evidence_id")
    path.write_text(json.dumps(payload), encoding="utf-8")
    assert any("evidence_id" in failure for failure in validate_environment_inventory(path))


@pytest.mark.skipif(os.name == "nt", reason="Windows CI does not guarantee symlink privilege")
def test_directory_tree_hash_freezes_internal_links_and_rejects_escapes(tmp_path: Path):
    root = tmp_path / "index"
    root.mkdir()
    (root / "target.fa").write_text(">A\nAC\n")
    (root / "alias.fa").symlink_to("target.fa")
    internal = directory_tree_identity(root)
    assert internal["failures"] == []
    assert internal["files"] == 2
    (root / "absolute_alias.fa").symlink_to((root / "target.fa").resolve())
    assert directory_tree_identity(root)["failures"] == []
    outside = tmp_path / "outside.fa"
    outside.write_text(">B\nGT\n")
    (root / "escape.fa").symlink_to(outside)
    assert any("escapes" in failure for failure in directory_tree_identity(root)["failures"])


def test_workflow_lock_detects_changes_and_masked_success(tmp_path: Path):
    project = tmp_path / "project"
    workflow = project / "workflow"
    (workflow / "modules").mkdir(parents=True)
    (workflow / "conf").mkdir()
    required = {
        "main.nf": "nextflow.enable.dsl = 2\n",
        "nextflow.config": "process.errorStrategy = 'terminate'\n",
        "conf/roihu_params.yaml": "hlahd_db: /app/hlahd.1.4.0\n",
        "conf/polysolver_wrapper_patch.json": "{}\n",
        "bin/detect_chr6_contig.sh": "#!/bin/bash\nprintf 'chr6\\n'\n",
        "bin/patch_polysolver_wrapper.py": "print('fixture')\n",
    }
    for module in ("arcashla", "bam_to_fastq", "hlahd", "kourami", "optitype",
                   "polysolver", "spechla", "t1k"):
        required[f"modules/{module}.nf"] = f"process {module.upper()} {{ script: 'exit 0' }}\n"
    for relative, content in required.items():
        target = workflow / relative
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text(content, encoding="utf-8")
    for relative in (
        "scripts/roihu_preflight.sh", "scripts/roihu_promote_preflight.sh",
        "scripts/roihu_assess_storage.sh", "scripts/roihu_promote_storage_gate.sh",
        "scripts/roihu_submit_wave.sh", "scripts/roihu_stage_inputs.sbatch",
        "scripts/roihu_run_sample.sbatch", "scripts/roihu_finalize_sample.sbatch",
        "src/champhla_confirmation/cli.py", "src/champhla_confirmation/io.py",
        "src/champhla_confirmation/manifests.py", "src/champhla_confirmation/panels.py",
        "src/champhla_confirmation/parsers.py", "src/champhla_confirmation/roihu.py",
        "src/champhla_confirmation/staging.py",
    ):
        target = project / relative
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text("# locked operational fixture\n", encoding="utf-8")
    lock = project / "lock.json"
    freeze_workflow_lock(project, workflow, lock)
    assert validate_workflow_lock(lock, project)["passed"] is True
    locked = json.loads(lock.read_text(encoding="utf-8"))["files"]
    assert "src/champhla_confirmation/roihu.py" in locked
    assert "src/champhla_confirmation/cli.py" in locked
    (workflow / "modules" / "t1k.nf").write_text("errorStrategy 'ignore'\n", encoding="utf-8")
    failures = validate_workflow_lock(lock, project)["failures"]
    assert any("hash mismatch" in failure for failure in failures)
    assert any("forbidden token" in failure for failure in failures)
