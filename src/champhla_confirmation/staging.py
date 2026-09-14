from __future__ import annotations

import hashlib
import os
import re
import subprocess
import time
from datetime import datetime, timezone
from pathlib import Path

from .io import read_tsv, write_tsv
from .roihu import HEX64, validate_run_identity


STAGE_LEDGER_FIELDS = (
    "run_id", "run_role", "cohort", "sample_id", "donor_id", "modality",
    "attempt", "scheduler_job_id", "state", "failure_class", "exit_code",
    "source_uri", "source_checksum_declared", "source_checksum_verified",
    "index_uri", "index_checksum_declared", "index_checksum_verified", "extracted_sha256",
    "runtime_seconds", "peak_storage_bytes", "superseded_attempt",
    "started_at_utc", "updated_at_utc",
)
STAGE_STATES = {
    "pending": {"submitted"},
    "submitted": {"running", "cancelled", "failed"},
    "running": {"validated", "failed"},
    "failed": {"superseded", "quarantined"},
    "validated": set(), "cancelled": set(), "superseded": set(), "quarantined": set(),
}
REMOTE_URI = re.compile(r"^(?:https?|ftp|s3)://", re.I)
DETERMINISTIC_PATTERNS = (
    r"\bHTTP[^\n]*(?:401|403|404)\b", r"\b(?:401|403|404)\b[^\n]*HTTP",
    r"reference.*(?:mismatch|not found|missing)", r"unknown reference",
    r"(?:contig|region).*(?:invalid|unknown|not found)", r"checksum.*(?:mismatch|failed)",
    r"(?:no such file|file not found|local corruption|truncated file)",
)
TRANSIENT_PATTERNS = (
    r"(?:CRAM|container header).*CRC(?:32)? failure", r"connection (?:reset|timed out)",
    r"operation timed out", r"temporary failure in name resolution",
    r"could not resolve host", r"\bHTTP[^\n]*429\b", r"\b429\b[^\n]*HTTP",
    r"\bHTTP[^\n]*5\d\d\b", r"\b5\d\d\b[^\n]*HTTP",
)


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat().replace("+00:00", "Z")


def classify_stage_failure(exit_code: int, stderr: str, source_uri: str) -> dict:
    """Classify conservatively: only an allowlisted remote transport error retries."""
    if exit_code == 0:
        return {"failure_class": "none", "retryable": False, "reason": "success"}
    if exit_code == 124 and REMOTE_URI.match(source_uri):
        return {"failure_class": "transient_transport", "retryable": True,
                "reason": "attempt timeout"}
    for pattern in DETERMINISTIC_PATTERNS:
        if re.search(pattern, stderr, re.I):
            return {"failure_class": "deterministic", "retryable": False,
                    "reason": pattern}
    if REMOTE_URI.match(source_uri):
        for pattern in TRANSIENT_PATTERNS:
            if re.search(pattern, stderr, re.I):
                return {"failure_class": "transient_transport", "retryable": True,
                        "reason": pattern}
    return {"failure_class": "deterministic", "retryable": False,
            "reason": f"unclassified exit {exit_code}"}


def initialize_stage_ledger(manifest_path: str | Path, output: str | Path,
                            run_id: str, run_role: str) -> list[dict[str, str]]:
    from .roihu import audit_run_manifest
    validate_run_identity(run_id, run_role)
    audit = audit_run_manifest(manifest_path)
    if not audit["passed"]:
        raise ValueError(f"run manifest failed: {audit['failures']}")
    rows = []
    for source in read_tsv(manifest_path):
        rows.append({
            "run_id": run_id, "run_role": run_role, "cohort": source["cohort"],
            "sample_id": source["sample_id"], "donor_id": source["donor_id"],
            "modality": source["modality"], "attempt": "1", "scheduler_job_id": "",
            "state": "pending", "failure_class": "", "exit_code": "",
            "source_uri": source["input_uri"],
            "source_checksum_declared": source["source_checksum"],
            "source_checksum_verified": "", "index_uri": source["index_uri"],
            "index_checksum_declared": source["index_checksum"],
            "index_checksum_verified": "", "extracted_sha256": "",
            "runtime_seconds": "", "peak_storage_bytes": "", "superseded_attempt": "",
            "started_at_utc": "", "updated_at_utc": "",
        })
    write_tsv(output, rows, list(STAGE_LEDGER_FIELDS))
    return rows


def _lock(path: Path) -> Path:
    lock = path.with_suffix(path.suffix + ".lockdir")
    deadline = time.monotonic() + 30
    while True:
        try:
            lock.mkdir()
            return lock
        except FileExistsError:
            if time.monotonic() >= deadline:
                raise TimeoutError(f"timed out acquiring stage-ledger lock: {lock}")
            time.sleep(.05)


def transition_stage_attempt(ledger_path: str | Path, cohort: str, sample_id: str,
                             modality: str, attempt: int, new_state: str,
                             append_attempt: bool = False, **updates: str) -> list[dict[str, str]]:
    path = Path(ledger_path)
    lock = _lock(path)
    try:
        rows = read_tsv(path)
        keys = (cohort, sample_id, modality, str(attempt))
        matches = [i for i, row in enumerate(rows) if
                   (row["cohort"], row["sample_id"], row["modality"], row["attempt"]) == keys]
        if append_attempt:
            if matches or attempt < 2:
                raise ValueError("new attempt must be unique and follow attempt 1")
            prior = next((r for r in rows if (r["cohort"], r["sample_id"], r["modality"], r["attempt"])
                          == (cohort, sample_id, modality, str(attempt - 1))), None)
            if not prior or prior["state"] not in {"failed", "superseded", "quarantined"}:
                raise ValueError("previous staging attempt has not failed")
            row = {field: "" for field in STAGE_LEDGER_FIELDS}
            for field in ("run_id", "run_role", "cohort", "sample_id", "donor_id", "modality",
                          "source_uri", "source_checksum_declared", "index_uri",
                          "index_checksum_declared"):
                row[field] = prior[field]
            row.update({"attempt": str(attempt), "state": "submitted",
                        "superseded_attempt": str(attempt - 1), "updated_at_utc": utc_now()})
            row.update({key: str(value) for key, value in updates.items()})
            rows.append(row)
        else:
            if len(matches) != 1:
                raise ValueError(f"expected one stage attempt, matched {len(matches)}")
            row = rows[matches[0]]
            if new_state not in STAGE_STATES.get(row["state"], set()):
                raise ValueError(f"invalid stage-ledger transition: {row['state']} -> {new_state}")
            row = {**row, **{key: str(value) for key, value in updates.items()},
                   "state": new_state, "updated_at_utc": utc_now()}
            if new_state == "running" and not row.get("started_at_utc"):
                row["started_at_utc"] = utc_now()
            if new_state == "validated":
                if row.get("exit_code") != "0" or not HEX64.fullmatch(row.get("extracted_sha256", "")):
                    raise ValueError("validated staging requires exit_code=0 and extracted SHA-256")
            if new_state == "failed" and not row.get("failure_class"):
                raise ValueError("failed staging requires failure_class")
            rows[matches[0]] = row
        rows.sort(key=lambda r: (r["cohort"], r["modality"], r["sample_id"], int(r["attempt"])))
        temporary = path.with_suffix(path.suffix + ".tmp")
        write_tsv(temporary, rows, list(STAGE_LEDGER_FIELDS))
        os.replace(temporary, path)
        return rows
    finally:
        lock.rmdir()


def validate_certificate_pair(identity: str | Path, certificate: str | Path,
                              now_epoch: int | None = None) -> dict:
    """Validate a user certificate without reading or reporting private key material."""
    identity, certificate = Path(identity), Path(certificate)
    if not identity.is_file() or not certificate.is_file():
        raise ValueError("identity and certificate files must exist")
    pub = subprocess.run(["ssh-keygen", "-y", "-f", str(identity)], capture_output=True,
                         text=True, check=True).stdout
    cert_pub = subprocess.run(["ssh-keygen", "-L", "-f", str(certificate)], capture_output=True,
                              text=True, check=True).stdout
    identity_fp = subprocess.run(["ssh-keygen", "-lf", "-", "-E", "sha256"], input=pub,
                                 capture_output=True, text=True, check=True).stdout.split()[1]
    cert_identity = re.search(r"Public key:.*?(SHA256:[A-Za-z0-9+/=]+)", cert_pub)
    validity = re.search(r"Valid: from (\S+) to (\S+)", cert_pub)
    if not cert_identity or not validity:
        raise ValueError("could not parse OpenSSH certificate")
    from datetime import datetime
    start = datetime.fromisoformat(validity.group(1)).timestamp()
    end = datetime.fromisoformat(validity.group(2)).timestamp()
    now = time.time() if now_epoch is None else now_epoch
    return {"passed": identity_fp == cert_identity.group(1) and start <= now <= end,
            "fingerprint_match": identity_fp == cert_identity.group(1),
            "valid_now": start <= now <= end, "valid_after_epoch": int(start),
            "valid_before_epoch": int(end),
            "certificate_sha256": hashlib.sha256(certificate.read_bytes()).hexdigest()}


def build_run_disposition(stage_ledgers: list[str | Path], caller_ledgers: list[str | Path],
                          output: str | Path) -> list[dict[str, str]]:
    """Render operational disposition only; intentionally exposes no accuracy fields."""
    fields = ("run_id", "run_role", "sample_id", "modality", "evidence_layer",
              "component", "attempt", "scheduler_job_id", "state", "failure_class",
              "exit_code", "superseded_attempt")
    rows: list[dict[str, str]] = []
    for path in stage_ledgers:
        for source in read_tsv(path):
            rows.append({field: value for field, value in {
                "run_id": source["run_id"], "run_role": source["run_role"],
                "sample_id": source["sample_id"], "modality": source["modality"],
                "evidence_layer": "staging", "component": "remote_cram_extraction",
                "attempt": source["attempt"], "scheduler_job_id": source["scheduler_job_id"],
                "state": source["state"], "failure_class": source["failure_class"],
                "exit_code": source["exit_code"],
                "superseded_attempt": source["superseded_attempt"],
            }.items() if field in fields})
    for path in caller_ledgers:
        for source in read_tsv(path):
            rows.append({field: value for field, value in {
                "run_id": source["run_id"], "run_role": source["run_role"],
                "sample_id": source["sample_id"], "modality": source["modality"],
                "evidence_layer": "caller", "component": source["caller"],
                "attempt": source["attempt"], "scheduler_job_id": source["job_id"],
                "state": source["state"], "failure_class": "",
                "exit_code": source["exit_code"], "superseded_attempt": "",
            }.items() if field in fields})
    rows.sort(key=lambda r: (r["run_id"], r["modality"], r["sample_id"],
                             r["evidence_layer"], r["component"], int(r["attempt"] or 0)))
    write_tsv(output, rows, list(fields))
    return rows


def finalize_sample(stage_ledger: str | Path, caller_ledger: str | Path,
                    cohort: str, sample_id: str, modality: str,
                    stage_scheduler_state: str, caller_scheduler_state: str,
                    completion_marker: str | Path) -> dict:
    """Reconcile terminal scheduler state and enforce the complete-sample contract."""
    from .roihu import reconcile_terminal_run_sample
    stage_rows = [row for row in read_tsv(stage_ledger)
                  if (row["cohort"], row["sample_id"], row["modality"]) ==
                  (cohort, sample_id, modality)]
    if not stage_rows:
        raise ValueError("missing staging ledger row")
    latest = max(stage_rows, key=lambda row: int(row["attempt"]))
    if latest["state"] in {"pending", "submitted", "running"}:
        target = ("cancelled" if latest["state"] == "submitted" and
                  stage_scheduler_state.upper().startswith("CANCELLED") else "failed")
        updates = {"exit_code": stage_scheduler_state or "terminal_unknown"}
        if target == "failed":
            updates["failure_class"] = "scheduler_terminal"
        transition_stage_attempt(stage_ledger, cohort, sample_id, modality,
                                 int(latest["attempt"]), target, **updates)
    caller_rows = [row for row in read_tsv(caller_ledger)
                   if (row["cohort"], row["sample_id"], row["modality"]) ==
                   (cohort, sample_id, modality)]
    active = {"pending", "submitted", "running", "succeeded", "resubmitted"}
    if any(row["state"] in active for row in caller_rows):
        reconcile_terminal_run_sample(caller_ledger, cohort, sample_id, modality,
                                      caller_scheduler_state or "terminal_unknown")
    latest = max(
        (row for row in read_tsv(stage_ledger)
         if (row["cohort"], row["sample_id"], row["modality"]) ==
         (cohort, sample_id, modality)), key=lambda row: int(row["attempt"]),
    )
    caller_rows = [row for row in read_tsv(caller_ledger)
                   if (row["cohort"], row["sample_id"], row["modality"]) ==
                   (cohort, sample_id, modality)]
    failures = []
    if latest["state"] != "validated":
        failures.append(f"stage={latest['state']}")
    if not caller_rows or any(row["state"] not in {"validated", "frozen"}
                              for row in caller_rows):
        failures.append("caller ledger incomplete")
    if not Path(completion_marker).is_file():
        failures.append("CALLERS_COMPLETE missing")
    return {"passed": not failures, "failures": failures,
            "stage_state": latest["state"], "caller_rows": len(caller_rows)}
