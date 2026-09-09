from __future__ import annotations

import math
import os
import re
import shutil
import subprocess
import hashlib
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
import time
from pathlib import Path
from urllib.parse import urlparse
from urllib.request import urlopen

from .io import read_json, read_tsv, reject_truth_columns, sha256, write_json, write_tsv
from .panels import PANELS
from .parsers import parse_caller_call
from .manifests import validate_caller_reference_attestation


RUN_MANIFEST_FIELDS = (
    "cohort", "sample_id", "donor_id", "modality", "input_type", "input_uri",
    "index_uri", "index_checksum", "source_checksum", "reference_build", "read_layout",
    "independence_stratum", "evidence_role",
)
INPUT_TYPES = {"bam", "cram", "fastq_pair"}
EVIDENCE_ROLES = {
    "development", "same_resource_confirmation", "independent_validation", "exploratory",
}
INDEPENDENCE_STRATA = {"new_library_overlap", "donor_independent", "not_applicable"}
EXPECTED_REFERENCE = {"wgs": "GRCh38DH", "wes": "GRCh38DH", "rnaseq": "GRCh37"}
HEX64 = re.compile(r"^[0-9a-fA-F]{64}$")
HEX32 = re.compile(r"^[0-9a-fA-F]{32}$")
WORKFLOW_FORBIDDEN = (
    "errorStrategy 'ignore'",
    "? 'retry' : 'ignore'",
    "|| true",
    "No results generated",
    "(STUB)",
    "stub:",
)

RUN_LEDGER_FIELDS = (
    "cohort", "sample_id", "donor_id", "modality", "caller", "state", "job_id",
    "git_commit", "manifest_sha256", "exit_code", "runtime_seconds",
    "peak_memory_bytes", "peak_disk_bytes", "retained_disk_bytes", "output_sha256", "attempt",
    "supersedes_job_id", "updated_at_utc",
)
STATE_TRANSITIONS = {
    "pending": {"submitted"},
    "submitted": {"running", "failed"},
    "running": {"succeeded", "failed"},
    "succeeded": {"validated", "failed"},
    "failed": {"resubmitted", "quarantined"},
    "resubmitted": {"running", "failed"},
    "validated": {"frozen", "quarantined"},
    "frozen": set(),
    "quarantined": set(),
}


def _split(value: str) -> list[str]:
    return [item.strip() for item in value.split(";") if item.strip()]


def _public_uri(value: str) -> bool:
    parsed = urlparse(value)
    return parsed.scheme in {"http", "https", "ftp", "s3"} or Path(value).is_absolute()


def _valid_checksum(value: str) -> bool:
    token = value.strip()
    if token.startswith("sha256:"):
        token = token[7:]
        return bool(HEX64.fullmatch(token))
    if token.startswith("md5:"):
        token = token[4:]
        return bool(HEX32.fullmatch(token))
    return bool(HEX64.fullmatch(token) or HEX32.fullmatch(token))


def audit_run_manifest(path: str | Path) -> dict:
    rows = read_tsv(path)
    failures: list[str] = []
    if not rows:
        failures.append("run manifest is empty")
    else:
        try:
            reject_truth_columns(rows, "truth-free run manifest")
        except ValueError as error:
            failures.append(str(error))
        missing = sorted(set(RUN_MANIFEST_FIELDS) - set(rows[0]))
        if missing:
            failures.append(f"missing columns: {missing}")
    seen = set()
    modality_samples: dict[str, set[str]] = {key: set() for key in EXPECTED_REFERENCE}
    for line, row in enumerate(rows, 2):
        prefix = f"row {line}"
        key = (row.get("cohort", ""), row.get("sample_id", ""), row.get("modality", ""))
        if key in seen:
            failures.append(f"{prefix}: duplicate cohort/sample/modality {key}")
        seen.add(key)
        modality = row.get("modality", "").strip().lower()
        input_type = row.get("input_type", "").strip().lower()
        if modality not in EXPECTED_REFERENCE:
            failures.append(f"{prefix}: unsupported modality {modality!r}")
            continue
        modality_samples[modality].add(row.get("sample_id", ""))
        if not all(row.get(field, "").strip() for field in ("cohort", "sample_id", "donor_id")):
            failures.append(f"{prefix}: cohort, sample_id, and donor_id are required")
        if input_type not in INPUT_TYPES:
            failures.append(f"{prefix}: unsupported input_type {input_type!r}")
        uris, checksums = _split(row.get("input_uri", "")), _split(row.get("source_checksum", ""))
        expected_parts = 2 if input_type == "fastq_pair" else 1
        if len(uris) != expected_parts or len(checksums) != expected_parts:
            failures.append(
                f"{prefix}: {input_type or 'input'} requires {expected_parts} URI/checksum value(s)"
            )
        if any(not _public_uri(uri) for uri in uris):
            failures.append(f"{prefix}: input_uri is not a supported public or absolute URI")
        if any(not _valid_checksum(value) for value in checksums):
            failures.append(f"{prefix}: source_checksum must contain MD5 or SHA-256 values")
        if input_type in {"bam", "cram"} and not row.get("index_uri", "").strip():
            failures.append(f"{prefix}: indexed alignment input requires index_uri")
        if input_type in {"bam", "cram"} and not _valid_checksum(row.get("index_checksum", "")):
            failures.append(f"{prefix}: indexed alignment input requires an MD5 or SHA-256 index_checksum")
        if input_type == "fastq_pair" and row.get("index_uri", "").strip():
            failures.append(f"{prefix}: FASTQ pairs must not provide index_uri")
        if input_type == "fastq_pair" and row.get("index_checksum", "").strip():
            failures.append(f"{prefix}: FASTQ pairs must not provide index_checksum")
        if row.get("reference_build", "").strip() != EXPECTED_REFERENCE[modality]:
            failures.append(
                f"{prefix}: expected {EXPECTED_REFERENCE[modality]} for {modality}"
            )
        if row.get("read_layout", "").strip().lower() != "paired":
            failures.append(f"{prefix}: only paired inputs are frozen")
        if row.get("evidence_role", "").strip() not in EVIDENCE_ROLES:
            failures.append(f"{prefix}: unsupported evidence_role")
        if row.get("independence_stratum", "").strip() not in INDEPENDENCE_STRATA:
            failures.append(f"{prefix}: unsupported independence_stratum")
    return {
        "schema_version": "champhla-truth-free-run-manifest-audit-1",
        "passed": not failures,
        "truth_blind": not any("truth-bearing" in failure for failure in failures),
        "rows": len(rows),
        "samples_by_modality": {
            key: len(value) for key, value in modality_samples.items() if value
        },
        "source_sha256": sha256(path),
        "failures": failures,
    }


def freeze_run_manifest(path: str | Path, output: str | Path,
                        expected_counts: dict[str, int] | None = None) -> dict:
    result = audit_run_manifest(path)
    expected = {str(key): int(value) for key, value in (expected_counts or {}).items()}
    if expected:
        unknown = sorted(set(expected) - set(EXPECTED_REFERENCE))
        if unknown:
            result["failures"].append(f"unknown expected-count modalities: {unknown}")
        for modality, count in expected.items():
            observed = result["samples_by_modality"].get(modality, 0)
            if observed != count:
                result["failures"].append(
                    f"{modality}: expected {count} unique samples, observed {observed}"
                )
        result["passed"] = not result["failures"]
    if not result["passed"]:
        raise ValueError(f"run manifest failed: {result['failures']}")
    caller_locus_records = sum(
        expected.get(modality, 0) * len(PANELS[modality]) * 3
        for modality in expected
    ) if expected else None
    payload = {
        **result,
        "schema_version": "champhla-truth-free-run-manifest-freeze-1",
        "frozen": True,
        "path": Path(path).name,
        "expected_samples_by_modality": expected,
        "expected_caller_locus_records": caller_locus_records,
        "expected_plurality_rows": sum(expected.values()) * 3 if expected else None,
    }
    write_json(output, payload)
    return payload


def initialize_run_ledger(manifest_path: str | Path, output: str | Path,
                          git_commit: str = "", job_id: str = "") -> list[dict[str, str]]:
    audit = audit_run_manifest(manifest_path)
    if not audit["passed"]:
        raise ValueError(f"run manifest failed: {audit['failures']}")
    commit = git_commit.strip()
    if not commit:
        commit = subprocess.run(
            ["git", "rev-parse", "HEAD"], capture_output=True, text=True, check=False,
        ).stdout.strip() or "UNRESOLVED"
    rows = []
    for source in read_tsv(manifest_path):
        modality = source["modality"].lower()
        for caller in PANELS[modality]:
            rows.append({
                "cohort": source["cohort"], "sample_id": source["sample_id"],
                "donor_id": source["donor_id"], "modality": modality, "caller": caller,
                "state": "submitted" if job_id else "pending", "job_id": job_id,
                "git_commit": commit,
                "manifest_sha256": audit["source_sha256"], "exit_code": "",
                "runtime_seconds": "", "peak_memory_bytes": "", "peak_disk_bytes": "",
                "retained_disk_bytes": "", "output_sha256": "", "attempt": "1", "supersedes_job_id": "",
                "updated_at_utc": "",
            })
    rows.sort(key=lambda row: (row["cohort"], row["modality"], row["sample_id"], row["caller"]))
    write_tsv(output, rows, list(RUN_LEDGER_FIELDS))
    return rows


def transition_run_sample(ledger_path: str | Path, cohort: str, sample_id: str,
                          modality: str, new_state: str, **updates: str) -> list[dict[str, str]]:
    """Atomically transition every caller for one array sample on Roihu."""
    path = Path(ledger_path)
    lock_path = path.with_suffix(path.suffix + ".lockdir")
    deadline = time.monotonic() + 30
    while True:
        try:
            lock_path.mkdir()
            break
        except FileExistsError:
            if time.monotonic() >= deadline:
                raise TimeoutError(f"timed out acquiring run-ledger lock: {lock_path}")
            time.sleep(0.05)
    try:
        rows = read_tsv(path)
        matched = 0
        now = datetime.now(timezone.utc).isoformat().replace("+00:00", "Z")
        for index, row in enumerate(rows):
            if (row.get("cohort"), row.get("sample_id"), row.get("modality")) != (
                    cohort, sample_id, modality):
                continue
            rows[index] = transition_run_record(
                row, new_state, updated_at_utc=now, **updates,
            )
            matched += 1
        expected = len(PANELS.get(modality, ()))
        if matched != expected:
            raise ValueError(f"expected {expected} ledger rows, matched {matched}")
        temporary = path.with_suffix(path.suffix + ".tmp")
        write_tsv(temporary, rows, list(RUN_LEDGER_FIELDS))
        os.replace(temporary, path)
        return rows
    finally:
        lock_path.rmdir()


def transition_run_record(row: dict[str, str], new_state: str, **updates: str) -> dict[str, str]:
    old_state = row.get("state", "")
    if new_state not in STATE_TRANSITIONS.get(old_state, set()):
        raise ValueError(f"invalid run-ledger transition: {old_state} -> {new_state}")
    updated = {**row, **{key: str(value) for key, value in updates.items()}, "state": new_state}
    if new_state in {"succeeded", "validated", "frozen"}:
        if updated.get("exit_code") != "0" or not HEX64.fullmatch(updated.get("output_sha256", "")):
            raise ValueError(f"{new_state} requires exit_code=0 and an output SHA-256")
    if new_state == "resubmitted" and not updated.get("supersedes_job_id"):
        raise ValueError("resubmitted requires supersedes_job_id")
    return updated


def _nearest_rank(values: list[int], quantile: float) -> int:
    if not values:
        raise ValueError("no pilot disk measurements")
    ordered = sorted(values)
    return ordered[max(0, math.ceil(quantile * len(ordered)) - 1)]


def assess_storage(pilot_ledger: str | Path, targets: dict[str, int],
                   available_bytes: int, output: str | Path | None = None,
                   availability_source: str = "project_allocation",
                   environment_inventory_sha256: str = "") -> dict:
    rows = read_tsv(pilot_ledger)
    by_modality: dict[str, list[int]] = {}
    retained_by_modality: dict[str, list[int]] = {}
    for row in rows:
        if row.get("state") not in {"validated", "frozen"}:
            continue
        value = int(row.get("peak_disk_bytes", "0") or 0)
        if value > 0:
            by_modality.setdefault(row["modality"], []).append(value)
        retained = int(row.get("retained_disk_bytes", "0") or 0)
        if retained > 0:
            retained_by_modality.setdefault(row["modality"], []).append(retained)
    missing = sorted(
        modality for modality in targets
        if not by_modality.get(modality) or not retained_by_modality.get(modality)
    )
    p95 = {modality: _nearest_rank(values, 0.95) for modality, values in by_modality.items()}
    retained_p95 = {
        modality: _nearest_rank(values, 0.95)
        for modality, values in retained_by_modality.items()
    }
    projected = sum(p95.get(modality, 0) * int(count) for modality, count in targets.items())
    retained_projected = sum(
        retained_p95.get(modality, 0) * int(count) for modality, count in targets.items()
    )
    reserve = 100 * 1024 ** 3
    required = math.ceil(1.25 * projected) + reserve
    transient_p95 = {
        modality: max(0, p95.get(modality, 0) - retained_p95.get(modality, 0))
        for modality in targets
    }
    sequential_required = (
        retained_projected + math.ceil(1.25 * max(transient_p95.values(), default=0)) + reserve
    )
    source_valid = availability_source == "project_allocation"
    full_scale_passed = source_valid and not missing and int(available_bytes) >= required
    sequential_passed = (
        source_valid and not missing and int(available_bytes) >= sequential_required
    )
    execution_mode = (
        "full_scale" if full_scale_passed
        else "sequential_low_storage" if sequential_passed
        else "blocked"
    )
    result = {
        "schema_version": "champhla-roihu-storage-gate-2",
        "passed": full_scale_passed or sequential_passed,
        "full_scale_passed": full_scale_passed,
        "sequential_low_storage_passed": sequential_passed,
        "execution_mode": execution_mode,
        "pilot_p95_bytes_by_modality": p95,
        "pilot_retained_p95_bytes_by_modality": retained_p95,
        "pilot_transient_p95_bytes_by_modality": transient_p95,
        "targets": targets,
        "projected_peak_bytes": projected,
        "reserve_bytes": reserve,
        "required_available_bytes": required,
        "sequential_required_available_bytes": sequential_required,
        "observed_available_bytes": int(available_bytes),
        "availability_source": availability_source,
        "environment_inventory_sha256": environment_inventory_sha256,
        "missing_pilot_modalities": missing,
        "failures": ([] if source_valid else ["availability must come from project allocation"]),
        "policy": (
            "full scale: available >= 1.25 * summed peak + 100 GiB; "
            "sequential: available >= projected retained + 1.25 * largest transient + 100 GiB"
        ),
    }
    if output:
        write_json(output, result)
    return result


def build_cleanup_plan(run_root: str | Path, freeze_validation: str | Path,
                       output: str | Path) -> dict:
    root = Path(run_root).resolve()
    freeze = read_json(freeze_validation)
    failures = []
    if not root.is_dir() or root == Path("/"):
        failures.append("run root is missing or unsafe")
    if not freeze.get("valid"):
        failures.append("prediction freeze is not valid")
    candidates = []
    for relative in ("work", "tmp"):
        path = root / relative
        if path.is_dir():
            candidates.append({
                "path": path.relative_to(root).as_posix(),
                "absolute_path": path.as_posix(),
                "bytes": sum(item.stat().st_size for item in path.rglob("*") if item.is_file()),
            })
    result = {
        "schema_version": "champhla-roihu-cleanup-plan-1",
        "dry_run": True,
        "executable": not failures,
        "run_root": root.as_posix(),
        "eligible_paths": candidates,
        "protected_classes": [
            "source_inputs", "extracted_hla_inputs", "native_caller_outputs", "logs",
            "predictions", "reference_manifests", "release_evidence",
        ],
        "failures": failures,
    }
    write_json(output, result)
    return result


def directory_tree_identity(path: str | Path) -> dict:
    """Hash a directory deterministically using POSIX relative paths and file hashes."""
    root = Path(path).resolve()
    if not root.is_dir():
        return {"present": False, "path": root.as_posix(), "files": 0,
                "bytes": 0, "tree_sha256": "", "failures": ["directory missing"]}
    digest = hashlib.sha256()
    total = 0
    count = 0
    failures = []
    for item in sorted(root.rglob("*"), key=lambda value: value.relative_to(root).as_posix()):
        relative = item.relative_to(root).as_posix()
        if item.is_symlink():
            link_target = Path(os.readlink(item))
            resolved = item.resolve()
            try:
                target_relative = resolved.relative_to(root).as_posix()
            except ValueError:
                failures.append(f"symlink escapes frozen directory: {relative}")
                continue
            if link_target.is_absolute() or not resolved.is_file():
                failures.append(f"symlink is absolute or broken: {relative}")
                continue
            digest.update(relative.encode("utf-8"))
            digest.update(b"\0SYMLINK\0")
            digest.update(target_relative.encode("utf-8"))
            count += 1
            continue
        if not item.is_file():
            continue
        file_hash = sha256(item)
        digest.update(relative.encode("utf-8"))
        digest.update(b"\0")
        digest.update(bytes.fromhex(file_hash))
        total += item.stat().st_size
        count += 1
    return {
        "present": True, "path": root.as_posix(), "files": count, "bytes": total,
        "tree_sha256": digest.hexdigest(), "failures": failures,
    }


def validate_workflow_lock(lock_path: str | Path, project_root: str | Path) -> dict:
    """Validate every vendored execution file and reject masked-success constructs."""
    path = Path(lock_path)
    project = Path(project_root).resolve()
    failures = []
    if not path.is_file():
        return {"schema_version": "champhla-workflow-lock-audit-1", "passed": False,
                "failures": ["workflow lock missing"], "files": 0}
    lock = read_json(path)
    if lock.get("schema_version") != "champhla-workflow-lock-1":
        failures.append("unsupported workflow lock schema")
    if lock.get("status") != "FROZEN":
        failures.append("workflow lock is not FROZEN")
    files = lock.get("files", {})
    if not isinstance(files, dict) or not files:
        failures.append("workflow lock has no files")
        files = {}
    for relative, expected in sorted(files.items()):
        candidate = Path(relative)
        if candidate.is_absolute() or ".." in candidate.parts:
            failures.append(f"unsafe workflow path: {relative}")
            continue
        source = (project / candidate).resolve()
        try:
            source.relative_to(project)
        except ValueError:
            failures.append(f"workflow path escapes project: {relative}")
            continue
        if not source.is_file():
            failures.append(f"workflow file missing: {relative}")
            continue
        observed = sha256(source)
        if not HEX64.fullmatch(str(expected)) or observed != expected:
            failures.append(f"workflow hash mismatch: {relative}")
        if source.suffix in {".nf", ".config"}:
            content = source.read_text(encoding="utf-8", errors="strict")
            for token in WORKFLOW_FORBIDDEN:
                if token in content:
                    failures.append(f"workflow forbidden token {token!r}: {relative}")
    required = {
        "workflow/main.nf", "workflow/nextflow.config", "workflow/conf/roihu_params.yaml",
        *{f"workflow/modules/{name}.nf" for name in (
            "arcashla", "bam_to_fastq", "hlahd", "kourami", "optitype",
            "polysolver", "spechla", "t1k",
        )},
    }
    missing_required = sorted(required - set(files))
    if missing_required:
        failures.append(f"required workflow files absent from lock: {missing_required}")
    return {
        "schema_version": "champhla-workflow-lock-audit-1",
        "passed": not failures, "lock_sha256": sha256(path),
        "files": len(files), "failures": failures,
    }


def freeze_workflow_lock(project_root: str | Path, workflow_root: str | Path,
                         output: str | Path) -> dict:
    """Create the immutable lock for the compact, repository-owned workflow."""
    project = Path(project_root).resolve()
    workflow = Path(workflow_root).resolve()
    workflow.relative_to(project)
    selected = sorted(
        item for item in workflow.rglob("*")
        if item.is_file() and item.suffix in {".nf", ".config", ".yaml", ".py"}
    )
    files = {item.relative_to(project).as_posix(): sha256(item) for item in selected}
    payload = {
        "schema_version": "champhla-workflow-lock-1", "status": "FROZEN",
        "workflow_root": workflow.relative_to(project).as_posix(), "files": files,
        "policy": "exact files; POSIX paths; masked-success constructs forbidden",
    }
    write_json(output, payload)
    result = validate_workflow_lock(output, project)
    if not result["passed"]:
        Path(output).unlink(missing_ok=True)
        raise ValueError(f"workflow cannot be frozen: {result['failures']}")
    return payload


def _version(command: str, args: list[str]) -> str:
    path = shutil.which(command)
    if not path:
        return "MISSING"
    completed = subprocess.run([path, *args], capture_output=True, text=True, check=False)
    text = (completed.stdout or completed.stderr).strip().splitlines()
    return text[0] if text else f"exit={completed.returncode}"


def _du_bytes(path: Path) -> int:
    completed = subprocess.run(
        ["du", "-s", "-B1", str(path)], capture_output=True, text=True, check=False,
    )
    if completed.returncode != 0 or not completed.stdout.strip():
        raise RuntimeError(f"project allocation usage unavailable: {completed.stderr.strip()}")
    return int(completed.stdout.split()[0])


def _repository_identity(project_root: Path) -> dict:
    commit = subprocess.run(
        ["git", "-C", str(project_root), "rev-parse", "HEAD"],
        capture_output=True, text=True, check=False,
    ).stdout.strip()
    status = subprocess.run(
        ["git", "-C", str(project_root), "status", "--porcelain"],
        capture_output=True, text=True, check=False,
    ).stdout.splitlines()
    tracked = subprocess.run(
        ["git", "-C", str(project_root), "ls-files", "-z"],
        capture_output=True, check=False,
    ).stdout.split(b"\0")
    digest = hashlib.sha256()
    missing = []
    for raw in sorted(item for item in tracked if item):
        relative = raw.decode("utf-8", errors="strict")
        path = project_root / relative
        if not path.is_file():
            missing.append(relative)
            continue
        digest.update(relative.replace("\\", "/").encode("utf-8"))
        digest.update(b"\0")
        digest.update(bytes.fromhex(sha256(path)))
    return {
        "commit": commit or "UNRESOLVED",
        "clean": not status,
        "status": status,
        "tracked_tree_sha256": digest.hexdigest(),
        "missing_tracked_files": missing,
    }


def inventory_environment(site_config: str | Path, output: str | Path) -> dict:
    config = read_json(site_config)
    roots = {name: Path(value).resolve() for name, value in config["roots"].items()}
    files = {}
    failures = []
    for group in ("references", "containers", "software"):
        files[group] = {}
        for name, raw_path in config.get(group, {}).items():
            path = Path(raw_path)
            present = path.is_file()
            files[group][name] = {
                "path": path.as_posix(), "present": present,
                "bytes": path.stat().st_size if present else 0,
                "sha256": sha256(path) if present else "",
            }
            if not present:
                failures.append(f"{group}:{name}:missing")
    files["reference_directories"] = {}
    for name, raw_path in config.get("reference_directories", {}).items():
        record = directory_tree_identity(raw_path)
        files["reference_directories"][name] = record
        failures.extend(f"reference_directories:{name}:{failure}"
                        for failure in record["failures"])
        if not record["present"]:
            failures.append(f"reference_directories:{name}:missing")
    project_storage_config = config.get("project_storage", {})
    allocation_root = Path(project_storage_config.get("allocation_root", "")).resolve()
    allocation_bytes = int(project_storage_config.get("allocation_bytes", 0) or 0)
    try:
        project_used = _du_bytes(allocation_root) if allocation_root.is_dir() else 0
    except RuntimeError as error:
        project_used = 0
        failures.append(str(error))
    if not allocation_root.is_dir():
        failures.append("project_storage:allocation_root_missing")
    if allocation_bytes <= 0:
        failures.append("project_storage:allocation_bytes_missing")
    if project_used > allocation_bytes:
        failures.append("project_storage:usage_exceeds_declared_allocation")
    project_free = max(0, allocation_bytes - project_used)
    tools = {
        "python": _version("python3", ["--version"]),
        "samtools": _version("samtools", ["--version"]),
        "nextflow": _version("nextflow", ["-version"]),
        "apptainer": _version("apptainer", ["--version"]),
        "sbatch": _version("sbatch", ["--version"]),
    }
    wgs_artifacts = {}
    for caller in PANELS["wgs"]:
        group = "software" if caller == "SpecHLA" else "containers"
        wgs_artifacts[caller] = files.get(group, {}).get(caller, {}).get("sha256", "")
    reference_record = files.get("references", {}).get("GRCh38DH", {})
    repository = _repository_identity(roots["project_root"])
    workflow_lock = validate_workflow_lock(
        config.get("workflow_lock", ""), roots["project_root"],
    )
    if not workflow_lock["passed"]:
        failures.extend(f"workflow_lock:{failure}" for failure in workflow_lock["failures"])
    for name, version in tools.items():
        if version == "MISSING" or version.startswith("exit="):
            failures.append(f"tool:{name}:missing_or_unusable")
    expected_python = str(config.get("production_python", ""))
    if expected_python and expected_python not in tools["python"]:
        failures.append(f"python:expected_{expected_python}:observed_{tools['python']}")
    if not repository["clean"] or repository["missing_tracked_files"]:
        failures.append("repository:not_clean_or_complete")
    attestation_path = config.get("caller_reference_attestation", "")
    if not Path(attestation_path).is_file():
        caller_reference_failures = ["caller-reference attestation missing"]
    else:
        caller_reference_failures = validate_caller_reference_attestation(
            attestation_path, require_ready=True,
        )
    failures.extend(f"caller_reference:{failure}" for failure in caller_reference_failures)
    designs = config.get("reference_designs", {})
    for modality, expected in EXPECTED_REFERENCE.items():
        design = designs.get(modality, {})
        if design.get("source_build") != expected:
            failures.append(f"reference_design:{modality}:expected_{expected}")
        runtime_name = design.get("runtime_fasta")
        if runtime_name and not files.get("references", {}).get(runtime_name, {}).get("present"):
            failures.append(f"reference_design:{modality}:runtime_fasta_missing")
    result = {
        "schema_version": "champhla-roihu-environment-inventory-1",
        "passed": not failures,
        "hostname": os.uname().nodename,
        "tools": tools,
        "roots": {name: path.as_posix() for name, path in roots.items()},
        "repository": repository,
        "project_storage": {
            "allocation_root": allocation_root.as_posix(),
            "allocation_bytes": allocation_bytes,
            "used_bytes": project_used,
            "free_bytes": project_free,
            "accounting": project_storage_config.get("accounting", ""),
            "filesystem_global_space_ignored": True,
        },
        "artifacts": files,
        "reference_build": "GRCh38DH",
        "reference_sha256": reference_record.get("sha256", ""),
        "caller_reference_attestation": {
            "path": str(attestation_path),
            "sha256": sha256(attestation_path) if Path(attestation_path).is_file() else "",
            "passed": not caller_reference_failures,
            "failures": caller_reference_failures,
        },
        "python_version": tools["python"],
        "samtools_version": tools["samtools"],
        "caller_artifacts": wgs_artifacts,
        "reference_designs": designs,
        "workflow_lock": workflow_lock,
        "failures": failures,
    }
    write_json(output, result)
    return result


CALLER_SUFFIX = {
    "HLA-HD": "hlahd", "Kourami": "kourami", "OptiType": "optitype",
    "POLYSOLVER": "polysolver", "SpecHLA": "spechla", "T1K": "t1k",
    "ArcasHLA": "arcashla",
}


def collect_run_outputs(caller_root: str | Path, manifest_path: str | Path,
                        output: str | Path, summary_output: str | Path) -> dict:
    manifest_audit = audit_run_manifest(manifest_path)
    if not manifest_audit["passed"]:
        raise ValueError(f"run manifest failed: {manifest_audit['failures']}")
    root = Path(caller_root)
    output_rows = []
    failures = []
    for sample in read_tsv(manifest_path):
        modality = sample["modality"]
        sample_id = sample["sample_id"]
        sample_root = root / modality / sample_id
        if not (sample_root / "CALLERS_COMPLETE").is_file():
            failures.append(f"{modality}:{sample_id}:completion_marker_missing")
            continue
        validation_path = sample_root / "caller_output_validation.tsv"
        validation = read_tsv(validation_path) if validation_path.is_file() else []
        for caller in PANELS[modality]:
            suffix = CALLER_SUFFIX[caller]
            matches = list((sample_root / "results").rglob(f"{sample_id}_{suffix}.txt"))
            if len(matches) != 1 or matches[0].stat().st_size == 0:
                failures.append(f"{modality}:{sample_id}:{caller}:native_output_count={len(matches)}")
                continue
            source = matches[0]
            source_hash = sha256(source)
            validation_rows = [row for row in validation if row.get("caller") == suffix]
            if len(validation_rows) != 1 or validation_rows[0].get("sha256") != source_hash:
                failures.append(f"{modality}:{sample_id}:{caller}:validation_hash_mismatch")
                continue
            native = source.read_text(encoding="utf-8", errors="replace")
            for gene in ("A", "B", "C"):
                try:
                    call = parse_caller_call(caller, native, gene)
                except ValueError as error:
                    failures.append(f"{modality}:{sample_id}:{caller}:{gene}:{error}")
                    continue
                output_rows.append({
                    "cohort": sample["cohort"], "subject": sample_id,
                    "donor": sample["donor_id"],
                    "independence_stratum": sample["independence_stratum"],
                    "evidence_role": sample["evidence_role"], "modality": modality,
                    "gene": gene, "caller": caller, "allele1_raw": call["allele1"],
                    "allele2_raw": call["allele2"], "allele1": call["allele1"],
                    "allele2": call["allele2"], "call_status": call["call_status"],
                    "source_path": source.resolve().as_posix(), "source_sha256": source_hash,
                })
    expected = sum(3 * len(PANELS[row["modality"]]) for row in read_tsv(manifest_path))
    passed = not failures and len(output_rows) == expected
    write_tsv(output, sorted(
        output_rows, key=lambda row: (row["cohort"], row["modality"], row["subject"],
                                      row["gene"], row["caller"]),
    ))
    summary = {
        "schema_version": "champhla-roihu-native-output-collection-1",
        "passed": passed, "truth_blind": True, "expected_records": expected,
        "observed_records": len(output_rows), "manifest_sha256": sha256(manifest_path),
        "output_sha256": sha256(output), "failures": failures,
        "partial_records": sum(row["call_status"] == "partial" for row in output_rows),
        "missing_records": sum(row["call_status"] == "missing" for row in output_rows),
    }
    write_json(summary_output, summary)
    return summary


def build_nci60_run_manifest(pilot_path: str | Path, ena_report_path: str | Path,
                             output: str | Path) -> dict:
    """Build the 11-subject exploratory lane from truth-free ENA metadata."""
    pilot = read_tsv(pilot_path)
    reject_truth_columns(pilot, "NCI-60 pilot roster")
    ready = {}
    for row in pilot:
        if (row.get("cohort") == "NCI60_PUBLIC_PILOT"
                and row.get("modality") == "rnaseq" and row.get("gate_status") == "ready"):
            ready[row["subject"]] = row
    if len(ready) != 11:
        raise ValueError(f"expected 11 strict NCI-60 RNA subjects, observed {len(ready)}")

    ena = {row.get("run_accession", ""): row for row in read_tsv(ena_report_path)}
    rows = []
    for subject, source in sorted(ready.items()):
        accession = source.get("source_accession", "")
        metadata = ena.get(accession)
        if not metadata:
            raise ValueError(f"{subject}: ENA checksum metadata missing for {accession}")
        urls = _split(metadata.get("fastq_ftp", ""))
        checksums = _split(metadata.get("fastq_md5", ""))
        if len(urls) != 2 or len(checksums) != 2:
            raise ValueError(f"{subject}: expected two ENA FASTQ URLs and MD5 values")
        urls = [url if "://" in url else f"https://{url}" for url in urls]
        rows.append({
            "cohort": "NCI60", "sample_id": accession, "donor_id": subject,
            "modality": "rnaseq", "input_type": "fastq_pair",
            "input_uri": ";".join(urls), "index_uri": "", "index_checksum": "",
            "source_checksum": ";".join(f"md5:{value}" for value in checksums),
            "reference_build": "GRCh37", "read_layout": "paired",
            "independence_stratum": "donor_independent", "evidence_role": "exploratory",
        })
    write_tsv(output, rows, list(RUN_MANIFEST_FIELDS))
    audit = audit_run_manifest(output)
    if not audit["passed"]:
        raise ValueError(f"generated NCI-60 manifest failed: {audit['failures']}")
    return audit


def resolve_same_resource_index_checksums(roster_path: str | Path,
                                          assay_manifest_path: str | Path,
                                          output: str | Path, workers: int = 8) -> dict:
    """Stream only missing public CRAI files and freeze both MD5 and SHA-256."""
    roster = read_tsv(roster_path)
    reject_truth_columns(roster, "same-resource truth-free roster")
    keys = {(row.get("subject", ""), row.get("modality", "").lower()) for row in roster}
    targets = []
    for row in read_tsv(assay_manifest_path):
        key = (row.get("sample_id", ""), row.get("modality", "").lower())
        if key not in keys or key[1] == "rnaseq" or row.get("index_md5", ""):
            continue
        uri = row.get("index_url", "")
        if not uri:
            raise ValueError(f"{key}: missing public index URI")
        targets.append((key, uri))

    def resolve(target: tuple[tuple[str, str], str]) -> dict[str, object]:
        key, uri = target
        md5 = hashlib.md5()
        sha = hashlib.sha256()
        total = 0
        with urlopen(uri, timeout=120) as response:
            while chunk := response.read(1024 * 1024):
                md5.update(chunk)
                sha.update(chunk)
                total += len(chunk)
        if total == 0:
            raise ValueError(f"{key}: empty public index response")
        return {
            "sample_id": key[0], "modality": key[1], "index_uri": uri,
            "bytes": total, "md5": md5.hexdigest(), "sha256": sha.hexdigest(),
        }

    if workers < 1 or workers > 16:
        raise ValueError("index checksum workers must be between 1 and 16")
    with ThreadPoolExecutor(max_workers=workers) as executor:
        records = list(executor.map(resolve, targets))
    records.sort(key=lambda row: (row["modality"], row["sample_id"]))
    write_tsv(output, records, ["sample_id", "modality", "index_uri", "bytes", "md5", "sha256"])
    return {"records": len(records), "output_sha256": sha256(output), "truth_blind": True}


def build_same_resource_run_manifest(roster_path: str | Path, assay_manifest_path: str | Path,
                                     ena_report_path: str | Path,
                                     index_checksums_path: str | Path, output: str | Path) -> dict:
    """Build the exact corrected 1000G run manifest without reading joined truth."""
    roster = read_tsv(roster_path)
    reject_truth_columns(roster, "same-resource truth-free roster")
    expected = {"wgs": 137, "wes": 130, "rnaseq": 107}
    observed = {
        modality: len({row.get("subject", "") for row in roster
                       if row.get("modality", "").lower() == modality})
        for modality in expected
    }
    if observed != expected:
        raise ValueError(f"same-resource roster counts differ: expected {expected}, observed {observed}")
    roster_keys = [(row.get("subject", ""), row.get("modality", "").lower()) for row in roster]
    if len(roster_keys) != len(set(roster_keys)):
        raise ValueError("same-resource roster contains duplicate subject/modality rows")

    assay_by_key: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
    for row in read_tsv(assay_manifest_path):
        assay_by_key[(row.get("sample_id", ""), row.get("modality", "").lower())].append(row)
    ena_by_run = {row.get("run_accession", ""): row for row in read_tsv(ena_report_path)}
    ena_by_subject: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in ena_by_run.values():
        uris = _split(row.get("submitted_ftp", ""))
        if uris:
            basename = Path(urlparse(uris[0]).path).name
            ena_by_subject[basename.split(".", 1)[0]].append(row)
    index_checksums = {
        (row.get("sample_id", ""), row.get("modality", "").lower()): row
        for row in read_tsv(index_checksums_path)
    }
    rows = []
    for source in sorted(roster, key=lambda row: (row["modality"], row["subject"])):
        subject = source["subject"]
        modality = source["modality"].lower()
        candidates = assay_by_key.get((subject, modality), [])
        if modality != "rnaseq" and len(candidates) != 1:
            raise ValueError(f"{subject}:{modality}: expected one official assay record, observed {len(candidates)}")
        assay = candidates[0] if candidates else {}
        if modality == "rnaseq":
            accession = assay.get("source_record", "")
            ena = ena_by_run.get(accession) if accession else None
            if not ena:
                subject_runs = ena_by_subject.get(subject, [])
                if len(subject_runs) != 1:
                    raise ValueError(f"{subject}: expected one ENA RNA run, observed {len(subject_runs)}")
                ena = subject_runs[0]
            uris = _split(ena.get("submitted_ftp", ""))
            checksums = _split(ena.get("submitted_md5", ""))
            if len(uris) != 2 or len(checksums) != 2:
                raise ValueError(f"{subject}: ENA submitted FASTQ pair/checksums incomplete")
            input_type = "fastq_pair"
            input_uri = ";".join(uri if "://" in uri else f"https://{uri}" for uri in uris)
            source_checksum = ";".join(f"md5:{value}" for value in checksums)
            index_uri = index_checksum = ""
        else:
            if assay.get("format", "").lower() != "cram":
                raise ValueError(f"{subject}:{modality}: expected full CRAM source")
            input_type = "cram"
            input_uri = assay.get("input_url", "")
            index_uri = assay.get("index_url", "")
            source_checksum = f"md5:{assay.get('input_md5', '')}"
            index_md5 = assay.get("index_md5", "")
            if not index_md5:
                resolved = index_checksums.get((subject, modality), {})
                if resolved.get("index_uri") != index_uri:
                    raise ValueError(f"{subject}:{modality}: frozen index checksum missing or URI differs")
                index_md5 = resolved.get("md5", "")
            index_checksum = f"md5:{index_md5}"
        rows.append({
            "cohort": "1000G_SAME_RESOURCE", "sample_id": subject, "donor_id": subject,
            "modality": modality, "input_type": input_type, "input_uri": input_uri,
            "index_uri": index_uri, "index_checksum": index_checksum,
            "source_checksum": source_checksum,
            "reference_build": EXPECTED_REFERENCE[modality], "read_layout": "paired",
            "independence_stratum": "new_library_overlap",
            "evidence_role": "same_resource_confirmation",
        })
    write_tsv(output, rows, list(RUN_MANIFEST_FIELDS))
    audit = audit_run_manifest(output)
    if not audit["passed"]:
        raise ValueError(f"generated same-resource manifest failed: {audit['failures']}")
    return audit


def build_hprc_run_manifest(roster_path: str | Path, output: str | Path) -> dict:
    """Convert a frozen 120-subject HPRC roster into a checksummed WGS manifest."""
    roster = read_tsv(roster_path)
    reject_truth_columns(roster, "HPRC short-read roster")
    if len(roster) != 120:
        raise ValueError(f"expected frozen 120-subject HPRC roster, observed {len(roster)}")
    rows = []
    for source in roster:
        sample = source.get("sample_id", "") or source.get("subject", "")
        reference = source.get("reference_build", "") or source.get("genome_build", "")
        if reference != "GRCh38DH":
            raise ValueError(f"{sample}: GRCh38DH header verification is not complete")
        rows.append({
            "cohort": "HPRC_R2", "sample_id": sample, "donor_id": sample,
            "modality": "wgs", "input_type": "cram",
            "input_uri": source.get("input_uri", "") or source.get("input_url", ""),
            "index_uri": source.get("index_uri", "") or source.get("index_url", ""),
            "source_checksum": source.get("source_checksum", "") or source.get("cram_md5", ""),
            "index_checksum": source.get("index_checksum", "") or source.get("index_md5", ""),
            "reference_build": reference, "read_layout": "paired",
            "independence_stratum": "donor_independent",
            "evidence_role": "independent_validation",
        })
    write_tsv(output, rows, list(RUN_MANIFEST_FIELDS))
    audit = audit_run_manifest(output)
    if not audit["passed"]:
        raise ValueError(f"generated HPRC manifest failed: {audit['failures']}")
    return audit
