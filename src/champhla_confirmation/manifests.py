from __future__ import annotations

import re

from .io import read_json
from .panels import MODALITIES, canonical_method
from .panels import PANELS


REQUIRED_COMPARATOR_FIELDS = {
    "method_id", "kind", "modalities", "version", "artifact_pin",
    "reference_build", "imgt_hla_version", "configuration",
    "deployable", "primary_family", "disposition_reason", "command",
    "threads", "memory", "artifact_sha256", "source_sha256",
    "reference_artifacts", "workflow_hash",
}
UNRESOLVED_MARKERS = {
    "unknown", "unresolved", "required", "verify", "tbd", "pending",
    "not established", "repository-wide",
}
HEX64 = re.compile(r"^[0-9a-f]{64}$")


def _has_unresolved(value: object) -> bool:
    text = str(value).strip().lower()
    return not text or any(marker in text for marker in UNRESOLVED_MARKERS)


def validate_comparator_manifest(path: str, require_frozen: bool = False) -> list[str]:
    """Validate neutral comparator inclusion and reproducibility metadata."""
    manifest = read_json(path)
    failures: list[str] = []
    if manifest.get("schema_version") != "champhla-comparator-manifest-2":
        failures.append("unsupported comparator manifest schema")
    if not manifest.get("selection_policy"):
        failures.append("comparator manifest lacks a neutral selection policy")
    rows = manifest.get("comparators", [])
    if not rows:
        failures.append("comparator manifest has no entries")
        return failures
    seen = set()
    for index, row in enumerate(rows, start=1):
        label = row.get("method_id") or f"row-{index}"
        missing = sorted(REQUIRED_COMPARATOR_FIELDS - set(row))
        if missing:
            failures.append(f"{label} missing fields: {missing}")
            continue
        method = canonical_method(str(row["method_id"]))
        if method in seen:
            failures.append(f"duplicate canonical comparator method: {method}")
        seen.add(method)
        modalities = row["modalities"]
        if not isinstance(modalities, list) or not modalities:
            failures.append(f"{label} has no modalities")
        elif set(modalities) - set(MODALITIES):
            failures.append(f"{label} has unsupported modalities")
        if row["primary_family"] and not row["deployable"]:
            failures.append(f"{label} is nondeployable but in the primary family")
        if not str(row["disposition_reason"]).strip():
            failures.append(f"{label} lacks a neutral inclusion/exclusion reason")
        if not isinstance(row["threads"], int) or row["threads"] < 0:
            failures.append(f"{label}.threads must be a non-negative integer")
        if not isinstance(row["reference_artifacts"], list):
            failures.append(f"{label}.reference_artifacts must be a list")
        if require_frozen and row["deployable"]:
            for field in ("version", "artifact_pin", "reference_build",
                          "imgt_hla_version", "configuration", "memory"):
                if _has_unresolved(row[field]):
                    failures.append(f"{label}.{field} is unresolved")
            if not str(row["command"]).strip():
                failures.append(f"{label}.command is empty")
            if row["threads"] < 1:
                failures.append(f"{label}.threads must be positive when deployable")
            for field in ("artifact_sha256", "source_sha256", "workflow_hash"):
                if not HEX64.fullmatch(str(row[field]).lower()):
                    failures.append(f"{label}.{field} is not an exact SHA-256")
            references = row["reference_artifacts"]
            if not references:
                failures.append(f"{label}.reference_artifacts is empty")
            for reference_index, reference in enumerate(references, start=1):
                if not isinstance(reference, dict):
                    failures.append(f"{label}.reference_artifacts[{reference_index}] is not an object")
                    continue
                if _has_unresolved(reference.get("artifact_id", "")):
                    failures.append(f"{label}.reference_artifacts[{reference_index}].artifact_id is unresolved")
                if not HEX64.fullmatch(str(reference.get("sha256", "")).lower()):
                    failures.append(f"{label}.reference_artifacts[{reference_index}].sha256 is not exact")
    if require_frozen and manifest.get("status") != "FROZEN":
        failures.append("comparator manifest status is not FROZEN")
    return sorted(set(failures))


def validate_caller_reference_attestation(path: str, require_ready: bool = True) -> list[str]:
    """Validate the caller/database evidence without substituting a global release."""
    data = read_json(path)
    failures: list[str] = []
    if data.get("schema_version") != "champhla-caller-reference-attestation-1":
        failures.append("unsupported caller-reference attestation schema")
    callers = data.get("callers", {})
    expected = {caller for panel in PANELS.values() for caller in panel}
    missing = sorted(expected - set(callers))
    extra = sorted(set(callers) - expected)
    if missing:
        failures.append(f"caller-reference attestation missing callers: {missing}")
    if extra:
        failures.append(f"caller-reference attestation has unknown callers: {extra}")
    for caller in sorted(expected & set(callers)):
        row = callers[caller]
        if not HEX64.fullmatch(str(row.get("artifact_sha256", "")).lower()):
            failures.append(f"{caller}: artifact SHA-256 is not exact")
        components = row.get("reference_components", [])
        if not components:
            failures.append(f"{caller}: no reference components")
        for index, component in enumerate(components, start=1):
            if not component.get("artifact_id"):
                failures.append(f"{caller}: reference component {index} has no identifier")
            if not HEX64.fullmatch(str(component.get("sha256", "")).lower()):
                failures.append(f"{caller}: reference component {index} has no exact SHA-256")
        if require_ready and not row.get("release_resolved"):
            failures.append(f"{caller}: database release unresolved")
    if require_ready and data.get("production_ready") is not True:
        failures.append("caller-reference attestation is not production ready")
    return sorted(set(failures))


def validate_hprc_truth_protocol(path: str, require_frozen: bool = True) -> list[str]:
    """Keep assembly truth execution blocked until every independent artifact is pinned."""
    data = read_json(path)
    failures: list[str] = []
    if data.get("schema_version") != "champhla-hprc-assembly-truth-protocol-2":
        failures.append("unsupported HPRC truth protocol schema")
    if data.get("target_subjects") != 120:
        failures.append("HPRC target must remain 120 subjects")
    access = data.get("assembly_access", {})
    if not str(access.get("official_archive", "")).startswith("https://"):
        failures.append("HPRC official archive is not HTTPS")
    if require_frozen and not HEX64.fullmatch(
            str(access.get("archive_artifact_sha256", "")).lower()):
        failures.append("HPRC AGC archive SHA-256 is unresolved")
    methods = data.get("methods", [])
    by_name = {row.get("method"): row for row in methods if isinstance(row, dict)}
    if set(by_name) != {"HLA-ASM", "Immuannot"}:
        failures.append("HPRC truth protocol must contain exactly HLA-ASM and Immuannot")
    for name in ("HLA-ASM", "Immuannot"):
        row = by_name.get(name, {})
        for field in ("version", "source_commit", "artifact_sha256", "wrapper_sha256"):
            value = row.get(field, "")
            if require_frozen and field == "version" and _has_unresolved(value):
                failures.append(f"{name}.{field} is unresolved")
            if require_frozen and field != "version" and not HEX64.fullmatch(str(value).lower()):
                failures.append(f"{name}.{field} is not an exact SHA-256")
        references = row.get("reference_artifacts", [])
        if not isinstance(references, list) or not references:
            failures.append(f"{name}.reference_artifacts is empty")
        elif require_frozen:
            for index, reference in enumerate(references, 1):
                if not reference.get("artifact_id"):
                    failures.append(f"{name}.reference_artifacts[{index}] has no identifier")
                if not HEX64.fullmatch(str(reference.get("sha256", "")).lower()):
                    failures.append(f"{name}.reference_artifacts[{index}] SHA-256 is unresolved")
    if require_frozen and data.get("status") != "FROZEN":
        failures.append("HPRC truth protocol is not FROZEN")
    return sorted(set(failures))
