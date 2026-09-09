from __future__ import annotations

import re
from pathlib import Path

from .imgt_release import (
    EVIDENCE_KINDS,
    NOT_RELEASE_BEARING,
    NO_MATCH,
    caller_release_summary,
)
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
    if require_frozen:
        failures.extend(_cross_check_attestation(manifest, path))
    return sorted(set(failures))


def _cross_check_attestation(manifest: dict, manifest_path: str) -> list[str]:
    """Confirm the manifest reports the release the attestation evidence derives.

    The two files are written at different times by different steps. Without this check a
    caller could be frozen in the manifest at one release while its component evidence says
    another, and nothing would notice.
    """
    reference = manifest.get("provenance_attestation")
    if not reference:
        return ["frozen comparator manifest has no provenance attestation reference"]
    attestation = Path(manifest_path).resolve().parent.parent / reference
    if not attestation.is_file():
        attestation = Path(reference)
    if not attestation.is_file():
        return [f"provenance attestation not found: {reference}"]
    derived = caller_release_map(str(attestation))
    failures = []
    for row in manifest.get("comparators", []):
        if row.get("kind") != "caller":
            continue
        caller = str(row["method_id"]).split(":", 1)[-1]
        if caller not in derived:
            failures.append(f"{row['method_id']} has no attestation record")
            continue
        if str(row.get("imgt_hla_version", "")).strip() != derived[caller]:
            failures.append(
                f"{row['method_id']}.imgt_hla_version {row.get('imgt_hla_version')!r} "
                f"disagrees with attestation evidence {derived[caller]!r}"
            )
    return failures


def validate_caller_reference_attestation(path: str, require_ready: bool = True) -> list[str]:
    """Validate caller/database evidence at the component level.

    Release identity belongs to a component, not to a caller. A deployed database may mix
    components built from different releases; that is a property of the deployment, not a
    failure to determine it, so such a caller is recorded as HETEROGENEOUS and is resolved.
    A component that carries no allele content -- an exon-split table, for example -- is
    marked NOT_RELEASE_BEARING and is excluded from the caller summary.

    Resolution is derived from the component evidence rather than read from a flag, so a
    caller cannot be declared resolved by editing a boolean.
    """
    data = read_json(path)
    failures: list[str] = []
    if data.get("schema_version") != "champhla-caller-reference-attestation-2":
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
        component_releases: list[str] = []
        for index, component in enumerate(components, start=1):
            if not component.get("artifact_id"):
                failures.append(f"{caller}: reference component {index} has no identifier")
            if not HEX64.fullmatch(str(component.get("sha256", "")).lower()):
                failures.append(f"{caller}: reference component {index} has no exact SHA-256")
            release = str(component.get("release", "")).strip()
            evidence = str(component.get("release_evidence", "")).strip()
            if not release:
                if require_ready:
                    failures.append(f"{caller}: component {index} has no release")
                continue
            component_releases.append(release)
            if release == NOT_RELEASE_BEARING:
                continue
            if release == NO_MATCH and require_ready:
                failures.append(f"{caller}: component {index} matched no release")
            if evidence not in EVIDENCE_KINDS:
                failures.append(
                    f"{caller}: component {index} release evidence is not one of {sorted(EVIDENCE_KINDS)}"
                )
        derived = caller_release_summary(component_releases)
        declared = str(row.get("database_release", "")).strip()
        if declared and declared != derived:
            failures.append(
                f"{caller}: database_release {declared!r} disagrees with component evidence {derived!r}"
            )
        resolved = derived not in {NO_MATCH, ""}
        if bool(row.get("release_resolved")) != resolved:
            failures.append(
                f"{caller}: release_resolved does not match the component evidence"
            )
        if require_ready and not resolved:
            failures.append(f"{caller}: database release unresolved")
    if require_ready and data.get("production_ready") is not True:
        failures.append("caller-reference attestation is not production ready")
    return sorted(set(failures))


def caller_release_map(path: str) -> dict[str, str]:
    """Return the derived release summary per caller, for cross-checking other manifests."""
    data = read_json(path)
    out: dict[str, str] = {}
    for caller, row in data.get("callers", {}).items():
        releases = [
            str(component.get("release", "")).strip()
            for component in row.get("reference_components", [])
        ]
        out[caller] = caller_release_summary([value for value in releases if value])
    return out


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
