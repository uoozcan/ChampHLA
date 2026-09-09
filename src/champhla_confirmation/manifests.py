from __future__ import annotations

from .io import read_json
from .panels import MODALITIES, canonical_method


REQUIRED_COMPARATOR_FIELDS = {
    "method_id", "kind", "modalities", "version", "artifact_pin",
    "reference_build", "imgt_hla_version", "configuration",
    "deployable", "primary_family", "disposition_reason",
}
UNRESOLVED_MARKERS = {"unknown", "required", "verify", "tbd", "pending"}


def validate_comparator_manifest(path: str, require_frozen: bool = False) -> list[str]:
    """Validate neutral comparator inclusion and reproducibility metadata."""
    manifest = read_json(path)
    failures: list[str] = []
    if manifest.get("schema_version") != "champhla-comparator-manifest-1":
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
        if require_frozen and row["deployable"]:
            for field in ("version", "artifact_pin", "reference_build",
                          "imgt_hla_version", "configuration"):
                value = str(row[field]).strip().lower()
                if not value or any(marker in value for marker in UNRESOLVED_MARKERS):
                    failures.append(f"{label}.{field} is unresolved")
    if require_frozen and manifest.get("status") != "FROZEN":
        failures.append("comparator manifest status is not FROZEN")
    return sorted(set(failures))
