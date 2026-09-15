from __future__ import annotations

import csv
import re
import struct
from pathlib import Path

from .io import sha256, write_json


REQUIRED_FIELDS = {
    "figure_id", "panel_id", "title", "manuscript_role", "script",
    "script_sha256", "source_data", "source_data_sha256", "upstream_artifact",
    "upstream_sha256", "output_paths", "output_sha256", "width_in", "height_in",
    "resolution", "evidence_role", "generation_commit", "status", "blocker",
}
REQUIRED_MAIN = {"MAIN-1", "MAIN-2", "MAIN-3", "MAIN-4", "MAIN-5"}
HEX64 = re.compile(r"^[0-9a-f]{64}$")


def _rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _png_dimensions(path: Path) -> tuple[int, int]:
    with path.open("rb") as handle:
        header = handle.read(24)
    if len(header) != 24 or header[:8] != b"\x89PNG\r\n\x1a\n" or header[12:16] != b"IHDR":
        raise ValueError("invalid PNG header")
    return struct.unpack(">II", header[16:24])


def validate_figure_manifest(root: str | Path, manifest_path: str | Path,
                             output: str | Path | None = None) -> dict:
    project = Path(root).resolve()
    manifest = Path(manifest_path)
    if not manifest.is_absolute():
        manifest = project / manifest
    failures: list[str] = []
    rows = _rows(manifest) if manifest.is_file() else []
    if not rows:
        failures.append("figure manifest is missing or empty")
    elif set(rows[0]) != REQUIRED_FIELDS:
        failures.append("figure manifest columns differ from the frozen schema")
    ids = {row.get("figure_id", "") for row in rows}
    missing = sorted(REQUIRED_MAIN - ids)
    if missing:
        failures.append(f"required main figures missing: {missing}")
    blocked: list[str] = []
    generated = 0
    for row in rows:
        label = row.get("figure_id", "unknown")
        status = row.get("status", "")
        if status == "BLOCKED_EVIDENCE":
            blocked.append(label)
            if not row.get("blocker"):
                failures.append(f"{label}: blocked figure has no reason")
            if row.get("output_paths") or row.get("output_sha256"):
                failures.append(f"{label}: blocked figure must not claim outputs")
            continue
        if status != "GENERATED":
            failures.append(f"{label}: unsupported figure status {status!r}")
            continue
        generated += 1
        for path_field, hash_field in (
            ("script", "script_sha256"),
            ("source_data", "source_data_sha256"),
            ("upstream_artifact", "upstream_sha256"),
        ):
            relative = Path(row[path_field])
            candidate = project / relative
            expected = row[hash_field]
            if relative.is_absolute() or ".." in relative.parts or not candidate.is_file():
                failures.append(f"{label}: missing or unsafe {path_field}")
            elif not HEX64.fullmatch(expected) or sha256(candidate) != expected:
                failures.append(f"{label}: {path_field} checksum drift")
        outputs = row["output_paths"].split(";") if row["output_paths"] else []
        hashes = row["output_sha256"].split(";") if row["output_sha256"] else []
        if len(outputs) != 3 or len(hashes) != 3:
            failures.append(f"{label}: SVG, PDF, and PNG outputs are required")
            continue
        suffixes = {Path(value).suffix.lower() for value in outputs}
        if suffixes != {".svg", ".pdf", ".png"}:
            failures.append(f"{label}: output formats must be SVG, PDF, and PNG")
        for relative_text, expected in zip(outputs, hashes):
            relative = Path(relative_text)
            candidate = project / relative
            if relative.is_absolute() or ".." in relative.parts or not candidate.is_file():
                failures.append(f"{label}: output missing or unsafe: {relative_text}")
                continue
            if candidate.stat().st_size == 0:
                failures.append(f"{label}: empty output: {relative_text}")
            if not HEX64.fullmatch(expected) or sha256(candidate) != expected:
                failures.append(f"{label}: output checksum drift: {relative_text}")
            if candidate.suffix.lower() == ".png":
                try:
                    width, height = _png_dimensions(candidate)
                    expected_width = round(float(row["width_in"]) * 300)
                    expected_height = round(float(row["height_in"]) * 300)
                    if (width, height) != (expected_width, expected_height):
                        failures.append(f"{label}: PNG dimensions are not the declared 300-dpi size")
                except (ValueError, OSError):
                    failures.append(f"{label}: invalid PNG")
            elif candidate.suffix.lower() == ".pdf":
                if candidate.read_bytes()[:5] != b"%PDF-":
                    failures.append(f"{label}: invalid PDF")
            elif candidate.suffix.lower() == ".svg":
                text = candidate.read_text(encoding="utf-8", errors="replace")
                if "<svg" not in text or "viewBox" not in text:
                    failures.append(f"{label}: invalid SVG")
        if row.get("resolution") != "300 dpi PNG; vector SVG/PDF":
            failures.append(f"{label}: resolution declaration is not frozen")
        if not row.get("evidence_role") or not row.get("generation_commit"):
            failures.append(f"{label}: evidence role or generation commit missing")
    result = {
        "schema_version": "champhla-figure-audit-1",
        "passed": bool(rows) and not failures and not blocked and ids == REQUIRED_MAIN,
        "generated_rows_valid": bool(rows) and not failures,
        "source_hashes_valid": not any("checksum" in failure for failure in failures),
        "outputs_valid": not any(
            token in failure for failure in failures
            for token in ("output", "PNG", "PDF", "SVG", "resolution")
        ),
        "registry_agreement": not blocked,
        "generated_figures": generated,
        "blocked_figures": blocked,
        "failures": failures,
    }
    if output:
        write_json(output, result)
    return result
