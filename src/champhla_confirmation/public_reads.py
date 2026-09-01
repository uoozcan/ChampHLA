from __future__ import annotations

import csv
import hashlib
import io
import json
import time
import urllib.parse
import urllib.request
from pathlib import Path

from .io import read_tsv, sha256, write_json, write_tsv


ENA_ENDPOINT = "https://www.ebi.ac.uk/ena/portal/api/search"
ENA_FIELDS = [
    "run_accession", "study_accession", "sample_accession", "experiment_accession",
    "experiment_title", "sample_title", "scientific_name", "tax_id", "library_strategy",
    "library_source", "library_layout", "instrument_platform", "base_count", "read_count",
    "fastq_ftp", "submitted_ftp",
]
ALLOWED_REGISTRY_FIELDS = ["ihw_number", "primary_name", "ancestry"]


def sanitize_ihwg_registry(source_tsv: str, output_tsv: str, manifest_json: str) -> dict:
    """Remove genotype-derived columns before public-read discovery."""
    rows = read_tsv(source_tsv)
    if not rows:
        raise ValueError("IHWG registry is empty")
    missing = [field for field in ALLOWED_REGISTRY_FIELDS if field not in rows[0]]
    if missing:
        raise ValueError(f"IHWG registry missing columns: {missing}")
    clean = []
    seen = set()
    for row in rows:
        ihw = row["ihw_number"].strip()
        name = row["primary_name"].strip()
        if not ihw or not name or ihw in seen:
            continue
        seen.add(ihw)
        clean.append({field: row.get(field, "").strip() for field in ALLOWED_REGISTRY_FIELDS})
    clean.sort(key=lambda row: row["ihw_number"])
    write_tsv(output_tsv, clean, ALLOWED_REGISTRY_FIELDS)
    manifest = {
        "schema_version": "ihwg-public-read-registry-1",
        "truth_blind": True,
        "allowed_columns": ALLOWED_REGISTRY_FIELDS,
        "discarded_columns": sorted(set(rows[0]) - set(ALLOWED_REGISTRY_FIELDS)),
        "source_sha256": sha256(source_tsv),
        "output_sha256": sha256(output_tsv),
        "subjects": len(clean),
    }
    write_json(manifest_json, manifest)
    return manifest


def _query_url(primary_name: str, limit: int) -> str:
    safe_name = primary_name.replace("\\", "\\\\").replace('"', '\\"')
    parameters = {
        "result": "read_run",
        "query": f'sample_title="{safe_name}"',
        "fields": ",".join(ENA_FIELDS),
        "format": "tsv",
        "limit": str(limit),
    }
    return ENA_ENDPOINT + "?" + urllib.parse.urlencode(parameters)


def _parse_ena_tsv(text: str) -> list[dict[str, str]]:
    if not text.strip():
        return []
    return list(csv.DictReader(io.StringIO(text), delimiter="\t"))


def _classify(row: dict[str, str], primary_name: str) -> dict[str, str]:
    strategy = row.get("library_strategy", "").strip().upper().replace("_", "-")
    source = row.get("library_source", "").strip().upper()
    layout = row.get("library_layout", "").strip().upper()
    platform = row.get("instrument_platform", "").strip().upper()
    organism = row.get("scientific_name", "").strip().lower()
    modality = ""
    if strategy in {"WXS", "WES"} and source == "GENOMIC":
        modality = "wes"
    elif strategy in {"RNA-SEQ", "RNASEQ"} and source == "TRANSCRIPTOMIC":
        modality = "rnaseq"
    technical_match = bool(
        modality and layout == "PAIRED" and platform == "ILLUMINA"
        and organism in {"homo sapiens", "human"}
    )
    exact_title = row.get("sample_title", "").strip().casefold() == primary_name.strip().casefold()
    return {
        "modality": modality,
        "exact_sample_title": "1" if exact_title else "0",
        "technical_match": "1" if technical_match else "0",
        "provenance_status": "requires_study_level_review" if technical_match and exact_title else "ineligible",
    }


def scout_ihwg_public_reads(
    registry_tsv: str,
    output_tsv: str,
    summary_json: str,
    cache_dir: str,
    sleep_seconds: float = 0.1,
    limit: int = 1000,
    offline: bool = False,
) -> dict:
    """Cross-match IHWG names against ENA without using HLA truth.

    Exact short-name matches are candidates only. Every technically plausible
    run remains locked as ``requires_study_level_review`` until its publication,
    BioSample and cell-line provenance are independently verified.
    """
    registry = read_tsv(registry_tsv)
    if not registry:
        raise ValueError("sanitized IHWG registry is empty")
    unexpected = sorted(set(registry[0]) - set(ALLOWED_REGISTRY_FIELDS))
    if unexpected:
        raise ValueError(f"registry is not sanitized; unexpected columns: {unexpected}")
    cache_root = Path(cache_dir)
    cache_root.mkdir(parents=True, exist_ok=True)
    output = []
    failures = []
    for index, subject in enumerate(registry):
        name = subject["primary_name"].strip()
        cache_key = hashlib.sha256(name.casefold().encode("utf-8")).hexdigest()
        cache_path = cache_root / f"{cache_key}.tsv"
        if cache_path.exists():
            text = cache_path.read_text(encoding="utf-8", errors="replace")
        elif offline:
            failures.append({"ihw_number": subject["ihw_number"], "reason": "cache_missing"})
            continue
        else:
            try:
                request = urllib.request.Request(
                    _query_url(name, limit), headers={"User-Agent": "ChampHLA-CC-confirmation/0.1"}
                )
                with urllib.request.urlopen(request, timeout=90) as response:
                    text = response.read().decode("utf-8", errors="replace")
                cache_path.write_text(text, encoding="utf-8")
            except Exception as error:  # the failure is recorded and the search remains resumable
                failures.append({"ihw_number": subject["ihw_number"], "reason": str(error)})
                continue
            if sleep_seconds and index + 1 < len(registry):
                time.sleep(sleep_seconds)
        for run in _parse_ena_tsv(text):
            classified = _classify(run, name)
            output.append({
                "ihw_number": subject["ihw_number"],
                "primary_name": name,
                "ancestry": subject.get("ancestry", ""),
                **{field: run.get(field, "") for field in ENA_FIELDS},
                **classified,
                "registry_truth_used": "0",
            })
    output.sort(key=lambda row: (row["ihw_number"], row.get("run_accession", "")))
    fields = ALLOWED_REGISTRY_FIELDS + ENA_FIELDS + [
        "modality", "exact_sample_title", "technical_match", "provenance_status",
        "registry_truth_used",
    ]
    write_tsv(output_tsv, output, fields)
    plausible = [row for row in output if row["provenance_status"] == "requires_study_level_review"]
    subjects_by_modality = {
        modality: len({row["ihw_number"] for row in plausible if row["modality"] == modality})
        for modality in ("wes", "rnaseq")
    }
    summary = {
        "schema_version": "ihwg-ena-scout-1",
        "truth_blind": True,
        "registry_sha256": sha256(registry_tsv),
        "registry_subjects": len(registry),
        "queried_or_cached_subjects": len(registry) - len(failures),
        "query_failures": failures,
        "all_runs": len(output),
        "runs_requiring_study_level_review": len(plausible),
        "candidate_subjects_by_modality": subjects_by_modality,
        "minimums": {"wes": 89, "rnaseq": 130},
        "availability_not_established_until_provenance_review": True,
        "endpoint": ENA_ENDPOINT,
        "output_sha256": sha256(output_tsv),
    }
    write_json(summary_json, summary)
    return summary


REVIEW_FIELDS = [
    "cell_line_provenance_verified", "orthogonal_dna_truth_available",
    "development_or_relative_overlap", "excluded_study", "evidence_url",
    "reviewer", "review_notes",
]


def build_ihwg_provenance_review_packet(matches_tsv: str, output_tsv: str,
                                        summary_json: str) -> dict:
    """Create a blank, fail-closed review sheet from technical name matches."""
    rows = read_tsv(matches_tsv)
    candidates = []
    for row in rows:
        if row.get("provenance_status") != "requires_study_level_review":
            continue
        record = dict(row)
        record.update({field: "" for field in REVIEW_FIELDS})
        candidates.append(record)
    candidates.sort(key=lambda row: (
        row.get("modality", ""), row.get("ihw_number", ""),
        -int(row.get("base_count", "0") or 0), row.get("run_accession", ""),
    ))
    fields = (list(rows[0]) if rows else ALLOWED_REGISTRY_FIELDS + ENA_FIELDS) + REVIEW_FIELDS
    write_tsv(output_tsv, candidates, fields)
    summary = {
        "schema_version": "ihwg-provenance-review-packet-1",
        "truth_blind": True,
        "source_sha256": sha256(matches_tsv),
        "candidate_runs": len(candidates),
        "candidate_subjects_by_modality": {
            modality: len({row["ihw_number"] for row in candidates if row["modality"] == modality})
            for modality in ("wes", "rnaseq")
        },
        "blank_review_fields": REVIEW_FIELDS,
        "fail_closed": True,
        "output_sha256": sha256(output_tsv),
    }
    write_json(summary_json, summary)
    return summary


def _yes(value: str) -> bool:
    return value.strip().lower() in {"1", "true", "yes", "y"}


def finalize_ihwg_read_roster(review_tsv: str, output_tsv: str, summary_json: str,
                              min_wes: int = 89, min_rnaseq: int = 130) -> dict:
    """Lock one provenance-verified run per line using truth-independent size."""
    rows = read_tsv(review_tsv)
    if not rows:
        raise ValueError("IHWG provenance review is empty")
    required = set(REVIEW_FIELDS) | {
        "ihw_number", "run_accession", "study_accession", "modality", "base_count",
        "technical_match", "exact_sample_title",
    }
    missing = sorted(required - set(rows[0]))
    if missing:
        raise ValueError(f"IHWG provenance review missing columns: {missing}")
    accepted = []
    for row in rows:
        complete_review = bool(row.get("reviewer", "").strip() and row.get("evidence_url", "").strip())
        eligible = (
            row.get("technical_match") == "1"
            and row.get("exact_sample_title") == "1"
            and _yes(row.get("cell_line_provenance_verified", ""))
            and _yes(row.get("orthogonal_dna_truth_available", ""))
            and not _yes(row.get("development_or_relative_overlap", ""))
            and not _yes(row.get("excluded_study", ""))
            and complete_review
            and row.get("modality") in {"wes", "rnaseq"}
        )
        if eligible:
            accepted.append(row)
    by_subject_modality: dict[tuple[str, str], list[dict[str, str]]] = {}
    for row in accepted:
        by_subject_modality.setdefault((row["ihw_number"], row["modality"]), []).append(row)
    selected = []
    for key, candidates in by_subject_modality.items():
        candidates.sort(key=lambda row: (
            -int(row.get("base_count", "0") or 0), row.get("run_accession", "")
        ))
        record = dict(candidates[0])
        record["selection_rule"] = "largest total post-QC library size; run accession tie-break"
        selected.append(record)
    selected.sort(key=lambda row: (row["modality"], row["ihw_number"]))
    fields = list(rows[0]) + ["selection_rule"]
    write_tsv(output_tsv, selected, fields)
    counts = {
        modality: len({row["ihw_number"] for row in selected if row["modality"] == modality})
        for modality in ("wes", "rnaseq")
    }
    summary = {
        "schema_version": "ihwg-public-read-roster-1",
        "truth_blind_selection": True,
        "review_sha256": sha256(review_tsv),
        "selected_subjects_by_modality": counts,
        "minimums": {"wes": min_wes, "rnaseq": min_rnaseq},
        "minimum_met": {"wes": counts["wes"] >= min_wes,
                        "rnaseq": counts["rnaseq"] >= min_rnaseq},
        "selection_rule": "largest total post-QC library size; run accession tie-break",
        "hla_depth_or_correctness_used": False,
        "output_sha256": sha256(output_tsv),
    }
    write_json(summary_json, summary)
    return summary
