from __future__ import annotations

from pathlib import Path

from .io import read_tsv, reject_truth_columns, sha256, write_json, write_tsv
from .panels import GENES, PANELS
from .parsers import parse_caller_pair


CALLER_SUFFIX = {
    "HLA-HD": "hlahd",
    "Kourami": "kourami",
    "OptiType": "optitype",
    "SpecHLA": "spechla",
    "T1K": "t1k",
}
PINNED_VERSION = {
    "HLA-HD": "1.4.0",
    "Kourami": "0.9.6",
    "OptiType": "pinned-container-see-environment-manifest",
    "SpecHLA": "pinned-local-wrapper-see-environment-manifest",
    "T1K": "1.0.9",
}


def collect_full_cram_wgs_pilot(caller_root: str, manifest_tsv: str,
                                samples: list[str], output_tsv: str,
                                summary_json: str) -> dict:
    """Collect five caller-native summaries only after the strict run marker."""
    manifest = read_tsv(manifest_tsv)
    reject_truth_columns(manifest, "WGS pilot manifest")
    manifest_samples = {row["sample_id"]: row for row in manifest}
    selected = samples or sorted(manifest_samples)
    unknown = sorted(set(selected) - set(manifest_samples))
    if unknown:
        raise ValueError(f"pilot samples absent from manifest: {unknown}")
    root = Path(caller_root)
    rows = []
    file_manifest = []
    for sample in selected:
        sample_root = root / sample
        marker = sample_root / "CALLERS_COMPLETE"
        validation = sample_root / "caller_output_validation.tsv"
        if not marker.is_file() or not validation.is_file():
            raise ValueError(f"strict five-caller completion marker absent for {sample}")
        validated = read_tsv(validation)
        if len(validated) != len(PANELS["wgs"]):
            raise ValueError(f"invalid caller-output validation record count for {sample}")
        for caller in PANELS["wgs"]:
            suffix = CALLER_SUFFIX[caller]
            matches = list((sample_root / "results").rglob(f"{sample}_{suffix}.txt"))
            if len(matches) != 1 or not matches[0].is_file() or matches[0].stat().st_size == 0:
                raise ValueError(f"expected one nonempty {caller} summary for {sample}, found {len(matches)}")
            source = matches[0]
            source_hash = sha256(source)
            validation_rows = [row for row in validated if row.get("caller") == suffix]
            if len(validation_rows) != 1 or validation_rows[0].get("sha256") != source_hash:
                raise ValueError(f"caller-output checksum mismatch for {sample} {caller}")
            text = source.read_text(encoding="utf-8", errors="replace")
            for gene in GENES:
                pair = parse_caller_pair(caller, text, gene)
                rows.append({
                    "cohort": "1000G-WGS-full-CRAM-pilot",
                    "subject": sample,
                    "modality": "wgs",
                    "gene": gene,
                    "caller": caller,
                    "allele1_raw": pair[0] if pair else "",
                    "allele2_raw": pair[1] if pair else "",
                    "allele1": pair[0] if pair else "",
                    "allele2": pair[1] if pair else "",
                    "call_status": "callable" if pair else "missing",
                    "caller_version": PINNED_VERSION[caller],
                    "source_path": str(source.resolve()),
                    "source_sha256": source_hash,
                })
            file_manifest.append({"sample": sample, "caller": caller,
                                  "path": str(source.resolve()), "sha256": source_hash})
    rows.sort(key=lambda row: (row["subject"], row["gene"], row["caller"]))
    write_tsv(output_tsv, rows)
    expected = len(selected) * len(GENES) * len(PANELS["wgs"])
    summary = {
        "schema_version": "full-cram-wgs-pilot-collection-1",
        "truth_blind": True,
        "samples": selected,
        "callers": list(PANELS["wgs"]),
        "expected_locus_caller_records": expected,
        "observed_locus_caller_records": len(rows),
        "explicit_status_for_every_record": len(rows) == expected,
        "callable_records": sum(row["call_status"] == "callable" for row in rows),
        "manifest_sha256": sha256(manifest_tsv),
        "caller_files": file_manifest,
        "output_sha256": sha256(output_tsv),
        "truth_or_correctness_used": False,
    }
    write_json(summary_json, summary)
    return summary
