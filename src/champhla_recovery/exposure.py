from __future__ import annotations

import re
from pathlib import Path

from .io import sha256, write_json, write_tsv


SUBJECT_RE = re.compile(r"(?<![A-Za-z0-9])(?:HG|NA)\d{5}(?![A-Za-z0-9])", re.I)


def _role(path: Path) -> str:
    name = str(path).lower()
    if any(token in name for token in ("train", "fold", "development", "label")):
        return "development_exposed"
    if "holdout" in name:
        return "holdout_exposed"
    if any(token in name for token in ("predict", "evaluation", "joined", "pilot")):
        return "evaluation_exposed"
    return "manifest_exposed"


def build_exposure_ledger(sources: list[str], output: str, summary: str) -> dict:
    records: dict[tuple[str, str], dict] = {}
    source_records = []
    for source in sorted({str(Path(item).resolve()) for item in sources}):
        path = Path(source)
        if not path.is_file():
            raise FileNotFoundError(source)
        digest = sha256(path)
        text = path.read_text(encoding="utf-8", errors="replace")
        subjects = sorted({match.upper() for match in SUBJECT_RE.findall(text)})
        role = _role(path)
        for subject in subjects:
            records[(subject, source)] = {
                "subject": subject,
                "canonical_subject": subject,
                "exposure_role": role,
                "source_path": source,
                "source_sha256": digest,
            }
        source_records.append({"path": source, "sha256": digest, "subjects": len(subjects)})
    rows = sorted(records.values(), key=lambda row: (row["subject"], row["source_path"]))
    write_tsv(output, rows, [
        "subject", "canonical_subject", "exposure_role", "source_path", "source_sha256",
    ])
    result = {
        "schema_version": "subject-exposure-ledger-1",
        "sources": source_records,
        "ledger_rows": len(rows),
        "unique_subjects": len({row["subject"] for row in rows}),
        "output_sha256": sha256(output),
    }
    write_json(summary, result)
    return result

