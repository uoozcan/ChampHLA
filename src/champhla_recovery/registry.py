from __future__ import annotations

from pathlib import Path

from .io import read_tsv, sha256


REQUIRED = {
    "result_id", "cohort", "modality", "method", "analysis_status", "validity",
    "subjects", "loci", "correct", "accuracy", "artifact_path", "artifact_sha256",
    "abstract_allowed", "claim_boundary",
}


def validate_registry(path: str, root: str | None = None) -> list[str]:
    rows = read_tsv(path)
    failures = []
    if not rows:
        return ["result registry is empty"]
    missing = sorted(REQUIRED - set(rows[0]))
    if missing:
        return [f"missing registry columns: {missing}"]
    seen = set()
    base = Path(root) if root else Path(path).resolve().parent
    for row in rows:
        rid = row["result_id"]
        if rid in seen:
            failures.append(f"duplicate result_id: {rid}")
        seen.add(rid)
        artifact = base / row["artifact_path"]
        if not artifact.is_file():
            failures.append(f"missing artifact for {rid}: {artifact}")
        elif sha256(artifact) != row["artifact_sha256"]:
            failures.append(f"artifact checksum mismatch for {rid}")
        if row["abstract_allowed"] == "1" and row["validity"] != "valid":
            failures.append(f"invalid result allowed in abstract: {rid}")
        if row["analysis_status"] == "discovery" and row["abstract_allowed"] == "1":
            failures.append(f"discovery result allowed in abstract: {rid}")
    return failures


def render_tables(registry_path: str, output_path: str) -> None:
    rows = read_tsv(registry_path)
    lines = [
        "# Registered ChampHLA results", "",
        "| Result | Cohort | Modality | Method | Status | Validity | Correct / loci | Accuracy | Abstract |",
        "|---|---|---|---|---|---|---:|---:|---|",
    ]
    for row in rows:
        count = f"{row['correct']} / {row['loci']}" if row["loci"] else "—"
        accuracy = row["accuracy"] or "—"
        lines.append(
            f"| {row['result_id']} | {row['cohort']} | {row['modality']} | {row['method']} | "
            f"{row['analysis_status']} | {row['validity']} | {count} | {accuracy} | "
            f"{'yes' if row['abstract_allowed'] == '1' else 'no'} |"
        )
    target = Path(output_path)
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text("\n".join(lines) + "\n", encoding="utf-8")

