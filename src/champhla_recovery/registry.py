from __future__ import annotations

from pathlib import Path

from .io import read_tsv, sha256


REQUIRED = {
    "result_id", "cohort", "modality", "method", "analysis_status", "validity",
    "subjects", "loci", "correct", "accuracy", "artifact_path", "artifact_sha256",
    "artifact_record_key", "source_artifact_path", "source_artifact_sha256",
    "evidence_role", "abstract_allowed", "claim_boundary",
}
EVIDENCE_ROLES = {
    "development", "same_resource_confirmation", "independent_validation",
    "exploratory", "invalid",
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
        elif artifact.suffix.lower() == ".tsv":
            selector = row["artifact_record_key"].split("=", 1)
            if len(selector) != 2 or not all(selector):
                failures.append(f"invalid artifact_record_key for {rid}")
            else:
                key, value = selector
                matches = [record for record in read_tsv(artifact) if record.get(key) == value]
                if len(matches) != 1:
                    failures.append(f"artifact row resolution for {rid} returned {len(matches)} rows")
                else:
                    record = matches[0]
                    for field in ("cohort", "modality", "method", "evidence_role",
                                  "analysis_status", "validity", "subjects", "loci",
                                  "correct", "accuracy"):
                        if record.get(field, "") != row.get(field, ""):
                            failures.append(f"artifact row mismatch for {rid}.{field}")
        source = base / row["source_artifact_path"]
        if not source.is_file():
            failures.append(f"missing source artifact for {rid}: {source}")
        elif sha256(source) != row["source_artifact_sha256"]:
            failures.append(f"source artifact checksum mismatch for {rid}")
        if row["evidence_role"] not in EVIDENCE_ROLES:
            failures.append(f"invalid evidence role for {rid}: {row['evidence_role']}")
        if row["validity"] not in {"valid", "invalid"}:
            failures.append(f"invalid validity for {rid}: {row['validity']}")
        for field in ("cohort", "modality", "method", "analysis_status", "claim_boundary"):
            if not row[field].strip():
                failures.append(f"empty {field} for {rid}")
        if row["loci"]:
            try:
                subjects = int(row["subjects"])
                loci = int(row["loci"])
                correct = int(row["correct"])
                accuracy = float(row["accuracy"])
                if subjects < 0 or loci <= 0 or not 0 <= correct <= loci:
                    failures.append(f"invalid counts for {rid}")
                elif abs(accuracy - correct / loci) > 1e-6:
                    failures.append(f"accuracy/count mismatch for {rid}")
            except ValueError:
                failures.append(f"nonnumeric counts for {rid}")
        elif row["correct"] or row["accuracy"]:
            failures.append(f"partial result counts for {rid}")
        if row["validity"] == "invalid" and row["evidence_role"] != "invalid":
            failures.append(f"invalid result lacks invalid evidence role: {rid}")
        if row["evidence_role"] == "invalid" and row["validity"] != "invalid":
            failures.append(f"invalid evidence role has valid result: {rid}")
        if row["abstract_allowed"] == "1" and row["validity"] != "valid":
            failures.append(f"invalid result allowed in abstract: {rid}")
        if row["evidence_role"] in {"development", "exploratory", "invalid"} and row["abstract_allowed"] == "1":
            failures.append(f"discovery result allowed in abstract: {rid}")
    return failures


def render_tables(registry_path: str, output_path: str) -> None:
    rows = read_tsv(registry_path)
    lines = [
        "# Registered ChampHLA results", "",
        "| Result | Cohort | Modality | Method | Evidence role | Validity | Correct / loci | Accuracy | Abstract |",
        "|---|---|---|---|---|---|---:|---:|---|",
    ]
    for row in rows:
        count = f"{row['correct']} / {row['loci']}" if row["loci"] else "—"
        accuracy = row["accuracy"] or "—"
        lines.append(
            f"| {row['result_id']} | {row['cohort']} | {row['modality']} | {row['method']} | "
            f"{row['evidence_role']} | {row['validity']} | {count} | {accuracy} | "
            f"{'yes' if row['abstract_allowed'] == '1' else 'no'} |"
        )
    lines.extend(["", "## Generated result sentences", ""])
    for row in rows:
        if not row["loci"] or row["validity"] != "valid":
            continue
        lines.append(
            f"- {row['method']} produced {row['correct']} correct genotypes among "
            f"{row['loci']} eligible {row['modality']} loci in {row['cohort']} "
            f"({row['evidence_role']}) [RESULT:{row['result_id']}]."
        )
    lines.extend(["", "## Invalid diagnostic records (supplement only)", ""])
    for row in rows:
        if not row["loci"] or row["validity"] != "invalid":
            continue
        lines.append(
            f"- INVALID DIAGNOSTIC: {row['method']} recorded {row['correct']} of "
            f"{row['loci']} loci; this cannot be used as performance evidence "
            f"[RESULT:{row['result_id']}]."
        )
    target = Path(output_path)
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text("\n".join(lines) + "\n", encoding="utf-8")
