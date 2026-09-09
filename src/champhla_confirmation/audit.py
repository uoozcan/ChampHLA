from __future__ import annotations

from collections import Counter, defaultdict
from pathlib import Path

from .io import canonical_allele, canonical_pair, normalize_gene, normalize_modality, read_json, read_tsv, sha256, write_json, write_tsv
from .panels import GENES, PANELS
from .parsers import parse_caller_call


def _safe_pair(a: str, b: str, gene: str):
    try:
        return canonical_pair(a, b, gene), ""
    except ValueError as exc:
        return None, str(exc)


def _raw_tokens_found(path: Path, alleles: tuple[str, ...] | None) -> bool:
    if not alleles:
        return True
    try:
        text = path.read_text(encoding="utf-8", errors="replace").upper().replace("HLA-", "")
    except OSError:
        return False
    tokens = []
    for allele in alleles:
        tokens.extend((allele.upper(), allele.split("*", 1)[-1].upper()))
    return all(any(token in text for token in pair) for pair in zip(tokens[0::2], tokens[1::2]))


def audit_wgs(harmonized_path: str | Path, output_dir: str | Path,
              manual_review_path: str | Path | None = None,
              environment_manifest_path: str | Path | None = None) -> dict:
    source_rows = read_tsv(harmonized_path)
    rows = [row for row in source_rows if normalize_modality(row.get("modality", "")) == "wgs"
            and row.get("tool", row.get("caller", "")) in PANELS["wgs"]]
    if not rows:
        raise ValueError("no eligible WGS rows")
    cohort_default = "1000G-development"
    observed = {}
    observed_subjects = set()
    populations = {}
    for row in rows:
        subject = row.get("subject", row.get("sample", "")).strip()
        gene = normalize_gene(row["gene"])
        caller = row.get("caller", row.get("tool", "")).strip()
        key = (row.get("cohort", cohort_default), subject, gene, caller)
        if key in observed:
            raise ValueError(f"duplicate WGS caller/locus row: {key}")
        observed[key] = row
        observed_subjects.add((key[0], subject))
        populations[subject] = row.get("superpopulation", row.get("population", ""))

    universe = {(cohort, subject, gene) for cohort, subject in observed_subjects for gene in GENES}

    file_hash_cache = {}
    audit_rows = []
    for cohort, subject, gene in sorted(universe):
        for caller in PANELS["wgs"]:
            key = (cohort, subject, gene, caller)
            row = observed.get(key)
            if row is None:
                audit_rows.append({
                    "cohort": cohort, "subject": subject, "superpopulation": populations.get(subject, ""),
                    "gene": gene, "caller": caller, "call_status": "missing_unrecorded",
                    "raw_allele1": "", "raw_allele2": "", "parsed_allele1": "", "parsed_allele2": "",
                    "raw_parsed_match": 0, "source_path": "", "source_exists": 0,
                    "source_sha256": "", "source_hash_match": 0, "raw_tokens_found": 0,
                    "audit_status": "fail_absent_expected_record", "audit_detail": "",
                })
                continue
            status = row.get("call_status", "").strip().lower() or (
                "callable" if row.get("allele1") and row.get("allele2") else "missing")
            raw_a1 = row.get("allele1_raw", row.get("raw_allele1", row.get("allele1", "")))
            raw_a2 = row.get("allele2_raw", row.get("raw_allele2", row.get("allele2", "")))
            parsed_a1, parsed_a2 = row.get("allele1", ""), row.get("allele2", "")
            raw_pair, raw_error = _safe_pair(raw_a1, raw_a2, gene) if status == "callable" else (None, "")
            parsed_pair, parsed_error = _safe_pair(parsed_a1, parsed_a2, gene) if status == "callable" else (None, "")
            partial_allele = None
            raw_partial_allele = None
            if status == "partial":
                try:
                    partial_allele = canonical_allele(parsed_a1, gene)
                    if parsed_a2:
                        raise ValueError("partial row has a second parsed allele")
                except ValueError as exc:
                    parsed_error = str(exc)
                try:
                    raw_partial_allele = canonical_allele(raw_a1, gene)
                    if raw_a2.strip() not in {"", "-"}:
                        raise ValueError("partial raw row has a second allele")
                except ValueError as exc:
                    raw_error = str(exc)
            pair_match = status not in {"callable", "partial"} or (
                bool(raw_pair and raw_pair == parsed_pair) if status == "callable"
                else bool(partial_allele and raw_partial_allele == partial_allele)
            )
            source_path = row.get("source_path", row.get("source_file", ""))
            path = Path(source_path) if source_path else None
            exists = bool(path and path.is_file())
            actual_hash = ""
            if exists:
                if source_path not in file_hash_cache:
                    file_hash_cache[source_path] = sha256(path)
                actual_hash = file_hash_cache[source_path]
            expected_hash = row.get("source_sha256", "")
            hash_match = bool(exists and (not expected_hash or expected_hash == actual_hash))
            audit_alleles = raw_pair if status == "callable" else ((partial_allele,) if partial_allele else None)
            tokens_found = bool(exists and _raw_tokens_found(path, audit_alleles))
            source_call = None
            if exists and status in {"callable", "partial"}:
                try:
                    source_call = parse_caller_call(
                        caller, path.read_text(encoding="utf-8", errors="replace"), gene)
                except (OSError, ValueError):
                    source_call = None
            source_pair_match = status not in {"callable", "partial"} or bool(
                source_call
                and source_call["call_status"] == status
                and source_call["allele1"] == parsed_a1
                and source_call["allele2"] == parsed_a2
            )
            failures = []
            if status in {"callable", "partial"} and not pair_match:
                failures.append("raw_parsed_mismatch")
            if status in {"callable", "partial"} and (raw_error or parsed_error):
                failures.append("unresolved_allele")
            if not exists:
                failures.append("source_unavailable")
            elif not hash_match:
                failures.append("source_hash_mismatch")
            if status in {"callable", "partial"} and exists and not tokens_found:
                failures.append("raw_tokens_not_found")
            if status in {"callable", "partial"} and exists and not source_pair_match:
                failures.append("source_pair_unparsed_or_mismatch")
            audit_rows.append({
                "cohort": cohort, "subject": subject, "superpopulation": populations.get(subject, ""),
                "gene": gene, "caller": caller, "call_status": status,
                "raw_allele1": raw_a1, "raw_allele2": raw_a2,
                "parsed_allele1": parsed_a1, "parsed_allele2": parsed_a2,
                "raw_parsed_match": int(pair_match), "source_path": source_path,
                "source_exists": int(exists), "source_sha256": actual_hash,
                "source_hash_match": int(hash_match), "raw_tokens_found": int(tokens_found),
                "source_pair_match": int(source_pair_match),
                "audit_status": "pass" if not failures else "fail_" + "+".join(failures),
                "audit_detail": raw_error or parsed_error,
            })

    existing_review = {}
    if manual_review_path:
        for row in read_tsv(manual_review_path):
            existing_review[(row["subject"], row["gene"], row["caller"])] = row
    by_caller = defaultdict(list)
    for row in audit_rows:
        if row["call_status"] == "missing":
            row["review_stratum"] = "missing"
        elif row["call_status"] in {"callable", "partial"}:
            homo = row["parsed_allele1"] == row["parsed_allele2"]
            high = row["raw_allele1"].count(":") > 1 or row["raw_allele2"].count(":") > 1
            row["review_stratum"] = ("partial" if row["call_status"] == "partial" else
                                     "high_field" if high else "homozygous" if homo else "heterozygous")
        else:
            continue
        by_caller[row["caller"]].append(row)
    manual_rows = []
    for caller in PANELS["wgs"]:
        candidates = by_caller[caller]
        selected = []
        for stratum in ("partial", "high_field", "homozygous", "heterozygous", "missing"):
            for row in candidates:
                if row["review_stratum"] == stratum and row not in selected:
                    selected.append(row)
                    if len([r for r in selected if r["review_stratum"] == stratum]) >= 3:
                        break
        for row in candidates:
            if len(selected) >= 10:
                break
            if row not in selected:
                selected.append(row)
        for row in selected[:10]:
            previous = existing_review.get((row["subject"], row["gene"], caller), {})
            manual_rows.append({
                "subject": row["subject"], "gene": row["gene"], "caller": caller,
                "review_stratum": row["review_stratum"], "source_path": row["source_path"],
                "raw_allele1": row["raw_allele1"], "raw_allele2": row["raw_allele2"],
                "parsed_allele1": row["parsed_allele1"], "parsed_allele2": row["parsed_allele2"],
                "review_status": previous.get("review_status", "pending"),
                "reviewer": previous.get("reviewer", ""), "review_note": previous.get("review_note", ""),
            })
    counts = Counter(row["audit_status"] for row in audit_rows)
    review_counts = Counter(row["review_status"] for row in manual_rows)
    automated_passed = bool(audit_rows and all(row["audit_status"] == "pass" for row in audit_rows))
    manual_review_complete = bool(
        len(manual_rows) == 10 * len(PANELS["wgs"])
        and all(row["review_status"] == "pass" for row in manual_rows)
        and all(row["reviewer"].strip() for row in manual_rows)
    )
    environment_failures = []
    environment = read_json(environment_manifest_path) if environment_manifest_path else {}
    required_environment = {
        "reference_build", "reference_sha256", "caller_reference_attestation",
        "python_version", "samtools_version", "caller_artifacts",
    }
    if required_environment - set(environment):
        environment_failures.append(
            f"missing environment fields: {sorted(required_environment - set(environment))}"
        )
    elif set(environment["caller_artifacts"]) != set(PANELS["wgs"]):
        environment_failures.append("caller_artifacts does not exactly match the WGS panel")
    elif not all(str(value).strip() for value in environment["caller_artifacts"].values()):
        environment_failures.append("one or more caller artifact hashes are empty")
    elif environment["caller_reference_attestation"].get("passed") is not True:
        environment_failures.append("caller/reference attestation is not production ready")
    environment_passed = not environment_failures
    passed = automated_passed and manual_review_complete and environment_passed
    out = Path(output_dir)
    write_tsv(out / "wgs_audit.tsv", audit_rows)
    write_tsv(out / "wgs_manual_review.tsv", manual_rows)
    summary = {
        "schema_version": "wgs-input-audit-1", "passed": passed,
        "automated_passed": automated_passed,
        "manual_review_complete": manual_review_complete,
        "human_review_passed": manual_review_complete,
        "named_reviewer": ",".join(sorted({row["reviewer"] for row in manual_rows if row["reviewer"]})),
        "environment_passed": environment_passed,
        "environment_failures": environment_failures,
        "environment_manifest": str(Path(environment_manifest_path).resolve()) if environment_manifest_path else "",
        "environment_manifest_sha256": sha256(environment_manifest_path) if environment_manifest_path else "",
        "eligible_panel": list(PANELS["wgs"]), "locus_caller_records": len(audit_rows),
        "expected_locus_caller_records": len(universe) * len(PANELS["wgs"]),
        "observed_input_locus_caller_records": len(observed),
        "expected_record_fraction": (len(observed) / (len(universe) * len(PANELS["wgs"]))
                                     if universe else 0.0),
        "native_output_checksums_passed": bool(audit_rows and all(
            row["source_hash_match"] for row in audit_rows)),
        "parser_round_trips_passed": bool(audit_rows and all(
            row.get("source_pair_match", 0) for row in audit_rows)),
        "audit_status_counts": dict(sorted(counts.items())),
        "manual_review_counts": dict(sorted(review_counts.items())),
        "source_harmonized": str(Path(harmonized_path).resolve()),
        "source_harmonized_sha256": sha256(harmonized_path),
        "failure_policy": "fail closed until every expected record, caller-native source, raw/parsed pair, checksum, environment/reference manifest, and 50-record named manual review passes",
    }
    write_json(out / "wgs_audit_summary.json", summary)
    return summary
