from __future__ import annotations

from .io import canonical_allele, canonical_pair, normalize_gene, normalize_modality
from .panels import GENES, PANELS, canonical_method

MISSING_ALLELE_TOKENS = {"", "-", "NOTTYPED", "NOT_TYPED", "NA", "N/A", "NONE", "UNKNOWN", "."}


def validate_production_caller_matrix(rows: list[dict[str, str]]) -> None:
    """Require one explicit frozen-panel record for every subject/modality/locus."""
    by_locus: dict[tuple[str, str, str, str], set[str]] = {}
    subject_modalities = set()
    for row in rows:
        key = (row.get("cohort", ""), row["subject"], row["modality"], row["gene"])
        caller = row["caller"]
        if caller not in PANELS[row["modality"]]:
            raise ValueError(f"assay-incompatible caller {caller} at {key}")
        callers = by_locus.setdefault(key, set())
        if caller in callers:
            raise ValueError(f"duplicate production caller row for {caller} at {key}")
        callers.add(caller)
        subject_modalities.add(key[:3])
    failures = []
    for cohort, subject, modality in sorted(subject_modalities):
        expected = set(PANELS[modality])
        for gene in GENES:
            key = (cohort, subject, modality, gene)
            observed = by_locus.get(key, set())
            if observed != expected:
                failures.append(
                    f"{key}: missing={sorted(expected - observed)} extra={sorted(observed - expected)}"
                )
    if failures:
        raise ValueError(f"incomplete production caller matrix: {failures[:10]}")


def normalize_status(value: str, has_pair: bool) -> str:
    token = value.strip().lower()
    if token in {"callable", "called", "call", "pass", "passed", "typed"}:
        return "callable" if has_pair else "missing"
    if token in {"", "unknown"}:
        return "callable" if has_pair else "missing"
    return token


def normalize_caller_row(row: dict[str, str], cohort_default: str = "") -> dict[str, str]:
    subject = row.get("subject", row.get("sample", "")).strip()
    caller = row.get("caller", row.get("tool", "")).strip()
    if not subject or not caller:
        raise ValueError("caller rows require subject/sample and caller/tool")
    modality = normalize_modality(row["modality"])
    gene = normalize_gene(row["gene"])
    raw_a1, raw_a2 = row.get("allele1", "").strip(), row.get("allele2", "").strip()
    present = [value for value in (raw_a1, raw_a2) if value.upper() not in MISSING_ALLELE_TOKENS]
    explicit_homozygous = row.get("call_status", "").strip().lower() in {
        "homozygous", "homozygote", "homo", "hom",
    }
    status = normalize_status(row.get("call_status", ""), len(present) == 2)
    if len(present) == 2:
        if status == "partial":
            raise ValueError("partial caller row contains two alleles")
        a1, a2 = canonical_pair(present[0], present[1], gene)
        status = "callable"
    elif len(present) == 1:
        if explicit_homozygous:
            a1, a2 = canonical_pair(present[0], present[0], gene)
            status = "callable"
        else:
            a1, a2 = canonical_allele(present[0], gene), ""
            status = "partial"
    else:
        a1 = a2 = ""
        status = "missing"
    return {
        "cohort": row.get("cohort", cohort_default), "subject": subject,
        "superpopulation": row.get("superpopulation", row.get("population", "")),
        "modality": modality, "gene": gene, "caller": caller,
        "allele1": a1, "allele2": a2, "call_status": status,
        "allele1_raw": row.get("allele1_raw", raw_a1),
        "allele2_raw": row.get("allele2_raw", raw_a2),
        "caller_version": row.get("caller_version", ""),
        "source_path": row.get("source_path", row.get("source_file", "")),
        "source_sha256": row.get("source_sha256", ""),
    }


def normalize_raw_cc_row(row: dict[str, str], cohort_default: str = "") -> dict[str, str]:
    subject = row.get("subject", row.get("sample", "")).strip()
    modality = normalize_modality(row["modality"])
    gene = normalize_gene(row["gene"])
    a1, a2 = row.get("allele1", "").strip(), row.get("allele2", "").strip()
    status = normalize_status(row.get("call_status", ""), bool(a1 and a2))
    if a1.upper() in MISSING_ALLELE_TOKENS or a2.upper() in MISSING_ALLELE_TOKENS:
        status = "missing"
    if status == "callable" and a1 and a2:
        a1, a2 = canonical_pair(a1, a2, gene)
    else:
        status, a1, a2 = "missing", "", ""
    return {
        "cohort": row.get("cohort", cohort_default), "subject": subject,
        "superpopulation": row.get("superpopulation", row.get("population", "")),
        "modality": modality, "gene": gene, "method": "ChampionChallenger",
        "allele1": a1, "allele2": a2, "call_status": status,
    }


def normalize_method_row(row: dict[str, str], cohort_default: str = "") -> dict[str, str]:
    normalized = normalize_raw_cc_row({**row, "method": "ChampionChallenger"}, cohort_default)
    normalized["method"] = canonical_method(row.get("method", "").strip())
    if not normalized["method"]:
        raise ValueError("secondary method row has no method")
    return normalized


def explode_candidate_rows(rows: list[dict[str, str]], cohort_default: str = "") -> list[dict[str, str]]:
    """Convert candidate/support rows into one exact call per supporting caller."""
    output = []
    seen = set()
    for row in rows:
        callers = row.get("supporting_callers", row.get("supporting_tools", "")).split(",")
        for caller in filter(None, (value.strip() for value in callers)):
            converted = normalize_caller_row({**row, "caller": caller}, cohort_default)
            key = (converted["cohort"], converted["subject"], converted["modality"],
                   converted["gene"], converted["caller"])
            if key in seen:
                raise ValueError(f"caller appears in multiple candidate pairs at {key}")
            seen.add(key)
            output.append(converted)
    return output
