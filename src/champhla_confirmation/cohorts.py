from __future__ import annotations

from collections import defaultdict

from .io import read_tsv, write_json, write_tsv


def _tokens(row: dict[str, str]) -> set[str]:
    fields = ("subject", "sample", "sample_id", "alias", "aliases", "coriell_id", "ihwg_id",
              "biosample", "biosample_accession", "individual_id")
    result = set()
    for field in fields:
        for token in row.get(field, "").replace(";", ",").split(","):
            token = token.strip().upper()
            if token:
                result.add(token)
    return result


def build_overlap_crosswalk(development_path: str, candidate_path: str,
                            output_tsv: str, summary_json: str) -> dict:
    development = read_tsv(development_path)
    candidates = read_tsv(candidate_path)
    dev_index = defaultdict(set)
    for row in development:
        canonical = row.get("subject", row.get("sample", row.get("sample_id", "")))
        for token in _tokens(row):
            dev_index[token].add(canonical)
    output = []
    for row in candidates:
        subject = row.get("subject", row.get("sample", row.get("sample_id", "")))
        matches = sorted({match for token in _tokens(row) for match in dev_index.get(token, set())})
        related = row.get("development_relative", "").strip().lower() in {"1", "true", "yes"}
        output.append({
            "cohort": row.get("cohort", ""), "subject": subject,
            "candidate_aliases": ",".join(sorted(_tokens(row))),
            "development_matches": ",".join(matches),
            "known_development_relative": int(related),
            "eligible_nonoverlap": int(not matches and not related),
            "exclusion_reason": "identifier_overlap" if matches else "known_relative" if related else "",
        })
    summary = {
        "schema_version": "subject-overlap-crosswalk-1", "candidate_subjects": len(output),
        "eligible_nonoverlap": sum(row["eligible_nonoverlap"] for row in output),
        "excluded_overlap_or_relative": sum(not row["eligible_nonoverlap"] for row in output),
    }
    write_tsv(output_tsv, output)
    write_json(summary_json, summary)
    return summary

