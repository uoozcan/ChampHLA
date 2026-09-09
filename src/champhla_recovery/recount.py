from __future__ import annotations

from collections import defaultdict
from pathlib import Path

from champhla_confirmation.io import canonical_pair, read_tsv, sha256, write_json, write_tsv
from champhla_confirmation.panels import canonical_method


def independent_recount(joined_path: str | Path, output_tsv: str | Path,
                        output_json: str | Path, method: str = "SimplePluralityLex") -> dict:
    """Recount primary exact matches without importing the production evaluator."""
    wanted = canonical_method(method)
    buckets = defaultdict(lambda: {"subjects": set(), "loci": 0, "correct": 0})
    seen = set()
    for row in read_tsv(joined_path):
        if canonical_method(row.get("method", "")) != wanted:
            continue
        key = (row.get("cohort", ""), row["subject"], row["modality"], row["gene"])
        if key in seen:
            raise ValueError(f"duplicate primary prediction row: {key}")
        seen.add(key)
        stratum = row.get("independence_stratum", "") or "unspecified"
        bucket = buckets[(row["modality"], stratum)]
        bucket["subjects"].add(row.get("donor", "") or row["subject"])
        bucket["loci"] += 1
        try:
            predicted = canonical_pair(row.get("allele1", ""), row.get("allele2", ""), row["gene"])
            truth = canonical_pair(row["truth_allele1"], row["truth_allele2"], row["gene"])
            bucket["correct"] += int(predicted == truth and row.get("truth_status", "resolved") == "resolved")
        except ValueError:
            pass
    rows = []
    for (modality, stratum), values in sorted(buckets.items()):
        rows.append({
            "modality": modality, "independence_stratum": stratum,
            "method": wanted, "subjects": len(values["subjects"]),
            "loci": values["loci"], "correct": values["correct"],
            "accuracy": values["correct"] / values["loci"] if values["loci"] else 0.0,
        })
    write_tsv(output_tsv, rows)
    result = {
        "schema_version": "champhla-independent-recount-1",
        "method": wanted, "rows": rows, "input_sha256": sha256(joined_path),
        "output_sha256": sha256(output_tsv), "imports_primary_evaluator": False,
    }
    write_json(output_json, result)
    return result
