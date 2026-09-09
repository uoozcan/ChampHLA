#!/usr/bin/env python3
from __future__ import annotations

import csv
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "artifacts/source_results/primary_modality_results.tsv"
OUTPUT = ROOT / "artifacts/verification/discovery_reproduction.json"
EXPECTED = {
    "wes": {"subjects": 130, "loci": 390, "baseline_correct": 328,
            "guarded_correct": 369, "plurality_correct": 365},
    "rnaseq": {"subjects": 107, "loci": 321, "baseline_correct": 296,
               "guarded_correct": 306, "plurality_correct": 306},
}


def main() -> int:
    with SOURCE.open(newline="", encoding="utf-8") as handle:
        rows = {row["modality"]: row for row in csv.DictReader(handle, delimiter="\t")}
    failures = []
    observed = {}
    for modality, expected in EXPECTED.items():
        observed[modality] = {}
        for field, value in expected.items():
            actual = int(float(rows[modality][field]))
            observed[modality][field] = actual
            if actual != value:
                failures.append(f"{modality}.{field}: expected {value}, observed {actual}")
    result = {
        "schema_version": "development-table-reconciliation-2",
        "passed": not failures,
        "scope": "table-level count reconciliation against locked expected values",
        "independent_locus_recalculation": False,
        "release_gate_satisfied": False,
        "observed": observed,
        "failures": failures,
        "wgs_status": "excluded_invalid_historical_input",
    }
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return 0 if result["passed"] else 2


if __name__ == "__main__":
    raise SystemExit(main())
