from __future__ import annotations

import random
from collections import defaultdict

from .io import read_json, read_tsv, write_json
from .panels import MODALITIES
from .statistics import exact_cluster_signflip, holm_adjust


def simulate(paired_path: str, output_json: str, targets: dict[str, int],
             iterations: int = 10000, seed: int = 20260831,
             wgs_audit_summary: str | None = None) -> dict:
    rows = read_tsv(paired_path)
    by_modality_subject = defaultdict(lambda: defaultdict(int))
    for row in rows:
        by_modality_subject[row["modality"]][row["subject"]] += (
            int(row["guarded_correct"]) - int(row["baseline_correct"]))
    missing = [m for m in MODALITIES if not by_modality_subject[m]]
    if missing:
        raise ValueError(f"power simulation requires discovery deltas for all modalities: {missing}")
    rng = random.Random(seed)
    modality_success = defaultdict(int)
    all_success = 0
    for _ in range(iterations):
        tests = []
        for modality in MODALITIES:
            observed = list(by_modality_subject[modality].values())
            sampled = [rng.choice(observed) for _ in range(targets[modality])]
            tests.append({"modality": modality, "delta": sum(sampled),
                          "p_value": exact_cluster_signflip(sampled)})
        tests = holm_adjust(tests)
        passed = []
        for test in tests:
            ok = test["delta"] > 0 and bool(test["holm_significant"])
            modality_success[test["modality"]] += int(ok)
            passed.append(ok)
        all_success += int(all(passed))
    wgs_valid = bool(wgs_audit_summary and read_json(wgs_audit_summary).get("passed"))
    payload = {
        "schema_version": "confirmation-power-simulation-1", "iterations": iterations,
        "seed": seed, "targets": targets,
        "modality_power": {m: modality_success[m] / iterations for m in MODALITIES},
        "all_three_power": all_success / iterations,
        "wgs_discovery_audit_pass": wgs_valid,
        "passed_80_percent": bool(wgs_valid and all_success / iterations >= 0.80),
        "status": "power_gate_passed" if wgs_valid and all_success / iterations >= 0.80
                  else "wgs_repair_required_before_valid_power_gate",
        "warning": "discovery-resampling projection only; external significance remains the decision gate",
    }
    write_json(output_json, payload)
    return payload
