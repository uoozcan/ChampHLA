from __future__ import annotations

import hashlib
import math
from collections import Counter, defaultdict

from .io import read_json, read_tsv, sha256, write_json, write_tsv


def _allocate(counts: Counter, target: int) -> dict[str, int]:
    total = sum(counts.values())
    if target > total:
        raise ValueError(f"target {target} exceeds eligible subjects {total}")
    exact = {key: target * value / total for key, value in counts.items()}
    allocation = {key: min(counts[key], math.floor(value)) for key, value in exact.items()}
    remaining = target - sum(allocation.values())
    order = sorted(counts, key=lambda key: (-(exact[key] - allocation[key]), key))
    while remaining:
        progressed = False
        for key in order:
            if allocation[key] < counts[key] and remaining:
                allocation[key] += 1
                remaining -= 1
                progressed = True
        if not progressed:
            raise RuntimeError("population allocation stalled")
    return allocation


def freeze_rosters(registry_path: str, assay_manifest_path: str, exposure_path: str,
                   config_path: str, output_path: str, summary_path: str) -> dict:
    registry = {row["subject"].upper(): row for row in read_tsv(registry_path)}
    exposed = {row["subject"].upper() for row in read_tsv(exposure_path)}
    config = read_json(config_path)
    seed = config["selection_seed"]
    candidates: dict[str, dict[str, dict]] = defaultdict(dict)
    for row in read_tsv(assay_manifest_path):
        subject = (row.get("sample_id") or row.get("subject") or "").upper()
        modality = row.get("modality", "").lower().replace("-", "")
        if modality == "rna":
            modality = "rnaseq"
        truth = registry.get(subject)
        primary_eligible = truth and truth.get("primary_eligible", truth.get("complete_truth_available")) == "1"
        if not primary_eligible or subject in exposed:
            continue
        if row.get("development_relative", "").lower() in {"1", "true", "yes"}:
            continue
        population = row.get("population") or truth.get("population") or "UNRECORDED"
        candidates[modality][subject] = {**row, "subject": subject, "population": population}

    output = []
    summaries = {}
    all_minimums_met = True
    for modality in ("wgs", "wes", "rnaseq"):
        target = int(config["initial_targets"][modality])
        by_population: dict[str, list[dict]] = defaultdict(list)
        for row in candidates.get(modality, {}).values():
            by_population[row["population"]].append(row)
        counts = Counter({key: len(value) for key, value in by_population.items()})
        minimum_met = sum(counts.values()) >= target
        all_minimums_met = all_minimums_met and minimum_met
        allocation = _allocate(counts, target) if minimum_met else dict(counts)
        selected = []
        for population in sorted(by_population):
            ranked = sorted(
                by_population[population],
                key=lambda row: hashlib.sha256(
                    f"{seed}|{modality}|{row['subject']}".encode()
                ).hexdigest(),
            )
            selected.extend(ranked[:allocation[population]])
        for row in sorted(selected, key=lambda item: item["subject"]):
            output.append({
                "cohort": "1000G-2014-unseen",
                "subject": row["subject"],
                "modality": modality,
                "population": row["population"],
                "input_url": row.get("input_url", ""),
                "index_url": row.get("index_url", ""),
                "input_md5": row.get("input_md5", ""),
                "genome_build": row.get("genome_build", ""),
                "selection_seed": seed,
                "primary_confirmation_eligible": int(minimum_met),
            })
        summaries[modality] = {
            "eligible": sum(counts.values()), "selected": len(selected),
            "target": target, "minimum_met": minimum_met,
            "status": "roster_frozen" if minimum_met else "infeasible_minimum_not_met",
            "population_allocation": dict(sorted(allocation.items())),
        }
    write_tsv(output_path, output, [
        "cohort", "subject", "modality", "population", "input_url", "index_url",
        "input_md5", "genome_build", "selection_seed", "primary_confirmation_eligible",
    ])
    result = {
        "schema_version": "frozen-extension-rosters-1",
        "truth_blind": True,
        "selection_seed": seed,
        "all_minimums_met": all_minimums_met,
        "modalities": summaries,
        "input_sha256": {
            "registry": sha256(registry_path), "assay_manifest": sha256(assay_manifest_path),
            "exposure": sha256(exposure_path), "config": sha256(config_path),
        },
        "output_sha256": sha256(output_path),
    }
    write_json(summary_path, result)
    return result
