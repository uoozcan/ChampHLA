from __future__ import annotations

from collections import defaultdict
from pathlib import Path

from .io import canonical_pair, is_callable, read_json, read_tsv, sha256, write_json, write_tsv
from .panels import METHOD_BASELINE, METHOD_GUARDED_CC, METHOD_PLURALITY, MODALITIES
from .statistics import cluster_bootstrap_ci, exact_cluster_signflip, holm_adjust, subject_deltas
from .freeze import validate_freeze


def prediction_key(row: dict[str, str]) -> tuple[str, str, str, str, str]:
    return (row.get("cohort", ""), row["subject"], row["modality"], row["gene"], row["method"])


def truth_key(row: dict[str, str]) -> tuple[str, str, str, str]:
    return (row.get("cohort", ""), row["subject"], row["modality"], row["gene"])


def join_truth(predictions_path: str, truth_path: str, freeze_manifest_path: str,
               output_path: str, join_manifest_path: str) -> dict:
    freeze = read_json(freeze_manifest_path)
    validation = validate_freeze(freeze_manifest_path)
    if not validation["valid"]:
        raise ValueError(f"confirmation freeze is invalid: {validation['failures']}")
    prediction_hash = sha256(predictions_path)
    frozen_prediction = freeze.get("inputs", {}).get("predictions")
    if not frozen_prediction or prediction_hash != frozen_prediction.get("sha256"):
        raise ValueError("prediction file is absent from the freeze or its checksum changed")
    predictions = read_tsv(predictions_path)
    truth_rows = read_tsv(truth_path)
    required = {"truth_allele1", "truth_allele2", "truth_status", "truth_source", "source_sha256"}
    if not truth_rows or not required.issubset(truth_rows[0]):
        raise ValueError(f"truth table missing required columns: {sorted(required - set(truth_rows[0] if truth_rows else {}))}")
    truth = {}
    for row in truth_rows:
        key = truth_key(row)
        if key in truth:
            raise ValueError(f"duplicate truth locus: {key}")
        truth[key] = row
    joined = []
    prediction_loci = {prediction_key(row)[:4] for row in predictions}
    if prediction_loci != set(truth):
        raise ValueError(f"prediction/truth locus mismatch predictions={len(prediction_loci)} truth={len(truth)}")
    for row in predictions:
        t = truth[prediction_key(row)[:4]]
        joined.append({**row, "truth_allele1": t["truth_allele1"],
                       "truth_allele2": t["truth_allele2"], "truth_status": t["truth_status"],
                       "truth_source": t["truth_source"], "truth_source_sha256": t["source_sha256"]})
    write_tsv(output_path, joined)
    payload = {
        "schema_version": "truth-join-1", "join_once": True, "loci": len(truth),
        "prediction_rows": len(predictions), "prediction_sha256": prediction_hash,
        "truth_sha256": sha256(truth_path), "joined_sha256": sha256(output_path),
        "freeze_manifest_sha256": sha256(freeze_manifest_path),
    }
    write_json(join_manifest_path, payload)
    return payload


def _evaluate_rows(joined: list[dict[str, str]]) -> list[dict]:
    output = []
    for row in joined:
        truth = canonical_pair(row["truth_allele1"], row["truth_allele2"], row["gene"])
        callable_ = is_callable(row)
        pair = canonical_pair(row["allele1"], row["allele2"], row["gene"]) if callable_ else None
        output.append({**row, "is_callable": int(callable_), "is_correct": int(pair == truth) if callable_ else 0})
    return output


def evaluate(joined_path: str, output_dir: str, bootstrap: int = 100000,
             seed: int = 20260831, allow_partial: bool = False,
             mode: str = "discovery", wgs_audit_summary: str | None = None) -> dict:
    evaluated = _evaluate_rows(read_tsv(joined_path))
    by_locus_method = {}
    for row in evaluated:
        key = prediction_key(row)
        if key in by_locus_method:
            raise ValueError(f"duplicate prediction method/locus row: {key}")
        by_locus_method[key] = row
    methods = {row["method"] for row in evaluated}
    required_methods = {METHOD_BASELINE, METHOD_GUARDED_CC, METHOD_PLURALITY}
    if not required_methods.issubset(methods):
        raise ValueError(f"missing required methods: {sorted(required_methods - methods)}")
    loci = sorted({key[:4] for key in by_locus_method})
    paired = []
    for locus in loci:
        rows = {method: by_locus_method.get((*locus, method)) for method in required_methods}
        if any(value is None for value in rows.values()):
            raise ValueError(f"incomplete required method set at {locus}")
        baseline, guarded, plurality = rows[METHOD_BASELINE], rows[METHOD_GUARDED_CC], rows[METHOD_PLURALITY]
        paired.append({
            "cohort": locus[0], "subject": locus[1], "modality": locus[2], "gene": locus[3],
            "superpopulation": guarded.get("superpopulation", ""),
            "truth_homozygous": int(guarded["truth_allele1"] == guarded["truth_allele2"]),
            "baseline_callable": baseline["is_callable"], "baseline_correct": baseline["is_correct"],
            "guarded_correct": guarded["is_correct"], "plurality_correct": plurality["is_correct"],
            "baseline_pair": f"{baseline['allele1']}+{baseline['allele2']}" if int(baseline["is_callable"]) else "",
            "guarded_pair": f"{guarded['allele1']}+{guarded['allele2']}",
            "plurality_pair": f"{plurality['allele1']}+{plurality['allele2']}",
            "guarded_decision_reason": guarded.get("decision_reason", ""),
        })
    present_modalities = sorted({row["modality"] for row in paired})
    if not allow_partial and set(present_modalities) != set(MODALITIES):
        raise ValueError(f"confirmation requires all modalities; present={present_modalities}")
    primary = []
    for index, modality in enumerate(MODALITIES):
        subset = [row for row in paired if row["modality"] == modality]
        if not subset:
            continue
        deltas = subject_deltas(subset, "baseline_correct", "guarded_correct")
        ci_lo, ci_hi = cluster_bootstrap_ci(subset, "baseline_correct", "guarded_correct",
                                            iterations=bootstrap, seed=seed + index)
        baseline_correct = sum(int(row["baseline_correct"]) for row in subset)
        guarded_correct = sum(int(row["guarded_correct"]) for row in subset)
        plurality_correct = sum(int(row["plurality_correct"]) for row in subset)
        primary.append({
            "modality": modality, "subjects": len({row["subject"] for row in subset}), "loci": len(subset),
            "baseline_callable": sum(int(row["baseline_callable"]) for row in subset),
            "baseline_correct": baseline_correct, "guarded_correct": guarded_correct,
            "plurality_correct": plurality_correct,
            "baseline_accuracy": baseline_correct / len(subset),
            "guarded_accuracy": guarded_correct / len(subset),
            "plurality_accuracy": plurality_correct / len(subset),
            "delta_points": 100 * (guarded_correct - baseline_correct) / len(subset),
            "plurality_delta_points": 100 * (guarded_correct - plurality_correct) / len(subset),
            "wins": sum(int(row["guarded_correct"]) > int(row["baseline_correct"]) for row in subset),
            "losses": sum(int(row["guarded_correct"]) < int(row["baseline_correct"]) for row in subset),
            "ties": sum(int(row["guarded_correct"]) == int(row["baseline_correct"]) for row in subset),
            "discordant_correctness_subjects": sum(value != 0 for value in deltas.values()),
            "ci_lo_points": 100 * ci_lo, "ci_hi_points": 100 * ci_hi,
            "p_value": exact_cluster_signflip(deltas.values()),
        })
    primary = holm_adjust(primary)
    gene_tests = []
    for modality in present_modalities:
        for gene in ("A", "B", "C"):
            subset = [row for row in paired if row["modality"] == modality and row["gene"] == gene]
            if not subset:
                continue
            deltas = subject_deltas(subset, "baseline_correct", "guarded_correct")
            delta = sum(deltas.values()) / len(subset)
            gene_tests.append({"modality": modality, "gene": gene, "loci": len(subset),
                               "delta_points": 100 * delta,
                               "p_value": exact_cluster_signflip(deltas.values())})
    gene_tests = holm_adjust(gene_tests)
    for row in gene_tests:
        row["significant_harm"] = int(row["delta_points"] < 0 and row["holm_significant"])
    all_three = set(present_modalities) == set(MODALITIES)
    minimums = {"wgs": 120, "wes": 89, "rnaseq": 130}
    cohort_sizes = {row["modality"]: row["subjects"] for row in primary}
    sample_size_pass = all(cohort_sizes.get(modality, 0) >= minimums[modality] for modality in MODALITIES)
    wgs_audit_pass = False
    if wgs_audit_summary:
        wgs_audit_pass = bool(read_json(wgs_audit_summary).get("passed"))
    if mode not in {"discovery", "external"}:
        raise ValueError(f"unsupported evaluation mode: {mode}")
    if mode == "external" and not wgs_audit_summary:
        raise ValueError("external evaluation requires --wgs-audit-summary")
    primary_pass = bool(all_three and all(
        row["delta_points"] > 0 and row["ci_lo_points"] > 0 and row["holm_significant"] for row in primary))
    plurality_retention = all(row["plurality_delta_points"] >= -2.0 for row in primary)
    no_gene_harm = not any(row["significant_harm"] for row in gene_tests)

    secondary_methods = sorted(methods - {METHOD_GUARDED_CC, METHOD_BASELINE})
    secondary_tests = []
    for modality in present_modalities:
        modality_loci = [locus for locus in loci if locus[2] == modality]
        guarded_map = {locus: int(by_locus_method[(*locus, METHOD_GUARDED_CC)]["is_correct"])
                       for locus in modality_loci}
        for method in secondary_methods:
            if not any((*locus, method) in by_locus_method for locus in modality_loci):
                continue
            subject_delta = defaultdict(int)
            method_correct = guarded_correct = wins = losses = 0
            for locus in modality_loci:
                candidate = by_locus_method.get((*locus, method))
                candidate_correct = int(candidate["is_correct"]) if candidate else 0
                selected_correct = guarded_map[locus]
                method_correct += candidate_correct
                guarded_correct += selected_correct
                wins += int(selected_correct > candidate_correct)
                losses += int(selected_correct < candidate_correct)
                subject_delta[locus[1]] += selected_correct - candidate_correct
            secondary_tests.append({
                "modality": modality, "comparator": method, "loci": len(modality_loci),
                "guarded_correct": guarded_correct, "comparator_correct": method_correct,
                "delta_points": 100 * (guarded_correct - method_correct) / len(modality_loci),
                "wins": wins, "losses": losses,
                "p_value": exact_cluster_signflip(subject_delta.values()),
            })
    secondary_tests = holm_adjust(secondary_tests)

    secondary_endpoints = []
    for modality in present_modalities:
        subset = [row for row in paired if row["modality"] == modality]
        called = [row for row in subset if int(row["baseline_callable"])]
        caller_methods = [method for method in methods if method.startswith("Caller:")]
        oracle = 0
        for row in subset:
            locus = (row["cohort"], row["subject"], row["modality"], row["gene"])
            oracle += int(any(by_locus_method.get((*locus, method), {}).get("is_correct") == 1
                              for method in caller_methods))
        secondary_endpoints.append({
            "modality": modality, "loci": len(subset), "baseline_callable": len(called),
            "baseline_callability": len(called) / len(subset),
            "baseline_called_only_accuracy": (sum(int(row["baseline_correct"]) for row in called) / len(called)
                                                if called else 0.0),
            "guarded_accuracy_at_baseline_coverage": (sum(int(row["guarded_correct"]) for row in called) / len(called)
                                                       if called else 0.0),
            "candidate_oracle_correct": oracle, "candidate_oracle_accuracy": oracle / len(subset),
            "protected_loci": sum(row["guarded_decision_reason"] == "protected_two_thirds_consensus"
                                  for row in subset),
            "challenged_loci": sum(row["guarded_decision_reason"] == "champion_challenger_resolved_disagreement"
                                   for row in subset),
        })

    stratum_rows = []
    for field in ("gene", "superpopulation", "truth_homozygous"):
        values = sorted({str(row[field]) for row in paired})
        for value in values:
            subset = [row for row in paired if str(row[field]) == value]
            if not subset:
                continue
            stratum_rows.append({
                "stratum_type": field, "stratum": value, "loci": len(subset),
                "baseline_correct": sum(int(row["baseline_correct"]) for row in subset),
                "guarded_correct": sum(int(row["guarded_correct"]) for row in subset),
                "delta_points": 100 * sum(int(row["guarded_correct"]) - int(row["baseline_correct"])
                                          for row in subset) / len(subset),
            })
    external_prerequisites = mode == "external" and wgs_audit_pass and sample_size_pass
    headline = bool(external_prerequisites and primary_pass and plurality_retention and no_gene_harm)
    final = {
        "schema_version": "three-modality-confirmation-evaluation-1",
        "mode": mode,
        "status": "three_modality_confirmed" if headline else
                  "discovery_only" if mode == "discovery" else "confirmation_not_achieved",
        "primary_pass": primary_pass, "plurality_retention_pass": plurality_retention,
        "per_gene_harm_pass": no_gene_harm,
        "wgs_audit_pass": wgs_audit_pass, "sample_size_pass": sample_size_pass,
        "cohort_subjects": cohort_sizes, "minimum_subjects": minimums,
        "external_prerequisites_pass": external_prerequisites,
        "headline_retained": headline,
        "claim_boundary": "resolution of two-thirds consensus failures at fixed denominator",
        "input_sha256": sha256(joined_path),
    }
    out = Path(output_dir)
    write_tsv(out / "paired_locus_audit.tsv", paired)
    write_tsv(out / "primary_modality_results.tsv", primary)
    write_tsv(out / "per_gene_results.tsv", gene_tests)
    write_tsv(out / "secondary_method_results.tsv", secondary_tests)
    write_tsv(out / "secondary_endpoints.tsv", secondary_endpoints)
    write_tsv(out / "stratified_results.tsv", stratum_rows)
    write_json(out / "confirmation_decision.json", final)
    return final


def capacity(predictions_path: str, output_json: str) -> dict:
    rows = read_tsv(predictions_path)
    by_key = {(row.get("cohort", ""), row["subject"], row["modality"], row["gene"], row["method"]): row
              for row in rows}
    loci = sorted({key[:4] for key in by_key})
    counts = defaultdict(set)
    for locus in loci:
        baseline = by_key.get((*locus, METHOD_BASELINE))
        guarded = by_key.get((*locus, METHOD_GUARDED_CC))
        if not baseline or not guarded:
            continue
        baseline_pair = (baseline.get("allele1", ""), baseline.get("allele2", ""), baseline.get("call_status", ""))
        guarded_pair = (guarded.get("allele1", ""), guarded.get("allele2", ""), guarded.get("call_status", ""))
        if baseline_pair != guarded_pair:
            counts[locus[2]].add(locus[1])
    values = {modality: len(counts[modality]) for modality in MODALITIES}
    result = {
        "schema_version": "truth-blind-capacity-gate-1", "discordant_subjects": values,
        "minimum_six_each": all(values[m] >= 6 for m in MODALITIES),
        "minimum_seven_in_two": sum(values[m] >= 7 for m in MODALITIES) >= 2,
    }
    result["passed"] = result["minimum_six_each"] and result["minimum_seven_in_two"]
    write_json(output_json, result)
    return result
