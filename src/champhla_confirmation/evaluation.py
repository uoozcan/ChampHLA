from __future__ import annotations

from collections import Counter, defaultdict
from pathlib import Path

from .io import canonical_allele, canonical_pair, is_callable, read_json, read_tsv, sha256, write_json, write_tsv
from .panels import (
    METHOD_BASELINE,
    METHOD_GUARDED_CC,
    METHOD_PLURALITY,
    MODALITIES,
    canonical_method,
)
from .statistics import (
    cluster_bootstrap_ci,
    exact_cluster_signflip,
    holm_adjust,
    simultaneous_cluster_bootstrap_ci,
    subject_deltas,
)
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
        predicted_alleles = []
        for field in ("allele1", "allele2"):
            if row.get(field, "").strip() and row[field].strip() not in {"-", "."}:
                predicted_alleles.append(canonical_allele(row[field], row["gene"]))
        allele_correct = sum((Counter(predicted_alleles) & Counter(truth)).values())
        cluster_id = row.get("donor", "").strip() or f"{row.get('cohort', '')}:{row['subject']}"
        output.append({**row, "source_method": row["method"],
                       "method": canonical_method(row["method"]),
                       "cluster_id": cluster_id,
                       "is_callable": int(callable_),
                       "is_correct": int(pair == truth) if callable_ else 0,
                       "alleles_correct": allele_correct})
    return output


def _index_evaluated(rows: list[dict]) -> dict[tuple[str, str, str, str, str], dict]:
    """Index predictions, accepting only an identical documented plurality alias."""
    indexed = {}
    for row in rows:
        key = prediction_key(row)
        if key in indexed:
            prior = indexed[key]
            same = all(prior.get(field, "") == row.get(field, "") for field in (
                "allele1", "allele2", "is_correct", "is_callable",
            ))
            alias_pair = {prior.get("source_method"), row.get("source_method")}
            if same and alias_pair == {METHOD_PLURALITY, "MajorityVote"}:
                if row.get("source_method") == METHOD_PLURALITY:
                    indexed[key] = row
                continue
            raise ValueError(f"duplicate prediction method/locus row: {key}")
        indexed[key] = row
    return indexed


def evaluate_legacy_guarded_cc(joined_path: str, output_dir: str, bootstrap: int = 100000,
             seed: int = 20260831, allow_partial: bool = False,
             mode: str = "discovery", wgs_audit_summary: str | None = None) -> dict:
    evaluated = _evaluate_rows(read_tsv(joined_path))
    by_locus_method = _index_evaluated(evaluated)
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
        "bootstrap_iterations": bootstrap,
        "random_seed": seed,
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


def _comparison_design(path: str | None, present_modalities: list[str]) -> dict:
    design = read_json(path) if path else {}
    valid = design.get("valid_modalities", present_modalities)
    invalid = design.get("invalid_modalities", {})
    if set(valid) & set(invalid):
        raise ValueError("evaluation design marks a modality both valid and invalid")
    return {
        "reference_method": canonical_method(design.get("reference_method", METHOD_PLURALITY)),
        "valid_modalities": [value for value in MODALITIES if value in valid],
        "invalid_modalities": invalid,
        "noninferiority_margin_points": float(design.get("noninferiority_margin_points", 2.0)),
        "minimum_subjects": design.get("minimum_subjects", {"wgs": 120, "wes": 89, "rnaseq": 130}),
        "prospective": bool(design.get("prospective", False)),
        "evidence_tier": design.get(
            "evidence_tier", "independent_validation" if design.get("prospective") else "development"
        ),
        "primary_independence_stratum": design.get("primary_independence_stratum", ""),
    }


def _manifest_families(path: str | None, methods: set[str], modalities: list[str],
                       reference_method: str) -> dict[str, list[str]]:
    if not path:
        candidates = sorted(methods - {reference_method})
        return {modality: candidates for modality in modalities}
    manifest = read_json(path)
    entries = manifest.get("comparators", [])
    if not entries:
        raise ValueError("comparator manifest has no comparators")
    families = {modality: [] for modality in modalities}
    seen = set()
    for entry in entries:
        method = canonical_method(entry["method_id"])
        if method in seen:
            raise ValueError(f"duplicate comparator manifest method: {method}")
        seen.add(method)
        if not entry.get("deployable", False) or not entry.get("primary_family", False):
            continue
        for modality in entry.get("modalities", []):
            if modality in families and method != reference_method:
                families[modality].append(method)
    return {key: sorted(value) for key, value in families.items()}


def evaluate(joined_path: str, output_dir: str, bootstrap: int = 100000,
             seed: int = 20260831, allow_partial: bool = False,
             mode: str = "discovery", wgs_audit_summary: str | None = None,
             evaluation_design: str | None = None,
             comparator_manifest: str | None = None) -> dict:
    """Evaluate pair-level plurality against frozen modality-specific families.

    This v2 evaluator is intentionally method-generic. The legacy Guarded-CC
    evaluator remains available as ``evaluate_legacy_guarded_cc`` for archived
    artifacts, but cannot authorize the plurality-centered manuscript.
    """
    if mode not in {"discovery", "external"}:
        raise ValueError(f"unsupported evaluation mode: {mode}")
    evaluated = _evaluate_rows(read_tsv(joined_path))
    by_locus_method = _index_evaluated(evaluated)
    methods = {row["method"] for row in evaluated}
    loci = sorted({key[:4] for key in by_locus_method})
    present_modalities = sorted({locus[2] for locus in loci})
    if not allow_partial and set(present_modalities) != set(MODALITIES):
        raise ValueError(f"confirmation requires all modalities; present={present_modalities}")
    design = _comparison_design(evaluation_design, present_modalities)
    reference_method = design["reference_method"]
    if reference_method not in methods:
        raise ValueError(f"missing reference method: {reference_method}")
    families = _manifest_families(
        comparator_manifest, methods, present_modalities, reference_method,
    )

    reference_rows = []
    for locus in loci:
        row = by_locus_method.get((*locus, reference_method))
        if row is None:
            raise ValueError(f"reference method missing at {locus}")
        reference_rows.append(row)

    strata = sorted({row.get("independence_stratum", "").strip() or "unspecified"
                     for row in reference_rows})
    primary_stratum = design["primary_independence_stratum"]
    if len(strata) > 1 and not primary_stratum:
        raise ValueError(
            "multiple independence strata cannot be pooled; set primary_independence_stratum"
        )
    if primary_stratum and primary_stratum not in strata:
        raise ValueError(f"primary independence stratum is absent: {primary_stratum}")
    primary_stratum = primary_stratum or strata[0]
    primary_reference_rows = [
        row for row in reference_rows
        if (row.get("independence_stratum", "").strip() or "unspecified") == primary_stratum
    ]
    primary_loci = {
        (row.get("cohort", ""), row["subject"], row["modality"], row["gene"])
        for row in primary_reference_rows
    }

    primary = []
    secondary_endpoints = []
    for modality in present_modalities:
        subset = [row for row in primary_reference_rows if row["modality"] == modality]
        if not subset:
            continue
        correct = sum(int(row["is_correct"]) for row in subset)
        alleles_correct = sum(int(row["alleles_correct"]) for row in subset)
        callable_count = sum(int(row["is_callable"]) for row in subset)
        partial_count = sum(row.get("call_status", "").lower() == "partial" for row in subset)
        subset_loci = [
            (row.get("cohort", ""), row["subject"], row["modality"], row["gene"])
            for row in subset
        ]
        caller_methods = sorted(method for method in methods if method.startswith("Caller:"))
        caller_totals = {
            method: sum(int(by_locus_method.get((*locus, method), {}).get("is_correct", 0))
                        for locus in subset_loci)
            for method in caller_methods
            if any((*locus, method) in by_locus_method for locus in subset_loci)
        }
        best_caller, best_caller_correct = (
            sorted(caller_totals.items(), key=lambda item: (-item[1], item[0]))[0]
            if caller_totals else ("", 0)
        )
        oracle_correct = sum(
            any(int(by_locus_method.get((*locus, method), {}).get("is_correct", 0))
                for method in caller_methods)
            for locus in subset_loci
        )
        primary.append({
            "modality": modality,
            "validity": "valid" if modality in design["valid_modalities"] else "invalid",
            "reference_method": reference_method,
            "independence_stratum": primary_stratum,
            "subjects": len({row["cluster_id"] for row in subset}),
            "loci": len(subset),
            "correct": correct,
            "accuracy": correct / len(subset),
            "alleles_correct": alleles_correct,
            "allele_accuracy": alleles_correct / (2 * len(subset)),
            "callable": callable_count,
            "call_rate": callable_count / len(subset),
            "partial": partial_count,
            "partial_rate": partial_count / len(subset),
        })
        secondary_endpoints.append({
            "modality": modality,
            "loci": len(subset),
            "reference_method": reference_method,
            "callable": callable_count,
            "partial": partial_count,
            "missing_or_no_evidence": len(subset) - callable_count - partial_count,
            "ties": sum(int(row.get("tie_at_top", 0) or 0) for row in subset),
            "called_only_accuracy": (correct / callable_count if callable_count else 0.0),
            "allele_accuracy": alleles_correct / (2 * len(subset)),
            "best_observed_individual_caller": best_caller,
            "best_observed_individual_caller_correct": best_caller_correct,
            "best_observed_individual_caller_accuracy": best_caller_correct / len(subset),
            "best_caller_selected_after_observing_results": 1,
            "candidate_oracle_correct": oracle_correct,
            "candidate_oracle_accuracy": oracle_correct / len(subset),
        })

    comparison_rows = []
    paired_rows = []
    incomplete = []
    noninferiority_by_modality = {}
    for modality_index, modality in enumerate(present_modalities):
        modality_loci = [locus for locus in sorted(primary_loci) if locus[2] == modality]
        family_pairs = {}
        raw_records = []
        for comparator in families.get(modality, []):
            rows = []
            missing = 0
            for locus in modality_loci:
                reference = by_locus_method[(*locus, reference_method)]
                candidate = by_locus_method.get((*locus, comparator))
                if candidate is None:
                    missing += 1
                    candidate_correct = candidate_callable = 0
                    candidate_status = "missing_prediction_row"
                else:
                    candidate_correct = int(candidate["is_correct"])
                    candidate_callable = int(candidate["is_callable"])
                    candidate_status = candidate.get("call_status", "")
                audit = {
                    "cohort": locus[0], "subject": locus[1], "modality": modality,
                    "gene": locus[3], "reference_method": reference_method,
                    "cluster_id": reference["cluster_id"],
                    "comparator": comparator,
                    "reference_correct": int(reference["is_correct"]),
                    "comparator_correct": candidate_correct,
                    "delta": candidate_correct - int(reference["is_correct"]),
                    "reference_callable": int(reference["is_callable"]),
                    "comparator_callable": candidate_callable,
                    "comparator_status": candidate_status,
                }
                rows.append(audit)
                paired_rows.append(audit)
            if missing:
                incomplete.append(f"{modality}:{comparator}:{missing}")
            family_pairs[comparator] = rows
        intervals = simultaneous_cluster_bootstrap_ci(
            family_pairs, "reference_correct", "comparator_correct",
            iterations=bootstrap, seed=seed + modality_index * 1000,
        )
        for comparator, rows in sorted(family_pairs.items()):
            reference_correct = sum(row["reference_correct"] for row in rows)
            comparator_correct = sum(row["comparator_correct"] for row in rows)
            deltas = subject_deltas(rows, "reference_correct", "comparator_correct")
            ci_lo, ci_hi = intervals[comparator]
            raw_records.append({
                "modality": modality,
                "reference_method": reference_method,
                "comparator": comparator,
                "loci": len(rows),
                "reference_correct": reference_correct,
                "comparator_correct": comparator_correct,
                "comparator_minus_reference_points": 100 * (comparator_correct - reference_correct) / len(rows),
                "simultaneous_ci_lo_points": 100 * ci_lo,
                "simultaneous_ci_hi_points": 100 * ci_hi,
                "wins": sum(row["delta"] > 0 for row in rows),
                "losses": sum(row["delta"] < 0 for row in rows),
                "missing_prediction_rows": missing,
                "comparison_valid": int(missing == 0),
                "p_value": exact_cluster_signflip(deltas.values()),
            })
        adjusted = holm_adjust(raw_records)
        margin = design["noninferiority_margin_points"]
        for record in adjusted:
            record["holm_significant_superiority"] = int(
                record["comparator_minus_reference_points"] > 0 and record["holm_significant"]
            )
            record["reference_noninferior_2pp"] = int(
                record["simultaneous_ci_hi_points"] < margin
            )
            comparison_rows.append(record)
        noninferiority_by_modality[modality] = bool(
            adjusted and all(row["comparison_valid"] and row["reference_noninferior_2pp"]
                             for row in adjusted)
            and not any(item.startswith(f"{modality}:") for item in incomplete)
        )

    valid_comparisons = [row for row in comparison_rows if row["modality"] in design["valid_modalities"]]
    no_comparator_holm_superior = bool(
        valid_comparisons
        and not any(item.split(":", 1)[0] in design["valid_modalities"] for item in incomplete)
        and not any(row["holm_significant_superiority"] for row in valid_comparisons)
    )
    consensus_noninferior = bool(
        design["valid_modalities"] and all(
            noninferiority_by_modality.get(modality, False)
            for modality in design["valid_modalities"]
        )
    )

    gene_results = []
    for modality in present_modalities:
        for gene in ("A", "B", "C"):
            subset = [row for row in primary_reference_rows
                      if row["modality"] == modality and row["gene"] == gene]
            if subset:
                correct = sum(int(row["is_correct"]) for row in subset)
                gene_results.append({
                    "modality": modality, "gene": gene, "reference_method": reference_method,
                    "loci": len(subset), "correct": correct, "accuracy": correct / len(subset),
                })

    stratum_rows = []
    for field in ("gene", "superpopulation"):
        for value in sorted({str(row.get(field, "")) for row in primary_reference_rows}):
            subset = [row for row in primary_reference_rows if str(row.get(field, "")) == value]
            if subset:
                correct = sum(int(row["is_correct"]) for row in subset)
                stratum_rows.append({
                    "stratum_type": field, "stratum": value, "loci": len(subset),
                    "reference_method": reference_method, "correct": correct,
                    "accuracy": correct / len(subset),
                })

    external_strata = []
    for modality in present_modalities:
        for stratum in strata:
            subset = [
                row for row in reference_rows
                if row["modality"] == modality
                and (row.get("independence_stratum", "").strip() or "unspecified") == stratum
            ]
            if subset:
                external_strata.append({
                    "modality": modality,
                    "independence_stratum": stratum,
                    "subjects": len({row["cluster_id"] for row in subset}),
                    "loci": len(subset),
                    "correct": sum(int(row["is_correct"]) for row in subset),
                    "accuracy": sum(int(row["is_correct"]) for row in subset) / len(subset),
                    "pooled_into_primary": int(stratum == primary_stratum),
                })

    caller_status = []
    caller_rows = [row for row in evaluated if row["method"].startswith("Caller:")]
    for modality in present_modalities:
        for method in sorted({row["method"] for row in caller_rows if row["modality"] == modality}):
            for gene in ("A", "B", "C"):
                subset = [row for row in caller_rows
                          if row["modality"] == modality and row["method"] == method
                          and row["gene"] == gene]
                if not subset:
                    continue
                caller_status.append({
                    "modality": modality, "gene": gene, "caller_method": method,
                    "records": len(subset),
                    "complete": sum(int(row["is_callable"]) for row in subset),
                    "partial": sum(row.get("call_status", "").lower() == "partial" for row in subset),
                    "missing": sum(row.get("call_status", "").lower() in {
                        "missing", "no_evidence", "no_consensus"} for row in subset),
                })

    cohort_sizes = {row["modality"]: row["subjects"] for row in primary}
    minimums = {key: int(value) for key, value in design["minimum_subjects"].items()}
    all_three_valid = set(design["valid_modalities"]) == set(MODALITIES)
    sample_size_pass = all(cohort_sizes.get(modality, 0) >= minimums[modality] for modality in MODALITIES)
    wgs_audit_pass = bool(wgs_audit_summary and read_json(wgs_audit_summary).get("passed"))
    same_resource_ready = bool(
        mode == "external" and design["evidence_tier"] == "same_resource_confirmation"
        and all_three_valid and sample_size_pass and wgs_audit_pass and not incomplete
    )
    three_modality_ready = bool(
        mode == "external" and design["evidence_tier"] == "independent_validation"
        and design["prospective"] and all_three_valid
        and sample_size_pass and wgs_audit_pass and no_comparator_holm_superior
        and all(noninferiority_by_modality.get(modality, False) for modality in MODALITIES)
        and not incomplete
    )
    final = {
        "schema_version": "plurality-consensus-evaluation-2",
        "mode": mode,
        "status": "three_modality_confirmed" if three_modality_ready else
                  "same_resource_benchmark_complete" if same_resource_ready else
                  "development_only" if mode == "discovery" else "confirmation_not_achieved",
        "reference_method": reference_method,
        "evidence_tier": design["evidence_tier"],
        "primary_independence_stratum": primary_stratum,
        "observed_independence_strata": strata,
        "valid_modalities": design["valid_modalities"],
        "invalid_modalities": design["invalid_modalities"],
        "no_comparator_holm_superior": no_comparator_holm_superior,
        "consensus_noninferior_2pp": consensus_noninferior,
        "noninferiority_margin_points": design["noninferiority_margin_points"],
        "noninferiority_by_modality": noninferiority_by_modality,
        "three_modality_claim_ready": three_modality_ready,
        "same_resource_benchmark_ready": same_resource_ready,
        "prospective_design": design["prospective"],
        "wgs_audit_pass": wgs_audit_pass,
        "sample_size_pass": sample_size_pass,
        "cohort_subjects": cohort_sizes,
        "minimum_subjects": minimums,
        "incomplete_comparator_rows": incomplete,
        "claim_boundary": (
            "No evaluated comparator significantly exceeded pair-level plurality; "
            "noninferiority requires the separate prospective two-point gate."
        ),
        "input_sha256": sha256(joined_path),
        "bootstrap_iterations": bootstrap,
        "random_seed": seed,
        "evaluation_design": evaluation_design or "built_in_defaults",
        "evaluation_design_sha256": sha256(evaluation_design) if evaluation_design else "",
        "comparator_manifest": comparator_manifest or "all_observed_methods",
        "comparator_manifest_sha256": sha256(comparator_manifest) if comparator_manifest else "",
    }
    out = Path(output_dir)
    write_tsv(out / "paired_locus_audit.tsv", paired_rows)
    write_tsv(out / "primary_modality_results.tsv", primary)
    write_tsv(out / "head_to_head.tsv", comparison_rows)
    write_tsv(out / "secondary_method_results.tsv", comparison_rows)
    write_tsv(out / "per_gene_results.tsv", gene_results)
    write_tsv(out / "secondary_endpoints.tsv", secondary_endpoints)
    write_tsv(out / "stratified_results.tsv", stratum_rows)
    write_tsv(out / "external_strata_results.tsv", external_strata)
    write_tsv(out / "caller_call_status.tsv", caller_status)
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
