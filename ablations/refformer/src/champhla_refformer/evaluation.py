"""Subject-grouped nested evaluation and predeclared RefFormer decision gates."""

from __future__ import annotations

import hashlib
import json
import math
import random
from collections import Counter, defaultdict
from pathlib import Path

import torch

from . import GENES, MODALITIES
from .data import build_loci, fit_runtime_context, read_tsv, write_tsv
from .training import fit_temperature, predict_loci, save_bundle, train_model


BASELINES = ("MajorityVote", "ChampionChallenger", "MVFlooredCC",
             "BestSingleTool_nestedCV", "MetaConsensus")


def stable_int(text: str) -> int:
    return int(hashlib.sha256(text.encode()).hexdigest()[:12], 16)


def assign_grouped_folds(rows: list[dict], folds: int, seed: int) -> dict[str, int]:
    populations = {}
    for row in rows:
        populations.setdefault(row["sample"], row.get("superpopulation") or "unknown")
    grouped = defaultdict(list)
    for subject, population in populations.items():
        grouped[population].append(subject)
    result = {}
    for population, subjects in sorted(grouped.items()):
        subjects.sort(key=lambda subject: stable_int(f"{seed}|{population}|{subject}"))
        for index, subject in enumerate(subjects):
            result[subject] = index % folds
    return result


def rate(rows: list[dict]) -> float:
    return sum(str(row.get("is_correct", "0")) == "1" for row in rows) / max(1, len(rows))


def confidence(row: dict) -> float:
    try:
        return max(0.0, min(1.0, float(row.get("confidence", 0.0))))
    except (TypeError, ValueError):
        return 0.0


def secondary_metrics(rows: list[dict]) -> dict:
    if not rows:
        return {"callability": 0.0, "aurc": 0.0, "brier": "", "ece": "",
                "accuracy_at_95pct_coverage": 0.0, "accuracy_at_90pct_coverage": 0.0,
                "accuracy_at_80pct_coverage": 0.0}
    ordered = sorted(rows, key=lambda row: (-confidence(row), row["sample"], row["gene"]))
    errors, risks = 0, []
    for index, row in enumerate(ordered, 1):
        errors += int(str(row.get("is_correct", "0")) != "1")
        risks.append(errors / index)
    output = {"callability": sum(str(row.get("is_callable", "1")) == "1" for row in rows) / len(rows),
              "aurc": sum(risks) / len(risks)}
    for coverage in (0.95, 0.90, 0.80):
        count = max(1, math.ceil(coverage * len(rows)))
        output[f"accuracy_at_{int(coverage * 100)}pct_coverage"] = rate(ordered[:count])
    probabilities = [confidence(row) for row in rows]
    labels = [int(str(row.get("is_correct", "0")) == "1") for row in rows]
    output["brier"] = sum((p - y) ** 2 for p, y in zip(probabilities, labels)) / len(rows)
    ece = 0.0
    for index in range(10):
        members = [j for j, value in enumerate(probabilities)
                   if index / 10 <= value < (index + 1) / 10 or (index == 9 and value == 1.0)]
        if members:
            ece += len(members) / len(rows) * abs(sum(probabilities[j] for j in members) / len(members) -
                                                   sum(labels[j] for j in members) / len(members))
    output["ece"] = ece
    return output


def call_map(rows: list[dict]) -> dict:
    return {(row["sample"], row["modality"], row["gene"]): int(str(row.get("is_correct", "0")) == "1")
            for row in rows}


def fixed_universe(calls: list[dict], rows: list[dict], modality: str, method: str) -> list[dict]:
    truth = {}
    for row in rows:
        if row.get("modality") == modality and row.get("gene") in GENES:
            truth[(row["sample"], row["gene"])] = (row.get("truth_allele1", ""), row.get("truth_allele2", ""))
    observed = {(row["sample"], row["gene"]): row for row in calls if row.get("modality") == modality}
    output = []
    for (sample, gene), truth_pair in sorted(truth.items()):
        if (sample, gene) in observed:
            output.append(observed[(sample, gene)])
        else:
            output.append({"sample": sample, "modality": modality, "gene": gene, "method": method,
                           "allele1": "", "allele2": "", "is_correct": "0", "is_correct_2field": "0",
                           "is_callable": "0", "call_status": "no_call", "confidence": 0.0,
                           "abstention_reason": "no_candidate"})
    return output


def cluster_bootstrap_delta(selected: list[dict], reference: list[dict], iterations: int, seed: int):
    first, second = call_map(selected), call_map(reference)
    subjects = sorted({key[0] for key in first} & {key[0] for key in second})
    per_subject = defaultdict(list)
    for key in set(first) & set(second):
        per_subject[key[0]].append(first[key] - second[key])
    rng, values = random.Random(seed), []
    for _ in range(iterations):
        drawn = [rng.choice(subjects) for _ in subjects]
        differences = [value for subject in drawn for value in per_subject[subject]]
        values.append(sum(differences) / max(1, len(differences)))
    values.sort()
    return values[int(0.025 * (len(values) - 1))], values[int(0.975 * (len(values) - 1))]


def cluster_permutation_p(selected: list[dict], reference: list[dict], iterations: int, seed: int):
    first, second = call_map(selected), call_map(reference)
    per_subject = defaultdict(int)
    for key in set(first) & set(second):
        per_subject[key[0]] += first[key] - second[key]
    values = list(per_subject.values())
    observed = abs(sum(values))
    if observed == 0:
        return 1.0
    rng, extreme = random.Random(seed), 0
    for _ in range(iterations):
        statistic = abs(sum(value if rng.random() < 0.5 else -value for value in values))
        extreme += int(statistic >= observed)
    return (extreme + 1) / (iterations + 1)


def mcnemar_exact(first: list[dict], second: list[dict]):
    a, b = call_map(first), call_map(second)
    better = sum(a[key] == 1 and b[key] == 0 for key in set(a) & set(b))
    worse = sum(a[key] == 0 and b[key] == 1 for key in set(a) & set(b))
    total = better + worse
    if total == 0:
        return better, worse, 1.0
    tail = sum(math.comb(total, index) for index in range(0, min(better, worse) + 1)) / (2 ** total)
    return better, worse, min(1.0, 2 * tail)


def holm(rows: list[dict]):
    ordered = sorted(enumerate(rows), key=lambda item: item[1]["mcnemar_p_raw"])
    adjusted, running = [1.0] * len(rows), 0.0
    for rank, (index, row) in enumerate(ordered):
        running = max(running, min(1.0, (len(rows) - rank) * row["mcnemar_p_raw"]))
        adjusted[index] = running
    for row, value in zip(rows, adjusted):
        row["mcnemar_p_holm"] = value


def inner_predictions(rows: list[dict], reference, store, config: dict, folds: int, seed: int,
                      no_imgt: bool):
    fold_of = assign_grouped_folds(rows, folds, seed)
    calls, raw, summaries = [], [], []
    for fold in range(folds):
        heldout = {subject for subject, value in fold_of.items() if value == fold}
        training_rows = [row for row in rows if row["sample"] not in heldout]
        test_rows = [row for row in rows if row["sample"] in heldout]
        context = fit_runtime_context(training_rows)
        train_loci = build_loci(training_rows, context, reference, include_labels=True)
        test_loci = build_loci(test_rows, context, reference, include_labels=True)
        model, summary = train_model(train_loci, test_loci, store, config, seed + fold, no_imgt)
        fold_calls, _ = predict_loci(model, test_loci, store)
        calls.extend(fold_calls)
        with torch.no_grad():
            for locus in test_loci:
                if locus.get("candidate_set_oracle"):
                    scores, _ = model.forward_locus(locus, store)
                    truth_index = next(index for index, item in enumerate(locus["candidates"]) if item.get("label") == 1)
                    raw.append((scores.detach().cpu().tolist(), truth_index))
        summaries.append(summary)
    return calls, raw, summaries


def nested_run(rows: list[dict], fold_assignments: dict[str, int], reference, store, config: dict,
               outer_folds: int, inner_folds: int, seed: int, outdir: Path,
               outer_fold_ids: list[int] | tuple[int, ...] | None = None):
    all_calls, all_ablation, all_audit, policies, selected_configs = [], [], [], [], []
    grid = config.get("grid") or [{key: config[key] for key in config if key != "grid"}]
    requested_folds = list(range(outer_folds)) if outer_fold_ids is None else list(outer_fold_ids)
    if not requested_folds:
        raise ValueError("at least one outer fold must be requested")
    if len(set(requested_folds)) != len(requested_folds):
        raise ValueError("outer fold IDs must be unique")
    invalid = [fold for fold in requested_folds if fold < 0 or fold >= outer_folds]
    if invalid:
        raise ValueError(f"outer fold IDs outside [0, {outer_folds}): {invalid}")
    for outer in requested_folds:
        test_subjects = {sample for sample, fold in fold_assignments.items() if fold == outer}
        train_rows = [row for row in rows if row["sample"] not in test_subjects]
        test_rows = [row for row in rows if row["sample"] in test_subjects]
        if not train_rows or not test_rows:
            continue
        best = None
        for grid_index, override in enumerate(grid):
            candidate_config = {**config, **override}
            inner_calls, inner_raw, summaries = inner_predictions(
                train_rows, reference, store, candidate_config, inner_folds, seed + outer * 100 + grid_index * 10, False
            )
            score = rate(inner_calls)
            record = (score, -grid_index, candidate_config, inner_calls, inner_raw, summaries)
            if best is None or record[:2] > best[:2]:
                best = record
        _, _, selected_config, inner_calls, inner_raw, summaries = best
        temperature = fit_temperature(inner_raw)
        context = fit_runtime_context(train_rows)
        train_loci = build_loci(train_rows, context, reference, include_labels=True)
        test_loci = build_loci(test_rows, context, reference, include_labels=True)
        epochs = [item["epochs_completed"] for item in summaries]
        final_config = dict(selected_config)
        final_config["epochs"] = max(1, round(sum(epochs) / max(1, len(epochs))))
        final_config["patience"] = final_config["epochs"] + 1
        model, train_summary = train_model(train_loci, [], store, final_config, seed + outer, False)
        fold_calls, fold_audit = predict_loci(model, test_loci, store, temperature)
        ablation_inner, ablation_raw, _ = inner_predictions(
            train_rows, reference, store, final_config, inner_folds, seed + outer * 1000 + 500, True
        )
        ablation_temperature = fit_temperature(ablation_raw)
        ablation_model, ablation_summary = train_model(train_loci, [], store, final_config, seed + 1000 + outer, True)
        ablation_calls, _ = predict_loci(ablation_model, test_loci, store, ablation_temperature)
        fold_calls = [row for modality in MODALITIES
                      for row in fixed_universe(fold_calls, test_rows, modality, "ChampHLA-RefFormer")]
        ablation_calls = [row for modality in MODALITIES
                          for row in fixed_universe(ablation_calls, test_rows, modality, "NoIMGTSetTransformer")]
        for row in fold_calls + fold_audit + ablation_calls:
            row["outer_fold"] = outer
        all_calls.extend(fold_calls); all_audit.extend(fold_audit); all_ablation.extend(ablation_calls)
        policies.append({"outer_fold": outer, "temperature": temperature,
                         "ablation_temperature": ablation_temperature,
                         "inner_accuracy": rate(inner_calls), "ablation_inner_accuracy": rate(ablation_inner),
                         "training_subjects": len({row['sample'] for row in train_rows}),
                         "test_subjects": len(test_subjects), "epochs": final_config["epochs"]})
        selected_configs.append(final_config)
        print(f"[RefFormer] outer fold {outer}: {len(test_subjects)} subjects", flush=True)
    write_tsv(outdir / "refformer_oof_calls.tsv", all_calls)
    write_tsv(outdir / "no_imgt_ablation_oof_calls.tsv", all_ablation)
    write_tsv(outdir / "candidate_audit.tsv", all_audit)
    write_tsv(outdir / "outer_fold_training.tsv", policies)
    return all_calls, all_ablation, all_audit, selected_configs


def summarize(refformer: list[dict], ablation: list[dict], audit: list[dict], baseline_rows: list[dict], bootstrap: int,
              permutation: int, seed: int):
    summaries, decisions, gene_tests = [], [], []
    for modality in MODALITIES:
        ref = [row for row in refformer if row["modality"] == modality]
        abl = [row for row in ablation if row["modality"] == modality]
        candidates = {method: [row for row in baseline_rows if row.get("modality") == modality and
                               row.get("method") == method] for method in BASELINES}
        baseline_rates = {method: rate(rows) for method, rows in candidates.items() if rows}
        strongest = sorted(baseline_rates, key=lambda method: (-baseline_rates[method], method))[0]
        base = candidates[strongest]
        delta = rate(ref) - rate(base)
        lo, hi = cluster_bootstrap_delta(ref, base, bootstrap, seed + len(decisions))
        pvalue = cluster_permutation_p(ref, base, permutation, seed + 100 + len(decisions))
        for method, rows in [("ChampHLA-RefFormer", ref), ("NoIMGTSetTransformer", abl),
                             *[(name, values) for name, values in candidates.items() if values]]:
            secondary = secondary_metrics(rows)
            summaries.append({"modality": modality, "method": method,
                              "n_subjects": len({row["sample"] for row in rows}), "n_gene_loci": len(rows),
                              "overall_correct_call_rate": rate(rows),
                              "callability": secondary["callability"], "aurc": secondary["aurc"],
                              "accuracy_at_95pct_coverage": secondary["accuracy_at_95pct_coverage"],
                              "accuracy_at_90pct_coverage": secondary["accuracy_at_90pct_coverage"],
                              "accuracy_at_80pct_coverage": secondary["accuracy_at_80pct_coverage"],
                              "brier_score": secondary["brier"] if method == "ChampHLA-RefFormer" else "",
                              "expected_calibration_error": secondary["ece"] if method == "ChampHLA-RefFormer" else "",
                              "candidate_set_oracle_ceiling": (sum(int(row.get("candidate_set_oracle", 0) or 0)
                                  for row in rows) / max(1, len(rows))) if method == "ChampHLA-RefFormer" else "",
                              "strongest_baseline": strongest if method == "ChampHLA-RefFormer" else "",
                              "delta_vs_strongest": delta if method == "ChampHLA-RefFormer" else "",
                              "delta_cluster_ci_lo": lo if method == "ChampHLA-RefFormer" else "",
                              "delta_cluster_ci_hi": hi if method == "ChampHLA-RefFormer" else "",
                              "cluster_permutation_p": pvalue if method == "ChampHLA-RefFormer" else ""})
        for gene in GENES:
            first = [row for row in ref if row["gene"] == gene]
            second = [row for row in base if row["gene"] == gene]
            better, worse, p = mcnemar_exact(first, second)
            gene_tests.append({"modality": modality, "gene": gene, "strongest_baseline": strongest,
                               "refformer_accuracy": rate(first), "baseline_accuracy": rate(second),
                               "delta": rate(first) - rate(second), "mcnemar_b": better,
                               "mcnemar_c": worse, "mcnemar_p_raw": p})
        decisions.append({"modality": modality, "strongest_baseline": strongest,
                          "refformer_accuracy": rate(ref), "baseline_accuracy": rate(base),
                          "delta": delta, "delta_ci_lo": lo, "delta_ci_hi": hi,
                          "meaningful_improvement": int(delta >= 0.02 and lo > 0),
                          "within_2pp_nonregression": int(delta >= -0.02),
                          "no_imgt_accuracy": rate(abl), "imgt_delta": rate(ref) - rate(abl),
                          "imgt_modality_nonregression": int(rate(ref) - rate(abl) >= -0.01)})
    holm(gene_tests)
    for decision in decisions:
        decision["significant_gene_harm"] = int(any(row["modality"] == decision["modality"] and
            row["delta"] < 0 and row["mcnemar_p_holm"] < 0.05 for row in gene_tests))
    pooled_imgt_delta = rate(refformer) - rate(ablation)
    passed = any(row["meaningful_improvement"] for row in decisions) and \
        all(row["within_2pp_nonregression"] for row in decisions) and \
        not any(row["significant_gene_harm"] for row in decisions) and \
        pooled_imgt_delta >= 0.01 and all(row["imgt_modality_nonregression"] for row in decisions)
    gate = {"gate_name": "ChampHLA-RefFormer nested candidate-reranker kill test", "passed": passed,
            "decision": "freeze_for_external_validation" if passed else "activate_standalone_reference_ranker_fallback",
            "pooled_imgt_delta_vs_no_imgt": pooled_imgt_delta, "modalities": decisions,
            "external_truth_changed_policy": False, "runtime_default_changed": False}
    topk = []
    for modality in MODALITIES:
        entries = [row for row in audit if row.get("modality") == modality and str(row.get("label", "")) == "1"]
        for k in (1, 2, 3, 5):
            topk.append({"modality": modality, "k": k, "oracle_positive_loci": len(entries),
                         "top_k_recall": sum(int(row.get("rank", 999)) <= k for row in entries) / max(1, len(entries))})
    return summaries, gene_tests, topk, gate
