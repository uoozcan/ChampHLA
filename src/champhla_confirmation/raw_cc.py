from __future__ import annotations

import math
from collections import defaultdict

from .io import is_callable, read_json, reject_truth_columns, sha256, write_json
from .panels import PANELS


FROZEN_POLICY = {
    "wgs": {
        "champions": {"A": "OptiType", "B": "OptiType", "C": "T1K"},
        "min_challenger_support_fraction": 0.65,
        "min_challenger_margin": 0.20,
        "min_supporting_tools": 2,
    },
    "wes": {
        "champions": {"A": "OptiType", "B": "OptiType", "C": "OptiType"},
        "min_challenger_support_fraction": 0.35,
        "min_challenger_margin": 0.0,
        "min_supporting_tools": 1,
    },
    "rnaseq": {
        "champions": {"A": "OptiType", "B": "ArcasHLA", "C": "ArcasHLA"},
        "min_challenger_support_fraction": 0.35,
        "min_challenger_margin": 0.0,
        "min_supporting_tools": 2,
    },
}


def _float(value: str | float | int | None) -> float | None:
    try:
        return float(value) if value not in (None, "") else None
    except (TypeError, ValueError):
        return None


def _calibration_stats(rows: list[dict[str, str]]) -> tuple[float | None, float | None]:
    scored = [(p, int(row.get("is_correct", "0") or 0)) for row in rows
              if (p := _float(row.get("calibrated_probability"))) is not None]
    if not scored:
        return None, None
    brier = sum((p - y) ** 2 for p, y in scored) / len(scored)
    bins: dict[int, list[tuple[float, int]]] = defaultdict(list)
    for p, y in scored:
        bins[min(4, int(max(0.0, min(1.0, p)) * 5))].append((p, y))
    ece = sum(abs(sum(p for p, _ in values) / len(values)
                  - sum(y for _, y in values) / len(values)) * len(values)
              for values in bins.values()) / len(scored)
    return brier, ece


def freeze_raw_cc_policy(development_rows: list[dict[str, str]], source_path: str,
                         output_json: str) -> dict:
    """Fit the existing reliability/confidence weight formula on development only."""
    required = {"tool", "modality", "gene", "is_callable", "is_correct"}
    if not development_rows or not required.issubset(development_rows[0]):
        raise ValueError(f"development rows missing policy-fit columns: {sorted(required)}")
    grouped: dict[tuple[str, str, str], list[dict[str, str]]] = defaultdict(list)
    for row in development_rows:
        modality = row.get("modality", "").lower()
        tool = row.get("tool", "")
        gene = row.get("gene", "").replace("HLA-", "")
        if modality in PANELS and tool in PANELS[modality] and gene in {"A", "B", "C"}:
            grouped[(modality, tool, gene)].append(row)
    weights: dict[str, dict[str, dict[str, float]]] = defaultdict(lambda: defaultdict(dict))
    diagnostics = []
    for (modality, tool, gene), rows in sorted(grouped.items()):
        reliability = sum(int(row.get("is_correct", "0") or 0) for row in rows) / len(rows)
        callable_rows = [row for row in rows if row.get("is_callable") == "1"]
        confidence_coverage = (sum(_float(row.get("confidence_score")) is not None
                                   for row in callable_rows) / len(callable_rows)
                               if callable_rows else 0.0)
        probabilities = [_float(row.get("calibrated_probability")) for row in callable_rows]
        probabilities = [value for value in probabilities if value is not None]
        mean_probability = sum(probabilities) / len(probabilities) if probabilities else None
        brier, ece = _calibration_stats(rows)
        if (mean_probability is None or confidence_coverage < 0.5 or brier is None or ece is None
                or max(brier, ece) > 0.35):
            effective = reliability
            guardrail = "reliability_fallback"
        else:
            shrink = max(0.0, 1.0 - max(brier, ece))
            effective = reliability + shrink * (mean_probability - reliability)
            guardrail = "calibrated_confidence_applied"
        weight = round(0.7 * reliability + 0.3 * effective, 4)
        weights[modality][tool][gene] = weight
        diagnostics.append({
            "modality": modality, "caller": tool, "gene": gene, "rows": len(rows),
            "reliability": round(reliability, 6), "confidence_coverage": round(confidence_coverage, 6),
            "mean_calibrated_probability": None if mean_probability is None else round(mean_probability, 6),
            "brier": None if brier is None else round(brier, 6),
            "ece": None if ece is None else round(ece, 6), "guardrail": guardrail,
            "final_weight": weight,
        })
    missing = [(modality, caller, gene) for modality, callers in PANELS.items()
               for caller in callers for gene in ("A", "B", "C")
               if gene not in weights.get(modality, {}).get(caller, {})]
    if missing:
        raise ValueError(f"cannot freeze raw CC policy; missing weights: {missing}")
    payload = {
        "schema_version": "truth-free-raw-cc-policy-1",
        "fit_scope": "development_only",
        "external_truth_used": False,
        "source_sha256": sha256(source_path),
        "weight_formula": "0.7*exact_reliability + 0.3*guarded_effective_calibrated_confidence",
        "confidence_guardrail": {"min_coverage": 0.5, "max_brier_or_ece": 0.35},
        "policies": FROZEN_POLICY,
        "weights": {m: {t: dict(g) for t, g in tools.items()} for m, tools in weights.items()},
        "diagnostics": diagnostics,
    }
    write_json(output_json, payload)
    return payload


def _pair(row: dict[str, str]) -> tuple[str, str] | None:
    return (row["allele1"], row["allele2"]) if is_callable(row) else None


def _rank(entries: list[dict[str, str]], weights: dict[str, dict[str, float]]):
    candidates: dict[tuple[str, str], dict] = defaultdict(lambda: {"weight": 0.0, "callers": []})
    total = 0.0
    for row in entries:
        pair = _pair(row)
        if not pair:
            continue
        weight = float(weights.get(row["caller"], {}).get(row["gene"], 0.0))
        if weight <= 0:
            continue
        candidates[pair]["weight"] += weight
        candidates[pair]["callers"].append(row["caller"])
        total += weight
    ranked = sorted(candidates.items(), key=lambda item: (-item[1]["weight"], item[0]))
    return ranked, total


def predict_raw_cc(calls: list[dict[str, str]], policy_json: str) -> list[dict]:
    reject_truth_columns(calls, "raw CC runtime caller calls")
    bundle = read_json(policy_json)
    if bundle.get("schema_version") != "truth-free-raw-cc-policy-1":
        raise ValueError("unsupported raw CC policy bundle")
    grouped: dict[tuple[str, str, str, str], list[dict[str, str]]] = defaultdict(list)
    for row in calls:
        grouped[(row.get("cohort", ""), row["subject"], row["modality"], row["gene"])].append(row)
    output = []
    for key, entries in sorted(grouped.items()):
        cohort, subject, modality, gene = key
        allowed = set(PANELS[modality])
        entries = [row for row in entries if row["caller"] in allowed]
        policy = bundle["policies"][modality]
        weights = bundle["weights"][modality]
        champion = policy["champions"][gene]
        champion_row = next((row for row in entries if row["caller"] == champion), None)
        champion_pair = _pair(champion_row) if champion_row else None
        ranked, total = _rank(entries, weights)
        chosen = None
        reason = "champion_missing_no_weighted_consensus"
        challenger_pair = None
        support = margin = 0.0
        supporting_callers: list[str] = []
        if champion_pair:
            chosen = champion_pair
            reason = "champion_retained"
            challenger = next(((pair, meta) for pair, meta in ranked if pair != champion_pair), None)
            if challenger and total:
                challenger_pair, meta = challenger
                support = meta["weight"] / total
                if ranked and ranked[0][0] == champion_pair:
                    comparator = ranked[1][1]["weight"] if len(ranked) > 1 else 0.0
                else:
                    comparator = ranked[0][1]["weight"] if ranked else 0.0
                margin = (meta["weight"] - comparator) / total
                supporting_callers = sorted(set(meta["callers"]))
                if (support >= policy["min_challenger_support_fraction"]
                        and margin >= policy["min_challenger_margin"]
                        and len(supporting_callers) >= policy["min_supporting_tools"]
                        and not any("/" in allele for allele in challenger_pair)):
                    chosen = challenger_pair
                    reason = "challenger_override"
        elif ranked and total:
            top_pair, top_meta = ranked[0]
            second = ranked[1][1]["weight"] if len(ranked) > 1 else 0.0
            support = top_meta["weight"] / total
            margin = (top_meta["weight"] - second) / total
            supporting_callers = sorted(set(top_meta["callers"]))
            tied = len(ranked) > 1 and math.isclose(top_meta["weight"], second, abs_tol=1e-12)
            if not tied and support >= 0.55 and margin >= 0.15:
                chosen = top_pair
                reason = "champion_missing_weighted_consensus"
        output.append({
            "cohort": cohort, "subject": subject,
            "superpopulation": entries[0].get("superpopulation", "") if entries else "",
            "modality": modality, "gene": gene, "method": "ChampionChallenger",
            "allele1": chosen[0] if chosen else "", "allele2": chosen[1] if chosen else "",
            "call_status": "callable" if chosen else "missing", "decision_reason": reason,
            "champion_caller": champion,
            "champion_pair": "+".join(champion_pair) if champion_pair else "",
            "challenger_pair": "+".join(challenger_pair) if challenger_pair else "",
            "challenger_support_fraction": round(support, 6),
            "challenger_support_margin": round(margin, 6),
            "supporting_callers": ",".join(supporting_callers),
            "policy_sha256": sha256(policy_json),
        })
    return output
