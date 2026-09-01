"""Truth-free candidate construction and fold-local runtime context."""

from __future__ import annotations

import csv
import math
from collections import Counter, defaultdict
from pathlib import Path

from . import GENES, MODALITIES
from .alleles import canonical_pair, pair_key


PANELS = {
    "wgs": ("HLA-HD", "Kourami", "OptiType", "SpecHLA", "T1K"),
    "wes": ("HLA-HD", "OptiType", "POLYSOLVER", "SpecHLA", "T1K"),
    "rnaseq": ("ArcasHLA", "HLA-HD", "OptiType", "T1K"),
}
TOOLS = tuple(sorted({tool for panel in PANELS.values() for tool in panel}))
FAMILIES = {
    "HLA-HD": "alignment", "Kourami": "graph", "OptiType": "optimization",
    "SpecHLA": "alignment", "T1K": "kmer", "POLYSOLVER": "alignment",
    "ArcasHLA": "transcript",
}
GLOBAL_FEATURE_NAMES = (
    "equal_support_fraction", "weighted_support_fraction", "vote_margin",
    "vote_entropy", "callable_fraction", "missing_fraction", "homozygous",
    "reference_resolved_fraction", "reference_partial", "reference_alias_or_deleted",
)
TOOL_FEATURE_NAMES = ("supports_candidate", "callable", "calibrated_confidence",
                      "confidence_missing", "reliability")


def read_tsv(path: str | Path) -> list[dict]:
    with Path(path).open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def exclude_subjects(rows: list[dict], subjects: set[str]) -> tuple[list[dict], set[str]]:
    """Remove reserved subjects globally across every modality and gene."""
    observed = {row.get("sample", "") for row in rows}
    overlap = observed & subjects
    return [row for row in rows if row.get("sample") not in subjects], overlap


def write_tsv(path: str | Path, rows: list[dict], fields=None):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    if fields is None:
        fields = []
        for row in rows:
            for key in row:
                if key not in fields:
                    fields.append(key)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def truth_pair(row: dict) -> tuple[str, str]:
    return canonical_pair(row.get("truth_allele1"), row.get("truth_allele2"), row.get("gene"))


def call_pair(row: dict) -> tuple[str, str]:
    return canonical_pair(row.get("allele1"), row.get("allele2"), row.get("gene"))


def is_callable(row: dict) -> bool:
    return str(row.get("is_callable", "")).strip() == "1" and call_pair(row) != ("", "")


def numeric_confidence(row: dict):
    for key in ("cv_calibrated_probability", "confidence_score", "raw_score_value", "confidence"):
        try:
            value = float(row.get(key, ""))
            if math.isfinite(value):
                return value
        except (TypeError, ValueError):
            continue
    return None


def _sigmoid(value: float) -> float:
    value = max(-30.0, min(30.0, value))
    return 1.0 / (1.0 + math.exp(-value))


def _fit_platt(points: list[tuple[float, int]], fallback: float) -> dict:
    if len(points) < 10 or len({label for _, label in points}) < 2:
        return {"slope": 0.0, "intercept": math.log(max(1e-6, fallback) / max(1e-6, 1 - fallback)),
                "mean": 0.0, "scale": 1.0, "n": len(points), "fallback": True}
    mean = sum(value for value, _ in points) / len(points)
    variance = sum((value - mean) ** 2 for value, _ in points) / len(points)
    scale = max(1e-6, math.sqrt(variance))
    slope, intercept = 0.0, math.log(max(1e-6, fallback) / max(1e-6, 1 - fallback))
    for _ in range(400):
        grad_slope = grad_intercept = 0.0
        for value, label in points:
            z = (value - mean) / scale
            error = _sigmoid(intercept + slope * z) - label
            grad_slope += error * z
            grad_intercept += error
        slope -= 0.05 * (grad_slope / len(points) + 0.01 * slope)
        intercept -= 0.05 * grad_intercept / len(points)
    return {"slope": slope, "intercept": intercept, "mean": mean, "scale": scale,
            "n": len(points), "fallback": False}


def fit_runtime_context(rows: list[dict]) -> dict:
    successes, totals, confidence = Counter(), Counter(), defaultdict(list)
    for row in rows:
        modality, gene, tool = row.get("modality", "").lower(), row.get("gene"), row.get("tool")
        if modality not in PANELS or gene not in GENES or tool not in PANELS[modality] or not is_callable(row):
            continue
        key = f"{modality}|{gene}|{tool}"
        correct = int(call_pair(row) == truth_pair(row) and truth_pair(row) != ("", ""))
        totals[key] += 1
        successes[key] += correct
        score = numeric_confidence(row)
        if score is not None:
            confidence[key].append((score, correct))
    reliability, calibration = {}, {}
    for modality in MODALITIES:
        for gene in GENES:
            for tool in PANELS[modality]:
                key = f"{modality}|{gene}|{tool}"
                rate = (successes[key] + 1.0) / (totals[key] + 2.0)
                reliability[key] = rate
                calibration[key] = _fit_platt(confidence[key], rate)
    return {"reliability": reliability, "calibration": calibration,
            "panels": {key: list(value) for key, value in PANELS.items()}}


def calibrated_confidence(context: dict, row: dict) -> tuple[float, float]:
    key = f"{row.get('modality', '').lower()}|{row.get('gene')}|{row.get('tool')}"
    reliability = float(context["reliability"].get(key, 0.5))
    value = numeric_confidence(row)
    if value is None:
        return reliability, 1.0
    model = context["calibration"].get(key, {})
    z = (value - float(model.get("mean", 0.0))) / max(1e-6, float(model.get("scale", 1.0)))
    return _sigmoid(float(model.get("intercept", 0.0)) + float(model.get("slope", 0.0)) * z), 0.0


def entropy(counts: list[int]) -> float:
    total = sum(counts)
    if total <= 1 or len(counts) <= 1:
        return 0.0
    raw = -sum((count / total) * math.log((count / total) + 1e-12) for count in counts if count)
    return raw / math.log(len(counts))


def build_loci(rows: list[dict], context: dict, reference, include_labels: bool = False) -> list[dict]:
    grouped = defaultdict(list)
    for row in rows:
        modality = row.get("modality", "").lower()
        if modality in PANELS and row.get("gene") in GENES and row.get("tool") in PANELS[modality]:
            grouped[(row.get("sample"), modality, row.get("gene"))].append(dict(row))
    loci = []
    for (sample, modality, gene), entries in sorted(grouped.items()):
        by_tool = {row.get("tool"): row for row in entries}
        proposals = []
        for tool in PANELS[modality]:
            row = by_tool.get(tool, {"tool": tool, "modality": modality, "gene": gene})
            pair = call_pair(row) if is_callable(row) else ("", "")
            calibrated, missing = calibrated_confidence(context, row)
            key = f"{modality}|{gene}|{tool}"
            proposals.append({"tool": tool, "family": FAMILIES[tool], "pair": pair,
                              "callable": float(pair != ("", "")), "confidence": calibrated,
                              "confidence_missing": missing,
                              "reliability": float(context["reliability"].get(key, 0.5))})
        candidates = sorted({proposal["pair"] for proposal in proposals if proposal["pair"] != ("", "")})
        if not candidates:
            continue
        counts = Counter(proposal["pair"] for proposal in proposals if proposal["pair"] != ("", ""))
        ordered_counts = sorted(counts.values(), reverse=True)
        top, second = ordered_counts[0], ordered_counts[1] if len(ordered_counts) > 1 else 0
        callable_count = sum(proposal["callable"] for proposal in proposals)
        truth = truth_pair(entries[0]) if include_labels else ("", "")
        candidate_rows = []
        for candidate in candidates:
            weighted_total = sum(proposal["reliability"] for proposal in proposals if proposal["callable"])
            weighted_support = sum(proposal["reliability"] for proposal in proposals if proposal["pair"] == candidate)
            flags = [reference.resolve(allele)[1] for allele in candidate]
            global_features = [
                counts[candidate] / max(1.0, callable_count),
                weighted_support / max(1e-9, weighted_total),
                (top - second) / max(1.0, callable_count),
                entropy(list(counts.values())),
                callable_count / len(PANELS[modality]),
                1.0 - callable_count / len(PANELS[modality]),
                float(candidate[0] == candidate[1]),
                sum(not item["unresolved"] for item in flags) / 2.0,
                float(any(item["partial"] for item in flags)),
                float(any(item["alias_resolved"] or item["deleted"] for item in flags)),
            ]
            tool_tokens = []
            for proposal in proposals:
                tool_tokens.append({**proposal, "features": [
                    float(proposal["pair"] == candidate), proposal["callable"], proposal["confidence"],
                    proposal["confidence_missing"], proposal["reliability"],
                ]})
            record = {"pair": candidate, "pair_key": pair_key(candidate), "global_features": global_features,
                      "tool_tokens": tool_tokens, "reference_flags": flags}
            if include_labels:
                record["label"] = int(candidate == truth)
            candidate_rows.append(record)
        locus = {"sample": sample, "modality": modality, "gene": gene, "candidates": candidate_rows,
                 "candidate_count": len(candidate_rows), "truth_pair": truth if include_labels else ("", ""),
                 "candidate_set_oracle": int(truth in candidates) if include_labels else None}
        loci.append(locus)
    return loci
