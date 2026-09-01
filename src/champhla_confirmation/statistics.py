from __future__ import annotations

import random
from collections import Counter, defaultdict


def exact_cluster_signflip(subject_deltas) -> float:
    """Two-sided exact sign-flip test using an aggregated integer distribution."""
    values = [int(value) for value in subject_deltas if int(value) != 0]
    if not values:
        return 1.0
    distribution = Counter({0: 1})
    for value in values:
        updated = Counter()
        for total, count in distribution.items():
            updated[total + value] += count
            updated[total - value] += count
        distribution = updated
    observed = abs(sum(values))
    numerator = sum(count for total, count in distribution.items() if abs(total) >= observed)
    return numerator / sum(distribution.values())


def cluster_bootstrap_ci(rows: list[dict], baseline_field: str, selected_field: str,
                         iterations: int = 100000, seed: int = 20260831) -> tuple[float, float]:
    """Subject bootstrap of the fixed-denominator ratio-estimator difference."""
    by_subject: dict[str, list[int]] = defaultdict(list)
    for row in rows:
        by_subject[row["subject"]].append(int(row[selected_field]) - int(row[baseline_field]))
    subjects = sorted(by_subject)
    if not subjects:
        return 0.0, 0.0
    rng = random.Random(seed)
    estimates = []
    for _ in range(iterations):
        numerator = denominator = 0
        for _sample_index in subjects:
            subject = rng.choice(subjects)
            values = by_subject[subject]
            numerator += sum(values)
            denominator += len(values)
        estimates.append(numerator / denominator if denominator else 0.0)
    estimates.sort()
    lo = estimates[int(0.025 * iterations)]
    hi = estimates[min(iterations - 1, int(0.975 * iterations))]
    return lo, hi


def holm_adjust(records: list[dict], p_field: str = "p_value") -> list[dict]:
    """Return records with monotone Holm-adjusted p-values and rejection flags."""
    ordered = sorted(enumerate(records), key=lambda item: float(item[1][p_field]))
    adjusted = [1.0] * len(records)
    running = 0.0
    total = len(records)
    for rank, (index, row) in enumerate(ordered):
        value = min(1.0, (total - rank) * float(row[p_field]))
        running = max(running, value)
        adjusted[index] = running
    result = []
    for index, row in enumerate(records):
        result.append({**row, "holm_adjusted_p": adjusted[index],
                       "holm_significant": int(adjusted[index] < 0.05)})
    return result


def subject_deltas(rows: list[dict], baseline_field: str, selected_field: str) -> dict[str, int]:
    totals: dict[str, int] = defaultdict(int)
    for row in rows:
        totals[row["subject"]] += int(row[selected_field]) - int(row[baseline_field])
    return dict(totals)

