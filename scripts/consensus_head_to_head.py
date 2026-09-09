#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path

from champhla_confirmation.evaluation import evaluate
from champhla_confirmation.io import read_json, read_tsv, write_json
from champhla_recovery.registry import validate_registry


def _render(primary: list[dict[str, str]], comparisons: list[dict[str, str]],
            endpoints: list[dict[str, str]], method_kinds: dict[str, str],
            integration_count: int) -> str:
    by_modality = {row["modality"]: row for row in primary}
    endpoints_by_modality = {row["modality"]: row for row in endpoints}
    lines = [
        "# Compact plurality head-to-head",
        "",
        f"The current comparator manifest contains {integration_count} integration or historical-ablation methods.",
        "",
        "| Modality | Plurality correct / loci | Plurality call rate | Comparator role | Comparator | Comparator correct | Comparator − plurality | Simultaneous 95% CI | Holm-adjusted p | Status |",
        "|---|---:|---:|---|---|---:|---:|---:|---:|---|",
    ]
    for modality in ("wgs", "wes", "rnaseq"):
        row = by_modality.get(modality)
        if not row:
            lines.append(f"| {modality} | — | — | — | — | — | — | — | — | awaiting corrected rerun |")
            continue
        if row.get("validity") != "valid":
            lines.append(f"| {modality} | — | — | — | — | — | — | — | — | invalid source excluded; awaiting corrected rerun |")
            continue
        candidates = [value for value in comparisons
                      if value["modality"] == modality
                      and int(value.get("comparison_valid", 1))]
        count = f"{row['correct']} / {row['loci']}"
        call_rate = f"{100 * float(row.get('call_rate', 0)):.1f}%"
        endpoint = endpoints_by_modality.get(modality, {})
        best_caller = endpoint.get("best_observed_individual_caller", "")
        integrations = [value for value in candidates
                        if method_kinds.get(value["comparator"]) == "integration"]
        selected = []
        if best_caller:
            match = next((value for value in candidates if value["comparator"] == best_caller), None)
            if match:
                selected.append(("best observed caller (descriptive)", match))
        if integrations:
            strongest = max(integrations, key=lambda value: (int(value["comparator_correct"]),
                                                              value["comparator"]))
            selected.append(("strongest integration", strongest))
        if not selected:
            lines.append(f"| {modality} | {count} | {call_rate} | — | — | — | — | — | — | comparator artifact unavailable |")
            continue
        status = "valid development result"
        for role, selected_row in selected:
            ci = (f"[{float(selected_row['simultaneous_ci_lo_points']):.2f}, "
                  f"{float(selected_row['simultaneous_ci_hi_points']):.2f}]")
            lines.append(
                f"| {modality} | {count} | {call_rate} | {role} | {selected_row['comparator']} | "
                f"{selected_row['comparator_correct']} | "
                f"{float(selected_row['comparator_minus_reference_points']):.2f} | {ci} | "
                f"{float(selected_row['holm_adjusted_p']):.4g} | {status} |"
            )
    return "\n".join(lines) + "\n"


def main() -> int:
    parser = argparse.ArgumentParser(description="Generate the registry-backed plurality main table")
    parser.add_argument("--registry", required=True)
    parser.add_argument("--comparators", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--summary", required=True)
    parser.add_argument("--joined")
    parser.add_argument("--evaluation-dir")
    parser.add_argument("--evaluation-design")
    parser.add_argument("--bootstrap", type=int, default=100000)
    args = parser.parse_args()
    failures = validate_registry(args.registry, str(Path(args.registry).resolve().parent))
    if failures:
        raise ValueError("; ".join(failures))
    if args.joined and not args.evaluation_dir:
        raise ValueError("--joined requires --evaluation-dir")
    if args.joined:
        evaluate(
            args.joined,
            args.evaluation_dir,
            bootstrap=args.bootstrap,
            allow_partial=True,
            mode="discovery",
            evaluation_design=args.evaluation_design,
            comparator_manifest=args.comparators,
        )
        primary = read_tsv(Path(args.evaluation_dir) / "primary_modality_results.tsv")
        comparisons = read_tsv(Path(args.evaluation_dir) / "head_to_head.tsv")
        endpoints = read_tsv(Path(args.evaluation_dir) / "secondary_endpoints.tsv")
        source = "fresh_evaluation"
    elif args.evaluation_dir:
        primary = read_tsv(Path(args.evaluation_dir) / "primary_modality_results.tsv")
        comparisons = read_tsv(Path(args.evaluation_dir) / "head_to_head.tsv")
        endpoints = read_tsv(Path(args.evaluation_dir) / "secondary_endpoints.tsv")
        existing_decision = read_json(Path(args.evaluation_dir) / "confirmation_decision.json")
        source = "existing_evaluation"
    else:
        rows = read_tsv(args.registry)
        primary = [
            {**row, "correct": row["correct"], "loci": row["loci"]}
            for row in rows if row["method"] == "SimplePluralityLex"
        ]
        comparisons = []
        endpoints = []
        source = "registry_only_comparators_pending"
    manifest = read_json(args.comparators)
    integration_count = sum(
        row.get("kind") in {"integration", "historical_ablation"}
        for row in manifest["comparators"]
    )
    method_kinds = {row["method_id"]: row["kind"] for row in manifest["comparators"]}
    target = Path(args.output)
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(
        _render(primary, comparisons, endpoints, method_kinds, integration_count),
        encoding="utf-8",
    )
    write_json(args.summary, {
        "schema_version": "plurality-head-to-head-render-1",
        "source": source,
        "bootstrap_iterations": (args.bootstrap if args.joined else
                                 existing_decision.get("bootstrap_iterations", 0)
                                 if args.evaluation_dir else 0),
        "evaluation_design": (args.evaluation_design or
                              existing_decision.get("evaluation_design", "")
                              if args.evaluation_dir else ""),
        "integration_method_count": integration_count,
        "comparison_rows": len(comparisons),
        "output": target.as_posix(),
    })
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
