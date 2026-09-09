from __future__ import annotations

from collections import defaultdict
from pathlib import Path

from .io import (
    canonical_pair, is_callable, normalize_gene, normalize_modality, read_json,
    reject_truth_columns, sha256,
)
from .panels import (
    METHOD_BASELINE,
    METHOD_GUARDED_CC,
    METHOD_MV_FLOOR,
    METHOD_PLURALITY,
    METHOD_RAW_CC,
    PANELS,
    PLURALITY_METHOD_VERSION,
)


def _key(row: dict[str, str]) -> tuple[str, str, str, str]:
    subject = row.get("subject", row.get("sample", "")).strip()
    if not subject:
        raise ValueError("caller row has no subject/sample")
    return (row.get("cohort", ""), subject, normalize_modality(row["modality"]),
            normalize_gene(row["gene"]))


def _caller(row: dict[str, str]) -> str:
    return row.get("caller", row.get("tool", "")).strip()


def build_consensus(caller_rows: list[dict[str, str]]) -> list[dict]:
    """Build the frozen plurality and two-thirds rules over intended-use callers."""
    reject_truth_columns(caller_rows, "caller prediction input")
    grouped: dict[tuple[str, str, str, str], list[dict[str, str]]] = defaultdict(list)
    for row in caller_rows:
        key = _key(row)
        caller = _caller(row)
        if caller not in PANELS[key[2]]:
            continue
        grouped[key].append(row)
    if not grouped:
        raise ValueError("no intended-use caller rows")

    output = []
    for key in sorted(grouped):
        cohort, subject, modality, gene = key
        tool_pair: dict[str, tuple[str, str]] = {}
        pair_tools: dict[tuple[str, str], set[str]] = defaultdict(set)
        observed_callers = set()
        partial_callers = set()
        source_hashes = set()
        for row in grouped[key]:
            caller = _caller(row)
            if caller in observed_callers:
                raise ValueError(f"duplicate caller row for {caller} at {key}")
            observed_callers.add(caller)
            if row.get("source_sha256", "").strip():
                source_hashes.add(f"{caller}:{row['source_sha256'].strip()}")
            if row.get("call_status", "").strip().lower() == "partial":
                partial_callers.add(caller)
            if not is_callable(row):
                continue
            pair = canonical_pair(row["allele1"], row["allele2"], gene)
            tool_pair[caller] = pair
            pair_tools[pair].add(caller)
        complete_tools = len(tool_pair)
        partial_tools = len(partial_callers)
        missing_tools = len(PANELS[modality]) - complete_tools - partial_tools
        if missing_tools < 0:
            raise ValueError(f"caller accounting exceeds frozen panel at {key}")
        population = grouped[key][0].get("superpopulation", grouped[key][0].get("population", ""))
        common_counts = {
            "cohort": cohort,
            "subject": subject,
            "superpopulation": population,
            "modality": modality,
            "gene": gene,
            "complete_tools": complete_tools,
            "partial_tools": partial_tools,
            "missing_tools": missing_tools,
            "callable_tools": complete_tools,
            "partial_callers": ",".join(sorted(partial_callers)),
            "source_hashes": ",".join(sorted(source_hashes)),
            "method_version": PLURALITY_METHOD_VERSION,
        }
        if not tool_pair:
            common = {
                **common_counts,
                "top_support": 0, "support_fraction": "0", "supporting_callers": "",
                "tie_at_top": 0, "tied_pairs": "",
            }
            output.append({**common, "method": METHOD_PLURALITY, "allele1": "", "allele2": "",
                           "call_status": "no_evidence", "decision_reason": "no_callable_intended_caller"})
            output.append({**common, "method": METHOD_BASELINE, "allele1": "", "allele2": "",
                           "call_status": "no_consensus", "decision_reason": "no_callable_intended_caller"})
            continue
        ranked = sorted(pair_tools, key=lambda pair: (-len(pair_tools[pair]), pair))
        top = ranked[0]
        top_support = len(pair_tools[top])
        tied_values = [pair for pair in ranked if len(pair_tools[pair]) == top_support]
        tied = len(tied_values) > 1
        support_fraction = top_support / complete_tools
        common = {
            **common_counts,
            "top_support": top_support, "support_fraction": f"{support_fraction:.12g}",
            "supporting_callers": ",".join(sorted(pair_tools[top])), "tie_at_top": int(tied),
            "tied_pairs": "|".join(f"{a}+{b}" for a, b in tied_values),
        }
        output.append({
            **common, "method": METHOD_PLURALITY, "allele1": top[0], "allele2": top[1],
            "call_status": "callable", "decision_reason": "max_pair_votes_then_lexicographic",
        })
        accepted = top_support * 3 >= complete_tools * 2
        output.append({
            **common, "method": METHOD_BASELINE,
            "allele1": top[0] if accepted else "", "allele2": top[1] if accepted else "",
            "call_status": "callable" if accepted else "no_consensus",
            "decision_reason": "pair_votes_at_least_two_thirds" if accepted else "below_two_thirds",
        })
    return output


def build_guarded_cc(caller_rows: list[dict[str, str]], raw_cc_rows: list[dict[str, str]]) -> list[dict]:
    """Protect two-thirds calls and invoke the frozen raw CC only below the threshold."""
    reject_truth_columns(caller_rows, "guarded CC caller input")
    reject_truth_columns(raw_cc_rows, "guarded CC raw method input")
    consensus_rows = build_consensus(caller_rows)
    baseline = {_key(r): r for r in consensus_rows if r["method"] == METHOD_BASELINE}
    plurality = {_key(r): r for r in consensus_rows if r["method"] == METHOD_PLURALITY}
    raw_cc = {}
    for row in raw_cc_rows:
        if row.get("method", METHOD_RAW_CC) != METHOD_RAW_CC:
            continue
        key = _key(row)
        if key in raw_cc:
            raise ValueError(f"duplicate raw CC row at {key}")
        raw_cc[key] = row
    unknown = set(raw_cc) - set(baseline)
    if unknown:
        raise ValueError(f"raw CC contains {len(unknown)} loci absent from caller universe")

    output = []
    for key in sorted(baseline):
        guard, top = baseline[key], plurality[key]
        cc = raw_cc.get(key)
        protected = guard["call_status"] == "callable"
        cc_callable = bool(cc and is_callable(cc))
        if protected:
            pair = canonical_pair(guard["allele1"], guard["allele2"], key[3])
            reason = "protected_two_thirds_consensus"
        elif cc_callable:
            pair = canonical_pair(cc["allele1"], cc["allele2"], key[3])
            reason = "champion_challenger_resolved_disagreement"
        elif top.get("call_status") == "callable":
            pair = canonical_pair(top["allele1"], top["allele2"], key[3])
            reason = "plurality_fallback_missing_cc"
        else:
            pair = None
            reason = "no_evidence_all_callers_and_raw_cc_missing"
        cc_pair = canonical_pair(cc["allele1"], cc["allele2"], key[3]) if cc_callable else None
        output.append({
            **{k: guard[k] for k in (
                "cohort", "subject", "superpopulation", "modality", "gene", "callable_tools",
                "complete_tools", "partial_tools", "missing_tools", "partial_callers",
                "top_support", "support_fraction", "supporting_callers", "tie_at_top",
                "tied_pairs", "source_hashes", "method_version")},
            "method": METHOD_GUARDED_CC, "allele1": pair[0] if pair else "",
            "allele2": pair[1] if pair else "",
            "call_status": "callable" if pair else "no_evidence", "decision_reason": reason,
            "raw_cc_same_call": int(cc_pair == pair) if cc_pair and pair else 0,
        })
    return consensus_rows + output


def build_mv_floored_cc(caller_rows: list[dict[str, str]], raw_cc_rows: list[dict[str, str]],
                        policy_path: str | Path) -> list[dict]:
    """Apply the frozen MV-floor routing policy without reading truth.

    The historical policy chooses between raw Champion–Challenger and plurality
    within a vote stratum. Missing raw-CC output always falls back to plurality;
    it never creates a synthetic success or an inferred genotype.
    """
    reject_truth_columns(caller_rows, "MV-floor caller input")
    reject_truth_columns(raw_cc_rows, "MV-floor raw method input")
    policy = read_json(policy_path)
    routes = policy.get("policy", {})
    if policy.get("default_route") not in {"mv", "cc"}:
        raise ValueError("MV-floor policy has no valid default_route")

    consensus_rows = build_consensus(caller_rows)
    plurality = {_key(row): row for row in consensus_rows if row["method"] == METHOD_PLURALITY}
    raw = {}
    for row in raw_cc_rows:
        if row.get("method", METHOD_RAW_CC) != METHOD_RAW_CC:
            continue
        key = _key(row)
        if key in raw:
            raise ValueError(f"duplicate raw CC row at {key}")
        raw[key] = row
    unknown = set(raw) - set(plurality)
    if unknown:
        raise ValueError(f"raw CC contains {len(unknown)} loci absent from caller universe")

    output = []
    audit_fields = (
        "cohort", "subject", "superpopulation", "modality", "gene", "callable_tools",
        "complete_tools", "partial_tools", "missing_tools", "partial_callers",
        "top_support", "support_fraction", "supporting_callers", "tie_at_top",
        "tied_pairs", "source_hashes",
    )
    for key in sorted(plurality):
        top = plurality[key]
        complete = int(top["complete_tools"])
        support = int(top["top_support"])
        if complete <= 1:
            stratum = "single_tool"
        elif support == complete:
            stratum = "unanimous"
        elif int(top["tie_at_top"]):
            stratum = "split"
        elif support * 2 > complete:
            stratum = "clear_majority"
        else:
            stratum = "split"
        route = routes.get(key[2], {}).get(stratum, policy["default_route"])
        if route not in {"mv", "cc"}:
            raise ValueError(f"invalid MV-floor route {route!r} for {key[2]}:{stratum}")
        cc = raw.get(key)
        cc_callable = bool(cc and is_callable(cc))
        selected = cc if route == "cc" and cc_callable else top
        callable_selected = is_callable(selected)
        output.append({
            **{field: top[field] for field in audit_fields},
            "method": METHOD_MV_FLOOR,
            "method_version": "mv-floor-frozen-policy-v1",
            "policy_sha256": sha256(policy_path),
            "allele1": selected["allele1"] if callable_selected else "",
            "allele2": selected["allele2"] if callable_selected else "",
            "call_status": "callable" if callable_selected else "no_evidence",
            "vote_stratum": stratum,
            "selected_route": route if route == "mv" or cc_callable else "mv_fallback",
            "decision_reason": (
                f"frozen_policy_{route}_{stratum}" if route == "mv" or cc_callable
                else f"raw_cc_missing_plurality_fallback_{stratum}"
            ),
        })
    return output
