from __future__ import annotations

from collections import defaultdict

from .io import canonical_pair, is_callable, normalize_gene, normalize_modality, reject_truth_columns
from .panels import METHOD_BASELINE, METHOD_GUARDED_CC, METHOD_PLURALITY, METHOD_RAW_CC, PANELS


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
        for row in grouped[key]:
            caller = _caller(row)
            if not is_callable(row):
                continue
            pair = canonical_pair(row["allele1"], row["allele2"], gene)
            if caller in tool_pair and tool_pair[caller] != pair:
                raise ValueError(f"caller {caller} proposes multiple pairs at {key}")
            tool_pair[caller] = pair
            pair_tools[pair].add(caller)
        if not tool_pair:
            population = grouped[key][0].get("superpopulation", grouped[key][0].get("population", ""))
            common = {
                "cohort": cohort, "subject": subject, "superpopulation": population,
                "modality": modality, "gene": gene, "callable_tools": 0,
                "top_support": 0, "support_fraction": "0", "supporting_callers": "",
                "tie_at_top": 0,
            }
            output.append({**common, "method": METHOD_PLURALITY, "allele1": "", "allele2": "",
                           "call_status": "no_evidence", "decision_reason": "no_callable_intended_caller"})
            output.append({**common, "method": METHOD_BASELINE, "allele1": "", "allele2": "",
                           "call_status": "no_consensus", "decision_reason": "no_callable_intended_caller"})
            continue
        ranked = sorted(pair_tools, key=lambda pair: (-len(pair_tools[pair]), pair))
        top = ranked[0]
        top_support = len(pair_tools[top])
        callable_tools = len(tool_pair)
        tied = sum(len(pair_tools[pair]) == top_support for pair in ranked) > 1
        support_fraction = top_support / callable_tools
        population = grouped[key][0].get("superpopulation", grouped[key][0].get("population", ""))
        common = {
            "cohort": cohort, "subject": subject, "superpopulation": population,
            "modality": modality, "gene": gene, "callable_tools": callable_tools,
            "top_support": top_support, "support_fraction": f"{support_fraction:.12g}",
            "supporting_callers": ",".join(sorted(pair_tools[top])), "tie_at_top": int(tied),
        }
        output.append({
            **common, "method": METHOD_PLURALITY, "allele1": top[0], "allele2": top[1],
            "call_status": "callable", "decision_reason": "max_pair_votes_then_lexicographic",
        })
        accepted = top_support * 3 >= callable_tools * 2
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
                "top_support", "support_fraction", "supporting_callers", "tie_at_top")},
            "method": METHOD_GUARDED_CC, "allele1": pair[0] if pair else "",
            "allele2": pair[1] if pair else "",
            "call_status": "callable" if pair else "no_evidence", "decision_reason": reason,
            "raw_cc_same_call": int(cc_pair == pair) if cc_pair and pair else 0,
        })
    return consensus_rows + output
