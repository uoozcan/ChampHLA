#!/usr/bin/env python3
"""Nested cross-validated evaluation of the Champion-Challenger (CC) consensus.

Motivation
----------
The headline WES/RNA CC accuracy in the manuscript (e.g. WES 0.9487) is the
*maximum* over a 49-point override-policy grid swept on the SAME cohort it is
scored on (``wes_champion_override_sweep.tsv``), with no held-out split. That is
in-sample operating-point selection and inflates the number. This script fixes
that: it selects the CC champions and override policy strictly OUT OF FOLD and
reports CC accuracy only on held-out samples, pooled across folds. Majority vote
(MV) needs no tuning and is scored on the identical held-out loci, so the paired
McNemar test is a fair CC-vs-MV comparison.

It reuses the pipeline's own scoring functions from ``hla_benchmark`` (the CC and
MV builders), so the re-analysis cannot silently diverge from the production
algorithm — it only changes *which loci the operating point is chosen on*.

Inputs are the already-produced harmonized rows + runtime weights from an
existing benchmark run (no HLA tools are re-run).
"""

import argparse
import json
import math
import random
import sys
from collections import defaultdict
from pathlib import Path

BIN_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(BIN_DIR))
import hla_benchmark as hb  # noqa: E402
import ciwd as ciwd_mod  # noqa: E402

# 1000G superpopulation → continental ancestry (for the per-ancestry robustness check).
ANCESTRY = {"CEU": "EUR", "FIN": "EUR", "GBR": "EUR", "TSI": "EUR", "IBS": "EUR",
            "YRI": "AFR", "LWK": "AFR", "CHB": "EAS", "JPT": "EAS"}

# The manuscript override-policy grid (Methods §2).
SUPPORT_GRID = [0.20, 0.35, 0.50, 0.65]
MARGIN_GRID = [0.0, 0.05, 0.10, 0.20]
TOOLS_GRID = [1, 2, 3]


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--harmonized", required=True, help="harmonized_benchmark_rows.tsv from an existing run")
    p.add_argument("--weights", required=True, help="consensus_runtime_weights.json from the same run")
    p.add_argument("--modality", required=True, help="modality label as it appears in the rows (e.g. wes, rnaseq, wgs)")
    p.add_argument("--genes", default="A,B,C")
    p.add_argument("--folds", type=int, default=10)
    p.add_argument("--seed", type=int, default=42)
    p.add_argument("--min-support", type=float, default=0.55, help="weighted-consensus fallback support threshold")
    p.add_argument("--min-margin", type=float, default=0.15, help="weighted-consensus fallback margin threshold")
    p.add_argument("--out", required=True, help="output directory for nested-CV tables")
    return p.parse_args()


def load_rows(path, genes, modality):
    rows = hb.read_table(Path(path))
    keep = []
    gene_set = set(genes)
    for r in rows:
        g = hb.normalize_gene(r.get("gene", ""))
        if g not in gene_set:
            continue
        if hb.clean_token(r.get("modality", "")).lower() != modality.lower():
            continue
        # only truth-backed loci are scorable
        if not hb.clean_token(r.get("truth_allele1", "")):
            continue
        r["gene"] = g
        keep.append(r)
    return keep


def stratified_folds(samples_by_pop, n_folds, seed):
    """Population-stratified k-fold assignment: within each superpopulation,
    shuffle deterministically and round-robin into folds so every fold keeps
    the cohort's ancestry balance."""
    rng = random.Random(seed)
    fold_of = {}
    for pop in sorted(samples_by_pop):
        members = sorted(samples_by_pop[pop])
        rng.shuffle(members)
        for i, s in enumerate(members):
            fold_of[s] = i % n_folds
    return fold_of


def learn_champions(train_rows, genes):
    """Champion per gene = single tool with the highest 2-field accuracy on the
    training fold (ties broken by callable count then name)."""
    stats = defaultdict(lambda: defaultdict(lambda: [0, 0, 0]))  # gene -> tool -> [correct, callable, total]
    for r in train_rows:
        g = r["gene"]
        tool = hb.clean_token(r.get("tool", ""))
        if not tool:
            continue
        s = stats[g][tool]
        s[2] += 1
        if r.get("is_callable") == "1":
            s[1] += 1
        if r.get("is_correct") == "1":
            s[0] += 1
    champions = {}
    for g in genes:
        tools = stats.get(g, {})
        if not tools:
            continue
        best = max(
            tools.items(),
            key=lambda kv: (kv[1][0] / kv[1][2] if kv[1][2] else 0.0, kv[1][1], kv[0]),
        )
        champions[g] = best[0]
    return champions


def make_config(champions, policy, min_support, min_margin, mode):
    return {
        "benchmark": {
            "mode": mode,
            "champion_challenger": {
                "enabled": True,
                "champion_by_gene": dict(champions),
                "fallback_method": "weighted_consensus",
                "override_policy": {
                    "min_challenger_support_fraction": policy[0],
                    "min_challenger_margin": policy[1],
                    "min_supporting_tools": policy[2],
                    "require_non_ambiguity_override": True,
                },
            },
            "consensus": {"min_support": min_support, "min_margin": min_margin},
        }
    }


def accuracy(call_rows):
    n = len(call_rows)
    correct = sum(1 for r in call_rows if r.get("is_correct") == "1")
    callable_n = sum(1 for r in call_rows if r.get("is_callable") == "1")
    return {
        "n": n,
        "correct": correct,
        "overall": correct / n if n else 0.0,
        "callable_rate": callable_n / n if n else 0.0,
        "acc_among_callable": correct / callable_n if callable_n else 0.0,
    }


def cc_calls(rows, weights, config, genes, mode):
    call_rows, *_ = hb.build_champion_challenger_outputs(rows, weights, config, mode, genes)
    return call_rows


def choose_policy(train_rows, weights, champions, genes, min_support, min_margin, mode):
    """Argmax over the grid on the TRAINING fold. Tie-breaks favour a more
    conservative gate (higher support, higher margin, fewer overrides)."""
    best = None
    for support in SUPPORT_GRID:
        for margin in MARGIN_GRID:
            for tools in TOOLS_GRID:
                policy = (support, margin, tools)
                cfg = make_config(champions, policy, min_support, min_margin, mode)
                calls = cc_calls(train_rows, weights, cfg, genes, mode)
                overrides = sum(1 for r in calls if r.get("override_triggered") == "1")
                acc = accuracy(calls)["overall"]
                key = (acc, support, margin, -overrides)
                if best is None or key > best[0]:
                    best = (key, policy)
    return best[1]


def wilson(k, n, z=1.96):
    if n == 0:
        return (0.0, 0.0)
    lo, hi = hb.wilson_ci(k, n, z)
    return (round(lo, 4), round(hi, 4))


def mcnemar_exact(pairs_a, pairs_b):
    """Exact two-sided binomial McNemar over paired boolean-correct dicts keyed
    by (sample, gene). b = A right/B wrong; c = A wrong/B right."""
    b = c = 0
    keys = set(pairs_a) & set(pairs_b)
    for k in keys:
        a_ok, b_ok = pairs_a[k], pairs_b[k]
        if a_ok and not b_ok:
            b += 1
        elif b_ok and not a_ok:
            c += 1
    n = b + c
    if n == 0:
        return {"b": b, "c": c, "n_discordant": 0, "p_value": 1.0}
    lo = min(b, c)
    # two-sided exact binomial p at p=0.5
    tail = sum(math.comb(n, i) for i in range(0, lo + 1)) / (2 ** n)
    p = min(1.0, 2 * tail)
    return {"b": b, "c": c, "n_discordant": n, "p_value": round(p, 5)}


def correct_map(call_rows):
    return {(r["sample"], r["gene"]): (r.get("is_correct") == "1") for r in call_rows}


def main():
    args = parse_args()
    genes = [hb.normalize_gene(g) for g in args.genes.split(",") if g.strip()]
    mode = "probabilistic_recalibrated"
    rows = load_rows(args.harmonized, genes, args.modality)
    weights = json.loads(Path(args.weights).read_text(encoding="utf-8"))
    if not rows:
        raise SystemExit("No scorable truth-backed rows for modality=%s genes=%s" % (args.modality, genes))

    samples_by_pop = defaultdict(set)
    for r in rows:
        samples_by_pop[hb.clean_token(r.get("superpopulation", "")) or "unknown"].add(r["sample"])
    all_samples = sorted({r["sample"] for r in rows})
    n_folds = min(args.folds, len(all_samples))
    fold_of = stratified_folds(samples_by_pop, n_folds, args.seed)

    by_sample = defaultdict(list)
    for r in rows:
        by_sample[r["sample"]].append(r)

    pooled_cc, pooled_mv = [], []
    fold_records = []
    for fold in range(n_folds):
        test_samples = [s for s in all_samples if fold_of.get(s, 0) == fold]
        train_samples = [s for s in all_samples if fold_of.get(s, 0) != fold]
        if not test_samples or not train_samples:
            continue
        train_rows = [r for s in train_samples for r in by_sample[s]]
        test_rows = [r for s in test_samples for r in by_sample[s]]

        champions = learn_champions(train_rows, genes)
        policy = choose_policy(train_rows, weights, champions, genes, args.min_support, args.min_margin, mode)

        cfg = make_config(champions, policy, args.min_support, args.min_margin, mode)
        cc = cc_calls(test_rows, weights, cfg, genes, mode)
        mv = hb.build_majority_vote_rows(test_rows)
        pooled_cc.extend(cc)
        pooled_mv.extend(mv)
        overrides = sum(1 for r in cc if r.get("override_triggered") == "1")
        fold_records.append({
            "fold": fold,
            "n_test_samples": len(test_samples),
            "champion_A": champions.get("A", ""),
            "champion_B": champions.get("B", ""),
            "champion_C": champions.get("C", ""),
            "min_challenger_support_fraction": policy[0],
            "min_challenger_margin": policy[1],
            "min_supporting_tools": policy[2],
            "holdout_overrides": overrides,
            "holdout_cc_overall": round(accuracy(cc)["overall"], 4),
            "holdout_mv_overall": round(accuracy(mv)["overall"], 4),
        })

    # ---- pooled held-out metrics ----
    cc_m = accuracy(pooled_cc)
    mv_m = accuracy(pooled_mv)

    # best single tool over the full cohort (untuned reference)
    tool_stats = defaultdict(lambda: [0, 0])  # tool -> [correct, n]
    for r in rows:
        t = hb.clean_token(r.get("tool", ""))
        tool_stats[t][1] += 1
        if r.get("is_correct") == "1":
            tool_stats[t][0] += 1
    best_tool, best_tool_acc = "", 0.0
    for t, (k, n) in tool_stats.items():
        acc = k / n if n else 0.0
        if acc > best_tool_acc:
            best_tool, best_tool_acc = t, acc

    mcnemar_cc_mv = mcnemar_exact(correct_map(pooled_cc), correct_map(pooled_mv))

    out = Path(args.out)
    (out).mkdir(parents=True, exist_ok=True)

    comp_rows = []
    for name, m in [("MajorityVote", mv_m), ("ChampionChallenger_nestedCV", cc_m)]:
        lo, hi = wilson(m["correct"], m["n"])
        comp_rows.append({
            "method": name,
            "modality": args.modality,
            "n_holdout_gene_rows": m["n"],
            "callable_rate": round(m["callable_rate"], 4),
            "accuracy_among_callable": round(m["acc_among_callable"], 4),
            "overall_correct_call_rate": round(m["overall"], 4),
            "ci_lo": lo,
            "ci_hi": hi,
        })
    comp_rows.append({
        "method": "BestSingleTool(%s)" % best_tool,
        "modality": args.modality,
        "n_holdout_gene_rows": tool_stats[best_tool][1],
        "callable_rate": "",
        "accuracy_among_callable": "",
        "overall_correct_call_rate": round(best_tool_acc, 4),
        "ci_lo": wilson(tool_stats[best_tool][0], tool_stats[best_tool][1])[0],
        "ci_hi": wilson(tool_stats[best_tool][0], tool_stats[best_tool][1])[1],
    })
    hb.write_tsv(out / "nested_cv_method_comparison.tsv", comp_rows,
                 ["method", "modality", "n_holdout_gene_rows", "callable_rate",
                  "accuracy_among_callable", "overall_correct_call_rate", "ci_lo", "ci_hi"])
    hb.write_tsv(out / "nested_cv_fold_policies.tsv", fold_records,
                 list(fold_records[0].keys()) if fold_records else ["fold"])
    mcnemar_row = dict(mcnemar_cc_mv)
    mcnemar_row.update({"comparison": "ChampionChallenger_vs_MajorityVote", "modality": args.modality,
                        "cc_overall": round(cc_m["overall"], 4), "mv_overall": round(mv_m["overall"], 4),
                        "delta_cc_minus_mv": round(cc_m["overall"] - mv_m["overall"], 4)})
    hb.write_tsv(out / "nested_cv_mcnemar.tsv", [mcnemar_row],
                 ["comparison", "modality", "cc_overall", "mv_overall", "delta_cc_minus_mv",
                  "b", "c", "n_discordant", "p_value"])

    # ---- per-gene held-out metrics (for Table 4 and Figure 3) ----
    per_gene_rows, per_gene_mcnemar, per_gene_gain = [], [], []
    for g in genes:
        cc_g = [r for r in pooled_cc if r.get("gene") == g]
        mv_g = [r for r in pooled_mv if r.get("gene") == g]
        if not cc_g or not mv_g:
            continue
        cc_gm, mv_gm = accuracy(cc_g), accuracy(mv_g)
        # best single tool for this gene (untuned reference)
        g_tool = defaultdict(lambda: [0, 0])
        for r in rows:
            if r.get("gene") != g:
                continue
            t = hb.clean_token(r.get("tool", ""))
            g_tool[t][1] += 1
            if r.get("is_correct") == "1":
                g_tool[t][0] += 1
        bt, bt_acc, bt_k, bt_n = "", 0.0, 0, 0
        for t, (k, n) in g_tool.items():
            a = k / n if n else 0.0
            if a > bt_acc:
                bt, bt_acc, bt_k, bt_n = t, a, k, n
        for name, m in [("MajorityVote", mv_gm), ("ChampionChallenger_nestedCV", cc_gm)]:
            lo, hi = wilson(m["correct"], m["n"])
            per_gene_rows.append({
                "gene": g, "method": name, "modality": args.modality,
                "n_gene_rows": m["n"], "callable_rate": round(m["callable_rate"], 4),
                "accuracy_among_callable": round(m["acc_among_callable"], 4),
                "overall_correct_call_rate": round(m["overall"], 4), "ci_lo": lo, "ci_hi": hi,
            })
        lo, hi = wilson(bt_k, bt_n)
        per_gene_rows.append({
            "gene": g, "method": "BestSingleTool(%s)" % bt, "modality": args.modality,
            "n_gene_rows": bt_n, "callable_rate": "", "accuracy_among_callable": "",
            "overall_correct_call_rate": round(bt_acc, 4), "ci_lo": lo, "ci_hi": hi,
        })
        mc_g = mcnemar_exact(correct_map(cc_g), correct_map(mv_g))
        mc_g.update({"gene": g, "modality": args.modality,
                     "cc_overall": round(cc_gm["overall"], 4), "mv_overall": round(mv_gm["overall"], 4),
                     "delta_cc_minus_mv": round(cc_gm["overall"] - mv_gm["overall"], 4),
                     "best_tool": bt, "best_tool_overall": round(bt_acc, 4)})
        per_gene_mcnemar.append(mc_g)
        per_gene_gain.append({"gene": g, "modality": args.modality,
                              "cc_gain_vs_mv": round(cc_gm["overall"] - mv_gm["overall"], 4),
                              "besttool_gain_vs_mv": round(bt_acc - mv_gm["overall"], 4)})
    hb.write_tsv(out / "nested_cv_per_gene.tsv", per_gene_rows,
                 ["gene", "method", "modality", "n_gene_rows", "callable_rate",
                  "accuracy_among_callable", "overall_correct_call_rate", "ci_lo", "ci_hi"])
    hb.write_tsv(out / "nested_cv_per_gene_mcnemar.tsv", per_gene_mcnemar,
                 ["gene", "modality", "cc_overall", "mv_overall", "delta_cc_minus_mv",
                  "best_tool", "best_tool_overall", "b", "c", "n_discordant", "p_value"])
    hb.write_tsv(out / "nested_cv_per_gene_gain.tsv", per_gene_gain,
                 ["gene", "modality", "cc_gain_vs_mv", "besttool_gain_vs_mv"])

    # ---- stratified held-out metrics: CIWD commonness + continental ancestry ----
    try:
        cat = ciwd_mod.load_ciwd()
    except Exception:
        cat = None

    def _strat_buckets(call_rows, stratifier):
        buckets = defaultdict(lambda: [0, 0])  # stratum -> [correct, n]
        for r in call_rows:
            b = buckets[stratifier(r)]
            b[1] += 1
            if r.get("is_correct") == "1":
                b[0] += 1
        return buckets

    def _ciwd_strat(r):
        if cat is None:
            return "unknown"
        return cat.genotype_stratum(r.get("truth_allele1", ""), r.get("truth_allele2", ""))

    def _anc_strat(r):
        return ANCESTRY.get(hb.clean_token(r.get("superpopulation", "")), "unknown")

    for stratifier, fname, colname in [
        (_ciwd_strat, "nested_cv_by_ciwd.tsv", "ciwd_stratum"),
        (_anc_strat, "nested_cv_by_ancestry.tsv", "ancestry"),
    ]:
        cc_b, mv_b = _strat_buckets(pooled_cc, stratifier), _strat_buckets(pooled_mv, stratifier)
        strat_rows = []
        for stratum in sorted(set(cc_b) | set(mv_b)):
            for method, b in [("MajorityVote", mv_b), ("ChampionChallenger_nestedCV", cc_b)]:
                k, n = b[stratum]
                lo, hi = wilson(k, n)
                strat_rows.append({colname: stratum, "method": method, "modality": args.modality,
                                   "n": n, "overall_correct_call_rate": round(k / n, 4) if n else "",
                                   "ci_lo": lo, "ci_hi": hi})
        hb.write_tsv(out / fname, strat_rows,
                     [colname, "method", "modality", "n", "overall_correct_call_rate", "ci_lo", "ci_hi"])

    print("== Nested-CV held-out results (%s, %d folds, n=%d samples) ==" % (args.modality, n_folds, len(all_samples)))
    print("  MajorityVote           overall = %.4f  [%.3f, %.3f]" % (mv_m["overall"], *wilson(mv_m["correct"], mv_m["n"])))
    print("  ChampionChallenger(CV) overall = %.4f  [%.3f, %.3f]" % (cc_m["overall"], *wilson(cc_m["correct"], cc_m["n"])))
    print("  BestSingleTool(%s)     overall = %.4f" % (best_tool, best_tool_acc))
    print("  delta (CC - MV) = %+.4f | McNemar b=%d c=%d p=%.4f" % (
        cc_m["overall"] - mv_m["overall"], mcnemar_cc_mv["b"], mcnemar_cc_mv["c"], mcnemar_cc_mv["p_value"]))
    print("  wrote:", out / "nested_cv_method_comparison.tsv")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
