#!/usr/bin/env python3
"""
Compute trimodal (WGS + WES + RNA-seq) joint consensus and compare against
per-modality and bimodal baselines.

Outputs:
  tables/trimodal_wgs_wes_rna_consensus.tsv
  tables/trimodal_accuracy_comparison.tsv

Usage:
  python3 compute_trimodal_consensus.py --tables-dir <dir>
"""
import argparse, json, math, sys
from pathlib import Path
from collections import defaultdict

# ── Wilson CI ────────────────────────────────────────────────────────────────
def wilson_ci(n_correct, n_total, z=1.96):
    if n_total == 0:
        return (0.0, 0.0)
    p = n_correct / n_total
    denom = 1 + z**2 / n_total
    centre = (p + z**2 / (2 * n_total)) / denom
    half = z * math.sqrt(p * (1 - p) / n_total + z**2 / (4 * n_total**2)) / denom
    return (max(0.0, centre - half), min(1.0, centre + half))

# ── CSV/TSV helpers ───────────────────────────────────────────────────────────
def read_tsv(path):
    path = Path(path)
    if not path.exists():
        return []
    with open(path) as f:
        lines = f.read().splitlines()
    if not lines:
        return []
    header = lines[0].split("\t")
    rows = []
    for line in lines[1:]:
        if not line.strip():
            continue
        parts = line.split("\t")
        rows.append(dict(zip(header, parts)))
    return rows

def write_tsv(path, rows, fieldnames):
    path = Path(path)
    with open(path, "w") as f:
        f.write("\t".join(fieldnames) + "\n")
        for row in rows:
            f.write("\t".join(str(row.get(k, "")) for k in fieldnames) + "\n")
    print(f"  Wrote {len(rows)} rows → {path}")

# ── Majority vote ─────────────────────────────────────────────────────────────
def majority_vote(calls):
    """calls: list of (allele1, allele2) tuples. Returns best pair or None."""
    tally = defaultdict(int)
    for a1, a2 in calls:
        key = tuple(sorted([a1, a2]))
        tally[key] += 1
    if not tally:
        return None
    return max(tally, key=tally.get)

# ── Weighted consensus ────────────────────────────────────────────────────────
MIN_SUPPORT = 0.55

def weighted_consensus(calls_weights):
    """calls_weights: list of (allele1, allele2, weight). Returns best pair or None."""
    tally = defaultdict(float)
    total_weight = 0.0
    for a1, a2, w in calls_weights:
        key = tuple(sorted([a1, a2]))
        tally[key] += w
        total_weight += w
    if not tally or total_weight == 0:
        return None, 0.0
    best_pair = max(tally, key=tally.get)
    support_frac = tally[best_pair] / total_weight
    if support_frac < MIN_SUPPORT:
        return None, support_frac  # abstain
    return best_pair, support_frac

# ── Accuracy summary ──────────────────────────────────────────────────────────
def accuracy_summary(rows, label, modality_tag, sample_ids):
    callable_rows = [r for r in rows if r.get("is_callable") == "1"]
    correct_rows  = [r for r in callable_rows if r.get("is_correct") == "1"]
    n_total = len(rows)
    n_call  = len(callable_rows)
    n_corr  = len(correct_rows)
    callable_rate = n_call / n_total if n_total else 0
    acc_callable  = n_corr / n_call  if n_call  else 0
    overall       = n_corr / n_total if n_total else 0
    ci_lo, ci_hi  = wilson_ci(n_corr, n_total)
    return {
        "comparison": label,
        "modality": modality_tag,
        "sample_count": len(sample_ids),
        "gene_rows": n_total,
        "callable_rate": round(callable_rate, 4),
        "accuracy_among_callable": round(acc_callable, 4),
        "overall_correct_call_rate": round(overall, 4),
        "overall_correct_call_rate_ci_lo": round(ci_lo, 4),
        "overall_correct_call_rate_ci_hi": round(ci_hi, 4),
    }

# ── Main ──────────────────────────────────────────────────────────────────────
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tables-dir", required=True)
    args = ap.parse_args()
    tables = Path(args.tables_dir)

    # Load harmonized rows (all tool calls)
    print("Loading harmonized_benchmark_rows.tsv …")
    bench_rows = read_tsv(tables / "harmonized_benchmark_rows.tsv")
    print(f"  {len(bench_rows)} rows")

    # Load runtime weights
    print("Loading consensus_runtime_weights.json …")
    with open(tables / "consensus_runtime_weights.json") as f:
        weight_payload = json.load(f)
    tool_weights = weight_payload.get("tool_weights", {})

    # Build callable index: {(sample, gene, modality)} → bool
    # Group callable calls by (sample, gene) → {modality: [(tool, a1, a2, weight)]}
    callable_by_key = defaultdict(lambda: defaultdict(list))

    for row in bench_rows:
        if row.get("is_callable") != "1":
            continue
        sample = row["sample"]
        gene   = row["gene"]
        mod    = row["modality"]
        tool   = row["tool"]
        a1     = row.get("allele1", "")
        a2     = row.get("allele2", "")
        if not a1 or not a2:
            continue
        # Get final_weight for this tool×modality
        tw = tool_weights.get(tool, {}).get(mod, {})
        weight = float(tw.get("final_weight", 0.5))
        callable_by_key[(sample, gene)][mod].append((tool, a1, a2, weight))

    # Trimodal eligibility: need ≥1 callable call in each of wgs, wes, rnaseq
    trimodal_keys = {
        k for k, mods in callable_by_key.items()
        if "wgs" in mods and "wes" in mods and "rnaseq" in mods
    }
    trimodal_samples = {k[0] for k in trimodal_keys}
    trimodal_genes   = sorted({k[1] for k in trimodal_keys})
    print(f"  Trimodal eligible: {len(trimodal_samples)} samples, "
          f"{len(trimodal_keys)} (sample,gene) pairs, genes: {trimodal_genes}")

    # Build truth lookup from bench_rows
    truth_lookup = {}
    for row in bench_rows:
        k = (row["sample"], row["gene"])
        if k not in truth_lookup:
            truth_lookup[k] = (row.get("truth_allele1", ""), row.get("truth_allele2", ""))

    # Generate trimodal consensus rows
    trimodal_rows = []
    for (sample, gene) in sorted(trimodal_keys):
        truth_a1, truth_a2 = truth_lookup.get((sample, gene), ("", ""))
        truth_pair = tuple(sorted([truth_a1, truth_a2]))

        # Pool all calls from all 3 modalities
        all_calls = []
        all_calls_w = []
        for mod in ("wgs", "wes", "rnaseq"):
            for (tool, a1, a2, w) in callable_by_key[(sample, gene)].get(mod, []):
                all_calls.append((a1, a2))
                all_calls_w.append((a1, a2, w))

        # ── TrimodalMajorityVote ──────────────────────────────────────────
        mv_pair = majority_vote(all_calls)
        if mv_pair:
            call_a1, call_a2 = mv_pair
            is_correct = int(tuple(sorted([call_a1, call_a2])) == truth_pair)
            trimodal_rows.append({
                "sample": sample, "modality": "wgs+wes+rnaseq", "gene": gene,
                "method": "TrimodalMajorityVote",
                "truth_allele1": truth_a1, "truth_allele2": truth_a2,
                "allele1": call_a1, "allele2": call_a2,
                "call_status": "called", "is_callable": "1", "is_correct": str(is_correct),
                "contributing_tools": str(len(all_calls)),
                "support_fraction": "",
            })
        else:
            trimodal_rows.append({
                "sample": sample, "modality": "wgs+wes+rnaseq", "gene": gene,
                "method": "TrimodalMajorityVote",
                "truth_allele1": truth_a1, "truth_allele2": truth_a2,
                "allele1": "", "allele2": "",
                "call_status": "no_call", "is_callable": "0", "is_correct": "0",
                "contributing_tools": "0", "support_fraction": "",
            })

        # ── TrimodalWeightedConsensus ────────────────────────────────────
        wc_pair, support = weighted_consensus(all_calls_w)
        if wc_pair:
            call_a1, call_a2 = wc_pair
            is_correct = int(tuple(sorted([call_a1, call_a2])) == truth_pair)
            trimodal_rows.append({
                "sample": sample, "modality": "wgs+wes+rnaseq", "gene": gene,
                "method": "TrimodalWeightedConsensus",
                "truth_allele1": truth_a1, "truth_allele2": truth_a2,
                "allele1": call_a1, "allele2": call_a2,
                "call_status": "called", "is_callable": "1", "is_correct": str(is_correct),
                "contributing_tools": str(len(all_calls_w)),
                "support_fraction": str(round(support, 4)),
            })
        else:
            trimodal_rows.append({
                "sample": sample, "modality": "wgs+wes+rnaseq", "gene": gene,
                "method": "TrimodalWeightedConsensus",
                "truth_allele1": truth_a1, "truth_allele2": truth_a2,
                "allele1": "", "allele2": "",
                "call_status": "low_confidence", "is_callable": "0", "is_correct": "0",
                "contributing_tools": str(len(all_calls_w)),
                "support_fraction": str(round(support, 4)) if support > 0 else "",
            })

    trimodal_fields = ["sample", "modality", "gene", "method",
                       "truth_allele1", "truth_allele2", "allele1", "allele2",
                       "call_status", "is_callable", "is_correct",
                       "contributing_tools", "support_fraction"]
    write_tsv(tables / "trimodal_wgs_wes_rna_consensus.tsv", trimodal_rows, trimodal_fields)

    # ── Accuracy comparison ───────────────────────────────────────────────────
    # Trimodal MV and WC rows
    tri_mv = [r for r in trimodal_rows if r["method"] == "TrimodalMajorityVote"]
    tri_wc = [r for r in trimodal_rows if r["method"] == "TrimodalWeightedConsensus"]

    # Per-modality baselines restricted to trimodal-eligible samples
    wgs_rows = [r for r in bench_rows
                if r["sample"] in trimodal_samples and r.get("modality") == "wgs"
                and r.get("method") in ("MajorityVote", "WeightedConsensus")]
    wes_rows = [r for r in bench_rows
                if r["sample"] in trimodal_samples and r.get("modality") == "wes"
                and r.get("method") in ("MajorityVote", "WeightedConsensus")]
    rna_rows = [r for r in bench_rows
                if r["sample"] in trimodal_samples and r.get("modality") == "rnaseq"
                and r.get("method") in ("MajorityVote", "WeightedConsensus")]

    # Load majority_vote_baseline for per-modality MV
    mv_rows = read_tsv(tables / "majority_vote_baseline.tsv")
    wgs_mv = [r for r in mv_rows if r["sample"] in trimodal_samples and r.get("modality") == "wgs"]
    wes_mv = [r for r in mv_rows if r["sample"] in trimodal_samples and r.get("modality") == "wes"]
    rna_mv = [r for r in mv_rows if r["sample"] in trimodal_samples and r.get("modality") == "rnaseq"]

    # Load weighted_consensus_calls for per-modality WC
    wc_rows = read_tsv(tables / "weighted_consensus_calls.tsv")
    wgs_wc = [r for r in wc_rows if r["sample"] in trimodal_samples and r.get("modality") == "wgs"]
    wes_wc = [r for r in wc_rows if r["sample"] in trimodal_samples and r.get("modality") == "wes"]
    rna_wc = [r for r in wc_rows if r["sample"] in trimodal_samples and r.get("modality") == "rnaseq"]

    # Load bimodal rows restricted to trimodal samples
    bi_rows = read_tsv(tables / "bimodal_wes_rna_consensus.tsv")
    bi_mv = [r for r in bi_rows
             if r["sample"] in trimodal_samples and r["method"] == "BimodalMajorityVote"]
    bi_wc = [r for r in bi_rows
             if r["sample"] in trimodal_samples and r["method"] == "BimodalWeightedConsensus"]

    comparisons = [
        ("WGS_MajorityVote",            wgs_mv,  "wgs"),
        ("WGS_WeightedConsensus",        wgs_wc,  "wgs"),
        ("WES_MajorityVote",             wes_mv,  "wes"),
        ("WES_WeightedConsensus",        wes_wc,  "wes"),
        ("RNA_MajorityVote",             rna_mv,  "rnaseq"),
        ("RNA_WeightedConsensus",        rna_wc,  "rnaseq"),
        ("Bimodal_MajorityVote",         bi_mv,   "wes+rnaseq"),
        ("Bimodal_WeightedConsensus",    bi_wc,   "wes+rnaseq"),
        ("Trimodal_MajorityVote",        tri_mv,  "wgs+wes+rnaseq"),
        ("Trimodal_WeightedConsensus",   tri_wc,  "wgs+wes+rnaseq"),
    ]

    comp_rows = []
    comp_fields = ["comparison", "modality", "sample_count", "gene_rows",
                   "callable_rate", "accuracy_among_callable",
                   "overall_correct_call_rate",
                   "overall_correct_call_rate_ci_lo", "overall_correct_call_rate_ci_hi"]

    for label, rows, mod_tag in comparisons:
        if not rows:
            print(f"  WARNING: no rows for {label}")
            continue
        samples_in = {r["sample"] for r in rows}
        comp_rows.append(accuracy_summary(rows, label, mod_tag, samples_in))

    write_tsv(tables / "trimodal_accuracy_comparison.tsv", comp_rows, comp_fields)

    # ── Print summary ─────────────────────────────────────────────────────────
    print("\n── Trimodal accuracy comparison ──────────────────────────────────────")
    print(f"{'Method':<35} {'Samples':>7} {'Gene rows':>9} {'Callable':>8} {'Acc(callable)':>13} {'Overall':>7}")
    for r in comp_rows:
        print(f"{r['comparison']:<35} {r['sample_count']:>7} {r['gene_rows']:>9} "
              f"{float(r['callable_rate']):.1%}  {float(r['accuracy_among_callable']):.2%}"
              f"         {float(r['overall_correct_call_rate']):.2%}")


if __name__ == "__main__":
    main()
