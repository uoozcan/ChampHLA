#!/usr/bin/env python3
"""validate_hla_a2_flow.py — validate computational HLA typing against flow-cytometry HLA-A2.

Flow cytometry (surface HLA-A2 antibody, typically BB7.2) is a protein-level assay and
serves as an orthogonal ground truth for the HLA-A2 serotype. This script derives an
HLA-A2 call from each donor's computational HLA-A alleles and compares it against the
flow-cytometry ground truth, reporting a 2x2 confusion matrix plus
sensitivity/specificity/PPV/NPV/accuracy (with Wilson 95% CIs).

HLA-A2 is a germline trait, so samples are matched at the DONOR level (across timepoints /
sequencing runs). Two serotype-mapping definitions are reported side by side:
  * strict:        A2-positive iff an HLA-A allele group is A*02
  * cross-reactive: A2-positive iff an HLA-A allele group is A*02, A*68 or A*69
                    (BB7.2 antibody cross-reactivity)

Usage:
    module load python-data/3.12
    python3 bin/validate_hla_a2_flow.py \
        --flow      /scratch/project_2008084/AB_HLA_A2_typing.xlsx \
        --calls     /scratch/project_2008084/pihla_local/fimm_results/fimm_hla_calls.tsv \
        --benchmark /scratch/project_2008084/pihla-publish/bin/hla_benchmark.py \
        --outdir    /scratch/project_2008084/pihla-publish/fimm_results/hla_a2_flow
"""

import argparse
import csv
import importlib.util
import os
import re
import sys
from collections import Counter, defaultdict

import pandas as pd


# ── Reuse hla_benchmark.py (normalize_allele, wilson_ci) ─────────────────────────
def load_benchmark(path):
    spec = importlib.util.spec_from_file_location("hla_benchmark", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


# ── Serotype mapping ─────────────────────────────────────────────────────────────
STRICT_GROUPS = {"A*02"}
CROSS_GROUPS = {"A*02", "A*68", "A*69"}
MAPPINGS = {"strict": STRICT_GROUPS, "cross": CROSS_GROUPS}

# Modality x tool columns present in fimm_hla_calls.tsv (gene A only)
MODALITY_TOOLS = [
    ("scrna", "arcashla"),
    ("scrna", "optitype"),
    ("bulkrna", "OptiType"),
    ("bulkrna", "arcasHLA"),
    ("bulkrna", "SpecHLA"),
    ("wes", "OptiType"),
    ("wes", "arcasHLA"),
    ("wes", "SpecHLA"),
]


def col_a(mod, tool):
    return f"{mod}_{tool}_A"


# ── Donor-key normalization ──────────────────────────────────────────────────────
def donor_key(raw):
    """Reduce a sample id (flow ID/Long ID or call sample_id) to a germline donor key.

    Examples:
      FH_10050_2                -> fh_10050
      FHRB_MV_55_41815          -> fhrb_mv_55
      vcp-4555_020615_2_D1      -> fhrb_4555   (matched via number; see resolve)
      FH_6753_14082017_9999_BM  -> fh_6753
      MV_Nov.14_Dg              -> mv_nov_14
    """
    s = str(raw).strip().lower()
    if not s:
        return ""
    s = s.replace("-", "_").replace(".", "_")
    if s.startswith("vcp_"):
        s = s[4:]
    parts = [p for p in s.split("_") if p != ""]
    if not parts:
        return ""
    # registry / cohort prefix tokens that precede the numeric donor id
    prefixes = {"fh", "fhrb", "ort", "mv", "nov"}
    key_tokens = []
    for tok in parts:
        if tok in prefixes:
            key_tokens.append(tok)
            continue
        # first non-prefix token: the donor number (may be like "14" or "55")
        key_tokens.append(tok)
        break
    return "_".join(key_tokens)


# ── Load flow ground truth ───────────────────────────────────────────────────────
def load_flow(path):
    df = pd.read_excel(path, engine="openpyxl", sheet_name=0)
    df.columns = [str(c).strip() for c in df.columns]
    a2_col = next(c for c in df.columns if c.lower().startswith("hla-a2"))
    rows = []
    for _, r in df.iterrows():
        sid = str(r.get("ID", "") or "").strip()
        long_id = str(r.get("Long ID", "") or "").strip()
        a2 = str(r.get(a2_col, "") or "").strip().lower()
        if not sid or a2 not in ("yes", "no"):
            continue
        # ID starting with a digit (e.g. 4555_020615_2) -> use Long ID for registry/number
        base = sid if not sid[0].isdigit() else long_id
        rows.append(
            {
                "flow_id": sid,
                "long_id": long_id,
                "donor": donor_key(base),
                "flow_a2": a2,
            }
        )
    return rows


# ── Load computational calls, grouped by donor ───────────────────────────────────
def load_calls(path):
    """Return {donor: {(mod,tool): set(allele_group_strings...)}, '_samples': [...]}."""
    by_donor = defaultdict(lambda: {"samples": [], "mt": defaultdict(list)})
    with open(path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            sid = row["sample_id"]
            dk = donor_key(sid)
            by_donor[dk]["samples"].append(sid)
            for mod, tool in MODALITY_TOOLS:
                cell = (row.get(col_a(mod, tool)) or "").strip()
                if cell and cell != "|":
                    by_donor[dk]["mt"][(mod, tool)].append(cell)
    return by_donor


def cell_groups(cell, normalize_allele):
    """Field-1 allele groups (e.g. {'A*02','A*03'}) from a 'A1|A2' cell."""
    groups = set()
    for part in cell.split("|"):
        g = normalize_allele(part, resolution=1)
        if g:
            groups.add(g)
    return groups


def predict_a2(groups, positive_groups):
    return "yes" if (groups & positive_groups) else "no"


# ── Metrics ──────────────────────────────────────────────────────────────────────
def confusion(pairs):
    """pairs: list of (truth, pred) in {'yes','no'}. Returns TP,FP,TN,FN."""
    tp = sum(1 for t, p in pairs if t == "yes" and p == "yes")
    fp = sum(1 for t, p in pairs if t == "no" and p == "yes")
    tn = sum(1 for t, p in pairs if t == "no" and p == "no")
    fn = sum(1 for t, p in pairs if t == "yes" and p == "no")
    return tp, fp, tn, fn


def ratio(k, n):
    return round(k / n, 4) if n else None


def fmt_ci(k, n, wilson_ci):
    if not n:
        return ""
    lo, hi = wilson_ci(k, n)
    return f"[{lo:.3f}-{hi:.3f}]"


def metrics_row(scope, mapping, pairs, wilson_ci):
    tp, fp, tn, fn = confusion(pairs)
    n = tp + fp + tn + fn
    return {
        "scope": scope,
        "mapping": mapping,
        "N": n,
        "TP": tp,
        "FP": fp,
        "TN": tn,
        "FN": fn,
        "sensitivity": ratio(tp, tp + fn),
        "sensitivity_95CI": fmt_ci(tp, tp + fn, wilson_ci),
        "specificity": ratio(tn, tn + fp),
        "specificity_95CI": fmt_ci(tn, tn + fp, wilson_ci),
        "PPV": ratio(tp, tp + fp),
        "NPV": ratio(tn, tn + fn),
        "accuracy": ratio(tp + tn, n),
        "accuracy_95CI": fmt_ci(tp + tn, n, wilson_ci),
    }


def consensus_call(preds):
    """Majority vote over a list of 'yes'/'no'. Returns (call, tie_flag)."""
    if not preds:
        return None, False
    c = Counter(preds)
    yes, no = c.get("yes", 0), c.get("no", 0)
    if yes == no:
        return "tie", True
    return ("yes" if yes > no else "no"), False


def write_tsv(path, rows, fieldnames):
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow(r)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--flow", required=True)
    ap.add_argument("--calls", required=True)
    ap.add_argument("--benchmark", required=True)
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()

    hb = load_benchmark(args.benchmark)
    normalize_allele, wilson_ci = hb.normalize_allele, hb.wilson_ci

    os.makedirs(args.outdir, exist_ok=True)

    flow = load_flow(args.flow)
    calls = load_calls(args.calls)

    # per-donor flow truth + single-valuedness assertion
    donor_truth = defaultdict(set)
    donor_flow_ids = defaultdict(list)
    for r in flow:
        donor_truth[r["donor"]].add(r["flow_a2"])
        donor_flow_ids[r["donor"]].append(r["flow_id"])
    conflicts = {d: v for d, v in donor_truth.items() if len(v) > 1}
    if conflicts:
        sys.exit(f"ERROR: flow HLA-A2 status is not single-valued for donor(s): {conflicts}")

    matched = sorted(d for d in donor_truth if d in calls and calls[d]["mt"])
    unmatched = sorted(d for d in donor_truth if d not in matched)

    # ── per-donor detail + prediction collection ────────────────────────────────
    per_donor_rows = []
    # pairs[mapping]['consensus'] and pairs[mapping][(mod,tool)]
    pairs = {m: defaultdict(list) for m in MAPPINGS}

    for d in matched:
        truth = next(iter(donor_truth[d]))
        row = {
            "donor": d,
            "flow_a2": truth,
            "flow_ids": ",".join(sorted(set(donor_flow_ids[d]))),
            "call_samples": ",".join(sorted(set(calls[d]["samples"]))),
        }
        # gather groups per (mod,tool) across all this donor's calls
        mt_groups = {}
        for mt, cells in calls[d]["mt"].items():
            g = set()
            for cell in cells:
                g |= cell_groups(cell, normalize_allele)
            mt_groups[mt] = g

        for mapping, pos_groups in MAPPINGS.items():
            mt_preds = []
            for mod, tool in MODALITY_TOOLS:
                mt = (mod, tool)
                if mt not in mt_groups:
                    continue
                pred = predict_a2(mt_groups[mt], pos_groups)
                mt_preds.append(pred)
                pairs[mapping][mt].append((truth, pred))
                row[f"{mod}_{tool}_{mapping}"] = pred
            cons, tie = consensus_call(mt_preds)
            row[f"consensus_{mapping}"] = cons
            # ties excluded from consensus metrics (cannot score a non-call)
            if cons in ("yes", "no"):
                pairs[mapping]["consensus"].append((truth, cons))
                row[f"consensus_{mapping}_agree"] = "agree" if cons == truth else "DISAGREE"
            else:
                row[f"consensus_{mapping}_agree"] = "tie"

        # record allele groups seen (for readability), strict-independent
        for mod, tool in MODALITY_TOOLS:
            mt = (mod, tool)
            if mt in mt_groups:
                row[f"{mod}_{tool}_alleles"] = "|".join(sorted(mt_groups[mt]))
        per_donor_rows.append(row)

    # ── metrics table ───────────────────────────────────────────────────────────
    metrics_rows = []
    for mapping in MAPPINGS:
        metrics_rows.append(metrics_row("consensus", mapping, pairs[mapping]["consensus"], wilson_ci))
    for mod, tool in MODALITY_TOOLS:
        for mapping in MAPPINGS:
            mt_pairs = pairs[mapping].get((mod, tool), [])
            if mt_pairs:
                metrics_rows.append(metrics_row(f"{mod}_{tool}", mapping, mt_pairs, wilson_ci))

    # ── unmatched (coverage gap) ────────────────────────────────────────────────
    unmatched_rows = []
    for r in flow:
        if r["donor"] in unmatched:
            unmatched_rows.append(
                {"donor": r["donor"], "flow_id": r["flow_id"], "flow_a2": r["flow_a2"]}
            )

    # ── write outputs ───────────────────────────────────────────────────────────
    detail_cols = (
        ["donor", "flow_a2", "flow_ids", "call_samples"]
        + [f"{m}_{t}_alleles" for m, t in MODALITY_TOOLS]
        + [f"{m}_{t}_{mp}" for m, t in MODALITY_TOOLS for mp in MAPPINGS]
        + [f"consensus_{mp}" for mp in MAPPINGS]
        + [f"consensus_{mp}_agree" for mp in MAPPINGS]
    )
    metric_cols = [
        "scope", "mapping", "N", "TP", "FP", "TN", "FN",
        "sensitivity", "sensitivity_95CI", "specificity", "specificity_95CI",
        "PPV", "NPV", "accuracy", "accuracy_95CI",
    ]
    p_detail = os.path.join(args.outdir, "hla_a2_flow_per_donor.tsv")
    p_metrics = os.path.join(args.outdir, "hla_a2_flow_metrics.tsv")
    p_unmatched = os.path.join(args.outdir, "hla_a2_flow_unmatched.tsv")
    write_tsv(p_detail, per_donor_rows, detail_cols)
    write_tsv(p_metrics, metrics_rows, metric_cols)
    write_tsv(p_unmatched, unmatched_rows, ["donor", "flow_id", "flow_a2"])

    # ── stdout summary ──────────────────────────────────────────────────────────
    n_pos = sum(1 for d in matched if next(iter(donor_truth[d])) == "yes")
    n_neg = len(matched) - n_pos
    print("=" * 72)
    print("HLA-A2 flow-cytometry vs computational typing — use-case validation")
    print("=" * 72)
    print(f"Flow donors total          : {len(donor_truth)}")
    print(f"Matched (have HLA-A calls) : {len(matched)}  (A2+ {n_pos}, A2- {n_neg})")
    print(f"Unmatched (coverage gap)   : {len(unmatched)} donors / "
          f"{len(unmatched_rows)} flow samples")
    print()
    print("Consensus (per-donor majority vote) metrics:")
    for m in metrics_rows:
        if m["scope"] != "consensus":
            continue
        print(f"  [{m['mapping']:5s}] N={m['N']:2d}  TP={m['TP']} FP={m['FP']} "
              f"TN={m['TN']} FN={m['FN']}  "
              f"sens={m['sensitivity']} {m['sensitivity_95CI']}  "
              f"spec={m['specificity']} {m['specificity_95CI']}  "
              f"acc={m['accuracy']} {m['accuracy_95CI']}")
    print()
    print("Outputs:")
    print(f"  {p_detail}")
    print(f"  {p_metrics}")
    print(f"  {p_unmatched}")


if __name__ == "__main__":
    main()
