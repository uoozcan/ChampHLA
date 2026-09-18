#!/usr/bin/env python3
"""champhla_detail_flow_cohort.py — detailed champHLA typing for the flow-cytometry cohort.

For every flow-cytometry donor that champHLA (pihla) has typed, this reports the FULL
HLA-A/B/C genotype across all modality x tool combinations (not just the binary HLA-A2
call), the majority consensus 2-field genotype per locus, cross-tool concordance,
zygosity, and — for locus A — the flow HLA-A2 truth and the A*02 subtype called.

Reuses donor_key / load_flow from validate_hla_a2_flow.py and
normalize_allele / sort_alleles from hla_benchmark.py.

Usage:
    module load python-data/3.12
    python3 bin/champhla_detail_flow_cohort.py \
        --flow      /scratch/project_2008084/AB_HLA_A2_typing.xlsx \
        --calls     /scratch/project_2008084/pihla_local/fimm_results/fimm_hla_calls.tsv \
        --benchmark /scratch/project_2008084/pihla-publish/bin/hla_benchmark.py \
        --validator /scratch/project_2008084/pihla-publish/bin/validate_hla_a2_flow.py \
        --outdir    /scratch/project_2008084/pihla-publish/fimm_results/hla_a2_flow
"""
import argparse
import csv
import importlib.util
import os
from collections import Counter, defaultdict

GENES = ["A", "B", "C"]
MODALITY_TOOLS = [
    ("scrna", "arcashla"), ("scrna", "optitype"),
    ("bulkrna", "arcasHLA"), ("bulkrna", "OptiType"), ("bulkrna", "SpecHLA"),
    ("wes", "arcasHLA"), ("wes", "OptiType"), ("wes", "SpecHLA"),
]
MT_DISPLAY = {
    ("scrna", "arcashla"): "scRNA_arcasHLA", ("scrna", "optitype"): "scRNA_OptiType",
    ("bulkrna", "arcasHLA"): "bulkRNA_arcasHLA", ("bulkrna", "OptiType"): "bulkRNA_OptiType",
    ("bulkrna", "SpecHLA"): "bulkRNA_SpecHLA",
    ("wes", "arcasHLA"): "WES_arcasHLA", ("wes", "OptiType"): "WES_OptiType",
    ("wes", "SpecHLA"): "WES_SpecHLA",
}


def load_module(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def genotype_2field(cell, sort_alleles):
    """Return a sorted 2-field ('A*02:01','A*03:01') tuple from a 'A1|A2' cell, or None."""
    cell = (cell or "").strip()
    if not cell or cell == "|":
        return None
    pair = sort_alleles(cell.split("|"), resolution=2)
    return pair if all(pair) else None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--flow", required=True)
    ap.add_argument("--calls", required=True)
    ap.add_argument("--benchmark", required=True)
    ap.add_argument("--validator", required=True)
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()

    hb = load_module(args.benchmark, "hla_benchmark")
    vd = load_module(args.validator, "validate_hla_a2_flow")
    donor_key, normalize_allele, sort_alleles = vd.donor_key, hb.normalize_allele, hb.sort_alleles

    # flow truth per donor
    flow = vd.load_flow(args.flow)
    flow_a2 = {}
    for r in flow:
        flow_a2[r["donor"]] = r["flow_a2"]

    # calls grouped by donor -> {gene: {sample: {(mod,tool): genotype}}}
    by_donor = defaultdict(lambda: {"samples": set(),
                                    "g": {g: defaultdict(dict) for g in GENES}})
    with open(args.calls) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            dk = donor_key(row["sample_id"])
            sid = row["sample_id"]
            by_donor[dk]["samples"].add(sid)
            for g in GENES:
                for mt in MODALITY_TOOLS:
                    gt = genotype_2field(row.get(f"{mt[0]}_{mt[1]}_{g}"), sort_alleles)
                    if gt:
                        by_donor[dk]["g"][g][sid][mt] = gt

    matched = sorted(d for d in flow_a2 if d in by_donor and
                     any(by_donor[d]["g"][g] for g in GENES))

    def irreconcilable(genos):
        """True if distinct genotypes share NO allele group (A*02 vs A*03/A*68 = swap),
        as opposed to mere sub-typing noise (A*03:04 vs A*03:16 share group A*03)."""
        grp_sets = [{normalize_allele(a, 1) for a in g} for g in genos]
        return len(genos) > 1 and not set.intersection(*grp_sets)

    by_locus_rows, wide_rows, concord_rows = [], [], []
    for d in matched:
        wide = {"donor": d, "flow_A2": flow_a2[d], "n_samples": len(by_donor[d]["samples"])}
        for g in GENES:
            per_sample = by_donor[d]["g"][g]           # sample -> {(mod,tool): geno}
            all_genos = [gt for mts in per_sample.values() for gt in mts.values()]
            tools_seen = {mt for mts in per_sample.values() for mt in mts}
            n_tools = len(tools_seen)
            distinct = set(all_genos)
            n_distinct = len(distinct)
            # sample-weighted consensus: per-sample modal genotype, then modal across samples
            sample_modal = [Counter(mts.values()).most_common(1)[0][0]
                            for mts in per_sample.values() if mts]
            consensus = Counter(sample_modal).most_common(1)[0][0] if sample_modal else ("", "")
            unanimous = n_distinct == 1 and n_tools > 0
            conflict = irreconcilable(distinct)
            # pairwise 2-field agreement across all (sample,tool) calls
            gl = all_genos
            pairs = [(i, j) for i in range(len(gl)) for j in range(i + 1, len(gl))]
            pw = round(sum(1 for i, j in pairs if gl[i] == gl[j]) / len(pairs), 3) if pairs else None
            zyg = "hom" if (consensus[0] and consensus[0] == consensus[1]) else \
                  ("het" if consensus[0] else "")
            cons_str = consensus[0] if zyg == "hom" else "/".join(consensus) if consensus[0] else ""
            # per (mod,tool) modal genotype across samples, for the display columns
            mt_geno = {}
            for mt in tools_seen:
                gts = [mts[mt] for mts in per_sample.values() if mt in mts]
                mt_geno[mt] = Counter(gts).most_common(1)[0][0]

            row = {"donor": d, "locus": g, "flow_A2": flow_a2[d] if g == "A" else "",
                   "consensus_2field": cons_str, "zygosity": zyg,
                   "n_tools": n_tools, "n_distinct_genotypes": n_distinct,
                   "unanimous": "yes" if unanimous else "no",
                   "within_donor_conflict": "SWAP?" if conflict else "no",
                   "distinct_genotypes": " ; ".join("/".join(x) for x in sorted(distinct)),
                   "pairwise_agree_2field": pw}
            if g == "A":
                a2sub = sorted({a for a in consensus if normalize_allele(a, 1) == "A*02"})
                cons_a2 = "yes" if a2sub else "no"
                row["A2_subtype"] = ",".join(a2sub)
                row["consensus_A2"] = cons_a2
                row["A2_matches_flow"] = "yes" if cons_a2 == flow_a2[d] else "NO"
            for mt in MODALITY_TOOLS:
                row[MT_DISPLAY[mt]] = "/".join(mt_geno[mt]) if mt in mt_geno else ""
            by_locus_rows.append(row)

            wide[f"{g}_consensus"] = cons_str
            wide[f"{g}_zygosity"] = zyg
            wide[f"{g}_n_tools"] = n_tools
            wide[f"{g}_unanimous"] = "yes" if unanimous else "no"
            wide[f"{g}_conflict"] = "SWAP?" if conflict else ""
            concord_rows.append({"donor": d, "locus": g, "unanimous": unanimous,
                                 "conflict": conflict, "n_tools": n_tools, "pw": pw})
        arow = next(r for r in by_locus_rows if r["donor"] == d and r["locus"] == "A")
        wide["A2_subtype"] = arow.get("A2_subtype", "")
        wide["A_consensus_A2"] = arow.get("consensus_A2", "")
        wide["A2_matches_flow"] = arow.get("A2_matches_flow", "")
        wide_rows.append(wide)

    # ── write TSVs ───────────────────────────────────────────────────────────────
    os.makedirs(args.outdir, exist_ok=True)
    locus_cols = (["donor", "locus", "flow_A2", "A2_subtype", "consensus_A2",
                   "A2_matches_flow", "consensus_2field", "zygosity", "n_tools",
                   "n_distinct_genotypes", "unanimous", "within_donor_conflict",
                   "distinct_genotypes", "pairwise_agree_2field"]
                  + [MT_DISPLAY[mt] for mt in MODALITY_TOOLS])
    wide_cols = (["donor", "flow_A2", "A2_subtype", "A_consensus_A2", "A2_matches_flow",
                  "n_samples"]
                 + [f"{g}_{k}" for g in GENES for k in ("consensus", "zygosity",
                                                        "n_tools", "unanimous", "conflict")])
    p_locus = os.path.join(args.outdir, "champhla_detail_by_locus.tsv")
    p_wide = os.path.join(args.outdir, "champhla_detail_wide.tsv")
    with open(p_locus, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=locus_cols, delimiter="\t", extrasaction="ignore")
        w.writeheader(); [w.writerow(r) for r in by_locus_rows]
    with open(p_wide, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=wide_cols, delimiter="\t", extrasaction="ignore")
        w.writeheader(); [w.writerow(r) for r in wide_rows]

    # ── per-locus concordance summary to stdout ──────────────────────────────────
    print("=" * 68)
    print("champHLA detailed typing — flow-cytometry cohort")
    print("=" * 68)
    print(f"Matched donors: {len(matched)}")
    print("\nCross-tool concordance (2-field, per locus, across matched donors):")
    for g in GENES:
        rows = [r for r in concord_rows if r["locus"] == g and r["n_tools"] > 1]
        if not rows:
            continue
        unan = sum(1 for r in rows if r["unanimous"])
        pw = [r["pw"] for r in rows if r["pw"] is not None]
        print(f"  HLA-{g}: {unan}/{len(rows)} donors fully concordant "
              f"({unan/len(rows)*100:.0f}%), mean pairwise 2-field agreement "
              f"{sum(pw)/len(pw):.2f}  (donors with >=2 tools)")

    conflicts = [r for r in by_locus_rows if r["within_donor_conflict"] == "SWAP?"]
    if conflicts:
        print("\nWithin-donor genotype CONFLICTS (irreconcilable calls — possible "
              "sample swap/contamination; flow adjudicates):")
        for r in conflicts:
            print(f"  {r['donor']} HLA-{r['locus']}: {r['distinct_genotypes']}"
                  + (f"   | flow A2={r['flow_A2']} -> flow-consistent genotype is the "
                     f"{'A*02' if r['flow_A2']=='yes' else 'non-A*02'} one" if r['locus']=='A' else ""))

    mism = [r for r in by_locus_rows if r["locus"] == "A" and r.get("A2_matches_flow") == "NO"]
    print(f"\nConsensus A2 vs flow: {len(matched)-len(mism)}/{len(matched)} donors match "
          f"(consensus derived from sample-weighted genotype).")
    for r in mism:
        print(f"  MISMATCH {r['donor']}: consensus {r['consensus_2field']} "
              f"(A2={r['consensus_A2']}) vs flow A2={r['flow_A2']}")
    print("\nOutputs:")
    print(f"  {p_locus}")
    print(f"  {p_wide}")


if __name__ == "__main__":
    main()
