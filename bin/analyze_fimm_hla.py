#!/usr/bin/env python3
"""
analyze_fimm_hla.py — HLA typing analysis for FIMM scRNA cohort.

Identifies FIMM samples from scRNA rows in the combined results Excel file,
matches them to BulkRNA/WES results by allele fingerprint, computes tool
concordance statistics, exports a clean table, and generates publication figures.

Usage:
    module load python-data
    cd /scratch/project_2008084/pihla-publish
    python3 bin/analyze_fimm_hla.py \
        --excel hla_typing_results_combined_v7_chatgpt.xlsx \
        --samples fimm_samples.txt \
        --outdir analysis/fimm_analysis/
"""

import argparse
import os
import sys
import math
import warnings
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

# ── Style constants (matching generate_figures_pub.py) ───────────────────────
DPI = 300
GRID_COLOR = "#E5E7EB"
WONG_PALETTE = ["#E69F00", "#56B4E9", "#009E73", "#F0E442",
                "#0072B2", "#D55E00", "#CC79A7", "#000000"]

TOOL_COLORS = {
    "arcasHLA":  "#0072B2",
    "OptiType":  "#D55E00",
    "SpecHLA":   "#009E73",
}
MODALITY_MARKERS = {
    "SCRNA":   "o",
    "BulkRNA": "s",
    "WES":     "^",
}
MODALITY_LABELS = {
    "SCRNA":   "scRNA",
    "BulkRNA": "Bulk RNA",
    "WES":     "WES",
}
GENES = ["A", "B", "C"]


# ── Allele normalisation ──────────────────────────────────────────────────────

def norm(allele):
    """Strip to 2-field resolution and standardise (None → '')."""
    if allele is None or (isinstance(allele, float) and math.isnan(allele)):
        return ""
    s = str(allele).strip()
    if s in ("", "nan", "None", "-", "NA"):
        return ""
    # Remove leading gene name if present (e.g. "HLA-A*02:01" → "A*02:01")
    if s.startswith("HLA-"):
        s = s[4:]
    # Keep only first two fields: A*02:01:01 → A*02:01
    parts = s.split(":")
    return ":".join(parts[:2])


def allele_pair(a1, a2):
    """Return a sorted tuple of two normalised alleles."""
    return tuple(sorted([norm(a1), norm(a2)]))


def row_fingerprint(row):
    """Return a dict of {gene: allele_pair} for a data row."""
    return {
        "A": allele_pair(row["HLA_A_1"], row["HLA_A_2"]),
        "B": allele_pair(row["HLA_B_1"], row["HLA_B_2"]),
        "C": allele_pair(row["HLA_C_1"], row["HLA_C_2"]),
    }


def fingerprint_score(fp1, fp2):
    """Count matching allele-pairs between two fingerprints (max 3)."""
    return sum(fp1[g] == fp2[g] for g in GENES)


def allele_agreement(pair1, pair2):
    """True if both allele pairs are non-empty and equal."""
    if not all(pair1) or not all(pair2):
        return None  # one side is missing
    return pair1 == pair2


# ── Wilson confidence interval ────────────────────────────────────────────────

def wilson_ci(k, n, z=1.96):
    """Wilson 95% CI for proportion k/n. Returns (lo, hi)."""
    if n == 0:
        return (0.0, 1.0)
    p = k / n
    denom = 1 + z**2 / n
    centre = (p + z**2 / (2 * n)) / denom
    spread = z * math.sqrt(p * (1 - p) / n + z**2 / (4 * n**2)) / denom
    return (max(0.0, centre - spread), min(1.0, centre + spread))


# ── Greedy 1-to-1 matching ────────────────────────────────────────────────────

def greedy_match(ref_fps, query_fps, min_score=2):
    """
    Match each ref row to the best unmatched query row.
    Returns a list of (ref_idx, query_idx | None, score).
    """
    scores = []
    for i, rfp in enumerate(ref_fps):
        for j, qfp in enumerate(query_fps):
            scores.append((fingerprint_score(rfp, qfp), i, j))
    scores.sort(reverse=True)

    used_ref = set()
    used_qry = set()
    matches = {}
    for score, i, j in scores:
        if i in used_ref or j in used_qry:
            continue
        if score >= min_score:
            matches[i] = (j, score)
            used_ref.add(i)
            used_qry.add(j)

    return [(i, matches.get(i, (None, 0))) for i in range(len(ref_fps))]


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--excel",   required=True,
                    help="Excel file path (relative to repo root or absolute)")
    ap.add_argument("--samples", default=None,
                    help="Optional: fimm_samples.txt with sample IDs (one per line)")
    ap.add_argument("--outdir",  default="analysis/fimm_analysis",
                    help="Output directory")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    # ── 1. Load Excel ─────────────────────────────────────────────────────────
    print(f"Loading {args.excel} …")
    df = pd.read_excel(args.excel, sheet_name="all_results_merged", engine="openpyxl")
    print(f"  {len(df)} rows, methods: {sorted(df.Method.unique())}")

    # ── 2. Split by modality ──────────────────────────────────────────────────
    scrna_arc  = df[(df.Seq_type == "SCRNA")  & (df.Method == "arcasHLA")].reset_index(drop=True)
    scrna_opt  = df[(df.Seq_type == "SCRNA")  & (df.Method == "OptiType")].reset_index(drop=True)
    bulk_all   = df[ df.Seq_type == "BulkRNA"].reset_index(drop=True)
    wes_all    = df[ df.Seq_type == "WES"    ].reset_index(drop=True)

    n_fimm = len(scrna_arc)
    print(f"  scRNA arcasHLA: {n_fimm} samples (reference)")
    print(f"  scRNA OptiType: {len(scrna_opt)} samples")

    # ── 3. Assign sample labels ───────────────────────────────────────────────
    sample_ids = []
    if args.samples and os.path.isfile(args.samples):
        with open(args.samples) as f:
            ids = [l.strip() for l in f if l.strip()]
        if len(ids) == n_fimm:
            sample_ids = ids
            print(f"  Using {len(ids)} sample IDs from {args.samples}")
        else:
            print(f"  WARNING: {args.samples} has {len(ids)} IDs but scRNA has {n_fimm} rows "
                  f"— using sequential labels")

    if not sample_ids:
        sample_ids = [f"fimm_{i+1:02d}" for i in range(n_fimm)]

    # Compute fingerprints for reference (arcasHLA scRNA)
    ref_fps = [row_fingerprint(scrna_arc.iloc[i]) for i in range(n_fimm)]

    # ── 4. Match OptiType → arcasHLA within scRNA ────────────────────────────
    opt_fps = [row_fingerprint(scrna_opt.iloc[i]) for i in range(len(scrna_opt))]
    matches_opt = greedy_match(ref_fps, opt_fps, min_score=2)
    # matches_opt[i] = (arcasHLA_idx, (optitype_idx | None, score))

    print("\nScRNA tool matching (arcasHLA → OptiType):")
    scrna_matched_rows = []
    for i, (j_info) in matches_opt:
        j, score = j_info
        sid = sample_ids[i]
        arc_fp = ref_fps[i]
        opt_fp = opt_fps[j] if j is not None else None
        row = {"sample_id": sid, "arc_idx": i, "opt_idx": j, "match_score": score}
        for g in GENES:
            row[f"scrna_arcashla_{g}"] = "|".join(arc_fp[g])
            row[f"scrna_optitype_{g}"] = "|".join(opt_fp[g]) if opt_fp else ""
            row[f"scrna_agree_{g}"] = (
                allele_agreement(arc_fp[g], opt_fp[g]) if opt_fp else None
            )
        scrna_matched_rows.append(row)

    scrna_df = pd.DataFrame(scrna_matched_rows)
    out_scrna = os.path.join(args.outdir, "fimm_scrna_matched.tsv")
    scrna_df.to_csv(out_scrna, sep="\t", index=False)
    print(f"  Saved: {out_scrna}")

    # ── 5. Cross-modality lookup ──────────────────────────────────────────────
    # For each FIMM sample fingerprint, find matching BulkRNA / WES rows
    def lookup_modality(ref_idx, modality_df, min_score=2):
        """Return best-matching rows per method for a given ref fingerprint."""
        rfp = ref_fps[ref_idx]
        results = {}
        for method in modality_df.Method.unique():
            sub = modality_df[modality_df.Method == method].reset_index(drop=True)
            best_score, best_fp = -1, None
            for k in range(len(sub)):
                fp = row_fingerprint(sub.iloc[k])
                s = fingerprint_score(rfp, fp)
                if s > best_score:
                    best_score, best_fp = s, fp
            if best_score >= min_score:
                results[method] = best_fp
        return results

    print("\nCross-modality lookup …")
    cross_rows = []
    for i in range(n_fimm):
        sid = sample_ids[i]
        bulk_match = lookup_modality(i, bulk_all)
        wes_match  = lookup_modality(i, wes_all)

        row = {"sample_id": sid}
        for g in GENES:
            arc_pair = ref_fps[i][g]
            row[f"scrna_arcashla_{g}"]  = "|".join(arc_pair)

            opt_j, _ = matches_opt[i][1]
            opt_pair = opt_fps[opt_j][g] if opt_j is not None else ("", "")
            row[f"scrna_optitype_{g}"]  = "|".join(opt_pair)

            for method, fp in bulk_match.items():
                row[f"bulkrna_{method}_{g}"] = "|".join(fp[g])
            for method, fp in wes_match.items():
                row[f"wes_{method}_{g}"] = "|".join(fp[g])
        cross_rows.append(row)

    cross_df = pd.DataFrame(cross_rows)
    out_cross = os.path.join(args.outdir, "fimm_hla_calls.tsv")
    cross_df.to_csv(out_cross, sep="\t", index=False)
    print(f"  Saved: {out_cross}")

    # ── 6. Concordance statistics ─────────────────────────────────────────────
    print("\nConcordance statistics:")
    stat_rows = []

    comparisons = [
        ("scRNA",    "arcasHLA", "OptiType",  "scrna_arcashla", "scrna_optitype"),
    ]
    # Add bulk and wes vs arcasHLA scRNA
    for method in bulk_all.Method.unique():
        comparisons.append(
            ("BulkRNA", "arcasHLA_scRNA", method,
             None, f"bulkrna_{method}")  # None = use ref_fps directly
        )
    for method in wes_all.Method.unique():
        comparisons.append(
            ("WES", "arcasHLA_scRNA", method,
             None, f"wes_{method}")
        )

    for modality, ref_label, query_label, ref_col_pfx, qry_col_pfx in comparisons:
        for g in GENES:
            agree = []
            for i in range(n_fimm):
                # Reference allele pair
                if ref_col_pfx is None:
                    rp = ref_fps[i][g]
                else:
                    val = scrna_df.loc[i, f"{ref_col_pfx}_{g}"]
                    rp = tuple(val.split("|")) if isinstance(val, str) and val.strip() else ("", "")

                # Query allele pair
                qry_key = f"{qry_col_pfx}_{g}"
                if qry_key in cross_df.columns:
                    val = cross_df.loc[i, qry_key]
                    qp = tuple(val.split("|")) if isinstance(val, str) and val.strip() else ("", "")
                else:
                    continue

                agr = allele_agreement(rp, qp)
                if agr is not None:
                    agree.append(agr)

            n = len(agree)
            k = sum(agree)
            lo, hi = wilson_ci(k, n)
            pct = k / n * 100 if n > 0 else float("nan")
            stat_rows.append({
                "modality": modality,
                "reference": ref_label,
                "query": query_label,
                "gene": g,
                "n_callable": n,
                "n_agree": k,
                "pct_agree": round(pct, 1),
                "ci_lo": round(lo * 100, 1),
                "ci_hi": round(hi * 100, 1),
            })
            print(f"  {modality:8s} {ref_label:20s} vs {query_label:12s} HLA-{g}: "
                  f"{k}/{n} = {pct:.1f}% [{lo*100:.1f}–{hi*100:.1f}]")

    stat_df = pd.DataFrame(stat_rows)
    out_stats = os.path.join(args.outdir, "fimm_concordance_stats.tsv")
    stat_df.to_csv(out_stats, sep="\t", index=False)
    print(f"  Saved: {out_stats}")

    # ── 7. Figures ────────────────────────────────────────────────────────────

    # Figure 1 — scRNA concordance heatmap
    print("\nGenerating Figure 1: scRNA concordance heatmap …")
    fig, ax = plt.subplots(figsize=(10, max(6, n_fimm * 0.28)), dpi=DPI)

    col_labels = [f"HLA-{g}{al}" for g in GENES for al in ("1", "2")]
    # Build matrix: 1=agree, 0=disagree, 0.5=missing OptiType
    mat = np.full((n_fimm, 6), 0.5)  # default: missing (grey)
    for i, row in scrna_df.iterrows():
        for gi, g in enumerate(GENES):
            agr = row.get(f"scrna_agree_{g}")
            if agr is True:
                mat[i, gi * 2]     = 1.0
                mat[i, gi * 2 + 1] = 1.0
            elif agr is False:
                mat[i, gi * 2]     = 0.0
                mat[i, gi * 2 + 1] = 0.0

    from matplotlib.colors import ListedColormap
    cmap = ListedColormap(["#D55E00", "#CCCCCC", "#009E73"])  # disagree, missing, agree
    im = ax.imshow(mat, aspect="auto", cmap=cmap, vmin=0, vmax=1,
                   interpolation="none")

    ax.set_xticks(range(6))
    ax.set_xticklabels(col_labels, fontsize=9)
    ax.set_yticks(range(n_fimm))
    ax.set_yticklabels(sample_ids, fontsize=7)
    ax.set_xlabel("HLA gene allele", fontsize=10)
    ax.set_ylabel("Sample", fontsize=10)
    ax.set_title("scRNA HLA typing: OptiType vs arcasHLA concordance\n(FIMM cohort)",
                 fontsize=11, fontweight="bold")

    legend_patches = [
        mpatches.Patch(color="#009E73", label="Agree"),
        mpatches.Patch(color="#D55E00", label="Disagree"),
        mpatches.Patch(color="#CCCCCC", label="OptiType missing"),
    ]
    ax.legend(handles=legend_patches, bbox_to_anchor=(1.01, 1), loc="upper left",
              fontsize=8, frameon=False)
    plt.tight_layout()
    for ext in ("pdf", "png"):
        p = os.path.join(args.outdir, f"fig_scrna_concordance_heatmap.{ext}")
        fig.savefig(p, dpi=DPI, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: fig_scrna_concordance_heatmap.pdf/.png")

    # Figure 2 — Per-gene concordance bar chart
    print("Generating Figure 2: per-gene concordance bar chart …")
    # Group: scRNA (OptiType vs arcasHLA), BulkRNA best tool, WES best tool
    plot_groups = [
        ("scRNA\nOptiType vs arcasHLA", "scRNA",    "arcasHLA", "OptiType"),
    ]
    for method in sorted(bulk_all.Method.unique()):
        plot_groups.append((f"BulkRNA\nvs {method}", "BulkRNA", "arcasHLA_scRNA", method))
    for method in sorted(wes_all.Method.unique()):
        plot_groups.append((f"WES\nvs {method}", "WES", "arcasHLA_scRNA", method))

    n_groups = len(plot_groups)
    fig, axes = plt.subplots(1, 3, figsize=(4 * n_groups * 0.6 + 2, 5), dpi=DPI,
                              sharey=True)
    if n_groups < 3:
        axes = list(axes)

    for gi, g in enumerate(GENES):
        ax = axes[gi]
        xs, ys, errs_lo, errs_hi, colors = [], [], [], [], []
        for xi, (label, mod, ref_lbl, qry_lbl) in enumerate(plot_groups):
            row = stat_df[
                (stat_df.modality == mod) &
                (stat_df.reference == ref_lbl) &
                (stat_df.query == qry_lbl) &
                (stat_df.gene == g)
            ]
            if row.empty:
                continue
            r = row.iloc[0]
            xs.append(xi)
            ys.append(r.pct_agree)
            errs_lo.append(r.pct_agree - r.ci_lo)
            errs_hi.append(r.ci_hi - r.pct_agree)
            colors.append(WONG_PALETTE[xi % len(WONG_PALETTE)])

        xlabels = [plot_groups[x][0] for x in xs]
        bars = ax.bar(range(len(xs)), ys, color=colors, edgecolor="white",
                      linewidth=0.5, width=0.7)
        ax.errorbar(range(len(xs)), ys,
                    yerr=[errs_lo, errs_hi],
                    fmt="none", color="black", capsize=3, linewidth=1)
        ax.set_xticks(range(len(xs)))
        ax.set_xticklabels(xlabels, fontsize=7, rotation=30, ha="right")
        ax.set_title(f"HLA-{g}", fontsize=10, fontweight="bold")
        ax.set_ylim(0, 105)
        ax.set_ylabel("% concordance" if gi == 0 else "", fontsize=9)
        ax.yaxis.grid(True, color=GRID_COLOR, linewidth=0.5)
        ax.set_axisbelow(True)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

    fig.suptitle("HLA typing concordance — FIMM cohort (vs arcasHLA scRNA reference)",
                 fontsize=10, fontweight="bold", y=1.01)
    plt.tight_layout()
    for ext in ("pdf", "png"):
        p = os.path.join(args.outdir, f"fig_concordance_by_gene.{ext}")
        fig.savefig(p, dpi=DPI, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: fig_concordance_by_gene.pdf/.png")

    # Figure 3 — Cross-modality allele calls dot plot (HLA-A, one allele per sample)
    print("Generating Figure 3: cross-modality allele comparison …")
    # One panel per gene; x=sample, y=allele (string); dots colored by tool
    fig, axes = plt.subplots(3, 1, figsize=(max(12, n_fimm * 0.35), 14), dpi=DPI)

    for gi, g in enumerate(GENES):
        ax = axes[gi]
        # Collect all allele calls per sample per (modality, tool)
        tracks = []
        # scRNA
        for tool_col, tool_label, color in [
            (f"scrna_arcashla_{g}", "scRNA arcasHLA", TOOL_COLORS["arcasHLA"]),
            (f"scrna_optitype_{g}", "scRNA OptiType",  TOOL_COLORS["OptiType"]),
        ]:
            for i, sid in enumerate(sample_ids):
                val = cross_df.loc[i, tool_col] if tool_col in cross_df.columns else ""
                if isinstance(val, str) and val.strip():
                    for allele in val.split("|"):
                        if allele:
                            tracks.append((i, allele, tool_label, color, "o"))
        # BulkRNA and WES
        for method in list(bulk_all.Method.unique()) + list(wes_all.Method.unique()):
            mod = "BulkRNA" if method in bulk_all.Method.unique() else "WES"
            col_key = f"{'bulkrna' if mod == 'BulkRNA' else 'wes'}_{method}_{g}"
            marker = MODALITY_MARKERS.get(mod, "D")
            color = TOOL_COLORS.get(method, WONG_PALETTE[4])
            label = f"{MODALITY_LABELS[mod]} {method}"
            for i, sid in enumerate(sample_ids):
                val = cross_df.loc[i, col_key] if col_key in cross_df.columns else ""
                if isinstance(val, str) and val.strip():
                    for allele in val.split("|"):
                        if allele:
                            tracks.append((i, allele, label, color, marker))

        if not tracks:
            ax.set_visible(False)
            continue

        # Get unique alleles sorted
        all_alleles = sorted(set(t[1] for t in tracks))
        allele_idx = {a: j for j, a in enumerate(all_alleles)}

        seen_labels = set()
        for xi, allele, label, color, marker in tracks:
            yi = allele_idx[allele]
            lbl = label if label not in seen_labels else None
            ax.scatter(xi, yi, color=color, marker=marker, s=25, alpha=0.75,
                       label=lbl, linewidths=0.3, edgecolors="white")
            seen_labels.add(label)

        ax.set_xticks(range(n_fimm))
        ax.set_xticklabels(sample_ids, fontsize=6, rotation=60, ha="right")
        ax.set_yticks(range(len(all_alleles)))
        ax.set_yticklabels(all_alleles, fontsize=7)
        ax.set_title(f"HLA-{g} allele calls across modalities", fontsize=10,
                     fontweight="bold")
        ax.yaxis.grid(True, color=GRID_COLOR, linewidth=0.4)
        ax.set_axisbelow(True)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        if gi == 0:
            ax.legend(bbox_to_anchor=(1.01, 1), loc="upper left",
                      fontsize=7, frameon=False, ncol=1)

    fig.suptitle("HLA allele calls per sample — FIMM cohort (all tools and modalities)",
                 fontsize=11, fontweight="bold")
    plt.tight_layout()
    for ext in ("pdf", "png"):
        p = os.path.join(args.outdir, f"fig_crossmodal_comparison.{ext}")
        fig.savefig(p, dpi=DPI, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: fig_crossmodal_comparison.pdf/.png")

    # ── Summary ───────────────────────────────────────────────────────────────
    print(f"\nDone. All outputs in: {args.outdir}/")
    print("  fimm_scrna_matched.tsv    — scRNA rows with tool match")
    print("  fimm_hla_calls.tsv        — wide-format table (sample × tool × modality)")
    print("  fimm_concordance_stats.tsv — per-gene concordance rates")
    print("  fig_scrna_concordance_heatmap.pdf/.png")
    print("  fig_concordance_by_gene.pdf/.png")
    print("  fig_crossmodal_comparison.pdf/.png")


if __name__ == "__main__":
    main()
