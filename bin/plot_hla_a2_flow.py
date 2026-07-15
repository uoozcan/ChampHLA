#!/usr/bin/env python3
"""plot_hla_a2_flow.py — figures for the flow-cytometry HLA-A2 validation.

Reads the TSVs produced by validate_hla_a2_flow.py and emits publication PNGs:
  fig1_confusion.png        consensus 2x2 confusion matrices (strict & cross-reactive)
  fig2_sens_spec.png        sensitivity & specificity by method (strict) with Wilson CIs
  fig3_donor_method.png     donor x method A2-call matrix vs flow truth (strict)
  fig4_mapping_effect.png   specificity: strict vs cross-reactive (the A*68/A*69 effect)

Usage:
    module load python-data/3.12
    python3 bin/plot_hla_a2_flow.py --indir <outdir-of-validate> [--outdir <indir>/figures]
"""
import argparse
import csv
import importlib.util
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

# ── palette (dataviz reference instance) ─────────────────────────────────────────
BLUE = "#2a78d6"      # categorical slot 1  (sensitivity)
AQUA = "#1baf7a"      # categorical slot 2  (specificity)
GOOD = "#0ca30c"      # status: correct
CRIT = "#d03b3b"      # status: error
GOOD_L = "#d5f0d5"    # correct, negative (light tint)
CRIT_L = "#f7d9d9"
GREY = "#e6e5e0"      # missing
INK = "#0b0b0b"
INK2 = "#52514e"
MUTED = "#898781"
GRID = "#e1e0d9"
SURF = "#fcfcfb"

plt.rcParams.update({
    "figure.facecolor": SURF, "axes.facecolor": SURF, "savefig.facecolor": SURF,
    "font.family": "sans-serif", "font.sans-serif": ["DejaVu Sans"],
    "text.color": INK, "axes.labelcolor": INK, "axes.edgecolor": "#c3c2b7",
    "xtick.color": INK2, "ytick.color": INK2, "axes.titlecolor": INK,
    "font.size": 10, "axes.titleweight": "bold",
})

METHOD_ORDER = [
    ("consensus", "Consensus"),
    ("scrna_arcashla", "scRNA arcasHLA"),
    ("scrna_optitype", "scRNA OptiType"),
    ("bulkrna_arcasHLA", "bulkRNA arcasHLA"),
    ("bulkrna_OptiType", "bulkRNA OptiType"),
    ("bulkrna_SpecHLA", "bulkRNA SpecHLA"),
    ("wes_arcasHLA", "WES arcasHLA"),
    ("wes_OptiType", "WES OptiType"),
    ("wes_SpecHLA", "WES SpecHLA"),
]
DETAIL_COLS = [  # (per_donor column stem, display) matching METHOD_ORDER minus consensus
    ("scrna_arcashla", "scRNA arcasHLA"),
    ("scrna_optitype", "scRNA OptiType"),
    ("bulkrna_arcasHLA", "bulkRNA arcasHLA"),
    ("bulkrna_OptiType", "bulkRNA OptiType"),
    ("bulkrna_SpecHLA", "bulkRNA SpecHLA"),
    ("wes_arcasHLA", "WES arcasHLA"),
    ("wes_OptiType", "WES OptiType"),
    ("wes_SpecHLA", "WES SpecHLA"),
]


def load_metrics(path):
    rows = {}
    with open(path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            rows[(r["scope"], r["mapping"])] = r
    return rows


def load_donors(path):
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def load_wilson(benchmark):
    spec = importlib.util.spec_from_file_location("hla_benchmark", benchmark)
    mod = importlib.util.module_from_spec(spec); spec.loader.exec_module(mod)
    return mod.wilson_ci


# ── Fig 1: confusion matrices ────────────────────────────────────────────────────
def fig_confusion(metrics, out):
    fig, axes = plt.subplots(1, 2, figsize=(8.2, 4.0))
    for ax, mapping, title in zip(axes, ("strict", "cross"),
                                  ("Strict  (A*02 only)", "Cross-reactive  (A*02 / A*68 / A*69)")):
        m = metrics[("consensus", mapping)]
        tp, fp, tn, fn = int(m["TP"]), int(m["FP"]), int(m["TN"]), int(m["FN"])
        # grid: rows = actual (flow) +/- ; cols = predicted +/-
        cells = [[(tp, "TP"), (fn, "FN")],
                 [(fp, "FP"), (tn, "TN")]]
        for i in range(2):
            for j in range(2):
                val, lab = cells[i][j]
                correct = lab in ("TP", "TN")
                if correct:
                    face = GOOD if lab == "TP" else GOOD_L
                    txt_c = "#ffffff" if lab == "TP" else INK
                elif val > 0:                      # a real error
                    face, txt_c = CRIT, "#ffffff"
                else:                              # empty error cell — stay neutral
                    face, txt_c = GREY, MUTED
                ax.add_patch(plt.Rectangle((j, 1 - i), 1, 1, facecolor=face,
                                           edgecolor=SURF, linewidth=3, zorder=1))
                ax.text(j + 0.5, 1 - i + 0.58, str(val), ha="center", va="center",
                        fontsize=22, fontweight="bold", color=txt_c, zorder=2)
                ax.text(j + 0.5, 1 - i + 0.24, lab, ha="center", va="center",
                        fontsize=9, color=txt_c, zorder=2)
        ax.set_xlim(0, 2); ax.set_ylim(0, 2)
        ax.set_xticks([0.5, 1.5]); ax.set_xticklabels(["A2 +", "A2 −"])
        ax.set_yticks([1.5, 0.5]); ax.set_yticklabels(["A2 +", "A2 −"])
        ax.set_xlabel("champHLA consensus call", color=INK2)
        ax.set_ylabel("Flow cytometry (truth)", color=INK2)
        acc = m["accuracy"]
        ax.set_title(f"{title}\nN={m['N']}  ·  accuracy {float(acc):.2f}", fontsize=10)
        for s in ax.spines.values():
            s.set_visible(False)
        ax.tick_params(length=0)
    fig.suptitle("HLA-A2: champHLA consensus vs flow cytometry", fontsize=12, fontweight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig(out, dpi=200); plt.close(fig)


# ── Fig 2: sensitivity & specificity by method (strict) ──────────────────────────
def _metric_ci(m, kind, wilson):
    if kind == "sens":
        k, n = int(m["TP"]), int(m["TP"]) + int(m["FN"])
    else:
        k, n = int(m["TN"]), int(m["TN"]) + int(m["FP"])
    val = k / n if n else float("nan")
    lo, hi = wilson(k, n) if n else (float("nan"), float("nan"))
    return val, lo, hi, n


def fig_sens_spec(metrics, wilson, out):
    labels = [d for _, d in METHOD_ORDER]
    y = list(range(len(labels)))[::-1]
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 5.2), sharey=True)
    for ax, kind, color, title in ((axes[0], "sens", BLUE, "Sensitivity  (detect A2 +)"),
                                   (axes[1], "spec", AQUA, "Specificity  (detect A2 −)")):
        for yi, (scope, _) in zip(y, METHOD_ORDER):
            m = metrics.get((scope, "strict"))
            if not m:
                continue
            val, lo, hi, n = _metric_ci(m, kind, wilson)
            ax.barh(yi, val, height=0.6, color=color, zorder=2,
                    edgecolor=SURF, linewidth=1)
            ax.plot([lo, hi], [yi, yi], color=INK2, lw=1.5, zorder=3)
            ax.plot([lo, lo, None, hi, hi], [yi - .12, yi + .12, None, yi - .12, yi + .12],
                    color=INK2, lw=1.5, zorder=3)
            ax.text(max(val, hi) + 0.04, yi, f"{val:.2f}  (n={n})", va="center",
                    ha="left", fontsize=8.5, color=INK)
        ax.set_xlim(0, 1.28); ax.set_xticks([0, .25, .5, .75, 1.0])
        ax.set_yticks(y); ax.set_yticklabels(labels, fontsize=9)
        ax.set_title(title, fontsize=10.5)
        ax.axvline(1.0, color=GRID, lw=1, zorder=1)
        ax.grid(axis="x", color=GRID, lw=0.8, zorder=0)
        ax.set_axisbelow(True)
        for s in ("top", "right", "left"):
            ax.spines[s].set_visible(False)
        ax.tick_params(length=0)
    fig.suptitle("champHLA HLA-A2 accuracy by method  (strict A*02, Wilson 95% CI)",
                 fontsize=12, fontweight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig(out, dpi=200); plt.close(fig)


# ── Fig 3: donor x method matrix ─────────────────────────────────────────────────
def fig_donor_method(donors, out):
    donors = sorted(donors, key=lambda r: (r["flow_a2"] != "yes", r["donor"]))
    rows = [r["donor"] for r in donors]
    ncol = len(DETAIL_COLS) + 1  # +1 flow truth column
    fig, ax = plt.subplots(figsize=(9.6, 0.42 * len(rows) + 2.2))
    for i, r in enumerate(donors):
        yv = len(rows) - 1 - i
        # flow truth column (x = -1)
        pos = r["flow_a2"] == "yes"
        ax.add_patch(plt.Rectangle((-1, yv), 1, 1, facecolor=INK2 if pos else "#d9d8d2",
                                   edgecolor=SURF, lw=2))
        ax.text(-0.5, yv + 0.5, "+" if pos else "−", ha="center", va="center",
                color="#ffffff" if pos else INK2, fontweight="bold", fontsize=11)
        for j, (stem, _) in enumerate(DETAIL_COLS):
            pred = (r.get(f"{stem}_strict") or "").strip()
            if pred == "":
                face, sym, sc = GREY, "·", MUTED
            else:
                pred_pos = pred == "yes"
                correct = pred == r["flow_a2"]
                if correct:
                    face = GOOD if pred_pos else GOOD_L
                    sc = "#ffffff" if pred_pos else INK
                else:
                    face = CRIT; sc = "#ffffff"
                sym = "+" if pred_pos else "−"
            ax.add_patch(plt.Rectangle((j, yv), 1, 1, facecolor=face,
                                       edgecolor=SURF, lw=2))
            ax.text(j + 0.5, yv + 0.5, sym, ha="center", va="center",
                    color=sc, fontweight="bold", fontsize=11)
    ax.set_xlim(-1, len(DETAIL_COLS)); ax.set_ylim(0, len(rows))
    ax.set_xticks([-0.5] + [j + 0.5 for j in range(len(DETAIL_COLS))])
    ax.set_xticklabels(["Flow\ntruth"] + [d for _, d in DETAIL_COLS], rotation=45,
                       ha="right", fontsize=8.5)
    ax.set_yticks([len(rows) - 1 - i + 0.5 for i in range(len(rows))])
    ax.set_yticklabels(rows, fontsize=9)
    ax.tick_params(length=0)
    for s in ax.spines.values():
        s.set_visible(False)
    ax.xaxis.set_ticks_position("top"); ax.xaxis.set_label_position("top")
    legend = [
        Patch(facecolor=GOOD, label="Correct, A2+ (TP)"),
        Patch(facecolor=GOOD_L, label="Correct, A2− (TN)"),
        Patch(facecolor=CRIT, label="Discordant (FP/FN)"),
        Patch(facecolor=GREY, label="No call"),
    ]
    ax.legend(handles=legend, loc="upper center", bbox_to_anchor=(0.5, -0.02),
              ncol=4, frameon=False, fontsize=8.5, handlelength=1.1)
    ax.set_title("Per-donor HLA-A2 call by champHLA method vs flow cytometry (strict A*02)",
                 fontsize=11.5, pad=34)
    fig.tight_layout()
    fig.savefig(out, dpi=200, bbox_inches="tight"); plt.close(fig)


# ── Fig 4: specificity strict vs cross ───────────────────────────────────────────
def fig_mapping_effect(metrics, wilson, out):
    labels = [d for _, d in METHOD_ORDER]
    y = list(range(len(labels)))[::-1]
    fig, ax = plt.subplots(figsize=(7.4, 5.0))
    h = 0.34
    for off, mapping, color, lab in ((h/2 + .02, "strict", AQUA, "strict A*02"),
                                     (-h/2 - .02, "cross", "#e08a2b", "+ A*68/A*69")):
        for yi, (scope, _) in zip(y, METHOD_ORDER):
            m = metrics.get((scope, mapping))
            if not m:
                continue
            val, lo, hi, n = _metric_ci(m, "spec", wilson)
            ax.barh(yi + off, val, height=h, color=color, edgecolor=SURF, lw=1, zorder=2,
                    label=lab if yi == y[0] else None)
            ax.text(val + 0.02, yi + off, f"{val:.2f}", va="center", ha="left",
                    fontsize=8, color=INK)
    ax.set_xlim(0, 1.2); ax.set_xticks([0, .25, .5, .75, 1.0])
    ax.set_yticks(y); ax.set_yticklabels(labels, fontsize=9)
    ax.grid(axis="x", color=GRID, lw=0.8); ax.set_axisbelow(True)
    for s in ("top", "right", "left"):
        ax.spines[s].set_visible(False)
    ax.tick_params(length=0)
    ax.legend(frameon=False, fontsize=9, loc="upper center",
              bbox_to_anchor=(0.5, -0.06), ncol=2)
    ax.set_title("Specificity drops when A*68/A*69 are counted as A2\n"
                 "(BB7.2 cross-reactivity) — champHLA methods",
                 fontsize=11, fontweight="bold")
    fig.tight_layout(rect=[0, 0.04, 1, 1])
    fig.savefig(out, dpi=200); plt.close(fig)


# ── Fig 5: champHLA cross-tool concordance (donor x locus) ───────────────────────
WARN = "#fab219"   # status: minor sub-typing difference


def fig_champhla_concordance(path_by_locus, out):
    if not os.path.exists(path_by_locus):
        return
    rows = load_donors(path_by_locus)
    donors = sorted({r["donor"] for r in rows},
                    key=lambda d: (next(x["flow_A2"] for x in rows
                                        if x["donor"] == d and x["locus"] == "A") != "yes", d))
    loci = ["A", "B", "C"]
    cell = {(r["donor"], r["locus"]): r for r in rows}

    def state(r):
        if r is None or int(r["n_tools"]) == 0:
            return GREY, "·", MUTED, "No call"
        if r["within_donor_conflict"] == "SWAP?":
            return CRIT, "✗", "#ffffff", "Swap / irreconcilable"
        if int(r["n_tools"]) == 1:
            return "#cfd8e6", "1", INK2, "Single tool (n/a)"
        if r["unanimous"] == "yes":
            return GOOD, "=", "#ffffff", "Unanimous (2-field)"
        return WARN, "~", INK, "Sub-typing difference"

    fig, ax = plt.subplots(figsize=(6.4, 0.42 * len(donors) + 2.0))
    seen = {}
    for i, d in enumerate(donors):
        yv = len(donors) - 1 - i
        flow = cell[(d, "A")]["flow_A2"]
        pos = flow == "yes"
        ax.add_patch(plt.Rectangle((-1, yv), 1, 1, facecolor=INK2 if pos else "#d9d8d2",
                                   edgecolor=SURF, lw=2))
        ax.text(-0.5, yv + 0.5, "+" if pos else "−", ha="center", va="center",
                color="#ffffff" if pos else INK2, fontweight="bold", fontsize=11)
        for j, g in enumerate(loci):
            face, sym, sc, lab = state(cell.get((d, g)))
            seen[lab] = face
            ax.add_patch(plt.Rectangle((j, yv), 1, 1, facecolor=face, edgecolor=SURF, lw=2))
            ax.text(j + 0.5, yv + 0.5, sym, ha="center", va="center", color=sc,
                    fontweight="bold", fontsize=11)
    ax.set_xlim(-1, len(loci)); ax.set_ylim(0, len(donors))
    ax.set_xticks([-0.5] + [j + 0.5 for j in range(len(loci))])
    ax.set_xticklabels(["Flow\nA2"] + [f"HLA-{g}" for g in loci], fontsize=9.5)
    ax.set_yticks([len(donors) - 1 - i + 0.5 for i in range(len(donors))])
    ax.set_yticklabels(donors, fontsize=9)
    ax.tick_params(length=0)
    for s in ax.spines.values():
        s.set_visible(False)
    ax.xaxis.set_ticks_position("top"); ax.xaxis.set_label_position("top")
    order = ["Unanimous (2-field)", "Sub-typing difference", "Swap / irreconcilable",
             "Single tool (n/a)", "No call"]
    handles = [Patch(facecolor=seen[l], label=l) for l in order if l in seen]
    ax.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, -0.02),
              ncol=3, frameon=False, fontsize=8.3, handlelength=1.1)
    ax.set_title("champHLA cross-tool genotype concordance (2-field)\n"
                 "per donor × locus, flow-cohort samples", fontsize=11, pad=30)
    fig.tight_layout()
    fig.savefig(out, dpi=200, bbox_inches="tight"); plt.close(fig)


# ── Figure 14 (manuscript): confusion + sens/spec by method ──────────────────────
def fig_manuscript_14(metrics, wilson, outdir, stem="figure_14_flow_hla_a2"):
    import matplotlib.gridspec as gridspec
    os.makedirs(outdir, exist_ok=True)
    fig = plt.figure(figsize=(11.0, 5.0))
    gs = gridspec.GridSpec(1, 2, width_ratios=[1.0, 1.9], wspace=0.28)

    # Panel A — consensus confusion matrix (strict)
    axA = fig.add_subplot(gs[0, 0])
    m = metrics[("consensus", "strict")]
    tp, fp, tn, fn = int(m["TP"]), int(m["FP"]), int(m["TN"]), int(m["FN"])
    cells = [[(tp, "TP"), (fn, "FN")], [(fp, "FP"), (tn, "TN")]]
    for i in range(2):
        for j in range(2):
            val, lab = cells[i][j]
            correct = lab in ("TP", "TN")
            if correct:
                face = GOOD if lab == "TP" else GOOD_L
                tc = "#ffffff" if lab == "TP" else INK
            elif val > 0:
                face, tc = CRIT, "#ffffff"
            else:
                face, tc = GREY, MUTED
            axA.add_patch(plt.Rectangle((j, 1 - i), 1, 1, facecolor=face,
                                        edgecolor=SURF, linewidth=3, zorder=1))
            axA.text(j + 0.5, 1 - i + 0.58, str(val), ha="center", va="center",
                     fontsize=20, fontweight="bold", color=tc, zorder=2)
            axA.text(j + 0.5, 1 - i + 0.26, lab, ha="center", va="center",
                     fontsize=8.5, color=tc, zorder=2)
    axA.set_xlim(0, 2); axA.set_ylim(0, 2)
    axA.set_xticks([0.5, 1.5]); axA.set_xticklabels(["A2 +", "A2 −"])
    axA.set_yticks([1.5, 0.5]); axA.set_yticklabels(["A2 +", "A2 −"])
    axA.set_xlabel("ChampHLA consensus", color=INK2)
    axA.set_ylabel("Flow cytometry (truth)", color=INK2)
    axA.set_title(f"A  Consensus vs flow (strict A*02)\nN={m['N']} · accuracy "
                  f"{float(m['accuracy']):.2f}", fontsize=10.5, loc="left")
    for s in axA.spines.values():
        s.set_visible(False)
    axA.tick_params(length=0)

    # Panel B — sensitivity & specificity by method (strict), grouped bars + CI
    axB = fig.add_subplot(gs[0, 1])
    labels = [d for _, d in METHOD_ORDER]
    y = list(range(len(labels)))[::-1]
    h = 0.36
    for off, kind, color, lab in ((h / 2 + .02, "sens", BLUE, "Sensitivity"),
                                  (-h / 2 - .02, "spec", AQUA, "Specificity")):
        for yi, (scope, _) in zip(y, METHOD_ORDER):
            mm = metrics.get((scope, "strict"))
            if not mm:
                continue
            val, lo, hi, n = _metric_ci(mm, kind, wilson)
            axB.barh(yi + off, val, height=h, color=color, edgecolor=SURF, lw=1,
                     zorder=2, label=lab if yi == y[0] else None)
            axB.plot([lo, hi], [yi + off, yi + off], color=INK2, lw=1.2, zorder=3)
            axB.text(max(val, hi) + 0.03, yi + off, f"{val:.2f}", va="center",
                     ha="left", fontsize=7.6, color=INK)
    axB.set_xlim(0, 1.24); axB.set_xticks([0, .25, .5, .75, 1.0])
    axB.set_yticks(y); axB.set_yticklabels(labels, fontsize=9)
    axB.grid(axis="x", color=GRID, lw=0.8); axB.set_axisbelow(True)
    for s in ("top", "right", "left"):
        axB.spines[s].set_visible(False)
    axB.tick_params(length=0)
    axB.legend(frameon=False, fontsize=8.5, loc="upper center",
               bbox_to_anchor=(0.5, -0.06), ncol=2)
    axB.set_title("B  Sensitivity & specificity by method (strict A*02)",
                  fontsize=10.5, loc="left")

    fig.tight_layout(rect=[0, 0.04, 1, 1])
    paths = []
    for ext in ("png", "pdf", "svg"):
        p = os.path.join(outdir, f"{stem}.{ext}")
        fig.savefig(p, dpi=200)
        paths.append(p)
    plt.close(fig)
    return paths


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--indir", required=True)
    ap.add_argument("--benchmark",
                    default="/scratch/project_2008084/pihla-publish/bin/hla_benchmark.py")
    ap.add_argument("--outdir", default=None)
    ap.add_argument("--manuscript-figdir", default=None,
                    help="if set, also write figure_14_flow_hla_a2.{png,pdf,svg} here")
    args = ap.parse_args()
    outdir = args.outdir or os.path.join(args.indir, "figures")
    os.makedirs(outdir, exist_ok=True)

    metrics = load_metrics(os.path.join(args.indir, "hla_a2_flow_metrics.tsv"))
    donors = load_donors(os.path.join(args.indir, "hla_a2_flow_per_donor.tsv"))
    wilson = load_wilson(args.benchmark)

    fig_confusion(metrics, os.path.join(outdir, "fig1_confusion.png"))
    fig_sens_spec(metrics, wilson, os.path.join(outdir, "fig2_sens_spec.png"))
    fig_donor_method(donors, os.path.join(outdir, "fig3_donor_method.png"))
    fig_mapping_effect(metrics, wilson, os.path.join(outdir, "fig4_mapping_effect.png"))
    fig_champhla_concordance(os.path.join(args.indir, "champhla_detail_by_locus.tsv"),
                             os.path.join(outdir, "fig5_champhla_concordance.png"))
    print("Figures written to", outdir)
    for f in ("fig1_confusion", "fig2_sens_spec", "fig3_donor_method",
              "fig4_mapping_effect", "fig5_champhla_concordance"):
        print("  ", os.path.join(outdir, f + ".png"))
    if args.manuscript_figdir:
        for p in fig_manuscript_14(metrics, wilson, args.manuscript_figdir):
            print("  ", p)


if __name__ == "__main__":
    main()
