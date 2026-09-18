#!/usr/bin/env python3.11
"""CIWD-stratified concordance — champHLA (MVHLA) manuscript auxiliary figure.

Reads `summary_ciwd_stratified.tsv` (written by bin/hla_benchmark.py) and plots two-field
concordance **by CIWD 3.0.0 allele-commonness stratum** (common → intermediate → well-documented →
not-CIWD → unknown/novel) for the ChampHLA consensus vs the baseline consensus methods, one panel per
modality. This isolates the reviewer question "does the consensus win only on easy common alleles, or
also on rare / well-documented ones?"

CIWD is an allele-classification catalogue, NOT a ground-truth cohort (assets/README_CIWD.md); this
is an auxiliary view over the existing ground-truth benchmark, not a new truth source.

Style (Arial, Wong palette, 300 DPI) comes from figure_hub/plot_style.py.
"""
import argparse
from collections import defaultdict
from pathlib import Path
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

SCRIPT_DIR = Path(__file__).resolve().parent
HUB_ROOT = SCRIPT_DIR.parent
sys.path.insert(0, str(HUB_ROOT))
from plot_style import apply_style, save_fig, legend_outside, MOD_COLORS, WONG  # noqa: E402

apply_style()

STRATA = ["common", "intermediate", "well_documented", "not_ciwd", "unknown"]
STRATA_LABEL = {
    "common": "Common", "intermediate": "Intermediate", "well_documented": "Well-doc.",
    "not_ciwd": "Not-CIWD", "unknown": "Novel/\nunseen",
}
# Methods to display (label -> colour); ChampHLA highlighted.
METHOD_STYLE = {
    "ChampionChallenger": ("ChampHLA (CC)", WONG["vermil"]),
    "WeightedConsensus": ("Weighted", WONG["blue"]),
    "MajorityVote": ("Majority", WONG["sky"]),
}


def read_tsv(path):
    rows = []
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        for ln in fh:
            rows.append(dict(zip(header, ln.rstrip("\n").split("\t"))))
    return rows


def load(path):
    """(modality, method, stratum) -> (concordance, n)."""
    data = {}
    modalities = []
    for r in read_tsv(path):
        mod = r["modality"]
        if mod not in modalities:
            modalities.append(mod)
        try:
            conc = float(r["overall_correct_call_rate"])
        except (ValueError, KeyError):
            conc = float("nan")
        data[(mod, r["method"], r["ciwd_stratum"])] = (conc, int(r["gene_rows"]))
    return data, modalities


def panel(ax, data, modality, methods):
    present = [s for s in STRATA if any((modality, m, s) in data for m in methods)]
    if not present:
        present = ["common"]
    x = np.arange(len(present))
    width = 0.8 / max(len(methods), 1)
    for i, method in enumerate(methods):
        label, color = METHOD_STYLE[method]
        vals = [data.get((modality, method, s), (float("nan"), 0))[0] for s in present]
        bars = ax.bar(x + (i - (len(methods) - 1) / 2) * width, vals, width,
                      label=label, color=color, edgecolor="white", linewidth=0.4)
        for rect, s in zip(bars, present):
            n = data.get((modality, method, s), (0, 0))[1]
            if n:
                ax.text(rect.get_x() + rect.get_width() / 2, 0.02, f"n={n}",
                        ha="center", va="bottom", fontsize=5.2, rotation=90, color="#333333")
    ax.set_xticks(x)
    ax.set_xticklabels([STRATA_LABEL.get(s, s) for s in present], fontsize=7)
    ax.set_ylim(0, 1.02)
    ax.set_ylabel("Two-field concordance")
    ax.set_title(MOD_COLORS.get(modality) and modality.upper() or modality.upper(),
                 color=MOD_COLORS.get(modality, WONG["black"]), fontsize=9)
    ax.axhline(1.0, color="#cccccc", lw=0.6, ls=":")


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--summary", required=True,
                   help="summary_ciwd_stratified.tsv from a benchmark run")
    p.add_argument("--out-dir", type=Path,
                   default=Path("/scratch/project_2008084/pihla-publish/analysis/figures_final"))
    p.add_argument("--stem", default="figure_ciwd_stratified")
    args = p.parse_args()

    data, modalities = load(args.summary)
    methods = [m for m in METHOD_STYLE if any(k[1] == m for k in data)]
    if not methods:
        raise SystemExit("No consensus methods found in %s" % args.summary)

    fig, axes = plt.subplots(1, len(modalities), figsize=(4.6 * len(modalities), 4.2), squeeze=False)
    for ax, modality in zip(axes[0], modalities):
        panel(ax, data, modality, methods)
    legend_outside(axes[0][-1], loc="lower center", ncol=1)
    fig.suptitle("Concordance by CIWD 3.0.0 allele-commonness stratum", fontsize=10)
    fig.tight_layout(rect=(0, 0.04, 1, 0.96))
    save_fig(fig, args.out_dir / args.stem)
    print(f"Wrote {args.out_dir / args.stem}.{{pdf,svg,png}}")


if __name__ == "__main__":
    main()
