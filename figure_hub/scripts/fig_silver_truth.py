#!/usr/bin/env python3.11
"""Figure 12 — Orthogonal silver-standard truth validation (HPRC vs IMGT-allele DB).

Two panels for the ChampHLA (MVHLA) manuscript §Results-6:
  (A) Per-locus + overall allele concordance vs 1000G gold (HPRC pangenome DB vs
      IPD-IMGT/HLA-allele DB), 30 samples, two-field; pilot overall annotated.
  (B) Calibration: concordance vs Locityper genotype-quality threshold (HPRC rises
      = calibrated; IMGT flat/declining = miscalibrated).

Real numbers are read from the reproducible summaries produced by
truth_bio_run/make_summaries.py:
  concordance_summary.tsv  (db, sample_set, profile, locus, n_alleles, conc_2field, conc_ggroup)
  calibration_summary.tsv  (db, gq_threshold, n_loci, frac_kept, conc_2field)
If a file is missing, documented constants from truth_bio_run/RESULTS.md are used.
Style (Arial, Wong palette, 300 DPI, legends outside axes) comes from figure_hub/plot_style.py.
"""
import argparse
from pathlib import Path
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

SCRIPT_DIR = Path(__file__).resolve().parent
HUB_ROOT = SCRIPT_DIR.parent
sys.path.insert(0, str(HUB_ROOT))
from plot_style import apply_style, save_fig, legend_outside, WONG  # noqa: E402

apply_style()

DB_COLOR = {"HPRC": WONG["blue"], "IMGT": WONG["vermil"]}
DB_LABEL = {"HPRC": "HPRC pangenome DB", "IMGT": "IPD-IMGT/HLA-allele DB"}

# Documented fallback constants (truth_bio_run/RESULTS.md) ----------------------
FALLBACK_CONC = {  # (db, sample_set) -> {locus: 2field}
    ("HPRC", "30samples"): {"HLA-A": .733, "HLA-B": .767, "HLA-C": .933, "ALL": .811},
    ("IMGT", "30samples"): {"HLA-A": .333, "HLA-B": .550, "HLA-C": .500, "ALL": .461},
}
FALLBACK_PILOT = {"HPRC": .944, "IMGT": .556}
FALLBACK_CAL = {  # db -> list of (gq, conc)
    "HPRC": [(0, .811), (3, .818), (5, .821), (10, .862), (15, .875), (20, .905)],
    "IMGT": [(0, .461), (3, .464), (5, .460), (10, .432), (15, .364), (20, .500)],
}


def read_tsv(path):
    if not Path(path).exists():
        return None
    rows = []
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        for ln in fh:
            rows.append(dict(zip(header, ln.rstrip("\n").split("\t"))))
    return rows


def load_concordance(path):
    rows = read_tsv(path)
    if not rows:
        return FALLBACK_CONC, FALLBACK_PILOT
    conc, pilot = {}, {}
    for r in rows:
        key = (r["db"], r["sample_set"])
        conc.setdefault(key, {})[r["locus"]] = float(r["conc_2field"])
        if r["sample_set"] == "3pilots" and r["locus"] == "ALL":
            pilot[r["db"]] = float(r["conc_2field"])
    return conc, (pilot or FALLBACK_PILOT)


def load_calibration(path):
    rows = read_tsv(path)
    if not rows:
        return FALLBACK_CAL
    cal = {}
    for r in rows:
        cal.setdefault(r["db"], []).append((float(r["gq_threshold"]), float(r["conc_2field"])))
    for db in cal:
        cal[db].sort()
    return cal


def panel_a(ax, conc, pilot):
    loci = ["HLA-A", "HLA-B", "HLA-C", "ALL"]
    labels = ["HLA-A", "HLA-B", "HLA-C", "Overall"]
    x = np.arange(len(loci))
    w = 0.38
    for i, db in enumerate(["HPRC", "IMGT"]):
        vals = [conc[(db, "30samples")][loc] for loc in loci]
        bars = ax.bar(x + (i - 0.5) * w, vals, w, color=DB_COLOR[db],
                      label=DB_LABEL[db], edgecolor="black", linewidth=0.4)
        for b, v in zip(bars, vals):
            ax.text(b.get_x() + b.get_width() / 2, v + 0.015, f"{v:.2f}",
                    ha="center", va="bottom", fontsize=6.5)
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.set_ylabel("Two-field concordance vs gold")
    ax.set_ylim(0, 1.08)
    ax.axhline(1.0, color="grey", lw=0.5, ls=":")
    ax.set_title("A  30-sample concordance (shared depth profile)", loc="left", fontsize=9)
    # pilot annotation
    ax.text(0.02, 0.97,
            f"3-sample per-sample-depth pilot (overall):\n"
            f"HPRC {pilot['HPRC']:.3f}   IMGT {pilot['IMGT']:.3f}",
            transform=ax.transAxes, va="top", ha="left", fontsize=6.8,
            bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="grey", lw=0.5))
    legend_outside(ax, loc="lower center", ncol=2)


def panel_b(ax, cal):
    for db in ["HPRC", "IMGT"]:
        xs = [g for g, _ in cal[db]]
        ys = [c for _, c in cal[db]]
        ax.plot(xs, ys, marker="o", ms=4, lw=1.6, color=DB_COLOR[db], label=DB_LABEL[db])
    ax.set_xlabel("Locityper genotype-quality threshold (GQ ≥)")
    ax.set_ylabel("Concordance among retained calls")
    ax.set_ylim(0.3, 1.0)
    ax.set_title("B  Confidence calibration (30 samples)", loc="left", fontsize=9)
    ax.annotate("calibrated\n(rises with GQ)", xy=(15, 0.875), xytext=(8.5, 0.96),
                fontsize=6.8, color=DB_COLOR["HPRC"],
                arrowprops=dict(arrowstyle="->", color=DB_COLOR["HPRC"], lw=0.8))
    ax.annotate("miscalibrated\n(flat / declining)", xy=(12, 0.43), xytext=(2.0, 0.34),
                fontsize=6.8, color=DB_COLOR["IMGT"],
                arrowprops=dict(arrowstyle="->", color=DB_COLOR["IMGT"], lw=0.8))
    legend_outside(ax, loc="lower center", ncol=1)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--conc", default="/scratch/project_2008084/truth_bio_run/concordance_summary.tsv")
    p.add_argument("--cal", default="/scratch/project_2008084/truth_bio_run/calibration_summary.tsv")
    p.add_argument("--out-dir", type=Path,
                   default=Path("/scratch/project_2008084/pihla-publish/analysis/figures_final"))
    p.add_argument("--stem", default="figure_12_silver_truth_hprc")
    args = p.parse_args()

    conc, pilot = load_concordance(args.conc)
    cal = load_calibration(args.cal)

    fig, axes = plt.subplots(1, 2, figsize=(9.6, 4.3))
    panel_a(axes[0], conc, pilot)
    panel_b(axes[1], cal)
    fig.tight_layout(rect=(0, 0.06, 1, 1))
    save_fig(fig, args.out_dir / args.stem)
    print(f"Wrote {args.out_dir / args.stem}.{{pdf,svg,png}}")


if __name__ == "__main__":
    main()
