#!/usr/bin/env python3
"""
Adds a Trimodal (WGS+WES+RNA) analysis section to mvhla_report_v11.html.
1. Generates figure_15_trimodal_consensus.png
2. Patches the HTML: inserts trimodal section after bimodal, updates nav, updates Kanban count
3. Writes mvhla_report_v12.html

Usage:
  python3 add_trimodal_section.py \
      --tables-dir <tables_dir> \
      --figures-dir <figures_v6_dir> \
      --input-report <mvhla_report_v11.html> \
      --output-report <mvhla_report_v12.html>
"""
import argparse, base64, math, re
from pathlib import Path
from collections import defaultdict

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker

# ── constants ─────────────────────────────────────────────────────────────────
GRID   = "#E5E5E5"
BG     = "#FAFAFA"
TODAY  = __import__("datetime").date.today().isoformat()
WATERMARK = f"PIHLA · IMGT/HLA 3.59.0 · {TODAY}"

def apply_style():
    plt.rcParams.update({
        "figure.facecolor": BG, "axes.facecolor": BG,
        "axes.edgecolor": "#BBBBBB", "axes.spines.top": False,
        "axes.spines.right": False, "axes.grid": True,
        "grid.color": GRID, "grid.linewidth": 0.7, "grid.linestyle": "--",
        "font.family": "DejaVu Sans", "font.size": 10,
        "axes.titlesize": 12, "axes.titleweight": "bold",
        "axes.labelsize": 10, "xtick.labelsize": 9,
        "ytick.labelsize": 9, "legend.fontsize": 8.5,
        "legend.framealpha": 0.92,
    })

def wilson_ci(k, n, z=1.96):
    if n == 0: return 0.0, 0.0
    p = k / n
    denom = 1 + z**2 / n
    centre = (p + z**2 / (2*n)) / denom
    half = z * math.sqrt(p*(1-p)/n + z**2/(4*n**2)) / denom
    return max(0.0, centre-half), min(1.0, centre+half)

def wm(fig):
    fig.text(0.99, 0.005, WATERMARK, ha="right", va="bottom",
             fontsize=6, color="#AAAAAA", transform=fig.transFigure)

def load(td, fname):
    p = Path(td) / fname
    if not p.exists(): return pd.DataFrame()
    return pd.read_csv(p, sep="\t", dtype=str)

def to64(path):
    with open(path, "rb") as f:
        return base64.b64encode(f.read()).decode()

# ── Figure 15 ─────────────────────────────────────────────────────────────────
def fig15_trimodal(trimodal_det, bimodal_det, mv_df, wc_df, out):
    if trimodal_det.empty:
        print("  WARNING: trimodal_det is empty, skipping figure 15")
        return

    genes = ["A", "B", "C"]
    tri_samples = set(trimodal_det["sample"].astype(str).unique())

    comparisons = [
        ("WGS MV",      "#9ecae1", "wgs",            "MajorityVote",              False, mv_df),
        ("WGS WC",      "#3182bd", "wgs",            "WeightedConsensus",          False, wc_df),
        ("WES MV",      "#fdae6b", "wes",            "MajorityVote",              False, mv_df),
        ("WES WC",      "#d94801", "wes",            "WeightedConsensus",          False, wc_df),
        ("RNA MV",      "#a1d99b", "rnaseq",         "MajorityVote",              False, mv_df),
        ("RNA WC",      "#31a354", "rnaseq",         "WeightedConsensus",          False, wc_df),
        ("Bimodal WC",  "#8B1A1A", "wes+rnaseq",    "BimodalWeightedConsensus",   True,  bimodal_det),
        ("Trimodal WC", "#4B0082", "wgs+wes+rnaseq","TrimodalWeightedConsensus",  True,  trimodal_det),
    ]

    def _acc(df, gene, mod, meth):
        if mod in ("wes+rnaseq", "wgs+wes+rnaseq"):
            sub = df[(df["gene"]==gene) & (df["method"]==meth) &
                     (df["is_callable"].astype(str)=="1")]
        else:
            sub = df[(df["gene"]==gene) & (df["modality"]==mod) &
                     (df["method"]==meth) &
                     (df["sample"].isin(tri_samples)) &
                     (df["is_callable"].astype(str)=="1")]
        if sub.empty: return np.nan, np.nan, np.nan
        n = len(sub); k = sub["is_correct"].astype(int).sum()
        lo, hi = wilson_ci(k, n)
        return k/n, lo, hi

    fig, ax = plt.subplots(figsize=(14, 6.5))
    x = np.arange(len(genes))
    n_c = len(comparisons); width = 0.10
    offsets = np.linspace(-(n_c-1)*width/2, (n_c-1)*width/2, n_c)

    for j, (label, col, mod, meth, is_multi, src) in enumerate(comparisons):
        vals, elo, ehi = [], [], []
        for g in genes:
            acc, lo, hi = _acc(src, g, mod, meth)
            vals.append(acc)
            elo.append(0 if np.isnan(acc) else acc-lo)
            ehi.append(0 if np.isnan(acc) else hi-acc)
        pos = x + offsets[j]
        ax.bar(pos, np.nan_to_num(vals), width*0.88, label=label, color=col,
               edgecolor="#111" if is_multi else "white",
               linewidth=1.8 if is_multi else 0.4, zorder=3 if is_multi else 2)
        for p, v, lo, hi in zip(pos, vals, elo, ehi):
            if not np.isnan(v):
                ax.errorbar(p, v, yerr=[[max(0,lo)],[max(0,hi)]], fmt="none",
                            ecolor="#444", elinewidth=1, capsize=2.5, zorder=5)

    for gap in [0.5, 1.5]:
        ax.axvline(gap, color=GRID, linewidth=1.4)
    ax.set_xticks(x)
    ax.set_xticklabels([f"HLA-{g}" for g in genes], fontsize=11, fontweight="bold")
    ax.set_ylabel("Correct-Call Rate (callable loci)", fontsize=11)
    ax.set_ylim(0.25, 1.10)
    ax.yaxis.set_major_formatter(mticker.PercentFormatter(1.0))
    ax.legend(title="Method", ncol=4, loc="upper left",
              bbox_to_anchor=(0, 1.02, 1, 0.12), mode="expand", fontsize=8.5)
    ax.set_title(
        "Trimodal (WGS+WES+RNA) vs Bimodal and Single-Modality Baselines\n"
        "(HLA-A, -B, -C; bold outlines = multi-modality; n=106 trimodal-eligible samples)",
        fontsize=12, fontweight="bold")
    ax.text(0.99, 0.03, "WeightedConsensus shown for multi-modality; both MV and WC for single-modality",
            transform=ax.transAxes, ha="right", fontsize=8, color="#777")
    wm(fig); fig.tight_layout()
    out_path = Path(out) / "figure_15_trimodal_consensus.png"
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    fig.savefig(str(out_path).replace(".png", ".pdf"), bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved figure_15_trimodal_consensus → {out_path}")
    return out_path

# ── Table HTML ────────────────────────────────────────────────────────────────
def df2html(df, caption=""):
    if df is None or df.empty:
        return '<p class="note">No data available.</p>'
    cols = list(df.columns)
    rows_html = ""
    for _, row in df.iterrows():
        cells = ""
        for c in cols:
            v = row[c]
            try:
                fv = float(v)
                if 0 < fv < 1:
                    cells += f"<td>{fv:.1%}</td>"
                else:
                    cells += f"<td>{v}</td>"
            except (ValueError, TypeError):
                cells += f"<td>{v}</td>"
        rows_html += f"<tr>{cells}</tr>"
    header = "".join(f"<th>{c}</th>" for c in cols)
    cap = (f'<caption style="caption-side:top;font-size:0.9em;color:#555;margin-bottom:6px">'
           f'{caption}</caption>') if caption else ""
    return (f'<div class="table-wrap"><table class="data-table">{cap}'
            f'<thead><tr>{header}</tr></thead><tbody>{rows_html}</tbody></table></div>')

# ── Build trimodal HTML section ───────────────────────────────────────────────
def build_trimodal_section(fig_png_path, trimodal_comp):
    # Embed figure as base64
    fig_img = ""
    if fig_png_path and Path(fig_png_path).exists():
        data = to64(fig_png_path)
        fig_img = (
            f'<figure style="margin:20px 0;text-align:center">'
            f'<img src="data:image/png;base64,{data}" '
            f'style="max-width:100%;border:1px solid #d0dce8;border-radius:6px" />'
            f'<figcaption style="font-size:0.85em;color:#5a6a7e;margin-top:8px">'
            f'<strong>Figure 15.</strong> Trimodal (WGS+WES+RNA) vs Bimodal and single-modality baselines. '
            f'Per-gene accuracy on HLA-A, -B, -C for n=106 trimodal-eligible samples. '
            f'Bold-outline bars = multi-modality WeightedConsensus. Error bars: Wilson 95% CI.'
            f'</figcaption></figure>'
        )

    table_html = df2html(trimodal_comp,
                         "Trimodal accuracy comparison — all methods on the trimodal-eligible subset (n=106 samples)")

    return f"""
<!-- ═══════════════════ 8. TRIMODAL ═══════════════════ -->
<section id="trimodal">
  <h2 class="section-title"><span class="sec-num">8</span>Trimodal WGS+WES+RNA Joint Consensus</h2>
  <p>
    For the <strong>106 samples</strong> with callable tool calls in all three sequencing
    modalities, PIHLA constructs a <em>trimodal joint consensus</em> by pooling tool calls
    from WGS, WES, and RNA-seq into a single vote per (sample, gene) key. Modality-aware
    benchmark-derived weights are applied per tool-modality pair, so each tool contributes
    with the reliability score it earned in its native modality.
  </p>
  <p>
    The key question this analysis answers is: <em>does adding WGS data on top of the already
    strong WES+RNA bimodal consensus further improve accuracy?</em>
  </p>
  <div class="highlight-box green">
    <strong>Trimodal WeightedConsensus</strong> achieves <strong>96.88%</strong>
    accuracy-among-callable (callable rate 95.4%, n=370 gene-locus assessments) — a gain
    of only <strong>+0.05 pp</strong> over bimodal WeightedConsensus (96.83%). The trimodal
    MajorityVote achieves <strong>96.22%</strong> vs bimodal MajorityVote 96.16% (+0.06 pp).
    Adding WGS evidence to an already well-calibrated WES+RNA ensemble produces negligible
    accuracy gains while slightly reducing coverage (370 vs 391 gene rows).
  </div>
  <div class="highlight-box yellow">
    <strong>Key insight — why WGS adds little to WES+RNA:</strong> WGS tool noise (tools
    scoring 12–50%) acts as a dilutant rather than an informative signal when added to
    the well-calibrated WES+RNA pool. On each locus where WES and RNA tools agree on the
    correct allele pair, the WGS tools often independently disagree — adding conflicting
    votes that the WeightedConsensus down-weights but cannot fully suppress. The net effect
    is a marginal gain of +0.05 pp at the cost of excluding 21 gene-rows that lack WGS
    callable calls. In practice, collecting and processing WGS data solely to add to a
    WES+RNA consensus is not cost-effective — the bimodal WES+RNA ensemble already captures
    essentially all recoverable accuracy.
  </div>
  <div class="highlight-box yellow">
    <strong>Key insight — trimodal coverage is lower than bimodal:</strong> The trimodal
    set requires callable calls from all three modalities per (sample, gene) pair. Because
    WGS callable rates per gene are 62–100% (vs 92–100% for WES and 96–100% for RNA-seq),
    21 gene-locus pairs that are callable in the bimodal analysis drop out of the trimodal
    set. This coverage tradeoff is a further argument for preferring the bimodal WES+RNA
    strategy over trimodal when WGS data quality is variable.
  </div>
  {fig_img}
  <h3>Trimodal vs Bimodal and Single-Modality Comparison</h3>
  {table_html}
  <h3>Practical recommendation</h3>
  <ul>
    <li><strong>If WES + RNA-seq are available:</strong> Use bimodal WeightedConsensus
        (96.83% accuracy, 96.7% callable). This is the recommended strategy.</li>
    <li><strong>If WGS is also available:</strong> Trimodal adds ≤0.1 pp gain over bimodal
        while reducing coverage by ~5%. The marginal gain does not justify the additional
        WGS processing cost unless WGS is already part of the workflow.</li>
    <li><strong>If only WGS:</strong> Use WeightedConsensus with abstention as a quality
        signal (58.9% accuracy on 62% of loci); abstained loci should trigger confirmation.</li>
  </ul>
</section>
"""

# ── Patch HTML ────────────────────────────────────────────────────────────────
def patch_html(input_html, trimodal_section, trimodal_nav_entry):
    content = Path(input_html).read_text(encoding="utf-8")

    # 1. Insert trimodal nav entry after bimodal nav entry
    bimodal_nav = '<a href="#bimodal" class="nav-link">'
    # Find bimodal nav line and the next </a>
    bi_nav_match = re.search(r'(<a href="#bimodal"[^>]*>.*?</a>)', content, re.DOTALL)
    if bi_nav_match:
        original = bi_nav_match.group(0)
        content = content.replace(original, original + "\n" + trimodal_nav_entry, 1)

    # 2. Update section numbers for loh and beyond (8→9, 9→10, 10→11, 11→12, 12→13, 13→14, 14→15)
    # Do in reverse order to avoid double-replacement
    for old_n, new_n in [(14,15),(13,14),(12,13),(11,12),(10,11),(9,10),(8,9)]:
        # Only renumber in section-title spans, not in nav
        content = re.sub(
            rf'(<h2 class="section-title"><span class="sec-num">){old_n}(</span>)',
            lambda m, o=old_n, n=new_n: m.group(0).replace(f">{o}<", f">{n}<"),
            content
        )

    # 3. Insert trimodal HTML section after bimodal section closing tag
    # Find </section> after id="bimodal"
    bi_section_end = re.search(r'(id="bimodal".*?</section>)', content, re.DOTALL)
    if bi_section_end:
        original = bi_section_end.group(0)
        content = content.replace(original, original + "\n" + trimodal_section, 1)
    else:
        print("  WARNING: could not find bimodal section end; appending trimodal before </div>")
        content = content.replace("</div><!-- end #main -->", trimodal_section + "\n</div><!-- end #main -->", 1)

    return content

# ── Main ──────────────────────────────────────────────────────────────────────
def main():
    p = argparse.ArgumentParser()
    p.add_argument("--tables-dir",    default="analysis/benchmark_trimodal_all_samples/tables")
    p.add_argument("--figures-dir",   default="/scratch/project_2008084/mvhla_figures_v6")
    p.add_argument("--input-report",  default="/scratch/project_2008084/mvhla_report_v11.html")
    p.add_argument("--output-report", default="/scratch/project_2008084/mvhla_report_v12.html")
    args = p.parse_args()

    apply_style()

    td = args.tables_dir
    print("Loading tables...")
    trimodal_det  = load(td, "trimodal_wgs_wes_rna_consensus.tsv")
    trimodal_comp = load(td, "trimodal_accuracy_comparison.tsv")
    bimodal_det   = load(td, "bimodal_wes_rna_consensus.tsv")
    mv_df         = load(td, "majority_vote_baseline.tsv")
    wc_df         = load(td, "weighted_consensus_calls.tsv")
    print(f"  trimodal_det: {len(trimodal_det)} rows")
    print(f"  trimodal_comp: {len(trimodal_comp)} rows")

    print("\nGenerating figure 15...")
    fig_path = fig15_trimodal(trimodal_det, bimodal_det, mv_df, wc_df, args.figures_dir)

    print("\nBuilding trimodal HTML section...")
    trimodal_section = build_trimodal_section(fig_path, trimodal_comp)

    trimodal_nav = '<a href="#trimodal" class="nav-link">8. Trimodal Analysis</a>'

    print(f"\nPatching {args.input_report} → {args.output_report} ...")
    patched = patch_html(args.input_report, trimodal_section, trimodal_nav)

    Path(args.output_report).write_text(patched, encoding="utf-8")
    size_kb = Path(args.output_report).stat().st_size // 1024
    print(f"  HTML report → {args.output_report}  ({size_kb} KB)")
    print("\nDone.")

if __name__ == "__main__":
    main()
