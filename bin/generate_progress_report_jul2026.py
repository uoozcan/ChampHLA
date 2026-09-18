#!/usr/bin/env python3
"""Build the champHLA July-2026 progress report as a single self-contained HTML file.

Reuses the June-2026 report's CSS (analysis/champhla_progress_report_jun2026.html), embeds every
figure as a base64 data: URI, and renders result tables from on-disk TSVs at build time, so the
output is fully offline (no external CSS/JS/font/image). Run after bin/build_ciwd_stratified.py and
figure_hub/scripts/fig_ciwd_stratified.py have (re)generated the CIWD assets.

Output: analysis/champhla_progress_report_jul2026.html
"""
import base64
import html
import re
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
JUNE = REPO / "analysis/champhla_progress_report_jun2026.html"
FIG = REPO / "analysis/figures_final"
CIWD_FIG = REPO / "figure_hub/outputs/fig_ciwd_stratified.png"
CIWD_TABLES = REPO / "analysis/ciwd_stratified/tables"
OUT = REPO / "analysis/champhla_progress_report_jul2026.html"

MOD_LABEL = {"rnaseq": "RNA-seq", "rna": "RNA-seq", "wes": "WES", "wgs": "WGS"}
MOD_BADGE = {"rnaseq": "rna", "rna": "rna", "wes": "wes", "wgs": "wgs"}
MOD_ORDER = {"rnaseq": 0, "rna": 0, "wes": 1, "wgs": 2}


# ---------------------------------------------------------------- helpers
def esc(x):
    return html.escape(str(x))


def data_uri(path):
    path = Path(path)
    b64 = base64.b64encode(path.read_bytes()).decode("ascii")
    return f"data:image/png;base64,{b64}"


def read_tsv(path):
    rows = []
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            rows.append(dict(zip(header, line.rstrip("\n").split("\t"))))
    return rows


def pct(x):
    try:
        return f"{float(x) * 100:.1f}%"
    except (TypeError, ValueError):
        return esc(x)


def extract_css():
    m = re.search(r"<style>.*?</style>", JUNE.read_text(), re.DOTALL)
    return m.group(0) if m else "<style></style>"


def figure_card(stem, title, caption, tag="new", tag_label="New"):
    src = data_uri(FIG / f"{stem}.png") if (FIG / f"{stem}.png").exists() else data_uri(stem)
    return f"""
      <div class="panel figure-card">
        <div class="figure-header">
          <h3>{esc(title)}</h3>
          <span class="fig-tag {tag}">{esc(tag_label)}</span>
        </div>
        <div class="figure-box"><img src="{src}" alt="{esc(title)}"></div>
        <div class="figure-copy"><p>{caption}</p></div>
      </div>"""


def table(headers, rows, caption=""):
    thead = "".join(f"<th>{esc(h)}</th>" for h in headers)
    body = ""
    for r in rows:
        body += "<tr>" + "".join(f"<td>{c}</td>" for c in r) + "</tr>"
    cap = f"<caption>{caption}</caption>" if caption else ""
    return f'<div class="table-wrap"><table>{cap}<thead><tr>{thead}</tr></thead><tbody>{body}</tbody></table></div>'


def badge(modality):
    return f'<span class="modality-badge {MOD_BADGE.get(modality, "")}">{esc(MOD_LABEL.get(modality, modality))}</span>'


# ---------------------------------------------------------------- data
# External validation (from analysis/EXTERNAL_VALIDATION_SUMMARY.md — curated, frozen policy).
EXT_VAL = [
    ("NCI-60", "Cancer lines", "rna", "42 / 81", "0.793", "0.817", "Adams 2005 SBT"),
    ("NCI-60", "Cancer lines", "wes", "42 / 82", "0.842", "0.842", "Adams 2005 SBT"),
    ("GIAB Ashkenazi trio", "Germline", "rna", "3 / 9", "0.889", "0.889", "Stanford SBT (gold)"),
    ("GIAB Ashkenazi trio", "Germline", "wes", "3 / 9", "1.000", "1.000", "Stanford SBT (gold)"),
    ("GIAB HG002", "Germline", "wgs", "1 / 3", "1.000", "1.000", "Chin 2020 clinical SBT"),
    ("IHWG MHC-ref cells", "Germline LCL", "wgs", "4 / 12", "0.833", "0.917", "IPD-IMGT IHIW multi-lab"),
]


# Coherent full-cohort sources (same samples across CC / Weighted / MV / per-tool), matching
# bin/build_ciwd_stratified.py. WGS uses the complete wave-2 truth-backed run.
BENCH_DIRS = {
    "rna": REPO / "analysis/1000g_realdata/benchmark_rna_truthbacked/run/tables",
    "wes": REPO / "analysis/1000g_realdata/benchmark_wes_truthbacked/run/tables",
    "wgs": REPO / "analysis/1000g_realdata/benchmark_wgs_wave2/tables",
}
METHOD_FILE = {
    "ChampionChallenger": "champion_challenger_calls.tsv",
    "WeightedConsensus": "weighted_consensus_calls.tsv",
    "MajorityVote": "majority_vote_baseline.tsv",
}


def _rate(rows):
    n = len(rows)
    if not n:
        return None, 0
    c = sum(1 for r in rows if r.get("is_correct", "0") == "1")
    return c / n, n


def core_benchmark_rows():
    """Per-modality, all columns from the SAME samples: best single tool, MV, Weighted, ChampHLA (CC)."""
    out = []
    for mod in ["rna", "wes", "wgs"]:
        tables = BENCH_DIRS[mod]
        method_rate = {}
        modality = mod
        for method, fname in METHOD_FILE.items():
            rows = read_tsv(tables / fname)
            if rows:
                modality = rows[0].get("modality", mod)
            method_rate[method], n = _rate(rows)
        # best single tool from the harmonized per-tool rows of the same cohort
        harm = read_tsv(tables / "harmonized_benchmark_rows.tsv")
        by_tool = {}
        for r in harm:
            by_tool.setdefault(r.get("tool", ""), []).append(r)
        best_tool, best_rate = max(
            ((t, _rate(rs)[0]) for t, rs in by_tool.items()),
            key=lambda kv: (kv[1] is not None, kv[1] or 0))
        out.append({
            "modality": modality,
            "n": n,
            "best_tool": best_tool,
            "best": f"{best_rate:.4f}" if best_rate is not None else "",
            "mv": f"{method_rate['MajorityVote']:.4f}" if method_rate['MajorityVote'] is not None else "",
            "wc": f"{method_rate['WeightedConsensus']:.4f}" if method_rate['WeightedConsensus'] is not None else "",
            "cc": f"{method_rate['ChampionChallenger']:.4f}" if method_rate['ChampionChallenger'] is not None else "",
        })
    return out


# ---------------------------------------------------------------- sections
def build():
    css = extract_css()

    # --- external validation table
    ext_rows = []
    for cohort, ctype, mod, n, cc, mv, truth in EXT_VAL:
        star = " ★" if cc == mv and float(cc) >= 0.99 else ""
        ext_rows.append([esc(cohort), esc(ctype), badge(mod), esc(n),
                         f"<strong>{pct(cc)}</strong>{star}", pct(mv), esc(truth)])
    ext_table = table(
        ["Cohort", "Type", "Modality", "n (samples / gene-rows)", "ChampHLA (CC)", "MajorityVote", "Truth"],
        ext_rows,
        "Frozen 1000G-derived Champion–Challenger policy applied unchanged to six independent, non-1000G cohorts.")

    # --- core benchmark table
    cb = core_benchmark_rows()
    cb_rows = []
    for r in cb:
        cb_rows.append([badge(r["modality"]), esc(r["n"]),
                        f'{esc(r["best_tool"])} ({pct(r["best"])})',
                        pct(r["mv"]), pct(r["wc"]), f'<strong>{pct(r["cc"])}</strong>'])
    core_table = table(
        ["Modality", "Gene-rows", "Best single tool", "MajorityVote", "WeightedConsensus", "ChampHLA (CC)"],
        cb_rows,
        "1000G truth-backed cohort, all methods on the same samples (two-field overall correct-call rate). "
        "CC = Champion–Challenger consensus. WGS = complete wave-2 truth-backed subset (smaller n).")

    # --- CIWD stratified table
    strat = read_tsv(CIWD_TABLES / "summary_ciwd_stratified.tsv")
    strat.sort(key=lambda r: (MOD_ORDER.get(r["modality"], 9), r["method"]))
    STRAT_LBL = {"common": "Common", "intermediate": "Intermediate", "well_documented": "Well-doc.",
                 "not_ciwd": "Not-CIWD", "unknown": "Novel/unseen"}
    ciwd_rows = []
    for r in strat:
        ciwd_rows.append([badge(r["modality"]), esc(r["method"]),
                          esc(STRAT_LBL.get(r["ciwd_stratum"], r["ciwd_stratum"])),
                          esc(r["gene_rows"]), pct(r["overall_correct_call_rate"])])
    ciwd_table = table(
        ["Modality", "Method", "CIWD stratum", "Gene-rows", "Concordance"],
        ciwd_rows,
        "Two-field concordance by CIWD 3.0.0 allele-commonness stratum (rarer truth allele).")

    # --- CIWD plausibility table (top implausible-rate tools)
    plaus = read_tsv(CIWD_TABLES / "summary_ciwd_plausibility.tsv")
    plaus.sort(key=lambda r: float(r["implausible_rate"] or 0), reverse=True)
    pl_rows = []
    for r in plaus[:10]:
        pl_rows.append([badge(r["modality"]), esc(r["tool"]), esc(r["callable_calls"]),
                        esc(r["implausible_calls"]), pct(r["implausible_rate"]),
                        pct(r["implausible_error_rate"]) if r["implausible_error_rate"] else "–"])
    plaus_table = table(
        ["Modality", "Tool", "Callable calls", "Implausible", "Implausible rate", "Error rate (of implausible)"],
        pl_rows,
        "QC view: calls that are not-CIWD / absent from the catalogue, and how error-enriched they are. Top 10 by implausible rate.")

    # --- figures
    fig_silver = figure_card("figure_12_silver_truth_hprc", "Silver-standard truth validation (HPRC)",
                             "Concordance against HPRC assembly-derived silver truth, supporting the external-validation generalisation claim.",
                             "updated", "Updated")
    fig_acc = figure_card("figure_2_accuracy_comparison", "Accuracy comparison across tools & consensus",
                          "Per-tool vs consensus two-field accuracy across modalities.", "updated", "Updated")
    fig_ciwd = figure_card(str(CIWD_FIG), "CIWD-stratified concordance (NEW)",
                           "ChampHLA (CC) vs Weighted vs Majority two-field concordance by allele-commonness stratum, one panel per modality. Consensus accuracy holds beyond common alleles; on RNA novel/unseen alleles CC reaches 100% where the weighted consensus drops to 50%.",
                           "new", "New")
    fig_workflow = figure_card("figure_1_workflow_architecture", "Pipeline workflow & architecture",
                               "End-to-end champHLA architecture: multi-tool typing → harmonisation → Champion–Challenger consensus.", "updated", "Updated")
    fig_ccmech = figure_card("figure_11a_cc_mechanism", "Champion–Challenger mechanism",
                             "How the modality-specific champion is overridden by challengers under the reliability-weighted gate.", "updated", "Updated")
    fig_mvcc = figure_card("figure_11b_mv_vs_cc", "Majority Voting vs Champion–Challenger",
                           "Head-to-head of the two consensus strategies.", "updated", "Updated")
    fig_pergene = figure_card("figure_3_per_gene_gains", "Per-gene consensus gains",
                              "Consensus accuracy gain over the best single tool, by HLA locus.", "updated", "Updated")
    fig_calib = figure_card("figure_4_confidence_calibration", "Confidence calibration",
                            "Calibration of the read_confidence_v2 scores used in weighting.", "updated", "Updated")
    fig_abst = figure_card("figure_5_abstention_tradeoff", "Abstention trade-off",
                           "Callable rate vs accuracy as the abstention threshold varies.", "updated", "Updated")
    fig_weights = figure_card("figure_7_confidence_weights", "Confidence weights (0.7/0.3)",
                              "The locked reliability/confidence weighting: 0.7×base_reliability + 0.3×confidence.", "updated", "Updated")
    fig_perf = figure_card("figure_8_computational_performance", "Computational performance",
                           "Runtime and memory per tool — the resource envelope carried into the Roihu migration.", "updated", "Updated")

    # --- navigation entries
    nav = [
        ("summary", "Executive summary"),
        ("external", "External validation"),
        ("ciwd", "CWD/CIWD stratification"),
        ("tools", "New tools"),
        ("benchmark", "Core benchmark results"),
        ("infra", "Puhti → Roihu migration"),
        ("manuscript", "Manuscript status"),
        ("next", "Next steps"),
    ]
    nav_html = "".join(f'<a href="#{sid}">{esc(lbl)}</a>' for sid, lbl in nav)

    doc = f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>champHLA — July 2026 Progress Report</title>
{css}
</head>
<body>
<div class="layout">
  <aside>
    <div class="brand">
      <p class="kicker">Progress Report</p>
      <h1>champHLA</h1>
      <p>Multi-tool HLA typing &amp; Champion–Challenger consensus</p>
    </div>
    <div class="sidebar-meta">
      <div class="sidebar-chip"><strong>Date</strong>13 July 2026</div>
      <div class="sidebar-chip"><strong>Since</strong>June 2026 report</div>
      <div class="sidebar-chip"><strong>Status</strong>Pre-submission (Bioinformatics)</div>
    </div>
    <nav class="nav">{nav_html}</nav>
  </aside>
  <main>
    <div class="hero">
      <p class="eyebrow">Update · July 2026</p>
      <h2>External validation, CWD/CIWD stratification &amp; the Roihu migration</h2>
      <p>Since June, champHLA completed a six-cohort external validation across RNA/WES/WGS on independent
      gold and consensus truth, added a CWD/CIWD allele-commonness stratification and QC layer, integrated
      two further tools (Locityper, Immuannot), and migrated the whole pipeline from CSC Puhti to CSC Roihu
      with a verified end-to-end run.</p>
    </div>

    <section id="summary">
      <div class="section-header"><h2>Executive summary</h2>
      <p>What changed since the June report, at a glance.</p></div>
      <div class="grid metrics-grid">
        <div class="metric"><div class="label">External cohorts</div><div class="value">6</div><div class="note">RNA + WES + WGS, all non-1000G</div></div>
        <div class="metric"><div class="label">Germline WES</div><div class="value">9/9</div><div class="note">ChampHLA perfect (GIAB gold truth)</div></div>
        <div class="metric"><div class="label">New tools</div><div class="value">+2</div><div class="note">Locityper, Immuannot</div></div>
        <div class="metric"><div class="label">CIWD catalogue</div><div class="value">3,249</div><div class="note">two-field alleles, annotation layer</div></div>
        <div class="metric"><div class="label">Roihu run</div><div class="value">~5 min</div><div class="note">verified end-to-end, 4 tools</div></div>
      </div>
      <div class="grid two-col" style="margin-top:18px;">
        <div class="callout success"><strong>Consistent finding.</strong> Across all six external arms and all three
        modalities, ChampHLA (Champion–Challenger) matches majority voting within overlapping CIs, and is perfect
        on germline WES (9/9) and near-perfect on germline RNA (8/9). Where it trails, the cause is the frozen
        modality-specific <em>champion</em> choice — recoverable — not the consensus logic.</div>
        <div class="callout warn"><strong>Infrastructure.</strong> Pipeline migrated Puhti → Roihu; OptiType,
        arcasHLA, T1K and native SpecHLA validated end-to-end. Remaining to enable all tools: transfer the HLA-HD
        database, fix the Kourami samtools wrapper, and stage input/reference data.</div>
      </div>
    </section>

    <section id="external">
      <div class="section-header"><h2>External validation</h2>
      <p>The frozen 1000G-derived Champion–Challenger policy + weights (no re-learning) applied to independent
      cohorts with experimental ground truth. ★ marks perfect germline arms.</p></div>
      {ext_table}
      <div class="grid two-col" style="margin-top:18px;">
        {fig_silver}
        {fig_acc}
      </div>
    </section>

    <section id="ciwd">
      <div class="section-header"><h2>CWD/CIWD allele-commonness stratification <span class="fig-tag new">New</span></h2>
      <p>We vendored the <strong>Common, Intermediate &amp; Well-Documented (CIWD) v3.0.0</strong> catalogue
      (Hurley et al. 2020, <em>HLA</em> 95:516–531; ~8M donor-registry typings across seven population groups)
      as an <em>annotation layer</em> — not a truth cohort (it carries no per-sample calls or reads). The older
      ASHI CWD 2.0.0 designation is retained in a column for lab/tool comparability. It answers the reviewer
      question: <em>does the consensus win only on easy common alleles?</em></p></div>
      {fig_ciwd}
      <div style="margin-top:18px;">{ciwd_table}</div>
      <div class="section-header" style="margin-top:26px;"><h2 style="font-size:20px;">Plausibility QC</h2>
      <p>Using CIWD as a biological-plausibility flag: a callable call that is not-CIWD or absent from the
      catalogue is a candidate typing error or novel allele. Such calls are strongly error-enriched, supporting
      the flag's use as a QC signal.</p></div>
      {plaus_table}
    </section>

    <section id="tools">
      <div class="section-header"><h2>New tools</h2>
      <p>Two further typing tools integrated end-to-end into the ensemble.</p></div>
      <div class="grid two-col">
        <div class="panel"><h3>Locityper</h3><p>Targeted genotyper wired into the pipeline
        (<code>--tools locityper</code>); requires a Locityper DB (<code>bin/build_locityper_db.sh</code>).
        Silver-truth generation supported via <code>--generate_truth</code>.</p></div>
        <div class="panel"><h3>Immuannot</h3><p>Contig-based HLA annotation integrated
        (<code>--tools immuannot</code>); needs <code>--immuannot_refdata</code> / <code>--contigs_dir</code>.
        Both are documented in the README overhaul.</p></div>
      </div>
    </section>

    <section id="benchmark">
      <div class="section-header"><h2>Core benchmark results</h2>
      <p>Full-cohort 1000G truth-backed two-field correct-call rates — the internal reference the external
      validation is measured against.</p></div>
      {core_table}
      <div class="grid two-col" style="margin-top:18px;">
        {fig_workflow}
        {fig_ccmech}
      </div>
      <div class="grid two-col" style="margin-top:18px;">
        {fig_mvcc}
        {fig_pergene}
      </div>
      <div class="grid two-col" style="margin-top:18px;">
        {fig_calib}
        {fig_abst}
      </div>
      <div class="grid two-col" style="margin-top:18px;">
        {fig_weights}
        {fig_perf}
      </div>
    </section>

    <section id="infra">
      <div class="section-header"><h2>Pipeline &amp; infrastructure — Puhti → Roihu migration</h2>
      <p>The whole pipeline was migrated from CSC Puhti to CSC Roihu with a new <code>roihu</code> Nextflow
      profile (<code>nextflow.config</code> + <code>conf/roihu_params.yaml</code>).</p></div>
      <div class="grid two-col">
        <div class="callout success"><strong>Verified end-to-end on Roihu.</strong> A real SLURM + Apptainer run
        on control sample HG00733 completed (Status SUCCESS, ~5 min). OptiType, arcasHLA, T1K and native SpecHLA
        all validated and cross-agree on HLA-A/B/C; SpecHLA produced full 4-field class I+II calls.</div>
        <div class="callout warn"><strong>Remaining to enable all tools.</strong> Transfer the HLA-HD database;
        fix the Kourami <code>samtools_wrapper.sh</code> (hardcoded Puhti path); stage input FASTQ/BAM and the
        reference/container assets not yet migrated.</div>
      </div>
      <div style="margin-top:18px;">
        {table(["Component", "Status on Roihu"], [
            ["<code>roihu</code> Nextflow profile", "✅ added &amp; verified"],
            ["OptiType / arcasHLA / T1K / native SpecHLA", "✅ validated end-to-end (~5 min)"],
            ["HLA-HD", "⚠️ DB not yet transferred"],
            ["Kourami", "⚠️ samtools wrapper path fix needed"],
            ["Locityper / Immuannot", "⚠️ wired; data assets pending"],
            ["Input / reference / container assets", "⚠️ staging in progress"],
        ], "Migration status snapshot (2026-07).")}
      </div>
    </section>

    <section id="manuscript">
      <div class="section-header"><h2>Manuscript status</h2>
      <p>Pre-submission.</p></div>
      <div class="panel"><ul>
        <li>Abstract, Methods prose, and figures integrated; target journal <strong>Bioinformatics</strong>.</li>
        <li>Weight formula locked at <strong>0.7 / 0.3</strong> (reliability / confidence); sensitivity analysis complete.</li>
        <li>CWD/CIWD stratification + plausibility QC wired into <code>docs/CHAMPHLA_MANUSCRIPT_V3.md</code>.</li>
        <li>External-validation section extends the NCI-60 results into a broader RNA+WES+WGS generalisation claim.</li>
      </ul></div>
    </section>

    <section id="next">
      <div class="section-header"><h2>Next steps</h2>
      <p>Open items to close the loop.</p></div>
      <div class="panel"><ul>
        <li>Close the Roihu gaps (HLA-HD DB, Kourami wrapper, data staging) to enable an all-tools run.</li>
        <li>Run the full all-tools benchmark on Roihu and confirm parity with Puhti.</li>
        <li>Finalise the external-validation manuscript section and supplementary CIWD tables.</li>
        <li>Regenerate the CIWD-stratified figure at full-cohort WGS scale once the WGS truth-backed run is staged.</li>
      </ul></div>
    </section>

    <p style="margin-top:34px; color:var(--muted); font-size:12px;">
      champHLA (internal: pihla) · generated 2026-07-13 · self-contained · figures embedded.
    </p>
  </main>
</div>
</body>
</html>"""

    OUT.write_text(doc, encoding="utf-8")
    kb = OUT.stat().st_size / 1024
    print(f"Wrote {OUT} ({kb:.0f} KB)")


if __name__ == "__main__":
    build()
