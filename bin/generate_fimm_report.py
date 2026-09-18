#!/usr/bin/env python3
"""generate_fimm_report.py — Detailed HTML QC & results report for FIMM HLA typing.

Usage:
    module load python-data
    python3 bin/generate_fimm_report.py \
        --calls  analysis/fimm_analysis/fimm_hla_calls.tsv \
        --scrna  analysis/fimm_analysis/fimm_scrna_matched.tsv \
        --stats  analysis/fimm_analysis/fimm_concordance_stats.tsv \
        --figdir analysis/fimm_analysis \
        --out    analysis/fimm_analysis/fimm_report.html
"""

import argparse
import base64
import os
import datetime
import math
from collections import defaultdict

import pandas as pd

GENES = ["A", "B", "C"]

# ── Colour scheme ──────────────────────────────────────────────────────────────
C_AGREE   = "#d4edda"   # light green
C_DISAGREE= "#f8d7da"   # light red
C_MISSING = "#f0f0f0"   # light grey
C_REF     = "#cfe2ff"   # light blue (scRNA arcasHLA reference)
C_WARN    = "#fff3cd"   # amber
C_INFO    = "#e8f4f8"   # pale blue

MODALITY_ORDER = [
    ("scrna",   "arcashla",  "scRNA arcasHLA"),
    ("scrna",   "optitype",  "scRNA OptiType"),
    ("bulkrna", "arcasHLA",  "BulkRNA arcasHLA"),
    ("bulkrna", "OptiType",  "BulkRNA OptiType"),
    ("bulkrna", "SpecHLA",   "BulkRNA SpecHLA"),
    ("wes",     "arcasHLA",  "WES arcasHLA"),
    ("wes",     "OptiType",  "WES OptiType"),
    ("wes",     "SpecHLA",   "WES SpecHLA"),
]


# ── Helpers ────────────────────────────────────────────────────────────────────

def patient_id(sample_id):
    """'FH_6088_2' → 'FH_6088'"""
    parts = sample_id.rsplit("_", 1)
    return parts[0] if len(parts) == 2 else sample_id


def cohort(sample_id):
    if sample_id.startswith("FH_") or sample_id.startswith("FHRB_"):
        return "FIMM"
    if sample_id.startswith("VX_"):
        return "VX"
    return "Other"


def allele_display(pair_str):
    """'A*02:01|A*03:01' → 'A*02:01 / A*03:01'; missing → '—'"""
    if not isinstance(pair_str, str) or not pair_str.strip():
        return "—"
    parts = [p.strip() for p in pair_str.split("|") if p.strip()]
    return " / ".join(parts) if parts else "—"


def alleles_set(pair_str):
    """Return frozenset of alleles from pipe-delimited string."""
    if not isinstance(pair_str, str) or not pair_str.strip():
        return frozenset()
    return frozenset(p.strip() for p in pair_str.split("|") if p.strip())


def col_key(mod, tool, gene):
    """Build column name matching fimm_hla_calls.tsv convention."""
    return f"{mod}_{tool}_{gene}"


def ref_val_for(row, gene):
    """scRNA arcasHLA value for a gene (the reference)."""
    return row.get(f"scrna_arcashla_{gene}", "")


def cell_bg(val_str, ref_str):
    """Background colour for a data cell compared to reference."""
    if not isinstance(val_str, str) or not val_str.strip():
        return C_MISSING
    val_alleles = alleles_set(val_str)
    ref_alleles = alleles_set(ref_str)
    if not val_alleles:
        return C_MISSING
    if not ref_alleles:
        return C_MISSING  # no reference → neutral grey
    return C_AGREE if val_alleles == ref_alleles else C_DISAGREE


def embed_png(path):
    """Base64-encode a PNG file → data URI string."""
    if not os.path.isfile(path):
        return None
    with open(path, "rb") as f:
        data = base64.b64encode(f.read()).decode("ascii")
    return f"data:image/png;base64,{data}"


def pct_bar(pct, width=80):
    """Simple inline SVG bar for concordance percentage."""
    if math.isnan(pct):
        return "N/A"
    filled = max(0, min(width, int(pct / 100 * width)))
    color = "#198754" if pct >= 90 else "#fd7e14" if pct >= 70 else "#dc3545"
    return (
        f'<span style="display:inline-block;vertical-align:middle;margin-right:4px">'
        f'<svg width="{width}" height="12">'
        f'<rect width="{width}" height="12" fill="#e9ecef" rx="3"/>'
        f'<rect width="{filled}" height="12" fill="{color}" rx="3"/>'
        f'</svg></span>'
        f'<b>{pct:.1f}%</b>'
    )


def bool_to_badge(val):
    if val is True or val == "True":
        return '<span style="color:#198754;font-weight:bold">✓ agree</span>'
    if val is False or val == "False":
        return '<span style="color:#dc3545;font-weight:bold">✗ discordant</span>'
    return '<span style="color:#6c757d">– missing</span>'


def h(text):
    """HTML-escape."""
    return str(text).replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")


# ── Data loading ───────────────────────────────────────────────────────────────

def load_data(calls_path, scrna_path, stats_path):
    calls_df = pd.read_csv(calls_path, sep="\t", dtype=str).fillna("")
    scrna_df = pd.read_csv(scrna_path, sep="\t", dtype=str).fillna("")
    stats_df = pd.read_csv(stats_path, sep="\t")
    scrna_df = scrna_df.set_index("sample_id", drop=False)
    calls_df = calls_df.set_index("sample_id", drop=False)
    return calls_df, scrna_df, stats_df


# ── QC issue detection ─────────────────────────────────────────────────────────

def detect_issues(calls_df, scrna_df):
    issues = []

    for sid, row in scrna_df.iterrows():
        opt_idx = row.get("opt_idx", "")
        score   = row.get("match_score", "")

        if not opt_idx or opt_idx == "nan":
            issues.append({"sample": sid, "severity": "WARNING",
                           "type": "No scRNA OptiType match",
                           "detail": "OptiType FASTQ run may have failed"})
        elif score and score != "nan":
            try:
                sc = int(float(score))
                if sc < 3:
                    issues.append({"sample": sid, "severity": "WARNING",
                                   "type": "Partial scRNA match",
                                   "detail": f"Match score {sc}/3 — only {sc} gene(s) matched"})
            except ValueError:
                pass

        for g in GENES:
            agr = row.get(f"scrna_agree_{g}", "")
            if agr in ("False", "false"):
                arc = allele_display(row.get(f"scrna_arcashla_{g}", ""))
                opt = allele_display(row.get(f"scrna_optitype_{g}", ""))
                issues.append({"sample": sid, "severity": "INFO",
                               "type": f"Discordant HLA-{g} in scRNA",
                               "detail": f"arcasHLA: {arc} | OptiType: {opt}"})

    for sid, row in calls_df.iterrows():
        has_bulk = any(row.get(f"bulkrna_{tool}_{g}", "").strip()
                       for mod, tool, _ in MODALITY_ORDER if mod == "bulkrna"
                       for g in GENES)
        has_wes  = any(row.get(f"wes_{tool}_{g}", "").strip()
                       for mod, tool, _ in MODALITY_ORDER if mod == "wes"
                       for g in GENES)
        if not has_bulk:
            issues.append({"sample": sid, "severity": "INFO",
                           "type": "No BulkRNA data", "detail": "No BulkRNA typing results found"})
        if not has_wes:
            issues.append({"sample": sid, "severity": "INFO",
                           "type": "No WES data", "detail": "No WES typing results found"})

    return issues


# ── HTML components ────────────────────────────────────────────────────────────

CSS = """
* { box-sizing: border-box; }
body { font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', sans-serif;
       margin: 0; padding: 0; background: #f8f9fa; color: #212529; font-size: 14px; }
.container { max-width: 1400px; margin: 0 auto; padding: 20px; }
h1 { font-size: 1.6rem; margin-bottom: 4px; }
h2 { font-size: 1.2rem; border-bottom: 2px solid #dee2e6; padding-bottom: 6px;
     margin-top: 32px; margin-bottom: 14px; }
h3 { font-size: 1.0rem; margin: 10px 0 6px 0; }
.meta { color: #6c757d; font-size: 12px; margin-bottom: 20px; }
.kpi-row { display: flex; gap: 16px; flex-wrap: wrap; margin-bottom: 20px; }
.kpi { background: #fff; border: 1px solid #dee2e6; border-radius: 8px;
       padding: 14px 20px; min-width: 160px; }
.kpi-val { font-size: 2rem; font-weight: bold; color: #0d6efd; }
.kpi-label { font-size: 12px; color: #6c757d; }
table { border-collapse: collapse; width: 100%; font-size: 13px; }
th { background: #343a40; color: #fff; padding: 6px 8px; text-align: left;
     white-space: nowrap; cursor: pointer; user-select: none; }
th:hover { background: #495057; }
td { padding: 5px 8px; border-bottom: 1px solid #dee2e6; vertical-align: middle; }
tr:hover td { filter: brightness(0.96); }
.table-wrap { overflow-x: auto; background: #fff; border: 1px solid #dee2e6;
              border-radius: 6px; margin-bottom: 20px; }
.badge { display: inline-block; padding: 2px 8px; border-radius: 12px;
         font-size: 11px; font-weight: bold; }
.badge-warn { background: #dc3545; color: #fff; }
.badge-info { background: #fd7e14; color: #fff; }
.badge-ok   { background: #198754; color: #fff; }
details { background: #fff; border: 1px solid #dee2e6; border-radius: 8px;
          margin-bottom: 10px; }
summary { padding: 12px 16px; cursor: pointer; font-weight: 600;
          list-style: none; display: flex; align-items: center; gap: 10px; }
summary::-webkit-details-marker { display: none; }
summary::before { content: '▶'; font-size: 10px; color: #6c757d;
                  transition: transform 0.2s; }
details[open] summary::before { transform: rotate(90deg); }
.patient-body { padding: 12px 16px 16px; border-top: 1px solid #dee2e6; }
.sample-section { margin-bottom: 18px; }
.sample-title { font-weight: 600; font-size: 13px; color: #0d6efd;
                margin-bottom: 6px; }
.hla-table { width: auto; min-width: 700px; }
.hla-table th { font-size: 12px; padding: 5px 7px; }
.hla-table td { font-size: 12px; padding: 4px 7px; white-space: nowrap; }
.hla-table td.ref { background: """ + C_REF + """ !important; font-weight: 600; }
.alert-box { padding: 6px 10px; border-radius: 4px; margin-bottom: 6px;
             font-size: 12px; }
.alert-warn { background: #fff3cd; border-left: 4px solid #ffc107; }
.alert-info { background: #e8f4f8; border-left: 4px solid #0dcaf0; }
.multi-compare { margin-top: 10px; }
.multi-compare table { width: auto; }
.toc { background: #fff; border: 1px solid #dee2e6; border-radius: 6px;
       padding: 12px 16px; display: inline-block; margin-bottom: 20px;
       min-width: 260px; }
.toc a { text-decoration: none; color: #0d6efd; font-size: 13px;
         display: block; padding: 2px 0; }
.toc a:hover { text-decoration: underline; }
.section-anchor { scroll-margin-top: 20px; }
.cohort-header { font-size: 11px; color: #6c757d; font-weight: normal; }
figure { margin: 0 0 20px; }
figure img { max-width: 100%; border: 1px solid #dee2e6; border-radius: 6px; }
figcaption { font-size: 12px; color: #6c757d; margin-top: 4px; }
.fig-row { display: flex; gap: 20px; flex-wrap: wrap; }
.fig-row figure { flex: 1 1 400px; }
.concordance-row-scrna { background: #f0fff0; }
.concordance-row-bulkrna { background: #f0f8ff; }
.concordance-row-wes { background: #fffaf0; }
"""

SORT_JS = """
function sortTable(tableId, col) {
    var tbl = document.getElementById(tableId);
    var rows = Array.from(tbl.querySelectorAll('tbody tr'));
    var asc = tbl.dataset.sortCol == col && tbl.dataset.sortDir == 'asc';
    rows.sort(function(a, b) {
        var va = a.cells[col] ? a.cells[col].innerText.trim() : '';
        var vb = b.cells[col] ? b.cells[col].innerText.trim() : '';
        var na = parseFloat(va), nb = parseFloat(vb);
        if (!isNaN(na) && !isNaN(nb)) return asc ? nb - na : na - nb;
        return asc ? vb.localeCompare(va) : va.localeCompare(vb);
    });
    var tbody = tbl.querySelector('tbody');
    rows.forEach(function(r) { tbody.appendChild(r); });
    tbl.dataset.sortCol = col;
    tbl.dataset.sortDir = asc ? 'desc' : 'asc';
}
"""


def render_hla_call_table(row, scrna_row, sample_id):
    """Render the per-gene HLA call table for one sample."""
    cols = [
        ("scrna_arcashla",  "scRNA<br>arcasHLA",  True),
        ("scrna_optitype",  "scRNA<br>OptiType",  False),
        ("bulkrna_arcasHLA","BulkRNA<br>arcasHLA",False),
        ("bulkrna_OptiType","BulkRNA<br>OptiType",False),
        ("bulkrna_SpecHLA", "BulkRNA<br>SpecHLA", False),
        ("wes_arcasHLA",    "WES<br>arcasHLA",    False),
        ("wes_OptiType",    "WES<br>OptiType",    False),
        ("wes_SpecHLA",     "WES<br>SpecHLA",     False),
    ]

    html = ['<table class="hla-table"><thead><tr>']
    html.append('<th>Gene</th>')
    for _, label, is_ref in cols:
        ref_mark = ' 🔵' if is_ref else ''
        html.append(f'<th>{label}{ref_mark}</th>')
    html.append('</tr></thead><tbody>')

    for g in GENES:
        ref_str = row.get(f"scrna_arcashla_{g}", "")
        html.append(f'<tr><td><b>HLA-{g}</b></td>')
        for col_prefix, _, is_ref in cols:
            key = f"{col_prefix}_{g}"
            val = row.get(key, "")
            display = allele_display(val)
            if is_ref:
                html.append(f'<td class="ref">{h(display)}</td>')
            else:
                bg = cell_bg(val, ref_str)
                html.append(f'<td style="background:{bg}">{h(display)}</td>')
        html.append('</tr>')

    html.append('</tbody></table>')
    return "".join(html)


def render_scrna_concordance_badges(scrna_row):
    """✓/✗/– badges for each gene."""
    parts = []
    for g in GENES:
        agr = scrna_row.get(f"scrna_agree_{g}", "") if scrna_row is not None else ""
        badge = bool_to_badge(agr if agr else None)
        parts.append(f'<span style="margin-right:12px"><b>HLA-{g}</b>: {badge}</span>')
    return "".join(parts)


def render_multi_sample_comparison(patient_samples, calls_df):
    """For patients with >1 sample, show a side-by-side allele comparison table."""
    html = ['<div class="multi-compare">']
    html.append('<h3>Cross-sample allele consistency</h3>')
    html.append('<table><thead><tr><th>Gene</th><th>Allele</th>')
    for sid in patient_samples:
        html.append(f'<th>{h(sid)}<br><span style="font-size:10px;font-weight:normal">(scRNA arcasHLA)</span></th>')
    html.append('</tr></thead><tbody>')

    for g in GENES:
        allele_sets = []
        for sid in patient_samples:
            row = calls_df.loc[sid] if sid in calls_df.index else None
            val = row.get(f"scrna_arcashla_{g}", "") if row is not None else ""
            allele_sets.append(alleles_set(val))

        for allele_label, allele_idx in [("allele 1", 0), ("allele 2", 1)]:
            html.append(f'<tr><td>{"HLA-" + g if allele_idx == 0 else ""}</td><td>{allele_label}</td>')
            all_vals = []
            for aset in allele_sets:
                sorted_a = sorted(aset)
                val = sorted_a[allele_idx] if len(sorted_a) > allele_idx else "—"
                all_vals.append(val)
            consistent = len(set(v for v in all_vals if v != "—")) <= 1
            bg = C_AGREE if consistent else C_DISAGREE
            for val in all_vals:
                html.append(f'<td style="background:{bg}">{h(val)}</td>')
            html.append('</tr>')

    html.append('</tbody></table></div>')
    return "".join(html)


def render_issues_table(issues):
    warn_count = sum(1 for i in issues if i["severity"] == "WARNING")
    info_count  = sum(1 for i in issues if i["severity"] == "INFO")

    html = [f'<p>Found <b>{warn_count}</b> warnings and <b>{info_count}</b> informational notes.</p>']
    html.append('<div class="table-wrap"><table id="tbl-issues"><thead><tr>')
    for i, hdr in enumerate(["#", "Severity", "Sample", "Issue", "Detail"]):
        html.append(f'<th onclick="sortTable(\'tbl-issues\',{i})">{hdr}</th>')
    html.append('</tr></thead><tbody>')

    severity_order = {"WARNING": 0, "INFO": 1}
    sorted_issues = sorted(issues, key=lambda x: (severity_order.get(x["severity"], 9), x["sample"]))

    for idx, iss in enumerate(sorted_issues, 1):
        sev = iss["severity"]
        badge_cls = "badge-warn" if sev == "WARNING" else "badge-info"
        html.append(
            f'<tr>'
            f'<td>{idx}</td>'
            f'<td><span class="badge {badge_cls}">{h(sev)}</span></td>'
            f'<td><code>{h(iss["sample"])}</code></td>'
            f'<td>{h(iss["type"])}</td>'
            f'<td style="color:#495057">{h(iss["detail"])}</td>'
            f'</tr>'
        )
    html.append('</tbody></table></div>')
    return "".join(html)


def render_concordance_stats(stats_df):
    html = ['<div class="table-wrap"><table id="tbl-stats"><thead><tr>']
    headers = ["Modality", "Reference", "Query", "Gene", "N callable", "N agree", "% Concordance", "95% CI"]
    for i, hdr in enumerate(headers):
        html.append(f'<th onclick="sortTable(\'tbl-stats\',{i})">{hdr}</th>')
    html.append('</tr></thead><tbody>')

    for _, row in stats_df.iterrows():
        mod = str(row["modality"]).lower()
        row_cls = f"concordance-row-{mod}"
        pct = float(row["pct_agree"]) if not math.isnan(float(row["pct_agree"])) else float("nan")
        ci_lo = float(row["ci_lo"])
        ci_hi = float(row["ci_hi"])
        ci_str = f"{ci_lo:.1f}–{ci_hi:.1f}%"
        html.append(
            f'<tr class="{row_cls}">'
            f'<td>{h(row["modality"])}</td>'
            f'<td>{h(row["reference"])}</td>'
            f'<td>{h(row["query"])}</td>'
            f'<td><b>HLA-{h(row["gene"])}</b></td>'
            f'<td>{int(row["n_callable"])}</td>'
            f'<td>{int(row["n_agree"])}</td>'
            f'<td>{pct_bar(pct)}</td>'
            f'<td style="font-size:11px;color:#6c757d">{ci_str}</td>'
            f'</tr>'
        )
    html.append('</tbody></table></div>')
    return "".join(html)


def render_full_data_table(calls_df, scrna_df):
    """Sortable full-data table with color coding."""
    all_cols = [c for c in calls_df.columns if c != "sample_id"]
    display_cols = ["sample_id"] + all_cols

    html = ['<div class="table-wrap"><table id="tbl-full"><thead><tr>']
    for i, col in enumerate(display_cols):
        label = col.replace("_", " ").replace("arcashla", "arcasHLA")
        html.append(f'<th onclick="sortTable(\'tbl-full\',{i})" style="font-size:11px">{h(label)}</th>')
    html.append('</tr></thead><tbody>')

    for sid, row in calls_df.iterrows():
        html.append('<tr>')
        for col in display_cols:
            val = row.get(col, "") if col != "sample_id" else sid
            if col == "sample_id":
                html.append(f'<td><code>{h(sid)}</code></td>')
            else:
                # Extract gene from column name
                parts = col.rsplit("_", 1)
                gene  = parts[-1] if parts[-1] in GENES else None
                if gene:
                    ref_str = row.get(f"scrna_arcashla_{gene}", "")
                    if col == f"scrna_arcashla_{gene}":
                        html.append(f'<td style="background:{C_REF};font-size:11px">{h(allele_display(val))}</td>')
                    else:
                        bg = cell_bg(val, ref_str)
                        html.append(f'<td style="background:{bg};font-size:11px">{h(allele_display(val))}</td>')
                else:
                    html.append(f'<td style="font-size:11px">{h(allele_display(val))}</td>')
        html.append('</tr>')

    html.append('</tbody></table></div>')
    return "".join(html)


def render_patient_cards(calls_df, scrna_df, issues_by_sample):
    """One <details> card per patient, sorted FIMM → VX."""
    # Group samples by patient
    patient_map = defaultdict(list)
    for sid in calls_df.index:
        patient_map[patient_id(sid)].append(sid)

    # Sort patients
    def patient_sort_key(pid):
        order = {"FIMM": 0, "VX": 1}
        c = cohort(pid + "_X")
        return (order.get(c, 9), pid)

    sorted_patients = sorted(patient_map.keys(), key=patient_sort_key)

    html = []
    for pid in sorted_patients:
        samples = sorted(patient_map[pid])
        c = cohort(samples[0])
        n_samples = len(samples)

        # Collect issues for this patient
        pat_issues = []
        for sid in samples:
            pat_issues.extend(issues_by_sample.get(sid, []))

        warn_count = sum(1 for i in pat_issues if i["severity"] == "WARNING")
        info_count  = sum(1 for i in pat_issues if i["severity"] == "INFO")

        # Build summary badges for header
        flag_html = ""
        if warn_count:
            flag_html += f' <span class="badge badge-warn">{warn_count} ⚠ WARNING</span>'
        if info_count:
            flag_html += f' <span class="badge badge-info">{info_count} ℹ INFO</span>'
        multi_badge = f' <span class="badge" style="background:#6f42c1;color:#fff">{n_samples} samples</span>' if n_samples > 1 else ""

        html.append(f'<details id="patient-{h(pid)}">')
        html.append(
            f'<summary>'
            f'<span style="font-size:14px">{h(pid)}</span>'
            f'<span class="cohort-header"> &nbsp;·&nbsp; cohort: {h(c)}</span>'
            f'{multi_badge}{flag_html}'
            f'</summary>'
        )
        html.append('<div class="patient-body">')

        # Per-sample inline alerts
        if pat_issues:
            for iss in pat_issues:
                cls = "alert-warn" if iss["severity"] == "WARNING" else "alert-info"
                html.append(
                    f'<div class="alert-box {cls}">'
                    f'<b>{h(iss["severity"])}</b> [{h(iss["sample"])}] '
                    f'{h(iss["type"])}'
                    + (f' — {h(iss["detail"])}' if iss["detail"] else "")
                    + '</div>'
                )

        # Each sample's HLA table
        for sid in samples:
            row = calls_df.loc[sid] if sid in calls_df.index else pd.Series()
            scrna_row = scrna_df.loc[sid] if sid in scrna_df.index else None

            html.append('<div class="sample-section">')
            html.append(f'<div class="sample-title">Sample: {h(sid)}</div>')

            if scrna_row is not None and isinstance(scrna_row, pd.Series):
                html.append('<div style="margin-bottom:6px">')
                html.append(render_scrna_concordance_badges(scrna_row))
                score = scrna_row.get("match_score", "")
                if score and score not in ("", "nan"):
                    html.append(f' &nbsp;<span style="font-size:11px;color:#6c757d">match score: {h(score)}/3</span>')
                html.append('</div>')

            html.append(render_hla_call_table(row, scrna_row, sid))
            html.append('</div>')

        # Multi-sample cross-comparison
        if n_samples > 1:
            html.append(render_multi_sample_comparison(samples, calls_df))

        html.append('</div>')  # patient-body
        html.append('</details>')

    return "".join(html)


def render_toc(calls_df):
    patient_map = defaultdict(list)
    for sid in calls_df.index:
        patient_map[patient_id(sid)].append(sid)

    def patient_sort_key(pid):
        order = {"FIMM": 0, "VX": 1}
        c = cohort(pid + "_X")
        return (order.get(c, 9), pid)

    sorted_patients = sorted(patient_map.keys(), key=patient_sort_key)

    html = ['<div class="toc"><b>Quick navigation</b><br>']
    html.append('<a href="#sec-summary">Cohort summary</a>')
    html.append('<a href="#sec-qc">QC issues</a>')
    html.append('<a href="#sec-concordance">Concordance statistics</a>')
    html.append('<a href="#sec-patients">Patient results</a>')
    html.append('<hr style="margin:4px 0">')
    for pid in sorted_patients:
        n = len(patient_map[pid])
        label = f"{pid}" + (f" ({n} samples)" if n > 1 else "")
        html.append(f'<a href="#patient-{pid}">{h(label)}</a>')
    html.append('<hr style="margin:4px 0">')
    html.append('<a href="#sec-fulltable">Full data table</a>')
    html.append('<a href="#sec-methods">Methods</a>')
    html.append('</div>')
    return "".join(html)


# ── Main ───────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--calls",  required=True)
    parser.add_argument("--scrna",  required=True)
    parser.add_argument("--stats",  required=True)
    parser.add_argument("--figdir", default=".")
    parser.add_argument("--out",    required=True)
    args = parser.parse_args()

    print("Loading data …")
    calls_df, scrna_df, stats_df = load_data(args.calls, args.scrna, args.stats)
    n_samples = len(calls_df)

    patient_map = defaultdict(list)
    for sid in calls_df.index:
        patient_map[patient_id(sid)].append(sid)
    n_patients = len(patient_map)

    cohort_counts = defaultdict(int)
    for pid in patient_map:
        cohort_counts[cohort(pid + "_X")] += 1

    print("Detecting QC issues …")
    issues = detect_issues(calls_df, scrna_df)
    issues_by_sample = defaultdict(list)
    for iss in issues:
        issues_by_sample[iss["sample"]].append(iss)

    print("Embedding figures …")
    fig_paths = {
        "heatmap": os.path.join(args.figdir, "fig_scrna_concordance_heatmap.png"),
        "barchart": os.path.join(args.figdir, "fig_concordance_by_gene.png"),
        "dotplot":  os.path.join(args.figdir, "fig_crossmodal_comparison.png"),
    }
    figs = {k: embed_png(v) for k, v in fig_paths.items()}

    warn_count = sum(1 for i in issues if i["severity"] == "WARNING")
    info_count  = sum(1 for i in issues if i["severity"] == "INFO")
    no_opt = sum(1 for i in issues if i["type"] == "No scRNA OptiType match")
    now = datetime.datetime.now().strftime("%Y-%m-%d %H:%M")

    print("Rendering patient cards …")
    patient_cards_html = render_patient_cards(calls_df, scrna_df, issues_by_sample)

    print("Composing HTML …")

    # KPI cards
    kpi_html = (
        f'<div class="kpi-row">'
        f'<div class="kpi"><div class="kpi-val">{n_patients}</div><div class="kpi-label">Patients</div></div>'
        f'<div class="kpi"><div class="kpi-val">{n_samples}</div><div class="kpi-label">Samples</div></div>'
        f'<div class="kpi"><div class="kpi-val">{cohort_counts.get("FIMM",0)}</div><div class="kpi-label">FIMM patients</div></div>'
        f'<div class="kpi"><div class="kpi-val">{cohort_counts.get("VX",0)}</div><div class="kpi-label">VX patients</div></div>'
        f'<div class="kpi"><div class="kpi-val" style="color:#dc3545">{warn_count}</div><div class="kpi-label">Warnings</div></div>'
        f'<div class="kpi"><div class="kpi-val" style="color:#fd7e14">{info_count}</div><div class="kpi-label">Info notes</div></div>'
        f'</div>'
        f'<p>{n_samples} samples from {n_patients} patients; '
        f'{no_opt} samples missing scRNA OptiType match; '
        f'{sum(1 for p,s in patient_map.items() if len(s)>1)} patients with multiple samples.</p>'
    )

    # Figures
    fig_html = '<div class="fig-row">'
    if figs["heatmap"]:
        fig_html += f'<figure><img src="{figs["heatmap"]}" alt="scRNA concordance heatmap"><figcaption>Fig. 1 — scRNA concordance heatmap (OptiType vs arcasHLA). Green = agree, red = disagree, grey = missing.</figcaption></figure>'
    if figs["barchart"]:
        fig_html += f'<figure><img src="{figs["barchart"]}" alt="Per-gene concordance bar chart"><figcaption>Fig. 2 — Per-gene concordance bar chart with 95% Wilson CI.</figcaption></figure>'
    fig_html += '</div>'
    if figs["dotplot"]:
        fig_html += f'<figure><img src="{figs["dotplot"]}" alt="Cross-modality allele comparison"><figcaption>Fig. 3 — Cross-modality allele calls per sample (all tools and modalities).</figcaption></figure>'

    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>FIMM HLA Typing — QC &amp; Results Report</title>
<style>
{CSS}
</style>
</head>
<body>
<div class="container">
  <h1>FIMM HLA Typing — QC &amp; Results Report</h1>
  <div class="meta">Generated: {now} &nbsp;|&nbsp; Tools: OptiType, arcasHLA, SpecHLA &nbsp;|&nbsp; Modalities: scRNA, BulkRNA, WES</div>

  {render_toc(calls_df)}

  <h2 id="sec-summary" class="section-anchor">1. Cohort Summary</h2>
  {kpi_html}

  <h2 id="sec-qc" class="section-anchor">2. QC Issue Flags</h2>
  {render_issues_table(issues)}

  <h2 id="sec-concordance" class="section-anchor">3. Concordance Statistics</h2>
  <p>Reference: scRNA arcasHLA calls. Concordance = both alleles match at 2-field resolution (e.g. A*02:01).
  Colour legend: <span style="background:{C_AGREE};padding:1px 6px;border-radius:3px">green=agree</span>
  <span style="background:{C_DISAGREE};padding:1px 6px;border-radius:3px">red=disagree</span>
  <span style="background:{C_MISSING};padding:1px 6px;border-radius:3px">grey=missing</span>
  <span style="background:{C_REF};padding:1px 6px;border-radius:3px">blue=reference</span>.
  </p>
  {render_concordance_stats(stats_df)}
  {fig_html}

  <h2 id="sec-patients" class="section-anchor">4. Patient-by-Patient Results</h2>
  <p>Each card shows all HLA typing results for a patient. Click to expand.
  For patients with multiple samples (timepoints), a cross-sample consistency table is shown at the bottom of the card.
  Cell colours compare each tool against the <b>scRNA arcasHLA</b> reference (blue column).</p>
  {patient_cards_html}

  <h2 id="sec-fulltable" class="section-anchor">5. Full Data Table</h2>
  <p>All 35 samples × all typing methods. Click column headers to sort.
  Allele pairs shown as <code>allele1 / allele2</code>. Dash (—) = not called.</p>
  {render_full_data_table(calls_df, scrna_df)}

  <h2 id="sec-methods" class="section-anchor">6. Methods</h2>
  <p>
  <b>HLA typing tools:</b> OptiType (read-mapping, 2-field resolution),
  arcasHLA (graph-based, 3-field resolution), SpecHLA (phased, 3-field resolution).<br>
  <b>Modalities:</b> scRNA-seq, bulk RNA-seq (BulkRNA), whole-exome sequencing (WES).<br>
  <b>Normalisation:</b> All alleles truncated to 2-field (e.g. <code>A*02:01:01 → A*02:01</code>).
  Allele pairs treated as unordered sets.<br>
  <b>Sample matching:</b> scRNA OptiType rows matched to arcasHLA rows by greedy 1:1 allele-fingerprint similarity
  (score = number of genes with matching allele pairs, max 3).
  Cross-modality rows matched by ≥2/3 gene-pair similarity.<br>
  <b>Concordance statistics:</b> Wilson 95% confidence intervals.
  </p>
  <hr>
  <p style="font-size:11px;color:#6c757d">Report generated by <code>generate_fimm_report.py</code> · PIHLA pipeline · project_2008084</p>
</div>
<script>
{SORT_JS}
</script>
</body>
</html>"""

    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    with open(args.out, "w", encoding="utf-8") as f:
        f.write(html)

    size_kb = os.path.getsize(args.out) // 1024
    print(f"\nReport saved: {args.out} ({size_kb} KB)")
    print(f"  {n_patients} patients, {n_samples} samples")
    print(f"  {warn_count} warnings, {info_count} info notes")


if __name__ == "__main__":
    main()
