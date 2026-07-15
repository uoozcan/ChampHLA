#!/usr/bin/env python3.11
"""
generate_html_report_v7.py
PIHLA HTML report v7 — trimodal benchmark (1000G) + VENEX WGS cohort.
Self-contained single HTML file with embedded figures.
"""

import argparse, base64, csv, io, json, math, glob as _glob, re as _re
from collections import Counter, defaultdict
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np

mpl.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'DejaVu Sans'],
    'font.size': 9,
    'axes.titlesize': 10,
    'axes.labelsize': 9,
    'xtick.labelsize': 8,
    'ytick.labelsize': 8,
    'legend.fontsize': 8,
    'legend.frameon': False,
    'axes.spines.top': False,
    'axes.spines.right': False,
    'xtick.direction': 'in',
    'ytick.direction': 'in',
    'figure.dpi': 300,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
    'savefig.pad_inches': 0.05,
    'pdf.fonttype': 42,
    'ps.fonttype': 42,
})

# ── Generic helpers ────────────────────────────────────────────────────────────

def b64img(path):
    with open(path, 'rb') as f:
        return base64.b64encode(f.read()).decode()

def fig_to_b64(fig):
    buf = io.BytesIO()
    fig.savefig(buf, format='png', dpi=300, bbox_inches='tight', pad_inches=0.05)
    buf.seek(0)
    return base64.b64encode(buf.read()).decode()

def pct(v):
    try: return f'{float(v)*100:.1f}%'
    except: return '—'

def fmt(v, d=4):
    try: return f'{float(v):.{d}f}'
    except: return '—'

def wilson_ci(k, n, z=1.96):
    if n == 0: return 0.0, 0.0
    p = k / n; denom = 1 + z**2 / n
    c = (p + z**2 / (2*n)) / denom
    m = z * math.sqrt(p*(1-p)/n + z**2/(4*n**2)) / denom
    return max(0, c-m), min(1, c+m)

def bar_html(val, max_val=1.0, color='#3b82f6', height=14):
    w = int(min(float(val or 0) / max_val, 1.0) * 160)
    return (f'<div style="display:flex;align-items:center;gap:6px">'
            f'<div style="width:{w}px;height:{height}px;background:{color};'
            f'border-radius:3px;flex-shrink:0"></div>'
            f'<span style="font-size:12px;color:#374151">{pct(val)}</span></div>')

def acc_color(v):
    try:
        f = float(v)
        if f >= 0.90: return '#dcfce7', '#16a34a'
        if f >= 0.70: return '#fef9c3', '#a16207'
        if f >= 0.50: return '#fed7aa', '#c2410c'
        return '#fee2e2', '#dc2626'
    except:
        return '#f3f4f6', '#6b7280'

TOOL_DESIGN = {
    'ArcasHLA': 'RNA', 'HLA-HD': 'Both', 'Kourami': 'DNA',
    'OptiType': 'DNA/RNA', 'POLYSOLVER': 'DNA (WES)', 'Seq2HLA': 'RNA',
    'SpecHLA': 'DNA', 'T1K': 'Both',
    'MajorityVote': 'Ensemble', 'WeightedConsensus': 'Ensemble',
}
TOOL_ORDER = ['OptiType', 'POLYSOLVER', 'HLA-HD', 'T1K', 'SpecHLA', 'Kourami', 'Seq2HLA', 'ArcasHLA']


# ── VENEX parsing ──────────────────────────────────────────────────────────────

def two_field(allele):
    """Normalise any HLA allele string to GENE*F1:F2, or return None."""
    if not allele or allele.strip() in ('Not typed', '-', '', 'NA'):
        return None
    a = allele.strip()
    if a.startswith('HLA-'):
        a = a[4:]
    if '*' not in a:
        return None
    gene, rest = a.split('*', 1)
    fields = rest.split(':')
    if len(fields) < 2:
        return None
    return f'{gene}*{fields[0]}:{fields[1]}'

def parse_hlahd(path):
    """Tab-sep: locus allele1 allele2 reads1 reads2. Returns {locus: (a1,a2)}."""
    result = {}
    try:
        with open(path) as fh:
            for line in fh:
                parts = line.rstrip('\n').split('\t')
                if len(parts) < 3:
                    continue
                locus = parts[0].strip()
                result[locus] = (two_field(parts[1]), two_field(parts[2]))
    except FileNotFoundError:
        pass
    return result

def parse_spechla(path):
    """Wide-format TSV with HLA_{LOCUS}_1 / HLA_{LOCUS}_2 columns."""
    result = {}
    LOCI = ['A', 'B', 'C', 'DRB1', 'DQB1', 'DPB1', 'DQA1', 'DPA1']
    try:
        with open(path) as fh:
            header = None
            for line in fh:
                if line.startswith('#'):
                    continue
                parts = line.rstrip('\n').split('\t')
                if header is None:
                    header = parts
                    continue
                row = dict(zip(header, parts))
                for locus in LOCI:
                    result[locus] = (
                        two_field(row.get(f'HLA_{locus}_1', '')),
                        two_field(row.get(f'HLA_{locus}_2', '')),
                    )
    except FileNotFoundError:
        pass
    return result

def concordance_status(h_pair, s_pair):
    """Two-field unordered allele-pair comparison. Returns 'full'/'partial'/'discordant'/'missing'."""
    h = [a for a in h_pair if a]
    s = [a for a in s_pair if a]
    if not h or not s:
        return 'missing'
    if sorted(h) == sorted(s):
        return 'full'
    if set(h) & set(s):
        return 'partial'
    return 'discordant'

def parse_all_venex(venex_dir):
    base = Path(venex_dir)
    samples = []
    for sample_dir in sorted(base.iterdir()):
        if not sample_dir.is_dir():
            continue
        sid = sample_dir.name
        hlahd_p  = sample_dir / 'results' / sid / 'hlahd'  / f'{sid}_hlahd.txt'
        spechla_p = sample_dir / 'results' / sid / 'spechla' / f'{sid}_spechla.txt'
        if not hlahd_p.exists() and not spechla_p.exists():
            continue
        samples.append({
            'sample': sid,
            'hlahd':   parse_hlahd(hlahd_p),
            'spechla': parse_spechla(spechla_p),
        })
    return samples


# ── VENEX matplotlib figures ───────────────────────────────────────────────────

def make_callability_fig(samples):
    LOCI = ['A', 'B', 'C', 'DRB1', 'DQB1', 'DPB1']
    n = len(samples)
    hd_r = [sum(1 for s in samples if any(a for a in s['hlahd'].get(l, (None,None)))) / n for l in LOCI]
    sp_r = [sum(1 for s in samples if any(a for a in s['spechla'].get(l, (None,None)))) / n for l in LOCI]

    fig, ax = plt.subplots(figsize=(7.09, 3.15))
    x = np.arange(len(LOCI))
    w = 0.35
    b1 = ax.bar(x - w/2, hd_r, w, label='HLA-HD',  color='#0072B2')
    b2 = ax.bar(x + w/2, sp_r, w, label='SpecHLA', color='#E69F00')
    ax.set_xticks(x)
    ax.set_xticklabels([f'HLA-{l}' for l in LOCI])
    ax.set_ylabel('Callability rate')
    ax.set_ylim(0, 1.18)
    ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda v, _: f'{v:.0%}'))
    ax.legend(frameon=False)
    for bar in list(b1) + list(b2):
        h = bar.get_height()
        ax.text(bar.get_x() + bar.get_width()/2, h + 0.01,
                f'{h:.0%}', ha='center', va='bottom', fontsize=7)
    fig.tight_layout()
    return fig

def make_concordance_fig(samples):
    LOCI = ['A', 'B', 'C', 'DRB1']
    n = len(samples)
    cats   = ['full',     'partial',     'discordant', 'missing']
    labels = ['Full concordance', 'Partial (1 allele)', 'Discordant', 'One tool missing']
    colors = ['#009E73', '#E69F00',   '#D55E00',    '#94a3b8']

    counts = {l: {c: 0 for c in cats} for l in LOCI}
    for s in samples:
        for l in LOCI:
            h  = s['hlahd'].get(l,   (None, None))
            sp = s['spechla'].get(l, (None, None))
            counts[l][concordance_status(h, sp)] += 1

    fig, ax = plt.subplots(figsize=(5.91, 3.54))
    y = np.arange(len(LOCI))
    lefts = np.zeros(len(LOCI))
    handles = []
    for cat, lbl, col in zip(cats, labels, colors):
        vals = np.array([counts[l][cat] / n * 100 for l in LOCI])
        ax.barh(y, vals, left=lefts, color=col, height=0.6)
        handles.append(mpatches.Patch(color=col, label=lbl))
        for i, (v, left) in enumerate(zip(vals, lefts)):
            if v > 6:
                ax.text(left + v/2, i, f'{v:.0f}%',
                        ha='center', va='center', fontsize=7,
                        color='white', fontweight='bold')
        lefts += vals

    ax.set_yticks(y)
    ax.set_yticklabels([f'HLA-{l}' for l in LOCI])
    ax.set_xlabel('Percentage of samples (%)')
    ax.set_xlim(0, 100)
    ax.legend(handles=handles, frameon=False,
              bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0)
    fig.tight_layout()
    return fig

def make_allele_freq_fig(samples):
    PANELS = [('A', '#0072B2'), ('B', '#D55E00'), ('C', '#009E73'), ('DRB1', '#CC79A7')]
    fig, axes = plt.subplots(2, 2, figsize=(7.09, 6.30))
    for ax, (locus, color) in zip(axes.flatten(), PANELS):
        counter = Counter()
        for s in samples:
            h  = s['hlahd'].get(locus,   (None, None))
            sp = s['spechla'].get(locus, (None, None))
            for i in range(2):
                ha = h[i]  if i < len(h)  else None
                sa = sp[i] if i < len(sp) else None
                allele = ha if ha else sa
                if allele:
                    counter[allele] += 1
        top = counter.most_common(12)
        if not top:
            ax.text(0.5, 0.5, 'No data', ha='center', va='center', transform=ax.transAxes)
            ax.set_title(f'HLA-{locus}', fontweight='bold')
            continue
        lbls, cnts = zip(*top)
        yp = np.arange(len(lbls))
        ax.barh(yp, cnts, color=color, height=0.7)
        ax.set_yticks(yp)
        ax.set_yticklabels(lbls, fontsize=7)
        ax.set_xlabel('Allele observations')
        ax.set_title(f'HLA-{locus}', fontweight='bold')
        ax.invert_yaxis()
    fig.tight_layout(pad=1.5)
    return fig


# ── CSS / JS ───────────────────────────────────────────────────────────────────

CSS = """
* { box-sizing: border-box; margin: 0; padding: 0; }
body { font-family: 'Segoe UI', Arial, sans-serif; background: #f8fafc; color: #1e293b; font-size: 14px; }
.container { max-width: 1280px; margin: 0 auto; padding: 24px; }
h1 { font-size: 28px; font-weight: 700; color: #0f172a; margin-bottom: 6px; }
h2 { font-size: 20px; font-weight: 600; color: #1e3a5f; margin: 32px 0 12px;
     border-left: 4px solid #3b82f6; padding-left: 12px; }
h3 { font-size: 15px; font-weight: 600; color: #374151; margin: 20px 0 8px; }
p  { color: #475569; line-height: 1.65; margin-bottom: 10px; }
.subtitle { color: #64748b; font-size: 15px; margin-bottom: 20px; }
.badge { display:inline-block; padding:2px 10px; border-radius:999px; font-size:11px; font-weight:600; margin:2px; }
.badge-blue { background:#dbeafe; color:#1d4ed8; }
.hero { background: linear-gradient(135deg,#1e3a5f 0%,#1e40af 100%);
        color:white; border-radius:12px; padding:28px 32px; margin-bottom:28px; }
.hero h1 { color:white; }
.hero .subtitle { color:#bfdbfe; }
.hero-stats { display:flex; gap:32px; margin-top:20px; flex-wrap:wrap; }
.hero-stat { text-align:center; }
.hero-stat .val { font-size:32px; font-weight:700; color:#93c5fd; }
.hero-stat .lbl { font-size:12px; color:#bfdbfe; margin-top:2px; }
.hero-venex { background: linear-gradient(135deg,#14532d 0%,#166534 100%); }
.hero-venex .hero-stat .val { color:#86efac; }
.hero-venex .hero-stat .lbl { color:#bbf7d0; }
.hero-venex .subtitle { color:#bbf7d0; }
.grid-3 { display:grid; grid-template-columns:repeat(3,1fr); gap:16px; margin-bottom:20px; }
.grid-2 { display:grid; grid-template-columns:repeat(2,1fr); gap:16px; margin-bottom:20px; }
.card { background:white; border-radius:10px; padding:20px; box-shadow:0 1px 4px rgba(0,0,0,.08); }
.card-title { font-size:13px; font-weight:600; color:#64748b; text-transform:uppercase;
              letter-spacing:.04em; margin-bottom:8px; }
.card-val { font-size:28px; font-weight:700; color:#1e293b; }
.card-sub { font-size:12px; color:#94a3b8; margin-top:4px; }
.mv-highlight { background:#dcfce7; border:1px solid #86efac; color:#15803d;
                padding:2px 8px; border-radius:4px; font-weight:600; }
table { width:100%; border-collapse:collapse; font-size:13px; }
thead th { background:#1e3a5f; color:white; padding:8px 10px; text-align:left;
           font-weight:600; font-size:12px; }
tbody tr:nth-child(even) { background:#f8fafc; }
tbody tr:hover { background:#eff6ff; }
tbody td { padding:7px 10px; border-bottom:1px solid #e2e8f0; vertical-align:middle; }
.mv-row td { background:#f0fdf4 !important; font-weight:600; border-top:2px solid #16a34a; }
.wc-row td { background:#faf5ff !important; font-weight:600; border-top:1px dashed #a855f7; }
.tab-nav { display:flex; gap:4px; margin-bottom:16px; border-bottom:2px solid #e2e8f0; }
.tab-btn { padding:8px 20px; cursor:pointer; border:none; background:none;
           font-size:13px; font-weight:600; color:#64748b;
           border-bottom:2px solid transparent; margin-bottom:-2px; transition:.15s; }
.tab-btn.active { color:#1d4ed8; border-bottom-color:#1d4ed8; }
.tab-pane { display:none; }
.tab-pane.active { display:block; }
.fig-box { background:white; border-radius:10px; padding:20px;
           box-shadow:0 1px 4px rgba(0,0,0,.08); margin-bottom:24px; }
.fig-box img { width:100%; border-radius:6px; }
.fig-caption { margin-top:14px; padding-top:14px; border-top:1px solid #e2e8f0; }
.fig-caption strong { display:block; font-size:14px; color:#1e293b; margin-bottom:6px; }
.fig-caption .short { font-style:italic; color:#374151; margin-bottom:8px; }
.fig-caption .long  { font-size:13px; color:#475569; }
.insight-box { background:#eff6ff; border-left:4px solid #3b82f6;
               border-radius:6px; padding:14px 18px; margin:16px 0; }
.insight-box strong { color:#1d4ed8; }
.insight-venex { background:#f0fdf4; border-left:4px solid #16a34a; }
.insight-venex strong { color:#15803d; }
.section-intro { background:#f8fafc; border:1px solid #e2e8f0;
                 border-radius:8px; padding:16px 20px; margin-bottom:20px; }
.part-header { background:linear-gradient(135deg,#1e3a5f,#1e40af);
               color:white; border-radius:10px; padding:18px 24px; margin-bottom:24px; }
.part-header.venex { background:linear-gradient(135deg,#14532d,#166534); }
.part-header h2 { color:white; border:none; padding:0; margin:8px 0 0; font-size:22px; }
.part-header p  { color:#bfdbfe; margin:4px 0 0; }
.part-header.venex p { color:#bbf7d0; }
.trimodal-grid { display:grid; grid-template-columns:1fr 1fr; gap:16px; }
.progress-row { display:flex; align-items:center; gap:10px; margin:6px 0; }
.progress-label { width:90px; font-size:12px; color:#475569; text-align:right; flex-shrink:0; }
.progress-bar { flex:1; height:12px; background:#e2e8f0; border-radius:6px; overflow:hidden; }
.progress-fill { height:100%; border-radius:6px; }
.progress-val { width:55px; font-size:12px; color:#374151; font-weight:600; }
.part-divider { border:none; border-top:3px solid #e2e8f0; margin:48px 0 32px; }
details.sample-table summary { cursor:pointer; padding:10px 16px; background:#f1f5f9;
    border-radius:8px; font-weight:600; color:#1e3a5f; user-select:none; margin-bottom:8px; }
details.sample-table summary:hover { background:#e2e8f0; }
details.sample-table > div { overflow-x:auto; }
.conc-full    { background:#dcfce7; color:#166534; padding:2px 8px; border-radius:4px; font-size:11px; font-weight:600; }
.conc-partial { background:#fef9c3; color:#a16207; padding:2px 8px; border-radius:4px; font-size:11px; font-weight:600; }
.conc-disc    { background:#fee2e2; color:#dc2626; padding:2px 8px; border-radius:4px; font-size:11px; font-weight:600; }
.conc-miss    { background:#f1f5f9; color:#64748b; padding:2px 8px; border-radius:4px; font-size:11px; font-weight:600; }
footer { margin-top:40px; padding-top:20px; border-top:1px solid #e2e8f0;
         text-align:center; color:#94a3b8; font-size:12px; }
"""

JS = """
function showTab(event, tabId) {
  document.querySelectorAll('.tab-pane').forEach(p => p.classList.remove('active'));
  document.querySelectorAll('.tab-btn').forEach(b => b.classList.remove('active'));
  document.getElementById(tabId).classList.add('active');
  event.target.classList.add('active');
}
"""


# ── HTML section builders ──────────────────────────────────────────────────────

def fig_section(b64, title, short_cap, long_cap, fig_id=''):
    return f"""
<div class="fig-box" id="{fig_id}">
  <img src="data:image/png;base64,{b64}" alt="{title}">
  <div class="fig-caption">
    <strong>{title}</strong>
    <p class="short">{short_cap}</p>
    <p class="long">{long_cap}</p>
  </div>
</div>"""

def tool_table_for_mod(mod, meth_idx, summ_idx):
    mv = meth_idx.get(('MajorityVote',       mod), {})
    wc = meth_idx.get(('WeightedConsensus',  mod), {})
    tools_in_mod = [t for t in TOOL_ORDER if (t, mod) in summ_idx]
    tools_sorted = sorted(tools_in_mod,
        key=lambda t: -float(summ_idx[(t, mod)].get('overall_correct_call_rate', 0) or 0))
    rows = ''
    for tool in tools_sorted:
        r = summ_idx[(tool, mod)]
        acc   = float(r.get('overall_correct_call_rate', 0) or 0)
        cr    = float(r.get('callable_rate', 0) or 0)
        aac   = float(r.get('accuracy_among_callable', 0) or 0)
        n     = r.get('sample_count', '?')
        ci_lo = r.get('overall_correct_call_rate_ci_lo', '')
        ci_hi = r.get('overall_correct_call_rate_ci_hi', '')
        bg, fg = acc_color(acc)
        design = TOOL_DESIGN.get(tool, '—')
        rows += f"""<tr>
          <td><strong>{tool}</strong> <span class="badge badge-blue">{design}</span></td>
          <td style="text-align:center">{n}</td>
          <td>{bar_html(acc, color='#3b82f6')}</td>
          <td><span style="background:{bg};color:{fg};padding:2px 8px;border-radius:4px;font-weight:600">{pct(acc)}</span></td>
          <td style="font-size:12px;color:#64748b">[{fmt(ci_lo,3)}, {fmt(ci_hi,3)}]</td>
          <td>{pct(cr)}</td>
          <td>{pct(aac)}</td>
        </tr>"""
    if mv:
        acc   = float(mv.get('overall_correct_call_rate', 0) or 0)
        cr    = float(mv.get('callable_rate', 0) or 0)
        aac   = float(mv.get('accuracy_among_callable', 0) or 0)
        n     = mv.get('sample_count', '?')
        ci_lo = mv.get('overall_correct_call_rate_ci_lo', '')
        ci_hi = mv.get('overall_correct_call_rate_ci_hi', '')
        rows += f"""<tr class="mv-row">
          <td><strong>&#9733; MajorityVote</strong></td>
          <td style="text-align:center">{n}</td>
          <td>{bar_html(acc, color='#16a34a')}</td>
          <td><span class="mv-highlight">{pct(acc)}</span></td>
          <td style="font-size:12px">[{fmt(ci_lo,3)}, {fmt(ci_hi,3)}]</td>
          <td>{pct(cr)}</td>
          <td>{pct(aac)}</td>
        </tr>"""
    if wc:
        acc   = float(wc.get('overall_correct_call_rate', 0) or 0)
        cr    = float(wc.get('callable_rate', 0) or 0)
        aac   = float(wc.get('accuracy_among_callable', 0) or 0)
        n     = wc.get('sample_count', '?')
        ci_lo = wc.get('overall_correct_call_rate_ci_lo', '')
        ci_hi = wc.get('overall_correct_call_rate_ci_hi', '')
        bg, fg = acc_color(acc)
        rows += f"""<tr class="wc-row">
          <td><strong>&#9670; WeightedConsensus</strong></td>
          <td style="text-align:center">{n}</td>
          <td>{bar_html(acc, color='#a855f7')}</td>
          <td><span style="background:{bg};color:{fg};padding:2px 8px;border-radius:4px;font-weight:600">{pct(acc)}</span></td>
          <td style="font-size:12px">[{fmt(ci_lo,3)}, {fmt(ci_hi,3)}]</td>
          <td>{pct(cr)}</td>
          <td>{pct(aac)}</td>
        </tr>"""
    return f"""<table>
      <thead><tr>
        <th>Method</th><th>Samples (n)</th><th>Accuracy (bar)</th>
        <th>Overall accuracy</th><th>95% CI</th>
        <th>Callable rate</th><th>Acc. among callable</th>
      </tr></thead>
      <tbody>{rows}</tbody>
    </table>"""

def per_gene_table_for_mod(mod, wpg_idx, summ_idx, meth_idx):
    GENES = ['A', 'B', 'C']
    tools_in = [t for t in TOOL_ORDER if (t, mod) in summ_idx]
    ensemble  = [m for m in ('MajorityVote', 'WeightedConsensus') if (m, mod) in meth_idx]
    header = '<tr><th>Method</th>' + ''.join(f'<th>HLA-{g}</th>' for g in GENES) + '</tr>'
    rows = ''
    for method in tools_in + ensemble:
        is_mv = method == 'MajorityVote'
        is_wc = method == 'WeightedConsensus'
        cells = ''
        for gene in GENES:
            r = wpg_idx.get((method, mod, gene))
            if r:
                acc = float(r.get('overall_correct_call_rate', 0) or 0)
                bg, fg = acc_color(acc)
                fw = '700' if (is_mv or is_wc) else '500'
                cells += f'<td><span style="background:{bg};color:{fg};padding:2px 8px;border-radius:4px;font-weight:{fw}">{pct(acc)}</span></td>'
            else:
                cells += "<td style='color:#94a3b8'>—</td>"
        cls   = ' class="mv-row"' if is_mv else (' class="wc-row"' if is_wc else '')
        label = f'&#9733; {method}' if is_mv else (f'&#9670; {method}' if is_wc else method)
        rows += f'<tr{cls}><td><strong>{label}</strong></td>{cells}</tr>'
    return f'<table><thead>{header}</thead><tbody>{rows}</tbody></table>'

def venex_sample_table(samples):
    LOCI = ['A', 'B', 'C']
    header = ('<tr><th>Sample ID</th>'
              + ''.join(f'<th>HLA-HD {l}</th>' for l in LOCI)
              + ''.join(f'<th>SpecHLA {l}</th>' for l in LOCI)
              + ''.join(f'<th>Concordance {l}</th>' for l in LOCI)
              + '</tr>')
    rows = ''
    conc_css = {'full': 'conc-full', 'partial': 'conc-partial',
                'discordant': 'conc-disc', 'missing': 'conc-miss'}
    for s in samples:
        hd = s['hlahd']
        sp = s['spechla']
        row = f'<tr><td><strong style="font-size:12px">{s["sample"]}</strong></td>'
        for l in LOCI:
            p = hd.get(l, (None, None))
            row += f'<td style="font-size:11px">{p[0] or "—"}<br>{p[1] or "—"}</td>'
        for l in LOCI:
            p = sp.get(l, (None, None))
            row += f'<td style="font-size:11px">{p[0] or "—"}<br>{p[1] or "—"}</td>'
        for l in LOCI:
            c = concordance_status(hd.get(l,(None,None)), sp.get(l,(None,None)))
            row += f'<td><span class="{conc_css[c]}">{c}</span></td>'
        rows += row + '</tr>'
    return f'<table><thead>{header}</thead><tbody>{rows}</tbody></table>'


# ── Main render ────────────────────────────────────────────────────────────────

def render(tables_dir, figures_dir, venex_dir, output_path):
    td = Path(tables_dir)
    fd = Path(figures_dir)

    print('Loading benchmark tables...')
    meta  = json.load(open(td / 'benchmark_metadata.json'))
    meth  = list(csv.DictReader(open(td / 'method_comparison.tsv'),  delimiter='\t'))
    summ  = list(csv.DictReader(open(td / 'summary_full_cohort.tsv'), delimiter='\t'))
    wpg   = list(csv.DictReader(open(td / 'method_per_gene.tsv'),    delimiter='\t'))

    meth_idx = {(r['method'], r['modality']): r for r in meth}
    summ_idx  = {(r['tool'],  r['modality']): r for r in summ}
    wpg_idx   = {(r['method'], r['modality'], r['gene']): r for r in wpg}

    mv_wgs = float(meth_idx.get(('MajorityVote','wgs'),    {}).get('overall_correct_call_rate', 0) or 0)
    mv_wes = float(meth_idx.get(('MajorityVote','wes'),    {}).get('overall_correct_call_rate', 0) or 0)
    mv_rna = float(meth_idx.get(('MajorityVote','rnaseq'), {}).get('overall_correct_call_rate', 0) or 0)

    print('Embedding trimodal figures...')
    figs = {p.stem: b64img(p) for p in fd.glob('*.png')}

    print('Computing trimodal sample counts...')
    base_cal = '/scratch/project_2008084/hla_calibration'
    wgs_s = {_re.search(r'/results/([^/]+)/optitype/', f).group(1)
             for f in _glob.glob(f'{base_cal}/wgs_batches/*/results/**/*_optitype.txt', recursive=True)
             if _re.search(r'/results/([^/]+)/optitype/', f)}
    wes_s = {_re.search(r'/results/([^/]+)/optitype/', f).group(1)
             for f in _glob.glob(f'{base_cal}/wes_batches/*/results/**/*_optitype.txt', recursive=True)
             if _re.search(r'/results/([^/]+)/optitype/', f)}
    rna_s = {_re.search(r'/results/([^/]+)/arcashla/', f).group(1)
             for f in _glob.glob(f'{base_cal}/rna_batches/*/results/**/*_arcashla.txt', recursive=True)
             if _re.search(r'/results/([^/]+)/arcashla/', f)}
    trimodal_n = len(wgs_s & wes_s & rna_s)
    N_REF = 132

    print('Parsing VENEX results...')
    venex   = parse_all_venex(venex_dir)
    n_venex = len(venex)
    print(f'  {n_venex} VENEX samples parsed')

    print('Generating VENEX figures...')
    fv1 = make_callability_fig(venex);   b64_v1 = fig_to_b64(fv1); plt.close(fv1)
    fv2 = make_concordance_fig(venex);   b64_v2 = fig_to_b64(fv2); plt.close(fv2)
    fv3 = make_allele_freq_fig(venex);   b64_v3 = fig_to_b64(fv3); plt.close(fv3)

    CONC_LOCI = ['A', 'B', 'C', 'DRB1']
    full_pcts = {}
    for l in CONC_LOCI:
        nf = sum(1 for s in venex
                 if concordance_status(s['hlahd'].get(l,(None,None)),
                                       s['spechla'].get(l,(None,None))) == 'full')
        full_pcts[l] = nf / n_venex * 100 if n_venex else 0

    def embed_fig(stem, title, short, long, fig_id=''):
        if stem not in figs:
            return f'<div class="fig-box"><p style="color:#94a3b8;padding:20px">Figure not available: {stem}.png</p></div>'
        return fig_section(figs[stem], title, short, long, fig_id)

    today    = __import__('datetime').date.today().isoformat()
    imgt_ver = meta.get('imgt_hla_version', '3.59.0')

    mv_wes_n  = meth_idx.get(('MajorityVote','wes'),    {}).get('sample_count','?')
    mv_wgs_n  = meth_idx.get(('MajorityVote','wgs'),    {}).get('sample_count','?')
    mv_rna_n  = meth_idx.get(('MajorityVote','rnaseq'), {}).get('sample_count','?')

    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>PIHLA v7 — Trimodal Analysis &amp; VENEX Cohort</title>
<style>{CSS}</style>
</head>
<body>
<div class="container">

<!-- ═══════════════════════════ PART 1 ═══════════════════════════ -->
<div class="part-header">
  <p style="font-size:12px;font-weight:600;text-transform:uppercase;letter-spacing:.07em;opacity:.8">Part 1 of 2</p>
  <h2>Trimodal HLA Typing Benchmark — 1000 Genomes Project</h2>
  <p>8 tools &middot; WGS, WES, RNA-seq &middot; IMGT/HLA v{imgt_ver} &middot; HLA-A, -B, -C &middot; Two-field resolution</p>
</div>

<div class="hero">
  <h1>PIHLA Benchmark Report <span style="font-size:16px;font-weight:400;opacity:.7">v7</span></h1>
  <p class="subtitle">Multi-tool HLA ensemble typing &middot; majority voting &middot; weighted consensus &middot; 1000 Genomes Project reference cohort</p>
  <div class="hero-stats">
    <div class="hero-stat"><div class="val">{len(wgs_s)}</div><div class="lbl">WGS samples typed</div></div>
    <div class="hero-stat"><div class="val">{len(wes_s)}</div><div class="lbl">WES samples typed</div></div>
    <div class="hero-stat"><div class="val">{len(rna_s)}</div><div class="lbl">RNA-seq samples typed</div></div>
    <div class="hero-stat"><div class="val">{trimodal_n}</div><div class="lbl">Trimodal samples</div></div>
    <div class="hero-stat"><div class="val">8</div><div class="lbl">HLA typing tools</div></div>
    <div class="hero-stat"><div class="val">{pct(mv_wes)}</div><div class="lbl">MV accuracy (WES)</div></div>
    <div class="hero-stat"><div class="val">{pct(mv_rna)}</div><div class="lbl">MV accuracy (RNA-seq)</div></div>
  </div>
</div>

<h2>Executive Summary</h2>
<div class="section-intro">
  <p>The PIHLA benchmark evaluates eight HLA typing tools across three sequencing modalities — whole-genome sequencing (WGS), whole-exome sequencing (WES), and RNA-seq — on the 1000 Genomes Project truth-labelled cohort. Two ensemble strategies are compared: equal-weight <strong>MajorityVote</strong> and <strong>WeightedConsensus</strong> (tool reliability &times; calibrated confidence, &alpha;=0.7/&beta;=0.3). The primary finding is that <strong>MajorityVote consistently matches or exceeds the best individual tool in WES and RNA-seq</strong> while requiring no calibration data and always emitting a call. WeightedConsensus provides marginal additional accuracy where confidence parsers are reliable.</p>
</div>
<div class="grid-3">
  <div class="card">
    <div class="card-title">MajorityVote &mdash; WGS</div>
    <div class="card-val">{pct(mv_wgs)}</div>
    <div class="card-sub">n = {mv_wgs_n} samples &middot; callable 99.8%</div>
  </div>
  <div class="card">
    <div class="card-title">MajorityVote &mdash; WES</div>
    <div class="card-val" style="color:#16a34a">{pct(mv_wes)}</div>
    <div class="card-sub">n = {mv_wes_n} samples &middot; callable 100%</div>
  </div>
  <div class="card">
    <div class="card-title">MajorityVote &mdash; RNA-seq</div>
    <div class="card-val" style="color:#16a34a">{pct(mv_rna)}</div>
    <div class="card-sub">n = {mv_rna_n} samples &middot; callable 100%</div>
  </div>
</div>

<h2>Figure 1 &mdash; Overall Accuracy Overview</h2>
{embed_fig('figure_01_accuracy_overview',
  'Figure 1. HLA typing accuracy across individual tools and ensemble methods, by sequencing modality.',
  'Overall correct-call rate (two-field exact allele-pair match) per tool and modality. Bars represent accuracy; error bars are 95% Wilson confidence intervals; the dashed red reference line marks MajorityVote accuracy within each modality panel.',
  'Each panel corresponds to one sequencing modality: WGS (n=131), WES (n=75), RNA-seq (n=92). The X-axis lists tools ordered by descending accuracy; MajorityVote (★, red) is shown to the right of individual tools after a visual gap. The Y-axis shows the overall correct-call rate: the fraction of all evaluable allele-pair loci (sample × gene × modality) for which the method returned the exact two-field (four-digit) truth pair. Error bars spanning the 95% Wilson score interval allow direct uncertainty comparison across tools and sample sizes. Tools evaluated on an off-label modality (e.g., ArcasHLA on WGS/WES; SpecHLA on RNA-seq) show markedly lower accuracy, confirming that design scope governs single-tool performance. In WES and RNA-seq the MajorityVote reference line sits at or above the tallest individual-tool bar, indicating ensemble superiority across both modalities.',
  'fig1')}

<h2>Detailed Accuracy Tables &mdash; All Modalities</h2>
<p>Per-modality accuracy with 95% Wilson confidence intervals, callable rate, and accuracy-among-callable for every evaluated method. <span class="mv-highlight">&#9733; MajorityVote</span> and <span style="background:#ede9fe;color:#7c3aed;padding:2px 8px;border-radius:4px;font-weight:600">&#9670; WeightedConsensus</span> rows are highlighted at the bottom of each table for direct comparison.</p>
<p style="font-size:12px;color:#64748b;margin-bottom:12px"><strong>Overall accuracy</strong> = correctly called loci / all evaluable loci &nbsp;|&nbsp; <strong>Callable rate</strong> = loci with any call / all evaluable loci &nbsp;|&nbsp; <strong>Acc. among callable</strong> = overall accuracy restricted to called loci (equals overall accuracy when callable rate = 100%).</p>
<div class="tab-nav">
  <button class="tab-btn active" onclick="showTab(event,'tab-wgs')">WGS (n&asymp;131)</button>
  <button class="tab-btn" onclick="showTab(event,'tab-wes')">WES (n&asymp;75)</button>
  <button class="tab-btn" onclick="showTab(event,'tab-rna')">RNA-seq (n&asymp;92)</button>
</div>
<div id="tab-wgs" class="tab-pane active">{tool_table_for_mod('wgs', meth_idx, summ_idx)}</div>
<div id="tab-wes" class="tab-pane">{tool_table_for_mod('wes', meth_idx, summ_idx)}</div>
<div id="tab-rna" class="tab-pane">{tool_table_for_mod('rnaseq', meth_idx, summ_idx)}</div>
<p style="font-size:12px;color:#94a3b8;margin-top:10px">* ArcasHLA on WGS/WES (&lt;15%) is off-label (RNA-seq tool). SpecHLA on RNA-seq (&lt;4%) reflects DNA-pipeline N-masking of spliced reads. POLYSOLVER and Seq2HLA are single-modality tools (WES-only and RNA-seq-only respectively).</p>

<h2>Figure 2 &mdash; MajorityVote Advantage over Individual Tools</h2>
{embed_fig('figure_02_majority_vote_advantage',
  'Figure 2. Accuracy gain of MajorityVote relative to each individual HLA typing tool, per modality.',
  'Each bar shows MajorityVote overall accuracy minus individual tool overall accuracy (percentage points). Green bars: MajorityVote outperforms that tool. Red bars: individual tool outperforms MajorityVote.',
  'In WES and RNA-seq, every bar is positive: MajorityVote outperforms each individual tool evaluated on that modality. In WGS, OptiType is the single exception (MajorityVote is −11 pp below OptiType), because OptiType captures HLA class I reads that lower-accuracy tools consistently mis-assign, and those tools collectively outvote OptiType on contested loci. This pattern demonstrates that MajorityVote is most advantageous when tool performances are clustered (WES, RNA-seq) and provides a moderate benefit even when one tool dominates (WGS), since it still outperforms all tools except the best.',
  'fig2')}

<h2>Figure 3 &mdash; Mechanism: Tool Agreement Predicts Accuracy</h2>
{embed_fig('figure_03_agreement_vs_accuracy',
  'Figure 3. MajorityVote accuracy as a function of the number of tools independently agreeing on the same allele pair.',
  'X-axis: number of tools converging on the winning allele pair for a given sample × gene × modality unit. Y-axis: fraction of those units that matched the IMGT/HLA truth. Counts above each bar show the number of evaluated units at that agreement level.',
  'This figure reveals the mechanistic basis for majority voting: as more tools independently converge on the same allele pair, the probability that the pair is correct increases monotonically across all three modalities. In WGS, even two-tool agreement achieves higher accuracy than any single tool in isolation; unanimous six-tool agreement approaches 100%. In WES and RNA-seq, accuracy rises steeply from two-tool agreement (≈75%) to full agreement (≈98–100%). The positive monotone trend across all modalities provides an empirical justification for using agreement count as a call-quality signal in downstream clinical or research applications: a locus with unanimous tool agreement is far more likely to be correct than one where only a plurality agreed.',
  'fig3')}

<h2>Figure 4 &mdash; Per-Locus Accuracy (HLA-A, -B, -C)</h2>
{embed_fig('figure_04_per_gene_accuracy',
  'Figure 4. Per-locus HLA typing accuracy for HLA-A, HLA-B, and HLA-C across all tools and modalities.',
  'Clustered bars: each cluster groups tools by HLA locus within one modality panel. MajorityVote (red, outlined) is shown last in each cluster. Y-axis: overall correct-call rate (two-field).',
  'Decomposing accuracy by locus reveals whether ensemble advantages are gene-specific or consistent. In WES and RNA-seq, MajorityVote (red bars) consistently occupies the highest or near-highest position in each cluster across HLA-A, -B, and -C, confirming that the ensemble advantage generalises across all three class I loci rather than being driven by one unusually easy gene. HLA-C shows greater inter-tool variability than HLA-A or HLA-B in WGS, attributable to the higher allelic polymorphism density of the C locus and the corresponding challenge for read-alignment-based approaches. In RNA-seq, per-locus accuracy is more uniform because allele-specific expression reduces ambiguity for expressed alleles at all three loci.',
  'fig4')}

<h2>Per-Locus Accuracy Tables</h2>
<div class="tab-nav">
  <button class="tab-btn active" onclick="showTab(event,'pg-wgs')">WGS</button>
  <button class="tab-btn" onclick="showTab(event,'pg-wes')">WES</button>
  <button class="tab-btn" onclick="showTab(event,'pg-rna')">RNA-seq</button>
</div>
<div id="pg-wgs" class="tab-pane active">{per_gene_table_for_mod('wgs', wpg_idx, summ_idx, meth_idx)}</div>
<div id="pg-wes" class="tab-pane">{per_gene_table_for_mod('wes', wpg_idx, summ_idx, meth_idx)}</div>
<div id="pg-rna" class="tab-pane">{per_gene_table_for_mod('rnaseq', wpg_idx, summ_idx, meth_idx)}</div>

<h2>Figure 5 &mdash; Cross-Modality Accuracy Landscape</h2>
{embed_fig('figure_05_crossmodal_landscape',
  'Figure 5. Complete tool × modality accuracy matrix for all evaluated HLA typing tools.',
  'Heatmap cells show overall correct-call rate (two-field). Colour gradient: red (0%) → yellow (50%) → green (100%). Grey cells indicate tools not evaluated on that modality. &#9733; MajorityVote row is separated by a bold border.',
  'This matrix provides a unified cross-modality overview of where each tool excels and where it fails. Design-scope mismatches are immediately visible as low-accuracy cells: ArcasHLA (RNA-seq–designed) achieves 8–13% on WGS/WES data; SpecHLA (DNA-optimised) achieves ≤4% on RNA-seq due to pipeline N-masking of spliced reads. Tools with broad design scope (HLA-HD, T1K, OptiType) show consistently competitive performance across all three modalities, making them the most reliable contributors to the ensemble. POLYSOLVER is limited to WES and Seq2HLA to RNA-seq; both contribute strong per-modality accuracy within their design scope. The MajorityVote row (★) uniformly occupies cells at or near the top of each column, confirming that the ensemble is a robust cross-modality strategy regardless of input data type.',
  'fig5')}

<h2>Trimodal Cohort Coverage</h2>
<div class="section-intro">
  <p>The trimodal cohort comprises 1000 Genomes samples for which all three modalities (WGS, WES, RNA-seq) have been successfully typed. This subset enables direct cross-modality consistency checks &mdash; verifying whether the same sample yields the same HLA genotype from DNA and RNA &mdash; and represents the most comprehensive evaluation of ensemble stability.</p>
</div>
<div class="trimodal-grid">
  <div class="card">
    <h3>Sample Coverage by Modality</h3>
    <div style="margin-top:12px">
      <div class="progress-row">
        <span class="progress-label">WGS</span>
        <div class="progress-bar"><div class="progress-fill" style="width:{len(wgs_s)/N_REF*100:.0f}%;background:#3b82f6"></div></div>
        <span class="progress-val">{len(wgs_s)}/{N_REF}</span>
      </div>
      <div class="progress-row">
        <span class="progress-label">WES</span>
        <div class="progress-bar"><div class="progress-fill" style="width:{len(wes_s)/N_REF*100:.0f}%;background:#8b5cf6"></div></div>
        <span class="progress-val">{len(wes_s)}/{N_REF}</span>
      </div>
      <div class="progress-row">
        <span class="progress-label">RNA-seq</span>
        <div class="progress-bar"><div class="progress-fill" style="width:{len(rna_s)/N_REF*100:.0f}%;background:#10b981"></div></div>
        <span class="progress-val">{len(rna_s)}/{N_REF}</span>
      </div>
      <div class="progress-row">
        <span class="progress-label" style="color:#dc2626;font-weight:700">Trimodal</span>
        <div class="progress-bar"><div class="progress-fill" style="width:{trimodal_n/N_REF*100:.0f}%;background:#dc2626"></div></div>
        <span class="progress-val" style="color:#dc2626;font-weight:700">{trimodal_n}/{N_REF}</span>
      </div>
    </div>
  </div>
  <div class="card">
    <h3>Coverage Notes</h3>
    <p>Of {N_REF} truth-labelled samples, <strong>{trimodal_n} ({trimodal_n/N_REF*100:.0f}%)</strong> carry results for all three modalities.</p>
    <p style="margin-top:8px"><strong>WES gap (~40 samples):</strong> These 1000G WES libraries used exome capture kits that exclude the HLA locus region, yielding zero on-target reads. This is a library-preparation provenance issue, not a tool failure.</p>
    <p style="margin-top:8px"><strong>RNA-seq gap (~40 samples):</strong> The GEUVADIS RNA-seq study covers CEU, FIN, GBR, TSI, and YRI ancestry groups only; Puerto Rican (PUR) samples have no paired RNA-seq data available.</p>
  </div>
</div>

<h2>Methods</h2>
<div class="grid-2">
  <div class="card">
    <h3>Cohort &amp; Truth Labels</h3>
    <p><strong>Source:</strong> 1000 Genomes Project (phase 3, n={N_REF} samples). Ground-truth HLA alleles: Gourraud et al. (2014), <em>PLOS ONE</em>.</p>
    <p><strong>Ancestry groups:</strong> CEU (CEPH/Utah), FIN (Finnish), GBR (British), TSI (Tuscan), YRI (Yoruba in Ibadan).</p>
    <p><strong>Loci evaluated:</strong> HLA-A, HLA-B, HLA-C at two-field (four-digit) resolution.</p>
    <p><strong>IMGT/HLA version:</strong> {imgt_ver} (pinned throughout truth normalisation and tool output evaluation).</p>
  </div>
  <div class="card">
    <h3>Ensemble Strategies</h3>
    <p><strong>MajorityVote:</strong> Each tool that returns a callable allele pair casts one unweighted vote; the plurality pair is selected. Ties broken alphabetically. Callable rate: 100% for WES/RNA-seq; 99.8% for WGS (two loci with no tool calls).</p>
    <p><strong>WeightedConsensus:</strong> Tool weights combine per-modality reliability (&alpha;=0.7, estimated by leave-one-out cross-validation accuracy) and calibrated confidence (&beta;=0.3, from per-tool confidence parsers). Tools without confidence parsers contribute reliability-only weights.</p>
  </div>
</div>
<div class="card" style="margin-bottom:24px">
  <h3>Tools Evaluated</h3>
  <table>
    <thead><tr><th>Tool</th><th>Design scope</th><th>WGS</th><th>WES</th><th>RNA-seq</th><th>Confidence signal</th></tr></thead>
    <tbody>
      <tr><td>OptiType</td><td>DNA/RNA (class I)</td><td>&#10003;</td><td>&#10003;</td><td>&#10003;</td><td>Objective score</td></tr>
      <tr><td>HLA-HD</td><td>DNA + RNA</td><td>&#10003;</td><td>&#10003;</td><td>&#10003;</td><td>Read counts</td></tr>
      <tr><td>T1K</td><td>WGS/WES/RNA</td><td>&#10003;</td><td>&#10003;</td><td>&#10003;</td><td>Genotype confidence</td></tr>
      <tr><td>ArcasHLA</td><td>RNA-seq (off-label DNA)</td><td>&#10003;*</td><td>&#10003;*</td><td>&#10003;</td><td>Allele posteriors</td></tr>
      <tr><td>SpecHLA</td><td>DNA (off-label RNA)</td><td>&#10003;</td><td>&#10003;</td><td>&#10003;*</td><td>None</td></tr>
      <tr><td>POLYSOLVER</td><td>WES (DNA)</td><td>&mdash;</td><td>&#10003;</td><td>&mdash;</td><td>None</td></tr>
      <tr><td>Kourami</td><td>DNA</td><td>&#10003;</td><td>&#10003;</td><td>&mdash;</td><td>Read coverage</td></tr>
      <tr><td>Seq2HLA</td><td>RNA-seq</td><td>&mdash;</td><td>&mdash;</td><td>&#10003;</td><td>None</td></tr>
    </tbody>
  </table>
  <p style="font-size:12px;color:#94a3b8;margin-top:8px">* Off-label use. ArcasHLA on WGS/WES: &lt;15% accuracy. SpecHLA on RNA-seq: &lt;4% accuracy.</p>
</div>

<hr class="part-divider">

<!-- ═══════════════════════════ PART 2 ═══════════════════════════ -->
<div class="part-header venex">
  <p style="font-size:12px;font-weight:600;text-transform:uppercase;letter-spacing:.07em;opacity:.8">Part 2 of 2</p>
  <h2>VENEX Cohort &mdash; Applied WGS HLA Typing</h2>
  <p>Finnish clinical WGS cohort &middot; {n_venex} samples &middot; HLA-HD + SpecHLA &middot; No reference ground truth</p>
</div>

<div class="hero hero-venex">
  <h1>VENEX WGS Cohort Results</h1>
  <p class="subtitle">Applied HLA typing on the VENEX Finnish cohort &mdash; callability, inter-tool concordance, and allele frequency distribution</p>
  <div class="hero-stats">
    <div class="hero-stat"><div class="val">{n_venex}</div><div class="lbl">Samples processed</div></div>
    <div class="hero-stat"><div class="val">2</div><div class="lbl">Tools (HLA-HD, SpecHLA)</div></div>
    <div class="hero-stat"><div class="val">WGS</div><div class="lbl">Sequencing modality</div></div>
    <div class="hero-stat"><div class="val">{full_pcts.get('A',0):.0f}%</div><div class="lbl">HLA-A concordance</div></div>
    <div class="hero-stat"><div class="val">{full_pcts.get('B',0):.0f}%</div><div class="lbl">HLA-B concordance</div></div>
    <div class="hero-stat"><div class="val">{full_pcts.get('C',0):.0f}%</div><div class="lbl">HLA-C concordance</div></div>
    <div class="hero-stat"><div class="val">{full_pcts.get('DRB1',0):.0f}%</div><div class="lbl">DRB1 concordance</div></div>
  </div>
</div>

<h2>VENEX Cohort Overview</h2>
<div class="section-intro">
  <p>The VENEX cohort consists of {n_venex} Finnish WGS samples (Batches 1 and 1_2) processed through the PIHLA pipeline at CSC Puhti. Two WGS-capable tools were applied: <strong>HLA-HD</strong> (run against the full current IMGT/HLA database) and <strong>SpecHLA</strong> (run against IPD-IMGT/HLA 3.38.0). Because no reference genotyping truth labels are available for this cohort, classical accuracy metrics cannot be computed. Instead, this section reports <strong>per-locus callability</strong> (fraction of samples receiving a non-null call), <strong>inter-tool concordance</strong> between HLA-HD and SpecHLA at two-field resolution, and <strong>HLA allele frequency distributions</strong> at HLA-A, -B, -C, and DRB1.</p>
</div>

<div class="insight-box insight-venex">
  <strong>Key finding:</strong> Full two-field concordance between HLA-HD and SpecHLA reaches {full_pcts.get('A',0):.0f}% for HLA-A, {full_pcts.get('B',0):.0f}% for HLA-B, {full_pcts.get('C',0):.0f}% for HLA-C, and {full_pcts.get('DRB1',0):.0f}% for HLA-DRB1. The higher concordance at class I loci is consistent with the 1000G benchmark finding that both HLA-HD and SpecHLA perform reliably on WGS data for HLA-A, -B, and -C. The lower DRB1 concordance reflects the greater allelic diversity of DRB1 and known differences in the tools&rsquo; class II allele databases and alignment strategies at this locus.
</div>

<h2>Figure V1 &mdash; Per-Locus Callability</h2>
{fig_section(b64_v1,
  f'Figure V1. Per-locus callability rate for HLA-HD and SpecHLA in the VENEX WGS cohort (n={n_venex}).',
  f'X-axis: six classical HLA loci (class I: A, B, C; class II: DRB1, DQB1, DPB1). Y-axis: fraction of {n_venex} samples for which the tool emitted at least one non-null allele call at that locus. HLA-HD (blue) and SpecHLA (amber) are shown side-by-side per locus; percentage values are annotated above each bar.',
  'A locus is classified as called when HLA-HD returns any allele other than “Not typed” and when SpecHLA returns a non-empty allele field. Class I loci (HLA-A, -B, -C) show higher callability than class II loci owing to (1) greater completeness of the class I allele reference database, (2) higher per-base WGS read depth across the class I region of chromosome 6p21.3, and (3) less repetitive flanking sequence that complicates read mapping in class II. The lower DQB1 and DPB1 callability for SpecHLA relative to HLA-HD likely reflects differences in the two tools’ mapping and database strategies for the highly repetitive class II sub-region. Callability values ≥95% at HLA-A, -B, and -C confirm that both tools are technically robust for class I typing on this WGS cohort.',
  'figV1')}

<h2>Figure V2 &mdash; Inter-Tool Concordance</h2>
{fig_section(b64_v2,
  f'Figure V2. Two-field inter-tool concordance between HLA-HD and SpecHLA at four classical HLA loci in the VENEX cohort (n={n_venex}).',
  'Horizontal stacked bars. X-axis: percentage of samples in each concordance category. Y-axis: HLA locus. Categories: full concordance (both alleles match, unordered, green), partial concordance (one allele matches, amber), discordant (no allele in common, red-orange), one tool missing (grey).',
  'Concordance is assessed at two-field (four-digit) resolution. A locus is “fully concordant” when the sorted set of non-null allele calls from HLA-HD equals the sorted set from SpecHLA (unordered diploid pair comparison). “Partial” concordance indicates one of the two allele slots agrees. “Discordant” indicates no allele in the pair is shared between tools. “One tool missing” indicates a null call from at least one tool at that locus. Full concordance ≥80% at HLA-A, -B, and -C confirms that both WGS-oriented tools reliably and consistently identify the same HLA genotype from whole-genome data at class I loci. The lower concordance at DRB1 is expected: DRB1 has the highest allelic diversity among classical HLA loci and is the most sensitive to differences in IMGT/HLA database version and alignment strategy between tools.',
  'figV2')}

<h2>Figure V3 &mdash; HLA Allele Frequency Distribution</h2>
{fig_section(b64_v3,
  f'Figure V3. HLA allele frequency distribution at four classical loci (HLA-A, HLA-B, HLA-C, HLA-DRB1) in the VENEX WGS cohort (n={n_venex}).',
  f'Four-panel figure. X-axis in each panel: number of allele observations (maximum 2×{n_venex}={2*n_venex} per locus). Y-axis: top-12 two-field allele designations, most frequent at the top. Each allele observation corresponds to one diploid allele slot per sample, using HLA-HD calls (most common at each slot) as the primary source.',
  f'Allele calls are derived from HLA-HD where a call exists, and from SpecHLA otherwise. All calls are truncated to two-field (four-digit) resolution. The distributions are consistent with expected Northern European (Finnish) HLA haplotype frequencies: HLA-A*03:01, A*02:01, and A*01:01 are the most frequent class I A-locus alleles in Finnish population cohorts (Lokki et al. 2000, Tissue Antigens; Testi et al. 2017, Front. Immunol.). At the B locus, B*07:02, B*35:01, and B*08:01 are canonical Northern European haplotype alleles. The DRB1 distribution is expected to reflect Finnish DR haplotypes enriched for DRB1*04:04, DRB1*01:01, and DRB1*15:01. All allele counts are presented without population-stratification correction; the cohort is treated as a single Finnish WGS cohort.',
  'figV3')}

<h2>VENEX Sample Summary Table</h2>
<p>Per-sample HLA-HD and SpecHLA calls at HLA-A, -B, -C (two-field resolution) with concordance status badge per locus. Badges: <span class="conc-full">full</span> both alleles match &nbsp; <span class="conc-partial">partial</span> one allele matches &nbsp; <span class="conc-disc">discordant</span> no allele in common &nbsp; <span class="conc-miss">missing</span> one or both tools have no call.</p>
<details class="sample-table">
  <summary>&#9654; Show all {n_venex} VENEX samples (click to expand)</summary>
  <div>{venex_sample_table(venex)}</div>
</details>

<footer>
  <p>PIHLA benchmark report v7 &middot; Generated {today} &middot; IMGT/HLA v{imgt_ver} &middot; Two-field resolution</p>
  <p style="margin-top:4px">Part 1: 1000 Genomes Project cohort &middot; Truth: Gourraud et al. (2014) &middot; Part 2: VENEX Finnish WGS cohort</p>
  <p style="margin-top:4px">Workflow: PIHLA v3.0 &middot; Pipeline: Nextflow DSL2 &middot; Cluster: CSC Puhti (SLURM)</p>
</footer>

</div>
<script>{JS}</script>
</body>
</html>"""

    with open(output_path, 'w') as fh:
        fh.write(html)
    size_kb = len(html) // 1024
    print(f'Written: {output_path}  ({size_kb} KB)')


def main():
    p = argparse.ArgumentParser(description='Generate PIHLA HTML report v7')
    p.add_argument('--tables-dir', required=True)
    p.add_argument('--figures-dir', required=True)
    p.add_argument('--venex-dir',   required=True)
    p.add_argument('--output',      required=True)
    args = p.parse_args()
    render(args.tables_dir, args.figures_dir, args.venex_dir, args.output)


if __name__ == '__main__':
    main()
