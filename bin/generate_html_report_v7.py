#!/usr/bin/env python3.11
"""
generate_html_report_v7.py
MVHLA HTML report v7 — trimodal benchmark (1000G) + VENEX WGS cohort.
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
:root {
  --bg: #f2f6fb;
  --surface: #ffffff;
  --surface-soft: #f7fbff;
  --line: #d9e3ef;
  --ink: #102033;
  --muted: #5b6674;
  --navy: #0f2744;
  --blue: #2a63b8;
  --blue-deep: #173d74;
  --green: #166534;
  --green-soft: #effcf4;
  --amber: #a16207;
  --amber-soft: #fff8e8;
  --rose: #c2410c;
  --shadow: 0 18px 48px rgba(15, 39, 68, 0.10);
  --radius: 22px;
  --sidebar-width: 300px;
}

* { box-sizing: border-box; margin: 0; padding: 0; }
html { scroll-behavior: smooth; }
body {
  font-family: "Segoe UI", Arial, sans-serif;
  background: linear-gradient(180deg, #edf3f9 0%, #f9fbfd 100%);
  color: var(--ink);
  font-size: 14px;
}

.layout {
  display: grid;
  grid-template-columns: var(--sidebar-width) minmax(0, 1fr);
  min-height: 100vh;
}

aside {
  position: sticky;
  top: 0;
  align-self: start;
  height: 100vh;
  padding: 28px 22px;
  background: linear-gradient(180deg, #102844 0%, #18467a 56%, #215fa8 100%);
  color: #fff;
  border-right: 1px solid rgba(255,255,255,0.10);
}

.brand {
  margin-bottom: 24px;
  padding-bottom: 18px;
  border-bottom: 1px solid rgba(255,255,255,0.14);
}

.kicker {
  margin: 0 0 10px;
  color: #c9ddf4;
  font-size: 12px;
  letter-spacing: 0.08em;
  text-transform: uppercase;
  font-weight: 700;
}

.brand h1 {
  margin: 0 0 10px;
  font-size: 28px;
  line-height: 1.08;
  color: #fff;
}

.brand p {
  color: #dcebfb;
  line-height: 1.5;
}

.sidebar-meta {
  display: grid;
  gap: 12px;
  margin: 18px 0 24px;
}

.sidebar-chip {
  padding: 10px 12px;
  border: 1px solid rgba(255,255,255,0.14);
  border-radius: 16px;
  background: rgba(255,255,255,0.07);
  color: #e6f0fb;
  font-size: 13px;
}

.sidebar-chip strong {
  display: block;
  margin-bottom: 4px;
  color: #fff;
  font-size: 12px;
  text-transform: uppercase;
  letter-spacing: 0.05em;
}

.nav { display: grid; gap: 8px; }

.nav a {
  color: #d8e7f7;
  text-decoration: none;
  padding: 10px 12px;
  border-radius: 14px;
  border: 1px solid transparent;
  font-size: 14px;
  line-height: 1.35;
}

.nav a:hover,
.nav a.active {
  background: rgba(255,255,255,0.10);
  border-color: rgba(255,255,255,0.16);
  color: #fff;
}

main { padding: 30px 32px 48px; }

.hero {
  background: linear-gradient(135deg, #102844 0%, #174576 55%, #2a63b8 100%);
  color: #fff;
  border-radius: 28px;
  box-shadow: var(--shadow);
  padding: 32px 34px;
  margin-bottom: 24px;
}

.hero h1 {
  font-size: 40px;
  line-height: 1.05;
  color: #fff;
  margin-bottom: 10px;
}

.hero p {
  color: #e3eef9;
  line-height: 1.65;
  max-width: 1040px;
}

.eyebrow {
  margin: 0 0 12px;
  font-size: 12px;
  letter-spacing: 0.08em;
  text-transform: uppercase;
  font-weight: 700;
  color: #cbe0f7;
}

.hero-stats {
  display: grid;
  grid-template-columns: repeat(auto-fit, minmax(140px, 1fr));
  gap: 14px;
  margin-top: 22px;
}

.hero-stat {
  padding: 14px 14px 12px;
  border-radius: 18px;
  background: rgba(255,255,255,0.08);
  border: 1px solid rgba(255,255,255,0.09);
}

.hero-stat .val {
  font-size: 28px;
  font-weight: 700;
  color: #9ec9ff;
  line-height: 1.1;
}

.hero-stat .lbl {
  margin-top: 6px;
  font-size: 12px;
  color: #d7e8fa;
}

section { margin-top: 24px; }

.section-header { margin-bottom: 16px; }

.section-header h2 {
  margin: 0 0 8px;
  font-size: 28px;
  color: var(--navy);
}

.section-header p {
  color: var(--muted);
  line-height: 1.65;
  max-width: 980px;
}

.grid { display: grid; gap: 18px; }
.grid-2 { grid-template-columns: repeat(2, minmax(0, 1fr)); }
.grid-3 { grid-template-columns: repeat(3, minmax(0, 1fr)); }
.metrics-grid { grid-template-columns: repeat(3, minmax(0, 1fr)); }
.two-col { grid-template-columns: 1.15fr 0.85fr; align-items: start; }
.trimodal-grid { grid-template-columns: repeat(2, minmax(0, 1fr)); }

.panel, .card, .metric, .fig-box {
  background: linear-gradient(180deg, #fbfdff 0%, #ffffff 100%);
  border: 1px solid var(--line);
  border-radius: var(--radius);
  box-shadow: 0 8px 24px rgba(15, 39, 68, 0.05);
}

.panel, .card, .metric, .fig-box { padding: 22px; }
.panel h3, .card h3 {
  margin: 0 0 12px;
  font-size: 20px;
  color: var(--navy);
}

p { color: var(--ink); line-height: 1.65; }
.subtitle { color: #dfeefe; margin-bottom: 0; }

.metric .card-title,
.metric .label,
.card-title {
  font-size: 12px;
  text-transform: uppercase;
  letter-spacing: 0.06em;
  color: var(--muted);
  font-weight: 700;
  margin-bottom: 8px;
}

.metric .card-val,
.metric .value,
.card-val {
  font-size: 30px;
  font-weight: 700;
  color: var(--navy);
  line-height: 1.1;
  margin-bottom: 8px;
}

.card-sub, .metric .note {
  font-size: 14px;
  line-height: 1.5;
  color: var(--muted);
}

.section-intro,
.callout {
  padding: 18px 20px;
  border-radius: 18px;
  border: 1px solid var(--line);
  background: var(--surface-soft);
}

.callout {
  border-left: 5px solid var(--blue);
  background: #f5f9ff;
}

.callout.venex {
  border-left-color: var(--green);
  background: var(--green-soft);
}

.callout.warn {
  border-left-color: var(--amber);
  background: var(--amber-soft);
}

.badge {
  display: inline-block;
  padding: 2px 10px;
  border-radius: 999px;
  font-size: 11px;
  font-weight: 600;
  margin: 2px;
}

.badge-blue { background: #dbeafe; color: #1d4ed8; }
.mv-highlight {
  background: #dcfce7;
  border: 1px solid #86efac;
  color: #15803d;
  padding: 2px 8px;
  border-radius: 4px;
  font-weight: 600;
}

.tab-nav {
  display: flex;
  flex-wrap: wrap;
  gap: 4px;
  margin-bottom: 16px;
  border-bottom: 2px solid #e2e8f0;
}

.tab-btn {
  padding: 8px 20px;
  cursor: pointer;
  border: none;
  background: none;
  font-size: 13px;
  font-weight: 600;
  color: #64748b;
  border-bottom: 2px solid transparent;
  margin-bottom: -2px;
  transition: .15s;
}

.tab-btn.active {
  color: #1d4ed8;
  border-bottom-color: #1d4ed8;
}

.tab-pane { display: none; }
.tab-pane.active { display: block; }

.table-wrap,
details.sample-table > div {
  overflow-x: auto;
  border: 1px solid var(--line);
  border-radius: 18px;
  background: #fff;
}

table { width: 100%; border-collapse: collapse; font-size: 13px; }
thead th {
  background: #eef5fd;
  color: var(--navy);
  padding: 10px 12px;
  text-align: left;
  font-weight: 700;
  font-size: 12px;
  text-transform: uppercase;
  letter-spacing: 0.05em;
}

tbody tr:nth-child(even) { background: #fbfdff; }
tbody tr:hover { background: #eff6ff; }
tbody td {
  padding: 8px 10px;
  border-bottom: 1px solid #e2e8f0;
  vertical-align: middle;
}

.mv-row td { background: #f0fdf4 !important; font-weight: 600; border-top: 2px solid #16a34a; }
.wc-row td { background: #faf5ff !important; font-weight: 600; border-top: 1px dashed #a855f7; }

.figure-grid { display: grid; gap: 20px; }
.fig-box img { width: 100%; border-radius: 10px; border: 1px solid #e6edf5; }
.fig-caption { margin-top: 14px; padding-top: 14px; border-top: 1px solid #e2e8f0; }
.fig-caption strong {
  display: block;
  font-size: 18px;
  color: var(--navy);
  margin-bottom: 6px;
}
.fig-caption .short { font-style: italic; color: #374151; margin-bottom: 8px; }
.fig-caption .long  { font-size: 13px; color: #475569; }

.meta-list { display: grid; gap: 10px; }
.meta-row {
  display: grid;
  grid-template-columns: 220px 1fr;
  gap: 14px;
  padding: 12px 0;
  border-bottom: 1px solid #e8eef5;
}
.meta-row:last-child { border-bottom: 0; }
.meta-row strong { color: var(--navy); }

.progress-row { display: flex; align-items: center; gap: 10px; margin: 8px 0; }
.progress-label { width: 90px; font-size: 12px; color: #475569; text-align: right; flex-shrink: 0; }
.progress-bar { flex: 1; height: 12px; background: #e2e8f0; border-radius: 6px; overflow: hidden; }
.progress-fill { height: 100%; border-radius: 6px; }
.progress-val { width: 62px; font-size: 12px; color: #374151; font-weight: 600; }

details.sample-table summary {
  cursor: pointer;
  padding: 10px 16px;
  background: #f1f5f9;
  border-radius: 8px;
  font-weight: 600;
  color: var(--navy);
  user-select: none;
  margin-bottom: 8px;
}
details.sample-table summary:hover { background: #e2e8f0; }

.conc-full    { background: #dcfce7; color: #166534; padding: 2px 8px; border-radius: 4px; font-size: 11px; font-weight: 600; }
.conc-partial { background: #fef9c3; color: #a16207; padding: 2px 8px; border-radius: 4px; font-size: 11px; font-weight: 600; }
.conc-disc    { background: #fee2e2; color: #dc2626; padding: 2px 8px; border-radius: 4px; font-size: 11px; font-weight: 600; }
.conc-miss    { background: #f1f5f9; color: #64748b; padding: 2px 8px; border-radius: 4px; font-size: 11px; font-weight: 600; }

footer {
  margin-top: 40px;
  padding-top: 20px;
  border-top: 1px solid #e2e8f0;
  text-align: center;
  color: #94a3b8;
  font-size: 12px;
}

@media (max-width: 1180px) {
  .layout { grid-template-columns: 1fr; }
  aside { position: relative; height: auto; }
  main { padding-top: 20px; }
  .grid-2, .grid-3, .metrics-grid, .two-col, .trimodal-grid { grid-template-columns: 1fr; }
}

@media (max-width: 720px) {
  main { padding: 16px; }
  aside { padding: 20px 16px; }
  .hero { padding: 24px 20px; }
  .hero h1 { font-size: 32px; }
  .section-header h2 { font-size: 24px; }
  .meta-row { grid-template-columns: 1fr; gap: 6px; }
}
"""

JS = """
function showTab(event, tabId) {
  document.querySelectorAll('.tab-pane').forEach(p => p.classList.remove('active'));
  document.querySelectorAll('.tab-btn').forEach(b => b.classList.remove('active'));
  document.getElementById(tabId).classList.add('active');
  event.target.classList.add('active');
}

const navLinks = Array.from(document.querySelectorAll('#sidebar-nav a'));
const sections = navLinks
  .map(link => document.querySelector(link.getAttribute('href')))
  .filter(Boolean);
const byId = new Map(navLinks.map(link => [link.getAttribute('href').slice(1), link]));
const observer = new IntersectionObserver((entries) => {
  entries.forEach(entry => {
    if (!entry.isIntersecting) return;
    navLinks.forEach(link => link.classList.remove('active'));
    const link = byId.get(entry.target.id);
    if (link) link.classList.add('active');
  });
}, { rootMargin: '-30% 0px -55% 0px', threshold: 0.01 });
sections.forEach(section => observer.observe(section));
if (navLinks.length) navLinks[0].classList.add('active');
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
<title>MVHLA v7 — Trimodal Analysis &amp; VENEX Cohort</title>
<style>{CSS}</style>
</head>
<body>
<div class="layout">
  <aside>
    <div class="brand">
      <p class="kicker">MVHLA report</p>
      <h1>MVHLA Benchmark Report</h1>
      <p>Trimodal 1000G benchmark plus VENEX WGS cohort in one standalone offline HTML report.</p>
    </div>
    <div class="sidebar-meta">
      <div class="sidebar-chip"><strong>Generated</strong>{today}</div>
      <div class="sidebar-chip"><strong>Benchmark figures</strong>{fd}</div>
      <div class="sidebar-chip"><strong>Benchmark tables</strong>{td}</div>
      <div class="sidebar-chip"><strong>VENEX cohort</strong>{venex_dir}</div>
    </div>
    <nav class="nav" id="sidebar-nav">
      <a href="#executive-summary">Executive Summary</a>
      <a href="#report-scope-and-provenance">Report Scope and Provenance</a>
      <a href="#trimodal-benchmark-context">Trimodal Benchmark Context</a>
      <a href="#benchmark-performance-tables">Benchmark Performance Tables</a>
      <a href="#benchmark-figure-gallery">Benchmark Figure Gallery</a>
      <a href="#methods-and-tool-coverage">Methods and Tool Coverage</a>
      <a href="#venex-cohort-overview">VENEX Cohort Overview</a>
      <a href="#venex-figure-gallery">VENEX Figure Gallery</a>
      <a href="#venex-sample-summary">VENEX Sample Summary</a>
      <a href="#reproducibility-metadata">Reproducibility Metadata</a>
      <a href="#appendix">Appendix</a>
    </nav>
  </aside>

  <main>
    <section class="hero">
      <p class="eyebrow">Offline standalone HTML benchmark report</p>
      <h1>MVHLA v7</h1>
      <p>Sidebar-based project report for the MVHLA (Majority Voting HLA) framework, covering the 1000 Genomes trimodal HLA typing benchmark and the VENEX Finnish WGS cohort. The report preserves the current v7 computations, embeds benchmark PNG figures directly, generates VENEX figures on the fly, and keeps provenance visible throughout.</p>
      <div class="hero-stats">
        <div class="hero-stat"><div class="val">{len(wgs_s)}</div><div class="lbl">WGS samples typed</div></div>
        <div class="hero-stat"><div class="val">{len(wes_s)}</div><div class="lbl">WES samples typed</div></div>
        <div class="hero-stat"><div class="val">{len(rna_s)}</div><div class="lbl">RNA-seq samples typed</div></div>
        <div class="hero-stat"><div class="val">{trimodal_n}</div><div class="lbl">Trimodal samples</div></div>
        <div class="hero-stat"><div class="val">{n_venex}</div><div class="lbl">VENEX samples</div></div>
        <div class="hero-stat"><div class="val">{pct(mv_wgs)}</div><div class="lbl">MV accuracy WGS</div></div>
        <div class="hero-stat"><div class="val">{pct(mv_wes)}</div><div class="lbl">MV accuracy WES</div></div>
        <div class="hero-stat"><div class="val">{pct(mv_rna)}</div><div class="lbl">MV accuracy RNA-seq</div></div>
      </div>
    </section>

    <section id="executive-summary">
      <div class="section-header">
        <h2>Executive Summary</h2>
        <p>This report combines the v7 benchmark interpretation for the truth-labelled 1000 Genomes trimodal cohort with an applied VENEX WGS cohort readout, using one navigable HTML document for MVHLA.</p>
      </div>
      <div class="section-intro">
        <p>The MVHLA benchmark evaluates eight HLA typing tools across whole-genome sequencing (WGS), whole-exome sequencing (WES), and RNA-seq. Two ensemble strategies are compared: equal-weight <strong>MajorityVote</strong> and <strong>WeightedConsensus</strong> (tool reliability &times; calibrated confidence, &alpha;=0.7/&beta;=0.3). The main benchmark finding remains that <strong>MajorityVote consistently matches or exceeds the best individual tool in WES and RNA-seq</strong> while requiring no calibration data and typically emitting a call. The VENEX section then shows how two WGS-oriented tools behave on a real Finnish cohort without external truth labels.</p>
      </div>
      <div class="grid metrics-grid">
        <div class="metric">
          <div class="label">MajorityVote WGS</div>
          <div class="value">{pct(mv_wgs)}</div>
          <div class="note">n = {mv_wgs_n} samples; callable rate about 99.8%.</div>
        </div>
        <div class="metric">
          <div class="label">MajorityVote WES</div>
          <div class="value">{pct(mv_wes)}</div>
          <div class="note">n = {mv_wes_n} samples; callable rate 100%.</div>
        </div>
        <div class="metric">
          <div class="label">MajorityVote RNA-seq</div>
          <div class="value">{pct(mv_rna)}</div>
          <div class="note">n = {mv_rna_n} samples; callable rate 100%.</div>
        </div>
      </div>
    </section>

    <section id="report-scope-and-provenance">
      <div class="section-header">
        <h2>Report Scope and Provenance</h2>
        <p>The report is intentionally explicit about where each component comes from so interpretation stays grounded in actual MVHLA v7 inputs.</p>
      </div>
      <div class="grid two-col">
        <div class="panel">
          <h3>Source boundaries</h3>
          <div class="meta-list">
            <div class="meta-row"><strong>Tables directory</strong><span>{td}</span></div>
            <div class="meta-row"><strong>Figures directory</strong><span>{fd}</span></div>
            <div class="meta-row"><strong>VENEX directory</strong><span>{venex_dir}</span></div>
            <div class="meta-row"><strong>Generated on</strong><span>{today}</span></div>
            <div class="meta-row"><strong>Benchmark truth context</strong><span>1000 Genomes Project with Gourraud et al. (2014) labels at two-field resolution.</span></div>
          </div>
        </div>
        <div class="panel">
          <h3>Interpretation notes</h3>
          <p>The benchmark portion is truth-labelled and reports classical accuracy metrics. The VENEX portion has no external truth labels, so it reports callability, inter-tool concordance, and allele frequency structure instead. Missing benchmark PNG assets are shown as placeholder panels rather than stopping report generation.</p>
        </div>
      </div>
    </section>

    <section id="trimodal-benchmark-context">
      <div class="section-header">
        <h2>Trimodal Benchmark Context</h2>
        <p>The trimodal cohort comprises 1000 Genomes samples with usable WGS, WES, and RNA-seq outputs, enabling modality-by-modality benchmarking and cross-modality stability checks.</p>
      </div>
      <div class="grid trimodal-grid">
        <div class="panel">
          <h3>Sample coverage by modality</h3>
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
        <div class="panel">
          <h3>Coverage notes</h3>
          <p>Of {N_REF} truth-labelled samples, <strong>{trimodal_n} ({trimodal_n/N_REF*100:.0f}%)</strong> carry results for all three modalities.</p>
          <p><strong>WES gap:</strong> some 1000G exome libraries exclude the HLA locus, yielding zero on-target reads.</p>
          <p><strong>RNA-seq gap:</strong> GEUVADIS coverage is population-limited, so paired RNA-seq does not exist for every truth-labelled sample.</p>
          <p><strong>Resolution:</strong> all benchmark summaries here are interpreted at two-field HLA-A, HLA-B, and HLA-C resolution under IMGT/HLA v{imgt_ver}.</p>
        </div>
      </div>
    </section>

    <section id="benchmark-performance-tables">
      <div class="section-header">
        <h2>Benchmark Performance Tables</h2>
        <p>Per-modality accuracy with Wilson confidence intervals, callable rate, and accuracy-among-callable for every evaluated method.</p>
      </div>
      <div class="callout">
        <p><strong>How to read these tables.</strong> <strong>Overall accuracy</strong> means correctly called loci divided by all evaluable loci. <strong>Callable rate</strong> means loci with any call divided by all evaluable loci. <strong>Accuracy among callable</strong> is overall accuracy restricted to called loci. <span class="mv-highlight">&#9733; MajorityVote</span> and the purple &#9670; <strong>WeightedConsensus</strong> rows are highlighted for quick comparison.</p>
      </div>
      <div class="panel" style="margin-top:18px">
        <div class="tab-nav">
          <button class="tab-btn active" onclick="showTab(event,'tab-wgs')">WGS (n&asymp;131)</button>
          <button class="tab-btn" onclick="showTab(event,'tab-wes')">WES (n&asymp;75)</button>
          <button class="tab-btn" onclick="showTab(event,'tab-rna')">RNA-seq (n&asymp;92)</button>
        </div>
        <div id="tab-wgs" class="tab-pane active"><div class="table-wrap">{tool_table_for_mod('wgs', meth_idx, summ_idx)}</div></div>
        <div id="tab-wes" class="tab-pane"><div class="table-wrap">{tool_table_for_mod('wes', meth_idx, summ_idx)}</div></div>
        <div id="tab-rna" class="tab-pane"><div class="table-wrap">{tool_table_for_mod('rnaseq', meth_idx, summ_idx)}</div></div>
      </div>
      <div class="panel" style="margin-top:18px">
        <h3>Per-locus accuracy tables</h3>
        <div class="tab-nav">
          <button class="tab-btn active" onclick="showTab(event,'pg-wgs')">WGS</button>
          <button class="tab-btn" onclick="showTab(event,'pg-wes')">WES</button>
          <button class="tab-btn" onclick="showTab(event,'pg-rna')">RNA-seq</button>
        </div>
        <div id="pg-wgs" class="tab-pane active"><div class="table-wrap">{per_gene_table_for_mod('wgs', wpg_idx, summ_idx, meth_idx)}</div></div>
        <div id="pg-wes" class="tab-pane"><div class="table-wrap">{per_gene_table_for_mod('wes', wpg_idx, summ_idx, meth_idx)}</div></div>
        <div id="pg-rna" class="tab-pane"><div class="table-wrap">{per_gene_table_for_mod('rnaseq', wpg_idx, summ_idx, meth_idx)}</div></div>
      </div>
      <p style="font-size:12px;color:#64748b;margin-top:12px">ArcasHLA on WGS/WES and SpecHLA on RNA-seq remain off-label comparisons and should be read as design-scope stress tests rather than intended production modalities.</p>
    </section>

    <section id="benchmark-figure-gallery">
      <div class="section-header">
        <h2>Benchmark Figure Gallery</h2>
        <p>The current MVHLA v7 benchmark narrative is preserved here as a navigable figure gallery. Each panel embeds a PNG directly when available and otherwise shows a graceful placeholder.</p>
      </div>
      <div class="callout">
        <p><strong>Benchmark reading guide.</strong> These figures focus on modality-specific accuracy, MajorityVote advantage, agreement-driven quality, gene-specific agreement behavior, per-gene performance, the cross-modality landscape of tool behavior, and a winner-focused summary of the strongest methods for HLA-A, HLA-B, and HLA-C in each sequencing setting.</p>
      </div>
      <div class="figure-grid" style="margin-top:18px">
        {embed_fig('figure_01_accuracy_overview',
          'Figure 1. HLA typing accuracy across individual tools and ensemble methods, by sequencing modality.',
          'Three-panel bar chart with one panel each for WGS, WES, and RNA-seq. The x-axis lists HLA typing tools plus MajorityVote; the y-axis shows overall correct-call rate, defined as the proportion of evaluable sample × gene loci for which the exact two-field allele pair matched the truth set. Vertical error bars denote 95% Wilson confidence intervals, and the dashed red horizontal line marks the MajorityVote accuracy within each modality.',
          'Bars are ordered within each modality from higher to lower accuracy among the individual tools, with MajorityVote shown after a visual separator so that it can be compared directly against the strongest single-tool baseline. Each bar summarizes all callable HLA-A, HLA-B, and HLA-C loci available for that method in that modality. This figure provides the primary benchmark comparison by showing both effect size and uncertainty, and it allows the reader to see whether MVHLA’s equal-weight voting rule lies within, above, or below the confidence range of the best single-method competitors.',
          'fig1')}
        {embed_fig('figure_02_majority_vote_advantage',
          'Figure 2. Accuracy gain of MajorityVote relative to each individual HLA typing tool, per modality.',
          'Three horizontal bar charts, one per modality. The y-axis lists individual tools, and the x-axis shows the difference in overall correct-call rate between MajorityVote and that tool (MajorityVote minus tool accuracy, in absolute proportion units). Green bars indicate that MajorityVote performed better, whereas red bars indicate that the individual tool performed better.',
          'Each data point is a modality-specific pairwise comparison between the MVHLA majority-vote baseline and one constituent tool, computed on the same benchmark summary tables used in Figure 1. The vertical zero line marks parity. This figure isolates ensemble benefit more clearly than raw accuracy plots because positive values directly quantify the gain achieved by combining tool outputs, whereas negative values reveal settings in which a single strong caller still retains an advantage.',
          'fig2')}
        {embed_fig('figure_03_agreement_vs_accuracy',
          'Figure 3. MajorityVote accuracy by tool agreement level, shown separately for HLA-A, HLA-B, and HLA-C.',
          'Nine-panel figure arranged as a 3 × 3 grid. Rows correspond to HLA-A, HLA-B, and HLA-C; columns correspond to WGS, WES, and RNA-seq. In each panel, the x-axis is the number of tools agreeing on the final majority-vote allele pair, and the y-axis is the empirical accuracy of those majority-vote calls, measured as the fraction of callable loci with the correct two-field truth pair. Bar labels report the number of benchmark observations contributing to each agreement level.',
          'Each bar therefore represents a subset of sample × gene observations grouped by consensus support count, and the dashed line shows the linear trend between support level and empirical correctness. This figure is mechanistically important because it tests whether agreement count is a valid confidence proxy for MVHLA. By separating HLA-A, HLA-B, and HLA-C, it also shows whether the agreement-to-accuracy relationship is consistent across all three class I genes or whether some loci require stronger tool consensus before the same level of reliability is achieved.',
          'fig3')}
        {embed_fig('figure_04_per_gene_accuracy',
          'Figure 4. Per-locus HLA typing accuracy for HLA-A, HLA-B, and HLA-C across all tools and modalities.',
          'Three-panel grouped bar chart with one panel per modality. Within each panel, the x-axis lists the three class I loci HLA-A, HLA-B, and HLA-C, and colored bars denote individual tools plus MajorityVote. The y-axis shows per-gene overall correct-call rate, defined on the subset of benchmark rows corresponding to that locus.',
          'Each bar is therefore a gene-specific summary rather than an aggregate across loci. This decomposition is important because a method that appears strong overall may derive that result from only one or two comparatively easy genes. The figure allows the reader to determine whether MVHLA improves performance consistently across HLA-A, HLA-B, and HLA-C or whether the ensemble advantage is concentrated in a specific locus, particularly the more difficult HLA-C setting.',
          'fig4')}
        {embed_fig('figure_05_crossmodal_landscape',
          'Figure 5. Complete tool × modality accuracy matrix for all evaluated HLA typing tools.',
          'Heatmap summarizing all evaluated methods across all sequencing modalities. Rows correspond to tools plus MajorityVote; columns correspond to WGS, WES, and RNA-seq. Cell color encodes overall correct-call rate from 0 to 1, numeric labels report the exact value in each populated cell, and grey cells indicate that a tool was not evaluated for that modality.',
          'This figure compresses the benchmark into a single cross-modality view and is particularly useful for identifying design-scope mismatches. Tools that are specialized for DNA or RNA data can be recognized immediately by strong performance in one column and weak performance in another, while MajorityVote can be evaluated as a modality-agnostic consensus strategy positioned against the full matrix of constituent methods.',
          'fig5')}
        {embed_fig('figure_06_best_tool_by_gene_modality',
          'Figure 6. Best-performing HLA typing methods for HLA-A, HLA-B, and HLA-C under WGS, WES, and RNA-seq conditions.',
          'Three-panel grouped bar chart with one panel each for WGS, WES, and RNA-seq. Within each panel, the x-axis is grouped by HLA-A, HLA-B, and HLA-C, and the y-axis shows per-gene overall correct-call rate. Bars correspond to the individual tools available in that modality plus MajorityVote; the highest-performing method within each gene cluster is highlighted and annotated.',
          'This panel is designed to answer a practical benchmark question more directly than the full per-gene comparison: which method is strongest for each class I gene in each sequencing context. Non-winning bars remain visible for context but are visually muted, while the winning method in each cluster is outlined and labeled with its accuracy value. The figure therefore acts as a concise winner-focused summary of empirical benchmark performance without replacing the broader method-level detail provided by Figure 4.',
          'fig6')}
      </div>
    </section>

    <section id="methods-and-tool-coverage">
      <div class="section-header">
        <h2>Methods and Tool Coverage</h2>
        <p>This section keeps the benchmark methodology and tool design scope explicit so readers can interpret both the tables and the figures safely.</p>
      </div>
      <div class="grid grid-2">
        <div class="panel">
          <h3>Cohort and truth labels</h3>
          <p><strong>Source:</strong> 1000 Genomes Project phase 3 (n={N_REF} truth-labelled samples).</p>
          <p><strong>Truth set:</strong> Gourraud et al. (2014), evaluated at two-field resolution for HLA-A, HLA-B, and HLA-C.</p>
          <p><strong>Ancestry groups:</strong> CEU, FIN, GBR, TSI, and YRI.</p>
          <p><strong>Database version:</strong> IMGT/HLA v{imgt_ver}.</p>
        </div>
        <div class="panel">
          <h3>Ensemble strategies</h3>
          <p><strong>MajorityVote:</strong> the core MVHLA rule. Each callable tool contributes one vote, and the plurality pair is selected. Ties are broken alphabetically.</p>
          <p><strong>WeightedConsensus:</strong> tool weights combine per-modality reliability and calibrated confidence. Tools without usable confidence signals contribute reliability-only weights.</p>
        </div>
      </div>
      <div class="panel" style="margin-top:18px">
        <h3>Tools evaluated</h3>
        <div class="table-wrap">
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
        </div>
      </div>
    </section>

    <section id="venex-cohort-overview">
      <div class="section-header">
        <h2>VENEX Cohort Overview</h2>
        <p>The VENEX section covers an applied Finnish WGS cohort where no external truth labels are available, so the focus shifts from benchmark accuracy to operational evidence.</p>
      </div>
      <div class="grid two-col">
        <div class="panel">
          <h3>Cohort summary</h3>
          <p>The VENEX cohort contains {n_venex} WGS samples processed through the MVHLA pipeline at CSC Puhti. Two WGS-oriented callers were compared: <strong>HLA-HD</strong> and <strong>SpecHLA</strong>.</p>
          <p>Without external truth labels, the key readouts are per-locus callability, two-field inter-tool concordance, and allele frequency distributions at HLA-A, HLA-B, HLA-C, and HLA-DRB1.</p>
        </div>
        <div class="panel">
          <h3>Key finding</h3>
          <p>Full two-field concordance reaches {full_pcts.get('A',0):.0f}% for HLA-A, {full_pcts.get('B',0):.0f}% for HLA-B, {full_pcts.get('C',0):.0f}% for HLA-C, and {full_pcts.get('DRB1',0):.0f}% for HLA-DRB1.</p>
          <p>Higher class I concordance is consistent with the benchmark observation that both HLA-HD and SpecHLA behave reliably on WGS data for HLA-A, HLA-B, and HLA-C, while DRB1 remains more sensitive to database and alignment differences.</p>
        </div>
      </div>
      <div class="callout venex" style="margin-top:18px">
        <p><strong>Interpretation note.</strong> VENEX results should be read as consistency evidence, not ground-truth accuracy. Agreement between independent tools increases confidence, but it does not replace orthogonal validation.</p>
      </div>
    </section>

    <section id="venex-figure-gallery">
      <div class="section-header">
        <h2>VENEX Figure Gallery</h2>
        <p>These figures are generated during report creation from parsed HLA-HD and SpecHLA outputs in the provided VENEX directory.</p>
      </div>
      <div class="figure-grid">
        {fig_section(b64_v1,
          f'Figure V1. Per-locus callability rate for HLA-HD and SpecHLA in the VENEX WGS cohort (n={n_venex}).',
          f'Side-by-side bar chart across six classical HLA loci. The x-axis lists HLA-A, HLA-B, HLA-C, HLA-DRB1, HLA-DQB1, and HLA-DPB1; the y-axis shows the callability rate, defined as the proportion of the {n_venex} VENEX samples for which the tool returned at least one non-null allele assignment at that locus. Bars are plotted separately for HLA-HD and SpecHLA, with percentage labels above each bar.',
          'Because no orthogonal truth labels are available for VENEX, this figure measures technical completeness rather than correctness. Each bar summarizes the fraction of samples that yielded a usable call for a given tool–locus combination. High class I callability supports the operational robustness of both WGS-oriented callers in HLA-A, HLA-B, and HLA-C, whereas lower values in selected class II loci highlight the greater alignment and reference-database challenges in those regions.',
          'figV1')}
        {fig_section(b64_v2,
          f'Figure V2. Two-field inter-tool concordance between HLA-HD and SpecHLA at four classical HLA loci in the VENEX cohort (n={n_venex}).',
          'Stacked horizontal bar chart with one bar per locus for HLA-A, HLA-B, HLA-C, and HLA-DRB1. The x-axis shows the percentage of VENEX samples, and each bar is partitioned into four concordance categories: full concordance, partial concordance, discordance, and missing-call comparison. Segment labels report the percentage contribution of each category where large enough to display.',
          'Each stacked bar summarizes an unordered pairwise comparison between HLA-HD and SpecHLA at two-field resolution. “Full concordance” indicates that both alleles matched between tools, “partial concordance” indicates that one allele matched, “discordance” indicates that neither allele matched, and “missing” indicates that at least one tool failed to return a call. In the absence of external truth data, strong class I concordance provides supporting evidence that the observed WGS genotype calls are stable across independent callers.',
          'figV2')}
        {fig_section(b64_v3,
          f'Figure V3. HLA allele frequency distribution at HLA-A, HLA-B, HLA-C, and HLA-DRB1 in the VENEX WGS cohort (n={n_venex}).',
          f'Four-panel horizontal bar chart showing the top 12 observed two-field alleles for HLA-A, HLA-B, HLA-C, and HLA-DRB1. Within each panel, the y-axis lists allele names ordered from most to less frequent, and the x-axis gives the number of observed allele copies, with a maximum possible count of {2*n_venex} per locus because each sample contributes up to two alleles.',
          'Counts are derived from the available HLA-HD call when present and fall back to the SpecHLA call when HLA-HD is missing, after truncation to two-field resolution. These panels do not estimate population frequencies with external normalization, but they do provide a cohort-level sanity check: allele count patterns that broadly match expected Northern European HLA structure support the biological plausibility of the combined VENEX callset and help reveal obvious anomalies that would warrant further validation.',
          'figV3')}
      </div>
    </section>

    <section id="venex-sample-summary">
      <div class="section-header">
        <h2>VENEX Sample Summary</h2>
        <p>Per-sample HLA-HD and SpecHLA calls at HLA-A, HLA-B, and HLA-C are shown below with a concordance badge for each locus.</p>
      </div>
      <div class="panel">
        <p>Badges: <span class="conc-full">full</span> both alleles match, <span class="conc-partial">partial</span> one allele matches, <span class="conc-disc">discordant</span> no allele in common, <span class="conc-miss">missing</span> one or both tools have no call.</p>
        <details class="sample-table" style="margin-top:14px">
          <summary>&#9654; Show all {n_venex} VENEX samples</summary>
          <div>{venex_sample_table(venex)}</div>
        </details>
      </div>
    </section>

    <section id="reproducibility-metadata">
      <div class="section-header">
        <h2>Reproducibility Metadata</h2>
        <p>Key benchmark configuration fields are surfaced here so the report can be interpreted and regenerated without hunting through external files.</p>
      </div>
      <div class="panel">
        <div class="meta-list">
          <div class="meta-row"><strong>IMGT/HLA version</strong><span>{meta.get('imgt_hla_version', '—')}</span></div>
          <div class="meta-row"><strong>Success definition</strong><span>{meta.get('success_definition', '—')}</span></div>
          <div class="meta-row"><strong>Benchmark splits</strong><span>{meta.get('benchmark_splits', {})}</span></div>
          <div class="meta-row"><strong>Confidence formula</strong><span>{meta.get('confidence_weighting', {}).get('formula', '—')}</span></div>
          <div class="meta-row"><strong>Confidence fallback</strong><span>{meta.get('confidence_weighting', {}).get('fallback', '—')}</span></div>
          <div class="meta-row"><strong>Runtime weight override</strong><span>{meta.get('confidence_weighting', {}).get('runtime_weight_override', '—')}</span></div>
        </div>
      </div>
    </section>

    <section id="appendix">
      <div class="section-header">
        <h2>Appendix</h2>
        <p>Concrete file paths and workflow context used by this generated report.</p>
      </div>
      <div class="grid grid-2">
        <div class="panel">
          <h3>Generation provenance</h3>
          <div class="meta-list">
            <div class="meta-row"><strong>Script</strong><span>/scratch/project_2008084/pihla-publish/bin/generate_html_report_v7.py</span></div>
            <div class="meta-row"><strong>Tables directory</strong><span>{td}</span></div>
            <div class="meta-row"><strong>Figures directory</strong><span>{fd}</span></div>
            <div class="meta-row"><strong>VENEX directory</strong><span>{venex_dir}</span></div>
            <div class="meta-row"><strong>Generated on</strong><span>{today}</span></div>
          </div>
        </div>
        <div class="panel">
          <h3>Workflow notes</h3>
          <p>Pipeline: MVHLA v3.0 on Nextflow DSL2.</p>
          <p>Cluster context: CSC Puhti with SLURM scheduling.</p>
          <p>Benchmark scope: 1000 Genomes truth-labelled cohort for trimodal evaluation plus a separate VENEX applied cohort without external truth labels.</p>
        </div>
      </div>
    </section>

    <footer>
      <p>MVHLA benchmark report v7 &middot; Generated {today} &middot; IMGT/HLA v{imgt_ver} &middot; Two-field resolution</p>
      <p style="margin-top:4px">1000 Genomes trimodal benchmark plus VENEX Finnish WGS cohort</p>
      <p style="margin-top:4px">Workflow: MVHLA v3.0 &middot; Pipeline: Nextflow DSL2 &middot; Cluster: CSC Puhti (SLURM)</p>
    </footer>
  </main>
</div>
<script>{JS}</script>
</body>
</html>"""

    with open(output_path, 'w') as fh:
        fh.write(html)
    size_kb = len(html) // 1024
    print(f'Written: {output_path}  ({size_kb} KB)')


def main():
    p = argparse.ArgumentParser(description='Generate MVHLA HTML report v7')
    p.add_argument('--tables-dir', required=True)
    p.add_argument('--figures-dir', required=True)
    p.add_argument('--venex-dir',   required=True)
    p.add_argument('--output',      required=True)
    args = p.parse_args()
    render(args.tables_dir, args.figures_dir, args.venex_dir, args.output)


if __name__ == '__main__':
    main()
