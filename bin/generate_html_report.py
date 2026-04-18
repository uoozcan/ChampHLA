#!/usr/bin/env python3
"""Generate a standalone PIHLA HTML benchmark report with embedded PNG figures."""

import argparse
import base64
import csv
import html
import json
from collections import defaultdict
from datetime import date, datetime, timezone
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

DEFAULT_FIGURES_DIR = Path("/scratch/project_2008084/pihla-publish/analysis/figures_v2")
DEFAULT_TABLES_DIR = Path("/scratch/project_2008084/pihla-publish/analysis/full_cohort_benchmark/tables")
DEFAULT_OUTPUT = Path("/users/ozcanumu/scratch/project_2008084/pihla-publish/analysis/pihla_benchmark_report_v2.html")

FIGURE_SPECS = [
    {
        "filename": "figure_01_cohort_overview.png",
        "label": "Figure 1",
        "title": "Cohort overview and modality coverage",
        "caption": "Overview of the benchmark cohort, truth support, and the distribution of usable modality-specific outputs that feed downstream evaluation.",
        "explanation": (
            "This figure summarizes the structure of the benchmark before any method-level interpretation is attempted. "
            "It is expected to show how the cohort is assembled, how many samples carry the required truth annotations, "
            "and how many WES, WGS, and RNA-seq observations survive into the analysis set used for benchmarking. "
            "Read it as the boundary condition for the entire report: every comparison in later figures is constrained by this coverage profile."
        ),
        "importance": (
            "This figure matters because benchmark credibility depends on cohort composition and data completeness. "
            "If coverage is sparse or uneven across modalities, that directly affects which tools can be compared fairly and how much confidence to place in apparent performance differences."
        ),
    },
    {
        "filename": "figure_02_method_performance.png",
        "label": "Figure 2",
        "title": "Method performance across modalities",
        "caption": "Comparison of single-tool baselines, majority vote, and weighted consensus across the evaluated sequencing modalities.",
        "explanation": (
            "This figure compares end-point benchmark performance for the major calling strategies. "
            "It should be interpreted by reading across modalities and then within each modality across method families: "
            "single tools, simple aggregation, and confidence-weighted consensus. "
            "The key analytical question is whether ensemble behavior improves call quality, callable rate, or both relative to the baselines available for the same modality."
        ),
        "importance": (
            "This is one of the central benchmark figures because it answers the headline question of whether PIHLA's ensemble strategy provides practical value beyond the constituent tools. "
            "It frames the rest of the report by showing where consensus helps, where it only matches the best baseline, and where modality-specific weaknesses remain."
        ),
    },
    {
        "filename": "figure_03_per_gene_heatmap.png",
        "label": "Figure 3",
        "title": "Per-gene performance heatmap",
        "caption": "Gene-level heatmap of comparative performance across tools and modalities for the benchmark loci A, B, and C.",
        "explanation": (
            "This figure breaks aggregate accuracy into gene-specific behavior so that improvements are not hidden by averaging. "
            "The heatmap structure is useful for detecting whether gains are broad-based across loci or concentrated in one gene or one modality. "
            "Patterns such as uniformly strong HLA-A results but weaker HLA-C behavior, or modality-specific instability, become much easier to see here than in the overall comparison panel."
        ),
        "importance": (
            "This figure is important because HLA benchmarking is rarely uniform across loci. "
            "Per-gene structure helps distinguish a genuinely robust ensemble from one that appears strong only because easy loci dominate the average, and it also indicates where PIHLA still needs calibration or parser improvements."
        ),
    },
    {
        "filename": "figure_04_calibration_diagrams.png",
        "label": "Figure 4",
        "title": "Confidence calibration diagrams",
        "caption": "Calibration plots comparing predicted confidence against observed correctness for tools and modalities with usable confidence signals.",
        "explanation": (
            "This figure examines whether tool-native confidence behaves like a trustworthy probability proxy. "
            "A well-calibrated tool should place high-confidence predictions in bins that are actually more often correct, while poorly calibrated tools will show overconfident or underconfident behavior. "
            "The panels therefore serve as a bridge between raw benchmark outcomes and the weighting logic used by the ensemble."
        ),
        "importance": (
            "Calibration is essential to PIHLA because weighted consensus should reward reliable confidence and down-weight misleading scores. "
            "If calibration is poor, confidence-aware weighting can amplify errors rather than reduce them, so this figure is a direct check on whether the ensemble weighting strategy is justified."
        ),
    },
    {
        "filename": "figure_05_abstention_tradeoff.png",
        "label": "Figure 5",
        "title": "Abstention versus accuracy tradeoff",
        "caption": "Tradeoff curve showing how stricter support thresholds affect callable fraction and correctness among emitted calls.",
        "explanation": (
            "This figure evaluates the operational cost of being selective. "
            "As the minimum support threshold rises, the system should emit fewer calls but ideally improve the accuracy of the calls it keeps. "
            "The shape of the curve indicates whether abstention is buying meaningful quality improvement or merely suppressing coverage without much payoff."
        ),
        "importance": (
            "This figure matters because real deployments must balance completeness against reliability. "
            "It shows whether PIHLA can be tuned for conservative reporting behavior, and it gives a quantitative basis for choosing thresholds appropriate for research versus higher-confidence downstream use."
        ),
    },
    {
        "filename": "figure_06_discordance_taxonomy.png",
        "label": "Figure 6",
        "title": "Discordance taxonomy",
        "caption": "Categorization of disagreement patterns observed across tools, modalities, and consensus outputs in the benchmark cohort.",
        "explanation": (
            "This figure moves beyond accuracy into error structure. "
            "Rather than treating all wrong calls equally, it groups disagreements into interpretable categories such as broad tool conflict, partial support, or modality-specific inconsistency. "
            "Those categories help explain why specific samples fail and whether the ensemble is facing random noise, systematic caller disagreement, or missing evidence."
        ),
        "importance": (
            "This figure is important because debugging and improving PIHLA depends on understanding the shape of failure, not only the count of failures. "
            "A taxonomy of discordance points toward the next engineering steps: parser fixes, confidence guardrails, modality-aware weighting, or targeted abstention rules."
        ),
    },
    {
        "filename": "figure_07_weight_heatmap.png",
        "label": "Figure 7",
        "title": "Confidence-weight heatmap",
        "caption": "Heatmap of learned tool weights across modalities, reflecting benchmark performance and usable confidence behavior.",
        "explanation": (
            "This figure visualizes the ensemble's learned weighting scheme. "
            "Each cell reflects how PIHLA translates benchmark evidence into modality-specific trust for a tool, incorporating both overall correctness and the quality of the confidence signal when available. "
            "The heatmap should be read as an interpretable summary of which tools the ensemble leans on most in each sequencing context."
        ),
        "importance": (
            "This figure matters because it makes the ensemble transparent. "
            "Instead of treating consensus as a black box, it shows how benchmark evidence shapes the final decision policy and whether the weight distribution aligns with the observed strengths and weaknesses in the performance tables."
        ),
    },
    {
        "filename": "figure_08_ensemble_advantage.png",
        "label": "Figure 8",
        "title": "Ensemble advantage summary",
        "caption": "Focused comparison of how weighted consensus performs relative to majority vote and the strongest available single-tool baselines.",
        "explanation": (
            "This figure isolates the central comparative claim of the report: when and how much the ensemble helps. "
            "Rather than listing all methods equally, it emphasizes the delta between weighted consensus, simpler aggregation, and the best individual tool performance available for the same setting. "
            "That makes it easier to judge whether PIHLA provides additive value or mainly recapitulates the best baseline."
        ),
        "importance": (
            "This figure is important because it translates a large set of benchmark tables into a simple strategic conclusion. "
            "If the ensemble advantage is clear, it supports PIHLA's design choice; if it is narrow or modality-specific, the report can say so directly and point to where additional model or parser work is needed."
        ),
    },
]

MODALITY_ORDER = ["wes", "wgs", "rnaseq"]
TOOL_ORDER = ["ArcasHLA", "HLA-HD", "Kourami", "OptiType", "POLYSOLVER", "Seq2HLA", "SpecHLA", "T1K"]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-html", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--figures-dir", type=Path, default=DEFAULT_FIGURES_DIR)
    parser.add_argument("--tables-dir", type=Path, default=DEFAULT_TABLES_DIR)
    parser.add_argument("--title", default="PIHLA Benchmark HTML Report")
    parser.add_argument(
        "--subtitle",
        default=(
            "Scientific benchmark report using the newest figures_v2 graphics and full-cohort benchmark tables"
        ),
    )
    parser.add_argument("--authors", default="Generated in Claude Code for the PIHLA project")
    return parser.parse_args()


def read_tsv(path: Path) -> List[Dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def load_json(path):
    with path.open(encoding="utf-8") as handle:
        return json.load(handle)


def as_data_uri(path: Path) -> str:
    raw = path.read_bytes()
    return f"data:image/png;base64,{base64.b64encode(raw).decode('ascii')}"


def pct(value) -> str:
    if value in (None, ""):
        return "-"
    number = float(value)
    return f"{number * 100:.1f}%"


def integer(value) -> str:
    if value in (None, ""):
        return "-"
    return f"{int(value):,}"


def safe(text: object) -> str:
    return html.escape(str(text), quote=True)


def section_id(title: str) -> str:
    return (
        title.lower()
        .replace(" ", "-")
        .replace("/", "-")
        .replace("and", "and")
        .replace("--", "-")
    )


def sort_rows(rows: Iterable[Dict[str, str]], keys: Tuple[str, ...]) -> List[Dict[str, str]]:
    def key_fn(row):
        parts = []
        for key in keys:
            value = row.get(key, "")
            if key == "modality":
                parts.append(MODALITY_ORDER.index(value) if value in MODALITY_ORDER else 999)
            elif key == "tool":
                parts.append(TOOL_ORDER.index(value) if value in TOOL_ORDER else 999)
            else:
                parts.append(value)
        return tuple(parts)

    return sorted(rows, key=key_fn)


def compute_performance_summary(method_rows):
    grouped = defaultdict(list)
    for row in method_rows:
        grouped[row["modality"]].append(row)

    summary_cards = []
    finding_lines = []
    for modality in MODALITY_ORDER:
        rows = grouped.get(modality, [])
        if not rows:
            continue
        single_tools = [r for r in rows if r["method_type"] == "single_tool"]
        single_best = max(single_tools, key=lambda r: float(r["overall_correct_call_rate"])) if single_tools else None
        ensemble = next((r for r in rows if r["method"] == "WeightedConsensus"), None)
        majority = next((r for r in rows if r["method"] == "MajorityVote"), None)
        summary_cards.append(
            {
                "modality": modality.upper(),
                "best_single": single_best["method"] if single_best else "-",
                "best_single_rate": pct(single_best["overall_correct_call_rate"]) if single_best else "-",
                "ensemble_rate": pct(ensemble["overall_correct_call_rate"]) if ensemble else "-",
                "majority_rate": pct(majority["overall_correct_call_rate"]) if majority else "-",
            }
        )
        if single_best and ensemble and majority:
            finding_lines.append(
                f"{modality.upper()}: best single-tool overall correct call rate is {single_best['method']} at {pct(single_best['overall_correct_call_rate'])}, "
                f"WeightedConsensus reaches {pct(ensemble['overall_correct_call_rate'])}, and MajorityVote reaches {pct(majority['overall_correct_call_rate'])}."
            )
    return summary_cards, finding_lines


def build_metric_cards(metadata, figure_count):
    splits = metadata.get("split_membership_summary", {})
    modalities = metadata.get("modalities", [])
    supported = metadata.get("supported_loci", [])
    return [
        {
            "label": "Final Tri-Modal Cohort",
            "value": integer(metadata.get("final_tri_modal_cohort_size")),
            "note": "Full benchmark cohort size across the selected mixed-source report inputs.",
        },
        {
            "label": "Benchmark Splits",
            "value": f"{integer(splits.get('training'))} / {integer(splits.get('validation'))} / {integer(splits.get('holdout'))}",
            "note": "Training, validation, and holdout sample counts used for weighting, tuning, and final reporting.",
        },
        {
            "label": "Modalities",
            "value": ", ".join(m.upper() for m in modalities) or "-",
            "note": "Sequencing modalities represented in the table-backed benchmark context.",
        },
        {
            "label": "Supported Loci",
            "value": ", ".join(supported) or "-",
            "note": f"Configured benchmark genes evaluated at primary resolution {metadata.get('resolution', '-')}.",
        },
        {
            "label": "Weight Version",
            "value": metadata.get("confidence_weighting", {}).get("weight_version", "-"),
            "note": "Current confidence-aware weighting version recorded in benchmark metadata.",
        },
        {
            "label": "Figures Embedded",
            "value": integer(figure_count),
            "note": "Newest figures_v2 PNG assets embedded directly into this offline HTML report.",
        },
    ]


def render_table(headers, rows, caption=None):
    head_html = "".join(f"<th>{safe(h)}</th>" for h in headers)
    body_rows = []
    for row in rows:
        body_rows.append("<tr>" + "".join(f"<td>{cell}</td>" for cell in row) + "</tr>")
    caption_html = f"<caption>{safe(caption)}</caption>" if caption else ""
    return (
        "<div class='table-wrap'><table>"
        f"{caption_html}<thead><tr>{head_html}</tr></thead>"
        f"<tbody>{''.join(body_rows)}</tbody></table></div>"
    )


def build_tool_status_rows(rows):
    table_rows = []
    for row in sort_rows(rows, ("modality", "tool", "status")):
        status = row["status"]
        badge = f"<span class='status {safe(status)}'>{safe(status.replace('_', ' '))}</span>"
        table_rows.append([
            safe(row["modality"].upper()),
            safe(row["tool"]),
            badge,
            safe(row["n_samples"]),
        ])
    return table_rows


def build_performance_rows(rows):
    sorted_rows = sorted(
        rows,
        key=lambda r: (
            MODALITY_ORDER.index(r["modality"]) if r["modality"] in MODALITY_ORDER else 999,
            {"ensemble": 0, "baseline": 1, "single_tool": 2}.get(r["method_type"], 9),
            -(float(r["overall_correct_call_rate"] or 0.0)),
            r["method"],
        ),
    )
    out = []
    for row in sorted_rows:
        out.append([
            safe(row["modality"].upper()),
            safe(row["method"]),
            safe(row["method_type"].replace("_", " ")),
            safe(row["sample_count"]),
            safe(row["gene_rows"]),
            pct(row["callable_rate"]),
            pct(row["accuracy_among_callable"]),
            pct(row["overall_correct_call_rate"]),
        ])
    return out


def figure_card(spec: Dict[str, str], figures_dir: Path) -> str:
    path = figures_dir / spec["filename"]
    if not path.exists():
        raise FileNotFoundError(f"Missing figure: {path}")
    image = as_data_uri(path)
    return f"""
    <article class="figure-card panel">
      <div class="figure-header">
        <div>
          <p class="eyebrow">{safe(spec['label'])}</p>
          <h3>{safe(spec['title'])}</h3>
        </div>
        <span class="figure-file">{safe(spec['filename'])}</span>
      </div>
      <div class="figure-box">
        <img src="{image}" alt="{safe(spec['title'])}" />
      </div>
      <div class="figure-copy">
        <p><strong>Caption.</strong> {safe(spec['caption'])}</p>
        <p><strong>What this figure shows.</strong> {safe(spec['explanation'])}</p>
        <p><strong>Meaning and importance.</strong> {safe(spec['importance'])}</p>
      </div>
    </article>
    """


def build_html(args: argparse.Namespace) -> str:
    figures_dir = args.figures_dir
    tables_dir = args.tables_dir

    metadata = load_json(tables_dir / "benchmark_metadata.json")
    method_rows = read_tsv(tables_dir / "method_comparison.tsv")
    tool_status_rows = read_tsv(tables_dir / "tool_availability_by_modality.tsv")
    cohort_rows = read_tsv(tables_dir / "cohort_overview.tsv")
    summary_rows = read_tsv(tables_dir / "summary_full_cohort.tsv")

    generated_at = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC")
    metric_cards = build_metric_cards(metadata, len(FIGURE_SPECS))
    performance_cards, finding_lines = compute_performance_summary(method_rows)

    sections = [
        "Executive Summary",
        "Report Scope and Provenance",
        "Cohort and Benchmark Context",
        "Tool Status by Modality",
        "Benchmark Performance Summary",
        "Figure Gallery",
        "Interpretation and Key Findings",
        "Limitations and Caveats",
        "Reproducibility Metadata",
        "Appendix",
    ]

    nav_html = "".join(
        f"<a href='#{section_id(section)}'>{safe(section)}</a>" for section in sections
    )

    metric_html = "".join(
        f"<div class='metric'><div class='label'>{safe(card['label'])}</div><div class='value'>{safe(card['value'])}</div><div class='note'>{safe(card['note'])}</div></div>"
        for card in metric_cards
    )

    performance_card_html = "".join(
        f"<div class='metric'><div class='label'>{safe(card['modality'])}</div><div class='value'>{safe(card['ensemble_rate'])}</div><div class='note'>WeightedConsensus overall correct call rate. Best single tool: {safe(card['best_single'])} at {safe(card['best_single_rate'])}. MajorityVote: {safe(card['majority_rate'])}.</div></div>"
        for card in performance_cards
    )

    tool_status_table = render_table(
        ["Modality", "Tool", "Status", "N samples"],
        build_tool_status_rows(tool_status_rows),
        "Tool availability status derived from full-cohort benchmark tables.",
    )

    performance_table = render_table(
        [
            "Modality",
            "Method",
            "Type",
            "Samples",
            "Gene rows",
            "Callable rate",
            "Accuracy among callable",
            "Overall correct call rate",
        ],
        build_performance_rows(method_rows),
        "Method comparison summary from full-cohort benchmark tables.",
    )

    cohort_table = render_table(
        ["Modality", "Tool", "Truth samples", "Completed tool samples", "Modality samples", "Shared gene count"],
        [
            [
                safe(row["modality"].upper()),
                safe(row["tool"]),
                safe(row["truth_samples"]),
                safe(row["completed_tool_samples"]),
                safe(row["modality_samples"]),
                safe(row["shared_gene_count"]),
            ]
            for row in sort_rows(cohort_rows, ("modality", "tool"))
        ],
        "Coverage and completion overview used to contextualize method-level comparisons.",
    )

    summary_table = render_table(
        ["Tool", "Modality", "Samples", "Gene rows", "Callable rate", "Accuracy among callable", "Overall correct call rate"],
        [
            [
                safe(row["tool"]),
                safe(row["modality"].upper()),
                safe(row["sample_count"]),
                safe(row["gene_rows"]),
                pct(row["callable_rate"]),
                pct(row["accuracy_among_callable"]),
                pct(row["overall_correct_call_rate"]),
            ]
            for row in sort_rows(summary_rows, ("modality", "tool"))
        ],
        "Summary full-cohort performance table used to ground headline comparisons.",
    )

    findings_html = "".join(f"<li>{safe(line)}</li>" for line in finding_lines)
    figure_html = "".join(figure_card(spec, figures_dir) for spec in FIGURE_SPECS)
    modalities = ", ".join(m.upper() for m in metadata.get("modalities", []))
    supported_loci = ", ".join(metadata.get("supported_loci", []))
    split_summary = metadata.get("split_membership_summary", {})
    confidence = metadata.get("confidence_weighting", {})

    appendix_rows = [
        [safe("Figures directory"), safe(str(figures_dir))],
        [safe("Tables directory"), safe(str(tables_dir))],
        [safe("Benchmark scope"), safe(metadata.get("tool_coverage_policy", "-"))],
        [safe("Weight version"), safe(confidence.get("weight_version", "-"))],
        [safe("Generated at"), safe(generated_at)],
        [safe("Generated on"), safe(str(date.today()))],
    ]
    appendix_table = render_table(["Field", "Value"], appendix_rows, "Report generation provenance.")

    return f"""<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8" />
  <meta name="viewport" content="width=device-width, initial-scale=1.0" />
  <title>{safe(args.title)}</title>
  <style>
    :root {{
      --bg: #f3f6fb;
      --surface: #ffffff;
      --surface-soft: #f8fbff;
      --ink: #102033;
      --muted: #5b6472;
      --line: #d7e1ee;
      --navy: #0f2744;
      --blue: #275dad;
      --blue-deep: #163c74;
      --teal: #0f766e;
      --amber: #b45309;
      --rose: #b91c1c;
      --shadow: 0 20px 44px rgba(15, 39, 68, 0.10);
      --radius: 22px;
      --sidebar-width: 300px;
    }}

    * {{ box-sizing: border-box; }}
    html {{ scroll-behavior: smooth; }}
    body {{
      margin: 0;
      font-family: Arial, Helvetica, sans-serif;
      color: var(--ink);
      background: linear-gradient(180deg, #eef3f9 0%, #f8fbfd 100%);
    }}

    .layout {{
      display: grid;
      grid-template-columns: var(--sidebar-width) minmax(0, 1fr);
      min-height: 100vh;
    }}

    aside {{
      position: sticky;
      top: 0;
      align-self: start;
      height: 100vh;
      padding: 28px 22px;
      background: linear-gradient(180deg, #102844 0%, #174576 55%, #1f5ea8 100%);
      color: #fff;
      border-right: 1px solid rgba(255,255,255,0.10);
    }}

    .brand {{
      margin-bottom: 22px;
      padding-bottom: 18px;
      border-bottom: 1px solid rgba(255,255,255,0.14);
    }}

    .brand .kicker {{
      margin: 0 0 10px;
      color: #c8dbf2;
      font-size: 12px;
      letter-spacing: 0.08em;
      text-transform: uppercase;
      font-weight: 700;
    }}

    .brand h1 {{
      margin: 0 0 10px;
      font-size: 28px;
      line-height: 1.08;
    }}

    .brand p {{
      margin: 0;
      color: #d7e6f7;
      font-size: 14px;
      line-height: 1.5;
    }}

    .sidebar-meta {{
      display: grid;
      gap: 12px;
      margin: 18px 0 24px;
    }}

    .sidebar-chip {{
      padding: 10px 12px;
      border: 1px solid rgba(255,255,255,0.14);
      border-radius: 16px;
      background: rgba(255,255,255,0.07);
      font-size: 13px;
      color: #e7f0fa;
    }}

    .sidebar-chip strong {{ display: block; margin-bottom: 4px; color: #fff; font-size: 12px; text-transform: uppercase; letter-spacing: 0.06em; }}

    .nav {{ display: grid; gap: 8px; }}

    .nav a {{
      color: #d9e8f7;
      text-decoration: none;
      padding: 10px 12px;
      border-radius: 14px;
      border: 1px solid transparent;
      font-size: 14px;
      line-height: 1.35;
    }}

    .nav a:hover,
    .nav a.active {{
      background: rgba(255,255,255,0.10);
      border-color: rgba(255,255,255,0.16);
      color: #fff;
    }}

    main {{ padding: 30px 32px 48px; }}

    .hero {{
      background: linear-gradient(135deg, #102844 0%, #174576 55%, #275dad 100%);
      color: #fff;
      border-radius: 28px;
      box-shadow: var(--shadow);
      padding: 32px 34px;
      margin-bottom: 24px;
    }}

    .hero .eyebrow,
    .eyebrow {{
      margin: 0 0 12px;
      font-size: 12px;
      letter-spacing: 0.08em;
      text-transform: uppercase;
      font-weight: 700;
      color: #c9dcf4;
    }}

    .hero h2 {{ margin: 0 0 10px; font-size: 38px; line-height: 1.08; }}
    .hero p {{ margin: 0; max-width: 1040px; font-size: 17px; line-height: 1.6; color: #e3eef9; }}

    section {{ margin-top: 22px; }}
    .section-header {{ margin-bottom: 16px; }}
    .section-header h2 {{ margin: 0 0 8px; font-size: 28px; color: var(--navy); }}
    .section-header p {{ margin: 0; color: var(--muted); font-size: 15px; line-height: 1.6; max-width: 960px; }}

    .grid {{ display: grid; gap: 18px; }}
    .metrics-grid {{ grid-template-columns: repeat(3, minmax(0, 1fr)); }}
    .two-col {{ grid-template-columns: 1.15fr 0.85fr; align-items: start; }}

    .panel, .metric {{
      background: linear-gradient(180deg, #fbfdff 0%, #ffffff 100%);
      border: 1px solid var(--line);
      border-radius: var(--radius);
      box-shadow: 0 8px 24px rgba(15, 39, 68, 0.05);
    }}

    .panel {{ padding: 22px; }}
    .panel h3 {{ margin: 0 0 12px; font-size: 20px; color: var(--navy); }}
    .panel p, .panel li {{ margin: 0; font-size: 15px; line-height: 1.65; color: var(--ink); }}
    .panel ul {{ margin: 0; padding-left: 18px; display: grid; gap: 10px; }}

    .metric {{ padding: 18px; }}
    .metric .label {{ font-size: 12px; text-transform: uppercase; letter-spacing: 0.06em; color: var(--muted); font-weight: 700; margin-bottom: 8px; }}
    .metric .value {{ font-size: 30px; font-weight: 700; color: var(--navy); line-height: 1.1; margin-bottom: 8px; }}
    .metric .note {{ font-size: 14px; line-height: 1.5; color: var(--muted); }}

    .callout {{
      padding: 18px 20px;
      border-left: 5px solid var(--blue);
      background: #f5f9ff;
      border-radius: 16px;
      color: var(--ink);
    }}

    .callout.warn {{ border-left-color: var(--amber); background: #fffaf2; }}

    .table-wrap {{ overflow-x: auto; border: 1px solid var(--line); border-radius: 18px; background: #fff; }}
    table {{ width: 100%; border-collapse: collapse; min-width: 720px; }}
    caption {{ text-align: left; padding: 16px 18px 0; font-size: 14px; color: var(--muted); }}
    thead th {{ position: sticky; top: 0; background: #eef5fd; color: var(--navy); font-size: 12px; text-transform: uppercase; letter-spacing: 0.05em; }}
    th, td {{ padding: 12px 14px; border-bottom: 1px solid #e7eef6; text-align: left; font-size: 14px; vertical-align: top; }}
    tbody tr:nth-child(even) {{ background: #fbfdff; }}

    .status {{ display: inline-block; padding: 5px 10px; border-radius: 999px; font-size: 12px; font-weight: 700; text-transform: uppercase; letter-spacing: 0.04em; }}
    .status.available {{ color: #065f46; background: #dff7ef; border: 1px solid #b6e7d7; }}
    .status.not_available {{ color: #9a3412; background: #fff1e8; border: 1px solid #f7d4c1; }}

    .figure-grid {{ display: grid; gap: 20px; }}
    .figure-card {{ padding: 22px; }}
    .figure-header {{ display: flex; justify-content: space-between; gap: 16px; align-items: flex-start; margin-bottom: 16px; }}
    .figure-header h3 {{ margin: 0; font-size: 22px; }}
    .figure-file {{ font-size: 12px; color: var(--muted); background: #f2f6fb; border: 1px solid var(--line); padding: 8px 10px; border-radius: 999px; }}
    .figure-box {{ border: 1px solid var(--line); border-radius: 20px; background: #f9fbfe; padding: 14px; display: flex; justify-content: center; align-items: center; min-height: 340px; margin-bottom: 16px; }}
    .figure-box img {{ display: block; max-width: 100%; width: 100%; height: auto; border-radius: 14px; }}
    .figure-copy {{ display: grid; gap: 12px; }}

    .meta-list {{ display: grid; gap: 10px; }}
    .meta-row {{ display: grid; grid-template-columns: 220px 1fr; gap: 14px; padding: 12px 0; border-bottom: 1px solid #e8eef5; }}
    .meta-row:last-child {{ border-bottom: 0; }}
    .meta-row strong {{ color: var(--navy); }}

    @media (max-width: 1180px) {{
      .layout {{ grid-template-columns: 1fr; }}
      aside {{ position: relative; height: auto; }}
      main {{ padding-top: 20px; }}
      .metrics-grid, .two-col {{ grid-template-columns: 1fr; }}
    }}

    @media (max-width: 720px) {{
      main {{ padding: 16px; }}
      aside {{ padding: 20px 16px; }}
      .hero {{ padding: 24px 20px; }}
      .hero h2 {{ font-size: 30px; }}
      .section-header h2 {{ font-size: 24px; }}
      .figure-header {{ flex-direction: column; }}
      .meta-row {{ grid-template-columns: 1fr; gap: 6px; }}
    }}
  </style>
</head>
<body>
  <div class="layout">
    <aside>
      <div class="brand">
        <p class="kicker">PIHLA report</p>
        <h1>{safe(args.title)}</h1>
        <p>{safe(args.subtitle)}</p>
      </div>
      <div class="sidebar-meta">
        <div class="sidebar-chip"><strong>Authors</strong>{safe(args.authors)}</div>
        <div class="sidebar-chip"><strong>Figures</strong>{safe(str(figures_dir))}</div>
        <div class="sidebar-chip"><strong>Tables</strong>{safe(str(tables_dir))}</div>
      </div>
      <nav class="nav" id="sidebar-nav">{nav_html}</nav>
    </aside>
    <main>
      <section class="hero">
        <p class="eyebrow">Offline standalone HTML benchmark report</p>
        <h2>{safe(args.title)}</h2>
        <p>{safe(args.subtitle)}. This report embeds all PNG figures directly into the document, uses the newest `figures_v2` graphics, and grounds its tables and context panels in the full-cohort benchmark outputs.</p>
      </section>

      <section id="{section_id('Executive Summary')}">
        <div class="section-header">
          <h2>Executive Summary</h2>
          <p>The report combines the newest PIHLA visualization set with full-cohort benchmark tables to produce an offline, scientific HTML report focused on modality coverage, tool availability, benchmark performance, and ensemble interpretation.</p>
        </div>
        <div class="grid metrics-grid">{metric_html}</div>
      </section>

      <section id="{section_id('Report Scope and Provenance')}">
        <div class="section-header">
          <h2>Report Scope and Provenance</h2>
          <p>This report is intentionally transparent about its mixed-source construction so readers can distinguish the newest figure assets from the benchmark tables that provide context and modality-level status information.</p>
        </div>
        <div class="grid two-col">
          <div class="panel">
            <h3>Source boundaries</h3>
            <div class="meta-list">
              <div class="meta-row"><strong>Figures</strong><span>{safe(str(figures_dir))}</span></div>
              <div class="meta-row"><strong>Tables</strong><span>{safe(str(tables_dir))}</span></div>
              <div class="meta-row"><strong>Generated at</strong><span>{safe(generated_at)}</span></div>
              <div class="meta-row"><strong>Mixed-source policy</strong><span>Figures come from `analysis/figures_v2`, while tabular summaries come from `analysis/full_cohort_benchmark/tables`.</span></div>
            </div>
          </div>
          <div class="panel">
            <h3>Why this pairing was chosen</h3>
            <ul>
              <li>The figure bundle under `figures_v2` is the newest visual asset set available in `/scratch`.</li>
              <li>The full-cohort tables provide the tri-modality coverage needed for a complete <em>Tool Status by Modality</em> section.</li>
              <li>The report keeps provenance explicit so visual interpretation and table-backed metrics are never conflated.</li>
            </ul>
          </div>
        </div>
      </section>

      <section id="{section_id('Cohort and Benchmark Context')}">
        <div class="section-header">
          <h2>Cohort and Benchmark Context</h2>
          <p>The benchmark metadata describes a full tri-modal cohort with explicit training, validation, and holdout splits, using loci {safe(supported_loci)} across modalities {safe(modalities)}.</p>
        </div>
        <div class="grid two-col">
          <div class="panel">
            <h3>Benchmark context</h3>
            <div class="meta-list">
              <div class="meta-row"><strong>Final tri-modal cohort size</strong><span>{integer(metadata.get('final_tri_modal_cohort_size'))}</span></div>
              <div class="meta-row"><strong>Split membership</strong><span>Training {integer(split_summary.get('training'))}, validation {integer(split_summary.get('validation'))}, holdout {integer(split_summary.get('holdout'))}</span></div>
              <div class="meta-row"><strong>Supported loci</strong><span>{safe(supported_loci)}</span></div>
              <div class="meta-row"><strong>Weight version</strong><span>{safe(confidence.get('weight_version', '-'))}</span></div>
              <div class="meta-row"><strong>Confidence rows</strong><span>{integer(confidence.get('harmonized_rows_with_confidence'))}</span></div>
              <div class="meta-row"><strong>Coverage policy</strong><span>{safe(metadata.get('tool_coverage_policy', '-'))}</span></div>
            </div>
          </div>
          <div class="panel">
            <h3>Coverage interpretation</h3>
            <ul>
              <li>The full-cohort benchmark preserves all three modalities in the metadata layer, which is critical for interpreting tool availability and fair comparison boundaries.</li>
              <li>Coverage is phase-gated rather than uniform, so the report distinguishes between theoretical modality membership and completed tool outputs.</li>
              <li>The cohort overview table below is the operational bridge between benchmark scope and the downstream performance panels.</li>
            </ul>
          </div>
        </div>
        <div class="panel">{cohort_table}</div>
      </section>

      <section id="{section_id('Tool Status by Modality')}">
        <div class="section-header">
          <h2>Tool Status by Modality</h2>
          <p>This table is built directly from `tool_availability_by_modality.tsv` and makes missingness visible instead of hiding it behind aggregate benchmark plots.</p>
        </div>
        <div class="panel">{tool_status_table}</div>
      </section>

      <section id="{section_id('Benchmark Performance Summary')}">
        <div class="section-header">
          <h2>Benchmark Performance Summary</h2>
          <p>Performance summaries below come from the full-cohort comparison tables and provide the quantitative context needed to interpret the newest figure set.</p>
        </div>
        <div class="grid metrics-grid">{performance_card_html}</div>
        <div class="panel" style="margin-top: 18px;">{performance_table}</div>
        <div class="panel" style="margin-top: 18px;">{summary_table}</div>
      </section>

      <section id="{section_id('Figure Gallery')}">
        <div class="section-header">
          <h2>Figure Gallery</h2>
          <p>Every figure is embedded directly into this HTML file and accompanied by a caption, a detailed explanation of what the figure is showing, and a statement of why the result matters for PIHLA.</p>
        </div>
        <div class="figure-grid">{figure_html}</div>
      </section>

      <section id="{section_id('Interpretation and Key Findings')}">
        <div class="section-header">
          <h2>Interpretation and Key Findings</h2>
          <p>The mixed-source report should be read as a scientific status update rather than a single monolithic benchmark export. The figure set is current, while the tables provide the most complete cross-modality grounding available.</p>
        </div>
        <div class="grid two-col">
          <div class="panel">
            <h3>Key findings</h3>
            <ul>{findings_html}</ul>
          </div>
          <div class="panel">
            <h3>How to read the report</h3>
            <ul>
              <li>Use the cohort and tool-status sections to understand what was actually available for evaluation in each modality.</li>
              <li>Use the performance tables and Figure 2 / Figure 8 together to judge whether weighted consensus adds value beyond majority vote and strong single-tool baselines.</li>
              <li>Use the calibration, abstention, and weight figures as explanation layers for why the ensemble behaves as it does.</li>
            </ul>
          </div>
        </div>
      </section>

      <section id="{section_id('Limitations and Caveats')}">
        <div class="section-header">
          <h2>Limitations and Caveats</h2>
          <p>The report deliberately prioritizes transparency over cosmetic simplification, especially where modality coverage and tool completion differ.</p>
        </div>
        <div class="grid two-col">
          <div class="callout warn">
            <strong>Mixed-source caveat.</strong> The visual narrative comes from `figures_v2`, while tabular evidence and modality availability come from the full-cohort benchmark directory. Readers should treat the provenance section as part of the interpretation, not as a footnote.
          </div>
          <div class="callout warn">
            <strong>Coverage caveat.</strong> Completed tool outputs are not uniform across modalities, so apparent differences in benchmark performance can reflect both model behavior and the subset of samples actually available for each tool.
          </div>
        </div>
      </section>

      <section id="{section_id('Reproducibility Metadata')}">
        <div class="section-header">
          <h2>Reproducibility Metadata</h2>
          <p>Key benchmark metadata fields are surfaced here so the report can be interpreted and regenerated without external context.</p>
        </div>
        <div class="panel">
          <div class="meta-list">
            <div class="meta-row"><strong>IMGT/HLA version</strong><span>{safe(metadata.get('imgt_hla_version', '-'))}</span></div>
            <div class="meta-row"><strong>Success definition</strong><span>{safe(metadata.get('success_definition', '-'))}</span></div>
            <div class="meta-row"><strong>Benchmark splits</strong><span>{safe(json.dumps(metadata.get('benchmark_splits', {}), sort_keys=True))}</span></div>
            <div class="meta-row"><strong>Confidence formula</strong><span>{safe(confidence.get('formula', '-'))}</span></div>
            <div class="meta-row"><strong>Confidence fallback</strong><span>{safe(confidence.get('fallback', '-'))}</span></div>
            <div class="meta-row"><strong>Runtime weight override</strong><span>{safe(confidence.get('runtime_weight_override', '-'))}</span></div>
          </div>
        </div>
      </section>

      <section id="{section_id('Appendix')}">
        <div class="section-header">
          <h2>Appendix</h2>
          <p>Appendix information captures the concrete file sources used by this generated report.</p>
        </div>
        <div class="panel">{appendix_table}</div>
      </section>
    </main>
  </div>
  <script>
    const navLinks = Array.from(document.querySelectorAll('#sidebar-nav a'));
    const sections = navLinks
      .map(link => document.querySelector(link.getAttribute('href')))
      .filter(Boolean);
    const byId = new Map(navLinks.map(link => [link.getAttribute('href').slice(1), link]));
    const observer = new IntersectionObserver((entries) => {{
      entries.forEach(entry => {{
        if (!entry.isIntersecting) return;
        navLinks.forEach(link => link.classList.remove('active'));
        const link = byId.get(entry.target.id);
        if (link) link.classList.add('active');
      }});
    }}, {{ rootMargin: '-30% 0px -55% 0px', threshold: 0.01 }});
    sections.forEach(section => observer.observe(section));
    if (navLinks.length) navLinks[0].classList.add('active');
  </script>
</body>
</html>
"""


def main() -> None:
    args = parse_args()
    args.output_html.parent.mkdir(parents=True, exist_ok=True)
    content = build_html(args)
    args.output_html.write_text(content, encoding="utf-8")
    print(args.output_html)


if __name__ == "__main__":
    main()
