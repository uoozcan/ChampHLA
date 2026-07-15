#!/usr/bin/env python3.11
"""Build the MVHLA presentation preparation PDF using matplotlib PdfPages."""

import argparse
import json
import math
import re
import textwrap
from datetime import datetime
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.image as mpimg
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages


PAGE_FOOTER = "MVHLA presentation prep"
FIGURE_ORDER = [
    "figure_1_workflow_architecture",
    "figure_2_accuracy_comparison",
    "figure_3_per_gene_gains",
    "figure_4_confidence_calibration",
    "figure_5_abstention_tradeoff",
    "figure_6_discordance_taxonomy",
    "figure_7_confidence_weights",
    "figure_09_bimodal_per_gene",
    "figure_10_trimodal_comparison",
    "figure_fimm_loh_wes",
]
OPTIONAL_DOC_KEYS = [
    "presentation_prompt",
    "manuscript",
    "evidence_map",
    "captions",
    "figures_readme",
]


def parse_args():
    repo_root = Path(__file__).resolve().parent.parent
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--content",
        type=Path,
        default=repo_root / "docs" / "PRESENTATION_PREP_CONTENT.md",
    )
    parser.add_argument(
        "--qa",
        type=Path,
        default=repo_root / "docs" / "PRESENTATION_QA_BANK.md",
    )
    parser.add_argument(
        "--figures-dir",
        type=Path,
        default=Path("/scratch/project_2008084/pihla-publish/analysis/figures_final"),
    )
    parser.add_argument(
        "--hub-outputs-dir",
        type=Path,
        default=repo_root / "figure_hub" / "outputs",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("/scratch/project_2008084/pihla-publish/analysis/presentation_prep/mvhla_presentation_prep.pdf"),
    )
    parser.add_argument(
        "--manifest",
        type=Path,
        default=Path("/scratch/project_2008084/pihla-publish/analysis/presentation_prep/mvhla_presentation_prep_manifest.json"),
    )
    parser.add_argument(
        "--fimm-figure",
        type=Path,
        default=Path("/scratch/project_2008084/pihla_local/analysis/figures_final/figure_fimm_loh_wes.png"),
    )
    parser.add_argument(
        "--fimm-survival-report",
        type=Path,
        default=Path("/scratch/project_2008084/pihla_local/fimm_survival_hed_report.html"),
    )
    parser.add_argument(
        "--page-size",
        choices=("letter", "a4"),
        default="letter",
    )
    return parser.parse_args()


def clean_text(value):
    if value is None:
        return ""
    value = str(value).strip()
    if value.startswith("`") and value.endswith("`"):
        return value[1:-1]
    return value


def load_markdown_sections(path):
    sections = []
    current = None
    subsection = None
    for raw_line in path.read_text(encoding="utf-8").splitlines():
        line = raw_line.rstrip()
        if line.startswith("## "):
            current = {"title": line[3:].strip(), "lines": [], "subsections": []}
            sections.append(current)
            subsection = None
            continue
        if line.startswith("### "):
            subsection = {"title": line[4:].strip(), "lines": []}
            if current is None:
                current = {"title": "Untitled", "lines": [], "subsections": []}
                sections.append(current)
            current["subsections"].append(subsection)
            continue
        target = subsection if subsection is not None else current
        if target is not None:
            target["lines"].append(line)
    return sections


def parse_bullet_map(lines):
    data = {}
    current_key = None
    for line in lines:
        stripped = line.strip()
        if not stripped:
            continue
        if stripped.startswith("- ") and ":" in stripped:
            key, value = stripped[2:].split(":", 1)
            key = key.strip()
            value = clean_text(value.strip())
            if value:
                data[key] = value
                current_key = None
            else:
                data[key] = []
                current_key = key
        elif stripped.startswith("- "):
            item = clean_text(stripped[2:].strip())
            if current_key is not None:
                data.setdefault(current_key, []).append(item)
            else:
                data.setdefault("items", []).append(item)
    return data


def section_map(sections):
    return {section["title"]: section for section in sections}


def subsection_map(section):
    return {sub["title"]: sub for sub in section.get("subsections", [])}


def load_qa_bank(path):
    blocks = load_markdown_sections(path)
    parsed = []
    for block in blocks:
        if not block["title"].startswith("QA-"):
            continue
        item = parse_bullet_map(block["lines"])
        item["question_id"] = block["title"]
        parsed.append(item)
    return parsed


def wrap_text_to_axes(ax, text, x, y, width=88, line_height=0.035, fontsize=10, color="#1f2937", bullet=False):
    lines = []
    for paragraph in text.split("\n"):
        if not paragraph.strip():
            lines.append("")
            continue
        prefix = "• " if bullet else ""
        wrapped = textwrap.wrap(paragraph, width=width) or [paragraph]
        for idx, segment in enumerate(wrapped):
            lines.append(prefix + segment if idx == 0 else ("  " + segment if bullet else segment))
    for idx, line in enumerate(lines):
        ax.text(x, y - idx * line_height, line, fontsize=fontsize, color=color, va="top", ha="left")
    return y - len(lines) * line_height


def add_footer(fig, page_no):
    fig.text(0.02, 0.015, PAGE_FOOTER, fontsize=8, color="#64748b", ha="left")
    fig.text(0.98, 0.015, f"Page {page_no}", fontsize=8, color="#64748b", ha="right")


def new_page(page_size):
    if page_size == "a4":
        return plt.figure(figsize=(8.27, 11.69), facecolor="#f8f6f0")
    return plt.figure(figsize=(8.5, 11.0), facecolor="#f8f6f0")


def render_cover_page(pdf, content, page_no, page_size):
    fig = new_page(page_size)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.axis("off")
    title = clean_text(content.get("project_title"))
    working = clean_text(content.get("working_title"))
    take_home = clean_text(content.get("one_sentence_take_home"))
    roles = clean_text(content.get("benchmark_roles"))
    ax.text(0.08, 0.88, title, fontsize=24, fontweight="bold", color="#16324f", va="top")
    ax.text(0.08, 0.80, "Presenter preparation pack", fontsize=14, color="#475569", va="top")
    ax.text(0.08, 0.72, "Working manuscript title", fontsize=11, fontweight="bold", color="#16324f", va="top")
    wrap_text_to_axes(ax, working, 0.08, 0.69, width=70, line_height=0.038, fontsize=12)
    ax.text(0.08, 0.56, "One-sentence take-home", fontsize=11, fontweight="bold", color="#16324f", va="top")
    wrap_text_to_axes(ax, take_home, 0.08, 0.53, width=78, line_height=0.038, fontsize=12)
    ax.text(0.08, 0.34, "Benchmark-role reminder", fontsize=11, fontweight="bold", color="#16324f", va="top")
    wrap_text_to_axes(ax, roles, 0.08, 0.31, width=78, line_height=0.04, fontsize=11)
    ax.text(0.08, 0.15, "Use this PDF to prepare how to speak about the work, not as a substitute scientific authority.", fontsize=10, color="#475569")
    add_footer(fig, page_no)
    pdf.savefig(fig)
    plt.close(fig)


def render_text_sections(pdf, title, lines, page_no_start, page_size, subtitle=None):
    body_lines = []
    for line in lines:
        stripped = line.strip()
        if not stripped:
            body_lines.append("")
        elif stripped.startswith("- "):
            body_lines.append("• " + clean_text(stripped[2:]))
        else:
            body_lines.append(clean_text(stripped))
    max_lines = 28
    pages = [body_lines[i:i + max_lines] for i in range(0, len(body_lines), max_lines)] or [[]]
    page_no = page_no_start
    for idx, chunk in enumerate(pages):
        fig = new_page(page_size)
        ax = fig.add_axes([0, 0, 1, 1])
        ax.axis("off")
        ax.text(0.08, 0.93, title, fontsize=20, fontweight="bold", color="#16324f", va="top")
        if subtitle:
            ax.text(0.08, 0.89, subtitle, fontsize=11, color="#475569", va="top")
        y = 0.84
        for line in chunk:
            if line == "":
                y -= 0.02
            else:
                y = wrap_text_to_axes(ax, line, 0.08, y, width=92, line_height=0.028, fontsize=10)
                y -= 0.008
        add_footer(fig, page_no)
        pdf.savefig(fig)
        plt.close(fig)
        page_no += 1
    return page_no


def render_summary_pages(pdf, sections_by_title, page_no, page_size):
    for key in ["How to Use This Packet", "One-Page Talk Summary", "Benchmark Hierarchy Cheat Sheet"]:
        section = sections_by_title[key]
        page_no = render_text_sections(pdf, section["title"], section["lines"], page_no, page_size)
    return page_no


def render_talk_flow_pages(pdf, section, page_no, page_size):
    for sub in section.get("subsections", []):
        info = parse_bullet_map(sub["lines"])
        fig = new_page(page_size)
        ax = fig.add_axes([0, 0, 1, 1])
        ax.axis("off")
        ax.text(0.08, 0.93, "Talk Flow Overview", fontsize=20, fontweight="bold", color="#16324f", va="top")
        ax.text(0.08, 0.88, sub["title"], fontsize=14, fontweight="bold", color="#9a5c1f", va="top")
        y = 0.82
        for label, key in [
            ("Objective", "objective"),
            ("Key line to say", "key_line_to_say"),
            ("Key line not to say", "key_line_not_to_say"),
            ("Transition", "transition"),
        ]:
            ax.text(0.08, y, label, fontsize=11, fontweight="bold", color="#16324f", va="top")
            y = wrap_text_to_axes(ax, clean_text(info.get(key, "")), 0.08, y - 0.03, width=90, line_height=0.032, fontsize=10)
            y -= 0.03
        add_footer(fig, page_no)
        pdf.savefig(fig)
        plt.close(fig)
        page_no += 1
    return page_no


def render_slide_note_pages(pdf, section, page_no, page_size):
    for sub in section.get("subsections", []):
        info = parse_bullet_map(sub["lines"])
        fig = new_page(page_size)
        ax = fig.add_axes([0, 0, 1, 1])
        ax.axis("off")
        ax.text(0.08, 0.93, "Slide-by-Slide Speaker Notes", fontsize=20, fontweight="bold", color="#16324f", va="top")
        ax.text(0.08, 0.88, sub["title"], fontsize=14, fontweight="bold", color="#9a5c1f", va="top")
        y = 0.82
        fields = [
            ("Objective", "objective"),
            ("30-second version", "short_version"),
            ("2-minute version", "long_version"),
            ("Key evidence", "key_evidence"),
            ("Common audience misunderstanding", "common_audience_misunderstanding"),
            ("Safe clarification", "safe_clarification"),
            ("Transition to next slide", "transition_to_next_slide"),
        ]
        for label, key in fields:
            ax.text(0.08, y, label, fontsize=10.5, fontweight="bold", color="#16324f", va="top")
            y = wrap_text_to_axes(ax, clean_text(info.get(key, "")), 0.08, y - 0.028, width=92, line_height=0.028, fontsize=9.6)
            y -= 0.02
        add_footer(fig, page_no)
        pdf.savefig(fig)
        plt.close(fig)
        page_no += 1
    return page_no


def load_captions_map(path):
    titles = {}
    current = None
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith("### Figure"):
            current = line[4:].strip()
            titles[current.split(".")[0].replace("Figure ", "").replace("S", "S")] = current
    return titles


def render_figure_with_notes(pdf, figure_path, heading, note_map, page_no, page_size):
    fig = new_page(page_size)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.axis("off")
    ax.text(0.08, 0.95, "Figure-by-Figure Guide", fontsize=20, fontweight="bold", color="#16324f", va="top")
    ax.text(0.08, 0.91, heading, fontsize=13, fontweight="bold", color="#9a5c1f", va="top")

    image_ax = fig.add_axes([0.08, 0.47, 0.84, 0.36])
    image_ax.axis("off")
    if figure_path.exists():
        image = mpimg.imread(str(figure_path))
        image_ax.imshow(image)
        image_ax.set_aspect("auto")
    else:
        image_ax.text(0.5, 0.5, f"Missing figure: {figure_path.name}", ha="center", va="center", fontsize=12, color="#b42318")

    y = 0.42
    for label, key in [
        ("Stem", "stem"),
        ("Visible title", "visible_title"),
        ("Benchmark role", "benchmark_role"),
        ("What the audience should notice first", "notice_first"),
        ("What to say", "what_to_say"),
        ("What not to say", "what_not_to_say"),
        ("Support if challenged", "support_if_challenged"),
    ]:
        ax.text(0.08, y, label, fontsize=10.5, fontweight="bold", color="#16324f", va="top")
        y = wrap_text_to_axes(ax, clean_text(note_map.get(key, "")), 0.08, y - 0.026, width=92, line_height=0.025, fontsize=9.4)
        y -= 0.014
    add_footer(fig, page_no)
    pdf.savefig(fig)
    plt.close(fig)
    return page_no + 1


def resolve_figure_image(stem, figures_dir, hub_outputs_dir=None):
    if hub_outputs_dir is not None:
        hub_png = hub_outputs_dir / f"{stem}.png"
        if hub_png.exists():
            return hub_png
    return figures_dir / f"{stem}.png"


def render_figure_pages(pdf, section, figures_dir, hub_outputs_dir, page_no, page_size):
    for sub in section.get("subsections", []):
        note_map = parse_bullet_map(sub["lines"])
        image_override = clean_text(note_map.get("image_path", ""))
        if image_override:
            image_path = Path(image_override)
        else:
            stem = clean_text(note_map.get("stem", "")).strip("`")
            image_path = resolve_figure_image(stem, figures_dir, hub_outputs_dir)
        page_no = render_figure_with_notes(pdf, image_path, sub["title"], note_map, page_no, page_size)
    return page_no


def render_qa_pages(pdf, qa_items, page_no, page_size):
    for item in qa_items:
        fig = new_page(page_size)
        ax = fig.add_axes([0, 0, 1, 1])
        ax.axis("off")
        ax.text(0.08, 0.93, "Likely Questions and Answers", fontsize=20, fontweight="bold", color="#16324f", va="top")
        ax.text(0.08, 0.89, item["question_id"], fontsize=12, fontweight="bold", color="#9a5c1f", va="top")
        y = 0.84
        fields = [
            ("Topic", "topic"),
            ("Why it will be asked", "why_it_will_be_asked"),
            ("Short answer", "short_answer"),
            ("Expanded answer", "expanded_answer"),
            ("Evidence root", "evidence_root"),
            ("Source files", "source_files"),
            ("Risk if misstated", "risk_if_misstated"),
            ("Forbidden answer patterns", "forbidden_answer_patterns"),
        ]
        for label, key in fields:
            ax.text(0.08, y, label, fontsize=10.5, fontweight="bold", color="#16324f", va="top")
            value = item.get(key, "")
            if isinstance(value, list):
                text = "\n".join(f"- {clean_text(v)}" for v in value)
            else:
                text = clean_text(value)
            y = wrap_text_to_axes(ax, text, 0.08, y - 0.026, width=92, line_height=0.026, fontsize=9.4)
            y -= 0.014
        add_footer(fig, page_no)
        pdf.savefig(fig)
        plt.close(fig)
        page_no += 1
    return page_no


def render_risk_pages(pdf, sections_by_title, page_no, page_size):
    for title in ["Dangerous Questions and Safe Framing", "Do-Not-Overclaim List", "Backup Data Points", "Final Checklist Before Speaking"]:
        page_no = render_text_sections(pdf, title, sections_by_title[title]["lines"], page_no, page_size)
    return page_no


def render_final_checklist_page(pdf, sections_by_title, page_no, page_size):
    section = sections_by_title["Final Checklist Before Speaking"]
    fig = new_page(page_size)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.axis("off")
    ax.text(0.08, 0.93, "Final Speaking Checklist", fontsize=20, fontweight="bold", color="#16324f", va="top")
    ax.text(0.08, 0.88, "Use this page in the final hour before the talk.", fontsize=11, color="#475569", va="top")
    y = 0.80
    for line in section["lines"]:
        stripped = line.strip()
        if stripped.startswith("- "):
            y = wrap_text_to_axes(ax, "• " + clean_text(stripped[2:]), 0.10, y, width=86, line_height=0.05, fontsize=12)
            y -= 0.012
    add_footer(fig, page_no)
    pdf.savefig(fig)
    plt.close(fig)
    return page_no + 1


def build_pdf(args, content_sections, qa_items, optional_docs):
    sections_by_title = section_map(content_sections)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    page_no = 1
    with PdfPages(args.output) as pdf:
        render_cover_page(pdf, parse_bullet_map(sections_by_title["Title Page"]["lines"]), page_no, args.page_size)
        page_no += 1
        page_no = render_summary_pages(pdf, sections_by_title, page_no, args.page_size)
        page_no = render_talk_flow_pages(pdf, sections_by_title["Talk Flow Overview"], page_no, args.page_size)
        page_no = render_slide_note_pages(pdf, sections_by_title["Slide-by-Slide Speaker Notes"], page_no, args.page_size)
        page_no = render_figure_pages(pdf, sections_by_title["Figure-by-Figure Interpretation Notes"], args.figures_dir, args.hub_outputs_dir, page_no, args.page_size)
        page_no = render_qa_pages(pdf, qa_items, page_no, args.page_size)
        page_no = render_risk_pages(pdf, sections_by_title, page_no, args.page_size)
        page_no = render_final_checklist_page(pdf, sections_by_title, page_no, args.page_size)

    manifest = {
        "generated_at": datetime.now().isoformat(),
        "content": str(args.content),
        "qa": str(args.qa),
        "figures_dir": str(args.figures_dir),
        "hub_outputs_dir": str(args.hub_outputs_dir),
        "fimm_figure": str(args.fimm_figure),
        "fimm_survival_report": str(args.fimm_survival_report),
        "figure_stems_used": FIGURE_ORDER,
        "benchmark_hierarchy": {
            "primary": "analysis/1000g_realdata/benchmark_wgs_wave2",
            "secondary": [
                "analysis/1000g_realdata/benchmark_wes_truthbacked/run",
                "analysis/1000g_realdata/benchmark_rna_truthbacked/run",
            ],
            "supplementary": "analysis/1000g_realdata/benchmark_trimodal_robustness/run",
        },
        "optional_docs_loaded": optional_docs,
        "output": str(args.output),
    }
    args.manifest.parent.mkdir(parents=True, exist_ok=True)
    args.manifest.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")


def main():
    args = parse_args()
    content_sections = load_markdown_sections(args.content)
    qa_items = load_qa_bank(args.qa)
    repo_root = Path(__file__).resolve().parent.parent
    optional_paths = {
        "presentation_prompt": repo_root / "docs" / "PRESENTATION_AI_PROMPT.md",
        "manuscript": repo_root / "docs" / "MANUSCRIPT_DRAFT_V1.md",
        "evidence_map": repo_root / "docs" / "MANUSCRIPT_EVIDENCE_MAP.md",
        "captions": args.figures_dir / "captions.md",
        "figures_readme": args.figures_dir / "README.md",
        "fimm_survival_report": args.fimm_survival_report,
    }
    optional_docs = {key: str(path) for key, path in optional_paths.items() if path.exists()}
    build_pdf(args, content_sections, qa_items, optional_docs)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
