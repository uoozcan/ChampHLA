#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import textwrap
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/champhla-matplotlib")
os.environ.setdefault("SOURCE_DATE_EPOCH", "946684800")
import matplotlib
matplotlib.use("Agg")
matplotlib.rcParams["svg.hashsalt"] = "champhla-publication"
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch


PALETTE = {
    "blue": "#0072B2", "orange": "#E69F00", "sky": "#56B4E9",
    "green": "#009E73", "yellow": "#F0E442", "navy": "#332288",
    "pink": "#CC79A7", "grey": "#777777", "light": "#F4F4F4",
}


def digest(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def write_tsv(path: Path, rows: list[dict], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def workflow_source(root: Path, path: Path) -> None:
    protocol = json.loads((root / "configs/confirmation_protocol.json").read_text())
    panel = " / ".join(protocol["modalities"]).upper()
    rows = [
        {"record_type": "node", "id": "inputs", "label": f"Truth-free inputs\n{panel}", "x": 0.10, "y": 0.68, "category": "input", "from": "", "to": "", "style": "solid"},
        {"record_type": "node", "id": "callers", "label": "Assay-compatible\ncaller panels", "x": 0.30, "y": 0.68, "category": "process", "from": "", "to": "", "style": "solid"},
        {"record_type": "node", "id": "native", "label": "Native outputs\nand source hashes", "x": 0.50, "y": 0.68, "category": "evidence", "from": "", "to": "", "style": "solid"},
        {"record_type": "node", "id": "harmonise", "label": "Complete / partial / missing\nharmonisation", "x": 0.70, "y": 0.68, "category": "process", "from": "", "to": "", "style": "solid"},
        {"record_type": "node", "id": "plurality", "label": "Pair-level plurality\nSimplePluralityLex", "x": 0.90, "y": 0.68, "category": "primary", "from": "", "to": "", "style": "solid"},
        {"record_type": "node", "id": "freeze", "label": "Prediction freeze\ncode + inputs + outputs", "x": 0.78, "y": 0.41, "category": "freeze", "from": "", "to": "", "style": "solid"},
        {"record_type": "node", "id": "truth", "label": "Segregated truth\nno prediction access", "x": 0.18, "y": 0.18, "category": "truth", "from": "", "to": "", "style": "solid"},
        {"record_type": "node", "id": "join", "label": "One-time non-overwriting\ntruth join", "x": 0.42, "y": 0.18, "category": "join", "from": "", "to": "", "style": "solid"},
        {"record_type": "node", "id": "evaluate", "label": "Donor-clustered evaluation\nand independent recount", "x": 0.66, "y": 0.18, "category": "evaluation", "from": "", "to": "", "style": "solid"},
        {"record_type": "node", "id": "registry", "label": "Result registry →\ntables, figures, prose", "x": 0.90, "y": 0.18, "category": "report", "from": "", "to": "", "style": "solid"},
    ]
    for source, target, style in [
        ("inputs", "callers", "solid"), ("callers", "native", "solid"),
        ("native", "harmonise", "solid"), ("harmonise", "plurality", "solid"),
        ("plurality", "freeze", "solid"), ("freeze", "join", "solid"),
        ("truth", "join", "solid"), ("join", "evaluate", "solid"),
        ("evaluate", "registry", "solid"), ("truth", "callers", "forbidden"),
    ]:
        rows.append({"record_type": "edge", "id": f"{source}_{target}", "label": "", "x": "", "y": "", "category": "flow", "from": source, "to": target, "style": style})
    write_tsv(path, rows, ["record_type", "id", "label", "x", "y", "category", "from", "to", "style"])


def cohort_source(root: Path, path: Path) -> None:
    roles = json.loads((root / "configs/dataset_roles.json").read_text())["datasets"]
    scales = {
        "1000G_DEVELOPMENT_WES_RNA": "WES 130; RNA 107 donors",
        "1000G_CORRECTED_WGS_WES_RNA": "WGS 137; WES 130; RNA 107 donors",
        "HPRC_RELEASE2_WGS": "target 120 donors",
        "NCI60_RNA": "11 strict subjects",
        "CCLE_NCI60_WES_OVERLAP": "3 overlapping cell lines",
        "AFGR_MKK_RNA_WGS": "RNA 166; WGS unresolved",
        "FNLCR_PBMC_RNA": "96 reported; controlled",
        "DICE_SORTED_CELL_RNA": "91 reported; controlled",
        "TRANSPLANT_GERMLINE_WES": "donor count unresolved",
        "REJECTED_OR_CIRCULAR_TRUTH": "excluded",
    }
    truth = {
        "1000G_DEVELOPMENT_WES_RNA": "2014 laboratory HLA truth",
        "1000G_CORRECTED_WGS_WES_RNA": "same frozen laboratory truth",
        "HPRC_RELEASE2_WGS": "dual-method phased-assembly truth pending",
        "NCI60_RNA": "independent sequence-based typing",
        "CCLE_NCI60_WES_OVERLAP": "complete A/B/C truth absent",
        "AFGR_MKK_RNA_WGS": "controlled MiSeq truth; mapping pending",
        "FNLCR_PBMC_RNA": "controlled original Sanger truth pending",
        "DICE_SORTED_CELL_RNA": "controlled DNA HLA file unverified",
        "TRANSPLANT_GERMLINE_WES": "match status is not allele truth",
        "REJECTED_OR_CIRCULAR_TRUTH": "rejected or circular",
    }
    rows = []
    for row in roles:
        rows.append({
            "dataset_id": row["dataset_id"], "role": row["role"],
            "independence_stratum": row["independence_stratum"],
            "scale": scales[row["dataset_id"]], "truth_source": truth[row["dataset_id"]],
            "accuracy_eligible": str(row["accuracy_eligible"]).lower(),
            "status": row["status"], "pooling": "reported separately; never pooled across strata",
        })
    write_tsv(path, rows, ["dataset_id", "role", "independence_stratum", "scale", "truth_source", "accuracy_eligible", "status", "pooling"])


def save(fig, stem: Path, width: float, height: float) -> list[Path]:
    stem.parent.mkdir(parents=True, exist_ok=True)
    paths = [stem.with_suffix(".svg"), stem.with_suffix(".pdf"), stem.with_suffix(".png")]
    fig.set_size_inches(width, height)
    fig.savefig(paths[0], format="svg", metadata={"Date": None})
    fig.savefig(paths[1], format="pdf", metadata={"CreationDate": None, "ModDate": None, "Creator": "ChampHLA"})
    fig.savefig(paths[2], format="png", dpi=300, metadata={"Software": "ChampHLA"})
    plt.close(fig)
    return paths


def draw_workflow(source: Path, output: Path) -> list[Path]:
    rows = list(csv.DictReader(source.open(), delimiter="\t"))
    nodes = {row["id"]: row for row in rows if row["record_type"] == "node"}
    fig, ax = plt.subplots()
    ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off")
    colors = {"input": PALETTE["sky"], "process": "#DDEBF7", "evidence": "#E8E8E8",
              "primary": PALETTE["orange"], "freeze": PALETTE["yellow"],
              "truth": PALETTE["pink"], "join": "#DDD3F3", "evaluation": "#CDEBDD",
              "report": "#C8D8F0"}
    for row in nodes.values():
        x, y = float(row["x"]), float(row["y"])
        box = FancyBboxPatch((x - .082, y - .065), .164, .13, boxstyle="round,pad=0.010",
                             linewidth=1.1, edgecolor="#333333", facecolor=colors[row["category"]])
        ax.add_patch(box)
        ax.text(x, y, row["label"], ha="center", va="center", fontsize=6.25)
    for row in (r for r in rows if r["record_type"] == "edge"):
        source_row, target_row = nodes[row["from"]], nodes[row["to"]]
        sx, sy = float(source_row["x"]), float(source_row["y"])
        tx, ty = float(target_row["x"]), float(target_row["y"])
        dx, dy = tx - sx, ty - sy
        if abs(dx) >= abs(dy):
            direction = 1 if dx > 0 else -1
            start, end = (sx + .084 * direction, sy), (tx - .084 * direction, ty)
        else:
            direction = 1 if dy > 0 else -1
            start, end = (sx, sy + .068 * direction), (tx, ty - .068 * direction)
        forbidden = row["style"] == "forbidden"
        arrow = FancyArrowPatch(start, end, arrowstyle="-|>", mutation_scale=10,
                                linewidth=1.2, linestyle="--" if forbidden else "-",
                                color=PALETTE["pink"] if forbidden else "#444444",
                                connectionstyle="arc3,rad=-0.12" if forbidden else "arc3,rad=0")
        ax.add_patch(arrow)
        if forbidden:
            ax.text(.20, .43, "forbidden truth flow", color=PALETTE["pink"], fontsize=6.4,
                    rotation=70, ha="center", va="center")
    ax.text(.02, .97, "ChampHLA workflow and truth firewall", fontsize=11, weight="bold", va="top")
    ax.text(.02, .92, "Truth cannot enter caller execution, harmonisation, consensus, or prediction freezing.", fontsize=7.5, va="top")
    return save(fig, output, 7.2, 4.5)


def draw_cohorts(source: Path, output: Path) -> list[Path]:
    rows = list(csv.DictReader(source.open(), delimiter="\t"))
    role_colors = {"development": PALETTE["blue"], "same_resource_confirmation": PALETTE["orange"],
                   "independent_validation": PALETTE["navy"], "exploratory": PALETTE["sky"],
                   "invalid": PALETTE["grey"]}
    fig, ax = plt.subplots()
    step = 1.28
    top = len(rows) * step
    ax.set_xlim(0, 1); ax.set_ylim(-1.0, top + 1.25); ax.axis("off")
    ax.text(.01, top + .95, "Cohort and evidence map", fontsize=11, weight="bold")
    ax.text(.01, top + .48, "Status is evidence state, not performance. Evidence strata are never pooled.", fontsize=7.5)
    ax.text(.07, top + .05, "Dataset, role, stratum, and scale", fontsize=6.2, color="#555555", weight="bold")
    ax.text(.49, top + .05, "Truth source and current lane state", fontsize=6.2, color="#555555", weight="bold")
    for index, row in enumerate(reversed(rows)):
        y = index * step
        ax.scatter(.035, y, s=48, marker="s", color=role_colors[row["role"]], edgecolor="#333333", linewidth=.5)
        ax.text(.07, y + .27, row["dataset_id"], fontsize=7.1, weight="bold", va="center")
        left = f"{row['role']} · {row['independence_stratum']} · {row['scale']}"
        ax.text(.07, y - .12, "\n".join(textwrap.wrap(left, 54)), fontsize=6.0, va="center")
        ax.text(.49, y + .28, "\n".join(textwrap.wrap(row["truth_source"], 48)), fontsize=6.2, va="center")
        ax.text(.49, y - .19, "\n".join(textwrap.wrap(row["status"], 62)), fontsize=5.55,
                color="#555555", va="center")
        ax.plot([.01, .99], [y - .58, y - .58], color="#DDDDDD", lw=.5)
    legend_y = -0.80
    x = .02
    for role in ("development", "same_resource_confirmation", "independent_validation", "exploratory", "invalid"):
        ax.scatter(x, legend_y, s=30, marker="s", color=role_colors[role])
        legend = {"same_resource_confirmation": "same-resource", "independent_validation": "independent"}.get(role, role)
        ax.text(x + .02, legend_y, legend.replace("_", " "), fontsize=5.7, va="center")
        x += .19
    return save(fig, output, 7.2, 7.4)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", default=".")
    parser.add_argument("--generation-commit", required=True)
    args = parser.parse_args()
    root = Path(args.root).resolve()
    base = root / "manuscripts/figures"
    source1 = base / "source_data/workflow_firewall.tsv"
    source2 = base / "source_data/cohort_evidence_map.tsv"
    workflow_source(root, source1)
    cohort_source(root, source2)
    script = Path(__file__).resolve()
    out1 = draw_workflow(source1, base / "main/Figure1_workflow_firewall")
    out2 = draw_cohorts(source2, base / "main/Figure2_cohort_evidence_map")
    rows = []
    for figure_id, title, source, upstream, outputs, width, height, role in [
        ("MAIN-1", "Workflow and truth firewall", source1, root / "configs/confirmation_protocol.json", out1, "7.2", "4.5", "methods"),
        ("MAIN-2", "Cohort and evidence map", source2, root / "configs/dataset_roles.json", out2, "7.2", "7.4", "evidence_architecture"),
    ]:
        rows.append({
            "figure_id": figure_id, "panel_id": "all", "title": title,
            "manuscript_role": role, "script": script.relative_to(root).as_posix(),
            "script_sha256": digest(script), "source_data": source.relative_to(root).as_posix(),
            "source_data_sha256": digest(source), "upstream_artifact": upstream.relative_to(root).as_posix(),
            "upstream_sha256": digest(upstream),
            "output_paths": ";".join(path.relative_to(root).as_posix() for path in outputs),
            "output_sha256": ";".join(digest(path) for path in outputs),
            "width_in": width, "height_in": height, "resolution": "300 dpi PNG; vector SVG/PDF",
            "evidence_role": "methodological" if figure_id == "MAIN-1" else "current_evidence_state",
            "generation_commit": args.generation_commit, "status": "GENERATED", "blocker": "",
        })
    for figure_id, title, role, blocker in [
        ("MAIN-3", "Corrected cross-assay performance", "primary_results", "corrected same-resource and independent evaluations are incomplete"),
        ("MAIN-4", "Paired method comparison and decision gates", "primary_inference", "frozen production evaluation and simultaneous intervals are incomplete"),
        ("MAIN-5", "Cross-assay conclusion and failure modes", "mechanism", "production call-state and tie/partial/missing evidence is incomplete"),
    ]:
        rows.append({
            "figure_id": figure_id, "panel_id": "all", "title": title,
            "manuscript_role": role, "script": "", "script_sha256": "",
            "source_data": "", "source_data_sha256": "", "upstream_artifact": "",
            "upstream_sha256": "", "output_paths": "", "output_sha256": "",
            "width_in": "", "height_in": "", "resolution": "", "evidence_role": "pending",
            "generation_commit": args.generation_commit, "status": "BLOCKED_EVIDENCE", "blocker": blocker,
        })
    write_tsv(base / "figure_manifest.tsv", rows, list(rows[0]))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
