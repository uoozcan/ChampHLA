#!/usr/bin/env python3
"""Generate, validate, archive, and promote the active publication figure package."""

import argparse
import shutil
import subprocess
from datetime import date
from pathlib import Path


EXPECTED_ACTIVE_STEMS = [
    "figure_1_workflow_architecture",
    "figure_2_accuracy_comparison",
    "figure_3_per_gene_gains",
    "figure_4_confidence_calibration",
    "figure_5_abstention_tradeoff",
    "figure_6_discordance_taxonomy",
    "figure_7_confidence_weights",
    "figure_09_bimodal_per_gene",
    "figure_10_trimodal_comparison",
    "figure_s1_per_gene_accuracy",
    "figure_s2_calibration_heatmap",
    "figure_s3_resolution_comparison",
]


def parse_args() -> argparse.Namespace:
    default_out = Path("/scratch/project_2008084/pihla-publish/analysis/figures_final")
    default_stage = Path("/scratch/project_2008084/pihla-publish/analysis/figures_final_staging")
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=default_out)
    parser.add_argument("--staging-dir", type=Path, default=default_stage)
    parser.add_argument("--keep-staging", action="store_true")
    return parser.parse_args()


def validate_stage(staging_dir: Path) -> None:
    for stem in EXPECTED_ACTIVE_STEMS:
        for ext in ("pdf", "png", "svg"):
            path = staging_dir / f"{stem}.{ext}"
            if not path.exists():
                raise FileNotFoundError(f"Missing staged asset: {path}")
    for name in ("captions.md", "README.md"):
        if not (staging_dir / name).exists():
            raise FileNotFoundError(f"Missing staged package file: {staging_dir / name}")


def main() -> int:
    args = parse_args()
    script_dir = Path(__file__).resolve().parent
    repo_root = script_dir.parent
    stage = args.staging_dir
    out_dir = args.output_dir
    if stage.exists():
        shutil.rmtree(stage)
    stage.mkdir(parents=True, exist_ok=True)

    hub_script = (repo_root / "figure_hub" / "scripts" / "fig_benchmark_main.py").resolve()
    subprocess.run(
        ["python3.11", str(hub_script), "--out-dir", str(stage)],
        check=True,
        cwd=str(repo_root),
    )

    validate_stage(stage)

    archive_dir = out_dir / f"legacy_pre_redesign_{date.today().isoformat()}"
    archive_dir.mkdir(parents=True, exist_ok=True)

    if out_dir.exists():
        for entry in out_dir.iterdir():
            if entry.name.startswith("legacy_pre_"):
                continue
            target = archive_dir / entry.name
            if target.exists():
                if target.is_dir():
                    shutil.rmtree(target)
                else:
                    target.unlink()
            shutil.move(str(entry), str(target))

    out_dir.mkdir(parents=True, exist_ok=True)
    for entry in stage.iterdir():
        shutil.move(str(entry), str(out_dir / entry.name))

    if not args.keep_staging and stage.exists():
        shutil.rmtree(stage)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
