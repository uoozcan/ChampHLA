#!/usr/bin/env python3
"""Run the 1000 Genomes benchmark and refreshed figure workflow."""

import argparse
import subprocess
from pathlib import Path


def parse_args():
    parser = argparse.ArgumentParser(description="Run the 1000G benchmark and generate refreshed figures.")
    parser.add_argument("--config", required=True, help="YAML config describing manifests, truth, and result globs.")
    parser.add_argument("--output-dir", required=True, help="Output directory for benchmark tables and figures.")
    parser.add_argument("--figures", default=None, help="Optional comma-separated benchmark figure keys to render.")
    parser.add_argument("--exploratory-figures", default="x1,x2,x3", help="Optional comma-separated exploratory figure keys to render.")
    return parser.parse_args()


def main():
    args = parse_args()
    script_dir = Path(__file__).resolve().parent
    repo_root = script_dir.parent
    config_path = Path(args.config).resolve()
    output_dir = Path(args.output_dir).resolve()
    subprocess.run([
        "python3",
        str((script_dir / "run_1000g_benchmark.py").resolve()),
        "--config",
        str(config_path),
        "--output-dir",
        str(output_dir),
    ], check=True, cwd=str(repo_root))

    cmd = [
        "python3.11",
        str((script_dir / "generate_figures.py").resolve()),
        "--tables-dir",
        str(output_dir / "tables"),
        "--figures-dir",
        str(output_dir / "figures"),
        "--source-dir",
        str(output_dir),
        "--include-exploratory",
        "always",
        "--exploratory-figures",
        args.exploratory_figures,
    ]
    if args.figures:
        cmd.extend(["--figures", args.figures])
    subprocess.run(cmd, check=True, cwd=str(repo_root))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
