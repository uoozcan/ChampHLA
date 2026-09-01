"""Console-entry wrappers for the isolated source-tree commands."""

from __future__ import annotations

import runpy
from pathlib import Path


def _dispatch(filename: str):
    script = Path(__file__).resolve().parents[2] / "scripts" / filename
    if not script.is_file():
        raise SystemExit(f"RefFormer source-tree script not found: {script}")
    runpy.run_path(str(script), run_name="__main__")


def train_refformer():
    _dispatch("train_refformer.py")


def predict_refformer():
    _dispatch("predict_refformer.py")


def nested_evaluation():
    _dispatch("evaluate_nested.py")
