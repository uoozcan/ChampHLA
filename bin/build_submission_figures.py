#!/usr/bin/env python3
"""Assemble a journal-ready figure set with caption-consistent filenames.

The manuscript keeps internal build-provenance annotations `(Source:
figures_final/<stem>)` in each caption; the on-disk figure files retain their
generator stems, so a display "Figure 5" may map to `figure_6_discordance...`.
This script parses the manuscript, resolves each display number to its source
file, and copies the figures into `submission_figures/` renamed to match the
captions (Figure_1.pdf … Figure_13.pdf, Supplementary_Figure_S1.pdf …), ready to
upload. The working repo and captions are left untouched.

Usage:  python3 bin/build_submission_figures.py
"""
import re
import shutil
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
MD = REPO / "docs" / "CHAMPHLA_MANUSCRIPT_V3.md"
ANALYSIS = REPO / "analysis"
OUT = ANALYSIS / "submission_figures"

CAP = re.compile(r"\*\*(Figure\s+(\d+)|Supplementary\s+Figure\s+(S\d+))")
# capture "(Source: `<figures_final|figures_final_candidate>/<stem>`)" in a caption
SRC = re.compile(r"\(Source:[^)]*?`(figures_final[a-z_]*)/([^`]+)`")


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    text = MD.read_text(encoding="utf-8")
    mapped, missing = [], []
    for line in text.splitlines():
        cap = CAP.search(line)
        src = SRC.search(line)
        if not (cap and src):
            continue
        label = ("Figure_" + cap.group(2)) if cap.group(2) else ("Supplementary_Figure_" + cap.group(3))
        srcdir, stem = src.group(1).strip(), src.group(2).strip()
        copied_any = False
        for ext in ("pdf", "png", "svg"):
            f = ANALYSIS / srcdir / f"{stem}.{ext}"
            if f.exists():
                shutil.copy2(f, OUT / f"{label}.{ext}")
                copied_any = True
        (mapped if copied_any else missing).append((label, stem))

    for label, stem in sorted(mapped, key=lambda x: (len(x[0]), x[0])):
        print(f"  {label:28s} <- {stem}")
    if missing:
        print("\nMISSING source files (regenerate figures first):")
        for label, stem in missing:
            print(f"  {label:28s} <- {stem}  (not found)")
    print(f"\n{len(mapped)} figures written to {OUT}")


if __name__ == "__main__":
    main()
