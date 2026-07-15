# MVHLA Figure Hub

Centralised location for all publication figure scripts and their outputs.

## Directory layout

```
figure_hub/
├── plot_style.py       ← import this in EVERY figure script
├── README.md           ← this file
├── scripts/            ← all figure-generating scripts
│   ├── fig_benchmark_main.py      — Figures 1–7, 9, 10, S1–S3 (benchmark)
│   ├── fig_loh_dropout_heatmap.py — Figure LOH: allele dropout heatmap (2-panel)
│   ├── fig_fimm_loh_wes.py        — Figure FIMM LOH WES (stacked bar)
│   └── fig_fimm_concordance.py    — FIMM concordance heatmaps + cross-modal plot
└── outputs/            ← symlinks to all final figures (PNG / PDF / SVG)
    └── README.md       — maps each figure to its generating script and source data
```

## Style standards (enforced by plot_style.py)

| Standard | Setting |
|---|---|
| Font | Arial (fallback: Helvetica, DejaVu Sans) |
| Legends | Always outside axes — use `legend_outside(ax)` or `legend_below(ax)` |
| Text overlap | Rotate x-tick labels ≥45° when > 6 labels; use `tight_with_legend()` |
| Save DPI | 300 DPI for PNG/PDF/SVG |
| Color palette | Wong 2011 colorblind-safe (`WONG` dict) |
| Formats | PDF + SVG + PNG via `save_fig(fig, path_stem)` |

## Quickstart — add a new figure script

```python
import sys
sys.path.insert(0, ".../figure_hub")
from plot_style import apply_style, legend_outside, save_fig, WONG, TOOL_COLORS

apply_style()               # call ONCE before any plt call

import matplotlib.pyplot as plt
import pandas as pd

# ... your figure code ...

ax.set_title("My Figure", fontweight="bold")
legend_outside(ax)          # legend to the right of axes
save_fig(fig, "figure_hub/outputs/figure_my_new_figure")
```

Preferred repo-local variant:

```python
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from plot_style import apply_style, legend_outside, save_fig, WONG
apply_style()
```

Then add a symlink row to `outputs/README.md`.

## FIMM sample IDs

All FIMM figures must use **de-identified files** (`*_deident.tsv`) as input.
Real ID → pseudo-ID mapping: `fimm_results/fimm_pseudo_id_map.tsv`

To regenerate mapping and de-identified files:
```bash
python3.11 bin/create_fimm_pseudo_ids.py
```

## Compatibility wrappers

| bin/ script | Hub equivalent |
|---|---|
| bin/generate_publication_figures.py | fig_benchmark_main.py |
| bin/analysis/figure_loh_allele_dropout_heatmap.py | fig_loh_dropout_heatmap.py |
| bin/analysis/analyze_fimm_loh_wes.py | fig_fimm_loh_wes.py (figure part only) |
| bin/analyze_fimm_hla.py | fig_fimm_concordance.py |

For new or modified figures, edit the hub copy in `scripts/`, not the bin/ originals.
Repo-tracked scripts outside `figure_hub/scripts/` should only remain as thin wrappers or orchestration entrypoints.
