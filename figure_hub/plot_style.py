"""
plot_style.py — shared matplotlib style module for all MVHLA figures.

Standards enforced:
  - Font: Arial (fallback: Helvetica, DejaVu Sans)
  - Legends: always outside axes (use legend_outside())
  - No overlapping text: use adjust_labels() or tight_layout() + subplots_adjust()
  - Resolution: 300 DPI for saved figures (PDF/SVG/PNG)
  - Color palette: Wong 2011 colorblind-safe (WONG dict)

Usage in any figure script:
    import sys
    sys.path.insert(0, "/scratch/project_2008084/pihla-publish/figure_hub")
    from plot_style import apply_style, legend_outside, save_fig, WONG
    apply_style()
"""

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
from pathlib import Path


# ---------------------------------------------------------------------------
# rcParams
# ---------------------------------------------------------------------------
ARIAL_RCPARAMS = {
    # Font
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 9,
    # Axes labels / titles
    "axes.titlesize": 10.5,
    "axes.titleweight": "bold",
    "axes.labelsize": 9,
    # Tick labels
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
    # Legend — always place outside via legend_outside()
    "legend.fontsize": 8,
    "legend.framealpha": 0.95,
    "legend.edgecolor": "#d8dee9",
    "legend.borderpad": 0.6,
    # Figure / save
    "figure.dpi": 150,
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
    "savefig.pad_inches": 0.08,
    # Spines & grid
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.grid": True,
    "axes.grid.axis": "y",
    "grid.alpha": 0.28,
    "grid.color": "#d7dde5",
    # Background
    "figure.facecolor": "#f8f6f0",
    "axes.facecolor": "#ffffff",
}


def apply_style() -> None:
    """Apply the shared MVHLA rcParams. Call once at script startup."""
    plt.rcParams.update(ARIAL_RCPARAMS)


# ---------------------------------------------------------------------------
# Legend helpers — keep legend outside axes so it never overlaps data
# ---------------------------------------------------------------------------

ROLE_COLORS = {
    "primary": "#16324f",
    "secondary_wes": "#2c5c9c",
    "secondary_rna": "#17806d",
    "supplementary": "#9a5c1f",
    "application": "#6d58c4",
}


METHOD_COLORS = {
    "OptiType": "#16324f",
    "WeightedConsensus": "#5e6b7c",
    "MajorityVote": "#b45309",
    "ChampionChallenger": "#2d728f",
    "GatedConsensus": "#8b5cf6",
    "LocusExpertConsensus": "#3b82f6",
    "BimodalMajorityVote": "#1d4ed8",
    "BimodalWeightedConsensus": "#60a5fa",
    "TrimodalMajorityVote": "#8b5cf6",
    "TrimodalWeightedConsensus": "#c084fc",
}


def legend_outside(
    ax,
    loc: str = "upper left",
    anchor: tuple = (1.02, 0.5),
    ncol: int = 1,
    title: str | None = None,
    **kwargs,
):
    """
    Place the legend outside the axes area.

    Call after all artists have been added to ax.
    The containing figure must be saved with bbox_inches='tight' (default in
    ARIAL_RCPARAMS / save_fig) so the legend is not clipped.

    Parameters
    ----------
    ax      : matplotlib Axes
    loc     : legend anchor point on the legend box
    bbox    : (x, y) in axes-fraction coordinates for the legend origin
    ncol    : number of legend columns
    title   : optional legend title
    """
    kwargs.setdefault("frameon", True)
    return ax.legend(
        loc=loc,
        bbox_to_anchor=anchor,
        borderaxespad=0,
        ncol=ncol,
        title=title,
        **kwargs,
    )


def legend_below(
    ax,
    ncol: int = 4,
    y: float = -0.18,
    title: str | None = None,
    **kwargs,
):
    """Place legend centered below the axes."""
    kwargs.setdefault("frameon", True)
    return ax.legend(
        loc="upper center",
        bbox_to_anchor=(0.5, y),
        borderaxespad=0,
        ncol=ncol,
        title=title,
        **kwargs,
    )


# ---------------------------------------------------------------------------
# Save helper
# ---------------------------------------------------------------------------

def save_fig(
    fig,
    path_stem: str | Path,
    dpi: int = 300,
    formats: tuple = ("pdf", "svg", "png"),
    close: bool = True,
) -> None:
    """
    Save figure in multiple formats and close it.

    Parameters
    ----------
    fig       : matplotlib Figure
    path_stem : path without extension, e.g. "analysis/figures_final/figure_2"
    formats   : tuple of format strings
    """
    path_stem = Path(path_stem)
    path_stem.parent.mkdir(parents=True, exist_ok=True)
    for fmt in formats:
        fig.savefig(str(path_stem) + f".{fmt}", format=fmt, dpi=dpi)
    if close:
        plt.close(fig)


# ---------------------------------------------------------------------------
# Wong 2011 colorblind-safe palette
# ---------------------------------------------------------------------------
WONG = {
    "black":   "#000000",
    "orange":  "#E69F00",
    "sky":     "#56B4E9",
    "green":   "#009E73",
    "yellow":  "#F0E442",
    "blue":    "#0072B2",
    "vermil":  "#D55E00",
    "purple":  "#CC79A7",
}

# Convenience list (same order as Wong 2011 Fig 2)
WONG_CYCLE = [
    WONG["orange"],
    WONG["sky"],
    WONG["green"],
    WONG["yellow"],
    WONG["blue"],
    WONG["vermil"],
    WONG["purple"],
    WONG["black"],
]

# Modality colours used throughout MVHLA
MOD_COLORS = {
    "wgs":    WONG["orange"],
    "wes":    WONG["blue"],
    "rna":    WONG["green"],
    "rnaseq": WONG["green"],
}

# Tool colours used throughout MVHLA
TOOL_COLORS = {
    "ArcasHLA":  WONG["orange"],
    "HLA-HD":    WONG["sky"],
    "Kourami":   WONG["green"],
    "OptiType":  WONG["yellow"],
    "POLYSOLVER": WONG["blue"],
    "Seq2HLA":   WONG["vermil"],
    "SpecHLA":   WONG["purple"],
    "T1K":       WONG["black"],
}


# ---------------------------------------------------------------------------
# Text overlap prevention helpers
# ---------------------------------------------------------------------------

def rotate_xticklabels(ax, angle: int = 20, ha: str = "right") -> None:
    """Rotate x-tick labels to prevent overlap."""
    ax.set_xticklabels(ax.get_xticklabels(), rotation=angle, ha=ha)


def avoid_xticklabel_overlap(ax, rotation: int = 45, ha: str = "right") -> None:
    rotate_xticklabels(ax, angle=rotation, ha=ha)


def set_figure_header(
    fig,
    title: str,
    role_text: str | None = None,
    role_color: str | None = None,
    title_y: float = 0.975,
    role_y: float = 0.938,
) -> None:
    fig.suptitle(title, fontsize=11, fontweight="bold", y=title_y)
    if role_text:
        fig.text(
            0.06,
            role_y,
            role_text,
            fontsize=8,
            fontweight="bold",
            color=role_color or ROLE_COLORS["primary"],
            ha="left",
            va="center",
        )


def set_figure_footer(fig, text: str, x: float = 0.06, y: float = 0.045) -> None:
    fig.text(x, y, text, fontsize=7, color="#475569", ha="left", va="bottom")


def tight_with_legend(
    fig,
    top: float = 0.90,
    bottom: float = 0.12,
    left: float = 0.08,
    right: float = 0.82,
) -> None:
    """
    Call instead of fig.tight_layout() when a right-side legend is present.
    Reserves space on the right so the legend is not clipped.
    """
    fig.tight_layout(rect=(left, bottom, right, top))


# ---------------------------------------------------------------------------
# Quick-start template (printed when script is run directly)
# ---------------------------------------------------------------------------
if __name__ == "__main__":
    print(__doc__)
    print(f"WONG palette entries: {list(WONG.keys())}")
    apply_style()
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.bar(list(WONG.keys()), range(len(WONG)), color=list(WONG.values()))
    legend_outside(ax)
    ax.set_title("plot_style.py — palette test")
    save_fig(fig, "/tmp/plot_style_test")
    print("Test figure saved to /tmp/plot_style_test.{pdf,svg,png}")
