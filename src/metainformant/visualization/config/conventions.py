"""Shared chart-convention enforcement utilities.

Centralizes the plotting conventions used by the visualization module and by
downstream consumers (e.g. rna/analysis atlas figures, campaign figures):

- Okabe-Ito color-blind-safe palette (single canonical definition).
- Hatch patterns paired with categorical fills so category identity survives
  grayscale printing and color-vision deficiency (hatch + text redundancy).
- Explicit zero marks on logarithmic axes (a log scale silently clamps zero
  to the axis floor; an unmarked floor invites misreading).
- Deterministic figure export: matplotlib PNG/PDF output embeds timestamps and
  host metadata by default; deterministic mode strips them so byte-identical
  inputs produce byte-identical files.

All statistics and labels rendered by consumers are DESCRIPTIVE ONLY; nothing
here computes inferential results.
"""

from __future__ import annotations

from typing import Dict, List, Sequence, Tuple

import matplotlib as mpl

from metainformant.core.utils import logging

logger = logging.get_logger(__name__)

# ---------------------------------------------------------------------------
# Okabe-Ito color-blind-safe palette (fixed order; deterministic).
# ---------------------------------------------------------------------------

OKABE_ITO: Tuple[str, ...] = (
    "#E69F00",  # orange
    "#56B4E9",  # sky blue
    "#009E73",  # bluish green
    "#F0E442",  # yellow
    "#0072B2",  # blue
    "#D55E00",  # vermillion
    "#CC79A7",  # reddish purple
    "#000000",  # black
)

# Hatch patterns paired with Okabe-Ito fills for redundant categorical
# encoding. Doubled characters render at higher density in publication-scale
# figures; order is fixed so legend construction is deterministic.
CATEGORY_HATCHES: Tuple[str, ...] = ("", "//", "xx", "\\\\", "..", "++", "oo", "**")


def okabe_ito(n: int) -> List[str]:
    """Return the first *n* Okabe-Ito colors (cycling if n > 8).

    Args:
        n: Number of colors needed (must be positive).

    Raises:
        ValueError: If n is not positive.
    """
    if n < 1:
        raise ValueError("n must be a positive integer")
    return [OKABE_ITO[i % len(OKABE_ITO)] for i in range(n)]


def category_style(labels: Sequence[str]) -> Dict[str, Tuple[str, str]]:
    """Map each label to a deterministic (color, hatch) style pair.

    Labels are assigned styles in sorted order over the unique labels, so the
    same label set always yields the same mapping.

    Args:
        labels: Category labels to style.

    Returns:
        Dict mapping label -> (Okabe-Ito color, hatch pattern).
    """
    unique = sorted(set(str(label) for label in labels))
    return {
        label: (OKABE_ITO[i % len(OKABE_ITO)], CATEGORY_HATCHES[i % len(CATEGORY_HATCHES)])
        for i, label in enumerate(unique)
    }


# ---------------------------------------------------------------------------
# Log-axis zero marks
# ---------------------------------------------------------------------------


def add_log_zero_mark(
    ax: "mpl.axes.Axes",
    axis: str = "y",
    *,
    label: str = "zero",
    color: str = "#666666",
) -> None:
    """Draw an explicit zero reference on a logarithmic axis.

    Matplotlib log scales clamp zero to the axis floor, so a dataset that
    includes zeros is silently clipped unless the axis floor is marked. This
    helper draws a dashed reference line at the current axis floor (where
    zero would be rendered) and appends an explanatory note to the axis
    label, satisfying the "explicit zero marks on logarithmic axes"
    convention.

    Args:
        ax: Axes whose log axis should be annotated.
        axis: Which axis is logarithmic: ``"y"`` or ``"x"``.
        label: Text describing the zero mark (appended to the axis label).
        color: Line color for the zero reference.

    Raises:
        ValueError: If the given axis is not log-scaled (refusing to draw a
            spurious zero mark on a linear axis) or axis is unknown.
    """
    if axis not in ("x", "y"):
        raise ValueError("axis must be 'x' or 'y'")
    getter = ax.get_xscale if axis == "x" else ax.get_yscale
    if getter() != "log":
        raise ValueError(f"{axis}-axis is not log-scaled; refusing to draw a spurious zero mark")

    if axis == "x":
        floor = ax.get_xlim()[0]
        ax.axvline(floor, color=color, linestyle="--", linewidth=0.8, zorder=0)
        existing = ax.get_xlabel()
    else:
        floor = ax.get_ylim()[0]
        ax.axhline(floor, color=color, linestyle="--", linewidth=0.8, zorder=0)
        existing = ax.get_ylabel()
    note = f"{label} at axis floor"
    if axis == "x":
        ax.set_xlabel(existing + (f" ({note})" if existing else note))
    else:
        ax.set_ylabel(existing + (f" ({note})" if existing else note))
    logger.debug(f"Drew zero mark on {axis}-axis floor {floor}")


# ---------------------------------------------------------------------------
# Deterministic figure export
# ---------------------------------------------------------------------------

#: rc-params that reduce nondeterministic content in figure output (font
#: embedding choices; per-file timestamps/producer metadata are stripped by
#: :func:`save_figure_deterministic` because matplotlib versions differ on
#: whether ``savefig.metadata`` is a valid rcParam).
DETERMINISTIC_RCPARAMS: Dict[str, object] = {
    "ps.useafm": True,
    "svg.fonttype": "none",
}


def apply_deterministic_rcparams() -> None:
    """Apply rc-params that strip nondeterministic metadata from figure output."""
    mpl.rcParams.update(DETERMINISTIC_RCPARAMS)
    logger.debug("Applied deterministic figure rc-params")


def save_figure_deterministic(fig: "mpl.figure.Figure", path, *, dpi: int = 300) -> None:
    """Save *fig* to *path* with reproducible bytes.

    Passes fixed (empty) metadata for Software/Creator/date fields so that
    identical figure state yields identical files across runs — required for
    regression tests on figure output and for reproducible publications.
    """
    metadata = {
        "Software": "",
        "Creator": "",
        "Creation Time": "",
        "Date": "",
    }
    fig.savefig(path, dpi=dpi, metadata=metadata)
