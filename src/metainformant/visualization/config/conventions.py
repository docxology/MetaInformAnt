"""Shared chart-convention enforcement utilities.

Centralizes the plotting conventions used by the visualization module and by
downstream consumers (e.g. rna/analysis atlas figures, campaign figures):

- Okabe-Ito color-blind-safe palette (single canonical definition).
- Hatch patterns paired with categorical fills so category identity survives
  grayscale printing and color-vision deficiency (hatch + text redundancy).
- Explicit zero marks on logarithmic axes (a log scale silently clamps zero
  to the axis floor; an unmarked floor invites misreading).
- Deterministic figure export: matplotlib PNG/PDF/SVG output embeds timestamps,
  host metadata and (SVG) random element ids by default; deterministic mode
  strips them so byte-identical inputs produce byte-identical files.

All statistics and labels rendered by consumers are DESCRIPTIVE ONLY; nothing
here computes inferential results.
"""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING, Any, Dict, List, Sequence, Tuple, cast

import matplotlib as mpl

from metainformant.core.utils import logging

if TYPE_CHECKING:
    from matplotlib.typing import RcKeyType

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
    mpl.rcParams.update(cast("Dict[RcKeyType, Any]", DETERMINISTIC_RCPARAMS))
    logger.debug("Applied deterministic figure rc-params")


# Metadata that blanks the producer and clock fields of each output format. The
# matplotlib backends do not share a vocabulary and each validates its own keys,
# so one mapping cannot serve them all: PNG takes free-form text chunks; PDF takes
# a fixed info dictionary (``None`` removes a predefined key) and only warns about
# other keys; SVG takes Dublin Core keys (``None`` drops the predefined ``Creator``
# and ``Date``) and raises on any other key; PS/EPS read only ``Creator``. Formats
# without an entry (JPEG, TIFF, WebP, GIF, AVIF, PGF, ...) accept no metadata at
# all, so none is passed for them.
_DETERMINISTIC_METADATA: Dict[str, Dict[str, str | None]] = {
    "png": {"Software": "", "Creator": "", "Creation Time": "", "Date": ""},
    "pdf": {"Creator": "", "Producer": "", "CreationDate": None, "ModDate": None},
    "svg": {"Creator": None, "Date": None},
    "svgz": {"Creator": None, "Date": None},
    "ps": {"Creator": ""},
    "eps": {"Creator": ""},
}

# Value of ``svg.hashsalt`` while an SVG is saved. matplotlib hashes the ids of
# clip paths, markers, glyphs and images with this salt and draws a fresh uuid4
# whenever it is ``None``, so the ids would differ on every save.
_SVG_HASH_SALT = "metainformant-deterministic-svg"


def save_figure_deterministic(
    fig: "mpl.figure.Figure",
    path: "str | Path",
    *,
    dpi: int = 300,
    bbox_inches: "str | None" = None,
) -> None:
    """Save *fig* to *path* with reproducible bytes.

    Blanks the producer and clock fields matplotlib otherwise stamps into the
    file, so identical figure state yields identical files across runs —
    required for regression tests on figure output and for reproducible
    publications. Each format receives only the metadata its backend accepts;
    the format is read from the suffix of *path*, as matplotlib does:

    - PNG: empty ``Software``, ``Creator``, ``Creation Time`` and ``Date``.
    - PDF: empty ``Creator`` and ``Producer``; no ``CreationDate``/``ModDate``.
    - SVG: no ``Creator`` and no ``Date``, and ``svg.hashsalt`` pinned for the
      duration of the save (matplotlib hashes element ids with a fresh uuid4
      otherwise) unless the caller already set one.
    - JPEG, TIFF, WebP, GIF and AVIF carry no metadata and are saved as is.
    - PS/EPS and SVGZ are only partly covered. PS/EPS get an empty ``Creator``,
      but matplotlib stamps ``%%CreationDate`` from the clock (or from
      ``SOURCE_DATE_EPOCH``) and no metadata key removes it, so they are
      reproducible only under a fixed ``SOURCE_DATE_EPOCH``. SVGZ gets the SVG
      treatment, but its gzip header carries the write time, which matplotlib
      does not expose.

    ``bbox_inches="tight"`` is accepted so pyplot-level
    ``savefig(..., bbox_inches="tight")`` call sites can route through this
    helper without changing output geometry.
    """
    fmt = Path(path).suffix.lstrip(".").lower() or fig.canvas.get_default_filetype()
    save_kwargs: Dict[str, Any] = {"dpi": dpi, "bbox_inches": bbox_inches}
    metadata = _DETERMINISTIC_METADATA.get(fmt)
    if metadata is not None:
        save_kwargs["metadata"] = dict(metadata)

    rc: Dict[RcKeyType, Any] = {}
    if fmt in ("svg", "svgz") and mpl.rcParams["svg.hashsalt"] is None:
        rc["svg.hashsalt"] = _SVG_HASH_SALT
    with mpl.rc_context(rc):
        fig.savefig(path, **save_kwargs)
