"""Color palette utilities for biological data visualization.

Provides curated palettes for chromosomes, expression data, significance
levels, and other common bioinformatics use cases. All palettes include
colorblind-safe options by default.
"""

from __future__ import annotations

from typing import Dict, List, Sequence

import numpy as np
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt

from metainformant.core import logging

logger = logging.get_logger(__name__)

# ---------------------------------------------------------------------------
# Chromosome palette (22 autosomes + X, Y, MT)
# ---------------------------------------------------------------------------

CHROMOSOME_COLORS: Dict[str, str] = {
    "1": "#E41A1C",
    "2": "#377EB8",
    "3": "#4DAF4A",
    "4": "#984EA3",
    "5": "#FF7F00",
    "6": "#A65628",
    "7": "#F781BF",
    "8": "#999999",
    "9": "#66C2A5",
    "10": "#FC8D62",
    "11": "#8DA0CB",
    "12": "#E78AC3",
    "13": "#A6D854",
    "14": "#FFD92F",
    "15": "#E5C494",
    "16": "#B3B3B3",
    "17": "#1B9E77",
    "18": "#D95F02",
    "19": "#7570B3",
    "20": "#E7298A",
    "21": "#66A61E",
    "22": "#E6AB02",
    "X": "#A6761D",
    "Y": "#666666",
    "MT": "#000000",
}


def chromosome_palette(chromosomes: Sequence[str] | None = None) -> List[str]:
    """Return hex colors for the given chromosomes.

    Args:
        chromosomes: Chromosome labels (e.g. ``["1", "2", "X"]``).
            If *None*, returns the full 25-color palette.
    """
    if chromosomes is None:
        return list(CHROMOSOME_COLORS.values())
    return [CHROMOSOME_COLORS.get(str(c).replace("chr", ""), "#888888") for c in chromosomes]


# ---------------------------------------------------------------------------
# Colorblind-safe categorical palettes (Wong 2011, Tol, IBM)
# ---------------------------------------------------------------------------

WONG: List[str] = [
    "#000000",  # black
    "#E69F00",  # orange
    "#56B4E9",  # sky blue
    "#009E73",  # bluish green
    "#F0E442",  # yellow
    "#0072B2",  # blue
    "#D55E00",  # vermilion
    "#CC79A7",  # reddish purple
]

TOL_BRIGHT: List[str] = [
    "#4477AA",  # blue
    "#EE6677",  # red
    "#228833",  # green
    "#CCBB44",  # yellow
    "#66CCEE",  # cyan
    "#AA3377",  # purple
    "#BBBBBB",  # grey
]

IBM_COLORBLIND: List[str] = [
    "#648FFF",  # ultramarine
    "#785EF0",  # purple
    "#DC267F",  # magenta
    "#FE6100",  # orange
    "#FFB000",  # gold
]


def categorical(n: int, *, palette: str = "wong") -> List[str]:
    """Return *n* colorblind-safe categorical colors.

    Args:
        n: Number of colors needed.
        palette: ``"wong"`` (default, 8 colors), ``"tol"`` (7 colors), or
            ``"ibm"`` (5 colors). Falls back to matplotlib tab20 for n > palette size.
    """
    palettes = {"wong": WONG, "tol": TOL_BRIGHT, "ibm": IBM_COLORBLIND}
    base = palettes.get(palette, WONG)
    if n <= len(base):
        return base[:n]
    # Fall back to tab20 for larger needs
    cmap = plt.cm.get_cmap("tab20", n)
    return [mcolors.to_hex(cmap(i)) for i in range(n)]


# ---------------------------------------------------------------------------
# Sequential / diverging palettes for biological data
# ---------------------------------------------------------------------------


def expression_gradient(n: int = 256, *, center_white: bool = True) -> mcolors.LinearSegmentedColormap:
    """Blue-white-red gradient for gene expression (log2 fold-change).

    Args:
        n: Number of color stops.
        center_white: Whether white is at the midpoint.
    """
    if center_white:
        colors = ["#2166AC", "#4393C3", "#D1E5F0", "#FFFFFF", "#FDDBC7", "#F4A582", "#D6604D", "#B2182B"]
    else:
        colors = ["#2166AC", "#4393C3", "#92C5DE", "#F4A582", "#D6604D", "#B2182B"]
    return mcolors.LinearSegmentedColormap.from_list("expression", colors, N=n)


def significance_palette() -> Dict[str, str]:
    """Standard colors for p-value significance levels."""
    return {
        "highly_significant": "#D32F2F",  # p < 0.001
        "significant": "#FF9800",  # p < 0.01
        "marginally_significant": "#FFC107",  # p < 0.05
        "not_significant": "#9E9E9E",  # p >= 0.05
    }


def significance_color(pvalue: float) -> str:
    """Return a color string for a given p-value."""
    if pvalue < 0.001:
        return "#D32F2F"
    if pvalue < 0.01:
        return "#FF9800"
    if pvalue < 0.05:
        return "#FFC107"
    return "#9E9E9E"


def heatmap_cmap(name: str = "expression") -> mcolors.Colormap:
    """Return a named colormap suitable for biological heatmaps.

    Args:
        name: ``"expression"`` (blue-white-red), ``"methylation"``
            (green-yellow-red), ``"quality"`` (red-yellow-green).
    """
    cmaps = {
        "expression": mcolors.LinearSegmentedColormap.from_list(
            "expression",
            ["#2166AC", "#FFFFFF", "#B2182B"],
            N=256,
        ),
        "methylation": mcolors.LinearSegmentedColormap.from_list(
            "methylation",
            ["#1A9850", "#FFFFBF", "#D73027"],
            N=256,
        ),
        "quality": mcolors.LinearSegmentedColormap.from_list(
            "quality",
            ["#D73027", "#FFFFBF", "#1A9850"],
            N=256,
        ),
    }
    if name not in cmaps:
        raise KeyError(f"Unknown heatmap colormap '{name}'. Available: {list(cmaps.keys())}")
    return cmaps[name]


def alternating_pair(n: int = 2) -> List[str]:
    """Return alternating color pairs for Manhattan-style plots.

    Args:
        n: Number of alternating color pairs.
    """
    base_pairs = [
        ("#4393C3", "#2166AC"),
        ("#F4A582", "#D6604D"),
        ("#A6D854", "#4DAF4A"),
        ("#FFD92F", "#E6AB02"),
    ]
    result = []
    for i in range(n):
        pair = base_pairs[i % len(base_pairs)]
        result.extend(pair)
    return result
