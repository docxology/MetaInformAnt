"""Theme and styling configuration for visualization.

Provides predefined themes for consistent, publication-ready figures across
the entire visualization module. Themes can be applied globally or via
a context manager for scoped styling.
"""

from __future__ import annotations

import contextlib
from typing import Any, Dict, Generator

import matplotlib as mpl
import matplotlib.pyplot as plt

from metainformant.core.utils import logging

logger = logging.get_logger(__name__)

# ---------------------------------------------------------------------------
# Theme definitions
# ---------------------------------------------------------------------------

_THEMES: Dict[str, Dict[str, Any]] = {
    "publication": {
        "figure.figsize": (7, 5),
        "figure.dpi": 300,
        "font.size": 10,
        "font.family": "sans-serif",
        "axes.titlesize": 12,
        "axes.labelsize": 10,
        "xtick.labelsize": 9,
        "ytick.labelsize": 9,
        "legend.fontsize": 9,
        "axes.linewidth": 1.0,
        "axes.grid": True,
        "grid.alpha": 0.3,
        "grid.linewidth": 0.5,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "savefig.dpi": 300,
        "savefig.bbox": "tight",
        "savefig.pad_inches": 0.1,
    },
    "presentation": {
        "figure.figsize": (10, 7),
        "figure.dpi": 150,
        "font.size": 14,
        "font.family": "sans-serif",
        "axes.titlesize": 18,
        "axes.labelsize": 14,
        "xtick.labelsize": 12,
        "ytick.labelsize": 12,
        "legend.fontsize": 12,
        "axes.linewidth": 1.5,
        "axes.grid": True,
        "grid.alpha": 0.25,
        "grid.linewidth": 0.8,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "lines.linewidth": 2.5,
        "lines.markersize": 8,
    },
    "poster": {
        "figure.figsize": (14, 10),
        "figure.dpi": 150,
        "font.size": 18,
        "font.family": "sans-serif",
        "axes.titlesize": 24,
        "axes.labelsize": 20,
        "xtick.labelsize": 16,
        "ytick.labelsize": 16,
        "legend.fontsize": 16,
        "axes.linewidth": 2.0,
        "axes.grid": True,
        "grid.alpha": 0.2,
        "grid.linewidth": 1.0,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "lines.linewidth": 3.0,
        "lines.markersize": 10,
    },
    "dark": {
        "figure.figsize": (8, 6),
        "figure.dpi": 150,
        "figure.facecolor": "#1e1e1e",
        "axes.facecolor": "#2d2d2d",
        "axes.edgecolor": "#555555",
        "axes.labelcolor": "#cccccc",
        "text.color": "#cccccc",
        "xtick.color": "#cccccc",
        "ytick.color": "#cccccc",
        "axes.titlesize": 12,
        "axes.labelsize": 10,
        "xtick.labelsize": 9,
        "ytick.labelsize": 9,
        "axes.grid": True,
        "grid.color": "#444444",
        "grid.alpha": 0.5,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "legend.facecolor": "#2d2d2d",
        "legend.edgecolor": "#555555",
        "savefig.facecolor": "#1e1e1e",
    },
    "minimal": {
        "figure.figsize": (7, 5),
        "figure.dpi": 150,
        "font.size": 10,
        "font.family": "sans-serif",
        "axes.titlesize": 11,
        "axes.labelsize": 10,
        "xtick.labelsize": 9,
        "ytick.labelsize": 9,
        "axes.linewidth": 0.8,
        "axes.grid": False,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.spines.left": True,
        "axes.spines.bottom": True,
    },
    "colorblind": {
        "figure.figsize": (7, 5),
        "figure.dpi": 150,
        "font.size": 10,
        "axes.prop_cycle": mpl.cycler(
            color=[
                "#0072B2",  # blue
                "#E69F00",  # orange
                "#009E73",  # green
                "#CC79A7",  # pink
                "#56B4E9",  # sky blue
                "#D55E00",  # vermilion
                "#F0E442",  # yellow
                "#000000",  # black
            ]
        ),
        "axes.grid": True,
        "grid.alpha": 0.3,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "lines.linewidth": 2.0,
    },
}


def list_themes() -> list[str]:
    """Return names of all available themes."""
    return list(_THEMES.keys())


def get_theme(name: str) -> Dict[str, Any]:
    """Return the rc-params dictionary for *name*.

    Raises:
        KeyError: If the theme name is not registered.
    """
    if name not in _THEMES:
        raise KeyError(f"Unknown theme '{name}'. Available: {list_themes()}")
    return dict(_THEMES[name])


def register_theme(name: str, params: Dict[str, Any]) -> None:
    """Register a custom theme.

    Args:
        name: Theme name.
        params: Dictionary of matplotlib rc-params.
    """
    _THEMES[name] = dict(params)
    logger.info(f"Registered custom theme '{name}'")


def apply_theme(name: str) -> None:
    """Apply a theme globally (persists until reset or another theme is applied).

    Args:
        name: Theme name to apply.
    """
    params = get_theme(name)
    mpl.rcParams.update(params)
    logger.info(f"Applied theme '{name}' globally")


def reset_theme() -> None:
    """Reset matplotlib rc-params to defaults."""
    mpl.rcdefaults()
    logger.info("Reset to default matplotlib theme")


@contextlib.contextmanager
def theme(name: str) -> Generator[None, None, None]:
    """Context manager that temporarily applies a theme.

    Usage::

        with theme("publication"):
            fig, ax = plt.subplots()
            ax.plot(x, y)
            fig.savefig("out.png")

    After the block, the previous rc-params are restored.
    """
    params = get_theme(name)
    with mpl.rc_context(params):
        logger.debug(f"Entering theme context: {name}")
        yield
    logger.debug(f"Exiting theme context: {name}")
