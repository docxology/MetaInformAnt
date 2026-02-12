"""GWAS visualization configuration and theming.

This module provides configurable plot styles, predefined themes, and utilities
for applying consistent visual settings across all GWAS visualization outputs.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, Union

from metainformant.core.utils import logging

logger = logging.get_logger(__name__)


@dataclass
class PlotStyle:
    """Configuration for GWAS plot visual style.

    Controls figure dimensions, font sizes, colors, and rendering options
    used across Manhattan plots, Q-Q plots, and regional association plots.

    Attributes:
        figsize: Figure width and height in inches.
        dpi: Resolution in dots per inch.
        font_size: Base font size for labels and ticks.
        title_size: Font size for plot titles.
        colormap: Matplotlib colormap name for chromosome coloring.
        significance_color: Color for genome-wide significance line.
        suggestive_color: Color for suggestive significance line.
        point_size: Marker size for scatter plots.
        alpha: Transparency for plot elements (0.0 to 1.0).
        theme: Theme name this style derives from.
        font_family: Font family for all text elements.
        grid: Whether to display grid lines.
        grid_alpha: Transparency for grid lines.
    """

    figsize: Tuple[float, float] = (12, 6)
    dpi: int = 300
    font_size: int = 12
    title_size: int = 14
    colormap: str = "viridis"
    significance_color: str = "red"
    suggestive_color: str = "blue"
    point_size: float = 8.0
    alpha: float = 0.7
    theme: str = "publication"
    font_family: str = "sans-serif"
    grid: bool = True
    grid_alpha: float = 0.3


THEMES: Dict[str, PlotStyle] = {
    "publication": PlotStyle(
        figsize=(10, 6),
        dpi=300,
        font_size=10,
        title_size=12,
    ),
    "presentation": PlotStyle(
        figsize=(14, 8),
        dpi=150,
        font_size=14,
        title_size=18,
    ),
    "poster": PlotStyle(
        figsize=(20, 12),
        dpi=300,
        font_size=18,
        title_size=24,
        point_size=12.0,
    ),
}


def get_style(theme: str = "publication", **overrides: Any) -> PlotStyle:
    """Get a PlotStyle for the given theme with optional overrides.

    If the requested theme is not found, falls back to ``"publication"``
    and emits a warning.

    Args:
        theme: Theme name (``"publication"``, ``"presentation"``, or ``"poster"``).
        **overrides: Keyword arguments to override individual style fields.

    Returns:
        A :class:`PlotStyle` instance with the theme defaults and any overrides applied.
    """
    if theme not in THEMES:
        logger.warning(
            "Unknown theme '%s', falling back to 'publication'. " "Available themes: %s",
            theme,
            list(THEMES.keys()),
        )
        theme = "publication"

    # Start from a copy of the theme's style
    base = THEMES[theme]
    style_kwargs: Dict[str, Any] = {f.name: getattr(base, f.name) for f in base.__dataclass_fields__.values()}

    # Apply caller overrides
    for key, value in overrides.items():
        if key in style_kwargs:
            style_kwargs[key] = value
        else:
            logger.warning("Ignoring unknown PlotStyle field: '%s'", key)

    return PlotStyle(**style_kwargs)


def apply_style(style: PlotStyle) -> None:
    """Apply a :class:`PlotStyle` to matplotlib's global rcParams.

    Sets ``font.size``, ``axes.titlesize``, ``figure.dpi``, ``figure.figsize``,
    ``font.family``, ``axes.grid``, and ``grid.alpha``.

    Args:
        style: The PlotStyle to apply.
    """
    try:
        import matplotlib as mpl

        mpl.rcParams["font.size"] = style.font_size
        mpl.rcParams["axes.titlesize"] = style.title_size
        mpl.rcParams["figure.dpi"] = style.dpi
        mpl.rcParams["figure.figsize"] = list(style.figsize)
        mpl.rcParams["font.family"] = style.font_family
        mpl.rcParams["axes.grid"] = style.grid
        mpl.rcParams["grid.alpha"] = style.grid_alpha
        logger.info("Applied PlotStyle (theme=%s) to matplotlib rcParams", style.theme)
    except ImportError:
        logger.warning("matplotlib is not installed; cannot apply PlotStyle to rcParams")


def style_from_config(config: Dict[str, Any]) -> PlotStyle:
    """Create a :class:`PlotStyle` from a nested configuration dictionary.

    Extracts values from ``config["output"]["visualization"]`` (matching the
    YAML config structure used by the GWAS workflow). Missing keys fall back
    to :class:`PlotStyle` defaults.

    Args:
        config: Configuration dictionary, typically loaded from YAML.

    Returns:
        A PlotStyle populated from the config values.
    """
    viz_config = config.get("output", {}).get("visualization", {})

    # Map config keys to PlotStyle fields with defaults from the dataclass
    defaults = PlotStyle()
    return PlotStyle(
        figsize=tuple(viz_config.get("figsize", defaults.figsize)),  # type: ignore[arg-type]
        dpi=viz_config.get("dpi", defaults.dpi),
        font_size=viz_config.get("font_size", defaults.font_size),
        title_size=viz_config.get("title_size", defaults.title_size),
        colormap=viz_config.get("colormap", defaults.colormap),
        significance_color=viz_config.get("significance_color", defaults.significance_color),
        suggestive_color=viz_config.get("suggestive_color", defaults.suggestive_color),
        point_size=viz_config.get("point_size", defaults.point_size),
        alpha=viz_config.get("alpha", defaults.alpha),
        theme=viz_config.get("theme", defaults.theme),
        font_family=viz_config.get("font_family", defaults.font_family),
        grid=viz_config.get("grid", defaults.grid),
        grid_alpha=viz_config.get("grid_alpha", defaults.grid_alpha),
    )
