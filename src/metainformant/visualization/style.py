"""Styling utilities for publication-quality figures.

This module provides style presets, color palettes, font management,
and figure size presets for creating publication-quality visualizations.
"""

from __future__ import annotations

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.style as mpl_style

# Use non-interactive backend by default for tests/headless
matplotlib.use("Agg", force=True)


# Publication-quality style presets
PUBLICATION_STYLE = {
    'figure.figsize': (8, 6),
    'figure.dpi': 300,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
    'font.size': 10,
    'axes.labelsize': 11,
    'axes.titlesize': 12,
    'xtick.labelsize': 9,
    'ytick.labelsize': 9,
    'legend.fontsize': 9,
    'figure.titlesize': 14,
    'axes.linewidth': 1.0,
    'grid.linewidth': 0.5,
    'lines.linewidth': 1.5,
    'patch.linewidth': 0.5,
    'xtick.major.width': 1.0,
    'ytick.major.width': 1.0,
    'xtick.minor.width': 0.5,
    'ytick.minor.width': 0.5,
    'axes.spines.left': True,
    'axes.spines.bottom': True,
    'axes.spines.top': False,
    'axes.spines.right': False,
}

# Scientific color palettes
SCIENTIFIC_COLORS = {
    'primary': ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b'],
    'accessible': ['#006BA4', '#FF800E', '#ABABAB', '#595959', '#5F9ED1', '#C85200'],
    'colorblind': ['#0173b2', '#de8f05', '#029e73', '#cc78bc', '#56b4e9', '#ece133'],
    'viridis_like': ['#440154', '#31688e', '#35b779', '#6ece58', '#fde725'],
}

# Default figure sizes
FIGURE_SIZES = {
    'small': (4, 3),
    'medium': (6, 4.5),
    'large': (8, 6),
    'wide': (10, 4),
    'tall': (6, 8),
    'square': (6, 6),
    'presentation': (10, 7.5),
    'poster': (12, 9),
}


def apply_publication_style(**overrides) -> None:
    """Apply publication-quality style settings.

    Args:
        **overrides: Style parameter overrides

    Example:
        >>> from metainformant.visualization.style import apply_publication_style
        >>> apply_publication_style(font.size=12)
    """
    style = PUBLICATION_STYLE.copy()
    style.update(overrides)
    plt.rcParams.update(style)


def get_color_palette(name: str = 'primary') -> list[str]:
    """Get a color palette by name.

    Args:
        name: Palette name ('primary', 'accessible', 'colorblind', 'viridis_like')

    Returns:
        List of color hex codes

    Example:
        >>> from metainformant.visualization.style import get_color_palette
        >>> colors = get_color_palette('colorblind')
    """
    return SCIENTIFIC_COLORS.get(name, SCIENTIFIC_COLORS['primary']).copy()


def get_figure_size(size: str = 'medium') -> tuple[float, float]:
    """Get a figure size by name.

    Args:
        size: Size name (see FIGURE_SIZES)

    Returns:
        Figure size tuple (width, height)

    Example:
        >>> from metainformant.visualization.style import get_figure_size
        >>> figsize = get_figure_size('large')
    """
    return FIGURE_SIZES.get(size, FIGURE_SIZES['medium'])


def set_font_family(family: str = 'sans-serif') -> None:
    """Set the default font family.

    Args:
        family: Font family ('serif', 'sans-serif', 'monospace')

    Example:
        >>> from metainformant.visualization.style import set_font_family
        >>> set_font_family('serif')
    """
    plt.rcParams['font.family'] = family


def set_font_size(size: float = 10) -> None:
    """Set the default font size.

    Args:
        size: Font size in points

    Example:
        >>> from metainformant.visualization.style import set_font_size
        >>> set_font_size(12)
    """
    plt.rcParams['font.size'] = size


def reset_style() -> None:
    """Reset matplotlib style to default.

    Example:
        >>> from metainformant.visualization.style import reset_style
        >>> reset_style()
    """
    plt.rcParams.update(plt.rcParamsDefault)


def apply_style(style_name: str) -> None:
    """Apply a named matplotlib style.

    Args:
        style_name: Style name (e.g., 'seaborn', 'ggplot', 'classic')

    Example:
        >>> from metainformant.visualization.style import apply_style
        >>> apply_style('seaborn-v0_8')
    """
    try:
        mpl_style.use(style_name)
    except OSError:
        # Fallback to default if style not found
        pass

