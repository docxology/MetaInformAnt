"""Export utilities for high-resolution figure saving.

This module provides utilities for exporting figures in multiple formats
with high-resolution settings and batch export capabilities.
"""

from __future__ import annotations

from pathlib import Path
from typing import Sequence

import matplotlib.pyplot as plt


def save_figure(
    fig: plt.Figure,
    path: str | Path,
    *,
    dpi: int = 300,
    bbox_inches: str = 'tight',
    format: str | None = None,
    **kwargs
) -> Path:
    """Save a figure with high-resolution settings.

    Args:
        fig: Matplotlib figure
        path: Output path
        dpi: Resolution in dots per inch
        bbox_inches: Bounding box mode ('tight', 'standard', etc.)
        format: Output format (if None, inferred from extension)
        **kwargs: Additional arguments for savefig

    Returns:
        Path to saved file

    Example:
        >>> from metainformant.visualization.export import save_figure
        >>> import matplotlib.pyplot as plt
        >>> fig, ax = plt.subplots()
        >>> ax.plot([1, 2, 3], [4, 5, 6])
        >>> path = save_figure(fig, 'output/plot.png', dpi=300)
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)

    # Determine format from extension if not provided
    if format is None:
        format = path.suffix[1:] if path.suffix else 'png'

    fig.savefig(
        path,
        dpi=dpi,
        bbox_inches=bbox_inches,
        format=format,
        **kwargs
    )

    return path


def save_figure_multiformat(
    fig: plt.Figure,
    base_path: str | Path,
    *,
    formats: Sequence[str] = ('png', 'pdf', 'svg'),
    dpi: int = 300,
    **kwargs
) -> list[Path]:
    """Save a figure in multiple formats.

    Args:
        fig: Matplotlib figure
        base_path: Base output path (without extension)
        formats: List of formats to save
        dpi: Resolution for raster formats
        **kwargs: Additional arguments for savefig

    Returns:
        List of saved file paths

    Example:
        >>> from metainformant.visualization.export import save_figure_multiformat
        >>> import matplotlib.pyplot as plt
        >>> fig, ax = plt.subplots()
        >>> paths = save_figure_multiformat(fig, 'output/plot', formats=['png', 'pdf'])
    """
    base_path = Path(base_path)
    base_path.parent.mkdir(parents=True, exist_ok=True)

    saved_paths = []
    for fmt in formats:
        path = base_path.with_suffix(f'.{fmt}')
        save_figure(fig, path, dpi=dpi if fmt in ('png', 'jpg', 'jpeg', 'tiff') else None,
                   format=fmt, **kwargs)
        saved_paths.append(path)

    return saved_paths


def batch_export_figures(
    figures: dict[str, plt.Figure],
    output_dir: str | Path,
    *,
    formats: Sequence[str] = ('png',),
    dpi: int = 300,
    **kwargs
) -> dict[str, list[Path]]:
    """Batch export multiple figures.

    Args:
        figures: Dictionary mapping names to figures
        output_dir: Output directory
        formats: List of formats to save
        dpi: Resolution for raster formats
        **kwargs: Additional arguments for savefig

    Returns:
        Dictionary mapping names to lists of saved paths

    Example:
        >>> from metainformant.visualization.export import batch_export_figures
        >>> import matplotlib.pyplot as plt
        >>> figs = {
        ...     'plot1': plt.figure(),
        ...     'plot2': plt.figure()
        ... }
        >>> paths = batch_export_figures(figs, 'output/plots')
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    saved_paths = {}
    for name, fig in figures.items():
        base_path = output_dir / name
        saved_paths[name] = save_figure_multiformat(fig, base_path, formats=formats,
                                                    dpi=dpi, **kwargs)

    return saved_paths


def get_supported_formats() -> list[str]:
    """Get list of supported export formats.

    Returns:
        List of format names

    Example:
        >>> from metainformant.visualization.export import get_supported_formats
        >>> formats = get_supported_formats()
    """
    return ['png', 'pdf', 'svg', 'eps', 'jpg', 'jpeg', 'tiff', 'ps', 'pgf']

