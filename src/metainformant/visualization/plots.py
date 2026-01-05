"""Enhanced plotting functions for biological data visualization.

This module provides a unified interface for common plotting tasks in biology,
including heatmaps, volcano plots, Manhattan plots, PCA plots, and more.
It re-exports and wraps functions from specialized visualization modules.

Main Functions:
    - correlation_heatmap: Plot correlation matrix as heatmap
    - volcano_plot: Plot log2FC vs p-value
    - expression_heatmap: Plot gene expression heatmap
    - pca_plot: Plot pre-computed PCA scatter
    - scatter_plot: Generic scatter plot
    - manhattan_plot: Plot genome-wide association results
    - qq_plot: Q-Q plot for QC

Example:
    >>> from metainformant.visualization import plots
    >>> import pandas as pd
    >>> data = pd.DataFrame({'A': [1, 2, 3], 'B': [4, 5, 6]})
    >>> ax = plots.correlation_heatmap(data)
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Optional

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.axes import Axes

from metainformant.core import logging, paths

logger = logging.get_logger(__name__)

# Optional imports
try:
    import seaborn as sns

    HAS_SEABORN = True
except ImportError:
    sns = None
    HAS_SEABORN = False


# Store original functions for wrapping
_original_correlation_heatmap = None
_original_qq_plot = None
_original_volcano_plot = None
_original_manhattan_plot = None
_original_scatter_plot = None


def _import_existing_functions() -> None:
    """Safely import functions from existing visualization modules."""
    global _original_correlation_heatmap, _original_qq_plot
    global _original_volcano_plot, _original_manhattan_plot
    global _original_scatter_plot

    try:
        from .statistical import correlation_heatmap as ch, qq_plot as qqp
        _original_correlation_heatmap = ch
        _original_qq_plot = qqp
    except ImportError:
        pass

    try:
        from .genomics import volcano_plot as vp, manhattan_plot as mp
        _original_volcano_plot = vp
        _original_manhattan_plot = mp
    except ImportError:
        pass

    try:
        from .basic import scatter_plot as sp
        _original_scatter_plot = sp
    except ImportError:
        pass


# Import on module load
_import_existing_functions()


def expression_heatmap(
    data: pd.DataFrame,
    *,
    ax: Optional[Axes] = None,
    output_path: Optional[str | Path] = None,
    **kwargs: Any,
) -> Axes:
    """Plot gene expression heatmap with hierarchical clustering.

    Creates a clustered heatmap of gene expression data showing samples
    and genes with color intensity representing expression levels.

    Args:
        data: Gene expression data (genes × samples or samples × genes)
        ax: matplotlib Axes object (creates new if None)
        output_path: Optional path to save figure
        **kwargs: Additional arguments for plotting:
            - cmap: colormap name (default: 'RdBu_r')
            - vmin/vmax: value limits for color scaling
            - figsize: figure size tuple
            - cbar_label: colorbar label

    Returns:
        matplotlib Axes object with the plot

    Raises:
        ValueError: If data is empty

    Example:
        >>> import pandas as pd
        >>> data = pd.DataFrame(
        ...     [[1, 2, 3], [4, 5, 6]],
        ...     index=['Gene1', 'Gene2'],
        ...     columns=['Sample1', 'Sample2', 'Sample3']
        ... )
        >>> ax = expression_heatmap(data)
        >>> assert ax is not None
    """
    if data.empty:
        raise ValueError("Expression data cannot be empty")

    # Create axes if needed
    if ax is None:
        figsize = kwargs.pop("figsize", (12, 8))
        fig, ax = plt.subplots(figsize=figsize)

    # Get plotting parameters
    cmap = kwargs.get("cmap", "RdBu_r")
    vmin = kwargs.get("vmin", None)
    vmax = kwargs.get("vmax", None)

    # Plot heatmap
    if HAS_SEABORN:
        sns.heatmap(
            data,
            ax=ax,
            cmap=cmap,
            vmin=vmin,
            vmax=vmax,
            cbar_kws={"label": kwargs.get("cbar_label", "Expression Level")},
        )
    else:
        # Fallback to basic matplotlib
        im = ax.imshow(data.values, aspect="auto", cmap=cmap, vmin=vmin, vmax=vmax)
        ax.set_xticks(range(len(data.columns)))
        ax.set_yticks(range(len(data.index)))
        ax.set_xticklabels(data.columns, rotation=45, ha="right")
        ax.set_yticklabels(data.index)
        cbar = plt.colorbar(im, ax=ax)
        cbar.set_label(kwargs.get("cbar_label", "Expression Level"))

    ax.set_title("Gene Expression Heatmap")

    # Save if requested
    if output_path:
        output_path = Path(output_path)
        paths.ensure_directory(output_path.parent)
        plt.savefig(output_path, dpi=300, bbox_inches="tight")
        logger.info(f"Expression heatmap saved to {output_path}")

    return ax


def pca_plot(
    data: pd.DataFrame,
    *,
    hue: Optional[str] = None,
    ax: Optional[Axes] = None,
    output_path: Optional[str | Path] = None,
    **kwargs: Any,
) -> Axes:
    """Plot pre-computed PCA scatter with variance information.

    Creates a scatter plot of pre-computed principal components with
    variance explained shown in axis labels.

    Args:
        data: DataFrame with columns:
            - 'PC1': First principal component
            - 'PC2': Second principal component
            - 'PC1_variance': Variance explained by PC1 (percentage)
            - 'PC2_variance': Variance explained by PC2 (percentage)
            - Optional: column specified by `hue` for color groups
        hue: Optional column name for coloring points (must exist in data)
        ax: matplotlib Axes object (creates new if None)
        output_path: Optional path to save figure
        **kwargs: Additional arguments:
            - figsize: figure size tuple (default: (8, 6))
            - alpha: point transparency (default: 0.6)
            - s: point size (default: 50)

    Returns:
        matplotlib Axes object with the plot

    Raises:
        ValueError: If required PC columns are missing or data is empty

    Example:
        >>> import pandas as pd
        >>> data = pd.DataFrame({
        ...     'PC1': [0.1, 0.2, 0.3],
        ...     'PC2': [0.4, 0.5, 0.6],
        ...     'PC1_variance': [25.5, 25.5, 25.5],
        ...     'PC2_variance': [18.3, 18.3, 18.3],
        ...     'cluster': ['A', 'B', 'A']
        ... })
        >>> ax = pca_plot(data, hue='cluster')
    """
    if data.empty:
        raise ValueError("PCA data cannot be empty")

    # Check for required columns
    required_cols = ["PC1", "PC2", "PC1_variance", "PC2_variance"]
    missing_cols = [col for col in required_cols if col not in data.columns]

    if missing_cols:
        raise ValueError(f"Missing required columns for PCA plot: {missing_cols}")

    # Create axes if needed
    if ax is None:
        figsize = kwargs.pop("figsize", (8, 6))
        fig, ax = plt.subplots(figsize=figsize)

    # Get plotting parameters
    alpha = kwargs.get("alpha", 0.6)
    s = kwargs.get("s", 50)
    c = None

    # Handle color by hue
    if hue:
        if hue not in data.columns:
            raise ValueError(f"Hue column '{hue}' not found in data")

        # Create color mapping
        unique_vals = data[hue].unique()
        colors = plt.cm.tab10(np.linspace(0, 1, len(unique_vals)))
        color_map = {val: colors[i] for i, val in enumerate(unique_vals)}
        c = [color_map[val] for val in data[hue]]

    # Plot scatter
    scatter = ax.scatter(
        data["PC1"],
        data["PC2"],
        c=c,
        alpha=alpha,
        s=s,
        edgecolors="black",
        linewidth=0.5,
    )

    # Add variance to labels
    pc1_var = data["PC1_variance"].iloc[0] if len(data) > 0 else 0
    pc2_var = data["PC2_variance"].iloc[0] if len(data) > 0 else 0

    ax.set_xlabel(f"PC1 ({pc1_var:.1f}% variance)")
    ax.set_ylabel(f"PC2 ({pc2_var:.1f}% variance)")
    ax.set_title("PCA Plot")
    ax.grid(True, alpha=0.3)

    # Add legend if using hue
    if hue and c is not None:
        unique_vals = data[hue].unique()
        colors = plt.cm.tab10(np.linspace(0, 1, len(unique_vals)))
        handles = [
            plt.Line2D([0], [0], marker="o", color="w", markerfacecolor=colors[i], markersize=8)
            for i in range(len(unique_vals))
        ]
        ax.legend(handles, unique_vals, title=hue, loc="best")

    # Save if requested
    if output_path:
        output_path = Path(output_path)
        paths.ensure_directory(output_path.parent)
        plt.savefig(output_path, dpi=300, bbox_inches="tight")
        logger.info(f"PCA plot saved to {output_path}")

    return ax


# Wrapper functions for re-exported modules to handle positional arguments
def correlation_heatmap(
    data: pd.DataFrame,
    *,
    ax: Optional[Axes] = None,
    output_path: Optional[str | Path] = None,
    **kwargs: Any,
) -> Axes:
    """Plot correlation matrix as heatmap.

    Wrapper around statistical.correlation_heatmap that provides consistent API.
    """
    if _original_correlation_heatmap is None:
        raise ImportError("correlation_heatmap not available in statistical module")

    # Remove annot from kwargs to avoid duplicate keyword argument error
    # Pass it explicitly if present
    annot = kwargs.pop("annot", True)

    return _original_correlation_heatmap(
        data, ax=ax, output_path=output_path, annot=annot, **kwargs
    )


def qq_plot(
    data: np.ndarray,
    *,
    ax: Optional[Axes] = None,
    output_path: Optional[str | Path] = None,
    **kwargs: Any,
) -> Axes:
    """Create Q-Q plot for QC analysis.

    Wrapper around statistical.qq_plot that provides consistent API.
    """
    if _original_qq_plot is None:
        raise ImportError("qq_plot not available in statistical module")

    return _original_qq_plot(data, ax=ax, output_path=output_path, **kwargs)


def volcano_plot(
    data: pd.DataFrame,
    log2fc_col: str = "log2FoldChange",
    pval_col: str = "padj",
    *,
    ax: Optional[Axes] = None,
    output_path: Optional[str | Path] = None,
    **kwargs: Any,
) -> Axes:
    """Plot log2 fold change vs p-value for differential expression.

    Wrapper around genomics.volcano_plot that supports both positional
    and keyword arguments for column names.
    """
    if _original_volcano_plot is None:
        raise ImportError("volcano_plot not available in genomics module")

    return _original_volcano_plot(
        data, log2fc_col=log2fc_col, pval_col=pval_col,
        ax=ax, output_path=output_path, **kwargs
    )


def manhattan_plot(
    data: pd.DataFrame,
    pos_col: str = "BP",
    pval_col: str = "P",
    *,
    ax: Optional[Axes] = None,
    output_path: Optional[str | Path] = None,
    **kwargs: Any,
) -> Axes:
    """Plot genome-wide association results.

    Wrapper around genomics.manhattan_plot that supports both positional
    and keyword arguments for column names.
    """
    if _original_manhattan_plot is None:
        raise ImportError("manhattan_plot not available in genomics module")

    # Auto-detect chromosome column if not provided
    chr_col = kwargs.pop("chr_col", None)
    if chr_col is None:
        # Look for common chromosome column names
        for col_name in ["chromosome", "chr", "CHR", "Chromosome"]:
            if col_name in data.columns:
                chr_col = col_name
                break
        else:
            # Default fallback
            chr_col = "CHR"

    return _original_manhattan_plot(
        data, chr_col=chr_col, pos_col=pos_col, pval_col=pval_col,
        ax=ax, output_path=output_path, **kwargs
    )


def scatter_plot(
    x: np.ndarray,
    y: np.ndarray,
    *,
    ax: Optional[Axes] = None,
    output_path: Optional[str | Path] = None,
    **kwargs: Any,
) -> Axes:
    """Generic scatter plot.

    Wrapper around basic.scatter_plot that provides consistent API.
    """
    if _original_scatter_plot is None:
        raise ImportError("scatter_plot not available in basic module")

    return _original_scatter_plot(x, y, ax=ax, output_path=output_path, **kwargs)


__all__ = [
    "correlation_heatmap",
    "volcano_plot",
    "manhattan_plot",
    "qq_plot",
    "scatter_plot",
    "expression_heatmap",
    "pca_plot",
]
