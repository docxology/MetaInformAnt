"""Dimensionality reduction visualization functions.

This module provides visualization functions for dimensionality reduction
techniques including PCA, UMAP, t-SNE, and related diagnostic plots.
"""

from __future__ import annotations

from typing import Sequence

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# Use non-interactive backend by default for tests/headless
matplotlib.use("Agg", force=True)

try:
    import seaborn as sns

    SEABORN_AVAILABLE = True
except ImportError:
    SEABORN_AVAILABLE = False
    sns = None


def pca_plot(
    data: pd.DataFrame,
    *,
    pc_x: int = 1,
    pc_y: int = 2,
    hue: str | None = None,
    ax: plt.Axes | None = None,
    **kwargs
) -> plt.Axes:
    """Create a PCA scatter plot.

    Args:
        data: DataFrame with PCA coordinates (PC1, PC2, etc.)
        pc_x: Principal component for x-axis
        pc_y: Principal component for y-axis
        hue: Column to color points by
        ax: Matplotlib axes (creates new if None)
        **kwargs: Additional arguments for scatter

    Returns:
        Matplotlib axes object

    Example:
        >>> from metainformant.visualization import pca_plot
        >>> import pandas as pd
        >>> import numpy as np
        >>> data = pd.DataFrame({
        ...     'PC1': np.random.normal(0, 1, 100),
        ...     'PC2': np.random.normal(0, 1, 100),
        ...     'group': ['A'] * 50 + ['B'] * 50
        ... })
        >>> ax = pca_plot(data, hue='group')
    """
    if ax is None:
        _, ax = plt.subplots()

    pc_x_col = f"PC{pc_x}"
    pc_y_col = f"PC{pc_y}"

    if pc_x_col not in data.columns or pc_y_col not in data.columns:
        raise ValueError(f"Columns {pc_x_col} and {pc_y_col} not found in data")

    if hue and hue in data.columns:
        # Color by hue column
        if SEABORN_AVAILABLE:
            sns.scatterplot(data=data, x=pc_x_col, y=pc_y_col, hue=hue, ax=ax, **kwargs)
        else:
            groups = data[hue].unique()
            colors = plt.cm.tab10(np.linspace(0, 1, len(groups)))
            for group, color in zip(groups, colors):
                group_data = data[data[hue] == group]
                ax.scatter(group_data[pc_x_col], group_data[pc_y_col],
                          label=group, c=[color], **kwargs)
            ax.legend()
    else:
        ax.scatter(data[pc_x_col], data[pc_y_col], **kwargs)

    variance_x = data.get(f'{pc_x_col}_variance', 'unknown')
    variance_y = data.get(f'{pc_y_col}_variance', 'unknown')

    ax.set_xlabel(f"{pc_x_col} ({variance_x}%)")
    ax.set_ylabel(f"{pc_y_col} ({variance_y}%)")
    ax.set_title(f"PCA Plot ({pc_x_col} vs {pc_y_col})")

    return ax


def umap_plot(
    data: pd.DataFrame | np.ndarray,
    *,
    x_col: str | int = 0,
    y_col: str | int = 1,
    hue: str | Sequence | None = None,
    ax: plt.Axes | None = None,
    title: str = "UMAP Plot",
    **kwargs
) -> plt.Axes:
    """Create a UMAP visualization plot.

    Args:
        data: DataFrame or array with UMAP coordinates
        x_col: Column name or index for x-axis
        y_col: Column name or index for y-axis
        hue: Column name or array to color points by
        ax: Matplotlib axes (creates new if None)
        title: Plot title
        **kwargs: Additional arguments for scatter

    Returns:
        Matplotlib axes object

    Example:
        >>> from metainformant.visualization import umap_plot
        >>> import pandas as pd
        >>> import numpy as np
        >>> data = pd.DataFrame({
        ...     'UMAP1': np.random.normal(0, 1, 100),
        ...     'UMAP2': np.random.normal(0, 1, 100),
        ...     'group': ['A'] * 50 + ['B'] * 50
        ... })
        >>> ax = umap_plot(data, x_col='UMAP1', y_col='UMAP2', hue='group')
    """
    if ax is None:
        _, ax = plt.subplots()

    if isinstance(data, pd.DataFrame):
        x_data = data[x_col] if isinstance(x_col, str) else data.iloc[:, x_col]
        y_data = data[y_col] if isinstance(y_col, str) else data.iloc[:, y_col]
        hue_data = data[hue] if hue and isinstance(hue, str) else hue
    else:
        data = np.array(data)
        x_data = data[:, x_col]
        y_data = data[:, y_col]
        hue_data = hue

    if hue_data is not None:
        if SEABORN_AVAILABLE and isinstance(hue_data, pd.Series):
            sns.scatterplot(x=x_data, y=y_data, hue=hue_data, ax=ax, **kwargs)
        else:
            # Manual coloring
            if isinstance(hue_data, (list, np.ndarray)):
                unique_values = list(set(hue_data))
                colors = plt.cm.tab10(np.linspace(0, 1, len(unique_values)))
                for val, color in zip(unique_values, colors):
                    mask = np.array(hue_data) == val
                    ax.scatter(x_data[mask], y_data[mask], label=str(val),
                              c=[color], **kwargs)
                ax.legend()
            else:
                ax.scatter(x_data, y_data, c=hue_data, **kwargs)
    else:
        ax.scatter(x_data, y_data, **kwargs)

    ax.set_xlabel("UMAP 1")
    ax.set_ylabel("UMAP 2")
    ax.set_title(title)

    return ax


def tsne_plot(
    data: pd.DataFrame | np.ndarray,
    *,
    x_col: str | int = 0,
    y_col: str | int = 1,
    hue: str | Sequence | None = None,
    ax: plt.Axes | None = None,
    title: str = "t-SNE Plot",
    **kwargs
) -> plt.Axes:
    """Create a t-SNE visualization plot.

    Args:
        data: DataFrame or array with t-SNE coordinates
        x_col: Column name or index for x-axis
        y_col: Column name or index for y-axis
        hue: Column name or array to color points by
        ax: Matplotlib axes (creates new if None)
        title: Plot title
        **kwargs: Additional arguments for scatter

    Returns:
        Matplotlib axes object

    Example:
        >>> from metainformant.visualization import tsne_plot
        >>> import pandas as pd
        >>> import numpy as np
        >>> data = pd.DataFrame({
        ...     'tSNE1': np.random.normal(0, 1, 100),
        ...     'tSNE2': np.random.normal(0, 1, 100),
        ...     'group': ['A'] * 50 + ['B'] * 50
        ... })
        >>> ax = tsne_plot(data, x_col='tSNE1', y_col='tSNE2', hue='group')
    """
    if ax is None:
        _, ax = plt.subplots()

    if isinstance(data, pd.DataFrame):
        x_data = data[x_col] if isinstance(x_col, str) else data.iloc[:, x_col]
        y_data = data[y_col] if isinstance(y_col, str) else data.iloc[:, y_col]
        hue_data = data[hue] if hue and isinstance(hue, str) else hue
    else:
        data = np.array(data)
        x_data = data[:, x_col]
        y_data = data[:, y_col]
        hue_data = hue

    if hue_data is not None:
        if SEABORN_AVAILABLE and isinstance(hue_data, pd.Series):
            sns.scatterplot(x=x_data, y=y_data, hue=hue_data, ax=ax, **kwargs)
        else:
            # Manual coloring
            if isinstance(hue_data, (list, np.ndarray)):
                unique_values = list(set(hue_data))
                colors = plt.cm.tab10(np.linspace(0, 1, len(unique_values)))
                for val, color in zip(unique_values, colors):
                    mask = np.array(hue_data) == val
                    ax.scatter(x_data[mask], y_data[mask], label=str(val),
                              c=[color], **kwargs)
                ax.legend()
            else:
                ax.scatter(x_data, y_data, c=hue_data, **kwargs)
    else:
        ax.scatter(x_data, y_data, **kwargs)

    ax.set_xlabel("t-SNE 1")
    ax.set_ylabel("t-SNE 2")
    ax.set_title(title)

    return ax


def pca_scree_plot(
    explained_variance: Sequence[float],
    *,
    n_components: int | None = None,
    ax: plt.Axes | None = None,
    title: str = "PCA Scree Plot",
    **kwargs
) -> plt.Axes:
    """Create a PCA scree plot showing variance explained.

    Args:
        explained_variance: Variance explained by each component
        n_components: Number of components to show (if None, shows all)
        ax: Matplotlib axes (creates new if None)
        title: Plot title
        **kwargs: Additional arguments

    Returns:
        Matplotlib axes object

    Example:
        >>> from metainformant.visualization import pca_scree_plot
        >>> import numpy as np
        >>> variance = np.array([0.4, 0.3, 0.2, 0.1])
        >>> ax = pca_scree_plot(variance)
    """
    if ax is None:
        _, ax = plt.subplots()

    explained_variance = np.array(explained_variance)
    if n_components:
        explained_variance = explained_variance[:n_components]

    components = range(1, len(explained_variance) + 1)

    ax.plot(components, explained_variance, 'o-', **kwargs)
    ax.set_xlabel("Principal Component")
    ax.set_ylabel("Variance Explained")
    ax.set_title(title)
    ax.grid(True, alpha=0.3)

    return ax


def pca_loadings_plot(
    loadings: pd.DataFrame | np.ndarray,
    pc_x: int = 0,
    pc_y: int = 1,
    *,
    feature_names: Sequence[str] | None = None,
    ax: plt.Axes | None = None,
    title: str = "PCA Loadings Plot",
    **kwargs
) -> plt.Axes:
    """Create a PCA loadings plot.

    Args:
        loadings: Loadings matrix (features x components)
        pc_x: Principal component for x-axis
        pc_y: Principal component for y-axis
        feature_names: Optional feature names
        ax: Matplotlib axes (creates new if None)
        title: Plot title
        **kwargs: Additional arguments

    Returns:
        Matplotlib axes object

    Example:
        >>> from metainformant.visualization import pca_loadings_plot
        >>> import numpy as np
        >>> loadings = np.random.random((10, 5))
        >>> ax = pca_loadings_plot(loadings, feature_names=[f'feature{i}' for i in range(10)])
    """
    if ax is None:
        _, ax = plt.subplots()

    if isinstance(loadings, pd.DataFrame):
        x_data = loadings.iloc[:, pc_x]
        y_data = loadings.iloc[:, pc_y]
        if feature_names is None:
            feature_names = loadings.index.tolist()
    else:
        loadings = np.array(loadings)
        x_data = loadings[:, pc_x]
        y_data = loadings[:, pc_y]

    ax.scatter(x_data, y_data, **kwargs)

    # Add feature names if provided
    if feature_names:
        for i, name in enumerate(feature_names):
            ax.annotate(name, (x_data[i], y_data[i]), fontsize=8, alpha=0.7)

    # Add axes lines
    ax.axhline(y=0, color='k', linestyle='--', alpha=0.5)
    ax.axvline(x=0, color='k', linestyle='--', alpha=0.5)

    ax.set_xlabel(f"PC{pc_x + 1}")
    ax.set_ylabel(f"PC{pc_y + 1}")
    ax.set_title(title)

    return ax


def biplot(
    scores: pd.DataFrame | np.ndarray,
    loadings: pd.DataFrame | np.ndarray,
    *,
    pc_x: int = 0,
    pc_y: int = 1,
    feature_names: Sequence[str] | None = None,
    sample_names: Sequence[str] | None = None,
    ax: plt.Axes | None = None,
    title: str = "PCA Biplot",
    **kwargs
) -> plt.Axes:
    """Create a PCA biplot showing both samples and loadings.

    Args:
        scores: PCA scores (samples x components)
        loadings: PCA loadings (features x components)
        pc_x: Principal component for x-axis
        pc_y: Principal component for y-axis
        feature_names: Optional feature names
        sample_names: Optional sample names
        ax: Matplotlib axes (creates new if None)
        title: Plot title
        **kwargs: Additional arguments

    Returns:
        Matplotlib axes object

    Example:
        >>> from metainformant.visualization import biplot
        >>> import numpy as np
        >>> scores = np.random.normal(0, 1, (20, 5))
        >>> loadings = np.random.random((10, 5))
        >>> ax = biplot(scores, loadings)
    """
    if ax is None:
        _, ax = plt.subplots()

    # Extract coordinates
    if isinstance(scores, pd.DataFrame):
        x_scores = scores.iloc[:, pc_x]
        y_scores = scores.iloc[:, pc_y]
    else:
        scores = np.array(scores)
        x_scores = scores[:, pc_x]
        y_scores = scores[:, pc_y]

    if isinstance(loadings, pd.DataFrame):
        x_loadings = loadings.iloc[:, pc_x]
        y_loadings = loadings.iloc[:, pc_y]
        if feature_names is None:
            feature_names = loadings.index.tolist()
    else:
        loadings = np.array(loadings)
        x_loadings = loadings[:, pc_x]
        y_loadings = loadings[:, pc_y]

    # Plot scores
    ax.scatter(x_scores, y_scores, alpha=0.6, label='Samples', **kwargs)

    # Plot loadings as arrows
    scale = max(abs(x_scores).max(), abs(y_scores).max()) / max(abs(x_loadings).max(), abs(y_loadings).max())
    x_loadings_scaled = x_loadings * scale
    y_loadings_scaled = y_loadings * scale

    for i, (x_ld, y_ld) in enumerate(zip(x_loadings_scaled, y_loadings_scaled)):
        ax.arrow(0, 0, x_ld, y_ld, head_width=0.05, head_length=0.05, fc='red', ec='red')
        if feature_names and i < len(feature_names):
            ax.text(x_ld * 1.1, y_ld * 1.1, feature_names[i], color='red', fontsize=8)

    # Add sample names if provided
    if sample_names:
        for i, name in enumerate(sample_names):
            if i < len(x_scores):
                ax.text(x_scores.iloc[i] if isinstance(x_scores, pd.Series) else x_scores[i],
                       y_scores.iloc[i] if isinstance(y_scores, pd.Series) else y_scores[i],
                       name, fontsize=6, alpha=0.7)

    ax.axhline(y=0, color='k', linestyle='--', alpha=0.5)
    ax.axvline(x=0, color='k', linestyle='--', alpha=0.5)
    ax.set_xlabel(f"PC{pc_x + 1}")
    ax.set_ylabel(f"PC{pc_y + 1}")
    ax.set_title(title)
    ax.legend()

    return ax

