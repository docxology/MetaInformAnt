"""Dimensionality reduction visualization functions.

This module provides specialized plotting functions for dimensionality reduction
techniques including PCA, UMAP, t-SNE, and associated diagnostic plots.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.axes import Axes

from metainformant.core.data import validation
from metainformant.core.io import paths
from metainformant.core.utils import logging
from metainformant.visualization.config.conventions import save_figure_deterministic

logger = logging.get_logger(__name__)

# Optional imports with graceful fallbacks
try:
    import seaborn as sns

    HAS_SEABORN = True
except ImportError:
    sns = None
    HAS_SEABORN = False

try:
    from sklearn.decomposition import PCA

    HAS_SKLEARN = True
except ImportError:
    PCA = None
    HAS_SKLEARN = False


def _save_plot(ax: Axes, output_path: str | Path, label: str) -> str:
    """Ensure the output directory exists and save the plotted figure deterministically.

    Consolidates the repeated ensure_directory / save_figure_deterministic / logger
    triple used by every plot function in this module. Saves ``ax.figure`` so the
    plotted figure is written even when a different pyplot figure is current.
    """
    paths.ensure_directory(Path(output_path).parent)
    save_figure_deterministic(ax.figure, output_path, dpi=300, bbox_inches="tight")
    logger.info(f"{label} saved to {output_path}")
    return str(output_path)


def _prepare_embedding_data(data: np.ndarray | pd.DataFrame, n_components: int) -> np.ndarray:
    """Validate embedding input data and component count (shared by PCA/UMAP/t-SNE plots)."""
    validation.validate_type(data, (np.ndarray, pd.DataFrame), "data")

    data_array = data.values if isinstance(data, pd.DataFrame) else data

    if data_array.ndim != 2:
        raise ValueError("Data must be 2D")

    if n_components not in (2, 3):
        raise ValueError("n_components must be 2 or 3")

    return data_array


def _plot_embedding(
    ax: Axes | None,
    result: np.ndarray,
    n_components: int,
    axis_labels: tuple[str, str, str],
    titles: tuple[str, str],
    **kwargs: Any,
) -> Axes:
    """Scatter a 2- or 3-column embedding onto ``ax`` (creating axes as needed).

    Shared scaffolding for the PCA/UMAP/t-SNE scatter plots; ``axis_labels`` are
    the x/y(/z) axis labels and ``titles`` the 2D/3D plot titles.
    """
    if n_components == 2:
        if ax is None:
            fig, ax = plt.subplots(figsize=kwargs.pop("figsize", (8, 6)))

        ax.scatter(result[:, 0], result[:, 1], **kwargs)
        ax.set_xlabel(axis_labels[0])
        ax.set_ylabel(axis_labels[1])
        ax.set_title(titles[0])
    else:  # 3D
        if ax is None:
            fig = plt.figure(figsize=kwargs.pop("figsize", (10, 8)))
            ax = fig.add_subplot(111, projection="3d")

        # mpl_toolkits.mplot3d exposes no mypy-visible stubs; the 3D axes are
        # dynamically typed here so scatter's zs slot and set_zlabel resolve.
        ax_3d: Any = ax
        ax_3d.scatter(result[:, 0], result[:, 1], result[:, 2], **kwargs)
        ax_3d.set_xlabel(axis_labels[0])
        ax_3d.set_ylabel(axis_labels[1])
        ax_3d.set_zlabel(axis_labels[2])
        ax_3d.set_title(titles[1])

    return ax


def plot_pca(
    data: np.ndarray | pd.DataFrame,
    *,
    n_components: int = 2,
    ax: Axes | None = None,
    output_path: str | Path | None = None,
    **kwargs: Any,
) -> Axes:
    """Create a PCA scatter plot.

    Args:
        data: Input data matrix or DataFrame
        n_components: Number of PCA components to show (2 or 3)
        ax: Optional matplotlib axes to plot on. Creates new if None.
        output_path: Optional path to save the figure.
        **kwargs: Additional arguments passed to matplotlib scatter().

    Returns:
        matplotlib Axes object

    Raises:
        ValueError: If n_components is not 2 or 3, or data is invalid
        ImportError: If scikit-learn is not available
    """
    if not HAS_SKLEARN:
        raise ImportError("scikit-learn required for PCA plotting")

    data_array = _prepare_embedding_data(data, n_components)

    # Perform PCA
    pca = PCA(n_components=n_components)
    pca_result = pca.fit_transform(data_array)
    evr = pca.explained_variance_ratio_

    ax = _plot_embedding(
        ax,
        pca_result,
        n_components,
        axis_labels=(
            f"PC1 ({evr[0]:.1%} variance)",
            f"PC2 ({evr[1]:.1%} variance)",
            f"PC3 ({evr[2]:.1%} variance)" if n_components == 3 else "",
        ),
        titles=("PCA Plot", "3D PCA Plot"),
        **kwargs,
    )

    if output_path:
        _save_plot(ax, output_path, "PCA plot")

    return ax


def plot_umap(
    data: np.ndarray | pd.DataFrame,
    *,
    n_components: int = 2,
    ax: Axes | None = None,
    output_path: str | Path | None = None,
    **kwargs: Any,
) -> Axes:
    """Create a UMAP scatter plot.

    Args:
        data: Input data matrix or DataFrame
        n_components: Number of UMAP components to show (2 or 3)
        ax: Optional matplotlib axes to plot on. Creates new if None.
        output_path: Optional path to save the figure.
        **kwargs: Additional arguments passed to matplotlib scatter().

    Returns:
        matplotlib Axes object

    Raises:
        ValueError: If n_components is not 2 or 3, or data is invalid
        ImportError: If umap-learn is not available
    """
    try:
        import umap
    except ImportError:
        raise ImportError("umap-learn required for UMAP plotting")

    data_array = _prepare_embedding_data(data, n_components)

    # Perform UMAP
    reducer = umap.UMAP(n_components=n_components, random_state=42)
    umap_result = reducer.fit_transform(data_array)

    ax = _plot_embedding(
        ax,
        umap_result,
        n_components,
        axis_labels=("UMAP 1", "UMAP 2", "UMAP 3"),
        titles=("UMAP Plot", "3D UMAP Plot"),
        **kwargs,
    )

    if output_path:
        _save_plot(ax, output_path, "UMAP plot")

    return ax


def plot_tsne(
    data: np.ndarray | pd.DataFrame,
    *,
    n_components: int = 2,
    ax: Axes | None = None,
    output_path: str | Path | None = None,
    **kwargs: Any,
) -> Axes:
    """Create a t-SNE scatter plot.

    Args:
        data: Input data matrix or DataFrame
        n_components: Number of t-SNE components to show (2 or 3)
        ax: Optional matplotlib axes to plot on. Creates new if None.
        output_path: Optional path to save the figure.
        **kwargs: Additional arguments passed to matplotlib scatter().

    Returns:
        matplotlib Axes object

    Raises:
        ValueError: If n_components is not 2 or 3, or data is invalid
        ImportError: If scikit-learn is not available
    """
    if not HAS_SKLEARN:
        raise ImportError("scikit-learn required for t-SNE plotting")

    try:
        from sklearn.manifold import TSNE
    except ImportError:
        raise ImportError("scikit-learn required for t-SNE plotting")

    data_array = _prepare_embedding_data(data, n_components)

    # Perform t-SNE (may be slow for large datasets)
    # Perplexity must be less than n_samples; default to min(30, n-1)
    perplexity = min(30.0, data_array.shape[0] - 1)
    tsne = TSNE(n_components=n_components, perplexity=perplexity, random_state=42)
    tsne_result = tsne.fit_transform(data_array)

    ax = _plot_embedding(
        ax,
        tsne_result,
        n_components,
        axis_labels=("t-SNE 1", "t-SNE 2", "t-SNE 3"),
        titles=("t-SNE Plot", "3D t-SNE Plot"),
        **kwargs,
    )

    if output_path:
        _save_plot(ax, output_path, "t-SNE plot")

    return ax


def plot_pca_loadings(
    pca_model: Any,
    *,
    n_components: int = 2,
    ax: Axes | None = None,
    output_path: str | Path | None = None,
    **kwargs: Any,
) -> Axes:
    """Create a PCA loadings plot.

    Args:
        pca_model: Fitted sklearn PCA model
        n_components: Number of components to show (must be 2)
        ax: Optional matplotlib axes to plot on. Creates new if None.
        output_path: Optional path to save the figure.
        **kwargs: Additional arguments passed to matplotlib scatter().

    Returns:
        matplotlib Axes object

    Raises:
        ValueError: If n_components is not 2 or PCA model is invalid
    """
    if not HAS_SKLEARN:
        raise ImportError("scikit-learn required for PCA loadings plotting")

    if n_components != 2:
        raise ValueError("PCA loadings plot only supports 2 components")

    # Validate PCA model
    if not hasattr(pca_model, "components_"):
        raise ValueError("pca_model must be a fitted sklearn PCA model")

    if ax is None:
        fig, ax = plt.subplots(figsize=kwargs.pop("figsize", (10, 8)))

    loadings = pca_model.components_[:2].T  # Shape: (n_features, 2)

    ax.scatter(loadings[:, 0], loadings[:, 1], **kwargs)

    # Add feature labels if available
    if hasattr(pca_model, "feature_names_in_"):
        feature_names = pca_model.feature_names_in_
        for i, name in enumerate(feature_names):
            ax.annotate(
                name, (loadings[i, 0], loadings[i, 1]), xytext=(5, 5), textcoords="offset points", fontsize=8, alpha=0.8
            )

    ax.axhline(y=0, color="k", linestyle="--", alpha=0.5)
    ax.axvline(x=0, color="k", linestyle="--", alpha=0.5)
    ax.set_xlabel(f"PC1 Loadings ({pca_model.explained_variance_ratio_[0]:.1%} variance)")
    ax.set_ylabel(f"PC2 Loadings ({pca_model.explained_variance_ratio_[1]:.1%} variance)")
    ax.set_title("PCA Loadings Plot")

    if output_path:
        _save_plot(ax, output_path, "PCA loadings plot")

    return ax


def biplot(
    data: np.ndarray | pd.DataFrame,
    pca_model: Any,
    *,
    ax: Axes | None = None,
    output_path: str | Path | None = None,
    **kwargs: Any,
) -> Axes:
    """Create a PCA biplot showing both samples and loadings.

    Args:
        data: Input data matrix or DataFrame
        pca_model: Fitted sklearn PCA model
        ax: Optional matplotlib axes to plot on. Creates new if None.
        output_path: Optional path to save the figure.
        **kwargs: Additional arguments passed to plotting functions.

    Returns:
        matplotlib Axes object

    Raises:
        ValueError: If data or PCA model is invalid
        ImportError: If scikit-learn is not available
    """
    if not HAS_SKLEARN:
        raise ImportError("scikit-learn required for biplot")

    validation.validate_type(data, (np.ndarray, pd.DataFrame), "data")

    if isinstance(data, pd.DataFrame):
        data_array = data.values
        feature_names = data.columns.tolist()
    else:
        data_array = data
        feature_names = [f"Feature_{i}" for i in range(data_array.shape[1])]

    # Validate PCA model
    if not hasattr(pca_model, "components_"):
        raise ValueError("pca_model must be a fitted sklearn PCA model")

    if ax is None:
        fig, ax = plt.subplots(figsize=kwargs.pop("figsize", (10, 8)))

    # Get PCA scores
    scores = pca_model.transform(data_array)

    # Plot samples
    ax.scatter(scores[:, 0], scores[:, 1], alpha=0.7, label="Samples", **kwargs)

    # Plot loadings as arrows
    loadings = pca_model.components_[:2].T
    scaling_factor = kwargs.get("scaling_factor", scores.std() / loadings.std() * 0.8)

    for i, (x, y) in enumerate(loadings):
        ax.arrow(
            0,
            0,
            x * scaling_factor,
            y * scaling_factor,
            head_width=0.05,
            head_length=0.05,
            fc="red",
            ec="red",
            alpha=0.7,
        )
        # Label features
        if len(feature_names) <= 20:  # Only label if not too many
            ax.text(
                x * scaling_factor * 1.1,
                y * scaling_factor * 1.1,
                feature_names[i],
                fontsize=8,
                ha="center",
                va="center",
            )

    ax.axhline(y=0, color="k", linestyle="--", alpha=0.3)
    ax.axvline(x=0, color="k", linestyle="--", alpha=0.3)
    ax.set_xlabel(f"PC1 ({pca_model.explained_variance_ratio_[0]:.1%} variance)")
    ax.set_ylabel(f"PC2 ({pca_model.explained_variance_ratio_[1]:.1%} variance)")
    ax.set_title("PCA Biplot")
    ax.legend()

    if output_path:
        _save_plot(ax, output_path, "PCA biplot")

    return ax
