"""Dimensionality reduction visualization functions.

This module provides specialized plotting functions for dimensionality reduction
techniques including PCA, UMAP, t-SNE, and associated diagnostic plots.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.axes import Axes

from metainformant.core import logging, paths, validation

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


def plot_pca(
    data: np.ndarray | pd.DataFrame,
    *,
    n_components: int = 2,
    ax: Axes | None = None,
    output_path: str | Path | None = None,
    **kwargs
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

    validation.validate_type(data, (np.ndarray, pd.DataFrame), "data")

    if isinstance(data, pd.DataFrame):
        data_array = data.values
    else:
        data_array = data

    if data_array.ndim != 2:
        raise ValueError("Data must be 2D")

    if n_components not in [2, 3]:
        raise ValueError("n_components must be 2 or 3")

    # Perform PCA
    pca = PCA(n_components=n_components)
    pca_result = pca.fit_transform(data_array)

    if n_components == 2:
        if ax is None:
            fig, ax = plt.subplots(figsize=kwargs.pop('figsize', (8, 6)))

        ax.scatter(
            pca_result[:, 0],
            pca_result[:, 1],
            **kwargs
        )
        ax.set_xlabel(f'PC1 ({pca.explained_variance_ratio_[0]:.1%} variance)')
        ax.set_ylabel(f'PC2 ({pca.explained_variance_ratio_[1]:.1%} variance)')
        ax.set_title('PCA Plot')
    else:  # 3D
        if ax is None:
            fig = plt.figure(figsize=kwargs.pop('figsize', (10, 8)))
            ax = fig.add_subplot(111, projection='3d')

        ax.scatter(
            pca_result[:, 0],
            pca_result[:, 1],
            pca_result[:, 2],
            **kwargs
        )
        ax.set_xlabel(f'PC1 ({pca.explained_variance_ratio_[0]:.1%} variance)')
        ax.set_ylabel(f'PC2 ({pca.explained_variance_ratio_[1]:.1%} variance)')
        ax.set_zlabel(f'PC3 ({pca.explained_variance_ratio_[2]:.1%} variance)')
        ax.set_title('3D PCA Plot')

    if output_path:
        output_path = paths.ensure_directory(Path(output_path).parent)
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        logger.info(f"PCA plot saved to {output_path}")

    return ax


def plot_umap(
    data: np.ndarray | pd.DataFrame,
    *,
    n_components: int = 2,
    ax: Axes | None = None,
    output_path: str | Path | None = None,
    **kwargs
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

    validation.validate_type(data, (np.ndarray, pd.DataFrame), "data")

    if isinstance(data, pd.DataFrame):
        data_array = data.values
    else:
        data_array = data

    if data_array.ndim != 2:
        raise ValueError("Data must be 2D")

    if n_components not in [2, 3]:
        raise ValueError("n_components must be 2 or 3")

    # Perform UMAP
    reducer = umap.UMAP(n_components=n_components, random_state=42)
    umap_result = reducer.fit_transform(data_array)

    if n_components == 2:
        if ax is None:
            fig, ax = plt.subplots(figsize=kwargs.pop('figsize', (8, 6)))

        ax.scatter(
            umap_result[:, 0],
            umap_result[:, 1],
            **kwargs
        )
        ax.set_xlabel('UMAP 1')
        ax.set_ylabel('UMAP 2')
        ax.set_title('UMAP Plot')
    else:  # 3D
        if ax is None:
            fig = plt.figure(figsize=kwargs.pop('figsize', (10, 8)))
            ax = fig.add_subplot(111, projection='3d')

        ax.scatter(
            umap_result[:, 0],
            umap_result[:, 1],
            umap_result[:, 2],
            **kwargs
        )
        ax.set_xlabel('UMAP 1')
        ax.set_ylabel('UMAP 2')
        ax.set_zlabel('UMAP 3')
        ax.set_title('3D UMAP Plot')

    if output_path:
        output_path = paths.ensure_directory(Path(output_path).parent)
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        logger.info(f"UMAP plot saved to {output_path}")

    return ax


def plot_tsne(
    data: np.ndarray | pd.DataFrame,
    *,
    n_components: int = 2,
    ax: Axes | None = None,
    output_path: str | Path | None = None,
    **kwargs
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

    validation.validate_type(data, (np.ndarray, pd.DataFrame), "data")

    if isinstance(data, pd.DataFrame):
        data_array = data.values
    else:
        data_array = data

    if data_array.ndim != 2:
        raise ValueError("Data must be 2D")

    if n_components not in [2, 3]:
        raise ValueError("n_components must be 2 or 3")

    # Perform t-SNE (may be slow for large datasets)
    tsne = TSNE(n_components=n_components, random_state=42)
    tsne_result = tsne.fit_transform(data_array)

    if n_components == 2:
        if ax is None:
            fig, ax = plt.subplots(figsize=kwargs.pop('figsize', (8, 6)))

        ax.scatter(
            tsne_result[:, 0],
            tsne_result[:, 1],
            **kwargs
        )
        ax.set_xlabel('t-SNE 1')
        ax.set_ylabel('t-SNE 2')
        ax.set_title('t-SNE Plot')
    else:  # 3D
        if ax is None:
            fig = plt.figure(figsize=kwargs.pop('figsize', (10, 8)))
            ax = fig.add_subplot(111, projection='3d')

        ax.scatter(
            tsne_result[:, 0],
            tsne_result[:, 1],
            tsne_result[:, 2],
            **kwargs
        )
        ax.set_xlabel('t-SNE 1')
        ax.set_ylabel('t-SNE 2')
        ax.set_zlabel('t-SNE 3')
        ax.set_title('3D t-SNE Plot')

    if output_path:
        output_path = paths.ensure_directory(Path(output_path).parent)
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        logger.info(f"t-SNE plot saved to {output_path}")

    return ax


def plot_pca_loadings(
    pca_model: Any,
    *,
    n_components: int = 2,
    ax: Axes | None = None,
    output_path: str | Path | None = None,
    **kwargs
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
    if not hasattr(pca_model, 'components_'):
        raise ValueError("pca_model must be a fitted sklearn PCA model")

    if ax is None:
        fig, ax = plt.subplots(figsize=kwargs.pop('figsize', (10, 8)))

    loadings = pca_model.components_[:2].T  # Shape: (n_features, 2)

    ax.scatter(
        loadings[:, 0],
        loadings[:, 1],
        **kwargs
    )

    # Add feature labels if available
    if hasattr(pca_model, '_feature_names_in'):
        feature_names = pca_model._feature_names_in
        for i, name in enumerate(feature_names):
            ax.annotate(name, (loadings[i, 0], loadings[i, 1]),
                       xytext=(5, 5), textcoords='offset points',
                       fontsize=8, alpha=0.8)

    ax.axhline(y=0, color='k', linestyle='--', alpha=0.5)
    ax.axvline(x=0, color='k', linestyle='--', alpha=0.5)
    ax.set_xlabel(f'PC1 Loadings ({pca_model.explained_variance_ratio_[0]:.1%} variance)')
    ax.set_ylabel(f'PC2 Loadings ({pca_model.explained_variance_ratio_[1]:.1%} variance)')
    ax.set_title('PCA Loadings Plot')

    if output_path:
        output_path = paths.ensure_directory(Path(output_path).parent)
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        logger.info(f"PCA loadings plot saved to {output_path}")

    return ax


def biplot(
    data: np.ndarray | pd.DataFrame,
    pca_model: Any,
    *,
    ax: Axes | None = None,
    output_path: str | Path | None = None,
    **kwargs
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
        sample_names = data.index.tolist()
        feature_names = data.columns.tolist()
    else:
        data_array = data
        sample_names = [f'Sample_{i}' for i in range(data_array.shape[0])]
        feature_names = [f'Feature_{i}' for i in range(data_array.shape[1])]

    # Validate PCA model
    if not hasattr(pca_model, 'components_'):
        raise ValueError("pca_model must be a fitted sklearn PCA model")

    if ax is None:
        fig, ax = plt.subplots(figsize=kwargs.pop('figsize', (10, 8)))

    # Get PCA scores
    scores = pca_model.transform(data_array)

    # Plot samples
    ax.scatter(
        scores[:, 0],
        scores[:, 1],
        alpha=0.7,
        label='Samples',
        **kwargs
    )

    # Plot loadings as arrows
    loadings = pca_model.components_[:2].T
    scaling_factor = kwargs.get('scaling_factor', scores.std() / loadings.std() * 0.8)

    for i, (x, y) in enumerate(loadings):
        ax.arrow(0, 0, x * scaling_factor, y * scaling_factor,
                head_width=0.05, head_length=0.05, fc='red', ec='red', alpha=0.7)
        # Label features
        if len(feature_names) <= 20:  # Only label if not too many
            ax.text(x * scaling_factor * 1.1, y * scaling_factor * 1.1,
                   feature_names[i], fontsize=8, ha='center', va='center')

    ax.axhline(y=0, color='k', linestyle='--', alpha=0.3)
    ax.axvline(x=0, color='k', linestyle='--', alpha=0.3)
    ax.set_xlabel(f'PC1 ({pca_model.explained_variance_ratio_[0]:.1%} variance)')
    ax.set_ylabel(f'PC2 ({pca_model.explained_variance_ratio_[1]:.1%} variance)')
    ax.set_title('PCA Biplot')
    ax.legend()

    if output_path:
        output_path = paths.ensure_directory(Path(output_path).parent)
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        logger.info(f"PCA biplot saved to {output_path}")

    return ax

