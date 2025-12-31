"""Single-cell data visualization functions.

This module provides plotting functions specifically designed for single-cell
RNA-seq data, including dimensionality reduction plots, trajectory visualizations,
and quality control plots.
"""

from __future__ import annotations

from typing import Dict, List, Optional, Any, Tuple
import numpy as np
import pandas as pd

from metainformant.core import logging, errors, validation

# Try to import matplotlib
try:
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches
    HAS_MATPLOTLIB = True
except ImportError:
    plt = None
    HAS_MATPLOTLIB = False
    logger.warning("matplotlib not available - visualization functions will return None")

logger = logging.get_logger(__name__)

# Import our SingleCellData
from .preprocessing import SingleCellData


def plot_umap(data: SingleCellData, color: Optional[str] = None,
             **kwargs) -> Any:
    """Create UMAP plot of single-cell data.

    Args:
        data: SingleCellData object with UMAP coordinates
        color: Column name to color points by (from obs)
        **kwargs: Additional arguments for scatter plot

    Returns:
        matplotlib Figure object

    Raises:
        ImportError: If matplotlib is not available
        TypeError: If data is not SingleCellData
        ValueError: If UMAP coordinates not found
    """
    if not HAS_MATPLOTLIB:
        raise ImportError("matplotlib required for plotting")

    validation.validate_type(data, SingleCellData, "data")

    # Find UMAP coordinates
    umap_cols = [col for col in (data.obs.columns if data.obs is not None else []) if col.upper().startswith('UMAP')]

    if len(umap_cols) < 2:
        raise errors.ValidationError("UMAP coordinates not found in data.obs. Run umap_reduction first.")

    x_col, y_col = umap_cols[0], umap_cols[1]

    fig, ax = plt.subplots(figsize=kwargs.get('figsize', (8, 6)))

    # Get coordinates
    x = data.obs[x_col].values
    y = data.obs[y_col].values

    # Color handling
    if color and data.obs is not None and color in data.obs.columns:
        colors = data.obs[color].values

        # Handle categorical colors
        if not np.issubdtype(colors.dtype, np.number):
            unique_colors = np.unique(colors)
            color_map = plt.cm.tab10(np.linspace(0, 1, len(unique_colors)))
            color_dict = dict(zip(unique_colors, color_map))

            for cat, cat_color in color_dict.items():
                mask = colors == cat
                ax.scatter(x[mask], y[mask], c=[cat_color], label=str(cat),
                          alpha=kwargs.get('alpha', 0.6), s=kwargs.get('s', 20))

            ax.legend(title=color, bbox_to_anchor=(1.05, 1), loc='upper left')
        else:
            scatter = ax.scatter(x, y, c=colors, cmap=kwargs.get('cmap', 'viridis'),
                               alpha=kwargs.get('alpha', 0.6), s=kwargs.get('s', 20))
            plt.colorbar(scatter, ax=ax, label=color)
    else:
        ax.scatter(x, y, alpha=kwargs.get('alpha', 0.6), s=kwargs.get('s', 20))

    ax.set_xlabel(x_col)
    ax.set_ylabel(y_col)
    ax.set_title('UMAP Projection')

    plt.tight_layout()
    return fig


def plot_tsne(data: SingleCellData, color: Optional[str] = None,
             **kwargs) -> Any:
    """Create t-SNE plot of single-cell data.

    Args:
        data: SingleCellData object with t-SNE coordinates
        color: Column name to color points by (from obs)
        **kwargs: Additional arguments for scatter plot

    Returns:
        matplotlib Figure object

    Raises:
        ImportError: If matplotlib is not available
        TypeError: If data is not SingleCellData
        ValueError: If t-SNE coordinates not found
    """
    if not HAS_MATPLOTLIB:
        raise ImportError("matplotlib required for plotting")

    validation.validate_type(data, SingleCellData, "data")

    # Find t-SNE coordinates
    tsne_cols = [col for col in (data.obs.columns if data.obs is not None else []) if 'tSNE' in col or 'TSNE' in col]

    if len(tsne_cols) < 2:
        raise errors.ValidationError("t-SNE coordinates not found in data.obs. Run tsne_reduction first.")

    x_col, y_col = tsne_cols[0], tsne_cols[1]

    fig, ax = plt.subplots(figsize=kwargs.get('figsize', (8, 6)))

    # Get coordinates
    x = data.obs[x_col].values
    y = data.obs[y_col].values

    # Color handling (same as UMAP)
    if color and data.obs is not None and color in data.obs.columns:
        colors = data.obs[color].values

        if not np.issubdtype(colors.dtype, np.number):
            unique_colors = np.unique(colors)
            color_map = plt.cm.tab10(np.linspace(0, 1, len(unique_colors)))
            color_dict = dict(zip(unique_colors, color_map))

            for cat, cat_color in color_dict.items():
                mask = colors == cat
                ax.scatter(x[mask], y[mask], c=[cat_color], label=str(cat),
                          alpha=kwargs.get('alpha', 0.6), s=kwargs.get('s', 20))

            ax.legend(title=color, bbox_to_anchor=(1.05, 1), loc='upper left')
        else:
            scatter = ax.scatter(x, y, c=colors, cmap=kwargs.get('cmap', 'viridis'),
                               alpha=kwargs.get('alpha', 0.6), s=kwargs.get('s', 20))
            plt.colorbar(scatter, ax=ax, label=color)
    else:
        ax.scatter(x, y, alpha=kwargs.get('alpha', 0.6), s=kwargs.get('s', 20))

    ax.set_xlabel(x_col)
    ax.set_ylabel(y_col)
    ax.set_title('t-SNE Projection')

    plt.tight_layout()
    return fig


def plot_pca(data: SingleCellData, color: Optional[str] = None,
            n_components: int = 2, **kwargs) -> Any:
    """Create PCA plot of single-cell data.

    Args:
        data: SingleCellData object with PCA coordinates
        color: Column name to color points by (from obs)
        n_components: Number of PCA components to plot (2 or 3)
        **kwargs: Additional arguments for plot

    Returns:
        matplotlib Figure object

    Raises:
        ImportError: If matplotlib is not available
        TypeError: If data is not SingleCellData
        ValueError: If PCA coordinates not found or n_components invalid
    """
    if not HAS_MATPLOTLIB:
        raise ImportError("matplotlib required for plotting")

    validation.validate_type(data, SingleCellData, "data")

    if n_components not in [2, 3]:
        raise errors.ValidationError("n_components must be 2 or 3")

    # Find PCA coordinates
    pca_cols = [col for col in (data.obs.columns if data.obs is not None else [])
               if col.upper().startswith('PC') and col[2:].isdigit()]

    if len(pca_cols) < n_components:
        raise errors.ValidationError(f"PCA coordinates not found in data.obs. Run pca_reduction first.")

    if n_components == 2:
        x_col, y_col = pca_cols[0], pca_cols[1]
        fig, ax = plt.subplots(figsize=kwargs.get('figsize', (8, 6)))

        x = data.obs[x_col].values
        y = data.obs[y_col].values

        if color and data.obs is not None and color in data.obs.columns:
            colors = data.obs[color].values

            if not np.issubdtype(colors.dtype, np.number):
                unique_colors = np.unique(colors)
                color_map = plt.cm.tab10(np.linspace(0, 1, len(unique_colors)))
                color_dict = dict(zip(unique_colors, color_map))

                for cat, cat_color in color_dict.items():
                    mask = colors == cat
                    ax.scatter(x[mask], y[mask], c=[cat_color], label=str(cat),
                              alpha=kwargs.get('alpha', 0.6), s=kwargs.get('s', 20))

                ax.legend(title=color, bbox_to_anchor=(1.05, 1), loc='upper left')
            else:
                scatter = ax.scatter(x, y, c=colors, cmap=kwargs.get('cmap', 'viridis'),
                                   alpha=kwargs.get('alpha', 0.6), s=kwargs.get('s', 20))
                plt.colorbar(scatter, ax=ax, label=color)
        else:
            ax.scatter(x, y, alpha=kwargs.get('alpha', 0.6), s=kwargs.get('s', 20))

        ax.set_xlabel(x_col)
        ax.set_ylabel(y_col)
        ax.set_title('PCA Projection')

    else:  # 3D plot
        from mpl_toolkits.mplot3d import Axes3D
        x_col, y_col, z_col = pca_cols[0], pca_cols[1], pca_cols[2]
        fig = plt.figure(figsize=kwargs.get('figsize', (10, 8)))
        ax = fig.add_subplot(111, projection='3d')

        x = data.obs[x_col].values
        y = data.obs[y_col].values
        z = data.obs[z_col].values

        if color and data.obs is not None and color in data.obs.columns:
            colors = data.obs[color].values
            scatter = ax.scatter(x, y, z, c=colors, cmap=kwargs.get('cmap', 'viridis'),
                               alpha=kwargs.get('alpha', 0.6), s=kwargs.get('s', 20))
            plt.colorbar(scatter, ax=ax, label=color)
        else:
            ax.scatter(x, y, z, alpha=kwargs.get('alpha', 0.6), s=kwargs.get('s', 20))

        ax.set_xlabel(x_col)
        ax.set_ylabel(y_col)
        ax.set_zlabel(z_col)
        ax.set_title('PCA Projection (3D)')

    plt.tight_layout()
    return fig


def plot_trajectory(data: SingleCellData, trajectory_key: str,
                   color_by_pseudotime: bool = True, **kwargs) -> Any:
    """Create trajectory plot showing pseudotime and cell ordering.

    Args:
        data: SingleCellData object with trajectory information
        trajectory_key: Key in uns containing trajectory data (e.g., 'dpt', 'paga')
        color_by_pseudotime: Whether to color by pseudotime
        **kwargs: Additional arguments for plot

    Returns:
        matplotlib Figure object

    Raises:
        ImportError: If matplotlib is not available
        TypeError: If data is not SingleCellData
        ValueError: If trajectory data not found
    """
    if not HAS_MATPLOTLIB:
        raise ImportError("matplotlib required for plotting")

    validation.validate_type(data, SingleCellData, "data")

    if data.uns is None or trajectory_key not in data.uns:
        raise errors.ValidationError(f"Trajectory data '{trajectory_key}' not found in data.uns")

    # Find embedding coordinates for plotting
    embedding_cols = None
    for method in ['umap', 'tsne', 'pca']:
        coords = [col for col in (data.obs.columns if data.obs is not None else [])
                 if method.upper() in col.upper()]
        if len(coords) >= 2:
            embedding_cols = coords[:2]
            break

    if embedding_cols is None:
        raise errors.ValidationError("No embedding coordinates found. Run dimensionality reduction first.")

    fig, ax = plt.subplots(figsize=kwargs.get('figsize', (10, 8)))

    x = data.obs[embedding_cols[0]].values
    y = data.obs[embedding_cols[1]].values

    # Color by pseudotime if available
    pseudotime_col = None
    if color_by_pseudotime:
        pseudotime_cols = [col for col in (data.obs.columns if data.obs is not None else [])
                          if 'pseudotime' in col.lower() or 'ptime' in col.lower()]
        if pseudotime_cols:
            pseudotime_col = pseudotime_cols[0]

    if pseudotime_col:
        colors = data.obs[pseudotime_col].values
        scatter = ax.scatter(x, y, c=colors, cmap=kwargs.get('cmap', 'viridis'),
                           alpha=kwargs.get('alpha', 0.7), s=kwargs.get('s', 30))
        plt.colorbar(scatter, ax=ax, label='Pseudotime')
    else:
        ax.scatter(x, y, alpha=kwargs.get('alpha', 0.7), s=kwargs.get('s', 30))

    # Add trajectory-specific elements
    trajectory_data = data.uns[trajectory_key]

    if trajectory_key == 'dpt':
        # Add root cell marker
        root_cell = trajectory_data.get('root_cell')
        if root_cell is not None and root_cell < len(x):
            ax.scatter([x[root_cell]], [y[root_cell]], c='red', s=100,
                      marker='*', edgecolors='black', linewidth=2, label='Root Cell')
            ax.legend()

    elif trajectory_key == 'paga':
        # Add connectivity information (simplified)
        # In practice, this would draw edges between connected clusters
        ax.set_title('PAGA Trajectory')

    ax.set_xlabel(embedding_cols[0])
    ax.set_ylabel(embedding_cols[1])
    ax.set_title(f'Trajectory Plot ({trajectory_key.upper()})')

    plt.tight_layout()
    return fig


def plot_marker_expression(data: SingleCellData, marker_genes: List[str],
                          method: str = 'dotplot', **kwargs) -> Any:
    """Create marker gene expression plots.

    Args:
        data: SingleCellData object
        marker_genes: List of genes to plot
        method: Plot type ('dotplot', 'heatmap', 'violin')
        **kwargs: Additional arguments for plot

    Returns:
        matplotlib Figure object

    Raises:
        ImportError: If matplotlib is not available
        TypeError: If data is not SingleCellData
        ValueError: If marker genes not found
    """
    if not HAS_MATPLOTLIB:
        raise ImportError("matplotlib required for plotting")

    validation.validate_type(data, SingleCellData, "data")

    # Check which marker genes are available
    available_genes = []
    if data.var is not None:
        gene_names = data.var.index.tolist()
        available_genes = [gene for gene in marker_genes if gene in gene_names]

    if not available_genes:
        raise errors.ValidationError(f"None of the marker genes {marker_genes} found in data")

    logger.info(f"Plotting {len(available_genes)} marker genes using {method}")

    # Get expression data
    X = data.X.toarray() if hasattr(data.X, 'toarray') else data.X
    gene_indices = [data.var.index.get_loc(gene) for gene in available_genes]

    if method == 'dotplot':
        return _plot_marker_dotplot(data, available_genes, gene_indices, X, **kwargs)
    elif method == 'heatmap':
        return _plot_marker_heatmap(data, available_genes, gene_indices, X, **kwargs)
    elif method == 'violin':
        return _plot_marker_violin(data, available_genes, gene_indices, X, **kwargs)
    else:
        raise errors.ValidationError(f"Unsupported plot method: {method}")


def _plot_marker_dotplot(data: SingleCellData, genes: List[str], gene_indices: List[int],
                        X: np.ndarray, **kwargs) -> Any:
    """Create dot plot for marker genes."""
    fig, ax = plt.subplots(figsize=kwargs.get('figsize', (len(genes) * 0.8, 6)))

    # Calculate fraction of cells expressing each gene and mean expression
    fractions = []
    means = []

    for gene_idx in gene_indices:
        expr = X[:, gene_idx]
        fraction = np.mean(expr > 0)
        mean_expr = np.mean(expr[expr > 0]) if fraction > 0 else 0

        fractions.append(fraction)
        means.append(mean_expr)

    # Create dot plot
    y_pos = np.arange(len(genes))

    # Size represents fraction expressing, color represents mean expression
    sizes = np.array(fractions) * 1000  # Scale for visibility
    sizes = np.clip(sizes, 10, 500)  # Reasonable size range

    scatter = ax.scatter(means, y_pos, s=sizes, c=means,
                        cmap=kwargs.get('cmap', 'Reds'), alpha=0.7, edgecolors='black')

    ax.set_yticks(y_pos)
    ax.set_yticklabels(genes)
    ax.set_xlabel('Mean Expression (among expressing cells)')
    ax.set_title('Marker Gene Expression (Dot Plot)')

    # Add colorbar
    plt.colorbar(scatter, ax=ax, label='Mean Expression')

    # Add size legend
    sizes_legend = [0.2, 0.5, 0.8]
    size_labels = ['20%', '50%', '80%']
    legend_elements = [plt.scatter([], [], s=size*1000, c='gray', alpha=0.7, edgecolors='black')
                      for size in sizes_legend]
    ax.legend(legend_elements, size_labels, title='Fraction Expressing',
             bbox_to_anchor=(1.05, 1), loc='upper left')

    plt.tight_layout()
    return fig


def _plot_marker_heatmap(data: SingleCellData, genes: List[str], gene_indices: List[int],
                        X: np.ndarray, **kwargs) -> Any:
    """Create heatmap for marker genes."""
    # Subsample cells for visualization if too many
    max_cells = kwargs.get('max_cells', 1000)
    n_cells = X.shape[0]

    if n_cells > max_cells:
        indices = np.random.choice(n_cells, max_cells, replace=False)
        X_plot = X[indices]
    else:
        X_plot = X
        indices = np.arange(n_cells)

    # Get expression for marker genes
    marker_expr = X_plot[:, gene_indices].T

    fig, ax = plt.subplots(figsize=kwargs.get('figsize', (12, len(genes) * 0.5)))

    # Create heatmap
    im = ax.imshow(marker_expr, aspect='auto', cmap=kwargs.get('cmap', 'viridis'),
                  interpolation='nearest')

    ax.set_yticks(np.arange(len(genes)))
    ax.set_yticklabels(genes)
    ax.set_xlabel('Cells')
    ax.set_title('Marker Gene Expression (Heatmap)')

    # Add colorbar
    plt.colorbar(im, ax=ax, label='Expression Level')

    plt.tight_layout()
    return fig


def _plot_marker_violin(data: SingleCellData, genes: List[str], gene_indices: List[int],
                       X: np.ndarray, **kwargs) -> Any:
    """Create violin plots for marker genes."""
    fig, ax = plt.subplots(figsize=kwargs.get('figsize', (len(genes) * 0.8, 6)))

    # Collect expression data for each gene
    expr_data = []
    for gene_idx in gene_indices:
        expr = X[:, gene_idx]
        # Only include non-zero expression for cleaner plots
        expr_data.append(expr[expr > 0])

    # Create violin plot
    parts = ax.violinplot(expr_data, showmeans=True, showextrema=True)

    # Customize violin plot
    for pc in parts['bodies']:
        pc.set_facecolor('lightblue')
        pc.set_edgecolor('black')
        pc.set_alpha(0.7)

    # Set labels
    ax.set_xticks(np.arange(1, len(genes) + 1))
    ax.set_xticklabels(genes, rotation=45, ha='right')
    ax.set_ylabel('Expression Level')
    ax.set_title('Marker Gene Expression (Violin Plot)')

    # Add grid
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    return fig


def plot_qc_metrics(data: SingleCellData, **kwargs) -> Any:
    """Create QC metrics visualization.

    Args:
        data: SingleCellData object with QC metrics
        **kwargs: Additional arguments for plot

    Returns:
        matplotlib Figure object

    Raises:
        ImportError: If matplotlib is not available
        TypeError: If data is not SingleCellData
    """
    if not HAS_MATPLOTLIB:
        raise ImportError("matplotlib required for plotting")

    validation.validate_type(data, SingleCellData, "data")

    if data.obs is None or 'n_counts' not in data.obs.columns:
        logger.warning("QC metrics not found. Run calculate_qc_metrics first.")
        data = data  # Would call calculate_qc_metrics in real implementation

    fig, axes = plt.subplots(2, 2, figsize=kwargs.get('figsize', (12, 10)))

    # Histogram of total counts
    if 'n_counts' in data.obs.columns:
        axes[0, 0].hist(data.obs['n_counts'], bins=50, alpha=0.7, color='blue', edgecolor='black')
        axes[0, 0].axvline(data.obs['n_counts'].median(), color='red', linestyle='--', label='Median')
        axes[0, 0].set_xlabel('Total Counts')
        axes[0, 0].set_ylabel('Number of Cells')
        axes[0, 0].set_title('Total Counts Distribution')
        axes[0, 0].legend()

    # Histogram of genes detected
    if 'n_genes' in data.obs.columns:
        axes[0, 1].hist(data.obs['n_genes'], bins=50, alpha=0.7, color='green', edgecolor='black')
        axes[0, 1].axvline(data.obs['n_genes'].median(), color='red', linestyle='--', label='Median')
        axes[0, 1].set_xlabel('Number of Genes')
        axes[0, 1].set_ylabel('Number of Cells')
        axes[0, 1].set_title('Genes Detected Distribution')
        axes[0, 1].legend()

    # Scatter plot: counts vs genes
    if 'n_counts' in data.obs.columns and 'n_genes' in data.obs.columns:
        axes[1, 0].scatter(data.obs['n_counts'], data.obs['n_genes'],
                          alpha=0.6, s=20, color='purple')
        axes[1, 0].set_xlabel('Total Counts')
        axes[1, 0].set_ylabel('Number of Genes')
        axes[1, 0].set_title('Counts vs Genes')

    # Mitochondrial percentage histogram
    if 'pct_mito' in data.obs.columns:
        axes[1, 1].hist(data.obs['pct_mito'], bins=50, alpha=0.7, color='red', edgecolor='black')
        axes[1, 1].axvline(data.obs['pct_mito'].median(), color='blue', linestyle='--', label='Median')
        axes[1, 1].set_xlabel('Mitochondrial Percentage')
        axes[1, 1].set_ylabel('Number of Cells')
        axes[1, 1].set_title('Mitochondrial Content')
        axes[1, 1].legend()

    plt.suptitle('Single-Cell QC Metrics', fontsize=16)
    plt.tight_layout()
    return fig


def plot_cluster_comparison(data: SingleCellData, cluster_cols: List[str],
                           embedding_cols: Optional[List[str]] = None, **kwargs) -> Any:
    """Compare different clustering results side by side.

    Args:
        data: SingleCellData object with multiple clusterings
        cluster_cols: List of cluster column names to compare
        embedding_cols: Columns for embedding coordinates (auto-detected if None)
        **kwargs: Additional arguments for plot

    Returns:
        matplotlib Figure object

    Raises:
        ImportError: If matplotlib is not available
        TypeError: If data is not SingleCellData
        ValueError: If cluster columns not found
    """
    if not HAS_MATPLOTLIB:
        raise ImportError("matplotlib required for plotting")

    validation.validate_type(data, SingleCellData, "data")

    if data.obs is None:
        raise errors.ValidationError("data.obs required for cluster comparison")

    # Check cluster columns
    missing_cols = [col for col in cluster_cols if col not in data.obs.columns]
    if missing_cols:
        raise errors.ValidationError(f"Cluster columns not found: {missing_cols}")

    # Find embedding coordinates
    if embedding_cols is None:
        for method in ['umap', 'tsne', 'pca']:
            coords = [col for col in data.obs.columns if method.upper() in col.upper()]
            if len(coords) >= 2:
                embedding_cols = coords[:2]
                break

    if embedding_cols is None or len(embedding_cols) < 2:
        raise errors.ValidationError("Embedding coordinates not found")

    n_clusters = len(cluster_cols)
    fig, axes = plt.subplots(1, n_clusters, figsize=kwargs.get('figsize', (6*n_clusters, 5)))

    if n_clusters == 1:
        axes = [axes]

    x = data.obs[embedding_cols[0]].values
    y = data.obs[embedding_cols[1]].values

    for i, cluster_col in enumerate(cluster_cols):
        ax = axes[i]

        clusters = data.obs[cluster_col].values
        unique_clusters = np.unique(clusters)

        # Color by cluster
        colors = plt.cm.tab10(np.linspace(0, 1, len(unique_clusters)))

        for j, cluster in enumerate(unique_clusters):
            mask = clusters == cluster
            ax.scatter(x[mask], y[mask], c=[colors[j]], label=str(cluster),
                      alpha=0.7, s=20, edgecolors='none')

        ax.set_xlabel(embedding_cols[0])
        ax.set_ylabel(embedding_cols[1])
        ax.set_title(f'{cluster_col.replace("_", " ").title()}')

        # Only show legend for first plot to avoid crowding
        if i == 0:
            ax.legend(title='Cluster', bbox_to_anchor=(1.05, 1), loc='upper left')

    plt.suptitle('Cluster Comparison', fontsize=16)
    plt.tight_layout()
    return fig
