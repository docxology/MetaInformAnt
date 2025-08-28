"""Visualization tools for single-cell data.

This module provides plotting functions for single-cell analysis results,
including quality control plots, dimensionality reduction visualizations,
and gene expression plots. All implementations use real plotting without mocking.
"""

from __future__ import annotations

from typing import Optional, List, Dict, Union, Tuple, Any
import warnings

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.colors import ListedColormap

from .preprocessing import SingleCellData


def plot_qc_metrics(
    data: SingleCellData,
    metrics: Optional[List[str]] = None,
    ncols: int = 3,
    figsize: Optional[Tuple[float, float]] = None,
) -> plt.Figure:
    """Plot quality control metrics.
    
    Args:
        data: SingleCellData object with QC metrics
        metrics: List of metrics to plot (default: standard QC metrics)
        ncols: Number of columns in subplot grid
        figsize: Figure size tuple
        
    Returns:
        Matplotlib figure object
    """
    if metrics is None:
        # Default QC metrics
        available_metrics = ['total_counts', 'n_genes', 'pct_mt', 'pct_ribo']
        metrics = [m for m in available_metrics if m in data.obs.columns]
    
    if not metrics:
        raise ValueError("No QC metrics found in data.obs")
    
    n_metrics = len(metrics)
    nrows = (n_metrics + ncols - 1) // ncols
    
    if figsize is None:
        figsize = (4 * ncols, 3 * nrows)
    
    fig, axes = plt.subplots(nrows, ncols, figsize=figsize)
    if nrows == 1 and ncols == 1:
        axes = [axes]
    elif nrows == 1 or ncols == 1:
        axes = axes.flatten()
    else:
        axes = axes.flatten()
    
    for i, metric in enumerate(metrics):
        ax = axes[i]
        
        values = data.obs[metric]
        
        # Create histogram
        ax.hist(values, bins=50, alpha=0.7, edgecolor='black', linewidth=0.5)
        ax.set_xlabel(metric.replace('_', ' ').title())
        ax.set_ylabel('Number of cells')
        ax.set_title(f'Distribution of {metric.replace("_", " ").title()}')
        
        # Add statistics text
        mean_val = np.mean(values)
        median_val = np.median(values)
        ax.axvline(mean_val, color='red', linestyle='--', alpha=0.7, label=f'Mean: {mean_val:.1f}')
        ax.axvline(median_val, color='blue', linestyle='--', alpha=0.7, label=f'Median: {median_val:.1f}')
        ax.legend()
        
        # Format axes
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
    
    # Hide unused subplots
    for i in range(n_metrics, len(axes)):
        axes[i].set_visible(False)
    
    plt.tight_layout()
    return fig


def plot_qc_scatter(
    data: SingleCellData,
    x: str = 'total_counts',
    y: str = 'n_genes',
    color: Optional[str] = 'pct_mt',
    figsize: Tuple[float, float] = (8, 6),
) -> plt.Figure:
    """Create scatter plot of QC metrics.
    
    Args:
        data: SingleCellData object
        x: X-axis metric
        y: Y-axis metric  
        color: Color-coding metric (optional)
        figsize: Figure size
        
    Returns:
        Matplotlib figure object
    """
    fig, ax = plt.subplots(figsize=figsize)
    
    # Get data
    x_vals = data.obs[x]
    y_vals = data.obs[y]
    
    if color and color in data.obs.columns:
        c_vals = data.obs[color]
        scatter = ax.scatter(x_vals, y_vals, c=c_vals, cmap='viridis', alpha=0.6, s=1)
        cbar = plt.colorbar(scatter, ax=ax)
        cbar.set_label(color.replace('_', ' ').title())
    else:
        ax.scatter(x_vals, y_vals, alpha=0.6, s=1)
    
    ax.set_xlabel(x.replace('_', ' ').title())
    ax.set_ylabel(y.replace('_', ' ').title())
    ax.set_title(f'{y.replace("_", " ").title()} vs {x.replace("_", " ").title()}')
    
    # Format axes
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    
    plt.tight_layout()
    return fig


def plot_dimensionality_reduction(
    data: SingleCellData,
    basis: str = 'umap',
    color: Optional[Union[str, List[str]]] = None,
    size: float = 1.0,
    alpha: float = 0.7,
    ncols: int = 2,
    figsize: Optional[Tuple[float, float]] = None,
    legend: bool = True,
) -> plt.Figure:
    """Plot dimensionality reduction (UMAP, t-SNE, PCA).
    
    Args:
        data: SingleCellData object
        basis: Embedding to plot ('umap', 'tsne', 'pca', 'diffmap')
        color: Variable(s) to color by
        size: Point size
        alpha: Point transparency
        ncols: Number of columns for multiple plots
        figsize: Figure size
        legend: Whether to show legend
        
    Returns:
        Matplotlib figure object
    """
    # Get coordinates
    coord_key = f'X_{basis}'
    if coord_key not in data.obsm:
        raise ValueError(f"No {basis.upper()} coordinates found. Run compute_{basis}() first.")
    
    coords = data.obsm[coord_key]
    
    # Handle color variables
    if color is None:
        color_vars = [None]
    elif isinstance(color, str):
        color_vars = [color]
    else:
        color_vars = color
    
    # Setup subplots
    n_plots = len(color_vars)
    nrows = (n_plots + ncols - 1) // ncols
    
    if figsize is None:
        figsize = (5 * ncols, 4 * nrows)
    
    fig, axes = plt.subplots(nrows, ncols, figsize=figsize)
    if n_plots == 1:
        axes = [axes] if not hasattr(axes, '__len__') else [axes]
    else:
        axes = axes.flatten() if hasattr(axes, 'flatten') else [axes]
    
    for i, color_var in enumerate(color_vars):
        ax = axes[i]
        
        x, y = coords[:, 0], coords[:, 1]
        
        if color_var is None:
            # Single color
            ax.scatter(x, y, s=size, alpha=alpha, c='lightblue', edgecolors='none')
        elif color_var in data.obs.columns:
            color_data = data.obs[color_var]
            
            if pd.api.types.is_categorical_dtype(color_data) or color_data.dtype == 'object':
                # Categorical coloring
                categories = color_data.unique()
                colors = plt.cm.tab20(np.linspace(0, 1, len(categories)))
                
                for j, cat in enumerate(categories):
                    mask = color_data == cat
                    ax.scatter(x[mask], y[mask], s=size, alpha=alpha, 
                              c=[colors[j]], label=str(cat), edgecolors='none')
                
                if legend and len(categories) < 20:  # Only show legend for reasonable number of categories
                    ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', frameon=False)
            else:
                # Continuous coloring
                scatter = ax.scatter(x, y, s=size, alpha=alpha, c=color_data, 
                                   cmap='viridis', edgecolors='none')
                cbar = plt.colorbar(scatter, ax=ax)
                cbar.set_label(color_var.replace('_', ' ').title())
        
        elif color_var in data.var.index or color_var in data.var.columns:
            # Gene expression coloring
            if color_var in data.var.index:
                gene_idx = data.var.index.get_loc(color_var)
            else:
                gene_idx = data.var[data.var['gene_name'] == color_var].index[0]
                gene_idx = data.var.index.get_loc(gene_idx)
            
            if hasattr(data.X, 'toarray'):
                expression = data.X[:, gene_idx].toarray().flatten()
            else:
                expression = data.X[:, gene_idx]
            
            scatter = ax.scatter(x, y, s=size, alpha=alpha, c=expression, 
                               cmap='Reds', edgecolors='none')
            cbar = plt.colorbar(scatter, ax=ax)
            cbar.set_label(f'{color_var} expression')
        
        else:
            # Default coloring
            ax.scatter(x, y, s=size, alpha=alpha, c='lightblue', edgecolors='none')
        
        # Format plot
        ax.set_xlabel(f'{basis.upper()}_1')
        ax.set_ylabel(f'{basis.upper()}_2')
        
        if color_var:
            ax.set_title(f'{basis.upper()} colored by {color_var}')
        else:
            ax.set_title(f'{basis.upper()} embedding')
        
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
    
    # Hide unused subplots
    for i in range(n_plots, len(axes)):
        axes[i].set_visible(False)
    
    plt.tight_layout()
    return fig


def plot_gene_expression(
    data: SingleCellData,
    genes: Union[str, List[str]],
    basis: str = 'umap',
    ncols: int = 3,
    size: float = 1.0,
    alpha: float = 0.7,
    figsize: Optional[Tuple[float, float]] = None,
    cmap: str = 'Reds',
) -> plt.Figure:
    """Plot gene expression on dimensionality reduction.
    
    Args:
        data: SingleCellData object
        genes: Gene name(s) to plot
        basis: Embedding basis
        ncols: Number of columns
        size: Point size
        alpha: Point transparency
        figsize: Figure size
        cmap: Colormap name
        
    Returns:
        Matplotlib figure object
    """
    if isinstance(genes, str):
        genes = [genes]
    
    # Get coordinates
    coord_key = f'X_{basis}'
    if coord_key not in data.obsm:
        raise ValueError(f"No {basis.upper()} coordinates found")
    
    coords = data.obsm[coord_key]
    x, y = coords[:, 0], coords[:, 1]
    
    # Setup subplots
    n_genes = len(genes)
    nrows = (n_genes + ncols - 1) // ncols
    
    if figsize is None:
        figsize = (4 * ncols, 3 * nrows)
    
    fig, axes = plt.subplots(nrows, ncols, figsize=figsize)
    if n_genes == 1:
        axes = [axes]
    else:
        axes = axes.flatten() if hasattr(axes, 'flatten') else [axes]
    
    for i, gene in enumerate(genes):
        ax = axes[i]
        
        # Find gene index
        try:
            if gene in data.var.index:
                gene_idx = data.var.index.get_loc(gene)
            elif 'gene_name' in data.var.columns and gene in data.var['gene_name'].values:
                gene_idx = data.var[data.var['gene_name'] == gene].index[0]
                gene_idx = data.var.index.get_loc(gene_idx)
            else:
                # Try partial matching
                matches = data.var.index[data.var.index.str.contains(gene, case=False)]
                if len(matches) == 0 and 'gene_name' in data.var.columns:
                    matches = data.var[data.var['gene_name'].str.contains(gene, case=False)].index
                
                if len(matches) == 0:
                    warnings.warn(f"Gene '{gene}' not found")
                    continue
                elif len(matches) > 1:
                    warnings.warn(f"Multiple matches for '{gene}': {matches[:5].tolist()}, using first")
                
                gene_idx = data.var.index.get_loc(matches[0])
            
            # Get expression values
            if hasattr(data.X, 'toarray'):
                expression = data.X[:, gene_idx].toarray().flatten()
            else:
                expression = data.X[:, gene_idx]
            
            # Plot
            scatter = ax.scatter(x, y, s=size, alpha=alpha, c=expression, 
                               cmap=cmap, edgecolors='none')
            cbar = plt.colorbar(scatter, ax=ax)
            cbar.set_label('Expression')
            
            ax.set_title(f'{gene}')
            
        except Exception as e:
            warnings.warn(f"Could not plot gene '{gene}': {e}")
            ax.text(0.5, 0.5, f"Gene '{gene}'\nnot found", 
                   ha='center', va='center', transform=ax.transAxes)
        
        ax.set_xlabel(f'{basis.upper()}_1')
        ax.set_ylabel(f'{basis.upper()}_2')
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
    
    # Hide unused subplots
    for i in range(n_genes, len(axes)):
        axes[i].set_visible(False)
    
    plt.tight_layout()
    return fig


def plot_clusters(
    data: SingleCellData,
    groupby: str,
    basis: str = 'umap',
    size: float = 1.0,
    alpha: float = 0.7,
    figsize: Tuple[float, float] = (8, 6),
    legend: bool = True,
) -> plt.Figure:
    """Plot clusters on dimensionality reduction.
    
    Args:
        data: SingleCellData object
        groupby: Clustering column name
        basis: Embedding basis
        size: Point size
        alpha: Point transparency
        figsize: Figure size
        legend: Whether to show legend
        
    Returns:
        Matplotlib figure object
    """
    if groupby not in data.obs.columns:
        raise ValueError(f"Column '{groupby}' not found in obs")
    
    # Get coordinates
    coord_key = f'X_{basis}'
    if coord_key not in data.obsm:
        raise ValueError(f"No {basis.upper()} coordinates found")
    
    coords = data.obsm[coord_key]
    x, y = coords[:, 0], coords[:, 1]
    
    fig, ax = plt.subplots(figsize=figsize)
    
    # Get cluster labels
    clusters = data.obs[groupby]
    unique_clusters = sorted(clusters.unique())
    
    # Generate colors
    if len(unique_clusters) <= 10:
        colors = plt.cm.tab10(np.linspace(0, 1, len(unique_clusters)))
    else:
        colors = plt.cm.tab20(np.linspace(0, 1, len(unique_clusters)))
    
    # Plot each cluster
    for i, cluster in enumerate(unique_clusters):
        mask = clusters == cluster
        ax.scatter(x[mask], y[mask], s=size, alpha=alpha, 
                  c=[colors[i]], label=f'Cluster {cluster}', edgecolors='none')
    
    ax.set_xlabel(f'{basis.upper()}_1')
    ax.set_ylabel(f'{basis.upper()}_2')
    ax.set_title(f'{basis.upper()} colored by {groupby}')
    
    if legend and len(unique_clusters) <= 20:
        ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', frameon=False)
    
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    
    plt.tight_layout()
    return fig


def plot_marker_genes_heatmap(
    data: SingleCellData,
    markers_df: pd.DataFrame,
    groupby: str,
    n_genes: int = 5,
    figsize: Optional[Tuple[float, float]] = None,
    cmap: str = 'RdBu_r',
) -> plt.Figure:
    """Plot heatmap of marker genes.
    
    Args:
        data: SingleCellData object
        markers_df: DataFrame from find_marker_genes()
        groupby: Clustering column
        n_genes: Number of top genes per cluster
        figsize: Figure size
        cmap: Colormap
        
    Returns:
        Matplotlib figure object
    """
    # Select top genes per group
    top_markers = markers_df.groupby('group').head(n_genes)
    gene_indices = top_markers['gene'].values
    gene_names = top_markers['gene_name'].values
    
    # Get expression data
    X = data.X
    if hasattr(X, 'toarray'):
        expr_data = X[:, gene_indices].toarray()
    else:
        expr_data = X[:, gene_indices]
    
    # Create DataFrame for plotting
    expr_df = pd.DataFrame(expr_data, columns=gene_names)
    expr_df[groupby] = data.obs[groupby].values
    
    # Calculate mean expression per cluster
    cluster_means = expr_df.groupby(groupby).mean()
    
    # Plot heatmap
    if figsize is None:
        figsize = (max(8, len(gene_names) * 0.4), max(6, len(cluster_means) * 0.3))
    
    fig, ax = plt.subplots(figsize=figsize)
    
    # Create heatmap
    im = ax.imshow(cluster_means.values, cmap=cmap, aspect='auto')
    
    # Set ticks and labels
    ax.set_xticks(range(len(gene_names)))
    ax.set_xticklabels(gene_names, rotation=45, ha='right')
    ax.set_yticks(range(len(cluster_means)))
    ax.set_yticklabels(cluster_means.index)
    
    ax.set_xlabel('Genes')
    ax.set_ylabel('Clusters')
    ax.set_title('Marker Gene Expression by Cluster')
    
    # Add colorbar
    cbar = plt.colorbar(im, ax=ax)
    cbar.set_label('Mean Expression')
    
    plt.tight_layout()
    return fig


def plot_pca_variance(
    data: SingleCellData,
    n_components: int = 50,
    figsize: Tuple[float, float] = (10, 4),
) -> plt.Figure:
    """Plot PCA variance explained.
    
    Args:
        data: SingleCellData object with PCA results
        n_components: Number of components to show
        figsize: Figure size
        
    Returns:
        Matplotlib figure object
    """
    if 'pca' not in data.uns:
        raise ValueError("No PCA results found. Run compute_pca() first.")
    
    pca_results = data.uns['pca']
    variance_ratio = pca_results['variance_ratio'][:n_components]
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=figsize)
    
    # Individual variance explained
    ax1.bar(range(1, len(variance_ratio) + 1), variance_ratio)
    ax1.set_xlabel('Principal Component')
    ax1.set_ylabel('Variance Explained')
    ax1.set_title('Variance Explained by Each PC')
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)
    
    # Cumulative variance explained
    cumulative_var = np.cumsum(variance_ratio)
    ax2.plot(range(1, len(cumulative_var) + 1), cumulative_var, 'o-')
    ax2.set_xlabel('Number of Components')
    ax2.set_ylabel('Cumulative Variance Explained')
    ax2.set_title('Cumulative Variance Explained')
    ax2.grid(True, alpha=0.3)
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)
    
    # Add text showing total variance explained
    total_var = cumulative_var[-1]
    ax2.text(0.02, 0.98, f'Total: {total_var:.1%}', 
             transform=ax2.transAxes, verticalalignment='top',
             bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    
    plt.tight_layout()
    return fig
