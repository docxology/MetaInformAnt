from __future__ import annotations

from typing import Iterable, Sequence

import matplotlib

# Use non-interactive backend by default for tests/headless
matplotlib.use("Agg", force=True)
import matplotlib.pyplot as plt  # noqa: E402
import pandas as pd  # noqa: E402

try:
    import seaborn as sns  # noqa: E402

    SEABORN_AVAILABLE = True
except ImportError:
    SEABORN_AVAILABLE = False
    sns = None


def lineplot(
    x: Sequence[float] | None,
    y: Sequence[float],
    *,
    label: str | None = None,
    ax: plt.Axes | None = None,
    style: str = "-",
    color: str | None = None,
) -> plt.Axes:
    """Simple line plot wrapper returning the Axes.

    If x is None, indices are used.
    """
    if ax is None:
        _, ax = plt.subplots()
    if x is None:
        x = list(range(len(y)))
    ax.plot(x, y, style, label=label, color=color)
    if label:
        ax.legend()
    return ax


def heatmap(
    data: Sequence[Sequence[float]] | pd.DataFrame,
    *,
    cmap: str = "viridis",
    cbar: bool = True,
    ax: plt.Axes | None = None,
    annot: bool = False,
) -> plt.Axes:
    """Heatmap using seaborn; accepts 2D sequences or DataFrame."""
    if ax is None:
        _, ax = plt.subplots()

    if SEABORN_AVAILABLE:
        sns.heatmap(data, cmap=cmap, cbar=cbar, ax=ax, annot=annot)
    else:
        # Fallback to matplotlib
        import numpy as np

        data_array = np.array(data)
        im = ax.imshow(data_array, cmap=cmap)
        if cbar:
            plt.colorbar(im, ax=ax)
        if annot:
            for i in range(data_array.shape[0]):
                for j in range(data_array.shape[1]):
                    ax.text(j, i, f"{data_array[i, j]:.2f}", ha="center", va="center")
    return ax


def pairplot_dataframe(df: pd.DataFrame, *, hue: str | None = None):
    """Pairplot for a tidy DataFrame, returns seaborn PairGrid or matplotlib figure."""
    if SEABORN_AVAILABLE:
        grid = sns.pairplot(df, hue=hue)
        return grid
    else:
        # Fallback: simple scatter plot matrix with matplotlib
        import numpy as np

        numeric_cols = df.select_dtypes(include=[np.number]).columns
        n_cols = len(numeric_cols)

        fig, axes = plt.subplots(n_cols, n_cols, figsize=(n_cols * 3, n_cols * 3))
        if n_cols == 1:
            axes = np.array([[axes]])
        elif n_cols == 2:
            axes = axes.reshape(2, 2)

        for i, col1 in enumerate(numeric_cols):
            for j, col2 in enumerate(numeric_cols):
                ax = axes[i, j]
                if i == j:
                    # Diagonal: histogram
                    ax.hist(df[col1], alpha=0.7)
                    ax.set_xlabel(col1)
                else:
                    # Off-diagonal: scatter plot
                    ax.scatter(df[col2], df[col1], alpha=0.6)
                    ax.set_xlabel(col2)
                    ax.set_ylabel(col1)

        plt.tight_layout()
        return fig


def scatter_plot(
    x: Sequence[float],
    y: Sequence[float],
    *,
    ax: plt.Axes | None = None,
    color: str | Sequence[str] | None = None,
    size: float | Sequence[float] = 20,
    alpha: float = 0.7,
    xlabel: str | None = None,
    ylabel: str | None = None,
    title: str | None = None,
) -> plt.Axes:
    """Create a scatter plot.

    Args:
        x: X-axis values
        y: Y-axis values
        ax: Matplotlib axes (creates new if None)
        color: Point colors
        size: Point sizes
        alpha: Point transparency
        xlabel: X-axis label
        ylabel: Y-axis label
        title: Plot title

    Returns:
        Matplotlib axes object
    """
    if ax is None:
        _, ax = plt.subplots()

    ax.scatter(x, y, c=color, s=size, alpha=alpha)

    if xlabel:
        ax.set_xlabel(xlabel)
    if ylabel:
        ax.set_ylabel(ylabel)
    if title:
        ax.set_title(title)

    return ax


def histogram(
    data: Sequence[float],
    *,
    bins: int | str | Sequence[float] = 30,
    ax: plt.Axes | None = None,
    density: bool = False,
    alpha: float = 0.7,
    color: str = "blue",
    xlabel: str | None = None,
    ylabel: str | None = None,
    title: str | None = None,
) -> plt.Axes:
    """Create a histogram.

    Args:
        data: Data to plot
        bins: Number of bins or bin edges
        ax: Matplotlib axes (creates new if None)
        density: Whether to normalize to density
        alpha: Bar transparency
        color: Bar color
        xlabel: X-axis label
        ylabel: Y-axis label
        title: Plot title

    Returns:
        Matplotlib axes object
    """
    if ax is None:
        _, ax = plt.subplots()

    ax.hist(data, bins=bins, density=density, alpha=alpha, color=color)

    if xlabel:
        ax.set_xlabel(xlabel)
    if ylabel:
        ax.set_ylabel(ylabel)
    if title:
        ax.set_title(title)

    return ax


def box_plot(
    data: Sequence[float] | Sequence[Sequence[float]],
    *,
    ax: plt.Axes | None = None,
    positions: Sequence[float] | None = None,
    labels: Sequence[str] | None = None,
    xlabel: str | None = None,
    ylabel: str | None = None,
    title: str | None = None,
) -> plt.Axes:
    """Create a box plot.

    Args:
        data: Data to plot (single array or list of arrays)
        ax: Matplotlib axes (creates new if None)
        positions: Box positions
        labels: Box labels
        xlabel: X-axis label
        ylabel: Y-axis label
        title: Plot title

    Returns:
        Matplotlib axes object
    """
    if ax is None:
        _, ax = plt.subplots()

    bp = ax.boxplot(data, positions=positions, labels=labels)

    if xlabel:
        ax.set_xlabel(xlabel)
    if ylabel:
        ax.set_ylabel(ylabel)
    if title:
        ax.set_title(title)

    return ax


def violin_plot(
    data: Sequence[float] | Sequence[Sequence[float]],
    *,
    ax: plt.Axes | None = None,
    positions: Sequence[float] | None = None,
    labels: Sequence[str] | None = None,
    xlabel: str | None = None,
    ylabel: str | None = None,
    title: str | None = None,
) -> plt.Axes:
    """Create a violin plot.

    Args:
        data: Data to plot (single array or list of arrays)
        ax: Matplotlib axes (creates new if None)
        positions: Violin positions
        labels: Violin labels
        xlabel: X-axis label
        ylabel: Y-axis label
        title: Plot title

    Returns:
        Matplotlib axes object
    """
    if ax is None:
        _, ax = plt.subplots()

    if SEABORN_AVAILABLE:
        if isinstance(data[0], (list, tuple)) if data else False:
            # Multiple datasets
            for i, dataset in enumerate(data):
                pos = positions[i] if positions else i + 1
                sns.violinplot(data=dataset, ax=ax, position=pos)
        else:
            # Single dataset
            sns.violinplot(data=data, ax=ax)
    else:
        # Fallback to boxplot
        return box_plot(data, ax=ax, positions=positions, labels=labels,
                       xlabel=xlabel, ylabel=ylabel, title=title)

    if labels:
        ax.set_xticklabels(labels)
    if xlabel:
        ax.set_xlabel(xlabel)
    if ylabel:
        ax.set_ylabel(ylabel)
    if title:
        ax.set_title(title)

    return ax


def correlation_heatmap(
    data: pd.DataFrame,
    *,
    method: str = "pearson",
    cmap: str = "RdBu_r",
    annot: bool = True,
    ax: plt.Axes | None = None,
    **kwargs
) -> plt.Axes:
    """Create a correlation heatmap.

    Args:
        data: DataFrame with numeric columns
        method: Correlation method ('pearson', 'spearman', 'kendall')
        cmap: Colormap for the heatmap
        annot: Whether to annotate cells with correlation values
        ax: Matplotlib axes (creates new if None)
        **kwargs: Additional arguments for seaborn.heatmap

    Returns:
        Matplotlib axes object
    """
    if ax is None:
        _, ax = plt.subplots()

    if SEABORN_AVAILABLE:
        corr_matrix = data.corr(method=method)
        sns.heatmap(corr_matrix, cmap=cmap, annot=annot, ax=ax, **kwargs)
    else:
        # Fallback to matplotlib
        import numpy as np
        corr_matrix = data.corr(method=method)
        im = ax.imshow(corr_matrix, cmap=cmap)
        if annot:
            for i in range(len(corr_matrix)):
                for j in range(len(corr_matrix)):
                    ax.text(j, i, f"{corr_matrix.iloc[i, j]:.2f}",
                           ha="center", va="center", color="white" if abs(corr_matrix.iloc[i, j]) > 0.5 else "black")

    return ax


def volcano_plot(
    data: pd.DataFrame,
    x_col: str,
    y_col: str,
    *,
    p_threshold: float = 0.05,
    fc_threshold: float = 1.0,
    ax: plt.Axes | None = None,
    **kwargs
) -> plt.Axes:
    """Create a volcano plot for differential expression analysis.

    Args:
        data: DataFrame with log fold change and p-value columns
        x_col: Column name for x-axis (usually log fold change)
        y_col: Column name for y-axis (usually -log10(p-value))
        p_threshold: P-value threshold for significance
        fc_threshold: Fold change threshold
        ax: Matplotlib axes (creates new if None)
        **kwargs: Additional arguments for scatter

    Returns:
        Matplotlib axes object
    """
    if ax is None:
        _, ax = plt.subplots()

    # Create significance mask
    significant = (data[y_col] > -np.log10(p_threshold)) & (abs(data[x_col]) > fc_threshold)

    # Plot non-significant points
    ax.scatter(data.loc[~significant, x_col], data.loc[~significant, y_col],
              color='gray', alpha=0.5, label='Not significant', **kwargs)

    # Plot significant points
    ax.scatter(data.loc[significant, x_col], data.loc[significant, y_col],
              color='red', alpha=0.8, label='Significant', **kwargs)

    # Add threshold lines
    ax.axvline(-fc_threshold, color='black', linestyle='--', alpha=0.7)
    ax.axvline(fc_threshold, color='black', linestyle='--', alpha=0.7)
    ax.axhline(-np.log10(p_threshold), color='black', linestyle='--', alpha=0.7)

    ax.set_xlabel(x_col)
    ax.set_ylabel(y_col)
    ax.legend()

    return ax


def expression_heatmap(
    data: pd.DataFrame,
    *,
    row_cluster: bool = True,
    col_cluster: bool = True,
    cmap: str = "RdYlBu_r",
    ax: plt.Axes | None = None,
    **kwargs
) -> plt.Axes:
    """Create an expression heatmap with clustering.

    Args:
        data: Expression data (genes x samples)
        row_cluster: Whether to cluster rows (genes)
        col_cluster: Whether to cluster columns (samples)
        cmap: Colormap for expression values
        ax: Matplotlib axes (creates new if None)
        **kwargs: Additional arguments for seaborn.clustermap

    Returns:
        Matplotlib axes object
    """
    if ax is None:
        _, ax = plt.subplots()

    if SEABORN_AVAILABLE:
        g = sns.clustermap(data, row_cluster=row_cluster, col_cluster=col_cluster,
                          cmap=cmap, **kwargs)
        return g.ax_heatmap
    else:
        # Fallback to basic heatmap
        return heatmap(data, cmap=cmap, ax=ax)


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
    """
    if ax is None:
        _, ax = plt.subplots()

    pc_x_col = f"PC{pc_x}"
    pc_y_col = f"PC{pc_y}"

    if pc_x_col not in data.columns or pc_y_col not in data.columns:
        raise ValueError(f"Columns {pc_x_col} and {pc_y_col} not found in data")

    scatter = ax.scatter(data[pc_x_col], data[pc_y_col], c=hue, **kwargs)

    ax.set_xlabel(f"{pc_x_col} ({data.get(f'{pc_x_col}_variance', 'unknown')}%)")
    ax.set_ylabel(f"{pc_y_col} ({data.get(f'{pc_y_col}_variance', 'unknown')}%)")
    ax.set_title(f"PCA Plot ({pc_x_col} vs {pc_y_col})")

    return ax


def manhattan_plot(
    data: pd.DataFrame,
    x_col: str,
    y_col: str,
    chromosome_col: str = "chromosome",
    *,
    p_threshold: float = 5e-8,
    highlight_color: str = "red",
    ax: plt.Axes | None = None,
    **kwargs
) -> plt.Axes:
    """Create a Manhattan plot for genome-wide association studies.

    Args:
        data: DataFrame with genomic positions and p-values
        x_col: Column name for x-axis (genomic position)
        y_col: Column name for y-axis (usually -log10(p-value))
        chromosome_col: Column name for chromosome information
        p_threshold: P-value threshold for significance
        highlight_color: Color for significant points
        ax: Matplotlib axes (creates new if None)
        **kwargs: Additional arguments for scatter

    Returns:
        Matplotlib axes object
    """
    if ax is None:
        _, ax = plt.subplots(figsize=(12, 6))

    # Create significance mask
    significant = data[y_col] > -np.log10(p_threshold)

    # Plot non-significant points
    if not significant.all():
        ax.scatter(data.loc[~significant, x_col], data.loc[~significant, y_col],
                  color='gray', alpha=0.6, s=10, **kwargs)

    # Plot significant points
    if significant.any():
        ax.scatter(data.loc[significant, x_col], data.loc[significant, y_col],
                  color=highlight_color, alpha=0.8, s=20, **kwargs)

    # Add threshold line
    ax.axhline(-np.log10(p_threshold), color='black', linestyle='--', alpha=0.7, label=f'P = {p_threshold}')

    # Group by chromosome for alternating colors
    chromosomes = data[chromosome_col].unique()
    colors = plt.cm.tab20(np.linspace(0, 1, len(chromosomes)))

    for i, chrom in enumerate(chromosomes):
        chrom_data = data[data[chromosome_col] == chrom]
        if not significant[chrom_data.index].all():
            ax.scatter(chrom_data.loc[~significant[chrom_data.index], x_col],
                      chrom_data.loc[~significant[chrom_data.index], y_col],
                      color=colors[i], alpha=0.6, s=10)

    ax.set_xlabel("Genomic Position")
    ax.set_ylabel("-log₁₀(P-value)")
    ax.set_title("Manhattan Plot")
    ax.legend()

    return ax


def qq_plot(
    p_values: list[float],
    *,
    ax: plt.Axes | None = None,
    title: str = "Q-Q Plot",
    **kwargs
) -> plt.Axes:
    """Create a Q-Q plot for p-value distribution analysis.

    Args:
        p_values: List of p-values to plot
        ax: Matplotlib axes (creates new if None)
        title: Plot title
        **kwargs: Additional arguments for scatter

    Returns:
        Matplotlib axes object
    """
    if ax is None:
        _, ax = plt.subplots()

    import scipy.stats as stats

    # Remove p-values of 0 (which would cause issues with log)
    p_values = [p for p in p_values if p > 0]

    if len(p_values) == 0:
        ax.text(0.5, 0.5, "No valid p-values", ha="center", va="center", transform=ax.transAxes)
        ax.set_title(title)
        return ax

    # Sort p-values
    p_sorted = np.sort(p_values)

    # Generate expected p-values under null hypothesis
    n = len(p_sorted)
    expected = np.array([(i + 0.5) / n for i in range(n)])

    # Create Q-Q plot
    ax.scatter(-np.log10(expected), -np.log10(p_sorted), alpha=0.6, **kwargs)

    # Add diagonal line (y = x)
    max_val = max(-np.log10(expected).max(), -np.log10(p_sorted).max())
    ax.plot([0, max_val], [0, max_val], 'r--', alpha=0.8, label='Expected under null')

    ax.set_xlabel("Expected -log₁₀(P)")
    ax.set_ylabel("Observed -log₁₀(P)")
    ax.set_title(title)
    ax.legend()

    return ax


def enrichment_plot(
    data: pd.DataFrame,
    x_col: str,
    y_col: str,
    *,
    p_threshold: float = 0.05,
    ax: plt.Axes | None = None,
    **kwargs
) -> plt.Axes:
    """Create an enrichment plot for pathway/gene set analysis.

    Args:
        data: DataFrame with enrichment results
        x_col: Column name for x-axis (usually enrichment ratio or odds ratio)
        y_col: Column name for y-axis (usually -log10(p-value))
        p_threshold: P-value threshold for significance
        ax: Matplotlib axes (creates new if None)
        **kwargs: Additional arguments for scatter

    Returns:
        Matplotlib axes object
    """
    if ax is None:
        _, ax = plt.subplots(figsize=(10, 8))

    # Create significance mask
    significant = data[y_col] > -np.log10(p_threshold)

    # Plot non-significant points
    if not significant.all():
        ax.scatter(data.loc[~significant, x_col], data.loc[~significant, y_col],
                  color='gray', alpha=0.5, label='Not significant', **kwargs)

    # Plot significant points
    if significant.any():
        ax.scatter(data.loc[significant, x_col], data.loc[significant, y_col],
                  color='red', alpha=0.8, label='Significant', **kwargs)

    # Add threshold line
    ax.axhline(-np.log10(p_threshold), color='black', linestyle='--', alpha=0.7)

    ax.set_xlabel(x_col)
    ax.set_ylabel(y_col)
    ax.set_title("Enrichment Analysis")
    ax.legend()

    return ax


def network_plot(
    nodes: list[str],
    edges: list[tuple[str, str]],
    *,
    node_sizes: list[float] | None = None,
    node_colors: list[str] | None = None,
    ax: plt.Axes | None = None,
    **kwargs
) -> plt.Axes:
    """Create a network graph visualization.

    Args:
        nodes: List of node names
        edges: List of (source, target) edge tuples
        node_sizes: Optional list of node sizes
        node_colors: Optional list of node colors
        ax: Matplotlib axes (creates new if None)
        **kwargs: Additional arguments for networkx.draw

    Returns:
        Matplotlib axes object
    """
    try:
        import networkx as nx
    except ImportError:
        # Fallback to simple text representation
        if ax is None:
            _, ax = plt.subplots()
        ax.text(0.5, 0.5, "NetworkX not available\nInstall with: uv pip install networkx",
               ha="center", va="center", transform=ax.transAxes)
        return ax

    if ax is None:
        _, ax = plt.subplots()

    # Create graph
    G = nx.Graph()
    G.add_nodes_from(nodes)
    G.add_edges_from(edges)

    # Set default node sizes
    if node_sizes is None:
        node_sizes = [300] * len(nodes)

    # Set default node colors
    if node_colors is None:
        node_colors = ['lightblue'] * len(nodes)

    # Draw network
    pos = nx.spring_layout(G, k=1, iterations=50)
    nx.draw(G, pos, ax=ax, node_size=node_sizes, node_color=node_colors,
            with_labels=True, font_size=8, **kwargs)

    ax.set_title("Network Graph")

    return ax
