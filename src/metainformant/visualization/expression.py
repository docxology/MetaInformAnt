"""Expression analysis visualization functions.

This module provides visualization functions for expression data including
expression heatmaps, enrichment plots, gene expression plots, differential
expression plots, and log fold change visualizations.
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

from .basic import heatmap


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

    Example:
        >>> from metainformant.visualization import expression_heatmap
        >>> import pandas as pd
        >>> import numpy as np
        >>> data = pd.DataFrame(np.random.random((10, 5)))
        >>> ax = expression_heatmap(data)
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

    Example:
        >>> from metainformant.visualization import enrichment_plot
        >>> import pandas as pd
        >>> data = pd.DataFrame({
        ...     'enrichment_ratio': [1.5, 2.0, 0.8, 1.2],
        ...     'pvalue': [0.001, 0.01, 0.5, 0.1]
        ... })
        >>> data['neg_log10_p'] = -np.log10(data['pvalue'])
        >>> ax = enrichment_plot(data, 'enrichment_ratio', 'neg_log10_p')
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


def gene_expression_plot(
    gene_name: str,
    expression_data: pd.DataFrame,
    sample_groups: Sequence[str] | None = None,
    *,
    ax: plt.Axes | None = None,
    title: str | None = None,
    **kwargs
) -> plt.Axes:
    """Plot expression levels for a single gene across samples.

    Args:
        gene_name: Name of the gene to plot
        expression_data: DataFrame with genes as rows and samples as columns
        sample_groups: Optional group labels for samples
        ax: Matplotlib axes (creates new if None)
        title: Plot title
        **kwargs: Additional arguments for plotting

    Returns:
        Matplotlib axes object

    Example:
        >>> from metainformant.visualization import gene_expression_plot
        >>> import pandas as pd
        >>> data = pd.DataFrame({
        ...     'sample1': [10, 20, 30],
        ...     'sample2': [15, 25, 35]
        ... }, index=['gene1', 'gene2', 'gene3'])
        >>> ax = gene_expression_plot('gene1', data)
    """
    if ax is None:
        _, ax = plt.subplots()

    if gene_name not in expression_data.index:
        raise ValueError(f"Gene {gene_name} not found in expression data")

    gene_data = expression_data.loc[gene_name]

    if sample_groups:
        # Group by sample groups
        groups = list(set(sample_groups))
        group_data = {group: [] for group in groups}
        for i, group in enumerate(sample_groups):
            if i < len(gene_data):
                group_data[group].append(gene_data.iloc[i])

        # Create box plot by group
        ax.boxplot([group_data[group] for group in groups], labels=groups, **kwargs)
        ax.set_ylabel("Expression Level")
        ax.set_title(title or f"Expression of {gene_name}")
    else:
        # Plot as bar chart
        ax.bar(range(len(gene_data)), gene_data.values, **kwargs)
        ax.set_xticks(range(len(gene_data)))
        ax.set_xticklabels(gene_data.index, rotation=45, ha='right')
        ax.set_ylabel("Expression Level")
        ax.set_title(title or f"Expression of {gene_name}")

    return ax


def differential_expression_plot(
    data: pd.DataFrame,
    gene_col: str = "gene",
    log2fc_col: str = "log2fc",
    pvalue_col: str = "pvalue",
    *,
    top_n: int = 20,
    ax: plt.Axes | None = None,
    **kwargs
) -> plt.Axes:
    """Plot top differentially expressed genes.

    Args:
        data: DataFrame with differential expression results
        gene_col: Column name for gene names
        log2fc_col: Column name for log2 fold change
        pvalue_col: Column name for p-values
        top_n: Number of top genes to plot
        ax: Matplotlib axes (creates new if None)
        **kwargs: Additional arguments

    Returns:
        Matplotlib axes object

    Example:
        >>> from metainformant.visualization import differential_expression_plot
        >>> import pandas as pd
        >>> data = pd.DataFrame({
        ...     'gene': ['gene1', 'gene2', 'gene3'],
        ...     'log2fc': [2.0, -1.5, 1.0],
        ...     'pvalue': [0.001, 0.01, 0.05]
        ... })
        >>> ax = differential_expression_plot(data, top_n=3)
    """
    if ax is None:
        _, ax = plt.subplots(figsize=(10, 6))

    # Sort by p-value and take top N
    data_sorted = data.sort_values(pvalue_col).head(top_n)

    # Create bar plot
    colors = ['red' if fc > 0 else 'blue' for fc in data_sorted[log2fc_col]]
    ax.barh(range(len(data_sorted)), data_sorted[log2fc_col], color=colors, **kwargs)

    ax.set_yticks(range(len(data_sorted)))
    ax.set_yticklabels(data_sorted[gene_col])
    ax.set_xlabel("Log2 Fold Change")
    ax.set_title(f"Top {top_n} Differentially Expressed Genes")
    ax.axvline(x=0, color='black', linestyle='-', linewidth=0.5)

    return ax


def log_fold_change_plot(
    data: pd.DataFrame,
    log2fc_col: str = "log2fc",
    group_col: str | None = None,
    *,
    ax: plt.Axes | None = None,
    title: str = "Log2 Fold Change Distribution",
    **kwargs
) -> plt.Axes:
    """Plot distribution of log fold changes.

    Args:
        data: DataFrame with log fold change values
        log2fc_col: Column name for log2 fold change
        group_col: Optional column for grouping
        ax: Matplotlib axes (creates new if None)
        title: Plot title
        **kwargs: Additional arguments

    Returns:
        Matplotlib axes object

    Example:
        >>> from metainformant.visualization import log_fold_change_plot
        >>> import pandas as pd
        >>> import numpy as np
        >>> data = pd.DataFrame({
        ...     'log2fc': np.random.normal(0, 1, 100)
        ... })
        >>> ax = log_fold_change_plot(data)
    """
    if ax is None:
        _, ax = plt.subplots()

    if group_col and group_col in data.columns:
        # Grouped histogram
        groups = data[group_col].unique()
        for group in groups:
            group_data = data[data[group_col] == group][log2fc_col]
            ax.hist(group_data, alpha=0.5, label=str(group), **kwargs)
        ax.legend()
    else:
        # Single histogram
        ax.hist(data[log2fc_col], **kwargs)

    ax.axvline(x=0, color='red', linestyle='--', alpha=0.7, label='No change')
    ax.set_xlabel("Log2 Fold Change")
    ax.set_ylabel("Frequency")
    ax.set_title(title)

    return ax

