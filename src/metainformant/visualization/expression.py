"""Gene expression analysis visualization functions.

This module provides specialized plotting functions for gene expression data
including expression heatmaps, pathway enrichment barplots, and differential
expression analysis visualizations.
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


def plot_expression_heatmap(
    data: pd.DataFrame,
    *,
    ax: Axes | None = None,
    output_path: str | Path | None = None,
    **kwargs
) -> Axes:
    """Create a gene expression heatmap.

    Args:
        data: DataFrame with genes as rows and samples as columns
        ax: Optional matplotlib axes to plot on. Creates new if None.
        output_path: Optional path to save the figure.
        **kwargs: Additional arguments passed to seaborn heatmap() or matplotlib imshow().

    Returns:
        matplotlib Axes object

    Raises:
        ValueError: If data is empty or contains non-numeric values
    """
    validation.validate_type(data, pd.DataFrame, "data")

    if data.empty:
        raise ValueError("Expression data DataFrame cannot be empty")

    # Check if data contains numeric values
    numeric_data = data.select_dtypes(include=[np.number])
    if numeric_data.empty:
        raise ValueError("Expression data must contain numeric columns")

    if ax is None:
        fig, ax = plt.subplots(figsize=kwargs.pop('figsize', (12, 8)))

    # Use seaborn if available, otherwise matplotlib
    if HAS_SEABORN:
        # Z-score normalize for better visualization
        z_scored = (numeric_data - numeric_data.mean()) / numeric_data.std()
        sns.heatmap(
            z_scored,
            ax=ax,
            cmap=kwargs.get('cmap', 'RdYlBu_r'),
            center=kwargs.get('center', 0),
            annot=kwargs.get('annot', False),
            fmt=kwargs.get('fmt', '.2f'),
            **kwargs
        )
        ax.set_title('Gene Expression Heatmap (Z-score normalized)')
    else:
        logger.warning("Seaborn not available, using basic heatmap")
        # Basic matplotlib heatmap
        im = ax.imshow(numeric_data.values, cmap=kwargs.get('cmap', 'viridis'), **kwargs)
        plt.colorbar(im, ax=ax)
        ax.set_title('Gene Expression Heatmap')

        # Set tick labels if not too many
        if len(data.index) <= 20:
            ax.set_yticks(range(len(data.index)))
            ax.set_yticklabels(data.index, fontsize=8)
        if len(data.columns) <= 20:
            ax.set_xticks(range(len(data.columns)))
            ax.set_xticklabels(data.columns, rotation=45, ha='right', fontsize=8)

    if output_path:
        output_path = paths.ensure_directory(Path(output_path).parent)
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        logger.info(f"Expression heatmap saved to {output_path}")

    return ax


def plot_enrichment_barplot(
    enrichment_results: pd.DataFrame,
    *,
    ax: Axes | None = None,
    output_path: str | Path | None = None,
    **kwargs
) -> Axes:
    """Create a pathway enrichment barplot.

    Args:
        enrichment_results: DataFrame with enrichment results containing
                           columns like 'Term', 'p.adjust', 'Count', 'Fold.Enrichment'
        ax: Optional matplotlib axes to plot on. Creates new if None.
        output_path: Optional path to save the figure.
        **kwargs: Additional arguments passed to matplotlib barh().

    Returns:
        matplotlib Axes object

    Raises:
        ValueError: If required columns are missing
    """
    validation.validate_type(enrichment_results, pd.DataFrame, "enrichment_results")

    if enrichment_results.empty:
        raise ValueError("Enrichment results DataFrame cannot be empty")

    # Check for required columns (flexible naming)
    term_col = None
    pval_col = None

    for col in enrichment_results.columns:
        col_lower = col.lower()
        if 'term' in col_lower or 'pathway' in col_lower or 'description' in col_lower:
            term_col = col
        elif 'padj' in col_lower or 'p.adjust' in col_lower or 'fdr' in col_lower:
            pval_col = col

    if term_col is None:
        raise ValueError("Enrichment results must contain a term/pathway column")
    if pval_col is None:
        raise ValueError("Enrichment results must contain an adjusted p-value column")

    if ax is None:
        fig, ax = plt.subplots(figsize=kwargs.pop('figsize', (10, 6)))

    # Sort by p-value and take top N
    sorted_results = enrichment_results.sort_values(pval_col).head(20)
    terms = sorted_results[term_col]
    neg_log_pvals = -np.log10(sorted_results[pval_col].clip(lower=1e-300))

    # Create horizontal bar plot
    y_pos = np.arange(len(terms))
    ax.barh(y_pos, neg_log_pvals, **kwargs)

    ax.set_yticks(y_pos)
    # Truncate long term names
    term_labels = [term[:50] + '...' if len(str(term)) > 50 else str(term) for term in terms]
    ax.set_yticklabels(term_labels, fontsize=8)
    ax.set_xlabel('-log₁₀(adjusted p-value)')
    ax.set_title('Pathway Enrichment Results')

    # Add count/size information if available
    count_col = None
    for col in enrichment_results.columns:
        if 'count' in col.lower() or 'size' in col.lower():
            count_col = col
            break

    if count_col and count_col in sorted_results.columns:
        # Add text annotations for counts
        for i, count in enumerate(sorted_results[count_col]):
            ax.text(neg_log_pvals.iloc[i] + 0.1, i, f'n={count}', va='center', fontsize=8)

    if output_path:
        output_path = paths.ensure_directory(Path(output_path).parent)
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        logger.info(f"Enrichment barplot saved to {output_path}")

    return ax


def plot_differential_expression(
    data: pd.DataFrame,
    *,
    ax: Axes | None = None,
    output_path: str | Path | None = None,
    **kwargs
) -> Axes:
    """Create a differential expression analysis plot.

    Args:
        data: DataFrame with differential expression results containing
              columns like 'log2FoldChange', 'padj', 'baseMean'
        ax: Optional matplotlib axes to plot on. Creates new if None.
        output_path: Optional path to save the figure.
        **kwargs: Additional arguments passed to matplotlib scatter().

    Returns:
        matplotlib Axes object

    Raises:
        ValueError: If required columns are missing
    """
    validation.validate_type(data, pd.DataFrame, "data")

    if data.empty:
        raise ValueError("Differential expression data DataFrame cannot be empty")

    # Check for required columns
    logfc_col = None
    pval_col = None
    mean_col = None

    for col in data.columns:
        col_lower = col.lower()
        if 'log2foldchange' in col_lower or 'logfc' in col_lower or 'log2fc' in col_lower:
            logfc_col = col
        elif 'padj' in col_lower or 'p.adjust' in col_lower or 'adj.p' in col_lower:
            pval_col = col
        elif 'basemean' in col_lower or 'mean' in col_lower or 'avg' in col_lower:
            mean_col = col

    if logfc_col is None:
        raise ValueError("Differential expression data must contain a log fold change column")
    if pval_col is None:
        raise ValueError("Differential expression data must contain an adjusted p-value column")

    if ax is None:
        fig, ax = plt.subplots(figsize=kwargs.pop('figsize', (10, 8)))

    # Prepare data
    plot_data = data.copy()
    plot_data['logP'] = -np.log10(plot_data[pval_col].clip(lower=1e-300))

    # Color points based on significance
    colors = []
    for _, row in plot_data.iterrows():
        if row['logP'] >= -np.log10(0.05) and abs(row[logfc_col]) >= 1:
            colors.append('red')  # Significant and large fold change
        elif row['logP'] >= -np.log10(0.05):
            colors.append('orange')  # Significant
        elif abs(row[logfc_col]) >= 1:
            colors.append('blue')  # Large fold change
        else:
            colors.append('gray')  # Not significant

    # Plot with size based on expression level if available
    sizes = kwargs.pop('s', 20)
    if mean_col and mean_col in plot_data.columns:
        # Scale point sizes by expression level
        mean_vals = plot_data[mean_col]
        sizes = 20 + (mean_vals - mean_vals.min()) / (mean_vals.max() - mean_vals.min()) * 80

    ax.scatter(
        plot_data[logfc_col],
        plot_data['logP'],
        c=colors,
        s=sizes,
        alpha=kwargs.get('alpha', 0.6),
        **kwargs
    )

    # Add threshold lines
    ax.axhline(y=-np.log10(0.05), color='black', linestyle='--', alpha=0.7)
    ax.axvline(x=1, color='black', linestyle='--', alpha=0.7)
    ax.axvline(x=-1, color='black', linestyle='--', alpha=0.7)

    ax.set_xlabel('log₂(Fold Change)')
    ax.set_ylabel('-log₁₀(adjusted p-value)')
    ax.set_title('Differential Expression Analysis')

    # Add legend
    legend_elements = [
        plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='red',
                  markersize=8, label='Significant & |FC| ≥ 2'),
        plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='orange',
                  markersize=8, label='Significant'),
        plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='blue',
                  markersize=8, label='|FC| ≥ 2'),
        plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='gray',
                  markersize=8, label='Not significant')
    ]
    ax.legend(handles=legend_elements, loc='best', fontsize=8)

    if output_path:
        output_path = paths.ensure_directory(Path(output_path).parent)
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        logger.info(f"Differential expression plot saved to {output_path}")

    return ax
