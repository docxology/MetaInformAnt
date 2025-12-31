"""Genomic data visualization functions for GWAS and genomic analysis.

This module provides specialized plotting functions for genomic data including
Manhattan plots, volcano plots, regional association plots, chromosome ideograms,
coverage plots, and variant visualizations.
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


def manhattan_plot(
    data: pd.DataFrame,
    *,
    chr_col: str = "CHR",
    pos_col: str = "BP",
    pval_col: str = "P",
    ax: Axes | None = None,
    output_path: str | Path | None = None,
    **kwargs
) -> Axes:
    """Create a Manhattan plot for GWAS results.

    Args:
        data: DataFrame with chromosome, position, and p-value columns
        chr_col: Column name for chromosome
        pos_col: Column name for position
        pval_col: Column name for p-value
        ax: Optional matplotlib axes to plot on. Creates new if None.
        output_path: Optional path to save the figure.
        **kwargs: Additional arguments passed to matplotlib scatter().

    Returns:
        matplotlib Axes object

    Raises:
        ValueError: If required columns are missing or data is invalid
    """
    validation.validate_type(data, pd.DataFrame, "data")

    required_cols = [chr_col, pos_col, pval_col]
    missing_cols = [col for col in required_cols if col not in data.columns]
    if missing_cols:
        raise ValueError(f"Missing required columns: {missing_cols}")

    if ax is None:
        fig, ax = plt.subplots(figsize=kwargs.pop('figsize', (12, 6)))

    # Prepare data
    plot_data = data.copy()
    plot_data['logP'] = -np.log10(plot_data[pval_col].clip(lower=1e-300))  # Avoid log(0)

    # Calculate cumulative positions for chromosomes
    chromosomes = sorted(plot_data[chr_col].unique())
    chr_starts = {}
    cumulative_pos = 0

    for chr_num in chromosomes:
        chr_mask = plot_data[chr_col] == chr_num
        chr_starts[chr_num] = cumulative_pos
        plot_data.loc[chr_mask, 'cumulative_pos'] = plot_data.loc[chr_mask, pos_col] + cumulative_pos
        cumulative_pos += plot_data.loc[chr_mask, pos_col].max() + 1e7  # Add gap between chromosomes

    # Plot points colored by chromosome
    colors = plt.cm.tab10(np.linspace(0, 1, len(chromosomes)))
    for i, chr_num in enumerate(chromosomes):
        chr_mask = plot_data[chr_col] == chr_num
        ax.scatter(
            plot_data.loc[chr_mask, 'cumulative_pos'],
            plot_data.loc[chr_mask, 'logP'],
            c=[colors[i]],
            s=kwargs.get('s', 1),
            alpha=kwargs.get('alpha', 0.8),
            **kwargs
        )

    # Add chromosome labels at centers
    chr_centers = []
    chr_labels = []
    for chr_num in chromosomes:
        chr_mask = plot_data[chr_col] == chr_num
        if chr_mask.any():
            center_pos = (plot_data.loc[chr_mask, 'cumulative_pos'].min() +
                         plot_data.loc[chr_mask, 'cumulative_pos'].max()) / 2
            chr_centers.append(center_pos)
            chr_labels.append(str(chr_num))

    ax.set_xticks(chr_centers)
    ax.set_xticklabels(chr_labels)
    ax.set_xlabel('Chromosome')
    ax.set_ylabel('-log₁₀(p-value)')
    ax.set_title('Manhattan Plot')
    ax.axhline(y=-np.log10(5e-8), color='red', linestyle='--', alpha=0.7, label='Genome-wide significance')

    # Add significance threshold legend if line is visible
    if plot_data['logP'].max() > -np.log10(5e-8):
        ax.legend()

    if output_path:
        output_path = paths.ensure_directory(Path(output_path).parent)
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        logger.info(f"Manhattan plot saved to {output_path}")

    return ax


def volcano_plot(
    data: pd.DataFrame,
    *,
    log2fc_col: str = "log2FoldChange",
    pval_col: str = "padj",
    ax: Axes | None = None,
    output_path: str | Path | None = None,
    **kwargs
) -> Axes:
    """Create a volcano plot for differential expression analysis.

    Args:
        data: DataFrame with log fold change and p-value columns
        log2fc_col: Column name for log2 fold change
        pval_col: Column name for adjusted p-value
        ax: Optional matplotlib axes to plot on. Creates new if None.
        output_path: Optional path to save the figure.
        **kwargs: Additional arguments passed to matplotlib scatter().

    Returns:
        matplotlib Axes object

    Raises:
        ValueError: If required columns are missing
    """
    validation.validate_type(data, pd.DataFrame, "data")

    required_cols = [log2fc_col, pval_col]
    missing_cols = [col for col in required_cols if col not in data.columns]
    if missing_cols:
        raise ValueError(f"Missing required columns: {missing_cols}")

    if ax is None:
        fig, ax = plt.subplots(figsize=kwargs.pop('figsize', (8, 6)))

    # Prepare data
    plot_data = data.copy()
    plot_data['logP'] = -np.log10(plot_data[pval_col].clip(lower=1e-300))

    # Color points based on significance and fold change
    colors = []
    for _, row in plot_data.iterrows():
        if row['logP'] >= -np.log10(0.05) and abs(row[log2fc_col]) >= 1:
            colors.append('red')  # Significant and large fold change
        elif row['logP'] >= -np.log10(0.05):
            colors.append('orange')  # Significant
        elif abs(row[log2fc_col]) >= 1:
            colors.append('blue')  # Large fold change
        else:
            colors.append('gray')  # Not significant

    ax.scatter(
        plot_data[log2fc_col],
        plot_data['logP'],
        c=colors,
        s=kwargs.get('s', 20),
        alpha=kwargs.get('alpha', 0.6),
        **kwargs
    )

    # Add threshold lines
    ax.axhline(y=-np.log10(0.05), color='black', linestyle='--', alpha=0.7)
    ax.axvline(x=1, color='black', linestyle='--', alpha=0.7)
    ax.axvline(x=-1, color='black', linestyle='--', alpha=0.7)

    ax.set_xlabel('log₂(Fold Change)')
    ax.set_ylabel('-log₁₀(adjusted p-value)')
    ax.set_title('Volcano Plot')

    if output_path:
        output_path = paths.ensure_directory(Path(output_path).parent)
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        logger.info(f"Volcano plot saved to {output_path}")

    return ax


def regional_plot(
    data: pd.DataFrame,
    *,
    chr: str,
    start: int,
    end: int,
    ax: Axes | None = None,
    output_path: str | Path | None = None,
    **kwargs
) -> Axes:
    """Create a regional association plot for a specific genomic region.

    Args:
        data: DataFrame with genomic data
        chr: Chromosome name
        start: Start position
        end: End position
        ax: Optional matplotlib axes to plot on. Creates new if None.
        output_path: Optional path to save the figure.
        **kwargs: Additional arguments passed to matplotlib scatter().

    Returns:
        matplotlib Axes object

    Raises:
        ValueError: If region parameters are invalid
    """
    validation.validate_type(data, pd.DataFrame, "data")

    if start >= end:
        raise ValueError("Start position must be less than end position")

    if ax is None:
        fig, ax = plt.subplots(figsize=kwargs.pop('figsize', (10, 6)))

    # Filter data for the specified region
    region_mask = (
        (data['CHR'].astype(str) == str(chr)) &
        (data['BP'] >= start) &
        (data['BP'] <= end)
    )
    region_data = data[region_mask].copy()

    if region_data.empty:
        logger.warning(f"No data found in region {chr}:{start}-{end}")
        ax.text(0.5, 0.5, f'No data in region {chr}:{start}-{end}',
                transform=ax.transAxes, ha='center', va='center')
        ax.set_xlim(start, end)
        ax.set_xlabel(f'Position on chromosome {chr}')
        ax.set_ylabel('-log₁₀(p-value)')
        ax.set_title(f'Regional Plot - {chr}:{start}-{end}')
        return ax

    # Calculate relative positions
    region_data['relative_pos'] = region_data['BP'] - start
    region_data['logP'] = -np.log10(region_data['P'].clip(lower=1e-300))

    # Plot association results
    ax.scatter(
        region_data['relative_pos'],
        region_data['logP'],
        s=kwargs.get('s', 30),
        alpha=kwargs.get('alpha', 0.7),
        **kwargs
    )

    # Add recombination rate if available (placeholder)
    # This would require additional data sources

    ax.set_xlabel(f'Position on chromosome {chr} (bp from {start:,})')
    ax.set_ylabel('-log₁₀(p-value)')
    ax.set_title(f'Regional Association Plot - {chr}:{start:,}-{end:,}')
    ax.axhline(y=-np.log10(5e-8), color='red', linestyle='--', alpha=0.7)

    if output_path:
        output_path = paths.ensure_directory(Path(output_path).parent)
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        logger.info(f"Regional plot saved to {output_path}")

    return ax


def circular_manhattan_plot(
    data: pd.DataFrame,
    *,
    ax: Axes | None = None,
    output_path: str | Path | None = None,
    **kwargs
) -> Axes:
    """Create a circular Manhattan plot.

    Args:
        data: DataFrame with chromosome, position, and p-value columns
        ax: Optional matplotlib axes to plot on. Creates new if None.
        output_path: Optional path to save the figure.
        **kwargs: Additional arguments passed to matplotlib scatter().

    Returns:
        matplotlib Axes object
    """
    validation.validate_type(data, pd.DataFrame, "data")

    if ax is None:
        fig, ax = plt.subplots(figsize=kwargs.pop('figsize', (8, 8)), subplot_kw={'projection': 'polar'})

    # Prepare data similar to linear Manhattan plot
    plot_data = data.copy()
    plot_data['logP'] = -np.log10(plot_data['P'].clip(lower=1e-300))

    # Calculate angles for chromosomes (divide circle evenly)
    chromosomes = sorted(plot_data['CHR'].unique())
    chr_angles = np.linspace(0, 2*np.pi, len(chromosomes) + 1)[:-1]

    # Plot each chromosome
    for i, (chr_num, angle) in enumerate(zip(chromosomes, chr_angles)):
        chr_mask = plot_data['CHR'] == chr_num
        chr_data = plot_data[chr_mask]

        # Normalize positions within chromosome to angle range
        chr_start = chr_data['BP'].min()
        chr_end = chr_data['BP'].max()
        chr_range = chr_end - chr_start

        if chr_range > 0:
            relative_angles = angle + (chr_data['BP'] - chr_start) / chr_range * (2*np.pi / len(chromosomes))

            ax.scatter(
                relative_angles,
                chr_data['logP'],
                s=kwargs.get('s', 2),
                alpha=kwargs.get('alpha', 0.7),
                label=f'Chr {chr_num}' if i < 5 else "",  # Label first few chromosomes
                **kwargs
            )

    # Add chromosome labels
    for angle, chr_num in zip(chr_angles, chromosomes):
        ax.text(angle, ax.get_ylim()[1] * 1.1, f'{chr_num}',
                ha='center', va='bottom', fontsize=8)

    ax.set_title('Circular Manhattan Plot')
    ax.set_rlabel_position(90)
    ax.grid(True, alpha=0.3)

    if output_path:
        output_path = paths.ensure_directory(Path(output_path).parent)
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        logger.info(f"Circular Manhattan plot saved to {output_path}")

    return ax


def chromosome_ideogram(
    *,
    ax: Axes | None = None,
    output_path: str | Path | None = None,
    **kwargs
) -> Axes:
    """Create a chromosome ideogram visualization.

    Args:
        ax: Optional matplotlib axes to plot on. Creates new if None.
        output_path: Optional path to save the figure.
        **kwargs: Additional arguments passed to plotting functions.

    Returns:
        matplotlib Axes object
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=kwargs.pop('figsize', (12, 4)))

    # Simplified chromosome ideogram (would need real cytogenetic data)
    chromosomes = list(range(1, 23)) + ['X', 'Y']
    chr_lengths = {
        1: 248956422, 2: 242193529, 3: 198295559, 4: 190214555, 5: 181538259,
        6: 170805979, 7: 159345973, 8: 145138636, 9: 138394717, 10: 133797422,
        11: 135086622, 12: 133275309, 13: 114364328, 14: 107043718, 15: 101991189,
        16: 90338345, 17: 83257441, 18: 80373285, 19: 58617616, 20: 64444167,
        21: 46709983, 22: 50818468, 'X': 156040895, 'Y': 57227415
    }

    # Plot chromosomes as horizontal bars
    y_positions = range(len(chromosomes))
    max_length = max(chr_lengths.values())

    for i, chr_name in enumerate(chromosomes):
        length = chr_lengths.get(chr_name, chr_lengths[1])  # fallback
        ax.barh(i, length / max_length, height=0.8, alpha=0.7,
                color='lightblue', edgecolor='black')

    ax.set_yticks(y_positions)
    ax.set_yticklabels([f'{c}' for c in chromosomes])
    ax.set_xlabel('Relative Chromosome Length')
    ax.set_title('Human Chromosome Ideogram')
    ax.grid(True, alpha=0.3)

    if output_path:
        output_path = paths.ensure_directory(Path(output_path).parent)
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        logger.info(f"Chromosome ideogram saved to {output_path}")

    return ax


def coverage_plot(
    coverage: np.ndarray,
    positions: np.ndarray,
    *,
    ax: Axes | None = None,
    output_path: str | Path | None = None,
    **kwargs
) -> Axes:
    """Create a sequencing coverage plot.

    Args:
        coverage: Array of coverage depths
        positions: Array of genomic positions
        ax: Optional matplotlib axes to plot on. Creates new if None.
        output_path: Optional path to save the figure.
        **kwargs: Additional arguments passed to matplotlib plot().

    Returns:
        matplotlib Axes object

    Raises:
        ValueError: If arrays have mismatched lengths
    """
    validation.validate_type(coverage, np.ndarray, "coverage")
    validation.validate_type(positions, np.ndarray, "positions")

    if len(coverage) != len(positions):
        raise ValueError("Coverage and positions arrays must have same length")

    if ax is None:
        fig, ax = plt.subplots(figsize=kwargs.pop('figsize', (12, 4)))

    ax.plot(positions, coverage, **kwargs)
    ax.fill_between(positions, coverage, alpha=0.3, **kwargs)

    ax.set_xlabel('Genomic Position')
    ax.set_ylabel('Coverage Depth')
    ax.set_title('Sequencing Coverage')
    ax.grid(True, alpha=0.3)

    if output_path:
        output_path = paths.ensure_directory(Path(output_path).parent)
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        logger.info(f"Coverage plot saved to {output_path}")

    return ax


def variant_plot(
    variants: pd.DataFrame,
    *,
    ax: Axes | None = None,
    output_path: str | Path | None = None,
    **kwargs
) -> Axes:
    """Create a variant visualization plot.

    Args:
        variants: DataFrame with variant information
        ax: Optional matplotlib axes to plot on. Creates new if None.
        output_path: Optional path to save the figure.
        **kwargs: Additional arguments passed to matplotlib scatter().

    Returns:
        matplotlib Axes object

    Raises:
        ValueError: If required columns are missing
    """
    validation.validate_type(variants, pd.DataFrame, "variants")

    if ax is None:
        fig, ax = plt.subplots(figsize=kwargs.pop('figsize', (12, 6)))

    # Assume variants has POS and some value column
    if 'POS' not in variants.columns:
        raise ValueError("Variants DataFrame must contain 'POS' column")

    # Simple variant visualization - plot positions
    positions = variants['POS'].values
    y_values = np.ones(len(positions))  # Default y-position

    # Add some jitter for visibility
    y_values += np.random.uniform(-0.1, 0.1, len(positions))

    ax.scatter(positions, y_values, s=kwargs.get('s', 10), alpha=0.7, **kwargs)

    ax.set_xlabel('Genomic Position')
    ax.set_ylabel('Variants')
    ax.set_title('Variant Positions')
    ax.set_yticks([])  # Hide y-axis ticks

    if output_path:
        output_path = paths.ensure_directory(Path(output_path).parent)
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        logger.info(f"Variant plot saved to {output_path}")

    return ax
