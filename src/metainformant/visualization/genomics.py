"""Genomic visualization functions for GWAS and sequence analysis.

This module provides genomic-specific plotting functions including Manhattan plots,
volcano plots, regional plots, circular Manhattan plots, chromosome ideograms,
coverage plots, and variant visualizations.
"""

from __future__ import annotations

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# Use non-interactive backend by default for tests/headless
matplotlib.use("Agg", force=True)


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
        y_col: Column name for y-axis (usually -log10(p-value), or raw p-values)
        chromosome_col: Column name for chromosome information
        p_threshold: P-value threshold for significance
        highlight_color: Color for significant points
        ax: Matplotlib axes (creates new if None)
        **kwargs: Additional arguments for scatter

    Returns:
        Matplotlib axes object

    Example:
        >>> from metainformant.visualization import manhattan_plot
        >>> import pandas as pd
        >>> data = pd.DataFrame({
        ...     'position': [1000, 2000, 3000],
        ...     'pvalue': [1e-6, 1e-9, 0.01],
        ...     'chromosome': ['chr1', 'chr1', 'chr2']
        ... })
        >>> data['neg_log10_p'] = -np.log10(data['pvalue'])
        >>> ax = manhattan_plot(data, 'position', 'neg_log10_p', 'chromosome')
    """
    if ax is None:
        _, ax = plt.subplots(figsize=(12, 6))

    # Detect if y_col contains raw p-values or -log10 transformed values
    y_values = data[y_col].values
    max_y = np.max(y_values[y_values > 0]) if np.any(y_values > 0) else 1.0
    
    if max_y <= 1.0:
        # Raw p-values: transform to -log10
        neg_log10_p = -np.log10(np.clip(y_values, 1e-300, 1.0))
        y_label = f"-log₁₀({y_col})"
    else:
        # Already -log10 transformed
        neg_log10_p = y_values
        y_label = "-log₁₀(P-value)"

    # Create significance mask (boolean array)
    significant = neg_log10_p > -np.log10(p_threshold)

    # Plot non-significant points
    if not significant.all():
        not_sig_mask = ~significant
        ax.scatter(data.loc[not_sig_mask, x_col].values, neg_log10_p[not_sig_mask],
                  color='gray', alpha=0.6, s=10, **kwargs)

    # Plot significant points
    if significant.any():
        ax.scatter(data.loc[significant, x_col].values, neg_log10_p[significant],
                  color=highlight_color, alpha=0.8, s=20, **kwargs)

    # Add threshold line
    ax.axhline(-np.log10(p_threshold), color='black', linestyle='--', alpha=0.7, label=f'P = {p_threshold}')

    # Group by chromosome for alternating colors (if chromosome_col exists)
    if chromosome_col in data.columns:
        chromosomes = data[chromosome_col].unique()
        colors = plt.cm.tab20(np.linspace(0, 1, len(chromosomes)))

        for i, chrom in enumerate(chromosomes):
            chrom_mask = data[chromosome_col] == chrom
            chrom_data = data[chrom_mask]
            chrom_significant = significant[chrom_mask]
            chrom_not_sig = ~chrom_significant
            if chrom_not_sig.any():
                ax.scatter(chrom_data.loc[chrom_not_sig, x_col].values,
                          neg_log10_p[chrom_mask][chrom_not_sig],
                          color=colors[i], alpha=0.6, s=10)

    ax.set_xlabel("Genomic Position")
    ax.set_ylabel(y_label)
    ax.set_title("Manhattan Plot")
    ax.legend()

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

    Example:
        >>> from metainformant.visualization import volcano_plot
        >>> import pandas as pd
        >>> data = pd.DataFrame({
        ...     'log2fc': [-2, 1, 0.5, -1.5],
        ...     'pvalue': [0.001, 0.01, 0.5, 0.001]
        ... })
        >>> data['neg_log10_p'] = -np.log10(data['pvalue'])
        >>> ax = volcano_plot(data, 'log2fc', 'neg_log10_p')
    """
    if ax is None:
        _, ax = plt.subplots()

    # Detect if y_col contains raw p-values or -log10 transformed values
    # If max value is <= 1, assume raw p-values; otherwise assume -log10 transformed
    y_values = data[y_col].values
    max_y = np.max(y_values[y_values > 0]) if np.any(y_values > 0) else 1.0
    
    if max_y <= 1.0:
        # Raw p-values: transform to -log10
        neg_log10_p = -np.log10(np.clip(y_values, 1e-300, 1.0))  # Clip to avoid log(0)
        significance_threshold = -np.log10(p_threshold)
    else:
        # Already -log10 transformed
        neg_log10_p = y_values
        significance_threshold = -np.log10(p_threshold)

    # Create significance mask
    significant = (neg_log10_p > significance_threshold) & (abs(data[x_col]) > fc_threshold)

    # Use transformed values for plotting
    if max_y <= 1.0:
        # Plot using transformed values
        plot_y = neg_log10_p
        y_label = f"-log₁₀({y_col})"
    else:
        # Already transformed
        plot_y = data[y_col].values
        y_label = y_col

    # Plot non-significant points
    ax.scatter(data.loc[~significant, x_col].values, plot_y[~significant],
              color='gray', alpha=0.5, label='Not significant', **kwargs)

    # Plot significant points
    ax.scatter(data.loc[significant, x_col].values, plot_y[significant],
              color='red', alpha=0.8, label='Significant', **kwargs)

    # Add threshold lines
    ax.axvline(-fc_threshold, color='black', linestyle='--', alpha=0.7)
    ax.axvline(fc_threshold, color='black', linestyle='--', alpha=0.7)
    ax.axhline(significance_threshold, color='black', linestyle='--', alpha=0.7)

    ax.set_xlabel(x_col)
    ax.set_ylabel(y_label)
    ax.legend()

    return ax


def regional_plot(
    data: pd.DataFrame,
    chromosome: str,
    start: int,
    end: int,
    x_col: str = "position",
    y_col: str = "neg_log10_p",
    *,
    p_threshold: float = 5e-8,
    ax: plt.Axes | None = None,
    **kwargs
) -> plt.Axes:
    """Create a regional plot for a specific genomic region.

    Args:
        data: DataFrame with genomic positions and p-values
        chromosome: Chromosome name
        start: Start position
        end: End position
        x_col: Column name for x-axis (genomic position)
        y_col: Column name for y-axis (usually -log10(p-value))
        p_threshold: P-value threshold for significance
        ax: Matplotlib axes (creates new if None)
        **kwargs: Additional arguments for scatter

    Returns:
        Matplotlib axes object

    Example:
        >>> from metainformant.visualization import regional_plot
        >>> import pandas as pd
        >>> data = pd.DataFrame({
        ...     'position': range(1000000, 1001000, 100),
        ...     'neg_log10_p': np.random.uniform(0, 5, 10)
        ... })
        >>> ax = regional_plot(data, 'chr1', 1000000, 1001000)
    """
    if ax is None:
        _, ax = plt.subplots(figsize=(10, 6))

    # Filter to region
    region_data = data[(data[x_col] >= start) & (data[x_col] <= end)].copy()

    if len(region_data) == 0:
        ax.text(0.5, 0.5, "No data in region", ha="center", va="center", transform=ax.transAxes)
        return ax

    # Create significance mask
    significant = region_data[y_col] > -np.log10(p_threshold)

    # Plot non-significant points
    if not significant.all():
        ax.scatter(region_data.loc[~significant, x_col], region_data.loc[~significant, y_col],
                  color='gray', alpha=0.6, s=20, **kwargs)

    # Plot significant points
    if significant.any():
        ax.scatter(region_data.loc[significant, x_col], region_data.loc[significant, y_col],
                  color='red', alpha=0.8, s=30, **kwargs)

    # Add threshold line
    ax.axhline(-np.log10(p_threshold), color='black', linestyle='--', alpha=0.7, label=f'P = {p_threshold}')

    ax.set_xlabel(f"Position on {chromosome}")
    ax.set_ylabel("-log₁₀(P-value)")
    ax.set_title(f"Regional Plot: {chromosome}:{start}-{end}")
    ax.legend()

    return ax


def circular_manhattan_plot(
    data: pd.DataFrame,
    x_col: str,
    y_col: str,
    chromosome_col: str = "chromosome",
    *,
    p_threshold: float = 5e-8,
    ax: plt.Axes | None = None,
    **kwargs
) -> plt.Axes:
    """Create a circular Manhattan plot for genome-wide visualization.

    Args:
        data: DataFrame with genomic positions and p-values
        x_col: Column name for x-axis (genomic position)
        y_col: Column name for y-axis (usually -log10(p-value))
        chromosome_col: Column name for chromosome information
        p_threshold: P-value threshold for significance
        ax: Matplotlib axes (creates new if None)
        **kwargs: Additional arguments

    Returns:
        Matplotlib axes object

    Example:
        >>> from metainformant.visualization import circular_manhattan_plot
        >>> import pandas as pd
        >>> data = pd.DataFrame({
        ...     'position': [1000, 2000, 3000],
        ...     'pvalue': [1e-6, 1e-9, 0.01],
        ...     'chromosome': ['chr1', 'chr1', 'chr2']
        ... })
        >>> data['neg_log10_p'] = -np.log10(data['pvalue'])
        >>> ax = circular_manhattan_plot(data, 'position', 'neg_log10_p', 'chromosome')
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=(10, 10), subplot_kw=dict(projection='polar'))

    # Group by chromosome
    chromosomes = sorted(data[chromosome_col].unique())
    n_chromosomes = len(chromosomes)
    
    # Calculate angles for each chromosome
    angle_per_chrom = 2 * np.pi / n_chromosomes
    
    # Plot each chromosome
    for i, chrom in enumerate(chromosomes):
        chrom_data = data[data[chromosome_col] == chrom]
        if len(chrom_data) == 0:
            continue
            
        # Normalize positions to 0-1 within chromosome
        positions = chrom_data[x_col].values
        pos_min, pos_max = positions.min(), positions.max()
        if pos_max > pos_min:
            normalized_pos = (positions - pos_min) / (pos_max - pos_min)
        else:
            normalized_pos = np.zeros(len(positions))
        
        # Calculate angles
        theta = i * angle_per_chrom + normalized_pos * angle_per_chrom
        r = chrom_data[y_col].values
        
        # Create significance mask
        significant = r > -np.log10(p_threshold)
        
        # Plot non-significant points
        if not significant.all():
            ax.scatter(theta[~significant], r[~significant], 
                      color='gray', alpha=0.6, s=10, **kwargs)
        
        # Plot significant points
        if significant.any():
            ax.scatter(theta[significant], r[significant],
                      color='red', alpha=0.8, s=20, **kwargs)
    
    # Set labels
    ax.set_xticks([i * angle_per_chrom + angle_per_chrom / 2 for i in range(n_chromosomes)])
    ax.set_xticklabels(chromosomes)
    ax.set_ylabel("-log₁₀(P-value)", labelpad=20)
    ax.set_title("Circular Manhattan Plot", pad=20)
    
    return ax


def chromosome_ideogram(
    chromosomes: list[str],
    positions: list[int],
    values: list[float],
    *,
    threshold: float | None = None,
    ax: plt.Axes | None = None,
    **kwargs
) -> plt.Axes:
    """Create a chromosome ideogram with marked positions.

    Args:
        chromosomes: List of chromosome names
        positions: List of positions on chromosomes
        values: Values to visualize (e.g., p-values, effect sizes)
        threshold: Optional threshold for highlighting
        ax: Matplotlib axes (creates new if None)
        **kwargs: Additional arguments

    Returns:
        Matplotlib axes object

    Example:
        >>> from metainformant.visualization import chromosome_ideogram
        >>> ax = chromosome_ideogram(
        ...     ['chr1', 'chr2'],
        ...     [1000000, 2000000],
        ...     [5.0, 6.0],
        ...     threshold=5.0
        ... )
    """
    if ax is None:
        _, ax = plt.subplots(figsize=(12, 8))

    # Group by chromosome
    chrom_dict = {}
    for chrom, pos, val in zip(chromosomes, positions, values):
        if chrom not in chrom_dict:
            chrom_dict[chrom] = {'positions': [], 'values': []}
        chrom_dict[chrom]['positions'].append(pos)
        chrom_dict[chrom]['values'].append(val)

    # Draw chromosomes
    y_pos = 0
    chrom_positions = {}
    for chrom in sorted(chrom_dict.keys()):
        chrom_positions[chrom] = y_pos
        # Draw chromosome rectangle
        ax.add_patch(plt.Rectangle((0, y_pos - 0.4), 1, 0.8, 
                                   fill=True, edgecolor='black', facecolor='lightgray'))
        ax.text(-0.05, y_pos, chrom, ha='right', va='center')
        y_pos += 1

    # Plot markers
    for chrom in sorted(chrom_dict.keys()):
        y_base = chrom_positions[chrom]
        pos_data = chrom_dict[chrom]['positions']
        val_data = chrom_dict[chrom]['values']
        
        # Normalize positions
        if len(pos_data) > 1:
            pos_min, pos_max = min(pos_data), max(pos_data)
            normalized_pos = [(p - pos_min) / (pos_max - pos_min) if pos_max > pos_min else 0.5 
                             for p in pos_data]
        else:
            normalized_pos = [0.5]
        
        # Plot markers
        for pos, val in zip(normalized_pos, val_data):
            color = 'red' if threshold and val > threshold else 'blue'
            ax.scatter([pos], [y_base], c=color, s=50, **kwargs)

    ax.set_xlim(-0.2, 1.0)
    ax.set_ylim(-0.5, y_pos - 0.5)
    ax.set_xlabel("Normalized Position")
    ax.set_ylabel("Chromosome")
    ax.set_title("Chromosome Ideogram")
    ax.axis('off')

    return ax


def coverage_plot(
    positions: Sequence[int],
    coverage: Sequence[float],
    *,
    ax: plt.Axes | None = None,
    title: str = "Coverage Plot",
    **kwargs
) -> plt.Axes:
    """Create a coverage plot for sequencing data.

    Args:
        positions: Genomic positions
        coverage: Coverage values
        ax: Matplotlib axes (creates new if None)
        title: Plot title
        **kwargs: Additional arguments

    Returns:
        Matplotlib axes object

    Example:
        >>> from metainformant.visualization import coverage_plot
        >>> import numpy as np
        >>> positions = np.arange(1000, 2000, 10)
        >>> coverage = np.random.poisson(50, len(positions))
        >>> ax = coverage_plot(positions, coverage)
    """
    if ax is None:
        _, ax = plt.subplots(figsize=(12, 4))

    ax.fill_between(positions, coverage, alpha=0.5, **kwargs)
    ax.plot(positions, coverage, linewidth=1, **kwargs)

    ax.set_xlabel("Genomic Position")
    ax.set_ylabel("Coverage")
    ax.set_title(title)

    return ax


def variant_plot(
    positions: Sequence[int],
    ref_alleles: Sequence[str],
    alt_alleles: Sequence[str],
    frequencies: Sequence[float] | None = None,
    *,
    ax: plt.Axes | None = None,
    title: str = "Variant Plot",
    **kwargs
) -> plt.Axes:
    """Create a variant visualization plot.

    Args:
        positions: Genomic positions
        ref_alleles: Reference alleles
        alt_alleles: Alternative alleles
        frequencies: Optional allele frequencies
        ax: Matplotlib axes (creates new if None)
        title: Plot title
        **kwargs: Additional arguments

    Returns:
        Matplotlib axes object

    Example:
        >>> from metainformant.visualization import variant_plot
        >>> ax = variant_plot(
        ...     [1000, 2000, 3000],
        ...     ['A', 'T', 'G'],
        ...     ['G', 'C', 'A'],
        ...     [0.1, 0.2, 0.05]
        ... )
    """
    if ax is None:
        _, ax = plt.subplots(figsize=(12, 6))

    # Create variant labels
    variants = [f"{ref}→{alt}" for ref, alt in zip(ref_alleles, alt_alleles)]

    # Plot variants
    y_pos = 0
    for i, (pos, variant) in enumerate(zip(positions, variants)):
        if frequencies:
            size = frequencies[i] * 1000 if frequencies[i] else 100
            color = plt.cm.viridis(frequencies[i])
        else:
            size = 100
            color = 'blue'
        
        ax.scatter([pos], [y_pos], s=size, c=[color], **kwargs)
        ax.text(pos, y_pos + 0.1, variant, ha='center', va='bottom', fontsize=8)
        y_pos += 1

    ax.set_xlabel("Genomic Position")
    ax.set_ylabel("Variant")
    ax.set_title(title)
    ax.set_yticks([])

    return ax

