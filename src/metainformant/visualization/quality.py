"""Quality control visualization functions.

This module provides visualization functions for quality control metrics including
QC metrics plots, quality score plots, per-base quality plots, adapter content
plots, and sequence length distributions.
"""

from __future__ import annotations

from typing import Sequence

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# Use non-interactive backend by default for tests/headless
matplotlib.use("Agg", force=True)


def qc_metrics_plot(
    metrics: dict[str, Sequence[float]],
    *,
    ncols: int = 3,
    figsize: tuple[int, int] | None = None,
    **kwargs
) -> plt.Figure:
    """Plot multiple quality control metrics.

    Args:
        metrics: Dictionary mapping metric names to values
        ncols: Number of columns in subplot grid
        figsize: Figure size tuple
        **kwargs: Additional arguments

    Returns:
        Matplotlib Figure object

    Example:
        >>> from metainformant.visualization import qc_metrics_plot
        >>> import numpy as np
        >>> metrics = {
        ...     'total_counts': np.random.poisson(1000, 100),
        ...     'n_genes': np.random.poisson(2000, 100),
        ...     'pct_mt': np.random.uniform(0, 10, 100)
        ... }
        >>> fig = qc_metrics_plot(metrics)
    """
    n_metrics = len(metrics)
    nrows = (n_metrics + ncols - 1) // ncols

    if figsize is None:
        figsize = (4 * ncols, 3 * nrows)

    fig, axes = plt.subplots(nrows, ncols, figsize=figsize)
    if nrows == 1 and ncols == 1:
        axes = [axes]
    elif nrows == 1 or ncols == 1:
        axes = axes.flatten() if hasattr(axes, 'flatten') else [axes]
    else:
        axes = axes.flatten()

    for i, (metric_name, values) in enumerate(metrics.items()):
        ax = axes[i]
        ax.hist(values, bins=50, alpha=0.7, edgecolor='black', linewidth=0.5, **kwargs)
        
        # Add statistics
        mean_val = np.mean(values)
        median_val = np.median(values)
        ax.axvline(mean_val, color='red', linestyle='--', alpha=0.7, label=f'Mean: {mean_val:.1f}')
        ax.axvline(median_val, color='blue', linestyle='--', alpha=0.7, label=f'Median: {median_val:.1f}')
        
        ax.set_xlabel(metric_name.replace('_', ' ').title())
        ax.set_ylabel("Frequency")
        ax.set_title(f'Distribution of {metric_name.replace("_", " ").title()}')
        ax.legend(fontsize=8)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

    # Hide unused subplots
    for i in range(n_metrics, len(axes)):
        axes[i].set_visible(False)

    plt.tight_layout()
    return fig


def quality_score_plot(
    quality_scores: Sequence[float],
    *,
    ax: plt.Axes | None = None,
    title: str = "Quality Score Distribution",
    **kwargs
) -> plt.Axes:
    """Plot distribution of quality scores.

    Args:
        quality_scores: Quality scores
        ax: Matplotlib axes (creates new if None)
        title: Plot title
        **kwargs: Additional arguments

    Returns:
        Matplotlib axes object

    Example:
        >>> from metainformant.visualization import quality_score_plot
        >>> import numpy as np
        >>> scores = np.random.uniform(0, 40, 1000)
        >>> ax = quality_score_plot(scores)
    """
    if ax is None:
        _, ax = plt.subplots(figsize=(10, 6))

    ax.hist(quality_scores, bins=50, alpha=0.7, edgecolor='black', **kwargs)
    
    # Add statistics
    mean_score = np.mean(quality_scores)
    median_score = np.median(quality_scores)
    ax.axvline(mean_score, color='red', linestyle='--', alpha=0.7, label=f'Mean: {mean_score:.1f}')
    ax.axvline(median_score, color='blue', linestyle='--', alpha=0.7, label=f'Median: {median_score:.1f}')
    
    ax.set_xlabel("Quality Score")
    ax.set_ylabel("Frequency")
    ax.set_title(title)
    ax.legend()
    ax.grid(True, alpha=0.3)

    return ax


def per_base_quality_plot(
    positions: Sequence[int],
    quality_scores: Sequence[float],
    *,
    ax: plt.Axes | None = None,
    title: str = "Per-Base Quality Scores",
    **kwargs
) -> plt.Axes:
    """Plot quality scores across read positions.

    Args:
        positions: Read positions
        quality_scores: Quality scores at each position
        ax: Matplotlib axes (creates new if None)
        title: Plot title
        **kwargs: Additional arguments

    Returns:
        Matplotlib axes object

    Example:
        >>> from metainformant.visualization import per_base_quality_plot
        >>> import numpy as np
        >>> positions = np.arange(1, 101)
        >>> quality = 30 - 0.1 * positions + np.random.normal(0, 2, 100)
        >>> ax = per_base_quality_plot(positions, quality)
    """
    if ax is None:
        _, ax = plt.subplots(figsize=(12, 6))

    positions = np.array(positions)
    quality_scores = np.array(quality_scores)

    # Plot mean quality
    ax.plot(positions, quality_scores, linewidth=2, **kwargs)
    
    # Add quality thresholds
    ax.axhline(y=28, color='green', linestyle='--', alpha=0.7, label='High quality (≥28)')
    ax.axhline(y=20, color='orange', linestyle='--', alpha=0.7, label='Medium quality (≥20)')
    ax.axhline(y=10, color='red', linestyle='--', alpha=0.7, label='Low quality (≥10)')
    
    ax.set_xlabel("Position in Read")
    ax.set_ylabel("Quality Score")
    ax.set_title(title)
    ax.legend()
    ax.grid(True, alpha=0.3)

    return ax


def adapter_content_plot(
    positions: Sequence[int],
    adapter_content: Sequence[float],
    *,
    threshold: float = 0.1,
    ax: plt.Axes | None = None,
    title: str = "Adapter Content",
    **kwargs
) -> plt.Axes:
    """Plot adapter content across read positions.

    Args:
        positions: Read positions
        adapter_content: Adapter content percentage at each position
        threshold: Threshold for adapter content warning
        ax: Matplotlib axes (creates new if None)
        title: Plot title
        **kwargs: Additional arguments

    Returns:
        Matplotlib axes object

    Example:
        >>> from metainformant.visualization import adapter_content_plot
        >>> import numpy as np
        >>> positions = np.arange(1, 101)
        >>> content = np.random.uniform(0, 0.2, 100)
        >>> ax = adapter_content_plot(positions, content)
    """
    if ax is None:
        _, ax = plt.subplots(figsize=(12, 6))

    positions = np.array(positions)
    adapter_content = np.array(adapter_content)

    ax.plot(positions, adapter_content * 100, linewidth=2, **kwargs)
    ax.axhline(y=threshold * 100, color='red', linestyle='--', alpha=0.7,
              label=f'Threshold ({threshold*100:.1f}%)')
    ax.fill_between(positions, 0, adapter_content * 100, alpha=0.3)

    ax.set_xlabel("Position in Read")
    ax.set_ylabel("Adapter Content (%)")
    ax.set_title(title)
    ax.legend()
    ax.grid(True, alpha=0.3)
    ax.set_ylim(0, max(adapter_content.max() * 100 * 1.1, threshold * 100 * 1.5))

    return ax


def sequence_length_distribution(
    lengths: Sequence[int],
    *,
    ax: plt.Axes | None = None,
    title: str = "Sequence Length Distribution",
    **kwargs
) -> plt.Axes:
    """Plot distribution of sequence lengths.

    Args:
        lengths: Sequence lengths
        ax: Matplotlib axes (creates new if None)
        title: Plot title
        **kwargs: Additional arguments

    Returns:
        Matplotlib axes object

    Example:
        >>> from metainformant.visualization import sequence_length_distribution
        >>> import numpy as np
        >>> lengths = np.random.normal(100, 10, 1000).astype(int)
        >>> ax = sequence_length_distribution(lengths)
    """
    if ax is None:
        _, ax = plt.subplots(figsize=(10, 6))

    lengths = np.array(lengths)

    ax.hist(lengths, bins=50, alpha=0.7, edgecolor='black', **kwargs)
    
    # Add statistics
    mean_length = np.mean(lengths)
    median_length = np.median(lengths)
    ax.axvline(mean_length, color='red', linestyle='--', alpha=0.7, label=f'Mean: {mean_length:.1f}')
    ax.axvline(median_length, color='blue', linestyle='--', alpha=0.7, label=f'Median: {median_length:.1f}')
    
    ax.set_xlabel("Sequence Length")
    ax.set_ylabel("Frequency")
    ax.set_title(title)
    ax.legend()
    ax.grid(True, alpha=0.3)

    return ax

