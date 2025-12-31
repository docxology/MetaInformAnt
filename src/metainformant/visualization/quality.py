"""Quality control data visualization functions.

This module provides specialized plotting functions for sequencing quality control
including quality metrics, adapter content, GC distribution, and read length analysis.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, List

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.axes import Axes

from metainformant.core import logging, paths, validation

logger = logging.get_logger(__name__)


def plot_quality_metrics(
    qc_data: Dict[str, Any],
    *,
    ax: Axes | None = None,
    output_path: str | Path | None = None,
    **kwargs
) -> Axes:
    """Create a comprehensive quality metrics visualization.

    Args:
        qc_data: Dictionary containing various QC metrics
        ax: Optional matplotlib axes to plot on. Creates new if None.
        output_path: Optional path to save the figure.
        **kwargs: Additional arguments passed to plotting.

    Returns:
        matplotlib Axes object

    Raises:
        ValueError: If qc_data is empty or malformed
    """
    validation.validate_type(qc_data, dict, "qc_data")

    if not qc_data:
        raise ValueError("QC data dictionary cannot be empty")

    if ax is None:
        fig, axes = plt.subplots(2, 2, figsize=kwargs.pop('figsize', (12, 10)))
        axes = axes.flatten()
    else:
        # If single ax provided, this is complex - just use it for a simple plot
        fig, axes = plt.subplots(2, 2, figsize=(12, 10))
        axes = axes.flatten()

    # Plot 1: Quality scores distribution
    if 'per_base_quality' in qc_data:
        qual_data = qc_data['per_base_quality']
        if 'positions' in qual_data and 'mean_qualities' in qual_data:
            axes[0].plot(qual_data['positions'], qual_data['mean_qualities'])
            axes[0].set_xlabel('Position in Read')
            axes[0].set_ylabel('Mean Quality Score')
            axes[0].set_title('Per-Base Quality')
            axes[0].grid(True, alpha=0.3)

    # Plot 2: GC content distribution
    if 'gc_content_distribution' in qc_data:
        gc_data = qc_data['gc_content_distribution']
        if 'bins' in gc_data and 'counts' in gc_data:
            # Convert bins to centers for plotting
            bin_centers = [(gc_data['bins'][i] + gc_data['bins'][i+1]) / 2
                          for i in range(len(gc_data['bins']) - 1)]
            axes[1].bar(bin_centers, gc_data['counts'], width=2, alpha=0.7)
            axes[1].set_xlabel('GC Content (%)')
            axes[1].set_ylabel('Frequency')
            axes[1].set_title('GC Content Distribution')
            axes[1].grid(True, alpha=0.3)

    # Plot 3: Read length distribution
    if 'sequence_length_distribution' in qc_data:
        length_data = qc_data['sequence_length_distribution']
        if 'lengths' in length_data and 'counts' in length_data:
            axes[2].bar(length_data['lengths'], length_data['counts'], alpha=0.7)
            axes[2].set_xlabel('Read Length')
            axes[2].set_ylabel('Count')
            axes[2].set_title('Read Length Distribution')
            axes[2].grid(True, alpha=0.3)

    # Plot 4: Basic statistics summary
    if 'basic_statistics' in qc_data:
        stats = qc_data['basic_statistics']
        stat_names = ['num_reads', 'total_bases', 'min_length', 'max_length', 'mean_length']
        stat_values = [stats.get(name, 0) for name in stat_names]

        bars = axes[3].bar(stat_names, stat_values, alpha=0.7)
        axes[3].set_ylabel('Value')
        axes[3].set_title('Basic Statistics')
        axes[3].tick_params(axis='x', rotation=45)

        # Add value labels on bars
        for bar, value in zip(bars, stat_values):
            axes[3].text(bar.get_x() + bar.get_width()/2, bar.get_height(),
                        f'{value:.0f}', ha='center', va='bottom', fontsize=8)

    plt.tight_layout()

    if output_path:
        output_path = paths.ensure_directory(Path(output_path).parent)
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        logger.info(f"Quality metrics plot saved to {output_path}")

    return axes[0]  # Return first axis for consistency


def plot_adapter_content(
    adapter_data: Dict[str, List[float]],
    *,
    ax: Axes | None = None,
    output_path: str | Path | None = None,
    **kwargs
) -> Axes:
    """Create an adapter content visualization.

    Args:
        adapter_data: Dictionary with adapter names as keys and content percentages as values
        ax: Optional matplotlib axes to plot on. Creates new if None.
        output_path: Optional path to save the figure.
        **kwargs: Additional arguments passed to plotting.

    Returns:
        matplotlib Axes object

    Raises:
        ValueError: If adapter_data is empty
    """
    validation.validate_type(adapter_data, dict, "adapter_data")

    if not adapter_data:
        raise ValueError("Adapter data dictionary cannot be empty")

    if ax is None:
        fig, ax = plt.subplots(figsize=kwargs.pop('figsize', (10, 6)))

    # Prepare data for plotting
    adapters = list(adapter_data.keys())
    # Assume each adapter has a list of percentages (could be by position)
    # For simplicity, take the maximum percentage for each adapter
    percentages = [max(values) if values else 0 for values in adapter_data.values()]

    # Sort by percentage
    sorted_indices = np.argsort(percentages)[::-1]
    adapters_sorted = [adapters[i] for i in sorted_indices]
    percentages_sorted = [percentages[i] for i in sorted_indices]

    bars = ax.bar(range(len(adapters_sorted)), percentages_sorted, alpha=0.7, **kwargs)
    ax.set_xlabel('Adapter Type')
    ax.set_ylabel('Maximum Content (%)')
    ax.set_title('Adapter Content Analysis')
    ax.set_xticks(range(len(adapters_sorted)))
    ax.set_xticklabels(adapters_sorted, rotation=45, ha='right')
    ax.grid(True, alpha=0.3)

    # Add percentage labels
    for bar, pct in zip(bars, percentages_sorted):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height(),
               f'{pct:.1f}%', ha='center', va='bottom', fontsize=8)

    if output_path:
        output_path = paths.ensure_directory(Path(output_path).parent)
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        logger.info(f"Adapter content plot saved to {output_path}")

    return ax


def plot_gc_distribution(
    gc_data: List[float],
    *,
    ax: Axes | None = None,
    output_path: str | Path | None = None,
    **kwargs
) -> Axes:
    """Create a GC content distribution plot.

    Args:
        gc_data: List of GC content percentages
        ax: Optional matplotlib axes to plot on. Creates new if None.
        output_path: Optional path to save the figure.
        **kwargs: Additional arguments passed to matplotlib hist().

    Returns:
        matplotlib Axes object

    Raises:
        ValueError: If gc_data is empty
    """
    validation.validate_type(gc_data, list, "gc_data")

    if not gc_data:
        raise ValueError("GC data list cannot be empty")

    if ax is None:
        fig, ax = plt.subplots(figsize=kwargs.pop('figsize', (8, 6)))

    # Plot histogram
    n, bins, patches = ax.hist(gc_data, bins=kwargs.get('bins', 20),
                              alpha=0.7, edgecolor='black', **kwargs)

    ax.set_xlabel('GC Content (%)')
    ax.set_ylabel('Frequency')
    ax.set_title('GC Content Distribution')
    ax.grid(True, alpha=0.3)

    # Add statistics as text
    mean_gc = np.mean(gc_data)
    median_gc = np.median(gc_data)
    std_gc = np.std(gc_data)

    stats_text = f'Mean: {mean_gc:.1f}%\nMedian: {median_gc:.1f}%\nStd: {std_gc:.1f}%'
    ax.text(0.02, 0.98, stats_text, transform=ax.transAxes,
           verticalalignment='top', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))

    if output_path:
        output_path = paths.ensure_directory(Path(output_path).parent)
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        logger.info(f"GC distribution plot saved to {output_path}")

    return ax


def plot_length_distribution(
    length_data: List[int],
    *,
    ax: Axes | None = None,
    output_path: str | Path | None = None,
    **kwargs
) -> Axes:
    """Create a read length distribution plot.

    Args:
        length_data: List of read lengths
        ax: Optional matplotlib axes to plot on. Creates new if None.
        output_path: Optional path to save the figure.
        **kwargs: Additional arguments passed to matplotlib hist().

    Returns:
        matplotlib Axes object

    Raises:
        ValueError: If length_data is empty
    """
    validation.validate_type(length_data, list, "length_data")

    if not length_data:
        raise ValueError("Length data list cannot be empty")

    if ax is None:
        fig, ax = plt.subplots(figsize=kwargs.pop('figsize', (8, 6)))

    # Plot histogram
    n, bins, patches = ax.hist(length_data, bins=kwargs.get('bins', 'auto'),
                              alpha=0.7, edgecolor='black', **kwargs)

    ax.set_xlabel('Read Length (bp)')
    ax.set_ylabel('Frequency')
    ax.set_title('Read Length Distribution')
    ax.grid(True, alpha=0.3)

    # Add statistics as text
    mean_len = np.mean(length_data)
    median_len = np.median(length_data)
    min_len = np.min(length_data)
    max_len = np.max(length_data)

    stats_text = f'Mean: {mean_len:.0f} bp\nMedian: {median_len:.0f} bp\nRange: {min_len}-{max_len} bp'
    ax.text(0.02, 0.98, stats_text, transform=ax.transAxes,
           verticalalignment='top', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))

    if output_path:
        output_path = paths.ensure_directory(Path(output_path).parent)
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        logger.info(f"Length distribution plot saved to {output_path}")

    return ax
