"""Genome-wide visualization functions for GWAS.

This module provides specialized plots for genomic data visualization.
"""

from __future__ import annotations

import math
from pathlib import Path
from typing import Any, Dict, List, Optional

import numpy as np

from metainformant.core import logging

logger = logging.get_logger(__name__)


def circular_manhattan_plot(results_df: Any, output_file: Optional[str | Path] = None,
                           significance_threshold: float = 5e-8,
                           title: str = "Circular Manhattan Plot") -> Optional[Any]:
    """Create a circular Manhattan plot for GWAS results.

    Args:
        results_df: DataFrame with columns 'CHR', 'BP', 'P'
        output_file: Optional output file path
        significance_threshold: P-value threshold for significance line
        title: Plot title

    Returns:
        Plot object if matplotlib available, None otherwise

    Example:
        >>> plot = circular_manhattan_plot(gwas_results)
    """
    try:
        import matplotlib.pyplot as plt
        from matplotlib.patches import Wedge
    except ImportError:
        logger.warning("matplotlib not available for circular Manhattan plot")
        return None

    # Validate input data
    if not hasattr(results_df, 'columns'):
        logger.error("Input data must be a DataFrame")
        return None

    required_cols = ['CHR', 'BP', 'P']
    missing_cols = [col for col in required_cols if col not in results_df.columns]
    if missing_cols:
        logger.error(f"Missing required columns: {missing_cols}")
        return None

    fig, ax = plt.subplots(1, 1, figsize=(10, 10))

    # Prepare data
    chroms = sorted(results_df['CHR'].unique())
    chrom_angles = {}

    # Calculate chromosome angles (equal spacing)
    angle_per_chrom = 360 / len(chroms)

    # Colors for chromosomes
    colors = plt.cm.tab20(np.linspace(0, 1, len(chroms)))

    # Plot each chromosome as a sector
    for i, chrom in enumerate(chroms):
        chrom_data = results_df[results_df['CHR'] == chrom].copy()
        if chrom_data.empty:
            continue

        start_angle = i * angle_per_chrom
        end_angle = (i + 1) * angle_per_chrom
        chrom_angles[chrom] = (start_angle, end_angle)

        # Convert p-values to radial positions
        chrom_data['neg_log_p'] = -np.log10(chrom_data['P'].clip(lower=1e-50))

        # Normalize positions within chromosome
        chrom_data['norm_pos'] = (chrom_data['BP'] - chrom_data['BP'].min()) / \
                                (chrom_data['BP'].max() - chrom_data['BP'].min())

        # Convert to polar coordinates
        angles = start_angle + chrom_data['norm_pos'] * angle_per_chrom
        radii = 1 + chrom_data['neg_log_p'] * 0.1  # Scale for visibility

        # Plot points
        ax.scatter(np.deg2rad(angles), radii, c=[colors[i]], s=1, alpha=0.6)

        # Add chromosome label
        label_angle = np.deg2rad(start_angle + angle_per_chrom / 2)
        label_radius = 1.2
        ax.text(label_angle, label_radius, str(chrom),
                ha='center', va='center', fontsize=10, fontweight='bold')

    # Add significance threshold ring
    threshold_radius = 1 + (-np.log10(significance_threshold)) * 0.1
    theta = np.linspace(0, 2*np.pi, 100)
    ax.plot(theta, [threshold_radius] * 100, 'r--', linewidth=2, alpha=0.8)

    # Add significance label
    ax.text(0, threshold_radius + 0.05, f'p = {significance_threshold}',
            ha='center', va='bottom', fontsize=10, color='red')

    # Configure plot
    ax.set_title(title, fontsize=14, pad=20)
    ax.set_rlabel_position(0)
    ax.set_rticks([])
    ax.set_thetagrids([])
    ax.grid(False)
    ax.spines['polar'].set_visible(False)

    # Set equal aspect ratio
    ax.set_aspect('equal')

    plt.tight_layout()

    if output_file:
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        logger.info(f"Saved circular Manhattan plot to {output_file}")

    return plt.gcf()


def chromosome_ideogram(chromosome_lengths: Dict[str, int],
                       highlighted_regions: Optional[List[Dict[str, Any]]] = None,
                       output_file: Optional[str | Path] = None,
                       title: str = "Chromosome Ideogram") -> Optional[Any]:
    """Create a chromosome ideogram showing chromosome structure.

    Args:
        chromosome_lengths: Dictionary mapping chromosome names to lengths
        highlighted_regions: Optional list of regions to highlight
        output_file: Optional output file path
        title: Plot title

    Returns:
        Plot object if matplotlib available, None otherwise
    """
    try:
        import matplotlib.pyplot as plt
        from matplotlib.patches import Rectangle
        import numpy as np
    except ImportError:
        logger.warning("matplotlib not available for chromosome ideogram")
        return None

    if not chromosome_lengths:
        logger.error("No chromosome lengths provided")
        return None

    # Sort chromosomes
    chroms = sorted(chromosome_lengths.keys(),
                   key=lambda x: int(x.replace('chr', '')) if x.replace('chr', '').isdigit() else float('inf'))

    fig, ax = plt.subplots(figsize=(12, 8))

    y_start = 0.1
    y_height = 0.8 / len(chroms) if chroms else 0.8
    colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b']

    max_length = max(chromosome_lengths.values()) if chromosome_lengths else 1

    for i, chrom in enumerate(chroms):
        length = chromosome_lengths[chrom]
        color = colors[i % len(colors)]

        # Draw chromosome body
        width = length / max_length * 0.8  # Scale to 80% of plot width
        rect = Rectangle((0.1, y_start + i * y_height), width, y_height * 0.8,
                        facecolor=color, alpha=0.7, edgecolor='black', linewidth=1)
        ax.add_patch(rect)

        # Add chromosome label
        ax.text(0.05, y_start + i * y_height + y_height * 0.4, chrom,
               ha='right', va='center', fontsize=10, fontweight='bold')

        # Add length label
        length_mb = length / 1e6
        ax.text(width + 0.12, y_start + i * y_height + y_height * 0.4,
               f'{length_mb:.1f} Mb', ha='left', va='center', fontsize=8)

        # Highlight regions if provided
        if highlighted_regions:
            for region in highlighted_regions:
                if region.get('chrom') == chrom:
                    region_start = region.get('start', 0)
                    region_end = region.get('end', 0)
                    if region_end > region_start:
                        # Convert to plot coordinates
                        x_start = 0.1 + (region_start / max_length) * 0.8
                        x_width = ((region_end - region_start) / max_length) * 0.8

                        highlight_color = region.get('color', 'red')
                        alpha = region.get('alpha', 0.5)

                        highlight_rect = Rectangle((x_start, y_start + i * y_height),
                                                 x_width, y_height * 0.8,
                                                 facecolor=highlight_color, alpha=alpha,
                                                 edgecolor='darkred', linewidth=1)
                        ax.add_patch(highlight_rect)

    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_aspect('equal')
    ax.axis('off')
    ax.set_title(title, fontsize=14, pad=20)

    plt.tight_layout()

    if output_file:
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        logger.info(f"Saved chromosome ideogram to {output_file}")

    return fig


def genome_wide_ld_heatmap(
    ld_data: List[Dict[str, Any]],
    output_file: Optional[str | Path] = None,
    figsize: tuple[int, int] = (12, 10),
    title: str = "Genome-wide LD Heatmap",
    chromosomes: Optional[List[str]] = None,
    ld_threshold: float = 0.8
) -> Dict[str, Any]:
    """Create a genome-wide LD heatmap visualization.

    Args:
        ld_data: List of LD data dictionaries with CHROM, POS1, POS2, r2 keys
        output_file: Optional output file path
        figsize: Figure size (width, height)
        title: Plot title
        chromosomes: List of chromosomes to include (None for all)
        ld_threshold: LD threshold for highlighting high LD regions

    Returns:
        Dictionary with status and metadata
    """
    try:
        import matplotlib.pyplot as plt
        import matplotlib.patches as patches
        import numpy as np
    except ImportError:
        logger.warning("matplotlib not available for genome-wide LD heatmap")
        return {"status": "error", "message": "matplotlib not available"}

    if not ld_data:
        logger.warning("No LD data provided")
        return {"status": "error", "message": "No LD data provided"}

    # Group data by chromosome
    chrom_data = {}
    for entry in ld_data:
        chrom = entry['CHROM']
        if chromosomes and chrom not in chromosomes:
            continue
        if chrom not in chrom_data:
            chrom_data[chrom] = []
        chrom_data[chrom].append(entry)

    if not chrom_data:
        logger.warning("No data found for specified chromosomes")
        return {"status": "error", "message": "No data for specified chromosomes"}

    # Create figure with subplots for each chromosome
    n_chroms = len(chrom_data)
    if n_chroms <= 3:
        nrows, ncols = 1, n_chroms
    elif n_chroms <= 6:
        nrows, ncols = 2, 3
    else:
        nrows = int(np.ceil(n_chroms / 4))
        ncols = 4

    fig, axes = plt.subplots(nrows, ncols, figsize=figsize,
                           squeeze=False, sharex=True, sharey=True)
    axes = axes.flatten()

    # Plot each chromosome
    for i, (chrom, data) in enumerate(sorted(chrom_data.items())):
        if i >= len(axes):
            break

        ax = axes[i]

        # Extract positions and LD values
        pos1 = [entry['POS1'] for entry in data]
        pos2 = [entry['POS2'] for entry in data]
        r2_values = [entry['r2'] for entry in data]

        # Create scatter plot colored by LD
        scatter = ax.scatter(pos1, pos2, c=r2_values, cmap='RdYlBu_r',
                           s=2, alpha=0.6, vmin=0, vmax=1)

        # Highlight high LD regions
        high_ld = [entry for entry in data if entry['r2'] >= ld_threshold]
        if high_ld:
            ax.scatter([entry['POS1'] for entry in high_ld],
                      [entry['POS2'] for entry in high_ld],
                      c='red', s=5, alpha=0.8, marker='x',
                      label=f'r² ≥ {ld_threshold}')

        ax.set_xlabel('Position 1 (bp)')
        ax.set_ylabel('Position 2 (bp)')
        ax.set_title(f'{chrom}')
        ax.grid(True, alpha=0.3)
        ax.legend(loc='upper right', fontsize=8)

        # Format axis labels
        def format_bp(x, pos):
            if x >= 1e6:
                return f'{x/1e6:.1f}M'
            elif x >= 1e3:
                return f'{x/1e3:.0f}K'
            else:
                return f'{x:.0f}'

        ax.xaxis.set_major_formatter(plt.FuncFormatter(format_bp))
        ax.yaxis.set_major_formatter(plt.FuncFormatter(format_bp))

    # Hide unused subplots
    for i in range(len(chrom_data), len(axes)):
        axes[i].set_visible(False)

    # Add colorbar
    cbar_ax = fig.add_axes([0.92, 0.15, 0.02, 0.7])
    cbar = fig.colorbar(scatter, cax=cbar_ax)
    cbar.set_label('LD (r²)')

    fig.suptitle(title, fontsize=14)

    plt.tight_layout(rect=[0, 0, 0.9, 0.95])

    if output_file:
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        logger.info(f"Saved genome-wide LD heatmap to {output_file}")

    return {
        "status": "success",
        "n_chromosomes": len(chrom_data),
        "total_ld_pairs": len(ld_data),
        "chromosomes": list(chrom_data.keys()),
        "ld_threshold": ld_threshold
    }
