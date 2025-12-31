"""GWAS effect size visualization utilities.

This module provides tools for visualizing and analyzing effect sizes
from GWAS results, including effect size distributions, comparisons,
and regional effect plots.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, List, Optional

import numpy as np

from metainformant.core import logging

logger = logging.get_logger(__name__)


def effect_size_plot(results: Any, output_file: Optional[str | Path] = None) -> Optional[Any]:
    """Create effect size distribution and analysis plot.

    Args:
        results: GWAS results DataFrame or dictionary
        output_file: Optional output file path

    Returns:
        Plot object if matplotlib available, None otherwise

    Example:
        >>> plot = effect_size_plot(gwas_results)
    """
    try:
        import matplotlib.pyplot as plt
        import seaborn as sns
        import pandas as pd
    except ImportError:
        logger.warning("matplotlib, seaborn, or pandas not available for plotting")
        return None

    # Convert to DataFrame if needed
    if isinstance(results, dict):
        results = pd.DataFrame(results)

    if not hasattr(results, 'columns') or 'BETA' not in results.columns:
        logger.warning("Results must contain BETA column for effect size plotting")
        return None

    fig, axes = plt.subplots(2, 2, figsize=(15, 12))
    fig.suptitle('GWAS Effect Size Analysis', fontsize=16)

    # 1. Effect size distribution
    ax1 = axes[0, 0]
    beta_values = results['BETA'].dropna()

    # Plot histogram
    ax1.hist(beta_values, bins=50, alpha=0.7, color='skyblue', edgecolor='black')
    ax1.axvline(x=0, color='red', linestyle='--', alpha=0.7, label='No effect')
    ax1.set_title('Effect Size Distribution')
    ax1.set_xlabel('Effect Size (Beta)')
    ax1.set_ylabel('Frequency')
    ax1.legend()
    ax1.grid(True, alpha=0.3)

    # 2. Effect size by chromosome
    ax2 = axes[0, 1]
    if 'CHR' in results.columns:
        # Group by chromosome
        chrom_effects = []
        chrom_labels = []

        for chrom in sorted(results['CHR'].unique()):
            chrom_data = results[results['CHR'] == chrom]['BETA'].dropna()
            if len(chrom_data) > 0:
                chrom_effects.append(chrom_data.values)
                chrom_labels.append(str(chrom))

        if chrom_effects:
            ax2.boxplot(chrom_effects, labels=chrom_labels)
            ax2.set_title('Effect Size by Chromosome')
            ax2.set_xlabel('Chromosome')
            ax2.set_ylabel('Effect Size (Beta)')
            ax2.grid(True, alpha=0.3)

    # 3. Effect size vs p-value
    ax3 = axes[1, 0]
    if 'P' in results.columns:
        p_values = results['P'].dropna()
        beta_subset = results['BETA'].dropna()

        # Align the data
        min_len = min(len(p_values), len(beta_subset))
        p_subset = p_values[:min_len]
        beta_aligned = beta_subset[:min_len]

        scatter = ax3.scatter(-np.log10(p_subset), beta_aligned,
                            alpha=0.6, s=2, c=beta_aligned,
                            cmap='RdYlBu_r')
        ax3.set_title('Effect Size vs Significance')
        ax3.set_xlabel('-log10(P-value)')
        ax3.set_ylabel('Effect Size (Beta)')
        plt.colorbar(scatter, ax=ax3, label='Effect Size')

    # 4. Effect size summary statistics
    ax4 = axes[1, 1]
    ax4.axis('off')

    # Calculate statistics
    stats = {
        'Mean Effect': '.4f',
        'Median Effect': '.4f',
        'Std Effect': '.4f',
        'Min Effect': '.4f',
        'Max Effect': '.4f',
        'Large Effects (>0.5)': (abs(beta_values) > 0.5).sum(),
        'Total SNPs': len(beta_values)
    }

    # Create table
    table_data = [[key, value] for key, value in stats.items()]
    table = ax4.table(cellText=table_data, colLabels=['Statistic', 'Value'],
                     loc='center', cellLoc='center')
    table.auto_set_font_size(False)
    table.set_fontsize(10)
    table.scale(1, 1.2)
    ax4.set_title('Effect Size Summary', pad=20)

    plt.tight_layout()

    if output_file:
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        logger.info(f"Saved effect size plot to {output_file}")

    return plt.gcf()


def compare_effect_sizes(study_results: Dict[str, Any],
                        output_file: Optional[str | Path] = None) -> Optional[Any]:
    """Compare effect sizes across different studies or methods.

    Args:
        study_results: Dictionary mapping study names to results
        output_file: Optional output file path

    Returns:
        Plot object if matplotlib available, None otherwise

    Example:
        >>> studies = {"Study1": results1, "Study2": results2}
        >>> plot = compare_effect_sizes(studies)
    """
    try:
        import matplotlib.pyplot as plt
        import seaborn as sns
        import pandas as pd
    except ImportError:
        logger.warning("matplotlib, seaborn, or pandas not available for plotting")
        return None

    if not study_results:
        return None

    fig, axes = plt.subplots(2, 2, figsize=(15, 12))
    fig.suptitle('Effect Size Comparison', fontsize=16)

    # Prepare data
    study_data = {}
    for study_name, results in study_results.items():
        if isinstance(results, dict):
            results = pd.DataFrame(results)

        if hasattr(results, 'columns') and 'BETA' in results.columns:
            study_data[study_name] = results['BETA'].dropna().values

    if not study_data:
        return None

    study_names = list(study_data.keys())

    # 1. Effect size distributions
    ax1 = axes[0, 0]
    for study_name, beta_values in study_data.items():
        ax1.hist(beta_values, bins=30, alpha=0.5, label=study_name, density=True)

    ax1.set_title('Effect Size Distributions')
    ax1.set_xlabel('Effect Size (Beta)')
    ax1.set_ylabel('Density')
    ax1.legend()
    ax1.grid(True, alpha=0.3)

    # 2. Box plot comparison
    ax2 = axes[0, 1]
    beta_lists = list(study_data.values())
    ax2.boxplot(beta_lists, labels=study_names)
    ax2.set_title('Effect Size Comparison')
    ax2.set_ylabel('Effect Size (Beta)')
    ax2.grid(True, alpha=0.3)

    # 3. Effect size correlation (if multiple studies)
    ax3 = axes[1, 0]
    if len(study_names) >= 2:
        # Compare first two studies
        study1_betas = study_data[study_names[0]]
        study2_betas = study_data[study_names[1]]

        # Align by SNP if possible (simplified)
        min_len = min(len(study1_betas), len(study2_betas))
        study1_subset = study1_betas[:min_len]
        study2_subset = study2_betas[:min_len]

        scatter = ax3.scatter(study1_subset, study2_subset, alpha=0.6, s=2)
        ax3.set_title(f'Effect Size Correlation: {study_names[0]} vs {study_names[1]}')
        ax3.set_xlabel(f'{study_names[0]} Effect Size')
        ax3.set_ylabel(f'{study_names[1]} Effect Size')

        # Add correlation line
        if len(study1_subset) > 1:
            corr = np.corrcoef(study1_subset, study2_subset)[0, 1]
            ax3.text(0.05, 0.95, f'r = {corr:.3f}',
                    transform=ax3.transAxes, fontsize=12,
                    verticalalignment='top',
                    bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

        # Add diagonal line
        min_val = min(min(study1_subset), min(study2_subset))
        max_val = max(max(study1_subset), max(study2_subset))
        ax3.plot([min_val, max_val], [min_val, max_val], 'r--', alpha=0.7)

    # 4. Effect size statistics comparison
    ax4 = axes[1, 1]
    ax4.axis('off')

    # Calculate statistics for each study
    stats_data = []
    for study_name, beta_values in study_data.items():
        stats_data.append([
            study_name,
            '.4f',
            '.4f',
            (abs(beta_values) > 0.5).sum(),
            len(beta_values)
        ])

    if stats_data:
        table = ax4.table(cellText=stats_data,
                         colLabels=['Study', 'Mean |Effect|', 'Max |Effect|', 'Large Effects', 'Total SNPs'],
                         loc='center', cellLoc='center')
        table.auto_set_font_size(False)
        table.set_fontsize(9)
        table.scale(1, 1.2)

    ax4.set_title('Effect Size Statistics', pad=20)

    plt.tight_layout()

    if output_file:
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        logger.info(f"Saved effect size comparison plot to {output_file}")

    return plt.gcf()


def effect_size_manhattan(results: Any, significance_threshold: float = 5e-8,
                         output_file: Optional[str | Path] = None) -> Optional[Any]:
    """Create Manhattan plot colored by effect size.

    Args:
        results: GWAS results DataFrame or dictionary
        significance_threshold: P-value threshold for significance
        output_file: Optional output file path

    Returns:
        Plot object if matplotlib available, None otherwise

    Example:
        >>> plot = effect_size_manhattan(gwas_results)
    """
    try:
        import matplotlib.pyplot as plt
        import pandas as pd
    except ImportError:
        logger.warning("matplotlib or pandas not available for plotting")
        return None

    # Convert to DataFrame if needed
    if isinstance(results, dict):
        results = pd.DataFrame(results)

    required_cols = ['CHR', 'BP', 'P', 'BETA']
    if not all(col in results.columns for col in required_cols):
        logger.warning(f"Results must contain columns: {required_cols}")
        return None

    fig, ax = plt.subplots(1, 1, figsize=(15, 8))

    # Prepare data
    df = results.dropna(subset=required_cols).copy()
    df['-logP'] = -np.log10(df['P'].clip(lower=1e-300))

    # Create chromosome positions
    chromosomes = sorted(df['CHR'].unique())
    chrom_starts = {}
    current_pos = 0

    for chrom in chromosomes:
        chrom_starts[chrom] = current_pos
        chrom_data = df[df['CHR'] == chrom]
        current_pos += len(chrom_data) + 1000  # Gap between chromosomes

    df['pos'] = df.apply(lambda row: chrom_starts[row['CHR']] +
                        df[df['CHR'] == row['CHR']].index.get_loc(row.name), axis=1)

    # Color by effect size
    effect_sizes = df['BETA'].abs()
    max_effect = effect_sizes.max()

    # Plot points colored by effect size
    scatter = ax.scatter(df['pos'], df['-logP'],
                        c=effect_sizes, cmap='RdYlBu_r',
                        alpha=0.8, s=3, vmin=0, vmax=max_effect)

    # Add significance threshold line
    ax.axhline(y=-np.log10(significance_threshold), color='red',
              linestyle='--', alpha=0.7, label='.0e')

    # Add chromosome labels
    for chrom in chromosomes:
        chrom_center = chrom_starts[chrom] + len(df[df['CHR'] == chrom]) / 2
        ax.text(chrom_center, ax.get_ylim()[0] - 0.5, str(chrom),
               ha='center', va='top', fontsize=10)

    ax.set_title('Manhattan Plot (Colored by Effect Size)', fontsize=14)
    ax.set_xlabel('Chromosome', fontsize=12)
    ax.set_ylabel('-log10(P-value)', fontsize=12)
    ax.legend()
    ax.grid(True, alpha=0.3)

    # Add colorbar
    cbar = plt.colorbar(scatter, ax=ax)
    cbar.set_label('Absolute Effect Size', fontsize=12)

    plt.tight_layout()

    if output_file:
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        logger.info(f"Saved effect size Manhattan plot to {output_file}")

    return plt.gcf()


def effect_size_qq(results: Any, output_file: Optional[str | Path] = None) -> Optional[Any]:
    """Create Q-Q plot with effect size information.

    Args:
        results: GWAS results DataFrame or dictionary
        output_file: Optional output file path

    Returns:
        Plot object if matplotlib available, None otherwise

    Example:
        >>> plot = effect_size_qq(gwas_results)
    """
    try:
        import matplotlib.pyplot as plt
        import pandas as pd
    except ImportError:
        logger.warning("matplotlib or pandas not available for plotting")
        return None

    # Convert to DataFrame if needed
    if isinstance(results, dict):
        results = pd.DataFrame(results)

    if 'P' not in results.columns or 'BETA' not in results.columns:
        logger.warning("Results must contain P and BETA columns")
        return None

    fig, ax = plt.subplots(1, 1, figsize=(10, 8))

    # Prepare Q-Q data
    p_values = results['P'].dropna().values
    p_values = np.sort(p_values)

    # Expected p-values under null
    n = len(p_values)
    expected = np.arange(1, n + 1) / (n + 1)

    # Convert to -log10
    observed_log = -np.log10(p_values.clip(lower=1e-300))
    expected_log = -np.log10(expected)

    # Color by effect size
    beta_values = results['BETA'].dropna().values
    beta_sorted = beta_values[np.argsort(results['P'].dropna().values)]
    effect_sizes = np.abs(beta_sorted)

    # Plot
    scatter = ax.scatter(expected_log, observed_log,
                        c=effect_sizes, cmap='RdYlBu_r',
                        alpha=0.7, s=4)

    # Add diagonal line
    ax.plot([0, max(expected_log.max(), observed_log.max())],
            [0, max(expected_log.max(), observed_log.max())],
            'k--', alpha=0.7, label='Expected')

    ax.set_title('Q-Q Plot (Colored by Effect Size)', fontsize=14)
    ax.set_xlabel('Expected -log10(P)', fontsize=12)
    ax.set_ylabel('Observed -log10(P)', fontsize=12)
    ax.legend()
    ax.grid(True, alpha=0.3)

    # Add colorbar
    cbar = plt.colorbar(scatter, ax=ax)
    cbar.set_label('Absolute Effect Size', fontsize=12)

    # Add inflation factor annotation
    lambda_gc = np.median(p_values) / 0.456
    ax.text(0.05, 0.95, f'λ = {lambda_gc:.3f}',
           transform=ax.transAxes, fontsize=12,
           verticalalignment='top',
           bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

    plt.tight_layout()

    if output_file:
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        logger.info(f"Saved effect size Q-Q plot to {output_file}")

    return plt.gcf()


