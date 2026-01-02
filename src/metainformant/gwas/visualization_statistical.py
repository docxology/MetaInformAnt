"""Statistical visualization functions for GWAS.

This module provides plots for statistical analysis and quality control.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, List, Optional

import numpy as np

from metainformant.core import logging

logger = logging.get_logger(__name__)


def lambda_gc_plot(lambda_values: List[float], output_file: Optional[str | Path] = None,
                  title: str = "Genomic Control Lambda Distribution") -> Optional[Any]:
    """Create a plot showing the distribution of genomic control lambda values.

    Args:
        lambda_values: List of lambda GC values from different analyses
        output_file: Optional output file path
        title: Plot title

    Returns:
        Plot object if matplotlib available, None otherwise

    Example:
        >>> lambdas = [1.05, 1.12, 1.08, 1.15]
        >>> plot = lambda_gc_plot(lambdas)
    """
    try:
        import matplotlib.pyplot as plt
        import seaborn as sns
    except ImportError:
        logger.warning("matplotlib or seaborn not available for lambda GC plot")
        return None

    if not lambda_values:
        logger.error("No lambda values provided")
        return None

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

    # Plot 1: Histogram of lambda values
    ax1.hist(lambda_values, bins=20, alpha=0.7, color='skyblue', edgecolor='black')
    ax1.axvline(x=np.mean(lambda_values), color='red', linestyle='--',
               linewidth=2, label=f'Mean: {np.mean(lambda_values):.3f}')
    ax1.axvline(x=1.0, color='green', linestyle='-', alpha=0.7,
               label='Expected (1.0)')
    ax1.set_xlabel('Genomic Control Lambda', fontsize=12)
    ax1.set_ylabel('Frequency', fontsize=12)
    ax1.set_title('Lambda GC Distribution', fontsize=14)
    ax1.legend()
    ax1.grid(True, alpha=0.3)

    # Plot 2: Q-Q plot against expected distribution
    expected_lambdas = np.random.normal(1.0, 0.1, len(lambda_values))
    expected_lambdas = np.sort(expected_lambdas)
    observed_lambdas = np.sort(lambda_values)

    ax2.scatter(expected_lambdas, observed_lambdas, alpha=0.7, color='orange')
    ax2.plot([min(expected_lambdas), max(expected_lambdas)],
             [min(expected_lambdas), max(expected_lambdas)],
             'k--', alpha=0.7, label='Expected')

    ax2.set_xlabel('Expected Lambda GC', fontsize=12)
    ax2.set_ylabel('Observed Lambda GC', fontsize=12)
    ax2.set_title('Q-Q Plot', fontsize=14)
    ax2.legend()
    ax2.grid(True, alpha=0.3)

    # Overall title
    fig.suptitle(title, fontsize=16, y=0.98)

    # Add summary statistics as text
    stats_text = f"""Summary Statistics:
Mean: {np.mean(lambda_values):.3f}
Median: {np.median(lambda_values):.3f}
SD: {np.std(lambda_values):.3f}
Range: {np.min(lambda_values):.3f} - {np.max(lambda_values):.3f}
Inflated: {np.mean(lambda_values) > 1.1}"""

    fig.text(0.02, 0.02, stats_text, fontsize=10,
             verticalalignment='bottom', fontfamily='monospace',
             bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))

    plt.tight_layout()

    if output_file:
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        logger.info(f"Saved lambda GC plot to {output_file}")

    return plt.gcf()


def power_plot(sample_sizes: List[int], effect_sizes: List[float],
              alpha: float = 0.05, output_file: Optional[str | Path] = None,
              title: str = "Statistical Power Analysis") -> Optional[Any]:
    """Create a power analysis plot showing statistical power vs sample size and effect size.

    Args:
        sample_sizes: List of sample sizes to evaluate
        effect_sizes: List of effect sizes to evaluate
        alpha: Significance level
        output_file: Optional output file path
        title: Plot title

    Returns:
        Plot object if matplotlib and scipy available, None otherwise
    """
    try:
        import matplotlib.pyplot as plt
        from scipy import stats
        import numpy as np
    except ImportError:
        logger.warning("matplotlib or scipy not available for power plot")
        return None

    if not sample_sizes or not effect_sizes:
        logger.error("No sample sizes or effect sizes provided")
        return None

    # Create meshgrid for power calculation
    n_grid, effect_grid = np.meshgrid(sample_sizes, effect_sizes)
    power_grid = np.zeros_like(n_grid, dtype=float)

    # Calculate power for each combination
    for i in range(len(effect_sizes)):
        for j in range(len(sample_sizes)):
            n = sample_sizes[j]
            effect = effect_sizes[i]

            # For simplicity, assume two-sample t-test power
            # Power = 1 - β, where β is probability of Type II error
            try:
                # Calculate non-centrality parameter
                ncp = effect * np.sqrt(n / 2)  # Approximation for equal sample sizes
                # Critical value for alpha
                t_crit = stats.t.ppf(1 - alpha/2, df=2*n-2)
                # Power calculation
                power = 1 - stats.nct.cdf(t_crit, df=2*n-2, nc=ncp)
                power_grid[i, j] = power
            except (ValueError, ZeroDivisionError):
                power_grid[i, j] = np.nan

    # Create plot
    fig, ax = plt.subplots(figsize=(10, 8))

    # Plot power surface
    cs = ax.contourf(n_grid, effect_grid, power_grid, levels=np.linspace(0, 1, 11),
                    cmap='RdYlBu_r', alpha=0.8)
    plt.colorbar(cs, ax=ax, label='Statistical Power')

    # Add contour lines
    contours = ax.contour(n_grid, effect_grid, power_grid, levels=[0.8, 0.9],
                         colors='black', linewidths=2)
    ax.clabel(contours, inline=True, fontsize=10)

    ax.set_xlabel('Sample Size', fontsize=12)
    ax.set_ylabel('Effect Size', fontsize=12)
    ax.set_title(title, fontsize=14, pad=20)
    ax.grid(True, alpha=0.3)

    # Add reference lines
    ax.axhline(y=0.2, color='red', linestyle='--', alpha=0.7, label='Small effect (d=0.2)')
    ax.axhline(y=0.5, color='orange', linestyle='--', alpha=0.7, label='Medium effect (d=0.5)')
    ax.axhline(y=0.8, color='green', linestyle='--', alpha=0.7, label='Large effect (d=0.8)')

    ax.legend()

    # Add text annotations for power thresholds
    ax.text(0.02, 0.98, 'Power Contours:\n0.8 (solid)\n0.9 (dashed)',
           transform=ax.transAxes, fontsize=10, verticalalignment='top',
           bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

    plt.tight_layout()

    if output_file:
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        logger.info(f"Saved power analysis plot to {output_file}")

    return fig
