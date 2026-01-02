"""Population structure visualization for GWAS.

This module provides plots for visualizing population stratification and structure.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, List, Optional

import numpy as np

from metainformant.core import logging

logger = logging.get_logger(__name__)


def pca_plot(pca_data: Dict[str, Any], output_file: Optional[str | Path] = None,
             title: str = "PCA Plot", pc1: int = 0, pc2: int = 1) -> Optional[Any]:
    """Create a PCA plot for population structure visualization.

    Args:
        pca_data: Dictionary containing PCA data with keys:
                 'pcs': 2D array of principal components
                 'explained_variance': array of explained variance ratios
                 'sample_names': list of sample names (optional)
                 'populations': list of population labels (optional)
        output_file: Optional output file path
        title: Plot title
        pc1: Index of first principal component to plot (0-based)
        pc2: Index of second principal component to plot (0-based)

    Returns:
        Plot object if matplotlib available, None otherwise

    Example:
        >>> data = {'pcs': pcs_array, 'explained_variance': var_ratios}
        >>> plot = pca_plot(data)
    """
    try:
        import matplotlib.pyplot as plt
        import seaborn as sns
    except ImportError:
        logger.warning("matplotlib or seaborn not available for PCA plot")
        return None

    # Validate input data
    required_keys = ['pcs', 'explained_variance']
    missing_keys = [key for key in required_keys if key not in pca_data]
    if missing_keys:
        logger.error(f"Missing required keys in pca_data: {missing_keys}")
        return None

    pcs = pca_data['pcs']
    explained_var = pca_data['explained_variance']

    if len(pcs.shape) != 2:
        logger.error("PCA data 'pcs' must be a 2D array")
        return None

    if pc1 >= pcs.shape[1] or pc2 >= pcs.shape[1]:
        logger.error(f"PC indices {pc1}, {pc2} out of range for {pcs.shape[1]} components")
        return None

    fig, ax = plt.subplots(1, 1, figsize=(10, 8))

    # Extract PC data
    x_data = pcs[:, pc1]
    y_data = pcs[:, pc2]

    # Check for population labels
    populations = pca_data.get('populations')
    sample_names = pca_data.get('sample_names')

    if populations is not None and len(populations) == len(x_data):
        # Color by population
        unique_pops = list(set(populations))
        colors = plt.cm.tab10(np.linspace(0, 1, len(unique_pops)))

        for i, pop in enumerate(unique_pops):
            mask = [p == pop for p in populations]
            ax.scatter(x_data[mask], y_data[mask],
                      c=[colors[i]], label=pop, alpha=0.7, s=50)
        ax.legend(title='Population', bbox_to_anchor=(1.05, 1), loc='upper left')
    else:
        # Simple scatter plot
        ax.scatter(x_data, y_data, alpha=0.7, s=50, c='blue')

    # Add sample labels if provided and not too many samples
    if sample_names is not None and len(sample_names) <= 50:
        for i, name in enumerate(sample_names):
            ax.annotate(name, (x_data[i], y_data[i]),
                       xytext=(5, 5), textcoords='offset points',
                       fontsize=8, alpha=0.8)

    # Labels and title
    var1 = explained_var[pc1] * 100
    var2 = explained_var[pc2] * 100

    ax.set_xlabel(f'PC{pc1+1} ({var1:.1f}% variance)', fontsize=12)
    ax.set_ylabel(f'PC{pc2+1} ({var2:.1f}% variance)', fontsize=12)
    ax.set_title(title, fontsize=14, pad=20)
    ax.grid(True, alpha=0.3)

    # Add center lines
    ax.axhline(y=0, color='gray', linestyle='--', alpha=0.5)
    ax.axvline(x=0, color='gray', linestyle='--', alpha=0.5)

    plt.tight_layout()

    if output_file:
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        logger.info(f"Saved PCA plot to {output_file}")

    return plt.gcf()


def pca_scree_plot(explained_variance: List[float], output_file: Optional[str | Path] = None,
                  title: str = "PCA Scree Plot") -> Optional[Any]:
    """Create a scree plot showing explained variance by principal component.

    Args:
        explained_variance: List of explained variance ratios (0-1)
        output_file: Optional output file path
        title: Plot title

    Returns:
        Plot object if matplotlib available, None otherwise
    """
    try:
        import matplotlib.pyplot as plt
        import numpy as np
    except ImportError:
        logger.warning("matplotlib not available for scree plot")
        return None

    if not explained_variance:
        logger.error("No explained variance data provided")
        return None

    # Convert to percentages if not already
    if max(explained_variance) <= 1.0:
        explained_variance_pct = [v * 100 for v in explained_variance]
    else:
        explained_variance_pct = explained_variance

    fig, ax = plt.subplots(figsize=(10, 6))

    # Create scree plot
    components = list(range(1, len(explained_variance) + 1))
    bars = ax.bar(components, explained_variance_pct, color='skyblue', edgecolor='navy', alpha=0.7)

    # Add cumulative variance line
    cumulative = np.cumsum(explained_variance_pct)
    ax.plot(components, cumulative, 'r-o', linewidth=2, markersize=4, label='Cumulative')

    # Add value labels on bars
    for bar, pct in zip(bars, explained_variance_pct):
        height = bar.get_height()
        ax.text(bar.get_x() + bar.get_width()/2., height + 0.5,
               f'{pct:.1f}%', ha='center', va='bottom', fontsize=8)

    ax.set_xlabel('Principal Component', fontsize=12)
    ax.set_ylabel('Explained Variance (%)', fontsize=12)
    ax.set_title(title, fontsize=14, pad=20)
    ax.set_xticks(components)
    ax.grid(True, alpha=0.3, axis='y')
    ax.legend()

    # Add elbow detection heuristic (find where additional variance drops below 5%)
    if len(explained_variance_pct) > 1:
        for i in range(1, len(explained_variance_pct)):
            if explained_variance_pct[i] < 5.0:
                ax.axvline(x=i + 0.5, color='red', linestyle='--', alpha=0.7,
                          label=f'Elbow at PC{i+1}')
                break

    plt.tight_layout()

    if output_file:
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        logger.info(f"Saved PCA scree plot to {output_file}")

    return fig
