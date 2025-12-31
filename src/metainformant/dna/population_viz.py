"""Population genetics visualization utilities.

This module provides plotting and visualization functions for population
genetics data including F_ST heatmaps, selection plots, and demographic visualizations.
"""

from __future__ import annotations

from typing import Dict, List, Optional, Tuple

import numpy as np

from metainformant.core import logging

logger = logging.get_logger(__name__)


def plot_fst_matrix(populations: Dict[str, List[str]], output_file: Optional[str] = None) -> Optional[any]:
    """Create F_ST matrix heatmap visualization.

    Args:
        populations: Dictionary mapping population names to sequence lists
        output_file: Optional output file path

    Returns:
        Plot object if matplotlib available, None otherwise

    Example:
        >>> pops = {"pop1": ["ATCG", "ATCG"], "pop2": ["GCTA", "GCTA"]}
        >>> plot = plot_fst_matrix(pops)
    """
    try:
        import matplotlib.pyplot as plt
        import seaborn as sns
    except ImportError:
        logger.warning("matplotlib and/or seaborn not available for plotting")
        return None

    pop_names = list(populations.keys())
    n_pops = len(pop_names)

    # Calculate F_ST matrix
    fst_matrix = np.zeros((n_pops, n_pops))

    for i in range(n_pops):
        for j in range(i + 1, n_pops):
            pop1_seqs = populations[pop_names[i]]
            pop2_seqs = populations[pop_names[j]]

            from .population_analysis import calculate_fst
            fst = calculate_fst(pop1_seqs, pop2_seqs)

            fst_matrix[i, j] = fst
            fst_matrix[j, i] = fst

    # Create heatmap
    plt.figure(figsize=(8, 6))
    sns.heatmap(fst_matrix, annot=True, cmap='YlOrRd', square=True,
                xticklabels=pop_names, yticklabels=pop_names)
    plt.title('F_ST Matrix')
    plt.tight_layout()

    if output_file:
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        logger.info(f"Saved F_ST matrix plot to {output_file}")

    return plt.gcf()


def plot_tajima_d_distribution(tajima_d_values: List[float], output_file: Optional[str] = None) -> Optional[any]:
    """Plot distribution of Tajima's D values.

    Args:
        tajima_d_values: List of Tajima's D statistics
        output_file: Optional output file path

    Returns:
        Plot object if matplotlib available, None otherwise

    Example:
        >>> tajima_values = [0.1, -0.2, 0.3, -0.1]
        >>> plot = plot_tajima_d_distribution(tajima_values)
    """
    try:
        import matplotlib.pyplot as plt
        import seaborn as sns
    except ImportError:
        logger.warning("matplotlib and/or seaborn not available for plotting")
        return None

    plt.figure(figsize=(10, 6))

    # Plot histogram
    plt.hist(tajima_d_values, bins=30, alpha=0.7, color='skyblue', edgecolor='black')
    plt.axvline(x=0, color='red', linestyle='--', alpha=0.7, label='Neutral expectation')

    plt.xlabel("Tajima's D")
    plt.ylabel("Frequency")
    plt.title("Distribution of Tajima's D Values")
    plt.legend()
    plt.grid(True, alpha=0.3)

    if output_file:
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        logger.info(f"Saved Tajima's D distribution plot to {output_file}")

    return plt.gcf()


def plot_selection_statistics(statistics: Dict[str, List[float]], output_file: Optional[str] = None) -> Optional[any]:
    """Plot multiple selection statistics.

    Args:
        statistics: Dictionary mapping statistic names to value lists
        output_file: Optional output file path

    Returns:
        Plot object if matplotlib available, None otherwise

    Example:
        >>> stats = {"tajima_d": [0.1, -0.2], "fu_li_d": [0.05, -0.1]}
        >>> plot = plot_selection_statistics(stats)
    """
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        logger.warning("matplotlib not available for plotting")
        return None

    plt.figure(figsize=(12, 8))

    stat_names = list(statistics.keys())
    n_stats = len(stat_names)

    for i, (stat_name, values) in enumerate(statistics.items()):
        plt.subplot(n_stats, 1, i + 1)
        plt.plot(values, 'o-', alpha=0.7)
        plt.axhline(y=0, color='red', linestyle='--', alpha=0.5)
        plt.ylabel(stat_name)
        plt.title(f"{stat_name} along sequence")
        plt.grid(True, alpha=0.3)

    plt.tight_layout()

    if output_file:
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        logger.info(f"Saved selection statistics plot to {output_file}")

    return plt.gcf()


def plot_population_diversity(diversity_data: Dict[str, float], output_file: Optional[str] = None) -> Optional[any]:
    """Plot population diversity comparison.

    Args:
        diversity_data: Dictionary mapping population names to diversity values
        output_file: Optional output file path

    Returns:
        Plot object if matplotlib available, None otherwise

    Example:
        >>> diversity = {"pop1": 0.01, "pop2": 0.015, "pop3": 0.008}
        >>> plot = plot_population_diversity(diversity)
    """
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        logger.warning("matplotlib not available for plotting")
        return None

    plt.figure(figsize=(10, 6))

    populations = list(diversity_data.keys())
    diversities = list(diversity_data.values())

    bars = plt.bar(populations, diversities, color='lightblue', alpha=0.7, edgecolor='black')
    plt.ylabel('Nucleotide Diversity (π)')
    plt.title('Population Diversity Comparison')
    plt.xticks(rotation=45)

    # Add value labels on bars
    for bar, value in zip(bars, diversities):
        plt.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.001,
                '.4f', ha='center', va='bottom')

    plt.grid(True, alpha=0.3, axis='y')

    if output_file:
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        logger.info(f"Saved population diversity plot to {output_file}")

    return plt.gcf()


def plot_ld_decay(ld_data: List[Tuple[int, float]], output_file: Optional[str] = None) -> Optional[any]:
    """Plot linkage disequilibrium decay with distance.

    Args:
        ld_data: List of tuples (distance, ld_value)
        output_file: Optional output file path

    Returns:
        Plot object if matplotlib available, None otherwise

    Example:
        >>> ld_data = [(1, 0.8), (2, 0.6), (3, 0.4)]
        >>> plot = plot_ld_decay(ld_data)
    """
    try:
        import matplotlib.pyplot as plt
        import numpy as np
    except ImportError:
        logger.warning("matplotlib not available for plotting")
        return None

    if not ld_data:
        return None

    distances, ld_values = zip(*ld_data)

    plt.figure(figsize=(10, 6))
    plt.plot(distances, ld_values, 'o-', color='darkgreen', alpha=0.7, markersize=4)

    # Add exponential decay fit if enough points
    if len(distances) > 3:
        try:
            from scipy.optimize import curve_fit

            def exp_decay(x, a, b):
                return a * np.exp(-b * x)

            popt, _ = curve_fit(exp_decay, distances, ld_values, p0=[1, 0.1])
            x_fit = np.linspace(min(distances), max(distances), 100)
            y_fit = exp_decay(x_fit, *popt)
            plt.plot(x_fit, y_fit, '--', color='red', alpha=0.7,
                    label='.3f')
        except ImportError:
            pass  # scipy not available

    plt.xlabel('Distance (bp)')
    plt.ylabel('Linkage Disequilibrium (r²)')
    plt.title('LD Decay with Distance')
    plt.grid(True, alpha=0.3)
    plt.legend()

    if output_file:
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        logger.info(f"Saved LD decay plot to {output_file}")

    return plt.gcf()


def plot_population_structure(pca_coords: np.ndarray, cluster_labels: Optional[List[int]] = None,
                            output_file: Optional[str] = None) -> Optional[any]:
    """Plot population structure using PCA coordinates.

    Args:
        pca_coords: PCA coordinates (n_samples, n_components)
        cluster_labels: Optional cluster assignments for coloring
        output_file: Optional output file path

    Returns:
        Plot object if matplotlib available, None otherwise

    Example:
        >>> import numpy as np
        >>> coords = np.random.rand(50, 2)
        >>> plot = plot_population_structure(coords)
    """
    try:
        import matplotlib.pyplot as plt
        from mpl_toolkits.mplot3d import Axes3D
    except ImportError:
        logger.warning("matplotlib not available for plotting")
        return None

    if pca_coords.shape[1] < 2:
        logger.warning("Need at least 2 PCA components for plotting")
        return None

    plt.figure(figsize=(10, 8))

    if pca_coords.shape[1] >= 3:
        # 3D plot
        ax = plt.axes(projection='3d')

        if cluster_labels is not None:
            scatter = ax.scatter(pca_coords[:, 0], pca_coords[:, 1], pca_coords[:, 2],
                               c=cluster_labels, cmap='tab10', alpha=0.7)
            plt.colorbar(scatter, label='Population Cluster')
        else:
            ax.scatter(pca_coords[:, 0], pca_coords[:, 1], pca_coords[:, 2],
                      alpha=0.7, color='blue')

        ax.set_xlabel('PC1')
        ax.set_ylabel('PC2')
        ax.set_zlabel('PC3')
        plt.title('Population Structure (3D PCA)')

    else:
        # 2D plot
        if cluster_labels is not None:
            scatter = plt.scatter(pca_coords[:, 0], pca_coords[:, 1],
                                c=cluster_labels, cmap='tab10', alpha=0.7)
            plt.colorbar(scatter, label='Population Cluster')
        else:
            plt.scatter(pca_coords[:, 0], pca_coords[:, 1],
                       alpha=0.7, color='blue')

        plt.xlabel('PC1')
        plt.ylabel('PC2')
        plt.title('Population Structure (2D PCA)')
        plt.grid(True, alpha=0.3)

    plt.tight_layout()

    if output_file:
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        logger.info(f"Saved population structure plot to {output_file}")

    return plt.gcf()


def plot_demographic_history(ne_estimates: List[float], generations: List[int],
                           output_file: Optional[str] = None) -> Optional[any]:
    """Plot demographic history (effective population size over time).

    Args:
        ne_estimates: Effective population size estimates
        generations: Generation numbers
        output_file: Optional output file path

    Returns:
        Plot object if matplotlib available, None otherwise

    Example:
        >>> ne_values = [1000, 1200, 800, 1500]
        >>> gens = [0, 1000, 2000, 3000]
        >>> plot = plot_demographic_history(ne_values, gens)
    """
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        logger.warning("matplotlib not available for plotting")
        return None

    plt.figure(figsize=(12, 6))

    plt.plot(generations, ne_estimates, 'o-', color='purple', linewidth=2, alpha=0.8)

    plt.xlabel('Generation')
    plt.ylabel('Effective Population Size (Ne)')
    plt.title('Demographic History')
    plt.grid(True, alpha=0.3)
    plt.yscale('log')  # Often more informative on log scale

    # Add generation markers
    for gen, ne in zip(generations, ne_estimates):
        plt.annotate('.0f', (gen, ne),
                    xytext=(5, 5), textcoords='offset points',
                    fontsize=8, alpha=0.8)

    if output_file:
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        logger.info(f"Saved demographic history plot to {output_file}")

    return plt.gcf()


def create_population_summary_plot(population_data: Dict[str, Dict[str, float]],
                                 output_file: Optional[str] = None) -> Optional[any]:
    """Create comprehensive population summary plot.

    Args:
        population_data: Dictionary with population statistics
        output_file: Optional output file path

    Returns:
        Plot object if matplotlib available, None otherwise

    Example:
        >>> data = {
        ...     "pop1": {"diversity": 0.01, "fst": 0.05},
        ...     "pop2": {"diversity": 0.015, "fst": 0.03}
        ... }
        >>> plot = create_population_summary_plot(data)
    """
    try:
        import matplotlib.pyplot as plt
        import numpy as np
    except ImportError:
        logger.warning("matplotlib not available for plotting")
        return None

    if not population_data:
        return None

    # Extract metrics
    populations = list(population_data.keys())
    metrics = {}

    # Get all available metrics
    all_metrics = set()
    for pop_data in population_data.values():
        all_metrics.update(pop_data.keys())

    # Prepare data for plotting
    for metric in all_metrics:
        metrics[metric] = [population_data[pop].get(metric, 0) for pop in populations]

    n_metrics = len(metrics)
    if n_metrics == 0:
        return None

    # Create subplots
    fig, axes = plt.subplots(n_metrics, 1, figsize=(10, 4 * n_metrics))
    if n_metrics == 1:
        axes = [axes]

    for i, (metric, values) in enumerate(metrics.items()):
        ax = axes[i]
        bars = ax.bar(populations, values, color='lightcoral', alpha=0.7, edgecolor='black')
        ax.set_ylabel(metric.replace('_', ' ').title())
        ax.set_title(f'{metric.replace("_", " ").title()} by Population')
        ax.grid(True, alpha=0.3, axis='y')

        # Add value labels
        for bar, value in zip(bars, values):
            ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + max(values) * 0.01,
                   '.4f', ha='center', va='bottom', fontsize=8)

    plt.tight_layout()

    if output_file:
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        logger.info(f"Saved population summary plot to {output_file}")

    return plt.gcf()


def plot_mutation_spectrum(mutation_counts: Dict[str, int], output_file: Optional[str] = None) -> Optional[any]:
    """Plot mutation spectrum (types of mutations).

    Args:
        mutation_counts: Dictionary mapping mutation types to counts
        output_file: Optional output file path

    Returns:
        Plot object if matplotlib available, None otherwise

    Example:
        >>> mutations = {"A>T": 10, "C>G": 15, "G>A": 8}
        >>> plot = plot_mutation_spectrum(mutations)
    """
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        logger.warning("matplotlib not available for plotting")
        return None

    if not mutation_counts:
        return None

    plt.figure(figsize=(10, 6))

    mutation_types = list(mutation_counts.keys())
    counts = list(mutation_counts.values())

    bars = plt.bar(mutation_types, counts, color='lightgreen', alpha=0.7, edgecolor='black')
    plt.ylabel('Count')
    plt.title('Mutation Spectrum')
    plt.xticks(rotation=45)

    # Add value labels
    for bar, count in zip(bars, counts):
        plt.text(bar.get_x() + bar.get_width()/2, bar.get_height() + max(counts) * 0.01,
                str(count), ha='center', va='bottom')

    plt.grid(True, alpha=0.3, axis='y')

    if output_file:
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        logger.info(f"Saved mutation spectrum plot to {output_file}")

    return plt.gcf()
