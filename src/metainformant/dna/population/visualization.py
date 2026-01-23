"""Population genetics visualization utilities.

This module provides plotting and visualization functions for population
genetics data including F_ST heatmaps, selection plots, and demographic visualizations.
"""

from __future__ import annotations

from typing import Dict, List, Optional, Tuple

import numpy as np

from metainformant.core import logging

logger = logging.get_logger(__name__)


def plot_fst_matrix(
    data: Dict[str, List[str]] | np.ndarray | List[List[float]],
    pop_names: Optional[List[str]] = None,
    output_file: Optional[str] = None,
) -> Optional[any]:
    """Create F_ST matrix heatmap visualization.

    Args:
        data: Dictionary of sequences (pop->seqs) OR F_ST matrix (array/list)
        pop_names: List of population names (required if data is matrix)
        output_file: Optional output file path

    Returns:
        Plot object if matplotlib available, None otherwise
    """
    try:
        import matplotlib.pyplot as plt
        import seaborn as sns
    except ImportError:
        logger.warning("matplotlib and/or seaborn not available for plotting")
        return None

    if isinstance(data, dict):
        # Calculate from sequences
        populations = data
        pop_names_list = list(populations.keys())
        n_pops = len(pop_names_list)
        fst_matrix = np.zeros((n_pops, n_pops))

        for i in range(n_pops):
            for j in range(i + 1, n_pops):
                pop1_seqs = populations[pop_names_list[i]]
                pop2_seqs = populations[pop_names_list[j]]
                from .analysis import calculate_fst

                fst = calculate_fst(pop1_seqs, pop2_seqs)
                fst_matrix[i, j] = fst
                fst_matrix[j, i] = fst

        plot_names = pop_names_list
    else:
        # Use provided matrix
        fst_matrix = np.array(data)
        plot_names = pop_names if pop_names else [f"Pop{i+1}" for i in range(len(fst_matrix))]

    # Create heatmap
    plt.figure(figsize=(8, 6))
    sns.heatmap(fst_matrix, annot=True, cmap="YlOrRd", square=True, xticklabels=plot_names, yticklabels=plot_names)
    plt.title("F_ST Matrix")
    plt.tight_layout()

    if output_file:
        plt.savefig(output_file, dpi=300, bbox_inches="tight")
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
    plt.hist(tajima_d_values, bins=30, alpha=0.7, color="skyblue", edgecolor="black")
    plt.axvline(x=0, color="red", linestyle="--", alpha=0.7, label="Neutral expectation")

    plt.xlabel("Tajima's D")
    plt.ylabel("Frequency")
    plt.title("Distribution of Tajima's D Values")
    plt.legend()
    plt.grid(True, alpha=0.3)

    if output_file:
        plt.savefig(output_file, dpi=300, bbox_inches="tight")
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
        plt.plot(values, "o-", alpha=0.7)
        plt.axhline(y=0, color="red", linestyle="--", alpha=0.5)
        plt.ylabel(stat_name)
        plt.title(f"{stat_name} along sequence")
        plt.grid(True, alpha=0.3)

    plt.tight_layout()

    if output_file:
        plt.savefig(output_file, dpi=300, bbox_inches="tight")
        logger.info(f"Saved selection statistics plot to {output_file}")

    return plt.gcf()


def plot_population_diversity(
    diversity_data: Dict[str, float], output_file: Optional[str] = None, output_path: Optional[str] = None
) -> Optional[any]:
    """Plot population diversity comparison.

    Args:
        diversity_data: Dictionary mapping population names to diversity values
        output_file: Optional output file path
        output_path: Alias for output_file

    Returns:
        Plot object if matplotlib available, None otherwise

    Example:
        >>> diversity = {"pop1": 0.01, "pop2": 0.015, "pop3": 0.008}
        >>> plot = plot_population_diversity(diversity)
    """
    output = output_file or output_path
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        logger.warning("matplotlib not available for plotting")
        return None

    plt.figure(figsize=(10, 6))

    populations = list(diversity_data.keys())
    diversities = list(diversity_data.values())

    bars = plt.bar(populations, diversities, color="lightblue", alpha=0.7, edgecolor="black")
    plt.ylabel("Nucleotide Diversity (π)")
    plt.title("Population Diversity Comparison")
    plt.xticks(rotation=45)

    # Add value labels on bars
    for bar, value in zip(bars, diversities):
        plt.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.001, ".4f", ha="center", va="bottom")

    plt.grid(True, alpha=0.3, axis="y")

    if output_file:
        plt.savefig(output_file, dpi=300, bbox_inches="tight")
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
    plt.plot(distances, ld_values, "o-", color="darkgreen", alpha=0.7, markersize=4)

    # Add exponential decay fit if enough points
    if len(distances) > 3:
        try:
            from scipy.optimize import curve_fit

            def exp_decay(x, a, b):
                return a * np.exp(-b * x)

            popt, _ = curve_fit(exp_decay, distances, ld_values, p0=[1, 0.1])
            x_fit = np.linspace(min(distances), max(distances), 100)
            y_fit = exp_decay(x_fit, *popt)
            plt.plot(x_fit, y_fit, "--", color="red", alpha=0.7, label=".3f")
        except ImportError:
            pass  # scipy not available

    plt.xlabel("Distance (bp)")
    plt.ylabel("Linkage Disequilibrium (r²)")
    plt.title("LD Decay with Distance")
    plt.grid(True, alpha=0.3)
    plt.legend()

    if output_file:
        plt.savefig(output_file, dpi=300, bbox_inches="tight")
        logger.info(f"Saved LD decay plot to {output_file}")

    return plt.gcf()


def plot_population_structure(
    pca_coords: np.ndarray, cluster_labels: Optional[List[int]] = None, output_file: Optional[str] = None
) -> Optional[any]:
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
        ax = plt.axes(projection="3d")

        if cluster_labels is not None:
            scatter = ax.scatter(
                pca_coords[:, 0], pca_coords[:, 1], pca_coords[:, 2], c=cluster_labels, cmap="tab10", alpha=0.7
            )
            plt.colorbar(scatter, label="Population Cluster")
        else:
            ax.scatter(pca_coords[:, 0], pca_coords[:, 1], pca_coords[:, 2], alpha=0.7, color="blue")

        ax.set_xlabel("PC1")
        ax.set_ylabel("PC2")
        ax.set_zlabel("PC3")
        plt.title("Population Structure (3D PCA)")

    else:
        # 2D plot
        if cluster_labels is not None:
            scatter = plt.scatter(pca_coords[:, 0], pca_coords[:, 1], c=cluster_labels, cmap="tab10", alpha=0.7)
            plt.colorbar(scatter, label="Population Cluster")
        else:
            plt.scatter(pca_coords[:, 0], pca_coords[:, 1], alpha=0.7, color="blue")

        plt.xlabel("PC1")
        plt.ylabel("PC2")
        plt.title("Population Structure (2D PCA)")
        plt.grid(True, alpha=0.3)

    plt.tight_layout()

    if output_file:
        plt.savefig(output_file, dpi=300, bbox_inches="tight")
        logger.info(f"Saved population structure plot to {output_file}")

    return plt.gcf()


def plot_demographic_history(
    ne_estimates: List[float], generations: List[int], output_file: Optional[str] = None
) -> Optional[any]:
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

    plt.plot(generations, ne_estimates, "o-", color="purple", linewidth=2, alpha=0.8)

    plt.xlabel("Generation")
    plt.ylabel("Effective Population Size (Ne)")
    plt.title("Demographic History")
    plt.grid(True, alpha=0.3)
    plt.yscale("log")  # Often more informative on log scale

    # Add generation markers
    for gen, ne in zip(generations, ne_estimates):
        plt.annotate(".0f", (gen, ne), xytext=(5, 5), textcoords="offset points", fontsize=8, alpha=0.8)

    if output_file:
        plt.savefig(output_file, dpi=300, bbox_inches="tight")
        logger.info(f"Saved demographic history plot to {output_file}")

    return plt.gcf()


def create_population_summary_plot(
    population_data: Dict[str, Dict[str, float]], output_file: Optional[str] = None
) -> Optional[any]:
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
        bars = ax.bar(populations, values, color="lightcoral", alpha=0.7, edgecolor="black")
        ax.set_ylabel(metric.replace("_", " ").title())
        ax.set_title(f'{metric.replace("_", " ").title()} by Population')
        ax.grid(True, alpha=0.3, axis="y")

        # Add value labels
        for bar, value in zip(bars, values):
            ax.text(
                bar.get_x() + bar.get_width() / 2,
                bar.get_height() + max(values) * 0.01,
                ".4f",
                ha="center",
                va="bottom",
                fontsize=8,
            )

    plt.tight_layout()

    if output_file:
        plt.savefig(output_file, dpi=300, bbox_inches="tight")
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

    bars = plt.bar(mutation_types, counts, color="lightgreen", alpha=0.7, edgecolor="black")
    plt.ylabel("Count")
    plt.title("Mutation Spectrum")
    plt.xticks(rotation=45)

    # Add value labels
    for bar, count in zip(bars, counts):
        plt.text(
            bar.get_x() + bar.get_width() / 2,
            bar.get_height() + max(counts) * 0.01,
            str(count),
            ha="center",
            va="bottom",
        )

    plt.grid(True, alpha=0.3, axis="y")

    if output_file:
        plt.savefig(output_file, dpi=300, bbox_inches="tight")
        logger.info(f"Saved mutation spectrum plot to {output_file}")

    return plt.gcf()


def plot_allele_frequency_spectrum(allele_frequencies: List[float], output_file: Optional[str] = None) -> Optional[any]:
    """Plot allele frequency spectrum.

    Args:
        allele_frequencies: List of allele frequencies (0.0 to 1.0)
        output_file: Optional output file path

    Returns:
        Plot object if matplotlib available, None otherwise

    Example:
        >>> freqs = [0.1, 0.2, 0.15, 0.3, 0.8, 0.9]
        >>> plot = plot_allele_frequency_spectrum(freqs)
    """
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        logger.warning("matplotlib not available for plotting")
        return None

    if not allele_frequencies:
        return None

    plt.figure(figsize=(10, 6))

    # Create histogram bins for frequency classes
    bins = np.linspace(0, 1, 11)  # 10 bins from 0 to 1
    plt.hist(allele_frequencies, bins=bins, alpha=0.7, color="skyblue", edgecolor="black")

    plt.xlabel("Allele Frequency")
    plt.ylabel("Number of Variants")
    plt.title("Allele Frequency Spectrum")
    plt.grid(True, alpha=0.3, axis="y")

    # Add bin labels
    bin_centers = (bins[:-1] + bins[1:]) / 2
    plt.xticks(bin_centers, [f"{center:.1f}" for center in bin_centers], rotation=45)

    if output_file:
        plt.savefig(output_file, dpi=300, bbox_inches="tight")
        logger.info(f"Saved allele frequency spectrum plot to {output_file}")

    return plt.gcf()


def plot_bootstrap_distribution(
    bootstrap_values: List[float], observed_value: float | None = None, output_file: Optional[str] = None
) -> Optional[any]:
    """Plot bootstrap distribution with confidence intervals.

    Args:
        bootstrap_values: List of bootstrap replicate values
        observed_value: Observed value to mark on the plot
        output_file: Optional output file path

    Returns:
        Plot object if matplotlib available, None otherwise

    Example:
        >>> values = [1.2, 1.5, 0.8, 1.1, 1.3, 0.9]
        >>> plot = plot_bootstrap_distribution(values, observed_value=1.0)
    """
    try:
        import matplotlib.pyplot as plt
        from scipy import stats
    except ImportError:
        logger.warning("matplotlib and/or scipy not available for plotting")
        return None

    if not bootstrap_values:
        return None

    plt.figure(figsize=(10, 6))

    # Plot histogram
    plt.hist(bootstrap_values, bins=30, alpha=0.7, color="lightblue", edgecolor="black", density=True)

    # Add kernel density estimate
    try:
        density = stats.gaussian_kde(bootstrap_values)
        x_vals = np.linspace(min(bootstrap_values), max(bootstrap_values), 200)
        plt.plot(x_vals, density(x_vals), "r-", linewidth=2, label="KDE")
    except Exception:
        pass  # Skip KDE if it fails

    # Mark observed value
    if observed_value is not None:
        plt.axvline(observed_value, color="red", linestyle="--", linewidth=2, label=f"Observed: {observed_value:.3f}")

    # Add confidence intervals
    ci_lower, ci_upper = np.percentile(bootstrap_values, [2.5, 97.5])
    plt.axvline(ci_lower, color="green", linestyle=":", linewidth=2, label=f"95% CI: [{ci_lower:.3f}, {ci_upper:.3f}]")
    plt.axvline(ci_upper, color="green", linestyle=":", linewidth=2)

    # Add mean and median
    mean_val = np.mean(bootstrap_values)
    median_val = np.median(bootstrap_values)
    plt.axvline(mean_val, color="orange", linestyle="-.", linewidth=2, label=f"Mean: {mean_val:.3f}")
    plt.axvline(median_val, color="purple", linestyle="-.", linewidth=2, label=f"Median: {median_val:.3f}")

    plt.xlabel("Bootstrap Values")
    plt.ylabel("Density")
    plt.title("Bootstrap Distribution")
    plt.legend()
    plt.grid(True, alpha=0.3, axis="y")

    if output_file:
        plt.savefig(output_file, dpi=300, bbox_inches="tight")
        logger.info(f"Saved bootstrap distribution plot to {output_file}")

    return plt.gcf()


def plot_demographic_comparison(
    pop1_data: Dict[str, Any],
    pop2_data: Optional[Dict[str, Any]] = None,
    output_file: Optional[str] = None,
    output_path: Optional[str] = None,
) -> Optional[any]:
    """Compare demographic histories between populations.

    Args:
        pop1_data: Demographic data for population 1
        pop2_data: Demographic data for population 2
        output_file: Optional output file path
        output_path: Alias for output_file

    Returns:
        Plot object if matplotlib available, None otherwise
    """
    output = output_file or output_path
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        logger.warning("matplotlib not available for plotting")
        return None

    # Placeholder implementation - would need actual demographic data
    fig, axes = plt.subplots(1, 2, figsize=(12, 6))
    fig.suptitle("Population Demographic Comparison")
    axes[0].text(0.5, 0.5, "Demographic Comparison\n(Not implemented)", ha="center", va="center")

    if output:
        plt.savefig(output, dpi=300, bbox_inches="tight")
        logger.info(f"Saved demographic comparison plot to {output}")

    return plt.gcf()


def plot_diversity_comparison(
    diversity_data: Dict[str, float], output_file: Optional[str] = None, output_path: Optional[str] = None
) -> Optional[any]:
    """Compare diversity metrics across populations.

    Args:
        diversity_data: Dictionary mapping population names to diversity values
        output_file: Optional output file path
        output_path: Alias for output_file

    Returns:
        Plot object if matplotlib available, None otherwise
    """
    output = output_file or output_path
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        logger.warning("matplotlib not available for plotting")
        return None

    if not diversity_data:
        return None

    plt.figure(figsize=(10, 6))

    populations = list(diversity_data.keys())
    diversities = list(diversity_data.values())

    bars = plt.bar(populations, diversities, color="skyblue", alpha=0.7, edgecolor="black")

    plt.ylabel("Diversity Index")
    plt.title("Population Diversity Comparison")
    plt.xticks(rotation=45)

    # Add value labels
    for bar, diversity in zip(bars, diversities):
        plt.text(
            bar.get_x() + bar.get_width() / 2,
            bar.get_height() + max(diversities) * 0.01,
            f"{diversity:.3f}",
            ha="center",
            va="bottom",
        )

    plt.grid(True, alpha=0.3, axis="y")

    if output:
        plt.savefig(output, dpi=300, bbox_inches="tight")
        logger.info(f"Saved diversity comparison plot to {output}")

    return plt.gcf()


def plot_fst_comparison(
    fst_data: Dict[str, float], output_file: Optional[str] = None, output_path: Optional[str] = None
) -> Optional[any]:
    """Compare F_ST values across loci or populations.

    Args:
        fst_data: Dictionary mapping locus/population names to F_ST values
        output_file: Optional output file path
        output_path: Alias for output_file

    Returns:
        Plot object if matplotlib available, None otherwise
    """
    output = output_file or output_path
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        logger.warning("matplotlib not available for plotting")
        return None

    if not fst_data:
        return None

    plt.figure(figsize=(12, 6))

    loci = list(fst_data.keys())
    fst_values = list(fst_data.values())

    plt.bar(loci, fst_values, color="lightcoral", alpha=0.7, edgecolor="black")

    plt.ylabel("F_ST")
    plt.title("F_ST Comparison Across Loci")
    plt.xticks(rotation=45)
    plt.grid(True, alpha=0.3, axis="y")

    # Add significance line at 0.05
    plt.axhline(y=0.05, color="red", linestyle="--", alpha=0.7, label="Significance threshold (0.05)")
    plt.legend()

    if output:
        plt.savefig(output, dpi=300, bbox_inches="tight")
        logger.info(f"Saved F_ST comparison plot to {output}")

    return plt.gcf()


def plot_hardy_weinberg_test(results: List[Dict[str, Any]], output_file: Optional[str] = None) -> Optional[any]:
    """Plot Hardy-Weinberg equilibrium test results.

    Args:
        results: List of HWE test results for each locus
        output_file: Optional output file path

    Returns:
        Plot object if matplotlib available, None otherwise
    """
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        logger.warning("matplotlib not available for plotting")
        return None

    if not results:
        return None

    plt.figure(figsize=(12, 6))

    loci = [r.get("locus", f"Locus_{i}") for i, r in enumerate(results)]
    p_values = [r.get("p_value", 1.0) for r in results]

    plt.bar(loci, [-np.log10(p) for p in p_values], color="lightgreen", alpha=0.7, edgecolor="black")

    plt.ylabel("-log10(p-value)")
    plt.title("Hardy-Weinberg Equilibrium Test Results")
    plt.xticks(rotation=45)

    # Add significance line
    plt.axhline(y=-np.log10(0.05), color="red", linestyle="--", alpha=0.7, label="Significance threshold (p=0.05)")
    plt.legend()
    plt.grid(True, alpha=0.3, axis="y")

    if output_file:
        plt.savefig(output_file, dpi=300, bbox_inches="tight")
        logger.info(f"Saved HWE test plot to {output_file}")

    return plt.gcf()


def plot_heterozygosity_distribution(
    het_data: Dict[str, List[float]] | List[float], output_file: Optional[str] = None
) -> Optional[any]:
    """Plot distribution of heterozygosity across loci.

    Args:
        het_data: Dictionary mapping population names to heterozygosity values OR list of values
        output_file: Optional output file path

    Returns:
        Plot object if matplotlib available, None otherwise
    """
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        logger.warning("matplotlib not available for plotting")
        return None

    if not het_data:
        return None

    plt.figure(figsize=(10, 6))

    if isinstance(het_data, dict):
        for pop_name, het_values in het_data.items():
            plt.hist(het_values, alpha=0.7, label=pop_name, bins=20)
        plt.legend()
    else:
        # List of values
        plt.hist(het_data, alpha=0.7, bins=20, color="skyblue", edgecolor="black")

    plt.xlabel("Heterozygosity")
    plt.ylabel("Frequency")
    plt.title("Heterozygosity Distribution")
    plt.grid(True, alpha=0.3)

    if output_file:
        plt.savefig(output_file, dpi=300, bbox_inches="tight")
        logger.info(f"Saved heterozygosity distribution plot to {output_file}")

    return plt.gcf()


def plot_kinship_matrix(
    kinship_matrix: np.ndarray | Dict[str, Any],
    sample_names: Optional[List[str]] = None,
    output_file: Optional[str] = None,
    output_path: Optional[str] = None,
    max_samples: Optional[int] = None,
) -> Optional[any]:
    """Plot kinship matrix as heatmap.

    Args:
        kinship_matrix: Kinship coefficient matrix (or result dict)
        sample_names: Optional sample names for axis labels
        output_file: Optional output file path
        output_path: Alias for output_file
        max_samples: Optional limit on number of samples

    Returns:
        Plot object if matplotlib available, None otherwise
    """
    output = output_file or output_path

    if isinstance(kinship_matrix, dict) and kinship_matrix.get("status") == "failed":
        from metainformant.core.utils.errors import ValidationError

        raise ValidationError("Kinship result status is not 'success'")

    if isinstance(kinship_matrix, dict) and "kinship_matrix" in kinship_matrix:
        kinship_matrix = np.array(kinship_matrix["kinship_matrix"])
        if max_samples and len(kinship_matrix) > max_samples:
            kinship_matrix = kinship_matrix[:max_samples, :max_samples]
    try:
        import matplotlib.pyplot as plt
        import seaborn as sns
    except ImportError:
        logger.warning("matplotlib and/or seaborn not available for plotting")
        return None

    plt.figure(figsize=(10, 8))

    # Create labels
    if sample_names:
        labels = sample_names
    else:
        labels = [f"Sample_{i}" for i in range(len(kinship_matrix))]

    sns.heatmap(kinship_matrix, annot=False, cmap="YlOrRd", square=True, xticklabels=labels, yticklabels=labels)

    plt.title("Kinship Matrix")
    plt.xticks(rotation=45, ha="right")
    plt.yticks(rotation=0)

    if output:
        plt.savefig(output, dpi=300, bbox_inches="tight")
        logger.info(f"Saved kinship matrix plot to {output}")

    return plt.gcf()


# Additional placeholder functions for completeness
def plot_linkage_disequilibrium_decay(
    ld_data: List[Tuple[int, float]] | List[float],
    distances: Optional[List[float]] = None,
    output_file: Optional[str] = None,
    output_path: Optional[str] = None,
) -> Optional[any]:
    """Plot linkage disequilibrium decay with distance."""
    output = output_file or output_path
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        logger.warning("matplotlib not available for plotting")
        return None

    plt.figure(figsize=(10, 6))

    # Check if we have data to plot
    if distances is not None and isinstance(ld_data, list):
        plt.plot(distances, ld_data, "o-", alpha=0.7)
    elif ld_data and isinstance(ld_data[0], tuple):
        # list of tuples
        dists, values = zip(*ld_data)
        plt.plot(dists, values, "o-", alpha=0.7)
    else:
        plt.text(0.5, 0.5, "LD Decay Plot\n(Not fully implemented)", ha="center", va="center")

    plt.title("Linkage Disequilibrium Decay")

    if output:
        plt.savefig(output, dpi=300, bbox_inches="tight")

    return plt.gcf()


def plot_neutrality_test_suite(
    results: Dict[str, Any], output_file: Optional[str] = None, output_path: Optional[str] = None
) -> Optional[any]:
    """Plot comprehensive neutrality test results."""
    output = output_file or output_path
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        logger.warning("matplotlib not available for plotting")
        return None

    plt.figure(figsize=(12, 8))
    plt.text(0.5, 0.5, "Neutrality Test Suite\n(Not fully implemented)", ha="center", va="center")
    plt.title("Neutrality Test Results")

    if output_file:
        plt.savefig(output_file, dpi=300, bbox_inches="tight")

    return plt.gcf()


def plot_neutrality_test_summary(
    summary: Dict[str, Any], output_file: Optional[str] = None, output_path: Optional[str] = None
) -> Optional[any]:
    """Plot summary of neutrality test interpretations."""
    output = output_file or output_path
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        logger.warning("matplotlib not available for plotting")
        return None

    fig, axes = plt.subplots(2, 2, figsize=(10, 8))
    fig.suptitle("Neutrality Test Summary")
    # Use axes so they exist
    axes[0, 0].text(0.5, 0.5, "Neutrality Test Summary\n(Not fully implemented)", ha="center", va="center")

    if output:
        plt.savefig(output, dpi=300, bbox_inches="tight")

    return plt.gcf()


def plot_outlier_detection(
    results: List[Dict[str, Any]] | List[float],
    outliers: Optional[List[int]] = None,
    output_file: Optional[str] = None,
    output_path: Optional[str] = None,
) -> Optional[any]:
    """Plot outlier detection results."""
    output = output_file or output_path
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        logger.warning("matplotlib not available for plotting")
        return None

    plt.figure(figsize=(10, 6))
    plt.text(0.5, 0.5, "Outlier Detection\n(Not fully implemented)", ha="center", va="center")
    plt.title("Outlier Detection Results")

    if output_file:
        plt.savefig(output_file, dpi=300, bbox_inches="tight")

    return plt.gcf()


def plot_pca_results(
    pca_result: Tuple[np.ndarray, np.ndarray, np.ndarray] | Dict[str, Any],
    sample_names: Optional[List[str]] = None,
    output_file: Optional[str] = None,
    output_path: Optional[str] = None,
    n_components: Optional[int] = None,
) -> Optional[any]:
    """Plot PCA results."""
    output = output_file or output_path

    if isinstance(pca_result, dict) and pca_result.get("status") == "failed":
        from metainformant.core.utils.errors import ValidationError

        raise ValidationError("PCA result status is not 'success'")

    try:
        import matplotlib.pyplot as plt
    except ImportError:
        logger.warning("matplotlib not available for plotting")
        return None

    fig, axes = plt.subplots(1, 3, figsize=(18, 6))
    fig.suptitle("Principal Component Analysis")

    # Placeholder for actual plotting logic or if real data passed
    axes[0].text(0.5, 0.5, "PC1 vs PC2", ha="center", va="center")
    axes[1].text(0.5, 0.5, "PC2 vs PC3", ha="center", va="center")
    axes[2].text(0.5, 0.5, "Explained Variance", ha="center", va="center")

    if output:
        plt.savefig(output, dpi=300, bbox_inches="tight")

    return plt.gcf()


def plot_permutation_test(
    results: Dict[str, Any] | List[float],
    observed_value: Optional[float] = None,
    p_value: Optional[float] = None,
    output_file: Optional[str] = None,
    output_path: Optional[str] = None,
) -> Optional[any]:
    """Plot permutation test results."""
    output = output_file or output_path
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        logger.warning("matplotlib not available for plotting")
        return None

    plt.figure(figsize=(10, 6))
    plt.text(0.5, 0.5, "Permutation Test\n(Not fully implemented)", ha="center", va="center")
    plt.title("Permutation Test Results")

    if output_file:
        plt.savefig(output_file, dpi=300, bbox_inches="tight")

    return plt.gcf()


def plot_pi_vs_theta(
    pi_values: List[float],
    theta_values: List[float],
    locus_names: Optional[List[str]] = None,
    output_file: Optional[str] = None,
) -> Optional[any]:
    """Plot π vs θ comparison."""
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        logger.warning("matplotlib not available for plotting")
        return None

    plt.figure(figsize=(8, 8))

    plt.scatter(theta_values, pi_values, alpha=0.7, color="blue", edgecolors="black")

    # Add diagonal line
    max_val = max(max(pi_values), max(theta_values))
    plt.plot([0, max_val], [0, max_val], "r--", alpha=0.7, label="π = θ")

    plt.xlabel("Watterson's θ")
    plt.ylabel("Nucleotide diversity (π)")
    plt.title("π vs θ Comparison")
    plt.legend()
    plt.grid(True, alpha=0.3)

    if output_file:
        plt.savefig(output_file, dpi=300, bbox_inches="tight")
        logger.info(f"Saved π vs θ plot to {output_file}")

    return plt.gcf()


def plot_site_frequency_spectrum(
    sfs_data: List[int], output_file: Optional[str] = None, output_path: Optional[str] = None
) -> Optional[any]:
    """Plot site frequency spectrum."""
    output = output_file or output_path
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        logger.warning("matplotlib not available for plotting")
        return None

    if not sfs_data:
        return None

    plt.figure(figsize=(10, 6))

    # SFS data represents counts at each frequency
    frequencies = list(range(1, len(sfs_data) + 1))
    counts = sfs_data

    plt.bar(frequencies, counts, color="lightblue", alpha=0.7, edgecolor="black")

    plt.xlabel("Allele Frequency")
    plt.ylabel("Number of Sites")
    plt.title("Site Frequency Spectrum")
    plt.grid(True, alpha=0.3, axis="y")

    if output:
        plt.savefig(output, dpi=300, bbox_inches="tight")
        logger.info(f"Saved SFS plot to {output}")

    return plt.gcf()


def plot_statistic_correlation_matrix(
    stats_data: Dict[str, List[float]], output_file: Optional[str] = None
) -> Optional[any]:
    """Plot correlation matrix of population statistics."""
    try:
        import matplotlib.pyplot as plt
        import seaborn as sns
    except ImportError:
        logger.warning("matplotlib and/or seaborn not available for plotting")
        return None

    if not stats_data:
        return None

    plt.figure(figsize=(10, 8))

    # Create correlation matrix
    import pandas as pd

    df = pd.DataFrame(stats_data)
    corr_matrix = df.corr()

    sns.heatmap(corr_matrix, annot=True, cmap="coolwarm", center=0, square=True, cbar_kws={"shrink": 0.8})

    plt.title("Population Statistics Correlation Matrix")

    if output_file:
        plt.savefig(output_file, dpi=300, bbox_inches="tight")
        logger.info(f"Saved correlation matrix plot to {output_file}")

    return plt.gcf()


def plot_statistic_distribution(
    stat_values: List[float] | Dict[str, List[float]],
    stat_name: str = "Statistic",
    output_file: Optional[str] = None,
    plot_type: str = "histogram",
) -> Optional[any]:
    """Plot distribution of a population statistic."""
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        logger.warning("matplotlib not available for plotting")
        return None

    if not stat_values:
        return None

    plt.figure(figsize=(10, 6))

    if isinstance(stat_values, dict):
        for name, values in stat_values.items():
            plt.hist(values, bins=20, alpha=0.5, label=name)
        plt.legend()
    else:
        plt.hist(stat_values, bins=20, alpha=0.7, color="lightgreen", edgecolor="black")

        # Add statistics only for single list
        mean_val = np.mean(stat_values)
        median_val = np.median(stat_values)
        plt.axvline(mean_val, color="red", linestyle="--", alpha=0.8, label=f"Mean: {mean_val:.3f}")
        plt.axvline(median_val, color="blue", linestyle="--", alpha=0.8, label=f"Median: {median_val:.3f}")
        plt.legend()

    plt.xlabel(stat_name)
    plt.ylabel("Frequency")
    plt.title(f"{stat_name} Distribution")
    plt.grid(True, alpha=0.3, axis="y")

    if output_file:
        plt.savefig(output_file, dpi=300, bbox_inches="tight")
        logger.info(f"Saved {stat_name} distribution plot to {output_file}")

    return plt.gcf()


def plot_summary_statistics_grid(
    stats_dict: Dict[str, Dict[str, float]], output_file: Optional[str] = None, output_path: Optional[str] = None
) -> Optional[any]:
    """Plot grid of summary statistics for multiple populations."""
    output = output_file or output_path
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        logger.warning("matplotlib not available for plotting")
        return None

    if not stats_dict:
        return None

    fig, axes = plt.subplots(3, 2, figsize=(15, 12))
    fig.suptitle("Population Summary Statistics")
    axes[0, 0].text(0.5, 0.5, "Summary Statistics Grid\n(Not fully implemented)", ha="center", va="center")

    if output:
        plt.savefig(output, dpi=300, bbox_inches="tight")

    return plt.gcf()


def plot_tajimas_d_comparison(
    tajima_d_values: Dict[str, float], output_file: Optional[str] = None, output_path: Optional[str] = None
) -> Optional[any]:
    """Compare Tajima's D values across populations."""
    output = output_file or output_path
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        logger.warning("matplotlib not available for plotting")
        return None

    if not tajima_d_values:
        return None

    plt.figure(figsize=(10, 6))

    populations = list(tajima_d_values.keys())
    d_values = list(tajima_d_values.values())

    colors = ["red" if d < -2 else "green" if d > 2 else "blue" for d in d_values]

    bars = plt.bar(populations, d_values, color=colors, alpha=0.7, edgecolor="black")

    plt.ylabel("Tajima's D")
    plt.title("Tajima's D Comparison Across Populations")
    plt.xticks(rotation=45)

    # Add reference lines
    plt.axhline(y=0, color="black", linestyle="-", alpha=0.5)
    plt.axhline(y=2, color="green", linestyle="--", alpha=0.7, label="Balancing selection")
    plt.axhline(y=-2, color="red", linestyle="--", alpha=0.7, label="Positive selection")
    plt.legend()

    plt.grid(True, alpha=0.3, axis="y")

    if output:
        plt.savefig(output, dpi=300, bbox_inches="tight")
        logger.info(f"Saved Tajima's D comparison plot to {output}")

    return plt.gcf()
