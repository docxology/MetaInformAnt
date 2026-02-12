"""Population genetics visualization utilities - statistical comparison plots.

This module provides statistical comparison and advanced visualization functions
for population genetics data including demographic comparisons, diversity
comparisons, F_ST comparisons, Hardy-Weinberg tests, heterozygosity, kinship
matrices, LD decay, neutrality tests, outlier detection, PCA results,
permutation tests, pi vs theta, site frequency spectrum, statistic correlation,
statistic distribution, summary statistics grid, and Tajima's D comparisons.
"""

from __future__ import annotations

from typing import Any, Dict, List, Optional, Tuple

import numpy as np

from metainformant.core.utils import logging

logger = logging.get_logger(__name__)


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
    """Plot pi vs theta comparison."""
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        logger.warning("matplotlib not available for plotting")
        return None

    plt.figure(figsize=(8, 8))

    plt.scatter(theta_values, pi_values, alpha=0.7, color="blue", edgecolors="black")

    # Add diagonal line
    max_val = max(max(pi_values), max(theta_values))
    plt.plot([0, max_val], [0, max_val], "r--", alpha=0.7, label="pi = theta")

    plt.xlabel("Watterson's theta")
    plt.ylabel("Nucleotide diversity (pi)")
    plt.title("pi vs theta Comparison")
    plt.legend()
    plt.grid(True, alpha=0.3)

    if output_file:
        plt.savefig(output_file, dpi=300, bbox_inches="tight")
        logger.info(f"Saved pi vs theta plot to {output_file}")

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
