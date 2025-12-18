"""Population genetics visualization functions."""

from __future__ import annotations

from collections.abc import Sequence
from pathlib import Path
from typing import Any

import matplotlib
import matplotlib.pyplot as plt
import numpy as np

matplotlib.use("Agg", force=True)

try:
    import seaborn as sns

    SEABORN_AVAILABLE = True
except ImportError:
    SEABORN_AVAILABLE = False
    sns = None


def plot_diversity_comparison(
    diversity_values: dict[str, float],
    *,
    output_path: str | Path | None = None,
    title: str = "Nucleotide Diversity Comparison",
    figsize: tuple[float, float] = (10, 6),
) -> plt.Figure:
    """Plot comparison of nucleotide diversity across scenarios.
    
    Args:
        diversity_values: Dictionary mapping scenario names to π values
        output_path: Optional path to save figure
        title: Plot title
        figsize: Figure size
    
    Returns:
        matplotlib Figure object
    """
    fig, ax = plt.subplots(figsize=figsize)
    
    scenarios = list(diversity_values.keys())
    values = list(diversity_values.values())
    
    # Create bar plot
    bars = ax.bar(range(len(scenarios)), values, color="steelblue", alpha=0.7)
    
    # Add value labels on bars
    for i, (bar, val) in enumerate(zip(bars, values)):
        height = bar.get_height()
        ax.text(
            bar.get_x() + bar.get_width() / 2.0,
            height,
            f"{val:.4f}",
            ha="center",
            va="bottom",
            fontsize=9,
        )
    
    ax.set_xlabel("Scenario", fontsize=12)
    ax.set_ylabel("Nucleotide Diversity (π)", fontsize=12)
    ax.set_title(title, fontsize=14, fontweight="bold")
    ax.set_xticks(range(len(scenarios)))
    ax.set_xticklabels(scenarios, rotation=45, ha="right")
    ax.grid(axis="y", alpha=0.3)
    
    plt.tight_layout()
    
    if output_path:
        fig.savefig(output_path, dpi=300, bbox_inches="tight")
    
    return fig


def plot_tajimas_d_comparison(
    tajimas_d_values: dict[str, float],
    *,
    output_path: str | Path | None = None,
    title: str = "Tajima's D Comparison",
    figsize: tuple[float, float] = (10, 6),
) -> plt.Figure:
    """Plot comparison of Tajima's D across scenarios.
    
    Args:
        tajimas_d_values: Dictionary mapping scenario names to Tajima's D values
        output_path: Optional path to save figure
        title: Plot title
        figsize: Figure size
    
    Returns:
        matplotlib Figure object
    """
    fig, ax = plt.subplots(figsize=figsize)
    
    scenarios = list(tajimas_d_values.keys())
    values = list(tajimas_d_values.values())
    
    # Color bars based on sign
    colors = ["red" if v < 0 else "blue" if v > 0 else "gray" for v in values]
    
    bars = ax.bar(range(len(scenarios)), values, color=colors, alpha=0.7)
    
    # Add value labels
    for i, (bar, val) in enumerate(zip(bars, values)):
        height = bar.get_height()
        ax.text(
            bar.get_x() + bar.get_width() / 2.0,
            height + (0.1 if height >= 0 else -0.15),
            f"{val:.3f}",
            ha="center",
            va="bottom" if height >= 0 else "top",
            fontsize=9,
        )
    
    # Add reference lines
    ax.axhline(y=0, color="black", linestyle="--", linewidth=1, alpha=0.5)
    ax.axhline(y=-2, color="red", linestyle=":", linewidth=1, alpha=0.5, label="Strong negative")
    ax.axhline(y=2, color="blue", linestyle=":", linewidth=1, alpha=0.5, label="Strong positive")
    
    ax.set_xlabel("Scenario", fontsize=12)
    ax.set_ylabel("Tajima's D", fontsize=12)
    ax.set_title(title, fontsize=14, fontweight="bold")
    ax.set_xticks(range(len(scenarios)))
    ax.set_xticklabels(scenarios, rotation=45, ha="right")
    ax.grid(axis="y", alpha=0.3)
    ax.legend(fontsize=9)
    
    plt.tight_layout()
    
    if output_path:
        fig.savefig(output_path, dpi=300, bbox_inches="tight")
    
    return fig


def plot_fst_comparison(
    fst_values: dict[str, float],
    *,
    output_path: str | Path | None = None,
    title: str = "Fst Comparison",
    figsize: tuple[float, float] = (8, 6),
) -> plt.Figure:
    """Plot comparison of Fst values.
    
    Args:
        fst_values: Dictionary mapping comparison names to Fst values
        output_path: Optional path to save figure
        title: Plot title
        figsize: Figure size
    
    Returns:
        matplotlib Figure object
    """
    fig, ax = plt.subplots(figsize=figsize)
    
    comparisons = list(fst_values.keys())
    values = list(fst_values.values())
    
    bars = ax.bar(range(len(comparisons)), values, color="darkgreen", alpha=0.7)
    
    # Add value labels
    for bar, val in zip(bars, values):
        height = bar.get_height()
        ax.text(
            bar.get_x() + bar.get_width() / 2.0,
            height,
            f"{val:.3f}",
            ha="center",
            va="bottom",
            fontsize=10,
        )
    
    # Add interpretation lines
    ax.axhline(y=0.05, color="orange", linestyle="--", linewidth=1, alpha=0.5, label="Low")
    ax.axhline(y=0.15, color="red", linestyle="--", linewidth=1, alpha=0.5, label="Moderate")
    ax.axhline(y=0.25, color="darkred", linestyle="--", linewidth=1, alpha=0.5, label="High")
    
    ax.set_xlabel("Comparison", fontsize=12)
    ax.set_ylabel("Fst", fontsize=12)
    ax.set_title(title, fontsize=14, fontweight="bold")
    ax.set_xticks(range(len(comparisons)))
    ax.set_xticklabels(comparisons, rotation=45, ha="right")
    ax.set_ylim(0, max(values) * 1.2 if values else 1.0)
    ax.grid(axis="y", alpha=0.3)
    ax.legend(fontsize=9)
    
    plt.tight_layout()
    
    if output_path:
        fig.savefig(output_path, dpi=300, bbox_inches="tight")
    
    return fig


def plot_pca_results(
    pca_result: dict[str, Any],
    *,
    output_path: str | Path | None = None,
    n_components: int = 3,
    title: str = "Principal Component Analysis",
    figsize: tuple[float, float] = (15, 5),
) -> plt.Figure:
    """Plot PCA results.
    
    Args:
        pca_result: Dictionary with PCA results (from compute_pca)
        output_path: Optional path to save figure
        n_components: Number of components to plot
        title: Plot title
        figsize: Figure size
    
    Returns:
        matplotlib Figure object
    """
    from ..core.errors import ValidationError
    if pca_result.get("status") != "success":
        raise ValidationError("PCA result status is not 'success'")
    
    pcs = np.array(pca_result["pcs"])
    explained_variance_ratio = pca_result.get("explained_variance_ratio", [])
    
    fig, axes = plt.subplots(1, 3, figsize=figsize)
    
    # Plot 1: PC1 vs PC2
    if pcs.shape[1] >= 2:
        ax = axes[0]
        ax.scatter(pcs[:, 0], pcs[:, 1], alpha=0.6, s=20)
        ax.set_xlabel(f"PC1 ({explained_variance_ratio[0]*100:.2f}% variance)" if explained_variance_ratio else "PC1")
        ax.set_ylabel(f"PC2 ({explained_variance_ratio[1]*100:.2f}% variance)" if explained_variance_ratio else "PC2")
        ax.set_title("PC1 vs PC2", fontweight="bold")
        ax.grid(alpha=0.3)
    
    # Plot 2: Explained variance
    if explained_variance_ratio:
        ax = axes[1]
        n_plot = min(n_components, len(explained_variance_ratio))
        ax.bar(range(1, n_plot + 1), explained_variance_ratio[:n_plot], color="steelblue", alpha=0.7)
        ax.set_xlabel("Principal Component", fontsize=11)
        ax.set_ylabel("Explained Variance Ratio", fontsize=11)
        ax.set_title("Explained Variance", fontweight="bold")
        ax.set_xticks(range(1, n_plot + 1))
        ax.grid(axis="y", alpha=0.3)
    
    # Plot 3: Cumulative explained variance
    if explained_variance_ratio:
        ax = axes[2]
        cumulative = np.cumsum(explained_variance_ratio[:n_components])
        ax.plot(range(1, len(cumulative) + 1), cumulative, marker="o", linewidth=2, markersize=6)
        ax.axhline(y=0.8, color="red", linestyle="--", alpha=0.5, label="80% variance")
        ax.set_xlabel("Number of Components", fontsize=11)
        ax.set_ylabel("Cumulative Explained Variance", fontsize=11)
        ax.set_title("Cumulative Variance", fontweight="bold")
        ax.grid(alpha=0.3)
        ax.legend()
    
    fig.suptitle(title, fontsize=14, fontweight="bold", y=1.02)
    plt.tight_layout()
    
    if output_path:
        fig.savefig(output_path, dpi=300, bbox_inches="tight")
    
    return fig


def plot_kinship_matrix(
    kinship_result: dict[str, Any],
    *,
    output_path: str | Path | None = None,
    max_samples: int = 100,
    title: str = "Kinship Matrix",
    figsize: tuple[float, float] = (10, 10),
) -> plt.Figure:
    """Plot kinship matrix as heatmap.
    
    Args:
        kinship_result: Dictionary with kinship results (from compute_kinship_matrix)
        output_path: Optional path to save figure
        max_samples: Maximum number of samples to plot (for large matrices)
        title: Plot title
        figsize: Figure size
    
    Returns:
        matplotlib Figure object
    """
    from ..core.errors import ValidationError
    if kinship_result.get("status") != "success":
        raise ValidationError("Kinship result status is not 'success'")
    
    kinship_matrix = np.array(kinship_result["kinship_matrix"])
    
    # Limit samples for visualization
    if kinship_matrix.shape[0] > max_samples:
        kinship_matrix = kinship_matrix[:max_samples, :max_samples]
    
    fig, ax = plt.subplots(figsize=figsize)
    
    im = ax.imshow(kinship_matrix, cmap="viridis", aspect="auto")
    ax.set_xlabel("Sample", fontsize=12)
    ax.set_ylabel("Sample", fontsize=12)
    ax.set_title(title, fontsize=14, fontweight="bold")
    
    plt.colorbar(im, ax=ax, label="Kinship Coefficient")
    
    plt.tight_layout()
    
    if output_path:
        fig.savefig(output_path, dpi=300, bbox_inches="tight")
    
    return fig


def plot_site_frequency_spectrum(
    sfs: Sequence[int],
    *,
    output_path: str | Path | None = None,
    title: str = "Site Frequency Spectrum",
    figsize: tuple[float, float] = (10, 6),
) -> plt.Figure:
    """Plot site frequency spectrum.
    
    Args:
        sfs: Site frequency spectrum (list of counts per frequency bin)
        output_path: Optional path to save figure
        title: Plot title
        figsize: Figure size
    
    Returns:
        matplotlib Figure object
    """
    fig, ax = plt.subplots(figsize=figsize)
    
    frequencies = list(range(1, len(sfs) + 1))
    counts = list(sfs)
    
    bars = ax.bar(frequencies, counts, color="steelblue", alpha=0.7)
    
    # Add value labels for non-zero bars
    for bar, count in zip(bars, counts):
        if count > 0:
            height = bar.get_height()
            ax.text(
                bar.get_x() + bar.get_width() / 2.0,
                height,
                str(count),
                ha="center",
                va="bottom",
                fontsize=8,
            )
    
    ax.set_xlabel("Minor Allele Frequency", fontsize=12)
    ax.set_ylabel("Number of Sites", fontsize=12)
    ax.set_title(title, fontsize=14, fontweight="bold")
    ax.grid(axis="y", alpha=0.3)
    
    plt.tight_layout()
    
    if output_path:
        fig.savefig(output_path, dpi=300, bbox_inches="tight")
    
    return fig


def plot_neutrality_test_summary(
    neutrality_results: dict[str, dict[str, Any]],
    *,
    output_path: str | Path | None = None,
    title: str = "Neutrality Test Summary",
    figsize: tuple[float, float] = (12, 8),
) -> plt.Figure:
    """Plot summary of neutrality tests across scenarios.
    
    Args:
        neutrality_results: Dictionary mapping scenario names to neutrality test results
        output_path: Optional path to save figure
        title: Plot title
        figsize: Figure size
    
    Returns:
        matplotlib Figure object
    """
    fig, axes = plt.subplots(2, 2, figsize=figsize)
    
    scenarios = list(neutrality_results.keys())
    
    # Plot 1: Tajima's D
    ax = axes[0, 0]
    tajimas_d_vals = [neutrality_results[s].get("tajimas_d", 0) for s in scenarios]
    colors = ["red" if v < 0 else "blue" if v > 0 else "gray" for v in tajimas_d_vals]
    ax.bar(range(len(scenarios)), tajimas_d_vals, color=colors, alpha=0.7)
    ax.axhline(y=0, color="black", linestyle="--", linewidth=1)
    ax.set_xticks(range(len(scenarios)))
    ax.set_xticklabels(scenarios, rotation=45, ha="right", fontsize=9)
    ax.set_ylabel("Tajima's D", fontsize=11)
    ax.set_title("Tajima's D by Scenario", fontweight="bold")
    ax.grid(axis="y", alpha=0.3)
    
    # Plot 2: π/θ ratio
    ax = axes[0, 1]
    pi_theta_ratios = [neutrality_results[s].get("pi_theta_ratio", 0) for s in scenarios]
    ax.bar(range(len(scenarios)), pi_theta_ratios, color="green", alpha=0.7)
    ax.axhline(y=1.0, color="red", linestyle="--", linewidth=1, label="Expected (neutral)")
    ax.set_xticks(range(len(scenarios)))
    ax.set_xticklabels(scenarios, rotation=45, ha="right", fontsize=9)
    ax.set_ylabel("π/θ Ratio", fontsize=11)
    ax.set_title("π/θ Ratio by Scenario", fontweight="bold")
    ax.legend(fontsize=8)
    ax.grid(axis="y", alpha=0.3)
    
    # Plot 3: Nucleotide diversity
    ax = axes[1, 0]
    pi_vals = [neutrality_results[s].get("nucleotide_diversity", 0) for s in scenarios]
    ax.bar(range(len(scenarios)), pi_vals, color="steelblue", alpha=0.7)
    ax.set_xticks(range(len(scenarios)))
    ax.set_xticklabels(scenarios, rotation=45, ha="right", fontsize=9)
    ax.set_ylabel("Nucleotide Diversity (π)", fontsize=11)
    ax.set_title("Nucleotide Diversity by Scenario", fontweight="bold")
    ax.grid(axis="y", alpha=0.3)
    
    # Plot 4: Segregating sites
    ax = axes[1, 1]
    seg_sites = [neutrality_results[s].get("segregating_sites", 0) for s in scenarios]
    ax.bar(range(len(scenarios)), seg_sites, color="purple", alpha=0.7)
    ax.set_xticks(range(len(scenarios)))
    ax.set_xticklabels(scenarios, rotation=45, ha="right", fontsize=9)
    ax.set_ylabel("Segregating Sites (S)", fontsize=11)
    ax.set_title("Segregating Sites by Scenario", fontweight="bold")
    ax.grid(axis="y", alpha=0.3)
    
    fig.suptitle(title, fontsize=14, fontweight="bold", y=0.995)
    plt.tight_layout()
    
    if output_path:
        fig.savefig(output_path, dpi=300, bbox_inches="tight")
    
    return fig


def plot_demographic_comparison(
    demographic_results: dict[str, dict[str, float]],
    *,
    output_path: str | Path | None = None,
    title: str = "Demographic Model Comparison",
    figsize: tuple[float, float] = (12, 6),
) -> plt.Figure:
    """Plot comparison of demographic models.
    
    Args:
        demographic_results: Dictionary with demographic model results
        output_path: Optional path to save figure
        title: Plot title
        figsize: Figure size
    
    Returns:
        matplotlib Figure object
    """
    fig, axes = plt.subplots(1, 2, figsize=figsize)
    
    models = list(demographic_results.keys())
    ne_values = [demographic_results[m].get("estimated_ne", 0) for m in models]
    diversity_values = [demographic_results[m].get("observed_diversity", 0) for m in models]
    
    # Plot 1: Effective population size
    ax = axes[0]
    ax.bar(models, ne_values, color="coral", alpha=0.7)
    ax.set_ylabel("Effective Population Size (Ne)", fontsize=12)
    ax.set_title("Estimated Ne by Model", fontweight="bold")
    ax.grid(axis="y", alpha=0.3)
    for i, (model, ne) in enumerate(zip(models, ne_values)):
        ax.text(i, ne, f"{ne:.1f}", ha="center", va="bottom", fontsize=9)
    
    # Plot 2: Observed diversity
    ax = axes[1]
    ax.bar(models, diversity_values, color="steelblue", alpha=0.7)
    ax.set_ylabel("Observed Diversity (π)", fontsize=12)
    ax.set_title("Observed Diversity by Model", fontweight="bold")
    ax.grid(axis="y", alpha=0.3)
    for i, (model, div) in enumerate(zip(models, diversity_values)):
        ax.text(i, div, f"{div:.4f}", ha="center", va="bottom", fontsize=9)
    
    fig.suptitle(title, fontsize=14, fontweight="bold", y=1.02)
    plt.tight_layout()
    
    if output_path:
        fig.savefig(output_path, dpi=300, bbox_inches="tight")
    
    return fig


def plot_summary_statistics_grid(
    summary_stats: dict[str, dict[str, Any]],
    *,
    output_path: str | Path | None = None,
    title: str = "Summary Statistics Grid",
    figsize: tuple[float, float] = (15, 10),
) -> plt.Figure:
    """Plot grid of summary statistics across scenarios.
    
    Args:
        summary_stats: Dictionary mapping scenario names to summary statistics
        output_path: Optional path to save figure
        title: Plot title
        figsize: Figure size
    
    Returns:
        matplotlib Figure object
    """
    scenarios = list(summary_stats.keys())
    
    fig, axes = plt.subplots(2, 3, figsize=figsize)
    axes = axes.flatten()
    
    # Plot 1: Nucleotide diversity
    ax = axes[0]
    pi_vals = [summary_stats[s].get("nucleotide_diversity", 0) for s in scenarios]
    ax.bar(range(len(scenarios)), pi_vals, color="steelblue", alpha=0.7)
    ax.set_xticks(range(len(scenarios)))
    ax.set_xticklabels(scenarios, rotation=45, ha="right", fontsize=8)
    ax.set_ylabel("π", fontsize=11)
    ax.set_title("Nucleotide Diversity", fontweight="bold")
    ax.grid(axis="y", alpha=0.3)
    
    # Plot 2: Segregating sites
    ax = axes[1]
    seg_vals = [summary_stats[s].get("segregating_sites", 0) for s in scenarios]
    ax.bar(range(len(scenarios)), seg_vals, color="purple", alpha=0.7)
    ax.set_xticks(range(len(scenarios)))
    ax.set_xticklabels(scenarios, rotation=45, ha="right", fontsize=8)
    ax.set_ylabel("S", fontsize=11)
    ax.set_title("Segregating Sites", fontweight="bold")
    ax.grid(axis="y", alpha=0.3)
    
    # Plot 3: Watterson's theta
    ax = axes[2]
    theta_vals = [summary_stats[s].get("wattersons_theta", 0) for s in scenarios]
    ax.bar(range(len(scenarios)), theta_vals, color="green", alpha=0.7)
    ax.set_xticks(range(len(scenarios)))
    ax.set_xticklabels(scenarios, rotation=45, ha="right", fontsize=8)
    ax.set_ylabel("θ_W", fontsize=11)
    ax.set_title("Watterson's Theta", fontweight="bold")
    ax.grid(axis="y", alpha=0.3)
    
    # Plot 4: Tajima's D
    ax = axes[3]
    d_vals = [summary_stats[s].get("tajimas_d", 0) for s in scenarios]
    colors = ["red" if v < 0 else "blue" if v > 0 else "gray" for v in d_vals]
    ax.bar(range(len(scenarios)), d_vals, color=colors, alpha=0.7)
    ax.axhline(y=0, color="black", linestyle="--", linewidth=1)
    ax.set_xticks(range(len(scenarios)))
    ax.set_xticklabels(scenarios, rotation=45, ha="right", fontsize=8)
    ax.set_ylabel("Tajima's D", fontsize=11)
    ax.set_title("Tajima's D", fontweight="bold")
    ax.grid(axis="y", alpha=0.3)
    
    # Plot 5: Sample size
    ax = axes[4]
    sample_sizes = [summary_stats[s].get("sample_size", 0) for s in scenarios]
    ax.bar(range(len(scenarios)), sample_sizes, color="orange", alpha=0.7)
    ax.set_xticks(range(len(scenarios)))
    ax.set_xticklabels(scenarios, rotation=45, ha="right", fontsize=8)
    ax.set_ylabel("Sample Size", fontsize=11)
    ax.set_title("Sample Size", fontweight="bold")
    ax.grid(axis="y", alpha=0.3)
    
    # Plot 6: Sequence length
    ax = axes[5]
    lengths = [summary_stats[s].get("sequence_length", 0) for s in scenarios]
    ax.bar(range(len(scenarios)), lengths, color="teal", alpha=0.7)
    ax.set_xticks(range(len(scenarios)))
    ax.set_xticklabels(scenarios, rotation=45, ha="right", fontsize=8)
    ax.set_ylabel("Sequence Length", fontsize=11)
    ax.set_title("Sequence Length", fontweight="bold")
    ax.grid(axis="y", alpha=0.3)
    
    fig.suptitle(title, fontsize=16, fontweight="bold", y=0.995)
    plt.tight_layout()
    
    if output_path:
        fig.savefig(output_path, dpi=300, bbox_inches="tight")
    
    return fig


def plot_linkage_disequilibrium_decay(
    ld_values: Sequence[float],
    distances: Sequence[float] | None = None,
    *,
    output_path: str | Path | None = None,
    title: str = "Linkage Disequilibrium Decay",
    figsize: tuple[float, float] = (10, 6),
) -> plt.Figure:
    """Plot linkage disequilibrium decay with distance.
    
    Args:
        ld_values: r² values at different distances
        distances: Optional distances (if None, uses indices)
        output_path: Optional path to save figure
        title: Plot title
        figsize: Figure size
    
    Returns:
        matplotlib Figure object
    """
    fig, ax = plt.subplots(figsize=figsize)
    
    if distances is None:
        distances = list(range(len(ld_values)))
    
    ax.plot(distances, ld_values, marker="o", linewidth=2, markersize=6, color="steelblue")
    ax.set_xlabel("Distance (sites)", fontsize=12)
    ax.set_ylabel("r² (Linkage Disequilibrium)", fontsize=12)
    ax.set_title(title, fontsize=14, fontweight="bold")
    ax.grid(alpha=0.3)
    ax.set_ylim(0, max(ld_values) * 1.1 if ld_values else 1.0)
    
    plt.tight_layout()
    
    if output_path:
        fig.savefig(output_path, dpi=300, bbox_inches="tight")
    
    return fig
def plot_allele_frequency_spectrum(
    sfs: Sequence[int],
    *,
    unfolded: bool = False,
    output_path: str | Path | None = None,
    title: str = "Site Frequency Spectrum",
    figsize: tuple[float, float] = (10, 6),
) -> plt.Figure:
    """Plot allele frequency spectrum (SFS).
    
    Args:
        sfs: Site frequency spectrum (counts of sites at each frequency)
        unfolded: Whether SFS is unfolded (requires ancestral state)
        output_path: Optional path to save figure
        title: Plot title
        figsize: Figure size
    
    Returns:
        matplotlib Figure object
    """
    fig, ax = plt.subplots(figsize=figsize)
    
    frequencies = list(range(1, len(sfs) + 1))
    ax.bar(frequencies, sfs, color="steelblue", alpha=0.7)
    ax.set_xlabel("Allele Frequency (minor allele count)", fontsize=12)
    ax.set_ylabel("Number of Sites", fontsize=12)
    ax.set_title(title, fontsize=14, fontweight="bold")
    ax.grid(axis="y", alpha=0.3)
    
    plt.tight_layout()
    
    if output_path:
        fig.savefig(output_path, dpi=300, bbox_inches="tight")
    
    return fig


def plot_pairwise_distance_distribution(
    pairwise_distances: Sequence[float],
    *,
    output_path: str | Path | None = None,
    title: str = "Pairwise Distance Distribution",
    figsize: tuple[float, float] = (10, 6),
    bins: int = 30,
) -> plt.Figure:
    """Plot distribution of pairwise nucleotide differences.
    
    Args:
        pairwise_distances: Pairwise distance values
        output_path: Optional path to save figure
        title: Plot title
        figsize: Figure size
        bins: Number of histogram bins
    
    Returns:
        matplotlib Figure object
    """
    fig, ax = plt.subplots(figsize=figsize)
    ax.hist(pairwise_distances, bins=bins, color="steelblue", alpha=0.7, edgecolor="black")
    ax.set_xlabel("Pairwise Distance", fontsize=12)
    ax.set_ylabel("Frequency", fontsize=12)
    ax.set_title(title, fontsize=14, fontweight="bold")
    ax.grid(axis="y", alpha=0.3)
    plt.tight_layout()
    if output_path:
        fig.savefig(output_path, dpi=300, bbox_inches="tight")
    return fig


def plot_heterozygosity_distribution(
    heterozygosity_values: Sequence[float],
    *,
    output_path: str | Path | None = None,
    title: str = "Heterozygosity Distribution",
    figsize: tuple[float, float] = (10, 6),
    bins: int = 30,
) -> plt.Figure:
    """Plot distribution of observed heterozygosity across sites."""
    fig, ax = plt.subplots(figsize=figsize)
    ax.hist(heterozygosity_values, bins=bins, color="steelblue", alpha=0.7, edgecolor="black")
    ax.set_xlabel("Heterozygosity", fontsize=12)
    ax.set_ylabel("Frequency", fontsize=12)
    ax.set_title(title, fontsize=14, fontweight="bold")
    ax.grid(axis="y", alpha=0.3)
    plt.tight_layout()
    if output_path:
        fig.savefig(output_path, dpi=300, bbox_inches="tight")
    return fig


def plot_statistic_distribution(
    statistic_values: dict[str, Sequence[float]],
    *,
    plot_type: str = "histogram",
    output_path: str | Path | None = None,
    title: str = "Statistic Distribution",
    figsize: tuple[float, float] = (10, 6),
) -> plt.Figure:
    """Plot distribution of statistics across scenarios."""
    fig, ax = plt.subplots(figsize=figsize)
    if plot_type == "box":
        data = [statistic_values[k] for k in statistic_values.keys()]
        labels = list(statistic_values.keys())
        ax.boxplot(data, labels=labels)
        ax.set_xticklabels(labels, rotation=45, ha="right")
    elif plot_type == "violin" and SEABORN_AVAILABLE:
        import pandas as pd
        df_data = []
        for scenario, values in statistic_values.items():
            for val in values:
                df_data.append({"scenario": scenario, "value": val})
        df = pd.DataFrame(df_data)
        sns.violinplot(data=df, x="scenario", y="value", ax=ax)
        ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha="right")
    else:
        for scenario, values in statistic_values.items():
            ax.hist(values, alpha=0.5, label=scenario, bins=20)
        ax.legend()
    ax.set_ylabel("Frequency", fontsize=12)
    ax.set_title(title, fontsize=14, fontweight="bold")
    ax.grid(axis="y", alpha=0.3)
    plt.tight_layout()
    if output_path:
        fig.savefig(output_path, dpi=300, bbox_inches="tight")
    return fig


def plot_pi_vs_theta(
    pi_values: Sequence[float],
    theta_values: Sequence[float],
    *,
    output_path: str | Path | None = None,
    title: str = "π vs θ_W Comparison",
    figsize: tuple[float, float] = (10, 6),
) -> plt.Figure:
    """Plot scatter of π vs θ_W with expected 1:1 line."""
    fig, ax = plt.subplots(figsize=figsize)
    ax.scatter(theta_values, pi_values, alpha=0.6, s=50, color="steelblue")
    max_val = max(max(pi_values), max(theta_values)) if pi_values and theta_values else 1.0
    ax.plot([0, max_val], [0, max_val], "r--", linewidth=2, label="Expected 1:1 (neutrality)")
    if len(pi_values) > 1 and len(theta_values) > 1:
        try:
            from scipy import stats
            slope, intercept, r_value, p_value, std_err = stats.linregress(theta_values, pi_values)
            x_line = np.array(theta_values)
            y_line = slope * x_line + intercept
            ax.plot(x_line, y_line, "g-", linewidth=2, label=f"Regression (R²={r_value**2:.3f})")
        except ImportError:
            pass
    ax.set_xlabel("θ_W (Watterson's Theta)", fontsize=12)
    ax.set_ylabel("π (Nucleotide Diversity)", fontsize=12)
    ax.set_title(title, fontsize=14, fontweight="bold")
    ax.legend()
    ax.grid(alpha=0.3)
    plt.tight_layout()
    if output_path:
        fig.savefig(output_path, dpi=300, bbox_inches="tight")
    return fig


def plot_statistic_correlation_matrix(
    statistics: dict[str, Sequence[float]],
    *,
    output_path: str | Path | None = None,
    title: str = "Statistic Correlation Matrix",
    figsize: tuple[float, float] = (10, 8),
) -> plt.Figure:
    """Plot correlation heatmap between statistics."""
    fig, ax = plt.subplots(figsize=figsize)
    stat_names = list(statistics.keys())
    n_stats = len(stat_names)
    corr_matrix = np.zeros((n_stats, n_stats))
    for i, stat1 in enumerate(stat_names):
        for j, stat2 in enumerate(stat_names):
            if i == j:
                corr_matrix[i, j] = 1.0
            else:
                vals1 = np.array(statistics[stat1])
                vals2 = np.array(statistics[stat2])
                if len(vals1) == len(vals2) and len(vals1) > 1:
                    corr = np.corrcoef(vals1, vals2)[0, 1]
                    corr_matrix[i, j] = corr if not np.isnan(corr) else 0.0
    im = ax.imshow(corr_matrix, cmap="coolwarm", vmin=-1, vmax=1, aspect="auto")
    cbar = plt.colorbar(im, ax=ax)
    cbar.set_label("Correlation", fontsize=12)
    ax.set_xticks(range(n_stats))
    ax.set_yticks(range(n_stats))
    ax.set_xticklabels(stat_names, rotation=45, ha="right")
    ax.set_yticklabels(stat_names)
    for i in range(n_stats):
        for j in range(n_stats):
            ax.text(j, i, f"{corr_matrix[i, j]:.2f}", ha="center", va="center", color="black", fontsize=8)
    ax.set_title(title, fontsize=14, fontweight="bold")
    plt.tight_layout()
    if output_path:
        fig.savefig(output_path, dpi=300, bbox_inches="tight")
    return fig


def plot_regression_analysis(
    x_values: Sequence[float],
    y_values: Sequence[float],
    *,
    x_label: str = "X",
    y_label: str = "Y",
    output_path: str | Path | None = None,
    title: str = "Regression Analysis",
    figsize: tuple[float, float] = (10, 6),
) -> plt.Figure:
    """Plot regression analysis with confidence bands."""
    fig, ax = plt.subplots(figsize=figsize)
    ax.scatter(x_values, y_values, alpha=0.6, s=50, color="steelblue")
    x_array = np.array(x_values)
    y_array = np.array(y_values)
    if len(x_array) > 1:
        try:
            from scipy import stats
            slope, intercept, r_value, p_value, std_err = stats.linregress(x_array, y_array)
            x_line = np.linspace(min(x_array), max(x_array), 100)
            y_line = slope * x_line + intercept
            ax.plot(x_line, y_line, "r-", linewidth=2, label=f"Regression (R²={r_value**2:.3f})")
            y_pred = slope * x_array + intercept
            residuals = y_array - y_pred
            mse = np.mean(residuals**2)
            std_dev = np.sqrt(mse)
            ax.fill_between(x_line, y_line - 2*std_dev, y_line + 2*std_dev, alpha=0.2, color="red")
        except ImportError:
            pass
    ax.set_xlabel(x_label, fontsize=12)
    ax.set_ylabel(y_label, fontsize=12)
    ax.set_title(title, fontsize=14, fontweight="bold")
    ax.legend()
    ax.grid(alpha=0.3)
    plt.tight_layout()
    if output_path:
        fig.savefig(output_path, dpi=300, bbox_inches="tight")
    return fig


def plot_statistic_relationships(
    statistics: dict[str, Sequence[float]],
    *,
    output_path: str | Path | None = None,
    title: str = "Statistic Relationships",
    figsize: tuple[float, float] = (12, 10),
) -> plt.Figure:
    """Plot multi-panel relationship plots between statistics."""
    stat_names = list(statistics.keys())
    n_stats = len(stat_names)
    if n_stats < 2:
        fig, ax = plt.subplots(figsize=figsize)
        ax.text(0.5, 0.5, "Need at least 2 statistics", ha="center", va="center")
        return fig
    fig, axes = plt.subplots(n_stats - 1, n_stats - 1, figsize=figsize)
    axes = axes.flatten() if n_stats > 2 else [axes]
    plot_idx = 0
    for i in range(n_stats):
        for j in range(i + 1, n_stats):
            if plot_idx < len(axes):
                ax = axes[plot_idx]
                x_vals = statistics[stat_names[i]]
                y_vals = statistics[stat_names[j]]
                if len(x_vals) == len(y_vals):
                    ax.scatter(x_vals, y_vals, alpha=0.6, s=30)
                    ax.set_xlabel(stat_names[i], fontsize=8)
                    ax.set_ylabel(stat_names[j], fontsize=8)
                    ax.grid(alpha=0.3)
                plot_idx += 1
    for idx in range(plot_idx, len(axes)):
        axes[idx].set_visible(False)
    fig.suptitle(title, fontsize=14, fontweight="bold")
    plt.tight_layout()
    if output_path:
        fig.savefig(output_path, dpi=300, bbox_inches="tight")
    return fig


def plot_fst_matrix(
    fst_matrix: Sequence[Sequence[float]] | dict[str, dict[str, float]],
    population_names: Sequence[str] | None = None,
    *,
    output_path: str | Path | None = None,
    title: str = "Fst Matrix",
    figsize: tuple[float, float] = (10, 8),
) -> plt.Figure:
    """Plot heatmap of pairwise Fst values."""
    fig, ax = plt.subplots(figsize=figsize)
    if isinstance(fst_matrix, dict):
        pop_names = list(fst_matrix.keys())
        matrix = np.zeros((len(pop_names), len(pop_names)))
        for i, pop1 in enumerate(pop_names):
            for j, pop2 in enumerate(pop_names):
                matrix[i, j] = fst_matrix[pop1].get(pop2, 0.0)
        fst_matrix = matrix.tolist()
        population_names = pop_names
    matrix_array = np.array(fst_matrix)
    if population_names is None:
        population_names = [f"Pop{i+1}" for i in range(len(matrix_array))]
    im = ax.imshow(matrix_array, cmap="YlOrRd", aspect="auto")
    cbar = plt.colorbar(im, ax=ax)
    cbar.set_label("Fst", fontsize=12)
    ax.set_xticks(range(len(population_names)))
    ax.set_yticks(range(len(population_names)))
    ax.set_xticklabels(population_names, rotation=45, ha="right")
    ax.set_yticklabels(population_names)
    for i in range(len(population_names)):
        for j in range(len(population_names)):
            ax.text(j, i, f"{matrix_array[i, j]:.3f}", ha="center", va="center", color="black", fontsize=8)
    ax.set_title(title, fontsize=14, fontweight="bold")
    plt.tight_layout()
    if output_path:
        fig.savefig(output_path, dpi=300, bbox_inches="tight")
    return fig


def plot_admixture_plot(
    ancestry_proportions: Sequence[Sequence[float]],
    population_names: Sequence[str] | None = None,
    *,
    output_path: str | Path | None = None,
    title: str = "Admixture Proportions",
    figsize: tuple[float, float] = (12, 6),
) -> plt.Figure:
    """Plot ancestry proportions as stacked bar plot."""
    fig, ax = plt.subplots(figsize=figsize)
    n_individuals = len(ancestry_proportions)
    n_components = len(ancestry_proportions[0]) if ancestry_proportions else 0
    if population_names is None:
        population_names = [f"K{i+1}" for i in range(n_components)]
    bottom = np.zeros(n_individuals)
    colors = plt.cm.Set3(np.linspace(0, 1, n_components))
    for i in range(n_components):
        proportions = [ancestry_proportions[j][i] for j in range(n_individuals)]
        ax.bar(range(n_individuals), proportions, bottom=bottom, label=population_names[i], color=colors[i])
        bottom += np.array(proportions)
    ax.set_xlabel("Individual", fontsize=12)
    ax.set_ylabel("Ancestry Proportion", fontsize=12)
    ax.set_title(title, fontsize=14, fontweight="bold")
    ax.legend()
    ax.set_ylim(0, 1)
    plt.tight_layout()
    if output_path:
        fig.savefig(output_path, dpi=300, bbox_inches="tight")
    return fig


def plot_population_structure_tree(
    fst_matrix: Sequence[Sequence[float]] | dict[str, dict[str, float]],
    population_names: Sequence[str] | None = None,
    *,
    output_path: str | Path | None = None,
    title: str = "Population Structure Tree",
    figsize: tuple[float, float] = (10, 8),
) -> plt.Figure:
    """Plot neighbor-joining tree based on Fst or genetic distances."""
    try:
        from scipy.cluster.hierarchy import dendrogram, linkage
        from scipy.spatial.distance import squareform
    except ImportError:
        fig, ax = plt.subplots(figsize=figsize)
        ax.text(0.5, 0.5, "scipy required for tree plotting", ha="center", va="center")
        if output_path:
            fig.savefig(output_path, dpi=300, bbox_inches="tight")
        return fig
    fig, ax = plt.subplots(figsize=figsize)
    if isinstance(fst_matrix, dict):
        pop_names = list(fst_matrix.keys())
        matrix = np.zeros((len(pop_names), len(pop_names)))
        for i, pop1 in enumerate(pop_names):
            for j, pop2 in enumerate(pop_names):
                matrix[i, j] = fst_matrix[pop1].get(pop2, 0.0)
        fst_matrix = matrix.tolist()
        population_names = pop_names
    matrix_array = np.array(fst_matrix)
    if population_names is None:
        population_names = [f"Pop{i+1}" for i in range(len(matrix_array))]
    dist_matrix = squareform(matrix_array)
    Z = linkage(dist_matrix, method="average")
    dendrogram(Z, labels=population_names, ax=ax, leaf_rotation=45, leaf_font_size=10)
    ax.set_title(title, fontsize=14, fontweight="bold")
    ax.set_ylabel("Distance", fontsize=12)
    plt.tight_layout()
    if output_path:
        fig.savefig(output_path, dpi=300, bbox_inches="tight")
    return fig


def plot_isolation_by_distance(
    fst_values: Sequence[float],
    geographic_distances: Sequence[float],
    *,
    output_path: str | Path | None = None,
    title: str = "Isolation by Distance",
    figsize: tuple[float, float] = (10, 6),
) -> plt.Figure:
    """Plot Fst vs geographic distance with Mantel test visualization."""
    fig, ax = plt.subplots(figsize=figsize)
    ax.scatter(geographic_distances, fst_values, alpha=0.6, s=50, color="steelblue")
    if len(fst_values) > 1 and len(geographic_distances) > 1:
        try:
            from scipy import stats
            slope, intercept, r_value, p_value, std_err = stats.linregress(geographic_distances, fst_values)
            x_line = np.linspace(min(geographic_distances), max(geographic_distances), 100)
            y_line = slope * x_line + intercept
            ax.plot(x_line, y_line, "r-", linewidth=2, label=f"Regression (R²={r_value**2:.3f}, p={p_value:.3f})")
        except ImportError:
            pass
    ax.set_xlabel("Geographic Distance", fontsize=12)
    ax.set_ylabel("Fst", fontsize=12)
    ax.set_title(title, fontsize=14, fontweight="bold")
    ax.legend()
    ax.grid(alpha=0.3)
    plt.tight_layout()
    if output_path:
        fig.savefig(output_path, dpi=300, bbox_inches="tight")
    return fig


def plot_three_population_f3(
    f3_values: dict[str, float],
    *,
    output_path: str | Path | None = None,
    title: str = "F3-Statistic (Admixture Test)",
    figsize: tuple[float, float] = (10, 6),
) -> plt.Figure:
    """Plot F3-statistic for admixture testing."""
    fig, ax = plt.subplots(figsize=figsize)
    labels = list(f3_values.keys())
    values = list(f3_values.values())
    colors = ["red" if v < 0 else "steelblue" for v in values]
    ax.barh(range(len(labels)), values, color=colors, alpha=0.7)
    ax.axvline(x=0, color="black", linestyle="--", linewidth=1)
    ax.set_yticks(range(len(labels)))
    ax.set_yticklabels(labels)
    ax.set_xlabel("F3-Statistic", fontsize=12)
    ax.set_title(title, fontsize=14, fontweight="bold")
    ax.grid(axis="x", alpha=0.3)
    plt.tight_layout()
    if output_path:
        fig.savefig(output_path, dpi=300, bbox_inches="tight")
    return fig


def plot_neutrality_test_suite(
    test_results: dict[str, dict[str, float]],
    *,
    output_path: str | Path | None = None,
    title: str = "Neutrality Test Suite",
    figsize: tuple[float, float] = (14, 10),
) -> plt.Figure:
    """Plot comprehensive panel of all neutrality tests.
    
    Args:
        test_results: Dictionary mapping scenario names to neutrality test results.
            Each result dict should contain test statistics like 'tajimas_d',
            'fu_and_li_d_star', 'fay_wu_h', etc.
        output_path: Optional path to save figure
        title: Plot title
        figsize: Figure size
    
    Returns:
        matplotlib Figure object
    """
    fig, axes = plt.subplots(2, 3, figsize=figsize)
    axes = axes.flatten()
    
    scenario_names = list(test_results.keys())
    
    # Collect all test statistics across scenarios
    test_stats = {}
    for scenario_name, result in test_results.items():
        # Extract Tajima's D
        if "tajimas_d" in result:
            if "tajimas_d" not in test_stats:
                test_stats["tajimas_d"] = []
            test_stats["tajimas_d"].append(result["tajimas_d"])
        # Extract Fu & Li's D*
        if "fu_and_li_d_star" in result:
            if "fu_and_li_d_star" not in test_stats:
                test_stats["fu_and_li_d_star"] = []
            test_stats["fu_and_li_d_star"].append(result["fu_and_li_d_star"])
        # Extract Fu & Li's F*
        if "fu_and_li_f_star" in result:
            if "fu_and_li_f_star" not in test_stats:
                test_stats["fu_and_li_f_star"] = []
            test_stats["fu_and_li_f_star"].append(result["fu_and_li_f_star"])
        # Extract Fay & Wu's H
        if "fay_wu_h" in result:
            if "fay_wu_h" not in test_stats:
                test_stats["fay_wu_h"] = []
            test_stats["fay_wu_h"].append(result["fay_wu_h"])
        # Extract π/θ ratio
        if "pi_theta_ratio" in result:
            if "pi_theta_ratio" not in test_stats:
                test_stats["pi_theta_ratio"] = []
            test_stats["pi_theta_ratio"].append(result["pi_theta_ratio"])
    
    # Plot 1: Tajima's D comparison
    if "tajimas_d" in test_stats and len(test_stats["tajimas_d"]) > 0:
        ax = axes[0]
        values = test_stats["tajimas_d"]
        colors = ["red" if v < -1 else "orange" if v < 0 else "green" if v < 1 else "blue" for v in values]
        ax.bar(range(len(values)), values, color=colors, alpha=0.7)
        ax.axhline(y=0, color="black", linestyle="-", linewidth=1)
        ax.axhline(y=-2, color="red", linestyle="--", linewidth=1, alpha=0.5, label="Strong selection")
        ax.axhline(y=2, color="blue", linestyle="--", linewidth=1, alpha=0.5, label="Balancing selection")
        ax.set_xticks(range(len(scenario_names)))
        ax.set_xticklabels([s.replace("_", " ").title() for s in scenario_names], rotation=45, ha="right", fontsize=8)
        ax.set_ylabel("Tajima's D", fontsize=10)
        ax.set_title("Tajima's D", fontsize=11, fontweight="bold")
        ax.legend(fontsize=7)
        ax.grid(axis="y", alpha=0.3)
    
    # Plot 2: π/θ ratio
    if "pi_theta_ratio" in test_stats and len(test_stats["pi_theta_ratio"]) > 0:
        ax = axes[1]
        values = test_stats["pi_theta_ratio"]
        ax.bar(range(len(values)), values, color="steelblue", alpha=0.7)
        ax.axhline(y=1.0, color="red", linestyle="--", linewidth=2, label="Neutral (π=θ)")
        ax.set_xticks(range(len(scenario_names)))
        ax.set_xticklabels([s.replace("_", " ").title() for s in scenario_names], rotation=45, ha="right", fontsize=8)
        ax.set_ylabel("π/θ Ratio", fontsize=10)
        ax.set_title("π vs θ_W Ratio", fontsize=11, fontweight="bold")
        ax.legend(fontsize=7)
        ax.grid(axis="y", alpha=0.3)
    
    # Plot 3: Fu & Li's D*
    if "fu_and_li_d_star" in test_stats and len(test_stats["fu_and_li_d_star"]) > 0:
        ax = axes[2]
        values = test_stats["fu_and_li_d_star"]
        colors = ["red" if v < -1 else "orange" if v < 0 else "green" if v < 1 else "blue" for v in values]
        ax.bar(range(len(values)), values, color=colors, alpha=0.7)
        ax.axhline(y=0, color="black", linestyle="-", linewidth=1)
        ax.set_xticks(range(len(values)))
        ax.set_xticklabels([s.replace("_", " ").title() for s in scenario_names[:len(values)]], rotation=45, ha="right", fontsize=8)
        ax.set_ylabel("Fu & Li D*", fontsize=10)
        ax.set_title("Fu & Li's D*", fontsize=11, fontweight="bold")
        ax.grid(axis="y", alpha=0.3)
    
    # Plot 4: Fu & Li's F*
    if "fu_and_li_f_star" in test_stats and len(test_stats["fu_and_li_f_star"]) > 0:
        ax = axes[3]
        values = test_stats["fu_and_li_f_star"]
        colors = ["red" if v < -1 else "orange" if v < 0 else "green" if v < 1 else "blue" for v in values]
        ax.bar(range(len(values)), values, color=colors, alpha=0.7)
        ax.axhline(y=0, color="black", linestyle="-", linewidth=1)
        ax.set_xticks(range(len(values)))
        ax.set_xticklabels([s.replace("_", " ").title() for s in scenario_names[:len(values)]], rotation=45, ha="right", fontsize=8)
        ax.set_ylabel("Fu & Li F*", fontsize=10)
        ax.set_title("Fu & Li's F*", fontsize=11, fontweight="bold")
        ax.grid(axis="y", alpha=0.3)
    
    # Plot 5: Fay & Wu's H
    if "fay_wu_h" in test_stats and len(test_stats["fay_wu_h"]) > 0:
        ax = axes[4]
        values = test_stats["fay_wu_h"]
        colors = ["red" if v < 0 else "steelblue" for v in values]
        ax.bar(range(len(values)), values, color=colors, alpha=0.7)
        ax.axhline(y=0, color="black", linestyle="-", linewidth=1)
        ax.set_xticks(range(len(values)))
        ax.set_xticklabels([s.replace("_", " ").title() for s in scenario_names[:len(values)]], rotation=45, ha="right", fontsize=8)
        ax.set_ylabel("Fay & Wu H", fontsize=10)
        ax.set_title("Fay & Wu's H", fontsize=11, fontweight="bold")
        ax.grid(axis="y", alpha=0.3)
    
    # Plot 6: Summary statistics comparison
    if "tajimas_d" in test_stats and len(test_stats["tajimas_d"]) > 0:
        ax = axes[5]
        # Create a summary plot showing all test statistics together
        x = np.arange(len(scenario_names))
        width = 0.2
        
        if "tajimas_d" in test_stats:
            ax.bar(x - 2*width, test_stats["tajimas_d"], width, label="Tajima's D", alpha=0.7)
        if "fu_and_li_d_star" in test_stats and len(test_stats["fu_and_li_d_star"]) == len(x):
            ax.bar(x - width, test_stats["fu_and_li_d_star"], width, label="Fu & Li D*", alpha=0.7)
        if "fu_and_li_f_star" in test_stats and len(test_stats["fu_and_li_f_star"]) == len(x):
            ax.bar(x, test_stats["fu_and_li_f_star"], width, label="Fu & Li F*", alpha=0.7)
        if "fay_wu_h" in test_stats and len(test_stats["fay_wu_h"]) == len(x):
            ax.bar(x + width, test_stats["fay_wu_h"], width, label="Fay & Wu H", alpha=0.7)
        
        ax.axhline(y=0, color="black", linestyle="-", linewidth=1)
        ax.set_xticks(x)
        ax.set_xticklabels([s.replace("_", " ").title() for s in scenario_names], rotation=45, ha="right", fontsize=7)
        ax.set_ylabel("Statistic Value", fontsize=10)
        ax.set_title("All Neutrality Tests Comparison", fontsize=11, fontweight="bold")
        ax.legend(fontsize=7, loc="best")
        ax.grid(axis="y", alpha=0.3)
    
    # Hide unused subplots
    for idx in range(len([s for s in [test_stats.get("tajimas_d"), test_stats.get("pi_theta_ratio"), 
                                      test_stats.get("fu_and_li_d_star"), test_stats.get("fu_and_li_f_star"),
                                      test_stats.get("fay_wu_h")] if s]), len(axes)):
        if idx < len(axes):
            axes[idx].set_visible(False)
    
    fig.suptitle(title, fontsize=14, fontweight="bold")
    plt.tight_layout()
    
    if output_path:
        fig.savefig(output_path, dpi=300, bbox_inches="tight")
    
    return fig


def plot_hardy_weinberg_test(
    hwe_results: dict[str, float] | Sequence[dict[str, float]],
    *,
    output_path: str | Path | None = None,
    title: str = "Hardy-Weinberg Equilibrium Test",
    figsize: tuple[float, float] = (10, 6),
) -> plt.Figure:
    """Plot Hardy-Weinberg test results."""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=figsize)
    if isinstance(hwe_results, dict):
        hwe_results = [hwe_results]
    chi_square_values = [r.get("chi_square", 0.0) for r in hwe_results]
    ax1.bar(range(len(chi_square_values)), chi_square_values, color="steelblue", alpha=0.7)
    ax1.set_xlabel("Site", fontsize=12)
    ax1.set_ylabel("Chi-Square", fontsize=12)
    ax1.set_title("Chi-Square Test Statistics", fontweight="bold")
    ax1.grid(axis="y", alpha=0.3)
    p_values = [r.get("p_value", 1.0) for r in hwe_results]
    ax2.bar(range(len(p_values)), p_values, color="coral", alpha=0.7)
    ax2.axhline(y=0.05, color="red", linestyle="--", linewidth=2, label="α=0.05")
    ax2.set_xlabel("Site", fontsize=12)
    ax2.set_ylabel("P-Value", fontsize=12)
    ax2.set_title("P-Values", fontweight="bold")
    ax2.set_ylim(0, 1.0)
    ax2.legend()
    ax2.grid(axis="y", alpha=0.3)
    fig.suptitle(title, fontsize=14, fontweight="bold")
    plt.tight_layout()
    if output_path:
        fig.savefig(output_path, dpi=300, bbox_inches="tight")
    return fig


def plot_bootstrap_distribution(
    bootstrap_values: Sequence[float],
    observed_value: float | None = None,
    *,
    confidence_level: float = 0.95,
    output_path: str | Path | None = None,
    title: str = "Bootstrap Distribution",
    figsize: tuple[float, float] = (10, 6),
) -> plt.Figure:
    """Plot bootstrap distribution with confidence intervals."""
    fig, ax = plt.subplots(figsize=figsize)
    ax.hist(bootstrap_values, bins=30, color="steelblue", alpha=0.7, edgecolor="black")
    alpha = 1.0 - confidence_level
    ci_lower = np.percentile(bootstrap_values, (alpha / 2.0) * 100.0)
    ci_upper = np.percentile(bootstrap_values, (1.0 - alpha / 2.0) * 100.0)
    ax.axvline(x=ci_lower, color="red", linestyle="--", linewidth=2, label=f"CI Lower ({ci_lower:.3f})")
    ax.axvline(x=ci_upper, color="red", linestyle="--", linewidth=2, label=f"CI Upper ({ci_upper:.3f})")
    if observed_value is not None:
        ax.axvline(x=observed_value, color="green", linestyle="-", linewidth=2, label=f"Observed ({observed_value:.3f})")
    ax.set_xlabel("Statistic Value", fontsize=12)
    ax.set_ylabel("Frequency", fontsize=12)
    ax.set_title(title, fontsize=14, fontweight="bold")
    ax.legend()
    ax.grid(axis="y", alpha=0.3)
    plt.tight_layout()
    if output_path:
        fig.savefig(output_path, dpi=300, bbox_inches="tight")
    return fig


def plot_permutation_test(
    permuted_values: Sequence[float],
    observed_value: float,
    *,
    p_value: float | None = None,
    output_path: str | Path | None = None,
    title: str = "Permutation Test",
    figsize: tuple[float, float] = (10, 6),
) -> plt.Figure:
    """Plot permuted vs observed distribution."""
    fig, ax = plt.subplots(figsize=figsize)
    ax.hist(permuted_values, bins=30, color="steelblue", alpha=0.7, edgecolor="black", label="Permuted Distribution")
    ax.axvline(x=observed_value, color="red", linestyle="-", linewidth=2, label=f"Observed ({observed_value:.3f})")
    if p_value is not None:
        ax.text(0.7, 0.9, f"P-value: {p_value:.4f}", transform=ax.transAxes, fontsize=12, 
                bbox=dict(boxstyle="round", facecolor="wheat", alpha=0.5))
    ax.set_xlabel("Statistic Value", fontsize=12)
    ax.set_ylabel("Frequency", fontsize=12)
    ax.set_title(title, fontsize=14, fontweight="bold")
    ax.legend()
    ax.grid(axis="y", alpha=0.3)
    plt.tight_layout()
    if output_path:
        fig.savefig(output_path, dpi=300, bbox_inches="tight")
    return fig


def plot_outlier_detection(
    statistic_values: Sequence[float],
    outlier_indices: Sequence[int] | None = None,
    *,
    positions: Sequence[float] | None = None,
    fdr_threshold: float | None = None,
    output_path: str | Path | None = None,
    title: str = "Outlier Detection",
    figsize: tuple[float, float] = (12, 6),
) -> plt.Figure:
    """Plot Manhattan plot of statistics with outliers highlighted."""
    fig, ax = plt.subplots(figsize=figsize)
    if positions is None:
        positions = list(range(len(statistic_values)))
    ax.scatter(positions, statistic_values, alpha=0.6, s=20, color="steelblue", label="All Loci")
    if outlier_indices:
        outlier_positions = [positions[i] for i in outlier_indices]
        outlier_values = [statistic_values[i] for i in outlier_indices]
        ax.scatter(outlier_positions, outlier_values, alpha=0.8, s=50, color="red", label="Outliers")
    if fdr_threshold is not None:
        ax.axhline(y=fdr_threshold, color="orange", linestyle="--", linewidth=2, label=f"FDR Threshold ({fdr_threshold:.3f})")
    ax.set_xlabel("Position", fontsize=12)
    ax.set_ylabel("Statistic Value", fontsize=12)
    ax.set_title(title, fontsize=14, fontweight="bold")
    ax.legend()
    ax.grid(alpha=0.3)
    plt.tight_layout()
    if output_path:
        fig.savefig(output_path, dpi=300, bbox_inches="tight")
    return fig
