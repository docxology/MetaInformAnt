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
    if pca_result.get("status") != "success":
        raise ValueError("PCA result status is not 'success'")
    
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
    if kinship_result.get("status") != "success":
        raise ValueError("Kinship result status is not 'success'")
    
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
