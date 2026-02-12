"""Mathematical biology and population genetics visualization functions.

This module provides comprehensive visualization capabilities for mathematical biology,
including population genetics models, evolutionary dynamics, selection analysis,
and theoretical biology simulations.
"""

from __future__ import annotations

import warnings
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple, Union

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.axes import Axes
from matplotlib.patches import Circle, FancyBboxPatch, Rectangle

from metainformant.core.data import validation
from metainformant.core.io import paths
from metainformant.core.utils import logging

logger = logging.get_logger(__name__)

try:
    import seaborn as sns

    HAS_SEABORN = True
except ImportError:
    HAS_SEABORN = False
    sns = None

try:
    import plotly.express as px
    import plotly.graph_objects as go

    HAS_PLOTLY = True
except ImportError:
    HAS_PLOTLY = False
    go = None
    px = None


def plot_allele_frequency_spectrum(
    allele_frequencies: np.ndarray,
    *,
    ax: Axes | None = None,
    output_path: str | Path | None = None,
    figsize: Tuple[float, float] = (10, 6),
    **kwargs,
) -> Axes:
    """Plot site frequency spectrum from allele frequencies.

    Args:
        allele_frequencies: Array of derived allele frequencies (0-1)
        ax: Optional matplotlib axes to plot on
        output_path: Optional path to save the figure
        figsize: Figure size as (width, height)
        **kwargs: Additional arguments for customization

    Returns:
        matplotlib Axes object
    """
    validation.validate_type(allele_frequencies, np.ndarray, "allele_frequencies")

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)

    # Create histogram of allele frequencies
    bins = np.linspace(0, 1, 21)  # 20 bins from 0 to 1
    hist, bin_edges = np.histogram(allele_frequencies, bins=bins)

    # Plot frequency spectrum
    ax.bar(bin_edges[:-1], hist, width=np.diff(bin_edges), alpha=0.7, color="skyblue", edgecolor="black")

    ax.set_xlabel("Derived Allele Frequency")
    ax.set_ylabel("Number of Sites")
    ax.set_title("Site Frequency Spectrum")
    ax.grid(True, alpha=0.3)

    if output_path:
        output_path = paths.ensure_directory(Path(output_path).parent)
        plt.savefig(output_path, dpi=300, bbox_inches="tight")
        logger.info(f"Allele frequency spectrum saved to {output_path}")

    return ax


def plot_population_genetics_summary(
    summary_stats: Dict[str, Any],
    *,
    ax: Axes | None = None,
    output_path: str | Path | None = None,
    figsize: Tuple[float, float] = (12, 8),
    **kwargs,
) -> Axes:
    """Plot comprehensive population genetics summary statistics.

    Args:
        summary_stats: Dictionary containing various population genetics statistics
        ax: Optional matplotlib axes to plot on
        output_path: Optional path to save the figure
        figsize: Figure size as (width, height)
        **kwargs: Additional arguments for customization

    Returns:
        matplotlib Axes object
    """
    validation.validate_type(summary_stats, dict, "summary_stats")

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)

    # Extract statistics
    stat_names = []
    stat_values = []
    stat_errors = []

    for key, value in summary_stats.items():
        if isinstance(value, dict) and "mean" in value:
            stat_names.append(key)
            stat_values.append(value["mean"])
            stat_errors.append(value.get("std", 0))
        elif isinstance(value, (int, float)):
            stat_names.append(key)
            stat_values.append(value)
            stat_errors.append(0)

    # Plot as bar chart with error bars
    x_positions = np.arange(len(stat_names))
    bars = ax.bar(x_positions, stat_values, yerr=stat_errors, capsize=5, alpha=0.7, color="lightcoral")

    ax.set_xlabel("Population Genetics Statistics")
    ax.set_ylabel("Value")
    ax.set_title("Population Genetics Summary")
    ax.set_xticks(x_positions)
    ax.set_xticklabels(stat_names, rotation=45, ha="right")
    ax.grid(True, alpha=0.3, axis="y")

    if output_path:
        output_path = paths.ensure_directory(Path(output_path).parent)
        plt.savefig(output_path, dpi=300, bbox_inches="tight")
        logger.info(f"Population genetics summary saved to {output_path}")

    return ax


def plot_evolutionary_trajectory(
    trajectory_data: np.ndarray,
    time_points: np.ndarray | None = None,
    *,
    ax: Axes | None = None,
    output_path: str | Path | None = None,
    figsize: Tuple[float, float] = (10, 6),
    **kwargs,
) -> Axes:
    """Plot evolutionary trajectory over time.

    Args:
        trajectory_data: Time series of trait values or frequencies
        time_points: Optional time points for x-axis
        ax: Optional matplotlib axes to plot on
        output_path: Optional path to save the figure
        figsize: Figure size as (width, height)
        **kwargs: Additional arguments for customization

    Returns:
        matplotlib Axes object
    """
    validation.validate_type(trajectory_data, np.ndarray, "trajectory_data")

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)

    x_data = time_points if time_points is not None else np.arange(len(trajectory_data))

    # Plot trajectory
    ax.plot(x_data, trajectory_data, "b-", linewidth=2, marker="o", markersize=4, alpha=0.8)

    # Add trend line if requested
    if kwargs.get("show_trend", False):
        from scipy.stats import linregress

        slope, intercept, r_value, p_value, std_err = linregress(x_data, trajectory_data)
        trend_line = slope * x_data + intercept
        ax.plot(x_data, trend_line, "r--", linewidth=1, alpha=0.7, label=f"Trend (R²={r_value**2:.3f})")

    ax.set_xlabel("Time/Generation")
    ax.set_ylabel("Trait Value/Frequency")
    ax.set_title("Evolutionary Trajectory")
    ax.legend()
    ax.grid(True, alpha=0.3)

    if output_path:
        output_path = paths.ensure_directory(Path(output_path).parent)
        plt.savefig(output_path, dpi=300, bbox_inches="tight")
        logger.info(f"Evolutionary trajectory saved to {output_path}")

    return ax


def plot_selection_coefficient_distribution(
    selection_coeffs: np.ndarray,
    *,
    ax: Axes | None = None,
    output_path: str | Path | None = None,
    figsize: Tuple[float, float] = (8, 6),
    **kwargs,
) -> Axes:
    """Plot distribution of selection coefficients.

    Args:
        selection_coeffs: Array of selection coefficient values
        ax: Optional matplotlib axes to plot on
        output_path: Optional path to save the figure
        figsize: Figure size as (width, height)
        **kwargs: Additional arguments for customization

    Returns:
        matplotlib Axes object
    """
    validation.validate_type(selection_coeffs, np.ndarray, "selection_coeffs")

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)

    if HAS_SEABORN:
        sns.histplot(selection_coeffs, kde=True, ax=ax, **kwargs)
    else:
        ax.hist(selection_coeffs, bins=30, alpha=0.7, density=True, **kwargs)

    ax.axvline(x=0, color="red", linestyle="--", alpha=0.7, label="Neutral (s=0)")
    ax.set_xlabel("Selection Coefficient (s)")
    ax.set_ylabel("Density")
    ax.set_title("Distribution of Selection Coefficients")
    ax.legend()
    ax.grid(True, alpha=0.3)

    if output_path:
        output_path = paths.ensure_directory(Path(output_path).parent)
        plt.savefig(output_path, dpi=300, bbox_inches="tight")
        logger.info(f"Selection coefficient distribution saved to {output_path}")

    return ax


def plot_fst_distribution(
    fst_values: np.ndarray,
    *,
    ax: Axes | None = None,
    output_path: str | Path | None = None,
    figsize: Tuple[float, float] = (8, 6),
    **kwargs,
) -> Axes:
    """Plot distribution of Fst values between populations.

    Args:
        fst_values: Array of Fst values between population pairs
        ax: Optional matplotlib axes to plot on
        output_path: Optional path to save the figure
        figsize: Figure size as (width, height)
        **kwargs: Additional arguments for customization

    Returns:
        matplotlib Axes object
    """
    validation.validate_type(fst_values, np.ndarray, "fst_values")

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)

    if HAS_SEABORN:
        sns.histplot(fst_values, kde=True, ax=ax, **kwargs)
    else:
        ax.hist(fst_values, bins=20, alpha=0.7, density=True, **kwargs)

    # Add reference lines
    ax.axvline(x=0, color="gray", linestyle="--", alpha=0.7, label="No differentiation")
    ax.axvline(x=0.05, color="orange", linestyle="--", alpha=0.7, label="Moderate differentiation")
    ax.axvline(x=0.15, color="red", linestyle="--", alpha=0.7, label="Strong differentiation")

    ax.set_xlabel("Fst Value")
    ax.set_ylabel("Density")
    ax.set_title("Fst Distribution Between Populations")
    ax.legend()
    ax.grid(True, alpha=0.3)

    if output_path:
        output_path = paths.ensure_directory(Path(output_path).parent)
        plt.savefig(output_path, dpi=300, bbox_inches="tight")
        logger.info(f"Fst distribution saved to {output_path}")

    return ax


def plot_coalescent_tree(
    tree_data: Any,
    *,
    ax: Axes | None = None,
    output_path: str | Path | None = None,
    figsize: Tuple[float, float] = (10, 6),
    **kwargs,
) -> Axes:
    """Plot coalescent tree from population genetics simulation.

    Args:
        tree_data: Coalescent tree data (NetworkX graph or similar)
        ax: Optional matplotlib axes to plot on
        output_path: Optional path to save the figure
        figsize: Figure size as (width, height)
        **kwargs: Additional arguments for customization

    Returns:
        matplotlib Axes object
    """
    try:
        import networkx as nx
    except ImportError:
        raise ImportError("networkx required for coalescent tree visualization")

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)

    # Plot tree
    if hasattr(tree_data, "nodes"):
        # NetworkX tree
        pos = nx.spring_layout(tree_data)
        nx.draw(tree_data, pos, with_labels=True, ax=ax, node_color="lightblue", node_size=300, font_size=8)
    else:
        # Simple tree representation
        ax.text(
            0.5,
            0.5,
            "Coalescent Tree\n(Not implemented for this format)",
            ha="center",
            va="center",
            transform=ax.transAxes,
        )

    ax.set_title("Coalescent Tree")
    ax.axis("off")

    if output_path:
        output_path = paths.ensure_directory(Path(output_path).parent)
        plt.savefig(output_path, dpi=300, bbox_inches="tight")
        logger.info(f"Coalescent tree saved to {output_path}")

    return ax


def plot_genetic_drift_simulation(
    frequency_trajectories: List[np.ndarray],
    *,
    ax: Axes | None = None,
    output_path: str | Path | None = None,
    figsize: Tuple[float, float] = (10, 6),
    **kwargs,
) -> Axes:
    """Plot genetic drift simulation trajectories.

    Args:
        frequency_trajectories: List of allele frequency trajectories over time
        ax: Optional matplotlib axes to plot on
        output_path: Optional path to save the figure
        figsize: Figure size as (width, height)
        **kwargs: Additional arguments for customization

    Returns:
        matplotlib Axes object
    """
    validation.validate_type(frequency_trajectories, list, "frequency_trajectories")

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)

    # Plot each trajectory
    colors = plt.cm.tab10(np.linspace(0, 1, len(frequency_trajectories)))

    for i, trajectory in enumerate(frequency_trajectories):
        generations = np.arange(len(trajectory))
        ax.plot(generations, trajectory, color=colors[i], alpha=0.7, linewidth=1)

    # Add starting and ending frequency lines
    ax.axhline(y=0.5, color="red", linestyle="--", alpha=0.5, label="Initial frequency")
    ax.axhline(y=0, color="gray", linestyle=":", alpha=0.5, label="Lost")
    ax.axhline(y=1, color="gray", linestyle=":", alpha=0.5, label="Fixed")

    ax.set_xlabel("Generation")
    ax.set_ylabel("Allele Frequency")
    ax.set_title("Genetic Drift Trajectories")
    ax.set_ylim(-0.05, 1.05)
    ax.legend()
    ax.grid(True, alpha=0.3)

    if output_path:
        output_path = paths.ensure_directory(Path(output_path).parent)
        plt.savefig(output_path, dpi=300, bbox_inches="tight")
        logger.info(f"Genetic drift simulation saved to {output_path}")

    return ax


def plot_moran_model_evolution(
    population_states: List[np.ndarray],
    *,
    ax: Axes | None = None,
    output_path: str | Path | None = None,
    figsize: Tuple[float, float] = (10, 6),
    **kwargs,
) -> Axes:
    """Plot evolution in the Moran model.

    Args:
        population_states: List of population states over time
        ax: Optional matplotlib axes to plot on
        output_path: Optional path to save the figure
        figsize: Figure size as (width, height)
        **kwargs: Additional arguments for customization

    Returns:
        matplotlib Axes object
    """
    validation.validate_type(population_states, list, "population_states")

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)

    # Convert to frequency of mutant alleles
    mutant_frequencies = [np.mean(state) for state in population_states]
    time_points = np.arange(len(mutant_frequencies))

    ax.plot(time_points, mutant_frequencies, "b-", linewidth=2, marker="o", markersize=3, alpha=0.8)

    ax.axhline(y=0, color="red", linestyle="--", alpha=0.7, label="All wild-type")
    ax.axhline(y=1, color="green", linestyle="--", alpha=0.7, label="All mutant")

    ax.set_xlabel("Time Step")
    ax.set_ylabel("Mutant Allele Frequency")
    ax.set_title("Moran Model Evolution")
    ax.set_ylim(-0.05, 1.05)
    ax.legend()
    ax.grid(True, alpha=0.3)

    if output_path:
        output_path = paths.ensure_directory(Path(output_path).parent)
        plt.savefig(output_path, dpi=300, bbox_inches="tight")
        logger.info(f"Moran model evolution saved to {output_path}")

    return ax


def plot_price_equation_components(
    price_components: Dict[str, np.ndarray],
    *,
    ax: Axes | None = None,
    output_path: str | Path | None = None,
    figsize: Tuple[float, float] = (12, 6),
    **kwargs,
) -> Axes:
    """Plot Price equation components over time.

    Args:
        price_components: Dictionary with Price equation component trajectories
        ax: Optional matplotlib axes to plot on
        output_path: Optional path to save the figure
        figsize: Figure size as (width, height)
        **kwargs: Additional arguments for customization

    Returns:
        matplotlib Axes object
    """
    validation.validate_type(price_components, dict, "price_components")

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)

    time_points = np.arange(len(next(iter(price_components.values()))))
    colors = plt.cm.Set1(np.linspace(0, 1, len(price_components)))

    for i, (component_name, values) in enumerate(price_components.items()):
        ax.plot(time_points, values, color=colors[i], linewidth=2, label=component_name, alpha=0.8)

    ax.set_xlabel("Generation")
    ax.set_ylabel("Component Value")
    ax.set_title("Price Equation Components")
    ax.legend()
    ax.grid(True, alpha=0.3)

    if output_path:
        output_path = paths.ensure_directory(Path(output_path).parent)
        plt.savefig(output_path, dpi=300, bbox_inches="tight")
        logger.info(f"Price equation components saved to {output_path}")

    return ax


def plot_epidemic_model_simulation(
    epidemic_data: Dict[str, np.ndarray],
    *,
    ax: Axes | None = None,
    output_path: str | Path | None = None,
    figsize: Tuple[float, float] = (10, 6),
    **kwargs,
) -> Axes:
    """Plot epidemic model simulation results.

    Args:
        epidemic_data: Dictionary with susceptible, infected, recovered trajectories
        ax: Optional matplotlib axes to plot on
        output_path: Optional path to save the figure
        figsize: Figure size as (width, height)
        **kwargs: Additional arguments for customization

    Returns:
        matplotlib Axes object
    """
    validation.validate_type(epidemic_data, dict, "epidemic_data")

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)

    time_points = np.arange(len(next(iter(epidemic_data.values()))))
    colors = {"susceptible": "blue", "infected": "red", "recovered": "green"}

    for compartment, values in epidemic_data.items():
        color = colors.get(compartment.lower(), "black")
        ax.plot(time_points, values, color=color, linewidth=2, label=compartment.capitalize(), alpha=0.8)

    ax.set_xlabel("Time")
    ax.set_ylabel("Population Size")
    ax.set_title("Epidemic Model Simulation")
    ax.legend()
    ax.grid(True, alpha=0.3)

    if output_path:
        output_path = paths.ensure_directory(Path(output_path).parent)
        plt.savefig(output_path, dpi=300, bbox_inches="tight")
        logger.info(f"Epidemic model simulation saved to {output_path}")

    return ax


def plot_population_structure_pca(
    pca_coords: np.ndarray,
    population_labels: np.ndarray | None = None,
    *,
    ax: Axes | None = None,
    output_path: str | Path | None = None,
    figsize: Tuple[float, float] = (8, 6),
    **kwargs,
) -> Axes:
    """Plot population structure from PCA of genetic data.

    Args:
        pca_coords: PCA coordinates (samples x components)
        population_labels: Optional population labels for coloring
        ax: Optional matplotlib axes to plot on
        output_path: Optional path to save the figure
        figsize: Figure size as (width, height)
        **kwargs: Additional arguments for customization

    Returns:
        matplotlib Axes object
    """
    validation.validate_type(pca_coords, np.ndarray, "pca_coords")

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)

    if population_labels is not None:
        unique_pops = np.unique(population_labels)
        colors = plt.cm.tab10(np.linspace(0, 1, len(unique_pops)))

        for i, pop in enumerate(unique_pops):
            mask = population_labels == pop
            ax.scatter(pca_coords[mask, 0], pca_coords[mask, 1], c=[colors[i]], label=str(pop), alpha=0.7, s=50)
    else:
        scatter = ax.scatter(
            pca_coords[:, 0], pca_coords[:, 1], c=range(len(pca_coords)), cmap="viridis", alpha=0.7, s=50
        )

    ax.set_xlabel("PC1")
    ax.set_ylabel("PC2")
    ax.set_title("Population Structure (PCA)")
    ax.grid(True, alpha=0.3)

    if population_labels is not None:
        ax.legend()
    else:
        # Add colorbar
        cbar = plt.colorbar(scatter, ax=ax)
        cbar.set_label("Sample Index")

    if output_path:
        output_path = paths.ensure_directory(Path(output_path).parent)
        plt.savefig(output_path, dpi=300, bbox_inches="tight")
        logger.info(f"Population structure PCA saved to {output_path}")

    return ax


def plot_linkage_disequilibrium_decay(
    ld_values: np.ndarray,
    distances: np.ndarray,
    *,
    ax: Axes | None = None,
    output_path: str | Path | None = None,
    figsize: Tuple[float, float] = (8, 6),
    **kwargs,
) -> Axes:
    """Plot linkage disequilibrium decay with distance.

    Args:
        ld_values: Linkage disequilibrium values (r²)
        distances: Physical or genetic distances
        ax: Optional matplotlib axes to plot on
        output_path: Optional path to save the figure
        figsize: Figure size as (width, height)
        **kwargs: Additional arguments for customization

    Returns:
        matplotlib Axes object
    """
    validation.validate_type(ld_values, np.ndarray, "ld_values")
    validation.validate_type(distances, np.ndarray, "distances")

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)

    ax.scatter(distances, ld_values, alpha=0.6, s=20, color="blue")

    # Add smoothed trend line
    if len(distances) > 10:
        # Bin the data for smoothing
        bins = np.logspace(np.log10(distances.min()), np.log10(distances.max()), 20)
        bin_centers = []
        bin_means = []

        for i in range(len(bins) - 1):
            mask = (distances >= bins[i]) & (distances < bins[i + 1])
            if np.any(mask):
                bin_centers.append((bins[i] + bins[i + 1]) / 2)
                bin_means.append(np.mean(ld_values[mask]))

        if bin_centers and bin_means:
            ax.plot(bin_centers, bin_means, "r-", linewidth=2, alpha=0.8, label="Trend")

    ax.set_xlabel("Distance")
    ax.set_ylabel("Linkage Disequilibrium (r²)")
    ax.set_title("LD Decay with Distance")
    ax.set_xscale("log")
    ax.legend()
    ax.grid(True, alpha=0.3)

    if output_path:
        output_path = paths.ensure_directory(Path(output_path).parent)
        plt.savefig(output_path, dpi=300, bbox_inches="tight")
        logger.info(f"LD decay plot saved to {output_path}")

    return ax


def create_interactive_population_genetics_dashboard(
    popgen_data: Dict[str, Any], *, output_path: str | Path | None = None, **kwargs
) -> Any:
    """Create an interactive population genetics dashboard using Plotly.

    Args:
        popgen_data: Dictionary containing population genetics data and statistics
        output_path: Optional path to save the HTML file
        **kwargs: Additional arguments for Plotly customization

    Returns:
        Plotly Figure object
    """
    if not HAS_PLOTLY:
        raise ImportError("Plotly required for interactive population genetics dashboard")

    validation.validate_type(popgen_data, dict, "popgen_data")

    # Create subplot figure
    fig = go.Figure()

    # Add allele frequency spectrum
    if "allele_frequencies" in popgen_data:
        freqs = popgen_data["allele_frequencies"]
        fig.add_trace(go.Histogram(x=freqs, nbinsx=20, name="Allele Frequencies"))

    # Add diversity statistics
    if "diversity_stats" in popgen_data:
        stats = popgen_data["diversity_stats"]
        fig.add_trace(go.Bar(x=list(stats.keys()), y=list(stats.values()), name="Diversity Statistics"))

    fig.update_layout(title="Interactive Population Genetics Dashboard", **kwargs)

    if output_path:
        output_path = paths.ensure_directory(Path(output_path).parent)
        html_path = Path(output_path).with_suffix(".html")
        fig.write_html(str(html_path))
        logger.info(f"Interactive population genetics dashboard saved to {html_path}")

    return fig
