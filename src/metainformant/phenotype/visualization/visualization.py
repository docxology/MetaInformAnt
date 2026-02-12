"""Phenotype data visualization functions.

This module provides comprehensive visualization capabilities for phenotypic data,
including trait distributions, life course analysis, morphological measurements,
and behavioral phenotype patterns.
"""

from __future__ import annotations

import warnings
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple, Union

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.axes import Axes
from matplotlib.patches import Circle, Rectangle

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


def plot_trait_distribution(
    trait_values: np.ndarray,
    trait_name: str = "Trait",
    *,
    ax: Axes | None = None,
    output_path: str | Path | None = None,
    figsize: Tuple[float, float] = (8, 6),
    **kwargs,
) -> Axes:
    """Plot distribution of phenotypic trait values.

    Args:
        trait_values: Array of trait measurements
        trait_name: Name of the trait for labeling
        ax: Optional matplotlib axes to plot on
        output_path: Optional path to save the figure
        figsize: Figure size as (width, height)
        **kwargs: Additional arguments for customization

    Returns:
        matplotlib Axes object
    """
    validation.validate_type(trait_values, np.ndarray, "trait_values")

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)

    # Plot histogram with density curve
    if HAS_SEABORN:
        sns.histplot(trait_values, kde=True, ax=ax, **kwargs)
    else:
        ax.hist(trait_values, bins=30, alpha=0.7, density=True, **kwargs)
        # Add simple density estimate
        from scipy.stats import gaussian_kde

        try:
            kde = gaussian_kde(trait_values)
            x_range = np.linspace(trait_values.min(), trait_values.max(), 100)
            ax.plot(x_range, kde(x_range), "r-", linewidth=2)
        except (ValueError, np.linalg.LinAlgError):
            pass  # Skip KDE if data is singular or insufficient

    ax.set_xlabel(f"{trait_name} Value")
    ax.set_ylabel("Density")
    ax.set_title(f"Distribution of {trait_name}")
    ax.grid(True, alpha=0.3)

    if output_path:
        output_path = paths.ensure_directory(Path(output_path).parent)
        plt.savefig(output_path, dpi=300, bbox_inches="tight")
        logger.info(f"Trait distribution plot saved to {output_path}")

    return ax


def plot_trait_correlation_matrix(
    trait_data: pd.DataFrame,
    *,
    ax: Axes | None = None,
    output_path: str | Path | None = None,
    figsize: Tuple[float, float] = (10, 8),
    **kwargs,
) -> Axes:
    """Plot correlation matrix between phenotypic traits.

    Args:
        trait_data: DataFrame with traits as columns
        ax: Optional matplotlib axes to plot on
        output_path: Optional path to save the figure
        figsize: Figure size as (width, height)
        **kwargs: Additional arguments for customization

    Returns:
        matplotlib Axes object
    """
    validation.validate_type(trait_data, pd.DataFrame, "trait_data")

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)

    # Calculate correlation matrix
    corr_matrix = trait_data.corr()

    # Plot heatmap
    if HAS_SEABORN:
        mask = np.triu(np.ones_like(corr_matrix, dtype=bool))
        sns.heatmap(corr_matrix, mask=mask, annot=True, fmt=".2f", cmap="coolwarm", center=0, ax=ax, **kwargs)
    else:
        im = ax.imshow(corr_matrix, cmap="coolwarm", aspect="equal", vmin=-1, vmax=1)
        ax.set_xticks(range(len(corr_matrix.columns)))
        ax.set_yticks(range(len(corr_matrix.index)))
        ax.set_xticklabels(corr_matrix.columns, rotation=45, ha="right")
        ax.set_yticklabels(corr_matrix.index)
        plt.colorbar(im, ax=ax)

    ax.set_title("Phenotypic Trait Correlations")

    if output_path:
        output_path = paths.ensure_directory(Path(output_path).parent)
        plt.savefig(output_path, dpi=300, bbox_inches="tight")
        logger.info(f"Trait correlation matrix saved to {output_path}")

    return ax


def plot_life_course_trajectory(
    life_events: List[Dict[str, Any]],
    individual_id: str | None = None,
    *,
    ax: Axes | None = None,
    output_path: str | Path | None = None,
    figsize: Tuple[float, float] = (12, 6),
    **kwargs,
) -> Axes:
    """Plot life course trajectory from event sequence.

    Args:
        life_events: List of life events with 'age' and 'event_type' keys
        individual_id: Optional individual identifier for labeling
        ax: Optional matplotlib axes to plot on
        output_path: Optional path to save the figure
        figsize: Figure size as (width, height)
        **kwargs: Additional arguments for customization

    Returns:
        matplotlib Axes object
    """
    validation.validate_type(life_events, list, "life_events")

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)

    # Extract ages and event types
    ages = []
    event_types = []
    descriptions = []

    for event in life_events:
        if "age" in event and "event_type" in event:
            ages.append(event["age"])
            event_types.append(event["event_type"])
            descriptions.append(event.get("description", ""))

    if not ages:
        raise ValueError("No valid life events with age information found")

    # Sort by age
    sorted_idx = np.argsort(ages)
    ages = np.array(ages)[sorted_idx]
    event_types = np.array(event_types)[sorted_idx]
    descriptions = np.array(descriptions)[sorted_idx]

    # Create timeline plot
    unique_events = np.unique(event_types)
    colors = plt.cm.tab10(np.linspace(0, 1, len(unique_events)))
    event_color_map = dict(zip(unique_events, colors))

    # Plot events as points on timeline
    for i, (age, event_type, desc) in enumerate(zip(ages, event_types, descriptions)):
        color = event_color_map[event_type]
        ax.scatter(age, 0, c=[color], s=100, alpha=0.8, edgecolors="black", zorder=5)

        # Add event label
        label = f"{event_type}"
        if desc:
            label += f": {desc}"
        ax.annotate(
            label,
            (age, 0),
            xytext=(0, 10 + i % 3 * 10),
            textcoords="offset points",
            ha="center",
            fontsize=8,
            bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.8),
        )

    # Add timeline
    ax.plot([ages.min(), ages.max()], [0, 0], "k-", linewidth=2, alpha=0.7)

    ax.set_xlabel("Age")
    ax.set_title(f'Life Course Trajectory{" - " + individual_id if individual_id else ""}')
    ax.set_yticks([])
    ax.set_xlim(ages.min() - 1, ages.max() + 1)
    ax.grid(True, alpha=0.3, axis="x")

    # Add legend
    legend_elements = [
        plt.scatter([], [], c=color, label=event, s=50, edgecolors="black") for event, color in event_color_map.items()
    ]
    ax.legend(handles=legend_elements, bbox_to_anchor=(1.05, 1), loc="upper left")

    if output_path:
        output_path = paths.ensure_directory(Path(output_path).parent)
        plt.savefig(output_path, dpi=300, bbox_inches="tight")
        logger.info(f"Life course trajectory saved to {output_path}")

    return ax


def plot_morphological_measurements(
    measurements: Dict[str, np.ndarray],
    *,
    ax: Axes | None = None,
    output_path: str | Path | None = None,
    figsize: Tuple[float, float] = (12, 8),
    **kwargs,
) -> Axes:
    """Plot multiple morphological measurements as box plots.

    Args:
        measurements: Dictionary mapping measurement names to value arrays
        ax: Optional matplotlib axes to plot on
        output_path: Optional path to save the figure
        figsize: Figure size as (width, height)
        **kwargs: Additional arguments for customization

    Returns:
        matplotlib Axes object
    """
    validation.validate_type(measurements, dict, "measurements")

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)

    # Prepare data for plotting
    measurement_names = list(measurements.keys())
    data = [measurements[name] for name in measurement_names]

    # Create box plot
    if HAS_SEABORN:
        df = pd.DataFrame({name: values for name, values in measurements.items()})
        df_melted = df.melt(var_name="Measurement", value_name="Value")
        sns.boxplot(data=df_melted, x="Measurement", y="Value", ax=ax, **kwargs)
        ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha="right")
    else:
        bp = ax.boxplot(data, labels=measurement_names, patch_artist=True, **kwargs)
        # Color boxes
        colors = plt.cm.Set3(np.linspace(0, 1, len(measurement_names)))
        for patch, color in zip(bp["boxes"], colors):
            patch.set_facecolor(color)

    ax.set_xlabel("Morphological Measurements")
    ax.set_ylabel("Value")
    ax.set_title("Morphological Trait Distributions")
    ax.grid(True, alpha=0.3, axis="y")

    if output_path:
        output_path = paths.ensure_directory(Path(output_path).parent)
        plt.savefig(output_path, dpi=300, bbox_inches="tight")
        logger.info(f"Morphological measurements plot saved to {output_path}")

    return ax


def plot_behavioral_patterns(
    behavioral_data: pd.DataFrame,
    time_column: str = "time",
    behavior_column: str = "behavior",
    *,
    ax: Axes | None = None,
    output_path: str | Path | None = None,
    figsize: Tuple[float, float] = (12, 6),
    **kwargs,
) -> Axes:
    """Plot behavioral patterns over time.

    Args:
        behavioral_data: DataFrame with time and behavior columns
        time_column: Name of time column
        behavior_column: Name of behavior column
        ax: Optional matplotlib axes to plot on
        output_path: Optional path to save the figure
        figsize: Figure size as (width, height)
        **kwargs: Additional arguments for customization

    Returns:
        matplotlib Axes object
    """
    validation.validate_type(behavioral_data, pd.DataFrame, "behavioral_data")

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)

    # Group by behavior and plot over time
    behaviors = behavioral_data[behavior_column].unique()
    colors = plt.cm.tab10(np.linspace(0, 1, len(behaviors)))

    for i, behavior in enumerate(behaviors):
        behavior_data = behavioral_data[behavioral_data[behavior_column] == behavior]
        times = behavior_data[time_column]

        # Create histogram of behavior timing
        ax.hist(times, bins=30, alpha=0.7, color=colors[i], label=str(behavior), density=True)

    ax.set_xlabel("Time")
    ax.set_ylabel("Behavior Frequency")
    ax.set_title("Behavioral Pattern Distribution")
    ax.legend()
    ax.grid(True, alpha=0.3)

    if output_path:
        output_path = paths.ensure_directory(Path(output_path).parent)
        plt.savefig(output_path, dpi=300, bbox_inches="tight")
        logger.info(f"Behavioral patterns plot saved to {output_path}")

    return ax


def plot_phenotype_pca(
    phenotype_data: np.ndarray,
    trait_names: List[str] | None = None,
    *,
    ax: Axes | None = None,
    output_path: str | Path | None = None,
    figsize: Tuple[float, float] = (8, 6),
    **kwargs,
) -> Axes:
    """Plot PCA of phenotypic trait data.

    Args:
        phenotype_data: Samples x traits matrix
        trait_names: Optional names for traits
        ax: Optional matplotlib axes to plot on
        output_path: Optional path to save the figure
        figsize: Figure size as (width, height)
        **kwargs: Additional arguments for customization

    Returns:
        matplotlib Axes object
    """
    validation.validate_type(phenotype_data, np.ndarray, "phenotype_data")

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)

    # Perform PCA
    from sklearn.decomposition import PCA

    pca = PCA(n_components=2)
    pca_coords = pca.fit_transform(phenotype_data)

    # Plot
    scatter = ax.scatter(
        pca_coords[:, 0], pca_coords[:, 1], alpha=0.7, c=range(len(pca_coords)), cmap="viridis", **kwargs
    )

    ax.set_xlabel(f"PC1 ({pca.explained_variance_ratio_[0]:.1%} variance)")
    ax.set_ylabel(f"PC2 ({pca.explained_variance_ratio_[1]:.1%} variance)")
    ax.set_title("Phenotype PCA")
    ax.grid(True, alpha=0.3)

    # Add colorbar
    cbar = plt.colorbar(scatter, ax=ax)
    cbar.set_label("Sample Index")

    if output_path:
        output_path = paths.ensure_directory(Path(output_path).parent)
        plt.savefig(output_path, dpi=300, bbox_inches="tight")
        logger.info(f"Phenotype PCA plot saved to {output_path}")

    return ax


def plot_trait_heritability(
    heritability_estimates: Dict[str, float],
    confidence_intervals: Dict[str, Tuple[float, float]] | None = None,
    *,
    ax: Axes | None = None,
    output_path: str | Path | None = None,
    figsize: Tuple[float, float] = (10, 6),
    **kwargs,
) -> Axes:
    """Plot trait heritability estimates.

    Args:
        heritability_estimates: Dictionary mapping trait names to heritability values
        confidence_intervals: Optional confidence intervals for each trait
        ax: Optional matplotlib axes to plot on
        output_path: Optional path to save the figure
        figsize: Figure size as (width, height)
        **kwargs: Additional arguments for customization

    Returns:
        matplotlib Axes object
    """
    validation.validate_type(heritability_estimates, dict, "heritability_estimates")

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)

    traits = list(heritability_estimates.keys())
    h2_values = list(heritability_estimates.values())

    # Plot heritability values
    bars = ax.bar(range(len(traits)), h2_values, color="skyblue", alpha=0.8, **kwargs)

    # Add confidence intervals if provided
    if confidence_intervals:
        for i, trait in enumerate(traits):
            if trait in confidence_intervals:
                lower, upper = confidence_intervals[trait]
                ax.errorbar(
                    i,
                    h2_values[i],
                    yerr=[[h2_values[i] - lower], [upper - h2_values[i]]],
                    fmt="none",
                    color="black",
                    capsize=3,
                    linewidth=1,
                )

    ax.set_xlabel("Traits")
    ax.set_ylabel("Heritability (h²)")
    ax.set_title("Trait Heritability Estimates")
    ax.set_xticks(range(len(traits)))
    ax.set_xticklabels(traits, rotation=45, ha="right")
    ax.set_ylim(0, 1)
    ax.grid(True, alpha=0.3, axis="y")

    # Add value labels on bars
    for bar, h2 in zip(bars, h2_values):
        ax.text(
            bar.get_x() + bar.get_width() / 2,
            bar.get_height() + 0.01,
            f"{h2:.2f}",
            ha="center",
            va="bottom",
            fontsize=8,
        )

    if output_path:
        output_path = paths.ensure_directory(Path(output_path).parent)
        plt.savefig(output_path, dpi=300, bbox_inches="tight")
        logger.info(f"Trait heritability plot saved to {output_path}")

    return ax


def plot_life_history_comparison(
    species_data: Dict[str, List[Dict[str, Any]]],
    *,
    ax: Axes | None = None,
    output_path: str | Path | None = None,
    figsize: Tuple[float, float] = (12, 8),
    **kwargs,
) -> Axes:
    """Compare life history patterns across species.

    Args:
        species_data: Dictionary mapping species names to life event lists
        ax: Optional matplotlib axes to plot on
        output_path: Optional path to save the figure
        figsize: Figure size as (width, height)
        **kwargs: Additional arguments for customization

    Returns:
        matplotlib Axes object
    """
    validation.validate_type(species_data, dict, "species_data")

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)

    species_names = list(species_data.keys())
    colors = plt.cm.Set1(np.linspace(0, 1, len(species_names)))

    # Plot life history curves for each species
    for i, (species, events) in enumerate(species_data.items()):
        ages = [event.get("age", 0) for event in events if "age" in event]
        if ages:
            # Create cumulative event plot
            sorted_ages = sorted(ages)
            cumulative = range(1, len(sorted_ages) + 1)
            ax.plot(sorted_ages, cumulative, "o-", color=colors[i], label=species, linewidth=2, markersize=4, alpha=0.8)

    ax.set_xlabel("Age")
    ax.set_ylabel("Cumulative Life Events")
    ax.set_title("Life History Comparison Across Species")
    ax.legend()
    ax.grid(True, alpha=0.3)

    if output_path:
        output_path = paths.ensure_directory(Path(output_path).parent)
        plt.savefig(output_path, dpi=300, bbox_inches="tight")
        logger.info(f"Life history comparison plot saved to {output_path}")

    return ax


def plot_phenotype_network(
    phenotype_correlations: np.ndarray,
    phenotype_names: List[str],
    *,
    ax: Axes | None = None,
    output_path: str | Path | None = None,
    figsize: Tuple[float, float] = (10, 8),
    **kwargs,
) -> Axes:
    """Plot phenotype correlation network.

    Args:
        phenotype_correlations: Correlation matrix between phenotypes
        phenotype_names: Names of phenotypes
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
        raise ImportError("networkx required for phenotype network visualization")

    validation.validate_type(phenotype_correlations, np.ndarray, "phenotype_correlations")

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)

    # Create network from correlation matrix
    threshold = kwargs.get("threshold", 0.5)
    G = nx.Graph()

    # Add nodes
    for i, name in enumerate(phenotype_names):
        G.add_node(i, name=name)

    # Add edges for strong correlations
    for i in range(len(phenotype_names)):
        for j in range(i + 1, len(phenotype_names)):
            if abs(phenotype_correlations[i, j]) > threshold:
                G.add_edge(i, j, weight=abs(phenotype_correlations[i, j]))

    # Plot network
    pos = nx.spring_layout(G, k=1, iterations=50)

    # Node colors based on degree
    degrees = dict(G.degree())
    node_colors = [degrees[node] for node in G.nodes()]

    nx.draw_networkx_nodes(G, pos, node_color=node_colors, node_size=300, cmap="viridis", ax=ax)
    nx.draw_networkx_edges(G, pos, alpha=0.5, ax=ax)
    nx.draw_networkx_labels(G, pos, {i: phenotype_names[i] for i in G.nodes()}, font_size=8, ax=ax)

    ax.set_title("Phenotype Correlation Network")
    ax.axis("off")

    if output_path:
        output_path = paths.ensure_directory(Path(output_path).parent)
        plt.savefig(output_path, dpi=300, bbox_inches="tight")
        logger.info(f"Phenotype network plot saved to {output_path}")

    return ax


def create_interactive_phenotype_browser(
    phenotype_data: pd.DataFrame, *, output_path: str | Path | None = None, **kwargs
) -> Any:
    """Create an interactive phenotype data browser using Plotly.

    Args:
        phenotype_data: DataFrame with phenotype measurements
        output_path: Optional path to save the HTML file
        **kwargs: Additional arguments for Plotly customization

    Returns:
        Plotly Figure object
    """
    if not HAS_PLOTLY:
        raise ImportError("Plotly required for interactive phenotype browser")

    validation.validate_type(phenotype_data, pd.DataFrame, "phenotype_data")

    # Create parallel coordinates plot
    fig = px.parallel_coordinates(phenotype_data, title="Interactive Phenotype Browser", **kwargs)

    if output_path:
        output_path = paths.ensure_directory(Path(output_path).parent)
        html_path = Path(output_path).with_suffix(".html")
        fig.write_html(str(html_path))
        logger.info(f"Interactive phenotype browser saved to {html_path}")

    return fig
