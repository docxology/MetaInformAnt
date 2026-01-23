"""Ecology and community analysis visualization functions.

This module provides comprehensive visualization capabilities for ecological data,
including biodiversity metrics, community structure, species abundance distributions,
and ecological network analysis.
"""

from __future__ import annotations

import warnings
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple, Union

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.axes import Axes
from matplotlib.patches import Rectangle, Circle
import pandas as pd

from metainformant.core import logging, paths, validation

logger = logging.get_logger(__name__)

try:
    import seaborn as sns

    HAS_SEABORN = True
except ImportError:
    HAS_SEABORN = False
    sns = None

try:
    import plotly.graph_objects as go
    import plotly.express as px

    HAS_PLOTLY = True
except ImportError:
    HAS_PLOTLY = False
    go = None
    px = None


def plot_species_abundance_distribution(
    abundance_data: np.ndarray,
    species_names: List[str] | None = None,
    *,
    ax: Axes | None = None,
    output_path: str | Path | None = None,
    figsize: Tuple[float, float] = (10, 6),
    **kwargs,
) -> Axes:
    """Plot species abundance distribution (rank-abundance curve).

    Args:
        abundance_data: Array of species abundances
        species_names: Optional names for species
        ax: Optional matplotlib axes to plot on
        output_path: Optional path to save the figure
        figsize: Figure size as (width, height)
        **kwargs: Additional arguments for customization

    Returns:
        matplotlib Axes object
    """
    validation.validate_type(abundance_data, np.ndarray, "abundance_data")

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)

    # Sort abundances in descending order
    sorted_abundances = np.sort(abundance_data)[::-1]
    ranks = np.arange(1, len(sorted_abundances) + 1)

    # Plot rank-abundance curve
    ax.plot(ranks, sorted_abundances, "b-o", linewidth=2, markersize=4, alpha=0.8)

    # Add log scale option
    if kwargs.get("log_scale", False):
        ax.set_yscale("log")
        ax.set_ylabel("Abundance (log scale)")
    else:
        ax.set_ylabel("Abundance")

    ax.set_xlabel("Species Rank")
    ax.set_title("Species Abundance Distribution")
    ax.grid(True, alpha=0.3)

    # Add species labels for top species
    if species_names and len(species_names) > 0:
        sorted_indices = np.argsort(abundance_data)[::-1]
        sorted_names = [species_names[i] for i in sorted_indices]

        # Label top 5 species
        for i in range(min(5, len(sorted_names))):
            ax.annotate(
                sorted_names[i],
                (i + 1, sorted_abundances[i]),
                xytext=(5, 5),
                textcoords="offset points",
                fontsize=8,
                bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.8),
            )

    if output_path:
        output_path = paths.ensure_directory(Path(output_path).parent)
        plt.savefig(output_path, dpi=300, bbox_inches="tight")
        logger.info(f"Species abundance distribution saved to {output_path}")

    return ax


def plot_diversity_accumulation_curve(
    diversity_data: List[Dict[str, Any]],
    *,
    ax: Axes | None = None,
    output_path: str | Path | None = None,
    figsize: Tuple[float, float] = (10, 6),
    **kwargs,
) -> Axes:
    """Plot species accumulation curve (rarefaction curve).

    Args:
        diversity_data: List of diversity measurements at different sample sizes
        ax: Optional matplotlib axes to plot on
        output_path: Optional path to save the figure
        figsize: Figure size as (width, height)
        **kwargs: Additional arguments for customization

    Returns:
        matplotlib Axes object
    """
    validation.validate_type(diversity_data, list, "diversity_data")

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)

    # Extract data
    sample_sizes = [d.get("sample_size", i + 1) for i, d in enumerate(diversity_data)]
    species_counts = [d.get("species_count", 0) for d in diversity_data]
    confidence_lower = [d.get("confidence_lower") for d in diversity_data]
    confidence_upper = [d.get("confidence_upper") for d in diversity_data]

    # Plot accumulation curve
    ax.plot(sample_sizes, species_counts, "b-", linewidth=2, label="Observed", alpha=0.8)

    # Add confidence intervals if available
    if confidence_lower[0] is not None and confidence_upper[0] is not None:
        ax.fill_between(sample_sizes, confidence_lower, confidence_upper, alpha=0.3, color="blue", label="95% CI")

    ax.set_xlabel("Number of Samples")
    ax.set_ylabel("Number of Species")
    ax.set_title("Species Accumulation Curve")
    ax.legend()
    ax.grid(True, alpha=0.3)

    if output_path:
        output_path = paths.ensure_directory(Path(output_path).parent)
        plt.savefig(output_path, dpi=300, bbox_inches="tight")
        logger.info(f"Diversity accumulation curve saved to {output_path}")

    return ax


def plot_community_composition(
    community_matrix: np.ndarray,
    species_names: List[str] | None = None,
    sample_names: List[str] | None = None,
    *,
    ax: Axes | None = None,
    output_path: str | Path | None = None,
    figsize: Tuple[float, float] = (12, 8),
    **kwargs,
) -> Axes:
    """Plot community composition as stacked bar chart.

    Args:
        community_matrix: Samples x species abundance matrix
        species_names: Optional names for species
        sample_names: Optional names for samples
        ax: Optional matplotlib axes to plot on
        output_path: Optional path to save the figure
        figsize: Figure size as (width, height)
        **kwargs: Additional arguments for customization

    Returns:
        matplotlib Axes object
    """
    validation.validate_type(community_matrix, np.ndarray, "community_matrix")

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)

    # Normalize to proportions
    row_sums = community_matrix.sum(axis=1, keepdims=True)
    proportions = community_matrix / row_sums

    # Plot stacked bars
    bottom = np.zeros(community_matrix.shape[0])
    colors = plt.cm.Set3(np.linspace(0, 1, community_matrix.shape[1]))

    for i in range(community_matrix.shape[1]):
        ax.bar(
            range(community_matrix.shape[0]),
            proportions[:, i],
            bottom=bottom,
            color=colors[i],
            alpha=0.8,
            label=species_names[i] if species_names else f"Species {i+1}",
        )
        bottom += proportions[:, i]

    ax.set_xlabel("Samples")
    ax.set_ylabel("Relative Abundance")
    ax.set_title("Community Composition")

    if sample_names:
        ax.set_xticks(range(len(sample_names)))
        ax.set_xticklabels(sample_names, rotation=45, ha="right")

    ax.legend(bbox_to_anchor=(1.05, 1), loc="upper left")
    ax.grid(True, alpha=0.3, axis="y")

    if output_path:
        output_path = paths.ensure_directory(Path(output_path).parent)
        plt.savefig(output_path, dpi=300, bbox_inches="tight")
        logger.info(f"Community composition plot saved to {output_path}")

    return ax


def plot_beta_diversity_ordination(
    ordination_coords: np.ndarray,
    sample_groups: np.ndarray | None = None,
    *,
    ax: Axes | None = None,
    output_path: str | Path | None = None,
    figsize: Tuple[float, float] = (8, 6),
    **kwargs,
) -> Axes:
    """Plot beta diversity ordination (PCoA/NMDS).

    Args:
        ordination_coords: Ordination coordinates (samples x dimensions)
        sample_groups: Optional grouping variable for coloring
        ax: Optional matplotlib axes to plot on
        output_path: Optional path to save the figure
        figsize: Figure size as (width, height)
        **kwargs: Additional arguments for customization

    Returns:
        matplotlib Axes object
    """
    validation.validate_type(ordination_coords, np.ndarray, "ordination_coords")

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)

    if ordination_coords.shape[1] < 2:
        raise ValueError("Ordination coordinates must have at least 2 dimensions")

    # Plot ordination
    if sample_groups is not None:
        unique_groups = np.unique(sample_groups)
        colors = plt.cm.tab10(np.linspace(0, 1, len(unique_groups)))

        for i, group in enumerate(unique_groups):
            mask = sample_groups == group
            ax.scatter(
                ordination_coords[mask, 0], ordination_coords[mask, 1], c=[colors[i]], label=str(group), alpha=0.7, s=50
            )
    else:
        scatter = ax.scatter(
            ordination_coords[:, 0],
            ordination_coords[:, 1],
            c=range(len(ordination_coords)),
            cmap="viridis",
            alpha=0.7,
            s=50,
        )

    ax.set_xlabel("Axis 1")
    ax.set_ylabel("Axis 2")
    ax.set_title("Beta Diversity Ordination")
    ax.grid(True, alpha=0.3)

    if sample_groups is not None:
        ax.legend()
    else:
        # Add colorbar
        cbar = plt.colorbar(scatter, ax=ax)
        cbar.set_label("Sample Index")

    if output_path:
        output_path = paths.ensure_directory(Path(output_path).parent)
        plt.savefig(output_path, dpi=300, bbox_inches="tight")
        logger.info(f"Beta diversity ordination saved to {output_path}")

    return ax


def plot_diversity_indices_comparison(
    diversity_indices: Dict[str, np.ndarray],
    index_names: List[str] | None = None,
    *,
    ax: Axes | None = None,
    output_path: str | Path | None = None,
    figsize: Tuple[float, float] = (10, 6),
    **kwargs,
) -> Axes:
    """Plot comparison of diversity indices across samples.

    Args:
        diversity_indices: Dictionary mapping index names to values
        index_names: Optional names for samples/sites
        ax: Optional matplotlib axes to plot on
        output_path: Optional path to save the figure
        figsize: Figure size as (width, height)
        **kwargs: Additional arguments for customization

    Returns:
        matplotlib Axes object
    """
    validation.validate_type(diversity_indices, dict, "diversity_indices")

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)

    # Prepare data
    index_types = list(diversity_indices.keys())
    n_samples = len(diversity_indices[index_types[0]])

    # Plot each diversity index
    x_positions = np.arange(n_samples)
    width = 0.8 / len(index_types)

    colors = plt.cm.Set1(np.linspace(0, 1, len(index_types)))

    for i, index_name in enumerate(index_types):
        values = diversity_indices[index_name]
        ax.bar(
            x_positions + i * width - width * len(index_types) / 2,
            values,
            width=width,
            color=colors[i],
            alpha=0.8,
            label=index_name,
        )

    ax.set_xlabel("Samples")
    ax.set_ylabel("Diversity Index Value")
    ax.set_title("Diversity Indices Comparison")

    if index_names:
        ax.set_xticks(x_positions)
        ax.set_xticklabels(index_names, rotation=45, ha="right")

    ax.legend()
    ax.grid(True, alpha=0.3, axis="y")

    if output_path:
        output_path = paths.ensure_directory(Path(output_path).parent)
        plt.savefig(output_path, dpi=300, bbox_inches="tight")
        logger.info(f"Diversity indices comparison saved to {output_path}")

    return ax


def plot_ecological_network(
    interaction_matrix: np.ndarray,
    species_names: List[str] | None = None,
    *,
    ax: Axes | None = None,
    output_path: str | Path | None = None,
    figsize: Tuple[float, float] = (10, 8),
    **kwargs,
) -> Axes:
    """Plot ecological interaction network.

    Args:
        interaction_matrix: Species x species interaction matrix
        species_names: Optional names for species
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
        raise ImportError("networkx required for ecological network visualization")

    validation.validate_type(interaction_matrix, np.ndarray, "interaction_matrix")

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)

    # Create network from interaction matrix
    G = nx.DiGraph()
    n_species = interaction_matrix.shape[0]

    # Add nodes
    for i in range(n_species):
        node_name = species_names[i] if species_names else f"Species {i+1}"
        G.add_node(i, name=node_name)

    # Add edges for interactions above threshold
    threshold = kwargs.get("threshold", 0.1)
    for i in range(n_species):
        for j in range(n_species):
            if abs(interaction_matrix[i, j]) > threshold and i != j:
                G.add_edge(i, j, weight=abs(interaction_matrix[i, j]))

    # Plot network
    pos = nx.spring_layout(G, k=1, iterations=50)

    # Node sizes based on degree
    degrees = dict(G.degree())
    node_sizes = [300 + degrees[node] * 50 for node in G.nodes()]

    nx.draw_networkx_nodes(G, pos, node_size=node_sizes, node_color="lightblue", ax=ax)
    nx.draw_networkx_edges(G, pos, edge_color="gray", alpha=0.5, arrows=True, arrowsize=10, ax=ax)

    # Add labels for important nodes only
    if len(G.nodes()) <= 20:
        labels = {i: species_names[i] if species_names else f"S{i+1}" for i in G.nodes()}
        nx.draw_networkx_labels(G, pos, labels, font_size=8, ax=ax)

    ax.set_title("Ecological Interaction Network")
    ax.axis("off")

    if output_path:
        output_path = paths.ensure_directory(Path(output_path).parent)
        plt.savefig(output_path, dpi=300, bbox_inches="tight")
        logger.info(f"Ecological network plot saved to {output_path}")

    return ax


def plot_rank_abundance_curve_comparison(
    abundance_datasets: Dict[str, np.ndarray],
    *,
    ax: Axes | None = None,
    output_path: str | Path | None = None,
    figsize: Tuple[float, float] = (10, 6),
    **kwargs,
) -> Axes:
    """Compare rank-abundance curves across different communities.

    Args:
        abundance_datasets: Dictionary mapping dataset names to abundance arrays
        ax: Optional matplotlib axes to plot on
        output_path: Optional path to save the figure
        figsize: Figure size as (width, height)
        **kwargs: Additional arguments for customization

    Returns:
        matplotlib Axes object
    """
    validation.validate_type(abundance_datasets, dict, "abundance_datasets")

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)

    colors = plt.cm.tab10(np.linspace(0, 1, len(abundance_datasets)))

    for i, (name, abundances) in enumerate(abundance_datasets.items()):
        # Sort abundances in descending order
        sorted_abundances = np.sort(abundances)[::-1]
        ranks = np.arange(1, len(sorted_abundances) + 1)

        ax.plot(ranks, sorted_abundances, "o-", color=colors[i], linewidth=2, markersize=3, alpha=0.8, label=name)

    ax.set_xlabel("Species Rank")
    ax.set_ylabel("Abundance")
    ax.set_title("Rank-Abundance Curve Comparison")
    ax.legend()
    ax.grid(True, alpha=0.3)

    # Set log scale if requested
    if kwargs.get("log_scale", False):
        ax.set_yscale("log")
        ax.set_ylabel("Abundance (log scale)")

    if output_path:
        output_path = paths.ensure_directory(Path(output_path).parent)
        plt.savefig(output_path, dpi=300, bbox_inches="tight")
        logger.info(f"Rank-abundance comparison saved to {output_path}")

    return ax


def plot_biodiversity_rarefaction(
    rarefaction_data: Dict[str, List[Dict[str, Any]]],
    *,
    ax: Axes | None = None,
    output_path: str | Path | None = None,
    figsize: Tuple[float, float] = (10, 6),
    **kwargs,
) -> Axes:
    """Plot biodiversity rarefaction curves for multiple communities.

    Args:
        rarefaction_data: Dictionary mapping community names to rarefaction data
        ax: Optional matplotlib axes to plot on
        output_path: Optional path to save the figure
        figsize: Figure size as (width, height)
        **kwargs: Additional arguments for customization

    Returns:
        matplotlib Axes object
    """
    validation.validate_type(rarefaction_data, dict, "rarefaction_data")

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)

    colors = plt.cm.Set1(np.linspace(0, 1, len(rarefaction_data)))

    for i, (community_name, data) in enumerate(rarefaction_data.items()):
        sample_sizes = [d.get("sample_size", 0) for d in data]
        species_counts = [d.get("species_count", 0) for d in data]

        ax.plot(
            sample_sizes,
            species_counts,
            "o-",
            color=colors[i],
            linewidth=2,
            markersize=4,
            alpha=0.8,
            label=community_name,
        )

    ax.set_xlabel("Number of Individuals Sampled")
    ax.set_ylabel("Number of Species")
    ax.set_title("Biodiversity Rarefaction Curves")
    ax.legend()
    ax.grid(True, alpha=0.3)

    if output_path:
        output_path = paths.ensure_directory(Path(output_path).parent)
        plt.savefig(output_path, dpi=300, bbox_inches="tight")
        logger.info(f"Biodiversity rarefaction plot saved to {output_path}")

    return ax


def plot_ecological_distance_heatmap(
    distance_matrix: np.ndarray,
    sample_names: List[str] | None = None,
    *,
    ax: Axes | None = None,
    output_path: str | Path | None = None,
    figsize: Tuple[float, float] = (10, 8),
    **kwargs,
) -> Axes:
    """Plot ecological distance/dissimilarity heatmap.

    Args:
        distance_matrix: Square distance matrix between samples
        sample_names: Optional names for samples
        ax: Optional matplotlib axes to plot on
        output_path: Optional path to save the figure
        figsize: Figure size as (width, height)
        **kwargs: Additional arguments for customization

    Returns:
        matplotlib Axes object
    """
    validation.validate_type(distance_matrix, np.ndarray, "distance_matrix")

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)

    # Plot heatmap
    im = ax.imshow(distance_matrix, cmap="YlOrRd", aspect="equal", origin="lower")

    ax.set_title("Ecological Distance Matrix")

    if sample_names and len(sample_names) <= 20:
        ax.set_xticks(range(len(sample_names)))
        ax.set_yticks(range(len(sample_names)))
        ax.set_xticklabels(sample_names, rotation=45, ha="right")
        ax.set_yticklabels(sample_names)

    # Add colorbar
    cbar = plt.colorbar(im, ax=ax)
    cbar.set_label("Ecological Distance")

    if output_path:
        output_path = paths.ensure_directory(Path(output_path).parent)
        plt.savefig(output_path, dpi=300, bbox_inches="tight")
        logger.info(f"Ecological distance heatmap saved to {output_path}")

    return ax


def create_interactive_ecology_dashboard(
    ecology_data: Dict[str, Any], *, output_path: str | Path | None = None, **kwargs
) -> Any:
    """Create an interactive ecology dashboard using Plotly.

    Args:
        ecology_data: Dictionary containing ecological data and metrics
        output_path: Optional path to save the HTML file
        **kwargs: Additional arguments for Plotly customization

    Returns:
        Plotly Figure object
    """
    if not HAS_PLOTLY:
        raise ImportError("Plotly required for interactive ecology dashboard")

    validation.validate_type(ecology_data, dict, "ecology_data")

    # Create subplot figure
    fig = go.Figure()

    # Add diversity indices subplot
    if "diversity_indices" in ecology_data:
        diversity_data = ecology_data["diversity_indices"]
        for i, (index_name, values) in enumerate(diversity_data.items()):
            fig.add_trace(go.Bar(name=index_name, x=list(range(len(values))), y=values, offsetgroup=i))

    fig.update_layout(title="Interactive Ecology Dashboard", barmode="group", **kwargs)

    if output_path:
        output_path = paths.ensure_directory(Path(output_path).parent)
        html_path = Path(output_path).with_suffix(".html")
        fig.write_html(str(html_path))
        logger.info(f"Interactive ecology dashboard saved to {html_path}")

    return fig
