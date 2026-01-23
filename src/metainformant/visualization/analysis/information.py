"""Information theory visualization functions.

This module provides specialized plotting functions for information-theoretic analysis
including entropy profiles, mutual information matrices, Rényi spectra,
information landscapes, and information flow networks.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, List

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.axes import Axes

from metainformant.core import logging, paths, validation

logger = logging.get_logger(__name__)

# Optional imports with graceful fallbacks
try:
    import networkx as nx

    HAS_NETWORKX = True
except ImportError:
    nx = None
    HAS_NETWORKX = False


def plot_entropy_profile(
    entropy_data: dict[str, List[float]], *, ax: Axes | None = None, output_path: str | Path | None = None, **kwargs
) -> Axes:
    """Create an entropy profile plot across different scales.

    Args:
        entropy_data: Dictionary with entropy values for different k-mer sizes or scales
        ax: Optional matplotlib axes to plot on. Creates new if None.
        output_path: Optional path to save the figure.
        **kwargs: Additional arguments passed to matplotlib plot().

    Returns:
        matplotlib Axes object

    Raises:
        ValueError: If entropy_data is empty or malformed
    """
    validation.validate_type(entropy_data, dict, "entropy_data")

    if not entropy_data:
        raise ValueError("Entropy data dictionary cannot be empty")

    if ax is None:
        fig, ax = plt.subplots(figsize=kwargs.pop("figsize", (10, 6)))

    # Plot entropy profiles for each dataset/sequence
    for label, entropy_values in entropy_data.items():
        if not isinstance(entropy_values, (list, np.ndarray)):
            raise ValueError(f"Entropy values for {label} must be a list or array")

        k_values = list(range(1, len(entropy_values) + 1))
        ax.plot(k_values, entropy_values, label=label, marker="o", **kwargs)

    ax.set_xlabel("k-mer size")
    ax.set_ylabel("Shannon Entropy (bits)")
    ax.set_title("Entropy Profile")
    ax.legend()
    ax.grid(True, alpha=0.3)

    if output_path:
        paths.ensure_directory(Path(output_path).parent)
        plt.savefig(output_path, dpi=300, bbox_inches="tight")
        logger.info(f"Entropy profile plot saved to {output_path}")

    return ax


def plot_mutual_information_matrix(
    mi_matrix: np.ndarray,
    labels: List[str] | None = None,
    *,
    ax: Axes | None = None,
    output_path: str | Path | None = None,
    **kwargs,
) -> Axes:
    """Create a mutual information matrix heatmap.

    Args:
        mi_matrix: Square matrix of mutual information values
        labels: Labels for matrix rows/columns
        ax: Optional matplotlib axes to plot on. Creates new if None.
        output_path: Optional path to save the figure.
        **kwargs: Additional arguments passed to matplotlib imshow().

    Returns:
        matplotlib Axes object

    Raises:
        ValueError: If mi_matrix is not square or labels don't match dimensions
    """
    validation.validate_type(mi_matrix, np.ndarray, "mi_matrix")

    if mi_matrix.ndim != 2 or mi_matrix.shape[0] != mi_matrix.shape[1]:
        raise ValueError("Mutual information matrix must be square")

    if labels is not None:
        validation.validate_type(labels, list, "labels")
        if len(labels) != mi_matrix.shape[0]:
            raise ValueError("Number of labels must match matrix dimensions")

    if ax is None:
        fig, ax = plt.subplots(figsize=kwargs.pop("figsize", (10, 8)))

    # Create heatmap
    im = ax.imshow(mi_matrix, cmap=kwargs.get("cmap", "viridis"), **kwargs)
    plt.colorbar(im, ax=ax)

    # Add labels if provided
    if labels:
        ax.set_xticks(np.arange(len(labels)))
        ax.set_yticks(np.arange(len(labels)))
        ax.set_xticklabels(labels, rotation=45, ha="right")
        ax.set_yticklabels(labels)

    ax.set_title("Mutual Information Matrix")
    ax.grid(False)  # Turn off grid for cleaner heatmap

    if output_path:
        paths.ensure_directory(Path(output_path).parent)
        plt.savefig(output_path, dpi=300, bbox_inches="tight")
        logger.info(f"Mutual information matrix plot saved to {output_path}")

    return ax


def plot_renyi_spectra(
    renyi_data: dict[str, List[float]],
    alpha_values: List[float],
    *,
    ax: Axes | None = None,
    output_path: str | Path | None = None,
    **kwargs,
) -> Axes:
    """Create a Rényi entropy spectra plot.

    Args:
        renyi_data: Dictionary with Rényi entropy values for different datasets
        alpha_values: List of alpha values used for Rényi entropy calculation
        ax: Optional matplotlib axes to plot on. Creates new if None.
        output_path: Optional path to save the figure.
        **kwargs: Additional arguments passed to matplotlib plot().

    Returns:
        matplotlib Axes object

    Raises:
        ValueError: If data dimensions don't match
    """
    validation.validate_type(renyi_data, dict, "renyi_data")
    validation.validate_type(alpha_values, list, "alpha_values")

    if not renyi_data:
        raise ValueError("Rényi data dictionary cannot be empty")

    if ax is None:
        fig, ax = plt.subplots(figsize=kwargs.pop("figsize", (10, 6)))

    # Plot Rényi spectra for each dataset
    for label, renyi_values in renyi_data.items():
        if len(renyi_values) != len(alpha_values):
            raise ValueError(f"Rényi values for {label} must match alpha values length")

        ax.plot(alpha_values, renyi_values, label=label, marker="o", **kwargs)

    ax.set_xlabel("α (Rényi parameter)")
    ax.set_ylabel("Rényi Entropy")
    ax.set_title("Rényi Entropy Spectra")
    ax.legend()
    ax.grid(True, alpha=0.3)

    # Highlight Shannon entropy (α=1)
    if 1.0 in alpha_values:
        ax.axvline(x=1.0, color="red", linestyle="--", alpha=0.7, label="Shannon entropy (α=1)")
        ax.legend()

    if output_path:
        paths.ensure_directory(Path(output_path).parent)
        plt.savefig(output_path, dpi=300, bbox_inches="tight")
        logger.info(f"Rényi spectra plot saved to {output_path}")

    return ax


def plot_information_landscape(
    landscape_data: np.ndarray,
    x_coords: np.ndarray | None = None,
    y_coords: np.ndarray | None = None,
    *,
    ax: Axes | None = None,
    output_path: str | Path | None = None,
    **kwargs,
) -> Axes:
    """Create an information landscape plot.

    Args:
        landscape_data: 2D array of information values
        x_coords: X coordinates for the landscape
        y_coords: Y coordinates for the landscape
        ax: Optional matplotlib axes to plot on. Creates new if None.
        output_path: Optional path to save the figure.
        **kwargs: Additional arguments passed to matplotlib contourf().

    Returns:
        matplotlib Axes object

    Raises:
        ValueError: If landscape_data dimensions are incorrect
    """
    validation.validate_type(landscape_data, np.ndarray, "landscape_data")

    if landscape_data.ndim != 2:
        raise ValueError("Information landscape data must be 2D")

    if ax is None:
        fig, ax = plt.subplots(figsize=kwargs.pop("figsize", (10, 8)))

    # Create coordinate grids if not provided
    if x_coords is None:
        x_coords = np.arange(landscape_data.shape[1])
    if y_coords is None:
        y_coords = np.arange(landscape_data.shape[0])

    X, Y = np.meshgrid(x_coords, y_coords)

    # Create contour plot
    cs = ax.contourf(
        X, Y, landscape_data, cmap=kwargs.get("cmap", "viridis"), levels=kwargs.get("levels", 20), **kwargs
    )
    plt.colorbar(cs, ax=ax)

    ax.set_xlabel("X coordinate")
    ax.set_ylabel("Y coordinate")
    ax.set_title("Information Landscape")

    if output_path:
        paths.ensure_directory(Path(output_path).parent)
        plt.savefig(output_path, dpi=300, bbox_inches="tight")
        logger.info(f"Information landscape plot saved to {output_path}")

    return ax


def plot_information_network(
    nodes: List[str],
    edges: List[tuple[str, str, float]],
    *,
    ax: Axes | None = None,
    output_path: str | Path | None = None,
    **kwargs,
) -> Axes:
    """Create an information flow network visualization.

    Args:
        nodes: List of node names
        edges: List of (source, target, weight) tuples for edges
        ax: Optional matplotlib axes to plot on. Creates new if None.
        output_path: Optional path to save the figure.
        **kwargs: Additional arguments passed to networkx draw().

    Returns:
        matplotlib Axes object

    Raises:
        ImportError: If NetworkX is not available
        ValueError: If nodes or edges are malformed
    """
    if not HAS_NETWORKX:
        raise ImportError("NetworkX required for information network plotting")

    validation.validate_type(nodes, list, "nodes")
    validation.validate_type(edges, list, "edges")

    if not nodes:
        raise ValueError("Nodes list cannot be empty")

    if ax is None:
        fig, ax = plt.subplots(figsize=kwargs.pop("figsize", (10, 8)))

    # Create network
    G = nx.DiGraph()
    G.add_nodes_from(nodes)

    # Add weighted edges
    for source, target, weight in edges:
        if source not in nodes or target not in nodes:
            raise ValueError(f"Edge connects unknown nodes: {source} -> {target}")
        G.add_edge(source, target, weight=weight)

    # Calculate node positions
    pos = nx.spring_layout(G, **kwargs.get("layout_kwargs", {}))

    # Draw network
    node_sizes = [G.degree(node) * 100 + 300 for node in G.nodes()]  # Size by degree
    edge_weights = [G[u][v]["weight"] * 2 for u, v in G.edges()]  # Width by weight

    nx.draw(
        G,
        pos=pos,
        ax=ax,
        with_labels=kwargs.get("with_labels", True),
        node_color=kwargs.get("node_color", "lightblue"),
        node_size=node_sizes,
        edge_color=kwargs.get("edge_color", "gray"),
        width=edge_weights,
        alpha=kwargs.get("alpha", 0.8),
        arrows=kwargs.get("arrows", True),
        **kwargs,
    )

    ax.set_title("Information Flow Network")

    if output_path:
        paths.ensure_directory(Path(output_path).parent)
        plt.savefig(output_path, dpi=300, bbox_inches="tight")
        logger.info(f"Information network plot saved to {output_path}")

    return ax
