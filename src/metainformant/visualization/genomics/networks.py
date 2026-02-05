"""Network visualization functions for graph analysis.

This module provides specialized plotting functions for biological networks
including basic layouts, circular arrangements, hierarchical structures,
force-directed layouts, and community-based visualizations.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict

import numpy as np
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


def plot_network_basic(G: Any, *, ax: Axes | None = None, output_path: str | Path | None = None, **kwargs) -> Axes:
    """Create a basic network visualization using spring layout.

    Args:
        G: NetworkX graph object
        ax: Optional matplotlib axes to plot on. Creates new if None.
        output_path: Optional path to save the figure.
        **kwargs: Additional arguments passed to networkx.draw().

    Returns:
        matplotlib Axes object

    Raises:
        ImportError: If NetworkX is not available
    """
    if not HAS_NETWORKX:
        raise ImportError("NetworkX required for network plotting")

    if ax is None:
        fig, ax = plt.subplots(figsize=kwargs.pop("figsize", (10, 8)))

    # Use spring layout for basic network visualization
    pos = nx.spring_layout(G, **kwargs.get("layout_kwargs", {}))

    nx.draw(
        G,
        pos=pos,
        ax=ax,
        with_labels=kwargs.get("with_labels", True),
        node_color=kwargs.get("node_color", "lightblue"),
        node_size=kwargs.get("node_size", 300),
        edge_color=kwargs.get("edge_color", "gray"),
        width=kwargs.get("width", 1),
        alpha=kwargs.get("alpha", 0.8),
        **kwargs,
    )

    ax.set_title("Network Visualization (Spring Layout)")

    if output_path:
        paths.ensure_directory(Path(output_path).parent)
        plt.savefig(output_path, dpi=300, bbox_inches="tight")
        logger.info(f"Basic network plot saved to {output_path}")

    return ax


def plot_network_circular(G: Any, *, ax: Axes | None = None, output_path: str | Path | None = None, **kwargs) -> Axes:
    """Create a circular network visualization.

    Args:
        G: NetworkX graph object
        ax: Optional matplotlib axes to plot on. Creates new if None.
        output_path: Optional path to save the figure.
        **kwargs: Additional arguments passed to networkx.draw().

    Returns:
        matplotlib Axes object

    Raises:
        ImportError: If NetworkX is not available
    """
    if not HAS_NETWORKX:
        raise ImportError("NetworkX required for network plotting")

    if ax is None:
        fig, ax = plt.subplots(figsize=kwargs.pop("figsize", (10, 10)))

    # Use circular layout
    pos = nx.circular_layout(G)

    nx.draw(
        G,
        pos=pos,
        ax=ax,
        with_labels=kwargs.get("with_labels", True),
        node_color=kwargs.get("node_color", "lightgreen"),
        node_size=kwargs.get("node_size", 400),
        edge_color=kwargs.get("edge_color", "gray"),
        width=kwargs.get("width", 1.5),
        alpha=kwargs.get("alpha", 0.8),
        **kwargs,
    )

    ax.set_title("Network Visualization (Circular Layout)")

    if output_path:
        paths.ensure_directory(Path(output_path).parent)
        plt.savefig(output_path, dpi=300, bbox_inches="tight")
        logger.info(f"Circular network plot saved to {output_path}")

    return ax


def plot_network_hierarchical(
    G: Any, *, ax: Axes | None = None, output_path: str | Path | None = None, **kwargs
) -> Axes:
    """Create a hierarchical network visualization.

    Args:
        G: NetworkX graph object (preferably a tree/DAG)
        ax: Optional matplotlib axes to plot on. Creates new if None.
        output_path: Optional path to save the figure.
        **kwargs: Additional arguments passed to networkx.draw().

    Returns:
        matplotlib Axes object

    Raises:
        ImportError: If NetworkX is not available
    """
    if not HAS_NETWORKX:
        raise ImportError("NetworkX required for network plotting")

    if ax is None:
        fig, ax = plt.subplots(figsize=kwargs.pop("figsize", (12, 8)))

    # Try to use hierarchical layout if available, otherwise use shell layout
    try:
        # For trees/DAGs, try to find a root and use a hierarchical approach
        if nx.is_directed(G):
            # Find root nodes (nodes with no incoming edges)
            roots = [n for n in G.nodes() if G.in_degree(n) == 0]
            if roots:
                # Use multipartite layout for hierarchical visualization
                # Assign levels based on distance from root
                root = roots[0]
                levels = {}
                for node in G.nodes():
                    try:
                        levels[node] = nx.shortest_path_length(G, root, node)
                    except nx.NetworkXNoPath:
                        levels[node] = len(G.nodes())  # Put unreachable nodes at the end

                pos = nx.multipartite_layout(G, subset_key=levels)
            else:
                pos = nx.shell_layout(G)
        else:
            # For undirected graphs, use shell layout as approximation
            pos = nx.shell_layout(G)
    except Exception:
        # Fallback to spring layout
        logger.warning("Could not create hierarchical layout, using spring layout")
        pos = nx.spring_layout(G)

    nx.draw(
        G,
        pos=pos,
        ax=ax,
        with_labels=kwargs.get("with_labels", True),
        node_color=kwargs.get("node_color", "lightcoral"),
        node_size=kwargs.get("node_size", 350),
        edge_color=kwargs.get("edge_color", "darkgray"),
        width=kwargs.get("width", 1.2),
        alpha=kwargs.get("alpha", 0.8),
        arrows=kwargs.get("arrows", True) if nx.is_directed(G) else False,
        **kwargs,
    )

    ax.set_title("Network Visualization (Hierarchical Layout)")

    if output_path:
        paths.ensure_directory(Path(output_path).parent)
        plt.savefig(output_path, dpi=300, bbox_inches="tight")
        logger.info(f"Hierarchical network plot saved to {output_path}")

    return ax


def plot_network_force_directed(
    G: Any, *, ax: Axes | None = None, output_path: str | Path | None = None, **kwargs
) -> Axes:
    """Create a force-directed network visualization.

    Args:
        G: NetworkX graph object
        ax: Optional matplotlib axes to plot on. Creates new if None.
        output_path: Optional path to save the figure.
        **kwargs: Additional arguments passed to networkx.draw().

    Returns:
        matplotlib Axes object

    Raises:
        ImportError: If NetworkX is not available
    """
    if not HAS_NETWORKX:
        raise ImportError("NetworkX required for network plotting")

    if ax is None:
        fig, ax = plt.subplots(figsize=kwargs.pop("figsize", (10, 8)))

    # Use Fruchterman-Reingold force-directed algorithm
    pos = nx.fruchterman_reingold_layout(G, **kwargs.get("layout_kwargs", {}))

    nx.draw(
        G,
        pos=pos,
        ax=ax,
        with_labels=kwargs.get("with_labels", False),  # Often too cluttered for force-directed
        node_color=kwargs.get("node_color", "lightskyblue"),
        node_size=kwargs.get("node_size", 250),
        edge_color=kwargs.get("edge_color", "gray"),
        width=kwargs.get("width", 0.8),
        alpha=kwargs.get("alpha", 0.8),
        **kwargs,
    )

    ax.set_title("Network Visualization (Force-Directed Layout)")

    if output_path:
        paths.ensure_directory(Path(output_path).parent)
        plt.savefig(output_path, dpi=300, bbox_inches="tight")
        logger.info(f"Force-directed network plot saved to {output_path}")

    return ax


def plot_community_network(
    G: Any, communities: Dict[str, int], *, ax: Axes | None = None, output_path: str | Path | None = None, **kwargs
) -> Axes:
    """Create a community-colored network visualization.

    Args:
        G: NetworkX graph object
        communities: Dictionary mapping node names to community IDs
        ax: Optional matplotlib axes to plot on. Creates new if None.
        output_path: Optional path to save the figure.
        **kwargs: Additional arguments passed to networkx.draw().

    Returns:
        matplotlib Axes object

    Raises:
        ImportError: If NetworkX is not available
        ValueError: If communities dict is invalid
    """
    if not HAS_NETWORKX:
        raise ImportError("NetworkX required for network plotting")

    validation.validate_type(communities, dict, "communities")

    if not communities:
        raise ValueError("Communities dictionary cannot be empty")

    if ax is None:
        fig, ax = plt.subplots(figsize=kwargs.pop("figsize", (10, 8)))

    # Create color map for communities
    unique_communities = set(communities.values())
    colors = plt.cm.tab10(np.linspace(0, 1, len(unique_communities)))
    community_colors = dict(zip(unique_communities, colors))

    # Assign colors to nodes
    node_colors = []
    for node in G.nodes():
        community_id = communities.get(node, -1)
        if community_id == -1:
            node_colors.append("gray")  # Unassigned nodes
        else:
            node_colors.append(community_colors[community_id])

    # Use spring layout for community visualization
    pos = nx.spring_layout(G, **kwargs.get("layout_kwargs", {}))

    nx.draw(
        G,
        pos=pos,
        ax=ax,
        with_labels=kwargs.get("with_labels", False),
        node_color=node_colors,
        node_size=kwargs.get("node_size", 300),
        edge_color=kwargs.get("edge_color", "lightgray"),
        width=kwargs.get("width", 1),
        alpha=kwargs.get("alpha", 0.8),
        **kwargs,
    )

    # Add legend
    legend_elements = [
        plt.Line2D([0], [0], marker="o", color="w", markerfacecolor=color, markersize=10, label=f"Community {comm_id}")
        for comm_id, color in community_colors.items()
    ]
    if any(isinstance(c, str) and c == "gray" for c in node_colors):
        legend_elements.append(
            plt.Line2D([0], [0], marker="o", color="w", markerfacecolor="gray", markersize=10, label="Unassigned")
        )

    ax.legend(handles=legend_elements, loc="best", fontsize=8)

    ax.set_title("Network Visualization (Community Colored)")

    if output_path:
        paths.ensure_directory(Path(output_path).parent)
        plt.savefig(output_path, dpi=300, bbox_inches="tight")
        logger.info(f"Community network plot saved to {output_path}")

    return ax
