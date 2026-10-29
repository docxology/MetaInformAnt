"""Network visualization functions for graph analysis.

This module provides specialized plotting functions for biological networks
including basic layouts, circular arrangements, hierarchical structures,
force-directed layouts, and community-based visualizations.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, List

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.axes import Axes

from metainformant.core.data import validation
from metainformant.core.io import paths
from metainformant.core.utils import logging
from metainformant.visualization.config.conventions import save_figure_deterministic
from metainformant.visualization.config.figures import root_figure

logger = logging.get_logger(__name__)

#: Seed of the force-directed layouts. ``spring_layout`` and ``fruchterman_reingold_layout`` start from random
#: positions, so without a seed the same graph was drawn differently on every call, which defeats the
#: byte-reproducible figure export; ``layout_kwargs={"seed": ...}`` replaces it.
DEFAULT_LAYOUT_SEED = 42

# Optional imports with graceful fallbacks
try:
    import networkx as nx

    HAS_NETWORKX = True
except ImportError:
    nx = None
    HAS_NETWORKX = False


def plot_network_basic(G: Any, *, ax: Axes | None = None, output_path: str | Path | None = None, **kwargs: Any) -> Axes:
    """Create a basic network visualization using spring layout.

    Args:
        G: NetworkX graph object
        ax: Optional matplotlib axes to plot on. Creates new if None.
        output_path: Optional path to save the figure.
        **kwargs: Additional arguments passed to networkx.draw(). ``layout_kwargs`` is passed to
            ``networkx.spring_layout`` on top of the default ``seed=42``, so the layout is reproducible.

    Returns:
        matplotlib Axes object

    Raises:
        ImportError: If NetworkX is not available
    """
    if not HAS_NETWORKX:
        raise ImportError("NetworkX required for network plotting")

    if ax is None:
        fig, ax = plt.subplots(figsize=kwargs.pop("figsize", (10, 8)))

    # Use spring layout for basic network visualization (seeded: see DEFAULT_LAYOUT_SEED)
    pos = nx.spring_layout(G, **{"seed": DEFAULT_LAYOUT_SEED, **kwargs.pop("layout_kwargs", {})})

    nx.draw(
        G,
        pos=pos,
        ax=ax,
        with_labels=kwargs.pop("with_labels", True),
        node_color=kwargs.pop("node_color", "lightblue"),
        node_size=kwargs.pop("node_size", 300),
        edge_color=kwargs.pop("edge_color", "gray"),
        width=kwargs.pop("width", 1),
        alpha=kwargs.pop("alpha", 0.8),
        **kwargs,
    )

    ax.set_title("Network Visualization (Spring Layout)")

    if output_path:
        paths.ensure_directory(Path(output_path).parent)
        save_figure_deterministic(root_figure(ax), output_path, dpi=300, bbox_inches="tight")
        logger.info(f"Basic network plot saved to {output_path}")

    return ax


def plot_network_circular(
    G: Any, *, ax: Axes | None = None, output_path: str | Path | None = None, **kwargs: Any
) -> Axes:
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
        with_labels=kwargs.pop("with_labels", True),
        node_color=kwargs.pop("node_color", "lightgreen"),
        node_size=kwargs.pop("node_size", 400),
        edge_color=kwargs.pop("edge_color", "gray"),
        width=kwargs.pop("width", 1.5),
        alpha=kwargs.pop("alpha", 0.8),
        **kwargs,
    )

    ax.set_title("Network Visualization (Circular Layout)")

    if output_path:
        paths.ensure_directory(Path(output_path).parent)
        save_figure_deterministic(root_figure(ax), output_path, dpi=300, bbox_inches="tight")
        logger.info(f"Circular network plot saved to {output_path}")

    return ax


def plot_network_hierarchical(
    G: Any, *, ax: Axes | None = None, output_path: str | Path | None = None, **kwargs: Any
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
                # Use multipartite layout for hierarchical visualization. A node's level is its distance from
                # the nearest root (every root starts level 0), so a graph with several roots is laid out
                # from all of them, not just the first
                layers: dict[int, list[Any]] = dict(enumerate(nx.bfs_layers(G, roots)))
                reached = {node for layer in layers.values() for node in layer}
                unreachable = [node for node in G.nodes() if node not in reached]
                if unreachable:
                    # Nodes no root reaches (e.g. a cycle with no root of its own) go in a final column
                    layers[len(layers)] = unreachable

                # multipartite_layout wants {level: nodes in that level}; handing it {node: level} made it
                # raise, and the blanket except below then always fell back to the spring layout
                pos = nx.multipartite_layout(G, subset_key=layers)
            else:
                pos = nx.shell_layout(G)
        else:
            # For undirected graphs, use shell layout as approximation
            pos = nx.shell_layout(G)
    except Exception:
        # Fallback to spring layout
        logger.warning("Could not create hierarchical layout, using spring layout")
        pos = nx.spring_layout(G, seed=DEFAULT_LAYOUT_SEED)

    nx.draw(
        G,
        pos=pos,
        ax=ax,
        with_labels=kwargs.pop("with_labels", True),
        node_color=kwargs.pop("node_color", "lightcoral"),
        node_size=kwargs.pop("node_size", 350),
        edge_color=kwargs.pop("edge_color", "darkgray"),
        width=kwargs.pop("width", 1.2),
        alpha=kwargs.pop("alpha", 0.8),
        arrows=kwargs.pop("arrows", True) if nx.is_directed(G) else False,
        **kwargs,
    )

    ax.set_title("Network Visualization (Hierarchical Layout)")

    if output_path:
        paths.ensure_directory(Path(output_path).parent)
        save_figure_deterministic(root_figure(ax), output_path, dpi=300, bbox_inches="tight")
        logger.info(f"Hierarchical network plot saved to {output_path}")

    return ax


def plot_network_force_directed(
    G: Any, *, ax: Axes | None = None, output_path: str | Path | None = None, **kwargs: Any
) -> Axes:
    """Create a force-directed network visualization.

    Args:
        G: NetworkX graph object
        ax: Optional matplotlib axes to plot on. Creates new if None.
        output_path: Optional path to save the figure.
        **kwargs: Additional arguments passed to networkx.draw(). ``layout_kwargs`` is passed to
            ``networkx.fruchterman_reingold_layout`` on top of the default ``seed=42``, so the layout is
            reproducible.

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
    pos = nx.fruchterman_reingold_layout(G, **{"seed": DEFAULT_LAYOUT_SEED, **kwargs.pop("layout_kwargs", {})})

    nx.draw(
        G,
        pos=pos,
        ax=ax,
        with_labels=kwargs.pop("with_labels", False),  # Often too cluttered for force-directed
        node_color=kwargs.pop("node_color", "lightskyblue"),
        node_size=kwargs.pop("node_size", 250),
        edge_color=kwargs.pop("edge_color", "gray"),
        width=kwargs.pop("width", 0.8),
        alpha=kwargs.pop("alpha", 0.8),
        **kwargs,
    )

    ax.set_title("Network Visualization (Force-Directed Layout)")

    if output_path:
        paths.ensure_directory(Path(output_path).parent)
        save_figure_deterministic(root_figure(ax), output_path, dpi=300, bbox_inches="tight")
        logger.info(f"Force-directed network plot saved to {output_path}")

    return ax


def plot_community_network(
    G: Any,
    communities: Dict[str, int],
    *,
    ax: Axes | None = None,
    output_path: str | Path | None = None,
    **kwargs: Any,
) -> Axes:
    """Create a community-colored network visualization.

    Communities are coloured in sorted order of their identifiers, so the same communities get the same colours
    whatever order the dictionary lists them in (and, for string identifiers, in every run). Nodes missing from
    ``communities`` are gray and listed as "Unassigned" in the legend.

    Args:
        G: NetworkX graph object
        communities: Dictionary mapping node names to community IDs
        ax: Optional matplotlib axes to plot on. Creates new if None.
        output_path: Optional path to save the figure.
        **kwargs: Additional arguments passed to networkx.draw(). ``layout_kwargs`` is passed to
            ``networkx.spring_layout`` on top of the default ``seed=42``, so the layout is reproducible.

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

    # Create color map for communities. A set of strings iterates in hash order, which changes with every
    # interpreter run, so the colours are assigned in sorted order instead.
    try:
        unique_communities = sorted(set(communities.values()))
    except TypeError:  # identifiers that cannot be compared with each other
        unique_communities = sorted(set(communities.values()), key=str)
    colors = plt.cm.tab10(np.linspace(0, 1, len(unique_communities)))
    community_colors: Dict[Any, Any] = dict(zip(unique_communities, colors))

    # Assign colors to nodes; a node the dictionary does not mention is unassigned (a real community may be called -1)
    node_colors: List[Any] = []
    for node in G.nodes():
        if node in communities:
            node_colors.append(community_colors[communities[node]])
        else:
            node_colors.append("gray")  # Unassigned nodes

    # Use spring layout for community visualization (seeded: see DEFAULT_LAYOUT_SEED)
    pos = nx.spring_layout(G, **{"seed": DEFAULT_LAYOUT_SEED, **kwargs.pop("layout_kwargs", {})})

    nx.draw(
        G,
        pos=pos,
        ax=ax,
        with_labels=kwargs.pop("with_labels", False),
        node_color=node_colors,
        node_size=kwargs.pop("node_size", 300),
        edge_color=kwargs.pop("edge_color", "lightgray"),
        width=kwargs.pop("width", 1),
        alpha=kwargs.pop("alpha", 0.8),
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
        save_figure_deterministic(root_figure(ax), output_path, dpi=300, bbox_inches="tight")
        logger.info(f"Community network plot saved to {output_path}")

    return ax
