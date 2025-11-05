"""Network visualization functions.

This module provides visualization functions for network graphs including
basic network plots, circular layouts, hierarchical layouts, force-directed
layouts, and community-based visualizations.
"""

from __future__ import annotations

from typing import Sequence

import matplotlib
import matplotlib.pyplot as plt
import numpy as np

# Use non-interactive backend by default for tests/headless
matplotlib.use("Agg", force=True)


def network_plot(
    nodes: list[str],
    edges: list[tuple[str, str]],
    *,
    node_sizes: list[float] | None = None,
    node_colors: list[str] | None = None,
    ax: plt.Axes | None = None,
    **kwargs
) -> plt.Axes:
    """Create a network graph visualization.

    Args:
        nodes: List of node names
        edges: List of (source, target) edge tuples
        node_sizes: Optional list of node sizes
        node_colors: Optional list of node colors
        ax: Matplotlib axes (creates new if None)
        **kwargs: Additional arguments for networkx.draw

    Returns:
        Matplotlib axes object

    Example:
        >>> from metainformant.visualization import network_plot
        >>> nodes = ['A', 'B', 'C']
        >>> edges = [('A', 'B'), ('B', 'C')]
        >>> ax = network_plot(nodes, edges)
    """
    try:
        import networkx as nx
    except ImportError:
        # Fallback to simple text representation
        if ax is None:
            _, ax = plt.subplots()
        ax.text(0.5, 0.5, "NetworkX not available\nInstall with: uv pip install networkx",
               ha="center", va="center", transform=ax.transAxes)
        return ax

    if ax is None:
        _, ax = plt.subplots()

    # Create graph
    G = nx.Graph()
    G.add_nodes_from(nodes)
    G.add_edges_from(edges)

    # Set default node sizes
    if node_sizes is None:
        node_sizes = [300] * len(nodes)

    # Set default node colors
    if node_colors is None:
        node_colors = ['lightblue'] * len(nodes)

    # Draw network
    pos = nx.spring_layout(G, k=1, iterations=50)
    nx.draw(G, pos, ax=ax, node_size=node_sizes, node_color=node_colors,
            with_labels=True, font_size=8, **kwargs)

    ax.set_title("Network Graph")

    return ax


def circular_network_plot(
    nodes: list[str],
    edges: list[tuple[str, str]],
    *,
    node_sizes: list[float] | None = None,
    node_colors: list[str] | None = None,
    ax: plt.Axes | None = None,
    **kwargs
) -> plt.Axes:
    """Create a circular network plot.

    Args:
        nodes: List of node names
        edges: List of (source, target) edge tuples
        node_sizes: Optional list of node sizes
        node_colors: Optional list of node colors
        ax: Matplotlib axes (creates new if None)
        **kwargs: Additional arguments

    Returns:
        Matplotlib axes object

    Example:
        >>> from metainformant.visualization import circular_network_plot
        >>> nodes = ['A', 'B', 'C', 'D']
        >>> edges = [('A', 'B'), ('B', 'C'), ('C', 'D'), ('D', 'A')]
        >>> ax = circular_network_plot(nodes, edges)
    """
    try:
        import networkx as nx
    except ImportError:
        if ax is None:
            _, ax = plt.subplots()
        ax.text(0.5, 0.5, "NetworkX not available", ha="center", va="center", transform=ax.transAxes)
        return ax

    if ax is None:
        _, ax = plt.subplots()

    G = nx.Graph()
    G.add_nodes_from(nodes)
    G.add_edges_from(edges)

    if node_sizes is None:
        node_sizes = [300] * len(nodes)
    if node_colors is None:
        node_colors = ['lightblue'] * len(nodes)

    pos = nx.circular_layout(G)
    nx.draw(G, pos, ax=ax, node_size=node_sizes, node_color=node_colors,
            with_labels=True, font_size=8, **kwargs)

    ax.set_title("Circular Network")
    ax.axis('off')

    return ax


def hierarchical_network_plot(
    nodes: list[str],
    edges: list[tuple[str, str]],
    root: str | None = None,
    *,
    node_sizes: list[float] | None = None,
    node_colors: list[str] | None = None,
    ax: plt.Axes | None = None,
    **kwargs
) -> plt.Axes:
    """Create a hierarchical network plot.

    Args:
        nodes: List of node names
        edges: List of (source, target) edge tuples
        root: Root node for hierarchy (if None, uses first node)
        node_sizes: Optional list of node sizes
        node_colors: Optional list of node colors
        ax: Matplotlib axes (creates new if None)
        **kwargs: Additional arguments

    Returns:
        Matplotlib axes object

    Example:
        >>> from metainformant.visualization import hierarchical_network_plot
        >>> nodes = ['A', 'B', 'C', 'D']
        >>> edges = [('A', 'B'), ('A', 'C'), ('C', 'D')]
        >>> ax = hierarchical_network_plot(nodes, edges, root='A')
    """
    try:
        import networkx as nx
    except ImportError:
        if ax is None:
            _, ax = plt.subplots()
        ax.text(0.5, 0.5, "NetworkX not available", ha="center", va="center", transform=ax.transAxes)
        return ax

    if ax is None:
        _, ax = plt.subplots()

    G = nx.DiGraph() if root else nx.Graph()
    G.add_nodes_from(nodes)
    G.add_edges_from(edges)

    if node_sizes is None:
        node_sizes = [300] * len(nodes)
    if node_colors is None:
        node_colors = ['lightblue'] * len(nodes)

    if root:
        pos = nx.spring_layout(G, k=2, iterations=50)
    else:
        pos = nx.spring_layout(G, k=1, iterations=50)

    nx.draw(G, pos, ax=ax, node_size=node_sizes, node_color=node_colors,
            with_labels=True, font_size=8, arrows=True, **kwargs)

    ax.set_title("Hierarchical Network")
    ax.axis('off')

    return ax


def force_directed_plot(
    nodes: list[str],
    edges: list[tuple[str, str]],
    *,
    node_sizes: list[float] | None = None,
    node_colors: list[str] | None = None,
    iterations: int = 50,
    ax: plt.Axes | None = None,
    **kwargs
) -> plt.Axes:
    """Create a force-directed network plot.

    Args:
        nodes: List of node names
        edges: List of (source, target) edge tuples
        node_sizes: Optional list of node sizes
        node_colors: Optional list of node colors
        iterations: Number of iterations for force-directed layout
        ax: Matplotlib axes (creates new if None)
        **kwargs: Additional arguments

    Returns:
        Matplotlib axes object

    Example:
        >>> from metainformant.visualization import force_directed_plot
        >>> nodes = ['A', 'B', 'C', 'D', 'E']
        >>> edges = [('A', 'B'), ('B', 'C'), ('C', 'D'), ('D', 'E'), ('E', 'A')]
        >>> ax = force_directed_plot(nodes, edges)
    """
    try:
        import networkx as nx
    except ImportError:
        if ax is None:
            _, ax = plt.subplots()
        ax.text(0.5, 0.5, "NetworkX not available", ha="center", va="center", transform=ax.transAxes)
        return ax

    if ax is None:
        _, ax = plt.subplots()

    G = nx.Graph()
    G.add_nodes_from(nodes)
    G.add_edges_from(edges)

    if node_sizes is None:
        node_sizes = [300] * len(nodes)
    if node_colors is None:
        node_colors = ['lightblue'] * len(nodes)

    pos = nx.spring_layout(G, k=1, iterations=iterations)
    nx.draw(G, pos, ax=ax, node_size=node_sizes, node_color=node_colors,
            with_labels=True, font_size=8, **kwargs)

    ax.set_title("Force-Directed Network")
    ax.axis('off')

    return ax


def community_network_plot(
    nodes: list[str],
    edges: list[tuple[str, str]],
    communities: dict[str, int] | None = None,
    *,
    node_sizes: list[float] | None = None,
    ax: plt.Axes | None = None,
    **kwargs
) -> plt.Axes:
    """Create a network plot colored by community membership.

    Args:
        nodes: List of node names
        edges: List of (source, target) edge tuples
        communities: Dictionary mapping node names to community IDs
        node_sizes: Optional list of node sizes
        ax: Matplotlib axes (creates new if None)
        **kwargs: Additional arguments

    Returns:
        Matplotlib axes object

    Example:
        >>> from metainformant.visualization import community_network_plot
        >>> nodes = ['A', 'B', 'C', 'D']
        >>> edges = [('A', 'B'), ('B', 'C'), ('C', 'D')]
        >>> communities = {'A': 0, 'B': 0, 'C': 1, 'D': 1}
        >>> ax = community_network_plot(nodes, edges, communities)
    """
    try:
        import networkx as nx
    except ImportError:
        if ax is None:
            _, ax = plt.subplots()
        ax.text(0.5, 0.5, "NetworkX not available", ha="center", va="center", transform=ax.transAxes)
        return ax

    if ax is None:
        _, ax = plt.subplots()

    G = nx.Graph()
    G.add_nodes_from(nodes)
    G.add_edges_from(edges)

    if node_sizes is None:
        node_sizes = [300] * len(nodes)

    # Assign colors based on communities
    if communities:
        unique_communities = sorted(set(communities.values()))
        colors = plt.cm.tab10(np.linspace(0, 1, len(unique_communities)))
        node_colors = [colors[communities.get(node, 0)] for node in nodes]
    else:
        node_colors = ['lightblue'] * len(nodes)

    pos = nx.spring_layout(G, k=1, iterations=50)
    nx.draw(G, pos, ax=ax, node_size=node_sizes, node_color=node_colors,
            with_labels=True, font_size=8, **kwargs)

    ax.set_title("Community Network")
    ax.axis('off')

    return ax

