"""Phylogenetic tree visualization functions.

This module provides specialized plotting functions for phylogenetic trees
including basic tree layouts, circular arrangements, unrooted trees,
tree comparisons, and annotated tree visualizations.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, List, Tuple

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.axes import Axes

from metainformant.core.data import validation
from metainformant.core.io import paths
from metainformant.core.utils import logging

logger = logging.get_logger(__name__)

# Optional imports with graceful fallbacks
try:
    import networkx as nx

    HAS_NETWORKX = True
except ImportError:
    nx = None
    HAS_NETWORKX = False


def plot_phylo_tree(tree: Any, *, ax: Axes | None = None, output_path: str | Path | None = None, **kwargs) -> Axes:
    """Create a basic phylogenetic tree visualization.

    Args:
        tree: Tree object (NetworkX DiGraph or similar tree structure)
        ax: Optional matplotlib axes to plot on. Creates new if None.
        output_path: Optional path to save the figure.
        **kwargs: Additional arguments passed to networkx draw().

    Returns:
        matplotlib Axes object

    Raises:
        ImportError: If NetworkX is not available
        ValueError: If tree structure is invalid
    """
    if not HAS_NETWORKX:
        raise ImportError("NetworkX required for phylogenetic tree plotting")

    if ax is None:
        fig, ax = plt.subplots(figsize=kwargs.pop("figsize", (10, 8)))

    # Convert tree to NetworkX if needed
    if isinstance(tree, nx.DiGraph):
        G = tree
    else:
        # Assume tree is a dict-based structure or similar
        G = _convert_tree_to_networkx(tree)

    # Use hierarchical layout for phylogenetic trees
    pos = _hierarchical_tree_layout(G)

    # Draw tree
    nx.draw(
        G,
        pos=pos,
        ax=ax,
        with_labels=kwargs.get("with_labels", True),
        node_color=kwargs.get("node_color", "lightgray"),
        node_size=kwargs.get("node_size", 200),
        edge_color=kwargs.get("edge_color", "black"),
        width=kwargs.get("width", 1),
        arrows=kwargs.get("arrows", False),  # Usually no arrows in trees
        **kwargs,
    )

    ax.set_title("Phylogenetic Tree")
    ax.set_xlabel("Branch Length")
    ax.set_ylabel("Taxa")

    if output_path:
        paths.ensure_directory(Path(output_path).parent)
        plt.savefig(output_path, dpi=300, bbox_inches="tight")
        logger.info(f"Phylogenetic tree plot saved to {output_path}")

    return ax


def circular_tree_plot(tree: Any, *, ax: Axes | None = None, output_path: str | Path | None = None, **kwargs) -> Axes:
    """Create a circular phylogenetic tree visualization.

    Args:
        tree: Tree object (NetworkX DiGraph or similar tree structure)
        ax: Optional matplotlib axes to plot on. Creates new if None.
        output_path: Optional path to save the figure.
        **kwargs: Additional arguments passed to plotting functions.

    Returns:
        matplotlib Axes object

    Raises:
        ImportError: If NetworkX is not available
    """
    if not HAS_NETWORKX:
        raise ImportError("NetworkX required for phylogenetic tree plotting")

    if ax is None:
        fig, ax = plt.subplots(figsize=kwargs.pop("figsize", (10, 10)), subplot_kw={"projection": "polar"})

    # Convert tree to NetworkX if needed
    if isinstance(tree, nx.DiGraph):
        G = tree
    else:
        G = _convert_tree_to_networkx(tree)

    # Create circular layout
    pos = _circular_tree_layout(G)

    # Draw circular tree
    nx.draw(
        G,
        pos=pos,
        ax=ax,
        with_labels=kwargs.get("with_labels", True),
        node_color=kwargs.get("node_color", "lightblue"),
        node_size=kwargs.get("node_size", 150),
        edge_color=kwargs.get("edge_color", "gray"),
        width=kwargs.get("width", 1),
        arrows=kwargs.get("arrows", False),
        **kwargs,
    )

    ax.set_title("Circular Phylogenetic Tree")
    ax.set_rlabel_position(0)  # Move radial labels

    if output_path:
        paths.ensure_directory(Path(output_path).parent)
        plt.savefig(output_path, dpi=300, bbox_inches="tight")
        logger.info(f"Circular phylogenetic tree plot saved to {output_path}")

    return ax


def unrooted_tree_plot(tree: Any, *, ax: Axes | None = None, output_path: str | Path | None = None, **kwargs) -> Axes:
    """Create an unrooted phylogenetic tree visualization.

    Args:
        tree: Tree object (NetworkX Graph or similar tree structure)
        ax: Optional matplotlib axes to plot on. Creates new if None.
        output_path: Optional path to save the figure.
        **kwargs: Additional arguments passed to networkx draw().

    Returns:
        matplotlib Axes object

    Raises:
        ImportError: If NetworkX is not available
    """
    if not HAS_NETWORKX:
        raise ImportError("NetworkX required for phylogenetic tree plotting")

    if ax is None:
        fig, ax = plt.subplots(figsize=kwargs.pop("figsize", (10, 10)))

    # Convert tree to NetworkX if needed
    if isinstance(tree, nx.Graph):
        G = tree
    else:
        G = _convert_tree_to_networkx(tree, directed=False)

    # Use spring layout for unrooted trees (approximates unrooted layout)
    pos = nx.spring_layout(G, **kwargs.get("layout_kwargs", {"seed": 42}))

    # Draw unrooted tree
    nx.draw(
        G,
        pos=pos,
        ax=ax,
        with_labels=kwargs.get("with_labels", True),
        node_color=kwargs.get("node_color", "lightgreen"),
        node_size=kwargs.get("node_size", 200),
        edge_color=kwargs.get("edge_color", "black"),
        width=kwargs.get("width", 1.5),
        **kwargs,
    )

    ax.set_title("Unrooted Phylogenetic Tree")

    if output_path:
        paths.ensure_directory(Path(output_path).parent)
        plt.savefig(output_path, dpi=300, bbox_inches="tight")
        logger.info(f"Unrooted phylogenetic tree plot saved to {output_path}")

    return ax


def tree_comparison_plot(
    tree1: Any, tree2: Any, *, ax: Axes | None = None, output_path: str | Path | None = None, **kwargs
) -> Axes:
    """Create a side-by-side comparison of two phylogenetic trees.

    Args:
        tree1: First tree object
        tree2: Second tree object
        ax: Optional matplotlib axes to plot on. Creates new if None.
        output_path: Optional path to save the figure.
        **kwargs: Additional arguments passed to plotting functions.

    Returns:
        matplotlib Axes object

    Raises:
        ImportError: If NetworkX is not available
    """
    if not HAS_NETWORKX:
        raise ImportError("NetworkX required for phylogenetic tree plotting")

    if ax is None:
        fig, axes = plt.subplots(1, 2, figsize=kwargs.pop("figsize", (15, 8)))
    else:
        # Assume ax is a single axes, create subplots anyway
        fig = ax.figure
        axes = [fig.add_subplot(1, 2, 1), fig.add_subplot(1, 2, 2)]

    # Convert trees to NetworkX
    G1 = _convert_tree_to_networkx(tree1) if not isinstance(tree1, nx.DiGraph) else tree1
    G2 = _convert_tree_to_networkx(tree2) if not isinstance(tree2, nx.DiGraph) else tree2

    # Plot first tree
    pos1 = _hierarchical_tree_layout(G1)
    nx.draw(G1, pos=pos1, ax=axes[0], with_labels=True, node_color="lightblue", node_size=150, arrows=False)
    axes[0].set_title("Tree 1")

    # Plot second tree
    pos2 = _hierarchical_tree_layout(G2)
    nx.draw(G2, pos=pos2, ax=axes[1], with_labels=True, node_color="lightcoral", node_size=150, arrows=False)
    axes[1].set_title("Tree 2")

    plt.tight_layout()

    if output_path:
        paths.ensure_directory(Path(output_path).parent)
        plt.savefig(output_path, dpi=300, bbox_inches="tight")
        logger.info(f"Tree comparison plot saved to {output_path}")

    return axes[0]  # Return first axes for consistency


def tree_annotation_plot(
    tree: Any, annotations: Dict[str, Any], *, ax: Axes | None = None, output_path: str | Path | None = None, **kwargs
) -> Axes:
    """Create an annotated phylogenetic tree visualization.

    Args:
        tree: Tree object
        annotations: Dictionary mapping node names to annotation data
        ax: Optional matplotlib axes to plot on. Creates new if None.
        output_path: Optional path to save the figure.
        **kwargs: Additional arguments passed to plotting functions.

    Returns:
        matplotlib Axes object

    Raises:
        ImportError: If NetworkX is not available
    """
    if not HAS_NETWORKX:
        raise ImportError("NetworkX required for phylogenetic tree plotting")

    if ax is None:
        fig, ax = plt.subplots(figsize=kwargs.pop("figsize", (12, 10)))

    # Convert tree to NetworkX
    G = _convert_tree_to_networkx(tree) if not isinstance(tree, nx.DiGraph) else tree

    # Create node colors based on annotations
    node_colors = []
    for node in G.nodes():
        if node in annotations:
            # Use annotation to determine color (simplified)
            annot = annotations[node]
            if isinstance(annot, dict) and "color" in annot:
                node_colors.append(annot["color"])
            else:
                node_colors.append("lightgray")
        else:
            node_colors.append("lightgray")

    # Use hierarchical layout
    pos = _hierarchical_tree_layout(G)

    # Draw annotated tree
    nx.draw(
        G,
        pos=pos,
        ax=ax,
        with_labels=kwargs.get("with_labels", True),
        node_color=node_colors,
        node_size=kwargs.get("node_size", 250),
        edge_color=kwargs.get("edge_color", "black"),
        width=kwargs.get("width", 1.2),
        arrows=kwargs.get("arrows", False),
        **kwargs,
    )

    # Add annotation legends or labels
    if annotations:
        # Create a simple legend for annotation types
        unique_annotations = set()
        for annot in annotations.values():
            if isinstance(annot, dict):
                for key, value in annot.items():
                    if key != "color":
                        unique_annotations.add(f"{key}: {value}")

        # Add text annotations for nodes
        for node, annot in annotations.items():
            if node in pos and isinstance(annot, dict):
                x, y = pos[node]
                if "label" in annot:
                    ax.text(
                        x + 0.1,
                        y,
                        annot["label"],
                        fontsize=8,
                        bbox=dict(boxstyle="round,pad=0.2", facecolor="white", alpha=0.8),
                    )

    ax.set_title("Annotated Phylogenetic Tree")

    if output_path:
        paths.ensure_directory(Path(output_path).parent)
        plt.savefig(output_path, dpi=300, bbox_inches="tight")
        logger.info(f"Annotated phylogenetic tree plot saved to {output_path}")

    return ax


def _convert_tree_to_networkx(tree: Any, directed: bool = True) -> nx.Graph:
    """Convert various tree formats to NetworkX graph.

    This is a helper function to handle different tree representations.
    """
    if directed:
        G = nx.DiGraph()
    else:
        G = nx.Graph()

    # Handle dict-based tree representation
    if isinstance(tree, dict):

        def add_edges(parent, children):
            for child, subtree in children.items():
                G.add_edge(parent, child)
                if isinstance(subtree, dict):
                    add_edges(child, subtree)
                elif isinstance(subtree, list):
                    for grandchild in subtree:
                        G.add_edge(child, grandchild)

        # Find root (node with no parent in edges)
        all_nodes = set()

        def collect_nodes(node, subtree):
            all_nodes.add(node)
            if isinstance(subtree, dict):
                for child, subsubtree in subtree.items():
                    collect_nodes(child, subsubtree)
            elif isinstance(subtree, list):
                for child in subtree:
                    all_nodes.add(child)

        for root, subtree in tree.items():
            collect_nodes(root, subtree)

        # Add all edges
        for root, subtree in tree.items():
            add_edges(root, subtree)
    else:
        # Handle Bio.Phylo tree objects
        try:
            from Bio import Phylo as BioPhylo

            if isinstance(tree, BioPhylo.BaseTree.Tree):
                # Convert Bio.Phylo tree to NetworkX graph
                node_counter = 0

                def _add_clade(clade, parent_name=None):
                    nonlocal node_counter
                    name = clade.name if clade.name else f"internal_{node_counter}"
                    node_counter += 1
                    G.add_node(name)
                    if parent_name is not None:
                        weight = clade.branch_length if clade.branch_length else 1.0
                        G.add_edge(parent_name, name, weight=weight)
                    for child in clade.clades:
                        _add_clade(child, name)

                _add_clade(tree.root)
                return G
        except ImportError:
            pass

        raise ValueError("Unsupported tree format. Please provide NetworkX graph, dict-based tree, or Bio.Phylo tree.")

    return G


def _hierarchical_tree_layout(G: nx.DiGraph) -> Dict[str, Tuple[float, float]]:
    """Create a hierarchical layout for phylogenetic trees."""
    # Simple hierarchical layout
    pos = {}

    # Find root (node with no incoming edges)
    roots = [n for n in G.nodes() if G.in_degree(n) == 0]
    if not roots:
        # If no clear root, pick an arbitrary node
        roots = [list(G.nodes())[0]]

    root = roots[0]

    # Assign levels based on distance from root
    levels = {}
    for node in nx.topological_sort(G):
        if node == root:
            levels[node] = 0
        else:
            # Find parent (assuming tree structure)
            parents = list(G.predecessors(node))
            if parents:
                levels[node] = levels[parents[0]] + 1
            else:
                levels[node] = 0

    # Assign positions
    max_level = max(levels.values()) if levels else 0
    nodes_per_level = {}
    for node, level in levels.items():
        if level not in nodes_per_level:
            nodes_per_level[level] = []
        nodes_per_level[level].append(node)

    for level, nodes in nodes_per_level.items():
        y_spacing = 1.0 / (len(nodes) + 1)
        for i, node in enumerate(nodes):
            x_pos = level * 1.0  # Branch length (simplified)
            y_pos = (i + 1) * y_spacing
            pos[node] = (x_pos, y_pos)

    return pos


def _circular_tree_layout(G: nx.DiGraph) -> Dict[str, Tuple[float, float]]:
    """Create a circular layout for phylogenetic trees."""
    pos = {}

    # Get leaves (nodes with no outgoing edges)
    leaves = [n for n in G.nodes() if G.out_degree(n) == 0]

    # Arrange leaves in a circle
    n_leaves = len(leaves)
    angles = [2 * 3.14159 * i / n_leaves for i in range(n_leaves)]

    # Position leaves on circle
    radius = 1.0
    for leaf, angle in zip(leaves, angles):
        pos[leaf] = (radius * np.cos(angle), radius * np.sin(angle))

    # Position internal nodes (simplified)
    # This is a basic implementation - real phylogenetic layouts are more complex
    internal_nodes = [n for n in G.nodes() if n not in leaves]
    for node in internal_nodes:
        # Position at centroid of children
        children = list(G.successors(node))
        if children:
            child_positions = [pos[child] for child in children if child in pos]
            if child_positions:
                avg_x = np.mean([p[0] for p in child_positions])
                avg_y = np.mean([p[1] for p in child_positions])
                # Move internal nodes inward
                dist_from_center = np.sqrt(avg_x**2 + avg_y**2)
                if dist_from_center > 0:
                    scale = 0.7  # Move closer to center
                    pos[node] = (avg_x * scale, avg_y * scale)
                else:
                    pos[node] = (avg_x, avg_y)

    return pos
