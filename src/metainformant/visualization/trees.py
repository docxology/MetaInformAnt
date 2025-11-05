"""Phylogenetic tree visualization functions.

This module provides visualization functions for phylogenetic trees including
basic tree plots, circular layouts, unrooted trees, tree comparisons, and
tree annotations.
"""

from __future__ import annotations

from typing import Any

import matplotlib
import matplotlib.pyplot as plt
import numpy as np

# Non-interactive for tests/headless
matplotlib.use("Agg", force=True)

try:  # optional Bio.Phylo ASCII exists already; we also expose matplotlib plot
    from Bio import Phylo
except Exception:  # pragma: no cover - optional runtime
    Phylo = None  # type: ignore


def plot_phylo_tree(tree: Any, *, ax: plt.Axes | None = None) -> plt.Axes:
    """Plot a Biopython Phylo tree to matplotlib Axes.

    Accepts any object compatible with Bio.Phylo.draw().

    Args:
        tree: Biopython Phylo tree object
        ax: Matplotlib axes (creates new if None)

    Returns:
        Matplotlib axes object

    Example:
        >>> from metainformant.visualization import plot_phylo_tree
        >>> from Bio import Phylo
        >>> tree = Phylo.read("tree.nwk", "newick")
        >>> ax = plot_phylo_tree(tree)
    """
    if Phylo is None:
        raise RuntimeError("Biopython Phylo not available")
    if ax is None:
        _, ax = plt.subplots()
    Phylo.draw(tree, do_show=False, axes=ax)
    return ax


def circular_tree_plot(tree: Any, *, ax: plt.Axes | None = None) -> plt.Axes:
    """Plot a phylogenetic tree in circular layout.

    Args:
        tree: Biopython Phylo tree object
        ax: Matplotlib axes (creates new if None, with polar projection)

    Returns:
        Matplotlib axes object

    Example:
        >>> from metainformant.visualization import circular_tree_plot
        >>> from Bio import Phylo
        >>> tree = Phylo.read("tree.nwk", "newick")
        >>> ax = circular_tree_plot(tree)
    """
    if Phylo is None:
        raise RuntimeError("Biopython Phylo not available")
    if ax is None:
        fig, ax = plt.subplots(figsize=(10, 10), subplot_kw=dict(projection='polar'))
    else:
        # Ensure polar projection
        if not hasattr(ax, 'set_theta_zero_location'):
            fig = ax.figure
            ax.remove()
            ax = fig.add_subplot(111, projection='polar')

    Phylo.draw(tree, do_show=False, axes=ax)
    return ax


def unrooted_tree_plot(tree: Any, *, ax: plt.Axes | None = None) -> plt.Axes:
    """Plot an unrooted phylogenetic tree.

    Args:
        tree: Biopython Phylo tree object
        ax: Matplotlib axes (creates new if None)

    Returns:
        Matplotlib axes object

    Example:
        >>> from metainformant.visualization import unrooted_tree_plot
        >>> from Bio import Phylo
        >>> tree = Phylo.read("tree.nwk", "newick")
        >>> ax = unrooted_tree_plot(tree)
    """
    if Phylo is None:
        raise RuntimeError("Biopython Phylo not available")
    if ax is None:
        _, ax = plt.subplots(figsize=(10, 10))

    # Use equal aspect for unrooted trees
    ax.set_aspect('equal')
    Phylo.draw(tree, do_show=False, axes=ax)
    ax.axis('off')
    return ax


def tree_comparison_plot(
    trees: list[Any],
    labels: list[str] | None = None,
    *,
    ncols: int = 2,
    figsize: tuple[int, int] | None = None,
) -> tuple[plt.Figure, np.ndarray]:
    """Plot multiple trees side by side for comparison.

    Args:
        trees: List of Biopython Phylo tree objects
        labels: Optional labels for each tree
        ncols: Number of columns in subplot grid
        figsize: Figure size tuple

    Returns:
        Tuple of (figure, axes array)

    Example:
        >>> from metainformant.visualization import tree_comparison_plot
        >>> from Bio import Phylo
        >>> tree1 = Phylo.read("tree1.nwk", "newick")
        >>> tree2 = Phylo.read("tree2.nwk", "newick")
        >>> fig, axes = tree_comparison_plot([tree1, tree2], labels=['Tree 1', 'Tree 2'])
    """
    if Phylo is None:
        raise RuntimeError("Biopython Phylo not available")

    n_trees = len(trees)
    nrows = (n_trees + ncols - 1) // ncols

    if figsize is None:
        figsize = (6 * ncols, 4 * nrows)

    fig, axes = plt.subplots(nrows, ncols, figsize=figsize)

    if nrows == 1 and ncols == 1:
        axes = np.array([[axes]])
    elif nrows == 1:
        axes = axes.reshape(1, -1)
    elif ncols == 1:
        axes = axes.reshape(-1, 1)

    for i, tree in enumerate(trees):
        row = i // ncols
        col = i % ncols
        ax = axes[row, col] if nrows > 1 else axes[col]
        Phylo.draw(tree, do_show=False, axes=ax)
        if labels and i < len(labels):
            ax.set_title(labels[i])

    # Hide unused subplots
    for i in range(n_trees, nrows * ncols):
        row = i // ncols
        col = i % ncols
        ax = axes[row, col] if nrows > 1 else axes[col]
        ax.set_visible(False)

    plt.tight_layout()
    return fig, axes


def tree_annotation_plot(
    tree: Any,
    annotations: dict[str, str] | None = None,
    *,
    ax: plt.Axes | None = None,
    **kwargs
) -> plt.Axes:
    """Plot a phylogenetic tree with annotations.

    Args:
        tree: Biopython Phylo tree object
        annotations: Dictionary mapping node names to annotation text
        ax: Matplotlib axes (creates new if None)
        **kwargs: Additional arguments

    Returns:
        Matplotlib axes object

    Example:
        >>> from metainformant.visualization import tree_annotation_plot
        >>> from Bio import Phylo
        >>> tree = Phylo.read("tree.nwk", "newick")
        >>> annotations = {'leaf1': 'Important', 'leaf2': 'Reference'}
        >>> ax = tree_annotation_plot(tree, annotations)
    """
    if Phylo is None:
        raise RuntimeError("Biopython Phylo not available")
    if ax is None:
        _, ax = plt.subplots()

    Phylo.draw(tree, do_show=False, axes=ax)

    # Add annotations if provided
    if annotations:
        # This is a simplified version - full implementation would
        # require tree traversal to find node positions
        for i, (name, text) in enumerate(annotations.items()):
            # Placeholder annotation - would need tree position calculation
            ax.text(0.1, 0.9 - i * 0.1, f"{name}: {text}",
                   transform=ax.transAxes, fontsize=8, **kwargs)

    return ax
