"""Phylogenetic tree visualization functions.

This module provides specialized plotting functions for phylogenetic trees
including basic tree layouts, circular arrangements, unrooted trees,
tree comparisons, and annotated tree visualizations.

Every plot function accepts the same tree representations:

* a :class:`networkx.DiGraph` of parent -> child edges (an edge attribute ``weight`` is read as the branch
  length, a node attribute ``support`` as a support value);
* a ``dna.phylogeny`` tree dictionary, as returned by ``neighbor_joining_tree``, ``upgma_tree``,
  ``nj_tree_from_kmer``, ``from_newick``, ``prune_tree`` and ``bootstrap_support``: every node is a key, a
  leaf maps to ``None`` and an internal node maps to ``{child: branch_length}``. A ``"bootstrap"`` entry of an
  internal node is its support value, not a child;
* a nested dictionary ``{parent: {child: subtree}}`` whose leaves are ``None``, ``{}``, a branch length or a
  list of leaf names;
* a Newick string such as ``"((A:0.1,B:0.2):0.3,C:0.4);"``;
* a Biopython ``Bio.Phylo`` tree or clade.
"""

from __future__ import annotations

import math
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, cast

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.axes import Axes

from metainformant.core.io import paths
from metainformant.core.utils import logging
from metainformant.core.utils.newick import from_newick
from metainformant.visualization.config.conventions import save_figure_deterministic
from metainformant.visualization.config.figures import root_figure

logger = logging.get_logger(__name__)

# Optional imports with graceful fallbacks
try:
    import networkx as nx

    HAS_NETWORKX = True
except ImportError:
    nx = None
    HAS_NETWORKX = False

#: Key under which ``dna.phylogeny.bootstrap_support`` stores a node's support value.
_SUPPORT_KEY = "bootstrap"

_UNSUPPORTED_TREE = (
    "Unsupported tree format. Please provide a NetworkX DiGraph, a dict-based tree (a dna.phylogeny tree or a "
    "nested dict), a Newick string, or a Bio.Phylo tree."
)

_ChildMap = Dict[Any, List[Tuple[Any, Optional[float]]]]


def _save_tree_figure(ax: Axes, output_path: str | Path | None, label: str) -> None:
    """Write the figure *ax* is drawn on to *output_path* when one is given."""
    if output_path:
        paths.ensure_directory(Path(output_path).parent)
        save_figure_deterministic(root_figure(ax), output_path, dpi=300, bbox_inches="tight")
        logger.info(f"{label} saved to {output_path}")


def _draw_graph(ax: Axes, graph: Any, pos: Dict[Any, Any], defaults: Dict[str, Any], overrides: Dict[str, Any]) -> None:
    """Draw *graph* with networkx; caller keyword arguments override the defaults instead of colliding with them."""
    nx.draw(graph, pos=pos, ax=ax, **{**defaults, **overrides})


def plot_phylo_tree(tree: Any, *, ax: Axes | None = None, output_path: str | Path | None = None, **kwargs: Any) -> Axes:
    """Create a basic phylogenetic tree visualization (a rectangular phylogram).

    Leaves occupy consecutive rows in depth-first order and each internal node sits at the mean row of its
    children. A node's horizontal position is its distance from the root: the sum of the branch lengths along
    its path when every branch has a length (and they are not all zero), otherwise the number of edges. Negative
    branch lengths, which neighbor joining can produce, are drawn as zero-length branches. When internal nodes
    carry bootstrap support values they are printed beside the node.

    Args:
        tree: Tree in any representation listed in the module docstring, such as the dictionary returned by
            ``dna.phylogeny.neighbor_joining_tree``.
        ax: Optional matplotlib axes to plot on. Creates new if None.
        output_path: Optional path to save the figure.
        **kwargs: Additional arguments passed to networkx draw(). ``figsize`` sizes a new figure and
            ``show_support=False`` hides the support values.

    Returns:
        matplotlib Axes object

    Raises:
        ImportError: If NetworkX is not available
        ValueError: If the tree is in an unsupported format or has an undirected graph without a root
    """
    if not HAS_NETWORKX:
        raise ImportError("NetworkX required for phylogenetic tree plotting")

    show_support = kwargs.pop("show_support", True)
    if ax is None:
        fig, ax = plt.subplots(figsize=kwargs.pop("figsize", (10, 8)))

    G = _convert_tree_to_networkx(tree)
    pos = _hierarchical_tree_layout(G)

    _draw_graph(
        ax,
        G,
        pos,
        {
            "with_labels": True,
            "node_color": "lightgray",
            "node_size": 200,
            "edge_color": "black",
            "width": 1,
            "arrows": False,  # Usually no arrows in trees
        },
        kwargs,
    )
    _label_phylogram_axes(ax, G)
    if show_support:
        _annotate_support(ax, G, pos)

    ax.set_title("Phylogenetic Tree")

    _save_tree_figure(ax, output_path, "Phylogenetic tree plot")

    return ax


def circular_tree_plot(
    tree: Any, *, ax: Axes | None = None, output_path: str | Path | None = None, **kwargs: Any
) -> Axes:
    """Create a circular phylogenetic tree visualization.

    Leaves are spaced evenly on the unit circle and internal nodes are pulled toward the centre (see
    ``_circular_tree_layout``). Branch lengths are not drawn to scale. The default axes are polar and the node
    coordinates are converted to (angle, radius), so the leaves lie on a circle; a caller-supplied Cartesian
    axes receives the same layout in (x, y) with an equal aspect ratio.

    Args:
        tree: Tree in any representation listed in the module docstring.
        ax: Optional matplotlib axes to plot on, polar or Cartesian. Creates a polar axes if None.
        output_path: Optional path to save the figure.
        **kwargs: Additional arguments passed to networkx draw(). ``figsize`` sizes a new figure and
            ``show_support=False`` hides bootstrap support values.

    Returns:
        matplotlib Axes object

    Raises:
        ImportError: If NetworkX is not available
        ValueError: If the tree is in an unsupported format
    """
    if not HAS_NETWORKX:
        raise ImportError("NetworkX required for phylogenetic tree plotting")

    show_support = kwargs.pop("show_support", True)
    if ax is None:
        fig, ax = plt.subplots(figsize=kwargs.pop("figsize", (10, 10)), subplot_kw={"projection": "polar"})

    G = _convert_tree_to_networkx(tree)

    # The layout is Cartesian (leaves on the unit circle). A polar axes reads data coordinates as
    # (angle, radius), so hand it the polar form of the same points.
    pos = _circular_tree_layout(G)
    polar = ax.name == "polar"
    if polar:
        pos = {node: (math.atan2(y, x), math.hypot(x, y)) for node, (x, y) in pos.items()}

    _draw_graph(
        ax,
        G,
        pos,
        {
            "with_labels": True,
            "node_color": "lightblue",
            "node_size": 150,
            "edge_color": "gray",
            "width": 1,
            "arrows": False,
        },
        kwargs,
    )
    if polar:
        # nx.draw's margins push the radial origin below zero; pin it so the root sits at the centre
        cast("Any", ax).set_rlim(0.0, 1.1)
        cast("Any", ax).set_rlabel_position(0)
    else:
        ax.set_aspect("equal")
    if show_support:
        _annotate_support(ax, G, pos)

    ax.set_title("Circular Phylogenetic Tree")

    _save_tree_figure(ax, output_path, "Circular phylogenetic tree plot")

    return ax


def unrooted_tree_plot(
    tree: Any, *, ax: Axes | None = None, output_path: str | Path | None = None, **kwargs: Any
) -> Axes:
    """Create an unrooted phylogenetic tree visualization.

    A rooted tree given as a dictionary, Newick string or Bio.Phylo tree is unrooted first: a root with exactly
    two children is removed and its children are joined by one edge whose length is the sum of the two root
    branches. The layout is a force-directed spring layout, so edge lengths are not proportional to branch
    lengths; it shows the topology only.

    Args:
        tree: Tree in any representation listed in the module docstring (a NetworkX graph is drawn as given,
            without its edge directions).
        ax: Optional matplotlib axes to plot on. Creates new if None.
        output_path: Optional path to save the figure.
        **kwargs: Additional arguments passed to networkx draw(). ``layout_kwargs`` is passed to
            ``networkx.spring_layout`` on top of the default ``seed=42``, so the layout stays reproducible
            unless the caller changes the seed. ``figsize`` sizes a new figure and ``show_support=False`` hides
            bootstrap support values.

    Returns:
        matplotlib Axes object

    Raises:
        ImportError: If NetworkX is not available
        ValueError: If the tree is in an unsupported format
    """
    if not HAS_NETWORKX:
        raise ImportError("NetworkX required for phylogenetic tree plotting")

    show_support = kwargs.pop("show_support", True)
    if ax is None:
        fig, ax = plt.subplots(figsize=kwargs.pop("figsize", (10, 10)))

    G = _convert_tree_to_networkx(tree, directed=False)

    # Use spring layout for unrooted trees (a topology-only approximation of an unrooted layout). The seed is a
    # default that layout_kwargs may override, not something a caller passing other layout arguments loses.
    layout_kwargs = {"seed": 42, **kwargs.pop("layout_kwargs", {})}
    pos = nx.spring_layout(G, **layout_kwargs)

    _draw_graph(
        ax,
        G,
        pos,
        {"with_labels": True, "node_color": "lightgreen", "node_size": 200, "edge_color": "black", "width": 1.5},
        kwargs,
    )
    if show_support:
        _annotate_support(ax, G, pos)

    ax.set_title("Unrooted Phylogenetic Tree")

    _save_tree_figure(ax, output_path, "Unrooted phylogenetic tree plot")

    return ax


def tree_comparison_plot(
    tree1: Any, tree2: Any, *, ax: Axes | None = None, output_path: str | Path | None = None, **kwargs: Any
) -> Axes:
    """Create a side-by-side comparison of two phylogenetic trees.

    Both trees are drawn as rectangular phylograms (see :func:`plot_phylo_tree`). With no ``ax`` a new
    one-by-two figure is created. A caller-supplied ``ax`` is split into two panels (child axes inset in the left
    and right halves of ``ax``) and its own frame is hidden, so the comparison is drawn on ``ax`` rather than on
    new subplots stacked over it.

    Args:
        tree1: First tree, in any representation listed in the module docstring
        tree2: Second tree
        ax: Optional matplotlib axes to plot on. Creates new if None.
        output_path: Optional path to save the figure.
        **kwargs: Additional arguments passed to networkx draw() for both panels. ``figsize`` sizes a new
            figure and ``show_support=False`` hides bootstrap support values.

    Returns:
        The axes of the first (left) panel

    Raises:
        ImportError: If NetworkX is not available
        ValueError: If a tree is in an unsupported format
    """
    if not HAS_NETWORKX:
        raise ImportError("NetworkX required for phylogenetic tree plotting")

    show_support = kwargs.pop("show_support", True)
    panels: List[Axes]
    if ax is None:
        fig, axes_pair = plt.subplots(1, 2, figsize=kwargs.pop("figsize", (15, 8)))
        panels = list(axes_pair)
    else:
        ax.set_axis_off()
        panels = [ax.inset_axes((0.0, 0.0, 0.48, 1.0)), ax.inset_axes((0.52, 0.0, 0.48, 1.0))]

    for panel, tree, title, color in zip(panels, (tree1, tree2), ("Tree 1", "Tree 2"), ("lightblue", "lightcoral")):
        G = _convert_tree_to_networkx(tree)
        pos = _hierarchical_tree_layout(G)
        _draw_graph(
            panel,
            G,
            pos,
            {"with_labels": True, "node_color": color, "node_size": 150, "arrows": False},
            kwargs,
        )
        _label_phylogram_axes(panel, G)
        if show_support:
            _annotate_support(panel, G, pos)
        panel.set_title(title)

    if ax is None:
        root_figure(panels[0]).tight_layout()

    _save_tree_figure(panels[0], output_path, "Tree comparison plot")

    return panels[0]  # Return first axes for consistency


def tree_annotation_plot(
    tree: Any,
    annotations: Dict[str, Any] | None = None,
    *,
    ax: Axes | None = None,
    output_path: str | Path | None = None,
    **kwargs: Any,
) -> Axes:
    """Create an annotated phylogenetic tree visualization.

    The tree is drawn as in :func:`plot_phylo_tree`. An annotation is a dictionary: its ``"color"`` colours the
    node and its ``"label"`` is printed beside it. Nodes without an annotation use the ``node_color`` keyword
    (light gray by default). Annotation keys that are not nodes of the tree, and annotations that are not
    dictionaries, are logged as warnings and otherwise ignored.

    Args:
        tree: Tree in any representation listed in the module docstring
        annotations: Dictionary mapping node names to annotation dictionaries (no annotations if None)
        ax: Optional matplotlib axes to plot on. Creates new if None.
        output_path: Optional path to save the figure.
        **kwargs: Additional arguments passed to networkx draw(). ``figsize`` sizes a new figure and
            ``show_support=False`` hides bootstrap support values.

    Returns:
        matplotlib Axes object

    Raises:
        ImportError: If NetworkX is not available
        ValueError: If the tree is in an unsupported format
    """
    if not HAS_NETWORKX:
        raise ImportError("NetworkX required for phylogenetic tree plotting")

    show_support = kwargs.pop("show_support", True)
    annotations = annotations or {}
    if ax is None:
        fig, ax = plt.subplots(figsize=kwargs.pop("figsize", (12, 10)))

    G = _convert_tree_to_networkx(tree)

    unknown = [node for node in annotations if node not in G]
    if unknown:
        logger.warning(f"Ignoring annotations for nodes that are not in the tree: {unknown}")
    not_dicts = [node for node, annot in annotations.items() if node in G and not isinstance(annot, dict)]
    if not_dicts:
        logger.warning(f"Ignoring annotations that are not dictionaries (node: {not_dicts})")

    # Create node colors based on annotations
    default_color = kwargs.pop("node_color", "lightgray")
    node_colors = []
    for node in G.nodes():
        annot = annotations.get(node)
        if isinstance(annot, dict) and "color" in annot:
            node_colors.append(annot["color"])
        else:
            node_colors.append(default_color)

    pos = _hierarchical_tree_layout(G)

    _draw_graph(
        ax,
        G,
        pos,
        {
            "with_labels": True,
            "node_color": node_colors,
            "node_size": 250,
            "edge_color": "black",
            "width": 1.2,
            "arrows": False,
        },
        kwargs,
    )
    _label_phylogram_axes(ax, G)
    if show_support:
        _annotate_support(ax, G, pos)

    # Text annotations beside the nodes, offset by a fraction of the tree's width
    x_values = [x for x, _ in pos.values()]
    offset = 0.03 * ((max(x_values) - min(x_values)) if x_values else 0.0) or 0.1
    for node, annot in annotations.items():
        if node in pos and isinstance(annot, dict) and "label" in annot:
            x, y = pos[node]
            ax.text(
                x + offset,
                y,
                annot["label"],
                fontsize=8,
                bbox=dict(boxstyle="round,pad=0.2", facecolor="white", alpha=0.8),
            )

    ax.set_title("Annotated Phylogenetic Tree")

    _save_tree_figure(ax, output_path, "Annotated phylogenetic tree plot")

    return ax


# ---------------------------------------------------------------------------
# Tree conversion
# ---------------------------------------------------------------------------


def _branch_length(value: Any) -> Optional[float]:
    """Return *value* as a finite float branch length, or None when it is not a number."""
    if isinstance(value, bool) or not isinstance(value, (int, float, np.integer, np.floating)):
        return None
    number = float(value)
    return number if math.isfinite(number) else None


def _dict_tree_children(tree: Dict[Any, Any]) -> Tuple[List[Any], _ChildMap, Dict[Any, float]]:
    """Read a dict-based tree into ``(roots, {node: [(child, branch_length)]}, {node: support})``.

    One reader serves the flat ``dna.phylogeny`` layout (``{leaf: None, internal: {child: length}}``, every node a
    key) and nested trees (``{parent: {child: subtree}}``): a node's children are merged from every place the
    node is a key, and the roots are the top-level keys that are nobody's child. A child's value is its subtree
    (a dict or list of leaf names), a branch length (a number) or a leaf marker (``None`` or ``{}``). A
    ``"bootstrap"`` entry with a numeric value is the support of the node holding it, unless a node of that name
    exists.
    """
    children: _ChildMap = {}
    supports: Dict[Any, float] = {}
    child_names: set[Any] = set()
    node_names = set(tree)

    work: List[Tuple[Any, Any]] = list(reversed(list(tree.items())))
    while work:
        node, value = work.pop()
        entries = children.setdefault(node, [])
        if isinstance(value, dict):
            pairs = list(value.items())
        elif isinstance(value, list):
            pairs = [(item, None) for item in value]
        else:
            pairs = []
        subtrees = []
        for child, sub in pairs:
            support = _branch_length(sub)
            if child == _SUPPORT_KEY and _SUPPORT_KEY not in node_names and support is not None:
                supports[node] = support
                continue
            entries.append((child, support))
            child_names.add(child)
            if isinstance(sub, (dict, list)):
                subtrees.append((child, sub))
        work.extend(reversed(subtrees))  # keep sibling order depth-first

    roots = [node for node in tree if node not in child_names]
    return roots, children, supports


def _biopython_tree_children(tree: Any) -> Tuple[List[Any], _ChildMap, Dict[Any, float]]:
    """Read a ``Bio.Phylo`` tree or clade into ``(roots, {node: [(child, branch_length)]}, {node: support})``.

    Unnamed clades are called ``internal_<n>``. A branch length of ``0.0`` stays zero (only a missing length is
    unknown), and a clade's ``confidence`` becomes its support value.

    Raises:
        ValueError: If two clades share a name, which would merge them into one node.
    """
    root_clade = tree.root if hasattr(tree, "root") else tree
    used: set[str] = {str(clade.name) for clade in root_clade.find_clades() if clade.name}
    counter = 0

    def fresh_name() -> str:
        nonlocal counter
        while f"internal_{counter}" in used:
            counter += 1
        name = f"internal_{counter}"
        used.add(name)
        return name

    children: _ChildMap = {}
    supports: Dict[Any, float] = {}
    seen: set[str] = set()
    names: Dict[int, str] = {}

    def name_of(clade: Any) -> str:
        if id(clade) not in names:
            name = str(clade.name) if clade.name else fresh_name()
            if name in seen:
                raise ValueError(f"Bio.Phylo tree has more than one clade named {name!r}; node names must be unique")
            seen.add(name)
            names[id(clade)] = name
        return names[id(clade)]

    root_name = name_of(root_clade)
    stack = [root_clade]
    while stack:
        clade = stack.pop()
        name = name_of(clade)
        children[name] = [(name_of(child), _branch_length(child.branch_length)) for child in clade.clades]
        confidence = _branch_length(clade.confidence)
        if clade.clades and confidence is not None:
            supports[name] = confidence
        stack.extend(reversed(clade.clades))
    return [root_name], children, supports


def _graph_from_children(roots: List[Any], children: _ChildMap, supports: Dict[Any, float]) -> Any:
    """Build a directed graph with the nodes added in depth-first order, so a layout sees the tree order.

    Edges carry a ``weight`` attribute (the branch length) when it is known, and nodes a ``support`` attribute.
    """
    graph = nx.DiGraph()

    def add_branch(parent: Any, child: Any, length: Optional[float]) -> None:
        if length is None:
            graph.add_edge(parent, child)
        else:
            graph.add_edge(parent, child, weight=length)

    # A node unreachable from a root (part of a cycle) is rooted last so that it is drawn and the layout reports it
    for root in [*roots, *children]:
        if root in graph:
            continue
        graph.add_node(root)
        stack = [(root, iter(children.get(root, ())))]
        while stack:
            parent, entries = stack[-1]
            for child, length in entries:
                is_new = child not in graph
                add_branch(parent, child, length)
                if is_new:
                    stack.append((child, iter(children.get(child, ()))))
                    break
            else:
                stack.pop()

    for node, value in supports.items():
        graph.nodes[node]["support"] = value
    return graph


def _unroot(graph: Any) -> Any:
    """Return an undirected copy of the rooted tree *graph* with a two-child root removed.

    The root's two branches are merged into one edge whose length is their sum (unknown if either is unknown).
    """
    roots = [node for node in graph.nodes() if graph.in_degree(node) == 0]
    undirected = graph.to_undirected()
    if len(roots) == 1 and graph.out_degree(roots[0]) == 2:
        root = roots[0]
        first, second = graph.successors(root)
        lengths = [_branch_length(graph[root][child].get("weight")) for child in (first, second)]
        undirected.remove_node(root)
        if None in lengths:
            undirected.add_edge(first, second)
        else:
            undirected.add_edge(first, second, weight=sum(cast("List[float]", lengths)))
    return undirected


def _convert_tree_to_networkx(tree: Any, directed: bool = True) -> nx.Graph:
    """Convert any supported tree representation to a NetworkX graph.

    Supported are NetworkX graphs, ``dna.phylogeny`` tree dictionaries, nested dictionaries, Newick strings and
    ``Bio.Phylo`` trees or clades (see the module docstring). The graph has a ``weight`` edge attribute wherever
    a branch length is known and a ``support`` node attribute wherever a support value is. Nodes are added in
    depth-first order.

    Args:
        tree: The tree to convert.
        directed: Return a directed parent -> child graph. When False the tree is returned undirected and, for a
            dictionary, Newick or Bio.Phylo tree, a root with two children is removed (an unrooted tree has none).

    Raises:
        ValueError: If the tree is in an unsupported format, a Newick string is malformed, or an undirected
            NetworkX graph is given where a directed one is needed.
    """
    if isinstance(tree, nx.Graph):
        if directed and not tree.is_directed():
            raise ValueError(
                "A rooted tree plot needs a directed graph of parent -> child edges; "
                "an undirected graph has no root (use unrooted_tree_plot)"
            )
        return tree if tree.is_directed() == directed else tree.to_undirected()

    if isinstance(tree, str):
        if not tree.strip().startswith("("):
            raise ValueError(f"{_UNSUPPORTED_TREE} A string must be a Newick tree starting with '('.")
        roots, children, supports = _dict_tree_children(from_newick(tree))
    elif isinstance(tree, dict):
        roots, children, supports = _dict_tree_children(tree)
    else:
        try:
            from Bio.Phylo import BaseTree as bio_base_tree
        except ImportError:
            raise ValueError(_UNSUPPORTED_TREE) from None
        if not isinstance(tree, (bio_base_tree.Tree, bio_base_tree.Clade)):
            raise ValueError(_UNSUPPORTED_TREE)
        roots, children, supports = _biopython_tree_children(tree)

    graph = _graph_from_children(roots, children, supports)
    return graph if directed else _unroot(graph)


# ---------------------------------------------------------------------------
# Layouts and annotation helpers
# ---------------------------------------------------------------------------


def _edge_branch_lengths(G: Any) -> Optional[Dict[Tuple[Any, Any], float]]:
    """Return ``{(parent, child): length}`` when *G* can be drawn to scale, else None.

    A tree can be drawn to scale only when every edge has a finite ``weight`` and they do not all vanish (the
    lengths of a Newick tree without branch lengths are all zero).
    """
    lengths: Dict[Tuple[Any, Any], float] = {}
    for parent, child, data in G.edges(data=True):
        length = _branch_length(data.get("weight"))
        if length is None:
            return None
        lengths[(parent, child)] = length
    if not lengths or sum(lengths.values()) <= 0:
        return None
    return lengths


def _hierarchical_tree_layout(G: nx.DiGraph) -> Dict[Any, Tuple[float, float]]:
    """Create a rectangular phylogram layout for a rooted tree.

    The horizontal coordinate is the distance from the root: the sum of the ``weight`` branch lengths on the path
    when :func:`_edge_branch_lengths` finds a length for every edge, otherwise the number of edges (negative
    lengths count as zero). The vertical coordinate is the leaf order: leaves take consecutive rows ``0, 1, ...``
    in depth-first order and every internal node sits at the mean row of its children. Several roots are laid out
    one below another. A node with more than one parent is placed under its first parent.

    Raises:
        networkx.NetworkXUnfeasible: If *G* contains a cycle (it is not a tree).
    """
    if G.number_of_nodes() == 0:
        return {}
    if not nx.is_directed_acyclic_graph(G):
        raise nx.NetworkXUnfeasible("Graph contains a cycle, so it is not a tree")

    lengths = _edge_branch_lengths(G)
    if lengths is not None and any(length < 0 for length in lengths.values()):
        logger.warning(
            f"{sum(length < 0 for length in lengths.values())} negative branch length(s) are drawn as zero length"
        )

    x: Dict[Any, float] = {}
    y: Dict[Any, float] = {}
    placed_children: Dict[Any, List[Any]] = {}
    row = 0
    for root in [node for node in G.nodes() if G.in_degree(node) == 0]:
        x[root] = 0.0
        placed_children[root] = []
        stack = [(root, iter(G.successors(root)))]
        while stack:
            node, successors = stack[-1]
            for child in successors:
                if child in x:  # a second parent: the node stays under its first one
                    continue
                step = 1.0 if lengths is None else max(lengths[(node, child)], 0.0)
                x[child] = x[node] + step
                placed_children[node].append(child)
                placed_children[child] = []
                stack.append((child, iter(G.successors(child))))
                break
            else:
                stack.pop()
                kids = placed_children[node]
                if kids:
                    y[node] = sum(y[kid] for kid in kids) / len(kids)
                else:
                    y[node] = float(row)
                    row += 1

    return {node: (x[node], y[node]) for node in G.nodes()}


def _label_phylogram_axes(ax: Axes, G: Any) -> None:
    """Show the horizontal axis of a phylogram layout with a label that says what it measures.

    ``nx.draw`` switches the axes off, which hides the scale and any axis label, so the horizontal axis is turned
    back on and the vertical one (the arbitrary leaf order) stays hidden.
    """
    ax.set_axis_on()
    ax.tick_params(axis="x", bottom=True, labelbottom=True)
    ax.yaxis.set_visible(False)
    for side in ("left", "right", "top"):
        ax.spines[side].set_visible(False)
    ax.set_xlabel("Branch Length" if _edge_branch_lengths(G) is not None else "Depth (edges from root)")


def _annotate_support(ax: Axes, G: Any, pos: Dict[Any, Tuple[float, float]]) -> None:
    """Print each node's ``support`` value (bootstrap percentage) beside the node."""
    if not pos:
        return
    x_values = [x for x, _ in pos.values()]
    y_values = [y for _, y in pos.values()]
    dx = 0.01 * ((max(x_values) - min(x_values)) or 1.0)
    dy = 0.02 * ((max(y_values) - min(y_values)) or 1.0)
    for node, support in G.nodes(data="support"):
        if support is None or node not in pos:
            continue
        x, y = pos[node]
        text = f"{support:g}" if isinstance(support, float) else str(support)
        ax.text(x + dx, y + dy, text, fontsize=7, color="firebrick", ha="left", va="bottom")


def _circular_tree_layout(G: nx.DiGraph) -> Dict[str, Tuple[float, float]]:
    """Create a circular layout for phylogenetic trees.

    Leaves (nodes without children) are spaced evenly on the unit circle. Every
    internal node is placed at 0.7 times the centroid of its children, which
    pulls it toward the centre. Nodes are visited children-first, so a node is
    positioned only once all of its children are, whatever order they were added
    to *G* in. This is a basic layout - real phylogenetic layouts are more complex.

    Raises:
        networkx.NetworkXUnfeasible: If *G* contains a cycle (it is not a tree).
    """
    pos: Dict[str, Tuple[float, float]] = {}

    # Get leaves (nodes with no outgoing edges) and arrange them on the circle
    leaves = [n for n in G.nodes() if G.out_degree(n) == 0]
    n_leaves = len(leaves)
    for i, leaf in enumerate(leaves):
        angle = 2 * np.pi * i / n_leaves
        pos[leaf] = (float(np.cos(angle)), float(np.sin(angle)))

    # Reverse topological order puts every child before its parent, so each
    # internal node finds all of its children already positioned.
    for node in reversed(list(nx.topological_sort(G))):
        if node in pos:
            continue
        child_positions = [pos[child] for child in G.successors(node)]
        avg_x = float(np.mean([p[0] for p in child_positions]))
        avg_y = float(np.mean([p[1] for p in child_positions]))
        # Move internal nodes inward
        pos[node] = (avg_x * 0.7, avg_y * 0.7)

    return pos
