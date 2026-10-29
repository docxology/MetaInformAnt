"""The tree plots must draw the trees ``dna.phylogeny`` builds, not only NetworkX graphs.

``neighbor_joining_tree`` and ``upgma_tree`` return a flat dictionary in which every node is a key, a leaf maps to
``None`` and an internal node maps to ``{child: branch_length}``; ``bootstrap_support`` adds a ``"bootstrap"`` entry
to each internal node. ``_convert_tree_to_networkx`` read every dictionary as a nested ``{parent: {child: subtree}}``
tree and called ``None.items()`` on the first leaf, so ``plot_phylo_tree`` and the four other tree plots raised
``AttributeError`` for every tree the package builds. The layout also placed nodes by the number of edges from the
root, ignoring branch lengths, although the plot's axis was labelled "Branch Length", and ``circular_tree_plot``
drew its Cartesian unit circle on a polar axes, which reads the same numbers as (angle, radius).

The trees are built from sequences whose p-distances are the exact small distance matrix of a known tree (every
branch mutates private sites), so neighbor joining must recover its distances. Every assertion reads the artists:
node offsets, label positions, edge segments and axis state. Real trees, real matplotlib Agg figures, no test
doubles.
"""

from __future__ import annotations

import itertools
import math
from io import StringIO
from typing import Any, Dict, Iterator, List, Tuple

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pytest

from metainformant.dna.phylogeny import (
    bootstrap_support,
    from_newick,
    neighbor_joining_tree,
    prune_tree,
    to_newick,
    upgma_tree,
)
from metainformant.visualization.genomics import trees

nx = pytest.importorskip("networkx")

matplotlib.use("Agg")

#: child -> (parent, number of private alignment sites mutated on the branch to it)
BRANCH_SITES: Dict[str, Tuple[str, int]] = {
    "n1": ("root", 4),
    "n2": ("root", 2),
    "Apis": ("n1", 3),
    "n3": ("n1", 2),
    "Bombus": ("n3", 2),
    "Camponotus": ("n3", 4),
    "Drosophila": ("n2", 3),
    "Eciton": ("n2", 5),
}
LEAVES = ["Apis", "Bombus", "Camponotus", "Drosophila", "Eciton"]
ALIGNMENT_LENGTH = 40


def _sequences() -> Dict[str, str]:
    """Sequences whose p-distance matrix is exactly the tree-metric matrix of BRANCH_SITES / ALIGNMENT_LENGTH."""
    private: Dict[str, List[int]] = {}
    next_site = 0
    for child, (_, count) in BRANCH_SITES.items():
        private[child] = list(range(next_site, next_site + count))
        next_site += count
    assert next_site <= ALIGNMENT_LENGTH

    def sequence(node: str) -> List[str]:
        if node == "root":
            return ["A"] * ALIGNMENT_LENGTH
        letters = sequence(BRANCH_SITES[node][0])
        for site in private[node]:
            letters[site] = "T"
        return letters

    return {leaf: "".join(sequence(leaf)) for leaf in LEAVES}


def _true_distance(a: str, b: str) -> float:
    def path(node: str) -> List[str]:
        out = [node]
        while node != "root":
            node = BRANCH_SITES[node][0]
            out.append(node)
        return out

    pa, pb = path(a), path(b)
    shared = set(pa) & set(pb)
    sites = sum(BRANCH_SITES[n][1] for n in pa if n not in shared) + sum(
        BRANCH_SITES[n][1] for n in pb if n not in shared
    )
    return sites / ALIGNMENT_LENGTH


@pytest.fixture(autouse=True)
def _close_figures() -> Iterator[None]:
    yield
    plt.close("all")


@pytest.fixture(scope="module")
def sequences() -> Dict[str, str]:
    return _sequences()


@pytest.fixture(scope="module")
def nj(sequences: Dict[str, str]) -> Dict[str, Any]:
    return neighbor_joining_tree(sequences)


@pytest.fixture(scope="module")
def upgma(sequences: Dict[str, str]) -> Dict[str, Any]:
    return upgma_tree(sequences)


# ---------------------------------------------------------------------------
# Independent readers of the dictionary tree and of the drawn artists
# ---------------------------------------------------------------------------


def _children(tree: Dict[str, Any], node: str) -> Dict[str, float]:
    entry = tree[node]
    return {} if entry is None else {child: length for child, length in entry.items() if child != "bootstrap"}


def _root(tree: Dict[str, Any]) -> str:
    kids = {child for node in tree for child in _children(tree, node)}
    (root,) = set(tree) - kids
    return root


def _root_to_tip(tree: Dict[str, Any]) -> Dict[str, float]:
    """Distance from the root to every leaf, summed along the dictionary's own branch lengths."""
    out: Dict[str, float] = {}
    stack = [(_root(tree), 0.0)]
    while stack:
        node, dist = stack.pop()
        kids = _children(tree, node)
        if not kids:
            out[node] = dist
        stack.extend((child, dist + max(length, 0.0)) for child, length in kids.items())
    return out


def _path_distances(tree: Dict[str, Any]) -> Dict[frozenset, float]:
    """Leaf-to-leaf path length for every pair of leaves."""
    parent = {child: node for node in tree for child in _children(tree, node)}
    blen = {child: length for node in tree for child, length in _children(tree, node).items()}

    def up(node: str) -> Dict[str, float]:
        out, total = {node: 0.0}, 0.0
        while node in parent:
            total += blen[node]
            node = parent[node]
            out[node] = total
        return out

    leaves = [n for n in tree if not _children(tree, n)]
    dist = {}
    for a, b in itertools.combinations(leaves, 2):
        ua, ub = up(a), up(b)
        meet = next(n for n in ub if n in ua)
        dist[frozenset((a, b))] = ua[meet] + ub[meet]
    return dist


def _node_offsets(ax: Any) -> np.ndarray:
    return np.asarray(ax.collections[0].get_offsets(), dtype=float)


def _key(point: Any) -> Tuple[float, float]:
    return (round(float(point[0]), 9), round(float(point[1]), 9))


def _drawn_tips(ax: Any) -> List[Tuple[float, float]]:
    """Positions of the drawn nodes that no edge starts from (the tips)."""
    starts = {_key(segment[0]) for segment in ax.collections[1].get_segments()}
    return [_key(point) for point in _node_offsets(ax) if _key(point) not in starts]


def _label_positions(ax: Any) -> Dict[str, Tuple[float, float]]:
    return {text.get_text(): _key(text.get_position()) for text in ax.texts}


# ---------------------------------------------------------------------------
# The builders' own output
# ---------------------------------------------------------------------------


def test_the_builders_return_the_flat_layout_the_plots_must_read(nj: Dict[str, Any], upgma: Dict[str, Any]) -> None:
    """Pin the premise: leaves are ``None`` keys, internal nodes ``{child: length}``, all nodes top-level."""
    for tree in (nj, upgma):
        assert [n for n, entry in tree.items() if entry is None] == LEAVES
        assert all(set(_children(tree, n)) <= set(tree) for n in tree)


def test_neighbor_joining_recovers_the_exact_distance_matrix(nj: Dict[str, Any]) -> None:
    """The fixture is a faithful small distance matrix: NJ's path lengths equal the true distances."""
    drawn = _path_distances(nj)
    for a, b in itertools.combinations(LEAVES, 2):
        assert drawn[frozenset((a, b))] == pytest.approx(_true_distance(a, b), abs=1e-9)


# ---------------------------------------------------------------------------
# plot_phylo_tree on neighbor-joining and UPGMA trees
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("builder", ["nj", "upgma"])
def test_plot_phylo_tree_draws_a_builder_tree_with_every_leaf_label_on_a_tip(
    builder: str, nj: Dict[str, Any], upgma: Dict[str, Any]
) -> None:
    tree = nj if builder == "nj" else upgma

    ax = trees.plot_phylo_tree(tree)

    assert len(_node_offsets(ax)) == len(tree)  # every node of the dictionary, none invented (no "bootstrap" node)
    labels = _label_positions(ax)
    assert set(labels) == set(tree)
    tips = _drawn_tips(ax)
    assert len(tips) == len(LEAVES) == 5
    assert sorted(labels[name] for name in LEAVES) == sorted(tips)  # each taxon label sits on its own tip


@pytest.mark.parametrize("builder", ["nj", "upgma"])
def test_plot_phylo_tree_places_each_tip_at_its_root_to_tip_branch_length(
    builder: str, nj: Dict[str, Any], upgma: Dict[str, Any]
) -> None:
    tree = nj if builder == "nj" else upgma

    ax = trees.plot_phylo_tree(tree)

    labels = _label_positions(ax)
    for leaf, distance in _root_to_tip(tree).items():
        assert labels[leaf][0] == pytest.approx(distance, abs=1e-8)
    assert ax.get_xlabel() == "Branch Length"
    assert ax.axison and ax.xaxis.get_visible()


def test_upgma_tips_all_sit_at_one_height(upgma: Dict[str, Any]) -> None:
    """UPGMA trees are ultrametric, so the drawn tips must line up on one vertical."""
    ax = trees.plot_phylo_tree(upgma)

    tip_x = {x for x, _ in _drawn_tips(ax)}
    assert len(tip_x) == 1
    assert next(iter(tip_x)) == pytest.approx(next(iter(_root_to_tip(upgma).values())), abs=1e-8)


@pytest.mark.parametrize("builder", ["nj", "upgma"])
def test_plot_phylo_tree_gives_every_internal_node_the_mean_row_of_its_children(
    builder: str, nj: Dict[str, Any], upgma: Dict[str, Any]
) -> None:
    tree = nj if builder == "nj" else upgma
    ax = trees.plot_phylo_tree(tree)
    labels = _label_positions(ax)

    assert len({labels[name][1] for name in LEAVES}) == len(LEAVES)  # one row per leaf
    for node in tree:
        kids = _children(tree, node)
        if kids:
            assert labels[node][1] == pytest.approx(np.mean([labels[child][1] for child in kids]), abs=1e-8)
            for child in kids:
                assert labels[child][0] >= labels[node][0] - 1e-8  # a child is never left of its parent


@pytest.mark.parametrize(
    "plot", ["plot_phylo_tree", "circular_tree_plot", "unrooted_tree_plot", "tree_annotation_plot"]
)
@pytest.mark.parametrize("builder", ["nj", "upgma"])
def test_every_tree_plot_draws_a_builder_tree(
    plot: str, builder: str, nj: Dict[str, Any], upgma: Dict[str, Any]
) -> None:
    """All of the single-tree plots used to raise ``AttributeError: 'NoneType' object has no attribute 'items'``."""
    tree = nj if builder == "nj" else upgma
    draw = getattr(trees, plot)

    ax = draw(tree, {}) if plot == "tree_annotation_plot" else draw(tree)

    expected_nodes = len(tree) - 1 if plot == "unrooted_tree_plot" else len(tree)  # an unrooted tree has no root
    assert len(_node_offsets(ax)) == expected_nodes
    assert set(LEAVES) <= {text.get_text() for text in ax.texts}


def test_tree_comparison_plot_draws_two_builder_trees(nj: Dict[str, Any], upgma: Dict[str, Any]) -> None:
    first = trees.tree_comparison_plot(nj, upgma)

    left, right = first.figure.axes
    assert len(_node_offsets(left)) == len(nj)
    assert len(_node_offsets(right)) == len(upgma)
    assert set(LEAVES) <= {text.get_text() for text in left.texts}
    assert set(LEAVES) <= {text.get_text() for text in right.texts}
    assert left.get_title() == "Tree 1" and right.get_title() == "Tree 2"


# ---------------------------------------------------------------------------
# Bootstrap-annotated trees: "bootstrap" is a support value, not a child
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("builder", ["nj", "upgma"])
def test_bootstrap_entries_are_support_values_not_nodes(
    builder: str, sequences: Dict[str, str], nj: Dict[str, Any], upgma: Dict[str, Any]
) -> None:
    tree = nj if builder == "nj" else upgma
    supported = bootstrap_support(tree, sequences, n_replicates=20, method=builder)
    internal = [n for n, entry in supported.items() if isinstance(entry, dict)]

    ax = trees.plot_phylo_tree(supported)

    assert len(_node_offsets(ax)) == len(tree)  # no node called "bootstrap", no edge into one
    assert "bootstrap" not in {text.get_text() for text in ax.texts}
    assert len(_drawn_tips(ax)) == len(LEAVES)
    shown = sorted(text.get_text() for text in ax.texts if text.get_text().isdigit())
    assert shown == sorted(str(supported[node]["bootstrap"]) for node in internal)


def test_show_support_false_hides_the_support_values(sequences: Dict[str, str], nj: Dict[str, Any]) -> None:
    supported = bootstrap_support(nj, sequences, n_replicates=10)

    ax = trees.plot_phylo_tree(supported, show_support=False)

    assert not [text for text in ax.texts if text.get_text().isdigit()]


# ---------------------------------------------------------------------------
# Other documented representations: Newick, nested dictionaries, Bio.Phylo, graphs
# ---------------------------------------------------------------------------


def test_plot_phylo_tree_draws_a_newick_string() -> None:
    ax = trees.plot_phylo_tree("((A:0.1,B:0.2):0.3,C:0.4);")

    labels = _label_positions(ax)
    assert {"A", "B", "C"} <= set(labels)
    assert len(_drawn_tips(ax)) == 3
    # the tips sit at their root-to-tip distances: A 0.3 + 0.1, B 0.3 + 0.2, C 0.4
    assert labels["A"][0] == pytest.approx(0.4)
    assert labels["B"][0] == pytest.approx(0.5)
    assert labels["C"][0] == pytest.approx(0.4)


def test_a_newick_tree_without_branch_lengths_is_drawn_by_depth() -> None:
    """``from_newick`` gives every edge length 0.0 when the string has none; that must not collapse the plot."""
    ax = trees.plot_phylo_tree("((A,B),C);")

    labels = _label_positions(ax)
    assert (labels["A"][0], labels["B"][0], labels["C"][0]) == (2.0, 2.0, 1.0)
    assert ax.get_xlabel() == "Depth (edges from root)"


def test_the_newick_written_for_a_builder_tree_draws_the_same_tips(nj: Dict[str, Any]) -> None:
    ax = trees.plot_phylo_tree(to_newick(nj))

    assert set(LEAVES) <= {text.get_text() for text in ax.texts}
    assert len(_drawn_tips(ax)) == len(LEAVES)


def test_a_pruned_builder_tree_still_draws(nj: Dict[str, Any]) -> None:
    pruned = prune_tree(nj, ["Apis", "Bombus", "Eciton"])

    ax = trees.plot_phylo_tree(pruned)

    assert {"Apis", "Bombus", "Eciton"} <= {text.get_text() for text in ax.texts}
    assert len(_drawn_tips(ax)) == 3


def test_a_parsed_newick_dictionary_draws(nj: Dict[str, Any]) -> None:
    ax = trees.plot_phylo_tree(from_newick(to_newick(nj)))

    assert len(_drawn_tips(ax)) == len(LEAVES)


def test_nested_dictionary_leaves_may_carry_branch_lengths() -> None:
    ax = trees.plot_phylo_tree({"root": {"A": 0.5, "B": 2.0}})

    labels = _label_positions(ax)
    assert (labels["A"][0], labels["B"][0]) == (0.5, 2.0)
    assert ax.get_xlabel() == "Branch Length"


def test_a_nested_dictionary_with_an_unmeasured_branch_is_drawn_by_depth() -> None:
    """The branch to an internal nested node has no length, so no branch is drawn to scale (no mixed scales)."""
    ax = trees.plot_phylo_tree({"root": {"A": {"A1": 0.5, "A2": 1.5}, "B": 2.0}})

    labels = _label_positions(ax)
    assert (labels["A1"][0], labels["A2"][0], labels["B"][0]) == (2.0, 2.0, 1.0)
    assert ax.get_xlabel() == "Depth (edges from root)"
    assert len(_drawn_tips(ax)) == 3


def test_nested_dictionary_leaf_lists_are_children() -> None:
    ax = trees.plot_phylo_tree({"root": {"A": ["a1", "a2"], "B": None}})

    assert {"a1", "a2", "B"} <= {text.get_text() for text in ax.texts}
    assert len(_drawn_tips(ax)) == 3


def test_bio_phylo_tree_keeps_a_zero_length_branch_at_zero() -> None:
    """The converter replaced a branch length of 0.0 with 1.0, so a zero-length tip was drawn a unit away."""
    phylo = pytest.importorskip("Bio.Phylo.TreeConstruction")
    matrix = phylo.DistanceMatrix(["A", "B", "C"], [[0], [0.1, 0], [0.2, 0.3, 0]])
    bio_tree = phylo.DistanceTreeConstructor().nj(matrix)  # A hangs off the root on a branch of length 0.0

    graph = trees._convert_tree_to_networkx(bio_tree)
    assert {(u, v): d["weight"] for u, v, d in graph.edges(data=True)}[("Inner1", "A")] == 0.0

    ax = trees.plot_phylo_tree(bio_tree)
    labels = _label_positions(ax)
    assert (labels["A"][0], labels["B"][0], labels["C"][0]) == pytest.approx((0.0, 0.1, 0.2))


def test_bio_phylo_clade_and_duplicate_names() -> None:
    bio = pytest.importorskip("Bio.Phylo")
    tree = bio.read(StringIO("((A:1,B:1):1,(C:1,D:1):1);"), "newick")
    ax = trees.plot_phylo_tree(tree.root)  # a clade, not a tree
    assert {"A", "B", "C", "D"} <= {text.get_text() for text in ax.texts}
    assert len(_drawn_tips(ax)) == 4

    duplicate = bio.read(StringIO("((A:1,A:1):1,C:1);"), "newick")
    with pytest.raises(ValueError, match="more than one clade named"):
        trees._convert_tree_to_networkx(duplicate)


def test_networkx_graph_edge_weights_are_branch_lengths() -> None:
    graph = nx.DiGraph()
    graph.add_edge("root", "A", weight=2.0)
    graph.add_edge("root", "B", weight=5.0)

    ax = trees.plot_phylo_tree(graph)

    labels = _label_positions(ax)
    assert (labels["A"][0], labels["B"][0]) == (2.0, 5.0)


def test_negative_branch_lengths_are_drawn_as_zero_length() -> None:
    graph = nx.DiGraph()
    graph.add_edge("root", "A", weight=-0.2)
    graph.add_edge("root", "B", weight=0.5)

    ax = trees.plot_phylo_tree(graph)

    assert _label_positions(ax)["A"][0] == 0.0


@pytest.mark.parametrize("bad", [42, 3.5, None, "not a tree", ["A", "B"], ""])
def test_unsupported_representations_are_refused(bad: Any) -> None:
    with pytest.raises(ValueError, match="Unsupported tree format"):
        trees._convert_tree_to_networkx(bad)


def test_an_undirected_graph_is_refused_by_the_rooted_plot_but_drawn_unrooted() -> None:
    graph = nx.Graph([("a", "b"), ("b", "c"), ("b", "d")])

    with pytest.raises(ValueError, match="directed"):
        trees.plot_phylo_tree(graph)
    assert len(_node_offsets(trees.unrooted_tree_plot(graph))) == 4


def test_a_cyclic_graph_is_refused() -> None:
    cyclic = nx.DiGraph([("root", "a"), ("a", "b"), ("b", "a")])

    with pytest.raises(nx.NetworkXUnfeasible):
        trees.plot_phylo_tree(cyclic)


def test_dictionary_node_order_follows_the_tree_not_the_builder() -> None:
    """Leaves come first in a builder dictionary; the graph must list them in depth-first (tree) order."""
    tree = {
        "a": None,
        "b": None,
        "c": None,
        "d": None,
        "x": {"a": 1.0, "c": 1.0},
        "y": {"b": 1.0, "d": 1.0},
        "r": {"x": 1.0, "y": 1.0},
    }

    graph = trees._convert_tree_to_networkx(tree)

    assert list(graph.nodes) == ["r", "x", "a", "c", "y", "b", "d"]


# ---------------------------------------------------------------------------
# circular_tree_plot geometry
# ---------------------------------------------------------------------------


def _display_polar(ax: Any) -> Dict[str, Tuple[float, float]]:
    """(angle in degrees, radius in pixels) of every labelled node, measured from the axes' polar centre."""
    ax.figure.canvas.draw()
    centre = ax.transData.transform((0.0, 0.0))
    out = {}
    for text in ax.texts:
        vector = ax.transData.transform(text.get_position()) - centre
        out[text.get_text()] = (math.degrees(math.atan2(vector[1], vector[0])) % 360.0, float(math.hypot(*vector)))
    return out


def test_circular_tree_plot_puts_the_leaves_on_a_circle_at_even_angles() -> None:
    """Cartesian (cos, sin) positions on a polar axes were read as (angle, radius): the leaves were not on a circle."""
    graph = nx.DiGraph([("root", "A"), ("root", "B"), ("root", "C"), ("root", "D")])

    ax = trees.circular_tree_plot(graph)

    assert ax.name == "polar"
    placed = _display_polar(ax)
    radii = [placed[leaf][1] for leaf in "ABCD"]
    assert radii == pytest.approx([radii[0]] * 4, rel=1e-6)
    angles = [placed[leaf][0] for leaf in "ABCD"]
    assert np.diff(angles) == pytest.approx([90.0, 90.0, 90.0], abs=1e-6)
    assert placed["root"][1] == pytest.approx(0.0, abs=1e-6)  # the root of a symmetric tree is at the centre


def test_circular_tree_plot_of_a_builder_tree_has_every_leaf_on_the_outer_circle(nj: Dict[str, Any]) -> None:
    ax = trees.circular_tree_plot(nj)

    placed = _display_polar(ax)
    leaf_radii = [placed[leaf][1] for leaf in LEAVES]
    assert leaf_radii == pytest.approx([leaf_radii[0]] * len(LEAVES), rel=1e-6)
    assert np.diff(sorted(placed[leaf][0] for leaf in LEAVES)) == pytest.approx([72.0] * 4, abs=1e-6)
    internal_radii = [r for name, (_, r) in placed.items() if name not in LEAVES]
    assert max(internal_radii) < leaf_radii[0]  # internal nodes are pulled inward


def test_circular_tree_plot_accepts_a_cartesian_axes() -> None:
    """A non-polar axes made ``set_rlabel_position`` raise ``AttributeError``."""
    graph = nx.DiGraph([("root", "A"), ("root", "B"), ("root", "C"), ("root", "D")])
    _, cartesian = plt.subplots()

    ax = trees.circular_tree_plot(graph, ax=cartesian)

    assert ax is cartesian
    labels = _label_positions(ax)
    for leaf, angle in zip("ABCD", (0.0, 90.0, 180.0, 270.0)):
        x, y = labels[leaf]
        assert math.hypot(x, y) == pytest.approx(1.0)
        assert math.degrees(math.atan2(y, x)) % 360.0 == pytest.approx(angle)
    assert ax.get_aspect() == 1.0


# ---------------------------------------------------------------------------
# unrooted_tree_plot
# ---------------------------------------------------------------------------


def test_unrooted_tree_plot_removes_the_root_and_joins_its_branches(nj: Dict[str, Any]) -> None:
    root = _root(nj)
    (first, first_length), (second, second_length) = _children(nj, root).items()

    graph = trees._convert_tree_to_networkx(nj, directed=False)

    assert not graph.is_directed()
    assert root not in graph
    assert graph.number_of_nodes() == len(nj) - 1
    assert graph.number_of_edges() == graph.number_of_nodes() - 1  # still a tree
    assert graph[first][second]["weight"] == pytest.approx(first_length + second_length)
    ax = trees.unrooted_tree_plot(nj)
    assert root not in {text.get_text() for text in ax.texts}


def test_unrooted_layout_is_reproducible_when_other_layout_arguments_are_given() -> None:
    """``layout_kwargs`` replaced the default ``{"seed": 42}``, so any other layout argument left the seed unset."""
    graph = nx.balanced_tree(2, 4)

    first = _node_offsets(trees.unrooted_tree_plot(graph, layout_kwargs={"iterations": 25}))
    second = _node_offsets(trees.unrooted_tree_plot(graph, layout_kwargs={"iterations": 25}))

    np.testing.assert_array_equal(first, second)
    other = _node_offsets(trees.unrooted_tree_plot(graph, layout_kwargs={"iterations": 25, "seed": 7}))
    assert not np.array_equal(first, other)  # an explicit seed still wins
    default = _node_offsets(trees.unrooted_tree_plot(graph))
    np.testing.assert_array_equal(default, _node_offsets(trees.unrooted_tree_plot(graph, layout_kwargs={"seed": 42})))


def test_unrooted_plot_of_a_directed_graph_has_no_arrows() -> None:
    ax = trees.unrooted_tree_plot(nx.DiGraph([("root", "A"), ("root", "B"), ("A", "A1")]))

    assert not ax.patches  # arrow heads would be FancyArrowPatch artists


# ---------------------------------------------------------------------------
# tree_comparison_plot and tree_annotation_plot
# ---------------------------------------------------------------------------


def test_tree_comparison_plot_draws_on_the_axes_it_is_given() -> None:
    """A supplied ``ax`` was ignored: two subplots were stacked on top of it and ``ax`` stayed empty."""
    figure, host = plt.subplots()
    first = nx.DiGraph([("r1", "a"), ("r1", "b")])
    second = nx.DiGraph([("r2", "c"), ("r2", "d"), ("d", "e")])

    panel = trees.tree_comparison_plot(first, second, ax=host)

    assert figure.axes == [host]  # no new subplots were added to the figure
    assert len(host.child_axes) == 2 and panel is host.child_axes[0]
    assert not host.axison  # the host is only a frame for the two panels
    assert len(_node_offsets(host.child_axes[0])) == 3
    assert len(_node_offsets(host.child_axes[1])) == 4


def test_tree_comparison_plot_forwards_drawing_keywords_to_both_panels() -> None:
    first = nx.DiGraph([("r1", "a"), ("r1", "b")])

    left = trees.tree_comparison_plot(first, first, node_size=33, with_labels=False)

    for panel in left.figure.axes:
        assert set(panel.collections[0].get_sizes()) == {33.0}
        assert not panel.texts


def test_tree_annotation_plot_colours_and_labels_annotated_nodes(nj: Dict[str, Any]) -> None:
    annotations = {"Apis": {"color": "#D32F2F", "label": "honeybee"}, "Eciton": {"label": "army ant"}}

    ax = trees.tree_annotation_plot(nj, annotations, node_color="#CCCCCC")

    texts = _label_positions(ax)
    assert {"honeybee", "army ant"} <= set(texts)
    names = list(trees._convert_tree_to_networkx(nj).nodes)
    colours = {name: tuple(np.round(rgba, 3)) for name, rgba in zip(names, ax.collections[0].get_facecolor())}
    assert colours["Apis"] == (0.827, 0.184, 0.184, 1.0)  # the annotated colour
    assert colours["Eciton"] == colours["Bombus"] == (0.8, 0.8, 0.8, 1.0)  # the node_color keyword, not a collision
    # the label sits just right of its node, offset by a fraction of the tree's width (0.1 data units used to be fixed)
    offset = texts["honeybee"][0] - texts["Apis"][0]
    width = max(x for x, _ in texts.values()) - min(x for x, _ in texts.values())
    assert 0 < offset < 0.1 * width


def test_tree_annotation_plot_without_annotations_is_a_plain_tree(nj: Dict[str, Any]) -> None:
    ax = trees.tree_annotation_plot(nj)

    assert len(_node_offsets(ax)) == len(nj)
