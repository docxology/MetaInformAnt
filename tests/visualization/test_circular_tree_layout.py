"""Regression tests: the circular tree layout must position every node of a multi-level tree.

``_circular_tree_layout`` put the leaves on the unit circle and then walked the internal nodes in the
order they were added to the graph, placing each one at the centroid of the children that already had
a position. A node visited before its internal children found no positioned child and was skipped, so
``nx.draw`` failed with ``NetworkXError: Node root has no position`` for every tree with more than one
level of internal nodes whose edges were added parent-first (root -> internal -> internal -> leaf).
A node with a mix of leaf children and not-yet-positioned internal children did not fail at all: it was
silently placed from its leaf children alone.

All trees are real ``networkx`` graphs, drawn onto real matplotlib Agg axes; no test doubles.
"""

from __future__ import annotations

import math
from typing import Any, Iterator

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pytest

from metainformant.visualization.genomics import trees

nx = pytest.importorskip("networkx")

matplotlib.use("Agg")

# Edges are listed parent-first, so internal nodes enter the graph before their internal children.
THREE_LEVEL = [("root", "A"), ("root", "B"), ("A", "A1"), ("A", "A2"), ("B", "B1"), ("B", "B2")]
FOUR_LEVEL_UNEVEN = [
    ("root", "A"),
    ("root", "B"),
    ("A", "A1"),
    ("A", "C"),
    ("C", "C1"),
    ("C", "C2"),
    ("C", "C3"),
    ("B", "B1"),
]
UNARY_CHAIN = [("root", "A"), ("A", "B"), ("B", "leaf")]

TREES = [
    pytest.param(THREE_LEVEL, id="three-level"),
    pytest.param(FOUR_LEVEL_UNEVEN, id="four-level-uneven"),
    pytest.param(UNARY_CHAIN, id="unary-chain"),
]


@pytest.fixture(autouse=True)
def _close_figures() -> Iterator[None]:
    yield
    plt.close("all")


def _leaves(graph: Any) -> list[Any]:
    return [node for node in graph.nodes if graph.out_degree(node) == 0]


@pytest.mark.parametrize("edges", TREES)
def test_circular_tree_plot_draws_every_node_of_a_multi_level_tree(edges: list[tuple[str, str]]) -> None:
    graph = nx.DiGraph(edges)

    ax = trees.circular_tree_plot(graph)

    drawn = ax.collections[0].get_offsets()
    assert drawn.shape == (graph.number_of_nodes(), 2)
    assert np.isfinite(np.asarray(drawn)).all()
    assert {text.get_text() for text in ax.texts} == set(graph.nodes)


def test_circular_tree_plot_accepts_a_nested_dict_three_levels_deep() -> None:
    tree = {"root": {"A": {"A1": None, "A2": None}, "B": {"B1": None, "B2": None}}}

    ax = trees.circular_tree_plot(tree)

    assert ax.collections[0].get_offsets().shape == (7, 2)


@pytest.mark.parametrize("edges", TREES)
def test_layout_places_leaves_on_the_unit_circle_and_internal_nodes_at_the_inward_child_centroid(
    edges: list[tuple[str, str]],
) -> None:
    graph = nx.DiGraph(edges)

    pos = trees._circular_tree_layout(graph)

    assert set(pos) == set(graph.nodes)
    leaves = _leaves(graph)
    for index, leaf in enumerate(leaves):
        angle = 2 * math.pi * index / len(leaves)
        assert pos[leaf] == pytest.approx((math.cos(angle), math.sin(angle)))
    for node in set(graph.nodes) - set(leaves):
        children = [pos[child] for child in graph.successors(node)]
        centroid_x = sum(x for x, _ in children) / len(children)
        centroid_y = sum(y for _, y in children) / len(children)
        assert pos[node] == pytest.approx((0.7 * centroid_x, 0.7 * centroid_y))
        assert math.hypot(*pos[node]) < 1.0


def test_layout_does_not_depend_on_the_order_internal_nodes_were_added() -> None:
    """Parent-first and children-first insertion of one tree must give one layout."""
    leaves = ["A1", "A2", "B1", "B2"]
    parent_first = nx.DiGraph()
    parent_first.add_nodes_from(leaves)  # fixes the leaf order so only the internal-node order differs
    parent_first.add_edges_from([("root", "A"), ("root", "B"), ("A", "A1"), ("A", "A2"), ("B", "B1"), ("B", "B2")])
    children_first = nx.DiGraph()
    children_first.add_nodes_from(leaves)
    children_first.add_edges_from([("A", "A1"), ("A", "A2"), ("B", "B1"), ("B", "B2"), ("root", "A"), ("root", "B")])

    first = trees._circular_tree_layout(parent_first)
    second = trees._circular_tree_layout(children_first)

    assert set(first) == set(second) == {"root", "A", "B", *leaves}
    for node, position in first.items():
        assert second[node] == pytest.approx(position)


def test_node_with_leaf_and_internal_children_is_placed_from_all_of_them() -> None:
    """The old layout visited ``root`` before ``A`` and placed it from ``L1`` alone, at x = 0.7."""
    graph = nx.DiGraph([("root", "L1"), ("root", "A"), ("A", "L2"), ("A", "L3")])

    pos = trees._circular_tree_layout(graph)

    # Leaves L1, L2, L3 sit at 0, 120 and 240 degrees: (1, 0), (-0.5, 0.866) and (-0.5, -0.866).
    # A = 0.7 * centroid(L2, L3) = 0.7 * (-0.5, 0); root = 0.7 * centroid(L1, A) = 0.7 * (0.325, 0).
    assert pos["A"] == pytest.approx((-0.35, 0.0), abs=1e-4)
    assert pos["root"] == pytest.approx((0.2275, 0.0), abs=1e-4)


def test_layout_of_degenerate_graphs() -> None:
    single = nx.DiGraph()
    single.add_node("solo")

    assert trees._circular_tree_layout(nx.DiGraph()) == {}
    assert trees._circular_tree_layout(single)["solo"] == pytest.approx((1.0, 0.0))


def test_graph_with_a_cycle_is_refused_instead_of_laid_out_partially() -> None:
    cyclic = nx.DiGraph([("a", "b"), ("b", "a"), ("a", "leaf")])

    with pytest.raises(nx.NetworkXUnfeasible):
        trees._circular_tree_layout(cyclic)
