"""``bootstrap_support`` must leave the caller's tree alone.

``bootstrap_support`` annotated ``tree.copy()``, a shallow copy that shares every internal-node dictionary with the
caller's tree, so the support values were also written into the input. ``basic_tree_stats`` then raised
``KeyError: 'bootstrap'`` on the input because ``_calculate_tree_height`` treated the stored support value as a child
node; it raised the same error for the annotated tree that ``bootstrap_support`` returns.

Real trees built from real sequences; no test doubles (real-implementation policy).
"""

from __future__ import annotations

import copy
from typing import Any, Dict

import pytest

from metainformant.dna.phylogeny.tree_analysis import (
    _find_root,
    basic_tree_stats,
    bootstrap_support,
    to_newick,
)
from metainformant.dna.phylogeny.tree_construction import neighbor_joining_tree, upgma_tree

SEQUENCES = {
    "A": "AAAAAAAACCCCCCCC",
    "B": "AAAAAAAACCCCCCCC",
    "C": "GGGGGGGGTTTTTTTT",
    "D": "GGGGGGGGTTTTTTTT",
}


def _internal_nodes(tree: Dict[str, Any]) -> list[str]:
    return [node for node, children in tree.items() if isinstance(children, dict)]


def _has_support_values(tree: Dict[str, Any]) -> bool:
    return any("bootstrap" in tree[node] for node in _internal_nodes(tree))


@pytest.fixture
def tree() -> Dict[str, Any]:
    return neighbor_joining_tree(SEQUENCES)


@pytest.mark.parametrize("method", ["nj", "upgma"])
def test_input_tree_is_unchanged(method: str) -> None:
    tree = neighbor_joining_tree(SEQUENCES) if method == "nj" else upgma_tree(SEQUENCES)
    before = copy.deepcopy(tree)

    supported = bootstrap_support(tree, SEQUENCES, n_replicates=10, method=method)

    assert _has_support_values(supported)  # the result really was annotated
    assert tree == before
    assert not _has_support_values(tree)


def test_input_tree_still_works_with_basic_tree_stats(tree: Dict[str, Any]) -> None:
    before = basic_tree_stats(tree)

    bootstrap_support(tree, SEQUENCES, n_replicates=10)

    assert basic_tree_stats(tree) == before  # raised KeyError: 'bootstrap' when the input had been annotated
    assert before["leaves"] == len(SEQUENCES)


def test_annotated_tree_works_with_basic_tree_stats(tree: Dict[str, Any]) -> None:
    supported = bootstrap_support(tree, SEQUENCES, n_replicates=10)

    # support values are annotations, not nodes: they change neither the counts nor the height
    assert basic_tree_stats(supported) == basic_tree_stats(tree)


def test_returned_tree_does_not_share_nodes_with_the_input(tree: Dict[str, Any]) -> None:
    before = copy.deepcopy(tree)
    supported = bootstrap_support(tree, SEQUENCES, n_replicates=10)

    for node in _internal_nodes(tree):
        assert supported[node] is not tree[node]

    # changing the result afterwards must not reach back into the input
    root = _find_root(supported)
    supported[root]["bootstrap"] = -1
    supported[root]["extra"] = 1.0
    assert tree == before


def test_repeated_calls_agree_and_each_leaves_the_input_unchanged(tree: Dict[str, Any]) -> None:
    before = copy.deepcopy(tree)

    first = bootstrap_support(tree, SEQUENCES, n_replicates=12)
    assert tree == before
    second = bootstrap_support(tree, SEQUENCES, n_replicates=12)

    assert first == second
    assert tree == before


def test_already_annotated_tree_is_recomputed_not_accumulated(tree: Dict[str, Any]) -> None:
    supported = bootstrap_support(tree, SEQUENCES, n_replicates=10)
    annotated_before = copy.deepcopy(supported)

    again = bootstrap_support(supported, SEQUENCES, n_replicates=10)

    assert again == annotated_before
    assert supported == annotated_before


@pytest.mark.parametrize("sequences", [{}, {"A": "", "B": "", "C": "", "D": ""}])
def test_early_returns_are_independent_copies_too(tree: Dict[str, Any], sequences: Dict[str, str]) -> None:
    before = copy.deepcopy(tree)

    result = bootstrap_support(tree, sequences, n_replicates=10)

    assert result == before
    root = _find_root(result)
    assert result[root] is not tree[root]
    result[root]["bootstrap"] = 100
    assert tree == before
    assert to_newick(tree) == to_newick(before)
