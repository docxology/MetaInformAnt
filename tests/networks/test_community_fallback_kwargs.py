"""Regression tests for the Louvain fallback of ``detect_communities(method="leiden")``.

When ``leidenalg`` / ``igraph`` cannot be imported, ``detect_communities`` falls back from Leiden to Louvain
(or to NetworkX greedy modularity when python-louvain is missing too). The keyword arguments of the call were handed to
the fallback unchanged, so a Leiden-only option such as ``n_iterations`` reached
``community_louvain.best_partition``, which does not accept it and raised ``TypeError`` instead of producing the
communities the caller asked for.

The fallback now keeps only the options Louvain understands. The signature-level check does not need python-louvain;
the end-to-end check runs only where python-louvain is installed and leidenalg is not (the combination that triggers the
fallback), and is skipped otherwise. The partition oracle is the known clique structure of the test graph.
No test doubles.
"""

from __future__ import annotations

import importlib.util

import networkx as nx
import pytest

from metainformant.networks.analysis import community as community_module
from metainformant.networks.analysis.community import detect_communities

CLIQUE_A = [f"a{i}" for i in range(6)]
CLIQUE_B = [f"b{i}" for i in range(6)]


def _two_cliques() -> nx.Graph:
    graph = nx.Graph()
    graph.add_edges_from((u, v) for i, u in enumerate(CLIQUE_A) for v in CLIQUE_A[i + 1 :])
    graph.add_edges_from((u, v) for i, u in enumerate(CLIQUE_B) for v in CLIQUE_B[i + 1 :])
    graph.add_edge("a0", "b0")  # one bridge
    return graph


def _groups(partition: dict[str, int]) -> set[frozenset[str]]:
    by_community: dict[int, set[str]] = {}
    for node, community_id in partition.items():
        by_community.setdefault(community_id, set()).add(node)
    return {frozenset(members) for members in by_community.values()}


def test_louvain_options_drop_what_leiden_alone_understands() -> None:
    kwargs = {
        "n_iterations": 5,  # Leiden only
        "resolution": 1.5,
        "random_state": 7,
        "weights": "weight",  # a leidenalg.find_partition option
        "initial_membership": [0, 1],  # a leidenalg.find_partition option
    }

    assert community_module._louvain_kwargs(kwargs) == {"resolution": 1.5, "random_state": 7}


def test_louvain_options_keep_every_option_louvain_accepts() -> None:
    kwargs = {"resolution": 0.8, "randomize": False, "random_state": 3, "weight": "w", "partition": {"a": 0}}

    assert community_module._louvain_kwargs(kwargs) == kwargs


def test_leiden_call_with_leiden_options_still_returns_the_cliques() -> None:
    # Works through whichever backend is installed: Leiden, else Louvain, else greedy modularity.
    partition = detect_communities(_two_cliques(), method="leiden", n_iterations=3, resolution=1.0, random_state=5)

    assert _groups(partition) == {frozenset(CLIQUE_A), frozenset(CLIQUE_B)}


@pytest.mark.skipif(
    importlib.util.find_spec("community") is None or importlib.util.find_spec("leidenalg") is not None,
    reason="needs python-louvain installed and leidenalg absent (the Leiden -> Louvain fallback)",
)
def test_leiden_to_louvain_fallback_does_not_raise_type_error() -> None:
    partition = detect_communities(_two_cliques(), method="leiden", n_iterations=3, random_state=5)

    assert _groups(partition) == {frozenset(CLIQUE_A), frozenset(CLIQUE_B)}
