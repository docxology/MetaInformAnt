"""Regression tests: ``louvain_communities`` must not seed the global random generators.

With a ``random_state`` the function ran ``random.seed(random_state)`` and ``np.random.seed(random_state)`` before
calling python-louvain: the caller's global ``random`` / ``numpy.random`` streams were silently reset, and with NumPy
missing the second call was a ``NameError`` on ``np``. python-louvain has a ``random_state`` argument of its own that
builds a local generator, so the seed is now handed to it and nothing global is touched. (python-louvain treats
``randomize`` as a deprecated alternative to ``random_state``, so with a seed only the seed is passed on.)

python-louvain is an optional dependency and is not installed in every environment, so there are two checks: a static
one that parses the source and finds no call to a ``.seed`` function in the module, and an end-to-end one that runs
only where python-louvain is available (it is skipped otherwise). No test doubles.
"""

from __future__ import annotations

import ast
import inspect
import random
import textwrap

import networkx as nx
import numpy as np
import pytest

from metainformant.networks.analysis import community as community_module
from metainformant.networks.analysis.community import detect_communities, louvain_communities


def _seed_calls(source: str) -> list[str]:
    """Every ``<something>.seed(...)`` call in the source, as written."""
    tree = ast.parse(textwrap.dedent(source))
    return [
        ast.unparse(node.func)
        for node in ast.walk(tree)
        if isinstance(node, ast.Call) and isinstance(node.func, ast.Attribute) and node.func.attr == "seed"
    ]


def test_louvain_communities_does_not_call_a_global_seed_function() -> None:
    assert _seed_calls(inspect.getsource(louvain_communities)) == []


def test_the_community_module_seeds_no_global_generator_anywhere() -> None:
    assert _seed_calls(inspect.getsource(community_module)) == []


def _two_cliques() -> nx.Graph:
    graph = nx.Graph()
    first, second = [f"a{i}" for i in range(6)], [f"b{i}" for i in range(6)]
    graph.add_edges_from((u, v) for i, u in enumerate(first) for v in first[i + 1 :])
    graph.add_edges_from((u, v) for i, u in enumerate(second) for v in second[i + 1 :])
    graph.add_edge("a0", "b0")
    return graph


def _louvain_installed() -> bool:
    try:
        import community as community_louvain  # type: ignore[import-not-found]
    except ImportError:
        return False
    return hasattr(community_louvain, "best_partition")


needs_python_louvain = pytest.mark.skipif(not _louvain_installed(), reason="python-louvain is not installed")


@needs_python_louvain
def test_a_seeded_run_is_reproducible_and_leaves_the_global_generators_alone() -> None:
    graph = _two_cliques()
    random.seed(20261027)
    np.random.seed(20261027)
    python_state, numpy_state = random.getstate(), np.random.get_state()

    first = louvain_communities(graph, random_state=3)
    second = louvain_communities(graph, random_state=3)

    assert first == second
    assert {frozenset(community) for community in first} == {
        frozenset(f"a{i}" for i in range(6)),
        frozenset(f"b{i}" for i in range(6)),
    }
    assert random.getstate() == python_state
    after = np.random.get_state()
    assert after[0] == numpy_state[0] and after[2:] == numpy_state[2:] and np.array_equal(after[1], numpy_state[1])


@needs_python_louvain
def test_detect_communities_passes_the_seed_through_to_louvain() -> None:
    graph = _two_cliques()

    assert detect_communities(graph, method="louvain", random_state=5) == detect_communities(
        graph, method="louvain", random_state=5
    )


@needs_python_louvain
def test_without_randomisation_the_result_does_not_depend_on_the_seed() -> None:
    graph = _two_cliques()

    assert louvain_communities(graph, randomize=False, random_state=1) == louvain_communities(
        graph, randomize=False, random_state=2
    )
