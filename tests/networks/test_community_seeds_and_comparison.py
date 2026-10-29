"""Community detection keeps the method and the seed it was given, and compares partitions without made-up numbers.

* ``community_stability`` handed ``random_state`` to every method. NetworkX's asynchronous label propagation, fluid
  communities, label propagation and Girvan-Newman do not accept that keyword, so each run raised ``TypeError``,
  ``detect_communities`` caught it and silently ran greedy modularity instead: the "stability" of asyn_lpa, fluid,
  label_propagation and girvan_newman was the (perfect) stability of a deterministic algorithm that was never asked
  for. Stochastic methods now take ``random_state`` (NetworkX's ``seed``); deterministic ones accept and ignore it.
* ``NetworkWorkflow`` never read ``CommunityDetectionConfig.random_state``, and sent ``n_communities`` as ``k`` to
  every method and ``resolution`` to methods without one, which ended in the same silent greedy fallback.
* ``compare_communities`` fell back, when scikit-learn raised, to the fraction of nodes with equal labels reported as
  both NMI and ARI, which scores a relabelled copy of a partition 0.0. It now computes both exactly.
* ``modularity`` swallowed NetworkX's ``NotAPartition`` and returned an unweighted edge-count formula for any list of
  node sets.

Oracles: NetworkX called directly with the same seeds, and scikit-learn's NMI / ARI. The no-scikit-learn case runs the
real module in a subprocess in which scikit-learn cannot be imported. No test doubles.
"""

from __future__ import annotations

import json
import os
import subprocess
import sys
from pathlib import Path

import networkx as nx
import numpy as np
import pytest
from networkx.algorithms import community as nx_community

import metainformant
from metainformant.networks.analysis.community import (
    asyn_lpa_communities,
    community_stability,
    compare_communities,
    detect_communities,
    fluid_communities,
    girvan_newman_communities,
    greedy_modularity_communities,
    label_propagation_communities,
    modularity,
)
from metainformant.networks.config.config import CommunityDetectionConfig, NetworkWorkflowConfig
from metainformant.networks.workflow.workflow import NetworkWorkflow


def _graph() -> nx.Graph:
    """A connected 30-node random graph with no clear communities: stochastic methods differ between seeds."""
    graph = nx.gnp_random_graph(30, 0.12, seed=1)
    assert nx.is_connected(graph)
    return nx.relabel_nodes(graph, {node: f"n{node}" for node in graph})


def _groups(partition: dict[str, int]) -> set[frozenset[str]]:
    by_community: dict[int, set[str]] = {}
    for node, community_id in partition.items():
        by_community.setdefault(community_id, set()).add(node)
    return {frozenset(members) for members in by_community.values()}


def _as_groups(communities) -> set[frozenset[str]]:
    return {frozenset(members) for members in communities}


def _subprocess_env() -> dict[str, str]:
    package_root = str(Path(metainformant.__file__).parent.parent)
    existing = os.environ.get("PYTHONPATH")
    return {**os.environ, "PYTHONPATH": os.pathsep.join([package_root] + ([existing] if existing else []))}


class TestSeededMethods:
    def test_the_test_graph_separates_the_methods(self) -> None:
        graph = _graph()
        greedy = _as_groups(greedy_modularity_communities(graph))

        assert _as_groups(label_propagation_communities(graph)) != greedy
        assert _as_groups(girvan_newman_communities(graph)) != greedy
        assert len({frozenset(_as_groups(asyn_lpa_communities(graph, random_state=s))) for s in range(5, 9)}) > 1

    @pytest.mark.parametrize("seed", [5, 6, 7, 8])
    def test_asyn_lpa_takes_the_seed_as_networkx_does(self, seed: int) -> None:
        graph = _graph()

        got = _as_groups(asyn_lpa_communities(graph, random_state=seed))

        assert got == _as_groups(nx_community.asyn_lpa_communities(graph.copy(), seed=seed))
        assert got == _as_groups(detect_communities_lists(graph, "asyn_lpa", random_state=seed))

    @pytest.mark.parametrize("seed", [5, 6, 7, 8])
    def test_fluid_communities_take_the_seed_as_networkx_does(self, seed: int) -> None:
        graph = _graph()

        got = _as_groups(fluid_communities(graph, 3, random_state=seed))

        assert got == _as_groups(nx_community.asyn_fluidc(graph.copy(), 3, seed=seed))

    def test_the_same_seed_gives_the_same_partition_and_leaves_global_random_state_alone(self) -> None:
        import random

        graph = _graph()
        before = (random.getstate(), np.random.get_state()[1].tobytes())

        first = detect_communities(graph, "asyn_lpa", random_state=11)
        second = detect_communities(graph, "asyn_lpa", random_state=11)

        assert first == second
        assert (random.getstate(), np.random.get_state()[1].tobytes()) == before

    @pytest.mark.parametrize("method", ["label_propagation", "girvan_newman", "greedy"])
    def test_deterministic_methods_accept_a_seed_and_still_run_that_method(self, method: str) -> None:
        graph = _graph()
        direct = {
            "label_propagation": label_propagation_communities,
            "girvan_newman": girvan_newman_communities,
            "greedy": greedy_modularity_communities,
        }[method](graph)

        seeded = detect_communities(graph, method, random_state=3)

        assert _groups(seeded) == _as_groups(direct)
        if method != "greedy":
            # the old TypeError -> greedy fallback returned the greedy partition here
            assert _groups(seeded) != _as_groups(greedy_modularity_communities(graph))


def detect_communities_lists(graph: nx.Graph, method: str, **kwargs) -> list[list[str]]:
    """The communities of ``detect_communities`` as lists of nodes."""
    partition = detect_communities(graph, method, **kwargs)
    grouped: dict[int, list[str]] = {}
    for node, community_id in partition.items():
        grouped.setdefault(community_id, []).append(node)
    return list(grouped.values())


class TestStability:
    def test_stability_of_a_stochastic_method_is_that_methods_stability(self) -> None:
        graph = _graph()
        seed, runs = 5, 4

        result = community_stability(graph, method="asyn_lpa", n_runs=runs, seed=seed)

        oracle = [list(nx_community.asyn_lpa_communities(graph.copy(), seed=seed + run)) for run in range(runs)]
        assert result["modularity_values"] == pytest.approx([nx_community.modularity(graph, part) for part in oracle])
        # greedy modularity (the silent replacement) gives one modularity four times; asyn_lpa does not
        assert len({round(value, 9) for value in result["modularity_values"]}) > 1
        assert result["stability_score"] < 1.0

    def test_stability_of_fluid_communities_uses_the_seeds_too(self) -> None:
        graph = _graph()

        result = community_stability(graph, method="fluid", n_runs=4, seed=5, k=3)

        oracle = [list(nx_community.asyn_fluidc(graph.copy(), 3, seed=5 + run)) for run in range(4)]
        assert result["modularity_values"] == pytest.approx([nx_community.modularity(graph, part) for part in oracle])

    @pytest.mark.parametrize("method", ["label_propagation", "girvan_newman"])
    def test_stability_of_a_deterministic_method_is_one_with_that_methods_modularity(self, method: str) -> None:
        graph = _graph()
        direct = {"label_propagation": label_propagation_communities, "girvan_newman": girvan_newman_communities}[
            method
        ]

        result = community_stability(graph, method=method, n_runs=3, seed=2)

        assert result["stability_score"] == pytest.approx(1.0)
        assert result["modularity_values"] == pytest.approx([modularity(graph, direct(graph))] * 3)
        greedy_value = modularity(graph, greedy_modularity_communities(graph))
        assert abs(result["modularity_values"][0] - greedy_value) > 1e-6  # not the greedy fallback's value


class TestWorkflowConfig:
    def _workflow(self, **community) -> NetworkWorkflow:
        config = NetworkWorkflowConfig(community=CommunityDetectionConfig(**community))
        edges = [(u, v) for u, v in _graph().edges]
        return NetworkWorkflow(config).build_network(edges=edges).detect_communities()

    def _run(self, **community) -> dict:
        return self._workflow(**community).results["communities"]

    @pytest.mark.parametrize("seed", [5, 6, 7, 8])
    def test_the_configured_seed_reaches_the_method(self, seed: int) -> None:
        workflow = self._workflow(method="asyn_lpa", random_state=seed)

        # the oracle runs on a copy of the very graph the workflow built (asyn_lpa depends on neighbour order, and
        # the community functions work on a copy, as NetworkX itself does for a graph's .copy())
        built = workflow.results["network"].graph
        assert _groups(workflow.results["communities"]) == _as_groups(
            nx_community.asyn_lpa_communities(built.copy(), seed=seed)
        )

    def test_the_target_number_of_communities_is_passed_under_each_methods_own_keyword(self) -> None:
        fluid = self._run(method="fluid", n_communities=3, random_state=2)
        girvan_newman = self._run(method="girvan_newman", n_communities=4)

        assert len(set(fluid.values())) == 3
        assert len(set(girvan_newman.values())) >= 4

    @pytest.mark.parametrize("method", ["louvain", "greedy", "label_propagation", "asyn_lpa"])
    def test_a_target_number_of_communities_the_method_cannot_honour_raises(self, method: str) -> None:
        with pytest.raises(ValueError, match="n_communities"):
            self._run(method=method, n_communities=3)

    @pytest.mark.parametrize("method", ["label_propagation", "girvan_newman", "asyn_lpa", "fluid"])
    def test_a_resolution_the_method_does_not_take_raises(self, method: str) -> None:
        with pytest.raises(ValueError, match="resolution"):
            self._run(method=method, resolution=1.5)

    def test_a_resolution_greedy_takes_is_passed_on(self) -> None:
        low = self._run(method="greedy", resolution=0.3)
        high = self._run(method="greedy", resolution=3.0)

        assert len(set(high.values())) > len(set(low.values()))


class TestCompareCommunities:
    def test_matches_scikit_learn_on_random_assignments(self) -> None:
        sklearn_metrics = pytest.importorskip("sklearn.metrics")
        rng = np.random.default_rng(20261028)
        for _ in range(200):
            n = int(rng.integers(1, 50))
            first = rng.integers(0, int(rng.integers(1, 7)), size=n)
            second = rng.integers(0, int(rng.integers(1, 7)), size=n)
            nodes = [f"n{i}" for i in range(n)]

            got = compare_communities(
                {node: int(label) for node, label in zip(nodes, first)},
                {node: int(label) for node, label in zip(nodes, second)},
            )

            assert got["normalized_mutual_information"] == pytest.approx(
                sklearn_metrics.normalized_mutual_info_score(first, second), abs=1e-12
            )
            assert got["adjusted_rand_index"] == pytest.approx(
                sklearn_metrics.adjusted_rand_score(first, second), abs=1e-12
            )

    def test_a_relabelled_copy_of_a_partition_is_identical(self) -> None:
        first = {"a": 0, "b": 0, "c": 1, "d": 1, "e": 2}
        second = {"a": 7, "b": 7, "c": 3, "d": 3, "e": 9}

        assert compare_communities(first, second) == {"normalized_mutual_information": 1.0, "adjusted_rand_index": 1.0}

    def test_only_the_shared_nodes_are_compared_and_disjoint_assignments_score_zero(self) -> None:
        first = {"a": 0, "b": 0, "c": 1, "x": 5}
        second = {"a": 0, "b": 0, "c": 1, "y": 6}

        assert compare_communities(first, second)["adjusted_rand_index"] == 1.0
        assert compare_communities({"a": 0}, {"b": 0}) == {
            "normalized_mutual_information": 0.0,
            "adjusted_rand_index": 0.0,
        }

    def test_independent_assignments_score_near_zero_not_the_label_agreement(self) -> None:
        first = {f"n{i}": i % 2 for i in range(40)}
        second = {f"n{i}": (i // 2) % 2 for i in range(40)}  # independent of the first, but half the labels agree

        got = compare_communities(first, second)

        assert got["normalized_mutual_information"] == pytest.approx(0.0, abs=1e-12)
        assert abs(got["adjusted_rand_index"]) < 0.05

    def test_it_does_not_need_scikit_learn(self) -> None:
        sklearn_metrics = pytest.importorskip("sklearn.metrics")
        script = """
import json, sys
sys.modules["sklearn"] = None  # any `import sklearn...` now raises ImportError

from metainformant.networks.analysis.community import compare_communities

relabelled = compare_communities({"a": 0, "b": 0, "c": 1, "d": 1}, {"a": 1, "b": 1, "c": 0, "d": 0})
partial = compare_communities({"a": 0, "b": 0, "c": 1, "d": 1}, {"a": 0, "b": 1, "c": 1, "d": 1})
print(json.dumps([relabelled, partial]))
"""
        completed = subprocess.run(
            [sys.executable, "-c", script],
            capture_output=True,
            text=True,
            env=_subprocess_env(),
            timeout=300,
            check=False,
        )
        assert completed.returncode == 0, completed.stderr
        relabelled, partial = json.loads(completed.stdout.strip().splitlines()[-1])

        # the old fallback reported the share of equal labels (0.0) for a relabelled copy
        assert relabelled == {"normalized_mutual_information": 1.0, "adjusted_rand_index": 1.0}
        first, second = [0, 0, 1, 1], [0, 1, 1, 1]  # contingency [[1, 1], [0, 2]]
        assert partial["adjusted_rand_index"] == pytest.approx(
            sklearn_metrics.adjusted_rand_score(first, second), abs=1e-12
        )
        assert partial["normalized_mutual_information"] == pytest.approx(
            sklearn_metrics.normalized_mutual_info_score(first, second), abs=1e-12
        )
        assert 0.0 < partial["normalized_mutual_information"] < 1.0


class TestModularity:
    def test_a_partition_gets_the_networkx_value_weights_included(self) -> None:
        graph = nx.Graph()
        graph.add_weighted_edges_from([("a", "b", 5.0), ("b", "c", 0.1), ("c", "d", 5.0), ("a", "c", 0.1)])
        parts = [["a", "b"], ["c", "d"]]

        assert modularity(graph, parts) == pytest.approx(nx_community.modularity(graph, parts))
        assert modularity(graph, parts) > 0.4

    def test_overlapping_communities_are_an_error_not_a_number(self) -> None:
        graph = _graph()
        nodes = list(graph)

        with pytest.raises(nx.NetworkXError):
            modularity(graph, [nodes[:20], nodes[10:]])  # nodes 10-19 are in both

    def test_unlisted_nodes_are_singletons_and_unknown_names_are_ignored(self) -> None:
        graph = _graph()
        nodes = list(graph)
        listed = [nodes[:8], nodes[8:16]]
        completed = listed + [[node] for node in nodes[16:]]

        # a partial assignment used to get an unweighted edge-count approximation of only the listed communities
        assert modularity(graph, listed) == pytest.approx(nx_community.modularity(graph, completed))
        assert modularity(graph, listed + [["not_a_node", "neither"]]) == pytest.approx(
            nx_community.modularity(graph, completed)
        )
        mapping = {node: 0 for node in nodes[:8]} | {node: 1 for node in nodes[8:16]} | {"ghost": 2}
        assert modularity(graph, mapping) == pytest.approx(nx_community.modularity(graph, completed))

    def test_an_empty_assignment_or_an_edgeless_graph_is_zero(self) -> None:
        assert modularity(_graph(), []) == 0.0
        assert modularity(nx.empty_graph(4), [[0, 1], [2, 3]]) == 0.0
