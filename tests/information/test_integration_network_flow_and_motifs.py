"""Information flow is per source and never silently empty or NaN; motif analysis is seeded and respects direction.

* ``information_flow(method="random_walk")`` recorded flows only when ``target_nodes`` was given, and then walked the
  mixture of all sources: every source got the same flow to a target, and with the default ``target_nodes=None`` (or
  with ``source_nodes=None``) the flow matrix was all zeros and the path-length entropy 0. Each source is now walked on
  its own and ``flow_matrix[s][t]`` is the sum over the steps of the probability of being at ``t``.
* A source node that is not in the graph was skipped silently; with no valid source the diffusion signal was all zeros
  and its renormalisation returned NaN. Unknown names (and empty lists) now raise ``ValueError``.
  ``method="diffusion"`` also ignored ``target_nodes``.
* ``network_motif_information`` built its null model with NetworkX's configuration model without a seed (the global
  ``random`` state: results changed from call to call and the caller's random state was consumed), subsampled big
  graphs with a hidden fixed seed after listing every combination, never said it had sampled, built the null of a
  directed graph as an undirected one, labelled directed subgraphs by edge count and total degree (the feed-forward
  loop and the three-cycle collided) and gave ``+inf`` to motifs below the null mean.

Oracles: the matrix power sum ``sum_k P^k`` computed in the test, and hand-built three-node motifs. No test doubles.
"""

from __future__ import annotations

import math
import random

import numpy as np
import pytest

nx = pytest.importorskip("networkx")

from metainformant.information.integration.networks import (  # noqa: E402
    _count_motifs,
    information_flow,
    network_motif_information,
)


def _visit_oracle(graph, steps: int) -> np.ndarray:
    """sum_{k=1..steps} P^k for the random walk on ``graph`` (rows indexed by ``list(graph.nodes)``)."""
    adjacency = nx.to_numpy_array(graph)
    degrees = adjacency.sum(axis=1)
    transition = adjacency / np.where(degrees == 0, 1, degrees)[:, None]
    power, total = np.eye(len(graph)), np.zeros((len(graph), len(graph)))
    for _ in range(steps):
        power = power @ transition
        total += power
    return total


class TestRandomWalkFlow:
    def test_each_source_has_its_own_flow(self) -> None:
        graph = nx.path_graph(5)
        oracle = _visit_oracle(graph, 6)

        result = information_flow(graph, source_nodes=[0, 2], target_nodes=[4], method="random_walk", steps=6)

        flows = result["flow_matrix"]
        assert flows[0][4] == pytest.approx(oracle[0, 4], abs=1e-12)
        assert flows[2][4] == pytest.approx(oracle[2, 4], abs=1e-12)
        assert flows[2][4] > flows[0][4] > 0.0  # node 2 is nearer to node 4 than node 0 is; the old rows were equal

    def test_the_default_arguments_give_the_whole_flow_matrix_not_zeros(self) -> None:
        graph = nx.karate_club_graph()
        oracle = _visit_oracle(graph, 5)

        result = information_flow(graph, steps=5)

        for source in (0, 7, 33):
            for target in (0, 2, 33):
                assert result["flow_matrix"][source][target] == pytest.approx(oracle[source, target], abs=1e-12)
        assert result["path_length_entropy"] > 1.0
        assert len(result["flow_matrix"]) == len(graph)

    def test_sources_without_targets_and_targets_without_sources(self) -> None:
        graph = nx.path_graph(5)
        oracle = _visit_oracle(graph, 5)

        only_sources = information_flow(graph, source_nodes=[0], steps=5)["flow_matrix"]
        only_targets = information_flow(graph, target_nodes=[2], steps=5)["flow_matrix"]

        assert only_sources[0] == pytest.approx({t: oracle[0, t] for t in range(5)}, abs=1e-12)
        assert {s: row[2] for s, row in only_targets.items()} == pytest.approx(
            {s: oracle[s, 2] for s in range(5)}, abs=1e-12
        )

    def test_the_adjacency_array_input_gives_the_same_flows(self) -> None:
        graph = nx.cycle_graph(6)

        from_graph = information_flow(graph, steps=4)["flow_matrix"]
        from_array = information_flow(nx.to_numpy_array(graph), steps=4)["flow_matrix"]

        for source, row in from_graph.items():
            assert from_array[source] == pytest.approx(row)

    def test_flow_is_expected_visits_so_it_adds_up_over_a_connected_graph(self) -> None:
        graph = nx.complete_graph(5)  # no isolated nodes: a walker is always somewhere
        steps = 7

        result = information_flow(graph, steps=steps)

        for row in result["flow_matrix"].values():
            assert sum(row.values()) == pytest.approx(steps, abs=1e-9)


class TestNodeValidation:
    @pytest.mark.parametrize("method", ["random_walk", "diffusion"])
    @pytest.mark.parametrize("role", ["source_nodes", "target_nodes"])
    def test_unknown_nodes_raise(self, method: str, role: str) -> None:
        with pytest.raises(ValueError, match="not in the graph"):
            information_flow(nx.path_graph(4), method=method, steps=3, **{role: ["nope"]})

    @pytest.mark.parametrize("method", ["random_walk", "diffusion"])
    def test_a_partly_unknown_list_raises_too(self, method: str) -> None:
        with pytest.raises(ValueError, match="'ghost'"):
            information_flow(nx.path_graph(4), source_nodes=[0, "ghost"], method=method, steps=3)

    @pytest.mark.parametrize("method", ["random_walk", "diffusion"])
    def test_an_empty_list_raises(self, method: str) -> None:
        with pytest.raises(ValueError, match="empty"):
            information_flow(nx.path_graph(4), source_nodes=[], method=method, steps=3)


class TestDiffusion:
    def test_target_nodes_select_what_is_reported(self) -> None:
        graph = nx.path_graph(6)

        everything = information_flow(graph, source_nodes=[0], method="diffusion", steps=8)
        selected = information_flow(graph, source_nodes=[0], target_nodes=[3, 5], method="diffusion", steps=8)

        assert set(everything["node_flows"]) == set(range(6))
        assert set(selected["node_flows"]) == {3, 5}
        assert selected["target_nodes"] == [3, 5]
        for node in (3, 5):
            assert selected["node_flows"][node] == pytest.approx(everything["node_flows"][node])

    def test_the_signal_is_a_distribution_and_never_nan(self) -> None:
        result = information_flow(nx.karate_club_graph(), source_nodes=[0, 33], method="diffusion", steps=20)

        finals = [flow["final_signal"] for flow in result["node_flows"].values()]
        assert all(math.isfinite(value) and value >= 0.0 for value in finals)
        assert sum(finals) == pytest.approx(1.0)


class TestMotifSeeding:
    GRAPH = nx.gnp_random_graph(14, 0.3, seed=1)

    def test_the_same_seed_gives_the_same_result_and_another_seed_differs(self) -> None:
        first = network_motif_information(self.GRAPH, motif_size=3, n_random=20, random_state=5)
        second = network_motif_information(self.GRAPH, motif_size=3, n_random=20, random_state=5)
        third = network_motif_information(self.GRAPH, motif_size=3, n_random=20, random_state=6)

        assert first == second
        assert first["z_scores"] != third["z_scores"]
        assert first["motifs_found"] == third["motifs_found"]  # the observed counts do not depend on the seed

    def test_global_random_state_is_left_alone(self) -> None:
        before = (random.getstate(), np.random.get_state()[1].tobytes())

        network_motif_information(self.GRAPH, motif_size=3, n_random=5)
        network_motif_information(self.GRAPH, motif_size=3, n_random=5, random_state=1)

        assert (random.getstate(), np.random.get_state()[1].tobytes()) == before

    def test_a_generator_can_be_passed(self) -> None:
        by_seed = network_motif_information(self.GRAPH, motif_size=3, n_random=10, random_state=9)
        by_generator = network_motif_information(
            self.GRAPH, motif_size=3, n_random=10, random_state=np.random.default_rng(9)
        )

        assert by_generator == by_seed

    def test_a_small_graph_is_counted_exactly_and_says_so(self) -> None:
        result = network_motif_information(nx.complete_graph(5), motif_size=3, n_random=3, random_state=0)

        assert result["sampled"] is False
        assert result["n_combinations_examined"] == math.comb(5, 3)
        assert result["motifs_found"] == {"n3_e3_d222": 10}  # ten triangles

    def test_sampling_is_reported_seeded_and_does_not_list_every_combination(self) -> None:
        graph = nx.gnp_random_graph(400, 0.02, seed=4)  # 10.6 million triples; the old code listed them all

        first = _count_motifs(graph, 3, np.random.default_rng(1), max_subgraphs=1000)
        again = _count_motifs(graph, 3, np.random.default_rng(1), max_subgraphs=1000)
        other = _count_motifs(graph, 3, np.random.default_rng(2), max_subgraphs=1000)

        assert first == again
        assert first[1] == 1000 and first[2] is True
        assert first[0] != other[0]  # a different seed examines different node triples

        result = network_motif_information(graph, motif_size=3, n_random=2, random_state=1, max_subgraphs=1000)
        assert result["sampled"] is True and result["n_combinations_examined"] == 1000


class TestDirectedMotifs:
    def test_the_feed_forward_loop_and_the_three_cycle_are_different_motifs(self) -> None:
        feed_forward = nx.DiGraph([(0, 1), (1, 2), (0, 2)])
        cycle = nx.DiGraph([(0, 1), (1, 2), (2, 0)])

        ffl = _count_motifs(feed_forward, 3)[0]
        three_cycle = _count_motifs(cycle, 3)[0]

        assert sum(ffl.values()) == 1 and sum(three_cycle.values()) == 1
        assert set(ffl) != set(three_cycle)  # both had the key n3_e3_d222

    def test_isomorphic_directed_subgraphs_share_one_class(self) -> None:
        graph = nx.DiGraph([(0, 1), (1, 2), (0, 2), (10, 11), (11, 12), (10, 12)])  # two feed-forward loops

        counts = _count_motifs(graph, 3)[0]

        assert sorted(counts.values()) == [2]

    def test_the_null_model_of_a_directed_graph_is_directed(self) -> None:
        graph = nx.gnp_random_graph(12, 0.25, seed=2, directed=True)

        result = network_motif_information(graph, motif_size=3, n_random=10, random_state=3)

        assert result["n_random_graphs"] == 10
        assert all("dir" in key for key in result["motifs_found"])
        # most classes of the observed graph also occur in the directed null, so their z-scores are finite (the
        # old undirected null gave inf to classes it could not produce); a class the null never produces is +inf
        z_values = list(result["z_scores"].values())
        assert not any(math.isnan(z) for z in z_values)
        assert sum(math.isfinite(z) for z in z_values) >= 0.75 * len(z_values)

    def test_undirected_motif_keys_are_unchanged(self) -> None:
        path = nx.path_graph(3)
        triangle = nx.complete_graph(3)

        assert _count_motifs(path, 3)[0] == {"n3_e2_d112": 1}
        assert _count_motifs(triangle, 3)[0] == {"n3_e3_d222": 1}
        assert _count_motifs(nx.star_graph(3), 4)[0] == {"n4_e3_d1113": 1}
