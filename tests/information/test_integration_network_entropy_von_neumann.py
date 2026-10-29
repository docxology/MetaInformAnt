"""``integration.networks.network_entropy(method="von_neumann")`` is the Von Neumann entropy of a density matrix.

The private ``_von_neumann_entropy`` took the eigenvalues of the row-normalised adjacency (random-walk transition)
matrix, kept the positive ones and summed ``-lambda * log2(lambda)``. That matrix is not a density matrix: a simple
graph has ``trace(P) = 0``, the top eigenvalue is always exactly 1 (contributing ``-1 * log2(1) = 0``) and the
remaining positive eigenvalues are not probabilities. Every complete graph and every star came out at 0 (to rounding
noise, sometimes slightly negative) instead of the known closed forms, and the 6-cycle and two triangles got 1 and 0.

The function now builds the Laplacian density matrix ``rho = L / trace(L)`` (Braunstein, Ghosh and Severini, 2006) by
reusing ``network_info.information_flow.network_entropy``. The oracles are the closed forms for complete graphs and
stars, and NetworkX's Laplacian fed to NumPy's symmetric eigensolver with a separate entropy sum. No test doubles.
"""

from __future__ import annotations

import math
import random

import numpy as np
import pytest

nx = pytest.importorskip("networkx")

from metainformant.information.integration.networks import network_entropy  # noqa: E402
from metainformant.information.network_info.information_flow import network_entropy as laplacian_entropy  # noqa: E402


def _bits(eigenvalues) -> float:
    positive = [float(value) for value in eigenvalues if value > 1e-14]
    return -sum(value * math.log2(value) for value in positive)


def _oracle(graph) -> float:
    """Von Neumann entropy of rho = L / trace(L) from networkx's Laplacian and NumPy's symmetric eigensolver."""
    laplacian = nx.laplacian_matrix(graph, weight="weight").toarray().astype(float)
    return _bits(np.linalg.eigvalsh(laplacian / np.trace(laplacian)))


def _star_closed_form(leaves: int) -> float:
    # Laplacian spectrum of a star with m leaves: 0, 1 (m - 1 times), m + 1; trace 2m
    small, large = 1 / (2 * leaves), (leaves + 1) / (2 * leaves)
    return -(leaves - 1) * small * math.log2(small) - large * math.log2(large)


@pytest.mark.parametrize("n", [3, 4, 5, 6, 8, 12, 20])
def test_a_complete_graph_has_entropy_log2_of_n_minus_one(n: int) -> None:
    assert network_entropy(nx.complete_graph(n), method="von_neumann") == pytest.approx(math.log2(n - 1), abs=1e-12)


@pytest.mark.parametrize("leaves", [2, 3, 4, 6, 9, 15])
def test_a_star_matches_its_closed_form(leaves: int) -> None:
    assert network_entropy(nx.star_graph(leaves), method="von_neumann") == pytest.approx(
        _star_closed_form(leaves), abs=1e-12
    )


def test_two_nodes_joined_by_one_edge_are_a_pure_state() -> None:
    # L has eigenvalues 0 and 2, so rho = diag(1, 0) and S = log2(2 - 1) = 0
    assert network_entropy(nx.path_graph(2), method="von_neumann") == pytest.approx(0.0, abs=1e-12)


def test_the_entropy_is_never_negative() -> None:
    # the old eigenvalue sum returned -3e-16 for K_8 and K_12
    for n in (3, 5, 8, 12, 30):
        assert network_entropy(nx.complete_graph(n), method="von_neumann") >= 0.0
        assert network_entropy(nx.star_graph(n), method="von_neumann") >= 0.0


def test_it_is_the_same_value_the_network_info_function_reports() -> None:
    graph = nx.gnp_random_graph(12, 0.4, seed=11)
    matrix = nx.to_numpy_array(graph).tolist()

    assert network_entropy(graph, method="von_neumann") == laplacian_entropy(matrix)["entropy"]


def test_random_weighted_graphs_match_the_networkx_oracle() -> None:
    rng = random.Random(20261028)
    checked = 0
    for _ in range(40):
        graph = nx.gnp_random_graph(rng.randint(3, 14), 0.4, seed=rng.randint(0, 10**6))
        for u, v in graph.edges:
            graph[u][v]["weight"] = rng.uniform(0.2, 3.0)
        if graph.number_of_edges() == 0:
            continue

        assert network_entropy(graph, method="von_neumann") == pytest.approx(_oracle(graph), abs=1e-10)
        checked += 1
    assert checked >= 30


def test_the_entropy_depends_on_the_structure_and_not_only_on_the_degrees() -> None:
    cycle = nx.cycle_graph(6)  # 2-regular
    triangles = nx.disjoint_union(nx.complete_graph(3), nx.complete_graph(3))  # 2-regular as well

    cycle_entropy = network_entropy(cycle, method="von_neumann")
    triangles_entropy = network_entropy(triangles, method="von_neumann")

    assert cycle_entropy == pytest.approx(_oracle(cycle), abs=1e-10)
    assert triangles_entropy == pytest.approx(2.0, abs=1e-12)  # four equal eigenvalues 1/4 of rho
    assert abs(cycle_entropy - triangles_entropy) > 0.1


def test_a_numpy_adjacency_array_gives_the_same_entropy_as_the_graph() -> None:
    graph = nx.gnp_random_graph(9, 0.5, seed=5)

    from_array = network_entropy(nx.to_numpy_array(graph), method="von_neumann")

    assert from_array == pytest.approx(network_entropy(graph, method="von_neumann"), abs=1e-14)
    assert network_entropy(nx.to_numpy_array(nx.complete_graph(7)), method="von_neumann") == pytest.approx(
        math.log2(6), abs=1e-12
    )


def test_the_edge_weights_are_used_and_only_their_relative_size_matters() -> None:
    graph = nx.gnp_random_graph(8, 0.5, seed=2)
    for u, v in graph.edges:
        graph[u][v]["weight"] = 1.0 + 0.5 * ((u + v) % 3)
    scaled = nx.Graph()
    scaled.add_weighted_edges_from((u, v, 7.5 * data["weight"]) for u, v, data in graph.edges(data=True))
    unweighted = nx.Graph(graph.edges)

    reference = network_entropy(graph, method="von_neumann")

    assert network_entropy(scaled, method="von_neumann") == pytest.approx(reference, abs=1e-12)
    assert reference == pytest.approx(_oracle(graph), abs=1e-10)
    assert abs(reference - network_entropy(unweighted, method="von_neumann")) > 1e-4


def test_parallel_edges_of_a_multigraph_add_their_weights() -> None:
    multi = nx.MultiGraph([(0, 1), (0, 1), (1, 2), (2, 3)])
    merged = nx.Graph()
    merged.add_weighted_edges_from([(0, 1, 2.0), (1, 2, 1.0), (2, 3, 1.0)])

    assert network_entropy(multi, method="von_neumann") == pytest.approx(_oracle(merged), abs=1e-12)


def test_self_loops_do_not_change_the_entropy() -> None:
    graph = nx.gnp_random_graph(8, 0.5, seed=3)
    with_loops = graph.copy()
    with_loops.add_edges_from((node, node) for node in graph.nodes)

    assert network_entropy(with_loops, method="von_neumann") == pytest.approx(
        network_entropy(graph, method="von_neumann"), abs=1e-12
    )


def test_relabelling_the_nodes_does_not_change_the_entropy() -> None:
    graph = nx.gnp_random_graph(9, 0.5, seed=4)
    order = list(range(9))
    random.Random(1).shuffle(order)
    relabelled = nx.relabel_nodes(graph, dict(enumerate(order)))

    assert network_entropy(graph, method="von_neumann") == pytest.approx(
        network_entropy(relabelled, method="von_neumann"), abs=1e-12
    )


@pytest.mark.parametrize("graph", [nx.Graph(), nx.empty_graph(1), nx.empty_graph(5)])
def test_graphs_without_edges_have_zero_entropy(graph) -> None:
    assert network_entropy(graph, method="von_neumann") == 0.0


def test_a_directed_graph_is_rejected_instead_of_being_given_a_meaningless_value() -> None:
    cycle = nx.DiGraph([(0, 1), (1, 2), (2, 0)])

    with pytest.raises(ValueError, match="undirected"):
        network_entropy(cycle, method="von_neumann")

    # the explicit conversion the message asks for analyses the undirected skeleton: a triangle, log2(3 - 1)
    assert network_entropy(cycle.to_undirected(), method="von_neumann") == pytest.approx(1.0, abs=1e-12)


def test_the_shannon_method_is_unchanged_and_unknown_methods_still_raise() -> None:
    star = nx.star_graph(4)  # degrees 4, 1, 1, 1, 1: P(4) = 1/5 and P(1) = 4/5
    expected = -(0.2 * math.log2(0.2) + 0.8 * math.log2(0.8))

    assert network_entropy(star) == pytest.approx(expected, abs=1e-12)
    with pytest.raises(ValueError, match="Unknown entropy method"):
        network_entropy(star, method="renyi")
