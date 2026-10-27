"""Regression tests: ``network_entropy`` computes the Von Neumann entropy its docstring describes.

The docstring promised the Von Neumann entropy of the graph's density matrix ("power iteration for eigenvalue
estimation"), but the code never computed an eigenvalue. It built ``rho = L / trace(L)`` and then used the *diagonal*
of ``rho`` (the normalised degrees) as if it were the spectrum, so the result was just the Shannon entropy of the
degree distribution: it did not depend on the structure of the network at all (a 6-cycle and two disjoint triangles,
both 2-regular, got the same value, and every complete graph got exactly ``log2(n)`` instead of ``log2(n - 1)``).

The function now does what it says: the density matrix ``rho = L / trace(L)`` of the weighted Laplacian (Braunstein,
Ghosh and Severini 2006) is eigendecomposed (``numpy.linalg.eigvalsh``, or a pure-Python Jacobi solver without NumPy)
and ``S = -sum(lambda log2 lambda)``. The oracles are the closed forms for complete graphs and stars, and
``networkx``'s Laplacian fed to NumPy's eigensolver and a separate entropy sum. No test doubles.
"""

from __future__ import annotations

import math
import random

import networkx as nx
import numpy as np
import pytest

from metainformant.information.network_info import information_flow
from metainformant.information.network_info.information_flow import network_entropy


def _bits(eigenvalues) -> float:
    positive = [float(value) for value in eigenvalues if value > 1e-14]
    return -sum(value * math.log2(value) for value in positive)


def _oracle(graph: nx.Graph) -> float:
    """Von Neumann entropy of rho = L / trace(L) from networkx's Laplacian and NumPy's symmetric eigensolver."""
    laplacian = nx.laplacian_matrix(graph, weight="weight").toarray().astype(float)
    return _bits(np.linalg.eigvalsh(laplacian / np.trace(laplacian)))


def _matrix(graph: nx.Graph) -> list[list[float]]:
    return nx.to_numpy_array(graph, weight="weight").tolist()


@pytest.mark.parametrize("n", [3, 4, 5, 6, 8, 12])
def test_a_complete_graph_has_entropy_log2_of_n_minus_one(n: int) -> None:
    result = network_entropy(_matrix(nx.complete_graph(n)))

    assert result["entropy"] == pytest.approx(math.log2(n - 1), abs=1e-12)
    assert result["normalized_entropy"] == pytest.approx(math.log2(n - 1) / math.log2(n), abs=1e-12)
    assert result["density"] == 1.0 and result["n_nodes"] == n


@pytest.mark.parametrize("n", [4, 5, 7, 10])
def test_a_star_matches_its_closed_form(n: int) -> None:
    # Laplacian spectrum of a star on n nodes: 0, 1 (n - 2 times), n; trace 2 (n - 1)
    small, large = 1 / (2 * (n - 1)), n / (2 * (n - 1))
    expected = -(n - 2) * small * math.log2(small) - large * math.log2(large)

    assert network_entropy(_matrix(nx.star_graph(n - 1)))["entropy"] == pytest.approx(expected, abs=1e-12)


def test_random_weighted_graphs_match_the_networkx_oracle() -> None:
    rng = random.Random(20261027)
    checked = 0
    for _ in range(40):
        graph = nx.gnp_random_graph(rng.randint(3, 14), 0.4, seed=rng.randint(0, 10**6))
        for u, v in graph.edges:
            graph[u][v]["weight"] = rng.uniform(0.2, 3.0)
        if graph.number_of_edges() == 0:
            continue

        assert network_entropy(_matrix(graph))["entropy"] == pytest.approx(_oracle(graph), abs=1e-10)
        checked += 1
    assert checked >= 30


def test_the_entropy_depends_on_the_structure_and_not_only_on_the_degrees() -> None:
    cycle = nx.cycle_graph(6)  # 2-regular
    triangles = nx.disjoint_union(nx.complete_graph(3), nx.complete_graph(3))  # 2-regular as well

    cycle_entropy = network_entropy(_matrix(cycle))["entropy"]
    triangles_entropy = network_entropy(_matrix(triangles))["entropy"]

    assert cycle_entropy == pytest.approx(_oracle(cycle), abs=1e-10)
    assert triangles_entropy == pytest.approx(2.0, abs=1e-12)  # four equal eigenvalues 1/4 of rho
    assert abs(cycle_entropy - triangles_entropy) > 0.1


def test_disconnected_graphs_match_the_oracle() -> None:
    graph = nx.disjoint_union_all([nx.path_graph(4), nx.complete_graph(3), nx.star_graph(3)])

    assert network_entropy(_matrix(graph))["entropy"] == pytest.approx(_oracle(graph), abs=1e-10)


def test_relabelling_the_nodes_does_not_change_the_entropy() -> None:
    graph = nx.gnp_random_graph(9, 0.5, seed=4)
    order = list(range(9))
    random.Random(1).shuffle(order)
    relabelled = nx.relabel_nodes(graph, dict(enumerate(order)))

    assert network_entropy(_matrix(graph))["entropy"] == pytest.approx(
        network_entropy(_matrix(relabelled))["entropy"], abs=1e-12
    )


def test_scale_sign_and_self_loops_of_the_weights() -> None:
    graph = nx.gnp_random_graph(8, 0.5, seed=2)
    for u, v in graph.edges:
        graph[u][v]["weight"] = 1.0 + 0.5 * ((u + v) % 3)
    matrix = np.array(_matrix(graph))
    reference = network_entropy(matrix.tolist())["entropy"]

    assert network_entropy((7.5 * matrix).tolist())["entropy"] == pytest.approx(reference, abs=1e-12)  # scale-free
    assert network_entropy((-matrix).tolist())["entropy"] == pytest.approx(reference, abs=1e-12)  # |weight| is used
    with_loops = matrix + np.diag(np.arange(1.0, 9.0))
    assert network_entropy(with_loops.tolist())["entropy"] == pytest.approx(reference, abs=1e-12)  # D - W drops loops


def test_networks_without_edges_have_zero_entropy() -> None:
    assert network_entropy([]) == {"entropy": 0.0, "normalized_entropy": 0.0, "n_nodes": 0, "density": 0.0}
    assert network_entropy([[0.0, 0.0], [0.0, 0.0]])["entropy"] == 0.0
    assert network_entropy([[0.0]])["entropy"] == 0.0
    assert network_entropy([[3.0]])["entropy"] == 0.0  # a lone self-loop is no edge


def test_non_square_and_asymmetric_matrices_are_rejected() -> None:
    with pytest.raises(ValueError, match="square"):
        network_entropy([[0.0, 1.0, 0.0], [1.0, 0.0]])
    with pytest.raises(ValueError, match="symmetric"):
        network_entropy([[0.0, 1.0], [0.0, 0.0]])


def test_the_pure_python_solver_agrees_with_numpy() -> None:
    rng = np.random.default_rng(3)
    for n in (2, 3, 7, 15, 25):
        symmetric = rng.normal(size=(n, n))
        symmetric = (symmetric + symmetric.T) / 2

        jacobi = np.sort(information_flow._jacobi_eigenvalues(symmetric.tolist()))

        assert jacobi == pytest.approx(np.linalg.eigvalsh(symmetric), abs=1e-10)


def test_the_entropy_from_the_pure_python_solver_equals_the_reported_entropy() -> None:
    graph = nx.gnp_random_graph(11, 0.4, seed=8)
    for u, v in graph.edges:
        graph[u][v]["weight"] = 0.5 + (u * v) % 4
    laplacian = nx.laplacian_matrix(graph, weight="weight").toarray().astype(float)
    rho = (laplacian / np.trace(laplacian)).tolist()

    assert _bits(information_flow._jacobi_eigenvalues(rho)) == pytest.approx(
        network_entropy(_matrix(graph))["entropy"], abs=1e-10
    )


def test_the_documentation_no_longer_promises_power_iteration() -> None:
    doc = " ".join((network_entropy.__doc__ or "").split())

    assert "power iteration" not in doc and "eigendecomposition" in doc and "log2(n - 1)" in doc
