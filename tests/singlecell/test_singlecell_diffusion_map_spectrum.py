"""The diffusion maps return eigenvectors of the Markov matrix, built from one weight per kNN edge.

``diffusion_map_reduction`` used to hand a non-symmetric matrix to ``numpy.linalg.eigh`` (which reads only its
lower triangle and so decomposes some other matrix) and divide the result by the eigenvalues, and
``compute_diffusion_map`` added the transpose of a directed kNN graph, which doubled the weight of every
edge that both cells list. These tests compare against an explicit eigendecomposition of the Markov matrix
built independently here.

Eigenvectors are only defined up to sign and scale (and, in a degenerate eigenspace, rotation), and a solver
may return either sign: every comparison is made up to a scalar factor, with tolerances, on spectra whose
leading eigenvalues are checked to be well separated first. Real data and real linear algebra only.
"""

from __future__ import annotations

import numpy as np
import pytest
from scipy import sparse
from scipy.spatial.distance import pdist, squareform

from metainformant.core.utils import errors
from metainformant.singlecell.analysis import nonlinear_methods
from metainformant.singlecell.analysis.nonlinear_methods import (
    compute_diffusion_map,
    compute_neighbors,
    diffusion_map_reduction,
)
from metainformant.singlecell.data.preprocessing import SingleCellData

N_COMPONENTS = 4


def _curve(n_cells: int = 40, seed: int = 0) -> SingleCellData:
    """Cells along a noisy curve: a smooth manifold gives a diffusion spectrum with distinct leading eigenvalues."""
    rng = np.random.default_rng(seed)
    t = np.linspace(0.0, 1.0, n_cells)
    X = np.column_stack([t, t**2, 0.02 * rng.standard_normal((n_cells, 4))])
    return SingleCellData(X)


def _explicit_markov_matrix(X: np.ndarray, sigma: float, alpha: float = 1.0) -> tuple[np.ndarray, np.ndarray]:
    """``P = D^-1 W`` for the Gaussian kernel of width sigma (as a lazy walk when alpha < 1), and the degrees."""
    distances = squareform(pdist(X))
    kernel = np.exp(-(distances**2) / (2 * sigma**2))
    degrees = kernel.sum(axis=1)
    markov = kernel / degrees[:, np.newaxis]
    return alpha * markov + (1 - alpha) * np.eye(len(X)), degrees


def _explicit_spectrum(markov: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Eigenvalues (largest first) and right eigenvectors from the general, non-symmetric eigensolver."""
    values, vectors = np.linalg.eig(markov)
    order = np.argsort(values.real)[::-1]
    return values.real[order], vectors.real[:, order]


def _directions_agree(a: np.ndarray, b: np.ndarray, tolerance: float = 1e-8) -> bool:
    """Whether two vectors are parallel (either sign, any scale)."""
    cosine = abs(float(a @ b)) / (np.linalg.norm(a) * np.linalg.norm(b))
    return cosine > 1.0 - tolerance


def _components(result: SingleCellData) -> np.ndarray:
    return result.obs[[f"DC{k + 1}" for k in range(N_COMPONENTS)]].to_numpy()


class TestDiffusionMapReductionIsTheMarkovSpectrum:
    @pytest.mark.parametrize("alpha", [1.0, 0.6])
    def test_the_components_are_right_eigenvectors_of_the_markov_matrix(self, alpha: float) -> None:
        data = _curve()

        result = diffusion_map_reduction(data, n_components=N_COMPONENTS, n_neighbors=5, alpha=alpha)

        info = result.uns["diffusion_map"]
        markov, _ = _explicit_markov_matrix(data.X, info["sigma"], alpha)
        values, vectors = _explicit_spectrum(markov)
        assert values[0] == pytest.approx(1.0, abs=1e-12)  # the trivial eigenvalue of a Markov matrix
        gaps = -np.diff(values[: N_COMPONENTS + 2])
        assert gaps.min() > 1e-3, "the test needs a spectrum with distinct leading eigenvalues"

        np.testing.assert_allclose(info["eigenvalues"], values[1 : N_COMPONENTS + 1], rtol=0, atol=1e-10)
        components = _components(result)
        for k in range(N_COMPONENTS):
            column = components[:, k]
            # P psi = lambda psi, whatever the sign or scale of psi
            np.testing.assert_allclose(markov @ column, values[k + 1] * column, atol=1e-10)
            assert _directions_agree(column, vectors[:, k + 1])

    def test_the_kernel_width_is_the_median_pairwise_distance(self) -> None:
        data = _curve()

        result = diffusion_map_reduction(data, n_components=N_COMPONENTS, n_neighbors=5)

        assert result.uns["diffusion_map"]["sigma"] == pytest.approx(float(np.median(squareform(pdist(data.X)))))

    def test_the_components_are_orthonormal_in_the_degree_weighted_inner_product(self) -> None:
        data = _curve()

        result = diffusion_map_reduction(data, n_components=N_COMPONENTS, n_neighbors=5)

        _, degrees = _explicit_markov_matrix(data.X, result.uns["diffusion_map"]["sigma"])
        components = _components(result)
        np.testing.assert_allclose(
            components.T @ (degrees[:, np.newaxis] * components), np.eye(N_COMPONENTS), atol=1e-10
        )
        # ... and orthogonal to the trivial constant eigenvector
        np.testing.assert_allclose(degrees @ components, np.zeros(N_COMPONENTS), atol=1e-10)

    def test_each_component_has_its_largest_entry_positive(self) -> None:
        components = _components(diffusion_map_reduction(_curve(), n_components=N_COMPONENTS, n_neighbors=5))

        pivot = np.argmax(np.abs(components), axis=0)
        assert np.all(components[pivot, np.arange(N_COMPONENTS)] > 0)

    def test_the_eigenvalues_are_decreasing_and_below_the_trivial_one(self) -> None:
        eigenvalues = np.array(
            diffusion_map_reduction(_curve(), n_components=N_COMPONENTS, n_neighbors=5).uns["diffusion_map"][
                "eigenvalues"
            ]
        )

        assert np.all(np.diff(eigenvalues) < 0)
        assert eigenvalues[0] < 1.0

    def test_repeated_runs_agree(self) -> None:
        first = _components(diffusion_map_reduction(_curve(), n_components=N_COMPONENTS, n_neighbors=5))
        second = _components(diffusion_map_reduction(_curve(), n_components=N_COMPONENTS, n_neighbors=5))

        np.testing.assert_allclose(first, second, rtol=1e-12, atol=1e-14)


class TestDiffusionMapReductionRejectsWhatItCannotDecompose:
    def test_more_components_than_cells_allow(self) -> None:
        data = _curve(n_cells=8)

        with pytest.raises(errors.ValidationError, match="less than the number of cells"):
            diffusion_map_reduction(data, n_components=8, n_neighbors=5)

    def test_identical_cells_have_no_kernel_width(self) -> None:
        data = SingleCellData(np.ones((12, 3)))

        with pytest.raises(errors.ValidationError, match="identical"):
            diffusion_map_reduction(data, n_components=2, n_neighbors=5)

    def test_alpha_zero_would_be_the_identity_operator(self) -> None:
        with pytest.raises(errors.ValidationError, match="alpha"):
            diffusion_map_reduction(_curve(), n_components=N_COMPONENTS, n_neighbors=5, alpha=0.0)

    def test_infinite_values_are_rejected_before_the_solver(self) -> None:
        data = _curve()
        data.X[3, 1] = np.inf

        with pytest.raises(errors.ValidationError, match="non-finite"):
            diffusion_map_reduction(data, n_components=N_COMPONENTS, n_neighbors=5)


def _neighbor_data(n_cells: int = 60, n_neighbors: int = 8) -> SingleCellData:
    return compute_neighbors(_curve(n_cells=n_cells, seed=1), n_neighbors=n_neighbors)


def _explicit_affinities(indices: np.ndarray, distances: np.ndarray) -> np.ndarray:
    """Dense symmetric affinities: one weight ``1 / (1 + d)`` per pair in which either cell lists the other."""
    n = indices.shape[0]
    affinity = np.zeros((n, n))
    for i in range(n):
        for j, distance in zip(indices[i], distances[i]):
            if i != j:
                affinity[i, j] = affinity[j, i] = 1.0 / (1.0 + distance)
    return affinity


class TestKnnAffinitiesAreNotDoubled:
    def test_mutual_and_one_way_edges_get_one_weight_each(self) -> None:
        # Cell 0 and 1 list each other (a mutual edge); cell 2 lists 1 but is not listed by it (a one-way edge)
        indices = np.array([[0, 1], [1, 0], [2, 1]])
        distances = np.array([[0.0, 1.0], [0.0, 1.0], [0.0, 3.0]])

        affinity = nonlinear_methods._knn_affinity_matrix(indices, distances).toarray()

        assert affinity[0, 1] == affinity[1, 0] == pytest.approx(1.0 / 2.0)  # mutual: 1/(1+1), not 2 * 1/(1+1)
        assert affinity[1, 2] == affinity[2, 1] == pytest.approx(1.0 / 4.0)  # one-way: 1/(1+3)
        assert affinity[0, 2] == 0.0
        assert np.all(np.diag(affinity) == 0.0)

    def test_the_real_graph_has_mutual_edges_and_each_weighs_one_over_one_plus_distance(self) -> None:
        graph = _neighbor_data().uns["neighbors"]
        indices, distances = graph["indices"], graph["distances"]
        mutual = [(i, j) for i in range(len(indices)) for j in indices[i] if i < j and i in indices[j]]
        assert mutual, "the check is only meaningful if cells list each other"

        affinity = nonlinear_methods._knn_affinity_matrix(indices, distances)

        np.testing.assert_allclose(affinity.toarray(), _explicit_affinities(indices, distances), rtol=1e-12)
        assert (affinity != affinity.T).nnz == 0
        assert affinity.max() <= 1.0  # 1 / (1 + d) never exceeds 1; a doubled mutual edge would exceed it at d < 1


class TestComputeDiffusionMapIsTheSymmetricNormalSpectrum:
    @pytest.mark.parametrize("n_cells, n_neighbors", [(60, 8), (30, 5)])
    def test_matches_an_explicit_decomposition_of_the_documented_operator(self, n_cells: int, n_neighbors: int) -> None:
        data = _neighbor_data(n_cells, n_neighbors)
        graph = data.uns["neighbors"]
        affinity = _explicit_affinities(graph["indices"], graph["distances"])
        degrees = affinity.sum(axis=1)
        operator = affinity / np.sqrt(np.outer(degrees, degrees))  # D^-1/2 A D^-1/2
        values, vectors = np.linalg.eigh(operator)
        values, vectors = values[::-1], vectors[:, ::-1]
        assert values[0] == pytest.approx(1.0, abs=1e-12)
        gaps = -np.diff(values[: N_COMPONENTS + 2])
        assert gaps.min() > 1e-4, "the test needs a spectrum with distinct leading eigenvalues"

        result = compute_diffusion_map(data, n_components=N_COMPONENTS)

        np.testing.assert_allclose(
            result.uns["diffmap"]["eigenvalues"], values[1 : N_COMPONENTS + 1], rtol=0, atol=1e-10
        )
        coordinates = result.obsm["X_diffmap"]
        assert coordinates.shape == (n_cells, N_COMPONENTS)
        for k in range(N_COMPONENTS):
            assert _directions_agree(coordinates[:, k], vectors[:, k + 1])
            np.testing.assert_allclose(operator @ coordinates[:, k], values[k + 1] * coordinates[:, k], atol=1e-10)

    def test_the_reported_eigenvalues_are_not_shifted(self) -> None:
        """A diagonal 'regularization' of 1e-6 used to be added, moving every reported eigenvalue."""
        data = _neighbor_data()
        graph = data.uns["neighbors"]
        affinity = _explicit_affinities(graph["indices"], graph["distances"])
        degrees = affinity.sum(axis=1)
        values = np.linalg.eigvalsh(affinity / np.sqrt(np.outer(degrees, degrees)))[::-1]

        result = compute_diffusion_map(data, n_components=N_COMPONENTS)

        assert np.max(np.abs(result.uns["diffmap"]["eigenvalues"] - values[1 : N_COMPONENTS + 1])) < 1e-10

    def test_each_component_has_its_largest_entry_positive(self) -> None:
        coordinates = compute_diffusion_map(_neighbor_data(), n_components=N_COMPONENTS).obsm["X_diffmap"]

        pivot = np.argmax(np.abs(coordinates), axis=0)
        assert np.all(coordinates[pivot, np.arange(N_COMPONENTS)] > 0)

    def test_a_stale_neighbor_graph_is_refused(self) -> None:
        data = _neighbor_data(60)
        smaller = SingleCellData(data.X[:50].copy(), uns=dict(data.uns))  # the 60-cell graph, 50 cells

        with pytest.raises(errors.ValidationError, match="60 cells but data has 50"):
            compute_diffusion_map(smaller, n_components=N_COMPONENTS)

    def test_a_graph_without_distances_is_refused(self) -> None:
        data = _neighbor_data(30, 5)
        del data.uns["neighbors"]["distances"]

        with pytest.raises(errors.ValidationError, match="distances"):
            compute_diffusion_map(data, n_components=N_COMPONENTS)


class TestTheSolverIsChosenBySizeNotByAnException:
    def test_arpack_and_the_dense_solver_agree_above_the_dense_limit(self) -> None:
        n = nonlinear_methods._DENSE_EIGEN_LIMIT + 200
        rng = np.random.default_rng(5)
        rows = np.repeat(np.arange(n), 6)
        cols = (rows + rng.integers(1, 40, size=rows.size)) % n  # a banded graph: distinct leading eigenvalues
        weights = rng.uniform(0.5, 1.0, size=rows.size)
        adjacency = sparse.csr_matrix((weights, (rows, cols)), shape=(n, n))
        adjacency = sparse.csr_matrix(adjacency.maximum(adjacency.T))
        degrees = np.asarray(adjacency.sum(axis=1)).ravel()
        scale = sparse.diags(1.0 / np.sqrt(degrees))
        operator = sparse.csr_matrix(scale @ adjacency @ scale)

        values, vectors = nonlinear_methods._leading_eigenpairs(operator, 4)

        expected_values, expected_vectors = np.linalg.eigh(operator.toarray())
        np.testing.assert_allclose(values, expected_values[::-1][:4], atol=1e-8)
        for k in range(4):
            assert _directions_agree(vectors[:, k], expected_vectors[:, ::-1][:, k], tolerance=1e-6)

    def test_a_single_cell_has_a_trivial_decomposition_without_any_exception_handling(self) -> None:
        values, vectors = nonlinear_methods._leading_eigenpairs(sparse.csr_matrix((1, 1)), 1)

        assert values.shape == (1,) and vectors.shape == (1, 1)

    @pytest.mark.parametrize(
        "column, expected",
        [
            ([0.2, -0.9, 0.1], [-0.2, 0.9, -0.1]),  # largest magnitude is negative: flip
            ([-0.5, -0.9, -0.1], [0.5, 0.9, 0.1]),
            ([-0.2, 0.9, 0.1], [-0.2, 0.9, 0.1]),  # largest magnitude already positive: unchanged
            ([-0.7, 0.7, 0.1], [0.7, -0.7, -0.1]),  # a tie takes the first entry
            ([0.0, 0.0, 0.0], [0.0, 0.0, 0.0]),
        ],
    )
    def test_the_sign_convention(self, column: list[float], expected: list[float]) -> None:
        flipped = nonlinear_methods._canonical_signs(np.array(column)[:, np.newaxis])

        np.testing.assert_array_equal(flipped[:, 0], np.array(expected))
