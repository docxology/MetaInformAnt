"""Diffusion pseudotime (Haghverdi et al., 2016) follows the cell graph, not the straight line to the root.

The pseudotime used to be the Euclidean distance from the root cell (a shortest path over the complete
distance graph is the direct edge), which turns back on any curved trajectory. These tests pin the real
algorithm down against a curved trajectory, scale invariance, the paper's brute-force definition and
the graph bookkeeping. Real data and real linear algebra only.
"""

from __future__ import annotations

import numpy as np
import pytest
from scipy.stats import spearmanr

from metainformant.core.utils import errors
from metainformant.singlecell.analysis.nonlinear_methods import compute_neighbors
from metainformant.singlecell.analysis.trajectory import compute_diffusion_pseudotime, dpt_trajectory
from metainformant.singlecell.data.preprocessing import SingleCellData


def _arc(n_cells: int = 60, sweep_degrees: float = 300.0) -> tuple[SingleCellData, np.ndarray]:
    """Cells ordered along a 300-degree arc: the true developmental order is the angle."""
    angle = np.radians(np.linspace(0.0, sweep_degrees, n_cells))
    X = np.zeros((n_cells, 5))
    X[:, 0], X[:, 1] = np.cos(angle), np.sin(angle)
    return SingleCellData(X), angle


def _euclidean_from_root(data: SingleCellData, root: int) -> np.ndarray:
    return np.linalg.norm(data.X - data.X[root], axis=1)


class TestCurvedTrajectory:
    def test_pseudotime_follows_the_arc(self) -> None:
        data, angle = _arc()

        pseudotime = compute_diffusion_pseudotime(data, root_cell=0).obs["dpt_pseudotime"].to_numpy()

        assert spearmanr(pseudotime, angle).statistic > 0.99
        assert np.all(np.diff(pseudotime) > 0)  # strictly increasing along the whole arc

    def test_it_is_not_the_straight_line_distance_from_the_root(self) -> None:
        data, angle = _arc()
        euclidean = _euclidean_from_root(data, 0)

        pseudotime = compute_diffusion_pseudotime(data, root_cell=0).obs["dpt_pseudotime"].to_numpy()

        # the straight-line distance rises to the far side of the arc and then comes back toward the root
        assert spearmanr(euclidean, angle).statistic < 0.7
        # compared on the same [0, 1] scale, the two orderings of the cells really are different
        assert np.max(np.abs(pseudotime / pseudotime.max() - euclidean / euclidean.max())) > 0.25

    def test_pseudotime_spans_zero_to_one_from_the_root(self) -> None:
        data, _ = _arc()

        result = compute_diffusion_pseudotime(data, root_cell=0)

        pseudotime = result.obs["dpt_pseudotime"].to_numpy()
        assert pseudotime[0] == pytest.approx(0.0)
        assert pseudotime.max() == pytest.approx(1.0)
        assert np.all(np.isfinite(pseudotime))
        assert result.uns["dpt"]["min_pseudotime"] == pytest.approx(0.0)
        assert result.uns["dpt"]["max_pseudotime"] == pytest.approx(1.0)

    def test_the_far_end_of_the_arc_can_be_the_root(self) -> None:
        data, angle = _arc()

        pseudotime = compute_diffusion_pseudotime(data, root_cell=data.n_obs - 1).obs["dpt_pseudotime"].to_numpy()

        assert pseudotime[-1] == pytest.approx(0.0)
        assert spearmanr(pseudotime, angle).statistic < -0.99


class TestScaleInvariance:
    def test_rescaling_the_expression_does_not_change_the_pseudotime(self) -> None:
        data, _ = _arc()
        scaled = SingleCellData(10.0 * data.X)

        original = compute_diffusion_pseudotime(data, root_cell=0).obs["dpt_pseudotime"].to_numpy()
        rescaled = compute_diffusion_pseudotime(scaled, root_cell=0).obs["dpt_pseudotime"].to_numpy()

        np.testing.assert_allclose(rescaled, original, atol=1e-8)


class TestAgainstTheDefinition:
    def test_matches_the_accumulated_transition_matrix(self) -> None:
        """With every non-trivial eigenvector the result equals the brute-force definition of the paper.

        M = sum_{t>=1} (T^t - steady state) = (I - T + steady state)^-1 - I, and the DPT distance between
        two cells is the Euclidean distance between their rows of M.
        """
        data, _ = _arc(n_cells=40)
        root = 3

        result = compute_diffusion_pseudotime(data, root_cell=root, n_components=data.n_obs - 1)

        affinities = result.uns["neighbors"]["connectivities"].toarray()
        degree = affinities.sum(axis=1)
        density_normalised = affinities / np.outer(degree, degree)
        inverse_root = 1.0 / np.sqrt(density_normalised.sum(axis=1))
        transition = density_normalised * np.outer(inverse_root, inverse_root)
        values, vectors = np.linalg.eigh(transition)
        steady_state = np.outer(vectors[:, -1], vectors[:, -1])
        accumulated = np.linalg.inv(np.eye(data.n_obs) - transition + steady_state) - np.eye(data.n_obs)
        expected = np.linalg.norm(accumulated - accumulated[root], axis=1)

        np.testing.assert_allclose(result.obs["dpt_pseudotime"], expected / expected.max(), atol=1e-9)
        assert values[-1] == pytest.approx(1.0)

    def test_fewer_components_change_the_distance_but_keep_the_order(self) -> None:
        data, angle = _arc()

        few = compute_diffusion_pseudotime(data, root_cell=0, n_components=2)
        many = compute_diffusion_pseudotime(data, root_cell=0, n_components=20)

        assert few.uns["dpt"]["n_components"] == 2
        assert many.uns["dpt"]["n_components"] == 20
        assert not np.allclose(few.obs["dpt_pseudotime"], many.obs["dpt_pseudotime"])
        assert spearmanr(few.obs["dpt_pseudotime"], angle).statistic > 0.99


class TestDisconnectedGraph:
    @staticmethod
    def _two_chains() -> SingleCellData:
        """Two chains of 30 cells a thousand units apart: with 3 neighbours they never connect."""
        X = np.zeros((60, 2))
        X[:30, 0] = np.arange(30)
        X[30:, 0] = 1000 + np.arange(30)
        return SingleCellData(X)

    def test_cells_in_other_components_get_infinite_pseudotime(self) -> None:
        result = compute_diffusion_pseudotime(self._two_chains(), root_cell=0, n_neighbors=3)

        pseudotime = result.obs["dpt_pseudotime"].to_numpy()
        assert np.all(np.isfinite(pseudotime[:30]))
        assert np.all(np.isinf(pseudotime[30:]))
        assert result.uns["dpt"]["n_unreachable"] == 30
        assert pseudotime[:30].max() == pytest.approx(1.0)
        assert result.uns["dpt"]["max_pseudotime"] == pytest.approx(1.0)

    def test_the_other_component_is_reachable_from_a_root_inside_it(self) -> None:
        result = compute_diffusion_pseudotime(self._two_chains(), root_cell=59, n_neighbors=3)

        pseudotime = result.obs["dpt_pseudotime"].to_numpy()
        assert np.all(np.isfinite(pseudotime[30:]))
        assert np.all(np.isinf(pseudotime[:30]))


class TestLargeGraphUsesTheSparseSolver:
    def test_matches_a_dense_eigendecomposition_of_the_same_graph(self) -> None:
        """Above 1000 cells the leading eigenpairs come from ARPACK; they must agree with a dense solve."""
        rng = np.random.default_rng(5)
        data = SingleCellData(rng.normal(size=(1100, 3)))

        result = compute_diffusion_pseudotime(data, root_cell=0, n_components=10)

        affinities = result.uns["neighbors"]["connectivities"].toarray()
        degree = affinities.sum(axis=1)
        density_normalised = affinities / np.outer(degree, degree)
        inverse_root = 1.0 / np.sqrt(density_normalised.sum(axis=1))
        values, vectors = np.linalg.eigh(density_normalised * np.outer(inverse_root, inverse_root))
        leading = np.argsort(values)[::-1][1:11]  # the 10 largest after the trivial eigenvalue 1
        weights = values[leading] / (1.0 - values[leading])
        expected = np.linalg.norm((vectors[:, leading] - vectors[0, leading]) * weights, axis=1)

        np.testing.assert_allclose(result.obs["dpt_pseudotime"], expected / expected.max(), atol=1e-6)
        assert result.uns["dpt"]["n_components"] == 10


class TestGraphAndRootBookkeeping:
    def test_an_existing_neighbor_graph_is_used(self) -> None:
        data, _ = _arc()
        with_graph = compute_neighbors(data, n_neighbors=10, use_rep="X")

        from_stored = compute_diffusion_pseudotime(with_graph, root_cell=0)
        built_inside = compute_diffusion_pseudotime(data, root_cell=0, n_neighbors=10)

        assert from_stored.uns["dpt"]["neighbors_source"] == "uns"
        assert built_inside.uns["dpt"]["neighbors_source"] == "X"
        assert from_stored.uns["dpt"]["n_neighbors"] == built_inside.uns["dpt"]["n_neighbors"] == 10
        np.testing.assert_allclose(from_stored.obs["dpt_pseudotime"], built_inside.obs["dpt_pseudotime"])

    def test_a_stale_neighbor_graph_is_rejected(self) -> None:
        data, _ = _arc()
        stale = SingleCellData(data.X[:40], uns={"neighbors": compute_neighbors(data, use_rep="X").uns["neighbors"]})

        with pytest.raises(errors.ValidationError, match="compute_neighbors"):
            compute_diffusion_pseudotime(stale, root_cell=0)

    def test_the_input_is_not_modified(self) -> None:
        data, _ = _arc()

        result = compute_diffusion_pseudotime(data, root_cell=0)

        assert "dpt_pseudotime" not in data.obs.columns
        assert "neighbors" not in data.uns and "dpt" not in data.uns
        assert "neighbors" in result.uns  # the graph that was used is kept on the result

    def test_the_root_that_was_used_is_recorded(self) -> None:
        data, _ = _arc()

        given = compute_diffusion_pseudotime(data, root_cell=7).uns["dpt"]
        guessed = compute_diffusion_pseudotime(data).uns["dpt"]

        assert (given["root_cell"], given["root_selection"]) == (7, "user")
        assert guessed["root_selection"] == "heuristic"
        assert isinstance(guessed["root_cell"], int) and 0 <= guessed["root_cell"] < data.n_obs

    def test_the_heuristic_root_has_pseudotime_zero(self) -> None:
        data, _ = _arc()

        result = compute_diffusion_pseudotime(data)

        assert result.obs["dpt_pseudotime"].to_numpy()[result.uns["dpt"]["root_cell"]] == pytest.approx(0.0)

    @pytest.mark.parametrize(
        "kwargs",
        [{"root_cell": -1}, {"root_cell": 60}, {"n_components": 0}, {"n_neighbors": 0}],
    )
    def test_invalid_arguments_are_rejected(self, kwargs: dict) -> None:
        data, _ = _arc()

        with pytest.raises(Exception):
            compute_diffusion_pseudotime(data, **kwargs)

    def test_a_single_cell_has_pseudotime_zero(self) -> None:
        result = compute_diffusion_pseudotime(SingleCellData(np.ones((1, 3))), root_cell=0)

        assert result.obs["dpt_pseudotime"].tolist() == [0.0]

    def test_dpt_trajectory_is_the_same_computation(self) -> None:
        data, _ = _arc()

        np.testing.assert_allclose(
            dpt_trajectory(data, root_cell=0).obs["dpt_pseudotime"],
            compute_diffusion_pseudotime(data, root_cell=0, n_components=10).obs["dpt_pseudotime"],
        )
