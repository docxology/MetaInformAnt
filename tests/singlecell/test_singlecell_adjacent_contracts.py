"""Contracts of functions whose docstrings promised more than the code did.

- ``compute_neighbors`` used ``data.X`` for any ``use_rep`` other than ``"X_pca"`` and recorded the one asked for;
- ``compute_dimensionality_metrics`` reported one number as both trustworthiness and continuity;
- ``bbknn_integration`` built no batch-balanced graph (and stored the shape of an adjacency that did not exist);
- ``scanorama_integration`` and ``mnn_integration`` stacked datasets whose genes differ;
- ``combat_integration`` accepted ``covariates`` "to preserve" and ignored them;
- ``detect_doublets`` replaced a failed SVD with raw gene columns;
- ``cell_type_composition`` turned an invalid confidence level into a 95% interval.

Real data and real library calls only; an independent numpy/scipy implementation is the oracle where one is needed.
"""

from __future__ import annotations

import ast
import inspect
import math
import textwrap

import numpy as np
import pandas as pd
import pytest
from scipy import stats as scipy_stats
from scipy.spatial.distance import cdist

from metainformant.core.utils import errors
from metainformant.singlecell.analysis import nonlinear_methods
from metainformant.singlecell.analysis.nonlinear_methods import (
    compute_diffusion_map,
    compute_dimensionality_metrics,
    compute_neighbors,
)
from metainformant.singlecell.celltyping import annotation
from metainformant.singlecell.celltyping.annotation import cell_type_composition
from metainformant.singlecell.data.integration import (
    bbknn_integration,
    combat_integration,
    mnn_integration,
    scanorama_integration,
)
from metainformant.singlecell.data.preprocessing import SingleCellData
from metainformant.singlecell.doublet import detection
from metainformant.singlecell.doublet.detection import detect_doublets


def _brute_force_neighbors(embedding: np.ndarray, k: int) -> np.ndarray:
    """Sorted index sets of the k nearest rows (itself included, at distance 0) of every row."""
    distances = cdist(embedding, embedding)
    return np.sort(np.argsort(distances, axis=1)[:, :k], axis=1)


class TestComputeNeighborsRepresentation:
    @staticmethod
    def _data() -> SingleCellData:
        rng = np.random.default_rng(0)
        data = SingleCellData(rng.standard_normal((40, 6)))
        data.obsm["X_pca"] = rng.standard_normal((40, 5))
        data.obsm["X_umap"] = rng.standard_normal((40, 2))
        return data

    def test_any_obsm_key_is_used_not_x(self) -> None:
        data = self._data()

        result = compute_neighbors(data, n_neighbors=6, use_rep="X_umap")

        graph = result.uns["neighbors"]
        np.testing.assert_array_equal(np.sort(graph["indices"], axis=1), _brute_force_neighbors(data.obsm["X_umap"], 6))
        assert graph["params"]["use_rep"] == "X_umap"
        assert not np.array_equal(
            np.sort(graph["indices"], axis=1), _brute_force_neighbors(data.X, 6)
        )  # it really is another space

    def test_x_pca_is_used_when_present_and_truncated_to_n_pcs(self) -> None:
        data = self._data()

        result = compute_neighbors(data, n_neighbors=6, n_pcs=3)

        np.testing.assert_array_equal(
            np.sort(result.uns["neighbors"]["indices"], axis=1), _brute_force_neighbors(data.obsm["X_pca"][:, :3], 6)
        )
        assert result.uns["neighbors"]["params"]["use_rep"] == "X_pca"

    def test_explicit_x_ignores_the_embeddings(self) -> None:
        data = self._data()

        result = compute_neighbors(data, n_neighbors=6, use_rep="X")

        np.testing.assert_array_equal(
            np.sort(result.uns["neighbors"]["indices"], axis=1), _brute_force_neighbors(data.X, 6)
        )
        assert result.uns["neighbors"]["params"]["use_rep"] == "X"

    def test_the_default_falls_back_to_x_without_a_pca_and_says_so(self) -> None:
        data = SingleCellData(np.random.default_rng(1).standard_normal((30, 5)))

        params = compute_neighbors(data, n_neighbors=5).uns["neighbors"]["params"]

        assert params["use_rep"] == "X"  # what was used
        assert params["use_rep_requested"] == "X_pca"  # what was asked for

    @pytest.mark.parametrize("name", ["X_tsne", "pca", "nope"])
    def test_an_unknown_representation_is_refused_not_replaced_by_x(self, name: str) -> None:
        with pytest.raises(errors.ValidationError, match=f"'{name}' not found") as raised:
            compute_neighbors(self._data(), n_neighbors=6, use_rep=name)

        assert isinstance(raised.value, ValueError)
        assert "X_umap" in str(raised.value)  # the message lists what exists

    def test_the_diffusion_map_can_be_built_on_a_chosen_embedding(self) -> None:
        data = compute_neighbors(self._data(), n_neighbors=8, use_rep="X_umap")

        result = compute_diffusion_map(data, n_components=3)

        assert result.obsm["X_diffmap"].shape == (40, 3)


class TestTrustworthinessAndContinuity:
    @staticmethod
    def _venna_kaski(original: np.ndarray, embedding: np.ndarray, k: int) -> tuple[float, float]:
        """Trustworthiness and continuity straight from the definition (rank penalties), by brute force."""
        n = len(original)

        def penalty(space_a: np.ndarray, space_b: np.ndarray) -> float:
            """Sum over cells of (rank in A - k) for the cells that are k-neighbors in B but not in A."""
            dist_a, dist_b = cdist(space_a, space_a), cdist(space_b, space_b)
            total = 0.0
            for i in range(n):
                order_a = [j for j in np.argsort(dist_a[i], kind="stable") if j != i]
                order_b = [j for j in np.argsort(dist_b[i], kind="stable") if j != i]
                rank_a = {j: r + 1 for r, j in enumerate(order_a)}
                total += sum(rank_a[j] - k for j in order_b[:k] if rank_a[j] > k)
            return 1.0 - 2.0 / (n * k * (2 * n - 3 * k - 1)) * total

        return penalty(original, embedding), penalty(embedding, original)

    def test_matches_the_definition_and_the_two_measures_differ(self) -> None:
        rng = np.random.default_rng(2)
        original = rng.standard_normal((30, 8))
        embedding = original[:, :2] + 0.3 * rng.standard_normal((30, 2))

        found = nonlinear_methods._compute_trustworthiness_continuity(original, embedding, k=5)

        trust, continuity = self._venna_kaski(original, embedding, 5)
        assert found["trustworthiness"] == pytest.approx(trust, rel=1e-12)
        assert found["continuity"] == pytest.approx(continuity, rel=1e-12)
        assert found["neighborhood_size"] == 5
        # One number reported under two names was the old behavior; the measures are different quantities
        assert abs(found["trustworthiness"] - found["continuity"]) > 1e-3

    def test_an_embedding_that_keeps_every_neighborhood_scores_one(self) -> None:
        original = np.random.default_rng(3).standard_normal((30, 4))

        found = nonlinear_methods._compute_trustworthiness_continuity(original, original.copy(), k=5)

        assert found["trustworthiness"] == pytest.approx(1.0) and found["continuity"] == pytest.approx(1.0)

    def test_k_is_lowered_for_a_small_dataset_and_absent_below_three_cells(self) -> None:
        rng = np.random.default_rng(4)
        small = nonlinear_methods._compute_trustworthiness_continuity(
            rng.standard_normal((8, 3)), rng.standard_normal((8, 2))
        )
        tiny = nonlinear_methods._compute_trustworthiness_continuity(
            rng.standard_normal((2, 3)), rng.standard_normal((2, 2))
        )

        assert small["neighborhood_size"] == 3  # (8 - 1) // 2: the measures need k < n / 2
        assert 0.0 <= small["trustworthiness"] <= 1.0 and 0.0 <= small["continuity"] <= 1.0
        assert math.isnan(tiny["trustworthiness"]) and math.isnan(tiny["continuity"])

    def test_the_public_metrics_report_both(self) -> None:
        rng = np.random.default_rng(5)
        data = SingleCellData(rng.standard_normal((40, 6)))
        data.obs["E1"], data.obs["E2"] = data.X[:, 0] + 0.5 * rng.standard_normal(40), rng.standard_normal(40)

        metrics = compute_dimensionality_metrics(data, ["E1", "E2"])

        assert metrics["neighborhood_size"] == 12
        assert abs(metrics["trustworthiness"] - metrics["continuity"]) > 1e-4

    def test_a_constant_embedding_has_no_distance_correlation(self) -> None:
        data = SingleCellData(np.random.default_rng(6).standard_normal((20, 4)))
        data.obs["E1"] = 1.0
        data.obs["E2"] = 2.0

        metrics = compute_dimensionality_metrics(data, ["E1", "E2"])

        assert math.isnan(metrics["distance_correlation"])  # was 0.0 / a RuntimeWarning


def _batched(sizes: tuple[int, ...] = (30, 12), n_genes: int = 10, seed: int = 0) -> SingleCellData:
    rng = np.random.default_rng(seed)
    X = np.vstack([rng.standard_normal((n, n_genes)) + 2.0 * i for i, n in enumerate(sizes)])
    labels = np.concatenate([[f"b{i}"] * n for i, n in enumerate(sizes)])
    obs = pd.DataFrame({"batch": labels}, index=[f"c{i}" for i in range(len(labels))])
    return SingleCellData(X, obs=obs, var=pd.DataFrame(index=[f"g{i}" for i in range(n_genes)]))


class TestBbknnBuildsABatchBalancedGraph:
    def test_every_cell_has_the_nearest_cells_of_each_batch(self) -> None:
        data = _batched((30, 12))

        result = bbknn_integration(data, batch_key="batch", n_neighbors=5)

        graph = result.uns["bbknn"]
        labels = data.obs["batch"].to_numpy()
        assert graph["indices"].shape == (42, 10)  # 5 from each of the two batches
        for row in graph["indices"]:
            assert (labels[row] == "b0").sum() == 5 and (labels[row] == "b1").sum() == 5
        assert (graph["distances"] >= 0).all()

    def test_a_batch_smaller_than_k_contributes_all_its_cells(self) -> None:
        data = _batched((30, 4))

        graph = bbknn_integration(data, batch_key="batch", n_neighbors=6).uns["bbknn"]

        labels = data.obs["batch"].to_numpy()
        assert graph["indices"].shape == (34, 10)  # 6 from b0, all 4 of b1
        assert all((labels[row] == "b1").sum() == 4 for row in graph["indices"])

    def test_the_neighbors_within_a_batch_are_exact(self) -> None:
        data = _batched((30, 12))
        result = bbknn_integration(data, batch_key="batch", n_neighbors=5)
        # The graph lives in the PCA space of the standardized data; with all 10 components kept PCA is a rotation,
        # which leaves the distances between the standardized rows unchanged
        scaled = (data.X - data.X.mean(axis=0)) / data.X.std(axis=0)
        graph = result.uns["bbknn"]
        labels = data.obs["batch"].to_numpy()

        distances = cdist(scaled, scaled)
        for batch in ("b0", "b1"):
            members = np.flatnonzero(labels == batch)
            for cell in (0, 17, 35):
                expected = set(members[np.argsort(distances[cell, members])[:5]].tolist())
                found = {int(j) for j in graph["indices"][cell] if labels[j] == batch}
                assert found == expected

    def test_the_graph_is_stored_where_the_consumers_read_it_and_the_data_is_unchanged(self) -> None:
        data = _batched((30, 12))

        result = bbknn_integration(data, batch_key="batch", n_neighbors=5)

        np.testing.assert_array_equal(result.X, data.X)
        graph = result.uns["bbknn"]
        assert result.uns["neighbors"]["params"]["method"] == "bbknn"
        np.testing.assert_array_equal(result.uns["neighbors"]["indices"], graph["indices"])
        assert graph["batch_adjacency_shape"] == graph["connectivities"].shape == (42, 42)
        assert (graph["connectivities"] != graph["connectivities"].T).nnz == 0
        diffusion = compute_diffusion_map(result, n_components=3)
        assert diffusion.obsm["X_diffmap"].shape == (42, 3)

    def test_a_missing_batch_key_is_refused(self) -> None:
        with pytest.raises(errors.ValidationError, match="not found"):
            bbknn_integration(_batched(), batch_key="missing")


class TestGenesMustAlign:
    @staticmethod
    def _pair(permute_names: bool = False, extra_gene: bool = False) -> list[SingleCellData]:
        first, second = _batched((25,), seed=1), _batched((25,), seed=2)
        if permute_names:
            second.var = pd.DataFrame(index=list(reversed(first.var.index)))
        if extra_gene:
            second = SingleCellData(
                np.hstack([second.X, np.zeros((25, 1))]), var=pd.DataFrame(index=[*second.var.index, "extra"])
            )
        return [first, second]

    @pytest.mark.parametrize("integrate", [scanorama_integration, mnn_integration])
    def test_datasets_with_other_gene_names_are_refused(self, integrate) -> None:
        with pytest.raises(errors.ValidationError, match="other gene names"):
            integrate(self._pair(permute_names=True), batch_key="batch")

    @pytest.mark.parametrize("integrate", [scanorama_integration, mnn_integration])
    def test_datasets_with_another_number_of_genes_are_refused(self, integrate) -> None:
        with pytest.raises(errors.ValidationError, match="11 genes but data_list\\[0\\] has 10"):
            integrate(self._pair(extra_gene=True), batch_key="batch")

    @pytest.mark.parametrize("integrate", [scanorama_integration, mnn_integration])
    def test_aligned_datasets_are_integrated_and_keep_every_cell(self, integrate) -> None:
        result = integrate(self._pair(), batch_key="batch")

        assert result.n_obs == 50
        assert set(result.obs["batch"]) == {"batch_0", "batch_1"}


class TestCombatDoesNotPretendToPreserveCovariates:

    def test_a_single_cell_reference_batch_is_a_clear_error_not_a_key_error(self) -> None:
        data = _batched((1, 20))

        with pytest.raises(errors.ValidationError, match="fewer than 2 cells"):
            combat_integration(data, batch_key="batch")

    def test_a_single_cell_batch_is_left_uncorrected(self) -> None:
        data = _batched((20, 1))

        result = combat_integration(data, batch_key="batch")

        # The reference batch is the baseline and the lone cell of the other batch has no estimate: nothing moves
        np.testing.assert_array_equal(result.X, data.X)


class TestDoubletDetectionInputs:
    @staticmethod
    def _counts() -> np.ndarray:
        return np.random.default_rng(0).poisson(3.0, (40, 25)).astype(float)

    def test_non_finite_counts_are_refused_before_the_svd(self) -> None:
        counts = self._counts()
        counts[3, 4] = np.nan

        with pytest.raises(errors.ValidationError, match="finite"):
            detect_doublets(counts)

    def test_negative_counts_are_refused(self) -> None:
        counts = self._counts()
        counts[0, 0] = -1.0

        with pytest.raises(errors.ValidationError, match="non-negative"):
            detect_doublets(counts)

    @pytest.mark.parametrize("shape", [(5,), (1, 10), (10, 0)])
    def test_wrong_shapes_are_refused(self, shape: tuple[int, ...]) -> None:
        with pytest.raises(errors.ValidationError, match="2D array"):
            detect_doublets(np.ones(shape))

    def test_no_exception_is_swallowed_around_the_svd(self) -> None:
        """A failed SVD used to be replaced by the first raw gene columns, which are not principal components."""
        tree = ast.parse(textwrap.dedent(inspect.getsource(detection.detect_doublets)))

        assert not [node for node in ast.walk(tree) if isinstance(node, ast.ExceptHandler)]

    def test_a_valid_matrix_still_works_and_is_reproducible(self) -> None:
        first = detect_doublets(self._counts(), random_state=3)
        again = detect_doublets(self._counts(), random_state=3)

        assert len(first.scores) == 40
        np.testing.assert_array_equal(first.scores, again.scores)


class TestConfidenceLevels:
    @pytest.mark.parametrize("level", [0.5, 0.8, 0.9, 0.95, 0.99, 0.999, 0.123])
    def test_the_z_score_is_the_exact_normal_quantile(self, level: float) -> None:
        assert annotation._z_score_for_confidence(level) == pytest.approx(
            float(scipy_stats.norm.ppf((1.0 + level) / 2.0)), rel=1e-12
        )

    @pytest.mark.parametrize("level", [0.0, 1.0, 1.5, -0.1, float("nan"), float("inf")])
    def test_an_invalid_level_is_refused_not_turned_into_95_percent(self, level: float) -> None:
        with pytest.raises(ValueError, match="strictly between 0 and 1"):
            annotation._z_score_for_confidence(level)

        with pytest.raises(ValueError, match="strictly between 0 and 1"):
            cell_type_composition(["a", "b", "a"], confidence_level=level)

    def test_the_wilson_interval_of_no_trials_is_undefined_not_zero(self) -> None:
        lower, upper = annotation._wilson_interval(0, 0, 1.96)

        assert math.isnan(lower) and math.isnan(upper)

    def test_intervals_widen_with_the_confidence_level(self) -> None:
        labels = ["a"] * 30 + ["b"] * 70

        narrow = cell_type_composition(labels, confidence_level=0.8)["confidence_intervals"]["all::a"]
        wide = cell_type_composition(labels, confidence_level=0.99)["confidence_intervals"]["all::a"]

        assert wide[0] < narrow[0] < 0.3 < narrow[1] < wide[1]
