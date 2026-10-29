"""Graph clustering honours ``use_weights`` and marker genes never get an invented p-value.

``leiden_clustering`` recorded ``use_weights`` but never handed the edge weights to leidenalg;
``louvain_clustering`` ignored ``use_weights=False`` (python-louvain reads the ``weight`` attribute whatever the
flag says); ``find_marker_genes`` turned a group without variance into a t statistic of +-1000 with p = 1e-10, an
exception into p = 1 and approximated the t distribution with a normal. leidenalg, python-igraph and
python-louvain are optional: the tests that need them are skipped without them, and the part that decides what
reaches them is tested through the pure helpers that build their input. Real data and real library calls only.
"""

from __future__ import annotations

import math
import warnings

import numpy as np
import pandas as pd
import pytest
from scipy import stats as scipy_stats
from scipy.spatial.distance import pdist, squareform

from metainformant.singlecell.analysis import clustering
from metainformant.singlecell.analysis.clustering import find_marker_genes, leiden_clustering, louvain_clustering
from metainformant.singlecell.data.preprocessing import SingleCellData


def _points(n_cells: int = 40, n_genes: int = 4, seed: int = 0) -> np.ndarray:
    return np.random.default_rng(seed).standard_normal((n_cells, n_genes))


def _brute_force_knn_edges(X: np.ndarray, k: int) -> dict[tuple[int, int], float]:
    """Undirected kNN edges (i < j) with weight 1 / (distance + 1e-10), by an independent construction."""
    distances = squareform(pdist(X))
    edges: dict[tuple[int, int], float] = {}
    for i in range(len(X)):
        nearest = [j for j in np.argsort(distances[i]) if j != i][:k]
        for j in nearest:
            edges[(min(i, int(j)), max(i, int(j)))] = 1.0 / (distances[i, j] + 1e-10)
    return edges


class TestLeidenReceivesTheEdgeWeights:
    def test_the_weights_are_the_inverse_distances_aligned_with_the_edges(self) -> None:
        X = _points()

        problem = clustering._leiden_problem(X, resolution=0.8, n_neighbors=5, random_state=3, use_weights=True)

        expected = _brute_force_knn_edges(X, 5)
        assert problem.weights is not None
        assert dict(zip(problem.edges, problem.weights)) == pytest.approx(expected, rel=1e-12)
        assert len(problem.edges) == len(set(problem.edges))  # each undirected edge once
        assert all(i < j for i, j in problem.edges)
        assert len(set(np.round(problem.weights, 8))) > 1  # real weights, not a constant
        assert (problem.n_vertices, problem.resolution, problem.seed) == (40, 0.8, 3)

    def test_without_use_weights_no_weights_are_passed(self) -> None:
        X = _points()

        weighted = clustering._leiden_problem(X, resolution=1.0, n_neighbors=5, random_state=0, use_weights=True)
        unweighted = clustering._leiden_problem(X, resolution=1.0, n_neighbors=5, random_state=0, use_weights=False)

        assert unweighted.weights is None
        assert unweighted.edges == weighted.edges

    def test_the_default_of_the_public_function_is_weighted(self) -> None:
        import inspect

        assert inspect.signature(leiden_clustering).parameters["use_weights"].default is True

    @pytest.mark.skipif(
        not (clustering.HAS_NETWORKX and clustering.HAS_IGRAPH and clustering.HAS_LEIDEN),
        reason="needs networkx, python-igraph and leidenalg",
    )
    def test_leidenalg_gets_exactly_the_problem_the_helper_describes(self) -> None:
        """With the Leiden stack installed, the public call equals a direct leidenalg call on the helper's graph."""
        X = _points(60, 5, seed=2)

        for use_weights in (True, False):
            problem = clustering._leiden_problem(
                X, resolution=1.0, n_neighbors=8, random_state=11, use_weights=use_weights
            )
            direct = clustering.leidenalg.find_partition(
                clustering.ig.Graph(n=problem.n_vertices, edges=problem.edges, directed=False),
                clustering.leidenalg.RBConfigurationVertexPartition,
                weights=problem.weights,
                resolution_parameter=problem.resolution,
                seed=problem.seed,
            )

            labels = leiden_clustering(
                feature_matrix=X, resolution=1.0, n_neighbors=8, random_state=11, use_weights=use_weights
            )

            np.testing.assert_array_equal(labels, np.array(direct.membership))


class TestLouvainHonoursUseWeights:
    def test_with_weights_the_edge_attribute_is_the_inverse_distance(self) -> None:
        X = _points()

        graph = clustering._louvain_graph(X, 5, use_weights=True)

        expected = _brute_force_knn_edges(X, 5)
        found = {(min(u, v), max(u, v)): data["weight"] for u, v, data in graph.edges(data=True)}
        assert found == pytest.approx(expected, rel=1e-12)

    def test_without_weights_every_edge_weighs_one(self) -> None:
        """python-louvain reads ``weight`` whatever the flag was, so the graph itself must carry unit weights."""
        X = _points()

        graph = clustering._louvain_graph(X, 5, use_weights=False)

        weights = [data["weight"] for _, _, data in graph.edges(data=True)]
        assert weights and set(weights) == {1.0}
        assert {(min(u, v), max(u, v)) for u, v in graph.edges()} == set(_brute_force_knn_edges(X, 5))

    @pytest.mark.skipif(
        not (clustering.HAS_NETWORKX and clustering.HAS_LOUVAIN), reason="needs networkx and python-louvain"
    )
    def test_python_louvain_partitions_the_graph_the_helper_builds(self) -> None:
        X = _points(60, 5, seed=2)

        for use_weights in (True, False):
            graph = clustering._louvain_graph(X, 8, use_weights=use_weights)
            direct = clustering.community_louvain.best_partition(graph, resolution=1.0, random_state=11)

            result = louvain_clustering(
                SingleCellData(X), resolution=1.0, n_neighbors=8, random_state=11, use_weights=use_weights
            )

            np.testing.assert_array_equal(result.obs["louvain_cluster"].to_numpy(), [direct[i] for i in range(60)])


def _two_group_data(shifts: list[float], n_per_group: int = 20, seed: int = 0) -> SingleCellData:
    """Cells x genes; group B is shifted by ``shifts[g]`` in gene g, on unit-variance noise around a mean of 10."""
    rng = np.random.default_rng(seed)
    X = 10.0 + rng.standard_normal((2 * n_per_group, len(shifts)))
    X[n_per_group:] += np.array(shifts)
    obs = pd.DataFrame({"group": ["A"] * n_per_group + ["B"] * n_per_group}, index=[f"c{i}" for i in range(len(X))])
    return SingleCellData(X, obs=obs, var=pd.DataFrame(index=[f"g{i}" for i in range(len(shifts))]))


def _row(markers: pd.DataFrame, gene: str, group: str) -> pd.Series:
    return markers[(markers["gene"] == gene) & (markers["group"] == group)].iloc[0]


class TestMarkerGeneTTest:
    def test_p_values_and_statistics_are_welchs_exact_test(self) -> None:
        data = _two_group_data([0.0, 0.5, 1.5, 3.0])

        markers = find_marker_genes(data, groupby="group", method="t-test", n_genes=10)

        in_group = data.X[20:]
        out_group = data.X[:20]
        for gene_index in range(4):
            expected = scipy_stats.ttest_ind(in_group[:, gene_index], out_group[:, gene_index], equal_var=False)
            row = _row(markers, f"g{gene_index}", "B")
            assert row["statistic"] == pytest.approx(float(expected.statistic), rel=1e-10)
            assert row["p_value"] == pytest.approx(float(expected.pvalue), rel=1e-9)
            assert row["test_status"] == "ok"

    def test_a_group_without_variance_is_tested_not_given_a_statistic_of_1000(self) -> None:
        """Only one group is constant: Welch's standard error comes from the other group, the test is defined."""
        data = _two_group_data([0.0, 0.0])
        data.X[20:, 0] = 2.0  # gene g0 is constant in group B, varied in group A

        markers = find_marker_genes(data, groupby="group", method="t-test", n_genes=10)

        row = _row(markers, "g0", "B")
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")  # scipy warns about precision for constant input
            expected = scipy_stats.ttest_ind(data.X[20:, 0], data.X[:20, 0], equal_var=False)
        assert abs(row["statistic"]) != 1000
        assert row["p_value"] != 1e-10
        assert row["statistic"] == pytest.approx(float(expected.statistic), rel=1e-9)
        assert row["p_value"] == pytest.approx(float(expected.pvalue), rel=1e-9)

    def test_a_gene_constant_in_both_groups_has_no_p_value_not_one(self) -> None:
        data = _two_group_data([0.0, 0.0])
        data.X[:, 1] = 3.0

        markers = find_marker_genes(data, groupby="group", method="t-test", n_genes=10)

        for group in ("A", "B"):
            row = _row(markers, "g1", group)
            assert math.isnan(row["p_value"]) and math.isnan(row["statistic"])
            assert row["test_status"].startswith("zero_variance")

    def test_a_group_of_one_cell_has_no_t_test(self) -> None:
        data = _two_group_data([1.0, 1.0], n_per_group=10)
        data.obs["group"] = ["A"] + ["B"] * 19

        markers = find_marker_genes(data, groupby="group", method="t-test", n_genes=10)

        for gene in ("g0", "g1"):
            row = _row(markers, gene, "A")
            assert math.isnan(row["p_value"])
            assert row["test_status"].startswith("insufficient_observations")
        # the other group is compared with the one-cell remainder, so it has none either
        assert _row(markers, "g0", "B")["test_status"].startswith("insufficient_observations")

    def test_a_single_group_has_nothing_to_compare_with(self) -> None:
        data = _two_group_data([1.0, 1.0])
        data.obs["group"] = ["A"] * data.n_obs

        for method in ("t-test", "wilcoxon"):
            markers = find_marker_genes(data, groupby="group", method=method, n_genes=10)

            assert len(markers) == 2
            assert markers["p_value"].isna().all()
            assert not markers["test_status"].eq("ok").any()

    def test_genes_without_a_p_value_are_ranked_after_the_others(self) -> None:
        data = _two_group_data([2.0, 0.0, 1.0])
        data.X[:, 1] = 5.0  # g1 cannot be tested

        markers = find_marker_genes(data, groupby="group", method="t-test", n_genes=10)

        group_b = markers[markers["group"] == "B"]
        assert list(group_b["gene"])[-1] == "g1"
        assert group_b["p_value"].iloc[:-1].notna().all()

    def test_p_values_that_underflow_to_zero_are_ranked_by_the_size_of_the_statistic(self) -> None:
        """With huge separations every p-value is 0.0; ranking them by gene order would pick arbitrary genes."""
        data = _two_group_data([20.0, 40.0, 60.0, 80.0, 100.0, 0.0], n_per_group=200)

        markers = find_marker_genes(data, groupby="group", method="t-test", n_genes=2)

        group_b = markers[markers["group"] == "B"]
        assert (group_b["p_value"] == 0.0).all()
        assert list(group_b["gene"]) == ["g4", "g3"]  # the two largest statistics, not g0 and g1


class TestMarkerGeneWilcoxon:
    def test_p_values_are_the_rank_sum_test(self) -> None:
        data = _two_group_data([0.0, 1.0, 3.0])

        markers = find_marker_genes(data, groupby="group", method="wilcoxon", n_genes=10)

        for gene_index in range(3):
            expected = scipy_stats.ranksums(data.X[20:, gene_index], data.X[:20, gene_index])
            row = _row(markers, f"g{gene_index}", "B")
            assert row["statistic"] == pytest.approx(float(expected.statistic))
            assert row["p_value"] == pytest.approx(float(expected.pvalue))
            assert row["test_status"] == "ok"

    def test_all_tied_values_have_no_p_value(self) -> None:
        data = _two_group_data([1.0, 1.0])
        data.X[:, 0] = 0.0  # a gene that is zero in every cell

        markers = find_marker_genes(data, groupby="group", method="wilcoxon", n_genes=10)

        for group in ("A", "B"):
            row = _row(markers, "g0", group)
            assert math.isnan(row["p_value"])
            assert row["test_status"].startswith("zero_variance")


class TestMarkerGeneInputs:
    def test_an_unsupported_method_is_refused_even_with_no_genes(self) -> None:
        from metainformant.core.utils.errors import ValidationError

        data = SingleCellData(np.zeros((6, 0)), obs=pd.DataFrame({"group": list("AAABBB")}))

        with pytest.raises(ValidationError, match="Unsupported method"):
            find_marker_genes(data, groupby="group", method="bogus")

    def test_a_negative_ratio_of_means_has_no_log_fold_change_instead_of_zero(self) -> None:
        data = _two_group_data([0.0])
        data.X[:20, 0] = -4.0  # scaled data can have negative means; the ratio of -4 to +4 has no logarithm
        data.X[20:, 0] = np.linspace(3.0, 5.0, 20)

        markers = find_marker_genes(data, groupby="group", method="t-test", n_genes=10)

        assert math.isnan(_row(markers, "g0", "B")["log_fold_change"])
        assert math.isnan(_row(markers, "g0", "A")["log_fold_change"])

    def test_log_fold_change_of_a_gene_absent_from_the_group_is_strongly_negative_not_zero(self) -> None:
        """The ratio of a zero mean to a positive one was set to a log fold change of exactly 0."""
        data = _two_group_data([0.0])
        data.X[:20, 0] = 5.0  # present in A ...
        data.X[20:, 0] = 0.0  # ... absent from B

        markers = find_marker_genes(data, groupby="group", method="wilcoxon", n_genes=10)

        assert _row(markers, "g0", "B")["log_fold_change"] < -20.0
        assert _row(markers, "g0", "A")["log_fold_change"] > 20.0
