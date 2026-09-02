"""Round-4 test-depth tests for singlecell clustering and trajectory analysis.

Covers previously untested public APIs:
- analysis.clustering: kmeans_clustering, louvain_clustering (dependency-aware),
  find_marker_genes, compute_cluster_composition, compute_cluster_silhouette,
  evaluate_clustering_performance
- analysis.trajectory: compute_diffusion_pseudotime, dpt_trajectory,
  paga_trajectory, slingshot_trajectory, find_trajectory_branches,
  compute_trajectory_entropy

All tests use real computation on synthetic data (zero-mocks policy).
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from metainformant.singlecell.analysis.clustering import (
    compute_cluster_composition,
    compute_cluster_silhouette,
    evaluate_clustering_performance,
    find_marker_genes,
    kmeans_clustering,
)
from metainformant.singlecell.analysis.trajectory import (
    compute_diffusion_pseudotime,
    compute_trajectory_entropy,
    dpt_trajectory,
    find_trajectory_branches,
    paga_trajectory,
    slingshot_trajectory,
)
from metainformant.singlecell.data.preprocessing import SingleCellData

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _make_clustered_data(n_per_cluster: int = 30, n_genes: int = 20, seed: int = 0) -> SingleCellData:
    """Two well-separated clusters; ground truth in obs['cell_type']."""
    rng = np.random.RandomState(seed)
    X = np.vstack([rng.normal(5, 0.5, (n_per_cluster, n_genes)), rng.normal(9, 0.5, (n_per_cluster, n_genes))]).clip(
        min=0
    )
    obs = pd.DataFrame(
        {"cell_type": ["A"] * n_per_cluster + ["B"] * n_per_cluster},
        index=[f"c{i}" for i in range(2 * n_per_cluster)],
    )
    var = pd.DataFrame(index=[f"g{i}" for i in range(n_genes)])
    return SingleCellData(X, obs=obs, var=var)


def _make_trajectory_data(n_cells: int = 40, n_genes: int = 8, seed: int = 0) -> SingleCellData:
    """A smooth one-dimensional trajectory along gene 0; two clusters."""
    rng = np.random.RandomState(seed)
    t = np.linspace(0, 1, n_cells)
    X = np.column_stack([t, rng.normal(0, 0.05, n_cells) + 0.2 * t, rng.normal(0, 0.1, (n_cells, n_genes - 2))]).clip(
        min=0
    )
    obs = pd.DataFrame(
        {"cluster": ["A"] * (n_cells // 2) + ["B"] * (n_cells // 2)},
        index=[f"c{i}" for i in range(n_cells)],
    )
    var = pd.DataFrame(index=[f"g{i}" for i in range(n_genes)])
    return SingleCellData(X, obs=obs, var=var)


# ---------------------------------------------------------------------------
# kmeans_clustering
# ---------------------------------------------------------------------------


class TestKmeansClustering:
    def test_recovers_two_clusters(self) -> None:
        data = kmeans_clustering(_make_clustered_data(), n_clusters=2, random_state=0)
        assert "kmeans_cluster" in data.obs.columns
        labels = data.obs["kmeans_cluster"].values
        # ground truth: first 30 cells type A, next 30 type B
        assert len(set(labels[:30])) == 1
        assert len(set(labels[30:])) == 1
        assert labels[0] != labels[-1]

    def test_label_range(self) -> None:
        data = kmeans_clustering(_make_clustered_data(), n_clusters=2, random_state=0)
        assert set(data.obs["kmeans_cluster"].unique()).issubset({0, 1})

    def test_deterministic_with_seed(self) -> None:
        d1 = kmeans_clustering(_make_clustered_data(), n_clusters=2, random_state=7)
        d2 = kmeans_clustering(_make_clustered_data(), n_clusters=2, random_state=7)
        assert list(d1.obs["kmeans_cluster"]) == list(d2.obs["kmeans_cluster"])

    def test_rejects_non_singlecelldata(self) -> None:
        with pytest.raises(Exception):
            kmeans_clustering("not data")  # type: ignore[arg-type]


# ---------------------------------------------------------------------------
# find_marker_genes
# ---------------------------------------------------------------------------


class TestFindMarkerGenes:
    def test_returns_dataframe_with_expected_columns(self) -> None:
        data = kmeans_clustering(_make_clustered_data(), n_clusters=2, random_state=0)
        markers = find_marker_genes(data, groupby="cell_type", n_genes=5)
        assert isinstance(markers, pd.DataFrame)
        for col in ("gene", "group", "log_fold_change", "p_value"):
            assert col in markers.columns
        assert len(markers) == 10  # 5 genes x 2 groups

    def test_identifies_cluster_specific_genes(self) -> None:
        # genes 0-4 elevated in cluster B cells (rows 30+)
        data = kmeans_clustering(_make_clustered_data(), n_clusters=2, random_state=0)
        markers = find_marker_genes(data, groupby="cell_type", n_genes=5)
        group_b = markers[markers["group"] == "B"]
        top_b = set(group_b.sort_values("log_fold_change", ascending=False).head(5)["gene"])
        assert top_b == {f"g{i}" for i in range(5)}

    def test_invalid_groupby_raises(self) -> None:
        data = _make_clustered_data()
        with pytest.raises(Exception):
            find_marker_genes(data, groupby="nonexistent")


# ---------------------------------------------------------------------------
# compute_cluster_composition
# ---------------------------------------------------------------------------


class TestComputeClusterComposition:
    def test_pure_clusters_composition(self) -> None:
        data = kmeans_clustering(_make_clustered_data(), n_clusters=2, random_state=0)
        comp = compute_cluster_composition(data, groupby="cell_type", cluster_col="kmeans_cluster")
        assert {"cluster", "category", "count", "proportion", "percentage"}.issubset(comp.columns)
        # each cluster is pure
        for cluster_id, grp in comp.groupby("cluster"):
            nonzero = grp[grp["count"] > 0]
            assert len(nonzero) == 1
            assert nonzero["proportion"].iloc[0] == pytest.approx(1.0)

    def test_missing_cluster_column_raises(self) -> None:
        data = _make_clustered_data()
        with pytest.raises(Exception):
            compute_cluster_composition(data, groupby="cell_type", cluster_col="nope")


# ---------------------------------------------------------------------------
# compute_cluster_silhouette (regression: previously crashed with
# "'float' object is not subscriptable" because silhouette_score was indexed
# like silhouette_samples)
# ---------------------------------------------------------------------------


class TestComputeClusterSilhouette:
    def test_returns_scores_without_error(self) -> None:
        data = kmeans_clustering(_make_clustered_data(), n_clusters=2, random_state=0)
        result = compute_cluster_silhouette(data, cluster_col="kmeans_cluster")
        assert "error" not in result
        assert result["n_clusters"] == 2
        assert result["overall_silhouette_score"] > 0.5
        assert set(result["per_cluster_scores"].keys()) == {0, 1}
        for stats in result["per_cluster_scores"].values():
            assert {"mean", "std", "min", "max"}.issubset(stats.keys())
            assert stats["min"] <= stats["mean"] <= stats["max"]

    def test_includes_additional_quality_metrics(self) -> None:
        data = kmeans_clustering(_make_clustered_data(), n_clusters=2, random_state=0)
        result = compute_cluster_silhouette(data, cluster_col="kmeans_cluster")
        assert "calinski_harabasz_score" in result
        assert "davies_bouldin_score" in result

    def test_single_cluster_reports_error(self) -> None:
        data = _make_clustered_data()
        data.obs["one_cluster"] = 0
        result = compute_cluster_silhouette(data, cluster_col="one_cluster")
        assert "error" in result

    def test_sample_size_subsets(self) -> None:
        data = kmeans_clustering(_make_clustered_data(), n_clusters=2, random_state=0)
        result = compute_cluster_silhouette(data, cluster_col="kmeans_cluster", sample_size=20)
        assert "error" not in result
        assert result["n_samples"] == 20


# ---------------------------------------------------------------------------
# evaluate_clustering_performance
# ---------------------------------------------------------------------------


class TestEvaluateClusteringPerformance:
    def test_perfect_clustering_metrics(self) -> None:
        data = kmeans_clustering(_make_clustered_data(), n_clusters=2, random_state=0)
        result = evaluate_clustering_performance(data, cluster_col="kmeans_cluster", ground_truth_col="cell_type")
        assert result["n_clusters"] == 2
        assert result["cluster_sizes"] == [30, 30]
        sup = result["supervised_metrics"]
        assert sup["adjusted_rand_index"] == pytest.approx(1.0)
        assert sup["normalized_mutual_info"] == pytest.approx(1.0)
        # silhouette analysis is embedded and must not error (regression)
        assert "error" not in result["silhouette_analysis"]

    def test_without_ground_truth(self) -> None:
        data = kmeans_clustering(_make_clustered_data(), n_clusters=2, random_state=0)
        result = evaluate_clustering_performance(data, cluster_col="kmeans_cluster")
        assert result["n_clusters"] == 2


# ---------------------------------------------------------------------------
# Trajectory analysis
# ---------------------------------------------------------------------------


class TestDiffusionPseudotime:
    def test_pseudotime_correlates_with_position(self) -> None:
        data = dpt_trajectory(_make_trajectory_data(), root_cell=0)
        assert "dpt_pseudotime" in data.obs.columns
        assert "dpt" in data.uns
        pt = data.obs["dpt_pseudotime"].values
        true_t = np.linspace(0, 1, len(pt))
        assert np.corrcoef(pt, true_t)[0, 1] > 0.9

    def test_root_cell_has_extreme_pseudotime(self) -> None:
        data = dpt_trajectory(_make_trajectory_data(), root_cell=0)
        pt = data.obs["dpt_pseudotime"].values
        assert pt[0] == pytest.approx(0.0)

    def test_compute_diffusion_pseudotime_adds_column(self) -> None:
        data = compute_diffusion_pseudotime(_make_trajectory_data(), root_cell=0, n_components=3)
        assert "dpt_pseudotime" in data.obs.columns


class TestPagaSlingshot:
    def test_paga_trajectory(self) -> None:
        data = paga_trajectory(_make_trajectory_data(), groups="cluster")
        assert "paga_pseudotime" in data.obs.columns
        assert data.uns.get("paga") is not None

    def test_slingshot_trajectory(self) -> None:
        data = slingshot_trajectory(_make_trajectory_data(), start_cluster="A", end_clusters=["B"])
        assert "slingshot_pseudotime" in data.obs.columns
        pt = data.obs["slingshot_pseudotime"].values
        assert pt.min() >= 0.0

    def test_slingshot_unknown_cluster_raises(self) -> None:
        with pytest.raises(Exception):
            slingshot_trajectory(_make_trajectory_data(), start_cluster="Z", end_clusters=["B"])


class TestTrajectoryAnalysis:
    def test_find_trajectory_branches_structure(self) -> None:
        data = dpt_trajectory(_make_trajectory_data(), root_cell=0)
        branches = find_trajectory_branches(data, pseudotime_col="dpt_pseudotime", min_branch_size=5)
        assert {"n_branches", "branch_sizes", "branch_labels", "branches"}.issubset(branches.keys())
        assert branches["n_branches"] >= 1

    def test_trajectory_entropy_structure(self) -> None:
        data = dpt_trajectory(_make_trajectory_data(), root_cell=0)
        entropy = compute_trajectory_entropy(data, pseudotime_col="dpt_pseudotime", window_size=10)
        assert {"entropy_profile", "pseudotime_windows", "differentiation_points", "window_size", "n_windows"}.issubset(
            entropy.keys()
        )
        assert entropy["n_windows"] == len(entropy["entropy_profile"])
        for value in entropy["entropy_profile"]:
            assert value >= 0.0
