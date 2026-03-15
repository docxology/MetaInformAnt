"""Tests for metainformant.spatial.analysis -- spatial statistics, clustering, deconvolution, neighborhood.

All tests use real implementations (NO mocking policy).
"""

from __future__ import annotations

import numpy as np
import pytest
from scipy import sparse as sp_sparse

from metainformant.spatial.analysis.autocorrelation import (
    GearyCResult,
    GetisOrdResult,
    LocalMoransResult,
    MoransIResult,
    VariogramResult,
    gearys_c,
    getis_ord_g,
    local_morans_i,
    morans_i,
    spatial_variogram,
    spatial_weights_matrix,
)
from metainformant.spatial.analysis.clustering import (
    SpatialClusterResult,
    build_spatial_graph,
    leiden_clustering,
    louvain_clustering,
    spatial_cluster,
)
from metainformant.spatial.analysis.deconvolution import (
    DeconvolutionResult,
    create_reference_profiles,
    deconvolve_spots,
    enrichment_score,
    estimate_cell_fractions,
    nnls_deconvolution,
)
from metainformant.spatial.analysis.neighborhood import (
    InteractionResult,
    NeighborhoodEnrichmentResult,
    NicheResult,
    RipleyKResult,
    compute_interaction_matrix,
    neighborhood_enrichment,
    niche_detection,
    ripley_k,
)

# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture()
def grid_coordinates() -> np.ndarray:
    """10x10 regular grid of (x, y) positions."""
    xs, ys = np.meshgrid(np.arange(10), np.arange(10))
    return np.column_stack([xs.ravel(), ys.ravel()]).astype(np.float64)


@pytest.fixture()
def spatially_autocorrelated_values(grid_coordinates: np.ndarray) -> np.ndarray:
    """Values with strong positive spatial autocorrelation (gradient)."""
    return grid_coordinates[:, 0] + grid_coordinates[:, 1]


@pytest.fixture()
def random_values() -> np.ndarray:
    """Spatially random values."""
    rng = np.random.RandomState(42)
    return rng.standard_normal(100)


@pytest.fixture()
def expression_matrix() -> np.ndarray:
    """Synthetic expression matrix (50 spots x 20 genes)."""
    rng = np.random.RandomState(42)
    return np.abs(rng.standard_normal((50, 20)))


@pytest.fixture()
def small_coordinates() -> np.ndarray:
    """50 random 2D positions."""
    rng = np.random.RandomState(42)
    return rng.uniform(0, 100, (50, 2))


# ---------------------------------------------------------------------------
# Spatial weights
# ---------------------------------------------------------------------------


class TestSpatialWeightsMatrix:
    def test_knn_weights(self, grid_coordinates: np.ndarray) -> None:
        W = spatial_weights_matrix(grid_coordinates, method="knn", k=4)
        assert W.shape == (100, 100)
        assert sp_sparse.issparse(W)
        assert W.nnz > 0

    def test_distance_weights(self, grid_coordinates: np.ndarray) -> None:
        W = spatial_weights_matrix(grid_coordinates, method="distance", bandwidth=1.5)
        assert W.shape == (100, 100)
        assert W.nnz > 0

    def test_binary_weights(self, grid_coordinates: np.ndarray) -> None:
        W = spatial_weights_matrix(grid_coordinates, method="binary", bandwidth=1.5)
        assert W.shape == (100, 100)
        assert W.nnz > 0

    def test_row_standardization(self, grid_coordinates: np.ndarray) -> None:
        W = spatial_weights_matrix(grid_coordinates, method="knn", k=4, row_standardize=True)
        row_sums = np.array(W.sum(axis=1)).flatten()
        # Non-zero rows should sum to approximately 1
        nonzero_rows = row_sums > 0
        np.testing.assert_allclose(row_sums[nonzero_rows], 1.0, atol=1e-10)


# ---------------------------------------------------------------------------
# Moran's I
# ---------------------------------------------------------------------------


class TestMoransI:
    def test_positive_autocorrelation(
        self,
        spatially_autocorrelated_values: np.ndarray,
        grid_coordinates: np.ndarray,
    ) -> None:
        W = spatial_weights_matrix(grid_coordinates, method="knn", k=4)
        result = morans_i(spatially_autocorrelated_values, W)
        assert isinstance(result, MoransIResult)
        assert result.I > 0.0, "Spatially autocorrelated data should yield positive I"
        assert result.n == 100

    def test_random_values_near_zero(
        self,
        random_values: np.ndarray,
        grid_coordinates: np.ndarray,
    ) -> None:
        W = spatial_weights_matrix(grid_coordinates, method="knn", k=4)
        result = morans_i(random_values, W)
        assert isinstance(result, MoransIResult)
        # Random data: I should be near the expected value (close to 0)
        assert abs(result.I) < 0.5

    def test_has_p_value(
        self,
        spatially_autocorrelated_values: np.ndarray,
        grid_coordinates: np.ndarray,
    ) -> None:
        W = spatial_weights_matrix(grid_coordinates, method="knn", k=4)
        result = morans_i(spatially_autocorrelated_values, W)
        assert 0.0 <= result.p_value <= 1.0


# ---------------------------------------------------------------------------
# Geary's C
# ---------------------------------------------------------------------------


class TestGearysC:
    def test_basic_operation(
        self,
        spatially_autocorrelated_values: np.ndarray,
        grid_coordinates: np.ndarray,
    ) -> None:
        W = spatial_weights_matrix(grid_coordinates, method="knn", k=4)
        result = gearys_c(spatially_autocorrelated_values, W)
        assert isinstance(result, GearyCResult)
        # Positive autocorrelation: C should be < 1
        assert result.C < 1.0
        assert result.n == 100
        assert 0.0 <= result.p_value <= 1.0


# ---------------------------------------------------------------------------
# Getis-Ord G*
# ---------------------------------------------------------------------------


class TestGetisOrdG:
    def test_hotspot_detection(self, grid_coordinates: np.ndarray) -> None:
        # Create a clear hotspot in the top-right corner
        values = np.zeros(100)
        values[grid_coordinates[:, 0] > 7] = 10.0
        W = spatial_weights_matrix(grid_coordinates, method="knn", k=4, row_standardize=False)
        result = getis_ord_g(values, W)
        assert isinstance(result, GetisOrdResult)
        assert result.g_star.shape == (100,)
        assert result.hot_spots.dtype == bool
        assert result.cold_spots.dtype == bool


# ---------------------------------------------------------------------------
# Local Moran's I (LISA)
# ---------------------------------------------------------------------------


class TestLocalMoransI:
    def test_identifies_local_clusters(
        self,
        spatially_autocorrelated_values: np.ndarray,
        grid_coordinates: np.ndarray,
    ) -> None:
        W = spatial_weights_matrix(grid_coordinates, method="knn", k=4)
        result = local_morans_i(spatially_autocorrelated_values, W)
        assert isinstance(result, LocalMoransResult)
        assert len(result.local_I) == 100
        assert len(result.cluster_labels) == 100
        valid_labels = {"HH", "LL", "HL", "LH", "NS"}
        assert all(l in valid_labels for l in result.cluster_labels)


# ---------------------------------------------------------------------------
# Variogram
# ---------------------------------------------------------------------------


class TestSpatialVariogram:
    def test_returns_variogram_result(
        self,
        spatially_autocorrelated_values: np.ndarray,
        grid_coordinates: np.ndarray,
    ) -> None:
        result = spatial_variogram(spatially_autocorrelated_values, grid_coordinates, n_bins=10)
        assert isinstance(result, VariogramResult)
        assert len(result.bin_centers) == 10
        assert len(result.semivariance) == 10
        assert result.sill >= 0.0


# ---------------------------------------------------------------------------
# Spatial graph & clustering
# ---------------------------------------------------------------------------


class TestBuildSpatialGraph:
    def test_knn_graph(self, small_coordinates: np.ndarray) -> None:
        adj = build_spatial_graph(small_coordinates, method="knn", n_neighbors=5)
        assert sp_sparse.issparse(adj)
        assert adj.shape == (50, 50)
        assert adj.nnz > 0

    def test_delaunay_graph(self, small_coordinates: np.ndarray) -> None:
        adj = build_spatial_graph(small_coordinates, method="delaunay")
        assert adj.shape == (50, 50)
        assert adj.nnz > 0

    def test_radius_graph(self, small_coordinates: np.ndarray) -> None:
        adj = build_spatial_graph(small_coordinates, method="radius", radius=30.0)
        assert adj.shape == (50, 50)


class TestSpatialCluster:
    def test_kmeans_clustering(self, expression_matrix: np.ndarray, small_coordinates: np.ndarray) -> None:
        result = spatial_cluster(
            expression_matrix,
            small_coordinates,
            n_clusters=3,
            method="kmeans",
            seed=42,
        )
        assert isinstance(result, SpatialClusterResult)
        assert len(result.labels) == 50
        assert result.n_clusters == 3
        assert result.method == "kmeans"

    def test_leiden_clustering(self, expression_matrix: np.ndarray, small_coordinates: np.ndarray) -> None:
        result = spatial_cluster(
            expression_matrix,
            small_coordinates,
            method="leiden",
            resolution=0.5,
            seed=42,
        )
        assert isinstance(result, SpatialClusterResult)
        assert len(result.labels) == 50
        assert result.n_clusters >= 1


class TestLeidenAndLouvain:
    def test_leiden_clustering_direct(self, small_coordinates: np.ndarray) -> None:
        adj = build_spatial_graph(small_coordinates, method="knn", n_neighbors=5)
        labels, modularity = leiden_clustering(adj, resolution=1.0, seed=42)
        assert len(labels) == 50
        assert isinstance(modularity, float)

    def test_louvain_clustering_direct(self, small_coordinates: np.ndarray) -> None:
        adj = build_spatial_graph(small_coordinates, method="knn", n_neighbors=5)
        labels, modularity = louvain_clustering(adj, resolution=1.0, seed=42)
        assert len(labels) == 50
        assert isinstance(modularity, float)


# ---------------------------------------------------------------------------
# Deconvolution
# ---------------------------------------------------------------------------


class TestDeconvolution:
    def test_deconvolve_spots_returns_result(self) -> None:
        rng = np.random.RandomState(42)
        # 3 cell types, 10 genes
        ref = np.abs(rng.standard_normal((3, 10)))
        # 20 spots, each a mixture of cell types
        fracs = rng.dirichlet([1, 1, 1], size=20)
        spatial_expr = fracs @ ref + rng.standard_normal((20, 10)) * 0.01
        spatial_expr = np.abs(spatial_expr)

        result = deconvolve_spots(spatial_expr, ref, method="nnls", cell_type_names=["A", "B", "C"])
        assert isinstance(result, DeconvolutionResult)
        assert result.n_spots == 20
        assert result.n_types == 3
        assert result.fractions.shape == (20, 3)

    def test_nnls_deconvolution_single_spot(self) -> None:
        rng = np.random.RandomState(42)
        ref = np.abs(rng.standard_normal((3, 10)))
        spot = ref[0, :] * 0.5 + ref[1, :] * 0.3 + ref[2, :] * 0.2
        weights, residual = nnls_deconvolution(spot, ref)
        assert weights.shape == (3,)
        assert isinstance(residual, float)

    def test_nnls_deconvolution_batch(self) -> None:
        rng = np.random.RandomState(42)
        ref = np.abs(rng.standard_normal((3, 10)))
        bulk = np.abs(rng.standard_normal((5, 10)))
        weights, residuals = nnls_deconvolution(bulk, ref)
        assert weights.shape == (5, 3)
        assert residuals.shape == (5,)

    def test_create_reference_profiles(self) -> None:
        rng = np.random.RandomState(42)
        data = np.abs(rng.standard_normal((30, 10)))
        labels = np.array(["A"] * 10 + ["B"] * 10 + ["C"] * 10)
        profiles, type_names, gene_names = create_reference_profiles(data, labels)
        assert profiles.shape == (3, 10)
        assert type_names == ["A", "B", "C"]
        assert len(gene_names) == 10

    def test_estimate_cell_fractions(self) -> None:
        weights = np.array([[3.0, 1.0, 0.0], [0.0, 2.0, 2.0]])
        fracs = estimate_cell_fractions(weights)
        np.testing.assert_allclose(fracs.sum(axis=1), 1.0)

    def test_enrichment_score(self) -> None:
        observed = np.array([0.5, 0.3, 0.2])
        expected = np.array([0.33, 0.33, 0.33])
        scores = enrichment_score(observed, expected)
        assert scores.shape == (3,)


# ---------------------------------------------------------------------------
# Neighborhood
# ---------------------------------------------------------------------------


class TestNeighborhoodAnalysis:
    def test_compute_interaction_matrix(self, small_coordinates: np.ndarray) -> None:
        adj = build_spatial_graph(small_coordinates, method="knn", n_neighbors=5)
        rng = np.random.RandomState(42)
        cell_types = rng.choice(["A", "B", "C"], size=50)
        result = compute_interaction_matrix(cell_types, adj)
        assert isinstance(result, InteractionResult)
        assert result.interaction_matrix.shape == (3, 3)

    def test_neighborhood_enrichment(self, small_coordinates: np.ndarray) -> None:
        rng = np.random.RandomState(42)
        cell_types = rng.choice(["A", "B"], size=50)
        result = neighborhood_enrichment(
            cell_types,
            small_coordinates,
            n_neighbors=5,
            n_permutations=50,
            seed=42,
        )
        assert isinstance(result, NeighborhoodEnrichmentResult)
        assert result.enrichment_matrix.shape == (2, 2)
        assert result.p_values.shape == (2, 2)

    def test_niche_detection(self, small_coordinates: np.ndarray) -> None:
        rng = np.random.RandomState(42)
        cell_types = rng.choice(["A", "B", "C"], size=50)
        result = niche_detection(
            cell_types,
            small_coordinates,
            n_niches=3,
            n_neighbors=5,
            seed=42,
        )
        assert isinstance(result, NicheResult)
        assert len(result.niche_labels) == 50
        assert result.n_niches == 3
        assert result.niche_compositions.shape == (3, 3)

    def test_ripley_k(self) -> None:
        rng = np.random.RandomState(42)
        points = rng.uniform(0, 100, (30, 2))
        radii = np.linspace(1, 30, 10)
        area = 100.0 * 100.0
        result = ripley_k(points, radii, area, n_simulations=10, seed=42)
        assert isinstance(result, RipleyKResult)
        assert len(result.k_values) == 10
        assert len(result.l_values) == 10
        assert result.n_points == 30
        assert result.area == area
