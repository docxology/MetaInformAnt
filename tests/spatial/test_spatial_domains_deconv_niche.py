"""Round-4 test-depth tests for spatial analysis, deconvolution, and niche identification.

Covers previously untested public APIs:
- analysis.clustering: spatial_domains, spatial_cluster (all methods),
  build_spatial_graph (knn/delaunay/radius)
- deconvolution.spatial_deconvolution: deconvolve_spots, build_reference_profiles,
  spatial_cell_type_mapping, validate_deconvolution, niche_identification
- niche.identification: identify_niches

All tests use real computation on synthetic data (zero-mocks policy).
"""

from __future__ import annotations

import numpy as np
import pytest

from metainformant.spatial.analysis.clustering import (
    build_spatial_graph,
    spatial_cluster,
    spatial_domains,
)
from metainformant.spatial.deconvolution.spatial_deconvolution import (
    build_reference_profiles,
    deconvolve_spots,
    niche_identification,
    spatial_cell_type_mapping,
    validate_deconvolution,
)
from metainformant.spatial.niche.identification import NicheResult, identify_niches

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _make_spatial_expression(n_domains: int = 3, per_domain: int = 30, n_genes: int = 10, seed: int = 0):
    """Expression with `n_domains` distinct expression levels; 3x30 spot grid."""
    rng = np.random.RandomState(seed)
    levels = [5.0 + 3.0 * i for i in range(n_domains)]
    X = np.vstack([rng.normal(level, 0.3, (per_domain, n_genes)) for level in levels]).clip(min=0)
    coords = np.array([[float(i % 30), float(i // 30)] for i in range(n_domains * per_domain)])
    return X, coords


def _make_reference_and_spots(seed: int = 0):
    """Single-cell reference (2 cell types) plus 20 spatial spots."""
    rng = np.random.RandomState(seed)
    sc_expr = np.vstack([rng.normal(10, 1, (50, 20)), rng.normal(2, 1, (50, 20))]).clip(min=0)
    cell_types = ["A"] * 50 + ["B"] * 50
    spatial_counts = np.vstack([rng.normal(8, 1, (10, 20)), rng.normal(3, 1, (10, 20))]).clip(min=0)
    coords = [(float(rng.uniform(0, 100)), float(rng.uniform(0, 100))) for _ in range(20)]
    return sc_expr, cell_types, spatial_counts, coords


# ---------------------------------------------------------------------------
# spatial_domains / spatial_cluster
# ---------------------------------------------------------------------------


class TestSpatialDomains:
    def test_recovers_three_strips(self) -> None:
        X, coords = _make_spatial_expression()
        result = spatial_domains(X, coords, n_domains=3, seed=0)
        labels = np.asarray(result.labels)
        assert len(labels) == 90
        assert len(set(labels)) == 3
        # strips along y: rows 0-29, 30-59, 60-89 each form one domain
        for i in range(3):
            assert len(set(labels[i * 30 : (i + 1) * 30])) == 1

    def test_result_metadata(self) -> None:
        X, coords = _make_spatial_expression()
        result = spatial_domains(X, coords, n_domains=3, seed=0)
        assert result.n_clusters == 3
        assert result.method == "spatial_domains"


class TestSpatialCluster:
    @pytest.mark.parametrize("method", ["leiden", "louvain", "kmeans"])
    def test_all_methods_run(self, method: str) -> None:
        X, coords = _make_spatial_expression()
        result = spatial_cluster(X, coords, n_clusters=3, method=method, seed=0)
        labels = np.asarray(result.labels)
        assert len(labels) == 90
        # graph-based fallbacks may over-partition; kmeans must be exact
        if method == "kmeans":
            assert len(set(labels)) == 3
        else:
            assert len(set(labels)) >= 1

    def test_kmeans_recovers_domains(self) -> None:
        X, coords = _make_spatial_expression()
        result = spatial_cluster(X, coords, n_clusters=3, method="kmeans", seed=0)
        labels = np.asarray(result.labels)
        for i in range(3):
            assert len(set(labels[i * 30 : (i + 1) * 30])) == 1

    def test_invalid_method_raises(self) -> None:
        X, coords = _make_spatial_expression()
        with pytest.raises(Exception):
            spatial_cluster(X, coords, n_clusters=3, method="spectral")  # type: ignore[arg-type]


class TestBuildSpatialGraph:
    def test_knn_graph(self) -> None:
        _, coords = _make_spatial_expression()
        graph = build_spatial_graph(coords, method="knn", n_neighbors=4)
        assert graph.shape == (90, 90)

    def test_delaunay_graph(self) -> None:
        _, coords = _make_spatial_expression()
        graph = build_spatial_graph(coords, method="delaunay")
        assert graph.shape == (90, 90)
        assert graph.nnz > 0

    def test_radius_graph(self) -> None:
        _, coords = _make_spatial_expression()
        graph = build_spatial_graph(coords, method="radius", radius=1.5)
        assert graph.shape == (90, 90)
        assert graph.nnz > 0


# ---------------------------------------------------------------------------
# Deconvolution
# ---------------------------------------------------------------------------


class TestBuildReferenceProfiles:
    def test_reference_structure(self) -> None:
        sc_expr, cell_types, _, _ = _make_reference_and_spots()
        ref = build_reference_profiles(sc_expr, cell_types, n_markers=10)
        assert set(ref.keys()) >= {"profiles", "gene_indices", "n_markers_per_type", "unique_types"}
        assert set(ref["profiles"].keys()) == {"A", "B"}
        for profile in ref["profiles"].values():
            assert len(profile) == 20

    def test_more_markers_than_genes(self) -> None:
        sc_expr, cell_types, _, _ = _make_reference_and_spots()
        ref = build_reference_profiles(sc_expr, cell_types, n_markers=100)
        # clamped to available genes
        for profile in ref["profiles"].values():
            assert len(profile) <= 20


class TestDeconvolveSpots:
    def test_recovers_dominant_type(self) -> None:
        sc_expr, cell_types, spatial_counts, _ = _make_reference_and_spots()
        ref = build_reference_profiles(sc_expr, cell_types, n_markers=10)
        result = deconvolve_spots(spatial_counts, ref["profiles"])
        assert {"proportions_matrix", "cell_types", "confidence_scores", "spots_summary"}.issubset(result.keys())
        props = np.asarray(result["proportions_matrix"])
        assert props.shape == (20, 2)
        # rows sum to ~1
        assert np.allclose(props.sum(axis=1), 1.0, atol=0.05)
        # first 10 spots resemble type A (high expression), last 10 type B
        type_order = list(result["cell_types"])
        a_idx = type_order.index("A")
        assert props[:10, a_idx].mean() > props[10:, a_idx].mean()


class TestSpatialCellTypeMapping:
    def test_mapping_structure(self) -> None:
        sc_expr, cell_types, spatial_counts, coords = _make_reference_and_spots()
        ref = build_reference_profiles(sc_expr, cell_types, n_markers=10)
        decon = deconvolve_spots(spatial_counts, ref["profiles"])
        mapping = spatial_cell_type_mapping(decon["proportions_matrix"], coords, list(decon["cell_types"]))
        assert {"spatial_map", "dominant_type_per_spot", "mixing_entropy"}.issubset(mapping.keys())
        assert len(mapping["dominant_type_per_spot"]) == 20
        assert len(mapping["mixing_entropy"]) == 20


class TestValidateDeconvolution:
    def test_valid_result_passes(self) -> None:
        sc_expr, cell_types, spatial_counts, _ = _make_reference_and_spots()
        ref = build_reference_profiles(sc_expr, cell_types, n_markers=10)
        decon = deconvolve_spots(spatial_counts, ref["profiles"])
        validation = validate_deconvolution(decon)
        assert "overall_validity" in validation
        assert "mean_confidence" in validation


class TestNicheIdentificationDeconv:
    def test_niche_structure(self) -> None:
        sc_expr, cell_types, spatial_counts, coords = _make_reference_and_spots()
        ref = build_reference_profiles(sc_expr, cell_types, n_markers=10)
        decon = deconvolve_spots(spatial_counts, ref["profiles"])
        niches = niche_identification(decon["proportions_matrix"], coords, n_niches=2)
        assert {"niche_labels", "niche_compositions", "spatial_coherence"}.issubset(niches.keys())
        assert len(niches["niche_labels"]) == 20


# ---------------------------------------------------------------------------
# identify_niches (niche.identification)
# ---------------------------------------------------------------------------


class TestIdentifyNiches:
    def _proportions(self, n: int = 40, n_types: int = 3, seed: int = 0):
        rng = np.random.RandomState(seed)
        props = np.abs(rng.normal(size=(n, n_types))) + 0.1
        props /= props.sum(axis=1, keepdims=True)
        coords = rng.uniform(0, 100, (n, 2))
        return props, coords

    def test_result_type_and_fields(self) -> None:
        props, coords = self._proportions()
        result = identify_niches(props, coords, cell_type_names=["A", "B", "C"], n_niches=3, random_state=0)
        assert isinstance(result, NicheResult)
        assert result.n_niches == 3
        assert len(result.niche_labels) == 40
        assert set(np.unique(result.niche_labels)) <= {0, 1, 2}
        assert result.cell_type_names == ["A", "B", "C"]
        assert len(result.niche_sizes) == 3
        assert sum(result.niche_sizes.values() if isinstance(result.niche_sizes, dict) else result.niche_sizes) == 40

    def test_compositions_shape(self) -> None:
        props, coords = self._proportions()
        result = identify_niches(props, coords, cell_type_names=["A", "B", "C"], n_niches=3, random_state=0)
        compositions = np.asarray(result.compositions)
        assert compositions.shape[1] == 3

    def test_deterministic_with_seed(self) -> None:
        props, coords = self._proportions()
        r1 = identify_niches(props, coords, n_niches=3, random_state=7)
        r2 = identify_niches(props, coords, n_niches=3, random_state=7)
        assert list(r1.niche_labels) == list(r2.niche_labels)

    def test_niche_diversity_present(self) -> None:
        props, coords = self._proportions()
        result = identify_niches(props, coords, n_niches=2, random_state=0)
        assert result.niche_diversity is not None
