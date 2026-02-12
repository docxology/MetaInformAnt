"""Tests for spatial neighborhood analysis: enrichment, interaction, L-R, niches, Ripley K.

NO MOCKING POLICY: All tests use real implementations.
"""
from __future__ import annotations

import numpy as np
import pytest

from metainformant.spatial.analysis.neighborhood import (
    InteractionResult,
    NeighborhoodEnrichmentResult,
    NicheResult,
    RipleyKResult,
    compute_interaction_matrix,
    ligand_receptor_spatial,
    neighborhood_enrichment,
    niche_detection,
    ripley_k,
)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _make_spatial_cells(n=50, n_types=3, seed=42):
    rng = np.random.RandomState(seed)
    coords = rng.rand(n, 2) * 100
    types = rng.choice(["typeA", "typeB", "typeC"][:n_types], size=n)
    return types, coords


def _make_adjacency(coords, k=5, seed=42):
    from scipy.spatial import KDTree
    from scipy import sparse

    tree = KDTree(coords)
    _, indices = tree.query(coords, k=k + 1)
    n = len(coords)
    rows, cols = [], []
    for i in range(n):
        for j_idx in range(1, k + 1):
            j = int(indices[i, j_idx])
            rows.extend([i, j])
            cols.extend([j, i])
    vals = [1.0] * len(rows)
    return sparse.csr_matrix((vals, (rows, cols)), shape=(n, n))


# ---------------------------------------------------------------------------
# neighborhood_enrichment
# ---------------------------------------------------------------------------


class TestNeighborhoodEnrichment:
    def test_basic_enrichment(self):
        types, coords = _make_spatial_cells(n=30, n_types=2)
        result = neighborhood_enrichment(types, coords, n_neighbors=4, n_permutations=50, seed=42)
        assert isinstance(result, NeighborhoodEnrichmentResult)
        assert result.enrichment_matrix.shape[0] == 2
        assert len(result.cell_type_names) == 2

    def test_radius_based(self):
        types, coords = _make_spatial_cells(n=30)
        result = neighborhood_enrichment(types, coords, radius=30.0, n_permutations=50, seed=42)
        assert isinstance(result, NeighborhoodEnrichmentResult)

    def test_p_values_range(self):
        types, coords = _make_spatial_cells(n=40, n_types=2)
        result = neighborhood_enrichment(types, coords, n_neighbors=4, n_permutations=100, seed=42)
        assert np.all(result.p_values >= 0.0)
        assert np.all(result.p_values <= 1.0)

    def test_count_matrix_non_negative(self):
        types, coords = _make_spatial_cells(n=30)
        result = neighborhood_enrichment(types, coords, n_neighbors=4, n_permutations=50, seed=42)
        assert np.all(result.count_matrix >= 0)

    def test_expected_matrix_non_negative(self):
        types, coords = _make_spatial_cells(n=30)
        result = neighborhood_enrichment(types, coords, n_neighbors=4, n_permutations=50, seed=42)
        assert np.all(result.expected_matrix >= 0)


# ---------------------------------------------------------------------------
# compute_interaction_matrix
# ---------------------------------------------------------------------------


class TestComputeInteractionMatrix:
    def test_basic_interaction(self):
        types, coords = _make_spatial_cells(n=40)
        adj = _make_adjacency(coords, k=4)
        result = compute_interaction_matrix(types, adj)
        assert isinstance(result, InteractionResult)
        assert result.interaction_matrix.shape[0] == 3

    def test_normalized_method(self):
        types, coords = _make_spatial_cells(n=40)
        adj = _make_adjacency(coords, k=4)
        result = compute_interaction_matrix(types, adj, normalize_by_type_frequency=True)
        assert result.method == "frequency_normalized"

    def test_raw_method(self):
        types, coords = _make_spatial_cells(n=40)
        adj = _make_adjacency(coords, k=4)
        result = compute_interaction_matrix(types, adj, normalize_by_type_frequency=False)
        assert result.method == "raw_count"

    def test_symmetric_matrix(self):
        types, coords = _make_spatial_cells(n=40)
        adj = _make_adjacency(coords, k=4)
        result = compute_interaction_matrix(types, adj)
        np.testing.assert_allclose(result.interaction_matrix, result.interaction_matrix.T, atol=1e-10)


# ---------------------------------------------------------------------------
# ligand_receptor_spatial
# ---------------------------------------------------------------------------


class TestLigandReceptorSpatial:
    def test_basic_lr(self):
        rng = np.random.RandomState(42)
        n_spots = 30
        expression = rng.rand(n_spots, 5)
        coords = rng.rand(n_spots, 2) * 100
        gene_names = ["ligandA", "receptorA", "ligandB", "receptorB", "other"]
        lr_pairs = [("ligandA", "receptorA"), ("ligandB", "receptorB")]

        result = ligand_receptor_spatial(expression, lr_pairs, coords, gene_names=gene_names, n_neighbors=4)
        assert result["n_pairs_tested"] == 2
        assert ("ligandA", "receptorA") in result["scores"]

    def test_missing_gene_skipped(self):
        rng = np.random.RandomState(42)
        expression = rng.rand(20, 3)
        coords = rng.rand(20, 2) * 100
        gene_names = ["geneA", "geneB", "geneC"]
        lr_pairs = [("geneA", "geneB"), ("missing1", "missing2")]

        result = ligand_receptor_spatial(expression, lr_pairs, coords, gene_names=gene_names, n_neighbors=3)
        assert result["n_pairs_tested"] == 1

    def test_radius_based_neighbors(self):
        rng = np.random.RandomState(42)
        expression = rng.rand(20, 4)
        coords = rng.rand(20, 2) * 100
        gene_names = ["L1", "R1", "L2", "R2"]
        lr_pairs = [("L1", "R1")]

        result = ligand_receptor_spatial(expression, lr_pairs, coords, gene_names=gene_names, radius=50.0)
        assert result["n_pairs_tested"] == 1

    def test_per_spot_scores_shape(self):
        rng = np.random.RandomState(42)
        n_spots = 25
        expression = rng.rand(n_spots, 3)
        coords = rng.rand(n_spots, 2) * 100
        gene_names = ["L", "R", "other"]
        lr_pairs = [("L", "R")]

        result = ligand_receptor_spatial(expression, lr_pairs, coords, gene_names=gene_names, n_neighbors=4)
        assert result["per_spot_scores"][("L", "R")].shape == (n_spots,)


# ---------------------------------------------------------------------------
# niche_detection
# ---------------------------------------------------------------------------


class TestNicheDetection:
    def test_basic_niches(self):
        types, coords = _make_spatial_cells(n=50, n_types=3)
        result = niche_detection(types, coords, n_niches=3, n_neighbors=5, seed=42)
        assert isinstance(result, NicheResult)
        assert result.n_niches == 3
        assert len(result.niche_labels) == 50

    def test_compositions_sum_to_one(self):
        types, coords = _make_spatial_cells(n=50, n_types=3)
        result = niche_detection(types, coords, n_niches=3, n_neighbors=5, seed=42)
        row_sums = result.niche_compositions.sum(axis=1)
        np.testing.assert_allclose(row_sums, 1.0, atol=1e-10)

    def test_niche_labels_valid(self):
        types, coords = _make_spatial_cells(n=50)
        result = niche_detection(types, coords, n_niches=4, n_neighbors=5, seed=42)
        assert np.all(result.niche_labels >= 0)
        assert np.all(result.niche_labels < result.n_niches)


# ---------------------------------------------------------------------------
# ripley_k
# ---------------------------------------------------------------------------


class TestRipleyK:
    def test_basic_ripley(self):
        rng = np.random.RandomState(42)
        points = rng.rand(30, 2) * 100
        radii = np.array([10.0, 20.0, 30.0])
        result = ripley_k(points, radii, area=10000.0, n_simulations=10, seed=42)
        assert isinstance(result, RipleyKResult)
        assert len(result.k_values) == 3
        assert len(result.l_values) == 3
        assert result.n_points == 30

    def test_single_point(self):
        points = np.array([[50.0, 50.0]])
        radii = np.array([10.0, 20.0])
        result = ripley_k(points, radii, area=10000.0, n_simulations=5, seed=42)
        np.testing.assert_allclose(result.k_values, 0.0)

    def test_csr_envelope_returned(self):
        rng = np.random.RandomState(42)
        points = rng.rand(30, 2) * 100
        radii = np.array([10.0, 20.0, 30.0])
        result = ripley_k(points, radii, area=10000.0, n_simulations=20, seed=42)
        assert result.csr_envelope_lower is not None
        assert result.csr_envelope_upper is not None
        assert len(result.csr_envelope_lower) == 3
