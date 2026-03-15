"""Tests for spatial deconvolution: NNLS, cell fractions, enrichment, deconvolve_spots.

NO MOCKING POLICY: All tests use real implementations.
"""
from __future__ import annotations

import numpy as np
import pytest

from metainformant.spatial.analysis.deconvolution import (
    DeconvolutionResult,
    create_reference_profiles,
    deconvolve_spots,
    enrichment_score,
    estimate_cell_fractions,
    nnls_deconvolution,
)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _make_reference_data(n_cells=100, n_genes=50, n_types=3, seed=42):
    rng = np.random.RandomState(seed)
    data = rng.rand(n_cells, n_genes)
    type_names = ["type_A", "type_B", "type_C"][:n_types]
    labels = np.array([type_names[i % n_types] for i in range(n_cells)])
    return data, labels


def _make_spatial_data(n_spots=20, n_genes=50, seed=42):
    rng = np.random.RandomState(seed)
    return rng.rand(n_spots, n_genes)


# ---------------------------------------------------------------------------
# create_reference_profiles
# ---------------------------------------------------------------------------


class TestCreateReferenceProfiles:
    def test_basic_profiles(self):
        data, labels = _make_reference_data()
        profiles, type_names, gene_names = create_reference_profiles(data, labels)
        assert profiles.shape == (3, 50)
        assert len(type_names) == 3
        assert len(gene_names) == 50

    def test_mean_method(self):
        data, labels = _make_reference_data()
        profiles, _, _ = create_reference_profiles(data, labels, method="mean")
        assert profiles.shape[0] == 3

    def test_median_method(self):
        data, labels = _make_reference_data()
        profiles, _, _ = create_reference_profiles(data, labels, method="median")
        assert profiles.shape[0] == 3

    def test_custom_gene_names(self):
        data, labels = _make_reference_data(n_genes=10)
        gene_names = [f"gene_{i}" for i in range(10)]
        _, _, returned_names = create_reference_profiles(data, labels, gene_names=gene_names)
        assert returned_names == gene_names


# ---------------------------------------------------------------------------
# nnls_deconvolution
# ---------------------------------------------------------------------------


class TestNNLSDeconvolution:
    def test_single_spot(self):
        ref = np.array([[1.0, 0.0, 0.5], [0.0, 1.0, 0.5]])  # 2 types, 3 genes
        bulk = np.array([0.5, 0.5, 0.5])
        weights, residual = nnls_deconvolution(bulk, ref)
        assert len(weights) == 2
        assert residual >= 0

    def test_batch_spots(self):
        ref = np.array([[1.0, 0.0], [0.0, 1.0]])  # 2 types, 2 genes
        bulk = np.array([[0.8, 0.2], [0.3, 0.7], [0.5, 0.5]])  # 3 spots
        weights, residuals = nnls_deconvolution(bulk, ref)
        assert weights.shape == (3, 2)
        assert len(residuals) == 3

    def test_non_negative_weights(self):
        rng = np.random.RandomState(42)
        ref = rng.rand(3, 20)
        bulk = rng.rand(10, 20)
        weights, _ = nnls_deconvolution(bulk, ref)
        assert np.all(weights >= 0)


# ---------------------------------------------------------------------------
# estimate_cell_fractions
# ---------------------------------------------------------------------------


class TestEstimateCellFractions:
    def test_single_vector(self):
        weights = np.array([2.0, 3.0, 5.0])
        fractions = estimate_cell_fractions(weights)
        assert fractions.sum() == pytest.approx(1.0)
        assert fractions[2] == pytest.approx(0.5)

    def test_matrix_rows_sum_to_one(self):
        weights = np.array([[1.0, 3.0], [2.0, 2.0], [0.5, 1.5]])
        fractions = estimate_cell_fractions(weights)
        row_sums = fractions.sum(axis=1)
        np.testing.assert_allclose(row_sums, 1.0, atol=1e-10)

    def test_zero_weights_handled(self):
        weights = np.array([0.0, 0.0, 0.0])
        fractions = estimate_cell_fractions(weights)
        assert np.all(np.isfinite(fractions))

    def test_deconvolution_result_input(self):
        dr = DeconvolutionResult(
            weights=np.array([[1.0, 2.0], [3.0, 1.0]]),
            fractions=np.array([[0.33, 0.67], [0.75, 0.25]]),
            cell_type_names=["A", "B"],
            residuals=np.array([0.1, 0.2]),
            method="nnls",
        )
        fractions = estimate_cell_fractions(dr)
        row_sums = fractions.sum(axis=1)
        np.testing.assert_allclose(row_sums, 1.0, atol=1e-10)


# ---------------------------------------------------------------------------
# enrichment_score
# ---------------------------------------------------------------------------


class TestEnrichmentScore:
    def test_equal_fractions_zero_enrichment(self):
        observed = np.array([0.25, 0.25, 0.25, 0.25])
        expected = np.array([0.25, 0.25, 0.25, 0.25])
        scores = enrichment_score(observed, expected)
        np.testing.assert_allclose(scores, 0.0, atol=0.01)

    def test_enriched_type(self):
        observed = np.array([0.8, 0.1, 0.1])
        expected = np.array([0.33, 0.33, 0.34])
        scores = enrichment_score(observed, expected)
        assert scores[0] > 0  # Type 0 is enriched

    def test_depleted_type(self):
        observed = np.array([0.05, 0.05, 0.9])
        expected = np.array([0.33, 0.33, 0.34])
        scores = enrichment_score(observed, expected)
        assert scores[0] < 0  # Type 0 is depleted

    def test_batch_enrichment(self):
        observed = np.array([[0.5, 0.5], [0.9, 0.1]])
        expected = np.array([0.5, 0.5])
        scores = enrichment_score(observed, expected)
        assert scores.shape == (2, 2)


# ---------------------------------------------------------------------------
# deconvolve_spots
# ---------------------------------------------------------------------------


class TestDeconvolveSpots:
    def test_nnls_method(self):
        data, labels = _make_reference_data(n_cells=60, n_genes=20)
        ref, type_names, gene_names = create_reference_profiles(data, labels)
        spatial = _make_spatial_data(n_spots=10, n_genes=20)

        result = deconvolve_spots(spatial, ref, method="nnls", cell_type_names=type_names)
        assert isinstance(result, DeconvolutionResult)
        assert result.n_spots == 10
        assert result.n_types == 3
        assert result.method == "nnls"

    def test_nmf_method(self):
        data, labels = _make_reference_data(n_cells=60, n_genes=20)
        ref, type_names, gene_names = create_reference_profiles(data, labels)
        spatial = _make_spatial_data(n_spots=10, n_genes=20)

        result = deconvolve_spots(spatial, ref, method="nmf", cell_type_names=type_names)
        assert result.method == "nmf"

    def test_fractions_sum_to_one(self):
        data, labels = _make_reference_data(n_cells=60, n_genes=20)
        ref, type_names, gene_names = create_reference_profiles(data, labels)
        spatial = _make_spatial_data(n_spots=10, n_genes=20)

        result = deconvolve_spots(spatial, ref, method="nnls", cell_type_names=type_names)
        row_sums = result.fractions.sum(axis=1)
        np.testing.assert_allclose(row_sums, 1.0, atol=1e-10)

    def test_auto_type_names(self):
        ref = np.random.RandomState(42).rand(3, 10)
        spatial = np.random.RandomState(42).rand(5, 10)
        result = deconvolve_spots(spatial, ref, method="nnls")
        assert len(result.cell_type_names) == 3

    def test_unknown_method_raises(self):
        ref = np.random.RandomState(42).rand(3, 10)
        spatial = np.random.RandomState(42).rand(5, 10)
        with pytest.raises(ValueError, match="Unknown deconvolution"):
            deconvolve_spots(spatial, ref, method="invalid")

    def test_gene_name_intersection(self):
        ref = np.random.RandomState(42).rand(2, 5)
        spatial = np.random.RandomState(42).rand(3, 4)
        shared = ["geneA", "geneB", "geneC"]
        spatial_genes = ["geneA", "geneB", "geneC", "geneD"]
        ref_genes = ["geneA", "geneB", "geneC", "geneE", "geneF"]

        result = deconvolve_spots(
            spatial,
            ref,
            method="nnls",
            spatial_gene_names=spatial_genes,
            reference_gene_names=ref_genes,
        )
        assert result.metadata["n_genes_used"] == 3
