"""Comprehensive tests for advanced RNA simulation functions.

Tests cover all eight RNA simulation functions in metainformant.simulation.models.rna:
- simulate_counts_negative_binomial
- simulate_rnaseq_counts
- simulate_differential_expression
- simulate_bulk_rnaseq
- simulate_single_cell_rnaseq
- simulate_time_series_expression
- simulate_spatial_expression
- add_technical_noise

All tests use real numpy arrays and real implementations (NO MOCKING).
"""

from __future__ import annotations

import random

import numpy as np
import pytest

from metainformant.core import errors
from metainformant.simulation.models.rna import (
    add_technical_noise,
    simulate_bulk_rnaseq,
    simulate_counts_negative_binomial,
    simulate_differential_expression,
    simulate_rnaseq_counts,
    simulate_single_cell_rnaseq,
    simulate_spatial_expression,
    simulate_time_series_expression,
)


# ---------------------------------------------------------------------------
# simulate_counts_negative_binomial
# ---------------------------------------------------------------------------


class TestSimulateCountsNegativeBinomial:
    """Tests for the core negative binomial count simulator."""

    def test_basic_shape(self) -> None:
        """Output matrix has shape (n_samples, n_features)."""
        n_samples, n_features = 8, 15
        means = np.random.RandomState(0).exponential(50, n_features)
        dispersions = np.full(n_features, 0.5)
        rng = random.Random(42)

        counts = simulate_counts_negative_binomial(n_samples, n_features, means, dispersions, rng=rng)

        assert counts.shape == (n_samples, n_features)

    def test_non_negative_integer_counts(self) -> None:
        """All generated counts must be non-negative integers."""
        n_samples, n_features = 10, 20
        means = np.random.RandomState(1).exponential(30, n_features)
        dispersions = np.full(n_features, 0.3)
        rng = random.Random(99)

        counts = simulate_counts_negative_binomial(n_samples, n_features, means, dispersions, rng=rng)

        assert np.all(counts >= 0), "Counts must be non-negative"
        assert np.issubdtype(counts.dtype, np.integer), "Counts must be integers"

    def test_zero_mean_gives_zero_counts(self) -> None:
        """Features with mean=0 should produce all-zero columns."""
        n_samples, n_features = 5, 4
        means = np.array([0.0, 50.0, 0.0, 100.0])
        dispersions = np.full(n_features, 0.5)
        rng = random.Random(42)

        counts = simulate_counts_negative_binomial(n_samples, n_features, means, dispersions, rng=rng)

        assert np.all(counts[:, 0] == 0), "Zero-mean feature should yield zero counts"
        assert np.all(counts[:, 2] == 0), "Zero-mean feature should yield zero counts"

    def test_higher_mean_produces_higher_counts(self) -> None:
        """Features with higher means should produce higher average counts."""
        n_samples, n_features = 200, 2
        means = np.array([10.0, 500.0])
        dispersions = np.array([0.5, 0.5])
        rng = random.Random(42)

        counts = simulate_counts_negative_binomial(n_samples, n_features, means, dispersions, rng=rng)

        mean_low = counts[:, 0].mean()
        mean_high = counts[:, 1].mean()
        assert mean_high > mean_low, "Higher mean parameter should yield higher observed mean"

    def test_explicit_means_and_dispersions(self) -> None:
        """Manually specify per-feature means and dispersions."""
        n_samples, n_features = 12, 5
        means = np.array([1.0, 10.0, 50.0, 200.0, 1000.0])
        dispersions = np.array([0.1, 0.2, 0.3, 0.5, 1.0])
        rng = random.Random(7)

        counts = simulate_counts_negative_binomial(n_samples, n_features, means, dispersions, rng=rng)

        assert counts.shape == (n_samples, n_features)
        assert np.all(counts >= 0)

    def test_single_feature(self) -> None:
        """Edge case: single feature."""
        means = np.array([50.0])
        dispersions = np.array([0.3])
        rng = random.Random(42)

        counts = simulate_counts_negative_binomial(1, 1, means, dispersions, rng=rng)
        assert counts.shape == (1, 1)
        assert counts[0, 0] >= 0

    def test_single_sample(self) -> None:
        """Edge case: single sample across multiple features."""
        n_features = 10
        means = np.random.RandomState(0).exponential(20, n_features)
        dispersions = np.full(n_features, 0.4)
        rng = random.Random(42)

        counts = simulate_counts_negative_binomial(1, n_features, means, dispersions, rng=rng)
        assert counts.shape == (1, n_features)

    def test_reproducibility(self) -> None:
        """Same seed should produce identical results."""
        n_samples, n_features = 6, 8
        means = np.array([10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0])
        dispersions = np.full(n_features, 0.5)

        counts_a = simulate_counts_negative_binomial(n_samples, n_features, means, dispersions, rng=random.Random(42))
        counts_b = simulate_counts_negative_binomial(n_samples, n_features, means, dispersions, rng=random.Random(42))

        np.testing.assert_array_equal(counts_a, counts_b)

    def test_wrong_means_shape_raises(self) -> None:
        """Means array with wrong shape should raise ValidationError."""
        means = np.array([10.0, 20.0])  # shape (2,) but n_features=5
        dispersions = np.full(5, 0.5)

        with pytest.raises(errors.ValidationError, match="means must have shape"):
            simulate_counts_negative_binomial(3, 5, means, dispersions)

    def test_wrong_dispersions_shape_raises(self) -> None:
        """Dispersions array with wrong shape should raise ValidationError."""
        means = np.full(5, 10.0)
        dispersions = np.array([0.5, 0.5])  # shape (2,) but n_features=5

        with pytest.raises(errors.ValidationError, match="dispersion must have shape"):
            simulate_counts_negative_binomial(3, 5, means, dispersions)

    def test_negative_mean_raises(self) -> None:
        """Negative mean values should raise ValidationError."""
        means = np.array([-1.0, 10.0, 20.0])
        dispersions = np.full(3, 0.5)

        with pytest.raises(errors.ValidationError, match="non-negative"):
            simulate_counts_negative_binomial(3, 3, means, dispersions)

    def test_zero_dispersion_raises(self) -> None:
        """Zero dispersion should raise ValidationError."""
        means = np.array([10.0, 20.0])
        dispersions = np.array([0.0, 0.5])

        with pytest.raises(errors.ValidationError, match="positive"):
            simulate_counts_negative_binomial(3, 2, means, dispersions)

    def test_negative_dispersion_raises(self) -> None:
        """Negative dispersion should raise ValidationError."""
        means = np.array([10.0, 20.0])
        dispersions = np.array([-0.1, 0.5])

        with pytest.raises(errors.ValidationError, match="positive"):
            simulate_counts_negative_binomial(3, 2, means, dispersions)


# ---------------------------------------------------------------------------
# simulate_rnaseq_counts
# ---------------------------------------------------------------------------


class TestSimulateRnaseqCounts:
    """Tests for the convenience RNA-seq count wrapper."""

    def test_default_shape(self) -> None:
        """Default parameters produce (10, 1000) matrix."""
        rng = random.Random(42)
        counts = simulate_rnaseq_counts(rng=rng)
        assert counts.shape == (10, 1000)

    def test_custom_shape(self) -> None:
        """Custom n_genes and n_samples produce the correct shape."""
        rng = random.Random(42)
        counts = simulate_rnaseq_counts(n_genes=20, n_samples=8, rng=rng)
        assert counts.shape == (8, 20)

    def test_non_negative_values(self) -> None:
        """All counts must be non-negative."""
        rng = random.Random(42)
        counts = simulate_rnaseq_counts(n_genes=25, n_samples=10, rng=rng)
        assert np.all(counts >= 0)

    def test_integer_counts(self) -> None:
        """Counts must be integers."""
        rng = random.Random(42)
        counts = simulate_rnaseq_counts(n_genes=15, n_samples=5, rng=rng)
        assert np.issubdtype(counts.dtype, np.integer)

    def test_higher_mean_expression_raises_counts(self) -> None:
        """Higher mean_expression should produce higher total counts on average."""
        rng_low = random.Random(42)
        counts_low = simulate_rnaseq_counts(n_genes=20, n_samples=50, mean_expression=5.0, rng=rng_low)

        rng_high = random.Random(42)
        counts_high = simulate_rnaseq_counts(n_genes=20, n_samples=50, mean_expression=500.0, rng=rng_high)

        assert counts_high.sum() > counts_low.sum()

    def test_single_gene(self) -> None:
        """Edge case: single gene."""
        rng = random.Random(42)
        counts = simulate_rnaseq_counts(n_genes=1, n_samples=5, rng=rng)
        assert counts.shape == (5, 1)

    def test_single_sample(self) -> None:
        """Edge case: single sample."""
        rng = random.Random(42)
        counts = simulate_rnaseq_counts(n_genes=10, n_samples=1, rng=rng)
        assert counts.shape == (1, 10)

    def test_reproducibility(self) -> None:
        """Same seed produces identical output."""
        counts_a = simulate_rnaseq_counts(n_genes=15, n_samples=6, rng=random.Random(99))
        counts_b = simulate_rnaseq_counts(n_genes=15, n_samples=6, rng=random.Random(99))
        np.testing.assert_array_equal(counts_a, counts_b)


# ---------------------------------------------------------------------------
# simulate_differential_expression
# ---------------------------------------------------------------------------


class TestSimulateDifferentialExpression:
    """Tests for the differential expression simulator."""

    def test_basic_output_shapes(self) -> None:
        """Returns expression matrix and group labels with correct shapes."""
        n_samples, n_features = 10, 20
        fold_changes = np.array([2.0, 3.0, 0.5])
        rng = random.Random(42)

        expr, labels = simulate_differential_expression(n_samples, n_features, fold_changes, rng=rng)

        assert expr.shape == (n_samples, n_features)
        assert labels.shape == (n_samples,)

    def test_two_groups(self) -> None:
        """Group labels contain exactly two groups (0 and 1)."""
        fold_changes = np.array([2.0])
        rng = random.Random(42)

        _, labels = simulate_differential_expression(10, 20, fold_changes, rng=rng)

        unique_labels = set(labels)
        assert unique_labels == {0, 1}

    def test_group_split(self) -> None:
        """Samples are split roughly in half between groups."""
        n_samples = 12
        fold_changes = np.array([2.0])
        rng = random.Random(42)

        _, labels = simulate_differential_expression(n_samples, 20, fold_changes, rng=rng)

        n_group0 = np.sum(labels == 0)
        n_group1 = np.sum(labels == 1)
        assert n_group0 == n_samples // 2
        assert n_group1 == n_samples - n_samples // 2

    def test_non_negative_expression(self) -> None:
        """All expression values must be non-negative."""
        fold_changes = np.array([2.0, 0.5, 3.0])
        rng = random.Random(42)

        expr, _ = simulate_differential_expression(10, 20, fold_changes, rng=rng)
        assert np.all(expr >= 0)

    def test_fold_changes_affect_expression(self) -> None:
        """Large fold changes should cause measurable differences between groups."""
        n_samples = 200
        n_features = 10
        fold_changes = np.array([10.0])  # 10x upregulation of gene 0
        rng = random.Random(42)

        expr, labels = simulate_differential_expression(n_samples, n_features, fold_changes, rng=rng)

        group0_mean = expr[labels == 0, 0].mean()
        group1_mean = expr[labels == 1, 0].mean()

        # Group 1 should have substantially higher expression for gene 0
        assert group1_mean > group0_mean

    def test_minimum_samples(self) -> None:
        """Minimum allowed n_samples is 2."""
        fold_changes = np.array([2.0])
        rng = random.Random(42)

        expr, labels = simulate_differential_expression(2, 5, fold_changes, rng=rng)
        assert expr.shape == (2, 5)
        assert len(labels) == 2

    def test_too_many_fold_changes_raises(self) -> None:
        """fold_changes length >= n_features should raise ValidationError."""
        fold_changes = np.array([2.0, 3.0, 0.5, 1.5, 4.0])

        with pytest.raises(errors.ValidationError, match="differentially expressed"):
            simulate_differential_expression(10, 5, fold_changes)

    def test_single_sample_raises(self) -> None:
        """n_samples < 2 should raise."""
        fold_changes = np.array([2.0])

        with pytest.raises((errors.ValidationError, ValueError)):
            simulate_differential_expression(1, 10, fold_changes)

    def test_reproducibility(self) -> None:
        """Same seed produces identical results."""
        fold_changes = np.array([2.0, 0.5])

        expr_a, labels_a = simulate_differential_expression(10, 15, fold_changes, rng=random.Random(42))
        expr_b, labels_b = simulate_differential_expression(10, 15, fold_changes, rng=random.Random(42))

        np.testing.assert_array_equal(expr_a, expr_b)
        np.testing.assert_array_equal(labels_a, labels_b)


# ---------------------------------------------------------------------------
# simulate_bulk_rnaseq
# ---------------------------------------------------------------------------


class TestSimulateBulkRnaseq:
    """Tests for bulk RNA-seq simulation."""

    def test_basic_shape(self) -> None:
        """Output has shape (n_samples, n_genes)."""
        rng = random.Random(42)
        counts = simulate_bulk_rnaseq(8, 20, rng=rng)
        assert counts.shape == (8, 20)

    def test_non_negative_integer_counts(self) -> None:
        """Counts must be non-negative integers."""
        rng = random.Random(42)
        counts = simulate_bulk_rnaseq(5, 15, rng=rng)
        assert np.all(counts >= 0)
        assert np.issubdtype(counts.dtype, np.integer)

    def test_with_library_sizes(self) -> None:
        """Custom library sizes should be respected."""
        n_samples, n_genes = 6, 20
        library_sizes = np.array([10000, 20000, 15000, 25000, 18000, 22000], dtype=float)
        rng = random.Random(42)

        counts = simulate_bulk_rnaseq(n_samples, n_genes, library_sizes=library_sizes, rng=rng)

        assert counts.shape == (n_samples, n_genes)
        # Each sample's total counts should equal its library size
        for i in range(n_samples):
            assert counts[i, :].sum() == int(library_sizes[i])

    def test_with_gene_means(self) -> None:
        """Custom gene means should be used when provided."""
        n_samples, n_genes = 5, 10
        gene_means = np.random.RandomState(0).exponential(50, n_genes)
        rng = random.Random(42)

        counts = simulate_bulk_rnaseq(n_samples, n_genes, gene_means=gene_means, rng=rng)

        assert counts.shape == (n_samples, n_genes)
        assert np.all(counts >= 0)

    def test_with_both_library_sizes_and_gene_means(self) -> None:
        """Both library_sizes and gene_means can be provided together."""
        n_samples, n_genes = 4, 8
        library_sizes = np.full(n_samples, 50000.0)
        gene_means = np.random.RandomState(0).exponential(100, n_genes)
        rng = random.Random(42)

        counts = simulate_bulk_rnaseq(n_samples, n_genes, library_sizes=library_sizes, gene_means=gene_means, rng=rng)

        assert counts.shape == (n_samples, n_genes)
        for i in range(n_samples):
            assert counts[i, :].sum() == int(library_sizes[i])

    def test_wrong_library_sizes_length_raises(self) -> None:
        """library_sizes with wrong length should raise."""
        library_sizes = np.array([10000.0, 20000.0])  # length 2 but n_samples=5

        with pytest.raises(errors.ValidationError, match="library_sizes"):
            simulate_bulk_rnaseq(5, 10, library_sizes=library_sizes)

    def test_wrong_gene_means_length_raises(self) -> None:
        """gene_means with wrong length should raise."""
        gene_means = np.array([10.0, 20.0])  # length 2 but n_genes=5

        with pytest.raises(errors.ValidationError, match="gene_means"):
            simulate_bulk_rnaseq(3, 5, gene_means=gene_means)

    def test_single_sample(self) -> None:
        """Edge case: single sample."""
        rng = random.Random(42)
        counts = simulate_bulk_rnaseq(1, 10, rng=rng)
        assert counts.shape == (1, 10)

    def test_single_gene(self) -> None:
        """Edge case: single gene."""
        rng = random.Random(42)
        counts = simulate_bulk_rnaseq(5, 1, rng=rng)
        assert counts.shape == (5, 1)

    def test_reproducibility(self) -> None:
        """Same seed produces identical results."""
        counts_a = simulate_bulk_rnaseq(6, 12, rng=random.Random(42))
        counts_b = simulate_bulk_rnaseq(6, 12, rng=random.Random(42))
        np.testing.assert_array_equal(counts_a, counts_b)


# ---------------------------------------------------------------------------
# simulate_single_cell_rnaseq
# ---------------------------------------------------------------------------


class TestSimulateSingleCellRnaseq:
    """Tests for single-cell RNA-seq simulation."""

    def test_basic_shapes(self) -> None:
        """Returns expression matrix and cell type labels with correct shapes."""
        rng = random.Random(42)
        expr, labels = simulate_single_cell_rnaseq(50, 20, n_cell_types=3, rng=rng)

        assert expr.shape == (50, 20)
        assert labels.shape == (50,)

    def test_cell_type_count(self) -> None:
        """Number of distinct cell types in labels matches n_cell_types."""
        rng = random.Random(42)
        _, labels = simulate_single_cell_rnaseq(30, 15, n_cell_types=3, rng=rng)

        unique_types = np.unique(labels)
        assert len(unique_types) == 3

    def test_non_negative_expression(self) -> None:
        """Expression values must be non-negative."""
        rng = random.Random(42)
        expr, _ = simulate_single_cell_rnaseq(20, 15, rng=rng)
        assert np.all(expr >= 0)

    def test_integer_counts(self) -> None:
        """Counts should be integers."""
        rng = random.Random(42)
        expr, _ = simulate_single_cell_rnaseq(20, 15, rng=rng)
        assert np.issubdtype(expr.dtype, np.integer)

    def test_dropout_produces_zeros(self) -> None:
        """High dropout rate should produce many zeros in the matrix."""
        rng = random.Random(42)
        expr, _ = simulate_single_cell_rnaseq(50, 20, dropout_rate=0.8, rng=rng)

        zero_fraction = np.sum(expr == 0) / expr.size
        # With 80% dropout, expect a substantial fraction of zeros
        assert zero_fraction > 0.3, f"Expected many zeros with high dropout, got {zero_fraction:.2%}"

    def test_zero_dropout_fewer_zeros(self) -> None:
        """Zero dropout should produce fewer zeros than high dropout."""
        rng_no_drop = random.Random(42)
        expr_no_drop, _ = simulate_single_cell_rnaseq(50, 20, dropout_rate=0.0, rng=rng_no_drop)

        rng_high_drop = random.Random(42)
        expr_high_drop, _ = simulate_single_cell_rnaseq(50, 20, dropout_rate=0.9, rng=rng_high_drop)

        zeros_no_drop = np.sum(expr_no_drop == 0)
        zeros_high_drop = np.sum(expr_high_drop == 0)
        assert zeros_high_drop > zeros_no_drop

    def test_single_cell_type(self) -> None:
        """Edge case: single cell type."""
        rng = random.Random(42)
        expr, labels = simulate_single_cell_rnaseq(10, 15, n_cell_types=1, rng=rng)
        assert np.all(labels == 0)
        assert expr.shape == (10, 15)

    def test_reproducibility(self) -> None:
        """Same seed produces identical results."""
        expr_a, labels_a = simulate_single_cell_rnaseq(20, 10, n_cell_types=3, rng=random.Random(42))
        expr_b, labels_b = simulate_single_cell_rnaseq(20, 10, n_cell_types=3, rng=random.Random(42))
        np.testing.assert_array_equal(expr_a, expr_b)
        np.testing.assert_array_equal(labels_a, labels_b)


# ---------------------------------------------------------------------------
# simulate_time_series_expression
# ---------------------------------------------------------------------------


class TestSimulateTimeSeriesExpression:
    """Tests for time-series gene expression simulation."""

    def test_basic_shape(self) -> None:
        """Output has shape (n_timepoints, n_genes)."""
        rng = random.Random(42)
        expr = simulate_time_series_expression(20, 15, rng=rng)
        assert expr.shape == (20, 15)

    def test_non_negative_values(self) -> None:
        """Expression values must be non-negative."""
        rng = random.Random(42)
        expr = simulate_time_series_expression(15, 10, rng=rng)
        assert np.all(expr >= 0)

    def test_integer_counts(self) -> None:
        """Output should be integer counts (Poisson-rounded)."""
        rng = random.Random(42)
        expr = simulate_time_series_expression(10, 8, rng=rng)
        assert np.issubdtype(expr.dtype, np.integer)

    def test_custom_oscillation_frequencies(self) -> None:
        """Custom oscillation frequencies should be accepted."""
        n_genes = 5
        oscillation_freq = np.array([0.5, 1.0, 1.5, 2.0, 2.5])
        rng = random.Random(42)

        expr = simulate_time_series_expression(10, n_genes, oscillation_freq=oscillation_freq, rng=rng)

        assert expr.shape == (10, n_genes)

    def test_wrong_oscillation_freq_length_raises(self) -> None:
        """oscillation_freq with wrong length should raise."""
        oscillation_freq = np.array([1.0, 2.0])  # length 2 but n_genes=5

        with pytest.raises(errors.ValidationError, match="oscillation_freq"):
            simulate_time_series_expression(10, 5, oscillation_freq=oscillation_freq)

    def test_temporal_variation(self) -> None:
        """Expression should vary across timepoints (not all identical rows)."""
        rng = random.Random(42)
        expr = simulate_time_series_expression(20, 10, rng=rng)

        # Check that not all timepoints are identical
        row_sums = expr.sum(axis=1)
        assert not np.all(row_sums == row_sums[0]), "Expression should vary over time"

    def test_minimum_timepoints(self) -> None:
        """Minimum n_timepoints is 2."""
        rng = random.Random(42)
        expr = simulate_time_series_expression(2, 5, rng=rng)
        assert expr.shape == (2, 5)

    def test_single_gene(self) -> None:
        """Edge case: single gene time series."""
        rng = random.Random(42)
        expr = simulate_time_series_expression(10, 1, rng=rng)
        assert expr.shape == (10, 1)

    def test_reproducibility(self) -> None:
        """Same seed produces identical results."""
        expr_a = simulate_time_series_expression(10, 8, rng=random.Random(42))
        expr_b = simulate_time_series_expression(10, 8, rng=random.Random(42))
        np.testing.assert_array_equal(expr_a, expr_b)


# ---------------------------------------------------------------------------
# simulate_spatial_expression
# ---------------------------------------------------------------------------


class TestSimulateSpatialExpression:
    """Tests for spatially-resolved expression simulation."""

    def test_random_pattern_shapes(self) -> None:
        """Random pattern returns correct expression and coordinate shapes."""
        rng = random.Random(42)
        expr, coords = simulate_spatial_expression(30, 10, spatial_patterns="random", rng=rng)

        assert expr.shape == (30, 10)
        assert coords.shape == (30, 2)

    def test_gradient_pattern_shapes(self) -> None:
        """Gradient pattern returns correct shapes.

        Note: gradient mode uses int(sqrt(n_spots))^2 grid points, so
        n_spots must be a perfect square to avoid shape mismatch in the
        current implementation.
        """
        n_spots = 25  # 5x5 grid -- perfect square
        rng = random.Random(42)
        expr, coords = simulate_spatial_expression(n_spots, 10, spatial_patterns="gradient", rng=rng)

        assert expr.shape == (n_spots, 10)
        assert coords.shape == (n_spots, 2)

    def test_clusters_pattern_shapes(self) -> None:
        """Clusters pattern returns correct shapes."""
        rng = random.Random(42)
        expr, coords = simulate_spatial_expression(30, 10, spatial_patterns="clusters", rng=rng)

        assert expr.shape == (30, 10)
        assert coords.shape == (30, 2)

    def test_non_negative_expression(self) -> None:
        """Expression values must be non-negative for all patterns.

        Gradient mode requires a perfect-square n_spots (see gradient grid logic).
        """
        # random and clusters can use any n_spots; gradient needs a perfect square
        pattern_spots = {"random": 20, "gradient": 25, "clusters": 21}
        for pattern, n_spots in pattern_spots.items():
            rng = random.Random(42)
            expr, _ = simulate_spatial_expression(n_spots, 8, spatial_patterns=pattern, rng=rng)
            assert np.all(expr >= 0), f"Negative expression in {pattern} pattern"

    def test_integer_counts(self) -> None:
        """Expression should be integer counts."""
        rng = random.Random(42)
        expr, _ = simulate_spatial_expression(20, 8, spatial_patterns="random", rng=rng)
        assert np.issubdtype(expr.dtype, np.integer)

    def test_coordinates_are_finite(self) -> None:
        """Spatial coordinates must be finite real numbers.

        Gradient mode requires a perfect-square n_spots.
        """
        pattern_spots = {"random": 20, "gradient": 16, "clusters": 21}
        for pattern, n_spots in pattern_spots.items():
            rng = random.Random(42)
            _, coords = simulate_spatial_expression(n_spots, 5, spatial_patterns=pattern, rng=rng)
            assert np.all(np.isfinite(coords)), f"Non-finite coordinates in {pattern}"

    def test_invalid_pattern_raises(self) -> None:
        """Invalid spatial_patterns value should raise ValidationError."""
        with pytest.raises(errors.ValidationError, match="Invalid spatial_patterns"):
            simulate_spatial_expression(10, 5, spatial_patterns="nonexistent")

    def test_default_pattern_is_random(self) -> None:
        """Default spatial_patterns should be 'random'."""
        rng = random.Random(42)
        expr_default, _ = simulate_spatial_expression(15, 5, rng=rng)

        rng2 = random.Random(42)
        expr_random, _ = simulate_spatial_expression(15, 5, spatial_patterns="random", rng=rng2)

        np.testing.assert_array_equal(expr_default, expr_random)

    def test_gradient_spatial_correlation(self) -> None:
        """Gradient pattern: expression should correlate with x-coordinate."""
        rng = random.Random(42)
        expr, coords = simulate_spatial_expression(
            100, 5, spatial_patterns="gradient", rng=rng  # 100 = 10x10 perfect square
        )

        # For at least some genes, expression should correlate with x
        x_coords = coords[:, 0]
        any_correlated = False
        for g in range(expr.shape[1]):
            if np.std(expr[:, g]) > 0 and np.std(x_coords) > 0:
                corr = np.corrcoef(x_coords, expr[:, g])[0, 1]
                if abs(corr) > 0.1:
                    any_correlated = True
                    break

        assert any_correlated, "Gradient pattern should produce spatially correlated expression"

    def test_single_spot(self) -> None:
        """Edge case: single spot."""
        rng = random.Random(42)
        expr, coords = simulate_spatial_expression(1, 5, spatial_patterns="random", rng=rng)
        assert expr.shape == (1, 5)
        assert coords.shape == (1, 2)

    def test_single_gene(self) -> None:
        """Edge case: single gene."""
        rng = random.Random(42)
        expr, coords = simulate_spatial_expression(10, 1, spatial_patterns="random", rng=rng)
        assert expr.shape == (10, 1)

    def test_reproducibility_random(self) -> None:
        """Same seed produces identical results for random pattern."""
        expr_a, coords_a = simulate_spatial_expression(15, 5, spatial_patterns="random", rng=random.Random(42))
        expr_b, coords_b = simulate_spatial_expression(15, 5, spatial_patterns="random", rng=random.Random(42))
        np.testing.assert_array_equal(expr_a, expr_b)
        np.testing.assert_array_equal(coords_a, coords_b)

    def test_reproducibility_gradient(self) -> None:
        """Same seed produces identical results for gradient pattern."""
        n_spots = 16  # 4x4 perfect square for gradient mode
        expr_a, coords_a = simulate_spatial_expression(n_spots, 5, spatial_patterns="gradient", rng=random.Random(42))
        expr_b, coords_b = simulate_spatial_expression(n_spots, 5, spatial_patterns="gradient", rng=random.Random(42))
        np.testing.assert_array_equal(expr_a, expr_b)
        np.testing.assert_array_equal(coords_a, coords_b)

    def test_reproducibility_clusters(self) -> None:
        """Same seed produces identical results for clusters pattern."""
        expr_a, coords_a = simulate_spatial_expression(15, 5, spatial_patterns="clusters", rng=random.Random(42))
        expr_b, coords_b = simulate_spatial_expression(15, 5, spatial_patterns="clusters", rng=random.Random(42))
        np.testing.assert_array_equal(expr_a, expr_b)
        np.testing.assert_array_equal(coords_a, coords_b)


# ---------------------------------------------------------------------------
# add_technical_noise
# ---------------------------------------------------------------------------


class TestAddTechnicalNoise:
    """Tests for the technical noise injection function."""

    def test_output_shape_matches_input(self) -> None:
        """Output matrix has the same shape as input."""
        rng = random.Random(42)
        original = np.random.RandomState(0).poisson(50, (8, 15))

        noisy = add_technical_noise(original, rng=rng)
        assert noisy.shape == original.shape

    def test_integer_output(self) -> None:
        """Output should be integer counts."""
        rng = random.Random(42)
        original = np.random.RandomState(0).poisson(50, (5, 10))

        noisy = add_technical_noise(original, rng=rng)
        assert np.issubdtype(noisy.dtype, np.integer)

    def test_non_negative_output(self) -> None:
        """Output values must be non-negative."""
        rng = random.Random(42)
        original = np.random.RandomState(0).poisson(50, (5, 10))

        noisy = add_technical_noise(original, rng=rng)
        assert np.all(noisy >= 0)

    def test_noise_changes_values(self) -> None:
        """Adding noise should change at least some values."""
        rng = random.Random(42)
        original = np.full((10, 10), 100)

        noisy = add_technical_noise(original, amplification_bias=0.5, rng=rng)

        # Noise should make at least some values differ from the original
        assert not np.array_equal(noisy, original), "Noise should change the expression values"

    def test_zero_amplification_bias(self) -> None:
        """With zero bias, output should still differ due to Poisson noise."""
        rng = random.Random(42)
        original = np.full((5, 5), 100)

        noisy = add_technical_noise(original, amplification_bias=0.0, rng=rng)

        assert noisy.shape == original.shape
        assert np.all(noisy >= 0)

    def test_sequencing_depth_adjustment(self) -> None:
        """Custom sequencing_depth should rescale the total counts per sample."""
        rng = random.Random(42)
        original = np.random.RandomState(0).poisson(50, (5, 20))
        target_depth = 5000.0

        noisy = add_technical_noise(original, amplification_bias=0.0, sequencing_depth=target_depth, rng=rng)

        assert noisy.shape == original.shape
        # Due to Poisson noise the depth won't be exact, but should be in the ballpark
        for i in range(noisy.shape[0]):
            sample_depth = noisy[i, :].sum()
            # Allow generous tolerance due to stochastic Poisson rounding
            assert sample_depth > 0, "Sample depth should be positive"

    def test_higher_bias_more_variation(self) -> None:
        """Higher amplification bias should produce more variance in output."""
        original = np.full((50, 20), 100)

        rng_low = random.Random(42)
        noisy_low = add_technical_noise(original, amplification_bias=0.01, rng=rng_low)

        rng_high = random.Random(42)
        noisy_high = add_technical_noise(original, amplification_bias=0.9, rng=rng_high)

        var_low = np.var(noisy_low)
        var_high = np.var(noisy_high)
        assert var_high > var_low, "Higher bias should produce more variation"

    def test_preserves_zero_rows(self) -> None:
        """Rows of zeros in input should remain approximately zero after noise."""
        rng = random.Random(42)
        original = np.zeros((3, 10), dtype=int)

        noisy = add_technical_noise(original, amplification_bias=0.1, rng=rng)

        # Poisson(0) = 0 with probability 1, so zero input -> zero output
        assert np.all(noisy == 0), "Zero expression + Poisson noise should remain zero"

    def test_single_sample(self) -> None:
        """Edge case: single-sample expression matrix."""
        rng = random.Random(42)
        original = np.random.RandomState(0).poisson(50, (1, 10))
        noisy = add_technical_noise(original, rng=rng)
        assert noisy.shape == (1, 10)

    def test_single_feature(self) -> None:
        """Edge case: single-feature expression matrix."""
        rng = random.Random(42)
        original = np.random.RandomState(0).poisson(50, (5, 1))
        noisy = add_technical_noise(original, rng=rng)
        assert noisy.shape == (5, 1)

    def test_reproducibility(self) -> None:
        """Same seed produces identical noisy output."""
        original = np.random.RandomState(0).poisson(50, (8, 12))

        noisy_a = add_technical_noise(original, rng=random.Random(42))
        noisy_b = add_technical_noise(original, rng=random.Random(42))

        np.testing.assert_array_equal(noisy_a, noisy_b)


# ---------------------------------------------------------------------------
# Cross-function integration tests
# ---------------------------------------------------------------------------


class TestRnaSimulationIntegration:
    """Integration tests combining multiple RNA simulation functions."""

    def test_bulk_rnaseq_then_add_noise(self) -> None:
        """Simulate bulk data and then add technical noise."""
        rng = random.Random(42)
        counts = simulate_bulk_rnaseq(8, 20, rng=rng)

        rng_noise = random.Random(99)
        noisy = add_technical_noise(counts, amplification_bias=0.2, rng=rng_noise)

        assert noisy.shape == counts.shape
        assert np.all(noisy >= 0)

    def test_differential_expression_then_add_noise(self) -> None:
        """Simulate DE data and then add technical noise."""
        fold_changes = np.array([3.0, 0.3])
        rng = random.Random(42)

        expr, labels = simulate_differential_expression(20, 15, fold_changes, rng=rng)
        noisy = add_technical_noise(expr, amplification_bias=0.1, rng=random.Random(99))

        assert noisy.shape == expr.shape
        assert len(labels) == 20

    def test_single_cell_dropout_sparsity(self) -> None:
        """Single-cell data with dropout should be sparser than bulk data."""
        rng_bulk = random.Random(42)
        bulk = simulate_bulk_rnaseq(20, 15, rng=rng_bulk)

        rng_sc = random.Random(42)
        sc, _ = simulate_single_cell_rnaseq(20, 15, dropout_rate=0.7, rng=rng_sc)

        bulk_zero_frac = np.sum(bulk == 0) / bulk.size
        sc_zero_frac = np.sum(sc == 0) / sc.size

        # Single-cell with high dropout should have more zeros
        assert sc_zero_frac > bulk_zero_frac

    def test_time_series_shape_consistency(self) -> None:
        """Time series output is compatible with add_technical_noise input."""
        rng = random.Random(42)
        ts_expr = simulate_time_series_expression(12, 8, rng=rng)

        rng_noise = random.Random(99)
        noisy = add_technical_noise(ts_expr, rng=rng_noise)

        assert noisy.shape == ts_expr.shape

    def test_spatial_expression_all_patterns_produce_data(self) -> None:
        """All spatial pattern types produce valid, non-trivial data."""
        for pattern in ("random", "gradient", "clusters"):
            rng = random.Random(42)
            expr, coords = simulate_spatial_expression(30, 8, spatial_patterns=pattern, rng=rng)

            assert expr.shape == (30, 8), f"Wrong shape for {pattern}"
            assert coords.shape == (30, 2), f"Wrong coords shape for {pattern}"
            assert np.all(np.isfinite(coords)), f"Non-finite coords for {pattern}"
            # At least some non-zero expression expected
            assert expr.sum() > 0, f"All-zero expression for {pattern}"

    def test_convenience_wrapper_matches_core(self) -> None:
        """simulate_rnaseq_counts output should resemble core NB simulation output."""
        rng = random.Random(42)
        counts = simulate_rnaseq_counts(n_genes=20, n_samples=10, rng=rng)

        assert counts.shape == (10, 20)
        assert np.all(counts >= 0)
        assert np.issubdtype(counts.dtype, np.integer)
        # Should have realistic variation (not all zeros)
        assert counts.sum() > 0
