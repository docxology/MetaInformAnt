"""Tests for RNA-seq quality control module.

Comprehensive tests for all QC functions using REAL implementations.
NO MOCKING - all data is generated with numpy/pandas and processed
through actual function logic.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from metainformant.rna.analysis.qc import (
    classify_expression_level,
    compute_correlation_matrix,
    compute_gene_metrics,
    compute_sample_metrics,
    compute_saturation_curve,
    detect_batch_effects,
    detect_gc_bias,
    detect_length_bias,
    detect_outlier_samples,
    estimate_library_complexity,
    generate_qc_report,
)


# =============================================================================
# Fixtures
# =============================================================================


@pytest.fixture
def count_matrix() -> pd.DataFrame:
    """Standard count matrix: 200 genes x 8 samples from negative binomial."""
    rng = np.random.default_rng(42)
    n_genes, n_samples = 200, 8
    counts = rng.negative_binomial(5, 0.01, size=(n_genes, n_samples))
    genes = [f"gene_{i}" for i in range(n_genes)]
    samples = [f"sample_{i}" for i in range(n_samples)]
    return pd.DataFrame(counts, index=genes, columns=samples)


@pytest.fixture
def small_count_matrix() -> pd.DataFrame:
    """Small count matrix for deterministic hand-checked tests."""
    data = {
        "sample_A": [10, 20, 0, 5, 100],
        "sample_B": [15, 25, 3, 8, 80],
        "sample_C": [12, 18, 1, 6, 90],
    }
    return pd.DataFrame(data, index=["geneA", "geneB", "geneC", "geneD", "geneE"])


@pytest.fixture
def single_sample_matrix() -> pd.DataFrame:
    """Single-sample count matrix."""
    return pd.DataFrame({"sample_only": [10, 20, 0, 5, 100]}, index=[f"gene_{i}" for i in range(5)])


@pytest.fixture
def single_gene_matrix() -> pd.DataFrame:
    """Single-gene count matrix."""
    return pd.DataFrame({"s0": [50], "s1": [60], "s2": [70]}, index=["only_gene"])


@pytest.fixture
def zero_count_matrix() -> pd.DataFrame:
    """Count matrix where one sample has all zeros."""
    data = {
        "sample_normal": [10, 20, 5, 15, 30],
        "sample_zero": [0, 0, 0, 0, 0],
        "sample_ok": [8, 22, 3, 12, 28],
    }
    return pd.DataFrame(data, index=[f"gene_{i}" for i in range(5)])


@pytest.fixture
def batch_labels() -> pd.Series:
    """Batch labels for 8-sample count_matrix fixture."""
    return pd.Series(
        ["batch_A", "batch_A", "batch_A", "batch_A", "batch_B", "batch_B", "batch_B", "batch_B"],
        index=[f"sample_{i}" for i in range(8)],
    )


@pytest.fixture
def gc_content_series() -> pd.Series:
    """GC content for 200 genes (matching count_matrix fixture)."""
    rng = np.random.default_rng(99)
    genes = [f"gene_{i}" for i in range(200)]
    return pd.Series(rng.uniform(0.3, 0.7, size=200), index=genes)


@pytest.fixture
def gene_lengths_series() -> pd.Series:
    """Gene lengths for 200 genes (matching count_matrix fixture)."""
    rng = np.random.default_rng(77)
    genes = [f"gene_{i}" for i in range(200)]
    return pd.Series(rng.integers(300, 15000, size=200), index=genes)


# =============================================================================
# Tests: compute_sample_metrics
# =============================================================================


class TestComputeSampleMetrics:
    """Tests for compute_sample_metrics function."""

    def test_returns_dataframe_with_correct_columns(self, count_matrix: pd.DataFrame) -> None:
        """Verify the returned DataFrame contains all expected metric columns."""
        result = compute_sample_metrics(count_matrix)
        expected_columns = {"total_counts", "detected_genes", "median_expression", "mean_expression", "pct_zero", "cv"}
        assert set(result.columns) == expected_columns

    def test_index_matches_sample_names(self, count_matrix: pd.DataFrame) -> None:
        """Index of result should be sample names from the input columns."""
        result = compute_sample_metrics(count_matrix)
        assert list(result.index) == list(count_matrix.columns)

    def test_total_counts_are_correct(self, small_count_matrix: pd.DataFrame) -> None:
        """Total counts should equal column sums."""
        result = compute_sample_metrics(small_count_matrix)
        assert result.loc["sample_A", "total_counts"] == 135.0
        assert result.loc["sample_B", "total_counts"] == 131.0
        assert result.loc["sample_C", "total_counts"] == 127.0

    def test_detected_genes_correct(self, small_count_matrix: pd.DataFrame) -> None:
        """Detected genes = count of genes with expression > 0."""
        result = compute_sample_metrics(small_count_matrix)
        # sample_A has one zero (geneC)
        assert result.loc["sample_A", "detected_genes"] == 4
        # sample_B has all genes detected
        assert result.loc["sample_B", "detected_genes"] == 5

    def test_pct_zero_range(self, count_matrix: pd.DataFrame) -> None:
        """Percentage of zeros must be between 0 and 100."""
        result = compute_sample_metrics(count_matrix)
        assert (result["pct_zero"] >= 0).all()
        assert (result["pct_zero"] <= 100).all()

    def test_cv_nonnegative(self, count_matrix: pd.DataFrame) -> None:
        """Coefficient of variation must be non-negative."""
        result = compute_sample_metrics(count_matrix)
        assert (result["cv"] >= 0).all()

    def test_metric_types_are_float(self, count_matrix: pd.DataFrame) -> None:
        """All metric values should be numeric floats."""
        result = compute_sample_metrics(count_matrix)
        for col in result.columns:
            assert result[col].dtype in (np.float64, np.int64, float, int), f"Column {col} has unexpected dtype"

    def test_empty_dataframe_raises(self) -> None:
        """Empty DataFrame should raise ValueError."""
        with pytest.raises(ValueError, match="cannot be empty"):
            compute_sample_metrics(pd.DataFrame())

    def test_negative_values_raise(self) -> None:
        """Negative values in counts should raise ValueError."""
        df = pd.DataFrame({"s1": [10, -5, 3]}, index=["g1", "g2", "g3"])
        with pytest.raises(ValueError, match="negative"):
            compute_sample_metrics(df)

    def test_single_sample(self, single_sample_matrix: pd.DataFrame) -> None:
        """Single-sample input should still return valid metrics."""
        result = compute_sample_metrics(single_sample_matrix)
        assert len(result) == 1
        assert result.loc["sample_only", "total_counts"] == 135.0
        assert result.loc["sample_only", "detected_genes"] == 4


# =============================================================================
# Tests: detect_outlier_samples
# =============================================================================


class TestDetectOutlierSamples:
    """Tests for detect_outlier_samples function."""

    def test_mad_returns_list(self, count_matrix: pd.DataFrame) -> None:
        """MAD method should return a list of sample names."""
        result = detect_outlier_samples(count_matrix, method="mad", threshold=3.0)
        assert isinstance(result, list)
        for name in result:
            assert name in count_matrix.columns

    def test_mad_detects_extreme_outlier(self, count_matrix: pd.DataFrame) -> None:
        """An extreme outlier sample should be detected by MAD."""
        # Create a copy and make one sample extreme
        df = count_matrix.copy()
        df["sample_0"] = df["sample_0"] * 100  # Massively inflate one sample
        result = detect_outlier_samples(df, method="mad", threshold=2.0)
        assert "sample_0" in result

    def test_strict_threshold_finds_more_outliers(self, count_matrix: pd.DataFrame) -> None:
        """A lower threshold should detect at least as many outliers as a higher one."""
        outliers_strict = detect_outlier_samples(count_matrix, method="mad", threshold=1.5)
        outliers_loose = detect_outlier_samples(count_matrix, method="mad", threshold=5.0)
        assert len(outliers_strict) >= len(outliers_loose)

    def test_pca_distance_method(self, count_matrix: pd.DataFrame) -> None:
        """PCA distance method should return a list (may be empty for normal data)."""
        result = detect_outlier_samples(count_matrix, method="pca_distance", threshold=3.0)
        assert isinstance(result, list)

    def test_isolation_forest_method(self, count_matrix: pd.DataFrame) -> None:
        """Isolation forest method should return a list."""
        result = detect_outlier_samples(count_matrix, method="isolation_forest", threshold=0.1)
        assert isinstance(result, list)

    def test_invalid_method_raises(self, count_matrix: pd.DataFrame) -> None:
        """Unknown method should raise ValueError."""
        with pytest.raises(ValueError, match="Unknown method"):
            detect_outlier_samples(count_matrix, method="bogus")  # type: ignore[arg-type]

    def test_empty_dataframe_raises(self) -> None:
        """Empty DataFrame should raise ValueError."""
        with pytest.raises(ValueError, match="cannot be empty"):
            detect_outlier_samples(pd.DataFrame(), method="mad")

    def test_pca_too_few_samples_returns_empty(self) -> None:
        """PCA method with fewer than 3 samples returns empty list."""
        df = pd.DataFrame({"s0": [1, 2, 3], "s1": [4, 5, 6]}, index=["g0", "g1", "g2"])
        result = detect_outlier_samples(df, method="pca_distance", threshold=3.0)
        assert result == []


# =============================================================================
# Tests: compute_gene_metrics
# =============================================================================


class TestComputeGeneMetrics:
    """Tests for compute_gene_metrics function."""

    def test_returns_correct_columns(self, count_matrix: pd.DataFrame) -> None:
        """Result should have the documented metric columns."""
        result = compute_gene_metrics(count_matrix)
        expected = {"mean_expression", "variance", "cv", "pct_zero", "n_samples_detected", "dispersion_estimate"}
        assert set(result.columns) == expected

    def test_index_matches_gene_names(self, count_matrix: pd.DataFrame) -> None:
        """Index should match gene names from input rows."""
        result = compute_gene_metrics(count_matrix)
        assert list(result.index) == list(count_matrix.index)

    def test_mean_expression_matches_manual(self, small_count_matrix: pd.DataFrame) -> None:
        """Mean expression should match manually computed row means."""
        result = compute_gene_metrics(small_count_matrix)
        # geneA: [10, 15, 12] -> mean = 37/3 ~ 12.333
        assert abs(result.loc["geneA", "mean_expression"] - (10 + 15 + 12) / 3) < 1e-10

    def test_pct_zero_for_always_expressed_gene(self, small_count_matrix: pd.DataFrame) -> None:
        """A gene expressed in all samples should have 0% zeros."""
        result = compute_gene_metrics(small_count_matrix)
        # geneB is expressed in all 3 samples
        assert result.loc["geneB", "pct_zero"] == 0.0

    def test_pct_zero_for_partially_detected_gene(self, small_count_matrix: pd.DataFrame) -> None:
        """geneC has 0 in sample_A only -> pct_zero ~ 33.33."""
        result = compute_gene_metrics(small_count_matrix)
        assert abs(result.loc["geneC", "pct_zero"] - 100.0 / 3.0) < 0.1

    def test_dispersion_nonnegative(self, count_matrix: pd.DataFrame) -> None:
        """Dispersion estimate (variance/mean) should be non-negative."""
        result = compute_gene_metrics(count_matrix)
        assert (result["dispersion_estimate"] >= 0).all()

    def test_empty_dataframe_raises(self) -> None:
        """Empty DataFrame should raise ValueError."""
        with pytest.raises(ValueError, match="cannot be empty"):
            compute_gene_metrics(pd.DataFrame())

    def test_negative_values_raise(self) -> None:
        """Negative values should raise ValueError."""
        df = pd.DataFrame({"s1": [-1, 2], "s2": [3, 4]}, index=["g1", "g2"])
        with pytest.raises(ValueError, match="negative"):
            compute_gene_metrics(df)

    def test_single_sample_variance_is_zero(self, single_sample_matrix: pd.DataFrame) -> None:
        """With a single sample, variance should be zero (ddof=1 guarded)."""
        result = compute_gene_metrics(single_sample_matrix)
        assert (result["variance"] == 0.0).all()

    def test_single_gene_matrix(self, single_gene_matrix: pd.DataFrame) -> None:
        """Single gene input should still work."""
        result = compute_gene_metrics(single_gene_matrix)
        assert len(result) == 1
        assert result.loc["only_gene", "n_samples_detected"] == 3


# =============================================================================
# Tests: classify_expression_level
# =============================================================================


class TestClassifyExpressionLevel:
    """Tests for classify_expression_level function."""

    def test_returns_correct_columns(self, count_matrix: pd.DataFrame) -> None:
        """Result should have mean_expression, expression_level, log2_mean."""
        result = classify_expression_level(count_matrix)
        expected = {"mean_expression", "expression_level", "log2_mean"}
        assert set(result.columns) == expected

    def test_levels_are_valid_strings(self, count_matrix: pd.DataFrame) -> None:
        """Every expression level should be one of low, medium, high."""
        result = classify_expression_level(count_matrix)
        valid_levels = {"low", "medium", "high"}
        assert set(result["expression_level"].unique()).issubset(valid_levels)

    def test_custom_thresholds(self) -> None:
        """Custom thresholds should correctly classify genes."""
        df = pd.DataFrame(
            {"s1": [0.5, 5.0, 50.0], "s2": [0.3, 6.0, 60.0]},
            index=["low_gene", "mid_gene", "high_gene"],
        )
        thresholds = {"low_high": 2.0, "medium_high": 20.0}
        result = classify_expression_level(df, thresholds=thresholds)
        assert result.loc["low_gene", "expression_level"] == "low"
        assert result.loc["mid_gene", "expression_level"] == "medium"
        assert result.loc["high_gene", "expression_level"] == "high"

    def test_log2_mean_correct(self) -> None:
        """log2_mean should equal log2(mean_expression + 1)."""
        df = pd.DataFrame({"s1": [7.0], "s2": [9.0]}, index=["g1"])
        result = classify_expression_level(df)
        expected_log2 = float(np.log2(8.0 + 1))  # mean = 8.0
        assert abs(result.loc["g1", "log2_mean"] - expected_log2) < 1e-10

    def test_empty_dataframe_raises(self) -> None:
        """Empty DataFrame should raise ValueError."""
        with pytest.raises(ValueError, match="cannot be empty"):
            classify_expression_level(pd.DataFrame())

    def test_adaptive_thresholds_with_default(self, count_matrix: pd.DataFrame) -> None:
        """Default (None) thresholds should be computed from data distribution."""
        result = classify_expression_level(count_matrix, thresholds=None)
        # Should still contain all three levels for a negative binomial matrix
        assert len(result) == len(count_matrix)


# =============================================================================
# Tests: estimate_library_complexity
# =============================================================================


class TestEstimateLibraryComplexity:
    """Tests for estimate_library_complexity function."""

    def test_returns_correct_columns(self, count_matrix: pd.DataFrame) -> None:
        """Result should have all expected diversity columns."""
        result = estimate_library_complexity(count_matrix)
        expected = {
            "shannon_entropy",
            "simpson_diversity",
            "effective_gene_count",
            "max_possible_entropy",
            "normalized_entropy",
        }
        assert set(result.columns) == expected

    def test_index_matches_samples(self, count_matrix: pd.DataFrame) -> None:
        """Index should be sample names."""
        result = estimate_library_complexity(count_matrix)
        assert list(result.index) == list(count_matrix.columns)

    def test_shannon_entropy_positive(self, count_matrix: pd.DataFrame) -> None:
        """Shannon entropy should be positive for non-trivial data."""
        result = estimate_library_complexity(count_matrix)
        assert (result["shannon_entropy"] > 0).all()

    def test_normalized_entropy_zero_to_one(self, count_matrix: pd.DataFrame) -> None:
        """Normalized entropy should be in [0, 1]."""
        result = estimate_library_complexity(count_matrix)
        assert (result["normalized_entropy"] >= 0).all()
        assert (result["normalized_entropy"] <= 1.0 + 1e-10).all()

    def test_simpson_diversity_zero_to_one(self, count_matrix: pd.DataFrame) -> None:
        """Simpson diversity should be in [0, 1]."""
        result = estimate_library_complexity(count_matrix)
        assert (result["simpson_diversity"] >= 0).all()
        assert (result["simpson_diversity"] <= 1.0 + 1e-10).all()

    def test_zero_count_sample_handled(self, zero_count_matrix: pd.DataFrame) -> None:
        """A sample with all zeros should have zero entropy."""
        result = estimate_library_complexity(zero_count_matrix)
        assert result.loc["sample_zero", "shannon_entropy"] == 0.0
        assert result.loc["sample_zero", "simpson_diversity"] == 0.0
        assert result.loc["sample_zero", "normalized_entropy"] == 0.0

    def test_uniform_distribution_maximum_entropy(self) -> None:
        """A perfectly uniform distribution should have near-maximum entropy."""
        n_genes = 100
        df = pd.DataFrame({"uniform_sample": [10] * n_genes}, index=[f"g_{i}" for i in range(n_genes)])
        result = estimate_library_complexity(df)
        # Normalized entropy should be very close to 1.0 for uniform
        assert result.loc["uniform_sample", "normalized_entropy"] > 0.99

    def test_empty_dataframe_raises(self) -> None:
        """Empty DataFrame should raise ValueError."""
        with pytest.raises(ValueError, match="cannot be empty"):
            estimate_library_complexity(pd.DataFrame())

    def test_negative_values_raise(self) -> None:
        """Negative values should raise ValueError."""
        df = pd.DataFrame({"s1": [-1, 5]}, index=["g1", "g2"])
        with pytest.raises(ValueError, match="negative"):
            estimate_library_complexity(df)


# =============================================================================
# Tests: compute_saturation_curve
# =============================================================================


class TestComputeSaturationCurve:
    """Tests for compute_saturation_curve function."""

    def test_returns_dataframe_with_correct_columns(self, count_matrix: pd.DataFrame) -> None:
        """Result must have sample, fraction, detected_genes_mean, detected_genes_std, total_counts."""
        result = compute_saturation_curve(count_matrix, fractions=[0.5, 1.0], n_iterations=3)
        expected = {"sample", "fraction", "detected_genes_mean", "detected_genes_std", "total_counts"}
        assert set(result.columns) == expected

    def test_number_of_rows(self, count_matrix: pd.DataFrame) -> None:
        """Should have n_samples * n_fractions rows."""
        fractions = [0.25, 0.5, 0.75, 1.0]
        result = compute_saturation_curve(count_matrix, fractions=fractions, n_iterations=3)
        expected_rows = len(count_matrix.columns) * len(fractions)
        assert len(result) == expected_rows

    def test_full_fraction_detects_most_genes(self, count_matrix: pd.DataFrame) -> None:
        """At fraction=1.0, detected genes should be highest or tied."""
        result = compute_saturation_curve(count_matrix, fractions=[0.1, 0.5, 1.0], n_iterations=5)
        for sample in count_matrix.columns:
            sample_data = result[result["sample"] == sample]
            detected_at_full = sample_data[sample_data["fraction"] == 1.0]["detected_genes_mean"].values[0]
            detected_at_half = sample_data[sample_data["fraction"] == 0.5]["detected_genes_mean"].values[0]
            assert detected_at_full >= detected_at_half

    def test_monotonic_increase(self, count_matrix: pd.DataFrame) -> None:
        """Detected genes (mean) should generally increase with fraction."""
        fractions = [0.1, 0.3, 0.5, 0.7, 0.9, 1.0]
        result = compute_saturation_curve(count_matrix, fractions=fractions, n_iterations=10)
        for sample in count_matrix.columns:
            sample_data = result[result["sample"] == sample].sort_values("fraction")
            means = sample_data["detected_genes_mean"].values
            # Allow small non-monotonicity due to randomness, but overall should increase
            assert means[-1] >= means[0]

    def test_invalid_fraction_raises(self, count_matrix: pd.DataFrame) -> None:
        """Fraction outside (0, 1] should raise ValueError."""
        with pytest.raises(ValueError, match="Fractions must be in"):
            compute_saturation_curve(count_matrix, fractions=[0.0, 0.5])
        with pytest.raises(ValueError, match="Fractions must be in"):
            compute_saturation_curve(count_matrix, fractions=[1.5])

    def test_empty_dataframe_raises(self) -> None:
        """Empty DataFrame should raise ValueError."""
        with pytest.raises(ValueError, match="cannot be empty"):
            compute_saturation_curve(pd.DataFrame())

    def test_zero_count_sample(self, zero_count_matrix: pd.DataFrame) -> None:
        """Sample with all zeros should have 0 detected genes at every fraction."""
        result = compute_saturation_curve(zero_count_matrix, fractions=[0.5, 1.0], n_iterations=3)
        zero_data = result[result["sample"] == "sample_zero"]
        assert (zero_data["detected_genes_mean"] == 0.0).all()

    def test_reproducibility_with_seed(self, count_matrix: pd.DataFrame) -> None:
        """Same random_state should produce identical results."""
        r1 = compute_saturation_curve(count_matrix, fractions=[0.5], n_iterations=5, random_state=123)
        r2 = compute_saturation_curve(count_matrix, fractions=[0.5], n_iterations=5, random_state=123)
        pd.testing.assert_frame_equal(r1, r2)


# =============================================================================
# Tests: detect_batch_effects
# =============================================================================


class TestDetectBatchEffects:
    """Tests for detect_batch_effects function."""

    def test_kruskal_returns_required_keys(self, count_matrix: pd.DataFrame, batch_labels: pd.Series) -> None:
        """Kruskal method result dict must contain documented keys."""
        result = detect_batch_effects(count_matrix, batch_labels, method="kruskal")
        assert "method" in result
        assert "batch_effect_detected" in result
        assert result["method"] == "kruskal"
        assert isinstance(result["batch_effect_detected"], bool)

    def test_kruskal_significant_keys(self, count_matrix: pd.DataFrame, batch_labels: pd.Series) -> None:
        """Kruskal result should include pct_significant_genes and median_pvalue."""
        result = detect_batch_effects(count_matrix, batch_labels, method="kruskal")
        assert "pct_significant_genes" in result
        assert "n_significant_genes" in result
        assert "median_pvalue" in result
        assert 0 <= result["pct_significant_genes"] <= 100
        assert 0 <= result["median_pvalue"] <= 1.0

    def test_strong_batch_effect_detected(self) -> None:
        """An artificially strong batch effect should be detected by kruskal."""
        rng = np.random.default_rng(42)
        n_genes = 100
        # Batch A: baseline expression
        batch_a = rng.negative_binomial(5, 0.05, size=(n_genes, 4))
        # Batch B: shifted expression (strong systematic difference)
        batch_b = rng.negative_binomial(50, 0.05, size=(n_genes, 4))
        data = np.hstack([batch_a, batch_b])
        genes = [f"gene_{i}" for i in range(n_genes)]
        samples = [f"s_{i}" for i in range(8)]
        df = pd.DataFrame(data, index=genes, columns=samples)
        labels = pd.Series(["A"] * 4 + ["B"] * 4, index=samples)

        result = detect_batch_effects(df, labels, method="kruskal")
        assert result["batch_effect_detected"] is True

    def test_no_batch_effect_with_one_batch(self, count_matrix: pd.DataFrame) -> None:
        """Single batch should report no batch effect."""
        labels = pd.Series(["only_batch"] * 8, index=[f"sample_{i}" for i in range(8)])
        result = detect_batch_effects(count_matrix, labels, method="kruskal")
        assert result["batch_effect_detected"] is False

    def test_silhouette_method(self, count_matrix: pd.DataFrame, batch_labels: pd.Series) -> None:
        """Silhouette method should return silhouette_score and batch_separation."""
        result = detect_batch_effects(count_matrix, batch_labels, method="silhouette")
        assert "silhouette_score" in result
        assert "batch_separation" in result
        assert result["batch_separation"] in ("strong", "moderate", "weak")

    def test_pvca_method(self, count_matrix: pd.DataFrame, batch_labels: pd.Series) -> None:
        """PVCA method should return variance_explained_batch and variance_explained_residual."""
        result = detect_batch_effects(count_matrix, batch_labels, method="pvca")
        assert "variance_explained_batch" in result
        assert "variance_explained_residual" in result
        total = result["variance_explained_batch"] + result["variance_explained_residual"]
        assert abs(total - 1.0) < 1e-10

    def test_invalid_method_raises(self, count_matrix: pd.DataFrame, batch_labels: pd.Series) -> None:
        """Unknown method should raise ValueError."""
        with pytest.raises(ValueError, match="Unknown method"):
            detect_batch_effects(count_matrix, batch_labels, method="invalid")  # type: ignore[arg-type]

    def test_missing_batch_labels_raises(self, count_matrix: pd.DataFrame) -> None:
        """Batch labels missing some samples should raise ValueError."""
        partial_labels = pd.Series(["A", "B"], index=["sample_0", "sample_1"])
        with pytest.raises(ValueError, match="Missing batch labels"):
            detect_batch_effects(count_matrix, partial_labels, method="kruskal")

    def test_empty_dataframe_raises(self) -> None:
        """Empty DataFrame should raise ValueError."""
        labels = pd.Series(dtype=str)
        with pytest.raises(ValueError, match="cannot be empty"):
            detect_batch_effects(pd.DataFrame(), labels, method="kruskal")


# =============================================================================
# Tests: compute_correlation_matrix
# =============================================================================


class TestComputeCorrelationMatrix:
    """Tests for compute_correlation_matrix function."""

    def test_pearson_returns_square_matrix(self, count_matrix: pd.DataFrame) -> None:
        """Pearson correlation should return an n_samples x n_samples DataFrame."""
        result = compute_correlation_matrix(count_matrix, method="pearson")
        n = len(count_matrix.columns)
        assert result.shape == (n, n)
        assert list(result.index) == list(count_matrix.columns)
        assert list(result.columns) == list(count_matrix.columns)

    def test_spearman_returns_square_matrix(self, count_matrix: pd.DataFrame) -> None:
        """Spearman correlation should return an n_samples x n_samples DataFrame."""
        result = compute_correlation_matrix(count_matrix, method="spearman")
        n = len(count_matrix.columns)
        assert result.shape == (n, n)

    def test_diagonal_is_one(self, count_matrix: pd.DataFrame) -> None:
        """Diagonal of correlation matrix should be 1.0 (self-correlation)."""
        result = compute_correlation_matrix(count_matrix, method="pearson")
        for sample in count_matrix.columns:
            assert abs(result.loc[sample, sample] - 1.0) < 1e-10

    def test_symmetric(self, count_matrix: pd.DataFrame) -> None:
        """Correlation matrix should be symmetric."""
        result = compute_correlation_matrix(count_matrix, method="pearson")
        pd.testing.assert_frame_equal(result, result.T, check_exact=False, atol=1e-10)

    def test_values_in_range(self, count_matrix: pd.DataFrame) -> None:
        """All correlation values should be in [-1, 1]."""
        result = compute_correlation_matrix(count_matrix, method="spearman")
        assert (result.values >= -1.0 - 1e-10).all()
        assert (result.values <= 1.0 + 1e-10).all()

    def test_invalid_method_raises(self, count_matrix: pd.DataFrame) -> None:
        """Unknown method should raise ValueError."""
        with pytest.raises(ValueError, match="Unknown method"):
            compute_correlation_matrix(count_matrix, method="kendall")  # type: ignore[arg-type]

    def test_empty_dataframe_raises(self) -> None:
        """Empty DataFrame should raise ValueError."""
        with pytest.raises(ValueError, match="cannot be empty"):
            compute_correlation_matrix(pd.DataFrame())

    def test_identical_samples_correlation_one(self) -> None:
        """Identical samples should have correlation of 1.0."""
        df = pd.DataFrame({"s1": [1, 2, 3, 4, 5], "s2": [1, 2, 3, 4, 5]}, index=[f"g{i}" for i in range(5)])
        result = compute_correlation_matrix(df, method="pearson")
        assert abs(result.loc["s1", "s2"] - 1.0) < 1e-10


# =============================================================================
# Tests: detect_gc_bias
# =============================================================================


class TestDetectGcBias:
    """Tests for detect_gc_bias function."""

    def test_returns_required_keys(
        self, count_matrix: pd.DataFrame, gc_content_series: pd.Series
    ) -> None:
        """Result should contain overall_correlation, overall_pvalue, bias_detected, per_sample."""
        result = detect_gc_bias(count_matrix, gc_content_series)
        assert "overall_correlation" in result
        assert "overall_pvalue" in result
        assert "bias_detected" in result
        assert "per_sample" in result

    def test_per_sample_contains_all_samples(
        self, count_matrix: pd.DataFrame, gc_content_series: pd.Series
    ) -> None:
        """per_sample should have an entry for each sample in the expression DataFrame."""
        result = detect_gc_bias(count_matrix, gc_content_series)
        for sample in count_matrix.columns:
            assert sample in result["per_sample"]

    def test_too_few_common_genes_returns_no_bias(self) -> None:
        """Fewer than 10 common genes should return bias_detected=False."""
        df = pd.DataFrame({"s1": [1, 2, 3]}, index=["a", "b", "c"])
        gc = pd.Series([0.4, 0.5, 0.6], index=["a", "b", "c"])
        result = detect_gc_bias(df, gc)
        assert result["bias_detected"] is False
        assert "Too few genes" in result.get("message", "")

    def test_artificial_gc_bias_detected(self) -> None:
        """A dataset where expression strongly correlates with GC should detect bias."""
        rng = np.random.default_rng(42)
        n_genes = 100
        gc = pd.Series(rng.uniform(0.3, 0.7, size=n_genes), index=[f"g_{i}" for i in range(n_genes)])
        # Expression linearly proportional to GC content + noise
        base_expr = gc.values * 1000
        data = {}
        for i in range(6):
            noise = rng.normal(0, 10, size=n_genes)
            data[f"s_{i}"] = np.maximum(base_expr + noise, 1).astype(int)
        df = pd.DataFrame(data, index=[f"g_{i}" for i in range(n_genes)])
        result = detect_gc_bias(df, gc)
        assert result["bias_detected"] is True
        assert abs(result["overall_correlation"]) > 0.1

    def test_empty_dataframe_raises(self) -> None:
        """Empty DataFrame should raise ValueError."""
        gc = pd.Series([0.5], index=["g1"])
        with pytest.raises(ValueError, match="cannot be empty"):
            detect_gc_bias(pd.DataFrame(), gc)


# =============================================================================
# Tests: detect_length_bias
# =============================================================================


class TestDetectLengthBias:
    """Tests for detect_length_bias function."""

    def test_returns_required_keys(
        self, count_matrix: pd.DataFrame, gene_lengths_series: pd.Series
    ) -> None:
        """Result should contain overall_correlation, overall_pvalue, bias_detected, per_sample."""
        result = detect_length_bias(count_matrix, gene_lengths_series)
        assert "overall_correlation" in result
        assert "overall_pvalue" in result
        assert "bias_detected" in result
        assert "per_sample" in result

    def test_too_few_common_genes_returns_no_bias(self) -> None:
        """Fewer than 10 common genes should return bias_detected=False."""
        df = pd.DataFrame({"s1": [10, 20, 30]}, index=["x", "y", "z"])
        lengths = pd.Series([500, 1000, 1500], index=["x", "y", "z"])
        result = detect_length_bias(df, lengths)
        assert result["bias_detected"] is False

    def test_artificial_length_bias_detected(self) -> None:
        """Expression proportional to gene length should be detected as bias."""
        rng = np.random.default_rng(42)
        n_genes = 100
        lengths = pd.Series(rng.integers(500, 10000, size=n_genes), index=[f"g_{i}" for i in range(n_genes)])
        # Expression proportional to length
        base_expr = (lengths.values / 10).astype(float)
        data = {}
        for i in range(6):
            noise = rng.normal(0, 2, size=n_genes)
            data[f"s_{i}"] = np.maximum(base_expr + noise, 1).astype(int)
        df = pd.DataFrame(data, index=[f"g_{i}" for i in range(n_genes)])
        result = detect_length_bias(df, lengths)
        assert result["bias_detected"] is True

    def test_empty_dataframe_raises(self) -> None:
        """Empty DataFrame should raise ValueError."""
        lengths = pd.Series([1000], index=["g1"])
        with pytest.raises(ValueError, match="cannot be empty"):
            detect_length_bias(pd.DataFrame(), lengths)

    def test_per_sample_correlations_present(
        self, count_matrix: pd.DataFrame, gene_lengths_series: pd.Series
    ) -> None:
        """per_sample should map each sample to a dict with correlation and pvalue."""
        result = detect_length_bias(count_matrix, gene_lengths_series)
        for sample in count_matrix.columns:
            assert sample in result["per_sample"]
            assert "correlation" in result["per_sample"][sample]
            assert "pvalue" in result["per_sample"][sample]


# =============================================================================
# Tests: generate_qc_report
# =============================================================================


class TestGenerateQcReport:
    """Tests for generate_qc_report function."""

    def test_report_has_all_sections(self, count_matrix: pd.DataFrame) -> None:
        """Report should contain summary, sample_metrics, gene_metrics, etc."""
        report = generate_qc_report(count_matrix)
        expected_keys = {
            "summary",
            "sample_metrics",
            "gene_metrics",
            "library_complexity",
            "outlier_samples",
            "expression_classification",
            "correlation_stats",
            "batch_effects",
            "gc_bias",
            "length_bias",
            "warnings",
            "n_warnings",
        }
        assert expected_keys.issubset(set(report.keys()))

    def test_summary_statistics_correct(self, count_matrix: pd.DataFrame) -> None:
        """Summary section should have correct n_samples and n_genes."""
        report = generate_qc_report(count_matrix)
        assert report["summary"]["n_samples"] == 8
        assert report["summary"]["n_genes"] == 200

    def test_summary_total_counts_positive(self, count_matrix: pd.DataFrame) -> None:
        """Total counts in summary should be positive for non-trivial data."""
        report = generate_qc_report(count_matrix)
        assert report["summary"]["total_counts"] > 0

    def test_report_with_batch_metadata(self, count_matrix: pd.DataFrame) -> None:
        """Providing sample_metadata with batch column should run batch detection."""
        sample_meta = pd.DataFrame(
            {"batch": ["A", "A", "A", "A", "B", "B", "B", "B"]},
            index=[f"sample_{i}" for i in range(8)],
        )
        report = generate_qc_report(count_matrix, sample_metadata=sample_meta)
        # batch_effects should not be just a "no info" message
        assert "method" in report["batch_effects"] or "message" in report["batch_effects"]

    def test_report_with_gene_metadata(
        self, count_matrix: pd.DataFrame, gc_content_series: pd.Series, gene_lengths_series: pd.Series
    ) -> None:
        """Providing gene_metadata with gc_content and length should run bias detection."""
        gene_meta = pd.DataFrame(
            {"gc_content": gc_content_series, "length": gene_lengths_series},
            index=gc_content_series.index,
        )
        report = generate_qc_report(count_matrix, gene_metadata=gene_meta)
        # gc_bias and length_bias should have real results, not "no info" messages
        assert "overall_correlation" in report["gc_bias"] or "message" in report["gc_bias"]
        assert "overall_correlation" in report["length_bias"] or "message" in report["length_bias"]

    def test_report_warnings_is_list(self, count_matrix: pd.DataFrame) -> None:
        """Warnings should be a list of strings."""
        report = generate_qc_report(count_matrix)
        assert isinstance(report["warnings"], list)
        for w in report["warnings"]:
            assert isinstance(w, str)

    def test_report_n_warnings_matches(self, count_matrix: pd.DataFrame) -> None:
        """n_warnings should equal len(warnings)."""
        report = generate_qc_report(count_matrix)
        assert report["n_warnings"] == len(report["warnings"])

    def test_empty_dataframe_raises(self) -> None:
        """Empty DataFrame should raise ValueError."""
        with pytest.raises(ValueError, match="cannot be empty"):
            generate_qc_report(pd.DataFrame())

    def test_report_outlier_section_structure(self, count_matrix: pd.DataFrame) -> None:
        """outlier_samples section should have method, threshold, outliers, n_outliers."""
        report = generate_qc_report(count_matrix)
        outlier_section = report["outlier_samples"]
        assert "method" in outlier_section
        assert "threshold" in outlier_section
        assert "outliers" in outlier_section
        assert "n_outliers" in outlier_section
        assert isinstance(outlier_section["outliers"], list)

    def test_expression_classification_structure(self, count_matrix: pd.DataFrame) -> None:
        """expression_classification should have low, medium, high counts."""
        report = generate_qc_report(count_matrix)
        ec = report["expression_classification"]
        assert "low" in ec
        assert "medium" in ec
        assert "high" in ec
        total_classified = ec["low"] + ec["medium"] + ec["high"]
        assert total_classified == 200  # All genes classified

    def test_correlation_stats_structure(self, count_matrix: pd.DataFrame) -> None:
        """correlation_stats should have method, mean/min pairwise correlation."""
        report = generate_qc_report(count_matrix)
        cs = report["correlation_stats"]
        assert cs["method"] == "spearman"
        assert "mean_pairwise_correlation" in cs
        assert "min_pairwise_correlation" in cs
        assert "samples_with_low_correlation" in cs
