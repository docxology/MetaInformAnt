"""Comprehensive tests for the RNA-seq expression analysis module.

Tests normalization, differential expression, PCA, gene filtering,
and visualization data preparation using REAL implementations with
realistic RNA-seq count matrices. NO MOCKING.
"""

from __future__ import annotations

from typing import List

import numpy as np
import pandas as pd
import pytest

from metainformant.rna.analysis.expression import (
    adjust_pvalues,
    compute_sample_distances,
    differential_expression,
    estimate_size_factors,
    filter_low_expression,
    get_highly_variable_genes,
    normalize_counts,
    pca_analysis,
    prepare_ma_data,
    prepare_volcano_data,
)


# =============================================================================
# Fixtures - Realistic RNA-seq data generated with numpy
# =============================================================================


@pytest.fixture
def count_matrix() -> pd.DataFrame:
    """Realistic gene expression count matrix using negative binomial distribution."""
    rng = np.random.default_rng(42)
    n_genes, n_samples = 100, 6
    counts = rng.negative_binomial(5, 0.01, size=(n_genes, n_samples))
    genes = [f"gene_{i}" for i in range(n_genes)]
    samples = [f"sample_{i}" for i in range(n_samples)]
    return pd.DataFrame(counts, index=genes, columns=samples)


@pytest.fixture
def small_count_matrix() -> pd.DataFrame:
    """Small deterministic count matrix for exact value verification."""
    return pd.DataFrame(
        {
            "ctrl_1": [100, 200, 50, 0, 500],
            "ctrl_2": [110, 180, 60, 0, 480],
            "ctrl_3": [90, 210, 55, 0, 520],
            "treat_1": [300, 100, 40, 0, 250],
            "treat_2": [280, 110, 35, 0, 260],
            "treat_3": [320, 90, 45, 0, 240],
        },
        index=["gene_A", "gene_B", "gene_C", "gene_zero", "gene_D"],
    )


@pytest.fixture
def conditions_6sample() -> List[str]:
    """Condition labels for a 6-sample experiment (3 control, 3 treatment)."""
    return ["control", "control", "control", "treatment", "treatment", "treatment"]


@pytest.fixture
def gene_lengths() -> pd.Series:
    """Gene lengths in base pairs for TPM/RPKM normalization."""
    rng = np.random.default_rng(99)
    n_genes = 100
    lengths = rng.integers(500, 10000, size=n_genes)
    return pd.Series(lengths, index=[f"gene_{i}" for i in range(n_genes)], dtype=float)


@pytest.fixture
def small_gene_lengths() -> pd.Series:
    """Gene lengths matching small_count_matrix."""
    return pd.Series(
        [1000, 2000, 500, 1500, 3000],
        index=["gene_A", "gene_B", "gene_C", "gene_zero", "gene_D"],
        dtype=float,
    )


@pytest.fixture
def de_results() -> pd.DataFrame:
    """Simulated differential expression results for visualization tests."""
    rng = np.random.default_rng(123)
    n_genes = 50
    return pd.DataFrame(
        {
            "gene": [f"gene_{i}" for i in range(n_genes)],
            "log2_fold_change": rng.normal(0, 2, n_genes),
            "p_value": rng.uniform(0, 1, n_genes),
            "adjusted_p_value": rng.uniform(0, 1, n_genes),
            "base_mean": rng.uniform(10, 10000, n_genes),
            "stat": rng.normal(0, 3, n_genes),
        }
    )


# =============================================================================
# Tests: normalize_counts
# =============================================================================


class TestNormalizeCPM:
    def test_cpm_sums_to_one_million(self, count_matrix: pd.DataFrame) -> None:
        result = normalize_counts(count_matrix, method="cpm")
        column_sums = result.sum(axis=0)
        for col_sum in column_sums:
            assert abs(col_sum - 1e6) < 1.0, f"CPM column sum {col_sum} not close to 1e6"

    def test_cpm_preserves_shape(self, count_matrix: pd.DataFrame) -> None:
        result = normalize_counts(count_matrix, method="cpm")
        assert result.shape == count_matrix.shape
        assert list(result.index) == list(count_matrix.index)
        assert list(result.columns) == list(count_matrix.columns)

    def test_cpm_values_nonnegative(self, count_matrix: pd.DataFrame) -> None:
        result = normalize_counts(count_matrix, method="cpm")
        assert (result.values >= 0).all()

    def test_cpm_proportional_to_counts(self, small_count_matrix: pd.DataFrame) -> None:
        result = normalize_counts(small_count_matrix, method="cpm")
        # gene_A in ctrl_1 has 100 counts out of 850 total => 100/850 * 1e6
        lib_size_ctrl1 = small_count_matrix["ctrl_1"].sum()
        expected_cpm = 100.0 / lib_size_ctrl1 * 1e6
        assert abs(result.loc["gene_A", "ctrl_1"] - expected_cpm) < 0.01


class TestNormalizeTPM:
    def test_tpm_requires_gene_lengths(self, count_matrix: pd.DataFrame) -> None:
        with pytest.raises(ValueError, match="gene_lengths required"):
            normalize_counts(count_matrix, method="tpm")

    def test_tpm_sums_to_one_million(self, count_matrix: pd.DataFrame, gene_lengths: pd.Series) -> None:
        result = normalize_counts(count_matrix, method="tpm", gene_lengths=gene_lengths)
        column_sums = result.sum(axis=0)
        for col_sum in column_sums:
            assert abs(col_sum - 1e6) < 1.0, f"TPM column sum {col_sum} not close to 1e6"

    def test_tpm_preserves_shape(self, count_matrix: pd.DataFrame, gene_lengths: pd.Series) -> None:
        result = normalize_counts(count_matrix, method="tpm", gene_lengths=gene_lengths)
        assert result.shape == count_matrix.shape

    def test_tpm_accounts_for_gene_length(self, small_count_matrix: pd.DataFrame) -> None:
        lengths = pd.Series(
            [1000, 2000, 1000, 1000, 1000],
            index=small_count_matrix.index,
            dtype=float,
        )
        result = normalize_counts(small_count_matrix, method="tpm", gene_lengths=lengths)
        # gene_B (length 2000) should have lower TPM relative to its counts
        # compared to gene_A (length 1000) in the same sample
        ratio_counts = small_count_matrix.loc["gene_B", "ctrl_1"] / small_count_matrix.loc["gene_A", "ctrl_1"]
        ratio_tpm = result.loc["gene_B", "ctrl_1"] / result.loc["gene_A", "ctrl_1"]
        assert ratio_tpm < ratio_counts, "TPM should reduce values for longer genes"

    def test_tpm_handles_missing_gene_lengths(self, count_matrix: pd.DataFrame) -> None:
        partial_lengths = pd.Series(
            [1000.0] * 50,
            index=[f"gene_{i}" for i in range(50)],
        )
        result = normalize_counts(count_matrix, method="tpm", gene_lengths=partial_lengths)
        assert result.shape == count_matrix.shape
        assert not result.isna().any().any()


class TestNormalizeRPKM:
    def test_rpkm_requires_gene_lengths(self, count_matrix: pd.DataFrame) -> None:
        with pytest.raises(ValueError, match="gene_lengths required"):
            normalize_counts(count_matrix, method="rpkm")

    def test_rpkm_preserves_shape(self, count_matrix: pd.DataFrame, gene_lengths: pd.Series) -> None:
        result = normalize_counts(count_matrix, method="rpkm", gene_lengths=gene_lengths)
        assert result.shape == count_matrix.shape

    def test_rpkm_values_nonnegative(self, count_matrix: pd.DataFrame, gene_lengths: pd.Series) -> None:
        result = normalize_counts(count_matrix, method="rpkm", gene_lengths=gene_lengths)
        assert (result.values >= 0).all()

    def test_rpkm_formula_correct(self, small_count_matrix: pd.DataFrame, small_gene_lengths: pd.Series) -> None:
        result = normalize_counts(small_count_matrix, method="rpkm", gene_lengths=small_gene_lengths)
        # RPKM = (count / (gene_length/1000)) / (library_size) * 1e6
        lib_size = small_count_matrix["ctrl_1"].sum()
        gene_len_kb = small_gene_lengths["gene_A"] / 1000.0
        expected = (100.0 / gene_len_kb) / lib_size * 1e6
        assert abs(result.loc["gene_A", "ctrl_1"] - expected) < 0.01


class TestNormalizeLog2CPM:
    def test_log2cpm_values_nonnegative(self, count_matrix: pd.DataFrame) -> None:
        result = normalize_counts(count_matrix, method="log2cpm")
        # log2(cpm + 1) is always >= 0 since cpm >= 0
        assert (result.values >= 0).all()

    def test_log2cpm_compresses_dynamic_range(self, count_matrix: pd.DataFrame) -> None:
        cpm = normalize_counts(count_matrix, method="cpm")
        log2cpm = normalize_counts(count_matrix, method="log2cpm")
        # log2 should compress the range
        cpm_range = cpm.values.max() - cpm.values.min()
        log_range = log2cpm.values.max() - log2cpm.values.min()
        assert log_range < cpm_range

    def test_log2cpm_pseudocount_prevents_negative_inf(self) -> None:
        counts = pd.DataFrame({"s1": [0, 10, 100], "s2": [0, 20, 200]}, index=["g1", "g2", "g3"])
        result = normalize_counts(counts, method="log2cpm")
        assert np.isfinite(result.values).all()


class TestNormalizeQuantile:
    def test_quantile_produces_equal_distributions(self, count_matrix: pd.DataFrame) -> None:
        result = normalize_counts(count_matrix, method="quantile")
        # After quantile normalization, sorted values should be very similar across samples
        sorted_vals = np.sort(result.values, axis=0)
        # Standard deviation across samples at each rank position should be small
        row_stds = sorted_vals.std(axis=1)
        mean_std = row_stds.mean()
        assert mean_std < 1.0, "Quantile normalization should produce nearly identical distributions"

    def test_quantile_preserves_shape(self, count_matrix: pd.DataFrame) -> None:
        result = normalize_counts(count_matrix, method="quantile")
        assert result.shape == count_matrix.shape


class TestNormalizeMedianRatio:
    def test_median_ratio_preserves_shape(self, count_matrix: pd.DataFrame) -> None:
        result = normalize_counts(count_matrix, method="median_ratio")
        assert result.shape == count_matrix.shape

    def test_median_ratio_adjusts_for_depth(self) -> None:
        # Sample 2 has exactly 2x the counts of sample 1
        counts = pd.DataFrame(
            {"s1": [100, 200, 300, 400], "s2": [200, 400, 600, 800]},
            index=["g1", "g2", "g3", "g4"],
        )
        result = normalize_counts(counts, method="median_ratio")
        # After normalization, values should be similar between samples
        ratio = result["s2"].mean() / result["s1"].mean()
        assert abs(ratio - 1.0) < 0.2, f"Median ratio normalization should equalize depth, got ratio {ratio}"


class TestNormalizeEdgeCases:
    def test_empty_dataframe_returns_empty(self) -> None:
        empty = pd.DataFrame()
        result = normalize_counts(empty, method="cpm")
        assert result.empty

    def test_negative_values_raise_error(self) -> None:
        counts = pd.DataFrame({"s1": [-1, 10], "s2": [5, 20]}, index=["g1", "g2"])
        with pytest.raises(ValueError, match="negative values"):
            normalize_counts(counts, method="cpm")

    def test_all_zero_counts(self) -> None:
        counts = pd.DataFrame({"s1": [0, 0, 0], "s2": [0, 0, 0]}, index=["g1", "g2", "g3"])
        result = normalize_counts(counts, method="cpm")
        assert result.shape == counts.shape
        # Division by zero handled: library_size = 0 replaced with 1
        assert np.isfinite(result.values).all()

    def test_single_gene(self) -> None:
        counts = pd.DataFrame({"s1": [500], "s2": [1000]}, index=["gene_only"])
        result = normalize_counts(counts, method="cpm")
        # Single gene is 100% of library => 1e6
        assert abs(result.loc["gene_only", "s1"] - 1e6) < 1.0
        assert abs(result.loc["gene_only", "s2"] - 1e6) < 1.0

    def test_single_sample(self) -> None:
        counts = pd.DataFrame({"s1": [100, 200, 300]}, index=["g1", "g2", "g3"])
        result = normalize_counts(counts, method="cpm")
        assert abs(result["s1"].sum() - 1e6) < 1.0

    def test_unknown_method_raises_error(self, count_matrix: pd.DataFrame) -> None:
        with pytest.raises(ValueError, match="Unknown normalization method"):
            normalize_counts(count_matrix, method="invalid_method")  # type: ignore[arg-type]


# =============================================================================
# Tests: estimate_size_factors
# =============================================================================


class TestEstimateSizeFactors:
    def test_size_factors_close_to_one_for_balanced(self) -> None:
        # Balanced library sizes should give size factors near 1
        counts = pd.DataFrame(
            {
                "s1": [100, 200, 300, 400],
                "s2": [105, 195, 310, 390],
                "s3": [98, 205, 295, 405],
            },
            index=["g1", "g2", "g3", "g4"],
        )
        sf = estimate_size_factors(counts)
        assert len(sf) == 3
        for val in sf.values:
            assert 0.8 < val < 1.2, f"Size factor {val} too far from 1.0 for balanced data"

    def test_size_factors_detect_depth_difference(self) -> None:
        counts = pd.DataFrame(
            {"s1": [100, 200, 300], "s2": [200, 400, 600]},
            index=["g1", "g2", "g3"],
        )
        sf = estimate_size_factors(counts)
        # s2 has 2x depth, so sf["s2"] should be about 2x sf["s1"]
        ratio = sf["s2"] / sf["s1"]
        assert 1.5 < ratio < 2.5, f"Size factor ratio {ratio} should be ~2.0"

    def test_size_factors_empty_input(self) -> None:
        sf = estimate_size_factors(pd.DataFrame())
        assert len(sf) == 0

    def test_size_factors_negative_values_raise(self) -> None:
        counts = pd.DataFrame({"s1": [-5, 10], "s2": [5, 20]}, index=["g1", "g2"])
        with pytest.raises(ValueError, match="negative values"):
            estimate_size_factors(counts)

    def test_size_factors_all_zeros_fallback(self) -> None:
        counts = pd.DataFrame({"s1": [0, 0], "s2": [0, 0]}, index=["g1", "g2"])
        sf = estimate_size_factors(counts)
        # Should fall back to library size normalization or return 1.0
        assert len(sf) == 2
        assert np.isfinite(sf.values).all()

    def test_size_factors_return_type(self, count_matrix: pd.DataFrame) -> None:
        sf = estimate_size_factors(count_matrix)
        assert isinstance(sf, pd.Series)
        assert len(sf) == count_matrix.shape[1]
        assert sf.index.tolist() == count_matrix.columns.tolist()


# =============================================================================
# Tests: differential_expression
# =============================================================================


class TestDifferentialExpressionDeseq2Like:
    def test_deseq2_like_returns_expected_columns(
        self, small_count_matrix: pd.DataFrame, conditions_6sample: List[str]
    ) -> None:
        result = differential_expression(small_count_matrix, conditions_6sample, method="deseq2_like")
        required_cols = {"gene", "log2_fold_change", "p_value", "adjusted_p_value", "base_mean", "stat"}
        assert required_cols.issubset(set(result.columns))

    def test_deseq2_like_detects_upregulation(
        self, small_count_matrix: pd.DataFrame, conditions_6sample: List[str]
    ) -> None:
        result = differential_expression(small_count_matrix, conditions_6sample, method="deseq2_like")
        # gene_A: ctrl ~100, treat ~300 => should be upregulated (positive log2FC)
        gene_a = result[result["gene"] == "gene_A"]
        if len(gene_a) > 0:
            assert gene_a["log2_fold_change"].values[0] > 0, "gene_A should be upregulated in treatment"

    def test_deseq2_like_detects_downregulation(
        self, small_count_matrix: pd.DataFrame, conditions_6sample: List[str]
    ) -> None:
        result = differential_expression(small_count_matrix, conditions_6sample, method="deseq2_like")
        # gene_B: ctrl ~200, treat ~100 => should be downregulated (negative log2FC)
        gene_b = result[result["gene"] == "gene_B"]
        if len(gene_b) > 0:
            assert gene_b["log2_fold_change"].values[0] < 0, "gene_B should be downregulated in treatment"

    def test_deseq2_like_pvalues_valid(self, small_count_matrix: pd.DataFrame, conditions_6sample: List[str]) -> None:
        result = differential_expression(small_count_matrix, conditions_6sample, method="deseq2_like")
        assert (result["p_value"] >= 0).all()
        assert (result["p_value"] <= 1).all()
        assert (result["adjusted_p_value"] >= 0).all()
        assert (result["adjusted_p_value"] <= 1).all()


class TestDifferentialExpressionTtest:
    def test_ttest_returns_expected_columns(
        self, small_count_matrix: pd.DataFrame, conditions_6sample: List[str]
    ) -> None:
        result = differential_expression(small_count_matrix, conditions_6sample, method="ttest")
        required_cols = {"gene", "log2_fold_change", "p_value", "adjusted_p_value", "base_mean", "stat"}
        assert required_cols.issubset(set(result.columns))

    def test_ttest_fold_change_direction(self, small_count_matrix: pd.DataFrame, conditions_6sample: List[str]) -> None:
        result = differential_expression(small_count_matrix, conditions_6sample, method="ttest")
        gene_a = result[result["gene"] == "gene_A"]
        if len(gene_a) > 0:
            assert gene_a["log2_fold_change"].values[0] > 0

    def test_ttest_sorted_by_adjusted_pvalue(
        self, small_count_matrix: pd.DataFrame, conditions_6sample: List[str]
    ) -> None:
        result = differential_expression(small_count_matrix, conditions_6sample, method="ttest")
        adj_pvals = result["adjusted_p_value"].values
        # Should be sorted ascending
        assert all(adj_pvals[i] <= adj_pvals[i + 1] for i in range(len(adj_pvals) - 1))


class TestDifferentialExpressionWilcoxon:
    def test_wilcoxon_returns_expected_columns(
        self, small_count_matrix: pd.DataFrame, conditions_6sample: List[str]
    ) -> None:
        result = differential_expression(small_count_matrix, conditions_6sample, method="wilcoxon")
        required_cols = {"gene", "log2_fold_change", "p_value", "adjusted_p_value", "base_mean", "stat"}
        assert required_cols.issubset(set(result.columns))

    def test_wilcoxon_nonparametric_robust(self) -> None:
        # Create data with an outlier; Wilcoxon should be robust
        counts = pd.DataFrame(
            {
                "ctrl_1": [100, 200],
                "ctrl_2": [110, 180],
                "ctrl_3": [9000, 210],
                "treat_1": [300, 50],
                "treat_2": [280, 60],
                "treat_3": [320, 55],
            },
            index=["gene_A", "gene_B"],
        )
        conditions = ["control", "control", "control", "treatment", "treatment", "treatment"]
        result = differential_expression(counts, conditions, method="wilcoxon")
        assert len(result) > 0
        assert (result["p_value"] >= 0).all()
        assert (result["p_value"] <= 1).all()


class TestDifferentialExpressionEdgeCases:
    def test_empty_count_matrix(self) -> None:
        result = differential_expression(pd.DataFrame(), ["a", "b"])
        assert result.empty
        expected_cols = {"gene", "log2_fold_change", "p_value", "adjusted_p_value", "base_mean", "stat"}
        assert expected_cols.issubset(set(result.columns))

    def test_wrong_number_of_conditions_raises(self, small_count_matrix: pd.DataFrame) -> None:
        three_conditions = ["a", "a", "b", "b", "c", "c"]
        with pytest.raises(ValueError, match="Expected exactly 2 conditions"):
            differential_expression(small_count_matrix, three_conditions)

    def test_mismatched_conditions_length_raises(self, small_count_matrix: pd.DataFrame) -> None:
        wrong_length = ["a", "b"]
        with pytest.raises(ValueError):
            differential_expression(small_count_matrix, wrong_length)

    def test_reference_condition_parameter(
        self, small_count_matrix: pd.DataFrame, conditions_6sample: List[str]
    ) -> None:
        result_default = differential_expression(small_count_matrix, conditions_6sample, method="ttest")
        result_flipped = differential_expression(
            small_count_matrix, conditions_6sample, method="ttest", reference="treatment"
        )
        # Fold changes should be opposite when reference is flipped
        gene_a_default = result_default[result_default["gene"] == "gene_A"]["log2_fold_change"].values
        gene_a_flipped = result_flipped[result_flipped["gene"] == "gene_A"]["log2_fold_change"].values
        if len(gene_a_default) > 0 and len(gene_a_flipped) > 0:
            assert gene_a_default[0] * gene_a_flipped[0] < 0, "Flipping reference should negate fold change"

    def test_min_count_filtering(self) -> None:
        counts = pd.DataFrame(
            {
                "ctrl_1": [5, 100],
                "ctrl_2": [3, 120],
                "treat_1": [4, 200],
                "treat_2": [6, 180],
            },
            index=["low_gene", "high_gene"],
        )
        conditions = ["control", "control", "treatment", "treatment"]
        result = differential_expression(counts, conditions, method="ttest", min_count=50)
        # low_gene total = 18, should be filtered out with min_count=50
        assert "low_gene" not in result["gene"].values

    def test_unknown_method_raises(self, small_count_matrix: pd.DataFrame, conditions_6sample: List[str]) -> None:
        with pytest.raises(ValueError, match="Unknown DE method"):
            differential_expression(small_count_matrix, conditions_6sample, method="invalid")  # type: ignore[arg-type]

    def test_conditions_as_series(self, small_count_matrix: pd.DataFrame) -> None:
        conditions = pd.Series(
            ["control", "control", "control", "treatment", "treatment", "treatment"],
            index=small_count_matrix.columns,
        )
        result = differential_expression(small_count_matrix, conditions, method="ttest")
        assert len(result) > 0


# =============================================================================
# Tests: adjust_pvalues
# =============================================================================


class TestAdjustPvalues:
    def test_bh_adjustment_increases_pvalues(self) -> None:
        pvals = np.array([0.001, 0.01, 0.03, 0.05, 0.1])
        adjusted = adjust_pvalues(pvals, method="bh")
        # BH adjusted p-values should be >= raw p-values
        assert (adjusted >= pvals - 1e-10).all()

    def test_bh_adjustment_capped_at_one(self) -> None:
        pvals = np.array([0.5, 0.6, 0.8, 0.9])
        adjusted = adjust_pvalues(pvals, method="bh")
        assert (adjusted <= 1.0).all()

    def test_bonferroni_multiplies_by_n(self) -> None:
        pvals = np.array([0.01, 0.02, 0.03])
        adjusted = adjust_pvalues(pvals, method="bonferroni")
        expected = np.clip(pvals * 3, 0, 1)
        np.testing.assert_allclose(adjusted, expected, atol=1e-10)

    def test_bonferroni_capped_at_one(self) -> None:
        pvals = np.array([0.5, 0.9])
        adjusted = adjust_pvalues(pvals, method="bonferroni")
        assert (adjusted <= 1.0).all()

    def test_fdr_alias_for_bh(self) -> None:
        pvals = np.array([0.001, 0.01, 0.05])
        bh_adjusted = adjust_pvalues(pvals, method="bh")
        fdr_adjusted = adjust_pvalues(pvals, method="fdr")
        np.testing.assert_array_equal(bh_adjusted, fdr_adjusted)

    def test_empty_pvalues(self) -> None:
        result = adjust_pvalues(np.array([]))
        assert len(result) == 0

    def test_single_pvalue(self) -> None:
        result = adjust_pvalues(np.array([0.05]), method="bh")
        assert len(result) == 1
        assert 0.0 <= result[0] <= 1.0

    def test_nan_handling(self) -> None:
        pvals = np.array([0.01, np.nan, 0.05])
        result = adjust_pvalues(pvals, method="bh")
        assert np.isnan(result[1])
        assert np.isfinite(result[0])
        assert np.isfinite(result[2])

    def test_unknown_method_raises(self) -> None:
        with pytest.raises(ValueError, match="Unknown p-value adjustment method"):
            adjust_pvalues(np.array([0.05]), method="invalid_method")  # type: ignore[arg-type]

    def test_bh_monotonicity(self) -> None:
        # BH adjusted p-values should be monotonically non-decreasing when sorted by raw p-value
        rng = np.random.default_rng(77)
        pvals = rng.uniform(0, 0.1, size=20)
        adjusted = adjust_pvalues(pvals, method="bh")
        sorted_indices = np.argsort(pvals)
        sorted_adjusted = adjusted[sorted_indices]
        # Monotonicity: each value >= previous
        for i in range(1, len(sorted_adjusted)):
            assert sorted_adjusted[i] >= sorted_adjusted[i - 1] - 1e-10


# =============================================================================
# Tests: pca_analysis
# =============================================================================


class TestPCAAnalysis:
    def test_pca_returns_correct_keys(self, count_matrix: pd.DataFrame) -> None:
        log_counts = np.log2(count_matrix + 1)
        result = pca_analysis(log_counts, n_components=2)
        assert "transformed" in result
        assert "explained_variance_ratio" in result
        assert "loadings" in result
        assert "components" in result

    def test_pca_transformed_shape(self, count_matrix: pd.DataFrame) -> None:
        log_counts = np.log2(count_matrix + 1)
        result = pca_analysis(log_counts, n_components=3)
        # transformed: samples x components
        assert result["transformed"].shape == (6, 3)

    def test_pca_loadings_shape(self, count_matrix: pd.DataFrame) -> None:
        log_counts = np.log2(count_matrix + 1)
        result = pca_analysis(log_counts, n_components=2)
        # loadings: genes x components
        assert result["loadings"].shape == (100, 2)

    def test_pca_explained_variance_sums_lte_one(self, count_matrix: pd.DataFrame) -> None:
        log_counts = np.log2(count_matrix + 1)
        result = pca_analysis(log_counts, n_components=2)
        total = result["explained_variance_ratio"].sum()
        assert total <= 1.0 + 1e-6, f"Explained variance {total} should not exceed 1.0"

    def test_pca_explained_variance_nonnegative(self, count_matrix: pd.DataFrame) -> None:
        log_counts = np.log2(count_matrix + 1)
        result = pca_analysis(log_counts, n_components=2)
        assert (result["explained_variance_ratio"] >= 0).all()

    def test_pca_pc1_explains_most_variance(self, count_matrix: pd.DataFrame) -> None:
        log_counts = np.log2(count_matrix + 1)
        result = pca_analysis(log_counts, n_components=3)
        evr = result["explained_variance_ratio"]
        # PC1 should explain more variance than PC2, PC2 more than PC3
        assert evr[0] >= evr[1], "PC1 should explain >= variance than PC2"

    def test_pca_column_names(self, count_matrix: pd.DataFrame) -> None:
        log_counts = np.log2(count_matrix + 1)
        result = pca_analysis(log_counts, n_components=2)
        assert list(result["transformed"].columns) == ["PC1", "PC2"]
        assert list(result["loadings"].columns) == ["PC1", "PC2"]

    def test_pca_sample_names_preserved(self, count_matrix: pd.DataFrame) -> None:
        log_counts = np.log2(count_matrix + 1)
        result = pca_analysis(log_counts, n_components=2)
        assert list(result["transformed"].index) == list(count_matrix.columns)

    def test_pca_gene_names_preserved(self, count_matrix: pd.DataFrame) -> None:
        log_counts = np.log2(count_matrix + 1)
        result = pca_analysis(log_counts, n_components=2)
        assert list(result["loadings"].index) == list(count_matrix.index)

    def test_pca_empty_input(self) -> None:
        result = pca_analysis(pd.DataFrame())
        assert result["transformed"].empty
        assert len(result["explained_variance_ratio"]) == 0

    def test_pca_caps_components_at_max(self) -> None:
        # 3 samples, 5 genes => max components = 3
        counts = pd.DataFrame(
            np.random.default_rng(42).integers(0, 100, (5, 3)).astype(float),
            index=[f"g{i}" for i in range(5)],
            columns=[f"s{i}" for i in range(3)],
        )
        result = pca_analysis(counts, n_components=10)
        n_components_actual = result["transformed"].shape[1]
        assert n_components_actual <= 3

    def test_pca_scale_parameter(self, count_matrix: pd.DataFrame) -> None:
        log_counts = np.log2(count_matrix + 1)
        result_scaled = pca_analysis(log_counts, n_components=2, scale=True)
        result_unscaled = pca_analysis(log_counts, n_components=2, scale=False)
        # Results should differ when scaling is toggled
        assert not np.allclose(
            result_scaled["transformed"].values,
            result_unscaled["transformed"].values,
        )


# =============================================================================
# Tests: compute_sample_distances
# =============================================================================


class TestComputeSampleDistances:
    def test_euclidean_distance_symmetric(self, count_matrix: pd.DataFrame) -> None:
        log_counts = np.log2(count_matrix + 1)
        dist = compute_sample_distances(log_counts, method="euclidean")
        assert dist.shape == (6, 6)
        # Symmetric: dist[i,j] == dist[j,i]
        np.testing.assert_allclose(dist.values, dist.values.T, atol=1e-10)

    def test_euclidean_distance_zero_diagonal(self, count_matrix: pd.DataFrame) -> None:
        log_counts = np.log2(count_matrix + 1)
        dist = compute_sample_distances(log_counts, method="euclidean")
        np.testing.assert_allclose(np.diag(dist.values), 0.0, atol=1e-10)

    def test_euclidean_distance_nonnegative(self, count_matrix: pd.DataFrame) -> None:
        log_counts = np.log2(count_matrix + 1)
        dist = compute_sample_distances(log_counts, method="euclidean")
        assert (dist.values >= -1e-10).all()

    def test_correlation_distance_range(self, count_matrix: pd.DataFrame) -> None:
        log_counts = np.log2(count_matrix + 1)
        dist = compute_sample_distances(log_counts, method="correlation")
        # Correlation distance is in [0, 2] (1 - r, where r in [-1, 1])
        assert (dist.values >= -1e-10).all()
        assert (dist.values <= 2.0 + 1e-10).all()

    def test_cosine_distance_range(self, count_matrix: pd.DataFrame) -> None:
        log_counts = np.log2(count_matrix + 1)
        dist = compute_sample_distances(log_counts, method="cosine")
        # Cosine distance is in [0, 2]
        assert (dist.values >= -1e-10).all()
        assert (dist.values <= 2.0 + 1e-10).all()

    def test_identical_samples_zero_distance(self) -> None:
        counts = pd.DataFrame(
            {"s1": [10, 20, 30], "s2": [10, 20, 30]},
            index=["g1", "g2", "g3"],
        )
        dist = compute_sample_distances(counts, method="euclidean")
        assert abs(dist.loc["s1", "s2"]) < 1e-10

    def test_column_and_index_names(self, count_matrix: pd.DataFrame) -> None:
        log_counts = np.log2(count_matrix + 1)
        dist = compute_sample_distances(log_counts, method="euclidean")
        assert list(dist.index) == list(count_matrix.columns)
        assert list(dist.columns) == list(count_matrix.columns)

    def test_empty_input(self) -> None:
        result = compute_sample_distances(pd.DataFrame())
        assert result.empty

    def test_unknown_method_raises(self, count_matrix: pd.DataFrame) -> None:
        with pytest.raises(ValueError, match="Unknown distance method"):
            compute_sample_distances(count_matrix, method="manhattan")  # type: ignore[arg-type]


# =============================================================================
# Tests: filter_low_expression
# =============================================================================


class TestFilterLowExpression:
    def test_removes_low_expression_genes(self, small_count_matrix: pd.DataFrame) -> None:
        result = filter_low_expression(small_count_matrix, min_count=10, min_samples=2)
        # gene_zero has all zeros => should be removed
        assert "gene_zero" not in result.index

    def test_keeps_high_expression_genes(self, small_count_matrix: pd.DataFrame) -> None:
        result = filter_low_expression(small_count_matrix, min_count=10, min_samples=2)
        assert "gene_A" in result.index
        assert "gene_B" in result.index
        assert "gene_D" in result.index

    def test_returns_fewer_genes(self, small_count_matrix: pd.DataFrame) -> None:
        result = filter_low_expression(small_count_matrix, min_count=10, min_samples=2)
        assert len(result) < len(small_count_matrix)

    def test_preserves_columns(self, small_count_matrix: pd.DataFrame) -> None:
        result = filter_low_expression(small_count_matrix, min_count=10, min_samples=2)
        assert list(result.columns) == list(small_count_matrix.columns)

    def test_min_count_threshold(self) -> None:
        counts = pd.DataFrame(
            {"s1": [5, 15, 100], "s2": [8, 20, 90], "s3": [3, 18, 110]},
            index=["low", "mid", "high"],
        )
        result = filter_low_expression(counts, min_count=10, min_samples=2)
        assert "low" not in result.index  # max is 8, never reaches 10
        assert "mid" in result.index  # 15, 20, 18 all >= 10
        assert "high" in result.index

    def test_min_samples_threshold(self) -> None:
        counts = pd.DataFrame(
            {"s1": [50, 50], "s2": [3, 50], "s3": [2, 50]},
            index=["sparse", "dense"],
        )
        # sparse: only 1 sample >= 10, needs 2
        result = filter_low_expression(counts, min_count=10, min_samples=2)
        assert "sparse" not in result.index
        assert "dense" in result.index

    def test_empty_input(self) -> None:
        result = filter_low_expression(pd.DataFrame())
        assert result.empty

    def test_all_genes_pass(self) -> None:
        counts = pd.DataFrame(
            {"s1": [100, 200], "s2": [150, 300]},
            index=["g1", "g2"],
        )
        result = filter_low_expression(counts, min_count=10, min_samples=1)
        assert len(result) == 2

    def test_all_genes_fail(self) -> None:
        counts = pd.DataFrame(
            {"s1": [1, 2], "s2": [3, 4]},
            index=["g1", "g2"],
        )
        result = filter_low_expression(counts, min_count=100, min_samples=1)
        assert len(result) == 0


# =============================================================================
# Tests: get_highly_variable_genes
# =============================================================================


class TestGetHighlyVariableGenes:
    def test_returns_correct_count(self, count_matrix: pd.DataFrame) -> None:
        log_counts = np.log2(count_matrix + 1)
        hvg = get_highly_variable_genes(log_counts, n_top=10)
        assert len(hvg) == 10

    def test_returns_gene_names(self, count_matrix: pd.DataFrame) -> None:
        log_counts = np.log2(count_matrix + 1)
        hvg = get_highly_variable_genes(log_counts, n_top=5)
        for gene in hvg:
            assert gene in count_matrix.index

    def test_cv_method_selects_variable_genes(self) -> None:
        # gene_variable has high CV, gene_constant has low CV
        counts = pd.DataFrame(
            {
                "s1": [10, 100],
                "s2": [100, 100],
                "s3": [1000, 100],
                "s4": [5, 100],
            },
            index=["gene_variable", "gene_constant"],
        )
        hvg = get_highly_variable_genes(counts, n_top=1, method="cv")
        assert hvg[0] == "gene_variable"

    def test_variance_method_selects_variable_genes(self) -> None:
        counts = pd.DataFrame(
            {
                "s1": [10, 100],
                "s2": [1000, 101],
                "s3": [500, 99],
                "s4": [50, 100],
            },
            index=["gene_variable", "gene_constant"],
        )
        hvg = get_highly_variable_genes(counts, n_top=1, method="variance")
        assert hvg[0] == "gene_variable"

    def test_n_top_exceeds_gene_count(self, count_matrix: pd.DataFrame) -> None:
        log_counts = np.log2(count_matrix + 1)
        hvg = get_highly_variable_genes(log_counts, n_top=9999)
        assert len(hvg) == len(count_matrix)

    def test_empty_input(self) -> None:
        result = get_highly_variable_genes(pd.DataFrame())
        assert result == []

    def test_unknown_method_raises(self, count_matrix: pd.DataFrame) -> None:
        with pytest.raises(ValueError, match="Unknown variance method"):
            get_highly_variable_genes(count_matrix, method="invalid")  # type: ignore[arg-type]


# =============================================================================
# Tests: prepare_volcano_data
# =============================================================================


class TestPrepareVolcanoData:
    def test_adds_regulation_column(self, de_results: pd.DataFrame) -> None:
        result = prepare_volcano_data(de_results, fc_threshold=1.0, pvalue_threshold=0.05)
        assert "regulation" in result.columns
        assert "neg_log10_pvalue" in result.columns

    def test_regulation_categories(self, de_results: pd.DataFrame) -> None:
        result = prepare_volcano_data(de_results, fc_threshold=1.0, pvalue_threshold=0.05)
        valid_labels = {"up", "down", "ns"}
        assert set(result["regulation"].unique()).issubset(valid_labels)

    def test_upregulated_criteria(self) -> None:
        de = pd.DataFrame(
            {
                "gene": ["up_gene"],
                "log2_fold_change": [2.0],
                "p_value": [0.001],
                "adjusted_p_value": [0.01],
                "base_mean": [500.0],
                "stat": [5.0],
            }
        )
        result = prepare_volcano_data(de, fc_threshold=1.0, pvalue_threshold=0.05)
        assert result["regulation"].values[0] == "up"

    def test_downregulated_criteria(self) -> None:
        de = pd.DataFrame(
            {
                "gene": ["down_gene"],
                "log2_fold_change": [-2.0],
                "p_value": [0.001],
                "adjusted_p_value": [0.01],
                "base_mean": [500.0],
                "stat": [-5.0],
            }
        )
        result = prepare_volcano_data(de, fc_threshold=1.0, pvalue_threshold=0.05)
        assert result["regulation"].values[0] == "down"

    def test_not_significant_criteria(self) -> None:
        de = pd.DataFrame(
            {
                "gene": ["ns_gene"],
                "log2_fold_change": [0.5],  # below fc_threshold
                "p_value": [0.001],
                "adjusted_p_value": [0.01],
                "base_mean": [500.0],
                "stat": [1.0],
            }
        )
        result = prepare_volcano_data(de, fc_threshold=1.0, pvalue_threshold=0.05)
        assert result["regulation"].values[0] == "ns"

    def test_ns_due_to_pvalue(self) -> None:
        de = pd.DataFrame(
            {
                "gene": ["ns_pval"],
                "log2_fold_change": [3.0],  # above fc_threshold
                "p_value": [0.5],
                "adjusted_p_value": [0.8],  # above pvalue_threshold
                "base_mean": [500.0],
                "stat": [3.0],
            }
        )
        result = prepare_volcano_data(de, fc_threshold=1.0, pvalue_threshold=0.05)
        assert result["regulation"].values[0] == "ns"

    def test_neg_log10_pvalue_computed(self) -> None:
        de = pd.DataFrame(
            {
                "gene": ["g1"],
                "log2_fold_change": [1.0],
                "p_value": [0.01],
                "adjusted_p_value": [0.05],
                "base_mean": [100.0],
                "stat": [2.0],
            }
        )
        result = prepare_volcano_data(de, use_adjusted=True)
        expected = -np.log10(0.05)
        assert abs(result["neg_log10_pvalue"].values[0] - expected) < 0.01

    def test_use_adjusted_false(self) -> None:
        de = pd.DataFrame(
            {
                "gene": ["g1"],
                "log2_fold_change": [2.0],
                "p_value": [0.01],
                "adjusted_p_value": [0.5],
                "base_mean": [100.0],
                "stat": [3.0],
            }
        )
        # With adjusted: p=0.5 > 0.05 => ns
        result_adj = prepare_volcano_data(de, fc_threshold=1.0, pvalue_threshold=0.05, use_adjusted=True)
        assert result_adj["regulation"].values[0] == "ns"
        # Without adjusted: p=0.01 < 0.05 => up
        result_raw = prepare_volcano_data(de, fc_threshold=1.0, pvalue_threshold=0.05, use_adjusted=False)
        assert result_raw["regulation"].values[0] == "up"

    def test_empty_input(self) -> None:
        result = prepare_volcano_data(pd.DataFrame())
        assert "regulation" in result.columns
        assert "neg_log10_pvalue" in result.columns

    def test_preserves_original_columns(self, de_results: pd.DataFrame) -> None:
        result = prepare_volcano_data(de_results)
        for col in de_results.columns:
            assert col in result.columns


# =============================================================================
# Tests: prepare_ma_data
# =============================================================================


class TestPrepareMAData:
    def test_adds_a_and_m_columns(self, de_results: pd.DataFrame) -> None:
        result = prepare_ma_data(de_results)
        assert "A" in result.columns
        assert "M" in result.columns

    def test_m_equals_log2fc(self, de_results: pd.DataFrame) -> None:
        result = prepare_ma_data(de_results)
        np.testing.assert_array_equal(result["M"].values, result["log2_fold_change"].values)

    def test_a_is_log2_base_mean(self) -> None:
        de = pd.DataFrame(
            {
                "gene": ["g1", "g2"],
                "log2_fold_change": [1.0, -1.0],
                "p_value": [0.01, 0.05],
                "adjusted_p_value": [0.05, 0.1],
                "base_mean": [100.0, 1000.0],
                "stat": [2.0, -2.0],
            }
        )
        result = prepare_ma_data(de)
        expected_a1 = np.log2(100.0 + 1)
        expected_a2 = np.log2(1000.0 + 1)
        assert abs(result["A"].values[0] - expected_a1) < 0.01
        assert abs(result["A"].values[1] - expected_a2) < 0.01

    def test_empty_input(self) -> None:
        result = prepare_ma_data(pd.DataFrame())
        assert "A" in result.columns
        assert "M" in result.columns

    def test_preserves_original_columns(self, de_results: pd.DataFrame) -> None:
        result = prepare_ma_data(de_results)
        for col in de_results.columns:
            assert col in result.columns


# =============================================================================
# Tests: Integration / End-to-End Workflows
# =============================================================================


class TestEndToEndWorkflow:
    def test_full_de_pipeline(self) -> None:
        """Test the complete DE analysis pipeline: normalize -> filter -> DE -> volcano."""
        n_genes = 60
        n_ctrl, n_treat = 4, 4

        # Build a deterministic count matrix with clear differential expression.
        # Control group: stable counts drawn from a tight range.
        # Treatment group: first 10 genes have 10x higher expression, next 10 have 10x lower.
        rng = np.random.default_rng(42)

        ctrl_base = rng.integers(400, 600, size=(n_genes, n_ctrl))
        treat_base = rng.integers(400, 600, size=(n_genes, n_treat))

        # Strong upregulation for genes 0-9 in treatment
        treat_base[:10, :] = rng.integers(4000, 6000, size=(10, n_treat))
        # Strong downregulation for genes 10-19 in treatment
        treat_base[10:20, :] = rng.integers(20, 60, size=(10, n_treat))

        counts_arr = np.hstack([ctrl_base, treat_base])
        genes = [f"gene_{i}" for i in range(n_genes)]
        samples = [f"ctrl_{i}" for i in range(n_ctrl)] + [f"treat_{i}" for i in range(n_treat)]
        counts = pd.DataFrame(counts_arr, index=genes, columns=samples)
        conditions = ["control"] * n_ctrl + ["treatment"] * n_treat

        # Step 1: Filter low expression
        filtered = filter_low_expression(counts, min_count=10, min_samples=2)
        assert len(filtered) <= n_genes

        # Step 2: Normalize
        normalized = normalize_counts(filtered, method="cpm")
        assert normalized.shape == filtered.shape

        # Step 3: Differential expression
        de_results = differential_expression(filtered, conditions, method="ttest")
        assert "gene" in de_results.columns
        assert "log2_fold_change" in de_results.columns
        assert "adjusted_p_value" in de_results.columns

        # Step 4: Volcano data
        volcano = prepare_volcano_data(de_results, fc_threshold=1.0, pvalue_threshold=0.05)
        assert "regulation" in volcano.columns
        regulation_counts = volcano["regulation"].value_counts()
        # With the engineered 10x differences and 4 replicates, we should see significant genes
        assert "up" in regulation_counts.index or "down" in regulation_counts.index

        # Step 5: MA data
        ma = prepare_ma_data(de_results)
        assert "A" in ma.columns
        assert "M" in ma.columns

    def test_normalize_then_pca(self, count_matrix: pd.DataFrame) -> None:
        """Test normalization followed by PCA analysis."""
        normalized = normalize_counts(count_matrix, method="log2cpm")
        pca_result = pca_analysis(normalized, n_components=2)
        assert pca_result["transformed"].shape == (6, 2)
        assert pca_result["explained_variance_ratio"].sum() <= 1.0 + 1e-6

    def test_normalize_then_distance(self, count_matrix: pd.DataFrame) -> None:
        """Test normalization followed by sample distance computation."""
        normalized = normalize_counts(count_matrix, method="log2cpm")
        distances = compute_sample_distances(normalized, method="correlation")
        assert distances.shape == (6, 6)
        np.testing.assert_allclose(np.diag(distances.values), 0.0, atol=1e-10)

    def test_filter_then_hvg(self, count_matrix: pd.DataFrame) -> None:
        """Test filtering low expression then selecting highly variable genes."""
        filtered = filter_low_expression(count_matrix, min_count=10, min_samples=2)
        normalized = normalize_counts(filtered, method="log2cpm")
        hvg = get_highly_variable_genes(normalized, n_top=10, method="cv")
        assert len(hvg) <= 10
        for gene in hvg:
            assert gene in filtered.index

    def test_all_normalization_methods_produce_valid_output(
        self, count_matrix: pd.DataFrame, gene_lengths: pd.Series
    ) -> None:
        """Verify every normalization method runs without error and produces finite output."""
        methods_no_lengths = ["cpm", "log2cpm", "quantile", "median_ratio"]
        methods_with_lengths = ["tpm", "rpkm"]

        for method in methods_no_lengths:
            result = normalize_counts(count_matrix, method=method)  # type: ignore[arg-type]
            assert result.shape == count_matrix.shape, f"{method} changed shape"
            assert np.isfinite(result.values).all(), f"{method} produced non-finite values"

        for method in methods_with_lengths:
            result = normalize_counts(count_matrix, method=method, gene_lengths=gene_lengths)  # type: ignore[arg-type]
            assert result.shape == count_matrix.shape, f"{method} changed shape"
            assert np.isfinite(result.values).all(), f"{method} produced non-finite values"

    def test_all_de_methods_produce_valid_output(
        self, small_count_matrix: pd.DataFrame, conditions_6sample: List[str]
    ) -> None:
        """Verify every DE method runs without error and produces valid results."""
        for method in ["deseq2_like", "ttest", "wilcoxon"]:
            result = differential_expression(small_count_matrix, conditions_6sample, method=method)  # type: ignore[arg-type]
            assert "gene" in result.columns
            assert "log2_fold_change" in result.columns
            assert "p_value" in result.columns
            assert "adjusted_p_value" in result.columns
            assert (result["p_value"] >= 0).all(), f"{method}: negative p-value"
            assert (result["p_value"] <= 1.0 + 1e-10).all(), f"{method}: p-value > 1"

    def test_all_distance_methods_produce_valid_output(self, count_matrix: pd.DataFrame) -> None:
        """Verify every distance method runs without error."""
        log_counts = np.log2(count_matrix + 1)
        for method in ["euclidean", "correlation", "cosine"]:
            result = compute_sample_distances(log_counts, method=method)  # type: ignore[arg-type]
            assert result.shape == (6, 6)
            np.testing.assert_allclose(result.values, result.values.T, atol=1e-10)
