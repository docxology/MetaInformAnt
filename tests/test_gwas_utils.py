"""Tests for GWAS shared utility functions.

Tests cover:
- R-squared computation (analysis/utils.py)
- P-value key detection (visualization/utils.py)
- Safe log10 p-value transform (visualization/utils.py)
"""

from __future__ import annotations

import math

import numpy as np
import pytest


class TestComputeRSquared:
    """Tests for the shared compute_r_squared function in analysis/utils.py."""

    def test_perfectly_correlated_genotypes(self) -> None:
        """Identical genotype vectors should yield r2 = 1.0."""
        from metainformant.gwas.analysis.utils import compute_r_squared

        geno = [0, 1, 2, 1, 0, 2, 1, 0, 2, 1]
        result = compute_r_squared(geno, list(geno))
        assert abs(result - 1.0) < 1e-6

    def test_uncorrelated_genotypes(self) -> None:
        """Genotypes constructed to be uncorrelated should yield r2 near 0."""
        from metainformant.gwas.analysis.utils import compute_r_squared

        # Construct two vectors with known zero correlation
        rng = np.random.RandomState(42)
        n = 500
        a = rng.choice([0, 1, 2], size=n).tolist()
        b = rng.choice([0, 1, 2], size=n).tolist()
        result = compute_r_squared(a, b)
        # With independent draws, r2 should be small (< 0.05)
        assert 0.0 <= result < 0.1

    def test_negatively_correlated_genotypes(self) -> None:
        """Anti-correlated vectors: r2 should still be positive (it's squared)."""
        from metainformant.gwas.analysis.utils import compute_r_squared

        a = [0, 0, 0, 2, 2, 2, 1, 1, 1, 0]
        b = [2, 2, 2, 0, 0, 0, 1, 1, 1, 2]
        result = compute_r_squared(a, b)
        assert result > 0.5  # Strong negative correlation -> high r2

    def test_missing_values_excluded(self) -> None:
        """Missing values (< 0) should be excluded from computation."""
        from metainformant.gwas.analysis.utils import compute_r_squared

        a = [0, 1, 2, -1, 1, 0, 2, 1, -1, 0]
        b = [0, 1, 2, 0, 1, 0, 2, 1, 1, 0]
        result = compute_r_squared(a, b)
        # The non-missing pairs are perfectly correlated
        assert result > 0.9

    def test_too_few_valid_pairs_returns_zero(self) -> None:
        """Fewer than 3 valid pairs should return 0.0."""
        from metainformant.gwas.analysis.utils import compute_r_squared

        # Only 2 valid pairs
        a = [1, 2, -1, -1, -1]
        b = [1, 2, -1, -1, -1]
        result = compute_r_squared(a, b)
        assert result == 0.0

    def test_constant_genotype_returns_zero(self) -> None:
        """Zero-variance (constant) genotype vectors should return 0.0."""
        from metainformant.gwas.analysis.utils import compute_r_squared

        a = [1, 1, 1, 1, 1, 1, 1, 1]
        b = [0, 1, 2, 0, 1, 2, 0, 1]
        result = compute_r_squared(a, b)
        assert result == 0.0

    def test_both_constant_returns_zero(self) -> None:
        """Both vectors constant should return 0.0."""
        from metainformant.gwas.analysis.utils import compute_r_squared

        a = [2, 2, 2, 2, 2]
        b = [1, 1, 1, 1, 1]
        result = compute_r_squared(a, b)
        assert result == 0.0

    def test_known_r_squared_value(self) -> None:
        """Verify against a hand-computed r-squared."""
        from metainformant.gwas.analysis.utils import compute_r_squared

        # Simple known values: a = [0,1,2,3], b = [0,2,4,6] -> r = 1.0
        a = [0, 1, 2, 1]
        b = [0, 2, 0, 2]
        result = compute_r_squared(a, b)
        # Compute expected r2 manually via numpy
        expected_r = float(np.corrcoef([0.0, 1.0, 2.0, 1.0], [0.0, 2.0, 0.0, 2.0])[0, 1])
        expected_r2 = expected_r * expected_r
        assert abs(result - expected_r2) < 1e-6

    def test_result_bounded_zero_to_one(self) -> None:
        """R-squared must always be in [0, 1]."""
        from metainformant.gwas.analysis.utils import compute_r_squared

        rng = np.random.RandomState(99)
        for _ in range(20):
            a = rng.choice([0, 1, 2], size=50).tolist()
            b = rng.choice([0, 1, 2], size=50).tolist()
            result = compute_r_squared(a, b)
            assert 0.0 <= result <= 1.0


class TestDetectPValueKey:
    """Tests for p-value key detection from result dicts."""

    def test_standard_p_value_key(self) -> None:
        """Detect the standard 'p_value' key."""
        from metainformant.gwas.visualization.utils import detect_p_value_key

        result = {"chrom": "1", "pos": 100, "p_value": 0.01}
        assert detect_p_value_key(result) == "p_value"

    def test_pvalue_no_underscore(self) -> None:
        """Detect 'pvalue' (no underscore)."""
        from metainformant.gwas.visualization.utils import detect_p_value_key

        result = {"variant": "rs123", "pvalue": 0.05}
        assert detect_p_value_key(result) == "pvalue"

    def test_pval_key(self) -> None:
        """Detect 'pval' key."""
        from metainformant.gwas.visualization.utils import detect_p_value_key

        result = {"pval": 0.001}
        assert detect_p_value_key(result) == "pval"

    def test_uppercase_P_key(self) -> None:
        """Detect uppercase 'P' key (common in PLINK output)."""
        from metainformant.gwas.visualization.utils import detect_p_value_key

        result = {"CHR": "1", "BP": 1000, "P": 1e-8}
        assert detect_p_value_key(result) == "P"

    def test_p_lowercase_key(self) -> None:
        """Detect lowercase 'p' key."""
        from metainformant.gwas.visualization.utils import detect_p_value_key

        result = {"snp": "rs456", "p": 0.5}
        assert detect_p_value_key(result) == "p"

    def test_no_p_value_key_returns_none(self) -> None:
        """Return None when no p-value key is found."""
        from metainformant.gwas.visualization.utils import detect_p_value_key

        result = {"chrom": "1", "pos": 100, "beta": 0.5}
        assert detect_p_value_key(result) is None

    def test_empty_dict_returns_none(self) -> None:
        """Empty dict should return None."""
        from metainformant.gwas.visualization.utils import detect_p_value_key

        assert detect_p_value_key({}) is None

    def test_priority_order_p_value_first(self) -> None:
        """When multiple keys exist, 'p_value' takes priority."""
        from metainformant.gwas.visualization.utils import detect_p_value_key

        result = {"p_value": 0.01, "P": 0.02, "pval": 0.03}
        assert detect_p_value_key(result) == "p_value"

    def test_priority_order_P_over_p(self) -> None:
        """'P' should have higher priority than 'p'."""
        from metainformant.gwas.visualization.utils import detect_p_value_key

        result = {"P": 0.001, "p": 0.002}
        assert detect_p_value_key(result) == "P"


class TestSafeLog10P:
    """Tests for safe_log10_p computation."""

    def test_normal_p_value(self) -> None:
        """Normal p-value should return -log10(p)."""
        from metainformant.gwas.visualization.utils import safe_log10_p

        result = safe_log10_p(0.001)
        assert abs(result - 3.0) < 1e-6

    def test_p_value_one(self) -> None:
        """p = 1.0 should return 0.0."""
        from metainformant.gwas.visualization.utils import safe_log10_p

        result = safe_log10_p(1.0)
        assert abs(result - 0.0) < 1e-6

    def test_p_value_zero_clamped(self) -> None:
        """p = 0 should be clamped (not produce -inf)."""
        from metainformant.gwas.visualization.utils import safe_log10_p

        result = safe_log10_p(0.0)
        assert math.isfinite(result)
        assert result > 0  # -log10(tiny) is positive
        assert result == 300.0  # -log10(1e-300) = 300

    def test_very_small_p_value(self) -> None:
        """Very small p-value should be handled without overflow."""
        from metainformant.gwas.visualization.utils import safe_log10_p

        result = safe_log10_p(1e-300)
        assert math.isfinite(result)
        assert abs(result - 300.0) < 1e-6

    def test_moderate_p_value(self) -> None:
        """Standard moderate p-value check."""
        from metainformant.gwas.visualization.utils import safe_log10_p

        result = safe_log10_p(0.05)
        expected = -math.log10(0.05)
        assert abs(result - expected) < 1e-6

    def test_genome_wide_significance(self) -> None:
        """p = 5e-8 should yield approximately 7.3."""
        from metainformant.gwas.visualization.utils import safe_log10_p

        result = safe_log10_p(5e-8)
        expected = -math.log10(5e-8)
        assert abs(result - expected) < 1e-6

    def test_negative_p_value_clamped(self) -> None:
        """Negative p-value should be clamped to the minimum."""
        from metainformant.gwas.visualization.utils import safe_log10_p

        result = safe_log10_p(-0.01)
        assert math.isfinite(result)
        assert result == 300.0  # Clamped to floor

    def test_return_type_is_float(self) -> None:
        """Return type should always be float."""
        from metainformant.gwas.visualization.utils import safe_log10_p

        assert isinstance(safe_log10_p(0.5), float)
        assert isinstance(safe_log10_p(0.0), float)
        assert isinstance(safe_log10_p(1e-200), float)


class TestImportsFromModules:
    """Verify that the refactored modules import correctly from utils."""

    def test_analysis_utils_importable(self) -> None:
        """analysis.utils module should be importable."""
        from metainformant.gwas.analysis.utils import compute_r_squared

        assert callable(compute_r_squared)

    def test_visualization_utils_importable(self) -> None:
        """visualization.utils module should be importable."""
        from metainformant.gwas.visualization.utils import detect_p_value_key, safe_log10_p

        assert callable(detect_p_value_key)
        assert callable(safe_log10_p)

    def test_analysis_init_exports(self) -> None:
        """analysis __init__ should export compute_r_squared."""
        from metainformant.gwas.analysis import compute_r_squared

        assert callable(compute_r_squared)

    def test_visualization_init_exports(self) -> None:
        """visualization __init__ should export detect_p_value_key and safe_log10_p."""
        from metainformant.gwas.visualization import detect_p_value_key, safe_log10_p

        assert callable(detect_p_value_key)
        assert callable(safe_log10_p)
