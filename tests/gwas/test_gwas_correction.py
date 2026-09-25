"""Tests for GWAS multiple testing correction."""

from __future__ import annotations

import math

import numpy as np
import pytest

from metainformant.gwas.analysis.correction import (
    EXPECTED_MEDIAN_CHI2_1DF,
    _estimate_pi0,
    adjust_p_values,
    bonferroni_correction,
    fdr_correction,
    genomic_control,
    lambda_gc_from_p_values,
    qvalue_estimation,
)


def test_bonferroni_correction_basic() -> None:
    """Test basic Bonferroni correction."""
    pvalues = [0.05, 0.01, 0.001, 0.0001, 0.5]
    n_tests = len(pvalues)

    result = bonferroni_correction(pvalues, alpha=0.05)

    assert result["status"] == "success"
    assert result["n_tests"] == n_tests
    assert result["corrected_alpha"] == 0.05 / n_tests
    assert result["corrected_alpha"] < 0.05
    assert result["significant_count"] >= 0
    assert result["significant_count"] <= n_tests


def test_bonferroni_correction_no_significant() -> None:
    """Test Bonferroni correction with no significant results."""
    pvalues = [0.1, 0.2, 0.3, 0.4, 0.5]

    result = bonferroni_correction(pvalues, alpha=0.05)

    assert result["status"] == "success"
    assert result["significant_count"] == 0


def test_bonferroni_correction_all_significant() -> None:
    """Test Bonferroni correction where all results are significant."""
    # Very small p-values
    pvalues = [1e-10, 1e-9, 1e-8, 1e-7, 1e-6]

    result = bonferroni_correction(pvalues, alpha=0.05)

    assert result["status"] == "success"
    # With very small p-values, some should still be significant after correction
    assert result["significant_count"] >= 0


def test_fdr_correction_basic() -> None:
    """Test basic FDR (Benjamini-Hochberg) correction."""
    pvalues = [0.05, 0.01, 0.001, 0.0001, 0.5]

    result = fdr_correction(pvalues, alpha=0.05)

    assert result["status"] == "success"
    assert "corrected_pvalues" in result
    assert len(result["corrected_pvalues"]) == len(pvalues)
    assert all(0 <= p <= 1 for p in result["corrected_pvalues"])
    assert result["significant_count"] >= 0


def test_fdr_correction_properties() -> None:
    """Test FDR correction mathematical properties."""
    pvalues = [0.001, 0.01, 0.05, 0.1, 0.5]

    result = fdr_correction(pvalues, alpha=0.05)

    assert result["status"] == "success"
    corrected = np.array(result["corrected_pvalues"])

    # FDR-corrected p-values should be >= original p-values
    assert all(corrected[i] >= pvalues[i] for i in range(len(pvalues)))

    # FDR-corrected p-values should be monotonic (non-decreasing after sorting)
    sorted_indices = np.argsort(pvalues)
    sorted_corrected = corrected[sorted_indices]
    assert all(
        sorted_corrected[i] <= sorted_corrected[i + 1]
        for i in range(len(sorted_corrected) - 1)
    )


def test_genomic_control_from_pvalues() -> None:
    """Test genomic control calculation from p-values."""
    # Create p-values that would give chi-square statistics
    pvalues = [0.001, 0.01, 0.05, 0.1, 0.5, 0.8, 0.9]

    result = genomic_control(pvalues=pvalues)

    assert result["status"] == "success"
    assert "lambda_gc" in result
    assert result["lambda_gc"] > 0
    assert "median_chi2" in result
    assert "n_tests" in result


def test_genomic_control_from_chi2() -> None:
    """Test genomic control calculation from chi-square statistics."""
    chi2_stats = [0.5, 1.0, 2.0, 3.0, 4.0, 5.0]

    result = genomic_control(chi2_stats=chi2_stats)

    assert result["status"] == "success"
    assert "lambda_gc" in result
    assert result["lambda_gc"] > 0


def test_genomic_control_no_data() -> None:
    """Test genomic control with no valid data."""
    result = genomic_control(pvalues=[])
    assert result["status"] == "failed"
    assert "error" in result


def test_genomic_control_inflation() -> None:
    """Test genomic control detects inflation."""
    # P-values that would indicate inflation (many very small p-values)
    # Under null, median chi2 ≈ 0.456
    # With inflation, median chi2 > 0.456, so lambda_GC > 1
    pvalues = [1e-10] * 50 + [0.5] * 50  # Many significant, many not

    result = genomic_control(pvalues=pvalues)

    assert result["status"] == "success"
    # Lambda_GC might be > 1 if there's real signal or inflation
    assert result["lambda_gc"] > 0


def test_lambda_gc_from_p_values_null_is_near_one() -> None:
    """Uniform p-values should give lambda close to 1 under the exact chi2 conversion."""
    np.random.seed(42)
    p_values = list(np.random.uniform(0, 1, 10000))

    lambda_gc = lambda_gc_from_p_values(p_values)

    assert lambda_gc is not None
    assert 0.8 < lambda_gc < 1.2


def test_lambda_gc_from_p_values_detects_inflation() -> None:
    """A mixture of strong signals must push the chi2 median far above the null median."""
    p_values = [1e-10] * 50 + [0.5] * 50

    lambda_gc = lambda_gc_from_p_values(p_values)

    assert lambda_gc is not None
    assert lambda_gc > 5


def test_lambda_gc_from_p_values_matches_genomic_control() -> None:
    """The shared helper and genomic_control must agree on the same p-value set."""
    p_values = [0.001, 0.01, 0.1, 0.5, 0.9]

    result = genomic_control(p_values=p_values)

    assert isinstance(result, dict)
    assert lambda_gc_from_p_values(p_values) == pytest.approx(result["lambda_gc"])


def test_lambda_gc_from_p_values_invalid_input() -> None:
    """Empty or all-invalid p-values return None; invalid entries are skipped."""
    assert lambda_gc_from_p_values([]) is None
    assert lambda_gc_from_p_values([0.0, -0.5, 1.5, float("nan")]) is None
    assert lambda_gc_from_p_values([0.5, "not-a-p-value", None]) == pytest.approx(
        1.0, abs=0.01
    )


def test_genomic_control_uses_one_df_chi_square_transform() -> None:
    """P-value based lambda GC should match the 1-df chi-square survival scale."""
    scipy_stats = pytest.importorskip("scipy.stats")
    chi2_stats = np.array([0.05, 0.20, 0.60, 1.40, 3.20, 8.00])
    pvalues = scipy_stats.chi2.sf(chi2_stats, 1).tolist()

    result = genomic_control(pvalues=pvalues)

    expected_median = float(np.median(chi2_stats))
    assert result["status"] == "success"
    assert result["median_chi2"] == pytest.approx(expected_median, rel=1e-10)
    assert result["lambda_gc"] == pytest.approx(
        expected_median / EXPECTED_MEDIAN_CHI2_1DF, rel=1e-10
    )
    assert all(0 <= p <= 1 for p in result["corrected_p_values"])


@pytest.mark.parametrize("method", ["bonferroni", "fdr", "genomic_control", "qvalue"])
def test_adjust_p_values_wrapper_returns_adjusted_values(method: str) -> None:
    """Generic p-value adjustment should use legacy tuple paths internally."""
    adjusted = adjust_p_values([0.001, 0.01, 0.2, 0.8], method=method)

    assert len(adjusted) == 4
    assert all(0 <= p <= 1 for p in adjusted)


def test_estimate_pi0_flat_null_is_one() -> None:
    """Uniform p-values (flat null histogram) must estimate pi0 near 1."""
    rng = np.random.default_rng(11)
    p_values = rng.uniform(0.0, 1.0, 4000).tolist()

    pi0 = _estimate_pi0(p_values)

    assert 0.0 < pi0 <= 1.0
    assert pi0 == pytest.approx(1.0, abs=0.05)


def test_estimate_pi0_mixture_below_one() -> None:
    """A block of strong signals must push the estimated pi0 below 1."""
    rng = np.random.default_rng(14)
    p_values = [1e-8] * 200 + rng.uniform(0.0, 1.0, 800).tolist()

    pi0 = _estimate_pi0(p_values)

    assert 0.0 < pi0 < 0.95
    assert pi0 == pytest.approx(1.0 - 200 / 1000, abs=0.15)


def test_estimate_pi0_bounds_and_empty() -> None:
    """pi0 must stay in (0, 1] for degenerate inputs."""
    assert _estimate_pi0([]) == 1.0
    assert _estimate_pi0([1e-300] * 50) > 0.0
    assert _estimate_pi0([1e-300] * 50) <= 1.0
    assert _estimate_pi0(["bad", None, float("nan")]) == 1.0


def test_qvalue_with_pi0_one_equals_bh() -> None:
    """With pi0 == 1 the q-value procedure is exactly Benjamini-Hochberg."""
    p_values = [0.001, 0.008, 0.039, 0.11, 0.34, 0.62, 0.9]

    q_values, pi0 = qvalue_estimation(p_values, pi0=1.0)
    _, adjusted = fdr_correction(p_values, return_dict=False)

    assert pi0 == 1.0
    assert q_values == adjusted


def test_qvalue_estimated_pi0_never_exceeds_bh() -> None:
    """Estimated pi0 <= 1 keeps q-values at or below the BH adjustment."""
    rng = np.random.default_rng(13)
    p_values = [1e-10] * 40 + rng.uniform(0.0, 1.0, 460).tolist()

    q_values, pi0 = qvalue_estimation(p_values)
    _, adjusted = fdr_correction(p_values, return_dict=False)

    assert 0.0 < pi0 <= 1.0
    assert pi0 < 0.95
    assert all(q <= b + 1e-12 for q, b in zip(q_values, adjusted))


def test_genomic_control_chi2_returns_corrected_p_values() -> None:
    """Chi-square input must return lambda-corrected 1-df p-values aligned to input."""
    scipy_stats = pytest.importorskip("scipy.stats")
    chi2_stats = [0.5, 1.0, 2.0, 3.0, 4.0, 5.0]

    result = genomic_control(chi2_stats=chi2_stats)

    assert result["status"] == "success"
    corrected = result["corrected_p_values"]
    assert len(corrected) == len(chi2_stats)
    assert result["n_tests"] == len(chi2_stats)
    lam = result["lambda_gc"]
    for stat, corrected_p in zip(chi2_stats, corrected):
        assert corrected_p == pytest.approx(
            float(scipy_stats.chi2.sf(stat / lam, 1)), rel=1e-9
        )
    assert all(0.0 < p <= 1.0 for p in corrected)

    # Uninterpretable entries stay aligned as NaN instead of shrinking output.
    mixed = genomic_control(chi2_stats=[1.0, -1.0, 2.0])
    assert len(mixed["corrected_p_values"]) == 3
    assert math.isnan(mixed["corrected_p_values"][1])
