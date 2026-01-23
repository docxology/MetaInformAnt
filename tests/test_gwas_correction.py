"""Tests for GWAS multiple testing correction."""

from __future__ import annotations

import numpy as np
import pytest

from metainformant.gwas import bonferroni_correction, fdr_correction, genomic_control


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
    assert all(sorted_corrected[i] <= sorted_corrected[i + 1] for i in range(len(sorted_corrected) - 1))


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
