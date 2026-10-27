"""Zero residual variance: a constant phenotype is a failed test; an exact fit has p = 0, never p = 1.

Zero-mocks: the real association functions run on small hand-built data.
"""

from __future__ import annotations

import math

import pytest

from metainformant.gwas.analysis.association import association_test_linear


def _assert_failed(result: dict) -> None:
    assert result["status"] == "error"
    assert "constant" in result["error"]
    for key in ("beta", "se", "t_stat", "p_value"):
        assert math.isnan(result[key]), key


def test_constant_phenotype_is_a_failed_test() -> None:
    result = association_test_linear([0, 1, 2, 1, 0, 2], [3.0] * 6)
    _assert_failed(result)


def test_exact_fit_has_an_infinite_statistic_and_p_zero() -> None:
    genotypes = [0, 1, 2, 1, 0, 2]
    result = association_test_linear(genotypes, [1.0 + 0.5 * g for g in genotypes])
    assert result["status"] == "success"
    assert result["beta"] == pytest.approx(0.5)
    assert result["t_stat"] == math.inf
    assert result["p_value"] == 0.0


def test_constant_phenotype_with_a_covariate_is_a_failed_test() -> None:
    genotypes = [0, 1, 2, 1, 0, 2, 1, 0]
    covariate = [0.3, -1.2, 0.8, 2.0, -0.5, 1.1, 0.0, -0.9]  # covariate-major: one list per covariate
    result = association_test_linear(genotypes, [4.0] * 8, covariates=[covariate])
    _assert_failed(result)


def test_an_ordinary_fit_is_unchanged() -> None:
    genotypes = [0, 1, 2, 1, 0, 2, 1, 0]
    phenotypes = [1.1, 1.9, 3.2, 2.1, 0.8, 2.7, 1.6, 1.3]
    result = association_test_linear(genotypes, phenotypes)
    assert result["status"] == "success"
    assert 0.0 < result["p_value"] < 0.01
    scipy_stats = pytest.importorskip("scipy.stats")
    reference = scipy_stats.linregress(genotypes, phenotypes)
    assert result["beta"] == pytest.approx(reference.slope, rel=1e-9)
    assert result["p_value"] == pytest.approx(reference.pvalue, rel=1e-6)
