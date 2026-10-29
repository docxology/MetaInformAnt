"""Regression tests: ``sv_association_test`` reports a test it cannot carry out as missing, and its p-values are exact.

A monomorphic SV (or covariates collinear with the genotype) made the design matrix singular. A ridge term of 1e-8
kept the solver from failing, so the test "ran": ``beta`` 0, ``odds_ratio`` 1 and ``p_value`` 1.0 for a monomorphic
SV, a finite p near 1 for collinear covariates, ``p_value`` 0.0 for a fit with no residual degrees of freedom, p 1.0
for a constant phenotype, and an unconverged logistic fit under separation (``beta`` 18, ``se`` 1,900, p 0.99) as if
it were an estimate. Those are now NaN with the reason in ``error``; the rank of the design matrix decides, never
whether the solver raised, because that depends on the BLAS build.

The p-values were also approximate: the t statistic of a linear fit went through the normal distribution (too small
at small n) and ``2 * (1 - cdf)`` underflows to exactly 0.0 beyond |z| = 8.3.

The oracles are independent: ``scipy.stats.linregress``, a least-squares solve with an explicit t-distribution tail,
and a general-purpose optimiser on the logistic log-likelihood.
"""

from __future__ import annotations

import json
import math
import os
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

import metainformant.structural_variants as structural_variants
from metainformant.structural_variants.population import sv_association_test

scipy_stats = pytest.importorskip("scipy.stats")
scipy_optimize = pytest.importorskip("scipy.optimize")


def _assert_untestable(result: dict[str, object], reason: str, binary: bool = False) -> None:
    for key in ("beta", "se", "p_value"):
        assert math.isnan(result[key])  # type: ignore[arg-type]
    if binary:
        assert isinstance(result["odds_ratio"], float) and math.isnan(result["odds_ratio"])
    else:
        assert result["odds_ratio"] is None
    assert isinstance(result["error"], str) and reason in result["error"]


class TestLinearRegressionAgainstTheOracles:
    def test_a_worked_example_by_hand(self) -> None:
        # g = 0,1,2,0,1,2 and y = 1,2,3,2,3,4: slope Sxy / Sxx = 4 / 4 = 1, intercept 1.5, residuals -0.5 and 0.5
        # (rss 1.5, 4 degrees of freedom, sigma^2 0.375), se = sqrt(0.375 / 4) = 0.30618621784789724.
        result = sv_association_test([0, 1, 2, 0, 1, 2], [1, 2, 3, 2, 3, 4])
        assert result["beta"] == pytest.approx(1.0, abs=1e-12)
        assert result["se"] == pytest.approx(math.sqrt(0.375 / 4), abs=1e-12)
        assert result["p_value"] == pytest.approx(
            scipy_stats.linregress([0, 1, 2, 0, 1, 2], [1, 2, 3, 2, 3, 4]).pvalue, rel=1e-9
        )
        assert result["error"] is None and result["odds_ratio"] is None and result["n_samples"] == 6

    def test_the_p_value_is_the_t_distribution_tail(self) -> None:
        # n = 12, t = 2.38: the t(10) tail is 0.0386, the normal tail (the old answer) 0.0173.
        rng = np.random.default_rng(1)
        genotypes = [0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2]
        phenotypes = [0.4 * g + rng.normal(0, 1) for g in genotypes]
        expected = scipy_stats.linregress(genotypes, phenotypes)
        assert abs(expected.slope / expected.stderr) == pytest.approx(2.38, abs=0.01)

        result = sv_association_test(genotypes, phenotypes)
        assert result["beta"] == pytest.approx(expected.slope, rel=1e-9)
        assert result["se"] == pytest.approx(expected.stderr, rel=1e-9)
        assert result["p_value"] == pytest.approx(expected.pvalue, rel=1e-9)
        assert result["p_value"] == pytest.approx(0.0386, abs=5e-4)
        assert result["p_value_method"] == "t"

    def test_covariates_match_a_least_squares_solve(self) -> None:
        rng = np.random.default_rng(2)
        n = 40
        genotypes = rng.integers(0, 3, n)
        covariates = rng.normal(size=(n, 2))
        phenotypes = 0.8 * genotypes + covariates @ np.array([0.5, -0.3]) + rng.normal(0, 1, n)

        design = np.column_stack([np.ones(n), genotypes, covariates])
        beta, *_ = np.linalg.lstsq(design, phenotypes, rcond=None)
        residuals = phenotypes - design @ beta
        df = n - design.shape[1]
        sigma2 = float(residuals @ residuals) / df
        se = np.sqrt(np.diag(np.linalg.inv(design.T @ design)) * sigma2)
        p = 2 * scipy_stats.t.sf(abs(beta[1] / se[1]), df)

        result = sv_association_test(genotypes.tolist(), phenotypes.tolist(), covariates=covariates.tolist())
        assert (result["beta"], result["se"], result["p_value"]) == pytest.approx((beta[1], se[1], p), rel=1e-9)

    def test_nearly_collinear_covariates_are_solved_accurately(self) -> None:
        # The covariate is 3 x genotype plus noise of 1e-6: full rank (the rank test passes, so the test runs) but
        # with a condition number near 1e7, which the normal equations square. The oracle is the SVD.
        rng = np.random.default_rng(4)
        n = 60
        genotypes = rng.integers(0, 3, n)
        covariate = 3.0 * genotypes + 1e-6 * rng.normal(size=n)
        phenotypes = 0.5 * genotypes + rng.normal(0, 1, n)

        design = np.column_stack([np.ones(n), genotypes, covariate])
        assert np.linalg.matrix_rank(design) == 3 and np.linalg.cond(design) > 1e6  # the premise
        u, singular, vt = np.linalg.svd(design, full_matrices=False)
        beta = vt.T @ ((u.T @ phenotypes) / singular)
        residuals = phenotypes - design @ beta
        sigma2 = float(residuals @ residuals) / (n - 3)
        se = np.sqrt(np.diag(vt.T @ np.diag(singular**-2.0) @ vt) * sigma2)

        result = sv_association_test(genotypes.tolist(), phenotypes.tolist(), covariates=covariate[:, None].tolist())
        # Both are huge (the genotype is nearly a multiple of the covariate); the old 1e-8 ridge returned 212 and 8,980.
        assert result["beta"] == pytest.approx(beta[1], rel=1e-6)
        assert result["se"] == pytest.approx(se[1], rel=1e-6)
        assert result["error"] is None

    def test_a_very_small_p_value_is_not_rounded_to_zero(self) -> None:
        rng = np.random.default_rng(3)
        genotypes = rng.integers(0, 3, 300)
        phenotypes = 2.0 * genotypes + rng.normal(0, 0.5, 300)
        expected = scipy_stats.linregress(genotypes, phenotypes).pvalue
        assert 0.0 < expected < 1e-100  # the premise: far beyond what 1 - cdf can represent

        result = sv_association_test(genotypes.tolist(), phenotypes.tolist())
        assert result["p_value"] > 0.0
        assert math.log10(result["p_value"]) == pytest.approx(math.log10(expected), abs=1e-6)

    def test_numpy_integer_input_is_accepted(self) -> None:
        genotypes = np.array([0, 1, 2, 0, 1, 2])
        phenotypes = np.array([1.0, 2.0, 3.0, 2.0, 3.0, 4.0])
        assert sv_association_test(genotypes, phenotypes)["beta"] == pytest.approx(1.0)  # type: ignore[arg-type]


class TestLogisticRegressionAgainstAnOptimiser:
    GENOTYPES = [0] * 10 + [1] * 10 + [2] * 10
    PHENOTYPES = [0.0] * 9 + [1.0] + [0.0] * 8 + [1.0] * 2 + [0.0] + [1.0] * 9

    def test_the_coefficient_and_wald_p_value_match(self) -> None:
        genotypes, phenotypes = np.array(self.GENOTYPES, float), np.array(self.PHENOTYPES)
        design = np.column_stack([np.ones(30), genotypes])

        def negative_log_likelihood(beta: np.ndarray) -> float:
            eta = design @ beta
            return float(np.sum(np.log1p(np.exp(eta)) - phenotypes * eta))

        fit = scipy_optimize.minimize(negative_log_likelihood, np.zeros(2), method="BFGS", options={"gtol": 1e-10})
        beta = fit.x
        mu = 1 / (1 + np.exp(-(design @ beta)))
        covariance = np.linalg.inv(design.T @ (design * (mu * (1 - mu))[:, None]))
        se = math.sqrt(covariance[1, 1])
        p = math.erfc(abs(beta[1] / se) / math.sqrt(2))

        result = sv_association_test(self.GENOTYPES, self.PHENOTYPES)
        assert result["beta"] == pytest.approx(beta[1], abs=1e-4)
        assert result["se"] == pytest.approx(se, abs=1e-4)
        assert result["p_value"] == pytest.approx(p, rel=1e-3)
        assert result["odds_ratio"] == pytest.approx(math.exp(beta[1]), rel=1e-3)
        assert result["p_value_method"] == "wald_normal" and result["error"] is None


SCIPY_FREE_SCRIPT = """
import json, sys
sys.modules["scipy"] = None  # import scipy raises ImportError, as in an install without it
sys.modules["scipy.stats"] = None
from metainformant.structural_variants.population import sv_association_test
r = sv_association_test([0, 1, 2, 0, 1, 2], [1, 2, 3, 2, 3, 4])
print(json.dumps({key: r[key] for key in ("beta", "se", "p_value", "p_value_method")}))
"""


class TestWithoutScipy:
    def test_the_normal_tail_stands_in_for_the_t_tail_and_the_result_says_so(self, tmp_path: Path) -> None:
        # The child must import the very copy of the package that this process is testing.
        package_root = Path(structural_variants.__file__).parents[2]
        env = dict(os.environ, PYTHONPATH=str(package_root))
        completed = subprocess.run(
            [sys.executable, "-c", SCIPY_FREE_SCRIPT],
            capture_output=True,
            text=True,
            cwd=tmp_path,
            env=env,
            timeout=300,
            check=False,
        )
        assert completed.returncode == 0, completed.stderr[-3000:]
        result = json.loads(completed.stdout.strip().splitlines()[-1])

        # By hand: slope 1, se sqrt(0.375 / 4), so z = 1 / 0.30618... = 3.266; the normal tail is erfc(z / sqrt 2).
        se = math.sqrt(0.375 / 4)
        assert result["beta"] == pytest.approx(1.0) and result["se"] == pytest.approx(se)
        assert result["p_value"] == pytest.approx(math.erfc((1.0 / se) / math.sqrt(2.0)), rel=1e-12)
        assert result["p_value_method"] == "normal_approximation_of_t"
        # It is smaller than the t(4) tail (0.0309) that SciPy gives, which is why the method is reported.
        assert result["p_value"] < scipy_stats.linregress([0, 1, 2, 0, 1, 2], [1, 2, 3, 2, 3, 4]).pvalue


class TestUntestableAssociations:
    def test_a_monomorphic_sv_is_missing_not_an_effect_of_zero(self) -> None:
        genotypes = [0] * 20
        # The old result: beta 0, se 60,781, p 1.0.
        _assert_untestable(sv_association_test(genotypes, [float(i) for i in range(20)]), "singular")
        _assert_untestable(sv_association_test(genotypes, [float(i % 2) for i in range(20)]), "singular", binary=True)

    def test_the_design_matrix_is_rank_deficient_by_the_rank_not_by_a_raised_error(self) -> None:
        genotypes = [0, 1, 2] * 8
        covariates = [[2.0 * g] for g in genotypes]  # exactly twice the genotype
        design = np.column_stack([np.ones(24), genotypes, covariates])
        assert np.linalg.matrix_rank(design) == 2 < design.shape[1]  # the premise, from the rank

        result = sv_association_test(genotypes, [float(i) for i in range(24)], covariates=covariates)
        _assert_untestable(result, "singular")  # the old result: beta 0.2, p 0.99999 from a 1e-8 ridge

    def test_a_constant_covariate_is_collinear_with_the_intercept(self) -> None:
        genotypes = [0, 1, 2] * 4
        result = sv_association_test(genotypes, [float(i) for i in range(12)], covariates=[[3.0]] * 12)
        _assert_untestable(result, "singular")

    def test_no_residual_degrees_of_freedom(self) -> None:
        # Two samples, two parameters: the old p-value was exactly 0.0.
        _assert_untestable(sv_association_test([0, 1], [1.0, 2.0]), "no residual degrees of freedom")
        # Four samples, intercept + genotype + two covariates.
        result = sv_association_test(
            [0, 1, 2, 1], [1.0, 3.0, 2.0, 5.0], covariates=[[0.5, 1.0], [2.0, 0.1], [1.5, 3.0], [0.2, 0.9]]
        )
        _assert_untestable(result, "no residual degrees of freedom")

    def test_a_constant_phenotype(self) -> None:
        # The old result: beta 0, p 1.0.
        _assert_untestable(sv_association_test([0, 1, 2, 0, 1, 2], [5.0] * 6), "constant")

    def test_an_exact_fit_has_no_residual_variance(self) -> None:
        # y = 1 + 2 g exactly, with 2 degrees of freedom: the old p-value was 0.0 from a clipped standard error.
        _assert_untestable(sv_association_test([0, 1, 2, 3], [1.0, 3.0, 5.0, 7.0]), "exact fit")

    @pytest.mark.parametrize(
        ("genotypes", "phenotypes"),
        [
            pytest.param([0] * 10 + [2] * 10, [0.0] * 10 + [1.0] * 10, id="complete-separation"),
            pytest.param([0] * 10 + [1] * 10 + [2] * 10, [0.0] * 15 + [1.0] * 15, id="quasi-complete-separation"),
            pytest.param([0, 1, 2] * 8, [1.0] * 24, id="all-cases"),
            pytest.param([0, 1, 2] * 8, [0.0] * 24, id="all-controls"),
        ],
    )
    def test_a_logistic_fit_that_has_no_estimate_is_missing(
        self, genotypes: list[int], phenotypes: list[float]
    ) -> None:
        # Under separation the old iteration stopped after 50 steps at beta ~ 18-21 with se in the thousands and a
        # p-value near 1, and reported it as an estimate.
        result = sv_association_test(genotypes, phenotypes)
        assert math.isnan(result["beta"]) and math.isnan(result["p_value"]) and math.isnan(result["odds_ratio"])
        assert isinstance(result["error"], str) and result["error"]

    def test_the_failure_names_its_cause(self) -> None:
        result = sv_association_test([0] * 10 + [2] * 10, [0.0] * 10 + [1.0] * 10)
        assert "separation" in result["error"]  # type: ignore[operator]


class TestInputs:
    def test_a_length_mismatch_still_raises(self) -> None:
        with pytest.raises(ValueError, match="Length mismatch"):
            sv_association_test([0, 1, 2], [1.0, 2.0])

    @pytest.mark.parametrize("bad", [float("nan"), float("inf")])
    def test_non_finite_values_are_refused(self, bad: float) -> None:
        with pytest.raises(ValueError, match="finite"):
            sv_association_test([0, 1, 2, 1], [1.0, 2.0, bad, 4.0])
        with pytest.raises(ValueError, match="finite"):
            sv_association_test([0, 1, bad, 1], [1.0, 2.0, 3.0, 4.0])
        with pytest.raises(ValueError, match="finite"):
            sv_association_test([0, 1, 2, 1, 0, 2], [1.0, 2.0, 3.0, 4.0, 2.0, 1.0], covariates=[[bad]] * 6)
