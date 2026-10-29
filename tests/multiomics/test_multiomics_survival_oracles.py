"""``cox_regression`` and ``log_rank_test`` agree with statsmodels, and refuse what they cannot fit.

Regression tests for ``metainformant.multiomics.survival.analysis``:

* ``cox_regression`` reported a log partial likelihood off by a constant per event (the stabilising shift of the risk
  scores was never added back: -112 instead of -578 on the data below), standard errors 3-4% too small (a Hessian
  whose off-diagonal terms were an approximation, plus ridge terms), Wald p-values that were wrong with them, and
  coefficients that were not the Breslow maximum when event times tie (the risk set shrank inside a tie). A singular
  Hessian was patched with a gradient step of ``0.01 * gradient`` and the iteration could stop unconverged without
  a word.
* ``log_rank_test`` summed ``(O - E)**2 / V`` over groups, which is the statistic for two groups only; with three or
  more groups the groups' scores are correlated and the chi-squared needs the full covariance matrix. A group that
  is never at risk when an event occurs kept a degree of freedom, and a comparison with nothing to compare returned
  chi-squared 0, p = 1.0.

The oracles are statsmodels' ``PHReg`` (``ties="breslow"``) and ``survdiff``, an O(n^2) log-likelihood written out
in the test, and the closed form of the score. Real computation; no test doubles (real-implementation policy).
"""

from __future__ import annotations

import math
from typing import Any

import numpy as np
import pytest

from metainformant.multiomics.survival.analysis import cox_regression, log_rank_test

hazard_regression = pytest.importorskip("statsmodels.duration.hazard_regression")
survfunc = pytest.importorskip("statsmodels.duration.survfunc")


def _cox_data(seed: int = 11, n: int = 200, ties: bool = False) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    rng = np.random.default_rng(seed)
    covariates = rng.standard_normal((n, 3))
    latent = rng.exponential(1.0 / np.exp(covariates @ np.array([0.8, -0.5, 0.3])))
    censoring = rng.exponential(2.0, n)
    time = np.minimum(latent, censoring)
    event = (latent <= censoring).astype(int)
    return (np.round(time, 1) if ties else time), event, covariates


def _breslow_log_likelihood(time: np.ndarray, event: np.ndarray, covariates: np.ndarray, beta: np.ndarray) -> float:
    """The textbook sum over events of eta_i - log(sum of exp(eta_j) over everyone with time_j >= time_i)."""
    eta = covariates @ beta
    return float(sum(eta[i] - np.log(np.exp(eta[time >= time[i]]).sum()) for i in np.flatnonzero(event == 1)))


class TestCoxAgainstStatsmodels:
    @pytest.mark.parametrize("ties", [False, True], ids=["untied", "tied_times"])
    @pytest.mark.parametrize("seed", [11, 12, 13])
    def test_coefficients_standard_errors_likelihood_and_p_values_match_phreg(self, seed: int, ties: bool) -> None:
        time, event, covariates = _cox_data(seed=seed, ties=ties)
        reference = hazard_regression.PHReg(time, covariates, status=event, ties="breslow").fit()

        result = cox_regression(time.tolist(), event.tolist(), covariates)

        np.testing.assert_allclose(result["coefficients"], reference.params, rtol=1e-7, atol=1e-9)
        np.testing.assert_allclose(result["se"], reference.bse, rtol=1e-7)
        np.testing.assert_allclose(result["p_values"], reference.pvalues, rtol=1e-6)
        assert result["log_likelihood"] == pytest.approx(reference.llf, abs=1e-8)
        np.testing.assert_allclose(result["hazard_ratios"], np.exp(reference.params), rtol=1e-7)

    def test_heavy_ties_use_the_whole_tied_risk_set(self) -> None:
        rng = np.random.default_rng(3)
        covariates = rng.standard_normal((120, 2))
        time = rng.integers(1, 6, 120).astype(float)  # five distinct times for 120 subjects
        event = rng.binomial(1, 0.7, 120)
        reference = hazard_regression.PHReg(time, covariates, status=event, ties="breslow").fit()

        result = cox_regression(time.tolist(), event.tolist(), covariates)

        np.testing.assert_allclose(result["coefficients"], reference.params, rtol=1e-6, atol=1e-8)
        np.testing.assert_allclose(result["se"], reference.bse, rtol=1e-6)
        assert result["log_likelihood"] == pytest.approx(reference.llf, abs=1e-7)

    def test_one_binary_covariate(self) -> None:
        rng = np.random.default_rng(4)
        group = rng.integers(0, 2, 90).astype(float)
        time = rng.exponential(1.0 / np.exp(0.7 * group)).round(2)
        event = rng.binomial(1, 0.8, 90)
        reference = hazard_regression.PHReg(time, group[:, None], status=event, ties="breslow").fit()

        result = cox_regression(time.tolist(), event.tolist(), group.tolist())  # a flat list is one covariate

        np.testing.assert_allclose(result["coefficients"], reference.params, rtol=1e-7)
        np.testing.assert_allclose(result["se"], reference.bse, rtol=1e-7)


class TestCoxSelfConsistency:
    """Checks that need no other library: the likelihood is the textbook sum and the fit is a maximum."""

    def test_the_reported_likelihood_is_the_textbook_breslow_sum_at_the_fitted_coefficients(self) -> None:
        time, event, covariates = _cox_data(seed=5, n=80, ties=True)

        result = cox_regression(time.tolist(), event.tolist(), covariates)

        direct = _breslow_log_likelihood(time, event, covariates, np.array(result["coefficients"]))
        assert result["log_likelihood"] == pytest.approx(direct, abs=1e-9)
        assert result["log_likelihood"] < -50  # the old value dropped a shift of max(eta) per event

    def test_the_fit_is_a_maximum_of_the_likelihood(self) -> None:
        time, event, covariates = _cox_data(seed=6, n=100, ties=True)
        result = cox_regression(time.tolist(), event.tolist(), covariates)
        beta = np.array(result["coefficients"])
        best = _breslow_log_likelihood(time, event, covariates, beta)

        for direction in np.vstack([np.eye(3), -np.eye(3)]):
            for step in (1e-3, 1e-2, 0.1):
                assert _breslow_log_likelihood(time, event, covariates, beta + step * direction) < best

        gradient = [
            (
                _breslow_log_likelihood(time, event, covariates, beta + 1e-5 * e)
                - _breslow_log_likelihood(time, event, covariates, beta - 1e-5 * e)
            )
            / 2e-5
            for e in np.eye(3)
        ]
        np.testing.assert_allclose(gradient, 0.0, atol=1e-5)

    def test_the_information_matrix_gives_the_standard_errors(self) -> None:
        time, event, covariates = _cox_data(seed=7, n=100)
        result = cox_regression(time.tolist(), event.tolist(), covariates)
        beta = np.array(result["coefficients"])

        def log_lik(b: np.ndarray) -> float:
            return _breslow_log_likelihood(time, event, covariates, b)

        step = 1e-4
        hessian = np.zeros((3, 3))
        for i in range(3):
            for j in range(3):
                ei, ej = np.eye(3)[i] * step, np.eye(3)[j] * step
                hessian[i, j] = (
                    log_lik(beta + ei + ej)
                    - log_lik(beta + ei - ej)
                    - log_lik(beta - ei + ej)
                    + log_lik(beta - ei - ej)
                ) / (4 * step**2)
        np.testing.assert_allclose(result["se"], np.sqrt(np.diag(np.linalg.inv(-hessian))), rtol=1e-4)

    def test_p_values_are_the_two_sided_normal_tail_of_beta_over_se(self) -> None:
        time, event, covariates = _cox_data(seed=8)
        result = cox_regression(time.tolist(), event.tolist(), covariates)

        for beta, se, p in zip(result["coefficients"], result["se"], result["p_values"]):
            assert p == pytest.approx(math.erfc(abs(beta / se) / math.sqrt(2.0)), rel=1e-12)

    def test_shuffling_the_subjects_changes_nothing(self) -> None:
        time, event, covariates = _cox_data(seed=9, n=60, ties=True)
        order = np.random.default_rng(1).permutation(60)

        first = cox_regression(time.tolist(), event.tolist(), covariates)
        shuffled = cox_regression(time[order].tolist(), event[order].tolist(), covariates[order])

        np.testing.assert_allclose(first["coefficients"], shuffled["coefficients"], rtol=1e-9)
        assert first["log_likelihood"] == pytest.approx(shuffled["log_likelihood"], abs=1e-9)

    def test_large_covariate_values_do_not_overflow_the_likelihood(self) -> None:
        time, event, covariates = _cox_data(seed=10, n=100)
        scaled = covariates * np.array([1.0, 1.0, 1.0]) + 400.0  # eta around 400 would overflow exp() unshifted

        result = cox_regression(time.tolist(), event.tolist(), scaled)
        centred = cox_regression(time.tolist(), event.tolist(), covariates)

        np.testing.assert_allclose(result["coefficients"], centred["coefficients"], rtol=1e-6)
        assert math.isfinite(result["log_likelihood"])


class TestCoxRefusesWhatItCannotFit:
    def test_collinear_covariates_raise_instead_of_being_patched_with_a_ridge(self) -> None:
        time, event, covariates = _cox_data(n=80)
        collinear = np.column_stack([covariates[:, 0], covariates[:, 1], covariates[:, 0] + covariates[:, 1]])

        with pytest.raises(ValueError, match=r"information matrix is singular"):
            cox_regression(time.tolist(), event.tolist(), collinear)

    def test_a_constant_covariate_raises(self) -> None:
        time, event, covariates = _cox_data(n=80)
        covariates[:, 1] = 3.0

        with pytest.raises(ValueError, match=r"singular"):
            cox_regression(time.tolist(), event.tolist(), covariates)

    def test_complete_separation_raises_instead_of_returning_an_unconverged_coefficient(self) -> None:
        # everyone with x = 1 fails before everyone with x = 0: the coefficient is infinite
        time = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0]
        event = [1, 1, 1, 1, 1, 1, 1, 1]
        x = [1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0]

        with pytest.raises(ValueError, match=r"did not converge"):
            cox_regression(time, event, x)

    def test_too_many_covariates_for_the_events_raise(self) -> None:
        # two events cannot determine six coefficients: the fit runs away (or the information is singular)
        rng = np.random.default_rng(2)
        with pytest.raises(ValueError, match=r"singular|did not converge"):
            cox_regression(rng.exponential(size=12).tolist(), [1, 1] + [0] * 10, rng.standard_normal((12, 6)))

    @pytest.mark.parametrize("bad", [2, -1, 0.5, True and 3])
    def test_event_indicators_must_be_zero_or_one(self, bad: Any) -> None:
        with pytest.raises(ValueError, match=r"only 0 \(censored\) and 1 \(event\)"):
            cox_regression([1.0, 2.0, 3.0, 4.0], [1, bad, 0, 1], np.arange(4.0))

    @pytest.mark.parametrize("bad", [float("nan"), float("inf")])
    def test_times_and_covariates_must_be_finite(self, bad: float) -> None:
        with pytest.raises(ValueError, match=r"time must be finite"):
            cox_regression([1.0, bad, 3.0, 4.0], [1, 1, 0, 1], np.arange(4.0))
        with pytest.raises(ValueError, match=r"covariates must be finite"):
            cox_regression([1.0, 2.0, 3.0, 4.0], [1, 1, 0, 1], np.array([0.0, bad, 1.0, 2.0]))


def _two_group_data(seed: int = 2, n: int = 90) -> tuple[list[float], list[int]]:
    rng = np.random.default_rng(seed)
    time = np.round(rng.exponential(5.0, n), 1)
    event = rng.binomial(1, 0.75, n)
    return time.tolist(), event.tolist()


class TestLogRankAgainstStatsmodels:
    @pytest.mark.parametrize("n_groups", [2, 3, 4, 5])
    @pytest.mark.parametrize("ties", [False, True], ids=["untied", "tied_times"])
    def test_the_chi_squared_and_p_value_match_survdiff(self, n_groups: int, ties: bool) -> None:
        rng = np.random.default_rng(20 + n_groups)
        time, event, _ = _cox_data(seed=n_groups, n=150, ties=ties)
        groups = rng.integers(0, n_groups, 150)
        reference_chi2, reference_p = survfunc.survdiff(time, event, groups)

        result = log_rank_test(time.tolist(), event.tolist(), groups.tolist())

        assert result["df"] == n_groups - 1
        assert result["chi2"] == pytest.approx(reference_chi2, rel=1e-9)
        assert result["p_value"] == pytest.approx(reference_p, rel=1e-8)

    def test_observed_and_expected_events_each_add_up_to_the_number_of_events(self) -> None:
        time, event, _ = _cox_data(seed=3, n=120, ties=True)
        groups = np.random.default_rng(3).integers(0, 3, 120)

        result = log_rank_test(time.tolist(), event.tolist(), groups.tolist())

        per_group = result["observed_expected_per_group"]
        assert sum(v["observed"] for v in per_group.values()) == int(event.sum())
        assert sum(v["expected"] for v in per_group.values()) == pytest.approx(float(event.sum()), rel=1e-12)

    def test_group_labels_may_be_any_sortable_values_and_their_order_does_not_matter(self) -> None:
        time, event = _two_group_data()
        numeric = [i % 3 for i in range(len(time))]
        labels = [("low", "mid", "high")[g] for g in numeric]

        by_number = log_rank_test(time, event, numeric)
        by_name = log_rank_test(time, event, labels)  # type: ignore[arg-type]
        shuffled = log_rank_test(time[::-1], event[::-1], numeric[::-1])

        assert by_name["chi2"] == pytest.approx(by_number["chi2"], rel=1e-12)
        assert shuffled["chi2"] == pytest.approx(by_number["chi2"], rel=1e-12)

    def test_a_group_never_at_risk_at_an_event_time_is_left_out_with_its_degree_of_freedom(self) -> None:
        # group 2 is censored at t = 0.5, before the first event (t = 1): it can say nothing about the others
        base_time, base_event = _two_group_data(seed=5, n=60)
        time = [t + 1.0 for t in base_time] + [0.5] * 6
        event = base_event + [0] * 6
        groups = [i % 2 for i in range(60)] + [2] * 6
        two_groups = log_rank_test(time[:60], event[:60], groups[:60])

        three_labels = log_rank_test(time, event, groups)

        assert three_labels["df"] == 1  # not 2
        assert three_labels["chi2"] == pytest.approx(two_groups["chi2"], rel=1e-12)
        assert three_labels["p_value"] == pytest.approx(two_groups["p_value"], rel=1e-12)
        assert three_labels["observed_expected_per_group"][2] == {"observed": 0.0, "expected": 0.0}

    def test_nothing_to_compare_raises_instead_of_reporting_p_equal_one(self) -> None:
        # the only event happens when only group 0 is still at risk
        with pytest.raises(ValueError, match=r"cannot be compared"):
            log_rank_test([1.0, 1.0, 5.0, 6.0], [0, 0, 1, 0], [1, 1, 0, 0])

    def test_two_groups_still_give_the_classical_statistic(self) -> None:
        time, event = _two_group_data(seed=8, n=70)
        groups = [i % 2 for i in range(70)]
        reference_chi2, reference_p = survfunc.survdiff(np.array(time), np.array(event), np.array(groups))

        result = log_rank_test(time, event, groups)

        assert result["chi2"] == pytest.approx(reference_chi2, rel=1e-9)
        assert result["p_value"] == pytest.approx(reference_p, rel=1e-8)


class TestSurvivalConventionsAreDocumented:
    def test_the_lasso_penalty_choice_is_stated(self) -> None:
        from metainformant.multiomics.survival.analysis import multi_omic_survival_model

        doc = " ".join((multi_omic_survival_model.__doc__ or "").split())

        assert "lambda_max / 10" in doc and "not tuned" in doc

    def test_the_no_pair_concordance_convention_is_stated_and_holds(self) -> None:
        from metainformant.multiomics.survival.analysis import compute_concordance_index

        doc = " ".join((compute_concordance_index.__doc__ or "").split())
        result = compute_concordance_index([1.0, 2.0, 3.0], [1.0, 2.0, 3.0], [0, 0, 0])  # no event, no comparable pair

        assert "reported as 0.5" in doc and "a convention, not an estimate" in doc
        assert result["c_index"] == 0.5 and result["se"] == 0.0


class TestRiskStratificationPropagatesTheComparisonError:
    def test_risk_groups_that_cannot_be_compared_raise_instead_of_reporting_p_one(self) -> None:
        # the low-risk subjects are all censored before the only events, so no event has two groups at risk
        from metainformant.multiomics.survival.analysis import risk_stratification

        scores = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0]  # the four lowest scores form the low-risk group
        time = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0]
        event = [0, 0, 0, 0, 1, 1]

        with pytest.raises(ValueError, match=r"cannot be compared"):
            risk_stratification(scores, time, event, n_groups=2)
