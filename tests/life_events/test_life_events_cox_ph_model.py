"""Regression tests for ``survival.time_to_event.cox_ph_model`` (Newton-Raphson on the Breslow likelihood).

The function was documented as Newton-Raphson on the partial likelihood but

* stepped every coefficient by ``gradient / |Hessian diagonal|`` (clipped to +-1), ignoring the off-diagonal
  terms, so correlated covariates converged to the wrong place or not at all and exactly collinear covariates
  produced nonsense coefficients (a constant covariate got a coefficient of about 3);
* treated tied event times sequentially (each tied event dropped the earlier tied events from its risk set)
  instead of Breslow's shared risk set;
* took the standard errors from the diagonal of the Hessian alone, which understates them whenever
  covariates are correlated;
* accepted NaN/inf values and event codes other than 0/1 silently.

It now runs on the NumPy Cox core shared with ``SurvivalPredictor``. The fit is checked against
statsmodels' ``PHReg`` (Breslow ties), against a brute-force partial likelihood, and against a tiny example whose
maximiser and curvature are derived by hand.
"""

from __future__ import annotations

import json
import math
from typing import Dict, List, Sequence, Tuple

import numpy as np
import pytest
from scipy import optimize

from metainformant.life_events.models.statistical_models import SurvivalPredictor
from metainformant.life_events.survival.time_to_event import cox_ph_model

Data = Tuple[np.ndarray, np.ndarray, np.ndarray]  # times, events, covariates

TRUE_BETA = np.array([0.7, -0.4, 0.2])


def simulate(n: int = 400, seed: int = 3, bin_width: float = 0.0, correlated: bool = False) -> Data:
    """Right-censored exponential data; ``bin_width`` rounds the times up to multiples of it, creating ties."""
    rng = np.random.default_rng(seed)
    first = rng.normal(size=n)
    second = 0.9 * first + 0.3 * rng.normal(size=n) if correlated else rng.normal(size=n)
    third = rng.poisson(1.0, size=n).astype(float)
    covariates = np.column_stack([first, second, third])
    event_times = rng.exponential(1.0 / (0.1 * np.exp(covariates @ TRUE_BETA)))
    censoring_times = rng.exponential(15.0, size=n)
    times = np.minimum(event_times, censoring_times)
    if bin_width:
        times = np.ceil(times / bin_width) * bin_width
    return times, (event_times <= censoring_times).astype(int), covariates


SCENARIOS: Dict[str, Data] = {
    "continuous": simulate(),
    "tied": simulate(bin_width=1.0),
    "heavily tied": simulate(bin_width=4.0),
    "correlated": simulate(correlated=True),
    "correlated and tied": simulate(correlated=True, bin_width=1.0),
}


def fit(data: Data):
    times, events, covariates = data
    return cox_ph_model(times.tolist(), events.tolist(), covariates.tolist())


def statsmodels_fit(data: Data):
    hazard = pytest.importorskip("statsmodels.duration.hazard_regression")
    times, events, covariates = data
    return hazard.PHReg(times, covariates, status=events, ties="breslow").fit(disp=False)


def breslow_negative_log_likelihood(beta: np.ndarray, data: Data) -> float:
    """The Breslow partial likelihood straight from its definition (every tied event shares the risk set)."""
    times, events, covariates = data
    eta = covariates @ beta
    total = 0.0
    for i in np.where(events == 1)[0]:
        total += eta[i] - np.log(np.exp(eta[times >= times[i]]).sum())
    return float(-total)


# --- against statsmodels' PHReg (Breslow ties) ------------------------------------------------------------------


@pytest.mark.parametrize("scenario", list(SCENARIOS))
def test_coefficients_match_statsmodels(scenario: str):
    data = SCENARIOS[scenario]
    if "tied" in scenario:
        assert len(np.unique(data[0][data[1] == 1])) < 0.35 * data[1].sum()  # many events share a time

    result = fit(data)

    np.testing.assert_allclose(result["coefficients"], statsmodels_fit(data).params, rtol=1e-6, atol=1e-7)


@pytest.mark.parametrize("scenario", list(SCENARIOS))
def test_standard_errors_match_statsmodels(scenario: str):
    data = SCENARIOS[scenario]

    result = fit(data)

    np.testing.assert_allclose(result["se"], statsmodels_fit(data).bse, rtol=1e-6)


@pytest.mark.parametrize("scenario", list(SCENARIOS))
def test_p_values_and_hazard_ratios_match_statsmodels(scenario: str):
    data = SCENARIOS[scenario]
    reference = statsmodels_fit(data)

    result = fit(data)

    np.testing.assert_allclose(result["p_values"], reference.pvalues, atol=1e-7)
    np.testing.assert_allclose(result["hazard_ratios"], np.exp(reference.params), rtol=1e-6)


@pytest.mark.parametrize("scenario", list(SCENARIOS))
def test_partial_log_likelihood_matches_statsmodels(scenario: str):
    data = SCENARIOS[scenario]

    result = fit(data)

    assert result["log_likelihood"] == pytest.approx(statsmodels_fit(data).llf, rel=1e-9)
    assert result["converged"] is True
    assert result["n_iterations"] <= 10  # Newton-Raphson converges quadratically


# --- independent oracles ------------------------------------------------------------------------------------------


@pytest.mark.parametrize("scenario", ["tied", "heavily tied", "correlated and tied"])
def test_coefficients_maximise_the_brute_force_breslow_likelihood(scenario: str):
    data = SCENARIOS[scenario]
    times, events, covariates = data
    subset = (times[:250], events[:250], covariates[:250])

    result = fit(subset)

    reference = optimize.minimize(
        breslow_negative_log_likelihood, np.zeros(3), args=(subset,), method="BFGS", options={"gtol": 1e-9}
    )
    np.testing.assert_allclose(result["coefficients"], reference.x, atol=1e-5)


def test_tied_event_times_share_one_risk_set():
    """Two events tie at t=1. Breslow's partial likelihood of this data is

    ``l(b) = 2b - 2 log 2 - 3 log(1 + e^b)``, maximal where ``e^b / (1 + e^b) = 2/3``, that is at ``b = ln 2``,
    with curvature ``-3 e^b / (1 + e^b)^2 = -2/3``; the old sequential tie rule gave 0.94 instead.
    """
    result = cox_ph_model([1.0, 1.0, 2.0, 3.0], [1, 1, 1, 0], [[1.0], [0.0], [1.0], [0.0]])

    assert result["coefficients"][0] == pytest.approx(math.log(2), abs=1e-9)
    assert result["se"][0] == pytest.approx(math.sqrt(1.5), abs=1e-9)
    assert result["log_likelihood"] == pytest.approx(-3 * math.log(3), abs=1e-9)
    assert result["hazard_ratios"][0] == pytest.approx(2.0, abs=1e-8)


def test_standard_errors_are_the_inverse_of_the_full_numerical_hessian():
    """Correlated covariates: the standard errors come from the whole information matrix, not its diagonal."""
    data = SCENARIOS["correlated and tied"]
    times, events, covariates = data
    subset = (times[:250], events[:250], covariates[:250])
    result = fit(subset)
    beta = np.array(result["coefficients"])

    step = 1e-4
    hessian = np.zeros((3, 3))
    for i in range(3):
        for j in range(3):
            shifts = [(step, step, 1), (step, -step, -1), (-step, step, -1), (-step, -step, 1)]
            for di, dj, sign in shifts:
                point = beta.copy()
                point[i] += di
                point[j] += dj
                hessian[i, j] += sign * breslow_negative_log_likelihood(point, subset)
    hessian /= 4 * step**2
    full = np.sqrt(np.diag(np.linalg.inv(hessian)))
    diagonal_only = 1.0 / np.sqrt(np.diag(hessian))

    np.testing.assert_allclose(result["se"], full, rtol=1e-3)
    assert np.abs(diagonal_only / full - 1).max() > 0.2  # the data are correlated enough for the difference to matter


# --- identifiability -----------------------------------------------------------------------------------------------


def test_collinear_covariates_are_not_identified_and_do_not_disturb_the_others():
    times, events, covariates = SCENARIOS["continuous"]
    first, third = covariates[:, 0], covariates[:, 2]
    duplicated = np.column_stack([first, first, third])  # the same covariate twice, plus an independent one
    reduced = statsmodels_fit((times, events, np.column_stack([first, third])))

    result = cox_ph_model(times.tolist(), events.tolist(), duplicated.tolist())

    # Only the sum of the two duplicated coefficients is determined: no finite standard error for either.
    assert result["se"][0] == math.inf and result["se"][1] == math.inf
    assert result["p_values"][0] == 1.0 and result["p_values"][1] == 1.0
    assert result["coefficients"][0] + result["coefficients"][1] == pytest.approx(reduced.params[0], abs=1e-6)
    # The independent covariate is estimated exactly as if the duplicate were not there.
    assert result["coefficients"][2] == pytest.approx(reduced.params[1], abs=1e-6)
    assert result["se"][2] == pytest.approx(reduced.bse[1], rel=1e-6)


def test_a_constant_covariate_is_not_identified():
    times, events, covariates = SCENARIOS["continuous"]
    with_constant = np.column_stack([covariates[:, 0], np.full(len(times), 5.0), covariates[:, 2]])
    reduced = statsmodels_fit((times, events, covariates[:, [0, 2]]))

    result = cox_ph_model(times.tolist(), events.tolist(), with_constant.tolist())

    assert result["coefficients"][1] == pytest.approx(0.0, abs=1e-9)  # the old fit reported about 2.9
    assert result["se"][1] == math.inf and result["p_values"][1] == 1.0
    assert result["coefficients"][0] == pytest.approx(reduced.params[0], abs=1e-6)
    assert result["coefficients"][2] == pytest.approx(reduced.params[1], abs=1e-6)


def test_data_without_events_has_nothing_to_estimate():
    result = cox_ph_model([1.0, 2.0, 3.0], [0, 0, 0], [[1.0], [2.0], [3.0]])

    assert result["coefficients"] == [0.0]
    assert result["se"] == [math.inf] and result["p_values"] == [1.0]
    assert result["concordance"] == 0.5
    assert result["n_events"] == 0 and result["n_subjects"] == 3


def test_a_covariate_that_orders_the_event_times_perfectly_has_a_huge_standard_error():
    """Monotone likelihood: the coefficient runs off until the likelihood stops improving."""
    result = cox_ph_model([10.0, 20.0, 30.0, 40.0], [1, 1, 1, 0], [[4.0], [3.0], [2.0], [1.0]])

    assert result["coefficients"][0] > 15
    assert result["se"][0] > 1e3
    assert result["p_values"][0] > 0.99
    assert result["concordance"] == 1.0


# --- the shared core, the return value and the inputs --------------------------------------------------------------


def test_the_fit_is_the_unpenalised_survival_predictor_cox_model():
    rng = np.random.default_rng(5)
    tokens = ["a", "b", "c", "d"]
    sequences = [[str(t) for t in rng.choice(tokens, size=int(rng.integers(1, 6)))] for _ in range(300)]
    counts = np.array([[sequence.count(token) for token in tokens] for sequence in sequences], dtype=float)
    event_times = rng.exponential(1 / (0.1 * np.exp(counts @ np.array([0.6, -0.3, 0.0, 0.2]))))
    censoring_times = rng.exponential(15.0, size=300)
    times, died = np.minimum(event_times, censoring_times), event_times <= censoring_times

    function = cox_ph_model(times.tolist(), died.astype(int).tolist(), counts.tolist(), covariate_names=tokens)
    predictor = SurvivalPredictor(method="cox", penalizer=0.0).fit(sequences, times, died)

    assert function["covariate_names"] == tokens
    np.testing.assert_allclose(function["coefficients"], [predictor.coefficients[t] for t in tokens], atol=1e-10)
    np.testing.assert_allclose(function["se"], [predictor.standard_errors[t] for t in tokens], atol=1e-10)


def test_return_contract():
    times, events, covariates = SCENARIOS["tied"]

    result = cox_ph_model(times.tolist(), events.tolist(), covariates.tolist(), covariate_names=["a", "b", "c"])

    documented = {
        "coefficients",
        "hazard_ratios",
        "se",
        "p_values",
        "concordance",
        "covariate_names",
        "n_subjects",
        "n_events",
    }
    assert documented <= set(result)
    assert set(result) - documented == {"log_likelihood", "converged", "n_iterations"}
    for key in ("coefficients", "hazard_ratios", "se", "p_values"):
        assert isinstance(result[key], list) and len(result[key]) == 3
        assert all(isinstance(value, float) for value in result[key])
    assert result["covariate_names"] == ["a", "b", "c"]
    assert result["n_subjects"] == len(times) and result["n_events"] == int(events.sum())
    assert isinstance(result["n_subjects"], int) and isinstance(result["n_events"], int)
    assert result["hazard_ratios"] == pytest.approx([math.exp(b) for b in result["coefficients"]])
    assert all(0.0 <= p <= 1.0 for p in result["p_values"])
    assert 0.5 < result["concordance"] <= 1.0
    assert json.loads(json.dumps(result)) == result  # plain Python values, so results can be written as JSON


def test_covariate_names_default_to_x_indices():
    result = cox_ph_model([1.0, 2.0, 3.0, 4.0], [1, 1, 0, 1], [[1.0, 0.0], [0.0, 1.0], [1.0, 1.0], [0.0, 0.0]])

    assert result["covariate_names"] == ["X0", "X1"]


def test_numpy_arrays_and_boolean_events_are_accepted():
    times, events, covariates = SCENARIOS["tied"]

    from_lists = cox_ph_model(times.tolist(), events.tolist(), covariates.tolist())
    from_arrays = cox_ph_model(times, events.astype(bool), covariates)

    assert from_arrays == from_lists


def test_concordance_matches_a_literal_double_loop():
    """Pairs are comparable when the subject with the event has the earlier time; tied risk scores are skipped."""
    rng = np.random.default_rng(2)
    n = 200
    covariates = rng.integers(0, 3, size=(n, 2)).astype(float)  # discrete, so risk scores tie often
    times = np.ceil(rng.exponential(5.0, size=n) / 2.0) * 2.0  # and event times tie as well
    events = (rng.random(n) < 0.7).astype(int)

    result = cox_ph_model(times.tolist(), events.tolist(), covariates.tolist())

    risk = covariates @ np.array(result["coefficients"])
    concordant = discordant = 0
    for i in range(n):
        if events[i] != 1:
            continue
        for j in range(n):
            if times[j] > times[i]:
                concordant += risk[i] > risk[j]
                discordant += risk[i] < risk[j]
    assert (risk[:, None] == risk[None, :]).any()  # the data really do contain tied risk scores
    assert result["concordance"] == pytest.approx(concordant / (concordant + discordant), abs=1e-12)


def test_the_concordance_of_an_informative_model_is_high():
    times, events, covariates = SCENARIOS["continuous"]

    result = fit((times, events, covariates))

    assert 0.62 < result["concordance"] < 0.9


@pytest.mark.parametrize(
    "times, events, covariates, message",
    [
        ([1.0, 2.0], [1], [[1.0], [2.0]], "same length"),
        ([1.0, 2.0], [1, 0], [[1.0]], "same length"),
        ([1.0, 2.0], [1, 0], [[], []], "At least one covariate"),
        ([], [], [], "At least one covariate"),
        ([1.0, float("nan"), 3.0], [1, 1, 0], [[1.0], [2.0], [3.0]], "finite"),
        ([1.0, 2.0, float("inf")], [1, 1, 0], [[1.0], [2.0], [3.0]], "finite"),
        ([1.0, 2.0, 3.0], [1, 1, 0], [[1.0], [float("nan")], [3.0]], "finite"),
        ([1.0, 2.0, 3.0], [1, 2, 0], [[1.0], [2.0], [3.0]], "only 0"),
        ([1.0, 2.0, 3.0], [1, 0.5, 0], [[1.0], [2.0], [3.0]], "only 0"),
        ([1.0, 2.0, 3.0], [1, 1, 0], [[1.0, 2.0], [2.0], [3.0, 1.0]], "same covariates"),
        ([1.0, 2.0, 3.0], [1, 1, 0], [1.0, 2.0, 3.0], "one covariate vector per subject"),
        (["a", "b", "c"], [1, 1, 0], [[1.0], [2.0], [3.0]], "numeric"),
    ],
)
def test_invalid_input_is_rejected(times: Sequence, events: Sequence, covariates: Sequence, message: str):
    with pytest.raises(ValueError, match=message):
        cox_ph_model(list(times), list(events), list(covariates))


def test_the_inputs_are_not_modified():
    times, events, covariates = [3.0, 1.0, 2.0, 5.0, 4.0], [1, 1, 0, 1, 0], [[1.0], [2.0], [0.5], [3.0], [1.5]]
    snapshot: List = [list(times), list(events), [list(row) for row in covariates]]

    cox_ph_model(times, events, covariates)

    assert [times, events, covariates] == snapshot


def test_a_large_cohort_is_fitted_in_a_few_newton_steps():
    """The NumPy fit and the vectorised concordance handle 20 000 subjects (the old fit looped over pairs in Python)."""
    rng = np.random.default_rng(1)
    n = 20_000
    covariates = rng.normal(size=(n, 5))
    event_times = rng.exponential(1 / (0.1 * np.exp(covariates @ np.array([0.5, -0.3, 0.2, 0.0, 0.1]))))
    censoring_times = rng.exponential(15.0, size=n)
    times, events = np.minimum(event_times, censoring_times), (event_times <= censoring_times).astype(int)

    result = cox_ph_model(times, events, covariates)

    assert result["converged"] and result["n_iterations"] <= 6
    np.testing.assert_allclose(result["coefficients"], [0.5, -0.3, 0.2, 0.0, 0.1], atol=0.06)
