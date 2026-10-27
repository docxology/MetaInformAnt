"""Regression tests for ``SurvivalPredictor`` (a real Cox proportional hazards fit in NumPy).

Defects covered:

* Without PyTorch, ``fit`` logged a warning, set ``is_fitted = True`` and returned without fitting
  anything, so the model looked trained although it never saw the data.
* The NumPy-only methods (exponential, Kaplan-Meier) were gated behind ``import torch`` as well.
* ``predict`` returned a constant 365.0 for every sequence whenever there was no model.
* The ``cox`` method trained a network with mean squared error on ``log(time + 1)``: not a Cox
  likelihood, and it ignored censoring.
* Unknown methods silently fell back to a Kaplan-Meier "baseline" whose survival curve was
  ``linspace(1, 0)``.

The Cox fit is checked against independent oracles: statsmodels' ``PHReg`` (Breslow ties), a brute-force
partial likelihood maximised by SciPy, the hand-computed Kaplan-Meier curve and the package's pure-Python
Kaplan-Meier estimator. ``lifelines`` is not a dependency of this project and is not installed in the
development environment, so it is not used as an oracle.
"""

from __future__ import annotations

import sys
from datetime import datetime
from typing import Any, Dict, List, Sequence, Tuple

import numpy as np
import pytest
from scipy import optimize
from scipy.integrate import trapezoid

from metainformant.life_events.core.events import Event, EventSequence
from metainformant.life_events.core.utils import generate_synthetic_life_events
from metainformant.life_events.models.statistical_models import SurvivalPredictor
from metainformant.life_events.survival.time_to_event import kaplan_meier_estimator

VOCABULARY = [
    "health:diagnosis",
    "occupation:job_change",
    "income:raise",
    "family:marriage",
    "education:degree",
]
# True log hazard ratio of one more occurrence of each event, and the hazard of a sequence without events.
TRUE_LOG_HAZARD: Dict[str, float] = {
    "health:diagnosis": 0.8,
    "occupation:job_change": -0.5,
    "income:raise": 0.0,
    "family:marriage": 0.3,
    "education:degree": 0.0,
}
BASE_RATE = 0.1
METHODS = ["cox", "exponential", "baseline"]


def simulate(
    n: int, seed: int, censoring_rate: float = 0.04, rounding: bool = False
) -> Tuple[List[List[str]], np.ndarray, np.ndarray]:
    """Cohort with exponential event times whose hazard is ``BASE_RATE * exp(sum of TRUE_LOG_HAZARD)``.

    Censoring times are independent exponentials. With ``rounding`` the times are rounded up to whole
    units, which creates many tied event times.
    """
    rng = np.random.default_rng(seed)
    sequences = [[str(token) for token in rng.choice(VOCABULARY, size=int(rng.integers(1, 7)))] for _ in range(n)]
    log_hazard = np.array([sum(TRUE_LOG_HAZARD[token] for token in sequence) for sequence in sequences])
    event_times = rng.exponential(1.0 / (BASE_RATE * np.exp(log_hazard)))
    censoring_times = rng.exponential(1.0 / censoring_rate, size=n)
    times = np.minimum(event_times, censoring_times)
    if rounding:
        times = np.ceil(times)
    return sequences, times, event_times <= censoring_times


def count_matrix(sequences: Sequence[Sequence[str]]) -> np.ndarray:
    return np.array([[sequence.count(token) for token in VOCABULARY] for sequence in sequences], dtype=float)


def coefficient_vector(model: SurvivalPredictor) -> np.ndarray:
    return np.array([model.coefficients.get(token, 0.0) for token in VOCABULARY])


def as_event_sequences(token_lists: List[List[str]]) -> List[EventSequence]:
    sequences = []
    for index, tokens in enumerate(token_lists):
        events = []
        for position, token in enumerate(tokens):
            domain, event_type = token.split(":")
            events.append(Event(event_type, datetime(2000 + position, 1, 1), domain))
        sequences.append(EventSequence(f"person_{index}", events))
    return sequences


def breslow_negative_log_likelihood(
    beta: np.ndarray, features: np.ndarray, times: np.ndarray, died: np.ndarray
) -> float:
    """The Breslow partial likelihood straight from its definition (every tied event shares the risk set)."""
    eta = features @ beta
    total = 0.0
    for i in np.where(died)[0]:
        total += eta[i] - np.log(np.exp(eta[times >= times[i]]).sum())
    return float(-total)


@pytest.fixture(scope="module")
def cohort() -> Tuple[List[List[str]], np.ndarray, np.ndarray]:
    return simulate(2000, seed=21)


# --- the Cox fit against independent oracles -----------------------------------------------------------------


@pytest.mark.parametrize("rounding", [False, True], ids=["continuous times", "tied times"])
def test_cox_coefficients_match_statsmodels(rounding: bool):
    statsmodels_hazard = pytest.importorskip("statsmodels.duration.hazard_regression")
    sequences, times, died = simulate(400, seed=5, rounding=rounding)
    if rounding:
        assert len(np.unique(times[died])) < 0.2 * died.sum()  # many ties, so Breslow's rule matters

    model = SurvivalPredictor(method="cox", penalizer=0.0).fit(sequences, times, died)

    reference = statsmodels_hazard.PHReg(times, count_matrix(sequences), status=died.astype(int), ties="breslow").fit(
        disp=False
    )
    np.testing.assert_allclose(coefficient_vector(model), reference.params, rtol=1e-7, atol=1e-9)
    np.testing.assert_allclose([model.standard_errors[token] for token in VOCABULARY], reference.bse, rtol=1e-6)


@pytest.mark.parametrize("rounding", [False, True], ids=["continuous times", "tied times"])
def test_cox_maximises_the_brute_force_partial_likelihood(rounding: bool):
    sequences, times, died = simulate(250, seed=8, censoring_rate=0.08, rounding=rounding)
    features = count_matrix(sequences)

    model = SurvivalPredictor(method="cox", penalizer=0.0).fit(sequences, times, died)

    reference = optimize.minimize(
        breslow_negative_log_likelihood,
        np.zeros(len(VOCABULARY)),
        args=(features, times, died),
        method="BFGS",
        options={"gtol": 1e-9},
    )
    np.testing.assert_allclose(coefficient_vector(model), reference.x, atol=1e-5)
    # No direction improves the likelihood: it is a maximum, not merely close to one.
    best = breslow_negative_log_likelihood(coefficient_vector(model), features, times, died)
    rng = np.random.default_rng(0)
    for _ in range(20):
        nudged = coefficient_vector(model) + rng.normal(0, 0.01, len(VOCABULARY))
        assert breslow_negative_log_likelihood(nudged, features, times, died) >= best


@pytest.mark.parametrize("rounding", [False, True], ids=["continuous times", "tied times"])
def test_survival_function_uses_the_breslow_baseline_hazard(rounding: bool):
    sequences, times, died = simulate(300, seed=3, rounding=rounding)
    features = count_matrix(sequences)
    model = SurvivalPredictor(method="cox", penalizer=0.0).fit(sequences, times, died)
    beta = coefficient_vector(model)

    # Breslow's estimate from its definition: sum over event times u <= t of deaths(u) / sum_{t_j >= u} exp(x_j beta).
    weights = np.exp(features @ beta)
    event_times = np.unique(times[died])
    increments = np.array([np.sum((times == u) & died) / weights[times >= u].sum() for u in event_times])
    grid = np.array([0.0, event_times[0], 3.0, 10.0, 25.0, event_times[-1], event_times[-1] + 100.0])
    probes = [
        ["health:diagnosis"],
        ["occupation:job_change", "occupation:job_change"],
        [],
        ["family:marriage", "income:raise"],
    ]

    expected = np.array(
        [
            [
                np.exp(-increments[event_times <= t].sum() * np.exp(sum(beta[VOCABULARY.index(tok)] for tok in probe)))
                for t in grid
            ]
            for probe in probes
        ]
    )
    np.testing.assert_allclose(model.predict_survival_function(probes, grid), expected, rtol=1e-9, atol=1e-12)


def test_newton_raphson_converges_in_a_few_iterations():
    sequences, times, died = simulate(400, seed=5)

    model = SurvivalPredictor(method="cox", penalizer=0.0).fit(sequences, times, died)

    assert model.model["converged"] is True
    assert model.model["n_iter"] <= 10  # quadratic convergence: the old gradient-descent stand-in needed hundreds


# --- recovering a known hazard ratio -------------------------------------------------------------------------


@pytest.mark.parametrize("method", ["cox", "exponential"])
def test_the_fit_recovers_the_known_log_hazard_ratios(method: str, cohort):
    sequences, times, died = cohort
    assert 0.2 < 1 - died.mean() < 0.5  # a realistic share of the cohort is censored

    model = SurvivalPredictor(method=method).fit(sequences, times, died)

    for token, truth in TRUE_LOG_HAZARD.items():
        assert model.coefficients[token] == pytest.approx(truth, abs=0.15), token


@pytest.mark.parametrize("method", ["cox", "exponential"])
def test_hazard_ratios_implied_by_the_survival_functions_match_the_truth(method: str, cohort):
    """The same check through the prediction interface alone: for a proportional hazards model

    ``log(-log S(t | token)) - log(-log S(t | no events))`` equals the log hazard ratio of the token at every t.
    """
    sequences, times, died = cohort
    model = SurvivalPredictor(method=method).fit(sequences, times, died)

    for time_point in (3.0, 8.0):
        grid = np.array([time_point])
        baseline = model.predict_survival_function([[]], grid)[0, 0]
        for token, truth in TRUE_LOG_HAZARD.items():
            with_token = model.predict_survival_function([[token]], grid)[0, 0]
            implied = np.log(np.log(with_token) / np.log(baseline))
            assert implied == pytest.approx(truth, abs=0.15), (token, time_point)


def test_cox_and_exponential_models_agree_on_exponential_data(cohort):
    sequences, times, died = cohort

    cox = SurvivalPredictor(method="cox").fit(sequences, times, died)
    exponential = SurvivalPredictor(method="exponential").fit(sequences, times, died)

    np.testing.assert_allclose(coefficient_vector(cox), coefficient_vector(exponential), atol=0.06)


def test_the_baseline_hazard_is_calibrated(cohort):
    sequences, times, died = cohort
    probes = [
        [],
        ["health:diagnosis"],
        ["occupation:job_change", "occupation:job_change"],
        ["family:marriage", "income:raise"],
    ]
    grid = np.array([2.0, 5.0, 10.0, 20.0])
    truth = np.array(
        [
            [np.exp(-BASE_RATE * np.exp(sum(TRUE_LOG_HAZARD[token] for token in probe)) * t) for t in grid]
            for probe in probes
        ]
    )

    for method in ("cox", "exponential"):
        survival = SurvivalPredictor(method=method).fit(sequences, times, died).predict_survival_function(probes, grid)
        assert np.abs(survival - truth).max() < 0.05, method


def test_the_fit_recovers_the_hazard_ratio_on_generated_event_sequences():
    """The package's own generator produces EventSequence objects with 30+ distinct events."""
    sequences, _ = generate_synthetic_life_events(
        n_sequences=1500, min_events_per_sequence=3, max_events_per_sequence=14, random_state=3
    )
    diagnoses = np.array([sum(event.event_type == "diagnosis" for event in seq.events) for seq in sequences])
    truth = 0.9
    rng = np.random.default_rng(0)
    event_times = rng.exponential(1 / (0.05 * np.exp(truth * diagnoses)))
    censoring_times = rng.exponential(1 / 0.02, size=len(sequences))
    times, died = np.minimum(event_times, censoring_times), event_times <= censoring_times

    for method in ("cox", "exponential"):
        model = SurvivalPredictor(method=method).fit(sequences, times, died)
        coefficients = model.coefficients
        assert model.num_events is not None and model.num_events > 20
        assert coefficients["health:diagnosis"] == pytest.approx(truth, abs=0.2), method
        assert max(coefficients, key=lambda token: abs(coefficients[token])) == "health:diagnosis"


def test_censoring_is_used_not_ignored():
    """With 60% censoring, treating censored follow-up as events is badly biased; the Cox fit is not."""
    sequences, times, died = simulate(2000, seed=5, censoring_rate=0.25)
    assert 1 - died.mean() > 0.5

    correct = coefficient_vector(SurvivalPredictor(method="cox").fit(sequences, times, died))
    all_events = coefficient_vector(SurvivalPredictor(method="cox").fit(sequences, times, np.ones(len(times), bool)))

    truth = np.array([TRUE_LOG_HAZARD[token] for token in VOCABULARY])
    assert np.abs(correct - truth).max() < 0.2
    assert np.abs(all_events - truth).max() > 0.3
    assert np.abs(correct - all_events).max() > 0.3


def test_the_penalizer_shrinks_coefficients_and_zero_is_unpenalised(cohort):
    sequences, times, died = cohort

    fits = {
        penalizer: SurvivalPredictor(method="cox", penalizer=penalizer).fit(sequences, times, died).coefficients
        for penalizer in (0.0, 0.1, 100.0, 1e6)
    }

    key = "health:diagnosis"
    assert abs(fits[0.0][key]) > abs(fits[0.1][key]) > abs(fits[100.0][key]) > abs(fits[1e6][key])
    assert abs(fits[1e6][key]) < 1e-3
    assert fits[0.1][key] == pytest.approx(fits[0.0][key], abs=0.01)  # the default penalty is mild


def test_perfectly_predictive_events_stay_finite_with_the_default_penalizer():
    sequences = [["early"]] * 20 + [["late"]] * 20
    times = np.concatenate([np.linspace(1, 2, 20), np.linspace(10, 20, 20)])
    died = np.ones(40, dtype=bool)

    penalised = SurvivalPredictor(method="cox").fit(sequences, times, died)
    unpenalised = SurvivalPredictor(method="cox", penalizer=0.0).fit(sequences, times, died)

    assert abs(penalised.coefficients["early"]) < 5
    assert abs(unpenalised.coefficients["early"]) > abs(penalised.coefficients["early"])
    assert np.isfinite(penalised.predict(sequences)).all()


def test_collinear_counts_do_not_blow_up_even_without_a_penalty():
    rng = np.random.default_rng(0)
    sequences = []
    for _ in range(80):
        tokens = (["a", "b"] if rng.random() < 0.5 else []) + (["c"] if rng.random() < 0.5 else [])
        sequences.append(tokens or ["c", "c"])  # "a" and "b" always occur together
    hazard = np.array([0.8 * ("a" in tokens) for tokens in sequences])
    times = rng.exponential(1 / (0.1 * np.exp(hazard)))
    died = np.ones(80, dtype=bool)

    free = SurvivalPredictor(method="cox", penalizer=0.0).fit(sequences, times, died)
    tiny = SurvivalPredictor(method="cox", penalizer=1e-8).fit(sequences, times, died)

    assert max(abs(value) for value in free.coefficients.values()) < 10
    np.testing.assert_allclose(free.predict([["a", "b"], ["c"]]), tiny.predict([["a", "b"], ["c"]]), rtol=1e-4)


# --- predictions ----------------------------------------------------------------------------------------------


RISK_ORDERED_PROBES = [
    ["health:diagnosis", "health:diagnosis"],  # log hazard 1.6
    ["health:diagnosis"],  # 0.8
    ["family:marriage"],  # 0.3
    [],  # 0
    ["occupation:job_change"],  # -0.5
    ["occupation:job_change", "occupation:job_change"],  # -1.0
]


@pytest.mark.parametrize("method", ["cox", "exponential"])
def test_predicted_survival_times_are_ordered_by_risk(method: str, cohort):
    sequences, times, died = cohort

    predicted = SurvivalPredictor(method=method).fit(sequences, times, died).predict(RISK_ORDERED_PROBES)

    assert np.isfinite(predicted).all() and (predicted > 0).all()
    assert (np.diff(predicted) > 0).all()  # more hazard, shorter predicted survival


@pytest.mark.parametrize("method", METHODS)
def test_predictions_depend_on_the_sequences_for_covariate_models_only(method: str, cohort):
    sequences, times, died = cohort

    predicted = SurvivalPredictor(method=method).fit(sequences, times, died).predict(RISK_ORDERED_PROBES)

    if method == "baseline":
        assert np.ptp(predicted) == 0  # Kaplan-Meier has no covariates: every sequence gets the same curve
    else:
        assert np.ptp(predicted) > 1.0


@pytest.mark.parametrize("method", METHODS)
def test_predicted_time_is_the_area_under_the_survival_function(method: str, cohort):
    sequences, times, died = cohort
    model = SurvivalPredictor(method=method).fit(sequences, times, died)
    horizon = times.max()  # the longest follow-up, which is the documented horizon
    grid = np.linspace(0.0, horizon, 400_001)

    curves = model.predict_survival_function(RISK_ORDERED_PROBES, grid)

    np.testing.assert_allclose(model.predict(RISK_ORDERED_PROBES), trapezoid(curves, grid, axis=1), rtol=2e-3)


@pytest.mark.parametrize("method", METHODS)
def test_survival_functions_are_proper_curves(method: str, cohort):
    sequences, times, died = cohort
    model = SurvivalPredictor(method=method).fit(sequences, times, died)
    grid = np.array([0.0, 1.0, 5.0, 10.0, 30.0, 1000.0])

    curves = model.predict_survival_function(RISK_ORDERED_PROBES, grid)

    assert curves.shape == (len(RISK_ORDERED_PROBES), len(grid))
    assert np.all((curves >= 0) & (curves <= 1))
    assert np.all(curves[:, 0] == 1.0)  # nobody has failed at time zero
    assert np.all(np.diff(curves, axis=1) <= 0)  # survival never increases
    assert curves[0, 3] < 0.99  # and the curve actually falls
    if method != "baseline":
        assert np.all(np.diff(curves[:, 1:], axis=0) > 0)  # riskier sequences survive less at every later time


def test_survival_function_accepts_a_list_of_times_and_validates_them(cohort):
    sequences, times, died = cohort
    model = SurvivalPredictor(method="cox").fit(sequences, times, died)

    assert model.predict_survival_function([["health:diagnosis"]], [1, 5, 10]).shape == (1, 3)
    for bad in ([-1.0, 2.0], [np.nan], [[1.0, 2.0]]):
        with pytest.raises(ValueError, match="times"):
            model.predict_survival_function([["health:diagnosis"]], bad)


def test_events_unseen_in_training_are_ignored(cohort):
    sequences, times, died = cohort
    model = SurvivalPredictor(method="cox").fit(sequences, times, died)

    known = model.predict([["health:diagnosis"], []])
    with_unseen = model.predict([["health:diagnosis", "never:seen"], ["another:unseen"]])

    np.testing.assert_allclose(known, with_unseen)


def test_event_sequences_and_token_lists_are_equivalent(cohort):
    sequences, times, died = cohort
    from_tokens = SurvivalPredictor(method="cox").fit(sequences, times, died)
    from_events = SurvivalPredictor(method="cox").fit(as_event_sequences(sequences), times, died.astype(int))

    assert from_tokens.coefficients == from_events.coefficients
    probes = RISK_ORDERED_PROBES
    np.testing.assert_array_equal(from_tokens.predict(probes), from_events.predict(as_event_sequences(probes)))


def test_the_cox_fit_is_deterministic_and_ignores_random_state(cohort):
    sequences, times, died = cohort

    first = SurvivalPredictor(method="cox", random_state=1).fit(sequences, times, died)
    second = SurvivalPredictor(method="cox", random_state=2).fit(sequences, times, died)

    assert first.coefficients == second.coefficients


# --- Kaplan-Meier ---------------------------------------------------------------------------------------------


def test_kaplan_meier_matches_hand_computed_values():
    sequences = [["a"]] * 5
    model = SurvivalPredictor(method="baseline").fit(sequences, [1, 2, 3, 4, 5], [1, 0, 1, 1, 0])

    #   t = 1: 4 of 5 survive; t = 2 is censored; t = 3: 2 of 3; t = 4: 1 of 2; t = 5 is censored.
    grid = [0.5, 1, 1.5, 2, 3, 3.5, 4, 5, 9]
    expected = [1.0, 0.8, 0.8, 0.8, 0.8 * 2 / 3, 0.8 * 2 / 3, 0.8 * 2 / 3 * 0.5, 0.8 * 2 / 3 * 0.5, 0.8 * 2 / 3 * 0.5]
    np.testing.assert_allclose(model.predict_survival_function([["a"], []], grid), [expected, expected])
    # Area under the curve up to the last follow-up (5): 1*1 + 0.8*2 + (8/15)*1 + (4/15)*1 = 3.4.
    np.testing.assert_allclose(model.predict([["a"], []]), [3.4, 3.4])
    assert model.method == "baseline" and model.coefficients == {}


def test_kaplan_meier_matches_the_package_reference_implementation():
    sequences, times, died = simulate(300, seed=2, censoring_rate=0.08, rounding=True)
    reference = kaplan_meier_estimator(list(times), [int(flag) for flag in died])

    model = SurvivalPredictor(method="kaplan_meier").fit(sequences, times, died)

    grid = np.array(reference["time_points"])
    np.testing.assert_allclose(
        model.predict_survival_function([[]], grid)[0], reference["survival_probability"], rtol=0, atol=1e-12
    )


# --- no PyTorch ----------------------------------------------------------------------------------------------


@pytest.fixture
def torch_is_missing(monkeypatch: pytest.MonkeyPatch) -> None:
    """Make ``import torch`` fail (Python's documented ``sys.modules[name] = None``), installed or not."""
    monkeypatch.setitem(sys.modules, "torch", None)


@pytest.mark.parametrize("method", METHODS)
def test_numpy_methods_fit_and_predict_without_torch(torch_is_missing, method: str, cohort):
    sequences, times, died = cohort

    model = SurvivalPredictor(method=method).fit(sequences, times, died)

    assert model.is_fitted and model.model is not None
    predicted = model.predict(RISK_ORDERED_PROBES)
    assert np.isfinite(predicted).all() and (predicted > 0).all()
    assert predicted.max() < times.max() + 1e-9  # restricted mean survival time is bounded by the horizon
    if method != "baseline":
        assert model.coefficients["health:diagnosis"] > 0.5  # the fit really used the data
        assert predicted[0] < predicted[-1]


@pytest.mark.parametrize("method", METHODS)
def test_predict_before_fit_is_an_error_not_a_constant(method: str):
    model = SurvivalPredictor(method=method)

    with pytest.raises(ValueError, match="must be fitted"):
        model.predict([["health:diagnosis"]])
    with pytest.raises(ValueError, match="must be fitted"):
        model.predict_survival_function([["health:diagnosis"]], [1.0])
    with pytest.raises(ValueError, match="must be fitted"):
        model.coefficients
    assert model.is_fitted is False and model.model is None


def test_a_failed_fit_does_not_leave_a_model_that_claims_to_be_fitted():
    model = SurvivalPredictor(method="cox")

    with pytest.raises(ValueError, match="At least one event"):
        model.fit([["a"], ["b"]], [1.0, 2.0], [False, False])

    assert model.is_fitted is False
    assert model.model is None
    assert model.event_vocab is None
    with pytest.raises(ValueError, match="must be fitted"):
        model.predict([["a"]])


def test_a_failed_refit_keeps_the_previous_model(cohort):
    sequences, times, died = cohort
    model = SurvivalPredictor(method="cox").fit(sequences, times, died)
    before = model.coefficients

    with pytest.raises(ValueError):
        model.fit([["a"], ["b"]], [1.0, 2.0], [False, False])

    assert model.is_fitted and model.coefficients == before


@pytest.mark.parametrize("method", ["random_survival_forest", "deepsurv", "COX", ""])
def test_unknown_methods_are_rejected_instead_of_falling_back_to_kaplan_meier(method: str):
    with pytest.raises(ValueError, match="Unknown survival method"):
        SurvivalPredictor(method=method)


def test_fit_returns_the_model_for_chaining(cohort):
    sequences, times, died = cohort
    model = SurvivalPredictor(method="cox")

    assert model.fit(sequences, times, died) is model
    assert model.is_fitted and model.num_events == len(VOCABULARY)
    assert model.event_vocab == {token: index for index, token in enumerate(sorted(VOCABULARY))}


# --- input validation -----------------------------------------------------------------------------------------


@pytest.mark.parametrize(
    "sequences, times, died, message",
    [
        ([], [], [], "cannot be empty"),
        ([["a"], ["b"]], [1.0], [True, True], "same length"),
        ([["a"], ["b"]], [1.0, 2.0], [True], "same length"),
        ([["a"], ["b"]], [1.0, -2.0], [True, True], "not negative"),
        ([["a"], ["b"]], [1.0, np.nan], [True, True], "finite"),
        ([["a"], ["b"]], [1.0, np.inf], [True, True], "finite"),
        ([["a"], ["b"]], [1.0, 2.0], [2, 1], "booleans or 0 and 1"),
        ([["a"], ["b"]], [1.0, 2.0], [0.5, 1.0], "booleans or 0 and 1"),
        ([["a"], ["b"]], [1.0, 2.0], [False, False], "At least one event"),
        ([["a"], ["b"]], [[1.0, 2.0]], [[True, True]], "one-dimensional"),
    ],
)
def test_invalid_survival_data_is_rejected(sequences: List[Any], times: Any, died: Any, message: str):
    with pytest.raises(ValueError, match=message):
        SurvivalPredictor(method="cox").fit(sequences, times, died)


def test_sequences_without_any_event_cannot_be_used_for_covariate_models():
    with pytest.raises(ValueError, match="no events"):
        SurvivalPredictor(method="cox").fit([[], []], [1.0, 2.0], [True, True])
    # Kaplan-Meier needs no covariates.
    model = SurvivalPredictor(method="baseline").fit([[], []], [1.0, 2.0], [True, True])
    np.testing.assert_allclose(model.predict_survival_function([[]], [0.5, 1.0, 2.0])[0], [1.0, 0.5, 0.0])


@pytest.mark.parametrize("options, message", [({"penalizer": -1.0}, "penalizer"), ({"max_iter": 0}, "max_iter")])
def test_invalid_options_are_rejected(options: Dict[str, Any], message: str):
    with pytest.raises(ValueError, match=message):
        SurvivalPredictor(method="cox", **options)


def test_not_converging_in_the_allowed_iterations_warns(cohort):
    sequences, times, died = cohort

    with pytest.warns(RuntimeWarning, match="did not converge"):
        model = SurvivalPredictor(method="cox", max_iter=1, penalizer=0.0).fit(sequences, times, died)

    assert model.model is not None and model.model["converged"] is False


def test_integer_and_boolean_event_indicators_are_equivalent(cohort):
    sequences, times, died = cohort

    as_bool = SurvivalPredictor(method="cox").fit(sequences, times, died)
    as_int = SurvivalPredictor(method="cox").fit(sequences, times, died.astype(int))
    as_float = SurvivalPredictor(method="cox").fit(sequences, times, died.astype(float))

    assert as_bool.coefficients == as_int.coefficients == as_float.coefficients
