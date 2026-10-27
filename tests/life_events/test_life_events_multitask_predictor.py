"""Regression tests for ``MultiTaskPredictor`` (a real multi-task model instead of stored means).

``fit`` only stored the mean outcome of every task and ``predict`` returned that constant for every
sequence (or zeros when the model was not fitted, or for an unknown task name). The model now learns one
ridge-regression or logistic-regression head per task on a shared representation of the sequences, so
the predictions depend on the events of each sequence.

Heads are checked against independent optimisers (closed-form ridge, SciPy L-BFGS for the logistic loss)
and the tasks against data in which one event drives each task.
"""

from __future__ import annotations

import logging
from datetime import datetime
from typing import Any, Dict, List, Tuple

import numpy as np
import pytest
from scipy import optimize
from scipy.special import expit

from metainformant.life_events.core.events import Event, EventSequence
from metainformant.life_events.core.utils import generate_synthetic_life_events
from metainformant.life_events.models.statistical_models import MultiTaskPredictor

VOCABULARY = [
    "health:diagnosis",
    "occupation:job_change",
    "income:raise",
    "family:marriage",
    "education:degree",
    "housing:move",
]
INCOME_DRIVER = "occupation:job_change"  # every occurrence adds INCOME_SLOPE to the income
HEALTH_DRIVER = "health:diagnosis"  # its presence makes the health flag 1
INCOME_SLOPE = 2.0
TASKS = {"income": "regression", "health": "classification"}

Cohort = Tuple[List[List[str]], Dict[str, np.ndarray]]


def make_cohort(n: int, seed: int) -> Cohort:
    rng = np.random.default_rng(seed)
    sequences = [[str(token) for token in rng.choice(VOCABULARY, size=int(rng.integers(2, 9)))] for _ in range(n)]
    income = np.array([INCOME_SLOPE * s.count(INCOME_DRIVER) + 10 + rng.normal(0, 0.3) for s in sequences])
    health = np.array([int(HEALTH_DRIVER in s) for s in sequences])
    return sequences, {"income": income, "health": health}


def make_model(**options: Any) -> MultiTaskPredictor:
    settings: Dict[str, Any] = {"embedding_dim": 16, "random_state": 0}
    settings.update(options)
    return MultiTaskPredictor(task_types=dict(TASKS), **settings)


@pytest.fixture(scope="module")
def training() -> Cohort:
    return make_cohort(300, seed=0)


@pytest.fixture(scope="module")
def held_out() -> Cohort:
    return make_cohort(200, seed=1)


@pytest.fixture(scope="module")
def fitted(training: Cohort) -> MultiTaskPredictor:
    sequences, outcomes = training
    return make_model().fit(sequences, outcomes)


def as_event_sequences(token_lists: List[List[str]]) -> List[EventSequence]:
    sequences = []
    for index, tokens in enumerate(token_lists):
        events = []
        for position, token in enumerate(tokens):
            domain, event_type = token.split(":")
            events.append(Event(event_type, datetime(2000 + position, 1, 1), domain))
        sequences.append(EventSequence(f"person_{index}", events))
    return sequences


# --- the predictions depend on the sequences ----------------------------------------------------------------


def test_predictions_vary_with_the_input(fitted: MultiTaskPredictor, held_out: Cohort):
    sequences, _ = held_out

    predictions = fitted.predict(sequences)

    for task in TASKS:
        assert len(predictions[task]) == len(sequences)
        assert len(set(np.round(predictions[task], 6))) > 20, f"{task} predictions are (nearly) constant"


def test_each_task_tracks_its_own_driver(fitted: MultiTaskPredictor, held_out: Cohort):
    sequences, outcomes = held_out
    income = np.array(fitted.predict(sequences, task_name="income"))
    health = np.array(fitted.predict(sequences, task_name="health"))
    job_changes = np.array([s.count(INCOME_DRIVER) for s in sequences])
    diagnosed = np.array([int(HEALTH_DRIVER in s) for s in sequences])

    assert np.corrcoef(income, job_changes)[0, 1] > 0.95
    assert np.corrcoef(health, diagnosed)[0, 1] > 0.95
    # Each task follows its own driver far more than the other task's driver.
    assert abs(np.corrcoef(income, diagnosed)[0, 1]) < 0.4
    assert abs(np.corrcoef(health, job_changes)[0, 1]) < 0.4
    # Predictions are good, not just correlated: R^2 of the regression task, accuracy of the classification task.
    residual = np.sum((income - outcomes["income"]) ** 2)
    assert 1 - residual / np.sum((outcomes["income"] - outcomes["income"].mean()) ** 2) > 0.9
    assert np.mean((health > 0.5) == (outcomes["health"] == 1)) > 0.95


def test_the_event_based_outcome_of_generated_sequences_is_learned():
    """EventSequence objects from the package's generator, with the generator's own documented outcome rule."""
    sequences, outcomes = generate_synthetic_life_events(
        n_sequences=500, generate_outcomes=True, outcome_relationship="event_based", random_state=3
    )
    labels = np.array(outcomes)
    assert 0.3 < labels.mean() < 0.9  # both classes occur

    model = MultiTaskPredictor({"flag": "classification"}, embedding_dim=16, random_state=0)
    model.fit(sequences[:350], {"flag": labels[:350]})
    probabilities = np.array(model.predict(sequences[350:], task_name="flag"))

    assert model.num_events is not None and model.num_events > 20
    assert np.mean((probabilities > 0.5) == (labels[350:] == 1)) > 0.9


def test_adding_an_event_moves_only_the_task_it_drives(fitted: MultiTaskPredictor, held_out: Cohort):
    sequences, _ = held_out
    base = [s for s in sequences if HEALTH_DRIVER not in s][:50]

    def change(task: str, event: str) -> np.ndarray:
        after = np.array(fitted.predict([s + [event] for s in base], task_name=task))
        return after - np.array(fitted.predict(base, task_name=task))

    # One more job change adds the true slope to the income and leaves the health probability alone ...
    assert change("income", INCOME_DRIVER).mean() == pytest.approx(INCOME_SLOPE, abs=0.15)
    assert np.abs(change("health", INCOME_DRIVER)).max() < 0.15
    # ... a diagnosis flips the health probability and leaves the income alone ...
    assert change("health", HEALTH_DRIVER).min() > 0.8
    assert np.abs(change("income", HEALTH_DRIVER)).max() < 0.15
    # ... and an unrelated event moves neither.
    assert np.abs(change("income", "income:raise")).max() < 0.15
    assert np.abs(change("health", "income:raise")).max() < 0.15


def test_regression_predictions_are_on_the_scale_of_the_outcomes(fitted: MultiTaskPredictor, training: Cohort):
    sequences, outcomes = training

    income = np.array(fitted.predict(sequences, task_name="income"))

    assert income.mean() == pytest.approx(outcomes["income"].mean(), abs=1e-6)  # ridge keeps the intercept unpenalised
    assert income.min() > 5 and income.max() < 30


def test_classification_predictions_are_probabilities_of_the_larger_label(held_out: Cohort):
    sequences, outcomes = make_cohort(300, seed=3)
    relabelled = {"income": outcomes["income"], "health": np.where(outcomes["health"] == 1, 20, 10)}

    model = make_model().fit(sequences, relabelled)
    probabilities = np.array(model.predict(held_out[0], task_name="health"))

    diagnosed = np.array([HEALTH_DRIVER in s for s in held_out[0]])
    assert np.all((probabilities >= 0) & (probabilities <= 1))
    assert probabilities[diagnosed].min() > 0.5 > probabilities[~diagnosed].max()  # label 20 is the positive class


def test_string_and_boolean_labels_are_supported(training: Cohort):
    sequences, outcomes = training
    labels = np.where(outcomes["health"] == 1, "sick", "healthy")
    model = MultiTaskPredictor({"flag": "classification"}, embedding_dim=8, random_state=0)
    model.fit(sequences, {"flag": labels})
    from_strings = np.array(model.predict(sequences, task_name="flag"))

    boolean = MultiTaskPredictor({"flag": "classification"}, embedding_dim=8, random_state=0)
    boolean.fit(sequences, {"flag": outcomes["health"] == 1})

    # "sick" is the larger of the two labels, as True is, so both are the probability of the same class.
    np.testing.assert_allclose(from_strings, np.array(boolean.predict(sequences, task_name="flag")), atol=1e-12)
    assert np.mean((from_strings > 0.5) == (outcomes["health"] == 1)) > 0.95


# --- the heads against independent optimisers --------------------------------------------------------------


def test_regression_head_is_the_closed_form_ridge_solution(training: Cohort):
    sequences, outcomes = training
    model = MultiTaskPredictor({"income": "regression"}, embedding_dim=16, random_state=0, alpha=3.0)
    model.fit(sequences, {"income": outcomes["income"]})

    features = model.shared_representation(sequences)
    centred = features - features.mean(axis=0)
    weights = np.linalg.solve(
        centred.T @ centred + 3.0 * np.eye(features.shape[1]),
        centred.T @ (outcomes["income"] - outcomes["income"].mean()),
    )
    expected = outcomes["income"].mean() + centred @ weights

    np.testing.assert_allclose(model.predict(sequences, task_name="income"), expected, rtol=1e-9, atol=1e-9)


def test_classification_head_is_the_optimum_of_the_penalised_logistic_loss(training: Cohort):
    sequences, outcomes = training
    # A noisy label (the diagnosis flag with 15% of the labels flipped) keeps the classes from being separable.
    rng = np.random.default_rng(4)
    labels = np.where(rng.random(len(sequences)) < 0.15, 1 - outcomes["health"], outcomes["health"])
    model = MultiTaskPredictor({"flag": "classification"}, embedding_dim=16, random_state=0, alpha=2.0)
    model.fit(sequences, {"flag": labels})
    features = model.shared_representation(sequences)

    def loss(params: np.ndarray) -> Tuple[float, np.ndarray]:
        """Penalised negative log-likelihood (0.5 * alpha * ||weights||^2, alpha = 2) and its gradient."""
        linear = features @ params[:-1] + params[-1]
        value = np.sum(np.logaddexp(0.0, linear) - labels * linear) + 1.0 * np.sum(params[:-1] ** 2)
        residual = expit(linear) - labels
        return float(value), np.append(features.T @ residual + 2.0 * params[:-1], residual.sum())

    head = model.models["flag"]
    fitted_params = np.append(head["weights"], head["bias"])
    # First-order optimality: the objective is strictly convex, so a vanishing gradient means this is the optimum.
    assert np.abs(loss(fitted_params)[1]).max() < 1e-5

    reference = optimize.minimize(
        loss,
        np.zeros(len(fitted_params)),
        jac=True,
        method="L-BFGS-B",
        options={"gtol": 1e-12, "ftol": 1e-16, "maxiter": 5000},
    )
    expected = expit(features @ reference.x[:-1] + reference.x[-1])
    np.testing.assert_allclose(model.predict(sequences, task_name="flag"), expected, atol=1e-5)
    assert 0.5 < np.mean((expected > 0.5) == (labels == 1)) < 0.95  # the flipped labels cannot all be fitted


def test_perfectly_separable_classes_give_finite_probabilities(training: Cohort):
    sequences, outcomes = training  # the health flag is a function of the sequence: perfectly separable

    probabilities = np.array(make_model().fit(sequences, outcomes).predict(sequences, task_name="health"))

    assert np.isfinite(probabilities).all()
    assert 0 < probabilities.min() and probabilities.max() < 1  # L2 keeps them away from 0 and 1


def test_alpha_shrinks_the_predictions_toward_the_mean(training: Cohort):
    sequences, outcomes = training
    spread = {}
    for alpha in (0.1, 100.0, 1e6):
        model = MultiTaskPredictor({"income": "regression"}, embedding_dim=16, random_state=0, alpha=alpha)
        model.fit(sequences, {"income": outcomes["income"]})
        spread[alpha] = np.std(model.predict(sequences, task_name="income"))

    assert spread[0.1] > spread[100.0] > spread[1e6]
    assert spread[1e6] < 0.01 * spread[0.1]


# --- seeds, embeddings and the hidden dimension ----------------------------------------------------------------


def test_same_seed_gives_identical_predictions(training: Cohort, held_out: Cohort):
    sequences, outcomes = training

    first = make_model(random_state=7).fit(sequences, outcomes).predict(held_out[0])
    second = make_model(random_state=7).fit(sequences, outcomes).predict(held_out[0])

    assert first == second


def test_the_seed_drives_the_embeddings_and_through_them_the_predictions(training: Cohort, held_out: Cohort):
    sequences, outcomes = training

    first = make_model(random_state=1).fit(sequences, outcomes)
    second = make_model(random_state=2).fit(sequences, outcomes)

    assert not np.array_equal(first.event_embeddings[HEALTH_DRIVER], second.event_embeddings[HEALTH_DRIVER])
    assert first.predict(held_out[0]) != second.predict(held_out[0])
    for model in (first, second):  # a different seed changes the details, not the quality
        assert (
            np.corrcoef(model.predict(held_out[0], "income"), [s.count(INCOME_DRIVER) for s in held_out[0]])[0, 1]
            > 0.95
        )


def test_pretrained_event_embeddings_replace_the_learned_ones(training: Cohort, held_out: Cohort):
    sequences, outcomes = training
    rng = np.random.default_rng(0)
    vectors = {token: rng.normal(size=5) for token in VOCABULARY[:-1]}  # the last event has no embedding

    first = make_model(random_state=1).fit(sequences, outcomes, event_embeddings=vectors)
    second = make_model(random_state=2).fit(sequences, outcomes, event_embeddings=vectors)

    assert set(first.event_embeddings) == set(vectors)
    assert first.predict(held_out[0]) == second.predict(held_out[0])  # no training, so no seed dependence
    assert np.corrcoef(first.predict(held_out[0], "income"), [s.count(INCOME_DRIVER) for s in held_out[0]])[0, 1] > 0.95


def test_embedding_dim_is_the_size_of_the_learned_event_embeddings(training: Cohort):
    sequences, outcomes = training

    for dimension in (4, 24):
        model = make_model(embedding_dim=dimension).fit(sequences, outcomes)
        assert set(model.event_embeddings) == set(VOCABULARY)
        assert all(vector.shape == (dimension,) for vector in model.event_embeddings.values())


def test_hidden_dim_bounds_the_embedding_part_of_the_shared_representation(training: Cohort, held_out: Cohort):
    sequences, outcomes = training
    n_events = len(VOCABULARY)

    widths = {}
    for hidden_dim in (1, 2, 1000):
        model = make_model(hidden_dim=hidden_dim).fit(sequences, outcomes)
        widths[hidden_dim] = model.shared_representation(held_out[0]).shape
    # The event counts are always part of the representation; hidden_dim caps the embedding components.
    assert widths[1] == (len(held_out[0]), n_events + 1)
    assert widths[2] == (len(held_out[0]), n_events + 2)
    assert n_events + 2 < widths[1000][1] <= n_events + 16


def test_a_tiny_hidden_dim_does_not_lose_the_drivers(training: Cohort, held_out: Cohort):
    sequences, outcomes = training

    model = make_model(hidden_dim=1).fit(sequences, outcomes)

    income = model.predict(held_out[0], task_name="income")
    assert np.corrcoef(income, [s.count(INCOME_DRIVER) for s in held_out[0]])[0, 1] > 0.95


# --- interface and state --------------------------------------------------------------------------------------


def test_predict_returns_the_documented_shapes(fitted: MultiTaskPredictor, held_out: Cohort):
    sequences, _ = held_out

    everything = fitted.predict(sequences)
    one_task = fitted.predict(sequences, task_name="income")

    assert isinstance(everything, dict) and set(everything) == set(TASKS)
    assert isinstance(one_task, list) and one_task == everything["income"]
    assert all(isinstance(value, float) for value in one_task)
    assert fitted.predict([]) == {"income": [], "health": []}
    assert fitted.predict([], task_name="health") == []


def test_fitted_state_is_exposed(fitted: MultiTaskPredictor):
    assert fitted.is_fitted
    assert fitted.num_events == len(VOCABULARY)
    assert fitted.event_vocab == {token: index for index, token in enumerate(sorted(VOCABULARY))}
    assert {task: head["type"] for task, head in fitted.models.items()} == {"income": "ridge", "health": "logistic"}
    assert "_baseline" not in fitted.models  # the stored-means stand-in is gone


def test_predict_before_fit_is_an_error_not_zeros():
    model = make_model()

    with pytest.raises(ValueError, match="must be fitted"):
        model.predict([["health:diagnosis"]])
    with pytest.raises(ValueError, match="must be fitted"):
        model.predict([["health:diagnosis"]], task_name="income")
    with pytest.raises(ValueError, match="must be fitted"):
        model.shared_representation([["health:diagnosis"]])
    assert model.is_fitted is False


def test_an_unknown_task_name_is_an_error_not_zeros(fitted: MultiTaskPredictor):
    with pytest.raises(ValueError, match="Unknown task 'wealth'"):
        fitted.predict([["health:diagnosis"]], task_name="wealth")


def test_events_unseen_in_training_are_ignored(fitted: MultiTaskPredictor):
    known = fitted.predict([[HEALTH_DRIVER, INCOME_DRIVER], [INCOME_DRIVER], []])
    with_unseen = fitted.predict(
        [[HEALTH_DRIVER, "never:seen", INCOME_DRIVER], [INCOME_DRIVER, "also:unseen", "and:another"], ["only:unseen"]]
    )

    # An event the model has never seen has no learned effect and does not dilute the mean embedding either.
    for task in TASKS:
        np.testing.assert_allclose(known[task], with_unseen[task], rtol=1e-12, atol=1e-12)
    assert with_unseen["health"][0] > 0.5 > with_unseen["health"][1]


def test_event_sequences_and_token_lists_are_equivalent(training: Cohort, held_out: Cohort):
    sequences, outcomes = training
    from_tokens = make_model().fit(sequences, outcomes)
    from_events = make_model().fit(as_event_sequences(sequences), outcomes)

    assert from_tokens.predict(held_out[0]) == from_events.predict(as_event_sequences(held_out[0]))


def test_fit_returns_the_model_for_chaining(training: Cohort):
    sequences, outcomes = training
    model = make_model()

    assert model.fit(sequences, outcomes) is model


def test_a_failed_fit_leaves_the_model_unfitted_and_a_failed_refit_keeps_the_old_one(training: Cohort):
    sequences, outcomes = training
    model = make_model()

    with pytest.raises(ValueError):
        model.fit(sequences, {"income": outcomes["income"]})  # the health task has no outcomes
    assert model.is_fitted is False and model.models == {} and model.event_vocab is None

    model.fit(sequences, outcomes)
    before = model.predict(sequences[:20])
    with pytest.raises(ValueError):
        model.fit(sequences, {**outcomes, "health": np.zeros(len(sequences))})  # a single class
    assert model.is_fitted and model.predict(sequences[:20]) == before


def test_unsupported_constructor_options_are_reported(caplog: pytest.LogCaptureFixture):
    with caplog.at_level(logging.WARNING):
        MultiTaskPredictor(TASKS, epochs=50, learning_rate=0.1)

    assert "ignores unsupported options" in caplog.text
    assert "epochs" in caplog.text and "learning_rate" in caplog.text


def test_construction_does_not_seed_the_global_random_generators():
    np.random.seed(5)
    expected = np.random.random()
    np.random.seed(5)

    MultiTaskPredictor(TASKS, random_state=123)

    assert np.random.random() == expected


# --- validation -----------------------------------------------------------------------------------------------


@pytest.mark.parametrize(
    "arguments, message",
    [
        ({"task_types": {}}, "must not be empty"),
        ({"task_types": {"a": "ranking"}}, "Unknown task types"),
        ({"task_types": TASKS, "embedding_dim": 0}, "must be positive"),
        ({"task_types": TASKS, "hidden_dim": 0}, "must be positive"),
        ({"task_types": TASKS, "alpha": 0.0}, "alpha must be positive"),
        ({"task_types": TASKS, "alpha": -1.0}, "alpha must be positive"),
    ],
)
def test_invalid_constructor_arguments_are_rejected(arguments: Dict[str, Any], message: str):
    with pytest.raises(ValueError, match=message):
        MultiTaskPredictor(**arguments)


def test_invalid_training_data_is_rejected(training: Cohort):
    sequences, outcomes = training
    short = {task: values[:10] for task, values in outcomes.items()}
    with_nan = {**outcomes, "income": np.where(np.arange(len(sequences)) == 3, np.nan, outcomes["income"])}
    three_classes = {**outcomes, "health": np.arange(len(sequences)) % 3}
    text_regression = {**outcomes, "income": np.array(["high"] * len(sequences))}
    two_dimensional = {**outcomes, "income": outcomes["income"].reshape(-1, 1)}

    cases: List[Tuple[Any, Dict[str, np.ndarray], str]] = [
        ([], {"income": [], "health": []}, "sequences cannot be empty"),
        (sequences, {"income": outcomes["income"]}, "outcomes are missing for the tasks"),
        (sequences, {**outcomes, "wealth": outcomes["income"]}, "not in task_types"),
        (sequences, short, "one value per sequence"),
        (sequences, with_nan, "must be finite"),
        (sequences, three_classes, "exactly two distinct labels"),
        (sequences, {**outcomes, "health": np.zeros(len(sequences))}, "exactly two distinct labels"),
        (sequences, text_regression, "must be numeric"),
        (sequences, two_dimensional, "one value per sequence"),
        ([[], []], {"income": [1.0, 2.0], "health": [0, 1]}, "no events"),
    ]
    for case_sequences, case_outcomes, message in cases:
        with pytest.raises(ValueError, match=message):
            make_model().fit(case_sequences, case_outcomes)


def test_invalid_pretrained_embeddings_are_rejected(training: Cohort):
    sequences, outcomes = training

    with pytest.raises(ValueError, match="event_embeddings"):
        make_model().fit(sequences, outcomes, event_embeddings={})
    with pytest.raises(ValueError, match="event_embeddings"):
        make_model().fit(sequences, outcomes, event_embeddings={VOCABULARY[0]: np.ones(3), VOCABULARY[1]: np.ones(4)})
    with pytest.raises(ValueError, match="event_embeddings"):
        make_model().fit(sequences, outcomes, event_embeddings={VOCABULARY[0]: np.ones((2, 2))})
