"""Regression tests: ML analysis routines take a seed and neither seed nor consume the global random state.

The routines below drew their randomness from the global generators. Several seeded those generators as a side effect
(``random.seed(random_state)`` / ``np.random.seed(random_state)``, and ``np.random.seed(None)``, which *re*seeds from
system entropy and throws away the caller's seed), and the Kernel SHAP and LIME explainers had no seed argument at
all:

* ``evaluation.validation``: ``permutation_importance_biological``, ``validate_model_stability``,
  ``bootstrap_validate``, ``learning_curve``
* ``models.regression``: ``analyze_prediction_uncertainty``
* ``automl.optimization``: ``random_search``, ``bayesian_optimization``
* ``interpretability.explainers``: ``compute_permutation_importance``, ``compute_shap_values_kernel``,
  ``compute_lime_explanation``
* ``interpretability.feature_selection``: ``boruta_selection`` (which also passed ``None`` to its random forest when
  ``random_state`` was ``0``, so that seed was not reproducible), ``stability_selection``

Every routine now draws from a generator of its own: the same seed reproduces the output, different seeds differ, and
``random.getstate()`` / ``numpy.random.get_state()`` are identical before and after a call. Oracles are exact where
the maths allows (cross-validated accuracy from scikit-learn, exact Shapley values of a linear model) and structural
otherwise. All models are real; the ones used here are deterministic given their inputs. No test doubles.
"""

from __future__ import annotations

import random

import numpy as np
import pytest
from sklearn.datasets import make_classification, make_regression
from sklearn.linear_model import LinearRegression, LogisticRegression
from sklearn.model_selection import StratifiedKFold, cross_val_score

from metainformant.ml import _numeric
from metainformant.ml.automl.optimization import bayesian_optimization, random_search
from metainformant.ml.evaluation.validation import (
    bootstrap_validate,
    learning_curve,
    permutation_importance_biological,
    validate_model_stability,
)
from metainformant.ml.interpretability.explainers import (
    compute_lime_explanation,
    compute_permutation_importance,
    compute_shap_values_kernel,
)
from metainformant.ml.interpretability.feature_selection import boruta_selection, stability_selection
from metainformant.ml.models.regression import analyze_prediction_uncertainty

X_CLS, Y_CLS = make_classification(
    n_samples=80, n_features=5, n_informative=3, n_redundant=0, class_sep=1.5, random_state=0
)
X_REG, Y_REG = make_regression(n_samples=60, n_features=4, noise=5.0, random_state=0)


class GlobalStateGuard:
    """Seeds both global generators to known states and checks that a block of code left them exactly there."""

    def __enter__(self) -> "GlobalStateGuard":
        random.seed(20261027)
        np.random.seed(20261027)
        self.python_state = random.getstate()
        self.numpy_state = np.random.get_state()
        return self

    def __exit__(self, *exc_info) -> None:
        if exc_info[0] is not None:
            return
        assert random.getstate() == self.python_state, "the global `random` state changed"
        after = np.random.get_state()
        assert after[0] == self.numpy_state[0] and after[2:] == self.numpy_state[2:]
        assert np.array_equal(after[1], self.numpy_state[1]), "the global NumPy random state changed"


def _linear() -> LinearRegression:
    return LinearRegression()


# ---------------------------------------------------------------------------
# evaluation.validation
# ---------------------------------------------------------------------------


def test_permutation_importance_biological_is_reproducible_and_leaves_global_state_alone() -> None:
    model = LogisticRegression(max_iter=500)
    with GlobalStateGuard():
        first = permutation_importance_biological(model, X_CLS, Y_CLS, n_repeats=2, random_state=3)
    second = permutation_importance_biological(model, X_CLS, Y_CLS, n_repeats=2, random_state=3)
    other = permutation_importance_biological(model, X_CLS, Y_CLS, n_repeats=2, random_state=4)

    assert first == second
    assert first["raw_scores"] != other["raw_scores"]


def test_permutation_importance_biological_first_draw_matches_a_hand_computed_score_drop() -> None:
    model = LogisticRegression(max_iter=500)
    folds = StratifiedKFold(n_splits=5, shuffle=True, random_state=3)
    baseline = cross_val_score(model, X_CLS, Y_CLS, cv=folds, scoring="accuracy").mean()
    permuted = X_CLS.copy()
    permuted[:, 0] = X_CLS[np.random.RandomState(3).permutation(len(X_CLS)), 0]  # the first permutation drawn
    expected_drop = baseline - cross_val_score(model, permuted, Y_CLS, cv=folds, scoring="accuracy").mean()

    result = permutation_importance_biological(model, X_CLS, Y_CLS, n_repeats=2, random_state=3)

    assert result["baseline_score"] == pytest.approx(float(baseline))
    assert result["raw_scores"][0][0] == pytest.approx(float(expected_drop))


def test_validate_model_stability_is_reproducible_and_leaves_global_state_alone() -> None:
    with GlobalStateGuard():
        first = validate_model_stability(
            X_CLS, Y_CLS, lambda: LogisticRegression(max_iter=500), n_bootstraps=6, random_state=4
        )
    second = validate_model_stability(
        X_CLS, Y_CLS, lambda: LogisticRegression(max_iter=500), n_bootstraps=6, random_state=4
    )
    other = validate_model_stability(
        X_CLS, Y_CLS, lambda: LogisticRegression(max_iter=500), n_bootstraps=6, random_state=5
    )

    assert first == second
    assert first["scores"]["values"] != other["scores"]["values"]


def _predictions(X_train, y_train, X_test, y_test):
    return LinearRegression().fit(X_train, y_train).predict(X_test)


def test_bootstrap_validate_is_reproducible_and_leaves_global_state_alone() -> None:
    with GlobalStateGuard():
        first = bootstrap_validate(X_REG, Y_REG, _predictions, n_bootstrap=5, random_state=5)
    second = bootstrap_validate(X_REG, Y_REG, _predictions, n_bootstrap=5, random_state=5)
    other = bootstrap_validate(X_REG, Y_REG, _predictions, n_bootstrap=5, random_state=6)

    assert first["mean_mse"] == second["mean_mse"]
    assert first["mean_mse"] != other["mean_mse"]


def test_learning_curve_is_reproducible_and_leaves_global_state_alone() -> None:
    sizes = np.array([0.3, 0.6, 1.0])
    with GlobalStateGuard():
        first = learning_curve(
            X_CLS, Y_CLS, lambda: LogisticRegression(max_iter=500), train_sizes=sizes, cv=3, random_state=6
        )
    second = learning_curve(
        X_CLS, Y_CLS, lambda: LogisticRegression(max_iter=500), train_sizes=sizes, cv=3, random_state=6
    )
    other = learning_curve(
        X_CLS, Y_CLS, lambda: LogisticRegression(max_iter=500), train_sizes=sizes, cv=3, random_state=7
    )

    assert first == second
    assert first["train_scores"] != other["train_scores"] or first["validation_scores"] != other["validation_scores"]


# ---------------------------------------------------------------------------
# models.regression
# ---------------------------------------------------------------------------


def test_prediction_uncertainty_is_reproducible_and_leaves_global_state_alone() -> None:
    X_test = X_REG[:5]
    with GlobalStateGuard():
        first = analyze_prediction_uncertainty(X_REG, Y_REG, X_test, method="ridge", n_bootstraps=8, random_state=7)
    second = analyze_prediction_uncertainty(X_REG, Y_REG, X_test, method="ridge", n_bootstraps=8, random_state=7)
    other = analyze_prediction_uncertainty(X_REG, Y_REG, X_test, method="ridge", n_bootstraps=8, random_state=8)

    assert first == second
    assert first["std_predictions"] != other["std_predictions"]


# ---------------------------------------------------------------------------
# automl.optimization
# ---------------------------------------------------------------------------

C_RANGE = {"C": {"low": 0.01, "high": 10.0, "log": True}}


def _logistic(C: float) -> LogisticRegression:
    return LogisticRegression(C=C, max_iter=500)


def test_random_search_is_reproducible_and_leaves_global_state_alone() -> None:
    with GlobalStateGuard():
        first = random_search(_logistic, C_RANGE, X_CLS, Y_CLS, n_iter=5, cv=3, random_state=11)
    second = random_search(_logistic, C_RANGE, X_CLS, Y_CLS, n_iter=5, cv=3, random_state=11)
    other = random_search(_logistic, C_RANGE, X_CLS, Y_CLS, n_iter=5, cv=3, random_state=12)

    assert first == second
    sampled = [entry["params"]["C"] for entry in first["all_results"]]
    assert all(0.01 <= value <= 10.0 for value in sampled)
    assert sampled != [entry["params"]["C"] for entry in other["all_results"]]


def test_random_search_without_a_seed_does_not_touch_the_global_state_either() -> None:
    with GlobalStateGuard():
        random_search(_logistic, C_RANGE, X_CLS, Y_CLS, n_iter=3, cv=3)


def _bowl(params: dict) -> float:
    return -((params["x"] - 0.3) ** 2) - (params["y"] - 0.7) ** 2


SPACE = {"x": {"low": 0.0, "high": 1.0}, "y": {"low": 0.0, "high": 1.0}}


def test_bayesian_optimization_is_reproducible_and_leaves_global_state_alone() -> None:
    with GlobalStateGuard():
        first = bayesian_optimization(_bowl, SPACE, n_iter=10, n_initial=4, random_state=13)
    second = bayesian_optimization(_bowl, SPACE, n_iter=10, n_initial=4, random_state=13)
    other = bayesian_optimization(_bowl, SPACE, n_iter=10, n_initial=4, random_state=14)

    assert first == second
    assert first["history"] != other["history"]
    assert first["best_score"] > -0.15, "ten evaluations of a smooth bowl should land close to its optimum"


def test_bayesian_optimization_without_a_seed_does_not_touch_the_global_state_either() -> None:
    with GlobalStateGuard():
        bayesian_optimization(_bowl, SPACE, n_iter=6, n_initial=3)


# ---------------------------------------------------------------------------
# interpretability.explainers
# ---------------------------------------------------------------------------


def test_permutation_importance_is_reproducible_ranks_the_signal_and_leaves_global_state_alone() -> None:
    model = LogisticRegression(max_iter=500).fit(X_CLS, Y_CLS)
    with GlobalStateGuard():
        first = compute_permutation_importance(model, X_CLS, Y_CLS, n_repeats=5, random_state=8)
    second = compute_permutation_importance(model, X_CLS, Y_CLS, n_repeats=5, random_state=8)
    other = compute_permutation_importance(model, X_CLS, Y_CLS, n_repeats=5, random_state=9)

    assert first == second
    assert first["importances_mean"] != other["importances_mean"]
    # the largest absolute coefficient belongs to the feature the model relies on most
    assert int(np.argmax(first["importances_mean"])) == int(np.argmax(np.abs(model.coef_[0])))


def test_kernel_shap_recovers_exact_shapley_values_of_a_linear_model_and_is_seedable() -> None:
    model = LinearRegression().fit(X_REG, Y_REG)
    instances = X_REG[:3]
    with GlobalStateGuard():
        first = compute_shap_values_kernel(model.predict, instances, n_samples=400, random_state=9)
    second = compute_shap_values_kernel(model.predict, instances, n_samples=400, random_state=9)
    other = compute_shap_values_kernel(model.predict, instances, n_samples=400, random_state=10)

    # the default background is the mean of the explained instances; for a linear model phi_j = coef_j * (x_j - mean_j)
    exact = model.coef_ * (instances - instances.mean(axis=0))
    assert np.array(first["shap_values"]) == pytest.approx(exact, rel=1e-5, abs=1e-5)
    assert first == second
    assert np.array(other["shap_values"]) == pytest.approx(exact, rel=1e-5, abs=1e-5)


def test_lime_explanation_is_reproducible_and_leaves_global_state_alone() -> None:
    model = LinearRegression().fit(X_REG, Y_REG)
    instance = [float(v) for v in X_REG[0]]
    names = [f"f{i}" for i in range(4)]
    with GlobalStateGuard():
        first = compute_lime_explanation(model.predict, instance, names, n_samples=300, random_state=10)
    second = compute_lime_explanation(model.predict, instance, names, n_samples=300, random_state=10)
    other = compute_lime_explanation(model.predict, instance, names, n_samples=300, random_state=11)

    assert first == second
    assert first["coefficients"] != other["coefficients"]


# ---------------------------------------------------------------------------
# interpretability.feature_selection
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("seed", [0, 21])
def test_boruta_is_reproducible_for_every_seed_including_zero(seed: int) -> None:
    with GlobalStateGuard():
        first = boruta_selection(X_CLS, Y_CLS, max_iter=4, random_state=seed)
    second = boruta_selection(X_CLS, Y_CLS, max_iter=4, random_state=seed)

    assert first == second
    assert first["importance_history"] == second["importance_history"]


def test_boruta_seeds_differ() -> None:
    first = boruta_selection(X_CLS, Y_CLS, max_iter=4, random_state=1)
    other = boruta_selection(X_CLS, Y_CLS, max_iter=4, random_state=2)

    assert first["importance_history"] != other["importance_history"]


def test_stability_selection_is_reproducible_and_leaves_global_state_alone() -> None:
    with GlobalStateGuard():
        first = stability_selection(X_REG, Y_REG, n_bootstrap=12, alpha=0.5, random_state=1)
    second = stability_selection(X_REG, Y_REG, n_bootstrap=12, alpha=0.5, random_state=1)
    other = stability_selection(X_REG, Y_REG, n_bootstrap=12, alpha=0.5, random_state=2)

    assert first == second
    assert first["selection_probabilities"] != other["selection_probabilities"] or first["selected_features"]


# ---------------------------------------------------------------------------
# the shared generator helper
# ---------------------------------------------------------------------------


def test_resolve_python_rng_accepts_every_documented_kind_of_seed() -> None:
    resolve = _numeric.resolve_python_rng

    assert resolve(5).random() == random.Random(5).random()
    assert resolve(np.int64(5)).random() == random.Random(5).random()
    own = random.Random(7)
    assert resolve(own) is own
    assert resolve(np.random.default_rng(3)).random() == resolve(np.random.default_rng(3)).random()
    assert resolve(None).random() != resolve(None).random()
