"""Regression tests for ``compare_validation_strategies`` and ``cross_validate_biological(method="ensemble")``.

* ``compare_validation_strategies`` handed a ``StratifiedKFold`` object to ``cross_validation_scores(cv=...)``, whose
  ``cv`` is the *number* of folds (it builds ``StratifiedKFold(n_splits=cv)`` itself). The default
  ``"stratified_5fold"`` strategy therefore raised on every call and was reported as ``{"error": ...}``.
* ``cross_validate_biological(method="ensemble", random_state=<int>)`` put ``random_state`` into the options it
  forwards to ``train_ensemble_classifier`` *and* passed it as a named argument, so Python raised ``TypeError: got
  multiple values for keyword argument 'random_state'`` for every seeded ensemble run.

Oracles are scikit-learn's own ``cross_val_score`` over the same folds (``StratifiedKFold(5, shuffle=True,
random_state)``) and a hand-assembled soft-voting ensemble with the same members, so the expected numbers do not come
from the code under test. All models are real and seeded; ``n_jobs=1`` keeps the tests independent of joblib process
pools. No test doubles.
"""

from __future__ import annotations

import numpy as np
import pytest
from sklearn.datasets import make_classification
from sklearn.ensemble import GradientBoostingClassifier, RandomForestClassifier, VotingClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold, cross_val_score

from metainformant.ml.evaluation.validation import compare_validation_strategies
from metainformant.ml.models.classification import cross_validate_biological

X, Y = make_classification(n_samples=120, n_features=6, n_informative=4, n_redundant=0, random_state=0)


def _logistic() -> LogisticRegression:
    return LogisticRegression(max_iter=1000)


def test_stratified_5fold_strategy_scores_instead_of_erroring() -> None:
    result = compare_validation_strategies(X, Y, _logistic, strategies=["stratified_5fold"], random_state=0)

    entry = result["comparison"]["stratified_5fold"]
    assert "error" not in entry, entry.get("error")
    oracle = cross_val_score(
        _logistic(), X, Y, cv=StratifiedKFold(n_splits=5, shuffle=True, random_state=0), scoring="accuracy"
    ).mean()
    assert entry["score"] == pytest.approx(oracle)
    assert result["best_strategy"] == "stratified_5fold"


def test_default_strategies_all_succeed() -> None:
    result = compare_validation_strategies(X, Y, _logistic, random_state=0)

    assert result["strategies_tested"] == ["holdout", "5fold_cv", "10fold_cv", "stratified_5fold"]
    failed = {name: entry["error"] for name, entry in result["comparison"].items() if "error" in entry}
    assert failed == {}
    assert all(0.0 <= entry["score"] <= 1.0 for entry in result["comparison"].values())
    assert result["best_strategy"] in result["comparison"]


def test_stratified_5fold_matches_the_stratified_folds_of_5fold_cv() -> None:
    # for a classification target both strategies are shuffled stratified 5-fold CV with the same seed
    result = compare_validation_strategies(X, Y, _logistic, strategies=["5fold_cv", "stratified_5fold"], random_state=3)

    assert result["comparison"]["stratified_5fold"]["score"] == pytest.approx(result["comparison"]["5fold_cv"]["score"])


def test_unknown_strategy_is_still_reported_as_an_error() -> None:
    result = compare_validation_strategies(X, Y, _logistic, strategies=["bogus", "stratified_5fold"], random_state=0)

    assert "Unknown strategy" in result["comparison"]["bogus"]["error"]
    assert "score" in result["comparison"]["stratified_5fold"]


def _soft_voting_ensemble(seed: int | None, n_estimators: int = 10) -> VotingClassifier:
    """The ensemble ``train_ensemble_classifier`` builds, assembled by hand."""
    return VotingClassifier(
        estimators=[
            ("rf", RandomForestClassifier(n_estimators=n_estimators, random_state=seed)),
            ("gb", GradientBoostingClassifier(n_estimators=n_estimators, random_state=seed)),
            ("lr", LogisticRegression(random_state=seed, max_iter=1000)),
        ],
        voting="soft",
        n_jobs=1,
    )


@pytest.mark.parametrize("seed", [7, 0])
def test_seeded_ensemble_cross_validation_runs_and_matches_the_hand_built_ensemble(seed: int) -> None:
    result = cross_validate_biological(
        X, Y, method="ensemble", cv_folds=3, random_state=seed, parallel_jobs=1, n_jobs=1
    )

    assert result["method"] == "ensemble"
    oracle = cross_val_score(
        _soft_voting_ensemble(seed),
        X,
        Y,
        cv=StratifiedKFold(n_splits=3, shuffle=True, random_state=seed),
        scoring="accuracy",
        n_jobs=1,
    )
    assert result["cross_validation"]["accuracy"]["scores"] == pytest.approx(oracle.tolist())
    assert result["accuracy"] == pytest.approx(float(oracle.mean()))


def test_unseeded_ensemble_cross_validation_still_runs() -> None:
    result = cross_validate_biological(X, Y, method="ensemble", cv_folds=3, parallel_jobs=1, n_jobs=1)

    scores = np.asarray(result["cross_validation"]["accuracy"]["scores"])
    assert scores.shape == (3,) and np.all((0.0 <= scores) & (scores <= 1.0))


def test_seeded_ensemble_cross_validation_is_reproducible() -> None:
    first = cross_validate_biological(X, Y, method="ensemble", cv_folds=3, random_state=5, parallel_jobs=1, n_jobs=1)
    second = cross_validate_biological(X, Y, method="ensemble", cv_folds=3, random_state=5, parallel_jobs=1, n_jobs=1)

    assert first == second
