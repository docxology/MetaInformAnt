"""Argument-contract tests for ``cross_validate`` and ``learning_curve``.

``cross_validate`` declares ``X`` and ``y`` with ``None`` defaults (they trail
the optional ``model`` parameter) and ``learning_curve`` documents
``train_sizes=None`` as "use the default fractions". These tests pin both
contracts with real NumPy data and real functions; no test doubles.
"""

from __future__ import annotations

import numpy as np
import pytest

from metainformant.ml.evaluation.validation import cross_validate, learning_curve


def _two_class_data(n: int = 60, seed: int = 11) -> tuple[np.ndarray, np.ndarray]:
    rng = np.random.RandomState(seed)
    X = rng.randn(n, 4)
    y = (X[:, 0] + 0.1 * rng.randn(n) > 0).astype(int)
    return X, y


def _majority_class_scorer(X_train: np.ndarray, y_train: np.ndarray, X_val: np.ndarray, y_val: np.ndarray) -> dict:
    """Real classifier function: predict the training majority class."""
    majority = int(np.bincount(y_train).argmax())
    return {"accuracy": float(np.mean(y_val == majority))}


class TestCrossValidateClassifierFunctionArguments:
    def test_missing_x_and_y_raises_a_clear_value_error(self) -> None:
        with pytest.raises(ValueError, match="X and y"):
            cross_validate(classifier_func=_majority_class_scorer, cv=3)

    def test_missing_y_raises_a_clear_value_error(self) -> None:
        X, _y = _two_class_data()

        with pytest.raises(ValueError, match="X and y"):
            cross_validate(X=X, classifier_func=_majority_class_scorer, cv=3)

    def test_missing_x_raises_a_clear_value_error(self) -> None:
        _X, y = _two_class_data()

        with pytest.raises(ValueError, match="X and y"):
            cross_validate(y=y, classifier_func=_majority_class_scorer, cv=3)

    def test_complete_arguments_still_cross_validate(self) -> None:
        X, y = _two_class_data()

        result = cross_validate(X=X, y=y, classifier_func=_majority_class_scorer, cv=3, random_state=0)

        assert result["cv_folds"] == 3
        assert len(result["fold_results"]) == 3
        assert 0.0 <= result["mean_accuracy"] <= 1.0


class TestLearningCurveDefaultTrainSizes:
    def test_omitting_train_sizes_uses_ten_fractions_from_ten_to_one_hundred_percent(self) -> None:
        X, y = _two_class_data(n=80)

        result = learning_curve(X, y, classifier_func=_majority_class_scorer, cv=3, random_state=0)

        assert result["train_sizes"] == pytest.approx(np.linspace(0.1, 1.0, 10).tolist())
        assert len(result["train_scores"]) == 10
        assert len(result["validation_scores"]) == 10
        assert all(len(fold_scores) == 3 for fold_scores in result["validation_scores"])

    def test_explicit_train_sizes_are_respected(self) -> None:
        X, y = _two_class_data(n=80)

        result = learning_curve(
            X,
            y,
            classifier_func=_majority_class_scorer,
            train_sizes=np.array([0.5, 1.0]),
            cv=3,
            random_state=0,
        )

        assert result["train_sizes"] == [0.5, 1.0]
