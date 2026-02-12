"""Tests for ML evaluation and validation utilities.

NO MOCKING POLICY: All tests use real implementations.
"""
from __future__ import annotations

import numpy as np
import pytest

from metainformant.ml.evaluation.validation import (
    biological_data_validator,
    bootstrap_validate,
    cross_validate,
    cross_validation_scores,
    k_fold_split,
    learning_curve,
    train_test_split_biological,
    validate_model_stability,
)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _make_classification_data(n=100, n_features=5, n_classes=2, seed=42):
    rng = np.random.RandomState(seed)
    X = rng.randn(n, n_features)
    y = rng.randint(0, n_classes, size=n)
    return X, y


def _make_regression_data(n=100, n_features=5, seed=42):
    rng = np.random.RandomState(seed)
    X = rng.randn(n, n_features)
    y = X @ rng.randn(n_features) + rng.randn(n) * 0.1
    return X, y


# ---------------------------------------------------------------------------
# train_test_split_biological
# ---------------------------------------------------------------------------


class TestTrainTestSplitBiological:
    def test_basic_split(self):
        X, y = _make_classification_data()
        X_train, X_test, y_train, y_test = train_test_split_biological(X, y, test_size=0.2, random_state=42)
        assert len(X_train) == 80
        assert len(X_test) == 20
        assert len(y_train) == 80
        assert len(y_test) == 20

    def test_auto_stratification(self):
        X, y = _make_classification_data(n=200, n_classes=2)
        X_train, X_test, y_train, y_test = train_test_split_biological(X, y, test_size=0.25, random_state=42)
        # Both classes should appear in both splits
        assert len(np.unique(y_train)) == 2
        assert len(np.unique(y_test)) == 2

    def test_custom_test_size(self):
        X, y = _make_classification_data(n=50)
        X_train, X_test, y_train, y_test = train_test_split_biological(X, y, test_size=0.4, random_state=42)
        assert len(X_test) == 20


# ---------------------------------------------------------------------------
# cross_validation_scores
# ---------------------------------------------------------------------------


class TestCrossValidationScores:
    def test_single_metric(self):
        from sklearn.tree import DecisionTreeClassifier

        X, y = _make_classification_data(n=100, n_features=4)
        model = DecisionTreeClassifier(random_state=42)
        result = cross_validation_scores(model, X, y, cv=3, scoring="accuracy", random_state=42)
        assert "accuracy" in result
        assert len(result["accuracy"]) == 3

    def test_multiple_metrics(self):
        from sklearn.tree import DecisionTreeClassifier

        X, y = _make_classification_data(n=100, n_features=4)
        model = DecisionTreeClassifier(random_state=42)
        result = cross_validation_scores(model, X, y, cv=3, scoring=["accuracy", "f1"], random_state=42)
        assert "accuracy" in result
        assert "f1" in result


# ---------------------------------------------------------------------------
# biological_data_validator
# ---------------------------------------------------------------------------


class TestBiologicalDataValidator:
    def test_clean_data_passes(self):
        X, y = _make_classification_data()
        result = biological_data_validator(X, y)
        assert result["passed"] is True
        assert result["n_samples"] == 100

    def test_missing_values_detected(self):
        X, y = _make_classification_data()
        X[0, 0] = np.nan
        result = biological_data_validator(X, y)
        assert result["checks"]["missing_values"]["passed"] == False
        assert result["passed"] == False

    def test_infinite_values_detected(self):
        X, y = _make_classification_data()
        X[5, 2] = np.inf
        result = biological_data_validator(X, y)
        assert result["checks"]["infinite_values"]["passed"] == False

    def test_constant_feature_detected(self):
        X, y = _make_classification_data()
        X[:, 0] = 1.0  # constant column
        result = biological_data_validator(X, y)
        assert result["checks"]["constant_features"]["passed"] is False
        assert 0 in result["checks"]["constant_features"]["constant_feature_indices"]

    def test_custom_checks(self):
        X, y = _make_classification_data()
        result = biological_data_validator(X, y, checks=["missing_values"])
        assert "missing_values" in result["checks"]
        assert "infinite_values" not in result["checks"]


# ---------------------------------------------------------------------------
# bootstrap_validate
# ---------------------------------------------------------------------------


class TestBootstrapValidate:
    def test_basic_bootstrap(self):
        X, y = _make_regression_data(n=50, n_features=3)

        def model_fn(X_train, y_train, X_test, y_test):
            from sklearn.linear_model import LinearRegression

            m = LinearRegression()
            m.fit(X_train, y_train)
            preds = m.predict(X_test)
            return preds

        result = bootstrap_validate(X, y, model_fn, n_bootstrap=5, random_state=42)
        assert "mean_mse" in result
        assert result["n_bootstrap"] == 5
        assert result["mean_mse"] >= 0


# ---------------------------------------------------------------------------
# cross_validate
# ---------------------------------------------------------------------------


class TestCrossValidate:
    def test_with_sklearn_model(self):
        from sklearn.tree import DecisionTreeClassifier

        X, y = _make_classification_data(n=100)
        result = cross_validate(model=DecisionTreeClassifier(random_state=42), X=X, y=y, cv=3)
        assert "mean_score" in result
        assert result["cv_folds"] == 3

    def test_with_classifier_func(self):
        from sklearn.linear_model import LogisticRegression

        X, y = _make_classification_data(n=100)

        def clf_fn(X_train, y_train, X_val, y_val):
            m = LogisticRegression(random_state=42, max_iter=200)
            m.fit(X_train, y_train)
            acc = float((m.predict(X_val) == y_val).mean())
            return {"accuracy": acc}

        result = cross_validate(X=X, y=y, classifier_func=clf_fn, cv=3)
        assert "mean_accuracy" in result
        assert result["cv_folds"] == 3


# ---------------------------------------------------------------------------
# k_fold_split
# ---------------------------------------------------------------------------


class TestKFoldSplit:
    def test_basic_split(self):
        X, y = _make_classification_data(n=50)
        folds = k_fold_split(X, y, k=5, random_state=42)
        assert len(folds) == 5
        for train_idx, test_idx in folds:
            assert len(train_idx) + len(test_idx) == 50

    def test_no_overlap(self):
        X, y = _make_classification_data(n=50)
        folds = k_fold_split(X, y, k=5, random_state=42)
        for train_idx, test_idx in folds:
            assert len(set(train_idx) & set(test_idx)) == 0


# ---------------------------------------------------------------------------
# learning_curve
# ---------------------------------------------------------------------------


class TestLearningCurve:
    def test_basic_learning_curve(self):
        from sklearn.tree import DecisionTreeClassifier

        X, y = _make_classification_data(n=80)
        result = learning_curve(
            X,
            y,
            model_factory=lambda: DecisionTreeClassifier(random_state=42),
            train_sizes=np.array([0.3, 0.6, 1.0]),
            cv=3,
            random_state=42,
        )
        assert "train_sizes" in result
        assert len(result["train_scores"]) == 3
        assert len(result["validation_scores"]) == 3
