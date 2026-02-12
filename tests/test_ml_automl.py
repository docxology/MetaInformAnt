"""Tests for ML AutoML optimization, model selection, and preprocessing.

Real implementation testing for automated machine learning methods.
No mocking used - all tests use real sklearn models and data.
"""

from __future__ import annotations

import math

import numpy as np
import pytest
from sklearn.datasets import make_classification, make_regression
from sklearn.ensemble import RandomForestClassifier
from sklearn.tree import DecisionTreeClassifier

from metainformant.ml.automl.optimization import (
    auto_preprocess,
    bayesian_optimization,
    grid_search,
    model_selection,
    random_search,
)

# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture()
def classification_data():
    """Small binary classification dataset."""
    X, y = make_classification(
        n_samples=100,
        n_features=8,
        n_informative=4,
        n_redundant=2,
        random_state=42,
    )
    return X, y


@pytest.fixture()
def regression_data():
    """Small regression dataset."""
    X, y = make_regression(
        n_samples=100,
        n_features=6,
        n_informative=3,
        noise=1.0,
        random_state=42,
    )
    return X, y


# ---------------------------------------------------------------------------
# Random Search
# ---------------------------------------------------------------------------


class TestRandomSearch:
    """Tests for random_search."""

    def test_basic_structure(self, classification_data):
        """Random search returns correct result keys."""
        X, y = classification_data

        def model_fn(n_estimators=10, max_depth=3):
            return RandomForestClassifier(
                n_estimators=int(n_estimators),
                max_depth=int(max_depth),
                random_state=42,
            )

        param_dists = {
            "n_estimators": [10, 20, 50],
            "max_depth": [2, 3, 5],
        }

        result = random_search(
            model_fn,
            param_dists,
            X,
            y,
            n_iter=6,
            cv=3,
            metric="accuracy",
            random_state=42,
        )

        assert "best_params" in result
        assert "best_score" in result
        assert "all_results" in result
        assert "n_evaluations" in result

        assert result["n_evaluations"] == 6
        assert len(result["all_results"]) == 6
        assert result["best_score"] > 0.0

    def test_continuous_param_distribution(self, classification_data):
        """Random search supports continuous parameter ranges."""
        X, y = classification_data

        def model_fn(max_depth=3):
            return DecisionTreeClassifier(
                max_depth=int(max_depth),
                random_state=42,
            )

        param_dists = {
            "max_depth": {"low": 2, "high": 10, "type": "int"},
        }

        result = random_search(
            model_fn,
            param_dists,
            X,
            y,
            n_iter=5,
            cv=3,
            random_state=42,
        )

        assert result["n_evaluations"] == 5
        # Best max_depth should be in range
        assert 2 <= result["best_params"]["max_depth"] <= 10

    def test_log_scale_param(self, classification_data):
        """Random search supports log-scale parameter sampling."""
        X, y = classification_data

        from sklearn.linear_model import LogisticRegression

        def model_fn(C=1.0):
            return LogisticRegression(C=C, max_iter=500, random_state=42)

        param_dists = {
            "C": {"low": 0.001, "high": 100.0, "log": True},
        }

        result = random_search(
            model_fn,
            param_dists,
            X,
            y,
            n_iter=8,
            cv=3,
            random_state=42,
        )
        assert result["best_params"]["C"] >= 0.001
        assert result["best_params"]["C"] <= 100.0

    def test_empty_param_dists_raises(self, classification_data):
        """Empty param_distributions raises ValueError."""
        X, y = classification_data

        def model_fn():
            return DecisionTreeClassifier()

        with pytest.raises(ValueError, match="must not be empty"):
            random_search(model_fn, {}, X, y)

    def test_reproducibility(self, classification_data):
        """Same random_state yields identical results."""
        X, y = classification_data

        def model_fn(max_depth=3):
            return DecisionTreeClassifier(max_depth=int(max_depth), random_state=42)

        dists = {"max_depth": [2, 3, 5, 7]}
        r1 = random_search(model_fn, dists, X, y, n_iter=4, random_state=99)
        r2 = random_search(model_fn, dists, X, y, n_iter=4, random_state=99)
        assert r1["best_params"] == r2["best_params"]
        assert r1["best_score"] == r2["best_score"]


# ---------------------------------------------------------------------------
# Bayesian Optimization
# ---------------------------------------------------------------------------


class TestBayesianOptimization:
    """Tests for bayesian_optimization."""

    def test_basic_structure(self):
        """Bayesian optimization returns correct structure."""

        def objective(params):
            # Simple quadratic with max at x=0.5
            return -((params["x"] - 0.5) ** 2)

        param_space = {"x": {"low": 0.0, "high": 1.0}}

        result = bayesian_optimization(
            objective,
            param_space,
            n_iter=15,
            n_initial=5,
            random_state=42,
        )

        assert "best_params" in result
        assert "best_score" in result
        assert "history" in result
        assert "surrogate_model" in result

        assert len(result["history"]) == 15
        # Should find something close to x=0.5
        assert abs(result["best_params"]["x"] - 0.5) < 0.35

    def test_multidimensional(self):
        """Bayesian optimization works in multiple dimensions."""

        def objective(params):
            return -((params["a"] - 0.3) ** 2) - (params["b"] - 0.7) ** 2

        param_space = {
            "a": {"low": 0.0, "high": 1.0},
            "b": {"low": 0.0, "high": 1.0},
        }

        result = bayesian_optimization(
            objective,
            param_space,
            n_iter=20,
            n_initial=5,
            random_state=42,
        )

        assert result["best_score"] > -0.5  # Reasonable optimization

    def test_log_scale_param(self):
        """Bayesian optimization supports log-scale parameters."""

        def objective(params):
            return -abs(math.log10(params["lr"]) + 2)  # optimal at lr=0.01

        param_space = {"lr": {"low": 0.0001, "high": 1.0, "log": True}}

        result = bayesian_optimization(
            objective,
            param_space,
            n_iter=15,
            n_initial=5,
            random_state=42,
        )
        assert 0.0001 <= result["best_params"]["lr"] <= 1.0

    def test_empty_param_space_raises(self):
        """Empty param_space raises ValueError."""
        with pytest.raises(ValueError, match="must not be empty"):
            bayesian_optimization(lambda p: 0.0, {})

    def test_n_initial_greater_than_n_iter_raises(self):
        """n_initial > n_iter raises ValueError."""
        with pytest.raises(ValueError, match="n_initial"):
            bayesian_optimization(
                lambda p: 0.0,
                {"x": {"low": 0, "high": 1}},
                n_iter=3,
                n_initial=10,
            )

    def test_surrogate_model_info(self):
        """Surrogate model info contains GP metadata."""

        def objective(params):
            return -params["x"] ** 2

        result = bayesian_optimization(
            objective,
            {"x": {"low": -1.0, "high": 1.0}},
            n_iter=10,
            n_initial=3,
            random_state=42,
        )

        sm = result["surrogate_model"]
        assert sm["type"] == "gaussian_process"
        assert sm["kernel"] == "rbf"
        assert sm["n_observations"] == 10

    def test_history_tracks_all_evaluations(self):
        """History should contain one entry per iteration."""
        counter = {"n": 0}

        def objective(params):
            counter["n"] += 1
            return -params["x"] ** 2

        result = bayesian_optimization(
            objective,
            {"x": {"low": 0.0, "high": 1.0}},
            n_iter=8,
            n_initial=3,
            random_state=42,
        )
        assert len(result["history"]) == 8
        assert counter["n"] == 8


# ---------------------------------------------------------------------------
# Grid Search
# ---------------------------------------------------------------------------


class TestGridSearch:
    """Tests for grid_search."""

    def test_basic_structure(self, classification_data):
        """Grid search returns correct structure."""
        X, y = classification_data

        def model_fn(max_depth=3, min_samples_leaf=1):
            return DecisionTreeClassifier(
                max_depth=max_depth,
                min_samples_leaf=min_samples_leaf,
                random_state=42,
            )

        param_grid = {
            "max_depth": [2, 3, 5],
            "min_samples_leaf": [1, 3],
        }

        result = grid_search(model_fn, param_grid, X, y, cv=3)

        assert "best_params" in result
        assert "best_score" in result
        assert "all_results" in result

        # 3 x 2 = 6 combinations
        assert len(result["all_results"]) == 6
        assert result["best_score"] > 0.0

    def test_single_param(self, classification_data):
        """Grid search works with single parameter."""
        X, y = classification_data

        def model_fn(max_depth=3):
            return DecisionTreeClassifier(max_depth=max_depth, random_state=42)

        result = grid_search(
            model_fn,
            {"max_depth": [2, 3, 5, 10]},
            X,
            y,
            cv=3,
        )
        assert len(result["all_results"]) == 4

    def test_empty_param_grid_raises(self, classification_data):
        """Empty param_grid raises ValueError."""
        X, y = classification_data
        with pytest.raises(ValueError, match="must not be empty"):
            grid_search(lambda: DecisionTreeClassifier(), {}, X, y)

    def test_best_score_is_max(self, classification_data):
        """best_score should be the maximum of all scores."""
        X, y = classification_data

        def model_fn(max_depth=3):
            return DecisionTreeClassifier(max_depth=max_depth, random_state=42)

        result = grid_search(
            model_fn,
            {"max_depth": [1, 2, 5, 10]},
            X,
            y,
            cv=3,
        )
        all_scores = [r["score"] for r in result["all_results"] if not math.isnan(r["score"])]
        assert result["best_score"] == max(all_scores)


# ---------------------------------------------------------------------------
# Model Selection
# ---------------------------------------------------------------------------


class TestModelSelection:
    """Tests for model_selection."""

    def test_classification(self, classification_data):
        """Model selection for classification returns rankings."""
        X, y = classification_data
        result = model_selection(X, y, task="classification", cv=3)

        assert "best_model_type" in result
        assert "rankings" in result
        assert "cv_results_per_model" in result

        # Should have tried several models
        assert len(result["rankings"]) >= 3
        # Rankings should be sorted descending by score
        scores = [s for _, s in result["rankings"]]
        assert scores == sorted(scores, reverse=True)
        # Best model should match first in rankings
        assert result["best_model_type"] == result["rankings"][0][0]

    def test_regression(self, regression_data):
        """Model selection for regression returns rankings."""
        X, y = regression_data
        result = model_selection(X, y, task="regression", cv=3)

        assert len(result["rankings"]) >= 3
        assert result["best_model_type"] in result["cv_results_per_model"]

    def test_invalid_task_raises(self, classification_data):
        """Invalid task string raises ValueError."""
        X, y = classification_data
        with pytest.raises(ValueError, match="task must be"):
            model_selection(X, y, task="clustering")

    def test_cv_results_contain_score(self, classification_data):
        """Each model in cv_results_per_model has a score."""
        X, y = classification_data
        result = model_selection(X, y, task="classification", cv=3)
        for name, info in result["cv_results_per_model"].items():
            assert "score" in info
            assert "model_type" in info


# ---------------------------------------------------------------------------
# Auto Preprocess
# ---------------------------------------------------------------------------


class TestAutoPreprocess:
    """Tests for auto_preprocess."""

    def test_basic_structure(self, classification_data):
        """auto_preprocess returns correct structure."""
        X, y = classification_data
        result = auto_preprocess(X, y)

        assert "X_processed" in result
        assert "transformations_applied" in result
        assert "feature_info" in result

        X_proc = result["X_processed"]
        assert len(X_proc) == X.shape[0]
        assert len(X_proc[0]) == X.shape[1]
        assert len(result["feature_info"]) == X.shape[1]

    def test_standard_scaling_applied(self, classification_data):
        """Continuous features should be standard-scaled."""
        X, _ = classification_data
        result = auto_preprocess(X)

        # Check that at least some features got standard_scaling
        scaling_applied = any("standard_scaling" in info.get("transformations", []) for info in result["feature_info"])
        assert scaling_applied

        # Scaled features should have approx 0 mean
        X_proc = np.array(result["X_processed"])
        for i in range(X_proc.shape[1]):
            col = X_proc[:, i]
            info = result["feature_info"][i]
            if "standard_scaling" in info.get("transformations", []):
                assert abs(col.mean()) < 0.1

    def test_missing_value_imputation(self):
        """NaN values should be imputed."""
        X = np.array(
            [
                [1.0, 2.0, 3.0],
                [4.0, float("nan"), 6.0],
                [7.0, 8.0, float("nan")],
                [10.0, 11.0, 12.0],
            ]
        )
        result = auto_preprocess(X)

        X_proc = result["X_processed"]
        # No NaN in output
        for row in X_proc:
            for val in row:
                assert not math.isnan(val)

        # Imputation should be recorded
        imputed = [info for info in result["feature_info"] if "mean_imputation" in info.get("transformations", [])]
        assert len(imputed) >= 1

    def test_empty_X_raises(self):
        """Empty X raises ValueError."""
        empty = np.array([]).reshape(0, 3)
        with pytest.raises(ValueError, match="must not be empty"):
            auto_preprocess(empty)

    def test_binary_features_not_scaled(self):
        """Binary features should not be standard-scaled."""
        X = np.array(
            [
                [0.0, 1.5, 3.0],
                [1.0, 2.5, 6.0],
                [0.0, 3.5, 9.0],
                [1.0, 4.5, 12.0],
                [0.0, 5.5, 15.0],
            ]
        )
        result = auto_preprocess(X)

        # Feature 0 is binary (0, 1) - should not be scaled
        info_0 = result["feature_info"][0]
        assert info_0["is_binary"]
        assert "standard_scaling" not in info_0.get("transformations", [])

    def test_feature_info_metadata(self, classification_data):
        """Feature info should contain metadata about each feature."""
        X, _ = classification_data
        result = auto_preprocess(X)

        for info in result["feature_info"]:
            assert "index" in info
            assert "n_unique" in info
            assert "has_missing" in info
            assert "is_categorical" in info
            assert "is_binary" in info

    def test_no_y_required(self, classification_data):
        """auto_preprocess works without y."""
        X, _ = classification_data
        result = auto_preprocess(X)
        assert len(result["X_processed"]) == X.shape[0]
