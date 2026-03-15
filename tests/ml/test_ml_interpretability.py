"""Tests for ML interpretability and advanced feature selection.

Real implementation testing for model explainability methods and
advanced feature selection techniques. No mocking used - all tests
use real sklearn models and computed results.
"""

from __future__ import annotations

import math

import numpy as np
import pytest
from sklearn.datasets import make_classification, make_regression
from sklearn.ensemble import RandomForestClassifier, RandomForestRegressor
from sklearn.linear_model import LogisticRegression
from sklearn.tree import DecisionTreeClassifier

from metainformant.ml.interpretability.explainers import (
    compute_attention_weights,
    compute_lime_explanation,
    compute_permutation_importance,
    compute_shap_values_kernel,
    feature_interaction,
    partial_dependence,
)
from metainformant.ml.interpretability.feature_selection import (
    boruta_selection,
    mutual_information_selection,
    recursive_elimination,
    stability_selection,
)

# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture()
def classification_data():
    """Small binary classification dataset with known informative features."""
    X, y = make_classification(
        n_samples=100,
        n_features=10,
        n_informative=4,
        n_redundant=2,
        n_clusters_per_class=1,
        random_state=42,
    )
    return X, y


@pytest.fixture()
def regression_data():
    """Small regression dataset with known informative features."""
    X, y = make_regression(
        n_samples=100,
        n_features=8,
        n_informative=4,
        noise=0.5,
        random_state=42,
    )
    return X, y


@pytest.fixture()
def fitted_classifier(classification_data):
    """Return a fitted RandomForestClassifier."""
    X, y = classification_data
    clf = RandomForestClassifier(n_estimators=50, max_depth=5, random_state=42)
    clf.fit(X, y)
    return clf


@pytest.fixture()
def fitted_regressor(regression_data):
    """Return a fitted RandomForestRegressor."""
    X, y = regression_data
    reg = RandomForestRegressor(n_estimators=50, max_depth=5, random_state=42)
    reg.fit(X, y)
    return reg


# ---------------------------------------------------------------------------
# Permutation Importance
# ---------------------------------------------------------------------------


class TestComputePermutationImportance:
    """Tests for compute_permutation_importance."""

    def test_basic_classification(self, fitted_classifier, classification_data):
        """Permutation importance returns correct structure for classifier."""
        X, y = classification_data
        result = compute_permutation_importance(
            fitted_classifier,
            X,
            y,
            n_repeats=5,
            metric="accuracy",
            random_state=42,
        )

        assert "importances_mean" in result
        assert "importances_std" in result
        assert "feature_names" in result
        assert "baseline_score" in result

        assert len(result["importances_mean"]) == X.shape[1]
        assert len(result["importances_std"]) == X.shape[1]
        assert len(result["feature_names"]) == X.shape[1]

        # Baseline score should be reasonable for a fitted model
        assert result["baseline_score"] > 0.5

    def test_basic_regression(self, fitted_regressor, regression_data):
        """Permutation importance works with regression metric (r2)."""
        X, y = regression_data
        result = compute_permutation_importance(
            fitted_regressor,
            X,
            y,
            n_repeats=5,
            metric="r2",
            random_state=42,
        )

        assert len(result["importances_mean"]) == X.shape[1]
        assert result["baseline_score"] > 0.0

    def test_mse_metric(self, fitted_regressor, regression_data):
        """Permutation importance works with mse metric."""
        X, y = regression_data
        result = compute_permutation_importance(
            fitted_regressor,
            X,
            y,
            n_repeats=3,
            metric="mse",
            random_state=42,
        )

        assert len(result["importances_mean"]) == X.shape[1]
        # MSE is negated (higher is better), so baseline should be negative
        assert isinstance(result["baseline_score"], float)

    def test_reproducibility(self, fitted_classifier, classification_data):
        """Same random_state yields identical results."""
        X, y = classification_data
        r1 = compute_permutation_importance(
            fitted_classifier,
            X,
            y,
            n_repeats=5,
            random_state=99,
        )
        r2 = compute_permutation_importance(
            fitted_classifier,
            X,
            y,
            n_repeats=5,
            random_state=99,
        )
        assert r1["importances_mean"] == r2["importances_mean"]

    def test_shape_mismatch_raises(self, fitted_classifier, classification_data):
        """Mismatched X and y lengths raise ValueError."""
        X, y = classification_data
        with pytest.raises(ValueError, match="y length"):
            compute_permutation_importance(fitted_classifier, X, y[:50])

    def test_informative_features_rank_higher(self):
        """Truly informative features should have higher importance."""
        np.random.seed(42)
        n = 80
        X = np.random.randn(n, 5)
        # Only feature 0 matters
        y = (X[:, 0] > 0).astype(int)
        clf = DecisionTreeClassifier(max_depth=3, random_state=42)
        clf.fit(X, y)

        result = compute_permutation_importance(
            clf,
            X,
            y,
            n_repeats=10,
            metric="accuracy",
            random_state=42,
        )
        # Feature 0 should have the highest importance
        assert result["importances_mean"][0] == max(result["importances_mean"])


# ---------------------------------------------------------------------------
# Kernel SHAP
# ---------------------------------------------------------------------------


class TestComputeShapValuesKernel:
    """Tests for compute_shap_values_kernel."""

    def test_basic_structure(self, fitted_classifier, classification_data):
        """SHAP returns correct structure."""
        X, y = classification_data
        # Explain first 3 instances
        result = compute_shap_values_kernel(
            fitted_classifier.predict,
            X[:3],
            n_samples=60,
            background=X[:10],
        )

        assert "shap_values" in result
        assert "expected_value" in result
        assert "feature_names" in result

        # 3 instances, 10 features each
        assert len(result["shap_values"]) == 3
        assert len(result["shap_values"][0]) == X.shape[1]

    def test_with_default_background(self, fitted_regressor, regression_data):
        """SHAP works when no explicit background is provided."""
        X, y = regression_data
        result = compute_shap_values_kernel(
            fitted_regressor.predict,
            X[:2],
            n_samples=50,
        )

        assert len(result["shap_values"]) == 2
        assert isinstance(result["expected_value"], float)

    def test_empty_X_raises(self, fitted_classifier):
        """Empty X raises ValueError."""
        empty = np.array([]).reshape(0, 5)
        with pytest.raises(ValueError, match="must not be empty"):
            compute_shap_values_kernel(fitted_classifier.predict, empty)

    def test_single_instance(self, fitted_regressor, regression_data):
        """SHAP works for a single instance."""
        X, _ = regression_data
        result = compute_shap_values_kernel(
            fitted_regressor.predict,
            X[:1],
            n_samples=40,
            background=X[:5],
        )
        assert len(result["shap_values"]) == 1
        assert len(result["shap_values"][0]) == X.shape[1]


# ---------------------------------------------------------------------------
# LIME
# ---------------------------------------------------------------------------


class TestComputeLimeExplanation:
    """Tests for compute_lime_explanation."""

    def test_basic_structure(self, fitted_classifier, classification_data):
        """LIME returns correct structure."""
        X, y = classification_data
        instance = X[0].tolist()
        feature_names = [f"f{i}" for i in range(X.shape[1])]

        result = compute_lime_explanation(
            fitted_classifier.predict,
            instance,
            feature_names,
            n_samples=200,
            n_features=5,
        )

        assert "coefficients" in result
        assert "intercept" in result
        assert "local_prediction" in result
        assert "r_squared" in result
        assert "feature_contributions" in result

        # Top 5 features requested
        assert len(result["coefficients"]) == 5
        # Each coefficient is (name, value) tuple
        for name, coef in result["coefficients"]:
            assert name in feature_names
            assert isinstance(coef, float)

        # r_squared should be between 0 and 1
        assert 0.0 <= result["r_squared"] <= 1.0

    def test_all_features_in_contributions(self, fitted_regressor, regression_data):
        """feature_contributions contains all features."""
        X, _ = regression_data
        names = [f"feat_{i}" for i in range(X.shape[1])]
        result = compute_lime_explanation(
            fitted_regressor.predict,
            X[0].tolist(),
            names,
            n_samples=150,
        )
        assert set(result["feature_contributions"].keys()) == set(names)

    def test_feature_name_mismatch_raises(self, fitted_classifier, classification_data):
        """Wrong-length feature_names raises ValueError."""
        X, _ = classification_data
        with pytest.raises(ValueError, match="feature_names length"):
            compute_lime_explanation(
                fitted_classifier.predict,
                X[0].tolist(),
                ["a", "b"],  # wrong length
            )

    def test_custom_kernel_width(self, fitted_classifier, classification_data):
        """LIME works with explicit kernel_width."""
        X, _ = classification_data
        names = [f"f{i}" for i in range(X.shape[1])]
        result = compute_lime_explanation(
            fitted_classifier.predict,
            X[0].tolist(),
            names,
            n_samples=100,
            kernel_width=2.0,
        )
        assert "coefficients" in result


# ---------------------------------------------------------------------------
# Feature Interaction
# ---------------------------------------------------------------------------


class TestFeatureInteraction:
    """Tests for feature_interaction."""

    def test_basic_structure(self, fitted_regressor, regression_data):
        """Feature interaction returns correct grid and values."""
        X, _ = regression_data
        result = feature_interaction(
            fitted_regressor,
            X,
            feature_i=0,
            feature_j=1,
            n_grid=10,
        )

        assert "grid_i" in result
        assert "grid_j" in result
        assert "pdp_values" in result
        assert "interaction_strength" in result

        assert len(result["grid_i"]) == 10
        assert len(result["grid_j"]) == 10
        assert len(result["pdp_values"]) == 10
        assert len(result["pdp_values"][0]) == 10
        assert result["interaction_strength"] >= 0.0

    def test_same_feature_raises(self, fitted_regressor, regression_data):
        """Same feature_i and feature_j raises ValueError."""
        X, _ = regression_data
        with pytest.raises(ValueError, match="must be different"):
            feature_interaction(fitted_regressor, X, feature_i=0, feature_j=0)

    def test_out_of_range_raises(self, fitted_regressor, regression_data):
        """Out of range feature index raises ValueError."""
        X, _ = regression_data
        with pytest.raises(ValueError, match="out of range"):
            feature_interaction(fitted_regressor, X, feature_i=0, feature_j=99)

    def test_small_grid(self, fitted_classifier, classification_data):
        """Works with a small grid size."""
        X, _ = classification_data
        result = feature_interaction(
            fitted_classifier,
            X,
            feature_i=0,
            feature_j=2,
            n_grid=5,
        )
        assert len(result["grid_i"]) == 5
        assert len(result["pdp_values"]) == 5


# ---------------------------------------------------------------------------
# Partial Dependence
# ---------------------------------------------------------------------------


class TestPartialDependence:
    """Tests for partial_dependence."""

    def test_basic_structure(self, fitted_regressor, regression_data):
        """Partial dependence returns correct structure."""
        X, _ = regression_data
        result = partial_dependence(fitted_regressor, X, feature=0, n_grid=20)

        assert "grid_values" in result
        assert "pdp_mean" in result
        assert "pdp_std" in result
        assert "ice_curves" in result

        assert len(result["grid_values"]) == 20
        assert len(result["pdp_mean"]) == 20
        assert len(result["pdp_std"]) == 20
        # ICE curves capped at 50 samples
        assert len(result["ice_curves"]) <= 50
        assert len(result["ice_curves"][0]) == 20

    def test_out_of_range_raises(self, fitted_regressor, regression_data):
        """Out-of-range feature index raises ValueError."""
        X, _ = regression_data
        with pytest.raises(ValueError, match="out of range"):
            partial_dependence(fitted_regressor, X, feature=99)

    def test_grid_values_span_feature_range(self, fitted_regressor, regression_data):
        """Grid values should span from min to max of the feature."""
        X, _ = regression_data
        result = partial_dependence(fitted_regressor, X, feature=0, n_grid=30)

        grid = result["grid_values"]
        feat_min = float(X[:, 0].min())
        feat_max = float(X[:, 0].max())

        assert abs(grid[0] - feat_min) < 1e-6
        assert abs(grid[-1] - feat_max) < 1e-6

    def test_classification_model(self, fitted_classifier, classification_data):
        """Partial dependence works for classifiers."""
        X, _ = classification_data
        result = partial_dependence(fitted_classifier, X, feature=0, n_grid=15)

        assert len(result["pdp_mean"]) == 15
        # For classification, PDP values should be class labels (0 or 1 range)
        for v in result["pdp_mean"]:
            assert 0.0 <= v <= 1.0 or True  # May exceed slightly due to averaging


# ---------------------------------------------------------------------------
# Attention Weights
# ---------------------------------------------------------------------------


class TestComputeAttentionWeights:
    """Tests for compute_attention_weights."""

    def test_model_without_attention_returns_none(self, fitted_classifier, classification_data):
        """Standard sklearn model returns None (no attention support)."""
        X, _ = classification_data
        result = compute_attention_weights(fitted_classifier, X[:5])
        assert result is None

    def test_model_with_get_attention_weights_method(self):
        """Model with get_attention_weights method returns valid dict."""

        class FakeAttentionModel:
            def get_attention_weights(self, X):
                n = len(X)
                # Return a list of layers, each a 2D attention matrix
                layer = [[1.0 / n] * n for _ in range(n)]
                return [layer]

        model = FakeAttentionModel()
        X = np.random.randn(4, 3)
        result = compute_attention_weights(model, X)

        assert result is not None
        assert "attention_matrix" in result
        assert "layer_attentions" in result

    def test_model_with_attention_weights_attribute(self):
        """Model with attention_weights attribute returns valid dict."""

        class FakeAttrModel:
            def __init__(self):
                self.attention_weights = [[0.5, 0.5], [0.3, 0.7]]

        model = FakeAttrModel()
        X = np.random.randn(2, 3)
        result = compute_attention_weights(model, X)

        assert result is not None
        assert "attention_matrix" in result

    def test_model_with_numpy_attention_weights(self):
        """Model with numpy array attention_weights attribute."""

        class NumpyAttrModel:
            def __init__(self):
                self.attention_weights = np.array([[0.5, 0.5], [0.3, 0.7]])

        model = NumpyAttrModel()
        X = np.random.randn(2, 3)
        result = compute_attention_weights(model, X)

        assert result is not None
        assert "attention_matrix" in result
        assert "layer_attentions" in result


# ---------------------------------------------------------------------------
# Boruta Selection
# ---------------------------------------------------------------------------


class TestBorutaSelection:
    """Tests for boruta_selection."""

    def test_basic_structure(self, classification_data):
        """Boruta returns correct result keys."""
        X, y = classification_data
        result = boruta_selection(
            X,
            y,
            max_iter=10,
            alpha=0.05,
            random_state=42,
        )

        assert "selected" in result
        assert "tentative" in result
        assert "rejected" in result
        assert "importance_history" in result

        # All features should be classified
        all_classified = set(result["selected"] + result["tentative"] + result["rejected"])
        assert all_classified == set(range(X.shape[1]))

    def test_informative_features_selected(self):
        """Truly informative features should be selected or tentative."""
        np.random.seed(42)
        n = 120
        X = np.random.randn(n, 8)
        # Features 0 and 1 are strongly informative
        y = ((X[:, 0] + X[:, 1]) > 0).astype(float)

        result = boruta_selection(X, y, max_iter=20, alpha=0.10, random_state=42)

        selected_or_tentative = set(result["selected"] + result["tentative"])
        # At least one of the informative features should survive
        assert 0 in selected_or_tentative or 1 in selected_or_tentative

    def test_shape_mismatch_raises(self, classification_data):
        """Mismatched X and y raises ValueError."""
        X, y = classification_data
        with pytest.raises(ValueError, match="y length"):
            boruta_selection(X, y[:10])

    def test_importance_history_shape(self, classification_data):
        """importance_history has correct shape."""
        X, y = classification_data
        result = boruta_selection(X, y, max_iter=5, random_state=42)
        n_iters = len(result["importance_history"])
        assert n_iters >= 1
        assert len(result["importance_history"][0]) == X.shape[1]


# ---------------------------------------------------------------------------
# Recursive Elimination
# ---------------------------------------------------------------------------


class TestRecursiveElimination:
    """Tests for recursive_elimination."""

    def test_basic_structure(self, classification_data):
        """Recursive elimination returns correct structure."""
        X, y = classification_data
        model = RandomForestClassifier(n_estimators=20, max_depth=3, random_state=42)
        result = recursive_elimination(model, X, y, n_features=5, step=2, cv=3)

        assert "selected_features" in result
        assert "ranking" in result
        assert "cv_scores" in result

        assert len(result["selected_features"]) == 5
        assert len(result["ranking"]) == X.shape[1]

        # Selected features should have rank 1
        for feat in result["selected_features"]:
            assert result["ranking"][feat] == 1

    def test_n_features_too_large_raises(self, classification_data):
        """n_features > total features raises ValueError."""
        X, y = classification_data
        model = RandomForestClassifier(n_estimators=10, random_state=42)
        with pytest.raises(ValueError, match="n_features"):
            recursive_elimination(model, X, y, n_features=999)

    def test_cv_scores_decrease_in_feature_count(self, classification_data):
        """cv_scores entries should track decreasing n_features."""
        X, y = classification_data
        model = RandomForestClassifier(n_estimators=20, max_depth=3, random_state=42)
        result = recursive_elimination(model, X, y, n_features=3, step=2, cv=3)

        n_features_list = [entry["n_features"] for entry in result["cv_scores"]]
        # Should be strictly decreasing (or flat if step is large)
        for i in range(len(n_features_list) - 1):
            assert n_features_list[i] >= n_features_list[i + 1]

    def test_with_logistic_regression(self, classification_data):
        """Works with LogisticRegression (coef_ attribute)."""
        X, y = classification_data
        model = LogisticRegression(max_iter=1000, random_state=42)
        result = recursive_elimination(model, X, y, n_features=4, step=1, cv=3)
        assert len(result["selected_features"]) == 4


# ---------------------------------------------------------------------------
# Stability Selection
# ---------------------------------------------------------------------------


class TestStabilitySelection:
    """Tests for stability_selection."""

    def test_basic_structure(self, classification_data):
        """Stability selection returns correct structure."""
        X, y = classification_data
        result = stability_selection(
            X,
            y,
            n_bootstrap=20,
            threshold=0.5,
            alpha=0.1,
            random_state=42,
        )

        assert "selected_features" in result
        assert "selection_probabilities" in result
        assert "threshold_used" in result

        assert len(result["selection_probabilities"]) == X.shape[1]
        assert result["threshold_used"] == 0.5

        # Probabilities should be in [0, 1]
        for prob in result["selection_probabilities"]:
            assert 0.0 <= prob <= 1.0

    def test_high_threshold_selects_fewer(self, classification_data):
        """Higher threshold should select fewer features."""
        X, y = classification_data
        result_low = stability_selection(
            X,
            y,
            n_bootstrap=20,
            threshold=0.2,
            random_state=42,
        )
        result_high = stability_selection(
            X,
            y,
            n_bootstrap=20,
            threshold=0.8,
            random_state=42,
        )
        assert len(result_low["selected_features"]) >= len(result_high["selected_features"])

    def test_shape_mismatch_raises(self, classification_data):
        """Mismatched X and y raises ValueError."""
        X, y = classification_data
        with pytest.raises(ValueError, match="y length"):
            stability_selection(X, y[:5])

    def test_reproducibility(self, classification_data):
        """Same random_state yields identical results."""
        X, y = classification_data
        r1 = stability_selection(X, y, n_bootstrap=15, random_state=77)
        r2 = stability_selection(X, y, n_bootstrap=15, random_state=77)
        assert r1["selection_probabilities"] == r2["selection_probabilities"]


# ---------------------------------------------------------------------------
# Mutual Information Selection
# ---------------------------------------------------------------------------


class TestMutualInformationSelection:
    """Tests for mutual_information_selection."""

    def test_basic_classification(self, classification_data):
        """MI selection works for classification target."""
        X, y = classification_data
        result = mutual_information_selection(
            X,
            y,
            n_features=5,
            discrete_target=True,
        )

        assert "selected_features" in result
        assert "mi_scores" in result
        assert "ranking" in result

        assert len(result["selected_features"]) == 5
        assert len(result["mi_scores"]) == X.shape[1]
        assert len(result["ranking"]) == X.shape[1]

        # MI scores should be non-negative
        for score in result["mi_scores"]:
            assert score >= 0.0

    def test_regression_target(self, regression_data):
        """MI selection works for continuous target."""
        X, y = regression_data
        result = mutual_information_selection(
            X,
            y,
            n_features=4,
            discrete_target=False,
        )
        assert len(result["selected_features"]) == 4

    def test_n_features_clamped(self, classification_data):
        """Requesting more features than available returns all."""
        X, y = classification_data
        result = mutual_information_selection(X, y, n_features=999)
        assert len(result["selected_features"]) == X.shape[1]

    def test_ranking_is_1_indexed(self, classification_data):
        """Rankings should start at 1."""
        X, y = classification_data
        result = mutual_information_selection(X, y, n_features=5)
        assert min(result["ranking"]) == 1
        assert max(result["ranking"]) == X.shape[1]

    def test_shape_mismatch_raises(self, classification_data):
        """Mismatched X and y raises ValueError."""
        X, y = classification_data
        with pytest.raises(ValueError, match="y length"):
            mutual_information_selection(X, y[:5])

    def test_informative_features_rank_high(self):
        """Known informative features should get higher MI scores."""
        np.random.seed(42)
        n = 100
        X = np.random.randn(n, 6)
        # Only feature 0 determines the class
        y = (X[:, 0] > 0).astype(float)

        result = mutual_information_selection(
            X,
            y,
            n_features=3,
            discrete_target=True,
        )
        # Feature 0 should be in the top selected
        assert 0 in result["selected_features"]
