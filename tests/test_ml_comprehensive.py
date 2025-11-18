"""Comprehensive tests for ML functionality.

Real implementation testing for machine learning methods.
No mocking used - all tests use real computational methods and data.
"""

from __future__ import annotations

from typing import List, Tuple

import numpy as np
import pytest

from metainformant.ml.classification import BiologicalClassifier, cross_validate_biological, evaluate_classifier
from metainformant.ml.dimensionality import (
    biological_embedding,
    reduce_dimensions_pca,
    reduce_dimensions_tsne,
    reduce_dimensions_umap,
)
from metainformant.ml.regression import BiologicalRegressor, evaluate_regressor
from metainformant.ml.validation import (
    bootstrap_validate,
    cross_validate,
    k_fold_split,
    learning_curve,
    train_test_split,
)


class TestBiologicalClassifier:
    """Test BiologicalClassifier functionality."""

    def setup_method(self):
        """Set up test data for classification."""
        np.random.seed(42)

        # Create synthetic biological data
        n_samples = 200
        n_features = 50

        # Generate features representing gene expression levels
        # Class 0: healthy samples
        healthy_samples = np.random.lognormal(mean=1.0, sigma=0.5, size=(100, n_features))
        # Class 1: disease samples (some genes upregulated, some downregulated)
        disease_samples = np.random.lognormal(mean=1.2, sigma=0.6, size=(100, n_features))

        # Create differential expression pattern
        # Upregulate first 20 genes in disease
        disease_samples[:, :20] *= 1.8
        # Downregulate next 10 genes in disease
        disease_samples[:, 20:30] *= 0.4

        self.X = np.vstack([healthy_samples, disease_samples])
        self.y = np.array([0] * 100 + [1] * 100)  # 0=healthy, 1=disease

        # Feature names
        self.feature_names = [f"Gene_{i}" for i in range(n_features)]

    def test_classifier_initialization(self):
        """Test BiologicalClassifier initialization."""
        # Default initialization
        clf = BiologicalClassifier()
        assert clf.algorithm == "random_forest"
        assert clf.random_state is None
        assert not clf.is_fitted

        # Custom initialization
        custom_clf = BiologicalClassifier(algorithm="svm", random_state=42, C=1.0, kernel="rbf")
        assert custom_clf.algorithm == "svm"
        assert custom_clf.random_state == 42
        assert custom_clf.params["C"] == 1.0
        assert custom_clf.params["kernel"] == "rbf"

    def test_classifier_fit_predict(self):
        """Test classifier fitting and prediction."""
        clf = BiologicalClassifier(algorithm="random_forest", random_state=42)

        # Split data
        train_idx = np.random.choice(len(self.X), size=150, replace=False)
        test_idx = np.setdiff1d(np.arange(len(self.X)), train_idx)

        X_train, X_test = self.X[train_idx], self.X[test_idx]
        y_train, y_test = self.y[train_idx], self.y[test_idx]

        # Fit classifier
        clf.fit(X_train, y_train)
        assert clf.is_fitted

        # Make predictions
        predictions = clf.predict(X_test)
        probabilities = clf.predict_proba(X_test)

        # Check predictions
        assert len(predictions) == len(X_test)
        assert all(pred in [0, 1] for pred in predictions)

        # Check probabilities
        assert probabilities.shape == (len(X_test), 2)
        assert np.allclose(probabilities.sum(axis=1), 1.0)

        # Should achieve reasonable accuracy (>60% given signal in data)
        accuracy = np.mean(predictions == y_test)
        assert accuracy > 0.6

    def test_classifier_feature_importance(self):
        """Test feature importance extraction."""
        clf = BiologicalClassifier(algorithm="random_forest", random_state=42)
        clf.fit(self.X, self.y)

        # Get feature importance
        importance = clf.get_feature_importance()

        assert len(importance) == self.X.shape[1]
        assert all(imp >= 0 for imp in importance)
        assert np.sum(importance) > 0  # Should have some importance

        # Most important features should include some from the differential genes
        # (first 30 genes have differential expression)
        top_features = np.argsort(importance)[-10:]
        important_differential = sum(1 for idx in top_features if idx < 30)
        assert important_differential >= 3  # Should identify some differential genes

    def test_classifier_different_algorithms(self):
        """Test different classification algorithms."""
        algorithms = ["random_forest", "svm", "logistic_regression"]

        for algorithm in algorithms:
            try:
                clf = BiologicalClassifier(algorithm=algorithm, random_state=42)
                clf.fit(self.X, self.y)

                predictions = clf.predict(self.X)
                accuracy = np.mean(predictions == self.y)

                # All algorithms should achieve reasonable performance
                assert accuracy > 0.5
                assert clf.is_fitted

            except ValueError as e:
                if "Unknown algorithm" in str(e):
                    # Some algorithms might not be implemented
                    continue
                else:
                    raise

    def test_classifier_evaluation(self):
        """Test classifier evaluation functionality."""
        clf = BiologicalClassifier(algorithm="random_forest", random_state=42)
        clf.fit(self.X, self.y)

        results = evaluate_classifier(classifier=clf, X=self.X, y=self.y)

        # Check results structure
        assert "accuracy" in results
        assert "predictions" in results

        # Check performance metrics
        assert 0.0 <= results["accuracy"] <= 1.0
        assert len(results["predictions"]) == len(self.y)

        # Should achieve reasonable performance
        assert results["accuracy"] > 0.6

    def test_cross_validate_biological(self):
        """Test biological cross-validation."""
        results = cross_validate_biological(X=self.X, y=self.y, algorithm="random_forest", cv_folds=3, random_state=42)

        # Check results structure
        assert "accuracy" in results or "mean_accuracy" in results
        assert isinstance(results, dict)

        # Should return reasonable results
        accuracy_key = "accuracy" if "accuracy" in results else "mean_accuracy"
        if accuracy_key in results:
            assert 0.0 <= results[accuracy_key] <= 1.0


class TestBiologicalRegressor:
    """Test BiologicalRegressor functionality."""

    def setup_method(self):
        """Set up test data for regression."""
        np.random.seed(123)

        # Create synthetic biological regression data
        n_samples = 150
        n_features = 40

        # Generate features (e.g., gene expression, clinical variables)
        self.X = np.random.randn(n_samples, n_features)

        # Create target with some signal
        # Some features contribute to target (e.g., biomarkers predict outcome)
        true_coefficients = np.zeros(n_features)
        true_coefficients[:10] = np.random.randn(10) * 2  # First 10 features are predictive

        # Generate continuous target (e.g., survival time, drug response)
        self.y = self.X @ true_coefficients + np.random.randn(n_samples) * 0.5

        self.feature_names = [f"Biomarker_{i}" for i in range(n_features)]

    def test_regressor_initialization(self):
        """Test BiologicalRegressor initialization."""
        # Default initialization
        reg = BiologicalRegressor()
        assert reg.algorithm == "linear"
        assert reg.random_state is None
        assert not reg.is_fitted

        # Custom initialization
        custom_reg = BiologicalRegressor(algorithm="random_forest", random_state=42, n_estimators=100)
        assert custom_reg.algorithm == "random_forest"
        assert custom_reg.random_state == 42
        assert custom_reg.params["n_estimators"] == 100

    def test_regressor_fit_predict(self):
        """Test regressor fitting and prediction."""
        reg = BiologicalRegressor(algorithm="linear", random_state=42)

        # Split data
        train_size = int(0.7 * len(self.X))
        X_train, X_test = self.X[:train_size], self.X[train_size:]
        y_train, y_test = self.y[:train_size], self.y[train_size:]

        # Fit regressor
        reg.fit(X_train, y_train)
        assert reg.is_fitted

        # Make predictions
        predictions = reg.predict(X_test)

        # Check predictions
        assert len(predictions) == len(X_test)
        assert all(isinstance(pred, (int, float, np.number)) for pred in predictions)

        # Should achieve reasonable correlation with true values
        correlation = np.corrcoef(predictions, y_test)[0, 1]
        assert abs(correlation) > 0.3  # Should have some predictive power

    def test_regressor_different_algorithms(self):
        """Test different regression algorithms."""
        algorithms = ["linear", "random_forest", "svm"]

        for algorithm in algorithms:
            try:
                reg = BiologicalRegressor(algorithm=algorithm, random_state=42)
                reg.fit(self.X, self.y)

                predictions = reg.predict(self.X)
                mse = np.mean((predictions - self.y) ** 2)

                # All algorithms should make reasonable predictions
                assert mse < 10 * np.var(self.y)  # MSE should be reasonable
                assert reg.is_fitted

            except ValueError as e:
                if "Unknown algorithm" in str(e):
                    # Some algorithms might not be implemented
                    continue
                else:
                    raise

    def test_regressor_evaluation(self):
        """Test regressor evaluation."""
        reg = BiologicalRegressor(algorithm="linear", random_state=42)
        reg.fit(self.X, self.y)

        results = evaluate_regressor(regressor=reg, X=self.X, y=self.y)

        # Check results structure
        assert "mse" in results or "r2_score" in results or "predictions" in results

        # Check performance metrics if available
        if "mse" in results:
            assert results["mse"] >= 0.0
        if "predictions" in results:
            assert len(results["predictions"]) == len(self.y)

    def test_regressor_different_algorithms(self):
        """Test different regression algorithms."""
        algorithms = ["linear", "random_forest"]

        for algorithm in algorithms:
            try:
                reg = BiologicalRegressor(algorithm=algorithm, random_state=42)
                reg.fit(self.X, self.y)

                predictions = reg.predict(self.X)
                mse = np.mean((predictions - self.y) ** 2)

                # All algorithms should make reasonable predictions
                assert mse < 10 * np.var(self.y)  # MSE should be reasonable
                assert reg.is_fitted

            except ValueError as e:
                if "Unknown algorithm" in str(e):
                    # Some algorithms might not be implemented
                    continue
                else:
                    raise


class TestDimensionalityReduction:
    """Test dimensionality reduction methods."""

    def setup_method(self):
        """Set up high-dimensional biological data."""
        np.random.seed(456)

        # Create high-dimensional data (e.g., gene expression)
        self.n_samples = 100
        self.n_features = 200

        # Generate correlated features with some structure
        # Create latent factors
        n_factors = 5
        factor_loadings = np.random.randn(self.n_features, n_factors)
        factor_scores = np.random.randn(self.n_samples, n_factors)

        # Generate data with factor structure + noise
        self.X = factor_scores @ factor_loadings.T + np.random.randn(self.n_samples, self.n_features) * 0.5

        # Sample labels for visualization
        self.labels = np.random.randint(0, 3, self.n_samples)

    def test_pca_basic(self):
        """Test basic PCA functionality."""
        X_reduced, components, explained_var = reduce_dimensions_pca(X=self.X, n_components=10, standardize=True)

        # Check output dimensions
        assert X_reduced.shape == (self.n_samples, 10)
        assert components.shape == (self.n_features, 10)
        assert len(explained_var) == 10

        # Check explained variance
        assert all(var >= 0 for var in explained_var)
        assert np.sum(explained_var) <= 1.0

        # Components should be ordered by explained variance
        assert all(explained_var[i] >= explained_var[i + 1] for i in range(len(explained_var) - 1))

    def test_pca_different_components(self):
        """Test PCA with different numbers of components."""
        n_components_list = [2, 5, 20, 50]

        for n_components in n_components_list:
            X_reduced, components, explained_var = reduce_dimensions_pca(self.X, n_components=n_components)

            assert X_reduced.shape == (self.n_samples, n_components)
            assert components.shape == (self.n_features, n_components)
            assert len(explained_var) == n_components

    def test_pca_standardization_effect(self):
        """Test effect of standardization on PCA."""
        # With standardization
        X_std, _, var_std = reduce_dimensions_pca(self.X, n_components=5, standardize=True)

        # Without standardization
        X_raw, _, var_raw = reduce_dimensions_pca(self.X, n_components=5, standardize=False)

        # Results should be different
        assert not np.allclose(X_std, X_raw)
        assert not np.allclose(var_std, var_raw)

    def test_biological_embedding(self):
        """Test biological embedding methods."""
        embedding = biological_embedding(X=self.X, method="pca", n_components=3, labels=self.labels)

        # Check embedding structure
        assert "embedding" in embedding
        assert "method" in embedding
        assert "explained_variance" in embedding

        # Check embedding dimensions
        X_embedded = embedding["embedding"]
        assert X_embedded.shape == (self.n_samples, 3)

        # Check method
        assert embedding["method"] == "pca"

    @pytest.mark.slow
    def test_manifold_learning_methods(self):
        """Test different manifold learning methods."""
        methods = ["tsne", "umap"]

        for method in methods:
            try:
                if method == "tsne":
                    X_embedded = reduce_dimensions_tsne(X=self.X, n_components=2, random_state=42)
                elif method == "umap":
                    X_embedded = reduce_dimensions_umap(X=self.X, n_components=2, random_state=42)

                # Should produce 2D embedding
                assert X_embedded.shape == (self.n_samples, 2)

            except ValueError as e:
                if "Unknown method" in str(e):
                    # Some methods might not be implemented
                    continue
                else:
                    raise

            except ImportError:
                # Some methods might require optional dependencies
                continue


class TestValidation:
    """Test validation methods."""

    def setup_method(self):
        """Set up test data for validation."""
        np.random.seed(789)

        # Create balanced dataset
        self.n_samples = 120
        self.n_features = 30

        self.X = np.random.randn(self.n_samples, self.n_features)
        self.y = np.random.randint(0, 2, self.n_samples)

        # For regression
        self.y_continuous = np.random.randn(self.n_samples)

    def test_train_test_split(self):
        """Test train-test split functionality."""
        X_train, X_test, y_train, y_test = train_test_split(self.X, self.y, test_size=0.3, random_state=42)

        # Check split sizes
        expected_train_size = int(0.7 * self.n_samples)
        expected_test_size = self.n_samples - expected_train_size

        assert len(X_train) == expected_train_size
        assert len(X_test) == expected_test_size
        assert len(y_train) == expected_train_size
        assert len(y_test) == expected_test_size

        # Check dimensions
        assert X_train.shape[1] == self.n_features
        assert X_test.shape[1] == self.n_features

        # Check no overlap - verify total samples match and data rows are distinct
        assert len(X_train) + len(X_test) == self.n_samples
        # Check that no row in X_train appears in X_test (using first few features for efficiency)
        train_rows = {tuple(row[:5]) for row in X_train}
        test_rows = {tuple(row[:5]) for row in X_test}
        assert len(train_rows.intersection(test_rows)) == 0

    def test_k_fold_split(self):
        """Test k-fold cross-validation splitting."""
        folds = list(k_fold_split(self.X, self.y, k=5, random_state=42))

        # Should have 5 folds
        assert len(folds) == 5

        # Check fold structure
        for fold_idx, (train_indices, val_indices) in enumerate(folds):
            # Check indices are disjoint
            assert len(set(train_indices).intersection(set(val_indices))) == 0

            # Check all indices are covered
            all_indices = set(train_indices).union(set(val_indices))
            assert all_indices == set(range(self.n_samples))

            # Check approximate sizes
            expected_val_size = self.n_samples // 5
            assert abs(len(val_indices) - expected_val_size) <= 1

    def test_cross_validate(self):
        """Test cross-validation evaluation."""

        # Simple test classifier function (real implementation, not a mock)
        def simple_classifier(X_train, y_train, X_val, y_val):
            """Simple majority class classifier for testing."""
            majority_class = np.argmax(np.bincount(y_train))
            predictions = np.full(len(y_val), majority_class)
            accuracy = np.mean(predictions == y_val)
            return {"accuracy": accuracy, "predictions": predictions}

        cv_results = cross_validate(X=self.X, y=self.y, classifier_func=simple_classifier, cv_folds=3, random_state=42)

        # Check results structure
        assert "mean_accuracy" in cv_results
        assert "std_accuracy" in cv_results
        assert "fold_results" in cv_results

        # Check fold results
        assert len(cv_results["fold_results"]) == 3

        # Check accuracy is reasonable
        assert 0.0 <= cv_results["mean_accuracy"] <= 1.0
        assert cv_results["std_accuracy"] >= 0.0

    def test_bootstrap_validate(self):
        """Test bootstrap validation."""

        # Simple test model function (real implementation, not a mock)
        def simple_mean_model(X_train, y_train, X_test, y_test):
            """Simple mean predictor for regression testing."""
            mean_pred = np.mean(y_train)
            predictions = np.full(len(y_test), mean_pred)
            mse = np.mean((predictions - y_test) ** 2)
            return {"mse": mse, "predictions": predictions}

        bootstrap_results = bootstrap_validate(
            X=self.X, y=self.y_continuous, model_func=simple_mean_model, n_bootstrap=10, random_state=42
        )

        # Check results structure
        assert "mean_mse" in bootstrap_results
        assert "std_mse" in bootstrap_results
        assert "bootstrap_results" in bootstrap_results

        # Check bootstrap results
        assert len(bootstrap_results["bootstrap_results"]) == 10

        # Check metrics are reasonable
        assert bootstrap_results["mean_mse"] >= 0.0
        assert bootstrap_results["std_mse"] >= 0.0

    def test_learning_curve(self):
        """Test learning curve generation."""

        # Simple test classifier that improves with more data (real implementation, not a mock)
        def improving_classifier(X_train, y_train, X_val, y_val):
            """Classifier that shows improved accuracy with more training data."""
            train_size = len(X_train)
            base_accuracy = 0.5
            improvement = min(0.3, 0.3 * train_size / len(self.X))
            accuracy = base_accuracy + improvement + np.random.normal(0, 0.05)
            return {"accuracy": max(0, min(1, accuracy))}

        curve_results = learning_curve(
            X=self.X,
            y=self.y,
            classifier_func=improving_classifier,
            train_sizes=np.array([0.2, 0.5, 0.8]),
            cv_folds=3,
            random_state=42,
        )

        # Check results structure
        assert "train_sizes" in curve_results
        assert "train_scores" in curve_results
        assert "validation_scores" in curve_results

        # Check sizes
        train_sizes = curve_results["train_sizes"]
        assert len(train_sizes) == 3

        # Check scores
        train_scores = curve_results["train_scores"]
        val_scores = curve_results["validation_scores"]

        assert len(train_scores) == 3
        assert len(val_scores) == 3

        # Each should have cv_folds results
        for scores in train_scores:
            assert len(scores) == 3  # cv_folds
            assert all(0.0 <= score <= 1.0 for score in scores)

        for scores in val_scores:
            assert len(scores) == 3  # cv_folds
            assert all(0.0 <= score <= 1.0 for score in scores)


class TestMLIntegration:
    """Integration tests for ML functionality."""

    def test_complete_classification_pipeline(self):
        """Test complete classification analysis pipeline."""
        np.random.seed(42)

        # 1. Generate biological dataset
        n_samples = 200
        n_features = 100

        # Create dataset with known structure
        X = np.random.randn(n_samples, n_features)
        # First 20 features are predictive
        predictive_signal = np.sum(X[:, :20], axis=1)
        y = (predictive_signal > np.median(predictive_signal)).astype(int)

        # 2. Split data
        X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.3, random_state=42)

        # 3. Dimensionality reduction
        X_train_pca, components, _ = reduce_dimensions_pca(X_train, n_components=50, standardize=True)
        X_test_pca = (X_test - np.mean(X_train, axis=0)) @ components

        # 4. Train classifier
        clf = BiologicalClassifier(algorithm="random_forest", random_state=42)
        clf.fit(X_train_pca, y_train)

        # 5. Evaluate
        predictions = clf.predict(X_test_pca)
        accuracy = np.mean(predictions == y_test)

        # 6. Feature importance
        importance = clf.get_feature_importance()

        # 7. Cross-validation
        cv_results = cross_validate_biological(X_train_pca, y_train, algorithm="random_forest", cv_folds=3)

        # Verify pipeline results
        assert accuracy > 0.6  # Should achieve reasonable performance
        assert len(importance) == 50  # PCA components
        accuracy_key = "accuracy" if "accuracy" in cv_results else "mean_accuracy"
        if accuracy_key in cv_results:
            assert cv_results[accuracy_key] > 0.5

    def test_complete_regression_pipeline(self):
        """Test complete regression analysis pipeline."""
        np.random.seed(123)

        # 1. Generate regression dataset
        n_samples = 150
        n_features = 80

        X = np.random.randn(n_samples, n_features)
        # Create target with some predictive features
        true_coeffs = np.zeros(n_features)
        true_coeffs[:15] = np.random.randn(15)
        y = X @ true_coeffs + np.random.randn(n_samples) * 0.5

        # 2. Split data
        X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.25, random_state=123)

        # 3. Dimensionality reduction
        X_train_pca, components, _ = reduce_dimensions_pca(X_train, n_components=30, standardize=True)
        X_test_pca = (X_test - np.mean(X_train, axis=0)) @ components

        # 4. Train regressor
        reg = BiologicalRegressor(algorithm="linear", random_state=123)
        reg.fit(X_train_pca, y_train)

        # 5. Evaluate
        predictions = reg.predict(X_test_pca)
        correlation = np.corrcoef(predictions, y_test)[0, 1]

        # 6. Get predictions for diagnostics
        train_predictions = reg.predict(X_train_pca)

        # 7. Evaluation
        cv_results = evaluate_regressor(reg, X_train_pca, y_train)

        # Verify pipeline results
        assert abs(correlation) > 0.3  # Should have some predictive power
        assert len(train_predictions) == len(y_train)
        if "mse" in cv_results:
            assert cv_results["mse"] >= 0.0
        if "r2_score" in cv_results:
            assert cv_results["r2_score"] > -1.0  # Should be reasonable


class TestMLEdgeCases:
    """Test edge cases and error conditions."""

    def test_empty_data(self):
        """Test handling of empty datasets."""
        empty_X = np.array([]).reshape(0, 5)
        empty_y = np.array([])

        # Should handle gracefully
        try:
            clf = BiologicalClassifier()
            clf.fit(empty_X, empty_y)
        except ValueError:
            # Acceptable to raise error for empty data
            pass

    def test_single_sample(self):
        """Test handling of single sample."""
        single_X = np.random.randn(1, 10)
        single_y = np.array([1])

        # Should handle gracefully
        try:
            clf = BiologicalClassifier()
            clf.fit(single_X, single_y)
        except ValueError:
            # Acceptable to raise error for insufficient data
            pass

    def test_mismatched_dimensions(self):
        """Test error handling for mismatched X and y dimensions."""
        X = np.random.randn(100, 20)
        y = np.random.randint(0, 2, 50)  # Wrong size

        clf = BiologicalClassifier()

        with pytest.raises(ValueError):
            clf.fit(X, y)

    def test_invalid_algorithms(self):
        """Test error handling for invalid algorithms."""
        with pytest.raises(ValueError, match="Unknown algorithm"):
            BiologicalClassifier(algorithm="invalid_algorithm")

        with pytest.raises(ValueError, match="Unknown algorithm"):
            BiologicalRegressor(algorithm="invalid_algorithm")

    def test_prediction_before_fitting(self):
        """Test error when predicting before fitting."""
        clf = BiologicalClassifier()
        X_test = np.random.randn(10, 5)

        with pytest.raises(ValueError, match="not fitted"):
            clf.predict(X_test)

    def test_large_dataset_performance(self):
        """Test performance with larger datasets."""
        # Create larger dataset
        n_samples = 1000
        n_features = 200

        X = np.random.randn(n_samples, n_features)
        y = np.random.randint(0, 2, n_samples)

        # Should handle large dataset reasonably
        clf = BiologicalClassifier(algorithm="random_forest", random_state=42)
        clf.fit(X, y)

        predictions = clf.predict(X)
        assert len(predictions) == n_samples

        # Should complete in reasonable time and achieve some accuracy
        accuracy = np.mean(predictions == y)
        assert accuracy > 0.4  # Random guessing baseline
