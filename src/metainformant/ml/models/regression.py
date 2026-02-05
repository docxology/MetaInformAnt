"""Machine learning regression utilities for METAINFORMANT.

This module provides regression models and evaluation tools specifically
designed for biological trait prediction and quantitative analysis.
"""

from __future__ import annotations

import numpy as np
from typing import Any, Dict, List, Optional, Tuple, Union

from metainformant.core.utils.logging import get_logger

logger = get_logger(__name__)

# Optional imports for ML functionality
try:
    from sklearn.base import BaseEstimator
    from sklearn.ensemble import RandomForestRegressor, GradientBoostingRegressor
    from sklearn.linear_model import LinearRegression, Ridge, Lasso, ElasticNet
    from sklearn.svm import SVR
    from sklearn.model_selection import cross_val_score, KFold
    from sklearn.metrics import (
        mean_squared_error,
        mean_absolute_error,
        r2_score,
        explained_variance_score,
        median_absolute_error,
    )

    HAS_SKLEARN = True
except ImportError:
    HAS_SKLEARN = False
    logger.warning("scikit-learn not available, ML regression disabled")


class BiologicalRegressor:
    """Wrapper for biological data regression with evaluation metrics.

    Can be initialized in two ways:
    1. Algorithm-based: BiologicalRegressor(algorithm="linear", random_state=42)
    2. Model-based: BiologicalRegressor(model=sklearn_model, model_type="rf")
    """

    _ALGORITHM_MAP = {
        "linear": lambda **kw: (
            LinearRegression(**{k: v for k, v in kw.items() if k != "random_state"}) if HAS_SKLEARN else None
        ),
        "random_forest": lambda **kw: RandomForestRegressor(**kw) if HAS_SKLEARN else None,
        "gradient_boosting": lambda **kw: GradientBoostingRegressor(**kw) if HAS_SKLEARN else None,
        "ridge": lambda **kw: Ridge(**kw) if HAS_SKLEARN else None,
        "lasso": lambda **kw: Lasso(**kw) if HAS_SKLEARN else None,
        "svm": lambda **kw: SVR(**{k: v for k, v in kw.items() if k != "random_state"}) if HAS_SKLEARN else None,
    }

    def __init__(
        self,
        model: Any = None,
        model_type: str = "unknown",
        algorithm: str | None = None,
        random_state: int | None = None,
        **params: Any,
    ):
        """Initialize biological regressor.

        Args:
            model: Pre-built sklearn model (model-based init)
            model_type: Type of model when using model-based init
            algorithm: Algorithm name for algorithm-based init
                       ("linear", "random_forest", "gradient_boosting", "ridge", "lasso", "svm")
            random_state: Random state for reproducibility
            **params: Additional parameters for the sklearn model
        """
        self.random_state = random_state
        self.params = params
        self.feature_names: Optional[List[str]] = None
        self._is_fitted = False

        if algorithm is not None:
            self.algorithm = algorithm
            self.model_type = algorithm
            if not HAS_SKLEARN:
                raise ImportError("scikit-learn required for regression")
            factory = self._ALGORITHM_MAP.get(algorithm)
            if factory is None:
                raise ValueError(f"Unknown algorithm: {algorithm}")
            model_params = {**params}
            if random_state is not None:
                model_params["random_state"] = random_state
            self.model = factory(**model_params)
        elif model is not None:
            self.algorithm = model_type
            self.model_type = model_type
            self.model = model
        else:
            # Default: linear regression
            self.algorithm = "linear"
            self.model_type = "linear"
            if HAS_SKLEARN:
                self.model = LinearRegression()
            else:
                self.model = None

    @property
    def is_fitted(self) -> bool:
        """Whether the model has been fitted."""
        return self._is_fitted

    @property
    def is_trained(self) -> bool:
        """Alias for is_fitted."""
        return self._is_fitted

    @is_trained.setter
    def is_trained(self, value: bool) -> None:
        """Setter for backward compatibility."""
        self._is_fitted = value

    def fit(self, X: np.ndarray, y: np.ndarray, feature_names: Optional[List[str]] = None) -> BiologicalRegressor:
        """Fit the regressor.

        Args:
            X: Feature matrix
            y: Target values
            feature_names: Optional feature names

        Returns:
            Self for chaining
        """
        if not HAS_SKLEARN:
            raise ImportError("scikit-learn required for regression")

        self.model.fit(X, y)
        self.feature_names = feature_names
        self._is_fitted = True

        logger.info(f"Trained {self.model_type} regressor with {X.shape[1]} features")
        return self

    def predict(self, X: np.ndarray) -> np.ndarray:
        """Make predictions.

        Args:
            X: Feature matrix

        Returns:
            Predicted values
        """
        if not self._is_fitted:
            raise ValueError("Model not fitted")
        return self.model.predict(X)

    def evaluate(self, X: np.ndarray, y: np.ndarray, detailed: bool = True) -> Dict[str, Any]:
        """Evaluate regressor performance.

        Args:
            X: Feature matrix
            y: True target values
            detailed: Whether to include detailed metrics

        Returns:
            Dictionary with evaluation metrics
        """
        if not self._is_fitted:
            raise ValueError("Model not fitted")

        y_pred = self.predict(X)

        results = {
            "r2": r2_score(y, y_pred),
            "mse": mean_squared_error(y, y_pred),
            "mae": mean_absolute_error(y, y_pred),
            "explained_variance": explained_variance_score(y, y_pred),
            "median_absolute_error": median_absolute_error(y, y_pred),
            "rmse": np.sqrt(mean_squared_error(y, y_pred)),
        }

        # Additional metrics
        if detailed:
            # Residual analysis
            residuals = y - y_pred
            results["residuals"] = {
                "mean": float(residuals.mean()),
                "std": float(residuals.std()),
                "min": float(residuals.min()),
                "max": float(residuals.max()),
            }

            # Prediction intervals (simplified)
            std_residuals = np.std(residuals)
            results["prediction_std"] = float(std_residuals)

            # Feature importance if available
            if hasattr(self.model, "feature_importances_"):
                results["feature_importances"] = self.model.feature_importances_.tolist()
            elif hasattr(self.model, "coef_"):
                results["coefficients"] = self.model.coef_.tolist()

        return results


def train_regressor(X: np.ndarray, y: np.ndarray, method: str = "rf", **kwargs: Any) -> BiologicalRegressor:
    """Train a regression model for biological traits.

    Args:
        X: Training feature matrix
        y: Training target values
        method: Regression method ('rf', 'gb', 'linear', 'ridge', 'lasso', 'svr')
        **kwargs: Additional parameters for the regressor

    Returns:
        Trained BiologicalRegressor instance

    Raises:
        ImportError: If scikit-learn not available
        ValueError: If invalid method specified
    """
    if not HAS_SKLEARN:
        raise ImportError("scikit-learn required for regression")

    # Select and configure regressor
    if method == "rf":
        model = RandomForestRegressor(
            n_estimators=kwargs.get("n_estimators", 100), random_state=kwargs.get("random_state", 42), **kwargs
        )
    elif method == "gb":
        model = GradientBoostingRegressor(
            n_estimators=kwargs.get("n_estimators", 100), random_state=kwargs.get("random_state", 42), **kwargs
        )
    elif method == "linear":
        model = LinearRegression(**kwargs)
    elif method == "ridge":
        model = Ridge(alpha=kwargs.get("alpha", 1.0), random_state=kwargs.get("random_state", 42), **kwargs)
    elif method == "lasso":
        model = Lasso(alpha=kwargs.get("alpha", 1.0), random_state=kwargs.get("random_state", 42), **kwargs)
    elif method == "elasticnet":
        model = ElasticNet(
            alpha=kwargs.get("alpha", 1.0),
            l1_ratio=kwargs.get("l1_ratio", 0.5),
            random_state=kwargs.get("random_state", 42),
            **kwargs,
        )
    elif method == "svr":
        model = SVR(kernel=kwargs.get("kernel", "rbf"), C=kwargs.get("C", 1.0), **kwargs)
    else:
        raise ValueError(f"Unknown regression method: {method}")

    # Train model
    regressor = BiologicalRegressor(model, method)
    regressor.fit(X, y)

    logger.info(f"Trained {method} regressor")
    return regressor


def evaluate_regressor(regressor: "BiologicalRegressor", X: np.ndarray, y: np.ndarray, **kwargs: Any) -> Dict[str, Any]:
    """Evaluate a trained regressor on data.

    If the regressor is already fitted, evaluates on the provided data directly.
    Otherwise fits on a train split and evaluates on a test split.

    Args:
        regressor: BiologicalRegressor (fitted or unfitted)
        X: Feature matrix
        y: Target values
        **kwargs: Optional test_size and random_state for unfitted models

    Returns:
        Dictionary with evaluation metrics and predictions
    """
    from sklearn.metrics import mean_squared_error, r2_score, mean_absolute_error

    if not regressor.is_fitted:
        from sklearn.model_selection import train_test_split

        test_size = kwargs.get("test_size", 0.2)
        random_state = kwargs.get("random_state", 42)
        X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=test_size, random_state=random_state)
        regressor.fit(X_train, y_train)
        y_pred = regressor.predict(X_test)
        y_eval = y_test
    else:
        y_pred = regressor.predict(X)
        y_eval = y

    mse = mean_squared_error(y_eval, y_pred)
    r2 = r2_score(y_eval, y_pred)
    mae = mean_absolute_error(y_eval, y_pred)

    results = {
        "mse": mse,
        "r2_score": r2,
        "mae": mae,
        "predictions": y_pred.tolist(),
    }

    logger.info(f"Regressor evaluation: MSE={mse:.4f}, R²={r2:.4f}")
    return results


def cross_validate_regressor(model: BaseEstimator, X: np.ndarray, y: np.ndarray, cv: int = 5) -> Dict[str, float]:
    """Cross-validate a regression model.

    Args:
        model: Trained sklearn regressor
        X: Feature matrix
        y: Target values
        cv: Number of cross-validation folds

    Returns:
        Dictionary with cross-validation metrics

    Raises:
        ImportError: If scikit-learn not available
    """
    if not HAS_SKLEARN:
        raise ImportError("scikit-learn required for cross-validation")

    # Define scoring metrics
    scoring = {
        "r2": "r2",
        "neg_mean_squared_error": "neg_mean_squared_error",
        "neg_mean_absolute_error": "neg_mean_absolute_error",
    }

    results = {}

    for metric_name, scoring_name in scoring.items():
        scores = cross_val_score(
            model, X, y, cv=KFold(n_splits=cv, shuffle=True, random_state=42), scoring=scoring_name
        )

        if "neg_" in metric_name:
            # Convert negative scores to positive
            scores = -scores
            metric_name = metric_name.replace("neg_", "")

        results[f"{metric_name}_mean"] = float(scores.mean())
        results[f"{metric_name}_std"] = float(scores.std())
        results[f"{metric_name}_scores"] = scores.tolist()

    logger.info(f"Cross-validated regressor: R² = {results['r2_mean']:.3f} ± {results['r2_std']:.3f}")
    return results


def create_ensemble_regressor(
    X: np.ndarray, y: np.ndarray, n_estimators: int = 10, random_state: int | None = None, **kwargs: Any
) -> BiologicalRegressor:
    """Create an ensemble regressor for improved prediction.

    Args:
        X: Training feature matrix
        y: Training target values
        n_estimators: Number of base estimators
        random_state: Random state for reproducibility
        **kwargs: Additional parameters

    Returns:
        Trained ensemble BiologicalRegressor

    Raises:
        ImportError: If scikit-learn not available
    """
    if not HAS_SKLEARN:
        raise ImportError("scikit-learn required for ensemble regression")

    # Create multiple regressors
    estimators = [
        (
            "rf",
            RandomForestRegressor(n_estimators=n_estimators, random_state=random_state, **kwargs.get("rf_params", {})),
        ),
        (
            "gb",
            GradientBoostingRegressor(
                n_estimators=n_estimators, random_state=random_state, **kwargs.get("gb_params", {})
            ),
        ),
        ("ridge", Ridge(random_state=random_state, **kwargs.get("ridge_params", {}))),
    ]

    # Simple averaging ensemble (could be improved with more sophisticated methods)
    class AveragingRegressor:
        def __init__(self, estimators):
            self.estimators = estimators
            self.is_fitted = False

        def fit(self, X, y):
            for name, estimator in self.estimators:
                estimator.fit(X, y)
            self.is_fitted = True
            return self

        def predict(self, X):
            if not self.is_fitted:
                raise ValueError("Model not fitted")
            predictions = [estimator.predict(X) for _, estimator in self.estimators]
            return np.mean(predictions, axis=0)

    ensemble = AveragingRegressor(estimators)
    regressor = BiologicalRegressor(ensemble, "ensemble")
    regressor.fit(X, y)

    logger.info(f"Trained ensemble regressor with {len(estimators)} base models")
    return regressor


def compare_regression_methods(
    X: np.ndarray,
    y: np.ndarray,
    methods: List[str] = None,
    cv_folds: int = 5,
    random_state: int | None = None,
) -> Dict[str, Any]:
    """Compare multiple regression methods.

    Args:
        X: Feature matrix
        y: Target values
        methods: List of methods to compare
        cv_folds: Number of cross-validation folds
        random_state: Random state for reproducibility

    Returns:
        Dictionary with comparison results
    """
    if methods is None:
        methods = ["rf", "gb", "ridge", "lasso"]

    results = {}

    for method in methods:
        try:
            # Train model
            regressor = train_regressor(X, y, method=method, random_state=random_state)

            # Cross-validate
            cv_results = cross_validate_regressor(regressor.model, X, y, cv=cv_folds)

            results[method] = {
                "regressor": regressor,
                "cv_results": cv_results,
                "r2_mean": cv_results["r2_mean"],
                "mse_mean": cv_results["mean_squared_error_mean"],
            }

        except Exception as e:
            logger.error(f"Failed to evaluate {method}: {e}")
            results[method] = {"error": str(e)}

    # Find best method
    valid_results = [
        (method, result["r2_mean"])
        for method, result in results.items()
        if isinstance(result, dict) and "r2_mean" in result
    ]

    if valid_results:
        best_method = max(valid_results, key=lambda x: x[1])[0]
    else:
        best_method = None

    return {
        "comparison": results,
        "best_method": best_method,
        "methods_tested": methods,
    }


def predict_phenotypic_traits(
    X_genetic: np.ndarray, y_phenotype: np.ndarray, X_predict: np.ndarray, method: str = "rf", **kwargs: Any
) -> Dict[str, Any]:
    """Predict phenotypic traits from genetic data.

    Args:
        X_genetic: Genetic feature matrix for training
        y_phenotype: Observed phenotypic values
        X_predict: Genetic features for prediction
        method: Regression method to use
        **kwargs: Additional parameters

    Returns:
        Dictionary with predictions and model evaluation

    Raises:
        ImportError: If scikit-learn not available
    """
    if not HAS_SKLEARN:
        raise ImportError("scikit-learn required for phenotypic prediction")

    # Train model
    regressor = train_regressor(X_genetic, y_phenotype, method=method, **kwargs)

    # Evaluate on training data
    train_eval = regressor.evaluate(X_genetic, y_phenotype)

    # Make predictions
    predictions = regressor.predict(X_predict)

    # Cross-validate
    cv_results = cross_validate_regressor(regressor.model, X_genetic, y_phenotype)

    return {
        "method": method,
        "regressor": regressor,
        "training_evaluation": train_eval,
        "cross_validation": cv_results,
        "predictions": predictions.tolist(),
        "prediction_stats": {
            "mean": float(predictions.mean()),
            "std": float(predictions.std()),
            "min": float(predictions.min()),
            "max": float(predictions.max()),
        },
    }


def analyze_prediction_uncertainty(
    X: np.ndarray,
    y: np.ndarray,
    X_test: np.ndarray,
    method: str = "rf",
    n_bootstraps: int = 100,
    random_state: int | None = None,
) -> Dict[str, Any]:
    """Analyze prediction uncertainty using bootstrapping.

    Args:
        X: Training feature matrix
        y: Training target values
        X_test: Test feature matrix for uncertainty analysis
        method: Regression method
        n_bootstraps: Number of bootstrap iterations
        random_state: Random state for reproducibility

    Returns:
        Dictionary with uncertainty analysis results

    Raises:
        ImportError: If scikit-learn not available
    """
    if not HAS_SKLEARN:
        raise ImportError("scikit-learn required for uncertainty analysis")

    np.random.seed(random_state)
    n_samples = len(X)

    # Store predictions from each bootstrap
    all_predictions = []

    for i in range(n_bootstraps):
        # Bootstrap sample
        indices = np.random.choice(n_samples, size=n_samples, replace=True)
        X_boot = X[indices]
        y_boot = y[indices]

        # Train model
        regressor = train_regressor(X_boot, y_boot, method=method, random_state=random_state)
        predictions = regressor.predict(X_test)
        all_predictions.append(predictions)

    # Convert to array
    all_predictions = np.array(all_predictions)

    # Calculate uncertainty metrics
    mean_predictions = np.mean(all_predictions, axis=0)
    std_predictions = np.std(all_predictions, axis=0)

    # Confidence intervals (95%)
    ci_lower = np.percentile(all_predictions, 2.5, axis=0)
    ci_upper = np.percentile(all_predictions, 97.5, axis=0)

    return {
        "n_bootstraps": n_bootstraps,
        "mean_predictions": mean_predictions.tolist(),
        "std_predictions": std_predictions.tolist(),
        "ci_lower_95": ci_lower.tolist(),
        "ci_upper_95": ci_upper.tolist(),
        "prediction_ranges": (ci_upper - ci_lower).tolist(),
    }
