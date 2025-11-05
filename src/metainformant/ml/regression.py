"""Regression methods for biological data analysis."""

from __future__ import annotations

import warnings
from typing import Any, Dict, List, Optional

import numpy as np


class BiologicalRegressor:
    """Regression wrapper for biological data analysis."""

    def __init__(self, algorithm: str = "linear", random_state: Optional[int] = None, **kwargs):
        """Initialize biological regressor.

        Args:
            algorithm: Regression algorithm to use
            random_state: Random seed for reproducibility
            **kwargs: Algorithm-specific parameters
        """
        # Validate algorithm
        supported_algorithms = ["linear", "ridge", "lasso", "random_forest"]
        if algorithm not in supported_algorithms:
            raise ValueError(f"Unknown algorithm: {algorithm}. Supported: {supported_algorithms}")
        
        self.algorithm = algorithm
        self.random_state = random_state
        self.params = kwargs
        self.is_fitted = False
        self.feature_importance_ = None

    def fit(self, X: np.ndarray, y: np.ndarray) -> "BiologicalRegressor":
        """Fit regressor to training data."""
        if self.random_state is not None:
            np.random.seed(self.random_state)

        self.X_train_ = X.copy()
        self.y_train_ = y.copy()

        # Calculate coefficients/importance
        if self.algorithm == "linear":
            self._fit_linear(X, y)
        elif self.algorithm == "ridge":
            self._fit_ridge(X, y)
        else:
            # Default linear
            self._fit_linear(X, y)

        self.is_fitted = True
        return self

    def predict(self, X: np.ndarray) -> np.ndarray:
        """Predict continuous values for samples."""
        if not self.is_fitted:
            raise ValueError("Regressor must be fitted before prediction")

        if hasattr(self, "coefficients_"):
            return X @ self.coefficients_ + self.intercept_
        else:
            # Fallback: mean prediction
            return np.full(X.shape[0], np.mean(self.y_train_))

    def _fit_linear(self, X: np.ndarray, y: np.ndarray) -> None:
        """Fit linear regression."""
        try:
            # Add bias term
            X_with_bias = np.column_stack([np.ones(X.shape[0]), X])

            # Normal equation with regularization for stability
            XtX = X_with_bias.T @ X_with_bias
            Xty = X_with_bias.T @ y

            # Add small regularization
            reg = 1e-6 * np.eye(X_with_bias.shape[1])
            coeffs = np.linalg.solve(XtX + reg, Xty)

            self.intercept_ = coeffs[0]
            self.coefficients_ = coeffs[1:]
            self.feature_importance_ = np.abs(self.coefficients_)

        except np.linalg.LinAlgError:
            # Fallback if matrix is singular
            self.intercept_ = np.mean(y)
            self.coefficients_ = np.zeros(X.shape[1])
            self.feature_importance_ = np.ones(X.shape[1]) / X.shape[1]

    def _fit_ridge(self, X: np.ndarray, y: np.ndarray) -> None:
        """Fit ridge regression."""
        alpha = self.params.get("alpha", 1.0)

        try:
            # Ridge regression with regularization
            XtX = X.T @ X
            Xty = X.T @ y

            # Add ridge penalty
            ridge_matrix = XtX + alpha * np.eye(X.shape[1])
            self.coefficients_ = np.linalg.solve(ridge_matrix, Xty)
            self.intercept_ = np.mean(y) - np.mean(X, axis=0) @ self.coefficients_
            self.feature_importance_ = np.abs(self.coefficients_)

        except np.linalg.LinAlgError:
            self._fit_linear(X, y)


def train_ensemble_regressor(
    X: np.ndarray, y: np.ndarray, algorithms: List[str] = None, random_state: Optional[int] = None
) -> Dict[str, BiologicalRegressor]:
    """Train ensemble of regressors."""
    if algorithms is None:
        algorithms = ["linear", "ridge"]

    ensemble = {}

    for algorithm in algorithms:
        try:
            regressor = BiologicalRegressor(algorithm=algorithm, random_state=random_state)
            regressor.fit(X, y)
            ensemble[algorithm] = regressor
        except Exception as e:
            warnings.warn(f"Failed to train {algorithm}: {e}")

    return ensemble


def evaluate_regressor(regressor: BiologicalRegressor, X_test: np.ndarray = None, y_test: np.ndarray = None, X: np.ndarray = None, y: np.ndarray = None) -> Dict[str, float]:
    """Evaluate regressor performance.

    Args:
        regressor: Fitted regressor
        X_test: Test features (can also use X as keyword)
        y_test: True test labels (can also use y as keyword)

    Returns:
        Dictionary of evaluation metrics
    """
    # Support both X_test/y_test and X/y parameter names for compatibility
    if X_test is None and X is not None:
        X_test = X
    if y_test is None and y is not None:
        y_test = y
    
    if X_test is None or y_test is None:
        raise ValueError("Must provide X_test/X and y_test/y parameters")
    
    predictions = regressor.predict(X_test)

    # Calculate metrics
    mse = np.mean((predictions - y_test) ** 2)
    mae = np.mean(np.abs(predictions - y_test))

    # R-squared
    ss_res = np.sum((y_test - predictions) ** 2)
    ss_tot = np.sum((y_test - np.mean(y_test)) ** 2)
    r2 = 1 - (ss_res / (ss_tot + 1e-10))

    return {"mse": mse, "rmse": np.sqrt(mse), "mae": mae, "r2": r2}
