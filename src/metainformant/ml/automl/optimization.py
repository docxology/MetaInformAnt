"""Automated hyperparameter optimization and model selection.

This module provides automated machine learning (AutoML) capabilities
including random search, Bayesian optimization with a Gaussian process
surrogate, exhaustive grid search, automatic model selection, and
automatic preprocessing pipelines. All methods include pure Python
fallbacks where possible.
"""

from __future__ import annotations

import math
import random
from collections import defaultdict
from itertools import product
from typing import Any, Callable

from metainformant.core.utils.logging import get_logger

logger = get_logger(__name__)

# Optional dependencies
try:
    import numpy as np

    HAS_NUMPY = True
except ImportError:
    HAS_NUMPY = False
    np = None  # type: ignore[assignment]

try:
    from sklearn.model_selection import cross_val_score
    from sklearn.preprocessing import LabelEncoder, StandardScaler

    HAS_SKLEARN = True
except ImportError:
    HAS_SKLEARN = False


def _to_2d_list(X: Any) -> list[list[float]]:
    """Convert input matrix to list of lists.

    Args:
        X: Input matrix.

    Returns:
        Matrix as list of lists.
    """
    if HAS_NUMPY and isinstance(X, np.ndarray):
        return [[float(X[i, j]) for j in range(X.shape[1])] for i in range(X.shape[0])]
    return [[float(v) for v in row] for row in X]


def _to_1d_list(y: Any) -> list[float]:
    """Convert input vector to list of floats.

    Args:
        y: Input vector.

    Returns:
        Vector as list of floats.
    """
    if HAS_NUMPY and isinstance(y, np.ndarray):
        return [float(v) for v in y.ravel()]
    return [float(v) for v in y]


def _get_shape(X: Any) -> tuple[int, int]:
    """Get shape of 2D matrix.

    Args:
        X: Input matrix.

    Returns:
        Tuple of (n_rows, n_cols).
    """
    if HAS_NUMPY and isinstance(X, np.ndarray):
        return int(X.shape[0]), int(X.shape[1])
    n_rows = len(X)
    n_cols = len(X[0]) if n_rows > 0 else 0
    return n_rows, n_cols


def _cross_validate(
    model: Any,
    X_list: list[list[float]],
    y_list: list[float],
    cv: int = 5,
    metric: str = "accuracy",
) -> float:
    """Cross-validate a model and return mean score.

    Args:
        model: Model with fit/predict interface.
        X_list: Feature matrix.
        y_list: Target values.
        cv: Number of folds.
        metric: Scoring metric.

    Returns:
        Mean CV score.
    """
    if HAS_SKLEARN and HAS_NUMPY:
        try:
            scoring = metric if metric != "accuracy" else "accuracy"
            scores = cross_val_score(model, np.array(X_list), np.array(y_list), cv=cv, scoring=scoring)
            return float(scores.mean())
        except Exception:
            pass

    # Pure Python fallback
    n = len(y_list)
    fold_size = n // cv
    scores: list[float] = []

    for fold in range(cv):
        start = fold * fold_size
        end = start + fold_size if fold < cv - 1 else n
        test_idx = set(range(start, end))

        train_X = [X_list[i] for i in range(n) if i not in test_idx]
        train_y = [y_list[i] for i in range(n) if i not in test_idx]
        test_X = [X_list[i] for i in range(start, end)]
        test_y = [y_list[i] for i in range(start, end)]

        if not train_X or not test_X:
            continue

        try:
            if HAS_NUMPY:
                model.fit(np.array(train_X), np.array(train_y))
                preds = model.predict(np.array(test_X))
            else:
                model.fit(train_X, train_y)
                preds = model.predict(test_X)

            preds_list = [float(p) for p in preds]

            if metric == "accuracy":
                correct = sum(1 for yt, yp in zip(test_y, preds_list) if round(yt) == round(yp))
                scores.append(correct / len(test_y))
            elif metric == "mse":
                mse = sum((yt - yp) ** 2 for yt, yp in zip(test_y, preds_list)) / len(test_y)
                scores.append(-mse)
            elif metric == "r2":
                mean_y = sum(test_y) / len(test_y)
                ss_tot = sum((yt - mean_y) ** 2 for yt in test_y)
                ss_res = sum((yt - yp) ** 2 for yt, yp in zip(test_y, preds_list))
                scores.append(1.0 - ss_res / ss_tot if ss_tot > 1e-15 else 0.0)
            else:
                correct = sum(1 for yt, yp in zip(test_y, preds_list) if round(yt) == round(yp))
                scores.append(correct / len(test_y))
        except Exception:
            scores.append(0.0)

    return sum(scores) / len(scores) if scores else 0.0


def _sample_from_distribution(name: str, dist: Any) -> Any:
    """Sample a value from a parameter distribution specification.

    Args:
        name: Parameter name.
        dist: Distribution specification. Can be:
            - list/tuple: uniform choice from values
            - dict with 'low', 'high': uniform float in range
            - dict with 'low', 'high', 'log': log-uniform float
            - dict with 'low', 'high', 'type': 'int': uniform int

    Returns:
        Sampled value.
    """
    if isinstance(dist, (list, tuple)):
        return random.choice(dist)

    if isinstance(dist, dict):
        low = dist.get("low", 0)
        high = dist.get("high", 1)

        if dist.get("log", False):
            log_low = math.log(max(low, 1e-15))
            log_high = math.log(max(high, 1e-15))
            return math.exp(random.uniform(log_low, log_high))

        if dist.get("type") == "int":
            return random.randint(int(low), int(high))

        return random.uniform(low, high)

    return dist


def random_search(
    model_fn: Any,
    param_distributions: dict,
    X: Any,
    y: Any,
    n_iter: int = 50,
    cv: int = 5,
    metric: str = "accuracy",
    random_state: int | None = None,
) -> dict:
    """Random search hyperparameter optimization with cross-validation.

    Samples random parameter combinations from specified distributions and
    evaluates each via cross-validation. More efficient than grid search
    for high-dimensional parameter spaces.

    Args:
        model_fn: Callable that takes keyword arguments and returns a model
            instance with fit/predict methods.
        param_distributions: Dict mapping parameter names to distributions.
            Each value can be a list (uniform choice), or a dict with
            'low'/'high' keys (uniform range), optionally 'log': True
            for log-uniform sampling.
        X: Feature matrix (n_samples x n_features).
        y: Target values (n_samples,).
        n_iter: Number of parameter combinations to try.
        cv: Number of cross-validation folds.
        metric: Scoring metric ('accuracy', 'mse', 'r2').
        random_state: Random seed for reproducibility.

    Returns:
        Dictionary containing:
            - best_params: Dict of best hyperparameters found.
            - best_score: Best cross-validation score.
            - all_results: List of dicts with params and score for each
              iteration.
            - n_evaluations: Total number of evaluations performed.

    Raises:
        ValueError: If param_distributions is empty.
    """
    if not param_distributions:
        raise ValueError("param_distributions must not be empty")

    if random_state is not None:
        random.seed(random_state)

    X_list = _to_2d_list(X)
    y_list = _to_1d_list(y)

    logger.info(
        "Random search: %d iterations, %d params, cv=%d, metric=%s",
        n_iter,
        len(param_distributions),
        cv,
        metric,
    )

    all_results: list[dict] = []
    best_score = float("-inf")
    best_params: dict = {}

    for i in range(n_iter):
        # Sample parameters
        params = {name: _sample_from_distribution(name, dist) for name, dist in param_distributions.items()}

        try:
            model = model_fn(**params)
            score = _cross_validate(model, X_list, y_list, cv=cv, metric=metric)

            all_results.append(
                {
                    "params": params,
                    "score": round(score, 6),
                    "iteration": i,
                }
            )

            if score > best_score:
                best_score = score
                best_params = params.copy()

        except Exception as exc:
            logger.warning("Iteration %d failed: %s", i, exc)
            all_results.append(
                {
                    "params": params,
                    "score": float("nan"),
                    "iteration": i,
                    "error": str(exc),
                }
            )

    logger.info("Random search complete: best_score=%.4f", best_score)

    return {
        "best_params": best_params,
        "best_score": round(best_score, 6),
        "all_results": all_results,
        "n_evaluations": len(all_results),
    }


def bayesian_optimization(
    objective_fn: Any,
    param_space: dict,
    n_iter: int = 30,
    n_initial: int = 5,
    random_state: int | None = None,
) -> dict:
    """Bayesian optimization using a Gaussian process surrogate.

    Uses a GP with RBF kernel as a surrogate model to efficiently search
    the parameter space. Balances exploration and exploitation using the
    Expected Improvement (EI) acquisition function.

    The algorithm:
    1. Evaluate n_initial random points to build initial GP model.
    2. For each subsequent iteration, fit GP to observed data.
    3. Maximize EI to select the next evaluation point.
    4. Evaluate the objective at the selected point.

    Args:
        objective_fn: Callable that takes a dict of parameters and returns
            a scalar score (higher is better).
        param_space: Dict mapping parameter names to bounds. Each value
            should be a dict with 'low' and 'high' keys, optionally
            'log': True for log-scale.
        n_iter: Total number of evaluations.
        n_initial: Number of initial random evaluations.
        random_state: Random seed for reproducibility.

    Returns:
        Dictionary containing:
            - best_params: Dict of best parameters found.
            - best_score: Best objective value.
            - history: List of dicts with params and score per iteration.
            - surrogate_model: Dict describing the GP state
              (for inspection/debugging).

    Raises:
        ValueError: If param_space is empty or n_initial > n_iter.
    """
    if not param_space:
        raise ValueError("param_space must not be empty")
    if n_initial > n_iter:
        raise ValueError(f"n_initial ({n_initial}) must be <= n_iter ({n_iter})")

    if random_state is not None:
        random.seed(random_state)

    logger.info(
        "Bayesian optimization: %d iterations, %d initial, %d params",
        n_iter,
        n_initial,
        len(param_space),
    )

    param_names = sorted(param_space.keys())
    n_params = len(param_names)

    # Convert param space to normalized bounds [0, 1]
    bounds_info: list[dict] = []
    for name in param_names:
        spec = param_space[name]
        if isinstance(spec, dict):
            bounds_info.append(
                {
                    "low": spec.get("low", 0),
                    "high": spec.get("high", 1),
                    "log": spec.get("log", False),
                }
            )
        else:
            bounds_info.append({"low": 0, "high": 1, "log": False})

    def _normalize(params_dict: dict) -> list[float]:
        normalized = []
        for i, name in enumerate(param_names):
            val = params_dict[name]
            info = bounds_info[i]
            if info["log"]:
                log_low = math.log(max(info["low"], 1e-15))
                log_high = math.log(max(info["high"], 1e-15))
                norm = (math.log(max(val, 1e-15)) - log_low) / max(log_high - log_low, 1e-15)
            else:
                norm = (val - info["low"]) / max(info["high"] - info["low"], 1e-15)
            normalized.append(max(0.0, min(1.0, norm)))
        return normalized

    def _denormalize(normalized: list[float]) -> dict:
        params = {}
        for i, name in enumerate(param_names):
            info = bounds_info[i]
            norm = max(0.0, min(1.0, normalized[i]))
            if info["log"]:
                log_low = math.log(max(info["low"], 1e-15))
                log_high = math.log(max(info["high"], 1e-15))
                params[name] = math.exp(log_low + norm * (log_high - log_low))
            else:
                params[name] = info["low"] + norm * (info["high"] - info["low"])
        return params

    history: list[dict] = []
    X_observed: list[list[float]] = []
    y_observed: list[float] = []
    best_score = float("-inf")
    best_params: dict = {}

    # Phase 1: Initial random evaluations
    for i in range(n_initial):
        norm_point = [random.random() for _ in range(n_params)]
        params = _denormalize(norm_point)

        try:
            score = float(objective_fn(params))
        except Exception as exc:
            logger.warning("Initial evaluation %d failed: %s", i, exc)
            score = float("-inf")

        X_observed.append(norm_point)
        y_observed.append(score)
        history.append({"params": params, "score": round(score, 6), "iteration": i})

        if score > best_score:
            best_score = score
            best_params = params.copy()

    # Phase 2: Bayesian optimization loop
    for i in range(n_initial, n_iter):
        # Fit GP surrogate and find point with maximum EI
        next_point = _maximize_expected_improvement(X_observed, y_observed, n_params, n_candidates=200)

        params = _denormalize(next_point)

        try:
            score = float(objective_fn(params))
        except Exception as exc:
            logger.warning("Iteration %d failed: %s", i, exc)
            score = float("-inf")

        X_observed.append(next_point)
        y_observed.append(score)
        history.append({"params": params, "score": round(score, 6), "iteration": i})

        if score > best_score:
            best_score = score
            best_params = params.copy()

    # Build surrogate model summary
    surrogate_model = {
        "type": "gaussian_process",
        "kernel": "rbf",
        "n_observations": len(X_observed),
        "y_mean": sum(y_observed) / len(y_observed) if y_observed else 0.0,
        "y_std": (
            math.sqrt(
                sum((y - sum(y_observed) / len(y_observed)) ** 2 for y in y_observed) / max(len(y_observed) - 1, 1)
            )
            if y_observed
            else 0.0
        ),
    }

    logger.info("Bayesian optimization complete: best_score=%.4f", best_score)

    return {
        "best_params": best_params,
        "best_score": round(best_score, 6),
        "history": history,
        "surrogate_model": surrogate_model,
    }


def _maximize_expected_improvement(
    X_obs: list[list[float]],
    y_obs: list[float],
    n_dims: int,
    n_candidates: int = 200,
    length_scale: float = 0.3,
) -> list[float]:
    """Find the point maximizing Expected Improvement using a GP surrogate.

    Uses a simple RBF kernel GP and evaluates EI on random candidate points.

    Args:
        X_obs: Observed points (normalized [0,1]).
        y_obs: Observed values.
        n_dims: Dimensionality.
        n_candidates: Number of random candidates to evaluate.
        length_scale: RBF kernel length scale.

    Returns:
        Point (normalized) with maximum EI.
    """
    if not X_obs:
        return [random.random() for _ in range(n_dims)]

    n = len(X_obs)
    best_y = max(y_obs)
    mean_y = sum(y_obs) / n

    # Compute kernel matrix K(X_obs, X_obs)
    K = [[0.0] * n for _ in range(n)]
    for i in range(n):
        for j in range(i, n):
            sq_dist = sum((X_obs[i][d] - X_obs[j][d]) ** 2 for d in range(n_dims))
            k_val = math.exp(-sq_dist / (2 * length_scale**2))
            K[i][j] = k_val
            K[j][i] = k_val
        K[i][i] += 1e-6  # Noise/regularization

    # Compute K_inv using Cholesky-like approach (or direct inversion for small n)
    K_inv = _invert_matrix(K)
    if K_inv is None:
        return [random.random() for _ in range(n_dims)]

    # K_inv @ y
    alpha = [sum(K_inv[i][j] * y_obs[j] for j in range(n)) for i in range(n)]

    best_ei = -1.0
    best_candidate = [random.random() for _ in range(n_dims)]

    for _ in range(n_candidates):
        candidate = [random.random() for _ in range(n_dims)]

        # Compute k(candidate, X_obs)
        k_star = [
            math.exp(-sum((candidate[d] - X_obs[j][d]) ** 2 for d in range(n_dims)) / (2 * length_scale**2))
            for j in range(n)
        ]

        # GP mean at candidate
        mu = sum(k_star[j] * alpha[j] for j in range(n))

        # GP variance at candidate
        k_ss = 1.0  # k(candidate, candidate) = 1 for RBF
        v = [sum(K_inv[i][j] * k_star[j] for j in range(n)) for i in range(n)]
        var = k_ss - sum(k_star[j] * v[j] for j in range(n))
        var = max(var, 1e-10)
        sigma = math.sqrt(var)

        # Expected Improvement
        if sigma < 1e-10:
            ei = 0.0
        else:
            z = (mu - best_y) / sigma
            # EI = (mu - best_y) * Phi(z) + sigma * phi(z)
            phi_z = math.exp(-0.5 * z * z) / math.sqrt(2 * math.pi)
            Phi_z = 0.5 * (1 + math.erf(z / math.sqrt(2)))
            ei = (mu - best_y) * Phi_z + sigma * phi_z

        if ei > best_ei:
            best_ei = ei
            best_candidate = candidate

    return best_candidate


def _invert_matrix(M: list[list[float]]) -> list[list[float]] | None:
    """Invert a square matrix via Gauss-Jordan elimination.

    Args:
        M: Square matrix.

    Returns:
        Inverse matrix, or None if singular.
    """
    n = len(M)
    # Augment with identity
    aug = [M[i][:] + [1.0 if j == i else 0.0 for j in range(n)] for i in range(n)]

    for col in range(n):
        # Pivot
        max_val = abs(aug[col][col])
        max_row = col
        for row in range(col + 1, n):
            if abs(aug[row][col]) > max_val:
                max_val = abs(aug[row][col])
                max_row = row

        if max_val < 1e-15:
            return None

        aug[col], aug[max_row] = aug[max_row], aug[col]

        pivot = aug[col][col]
        for j in range(2 * n):
            aug[col][j] /= pivot

        for row in range(n):
            if row == col:
                continue
            factor = aug[row][col]
            for j in range(2 * n):
                aug[row][j] -= factor * aug[col][j]

    return [aug[i][n:] for i in range(n)]


def grid_search(
    model_fn: Any,
    param_grid: dict,
    X: Any,
    y: Any,
    cv: int = 5,
    metric: str = "accuracy",
) -> dict:
    """Exhaustive grid search over hyperparameter combinations.

    Evaluates every combination in the parameter grid via cross-validation.
    Suitable for small parameter spaces where exhaustive coverage is desired.

    Args:
        model_fn: Callable that takes keyword arguments and returns a model.
        param_grid: Dict mapping parameter names to lists of values to try.
        X: Feature matrix (n_samples x n_features).
        y: Target values (n_samples,).
        cv: Number of cross-validation folds.
        metric: Scoring metric.

    Returns:
        Dictionary containing:
            - best_params: Dict of best hyperparameters found.
            - best_score: Best cross-validation score.
            - all_results: List of dicts with params and score per
              combination.

    Raises:
        ValueError: If param_grid is empty.
    """
    if not param_grid:
        raise ValueError("param_grid must not be empty")

    X_list = _to_2d_list(X)
    y_list = _to_1d_list(y)

    param_names = sorted(param_grid.keys())
    param_values = [param_grid[name] for name in param_names]
    total_combos = 1
    for vals in param_values:
        total_combos *= len(vals)

    logger.info("Grid search: %d combinations, %d params, cv=%d", total_combos, len(param_names), cv)

    all_results: list[dict] = []
    best_score = float("-inf")
    best_params: dict = {}

    for combo in product(*param_values):
        params = {name: val for name, val in zip(param_names, combo)}

        try:
            model = model_fn(**params)
            score = _cross_validate(model, X_list, y_list, cv=cv, metric=metric)

            all_results.append({"params": params, "score": round(score, 6)})

            if score > best_score:
                best_score = score
                best_params = params.copy()
        except Exception as exc:
            logger.warning("Grid search combo failed: %s, error: %s", params, exc)
            all_results.append(
                {
                    "params": params,
                    "score": float("nan"),
                    "error": str(exc),
                }
            )

    logger.info("Grid search complete: best_score=%.4f", best_score)

    return {
        "best_params": best_params,
        "best_score": round(best_score, 6),
        "all_results": all_results,
    }


def model_selection(
    X: Any,
    y: Any,
    task: str = "classification",
    cv: int = 5,
) -> dict:
    """Automatic model selection across multiple model types.

    Tries several model types (linear, tree-based, ensemble, KNN) and
    ranks them by cross-validation score. For classification tasks, uses
    accuracy; for regression, uses R-squared.

    Args:
        X: Feature matrix (n_samples x n_features).
        y: Target values (n_samples,).
        task: Either 'classification' or 'regression'.
        cv: Number of cross-validation folds.

    Returns:
        Dictionary containing:
            - best_model_type: Name of the best-performing model type.
            - rankings: List of (model_type, score) tuples sorted by score.
            - cv_results_per_model: Dict mapping model type to its CV
              score and metadata.

    Raises:
        ValueError: If task is not 'classification' or 'regression'.
    """
    if task not in ("classification", "regression"):
        raise ValueError(f"task must be 'classification' or 'regression', got '{task}'")

    X_list = _to_2d_list(X)
    y_list = _to_1d_list(y)
    n_samples, n_features = len(X_list), len(X_list[0]) if X_list else 0

    metric = "accuracy" if task == "classification" else "r2"

    logger.info(
        "Model selection: task=%s, %d samples, %d features, cv=%d",
        task,
        n_samples,
        n_features,
        cv,
    )

    model_configs: list[dict] = []

    if HAS_SKLEARN:
        if task == "classification":
            from sklearn.ensemble import GradientBoostingClassifier, RandomForestClassifier
            from sklearn.linear_model import LogisticRegression
            from sklearn.neighbors import KNeighborsClassifier
            from sklearn.svm import SVC
            from sklearn.tree import DecisionTreeClassifier

            model_configs = [
                {"name": "logistic_regression", "model": LogisticRegression(max_iter=1000, random_state=42)},
                {"name": "decision_tree", "model": DecisionTreeClassifier(max_depth=10, random_state=42)},
                {
                    "name": "random_forest",
                    "model": RandomForestClassifier(n_estimators=100, max_depth=10, random_state=42),
                },
                {
                    "name": "gradient_boosting",
                    "model": GradientBoostingClassifier(n_estimators=100, max_depth=5, random_state=42),
                },
                {"name": "knn", "model": KNeighborsClassifier(n_neighbors=min(5, n_samples - 1))},
            ]
            # Only add SVC if dataset is not too large
            if n_samples <= 10000:
                model_configs.append({"name": "svm", "model": SVC(kernel="rbf", random_state=42)})
        else:
            from sklearn.ensemble import GradientBoostingRegressor, RandomForestRegressor
            from sklearn.linear_model import Lasso, LinearRegression, Ridge
            from sklearn.neighbors import KNeighborsRegressor
            from sklearn.tree import DecisionTreeRegressor

            model_configs = [
                {"name": "linear_regression", "model": LinearRegression()},
                {"name": "ridge", "model": Ridge(alpha=1.0)},
                {"name": "lasso", "model": Lasso(alpha=0.1, max_iter=1000)},
                {"name": "decision_tree", "model": DecisionTreeRegressor(max_depth=10, random_state=42)},
                {
                    "name": "random_forest",
                    "model": RandomForestRegressor(n_estimators=100, max_depth=10, random_state=42),
                },
                {
                    "name": "gradient_boosting",
                    "model": GradientBoostingRegressor(n_estimators=100, max_depth=5, random_state=42),
                },
                {"name": "knn", "model": KNeighborsRegressor(n_neighbors=min(5, n_samples - 1))},
            ]
    else:
        logger.warning("scikit-learn not available; model selection limited to pure Python models")
        model_configs = []

    cv_results: dict[str, dict] = {}

    for config in model_configs:
        model_name = config["name"]
        model = config["model"]

        try:
            score = _cross_validate(model, X_list, y_list, cv=cv, metric=metric)
            cv_results[model_name] = {
                "score": round(score, 6),
                "model_type": model_name,
                "status": "success",
            }
            logger.info("Model %s: CV score = %.4f", model_name, score)
        except Exception as exc:
            cv_results[model_name] = {
                "score": 0.0,
                "model_type": model_name,
                "status": f"failed: {exc}",
            }
            logger.warning("Model %s failed: %s", model_name, exc)

    # Rank by score
    rankings = sorted(
        [(name, info["score"]) for name, info in cv_results.items()],
        key=lambda t: t[1],
        reverse=True,
    )

    best_model_type = rankings[0][0] if rankings else "none"

    logger.info("Model selection complete: best=%s (score=%.4f)", best_model_type, rankings[0][1] if rankings else 0.0)

    return {
        "best_model_type": best_model_type,
        "rankings": rankings,
        "cv_results_per_model": cv_results,
    }


def auto_preprocess(
    X: Any,
    y: Any | None = None,
) -> dict:
    """Automatic preprocessing pipeline.

    Analyzes the input data and automatically applies appropriate
    preprocessing steps: detects column data types, imputes missing values,
    scales numeric features, and encodes categorical features.

    Args:
        X: Feature matrix (n_samples x n_features). Can contain mixed types.
        y: Optional target values (used for supervised imputation hints).

    Returns:
        Dictionary containing:
            - X_processed: Preprocessed feature matrix as list of lists.
            - transformations_applied: List of transformation descriptions.
            - feature_info: List of dicts describing each feature's type
              and transformation.

    Raises:
        ValueError: If X is empty.
    """
    n_samples, n_features = _get_shape(X)
    if n_samples == 0:
        raise ValueError("X must not be empty")

    logger.info("Auto-preprocessing: %d samples, %d features", n_samples, n_features)

    X_list = _to_2d_list(X)
    transformations: list[str] = []
    feature_info: list[dict] = []

    # Analyze each feature
    processed_columns: list[list[float]] = []

    for feat in range(n_features):
        col = [X_list[i][feat] for i in range(n_samples)]

        # Detect missing values (NaN)
        has_missing = any(math.isnan(v) if isinstance(v, float) and not math.isinf(v) else False for v in col)

        # Detect if column is categorical-like (few unique values relative to n)
        unique_vals = set(col)
        unique_ratio = len(unique_vals) / n_samples
        is_categorical = unique_ratio < 0.05 and len(unique_vals) <= 20

        # Detect if binary
        is_binary = len(unique_vals) <= 2

        info: dict[str, Any] = {
            "index": feat,
            "n_unique": len(unique_vals),
            "has_missing": has_missing,
            "is_categorical": is_categorical,
            "is_binary": is_binary,
            "transformations": [],
        }

        # Step 1: Impute missing values
        if has_missing:
            valid_vals = [v for v in col if not (isinstance(v, float) and math.isnan(v))]
            if valid_vals:
                impute_val = sum(valid_vals) / len(valid_vals)
            else:
                impute_val = 0.0

            col = [impute_val if (isinstance(v, float) and math.isnan(v)) else v for v in col]
            info["transformations"].append("mean_imputation")
            transformations.append(f"Feature {feat}: mean imputation")

        # Step 2: Scale numeric (non-binary, non-categorical) features
        if not is_binary and not is_categorical:
            mean_val = sum(col) / n_samples
            var_val = sum((v - mean_val) ** 2 for v in col) / n_samples
            std_val = math.sqrt(var_val) if var_val > 1e-15 else 1.0

            col = [(v - mean_val) / std_val for v in col]
            info["transformations"].append("standard_scaling")
            transformations.append(f"Feature {feat}: standard scaling")

        processed_columns.append(col)
        feature_info.append(info)

    # Reconstruct X_processed
    X_processed = [[processed_columns[feat][i] for feat in range(n_features)] for i in range(n_samples)]

    logger.info(
        "Preprocessing complete: %d transformations applied",
        len(transformations),
    )

    return {
        "X_processed": X_processed,
        "transformations_applied": transformations,
        "feature_info": feature_info,
    }
