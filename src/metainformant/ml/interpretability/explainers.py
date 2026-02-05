"""Model interpretability and explainability methods.

This module provides implementations of popular model explanation techniques
including permutation importance, Kernel SHAP, LIME, partial dependence, and
feature interaction analysis. All methods work with any model that exposes a
predict or predict_proba callable, with pure Python fallbacks where possible.
"""

from __future__ import annotations

import math
import random
from collections import defaultdict
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
    from scipy import stats as scipy_stats

    HAS_SCIPY = True
except ImportError:
    HAS_SCIPY = False
    scipy_stats = None  # type: ignore[assignment]


def _to_2d_list(X: Any) -> list[list[float]]:
    """Convert input matrix to a list of lists of floats.

    Args:
        X: Input matrix (numpy array or list of lists).

    Returns:
        Matrix as list of lists.
    """
    if HAS_NUMPY and isinstance(X, np.ndarray):
        return [[float(X[i, j]) for j in range(X.shape[1])] for i in range(X.shape[0])]
    return [[float(v) for v in row] for row in X]


def _to_1d_list(y: Any) -> list[float]:
    """Convert input vector to a list of floats.

    Args:
        y: Input vector (numpy array or list).

    Returns:
        Vector as list of floats.
    """
    if HAS_NUMPY and isinstance(y, np.ndarray):
        return [float(v) for v in y.ravel()]
    return [float(v) for v in y]


def _get_shape(X: Any) -> tuple[int, int]:
    """Get shape of a 2D matrix.

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


def _predict_helper(model: Any, X: list[list[float]]) -> list[float]:
    """Call model.predict on data, handling numpy conversion.

    Args:
        model: Model with a predict method.
        X: Input data as list of lists.

    Returns:
        Predictions as list of floats.
    """
    if HAS_NUMPY:
        arr = np.array(X)
        preds = model.predict(arr)
        return [float(v) for v in preds]
    return [float(v) for v in model.predict(X)]


def _score_metric(y_true: list[float], y_pred: list[float], metric: str) -> float:
    """Compute a scoring metric.

    Args:
        y_true: True labels.
        y_pred: Predicted labels/values.
        metric: Metric name ('accuracy', 'mse', 'r2').

    Returns:
        Metric value.
    """
    n = len(y_true)
    if n == 0:
        return 0.0

    if metric == "accuracy":
        correct = sum(1 for yt, yp in zip(y_true, y_pred) if round(yt) == round(yp))
        return correct / n

    if metric == "mse":
        mse = sum((yt - yp) ** 2 for yt, yp in zip(y_true, y_pred)) / n
        return -mse  # Negative so higher is better

    if metric == "r2":
        mean_y = sum(y_true) / n
        ss_total = sum((yt - mean_y) ** 2 for yt in y_true)
        ss_res = sum((yt - yp) ** 2 for yt, yp in zip(y_true, y_pred))
        return 1.0 - ss_res / ss_total if ss_total > 1e-15 else 0.0

    # Default to accuracy
    correct = sum(1 for yt, yp in zip(y_true, y_pred) if round(yt) == round(yp))
    return correct / n


def compute_permutation_importance(
    model: Any,
    X: Any,
    y: Any,
    n_repeats: int = 10,
    metric: str = "accuracy",
    random_state: int | None = None,
) -> dict:
    """Compute permutation feature importance.

    Measures the decrease in model performance when each feature is randomly
    shuffled, breaking the relationship between the feature and the target.
    Features whose shuffling causes a large performance drop are considered
    important.

    Args:
        model: Fitted model with a predict method.
        X: Feature matrix (n_samples x n_features).
        y: True target values (n_samples,).
        n_repeats: Number of times to shuffle each feature.
        metric: Scoring metric ('accuracy', 'mse', 'r2').
        random_state: Random seed for reproducibility.

    Returns:
        Dictionary containing:
            - importances_mean: List of mean importance per feature.
            - importances_std: List of std of importance per feature.
            - feature_names: List of feature indices as strings.
            - baseline_score: Model score without permutation.

    Raises:
        ValueError: If X and y have incompatible shapes.
    """
    n_samples, n_features = _get_shape(X)
    y_list = _to_1d_list(y)
    X_list = _to_2d_list(X)

    if len(y_list) != n_samples:
        raise ValueError(f"y length ({len(y_list)}) must match X rows ({n_samples})")

    if random_state is not None:
        random.seed(random_state)

    logger.info(
        "Computing permutation importance: %d features, %d repeats, metric=%s",
        n_features,
        n_repeats,
        metric,
    )

    # Baseline score
    baseline_preds = _predict_helper(model, X_list)
    baseline_score = _score_metric(y_list, baseline_preds, metric)

    importances_all: list[list[float]] = []

    for feat_idx in range(n_features):
        feat_importances: list[float] = []

        for _ in range(n_repeats):
            # Create shuffled copy
            X_shuffled = [row[:] for row in X_list]
            col_values = [X_shuffled[i][feat_idx] for i in range(n_samples)]
            random.shuffle(col_values)
            for i in range(n_samples):
                X_shuffled[i][feat_idx] = col_values[i]

            shuffled_preds = _predict_helper(model, X_shuffled)
            shuffled_score = _score_metric(y_list, shuffled_preds, metric)
            feat_importances.append(baseline_score - shuffled_score)

        importances_all.append(feat_importances)

    # Compute mean and std
    importances_mean: list[float] = []
    importances_std: list[float] = []

    for feat_imps in importances_all:
        mean_val = sum(feat_imps) / len(feat_imps)
        var_val = sum((v - mean_val) ** 2 for v in feat_imps) / max(len(feat_imps) - 1, 1)
        importances_mean.append(round(mean_val, 8))
        importances_std.append(round(math.sqrt(var_val), 8))

    feature_names = [f"feature_{i}" for i in range(n_features)]

    logger.info("Permutation importance computed. Baseline score: %.4f", baseline_score)

    return {
        "importances_mean": importances_mean,
        "importances_std": importances_std,
        "feature_names": feature_names,
        "baseline_score": round(baseline_score, 8),
    }


def compute_shap_values_kernel(
    predict_fn: Any,
    X: Any,
    n_samples: int = 100,
    background: Any | None = None,
) -> dict:
    """Kernel SHAP approximation for model explanations.

    Approximates Shapley values by sampling coalition vectors (binary masks
    indicating which features are present), weighting each coalition by the
    Shapley kernel weight, and solving a weighted linear regression for each
    instance to determine feature contributions.

    Args:
        predict_fn: Callable that takes a 2D array and returns predictions.
            Can be model.predict or model.predict_proba.
        X: Instances to explain (n_instances x n_features).
        n_samples: Number of coalition samples per instance.
        background: Background dataset for computing baseline predictions.
            If None, uses the mean of X.

    Returns:
        Dictionary containing:
            - shap_values: Matrix of SHAP values (n_instances x n_features)
              as list of lists.
            - expected_value: Baseline prediction (mean prediction on
              background).
            - feature_names: List of feature index strings.

    Raises:
        ValueError: If X is empty.
    """
    n_instances, n_features = _get_shape(X)
    if n_instances == 0:
        raise ValueError("X must not be empty")

    X_list = _to_2d_list(X)

    logger.info(
        "Computing Kernel SHAP: %d instances, %d features, %d coalition samples",
        n_instances,
        n_features,
        n_samples,
    )

    # Compute background (mean of X or provided background)
    if background is not None:
        bg_list = _to_2d_list(background)
    else:
        bg_list = [X_list[0][:]]  # Use first instance as simple background
        for j in range(n_features):
            bg_list[0][j] = sum(X_list[i][j] for i in range(n_instances)) / n_instances

    # Get baseline prediction
    if HAS_NUMPY:
        bg_pred = predict_fn(np.array(bg_list))
        expected_value = float(bg_pred[0]) if hasattr(bg_pred, "__len__") else float(bg_pred)
    else:
        bg_pred = predict_fn(bg_list)
        expected_value = float(bg_pred[0]) if hasattr(bg_pred, "__len__") else float(bg_pred)

    shap_values: list[list[float]] = []

    for inst_idx in range(n_instances):
        instance = X_list[inst_idx]

        # Generate coalition vectors and compute weights
        coalitions: list[list[int]] = []
        weights: list[float] = []
        predictions: list[float] = []

        for _ in range(n_samples):
            # Sample a random coalition (binary mask)
            n_present = random.randint(1, n_features - 1)
            coalition = [0] * n_features
            indices = random.sample(range(n_features), n_present)
            for idx in indices:
                coalition[idx] = 1
            coalitions.append(coalition)

            # Shapley kernel weight: (M-1) / (C(M, |S|) * |S| * (M - |S|))
            s = sum(coalition)
            if 0 < s < n_features:
                # Use log to avoid overflow
                log_weight = -(_log_comb(n_features, s) + math.log(s) + math.log(n_features - s))
                weight = math.exp(log_weight) * (n_features - 1)
            else:
                weight = 1e6  # Very high weight for empty/full coalitions

            weights.append(weight)

            # Create mixed instance: use instance values where coalition=1,
            # background values where coalition=0
            mixed = [instance[j] if coalition[j] == 1 else bg_list[0][j] for j in range(n_features)]

            if HAS_NUMPY:
                pred = predict_fn(np.array([mixed]))
            else:
                pred = predict_fn([mixed])

            pred_val = float(pred[0]) if hasattr(pred, "__len__") else float(pred)
            predictions.append(pred_val - expected_value)

        # Solve weighted linear regression: coalition * shap = predictions
        feature_shaps = _weighted_linear_regression(coalitions, predictions, weights, n_features)
        shap_values.append([round(v, 8) for v in feature_shaps])

    feature_names = [f"feature_{i}" for i in range(n_features)]

    logger.info("Kernel SHAP complete for %d instances", n_instances)

    return {
        "shap_values": shap_values,
        "expected_value": round(expected_value, 8),
        "feature_names": feature_names,
    }


def _log_comb(n: int, k: int) -> float:
    """Compute log of binomial coefficient C(n, k).

    Args:
        n: Total items.
        k: Items chosen.

    Returns:
        Natural log of C(n, k).
    """
    if k < 0 or k > n:
        return -math.inf
    if k == 0 or k == n:
        return 0.0
    # Use lgamma: log(C(n,k)) = lgamma(n+1) - lgamma(k+1) - lgamma(n-k+1)
    return math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1)


def _weighted_linear_regression(
    X: list[list[int]],
    y: list[float],
    weights: list[float],
    n_features: int,
) -> list[float]:
    """Solve weighted linear regression using normal equations.

    Solves X^T W X beta = X^T W y where W is diagonal weight matrix.

    Args:
        X: Design matrix (n_samples x n_features), binary coalitions.
        y: Target values.
        weights: Sample weights.
        n_features: Number of features.

    Returns:
        Coefficient vector.
    """
    n = len(y)
    if n == 0:
        return [0.0] * n_features

    # X^T W X (n_features x n_features)
    xtx: list[list[float]] = [[0.0] * n_features for _ in range(n_features)]
    xty: list[float] = [0.0] * n_features

    for i in range(n):
        w = weights[i]
        for j in range(n_features):
            if X[i][j] == 0:
                continue
            xty[j] += w * X[i][j] * y[i]
            for k in range(j, n_features):
                if X[i][k] == 0:
                    continue
                val = w * X[i][j] * X[i][k]
                xtx[j][k] += val
                if j != k:
                    xtx[k][j] += val

    # Add small regularization for numerical stability
    for j in range(n_features):
        xtx[j][j] += 1e-8

    # Solve via Gauss elimination
    return _solve_linear_system(xtx, xty)


def _solve_linear_system(A: list[list[float]], b: list[float]) -> list[float]:
    """Solve linear system Ax = b via Gaussian elimination with partial pivoting.

    Args:
        A: Square matrix (n x n).
        b: Right-hand side vector (n,).

    Returns:
        Solution vector x.
    """
    n = len(b)
    # Augmented matrix
    aug = [A[i][:] + [b[i]] for i in range(n)]

    # Forward elimination with partial pivoting
    for col in range(n):
        # Find pivot
        max_val = abs(aug[col][col])
        max_row = col
        for row in range(col + 1, n):
            if abs(aug[row][col]) > max_val:
                max_val = abs(aug[row][col])
                max_row = row
        aug[col], aug[max_row] = aug[max_row], aug[col]

        if abs(aug[col][col]) < 1e-15:
            continue

        for row in range(col + 1, n):
            factor = aug[row][col] / aug[col][col]
            for j in range(col, n + 1):
                aug[row][j] -= factor * aug[col][j]

    # Back substitution
    x = [0.0] * n
    for i in range(n - 1, -1, -1):
        if abs(aug[i][i]) < 1e-15:
            x[i] = 0.0
            continue
        x[i] = aug[i][n]
        for j in range(i + 1, n):
            x[i] -= aug[i][j] * x[j]
        x[i] /= aug[i][i]

    return x


def compute_lime_explanation(
    predict_fn: Any,
    instance: list[float],
    feature_names: list[str],
    n_samples: int = 1000,
    n_features: int = 10,
    kernel_width: float | None = None,
) -> dict:
    """Generate a LIME (Local Interpretable Model-agnostic Explanations) explanation.

    Generates perturbed samples around the instance, weights them by
    proximity to the original instance, and fits a local weighted linear
    model to approximate the black-box model locally.

    Args:
        predict_fn: Callable that takes a 2D array and returns predictions.
        instance: Single instance to explain (1D list of features).
        feature_names: Names for each feature.
        n_samples: Number of perturbed samples to generate.
        n_features: Number of top features to include in the explanation.
        kernel_width: Width of the exponential kernel for proximity
            weighting. If None, uses sqrt(n_features) * 0.75.

    Returns:
        Dictionary containing:
            - coefficients: List of (feature_name, coefficient) tuples for
              top features.
            - intercept: Intercept of the local linear model.
            - local_prediction: Prediction of the local model at the
              original instance.
            - r_squared: R-squared of the local model fit.
            - feature_contributions: Dict mapping feature names to their
              contribution values.

    Raises:
        ValueError: If instance length does not match feature_names length.
    """
    n_feat = len(instance)
    if len(feature_names) != n_feat:
        raise ValueError(f"feature_names length ({len(feature_names)}) must match " f"instance length ({n_feat})")

    if kernel_width is None:
        kernel_width = math.sqrt(n_feat) * 0.75

    logger.info(
        "Computing LIME explanation: %d features, %d samples, kernel_width=%.3f",
        n_feat,
        n_samples,
        kernel_width,
    )

    # Compute feature statistics from instance (use instance +/- noise)
    feature_std = [max(abs(v) * 0.3, 0.01) for v in instance]

    # Generate perturbed samples
    perturbed: list[list[float]] = []
    binary_repr: list[list[int]] = []  # Whether each feature is "on" or "off"

    for _ in range(n_samples):
        sample = []
        binary = []
        for j in range(n_feat):
            # Randomly decide whether to perturb this feature
            if random.random() < 0.5:
                # Keep original value
                sample.append(instance[j])
                binary.append(1)
            else:
                # Perturb with Gaussian noise
                noise = random.gauss(0, feature_std[j])
                sample.append(instance[j] + noise)
                binary.append(0)
        perturbed.append(sample)
        binary_repr.append(binary)

    # Get predictions for perturbed samples
    if HAS_NUMPY:
        predictions_raw = predict_fn(np.array(perturbed))
    else:
        predictions_raw = predict_fn(perturbed)

    predictions = [
        float(predictions_raw[i]) if hasattr(predictions_raw, "__getitem__") else float(predictions_raw)
        for i in range(n_samples)
    ]

    # Compute proximity weights using exponential kernel
    weights: list[float] = []
    for i in range(n_samples):
        distance_sq = sum(((perturbed[i][j] - instance[j]) / max(feature_std[j], 1e-10)) ** 2 for j in range(n_feat))
        distance = math.sqrt(distance_sq)
        weight = math.exp(-(distance**2) / (kernel_width**2))
        weights.append(weight)

    # Fit weighted linear regression on binary representation
    # This gives us feature importance in terms of "feature on/off"
    coefficients_raw = _weighted_linear_regression(binary_repr, predictions, weights, n_feat)

    # Compute intercept (mean prediction when all features are "off")
    mean_pred = sum(w * p for w, p in zip(weights, predictions)) / max(sum(weights), 1e-15)
    intercept = mean_pred - sum(
        coefficients_raw[j] * sum(weights[i] * binary_repr[i][j] for i in range(n_samples)) / max(sum(weights), 1e-15)
        for j in range(n_feat)
    )

    # Local prediction at original instance (all features "on")
    local_prediction = intercept + sum(coefficients_raw)

    # Compute R-squared
    pred_local = [
        intercept + sum(coefficients_raw[j] * binary_repr[i][j] for j in range(n_feat)) for i in range(n_samples)
    ]
    ss_res = sum(weights[i] * (predictions[i] - pred_local[i]) ** 2 for i in range(n_samples))
    ss_tot = sum(weights[i] * (predictions[i] - mean_pred) ** 2 for i in range(n_samples))
    r_squared = 1.0 - ss_res / ss_tot if ss_tot > 1e-15 else 0.0

    # Select top features by absolute coefficient
    indexed_coefs = [(feature_names[j], coefficients_raw[j]) for j in range(n_feat)]
    indexed_coefs.sort(key=lambda x: abs(x[1]), reverse=True)
    top_coefficients = indexed_coefs[:n_features]

    feature_contributions = {feature_names[j]: round(coefficients_raw[j], 8) for j in range(n_feat)}

    logger.info(
        "LIME explanation: R^2=%.4f, top feature=%s (coef=%.4f)",
        r_squared,
        top_coefficients[0][0] if top_coefficients else "none",
        top_coefficients[0][1] if top_coefficients else 0.0,
    )

    return {
        "coefficients": [(name, round(coef, 8)) for name, coef in top_coefficients],
        "intercept": round(intercept, 8),
        "local_prediction": round(local_prediction, 8),
        "r_squared": round(max(0.0, r_squared), 6),
        "feature_contributions": feature_contributions,
    }


def feature_interaction(
    model: Any,
    X: Any,
    feature_i: int,
    feature_j: int,
    n_grid: int = 20,
) -> dict:
    """Compute 2D partial dependence and feature interaction effects.

    Evaluates the model prediction surface over a grid of values for two
    features, while marginalizing over the remaining features. The interaction
    effect is the deviation from the sum of individual partial dependences.

    Args:
        model: Fitted model with a predict method.
        X: Feature matrix (n_samples x n_features).
        feature_i: Index of the first feature.
        feature_j: Index of the second feature.
        n_grid: Number of grid points per feature dimension.

    Returns:
        Dictionary containing:
            - grid_i: List of grid values for feature_i.
            - grid_j: List of grid values for feature_j.
            - pdp_values: 2D matrix (n_grid x n_grid) of partial
              dependence values.
            - interaction_strength: Scalar measure of interaction strength
              (variance of the interaction term).

    Raises:
        ValueError: If feature indices are out of range or equal.
    """
    n_samples, n_features = _get_shape(X)
    X_list = _to_2d_list(X)

    if feature_i == feature_j:
        raise ValueError("feature_i and feature_j must be different")
    if feature_i < 0 or feature_i >= n_features:
        raise ValueError(f"feature_i ({feature_i}) out of range [0, {n_features})")
    if feature_j < 0 or feature_j >= n_features:
        raise ValueError(f"feature_j ({feature_j}) out of range [0, {n_features})")

    logger.info(
        "Computing feature interaction: features %d x %d, grid=%d",
        feature_i,
        feature_j,
        n_grid,
    )

    # Create grid values
    vals_i = [X_list[s][feature_i] for s in range(n_samples)]
    vals_j = [X_list[s][feature_j] for s in range(n_samples)]

    min_i, max_i = min(vals_i), max(vals_i)
    min_j, max_j = min(vals_j), max(vals_j)

    range_i = max_i - min_i if max_i > min_i else 1.0
    range_j = max_j - min_j if max_j > min_j else 1.0

    grid_i = [min_i + k * range_i / (n_grid - 1) for k in range(n_grid)]
    grid_j = [min_j + k * range_j / (n_grid - 1) for k in range(n_grid)]

    # Compute 2D partial dependence
    pdp_values: list[list[float]] = []

    for gi in range(n_grid):
        row: list[float] = []
        for gj in range(n_grid):
            # Create modified dataset with feature_i = grid_i[gi], feature_j = grid_j[gj]
            modified = [r[:] for r in X_list]
            for s in range(n_samples):
                modified[s][feature_i] = grid_i[gi]
                modified[s][feature_j] = grid_j[gj]

            preds = _predict_helper(model, modified)
            mean_pred = sum(preds) / len(preds)
            row.append(round(mean_pred, 8))
        pdp_values.append(row)

    # Compute 1D partial dependences for interaction strength
    pdp_i: list[float] = []
    for gi in range(n_grid):
        modified = [r[:] for r in X_list]
        for s in range(n_samples):
            modified[s][feature_i] = grid_i[gi]
        preds = _predict_helper(model, modified)
        pdp_i.append(sum(preds) / len(preds))

    pdp_j: list[float] = []
    for gj in range(n_grid):
        modified = [r[:] for r in X_list]
        for s in range(n_samples):
            modified[s][feature_j] = grid_j[gj]
        preds = _predict_helper(model, modified)
        pdp_j.append(sum(preds) / len(preds))

    # Overall mean
    overall_preds = _predict_helper(model, X_list)
    overall_mean = sum(overall_preds) / len(overall_preds)

    # Interaction = pdp_ij - pdp_i - pdp_j + overall_mean
    interaction_values: list[float] = []
    for gi in range(n_grid):
        for gj in range(n_grid):
            interaction_val = pdp_values[gi][gj] - pdp_i[gi] - pdp_j[gj] + overall_mean
            interaction_values.append(interaction_val)

    # Interaction strength = variance of interaction term
    mean_interaction = sum(interaction_values) / len(interaction_values)
    interaction_strength = sum((v - mean_interaction) ** 2 for v in interaction_values) / len(interaction_values)

    logger.info("Feature interaction strength: %.6f", interaction_strength)

    return {
        "grid_i": [round(v, 8) for v in grid_i],
        "grid_j": [round(v, 8) for v in grid_j],
        "pdp_values": pdp_values,
        "interaction_strength": round(interaction_strength, 8),
    }


def partial_dependence(
    model: Any,
    X: Any,
    feature: int,
    n_grid: int = 50,
) -> dict:
    """Compute 1D partial dependence plot data.

    Shows the marginal effect of a single feature on the predicted outcome,
    averaging over the values of all other features (i.e., marginalizing
    over the dataset).

    Args:
        model: Fitted model with a predict method.
        X: Feature matrix (n_samples x n_features).
        feature: Index of the feature to compute PDP for.
        n_grid: Number of grid points.

    Returns:
        Dictionary containing:
            - grid_values: List of feature values on the grid.
            - pdp_mean: Mean partial dependence at each grid point.
            - pdp_std: Std of predictions at each grid point.
            - ice_curves: Individual Conditional Expectation curves, a list
              of n_samples lists of length n_grid (capped at 50 samples).

    Raises:
        ValueError: If feature index is out of range.
    """
    n_samples, n_features = _get_shape(X)
    X_list = _to_2d_list(X)

    if feature < 0 or feature >= n_features:
        raise ValueError(f"feature ({feature}) out of range [0, {n_features})")

    logger.info("Computing partial dependence for feature %d, grid=%d", feature, n_grid)

    # Create grid
    vals = [X_list[s][feature] for s in range(n_samples)]
    min_val, max_val = min(vals), max(vals)
    val_range = max_val - min_val if max_val > min_val else 1.0
    grid_values = [min_val + k * val_range / (n_grid - 1) for k in range(n_grid)]

    # Compute ICE curves and PDP
    max_ice_samples = min(n_samples, 50)
    ice_curves: list[list[float]] = [[] for _ in range(max_ice_samples)]
    pdp_mean: list[float] = []
    pdp_std: list[float] = []

    for gi, grid_val in enumerate(grid_values):
        modified = [r[:] for r in X_list]
        for s in range(n_samples):
            modified[s][feature] = grid_val

        preds = _predict_helper(model, modified)

        mean_val = sum(preds) / len(preds)
        var_val = sum((p - mean_val) ** 2 for p in preds) / max(len(preds) - 1, 1)

        pdp_mean.append(round(mean_val, 8))
        pdp_std.append(round(math.sqrt(var_val), 8))

        for s in range(max_ice_samples):
            ice_curves[s].append(round(preds[s], 8))

    logger.info("Partial dependence computed for feature %d", feature)

    return {
        "grid_values": [round(v, 8) for v in grid_values],
        "pdp_mean": pdp_mean,
        "pdp_std": pdp_std,
        "ice_curves": ice_curves,
    }


def compute_attention_weights(
    model: Any,
    X: Any,
) -> dict | None:
    """Extract attention weights from a model if it supports them.

    Checks whether the model has an attention mechanism (common in
    transformer-based architectures) and extracts the attention weight
    matrices for the given input.

    Args:
        model: Model that may have attention weights. Must have a
            get_attention_weights method or an attention_weights attribute.
        X: Input data (n_instances x n_features).

    Returns:
        Dictionary containing:
            - attention_matrix: Attention weights for the last layer
              (list of lists).
            - layer_attentions: List of attention matrices per layer
              (if multi-layer).
        Returns None if the model does not support attention extraction.
    """
    logger.info("Attempting to extract attention weights from model")

    # Check for get_attention_weights method
    if hasattr(model, "get_attention_weights"):
        try:
            if HAS_NUMPY:
                weights = model.get_attention_weights(np.array(_to_2d_list(X)))
            else:
                weights = model.get_attention_weights(_to_2d_list(X))

            if isinstance(weights, dict):
                return weights

            # Convert to standard format
            if isinstance(weights, (list, tuple)):
                layer_attentions = []
                for layer_w in weights:
                    if HAS_NUMPY and isinstance(layer_w, np.ndarray):
                        layer_attentions.append(layer_w.tolist())
                    elif isinstance(layer_w, list):
                        layer_attentions.append(layer_w)

                return {
                    "attention_matrix": layer_attentions[-1] if layer_attentions else [],
                    "layer_attentions": layer_attentions,
                }

            logger.warning("Unexpected attention weight format: %s", type(weights))
            return None

        except (AttributeError, TypeError, RuntimeError) as exc:
            logger.warning("Failed to extract attention weights: %s", exc)
            return None

    # Check for attention_weights attribute
    if hasattr(model, "attention_weights"):
        try:
            weights = model.attention_weights
            if HAS_NUMPY and isinstance(weights, np.ndarray):
                return {
                    "attention_matrix": weights.tolist(),
                    "layer_attentions": [weights.tolist()],
                }
            if isinstance(weights, list):
                return {
                    "attention_matrix": weights,
                    "layer_attentions": [weights],
                }
        except (AttributeError, TypeError) as exc:
            logger.warning("Failed to access attention_weights: %s", exc)
            return None

    logger.info("Model does not support attention weight extraction")
    return None
