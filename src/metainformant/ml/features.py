"""Feature selection methods for biological data analysis."""

from __future__ import annotations

import warnings
from collections import defaultdict
from typing import Dict, List, Optional, Tuple, Union

import numpy as np


def select_features_univariate(
    X: np.ndarray, y: np.ndarray, method: str = "f_score", k: int = 100, p_threshold: float = 0.05
) -> Tuple[np.ndarray, List[int]]:
    """Select features using univariate statistical tests.

    Args:
        X: Feature matrix (samples x features)
        y: Target vector
        method: Statistical test ("f_score", "chi2", "mutual_info")
        k: Number of features to select
        p_threshold: P-value threshold for significance

    Returns:
        Tuple of (selected_features_matrix, selected_indices)
    """
    if len(X.shape) != 2:
        raise ValueError("X must be 2D array (samples x features)")
    if len(y) != X.shape[0]:
        raise ValueError("y length must match X samples")

    n_samples, n_features = X.shape
    k = min(k, n_features)

    if method == "f_score":
        scores, p_values = _f_score_test(X, y)
    elif method == "chi2":
        scores, p_values = _chi2_test(X, y)
    elif method == "mutual_info":
        scores = _mutual_info_score(X, y)
        p_values = np.ones(len(scores))  # MI doesn't have p-values
    else:
        raise ValueError(f"Unknown method: {method}")

    # Select by significance and/or top-k
    if method != "mutual_info":
        significant = p_values < p_threshold
        scores = scores * significant  # Zero out non-significant

    # Get top k features
    selected_indices = np.argsort(scores)[-k:][::-1]
    selected_X = X[:, selected_indices]

    return selected_X, selected_indices.tolist()


def select_features_recursive(
    X: np.ndarray,
    y: np.ndarray,
    estimator_type: str = "random_forest",
    n_features: int = 50,
    step: float = 0.1,
    cv_folds: int = 5,
    random_state: Optional[int] = None,
) -> Tuple[np.ndarray, List[int]]:
    """Recursive feature elimination with cross-validation.

    Args:
        X: Feature matrix
        y: Target vector
        estimator_type: Base estimator ("random_forest", "linear", "svm")
        n_features: Number of features to select
        step: Fraction of features to eliminate each step
        cv_folds: Cross-validation folds
        random_state: Random seed

    Returns:
        Tuple of (selected_features_matrix, selected_indices)
    """
    if random_state is not None:
        np.random.seed(random_state)

    n_samples, n_total_features = X.shape
    n_features = min(n_features, n_total_features)

    # Initialize with all features
    current_features = list(range(n_total_features))

    while len(current_features) > n_features:
        # Number of features to remove this step
        n_remove = max(1, int(len(current_features) * step))
        if len(current_features) - n_remove < n_features:
            n_remove = len(current_features) - n_features

        # Train estimator on current features
        X_current = X[:, current_features]
        feature_importance = _get_feature_importance(X_current, y, estimator_type, random_state)

        # Remove least important features
        worst_features = np.argsort(feature_importance)[:n_remove]
        for idx in sorted(worst_features, reverse=True):
            current_features.pop(idx)

    selected_X = X[:, current_features]
    return selected_X, current_features


def select_features_stability(
    X: np.ndarray,
    y: np.ndarray,
    method: str = "random_forest",
    n_bootstrap: int = 100,
    subsample_ratio: float = 0.8,
    stability_threshold: float = 0.6,
    random_state: Optional[int] = None,
) -> Tuple[np.ndarray, List[int]]:
    """Select features based on stability across bootstrap samples.

    Args:
        X: Feature matrix
        y: Target vector
        method: Feature selection method
        n_bootstrap: Number of bootstrap samples
        subsample_ratio: Fraction of samples in each bootstrap
        stability_threshold: Minimum selection frequency
        random_state: Random seed

    Returns:
        Tuple of (selected_features_matrix, selected_indices)
    """
    if random_state is not None:
        np.random.seed(random_state)

    n_samples, n_features = X.shape
    subsample_size = int(n_samples * subsample_ratio)

    # Track feature selection frequency
    selection_counts = np.zeros(n_features)

    for i in range(n_bootstrap):
        # Bootstrap sample
        bootstrap_indices = np.random.choice(n_samples, size=subsample_size, replace=True)
        X_boot = X[bootstrap_indices]
        y_boot = y[bootstrap_indices]

        # Select features for this bootstrap
        if method == "univariate":
            _, selected_features = select_features_univariate(X_boot, y_boot, k=n_features // 2)
        elif method == "random_forest":
            importance = _get_feature_importance(X_boot, y_boot, "random_forest", random_state)
            # Select top half features
            selected_features = np.argsort(importance)[-(n_features // 2) :]
        else:
            raise ValueError(f"Unknown method: {method}")

        # Update selection counts
        selection_counts[selected_features] += 1

    # Calculate selection frequencies
    selection_frequencies = selection_counts / n_bootstrap

    # Select stable features
    stable_features = np.where(selection_frequencies >= stability_threshold)[0]

    if len(stable_features) == 0:
        warnings.warn("No features meet stability threshold, selecting top features")
        stable_features = np.argsort(selection_frequencies)[-min(50, n_features) :]

    selected_X = X[:, stable_features]
    return selected_X, stable_features.tolist()


def biological_feature_ranking(
    X: np.ndarray,
    y: np.ndarray,
    feature_names: Optional[List[str]] = None,
    method: str = "combined",
    gene_weights: Optional[Dict[str, float]] = None,
) -> List[Tuple[int, str, float]]:
    """Rank features considering biological relevance.

    Args:
        X: Feature matrix
        y: Target vector
        feature_names: Names of features (genes, proteins, etc.)
        method: Ranking method ("combined", "statistical", "biological")
        gene_weights: Prior biological weights for features

    Returns:
        List of (index, name, score) tuples sorted by importance
    """
    n_features = X.shape[1]

    if feature_names is None:
        feature_names = [f"feature_{i}" for i in range(n_features)]

    if len(feature_names) != n_features:
        raise ValueError("feature_names length must match X columns")

    # Statistical scores
    stat_scores, _ = _f_score_test(X, y)
    stat_scores = (stat_scores - np.min(stat_scores)) / (np.max(stat_scores) - np.min(stat_scores) + 1e-10)

    # Biological weights
    if gene_weights is None:
        bio_scores = np.ones(n_features)
    else:
        bio_scores = np.array([gene_weights.get(name, 1.0) for name in feature_names])

    # Combined scoring
    if method == "statistical":
        final_scores = stat_scores
    elif method == "biological":
        final_scores = bio_scores
    elif method == "combined":
        # Weighted combination
        final_scores = 0.7 * stat_scores + 0.3 * bio_scores
    else:
        raise ValueError(f"Unknown method: {method}")

    # Create ranked list
    ranked_features = []
    for i, score in enumerate(final_scores):
        ranked_features.append((i, feature_names[i], score))

    # Sort by score (descending)
    ranked_features.sort(key=lambda x: x[2], reverse=True)

    return ranked_features


def _f_score_test(X: np.ndarray, y: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """Compute F-scores for features."""
    n_samples, n_features = X.shape

    # Get unique classes
    classes = np.unique(y)
    n_classes = len(classes)

    if n_classes < 2:
        raise ValueError("Need at least 2 classes for F-score test")

    f_scores = np.zeros(n_features)
    p_values = np.ones(n_features)

    for j in range(n_features):
        feature = X[:, j]

        # Group means
        group_means = []
        group_vars = []
        group_sizes = []

        for class_val in classes:
            class_mask = y == class_val
            class_data = feature[class_mask]

            if len(class_data) > 1:
                group_means.append(np.mean(class_data))
                group_vars.append(np.var(class_data, ddof=1))
                group_sizes.append(len(class_data))
            else:
                group_means.append(np.mean(class_data) if len(class_data) > 0 else 0)
                group_vars.append(0)
                group_sizes.append(len(class_data))

        # Overall mean
        overall_mean = np.mean(feature)

        # Between-group sum of squares
        ss_between = sum(n * (mean - overall_mean) ** 2 for n, mean in zip(group_sizes, group_means))

        # Within-group sum of squares
        ss_within = sum(n * var for n, var in zip(group_sizes, group_vars))

        # Degrees of freedom
        df_between = n_classes - 1
        df_within = n_samples - n_classes

        if df_within > 0 and ss_within > 0:
            f_score = (ss_between / df_between) / (ss_within / df_within)
            f_scores[j] = f_score

            # Approximate p-value (simplified)
            # In practice, would use F-distribution CDF
            p_values[j] = max(0.001, min(0.999, 1.0 / (1.0 + f_score)))

    return f_scores, p_values


def _chi2_test(X: np.ndarray, y: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """Chi-squared test for categorical features."""
    # Simplified implementation - assumes non-negative features
    X_pos = np.maximum(X, 0)  # Ensure non-negative

    n_features = X_pos.shape[1]
    chi2_scores = np.zeros(n_features)
    p_values = np.ones(n_features)

    classes = np.unique(y)

    for j in range(n_features):
        feature = X_pos[:, j]

        # Create contingency table (simplified)
        contingency = []
        for class_val in classes:
            class_mask = y == class_val
            class_feature_sum = np.sum(feature[class_mask])
            contingency.append(class_feature_sum)

        # Chi-squared statistic (simplified)
        expected = np.mean(contingency)
        if expected > 0:
            chi2_score = sum((observed - expected) ** 2 / expected for observed in contingency)
            chi2_scores[j] = chi2_score
            p_values[j] = max(0.001, min(0.999, 1.0 / (1.0 + chi2_score)))

    return chi2_scores, p_values


def _mutual_info_score(X: np.ndarray, y: np.ndarray) -> np.ndarray:
    """Simplified mutual information calculation."""
    n_features = X.shape[1]
    mi_scores = np.zeros(n_features)

    for j in range(n_features):
        feature = X[:, j]

        # Discretize continuous features (simplified)
        n_bins = min(10, len(np.unique(feature)))
        if n_bins > 1:
            try:
                feature_binned = np.digitize(feature, np.linspace(np.min(feature), np.max(feature), n_bins))
                mi_scores[j] = _calculate_mi(feature_binned, y)
            except:
                mi_scores[j] = 0.0

    return mi_scores


def _calculate_mi(x: np.ndarray, y: np.ndarray) -> float:
    """Calculate mutual information between discrete variables."""
    # Joint probability
    joint_counts = defaultdict(int)
    for xi, yi in zip(x, y):
        joint_counts[(xi, yi)] += 1

    n_samples = len(x)
    joint_probs = {k: v / n_samples for k, v in joint_counts.items()}

    # Marginal probabilities
    x_counts = defaultdict(int)
    y_counts = defaultdict(int)
    for xi, yi in zip(x, y):
        x_counts[xi] += 1
        y_counts[yi] += 1

    x_probs = {k: v / n_samples for k, v in x_counts.items()}
    y_probs = {k: v / n_samples for k, v in y_counts.items()}

    # Mutual information
    mi = 0.0
    for (xi, yi), joint_prob in joint_probs.items():
        if joint_prob > 0:
            marginal_prob = x_probs[xi] * y_probs[yi]
            if marginal_prob > 0:
                mi += joint_prob * np.log2(joint_prob / marginal_prob)

    return max(0.0, mi)


def _get_feature_importance(
    X: np.ndarray, y: np.ndarray, method: str, random_state: Optional[int] = None
) -> np.ndarray:
    """Get feature importance using specified method."""
    if method == "random_forest":
        return _random_forest_importance(X, y, random_state)
    elif method == "linear":
        return _linear_importance(X, y)
    elif method == "svm":
        return _svm_importance(X, y, random_state)
    else:
        raise ValueError(f"Unknown importance method: {method}")


def _random_forest_importance(X: np.ndarray, y: np.ndarray, random_state: Optional[int] = None) -> np.ndarray:
    """Simplified random forest feature importance."""
    # This is a very simplified version
    # In practice, would use proper random forest implementation

    if random_state is not None:
        np.random.seed(random_state)

    n_features = X.shape[1]
    importance = np.zeros(n_features)

    # Use correlation as proxy for importance
    for j in range(n_features):
        correlation = np.corrcoef(X[:, j], y)[0, 1]
        importance[j] = abs(correlation) if not np.isnan(correlation) else 0.0

    return importance


def _linear_importance(X: np.ndarray, y: np.ndarray) -> np.ndarray:
    """Linear regression coefficients as importance."""
    try:
        # Simple linear regression coefficients
        X_centered = X - np.mean(X, axis=0)
        y_centered = y - np.mean(y)

        # Normal equation (simplified, no regularization)
        XtX = X_centered.T @ X_centered
        Xty = X_centered.T @ y_centered

        # Add small regularization for stability
        regularization = 1e-6 * np.eye(X_centered.shape[1])
        coefficients = np.linalg.solve(XtX + regularization, Xty)

        return np.abs(coefficients)
    except:
        # Fallback to correlation
        return np.abs(
            [
                np.corrcoef(X[:, j], y)[0, 1] if not np.isnan(np.corrcoef(X[:, j], y)[0, 1]) else 0.0
                for j in range(X.shape[1])
            ]
        )


def _svm_importance(X: np.ndarray, y: np.ndarray, random_state: Optional[int] = None) -> np.ndarray:
    """SVM-based feature importance (simplified)."""
    # For simplicity, use linear coefficients as proxy
    return _linear_importance(X, y)
