"""Advanced feature selection methods for bioinformatics ML pipelines.

This module provides sophisticated feature selection techniques beyond
simple univariate methods, including Boruta all-relevant selection,
recursive feature elimination, stability selection, and mutual
information-based selection. All methods include pure Python fallbacks.
"""

from __future__ import annotations

import math
import random
from collections import defaultdict
from typing import Any

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
    from sklearn.ensemble import RandomForestClassifier, RandomForestRegressor
    from sklearn.model_selection import cross_val_score

    HAS_SKLEARN = True
except ImportError:
    HAS_SKLEARN = False


def _to_2d_list(X: Any) -> list[list[float]]:
    """Convert input matrix to list of lists.

    Args:
        X: Input matrix (numpy array or list of lists).

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


def _compute_feature_importances_rf(
    X_list: list[list[float]],
    y_list: list[float],
    n_samples: int,
    n_features: int,
    random_state: int | None = None,
) -> list[float]:
    """Compute feature importances using a simple decision stump ensemble.

    Pure Python random forest approximation that builds shallow decision
    stumps and measures feature importance by variance reduction.

    Args:
        X_list: Feature matrix as list of lists.
        y_list: Target values as list.
        n_samples: Number of samples.
        n_features: Number of features.
        random_state: Random seed.

    Returns:
        List of importance scores per feature.
    """
    if HAS_SKLEARN:
        X_arr = np.array(X_list)
        y_arr = np.array(y_list)
        rf = RandomForestClassifier(
            n_estimators=100,
            max_depth=5,
            random_state=random_state,
            n_jobs=-1,
        )
        try:
            rf.fit(X_arr, y_arr)
            return [float(v) for v in rf.feature_importances_]
        except Exception:
            pass

    # Pure Python fallback: variance reduction with random splits
    if random_state is not None:
        random.seed(random_state)

    importances = [0.0] * n_features
    n_trees = 50

    total_var = _variance(y_list)
    if total_var < 1e-15:
        return importances

    for _ in range(n_trees):
        # Bootstrap sample
        indices = [random.randint(0, n_samples - 1) for _ in range(n_samples)]
        boot_X = [X_list[i] for i in indices]
        boot_y = [y_list[i] for i in indices]

        # For each feature, find best split and measure variance reduction
        for feat in range(n_features):
            vals = [(boot_X[i][feat], boot_y[i]) for i in range(n_samples)]
            vals.sort(key=lambda t: t[0])

            best_reduction = 0.0
            for split_idx in range(1, n_samples):
                if vals[split_idx][0] == vals[split_idx - 1][0]:
                    continue

                left_y = [v[1] for v in vals[:split_idx]]
                right_y = [v[1] for v in vals[split_idx:]]

                left_var = _variance(left_y)
                right_var = _variance(right_y)

                weighted_var = (len(left_y) * left_var + len(right_y) * right_var) / n_samples

                reduction = total_var - weighted_var
                if reduction > best_reduction:
                    best_reduction = reduction

            importances[feat] += best_reduction / n_trees

    # Normalize
    total_imp = sum(importances)
    if total_imp > 1e-15:
        importances = [imp / total_imp for imp in importances]

    return importances


def _variance(values: list[float]) -> float:
    """Compute variance of a list of values.

    Args:
        values: List of numeric values.

    Returns:
        Variance.
    """
    n = len(values)
    if n < 2:
        return 0.0
    mean = sum(values) / n
    return sum((v - mean) ** 2 for v in values) / n


def boruta_selection(
    X: Any,
    y: Any,
    max_iter: int = 100,
    alpha: float = 0.05,
    random_state: int | None = None,
) -> dict:
    """Boruta all-relevant feature selection.

    Identifies all features that are relevant to the outcome by comparing
    real feature importances against shadow features (random permutations
    of real features). A feature is selected if its importance consistently
    exceeds the maximum shadow feature importance.

    The algorithm:
    1. Create shadow features by shuffling each real feature column.
    2. Train a model on both real and shadow features.
    3. Compare each real feature's importance to the max shadow importance.
    4. Use binomial test to determine if the feature is significantly
       better than random.

    Args:
        X: Feature matrix (n_samples x n_features).
        y: Target values (n_samples,).
        max_iter: Maximum number of Boruta iterations.
        alpha: Significance threshold for the binomial test.
        random_state: Random seed for reproducibility.

    Returns:
        Dictionary containing:
            - selected: List of indices of confirmed important features.
            - tentative: List of indices of tentatively important features.
            - rejected: List of indices of confirmed unimportant features.
            - importance_history: List of lists tracking importance scores
              across iterations (n_iter x n_features).

    Raises:
        ValueError: If X and y have incompatible shapes.
    """
    n_samples, n_features = _get_shape(X)
    X_list = _to_2d_list(X)
    y_list = _to_1d_list(y)

    if len(y_list) != n_samples:
        raise ValueError(f"y length ({len(y_list)}) must match X rows ({n_samples})")

    if random_state is not None:
        random.seed(random_state)

    logger.info(
        "Boruta selection: %d features, %d samples, max_iter=%d, alpha=%.4f",
        n_features,
        n_samples,
        max_iter,
        alpha,
    )

    # Track hits: how many times each feature beats max shadow
    hits = [0] * n_features
    importance_history: list[list[float]] = []
    n_completed = 0

    for iteration in range(max_iter):
        # Create shadow features by shuffling columns
        shadow_features: list[list[float]] = []
        for feat in range(n_features):
            col = [X_list[i][feat] for i in range(n_samples)]
            random.shuffle(col)
            shadow_features.append(col)

        # Combine real + shadow features
        combined_X: list[list[float]] = []
        for i in range(n_samples):
            row = X_list[i][:] + [shadow_features[f][i] for f in range(n_features)]
            combined_X.append(row)

        # Get feature importances
        combined_importances = _compute_feature_importances_rf(
            combined_X,
            y_list,
            n_samples,
            2 * n_features,
            random_state=(random_state + iteration) if random_state else None,
        )

        real_importances = combined_importances[:n_features]
        shadow_importances = combined_importances[n_features:]

        importance_history.append(real_importances)

        # Max shadow importance
        max_shadow = max(shadow_importances) if shadow_importances else 0.0

        # Count hits
        for feat in range(n_features):
            if real_importances[feat] > max_shadow:
                hits[feat] += 1

        n_completed = iteration + 1

        # Early stopping: check if all features are decided
        all_decided = True
        for feat in range(n_features):
            p_selected = _binomial_test_pvalue(hits[feat], n_completed, 0.5)
            p_rejected = _binomial_test_pvalue(n_completed - hits[feat], n_completed, 0.5)
            if p_selected > alpha and p_rejected > alpha:
                all_decided = False
                break

        if all_decided and n_completed >= 10:
            logger.info("Boruta converged at iteration %d", n_completed)
            break

    # Classify features
    selected: list[int] = []
    tentative: list[int] = []
    rejected: list[int] = []

    for feat in range(n_features):
        p_selected = _binomial_test_pvalue(hits[feat], n_completed, 0.5)
        if p_selected <= alpha:
            selected.append(feat)
        elif hits[feat] / max(n_completed, 1) > 0.3:
            tentative.append(feat)
        else:
            rejected.append(feat)

    logger.info(
        "Boruta complete: %d selected, %d tentative, %d rejected",
        len(selected),
        len(tentative),
        len(rejected),
    )

    return {
        "selected": selected,
        "tentative": tentative,
        "rejected": rejected,
        "importance_history": importance_history,
    }


def _binomial_test_pvalue(k: int, n: int, p: float) -> float:
    """Compute one-sided binomial test p-value (upper tail).

    Args:
        k: Number of successes.
        n: Number of trials.
        p: Null hypothesis probability.

    Returns:
        P-value.
    """
    if n == 0:
        return 1.0

    mean = n * p
    variance = n * p * (1 - p)

    if variance < 1e-15:
        return 0.0 if k > mean else 1.0

    z = (k - mean) / math.sqrt(variance)
    p_val = 0.5 * math.erfc(z / math.sqrt(2.0))
    return max(0.0, min(1.0, p_val))


def recursive_elimination(
    model: Any,
    X: Any,
    y: Any,
    n_features: int = 10,
    step: int = 1,
    cv: int = 5,
) -> dict:
    """Recursive feature elimination with cross-validation scoring.

    Iteratively removes the least important features, re-fitting the model
    at each step and evaluating via cross-validation. Tracks performance
    at each feature count to identify the optimal number of features.

    Args:
        model: Model with fit and predict methods. Should also provide
            feature_importances_ or coef_ after fitting.
        X: Feature matrix (n_samples x n_features).
        y: Target values (n_samples,).
        n_features: Target number of features to select.
        step: Number of features to remove per iteration.
        cv: Number of cross-validation folds.

    Returns:
        Dictionary containing:
            - selected_features: List of indices of selected features.
            - ranking: List of feature rankings (1 = selected, higher =
              eliminated earlier).
            - cv_scores: List of dicts with n_features and cv_score at
              each elimination step.

    Raises:
        ValueError: If n_features is larger than total features.
    """
    n_samples, total_features = _get_shape(X)
    X_list = _to_2d_list(X)
    y_list = _to_1d_list(y)

    if n_features > total_features:
        raise ValueError(f"n_features ({n_features}) must be <= total features ({total_features})")

    logger.info(
        "Recursive elimination: %d -> %d features, step=%d, cv=%d",
        total_features,
        n_features,
        step,
        cv,
    )

    # Track active features and their rankings
    active_features = list(range(total_features))
    ranking = [0] * total_features
    current_rank = total_features
    cv_scores: list[dict] = []

    while len(active_features) > n_features:
        # Build subset matrix
        X_subset = [[X_list[i][f] for f in active_features] for i in range(n_samples)]

        # Cross-validation score
        cv_score = _cross_validate_simple(model, X_subset, y_list, cv=cv)
        cv_scores.append({"n_features": len(active_features), "cv_score": round(cv_score, 6)})

        # Fit model and get importances
        importances = _get_model_importances(model, X_subset, y_list, len(active_features))

        # Remove least important features
        n_to_remove = min(step, len(active_features) - n_features)
        indexed = [(imp, idx) for idx, imp in enumerate(importances)]
        indexed.sort(key=lambda t: t[0])

        features_to_remove = [indexed[i][1] for i in range(n_to_remove)]
        features_to_remove_global = [active_features[i] for i in features_to_remove]

        for feat in features_to_remove_global:
            ranking[feat] = current_rank
            current_rank -= 1

        active_features = [f for i, f in enumerate(active_features) if i not in set(features_to_remove)]

    # Final CV score
    X_final = [[X_list[i][f] for f in active_features] for i in range(n_samples)]
    final_cv_score = _cross_validate_simple(model, X_final, y_list, cv=cv)
    cv_scores.append({"n_features": len(active_features), "cv_score": round(final_cv_score, 6)})

    # Assign rank 1 to remaining features
    for feat in active_features:
        ranking[feat] = 1

    logger.info(
        "RFE complete: selected %d features, best CV=%.4f",
        len(active_features),
        final_cv_score,
    )

    return {
        "selected_features": active_features,
        "ranking": ranking,
        "cv_scores": cv_scores,
    }


def _get_model_importances(
    model: Any,
    X_list: list[list[float]],
    y_list: list[float],
    n_features: int,
) -> list[float]:
    """Get feature importances from a fitted model.

    Args:
        model: Model to fit.
        X_list: Feature matrix as list of lists.
        y_list: Target values.
        n_features: Number of features.

    Returns:
        List of importance scores.
    """
    n_samples = len(y_list)

    try:
        if HAS_NUMPY:
            model.fit(np.array(X_list), np.array(y_list))
        else:
            model.fit(X_list, y_list)

        if hasattr(model, "feature_importances_"):
            return [float(v) for v in model.feature_importances_]
        if hasattr(model, "coef_"):
            coef = model.coef_
            if HAS_NUMPY and isinstance(coef, np.ndarray):
                if coef.ndim > 1:
                    return [float(abs(v)) for v in coef[0]]
                return [float(abs(v)) for v in coef]
            return [float(abs(v)) for v in coef]
    except Exception as exc:
        logger.warning("Model fitting failed: %s, using variance-based importance", exc)

    # Fallback: use correlation with target
    importances = []
    for feat in range(n_features):
        vals = [X_list[i][feat] for i in range(n_samples)]
        corr = _abs_correlation(vals, y_list)
        importances.append(corr)

    return importances


def _abs_correlation(x: list[float], y: list[float]) -> float:
    """Compute absolute Pearson correlation.

    Args:
        x: First variable.
        y: Second variable.

    Returns:
        Absolute correlation value.
    """
    n = len(x)
    if n < 3:
        return 0.0

    mean_x = sum(x) / n
    mean_y = sum(y) / n

    cov = sum((xi - mean_x) * (yi - mean_y) for xi, yi in zip(x, y))
    var_x = sum((xi - mean_x) ** 2 for xi in x)
    var_y = sum((yi - mean_y) ** 2 for yi in y)

    denom = math.sqrt(var_x * var_y)
    if denom < 1e-15:
        return 0.0

    return abs(cov / denom)


def _cross_validate_simple(
    model: Any,
    X_list: list[list[float]],
    y_list: list[float],
    cv: int = 5,
) -> float:
    """Simple cross-validation scoring.

    Args:
        model: Model with fit and predict methods.
        X_list: Feature matrix as list of lists.
        y_list: Target values.
        cv: Number of folds.

    Returns:
        Mean CV accuracy.
    """
    n = len(y_list)
    fold_size = n // cv
    scores: list[float] = []

    for fold in range(cv):
        start = fold * fold_size
        end = start + fold_size if fold < cv - 1 else n

        test_indices = set(range(start, end))
        train_X = [X_list[i] for i in range(n) if i not in test_indices]
        train_y = [y_list[i] for i in range(n) if i not in test_indices]
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
            correct = sum(1 for yt, yp in zip(test_y, preds_list) if round(yt) == round(yp))
            scores.append(correct / len(test_y))
        except Exception:
            scores.append(0.0)

    return sum(scores) / len(scores) if scores else 0.0


def stability_selection(
    X: Any,
    y: Any,
    n_bootstrap: int = 100,
    threshold: float = 0.6,
    alpha: float = 0.1,
    random_state: int | None = None,
) -> dict:
    """Stability selection via L1-regularized models on bootstrap subsamples.

    Runs L1-penalized regression on many random subsamples of the data and
    tracks how frequently each feature is selected (has a non-zero
    coefficient). Features that are consistently selected across subsamples
    are considered stable and important.

    Args:
        X: Feature matrix (n_samples x n_features).
        y: Target values (n_samples,).
        n_bootstrap: Number of bootstrap iterations.
        threshold: Minimum selection probability to include a feature.
        alpha: L1 regularization strength.
        random_state: Random seed for reproducibility.

    Returns:
        Dictionary containing:
            - selected_features: List of indices of stably selected features.
            - selection_probabilities: List of selection probability per
              feature.
            - threshold_used: The threshold that was applied.

    Raises:
        ValueError: If X and y have incompatible shapes.
    """
    n_samples, n_features = _get_shape(X)
    X_list = _to_2d_list(X)
    y_list = _to_1d_list(y)

    if len(y_list) != n_samples:
        raise ValueError(f"y length ({len(y_list)}) must match X rows ({n_samples})")

    if random_state is not None:
        random.seed(random_state)

    logger.info(
        "Stability selection: %d features, %d bootstrap, threshold=%.2f, alpha=%.4f",
        n_features,
        n_bootstrap,
        threshold,
        alpha,
    )

    selection_counts = [0] * n_features
    subsample_size = max(n_samples // 2, 2)

    for iteration in range(n_bootstrap):
        # Random subsample (without replacement)
        indices = random.sample(range(n_samples), subsample_size)
        sub_X = [X_list[i] for i in indices]
        sub_y = [y_list[i] for i in indices]

        # Fit L1-penalized model
        coefs = _fit_lasso_pure(sub_X, sub_y, alpha=alpha, max_iter=500)

        # Count non-zero features
        for feat in range(n_features):
            if abs(coefs[feat]) > 1e-8:
                selection_counts[feat] += 1

    # Compute selection probabilities
    selection_probabilities = [round(count / n_bootstrap, 4) for count in selection_counts]

    # Select features above threshold
    selected_features = [feat for feat in range(n_features) if selection_probabilities[feat] >= threshold]

    logger.info(
        "Stability selection: %d features selected (threshold=%.2f)",
        len(selected_features),
        threshold,
    )

    return {
        "selected_features": selected_features,
        "selection_probabilities": selection_probabilities,
        "threshold_used": threshold,
    }


def _fit_lasso_pure(
    X_list: list[list[float]],
    y_list: list[float],
    alpha: float = 0.1,
    max_iter: int = 500,
    tol: float = 1e-5,
) -> list[float]:
    """Fit lasso regression via coordinate descent (pure Python).

    Args:
        X_list: Feature matrix as list of lists.
        y_list: Target values.
        alpha: L1 regularization strength.
        max_iter: Maximum iterations.
        tol: Convergence tolerance.

    Returns:
        List of coefficients.
    """
    n = len(y_list)
    if n == 0 or not X_list or not X_list[0]:
        return []

    p = len(X_list[0])
    coefs = [0.0] * p

    # Center features and target
    means_x = [sum(X_list[i][j] for i in range(n)) / n for j in range(p)]
    mean_y = sum(y_list) / n

    x_centered = [[X_list[i][j] - means_x[j] for j in range(p)] for i in range(n)]
    y_centered = [y_list[i] - mean_y for i in range(n)]

    col_norms = [sum(x_centered[i][j] ** 2 for i in range(n)) for j in range(p)]

    for _iteration in range(max_iter):
        max_change = 0.0

        for j in range(p):
            residual = [y_centered[i] - sum(coefs[k] * x_centered[i][k] for k in range(p) if k != j) for i in range(n)]

            rho = sum(x_centered[i][j] * residual[i] for i in range(n))

            if col_norms[j] < 1e-15:
                new_coef = 0.0
            else:
                # Soft thresholding
                if rho > alpha * n:
                    new_coef = (rho - alpha * n) / col_norms[j]
                elif rho < -alpha * n:
                    new_coef = (rho + alpha * n) / col_norms[j]
                else:
                    new_coef = 0.0

            change = abs(new_coef - coefs[j])
            if change > max_change:
                max_change = change
            coefs[j] = new_coef

        if max_change < tol:
            break

    return coefs


def mutual_information_selection(
    X: Any,
    y: Any,
    n_features: int = 20,
    discrete_target: bool = True,
    n_bins: int = 10,
) -> dict:
    """Select features by mutual information with the target variable.

    Computes the mutual information between each feature and the target,
    then selects the top features with the highest MI scores.

    Args:
        X: Feature matrix (n_samples x n_features).
        y: Target values (n_samples,).
        n_features: Number of top features to select.
        discrete_target: Whether the target is discrete (classification)
            or continuous (regression).
        n_bins: Number of bins for MI estimation when using continuous
            variables.

    Returns:
        Dictionary containing:
            - selected_features: List of indices of selected features.
            - mi_scores: List of MI scores per feature.
            - ranking: List of feature rankings (1 = highest MI).

    Raises:
        ValueError: If X and y have incompatible shapes.
    """
    n_samples, total_features = _get_shape(X)
    X_list = _to_2d_list(X)
    y_list = _to_1d_list(y)

    if len(y_list) != n_samples:
        raise ValueError(f"y length ({len(y_list)}) must match X rows ({n_samples})")

    n_features = min(n_features, total_features)

    logger.info(
        "MI feature selection: %d features, selecting top %d, discrete=%s",
        total_features,
        n_features,
        discrete_target,
    )

    mi_scores: list[float] = []

    for feat in range(total_features):
        feat_values = [X_list[i][feat] for i in range(n_samples)]

        if discrete_target:
            mi = _mi_discrete_continuous(y_list, feat_values, n_bins=n_bins)
        else:
            mi = _mi_continuous(feat_values, y_list, n_bins=n_bins)

        mi_scores.append(round(mi, 8))

    # Rank features by MI (descending)
    indexed_scores = [(mi, idx) for idx, mi in enumerate(mi_scores)]
    indexed_scores.sort(key=lambda t: t[0], reverse=True)

    ranking = [0] * total_features
    for rank, (_, idx) in enumerate(indexed_scores):
        ranking[idx] = rank + 1

    selected_features = [idx for _, idx in indexed_scores[:n_features]]

    logger.info(
        "MI selection: top feature index %d (MI=%.6f)",
        selected_features[0] if selected_features else -1,
        indexed_scores[0][0] if indexed_scores else 0.0,
    )

    return {
        "selected_features": selected_features,
        "mi_scores": mi_scores,
        "ranking": ranking,
    }


def _mi_continuous(x: list[float], y: list[float], n_bins: int = 10) -> float:
    """Compute mutual information between two continuous variables.

    Args:
        x: First variable.
        y: Second variable.
        n_bins: Number of bins for discretization.

    Returns:
        Mutual information in nats.
    """
    n = len(x)
    if n == 0:
        return 0.0

    min_x, max_x = min(x), max(x)
    min_y, max_y = min(y), max(y)
    range_x = max_x - min_x if max_x > min_x else 1.0
    range_y = max_y - min_y if max_y > min_y else 1.0

    def _bin(val: float, vmin: float, vrange: float) -> int:
        idx = int((val - vmin) / vrange * n_bins)
        return max(0, min(n_bins - 1, idx))

    joint: dict[tuple[int, int], int] = defaultdict(int)
    mx: dict[int, int] = defaultdict(int)
    my: dict[int, int] = defaultdict(int)

    for xi, yi in zip(x, y):
        bx = _bin(xi, min_x, range_x)
        by = _bin(yi, min_y, range_y)
        joint[(bx, by)] += 1
        mx[bx] += 1
        my[by] += 1

    mi = 0.0
    for (bx, by), count in joint.items():
        p_xy = count / n
        p_x = mx[bx] / n
        p_y = my[by] / n
        if p_xy > 0 and p_x > 0 and p_y > 0:
            mi += p_xy * math.log(p_xy / (p_x * p_y))

    return max(0.0, mi)


def _mi_discrete_continuous(
    discrete: list[float],
    continuous: list[float],
    n_bins: int = 10,
) -> float:
    """Compute MI between a discrete variable and a continuous variable.

    Args:
        discrete: Discrete variable (class labels).
        continuous: Continuous variable (feature values).
        n_bins: Number of bins for continuous variable.

    Returns:
        Mutual information in nats.
    """
    n = len(discrete)
    if n == 0:
        return 0.0

    # Group continuous values by class
    classes: dict[float, list[float]] = defaultdict(list)
    for d, c in zip(discrete, continuous):
        classes[d].append(c)

    # Bin the continuous variable
    min_c, max_c = min(continuous), max(continuous)
    range_c = max_c - min_c if max_c > min_c else 1.0

    def _bin(val: float) -> int:
        idx = int((val - min_c) / range_c * n_bins)
        return max(0, min(n_bins - 1, idx))

    # Compute joint and marginal distributions
    joint: dict[tuple[float, int], int] = defaultdict(int)
    margin_class: dict[float, int] = defaultdict(int)
    margin_bin: dict[int, int] = defaultdict(int)

    for d, c in zip(discrete, continuous):
        b = _bin(c)
        joint[(d, b)] += 1
        margin_class[d] += 1
        margin_bin[b] += 1

    mi = 0.0
    for (d, b), count in joint.items():
        p_xy = count / n
        p_x = margin_class[d] / n
        p_y = margin_bin[b] / n
        if p_xy > 0 and p_x > 0 and p_y > 0:
            mi += p_xy * math.log(p_xy / (p_x * p_y))

    return max(0.0, mi)
