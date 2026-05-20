"""Feature selection utilities for METAINFORMANT.

This module provides comprehensive feature selection methods specifically
designed for biological data, including univariate selection, recursive
elimination, and stability-based selection.
"""

from __future__ import annotations

import warnings
from typing import Any, Dict, List, Tuple

import numpy as np

from metainformant.core.utils.logging import get_logger

logger = get_logger(__name__)

# Optional imports for ML functionality
try:
    from sklearn.ensemble import RandomForestClassifier
    from sklearn.feature_selection import (
        RFE,
        RFECV,
        SelectFromModel,
        SelectKBest,
        SelectPercentile,
    )
    from sklearn.linear_model import LogisticRegression
    from sklearn.model_selection import StratifiedKFold, cross_val_score

    HAS_SKLEARN = True
except ImportError:
    HAS_SKLEARN = False
    logger.warning("scikit-learn not available, feature selection disabled")


def select_features_univariate(
    X: np.ndarray, y: np.ndarray, method: str = "f_classif", k: int | str = "all", **kwargs: Any
) -> Tuple[np.ndarray, list[int]]:
    """Select features using univariate statistical tests.

    Args:
        X: Feature matrix (n_samples, n_features)
        y: Target labels (n_samples,)
        method: Statistical test method
            - 'f_classif': ANOVA F-test for classification
            - 'chi2': Chi-squared test
            - 'mutual_info_classif': Mutual information for classification
            - 'f_regression': F-test for regression
            - 'mutual_info_regression': Mutual information for regression
        k: Number of features to select ('all' for no selection)
        **kwargs: Additional parameters for the selector

    Returns:
        Tuple of (X_selected, selected_indices)

    Raises:
        ImportError: If scikit-learn not available
        ValueError: If invalid method or parameters
    """
    if not HAS_SKLEARN:
        raise ImportError("scikit-learn required for feature selection")

    X = np.asarray(X)
    y = np.asarray(y)
    if X.ndim != 2:
        raise ValueError("X must be 2D array")

    # Validate input dimensions
    if len(y) != X.shape[0]:
        raise ValueError(f"y length must match X samples: {len(y)} != {X.shape[0]}")

    # Import statistical tests
    if method in ("f_classif", "f_score"):
        from sklearn.feature_selection import f_classif

        score_func = f_classif
    elif method in ("chi2", "chi_squared"):
        from sklearn.feature_selection import chi2

        score_func = chi2
    elif method in ("mutual_info_classif", "mutual_info"):
        from sklearn.feature_selection import mutual_info_classif

        score_func = mutual_info_classif
    elif method in ("f_regression", "f_score_regression"):
        from sklearn.feature_selection import f_regression

        score_func = f_regression
    elif method == "mutual_info_regression":
        from sklearn.feature_selection import mutual_info_regression

        score_func = mutual_info_regression
    else:
        raise ValueError(f"Unknown univariate method: {method}")

    selector_kwargs = dict(kwargs)
    selector_kwargs.pop("p_threshold", None)

    # Create selector
    if k == "all":
        selector = SelectKBest(score_func=score_func, k="all", **selector_kwargs)
    elif isinstance(k, int):
        selector = SelectKBest(score_func=score_func, k=min(k, X.shape[1]), **selector_kwargs)
    elif isinstance(k, str) and k.endswith("%"):
        # Percentage-based selection
        percentile = int(k[:-1])
        selector = SelectPercentile(score_func=score_func, percentile=percentile, **selector_kwargs)
    else:
        raise ValueError(f"Invalid k parameter: {k}")

    # Fit and transform
    X_selected = selector.fit_transform(X, y)
    selected_indices = selector.get_support(indices=True).tolist()

    logger.info(f"Selected {len(selected_indices)} features using {method} " f"(k={k} from {X.shape[1]} total)")

    return X_selected, selected_indices


def select_features_recursive(
    X: np.ndarray,
    y: np.ndarray,
    estimator: Any = None,
    n_features_to_select: int | None = None,
    step: float = 0.1,
    cv: int = 5,
    **kwargs: Any,
) -> Tuple[np.ndarray, list[int]]:
    """Select features using recursive elimination.

    Args:
        X: Feature matrix (n_samples, n_features)
        y: Target labels (n_samples,)
        estimator: Base estimator (default: RandomForest)
        n_features_to_select: Number of features to select
        step: Step size for elimination (int or float)
        cv: Cross-validation folds for RFECV
        **kwargs: Additional parameters

    Returns:
        Tuple of (X_selected, selected_indices)

    Raises:
        ImportError: If scikit-learn not available
    """
    if not HAS_SKLEARN:
        raise ImportError("scikit-learn required for recursive feature selection")

    X = np.asarray(X)
    y = np.asarray(y)
    if X.ndim != 2:
        raise ValueError("X must be 2D array")
    if len(y) != X.shape[0]:
        raise ValueError(f"y length must match X samples: {len(y)} != {X.shape[0]}")

    estimator_type = kwargs.pop("estimator_type", None)
    if "n_features" in kwargs:
        n_features_to_select = kwargs.pop("n_features")
    random_state = kwargs.pop("random_state", 42)

    # Default estimator
    if estimator is None:
        if estimator_type in (None, "random_forest", "rf"):
            estimator = RandomForestClassifier(n_estimators=100, random_state=random_state)
        elif estimator_type in ("linear", "logistic"):
            estimator = LogisticRegression(max_iter=1000, solver="liblinear", random_state=random_state)
        elif estimator_type == "svm":
            from sklearn.svm import LinearSVC

            estimator = LinearSVC(dual=False, max_iter=5000, random_state=random_state)
        else:
            raise ValueError(f"Unknown importance method: {estimator_type}")

    if n_features_to_select is not None:
        n_features_to_select = min(int(n_features_to_select), X.shape[1])

    # Use RFECV if no specific number requested
    if n_features_to_select is None:
        selector = RFECV(estimator=estimator, step=step, cv=cv, min_features_to_select=1, **kwargs)
    else:
        selector = RFE(estimator=estimator, n_features_to_select=n_features_to_select, step=step, **kwargs)

    # Fit and transform
    X_selected = selector.fit_transform(X, y)
    selected_indices = selector.get_support(indices=True).tolist()

    logger.info(f"Selected {len(selected_indices)} features using recursive elimination " f"(from {X.shape[1]} total)")

    return X_selected, selected_indices


def select_features_stability(
    X: np.ndarray,
    y: np.ndarray,
    method: str = "rf",
    n_bootstraps: int = 100,
    threshold: float = 0.5,
    random_state: int | None = None,
    **kwargs: Any,
) -> Tuple[np.ndarray, list[int]]:
    """Select features using stability-based selection.

    This method runs feature selection multiple times on bootstrapped samples
    and selects features that are consistently selected above a threshold.

    Args:
        X: Feature matrix (n_samples, n_features)
        y: Target labels (n_samples,)
        method: Base selection method ('rf', 'lasso', 'univariate')
        n_bootstraps: Number of bootstrap iterations
        threshold: Selection frequency threshold (0-1)
        random_state: Random state for reproducibility
        **kwargs: Additional parameters for base method

    Returns:
        Tuple of (X_selected, selected_indices)

    Raises:
        ImportError: If scikit-learn not available
    """
    if not HAS_SKLEARN:
        raise ImportError("scikit-learn required for stability selection")

    if "n_bootstrap" in kwargs:
        n_bootstraps = int(kwargs.pop("n_bootstrap"))
    if "stability_threshold" in kwargs:
        threshold = float(kwargs.pop("stability_threshold"))
    subsample_ratio = float(kwargs.pop("subsample_ratio", 1.0))

    X = np.asarray(X)
    y = np.asarray(y)
    if X.ndim != 2:
        raise ValueError("X must be 2D array")
    if len(y) != X.shape[0]:
        raise ValueError(f"y length must match X samples: {len(y)} != {X.shape[0]}")

    rng = np.random.default_rng(random_state)
    n_samples, n_features = X.shape
    bootstrap_size = max(1, min(n_samples, int(round(n_samples * subsample_ratio))))

    # Track feature selection frequency
    selection_counts = np.zeros(n_features)

    # Perform bootstrap iterations
    for i in range(n_bootstraps):
        # Bootstrap sample
        indices = rng.choice(n_samples, size=bootstrap_size, replace=True)
        X_boot = X[indices]
        y_boot = y[indices]

        try:
            # Apply base selection method
            if method in ("rf", "random_forest"):
                selector = SelectFromModel(
                    RandomForestClassifier(n_estimators=50, random_state=random_state, **kwargs), threshold="median"
                )
            elif method == "lasso":
                selector = SelectFromModel(
                    LogisticRegression(penalty="l1", solver="liblinear", random_state=random_state, **kwargs)
                )
            elif method == "univariate":
                from sklearn.feature_selection import f_classif

                selector = SelectKBest(
                    score_func=f_classif,
                    k=max(1, int(n_features * 0.5)),
                )
            else:
                raise ValueError(f"Unknown stability method: {method}")

            selector.fit(X_boot, y_boot)
            selected = selector.get_support()
            selection_counts += selected.astype(int)

        except Exception as e:
            logger.warning(f"Bootstrap iteration {i} failed: {e}")
            continue

    # Calculate selection frequency
    selection_freq = selection_counts / n_bootstraps

    # Select features above threshold
    selected_indices = np.where(selection_freq >= threshold)[0]
    if len(selected_indices) == 0:
        warnings.warn(
            "No features meet stability threshold; returning top-ranked features instead",
            UserWarning,
            stacklevel=2,
        )
        n_fallback = max(1, min(n_features, int(np.ceil(n_features * 0.2))))
        selected_indices = np.argsort(selection_freq)[::-1][:n_fallback]
        selected_indices = np.sort(selected_indices)

    X_selected = X[:, selected_indices]
    selected_indices_list = selected_indices.tolist()

    logger.info(
        f"Selected {len(selected_indices)} stable features " f"(threshold={threshold}, bootstraps={n_bootstraps})"
    )

    return X_selected, selected_indices_list


def biological_feature_ranking(
    X: np.ndarray, y: np.ndarray, feature_names: List[str] | None = None, method: str = "importance", **kwargs: Any
) -> Dict[str, Any] | list[tuple[int, str, float]]:
    """Rank features for biological interpretation.

    Args:
        X: Feature matrix (n_samples, n_features)
        y: Target labels (n_samples,)
        feature_names: Optional feature names
        method: Ranking method ('importance', 'univariate', 'stability')
        **kwargs: Additional parameters

    Returns:
        Dictionary with ranking results

    Raises:
        ImportError: If scikit-learn not available
    """
    if not HAS_SKLEARN:
        raise ImportError("scikit-learn required for feature ranking")

    X = np.asarray(X)
    y = np.asarray(y)
    if X.ndim != 2:
        raise ValueError("X must be 2D array")

    if feature_names is None:
        feature_names = [f"feature_{i}" for i in range(X.shape[1])]
    if len(feature_names) != X.shape[1]:
        raise ValueError("feature_names length must match X columns")

    if method in ("statistical", "biological", "combined"):
        from sklearn.feature_selection import f_classif

        f_scores, _ = f_classif(X, y)
        statistical_scores = np.nan_to_num(f_scores, nan=0.0, posinf=0.0, neginf=0.0)
        gene_weights = kwargs.get("gene_weights") or {}
        biological_scores = np.array([float(gene_weights.get(name, 1.0)) for name in feature_names])

        if method == "statistical":
            scores = statistical_scores
        elif method == "biological":
            scores = biological_scores
        else:
            stat_max = float(np.max(statistical_scores)) if len(statistical_scores) else 0.0
            bio_max = float(np.max(biological_scores)) if len(biological_scores) else 0.0
            stat_norm = statistical_scores / stat_max if stat_max > 0 else statistical_scores
            bio_norm = biological_scores / bio_max if bio_max > 0 else biological_scores
            scores = stat_norm + bio_norm

        ranked_indices = np.argsort(scores)[::-1]
        return [(int(i), feature_names[int(i)], float(scores[int(i)])) for i in ranked_indices]

    results = {
        "method": method,
        "n_features": X.shape[1],
        "feature_names": feature_names,
    }

    if method == "importance":
        # Use Random Forest feature importance
        rf = RandomForestClassifier(n_estimators=100, random_state=42, **kwargs)
        rf.fit(X, y)

        importances = rf.feature_importances_
        sorted_indices = np.argsort(importances)[::-1]

        results["rankings"] = [
            {"feature": feature_names[i], "importance": float(importances[i]), "rank": rank + 1}
            for rank, i in enumerate(sorted_indices)
        ]

    elif method == "univariate":
        # Use F-test for classification
        from sklearn.feature_selection import f_classif

        f_scores, p_values = f_classif(X, y)

        sorted_indices = np.argsort(f_scores)[::-1]

        results["rankings"] = [
            {
                "feature": feature_names[i],
                "f_score": float(f_scores[i]),
                "p_value": float(p_values[i]),
                "rank": rank + 1,
            }
            for rank, i in enumerate(sorted_indices)
        ]

    elif method == "stability":
        # Use stability selection
        _, selected_indices = select_features_stability(X, y, **kwargs)

        results["rankings"] = [
            {"feature": feature_names[i], "selected": True, "rank": rank + 1} for rank, i in enumerate(selected_indices)
        ]

    else:
        raise ValueError(f"Unknown method: {method}")

    logger.info(f"Ranked {len(results['rankings'])} features using {method} method")
    return results


def select_features_biological(
    X: np.ndarray, y: np.ndarray, methods: List[str] = None, feature_names: List[str] | None = None, **kwargs: Any
) -> Dict[str, Any]:
    """Comprehensive biological feature selection.

    This function applies multiple feature selection methods and
    provides a consensus ranking.

    Args:
        X: Feature matrix
        y: Target labels
        methods: List of methods to apply
        feature_names: Optional feature names
        **kwargs: Parameters for individual methods

    Returns:
        Dictionary with selection results from all methods
    """
    if methods is None:
        methods = ["univariate", "recursive", "stability"]

    if feature_names is None:
        feature_names = [f"feature_{i}" for i in range(X.shape[1])]

    results = {
        "input_shape": X.shape,
        "methods": {},
        "consensus": {},
    }

    # Apply each method
    for method in methods:
        try:
            if method == "univariate":
                X_sel, indices = select_features_univariate(X, y, k="all", **kwargs.get("univariate", {}))
                selected_features = [feature_names[i] for i in indices]

            elif method == "recursive":
                X_sel, indices = select_features_recursive(X, y, **kwargs.get("recursive", {}))
                selected_features = [feature_names[i] for i in indices]

            elif method == "stability":
                X_sel, indices = select_features_stability(X, y, **kwargs.get("stability", {}))
                selected_features = [feature_names[i] for i in indices]

            else:
                raise ValueError(f"Unknown method: {method}")

            results["methods"][method] = {
                "selected_features": selected_features,
                "n_selected": len(selected_features),
                "selected_indices": indices if isinstance(indices, list) else indices.tolist(),
            }

        except Exception as e:
            logger.error(f"Feature selection method {method} failed: {e}")
            results["methods"][method] = {"error": str(e)}

    # Calculate consensus (features selected by multiple methods)
    if len(methods) > 1:
        all_selected = []
        for method_result in results["methods"].values():
            if "selected_features" in method_result:
                all_selected.extend(method_result["selected_features"])

        from collections import Counter

        consensus_counts = Counter(all_selected)

        consensus_features = [
            feature for feature, count in consensus_counts.items() if count >= 2  # Selected by at least 2 methods
        ]

        results["consensus"] = {
            "features": consensus_features,
            "n_features": len(consensus_features),
            "selection_counts": dict(consensus_counts),
        }

    return results


def evaluate_feature_selection(
    X: np.ndarray,
    y: np.ndarray,
    selected_indices: np.ndarray,
    cv_folds: int = 5,
    random_state: int | None = None,
) -> Dict[str, Any]:
    """Evaluate the quality of feature selection.

    Args:
        X: Full feature matrix
        y: Target labels
        selected_indices: Indices of selected features
        cv_folds: Number of cross-validation folds
        random_state: Random state for reproducibility

    Returns:
        Dictionary with evaluation results
    """
    if not HAS_SKLEARN:
        raise ImportError("scikit-learn required for evaluation")

    X_selected = X[:, selected_indices]

    # Train classifier on selected features
    clf = RandomForestClassifier(n_estimators=100, random_state=random_state)

    # Cross-validation
    cv = StratifiedKFold(n_splits=cv_folds, shuffle=True, random_state=random_state)
    scores = cross_val_score(clf, X_selected, y, cv=cv, scoring="accuracy")

    return {
        "n_features_selected": len(selected_indices),
        "n_features_total": X.shape[1],
        "selection_ratio": len(selected_indices) / X.shape[1],
        "cv_accuracy_mean": float(scores.mean()),
        "cv_accuracy_std": float(scores.std()),
        "cv_scores": scores.tolist(),
    }
