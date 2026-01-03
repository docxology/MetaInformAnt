"""Feature selection utilities for METAINFORMANT.

This module provides comprehensive feature selection methods specifically
designed for biological data, including univariate selection, recursive
elimination, and stability-based selection.
"""

from __future__ import annotations

import numpy as np
from typing import Any, Dict, List, Optional, Tuple, Union

from metainformant.core.logging import get_logger

logger = get_logger(__name__)

# Optional imports for ML functionality
try:
    from sklearn.base import BaseEstimator
    from sklearn.feature_selection import (
        SelectKBest, SelectPercentile, SelectFdr, SelectFpr,
        RFE, RFECV, SelectFromModel
    )
    from sklearn.ensemble import RandomForestClassifier, RandomForestRegressor
    from sklearn.linear_model import Lasso, LogisticRegression
    from sklearn.preprocessing import StandardScaler
    from sklearn.model_selection import cross_val_score, StratifiedKFold
    HAS_SKLEARN = True
except ImportError:
    HAS_SKLEARN = False
    logger.warning("scikit-learn not available, feature selection disabled")


def select_features_univariate(
    X: np.ndarray,
    y: np.ndarray,
    method: str = "f_classif",
    k: int | str = "all",
    **kwargs: Any
) -> Tuple[np.ndarray, np.ndarray]:
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

    # Import statistical tests
    if method == "f_classif":
        from sklearn.feature_selection import f_classif
        score_func = f_classif
    elif method == "chi2":
        from sklearn.feature_selection import chi2
        score_func = chi2
    elif method == "mutual_info_classif":
        from sklearn.feature_selection import mutual_info_classif
        score_func = mutual_info_classif
    elif method == "f_regression":
        from sklearn.feature_selection import f_regression
        score_func = f_regression
    elif method == "mutual_info_regression":
        from sklearn.feature_selection import mutual_info_regression
        score_func = mutual_info_regression
    else:
        raise ValueError(f"Unknown univariate method: {method}")

    # Create selector
    if k == "all":
        selector = SelectKBest(score_func=score_func, k="all", **kwargs)
    elif isinstance(k, int):
        selector = SelectKBest(score_func=score_func, k=k, **kwargs)
    elif isinstance(k, str) and k.endswith('%'):
        # Percentage-based selection
        percentile = int(k[:-1])
        selector = SelectPercentile(score_func=score_func, percentile=percentile, **kwargs)
    else:
        raise ValueError(f"Invalid k parameter: {k}")

    # Fit and transform
    X_selected = selector.fit_transform(X, y)
    selected_indices = selector.get_support(indices=True)

    logger.info(
        f"Selected {len(selected_indices)} features using {method} "
        f"(k={k} from {X.shape[1]} total)"
    )

    return X_selected, selected_indices


def select_features_recursive(
    X: np.ndarray,
    y: np.ndarray,
    estimator: Any = None,
    n_features_to_select: int | None = None,
    step: float = 0.1,
    cv: int = 5,
    **kwargs: Any
) -> Tuple[np.ndarray, np.ndarray]:
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

    # Default estimator
    if estimator is None:
        estimator = RandomForestClassifier(n_estimators=100, random_state=42)

    # Use RFECV if no specific number requested
    if n_features_to_select is None:
        selector = RFECV(
            estimator=estimator,
            step=step,
            cv=cv,
            min_features_to_select=1,
            **kwargs
        )
    else:
        selector = RFE(
            estimator=estimator,
            n_features_to_select=n_features_to_select,
            step=step,
            **kwargs
        )

    # Fit and transform
    X_selected = selector.fit_transform(X, y)
    selected_indices = selector.get_support(indices=True)

    logger.info(
        f"Selected {len(selected_indices)} features using recursive elimination "
        f"(from {X.shape[1]} total)"
    )

    return X_selected, selected_indices


def select_features_stability(
    X: np.ndarray,
    y: np.ndarray,
    method: str = "rf",
    n_bootstraps: int = 100,
    threshold: float = 0.5,
    random_state: int | None = None,
    **kwargs: Any
) -> Tuple[np.ndarray, np.ndarray]:
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

    np.random.seed(random_state)
    n_samples, n_features = X.shape

    # Track feature selection frequency
    selection_counts = np.zeros(n_features)

    # Perform bootstrap iterations
    for i in range(n_bootstraps):
        # Bootstrap sample
        indices = np.random.choice(n_samples, size=n_samples, replace=True)
        X_boot = X[indices]
        y_boot = y[indices]

        try:
            # Apply base selection method
            if method == "rf":
                selector = SelectFromModel(
                    RandomForestClassifier(
                        n_estimators=50,
                        random_state=random_state,
                        **kwargs
                    ),
                    threshold="median"
                )
            elif method == "lasso":
                selector = SelectFromModel(
                    LogisticRegression(
                        penalty='l1',
                        solver='liblinear',
                        random_state=random_state,
                        **kwargs
                    )
                )
            elif method == "univariate":
                selector = SelectKBest(
                    k=int(n_features * 0.5),  # Select 50% in each iteration
                    **kwargs
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
    X_selected = X[:, selected_indices]

    logger.info(
        f"Selected {len(selected_indices)} stable features "
        f"(threshold={threshold}, bootstraps={n_bootstraps})"
    )

    return X_selected, selected_indices


def biological_feature_ranking(
    X: np.ndarray,
    y: np.ndarray,
    feature_names: List[str] | None = None,
    method: str = "importance",
    **kwargs: Any
) -> Dict[str, Any]:
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

    if feature_names is None:
        feature_names = [f"feature_{i}" for i in range(X.shape[1])]

    results = {
        'method': method,
        'n_features': X.shape[1],
        'feature_names': feature_names,
    }

    if method == "importance":
        # Use Random Forest feature importance
        rf = RandomForestClassifier(
            n_estimators=100,
            random_state=42,
            **kwargs
        )
        rf.fit(X, y)

        importances = rf.feature_importances_
        sorted_indices = np.argsort(importances)[::-1]

        results['rankings'] = [
            {
                'feature': feature_names[i],
                'importance': float(importances[i]),
                'rank': rank + 1
            }
            for rank, i in enumerate(sorted_indices)
        ]

    elif method == "univariate":
        # Use F-test for classification
        from sklearn.feature_selection import f_classif
        f_scores, p_values = f_classif(X, y)

        sorted_indices = np.argsort(f_scores)[::-1]

        results['rankings'] = [
            {
                'feature': feature_names[i],
                'f_score': float(f_scores[i]),
                'p_value': float(p_values[i]),
                'rank': rank + 1
            }
            for rank, i in enumerate(sorted_indices)
        ]

    elif method == "stability":
        # Use stability selection
        _, selected_indices = select_features_stability(X, y, **kwargs)

        results['rankings'] = [
            {
                'feature': feature_names[i],
                'selected': True,
                'rank': rank + 1
            }
            for rank, i in enumerate(selected_indices)
        ]

    else:
        raise ValueError(f"Unknown ranking method: {method}")

    logger.info(f"Ranked {len(results['rankings'])} features using {method} method")
    return results


def select_features_biological(
    X: np.ndarray,
    y: np.ndarray,
    methods: List[str] = None,
    feature_names: List[str] | None = None,
    **kwargs: Any
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
        methods = ['univariate', 'recursive', 'stability']

    if feature_names is None:
        feature_names = [f"feature_{i}" for i in range(X.shape[1])]

    results = {
        'input_shape': X.shape,
        'methods': {},
        'consensus': {},
    }

    # Apply each method
    for method in methods:
        try:
            if method == 'univariate':
                X_sel, indices = select_features_univariate(
                    X, y, k='all', **kwargs.get('univariate', {})
                )
                selected_features = [feature_names[i] for i in indices]

            elif method == 'recursive':
                X_sel, indices = select_features_recursive(
                    X, y, **kwargs.get('recursive', {})
                )
                selected_features = [feature_names[i] for i in indices]

            elif method == 'stability':
                X_sel, indices = select_features_stability(
                    X, y, **kwargs.get('stability', {})
                )
                selected_features = [feature_names[i] for i in indices]

            else:
                raise ValueError(f"Unknown method: {method}")

            results['methods'][method] = {
                'selected_features': selected_features,
                'n_selected': len(selected_features),
                'selected_indices': indices.tolist(),
            }

        except Exception as e:
            logger.error(f"Feature selection method {method} failed: {e}")
            results['methods'][method] = {'error': str(e)}

    # Calculate consensus (features selected by multiple methods)
    if len(methods) > 1:
        all_selected = []
        for method_result in results['methods'].values():
            if 'selected_features' in method_result:
                all_selected.extend(method_result['selected_features'])

        from collections import Counter
        consensus_counts = Counter(all_selected)

        consensus_features = [
            feature for feature, count in consensus_counts.items()
            if count >= 2  # Selected by at least 2 methods
        ]

        results['consensus'] = {
            'features': consensus_features,
            'n_features': len(consensus_features),
            'selection_counts': dict(consensus_counts),
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
    clf = RandomForestClassifier(
        n_estimators=100,
        random_state=random_state
    )

    # Cross-validation
    cv = StratifiedKFold(n_splits=cv_folds, shuffle=True, random_state=random_state)
    scores = cross_val_score(clf, X_selected, y, cv=cv, scoring='accuracy')

    return {
        'n_features_selected': len(selected_indices),
        'n_features_total': X.shape[1],
        'selection_ratio': len(selected_indices) / X.shape[1],
        'cv_accuracy_mean': float(scores.mean()),
        'cv_accuracy_std': float(scores.std()),
        'cv_scores': scores.tolist(),
    }





