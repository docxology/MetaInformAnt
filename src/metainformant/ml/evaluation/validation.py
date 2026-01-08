"""Model validation and cross-validation utilities for METAINFORMANT.

This module provides comprehensive validation tools specifically designed
for biological data, including stratified splitting, cross-validation,
and permutation importance analysis.
"""

from __future__ import annotations

import numpy as np
from typing import Any, Dict, List, Optional, Tuple, Union

from metainformant.core.utils.logging import get_logger

logger = get_logger(__name__)

# Optional imports for ML functionality
try:
    from sklearn.base import BaseEstimator
    from sklearn.model_selection import (
        train_test_split, StratifiedKFold, KFold,
        cross_val_score, cross_validate
    )
    from sklearn.metrics import (
        accuracy_score, precision_score, recall_score, f1_score,
        roc_auc_score, mean_squared_error, r2_score
    )
    from sklearn.utils import resample
    HAS_SKLEARN = True
except ImportError:
    HAS_SKLEARN = False
    logger.warning("scikit-learn not available, ML validation disabled")


def train_test_split_biological(
    X: np.ndarray,
    y: np.ndarray,
    test_size: float = 0.2,
    stratify: np.ndarray | None = None,
    random_state: int | None = None,
    **kwargs: Any
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Split biological data into train/test sets with stratification.

    This function ensures balanced class distribution in train/test splits,
    which is crucial for biological classification tasks where classes
    may be imbalanced.

    Args:
        X: Feature matrix
        y: Target labels
        test_size: Proportion of data for testing (0-1)
        stratify: Array to stratify by (usually y for classification)
        random_state: Random state for reproducibility
        **kwargs: Additional parameters for train_test_split

    Returns:
        Tuple of (X_train, X_test, y_train, y_test)

    Raises:
        ImportError: If scikit-learn not available
    """
    if not HAS_SKLEARN:
        raise ImportError("scikit-learn required for data splitting")

    # Default stratification for biological data
    if stratify is None and len(np.unique(y)) < 20:  # Classification-like
        stratify = y

    X_train, X_test, y_train, y_test = train_test_split(
        X, y,
        test_size=test_size,
        stratify=stratify,
        random_state=random_state,
        **kwargs
    )

    logger.info(
        f"Split biological data: {len(X_train)} train, {len(X_test)} test samples "
        f"({len(np.unique(y_train))} classes in train, {len(np.unique(y_test))} in test)"
    )

    return X_train, X_test, y_train, y_test


def cross_validation_scores(
    model: Any,
    X: np.ndarray,
    y: np.ndarray,
    cv: int = 5,
    scoring: str | List[str] = "accuracy",
    stratified: bool = True,
    random_state: int | None = None,
) -> Dict[str, np.ndarray]:
    """Perform cross-validation with multiple scoring metrics.

    Args:
        model: Trained sklearn model
        X: Feature matrix
        y: Target labels/values
        cv: Number of cross-validation folds
        scoring: Scoring metric(s) to use
        stratified: Whether to use stratified folds for classification
        random_state: Random state for reproducibility

    Returns:
        Dictionary mapping metric names to score arrays

    Raises:
        ImportError: If scikit-learn not available
    """
    if not HAS_SKLEARN:
        raise ImportError("scikit-learn required for cross-validation")

    # Determine cross-validation strategy
    if stratified and len(np.unique(y)) < 20:  # Classification
        cv_strategy = StratifiedKFold(n_splits=cv, shuffle=True, random_state=random_state)
    else:
        cv_strategy = KFold(n_splits=cv, shuffle=True, random_state=random_state)

    # Handle single or multiple scoring metrics
    if isinstance(scoring, str):
        scoring = [scoring]

    results = {}

    for metric in scoring:
        try:
            scores = cross_val_score(
                model, X, y,
                cv=cv_strategy,
                scoring=metric
            )
            results[metric] = scores

            logger.debug(
                f"{metric}: {scores.mean():.3f} ± {scores.std():.3f} "
                f"({cv}-fold CV)"
            )

        except Exception as e:
            logger.warning(f"Could not compute {metric}: {e}")
            results[metric] = np.array([])

    return results


def permutation_importance_biological(
    model: Any,
    X: np.ndarray,
    y: np.ndarray,
    n_repeats: int = 10,
    random_state: int | None = None,
    scoring: str = "accuracy",
) -> Dict[str, Any]:
    """Calculate permutation feature importance for biological models.

    This method assesses feature importance by measuring how much
    model performance degrades when each feature is randomly permuted.

    Args:
        model: Trained sklearn model
        X: Feature matrix
        y: Target labels/values
        n_repeats: Number of permutation repeats per feature
        random_state: Random state for reproducibility
        scoring: Scoring metric to use

    Returns:
        Dictionary with importance results

    Raises:
        ImportError: If scikit-learn not available
    """
    if not HAS_SKLEARN:
        raise ImportError("scikit-learn required for permutation importance")

    np.random.seed(random_state)

    # Get baseline score
    baseline_score = cross_val_score(
        model, X, y,
        cv=StratifiedKFold(n_splits=5, shuffle=True, random_state=random_state),
        scoring=scoring
    ).mean()

    n_features = X.shape[1]
    importance_scores = np.zeros((n_features, n_repeats))
    feature_names = getattr(X, 'columns', None)

    if feature_names is None:
        feature_names = [f"feature_{i}" for i in range(n_features)]

    # Calculate importance for each feature
    for i in range(n_features):
        for j in range(n_repeats):
            # Create permuted version of this feature
            X_permuted = X.copy()
            perm_indices = np.random.permutation(len(X))
            X_permuted[:, i] = X[perm_indices, i]

            # Calculate score with permuted feature
            permuted_score = cross_val_score(
                model, X_permuted, y,
                cv=StratifiedKFold(n_splits=5, shuffle=True, random_state=random_state),
                scoring=scoring
            ).mean()

            # Importance is baseline - permuted score
            importance_scores[i, j] = baseline_score - permuted_score

    # Calculate statistics
    mean_importances = np.mean(importance_scores, axis=1)
    std_importances = np.std(importance_scores, axis=1)

    # Sort by importance
    sorted_indices = np.argsort(mean_importances)[::-1]

    results = {
        'baseline_score': float(baseline_score),
        'scoring_metric': scoring,
        'n_repeats': n_repeats,
        'feature_importances': [
            {
                'feature': feature_names[idx],
                'importance_mean': float(mean_importances[idx]),
                'importance_std': float(std_importances[idx]),
                'rank': rank + 1,
            }
            for rank, idx in enumerate(sorted_indices)
        ],
        'raw_scores': importance_scores.tolist(),
    }

    logger.info(
        f"Calculated permutation importance for {n_features} features "
        f"({n_repeats} repeats, {scoring} scoring)"
    )

    return results


def validate_model_stability(
    X: np.ndarray,
    y: np.ndarray,
    model_factory: callable,
    n_bootstraps: int = 50,
    test_size: float = 0.2,
    random_state: int | None = None,
    **kwargs: Any
) -> Dict[str, Any]:
    """Validate model stability using bootstrapping.

    This function assesses how stable model performance is across
    different bootstrap samples of the training data.

    Args:
        X: Feature matrix
        y: Target labels/values
        model_factory: Function that returns a new model instance
        n_bootstraps: Number of bootstrap iterations
        test_size: Proportion of data for testing
        random_state: Random state for reproducibility
        **kwargs: Additional parameters

    Returns:
        Dictionary with stability analysis results

    Raises:
        ImportError: If scikit-learn not available
    """
    if not HAS_SKLEARN:
        raise ImportError("scikit-learn required for stability validation")

    np.random.seed(random_state)
    n_samples = len(X)

    scores = []
    feature_importances = []

    for i in range(n_bootstraps):
        # Bootstrap sample
        indices = np.random.choice(n_samples, size=n_samples, replace=True)
        X_boot = X[indices]
        y_boot = y[indices]

        # Split into train/test
        X_train, X_test, y_train, y_test = train_test_split_biological(
            X_boot, y_boot, test_size=test_size, random_state=random_state
        )

        # Train model
        model = model_factory(**kwargs)
        model.fit(X_train, y_train)

        # Evaluate
        y_pred = model.predict(X_test)

        if len(np.unique(y)) < 20:  # Classification
            score = accuracy_score(y_test, y_pred)
        else:  # Regression
            score = r2_score(y_test, y_pred)

        scores.append(score)

        # Collect feature importances if available
        if hasattr(model, 'feature_importances_'):
            feature_importances.append(model.feature_importances_)
        elif hasattr(model, 'coef_'):
            feature_importances.append(np.abs(model.coef_))

    # Calculate stability metrics
    scores_array = np.array(scores)

    results = {
        'n_bootstraps': n_bootstraps,
        'scores': {
            'mean': float(scores_array.mean()),
            'std': float(scores_array.std()),
            'min': float(scores_array.min()),
            'max': float(scores_array.max()),
            'values': scores_array.tolist(),
        },
        'stability_metrics': {
            'coefficient_of_variation': float(scores_array.std() / scores_array.mean()),
            'stability_score': 1.0 - float(scores_array.std() / scores_array.mean()),
        },
    }

    # Feature importance stability if available
    if feature_importances:
        importance_array = np.array(feature_importances)
        results['feature_importance_stability'] = {
            'mean_importances': importance_array.mean(axis=0).tolist(),
            'std_importances': importance_array.std(axis=0).tolist(),
            'cv_importances': (importance_array.std(axis=0) / importance_array.mean(axis=0)).tolist(),
        }

    logger.info(
        f"Validated model stability: {results['scores']['mean']:.3f} ± "
        f"{results['scores']['std']:.3f} ({n_bootstraps} bootstraps)"
    )

    return results


def compare_validation_strategies(
    X: np.ndarray,
    y: np.ndarray,
    model_factory: callable,
    strategies: List[str] = None,
    random_state: int | None = None,
) -> Dict[str, Any]:
    """Compare different validation strategies.

    Args:
        X: Feature matrix
        y: Target labels/values
        model_factory: Function that returns a new model instance
        strategies: List of validation strategies to compare
        random_state: Random state for reproducibility

    Returns:
        Dictionary with comparison results
    """
    if strategies is None:
        strategies = ['holdout', '5fold_cv', '10fold_cv', 'stratified_5fold']

    results = {}

    for strategy in strategies:
        try:
            if strategy == 'holdout':
                X_train, X_test, y_train, y_test = train_test_split_biological(
                    X, y, test_size=0.2, random_state=random_state
                )
                model = model_factory()
                model.fit(X_train, y_train)
                y_pred = model.predict(X_test)

                if len(np.unique(y)) < 20:  # Classification
                    score = accuracy_score(y_test, y_pred)
                else:  # Regression
                    score = r2_score(y_test, y_pred)

            elif strategy.endswith('fold_cv'):
                cv_folds = int(strategy.split('_')[0])
                scores = cross_validation_scores(
                    model_factory(), X, y, cv=cv_folds, random_state=random_state
                )
                score = scores.get('accuracy', scores.get('r2', np.array([0]))).mean()

            elif strategy == 'stratified_5fold':
                scores = cross_validation_scores(
                    model_factory(), X, y,
                    cv=StratifiedKFold(n_splits=5, shuffle=True, random_state=random_state),
                    random_state=random_state
                )
                score = scores.get('accuracy', np.array([0])).mean()

            else:
                raise ValueError(f"Unknown strategy: {strategy}")

            results[strategy] = {
                'score': float(score),
                'strategy': strategy,
            }

        except Exception as e:
            logger.error(f"Validation strategy {strategy} failed: {e}")
            results[strategy] = {'error': str(e)}

    # Find best strategy
    valid_results = [
        (strategy, result.get('score', 0))
        for strategy, result in results.items()
        if isinstance(result, dict) and 'score' in result
    ]

    if valid_results:
        best_strategy = max(valid_results, key=lambda x: x[1])[0]
    else:
        best_strategy = None

    return {
        'comparison': results,
        'best_strategy': best_strategy,
        'strategies_tested': strategies,
    }


def biological_data_validator(
    X: np.ndarray,
    y: np.ndarray,
    checks: List[str] = None,
) -> Dict[str, Any]:
    """Validate biological data for machine learning.

    Args:
        X: Feature matrix
        y: Target labels/values
        checks: List of validation checks to perform

    Returns:
        Dictionary with validation results
    """
    if checks is None:
        checks = ['missing_values', 'infinite_values', 'constant_features', 'class_balance']

    results = {
        'n_samples': len(X),
        'n_features': X.shape[1],
        'checks': {},
        'passed': True,
    }

    # Check for missing values
    if 'missing_values' in checks:
        x_missing = np.isnan(X).sum()
        y_missing = np.isnan(y).sum() if np.issubdtype(y.dtype, np.floating) else 0

        results['checks']['missing_values'] = {
            'X_missing': int(x_missing),
            'y_missing': int(y_missing),
            'total_missing': int(x_missing + y_missing),
            'passed': x_missing == 0 and y_missing == 0,
        }

        if not results['checks']['missing_values']['passed']:
            results['passed'] = False

    # Check for infinite values
    if 'infinite_values' in checks:
        x_inf = np.isinf(X).sum()
        y_inf = np.isinf(y).sum() if np.issubdtype(y.dtype, np.floating) else 0

        results['checks']['infinite_values'] = {
            'X_infinite': int(x_inf),
            'y_infinite': int(y_inf),
            'total_infinite': int(x_inf + y_inf),
            'passed': x_inf == 0 and y_inf == 0,
        }

        if not results['checks']['infinite_values']['passed']:
            results['passed'] = False

    # Check for constant features
    if 'constant_features' in checks:
        constant_features = []
        for i in range(X.shape[1]):
            if np.unique(X[:, i]).size == 1:
                constant_features.append(i)

        results['checks']['constant_features'] = {
            'constant_feature_indices': constant_features,
            'n_constant_features': len(constant_features),
            'passed': len(constant_features) == 0,
        }

        if not results['checks']['constant_features']['passed']:
            results['passed'] = False

    # Check class balance for classification
    if 'class_balance' in checks and len(np.unique(y)) < 20:
        unique, counts = np.unique(y, return_counts=True)
        balance_ratio = counts.min() / counts.max()

        results['checks']['class_balance'] = {
            'class_counts': counts.tolist(),
            'class_labels': unique.tolist(),
            'balance_ratio': float(balance_ratio),
            'passed': balance_ratio > 0.1,  # At least 10% of majority class
        }

        if not results['checks']['class_balance']['passed']:
            results['passed'] = False

    logger.info(f"Biological data validation: {'PASSED' if results['passed'] else 'FAILED'}")
    return results


def bootstrap_validate(X: np.ndarray, y: np.ndarray, model_func: callable,
                      n_bootstrap: int = 100, test_size: float = 0.2,
                      random_state: int = 42) -> Dict[str, Any]:
    """Perform bootstrap validation of a model.

    Args:
        X: Feature matrix
        y: Target values
        model_func: Function that takes (X_train, y_train, X_test, y_test) and returns predictions
        n_bootstrap: Number of bootstrap iterations
        test_size: Fraction of data for testing
        random_state: Random seed

    Returns:
        Dictionary with bootstrap validation results
    """
    if not HAS_SKLEARN:
        raise ImportError("scikit-learn required for bootstrap validation")

    np.random.seed(random_state)
    n_samples = len(X)

    scores = []
    predictions = []

    for i in range(n_bootstrap):
        # Bootstrap sample
        bootstrap_indices = np.random.choice(n_samples, size=n_samples, replace=True)
        X_boot = X[bootstrap_indices]
        y_boot = y[bootstrap_indices]

        # Split into train/test
        X_train, X_test, y_train, y_test = train_test_split(
            X_boot, y_boot, test_size=test_size, random_state=random_state + i
        )

        # Train and predict
        y_pred = model_func(X_train, y_train, X_test, y_test)
        predictions.append(y_pred)

        # Calculate score (MSE for regression, accuracy for classification)
        if len(np.unique(y)) > 20:  # Regression
            score = -mean_squared_error(y_test, y_pred)  # Negative MSE
        else:  # Classification
            score = accuracy_score(y_test, y_pred)

        scores.append(score)

    scores_array = np.array(scores)

    results = {
        "mean_score": float(scores_array.mean()),
        "std_score": float(scores_array.std()),
        "scores": scores_array.tolist(),
        "n_bootstrap": n_bootstrap,
        "test_size": test_size,
        "mean_mse": float(-scores_array.mean()) if len(np.unique(y)) > 20 else None,
    }

    logger.info(f"Bootstrap validation: {results['mean_score']:.3f} ± {results['std_score']:.3f}")
    return results


def cross_validate(model: Any, X: np.ndarray, y: np.ndarray, cv: int = 5,
                  scoring: str = "accuracy") -> Dict[str, Any]:
    """Perform cross-validation with specified scoring metric.

    Args:
        model: sklearn model or model factory function
        X: Feature matrix
        y: Target values
        cv: Number of cross-validation folds
        scoring: Scoring metric

    Returns:
        Dictionary with cross-validation results
    """
    if not HAS_SKLEARN:
        raise ImportError("scikit-learn required for cross-validation")

    # Handle model factory vs trained model
    if callable(model):
        model_factory = model
    else:
        model_factory = lambda: model

    scores = cross_val_score(model_factory(), X, y, cv=cv, scoring=scoring)

    results = {
        "scores": scores.tolist(),
        "mean_score": float(scores.mean()),
        "std_score": float(scores.std()),
        "cv_folds": cv,
        "scoring": scoring,
    }

    logger.info(f"Cross-validation ({cv} folds, {scoring}): {results['mean_score']:.3f} ± {results['std_score']:.3f}")
    return results


def k_fold_split(X: np.ndarray, y: np.ndarray, k: int = 5,
                stratify: bool = True, random_state: int = 42) -> List[Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]]:
    """Split data into k folds for cross-validation.

    Args:
        X: Feature matrix
        y: Target values
        k: Number of folds
        stratify: Whether to stratify folds
        random_state: Random seed

    Returns:
        List of (X_train, X_test, y_train, y_test) tuples for each fold
    """
    if not HAS_SKLEARN:
        raise ImportError("scikit-learn required for k-fold splitting")

    if stratify and len(np.unique(y)) < 20:
        kf = StratifiedKFold(n_splits=k, shuffle=True, random_state=random_state)
    else:
        kf = KFold(n_splits=k, shuffle=True, random_state=random_state)

    folds = []
    for train_idx, test_idx in kf.split(X, y):
        X_train, X_test = X[train_idx], X[test_idx]
        y_train, y_test = y[train_idx], y[test_idx]
        folds.append((X_train, X_test, y_train, y_test))

    logger.info(f"Created {k}-fold split with {len(folds)} folds")
    return folds


def learning_curve(X: np.ndarray, y: np.ndarray, model_factory: callable,
                  train_sizes: np.ndarray = None, cv: int = 5,
                  random_state: int = 42) -> Dict[str, Any]:
    """Generate learning curves for model evaluation.

    Args:
        X: Feature matrix
        y: Target values
        model_factory: Function that returns a new model instance
        train_sizes: Array of training set sizes (fractions)
        cv: Number of cross-validation folds
        random_state: Random seed

    Returns:
        Dictionary with learning curve results
    """
    if not HAS_SKLEARN:
        raise ImportError("scikit-learn required for learning curves")

    if train_sizes is None:
        train_sizes = np.linspace(0.1, 1.0, 10)

    train_scores = []
    val_scores = []

    for train_size in train_sizes:
        fold_scores_train = []
        fold_scores_val = []

        # Create stratified folds
        if len(np.unique(y)) < 20:
            kf = StratifiedKFold(n_splits=cv, shuffle=True, random_state=random_state)
        else:
            kf = KFold(n_splits=cv, shuffle=True, random_state=random_state)

        for train_idx, val_idx in kf.split(X, y):
            # Subsample training data
            n_train = int(len(train_idx) * train_size)
            train_subset = np.random.choice(train_idx, size=n_train, replace=False)

            X_train_fold = X[train_subset]
            y_train_fold = y[train_subset]
            X_val_fold = X[val_idx]
            y_val_fold = y[val_idx]

            # Train model
            model = model_factory()
            model.fit(X_train_fold, y_train_fold)

            # Score on train and validation
            if len(np.unique(y)) < 20:  # Classification
                train_score = accuracy_score(y_train_fold, model.predict(X_train_fold))
                val_score = accuracy_score(y_val_fold, model.predict(X_val_fold))
            else:  # Regression
                train_score = r2_score(y_train_fold, model.predict(X_train_fold))
                val_score = r2_score(y_val_fold, model.predict(X_val_fold))

            fold_scores_train.append(train_score)
            fold_scores_val.append(val_score)

        train_scores.append(np.mean(fold_scores_train))
        val_scores.append(np.mean(fold_scores_val))

    results = {
        "train_sizes": train_sizes.tolist(),
        "train_scores": train_scores,
        "validation_scores": val_scores,
        "cv_folds": cv,
    }

    logger.info(f"Generated learning curve with {len(train_sizes)} training sizes")
    return results



