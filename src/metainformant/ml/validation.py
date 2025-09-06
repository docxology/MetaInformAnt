"""Model validation and cross-validation methods."""

from __future__ import annotations

import warnings
from typing import Any, Callable, Dict, List, Optional, Tuple

import numpy as np


def train_test_split(
    X: np.ndarray, y: np.ndarray, test_size: float = 0.2, random_state: Optional[int] = None, stratify: bool = False
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Split data into training and testing sets.

    Args:
        X: Feature matrix
        y: Target vector
        test_size: Fraction of data for testing
        random_state: Random seed
        stratify: Whether to stratify split by target

    Returns:
        Tuple of (X_train, X_test, y_train, y_test)
    """
    if random_state is not None:
        np.random.seed(random_state)

    n_samples = X.shape[0]
    n_test = int(n_samples * test_size)

    if stratify:
        # Group samples by class
        unique_classes = np.unique(y)
        train_indices = []
        test_indices = []

        for cls in unique_classes:
            cls_indices = np.where(y == cls)[0]
            n_cls_test = max(1, int(len(cls_indices) * test_size))

            cls_shuffled = np.random.permutation(cls_indices)
            test_indices.extend(cls_shuffled[:n_cls_test])
            train_indices.extend(cls_shuffled[n_cls_test:])

        train_indices = np.array(train_indices)
        test_indices = np.array(test_indices)

    else:
        indices = np.random.permutation(n_samples)
        test_indices = indices[:n_test]
        train_indices = indices[n_test:]

    return X[train_indices], X[test_indices], y[train_indices], y[test_indices]


def k_fold_split(
    X: np.ndarray, y: np.ndarray, k: int = 5, random_state: Optional[int] = None, stratified: bool = False
) -> List[Tuple[np.ndarray, np.ndarray]]:
    """Generate k-fold cross-validation splits.

    Args:
        X: Feature matrix
        y: Target vector
        k: Number of folds
        random_state: Random seed
        stratified: Whether to use stratified k-fold

    Returns:
        List of (train_indices, val_indices) tuples
    """
    if random_state is not None:
        np.random.seed(random_state)

    n_samples = X.shape[0]

    if stratified:
        # Stratified k-fold
        unique_classes = np.unique(y)
        fold_indices = [[] for _ in range(k)]

        for cls in unique_classes:
            cls_indices = np.where(y == cls)[0]
            np.random.shuffle(cls_indices)

            # Distribute class samples across folds
            for i, idx in enumerate(cls_indices):
                fold_indices[i % k].append(idx)

    else:
        # Regular k-fold
        indices = np.random.permutation(n_samples)
        fold_size = n_samples // k
        fold_indices = []

        for i in range(k):
            start = i * fold_size
            end = start + fold_size if i < k - 1 else n_samples
            fold_indices.append(indices[start:end])

    # Generate train/validation splits
    splits = []
    for i in range(k):
        val_indices = np.array(fold_indices[i])
        train_indices = np.concatenate([fold_indices[j] for j in range(k) if j != i])
        splits.append((train_indices, val_indices))

    return splits


def cross_validate(
    model_func: Callable,
    X: np.ndarray,
    y: np.ndarray,
    cv: int = 5,
    scoring: str = "accuracy",
    random_state: Optional[int] = None,
    **model_kwargs,
) -> Dict[str, Any]:
    """Perform cross-validation.

    Args:
        model_func: Function that returns trained model
        X: Feature matrix
        y: Target vector
        cv: Number of cross-validation folds
        scoring: Scoring metric
        random_state: Random seed
        **model_kwargs: Parameters for model_func

    Returns:
        Dictionary with validation scores and statistics
    """
    # Determine if this is classification or regression
    is_classification = len(np.unique(y)) < X.shape[0] * 0.1

    # Generate cross-validation splits
    splits = k_fold_split(X, y, k=cv, random_state=random_state, stratified=is_classification)

    scores = []
    predictions_all = []
    true_labels_all = []

    for fold_idx, (train_idx, val_idx) in enumerate(splits):
        X_train, X_val = X[train_idx], X[val_idx]
        y_train, y_val = y[train_idx], y[val_idx]

        # Train model
        model = model_func(X_train, y_train, **model_kwargs)

        # Make predictions
        if hasattr(model, "predict"):
            y_pred = model.predict(X_val)
        elif callable(model):
            y_pred = model(X_val)
        else:
            raise ValueError("Model must have 'predict' method or be callable")

        # Calculate score
        if scoring == "accuracy" and is_classification:
            score = np.mean(y_pred == y_val)
        elif scoring == "mse":
            score = np.mean((y_pred - y_val) ** 2)
        elif scoring == "rmse":
            score = np.sqrt(np.mean((y_pred - y_val) ** 2))
        elif scoring == "r2":
            ss_res = np.sum((y_val - y_pred) ** 2)
            ss_tot = np.sum((y_val - np.mean(y_val)) ** 2)
            score = 1 - (ss_res / (ss_tot + 1e-8))
        else:
            raise ValueError(f"Unknown scoring method: {scoring}")

        scores.append(score)
        predictions_all.extend(y_pred)
        true_labels_all.extend(y_val)

    predictions_all = np.array(predictions_all)
    true_labels_all = np.array(true_labels_all)

    return {
        "scores": np.array(scores),
        "mean_score": np.mean(scores),
        "std_score": np.std(scores),
        "predictions": predictions_all,
        "true_labels": true_labels_all,
        "scoring_method": scoring,
        "n_folds": cv,
    }


def bootstrap_validate(
    model_func: Callable,
    X: np.ndarray,
    y: np.ndarray,
    n_bootstrap: int = 100,
    sample_size: Optional[int] = None,
    scoring: str = "accuracy",
    random_state: Optional[int] = None,
    **model_kwargs,
) -> Dict[str, Any]:
    """Perform bootstrap validation.

    Args:
        model_func: Function that returns trained model
        X: Feature matrix
        y: Target vector
        n_bootstrap: Number of bootstrap iterations
        sample_size: Size of bootstrap samples
        scoring: Scoring metric
        random_state: Random seed
        **model_kwargs: Parameters for model_func

    Returns:
        Dictionary with bootstrap validation results
    """
    if random_state is not None:
        np.random.seed(random_state)

    n_samples = X.shape[0]
    if sample_size is None:
        sample_size = n_samples

    scores = []

    for i in range(n_bootstrap):
        # Bootstrap sample
        bootstrap_idx = np.random.choice(n_samples, size=sample_size, replace=True)
        out_of_bag_idx = np.setdiff1d(np.arange(n_samples), np.unique(bootstrap_idx))

        if len(out_of_bag_idx) == 0:
            continue

        X_train, y_train = X[bootstrap_idx], y[bootstrap_idx]
        X_test, y_test = X[out_of_bag_idx], y[out_of_bag_idx]

        # Train model
        model = model_func(X_train, y_train, **model_kwargs)

        # Make predictions
        if hasattr(model, "predict"):
            y_pred = model.predict(X_test)
        elif callable(model):
            y_pred = model(X_test)
        else:
            raise ValueError("Model must have 'predict' method or be callable")

        # Calculate score
        if scoring == "accuracy":
            score = np.mean(y_pred == y_test)
        elif scoring == "mse":
            score = np.mean((y_pred - y_test) ** 2)
        elif scoring == "rmse":
            score = np.sqrt(np.mean((y_pred - y_test) ** 2))
        elif scoring == "r2":
            ss_res = np.sum((y_test - y_pred) ** 2)
            ss_tot = np.sum((y_test - np.mean(y_test)) ** 2)
            score = 1 - (ss_res / (ss_tot + 1e-8))
        else:
            raise ValueError(f"Unknown scoring method: {scoring}")

        scores.append(score)

    scores = np.array(scores)

    return {
        "scores": scores,
        "mean_score": np.mean(scores),
        "std_score": np.std(scores),
        "confidence_interval_95": np.percentile(scores, [2.5, 97.5]),
        "n_bootstrap": len(scores),
        "scoring_method": scoring,
    }


def learning_curve(
    model_func: Callable,
    X: np.ndarray,
    y: np.ndarray,
    train_sizes: Optional[np.ndarray] = None,
    cv: int = 5,
    scoring: str = "accuracy",
    random_state: Optional[int] = None,
    **model_kwargs,
) -> Dict[str, Any]:
    """Generate learning curves.

    Args:
        model_func: Function that returns trained model
        X: Feature matrix
        y: Target vector
        train_sizes: Training set sizes to evaluate
        cv: Number of cross-validation folds
        scoring: Scoring metric
        random_state: Random seed
        **model_kwargs: Parameters for model_func

    Returns:
        Dictionary with learning curve data
    """
    if train_sizes is None:
        train_sizes = np.linspace(0.1, 1.0, 10)

    n_samples = X.shape[0]
    train_sizes_abs = (train_sizes * n_samples).astype(int)

    train_scores_all = []
    val_scores_all = []

    for train_size in train_sizes_abs:
        train_scores = []
        val_scores = []

        # Generate CV splits
        splits = k_fold_split(X, y, k=cv, random_state=random_state)

        for train_idx, val_idx in splits:
            # Subsample training data
            if len(train_idx) > train_size:
                train_idx = np.random.choice(train_idx, size=train_size, replace=False)

            X_train, y_train = X[train_idx], y[train_idx]
            X_val, y_val = X[val_idx], y[val_idx]

            # Train model
            model = model_func(X_train, y_train, **model_kwargs)

            # Predictions and scores
            if hasattr(model, "predict"):
                y_train_pred = model.predict(X_train)
                y_val_pred = model.predict(X_val)
            elif callable(model):
                y_train_pred = model(X_train)
                y_val_pred = model(X_val)
            else:
                raise ValueError("Model must have 'predict' method or be callable")

            # Calculate scores
            if scoring == "accuracy":
                train_score = np.mean(y_train_pred == y_train)
                val_score = np.mean(y_val_pred == y_val)
            elif scoring == "mse":
                train_score = np.mean((y_train_pred - y_train) ** 2)
                val_score = np.mean((y_val_pred - y_val) ** 2)
            elif scoring == "r2":
                ss_res_train = np.sum((y_train - y_train_pred) ** 2)
                ss_tot_train = np.sum((y_train - np.mean(y_train)) ** 2)
                train_score = 1 - (ss_res_train / (ss_tot_train + 1e-8))

                ss_res_val = np.sum((y_val - y_val_pred) ** 2)
                ss_tot_val = np.sum((y_val - np.mean(y_val)) ** 2)
                val_score = 1 - (ss_res_val / (ss_tot_val + 1e-8))
            else:
                raise ValueError(f"Unknown scoring method: {scoring}")

            train_scores.append(train_score)
            val_scores.append(val_score)

        train_scores_all.append(train_scores)
        val_scores_all.append(val_scores)

    train_scores_all = np.array(train_scores_all)
    val_scores_all = np.array(val_scores_all)

    return {
        "train_sizes": train_sizes_abs,
        "train_scores": train_scores_all,
        "val_scores": val_scores_all,
        "train_scores_mean": np.mean(train_scores_all, axis=1),
        "train_scores_std": np.std(train_scores_all, axis=1),
        "val_scores_mean": np.mean(val_scores_all, axis=1),
        "val_scores_std": np.std(val_scores_all, axis=1),
        "scoring_method": scoring,
    }
