"""Machine learning classification utilities for METAINFORMANT.

This module provides ensemble classifiers and evaluation tools specifically
designed for biological data analysis, with support for cross-validation
and biological data handling.
"""

from __future__ import annotations

import numpy as np
from typing import Any, Dict, List, Optional, Tuple, Union

from metainformant.core.logging import get_logger

logger = get_logger(__name__)

# Optional imports for ML functionality
try:
    from sklearn.base import BaseEstimator
    from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier, VotingClassifier
    from sklearn.linear_model import LogisticRegression
    from sklearn.model_selection import cross_val_score, StratifiedKFold
    from sklearn.metrics import (
        accuracy_score, precision_score, recall_score, f1_score,
        roc_auc_score, confusion_matrix, classification_report
    )
    HAS_SKLEARN = True
except ImportError:
    HAS_SKLEARN = False
    logger.warning("scikit-learn not available, ML classification disabled")


class BiologicalClassifier:
    """Wrapper for biological data classification with evaluation metrics."""

    def __init__(self, model: Any, model_type: str = "unknown"):
        """Initialize biological classifier.

        Args:
            model: Trained sklearn model
            model_type: Type of model (rf, gb, lr, etc.)
        """
        self.model = model
        self.model_type = model_type
        self.feature_names: Optional[List[str]] = None
        self.is_trained = False

    def fit(self, X: np.ndarray, y: np.ndarray, feature_names: Optional[List[str]] = None) -> BiologicalClassifier:
        """Fit the classifier.

        Args:
            X: Feature matrix
            y: Target labels
            feature_names: Optional feature names

        Returns:
            Self for chaining
        """
        if not HAS_SKLEARN:
            raise ImportError("scikit-learn required for classification")

        self.model.fit(X, y)
        self.feature_names = feature_names
        self.is_trained = True

        logger.info(f"Trained {self.model_type} classifier with {X.shape[1]} features")
        return self

    def predict(self, X: np.ndarray) -> np.ndarray:
        """Make predictions.

        Args:
            X: Feature matrix

        Returns:
            Predicted labels
        """
        if not self.is_trained:
            raise ValueError("Model not trained")
        return self.model.predict(X)

    def predict_proba(self, X: np.ndarray) -> np.ndarray:
        """Predict class probabilities.

        Args:
            X: Feature matrix

        Returns:
            Class probabilities
        """
        if not self.is_trained:
            raise ValueError("Model not trained")
        if hasattr(self.model, 'predict_proba'):
            return self.model.predict_proba(X)
        else:
            raise AttributeError(f"Model {self.model_type} does not support probability prediction")

    def evaluate(
        self,
        X: np.ndarray,
        y: np.ndarray,
        detailed: bool = True
    ) -> Dict[str, Any]:
        """Evaluate classifier performance.

        Args:
            X: Feature matrix
            y: True labels
            detailed: Whether to include detailed metrics

        Returns:
            Dictionary with evaluation metrics
        """
        if not self.is_trained:
            raise ValueError("Model not trained")

        y_pred = self.predict(X)

        results = {
            'accuracy': accuracy_score(y, y_pred),
            'precision': precision_score(y, y_pred, average='weighted', zero_division=0),
            'recall': recall_score(y, y_pred, average='weighted', zero_division=0),
            'f1': f1_score(y, y_pred, average='weighted', zero_division=0),
        }

        # Add probability-based metrics if available
        if hasattr(self.model, 'predict_proba'):
            try:
                y_proba = self.predict_proba(X)
                if len(np.unique(y)) == 2:  # Binary classification
                    results['roc_auc'] = roc_auc_score(y, y_proba[:, 1])
                else:  # Multi-class
                    results['roc_auc'] = roc_auc_score(y, y_proba, multi_class='ovr', average='weighted')
            except Exception as e:
                logger.warning(f"Could not compute ROC-AUC: {e}")

        if detailed:
            results['confusion_matrix'] = confusion_matrix(y, y_pred).tolist()
            results['classification_report'] = classification_report(y, y_pred, output_dict=True, zero_division=0)

        return results


def train_ensemble_classifier(
    X_train: np.ndarray,
    y_train: np.ndarray,
    n_estimators: int = 10,
    random_state: int | None = None,
    **kwargs: Any
) -> BiologicalClassifier:
    """Train an ensemble classifier for biological data.

    Args:
        X_train: Training feature matrix
        y_train: Training target labels
        n_estimators: Number of estimators in ensemble
        random_state: Random state for reproducibility
        **kwargs: Additional arguments for ensemble

    Returns:
        Trained BiologicalClassifier instance

    Raises:
        ImportError: If scikit-learn not available
    """
    if not HAS_SKLEARN:
        raise ImportError("scikit-learn required for ensemble classification")

    # Create ensemble with multiple algorithms
    estimators = [
        ('rf', RandomForestClassifier(
            n_estimators=n_estimators,
            random_state=random_state,
            **kwargs.get('rf_params', {})
        )),
        ('gb', GradientBoostingClassifier(
            n_estimators=n_estimators,
            random_state=random_state,
            **kwargs.get('gb_params', {})
        )),
        ('lr', LogisticRegression(
            random_state=random_state,
            max_iter=1000,
            **kwargs.get('lr_params', {})
        )),
    ]

    # Voting classifier
    ensemble = VotingClassifier(
        estimators=estimators,
        voting='soft',  # Use probability voting
        n_jobs=-1,  # Use all available cores
    )

    # Wrap in BiologicalClassifier
    classifier = BiologicalClassifier(ensemble, "ensemble")
    classifier.fit(X_train, y_train)

    logger.info(f"Trained ensemble classifier with {n_estimators} estimators each")
    return classifier


def evaluate_classifier(
    classifier: BiologicalClassifier,
    X_test: np.ndarray | None = None,
    y_test: np.ndarray | None = None,
    X: np.ndarray | None = None,
    y: np.ndarray | None = None,
) -> Dict[str, Any]:
    """Evaluate a trained classifier.

    Args:
        classifier: Trained BiologicalClassifier
        X_test: Test feature matrix (alternative to X)
        y_test: Test target labels (alternative to y)
        X: Feature matrix (alternative to X_test)
        y: Target labels (alternative to y_test)

    Returns:
        Dictionary with evaluation results

    Raises:
        ValueError: If insufficient test data provided
    """
    # Handle alternative parameter names
    if X_test is not None and y_test is not None:
        X_eval, y_eval = X_test, y_test
    elif X is not None and y is not None:
        X_eval, y_eval = X, y
    else:
        raise ValueError("Must provide either (X_test, y_test) or (X, y)")

    return classifier.evaluate(X_eval, y_eval, detailed=True)


def cross_validate_biological(
    X: np.ndarray,
    y: np.ndarray,
    method: str = "rf",
    cv_folds: int = 5,
    random_state: int | None = None,
    **kwargs: Any
) -> Dict[str, Any]:
    """Perform cross-validation for biological data classification.

    This function uses stratified k-fold cross-validation to ensure
    balanced class distribution in each fold, which is important for
    biological classification tasks.

    Args:
        X: Feature matrix
        y: Target labels
        method: Classification method ('rf', 'gb', 'lr', 'ensemble')
        cv_folds: Number of cross-validation folds
        random_state: Random state for reproducibility
        **kwargs: Additional parameters for classifier

    Returns:
        Dictionary with cross-validation results

    Raises:
        ImportError: If scikit-learn not available
        ValueError: If invalid method specified
    """
    if not HAS_SKLEARN:
        raise ImportError("scikit-learn required for cross-validation")

    # Select classifier
    if method == "rf":
        classifier = RandomForestClassifier(
            random_state=random_state,
            n_estimators=kwargs.get('n_estimators', 100),
            **kwargs
        )
    elif method == "gb":
        classifier = GradientBoostingClassifier(
            random_state=random_state,
            n_estimators=kwargs.get('n_estimators', 100),
            **kwargs
        )
    elif method == "lr":
        classifier = LogisticRegression(
            random_state=random_state,
            max_iter=1000,
            **kwargs
        )
    elif method == "ensemble":
        # Use ensemble classifier
        ensemble_classifier = train_ensemble_classifier(
            X, y, random_state=random_state, **kwargs
        )
        classifier = ensemble_classifier.model
    else:
        raise ValueError(f"Unknown classification method: {method}")

    # Perform stratified cross-validation
    cv = StratifiedKFold(n_splits=cv_folds, shuffle=True, random_state=random_state)

    # Evaluate multiple metrics
    scoring = ['accuracy', 'precision_weighted', 'recall_weighted', 'f1_weighted']

    cv_results = {}
    for metric in scoring:
        try:
            scores = cross_val_score(
                classifier, X, y,
                cv=cv,
                scoring=metric,
                n_jobs=-1
            )
            cv_results[metric] = {
                'mean': float(scores.mean()),
                'std': float(scores.std()),
                'scores': scores.tolist(),
            }
        except Exception as e:
            logger.warning(f"Could not compute {metric}: {e}")
            cv_results[metric] = None

    # Try ROC-AUC if binary classification
    if len(np.unique(y)) == 2:
        try:
            auc_scores = cross_val_score(
                classifier, X, y,
                cv=cv,
                scoring='roc_auc',
                n_jobs=-1
            )
            cv_results['roc_auc'] = {
                'mean': float(auc_scores.mean()),
                'std': float(auc_scores.std()),
                'scores': auc_scores.tolist(),
            }
        except Exception as e:
            logger.warning(f"Could not compute ROC-AUC: {e}")

    results = {
        'method': method,
        'cv_folds': cv_folds,
        'n_samples': len(X),
        'n_features': X.shape[1],
        'n_classes': len(np.unique(y)),
        'class_distribution': np.bincount(y).tolist(),
        'cross_validation': cv_results,
    }

    logger.info(
        f"Completed {cv_folds}-fold CV for {method} classifier: "
        f"accuracy={cv_results.get('accuracy', {}).get('mean', 'N/A'):.3f}"
    )

    return results


def create_biological_classifier(
    method: str = "rf",
    **kwargs: Any
) -> BiologicalClassifier:
    """Create a biological classifier with specified method.

    Args:
        method: Classification method ('rf', 'gb', 'lr')
        **kwargs: Parameters for the classifier

    Returns:
        BiologicalClassifier instance

    Raises:
        ImportError: If scikit-learn not available
        ValueError: If invalid method specified
    """
    if not HAS_SKLEARN:
        raise ImportError("scikit-learn required for classification")

    if method == "rf":
        model = RandomForestClassifier(**kwargs)
    elif method == "gb":
        model = GradientBoostingClassifier(**kwargs)
    elif method == "lr":
        model = LogisticRegression(**kwargs)
    else:
        raise ValueError(f"Unknown classification method: {method}")

    return BiologicalClassifier(model, method)


def compare_classifiers(
    X: np.ndarray,
    y: np.ndarray,
    methods: List[str] = None,
    cv_folds: int = 5,
    random_state: int | None = None,
) -> Dict[str, Any]:
    """Compare multiple classification methods.

    Args:
        X: Feature matrix
        y: Target labels
        methods: List of methods to compare (default: ['rf', 'gb', 'lr'])
        cv_folds: Number of cross-validation folds
        random_state: Random state for reproducibility

    Returns:
        Dictionary with comparison results
    """
    if methods is None:
        methods = ['rf', 'gb', 'lr']

    results = {}

    for method in methods:
        try:
            cv_result = cross_validate_biological(
                X, y, method=method,
                cv_folds=cv_folds,
                random_state=random_state
            )
            results[method] = cv_result
        except Exception as e:
            logger.error(f"Failed to evaluate {method}: {e}")
            results[method] = {'error': str(e)}

    return {
        'comparison': results,
        'best_method': max(
            [(m, r.get('cross_validation', {}).get('accuracy', {}).get('mean', 0))
             for m, r in results.items() if isinstance(r, dict)],
            key=lambda x: x[1],
            default=(None, 0)
        )[0],
    }






