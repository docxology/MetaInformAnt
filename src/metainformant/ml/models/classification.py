"""Machine learning classification utilities for METAINFORMANT.

This module provides ensemble classifiers and evaluation tools specifically
designed for biological data analysis, with support for cross-validation
and biological data handling.
"""

from __future__ import annotations

from typing import Any, Dict, List, Optional, cast

import numpy as np

from metainformant.core.utils.logging import get_logger

logger = get_logger(__name__)

# Optional imports for ML functionality
try:
    from sklearn.ensemble import GradientBoostingClassifier, RandomForestClassifier, VotingClassifier
    from sklearn.linear_model import LogisticRegression
    from sklearn.metrics import (
        accuracy_score,
        classification_report,
        confusion_matrix,
        f1_score,
        precision_score,
        recall_score,
        roc_auc_score,
    )
    from sklearn.model_selection import StratifiedKFold, cross_val_score

    HAS_SKLEARN = True
except ImportError:
    HAS_SKLEARN = False
    logger.warning("scikit-learn not available, ML classification disabled")


def _calibrated_svm(**kwargs: Any) -> Any:
    """Build an SVM with calibrated probabilities without deprecated SVC flags.

    Recent scikit-learn releases are deprecating the internal probability
    fitting path behind ``SVC(probability=True)``.  Calibration preserves the
    public ``predict_proba`` contract while making the probabilistic behavior
    explicit and version-stable.
    """
    from sklearn.calibration import CalibratedClassifierCV
    from sklearn.svm import SVC

    return CalibratedClassifierCV(estimator=SVC(**kwargs), cv=3)


class BiologicalClassifier:
    """Wrapper for biological data classification with evaluation metrics.

    Can be initialized in two ways:
    1. Algorithm-based: BiologicalClassifier(algorithm="random_forest", random_state=42)
    2. Model-based: BiologicalClassifier(model=sklearn_model, model_type="rf")
    """

    _ALGORITHM_MAP = {
        "random_forest": lambda **kw: RandomForestClassifier(**kw) if HAS_SKLEARN else None,
        "svm": lambda **kw: (_calibrated_svm(**kw) if HAS_SKLEARN else None),
        "logistic_regression": lambda **kw: LogisticRegression(max_iter=1000, **kw) if HAS_SKLEARN else None,
        "gradient_boosting": lambda **kw: GradientBoostingClassifier(**kw) if HAS_SKLEARN else None,
    }

    def __init__(
        self,
        model: Any = None,
        model_type: str = "unknown",
        algorithm: str | None = None,
        random_state: int | None = None,
        **params: Any,
    ):
        """Initialize biological classifier.

        Args:
            model: Pre-built sklearn model (model-based init)
            model_type: Type of model when using model-based init
            algorithm: Algorithm name for algorithm-based init
                       ("random_forest", "svm", "logistic_regression", "gradient_boosting")
            random_state: Random state for reproducibility
            **params: Additional parameters for the sklearn model
        """
        self.random_state = random_state
        self.params = params
        self.feature_names: Optional[List[str]] = None
        self._is_fitted = False

        if algorithm is not None:
            # Algorithm-based initialization
            self.algorithm = algorithm
            self.model_type = algorithm
            if not HAS_SKLEARN:
                raise ImportError("scikit-learn required for classification")
            factory = self._ALGORITHM_MAP.get(algorithm)
            if factory is None:
                raise ValueError(f"Unknown algorithm: {algorithm}")
            model_params = {**params}
            if random_state is not None:
                model_params["random_state"] = random_state
            self.model: Any = factory(**model_params)
        elif model is not None:
            # Model-based initialization
            self.algorithm = model_type
            self.model_type = model_type
            self.model = model
        else:
            # Default: random forest
            self.algorithm = "random_forest"
            self.model_type = "random_forest"
            if HAS_SKLEARN:
                rf_params = {**params}
                if random_state is not None:
                    rf_params["random_state"] = random_state
                self.model = RandomForestClassifier(**rf_params)
            else:
                self.model = None

    @property
    def is_fitted(self) -> bool:
        """Whether the model has been fitted."""
        return self._is_fitted

    @property
    def is_trained(self) -> bool:
        """Alias for is_fitted."""
        return self._is_fitted

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
        self.classes_ = getattr(self.model, "classes_", np.unique(y))
        self._is_fitted = True

        logger.info(f"Trained {self.model_type} classifier with {X.shape[1]} features")
        return self

    def predict(self, X: np.ndarray) -> np.ndarray:
        """Make predictions.

        Args:
            X: Feature matrix

        Returns:
            Predicted labels
        """
        if not self._is_fitted:
            raise ValueError("Model not fitted")
        return cast("np.ndarray", self.model.predict(X))

    def predict_proba(self, X: np.ndarray) -> np.ndarray:
        """Predict class probabilities.

        Args:
            X: Feature matrix

        Returns:
            Class probabilities
        """
        if not self._is_fitted:
            raise ValueError("Model not fitted")
        if hasattr(self.model, "predict_proba"):
            return cast("np.ndarray", self.model.predict_proba(X))
        else:
            raise AttributeError(f"Model {self.model_type} does not support probability prediction")

    def get_feature_importance(self) -> np.ndarray:
        """Get feature importance from the fitted model.

        Returns:
            Array of feature importance values
        """
        if not self._is_fitted:
            raise ValueError("Model not fitted")
        if hasattr(self.model, "feature_importances_"):
            return cast("np.ndarray", self.model.feature_importances_)
        elif hasattr(self.model, "coef_"):
            return cast(
                "np.ndarray",
                np.abs(self.model.coef_).mean(axis=0) if self.model.coef_.ndim > 1 else np.abs(self.model.coef_),
            )
        else:
            raise AttributeError(f"Model {self.model_type} does not support feature importance")

    def evaluate(self, X: np.ndarray, y: np.ndarray, detailed: bool = True) -> Dict[str, Any]:
        """Evaluate classifier performance.

        Args:
            X: Feature matrix
            y: True labels
            detailed: Whether to include detailed metrics

        Returns:
            Dictionary with evaluation metrics
        """
        if not self._is_fitted:
            raise ValueError("Model not fitted")

        y_pred = self.predict(X)

        results = {
            "accuracy": accuracy_score(y, y_pred),
            "precision": precision_score(y, y_pred, average="weighted", zero_division=0),
            "recall": recall_score(y, y_pred, average="weighted", zero_division=0),
            "f1": f1_score(y, y_pred, average="weighted", zero_division=0),
        }

        # Add probability-based metrics if available
        if hasattr(self.model, "predict_proba"):
            try:
                y_proba = self.predict_proba(X)
                if len(np.unique(y)) == 2:  # Binary classification
                    results["roc_auc"] = roc_auc_score(y, y_proba[:, 1])
                else:  # Multi-class
                    results["roc_auc"] = roc_auc_score(y, y_proba, multi_class="ovr", average="weighted")
            except Exception as e:
                logger.warning(f"Could not compute ROC-AUC: {e}")

        if detailed:
            results["confusion_matrix"] = confusion_matrix(y, y_pred).tolist()
            results["classification_report"] = classification_report(y, y_pred, output_dict=True, zero_division=0)

        return results


def train_ensemble_classifier(
    X_train: np.ndarray,
    y_train: np.ndarray,
    n_estimators: int = 10,
    random_state: int | None = None,
    n_jobs: int | None = -1,
    **kwargs: Any,
) -> BiologicalClassifier:
    """Train an ensemble classifier for biological data.

    Args:
        X_train: Training feature matrix
        y_train: Training target labels
        n_estimators: Number of estimators in ensemble
        random_state: Random state for reproducibility
        n_jobs: Number of parallel jobs used by the voting ensemble. Use 1
            for deterministic single-process examples and restricted runners.
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
        (
            "rf",
            RandomForestClassifier(n_estimators=n_estimators, random_state=random_state, **kwargs.get("rf_params", {})),
        ),
        (
            "gb",
            GradientBoostingClassifier(
                n_estimators=n_estimators, random_state=random_state, **kwargs.get("gb_params", {})
            ),
        ),
        ("lr", LogisticRegression(random_state=random_state, max_iter=1000, **kwargs.get("lr_params", {}))),
    ]

    # Voting classifier
    ensemble = VotingClassifier(
        estimators=estimators,
        voting="soft",  # Use probability voting
        n_jobs=n_jobs,
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

    results = classifier.evaluate(X_eval, y_eval, detailed=True)
    results["predictions"] = classifier.predict(X_eval).tolist()
    return results


def _cross_validated_metric(
    estimator: Any, X: np.ndarray, y: np.ndarray, cv: Any, metric: str, n_jobs: int | None
) -> Dict[str, Any]:
    """Cross-validate one scoring metric and summarize its fold scores.

    Failures propagate. ``error_score="raise"`` stops scikit-learn from replacing a
    failed fold with NaN, and a metric that is undefined without raising (non-finite
    fold scores) is rejected too, so a returned summary always holds real numbers.
    """
    try:
        scores = cross_val_score(estimator, X, y, cv=cv, scoring=metric, n_jobs=n_jobs, error_score="raise")
    except Exception as exc:
        exc.add_note(f"cross_validate_biological: scoring metric {metric!r} failed")
        raise
    if not np.all(np.isfinite(scores)):
        raise ValueError(
            f"cross_validate_biological: metric {metric!r} is not finite for every fold: {scores.tolist()}"
        )
    return {"mean": float(scores.mean()), "std": float(scores.std()), "scores": scores.tolist()}


def cross_validate_biological(
    X: np.ndarray,
    y: np.ndarray,
    method: str = "rf",
    cv_folds: int = 5,
    random_state: int | None = None,
    algorithm: str | None = None,
    parallel_jobs: int | None = -1,
    **kwargs: Any,
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
        parallel_jobs: Number of parallel jobs used for cross-validation and, for ``method="ensemble"``, by the
            voting ensemble itself (an explicit ``n_jobs`` keyword overrides the latter). Use 1 in restricted or
            resource-sensitive environments: it then starts no worker processes at all.
        **kwargs: Additional parameters for classifier

    Returns:
        Dictionary with cross-validation results. Every entry of
        ``cross_validation`` holds the finite per-fold scores of one metric
        (``accuracy``, ``precision_weighted``, ``recall_weighted``,
        ``f1_weighted``, plus ``roc_auc`` for two classes).

    Raises:
        ImportError: If scikit-learn not available
        ValueError: If invalid method specified, or if a metric cannot be
            computed (for example too few samples per class for ``cv_folds``,
            non-finite features, or ROC-AUC undefined in a fold). Failures are
            never converted into null or NaN metrics; the original exception
            propagates with the failing metric named in a note.
    """
    if not HAS_SKLEARN:
        raise ImportError("scikit-learn required for cross-validation")

    # Allow algorithm as alias for method (map full names to short codes)
    if algorithm is not None:
        _algo_to_method = {"random_forest": "rf", "gradient_boosting": "gb", "logistic_regression": "lr"}
        method = _algo_to_method.get(algorithm, algorithm)

    # Select classifier. Build each estimator's options once so callers can
    # override defaults without passing duplicate keyword arguments to
    # scikit-learn (which otherwise raises at runtime).
    classifier_kwargs = dict(kwargs)
    if random_state is not None:
        classifier_kwargs["random_state"] = random_state

    if method == "rf":
        classifier_kwargs.setdefault("n_estimators", 100)
        classifier = RandomForestClassifier(**classifier_kwargs)
    elif method == "gb":
        classifier_kwargs.setdefault("n_estimators", 100)
        classifier = GradientBoostingClassifier(**classifier_kwargs)
    elif method == "lr":
        classifier_kwargs.setdefault("max_iter", 1000)
        classifier = LogisticRegression(**classifier_kwargs)
    elif method == "ensemble":
        # Use ensemble classifier; train_ensemble_classifier takes random_state as its own argument, so it must not
        # also arrive through **classifier_kwargs (that is a "multiple values for keyword argument" TypeError)
        ensemble_jobs = classifier_kwargs.pop("n_jobs", parallel_jobs)
        classifier_kwargs.pop("random_state", None)
        ensemble_classifier = train_ensemble_classifier(
            X,
            y,
            random_state=random_state,
            n_jobs=ensemble_jobs,
            **classifier_kwargs,
        )
        classifier = ensemble_classifier.model
    else:
        raise ValueError(f"Unknown classification method: {method}")

    # Perform stratified cross-validation
    cv = StratifiedKFold(n_splits=cv_folds, shuffle=True, random_state=random_state)

    # Evaluate multiple metrics; a metric that cannot be computed raises (see _cross_validated_metric)
    scoring = ["accuracy", "precision_weighted", "recall_weighted", "f1_weighted"]
    if len(np.unique(y)) == 2:
        scoring.append("roc_auc")  # binary classification only

    cv_results: Dict[str, Any] = {
        metric: _cross_validated_metric(classifier, X, y, cv, metric, parallel_jobs) for metric in scoring
    }

    results: Dict[str, Any] = {
        "method": method,
        "cv_folds": cv_folds,
        "n_samples": len(X),
        "n_features": X.shape[1],
        "n_classes": len(np.unique(y)),
        "class_distribution": np.unique(y, return_counts=True)[1].tolist(),
        "cross_validation": cv_results,
    }

    # Add top-level accuracy keys for convenience
    accuracy_mean = cv_results["accuracy"]["mean"]
    results["accuracy"] = accuracy_mean
    results["mean_accuracy"] = accuracy_mean

    logger.info(f"Completed {cv_folds}-fold CV for {method} classifier: accuracy={accuracy_mean:.3f}")

    return results


def create_biological_classifier(method: str = "rf", **kwargs: Any) -> BiologicalClassifier:
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
    methods: List[str] | None = None,
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
        Dictionary with ``comparison`` (per-method cross-validation results, or
        ``{"error": message}`` for a method whose evaluation failed) and
        ``best_method`` (the successful method with the highest mean accuracy,
        or None when every method failed)
    """
    if methods is None:
        methods = ["rf", "gb", "lr"]

    results = {}

    for method in methods:
        try:
            cv_result = cross_validate_biological(X, y, method=method, cv_folds=cv_folds, random_state=random_state)
            results[method] = cv_result
        except Exception as e:
            logger.error(f"Failed to evaluate {method}: {e}")
            results[method] = {"error": str(e)}

    # Only methods that completed cross-validation compete; a failed method is reported under
    # ``comparison[method]["error"]`` and never wins. ``best_method`` is None when every method failed.
    scored = [(m, r["cross_validation"]["accuracy"]["mean"]) for m, r in results.items() if "error" not in r]
    best_method = max(scored, key=lambda item: item[1])[0] if scored else None

    return {
        "comparison": results,
        "best_method": best_method,
    }
