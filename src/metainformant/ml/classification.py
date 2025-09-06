"""Classification methods for biological data analysis."""

from __future__ import annotations

import warnings
from collections import defaultdict
from typing import Any, Dict, List, Optional, Tuple, Union

import numpy as np


class BiologicalClassifier:
    """Classification wrapper for biological data analysis."""

    def __init__(self, algorithm: str = "random_forest", random_state: Optional[int] = None, **kwargs):
        """Initialize biological classifier.

        Args:
            algorithm: Classification algorithm to use
            random_state: Random seed for reproducibility
            **kwargs: Algorithm-specific parameters
        """
        self.algorithm = algorithm
        self.random_state = random_state
        self.params = kwargs
        self.is_fitted = False
        self.classes_ = None
        self.feature_importance_ = None

    def fit(self, X: np.ndarray, y: np.ndarray) -> "BiologicalClassifier":
        """Fit classifier to training data.

        Args:
            X: Training feature matrix (samples x features)
            y: Training labels

        Returns:
            Fitted classifier
        """
        if self.random_state is not None:
            np.random.seed(self.random_state)

        self.classes_ = np.unique(y)
        n_classes = len(self.classes_)

        if n_classes < 2:
            raise ValueError("Need at least 2 classes for classification")

        # Store training data for simple algorithms
        self.X_train_ = X.copy()
        self.y_train_ = y.copy()

        # Calculate feature importance (simplified)
        if self.algorithm == "random_forest":
            self.feature_importance_ = self._rf_feature_importance(X, y)
        elif self.algorithm == "linear":
            self.feature_importance_ = self._linear_feature_importance(X, y)
        else:
            self.feature_importance_ = np.ones(X.shape[1]) / X.shape[1]

        self.is_fitted = True
        return self

    def predict(self, X: np.ndarray) -> np.ndarray:
        """Predict class labels for samples.

        Args:
            X: Test feature matrix

        Returns:
            Predicted class labels
        """
        if not self.is_fitted:
            raise ValueError("Classifier must be fitted before prediction")

        if self.algorithm == "knn":
            return self._knn_predict(X)
        elif self.algorithm == "naive_bayes":
            return self._nb_predict(X)
        elif self.algorithm == "random_forest":
            return self._rf_predict(X)
        elif self.algorithm == "linear":
            return self._linear_predict(X)
        else:
            # Default: nearest centroid
            return self._centroid_predict(X)

    def predict_proba(self, X: np.ndarray) -> np.ndarray:
        """Predict class probabilities for samples.

        Args:
            X: Test feature matrix

        Returns:
            Class probability matrix (samples x classes)
        """
        if not self.is_fitted:
            raise ValueError("Classifier must be fitted before prediction")

        n_samples = X.shape[0]
        n_classes = len(self.classes_)
        probas = np.zeros((n_samples, n_classes))

        # Simple probability estimation based on distances
        for i in range(n_samples):
            # Distance-based probabilities
            distances = []
            for class_idx, class_label in enumerate(self.classes_):
                class_mask = self.y_train_ == class_label
                class_samples = self.X_train_[class_mask]

                if len(class_samples) > 0:
                    # Average distance to class samples
                    dists = [np.linalg.norm(X[i] - sample) for sample in class_samples]
                    avg_dist = np.mean(dists)
                    distances.append(avg_dist)
                else:
                    distances.append(float("inf"))

            # Convert distances to probabilities (inverse relationship)
            distances = np.array(distances)
            if np.any(np.isfinite(distances)):
                # Avoid division by zero
                distances = distances + 1e-10
                inv_distances = 1.0 / distances
                probas[i] = inv_distances / np.sum(inv_distances)
            else:
                # Uniform probabilities if all distances are infinite
                probas[i] = np.ones(n_classes) / n_classes

        return probas

    def _rf_feature_importance(self, X: np.ndarray, y: np.ndarray) -> np.ndarray:
        """Calculate random forest-like feature importance."""
        n_features = X.shape[1]
        importance = np.zeros(n_features)

        # Simple permutation-based importance
        baseline_accuracy = self._calculate_accuracy(X, y, X, y)

        for feature_idx in range(n_features):
            # Permute this feature
            X_permuted = X.copy()
            np.random.shuffle(X_permuted[:, feature_idx])

            # Calculate accuracy with permuted feature
            permuted_accuracy = self._calculate_accuracy(X, y, X_permuted, y)

            # Importance is the drop in accuracy
            importance[feature_idx] = max(0, baseline_accuracy - permuted_accuracy)

        # Normalize
        importance = importance / (np.sum(importance) + 1e-10)
        return importance

    def _linear_feature_importance(self, X: np.ndarray, y: np.ndarray) -> np.ndarray:
        """Calculate linear classifier feature importance."""
        # Simplified: use correlation with target
        n_features = X.shape[1]
        importance = np.zeros(n_features)

        for feature_idx in range(n_features):
            corr = np.corrcoef(X[:, feature_idx], y)[0, 1]
            importance[feature_idx] = abs(corr) if not np.isnan(corr) else 0.0

        # Normalize
        importance = importance / (np.sum(importance) + 1e-10)
        return importance

    def _calculate_accuracy(
        self, X_train: np.ndarray, y_train: np.ndarray, X_test: np.ndarray, y_test: np.ndarray
    ) -> float:
        """Calculate classification accuracy."""
        # Simple centroid-based classification for accuracy calculation
        predictions = self._centroid_predict_with_data(X_test, X_train, y_train)
        accuracy = np.mean(predictions == y_test)
        return accuracy

    def _knn_predict(self, X: np.ndarray, k: int = 3) -> np.ndarray:
        """K-nearest neighbors prediction."""
        k = min(k, len(self.X_train_))
        predictions = []

        for sample in X:
            # Calculate distances to all training samples
            distances = [np.linalg.norm(sample - train_sample) for train_sample in self.X_train_]

            # Get k nearest neighbors
            nearest_indices = np.argsort(distances)[:k]
            nearest_labels = self.y_train_[nearest_indices]

            # Majority vote
            unique_labels, counts = np.unique(nearest_labels, return_counts=True)
            prediction = unique_labels[np.argmax(counts)]
            predictions.append(prediction)

        return np.array(predictions)

    def _nb_predict(self, X: np.ndarray) -> np.ndarray:
        """Naive Bayes prediction (simplified)."""
        predictions = []

        # Calculate class priors and feature statistics
        class_priors = {}
        class_means = {}
        class_stds = {}

        for class_label in self.classes_:
            class_mask = self.y_train_ == class_label
            class_samples = self.X_train_[class_mask]

            class_priors[class_label] = np.sum(class_mask) / len(self.y_train_)
            class_means[class_label] = np.mean(class_samples, axis=0)
            class_stds[class_label] = np.std(class_samples, axis=0) + 1e-8

        # Predict each sample
        for sample in X:
            log_probabilities = {}

            for class_label in self.classes_:
                log_prob = np.log(class_priors[class_label])

                # Gaussian likelihood (assuming independence)
                means = class_means[class_label]
                stds = class_stds[class_label]

                log_likelihood = np.sum(-0.5 * np.log(2 * np.pi * stds**2) - 0.5 * ((sample - means) / stds) ** 2)

                log_probabilities[class_label] = log_prob + log_likelihood

            # Predict class with highest log probability
            prediction = max(log_probabilities.keys(), key=lambda k: log_probabilities[k])
            predictions.append(prediction)

        return np.array(predictions)

    def _rf_predict(self, X: np.ndarray) -> np.ndarray:
        """Random forest prediction (simplified ensemble)."""
        # Create multiple "trees" with different feature subsets
        n_trees = self.params.get("n_trees", 10)
        n_features_per_tree = max(1, int(np.sqrt(X.shape[1])))

        tree_predictions = []

        for _ in range(n_trees):
            # Random feature subset
            feature_indices = np.random.choice(X.shape[1], size=n_features_per_tree, replace=False)
            X_subset = X[:, feature_indices]
            X_train_subset = self.X_train_[:, feature_indices]

            # Simple tree: centroid-based prediction on feature subset
            tree_pred = self._centroid_predict_with_data(X_subset, X_train_subset, self.y_train_)
            tree_predictions.append(tree_pred)

        # Ensemble voting
        tree_predictions = np.array(tree_predictions)
        final_predictions = []

        for sample_idx in range(X.shape[0]):
            sample_votes = tree_predictions[:, sample_idx]
            unique_votes, counts = np.unique(sample_votes, return_counts=True)
            final_prediction = unique_votes[np.argmax(counts)]
            final_predictions.append(final_prediction)

        return np.array(final_predictions)

    def _linear_predict(self, X: np.ndarray) -> np.ndarray:
        """Linear classifier prediction."""
        # For multi-class, use one-vs-rest approach
        if len(self.classes_) == 2:
            # Binary classification
            return self._binary_linear_predict(X)
        else:
            # Multi-class: one-vs-rest
            predictions = []

            for sample in X:
                class_scores = {}

                for class_label in self.classes_:
                    # Create binary problem: this class vs all others
                    binary_labels = (self.y_train_ == class_label).astype(int)

                    # Simple linear score (correlation-based)
                    score = 0
                    for feature_idx in range(len(sample)):
                        feature_corr = np.corrcoef(self.X_train_[:, feature_idx], binary_labels)[0, 1]
                        if not np.isnan(feature_corr):
                            score += feature_corr * sample[feature_idx]

                    class_scores[class_label] = score

                # Predict class with highest score
                prediction = max(class_scores.keys(), key=lambda k: class_scores[k])
                predictions.append(prediction)

            return np.array(predictions)

    def _binary_linear_predict(self, X: np.ndarray) -> np.ndarray:
        """Binary linear classifier."""
        # Simple linear boundary based on feature correlations
        predictions = []

        # Calculate decision boundary
        class_0_samples = self.X_train_[self.y_train_ == self.classes_[0]]
        class_1_samples = self.X_train_[self.y_train_ == self.classes_[1]]

        centroid_0 = np.mean(class_0_samples, axis=0)
        centroid_1 = np.mean(class_1_samples, axis=0)

        for sample in X:
            # Distance to each class centroid
            dist_0 = np.linalg.norm(sample - centroid_0)
            dist_1 = np.linalg.norm(sample - centroid_1)

            prediction = self.classes_[0] if dist_0 < dist_1 else self.classes_[1]
            predictions.append(prediction)

        return np.array(predictions)

    def _centroid_predict(self, X: np.ndarray) -> np.ndarray:
        """Centroid-based prediction."""
        return self._centroid_predict_with_data(X, self.X_train_, self.y_train_)

    def _centroid_predict_with_data(self, X: np.ndarray, X_train: np.ndarray, y_train: np.ndarray) -> np.ndarray:
        """Centroid-based prediction with specific training data."""
        predictions = []

        # Calculate class centroids
        centroids = {}
        for class_label in self.classes_:
            class_mask = y_train == class_label
            class_samples = X_train[class_mask]
            if len(class_samples) > 0:
                centroids[class_label] = np.mean(class_samples, axis=0)

        # Predict based on closest centroid
        for sample in X:
            distances = {}
            for class_label, centroid in centroids.items():
                distances[class_label] = np.linalg.norm(sample - centroid)

            if distances:
                prediction = min(distances.keys(), key=lambda k: distances[k])
            else:
                prediction = self.classes_[0]  # Default to first class
            predictions.append(prediction)

        return np.array(predictions)


def train_ensemble_classifier(
    X: np.ndarray, y: np.ndarray, algorithms: List[str] = None, random_state: Optional[int] = None
) -> Dict[str, BiologicalClassifier]:
    """Train ensemble of classifiers.

    Args:
        X: Training feature matrix
        y: Training labels
        algorithms: List of algorithms to include in ensemble
        random_state: Random seed

    Returns:
        Dictionary of algorithm_name -> fitted_classifier
    """
    if algorithms is None:
        algorithms = ["knn", "naive_bayes", "random_forest", "linear"]

    ensemble = {}

    for algorithm in algorithms:
        try:
            classifier = BiologicalClassifier(algorithm=algorithm, random_state=random_state)
            classifier.fit(X, y)
            ensemble[algorithm] = classifier
        except Exception as e:
            warnings.warn(f"Failed to train {algorithm}: {e}")

    return ensemble


def evaluate_classifier(classifier: BiologicalClassifier, X_test: np.ndarray, y_test: np.ndarray) -> Dict[str, float]:
    """Evaluate classifier performance.

    Args:
        classifier: Fitted classifier
        X_test: Test features
        y_test: True test labels

    Returns:
        Dictionary of evaluation metrics
    """
    predictions = classifier.predict(X_test)
    probabilities = classifier.predict_proba(X_test)

    # Basic metrics
    accuracy = np.mean(predictions == y_test)

    # Per-class metrics
    classes = classifier.classes_
    class_metrics = {}

    for i, class_label in enumerate(classes):
        class_mask = y_test == class_label

        if np.sum(class_mask) > 0:
            # Precision, recall, F1
            true_positives = np.sum((predictions == class_label) & class_mask)
            false_positives = np.sum((predictions == class_label) & ~class_mask)
            false_negatives = np.sum((predictions != class_label) & class_mask)

            precision = true_positives / (true_positives + false_positives + 1e-10)
            recall = true_positives / (true_positives + false_negatives + 1e-10)
            f1 = 2 * precision * recall / (precision + recall + 1e-10)

            class_metrics[f"precision_{class_label}"] = precision
            class_metrics[f"recall_{class_label}"] = recall
            class_metrics[f"f1_{class_label}"] = f1

    # Overall metrics
    precisions = [class_metrics.get(f"precision_{c}", 0) for c in classes]
    recalls = [class_metrics.get(f"recall_{c}", 0) for c in classes]
    f1_scores = [class_metrics.get(f"f1_{c}", 0) for c in classes]

    metrics = {
        "accuracy": accuracy,
        "macro_precision": np.mean(precisions),
        "macro_recall": np.mean(recalls),
        "macro_f1": np.mean(f1_scores),
        **class_metrics,
    }

    return metrics


def cross_validate_biological(
    X: np.ndarray,
    y: np.ndarray,
    algorithm: str = "random_forest",
    cv_folds: int = 5,
    random_state: Optional[int] = None,
) -> Dict[str, List[float]]:
    """Cross-validation for biological classifier.

    Args:
        X: Feature matrix
        y: Labels
        algorithm: Classification algorithm
        cv_folds: Number of CV folds
        random_state: Random seed

    Returns:
        Dictionary of metric_name -> list of fold scores
    """
    if random_state is not None:
        np.random.seed(random_state)

    n_samples = X.shape[0]
    fold_size = n_samples // cv_folds

    # Create fold indices
    indices = np.random.permutation(n_samples)
    folds = []

    for i in range(cv_folds):
        start_idx = i * fold_size
        end_idx = start_idx + fold_size if i < cv_folds - 1 else n_samples
        fold_indices = indices[start_idx:end_idx]
        folds.append(fold_indices)

    # Cross-validation
    cv_results = defaultdict(list)

    for fold_idx, test_indices in enumerate(folds):
        # Create train/test split
        train_indices = np.concatenate([fold for i, fold in enumerate(folds) if i != fold_idx])

        X_train, X_test = X[train_indices], X[test_indices]
        y_train, y_test = y[train_indices], y[test_indices]

        # Train and evaluate
        classifier = BiologicalClassifier(algorithm=algorithm, random_state=random_state)
        classifier.fit(X_train, y_train)

        metrics = evaluate_classifier(classifier, X_test, y_test)

        # Store metrics
        for metric_name, metric_value in metrics.items():
            cv_results[metric_name].append(metric_value)

    return dict(cv_results)
