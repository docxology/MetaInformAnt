"""Embedding-based prediction models for life event sequences.

This module provides the core prediction models for life events:
- EventSequencePredictor: Main predictor using embedding or simple approaches
- EnsemblePredictor: Combines multiple predictors via weighted average
- biological_embedding: Learn Word2Vec-style embeddings for events
- domain_specific_embeddings: Learn per-domain embeddings
"""

from __future__ import annotations

import json
import pickle
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, Union

import numpy as np

from metainformant.core.utils import logging

logger = logging.get_logger(__name__)


class EventSequencePredictor:
    """Predictor for life event sequence outcomes using various ML approaches.

    Supports multiple model types: embedding-based, simple statistical, and LSTM models.
    """

    def __init__(
        self,
        model_type: str = "embedding",
        embedding_dim: int = 100,
        random_seed: int = 42,
        # Parameter aliases for compatibility
        task_type: str | None = None,
        embedding: int | None = None,
        random_state: int | None = None,
        **kwargs: Any,
    ):
        """Initialize predictor.

        Args:
            model_type: Type of model ('embedding', 'simple', 'lstm')
            embedding_dim: Dimensionality for embeddings
            random_seed: Random seed for reproducibility
            task_type: Task type ('classification' or 'regression')
            embedding: Alias for embedding_dim (backward compatibility)
            random_state: Alias for random_seed (sklearn convention)
            **kwargs: Additional parameters (ignored)
        """
        if embedding is not None:
            embedding_dim = embedding
        if random_state is not None:
            random_seed = random_state

        # Store task_type separately — do NOT conflate with model_type
        self.task_type = task_type  # 'classification', 'regression', or None
        self.model_type = model_type
        self.embedding_dim = embedding_dim
        self.random_seed = random_seed
        self.is_fitted = False

        # Model components
        self.embeddings: Dict[str, np.ndarray] = {}
        self.classifier: Optional[Dict[str, np.ndarray]] = None
        self.regressor: Optional[Dict[str, np.ndarray]] = None
        self.event_vocab: Dict[str, int] = {}
        self.reverse_vocab: Dict[int, str] = {}
        self.classes_: Optional[np.ndarray] = None
        self.sequence_length_model: Optional[Dict[int, float]] = None

        np.random.seed(random_seed)

    def _extract_events(self, seq: Any) -> List[str]:
        """Extract event strings from a sequence."""
        if isinstance(seq, list):
            return seq
        elif hasattr(seq, "events"):
            return [e.event_type if hasattr(e, "event_type") else str(e) for e in seq.events]
        else:
            return [str(e) for e in seq]

    def fit(
        self,
        sequences: List[Any],
        outcomes: Union[List[int], List[float]],
        task: str | None = None,
        event_embeddings: Optional[Dict[str, np.ndarray]] = None,
    ) -> EventSequencePredictor:
        """Fit the model to training data.

        Args:
            sequences: List of EventSequence objects or lists of event strings
            outcomes: Target outcomes (classes for classification, values for regression)
            task: Task type ('classification' or 'regression'). Overrides task_type.
            event_embeddings: Optional pre-computed event embeddings.

        Returns:
            Self for method chaining
        """
        # Determine task type from arguments or init parameter
        if task is None:
            task = self.task_type or "classification"
        self.task_type = task

        if self.model_type == "embedding":
            self._fit_embedding_model(sequences, outcomes, task, event_embeddings=event_embeddings)
        elif self.model_type == "simple":
            self._fit_simple_model(sequences, outcomes, task)
        else:
            raise ValueError(f"Unsupported model type: {self.model_type}")

        self.is_fitted = True
        # Store unique classes for classification tasks (sklearn convention)
        if task == "classification":
            self.classes_ = np.unique(outcomes)
        logger.info(f"Fitted {self.model_type} model for {task}")
        return self

    def predict(self, sequences: List[Any]) -> Union[np.ndarray, List[float]]:
        """Make predictions for new sequences.

        For classification, returns class labels. For regression, returns values.
        """
        if not self.is_fitted:
            raise ValueError("Model must be fitted before prediction")

        if self.model_type == "embedding":
            raw = self._predict_embedding_raw(sequences)
        elif self.model_type == "simple":
            raw = self._predict_simple(sequences)
        else:
            raise ValueError(f"Unsupported model type: {self.model_type}")

        # For classification, convert raw scores to class labels
        if self.task_type == "classification" and self.classes_ is not None and self.classifier is not None:
            probs = 1.0 / (1.0 + np.exp(-np.clip(raw, -500, 500)))
            predictions = np.where(probs >= 0.5, self.classes_[-1], self.classes_[0])
            return predictions

        return raw

    def predict_proba(self, sequences: List[Any]) -> np.ndarray:
        """Predict class probabilities (classification only).

        Args:
            sequences: List of EventSequence objects

        Returns:
            Probability array

        Raises:
            ValueError: If task is not classification
        """
        if self.task_type != "classification":
            raise ValueError("predict_proba is only available for classification tasks")
        if not self.is_fitted or self.classifier is None:
            raise ValueError("Model must be fitted for classification before calling predict_proba")

        raw = self._predict_embedding_raw(sequences)
        probs = 1.0 / (1.0 + np.exp(-np.clip(raw, -500, 500)))
        return np.column_stack([1.0 - probs, probs])

    def save_model(self, path: str | Path) -> None:
        """Save model to disk.

        Raises:
            ValueError: If model is not fitted
        """
        if not self.is_fitted:
            raise ValueError("Model must be fitted before saving")

        embeddings_serializable = {k: v.tolist() if hasattr(v, "tolist") else v for k, v in self.embeddings.items()}

        model_data: Dict[str, Any] = {
            "model_type": self.model_type,
            "task_type": self.task_type or "classification",
            "embedding_dim": self.embedding_dim,
            "random_seed": self.random_seed,
            "is_fitted": self.is_fitted,
            "event_embeddings": embeddings_serializable,
            "vocab": self.event_vocab,
            "event_vocab": self.event_vocab,
            "reverse_vocab": {str(k): v for k, v in self.reverse_vocab.items()},
        }

        if self.classes_ is not None:
            model_data["classes_"] = self.classes_.tolist()

        if self.classifier is not None:
            model_data["classifier"] = {
                k: v.tolist() if hasattr(v, "tolist") else v for k, v in self.classifier.items()
            }
        if self.regressor is not None:
            model_data["regressor"] = {
                k: v.tolist() if hasattr(v, "tolist") else v for k, v in self.regressor.items()
            }

        with open(path, "w") as f:
            json.dump(model_data, f, indent=2)

        logger.info(f"Saved model to {path}")

    @classmethod
    def load_model(cls, path: str | Path) -> EventSequencePredictor:
        """Load model from disk.

        Raises:
            FileNotFoundError: If path doesn't exist
            ValueError: If JSON is invalid or required fields are missing
        """
        path = Path(path)
        if not path.exists():
            raise FileNotFoundError(f"Model file not found: {path}")

        try:
            with open(path, "r") as f:
                model_data = json.load(f)
        except json.JSONDecodeError as e:
            raise ValueError(f"Failed to parse model file: {e}")

        # Validate required fields
        required_fields = {"model_type", "is_fitted"}
        missing = required_fields - set(model_data.keys())
        if missing:
            raise ValueError(f"Model file missing fields: {missing}")

        if not model_data.get("is_fitted", False):
            raise ValueError("Cannot load unfitted model")

        model = cls(
            model_type=model_data["model_type"],
            embedding_dim=model_data.get("embedding_dim", 100),
            random_seed=model_data.get("random_seed", 42),
            task_type=model_data.get("task_type"),
        )

        model.is_fitted = model_data["is_fitted"]
        model.task_type = model_data.get("task_type", "classification")

        # Load embeddings
        embeddings_data = model_data.get("event_embeddings", model_data.get("embeddings", {}))
        model.embeddings = {k: np.array(v) if isinstance(v, list) else v for k, v in embeddings_data.items()}

        model.event_vocab = model_data.get("vocab", model_data.get("event_vocab", {}))
        model.reverse_vocab = {int(k) if k.isdigit() else k: v for k, v in model_data.get("reverse_vocab", {}).items()}

        if "classes_" in model_data:
            model.classes_ = np.array(model_data["classes_"])

        # Load classifier/regressor
        if "classifier" in model_data:
            clf_data = model_data["classifier"]
            if isinstance(clf_data, dict):
                model.classifier = {k: np.array(v) if isinstance(v, list) else v for k, v in clf_data.items()}
            else:
                model.classifier = pickle.loads(bytes.fromhex(clf_data))
        if "regressor" in model_data:
            reg_data = model_data["regressor"]
            if isinstance(reg_data, dict):
                model.regressor = {k: np.array(v) if isinstance(v, list) else v for k, v in reg_data.items()}
            else:
                model.regressor = pickle.loads(bytes.fromhex(reg_data))

        logger.info(f"Loaded model from {path}")
        return model

    def _fit_embedding_model(
        self,
        sequences: List[Any],
        outcomes: Union[List[int], List[float]],
        task: str,
        event_embeddings: Optional[Dict[str, np.ndarray]] = None,
    ) -> None:
        """Fit embedding-based model."""
        all_events = set()
        for seq in sequences:
            events = self._extract_events(seq)
            for event in events:
                all_events.add(event)

        self.event_vocab = {event: i for i, event in enumerate(all_events)}
        self.reverse_vocab = {i: event for event, i in self.event_vocab.items()}

        # Use pre-computed or randomly initialized embeddings
        if event_embeddings:
            for event in all_events:
                if event in event_embeddings:
                    self.embeddings[event] = event_embeddings[event]
                else:
                    self.embeddings[event] = np.random.normal(0, 0.1, self.embedding_dim)
        else:
            for event in all_events:
                self.embeddings[event] = np.random.normal(0, 0.1, self.embedding_dim)

        # Compute sequence features (average embeddings)
        X = []
        for seq in sequences:
            events = self._extract_events(seq)
            if events:
                event_embs = [self.embeddings.get(e, np.zeros(self.embedding_dim)) for e in events]
                avg_embedding = np.mean(event_embs, axis=0)
            else:
                avg_embedding = np.zeros(self.embedding_dim)
            X.append(avg_embedding)
        X = np.array(X)
        y = np.array(outcomes)

        if task == "classification":
            self.classifier = self._train_logistic_regression(X, y)
        else:
            self.regressor = self._train_linear_regression(X, y)

    def _train_logistic_regression(
        self, X: np.ndarray, y: np.ndarray, lr: float = 0.01, epochs: int = 100
    ) -> Dict[str, np.ndarray]:
        """Train logistic regression model using gradient descent."""
        n_samples, n_features = X.shape
        weights = np.zeros(n_features)
        bias = 0.0

        for _ in range(epochs):
            linear = np.dot(X, weights) + bias
            predictions = 1.0 / (1.0 + np.exp(-np.clip(linear, -500, 500)))
            error = predictions - y
            dw = (1.0 / n_samples) * np.dot(X.T, error)
            db = (1.0 / n_samples) * np.sum(error)
            weights -= lr * dw
            bias -= lr * db

        return {"weights": weights, "bias": np.float64(bias)}

    def _train_linear_regression(self, X: np.ndarray, y: np.ndarray) -> Dict[str, np.ndarray]:
        """Train linear regression model using closed-form solution."""
        X_bias = np.column_stack([np.ones(X.shape[0]), X])
        try:
            XTX = np.dot(X_bias.T, X_bias)
            XTy = np.dot(X_bias.T, y)
            XTX += 1e-6 * np.eye(XTX.shape[0])
            weights_full = np.linalg.solve(XTX, XTy)
        except np.linalg.LinAlgError:
            weights_full = np.dot(np.linalg.pinv(X_bias), y)

        return {"weights": weights_full[1:], "bias": weights_full[0]}

    def _fit_simple_model(self, sequences: List[Any], outcomes: Union[List[int], List[float]], task: str) -> None:
        """Fit simple statistical model."""
        # Build vocab for simple model
        all_events = set()
        for seq in sequences:
            events = self._extract_events(seq)
            for event in events:
                all_events.add(event)
        self.event_vocab = {event: i for i, event in enumerate(all_events)}

        self.sequence_length_model = {}
        for seq, outcome in zip(sequences, outcomes):
            events = self._extract_events(seq)
            length = len(events)
            if length not in self.sequence_length_model:
                self.sequence_length_model[length] = []
            self.sequence_length_model[length].append(outcome)

        for length in self.sequence_length_model:
            self.sequence_length_model[length] = np.mean(self.sequence_length_model[length])

        # For simple classification, store a basic classifier reference
        if task == "classification":
            self.classifier = {"type": "simple", "model": self.sequence_length_model}

    def _predict_embedding_raw(self, sequences: List[Any]) -> np.ndarray:
        """Make raw (un-thresholded) predictions using trained embedding model."""
        X = []
        for seq in sequences:
            events = self._extract_events(seq)
            if events:
                event_embs = [self.embeddings.get(e, np.zeros(self.embedding_dim)) for e in events]
                avg_embedding = np.mean(event_embs, axis=0)
            else:
                avg_embedding = np.zeros(self.embedding_dim)
            X.append(avg_embedding)
        X = np.array(X)

        if self.classifier is not None and isinstance(self.classifier, dict) and "weights" in self.classifier:
            linear = np.dot(X, self.classifier["weights"]) + self.classifier["bias"]
            return linear
        elif self.regressor is not None and isinstance(self.regressor, dict):
            predictions = np.dot(X, self.regressor["weights"]) + self.regressor["bias"]
            return predictions
        else:
            return np.array([np.linalg.norm(x) for x in X])

    def _predict_simple(self, sequences: List[Any]) -> np.ndarray:
        """Make predictions using simple model."""
        predictions = []
        for seq in sequences:
            events = self._extract_events(seq)
            length = len(events)
            if self.sequence_length_model:
                prediction = self.sequence_length_model.get(length, 0.0)
            else:
                prediction = 0.0
            predictions.append(prediction)
        return np.array(predictions)


class EnsemblePredictor:
    """Ensemble predictor combining multiple models for improved performance."""

    def __init__(
        self,
        models: List[EventSequencePredictor],
        weights: Optional[List[float]] = None,
        task_type: str | None = None,
    ):
        """Initialize ensemble.

        Args:
            models: List of trained EventSequencePredictor models
            weights: Weights for each model (uniform if None)
            task_type: Task type for the ensemble
        """
        self.models = models
        self.weights = weights or [1.0 / len(models)] * len(models)
        self.is_fitted = all(model.is_fitted for model in models)
        self.task_type = task_type

    def predict(self, sequences: List[Any]) -> np.ndarray:
        """Make ensemble predictions."""
        if not self.is_fitted:
            raise ValueError("All models must be fitted")

        all_predictions = []
        for model in self.models:
            preds = model.predict(sequences)
            all_predictions.append(preds)

        all_predictions = np.array(all_predictions)
        weights = np.array(self.weights)

        return np.average(all_predictions, axis=0, weights=weights)
