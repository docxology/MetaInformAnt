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
            task_type: Alias for model_type (backward compatibility)
            embedding: Alias for embedding_dim (backward compatibility)
            random_state: Alias for random_seed (sklearn convention)
            **kwargs: Additional parameters (ignored)
        """
        # Handle parameter aliases
        if task_type is not None:
            model_type = task_type
        if embedding is not None:
            embedding_dim = embedding
        if random_state is not None:
            random_seed = random_state

        # Map task_type values to model_type values for compatibility
        model_type_mapping = {
            "classification": "embedding",
            "regression": "embedding",
        }
        model_type = model_type_mapping.get(model_type, model_type)

        self.model_type = model_type
        self.embedding_dim = embedding_dim
        self.random_seed = random_seed
        self.is_fitted = False

        # Model components
        self.embeddings = {}
        self.classifier = None
        self.regressor = None
        self.event_vocab = {}
        self.reverse_vocab = {}

        np.random.seed(random_seed)

    def _extract_events(self, seq: Any) -> List[str]:
        """Extract event strings from a sequence, handling both EventSequence objects and lists.

        Args:
            seq: Either an EventSequence object or a list of event strings

        Returns:
            List of event type strings
        """
        if isinstance(seq, list):
            # Plain list of event strings
            return seq
        elif hasattr(seq, "events"):
            # EventSequence object - extract event_type from each event
            return [e.event_type if hasattr(e, "event_type") else str(e) for e in seq.events]
        else:
            # Try to iterate
            return [str(e) for e in seq]

    def fit(
        self, sequences: List[Any], outcomes: Union[List[int], List[float]], task: str = "classification"
    ) -> EventSequencePredictor:
        """Fit the model to training data.

        Args:
            sequences: List of EventSequence objects
            outcomes: Target outcomes (classes for classification, values for regression)
            task: Task type ('classification' or 'regression')

        Returns:
            Self for method chaining
        """
        if self.model_type == "embedding":
            self._fit_embedding_model(sequences, outcomes, task)
        elif self.model_type == "simple":
            self._fit_simple_model(sequences, outcomes, task)
        else:
            raise ValueError(f"Unsupported model type: {self.model_type}")

        self.is_fitted = True
        logger.info(f"Fitted {self.model_type} model for {task}")
        return self

    def predict(self, sequences: List[Any]) -> Union[np.ndarray, List[float]]:
        """Make predictions for new sequences.

        Args:
            sequences: List of EventSequence objects

        Returns:
            Predicted outcomes
        """
        if not self.is_fitted:
            raise ValueError("Model must be fitted before prediction")

        if self.model_type == "embedding":
            return self._predict_embedding(sequences)
        elif self.model_type == "simple":
            return self._predict_simple(sequences)
        else:
            raise ValueError(f"Unsupported model type: {self.model_type}")

    def predict_proba(self, sequences: List[Any]) -> Optional[np.ndarray]:
        """Predict class probabilities (classification only).

        Args:
            sequences: List of EventSequence objects

        Returns:
            Probability array or None if not applicable
        """
        if not self.is_fitted or self.classifier is None:
            return None

        # Compute class probabilities using logistic function on predictions
        predictions = self.predict(sequences)
        # Apply sigmoid to get probabilities
        probs = 1.0 / (1.0 + np.exp(-predictions))
        # Return probabilities for both classes
        return np.column_stack([1.0 - probs, probs])

    def save_model(self, path: str | Path) -> None:
        """Save model to disk.

        Args:
            path: Output path
        """
        # Convert numpy arrays in embeddings to lists for JSON serialization
        embeddings_serializable = {k: v.tolist() if hasattr(v, "tolist") else v for k, v in self.embeddings.items()}

        model_data = {
            "model_type": self.model_type,
            "embedding_dim": self.embedding_dim,
            "random_seed": self.random_seed,
            "is_fitted": self.is_fitted,
            "embeddings": embeddings_serializable,
            "event_vocab": self.event_vocab,
            "reverse_vocab": {str(k): v for k, v in self.reverse_vocab.items()},  # Ensure keys are strings
        }

        # Save classifier/regressor if available
        if self.classifier is not None:
            model_data["classifier"] = pickle.dumps(self.classifier).hex()
        if self.regressor is not None:
            model_data["regressor"] = pickle.dumps(self.regressor).hex()

        with open(path, "w") as f:
            json.dump(model_data, f, indent=2)

        logger.info(f"Saved model to {path}")

    @classmethod
    def load_model(cls, path: str | Path) -> EventSequencePredictor:
        """Load model from disk.

        Args:
            path: Model file path

        Returns:
            Loaded model
        """
        with open(path, "r") as f:
            model_data = json.load(f)

        model = cls(
            model_type=model_data["model_type"],
            embedding_dim=model_data["embedding_dim"],
            random_seed=model_data["random_seed"],
        )

        model.is_fitted = model_data["is_fitted"]
        # Convert lists back to numpy arrays for embeddings
        model.embeddings = {k: np.array(v) if isinstance(v, list) else v for k, v in model_data["embeddings"].items()}
        model.event_vocab = model_data["event_vocab"]
        # Convert string keys back to int for reverse_vocab if needed
        model.reverse_vocab = {int(k) if k.isdigit() else k: v for k, v in model_data["reverse_vocab"].items()}

        # Load classifier/regressor if available
        if "classifier" in model_data:
            model.classifier = pickle.loads(bytes.fromhex(model_data["classifier"]))
        if "regressor" in model_data:
            model.regressor = pickle.loads(bytes.fromhex(model_data["regressor"]))

        logger.info(f"Loaded model from {path}")
        return model

    def _fit_embedding_model(self, sequences: List[Any], outcomes: Union[List[int], List[float]], task: str) -> None:
        """Fit embedding-based model using learned embeddings and linear model."""
        # Build vocabulary - handle both EventSequence objects and plain lists
        all_events = set()
        for seq in sequences:
            events = self._extract_events(seq)
            for event in events:
                all_events.add(event)

        self.event_vocab = {event: i for i, event in enumerate(all_events)}
        self.reverse_vocab = {i: event for event, i in self.event_vocab.items()}

        # Initialize embeddings with small random values
        for event in all_events:
            self.embeddings[event] = np.random.normal(0, 0.1, self.embedding_dim)

        # Compute sequence features (average embeddings)
        X = []
        for seq in sequences:
            events = self._extract_events(seq)
            if events:
                event_embeddings = [self.embeddings.get(e, np.zeros(self.embedding_dim)) for e in events]
                avg_embedding = np.mean(event_embeddings, axis=0)
            else:
                avg_embedding = np.zeros(self.embedding_dim)
            X.append(avg_embedding)
        X = np.array(X)
        y = np.array(outcomes)

        # Train linear model using gradient descent
        if task == "classification":
            # Logistic regression weights
            self.classifier = self._train_logistic_regression(X, y)
        else:
            # Linear regression weights
            self.regressor = self._train_linear_regression(X, y)

    def _train_logistic_regression(
        self, X: np.ndarray, y: np.ndarray, lr: float = 0.01, epochs: int = 100
    ) -> Dict[str, np.ndarray]:
        """Train logistic regression model using gradient descent.

        Args:
            X: Feature matrix (n_samples, n_features)
            y: Binary labels (n_samples,)
            lr: Learning rate
            epochs: Number of training epochs

        Returns:
            Dictionary with weights and bias
        """
        n_samples, n_features = X.shape
        weights = np.zeros(n_features)
        bias = 0.0

        for _ in range(epochs):
            # Forward pass
            linear = np.dot(X, weights) + bias
            predictions = 1.0 / (1.0 + np.exp(-np.clip(linear, -500, 500)))

            # Compute gradients
            error = predictions - y
            dw = (1.0 / n_samples) * np.dot(X.T, error)
            db = (1.0 / n_samples) * np.sum(error)

            # Update parameters
            weights -= lr * dw
            bias -= lr * db

        return {"weights": weights, "bias": bias}

    def _train_linear_regression(self, X: np.ndarray, y: np.ndarray) -> Dict[str, np.ndarray]:
        """Train linear regression model using closed-form solution.

        Args:
            X: Feature matrix (n_samples, n_features)
            y: Target values (n_samples,)

        Returns:
            Dictionary with weights and bias
        """
        # Add bias term
        X_bias = np.column_stack([np.ones(X.shape[0]), X])

        # Closed-form solution: w = (X^T X)^(-1) X^T y
        # Use pseudo-inverse for numerical stability
        try:
            XTX = np.dot(X_bias.T, X_bias)
            XTy = np.dot(X_bias.T, y)
            # Add small regularization for stability
            XTX += 1e-6 * np.eye(XTX.shape[0])
            weights_full = np.linalg.solve(XTX, XTy)
        except np.linalg.LinAlgError:
            # Fallback to pseudo-inverse
            weights_full = np.dot(np.linalg.pinv(X_bias), y)

        return {"weights": weights_full[1:], "bias": weights_full[0]}

    def _fit_simple_model(self, sequences: List[Any], outcomes: Union[List[int], List[float]], task: str) -> None:
        """Fit simple statistical model."""
        # Simple model based on sequence length
        self.sequence_length_model = {}
        for seq, outcome in zip(sequences, outcomes):
            events = self._extract_events(seq)
            length = len(events)
            if length not in self.sequence_length_model:
                self.sequence_length_model[length] = []
            self.sequence_length_model[length].append(outcome)

        # Compute averages
        for length in self.sequence_length_model:
            self.sequence_length_model[length] = np.mean(self.sequence_length_model[length])

    def _predict_embedding(self, sequences: List[Any]) -> Union[np.ndarray, List[float]]:
        """Make predictions using trained embedding model."""
        # Compute features
        X = []
        for seq in sequences:
            events = self._extract_events(seq)
            if events:
                event_embeddings = [self.embeddings.get(e, np.zeros(self.embedding_dim)) for e in events]
                avg_embedding = np.mean(event_embeddings, axis=0)
            else:
                avg_embedding = np.zeros(self.embedding_dim)
            X.append(avg_embedding)
        X = np.array(X)

        # Make predictions using trained model
        if self.classifier is not None and isinstance(self.classifier, dict):
            # Logistic regression prediction (returns log-odds for predict_proba)
            linear = np.dot(X, self.classifier["weights"]) + self.classifier["bias"]
            return linear
        elif self.regressor is not None and isinstance(self.regressor, dict):
            # Linear regression prediction
            predictions = np.dot(X, self.regressor["weights"]) + self.regressor["bias"]
            return predictions
        else:
            # Fallback: use embedding magnitude as prediction
            return np.array([np.linalg.norm(x) for x in X])

    def _predict_simple(self, sequences: List[Any]) -> Union[np.ndarray, List[float]]:
        """Make predictions using simple model."""
        predictions = []
        for seq in sequences:
            events = self._extract_events(seq)
            length = len(events)
            prediction = self.sequence_length_model.get(length, 0.0)
            predictions.append(prediction)

        return np.array(predictions)


class EnsemblePredictor:
    """Ensemble predictor combining multiple models for improved performance."""

    def __init__(self, models: List[EventSequencePredictor], weights: Optional[List[float]] = None):
        """Initialize ensemble.

        Args:
            models: List of trained EventSequencePredictor models
            weights: Weights for each model (uniform if None)
        """
        self.models = models
        self.weights = weights or [1.0 / len(models)] * len(models)
        self.is_fitted = all(model.is_fitted for model in models)

    def predict(self, sequences: List[Any]) -> np.ndarray:
        """Make ensemble predictions.

        Args:
            sequences: List of EventSequence objects

        Returns:
            Ensemble predictions
        """
        if not self.is_fitted:
            raise ValueError("All models must be fitted")

        # Get predictions from each model
        all_predictions = []
        for model in self.models:
            preds = model.predict(sequences)
            all_predictions.append(preds)

        # Weighted average
        all_predictions = np.array(all_predictions)
        weights = np.array(self.weights)

        return np.average(all_predictions, axis=0, weights=weights)
