"""Embedding-based prediction models for life event sequences.

This module provides the core outcome-prediction models for life events:
- EventSequencePredictor: Main predictor using embedding or simple (sequence-length) approaches
- EnsemblePredictor: Combines multiple predictors via weighted average

Event embeddings are learned in ``metainformant.life_events.models.embeddings``. Recurrent
next-event models (LSTM/GRU, which require PyTorch) live in
``metainformant.life_events.models.sequence_models``; they do not predict outcomes and are not
available through :class:`EventSequencePredictor`.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Dict, List, Optional, Union

import numpy as np

from metainformant.core.utils import logging

logger = logging.get_logger(__name__)


def _unsupported_model_payload(component: str) -> str:
    """Message for a model file whose ``component`` is not a JSON object of arrays and numbers."""
    return (
        f"The model file's {component!r} is not a JSON object. Models are saved as JSON by save_model; a "
        "serialised (pickled) payload is not loaded, because unpickling a file runs whatever code it contains."
    )


def _unsupported_model_type(model_type: str) -> ValueError:
    """Build the error raised for a model type that EventSequencePredictor does not implement."""
    message = f"Unsupported model type: {model_type}. EventSequencePredictor supports 'embedding' and 'simple'."
    if model_type in {"lstm", "gru"}:
        message += (
            " Recurrent models predict the next event rather than an outcome; use LSTMSequenceModel or "
            "GRUSequenceModel from metainformant.life_events.models.sequence_models (requires PyTorch)."
        )
    return ValueError(message)


class EventSequencePredictor:
    """Predictor for life event sequence outcomes using various ML approaches.

    Supports two model types:

    - ``"embedding"``: averages the embeddings of a sequence's events and fits a logistic
      (classification, binary outcomes) or linear (regression) model on them.
    - ``"simple"``: predicts the mean training outcome of sequences with the same number of events
      (the nearest observed length for unseen lengths). For classification the outcome is the
      indicator of the positive class (the largest label), so the model is binary and its score is
      the probability of that class.

    Recurrent (LSTM/GRU) models are not outcome predictors; see
    ``metainformant.life_events.models.sequence_models``.
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
            model_type: Type of model ('embedding' or 'simple')
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
        self.classifier: Optional[Dict[str, Any]] = None
        self.regressor: Optional[Dict[str, Any]] = None
        self.event_vocab: Dict[str, int] = {}
        self.reverse_vocab: Dict[int, str] = {}
        self.classes_: Optional[np.ndarray] = None
        self.sequence_length_model: Optional[Dict[int, Any]] = None

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
        outcomes: Union[List[int], List[float], List[str]],
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

        if self.model_type not in {"embedding", "simple"}:
            raise _unsupported_model_type(self.model_type)
        if len(sequences) == 0:
            raise ValueError("sequences cannot be empty")
        if len(sequences) != len(outcomes):
            raise ValueError(f"sequences and outcomes must have the same length ({len(sequences)} != {len(outcomes)})")

        if self.model_type == "embedding":
            self._fit_embedding_model(sequences, outcomes, task, event_embeddings=event_embeddings)
        else:
            self._fit_simple_model(sequences, outcomes, task)

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
            raise _unsupported_model_type(self.model_type)

        # For classification, convert the raw scores to the probability of the positive class and
        # threshold it at 0.5 (the largest label is the positive class).
        if self.task_type == "classification" and self.classes_ is not None and self.classifier is not None:
            probs = self._positive_class_probability(raw)
            predictions = np.where(probs >= 0.5, self.classes_[-1], self.classes_[0])
            return np.asarray(predictions)

        return raw

    def _positive_class_probability(self, raw: np.ndarray) -> np.ndarray:
        """Convert raw model scores to the probability of the positive class.

        The embedding model's scores are logits and pass through the logistic function. The simple
        model's scores already are probabilities (mean positive-class indicator per sequence length)
        and must not be squashed again, which would map every non-negative score above 0.5.
        """
        if self.model_type == "simple":
            return np.asarray(np.clip(raw, 0.0, 1.0))
        return np.asarray(1.0 / (1.0 + np.exp(-np.clip(raw, -500, 500))))

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

        if self.model_type == "simple":
            raw = self._predict_simple(sequences)
        else:
            raw = self._predict_embedding_raw(sequences)
        probs = self._positive_class_probability(raw)
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

        if self.sequence_length_model is not None:
            # JSON object keys are strings; load_model converts them back to integer lengths.
            model_data["sequence_length_model"] = {
                str(length): float(value) for length, value in self.sequence_length_model.items()
            }

        if self.classes_ is not None:
            model_data["classes_"] = self.classes_.tolist()

        if self.classifier is not None:
            model_data["classifier"] = {
                k: v.tolist() if hasattr(v, "tolist") else v for k, v in self.classifier.items()
            }
        if self.regressor is not None:
            model_data["regressor"] = {k: v.tolist() if hasattr(v, "tolist") else v for k, v in self.regressor.items()}

        with open(path, "w") as f:
            json.dump(model_data, f, indent=2)

        logger.info(f"Saved model to {path}")

    @classmethod
    def load_model(cls, path: str | Path) -> EventSequencePredictor:
        """Load model from disk.

        Raises:
            FileNotFoundError: If path doesn't exist
            ValueError: If JSON is invalid, required fields are missing, or the classifier / regressor is a
                serialised (pickled) payload instead of the JSON object ``save_model`` writes
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
            if not isinstance(clf_data, dict):
                raise ValueError(_unsupported_model_payload("classifier"))
            model.classifier = {k: np.array(v) if isinstance(v, list) else v for k, v in clf_data.items()}
        if "regressor" in model_data:
            reg_data = model_data["regressor"]
            if not isinstance(reg_data, dict):
                raise ValueError(_unsupported_model_payload("regressor"))
            model.regressor = {k: np.array(v) if isinstance(v, list) else v for k, v in reg_data.items()}

        if model.model_type == "simple":
            model.sequence_length_model = cls._read_sequence_length_model(model_data)
            if model.classifier is not None:
                model.classifier = {"type": "simple", "model": model.sequence_length_model}

        logger.info(f"Loaded model from {path}")
        return model

    @staticmethod
    def _read_sequence_length_model(model_data: Dict[str, Any]) -> Dict[int, float]:
        """Read the simple model's length -> mean outcome table from a saved model.

        Raises:
            ValueError: If the file holds no table, so a model that would silently predict a
                constant is never returned.
        """
        table = model_data.get("sequence_length_model")
        if table is None:
            # Classification files written before the table was saved separately keep it here.
            classifier = model_data.get("classifier")
            if isinstance(classifier, dict) and isinstance(classifier.get("model"), dict):
                table = classifier["model"]
        if not isinstance(table, dict) or not table:
            raise ValueError(
                "Simple model file contains no sequence-length table (it was written by a version that did "
                "not save it); re-train the model and save it again"
            )
        return {int(length): float(value) for length, value in table.items()}

    def _fit_embedding_model(
        self,
        sequences: List[Any],
        outcomes: Union[List[int], List[float], List[str]],
        task: str,
        event_embeddings: Optional[Dict[str, np.ndarray]] = None,
    ) -> None:
        """Fit embedding-based model."""
        unique_events = set()
        for seq in sequences:
            events = self._extract_events(seq)
            for event in events:
                unique_events.add(event)

        # Sorted, not in set order: the order of a set of strings changes with the interpreter's hash seed, and
        # it decided both the vocabulary indices and which random vector each event received
        all_events = sorted(unique_events)
        self.event_vocab = {event: i for i, event in enumerate(all_events)}
        self.reverse_vocab = {i: event for event, i in self.event_vocab.items()}

        # Use pre-computed or randomly initialized embeddings. The generator is private and seeded from
        # random_seed at every fit, so the same seed gives the same embeddings however many predictors were built
        # before (the module-level numpy generator was reseeded in __init__ and drawn from here).
        rng = np.random.default_rng(self.random_seed)
        for event in all_events:
            if event_embeddings and event in event_embeddings:
                self.embeddings[event] = event_embeddings[event]
            else:
                self.embeddings[event] = rng.normal(0, 0.1, self.embedding_dim)

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
        X_arr = np.asarray(X, dtype=float)
        y = np.asarray(outcomes)

        if task == "classification":
            self.classifier = self._train_logistic_regression(X_arr, y)
        else:
            self.regressor = self._train_linear_regression(X_arr, y)

    def _train_logistic_regression(
        self, X: np.ndarray, y: np.ndarray, lr: float = 0.01, epochs: int = 100
    ) -> Dict[str, Any]:
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

        return {"weights": weights, "bias": float(bias)}

    def _train_linear_regression(self, X: np.ndarray, y: np.ndarray) -> Dict[str, Any]:
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

    def _fit_simple_model(
        self, sequences: List[Any], outcomes: Union[List[int], List[float], List[str]], task: str
    ) -> None:
        """Fit the sequence-length model: the mean outcome of the training sequences of each length.

        For classification the outcome is the indicator of the positive class (the largest label), so
        the stored value is the probability of that class for sequences of that length.
        """
        # Build vocab for simple model
        all_events = set()
        for seq in sequences:
            events = self._extract_events(seq)
            for event in events:
                all_events.add(event)
        self.event_vocab = {event: i for i, event in enumerate(all_events)}

        outcome_array = np.asarray(outcomes)
        if task == "classification":
            classes = np.unique(outcome_array)
            if classes.size > 2:
                raise ValueError(
                    f"The simple model supports binary classification only, got {classes.size} classes; "
                    "use model_type='embedding' or regression"
                )
            targets = (outcome_array == classes[-1]).astype(float)
        else:
            targets = outcome_array.astype(float)

        outcomes_by_length: Dict[int, List[float]] = {}
        for seq, target in zip(sequences, targets):
            outcomes_by_length.setdefault(len(self._extract_events(seq)), []).append(float(target))
        self.sequence_length_model = {length: float(np.mean(values)) for length, values in outcomes_by_length.items()}

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
        X_arr = np.asarray(X, dtype=float)

        if self.classifier is not None and isinstance(self.classifier, dict) and "weights" in self.classifier:
            linear = np.dot(X_arr, self.classifier["weights"]) + self.classifier["bias"]
            return np.asarray(linear, dtype=float)
        elif self.regressor is not None and isinstance(self.regressor, dict):
            predictions = np.dot(X_arr, self.regressor["weights"]) + self.regressor["bias"]
            return np.asarray(predictions, dtype=float)
        else:
            return np.asarray([np.linalg.norm(x) for x in X_arr], dtype=float)

    def _predict_simple(self, sequences: List[Any]) -> np.ndarray:
        """Predict the mean training outcome for each sequence's length.

        Lengths never seen in training use the nearest observed length (the shorter one on a tie).
        """
        table = self.sequence_length_model
        if not table:
            raise ValueError("Simple model has no sequence-length table; fit the model or load a fitted one")
        known_lengths = np.array(sorted(table))
        predictions = []
        for seq in sequences:
            length = len(self._extract_events(seq))
            if length not in table:
                length = int(known_lengths[np.argmin(np.abs(known_lengths - length))])
            predictions.append(table[length])
        return np.array(predictions, dtype=float)


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

        all_predictions_arr = np.asarray([np.asarray(model.predict(sequences), dtype=float) for model in self.models])
        weights = np.asarray(self.weights, dtype=float)

        return np.average(all_predictions_arr, axis=0, weights=weights)
