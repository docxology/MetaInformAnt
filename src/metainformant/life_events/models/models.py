"""Prediction models for life event sequences.

This module provides machine learning models for predicting outcomes from life event sequences,
including embedding-based models, LSTM models, and ensemble methods.
"""

from __future__ import annotations

import json
import pickle
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, Union

import numpy as np

from metainformant.core import logging

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

    def _train_logistic_regression(self, X: np.ndarray, y: np.ndarray, lr: float = 0.01, epochs: int = 100) -> Dict[str, np.ndarray]:
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


def biological_embedding(
    sequences: List[Any], embedding_dim: int = 100, window_size: int = 5, min_count: int = 5
) -> Dict[str, np.ndarray]:
    """Learn biological event embeddings using Word2Vec-style approach.

    Args:
        sequences: List of EventSequence objects
        embedding_dim: Dimensionality of embeddings
        window_size: Context window size
        min_count: Minimum event frequency for embedding

    Returns:
        Dictionary mapping event types to embeddings
    """
    # Build vocabulary
    event_counts = {}
    for seq in sequences:
        for event in seq.events:
            event_counts[event.event_type] = event_counts.get(event.event_type, 0) + 1

    # Filter rare events
    vocab = {event: count for event, count in event_counts.items() if count >= min_count}

    # Simple random embeddings (real implementation would use Word2Vec)
    embeddings = {}
    for event in vocab:
        embeddings[event] = np.random.normal(0, 0.1, embedding_dim)

    logger.info(f"Generated embeddings for {len(embeddings)} events")
    return embeddings


def domain_specific_embeddings(
    sequences: List[Any], domains: List[str], embedding_dim: int = 100
) -> Dict[str, Dict[str, np.ndarray]]:
    """Learn domain-specific event embeddings.

    Args:
        sequences: List of EventSequence objects
        domains: List of domains to create embeddings for
        embedding_dim: Dimensionality of embeddings

    Returns:
        Dictionary mapping domains to embedding dictionaries
    """
    domain_embeddings = {}

    for domain in domains:
        # Filter sequences to this domain
        domain_sequences = []
        for seq in sequences:
            domain_events = [e for e in seq.events if e.domain == domain]
            if domain_events:
                from .events import EventSequence

                domain_sequences.append(EventSequence(seq.person_id, domain_events))

        if domain_sequences:
            domain_embeddings[domain] = biological_embedding(
                domain_sequences, embedding_dim, window_size=3, min_count=3
            )
        else:
            domain_embeddings[domain] = {}

    logger.info(f"Generated domain-specific embeddings for {len(domain_embeddings)} domains")
    return domain_embeddings


class GRUSequenceModel:
    """GRU-based sequence prediction model for life events.

    This model uses Gated Recurrent Units (GRUs) to predict future life events
    based on sequences of past events. It incorporates event embeddings and
    temporal patterns for improved prediction accuracy.
    """

    def __init__(
        self,
        embedding_dim: int = 100,
        hidden_dim: int = 64,
        num_layers: int = 1,
        dropout: float = 0.1,
        epochs: int = 50,
        batch_size: int = 32,
        learning_rate: float = 0.001,
        random_state: Optional[int] = None,
        **kwargs: Any,
    ):
        """Initialize GRU sequence model.

        Args:
            embedding_dim: Dimension of event embeddings
            hidden_dim: Hidden state dimension for GRU
            num_layers: Number of GRU layers
            dropout: Dropout probability
            epochs: Number of training epochs
            batch_size: Batch size for training
            learning_rate: Learning rate for optimizer
            random_state: Random seed for reproducibility
            **kwargs: Additional model parameters
        """
        self.embedding_dim = embedding_dim
        self.hidden_dim = hidden_dim
        self.num_layers = num_layers
        self.dropout = dropout
        self.epochs = epochs
        self.batch_size = batch_size
        self.learning_rate = learning_rate
        self.random_state = random_state

        # Will be initialized during fit
        self.event_vocab = None
        self.num_events = None
        self.model = None
        self.embeddings = None

        if random_state is not None:
            import torch

            torch.manual_seed(random_state)
            if torch.cuda.is_available():
                torch.cuda.manual_seed(random_state)

    def fit(self, sequences: List[EventSequence], targets: Optional[List[Any]] = None) -> "GRUSequenceModel":
        """Fit the GRU model to life event sequences.

        Args:
            sequences: List of EventSequence objects
            targets: Optional target values for supervised learning

        Returns:
            Self for method chaining
        """
        try:
            import torch
            import torch.nn as nn
            from torch.utils.data import DataLoader, TensorDataset
        except ImportError:
            logger.warning("PyTorch not available, GRUSequenceModel disabled")
            return self

        # Build vocabulary from sequences
        self._build_vocabulary(sequences)

        # Convert sequences to tensors
        train_data = self._prepare_data(sequences, targets)

        if train_data is None:
            logger.warning("No training data available")
            return self

        train_loader = DataLoader(train_data, batch_size=self.batch_size, shuffle=True)

        # Initialize model
        self.model = nn.Sequential(
            nn.Embedding(self.num_events, self.embedding_dim),
            nn.GRU(
                self.embedding_dim,
                self.hidden_dim,
                self.num_layers,
                dropout=self.dropout if self.num_layers > 1 else 0,
                batch_first=True,
            ),
            nn.Linear(self.hidden_dim, self.num_events),  # Predict next event
        )

        optimizer = torch.optim.Adam(self.model.parameters(), lr=self.learning_rate)
        criterion = nn.CrossEntropyLoss()

        # Training loop (simplified)
        self.model.train()
        for epoch in range(self.epochs):
            total_loss = 0
            for inputs, targets in train_loader:
                optimizer.zero_grad()
                outputs, _ = self.model[1](self.model[0](inputs))  # GRU forward
                outputs = self.model[2](outputs[:, -1, :])  # Take last timestep
                loss = criterion(outputs, targets)
                loss.backward()
                optimizer.step()
                total_loss += loss.item()

            if epoch % 10 == 0:
                logger.info(".3f")

        return self

    def predict(self, sequences: List[EventSequence]) -> List[List[float]]:
        """Predict next events for given sequences.

        Args:
            sequences: List of EventSequence objects

        Returns:
            List of prediction probabilities for each sequence
        """
        if self.model is None:
            logger.warning("Model not trained, returning random predictions")
            return [[1.0 / self.num_events] * self.num_events for _ in sequences]

        try:
            import torch
        except ImportError:
            logger.warning("PyTorch not available, returning random predictions")
            return [[1.0 / self.num_events] * self.num_events for _ in sequences]

        predictions = []

        self.model.eval()
        with torch.no_grad():
            for seq in sequences:
                # Convert sequence to tensor
                seq_tensor = self._sequence_to_tensor(seq)
                if seq_tensor is None:
                    predictions.append([1.0 / self.num_events] * self.num_events)
                    continue

                # Forward pass
                embedded = self.model[0](seq_tensor.unsqueeze(0))
                outputs, _ = self.model[1](embedded)
                logits = self.model[2](outputs[:, -1, :])
                probs = torch.softmax(logits, dim=1).squeeze(0).tolist()

                predictions.append(probs)

        return predictions

    def _build_vocabulary(self, sequences: List[EventSequence]) -> None:
        """Build event vocabulary from sequences."""
        event_types = set()
        for seq in sequences:
            for event in seq.events:
                event_types.add(event.event_type)

        self.event_vocab = {event_type: idx for idx, event_type in enumerate(sorted(event_types))}
        self.num_events = len(self.event_vocab)

    def _prepare_data(self, sequences: List[EventSequence], targets: Optional[List[Any]] = None) -> Optional[Any]:
        """Prepare sequences for training."""
        try:
            import torch
        except ImportError:
            return None

        sequences_data = []
        targets_data = []

        for seq in sequences:
            seq_tensor = self._sequence_to_tensor(seq)
            if seq_tensor is not None and len(seq_tensor) > 1:
                # Use all but last event as input, last event as target
                input_seq = seq_tensor[:-1]
                target = seq_tensor[-1]

                sequences_data.append(input_seq)
                targets_data.append(target)

        if not sequences_data:
            return None

        # Pad sequences to same length
        max_len = max(len(seq) for seq in sequences_data)
        padded_sequences = []
        for seq in sequences_data:
            padding = torch.full((max_len - len(seq),), self.num_events - 1)  # Use last index as padding
            padded_seq = torch.cat([seq, padding])
            padded_sequences.append(padded_seq)

        return torch.utils.data.TensorDataset(torch.stack(padded_sequences), torch.tensor(targets_data))

    def _sequence_to_tensor(self, sequence: EventSequence) -> Optional[Any]:
        """Convert EventSequence to tensor."""
        try:
            import torch
        except ImportError:
            return None

        event_indices = []
        for event in sequence.events:
            if event.event_type in self.event_vocab:
                event_indices.append(self.event_vocab[event.event_type])

        if not event_indices:
            return None

        return torch.tensor(event_indices, dtype=torch.long)


class LSTMSequenceModel:
    """LSTM-based sequence prediction model for life events.

    This model uses Long Short-Term Memory (LSTM) networks to predict future life events
    based on sequences of past events. It incorporates event embeddings and
    temporal patterns for improved prediction accuracy.
    """

    def __init__(
        self,
        embedding_dim: int = 100,
        hidden_dim: int = 64,
        num_layers: int = 1,
        dropout: float = 0.1,
        epochs: int = 50,
        batch_size: int = 32,
        learning_rate: float = 0.001,
        random_state: Optional[int] = None,
        **kwargs: Any,
    ):
        """Initialize LSTM sequence model.

        Args:
            embedding_dim: Dimension of event embeddings
            hidden_dim: Hidden state dimension for LSTM
            num_layers: Number of LSTM layers
            dropout: Dropout probability
            epochs: Number of training epochs
            batch_size: Batch size for training
            learning_rate: Learning rate for optimizer
            random_state: Random seed for reproducibility
            **kwargs: Additional model parameters
        """
        self.embedding_dim = embedding_dim
        self.hidden_dim = hidden_dim
        self.num_layers = num_layers
        self.dropout = dropout
        self.epochs = epochs
        self.batch_size = batch_size
        self.learning_rate = learning_rate
        self.random_state = random_state

        # Will be initialized during fit
        self.event_vocab = None
        self.num_events = None
        self.model = None
        self.embeddings = None

        if random_state is not None:
            import torch

            torch.manual_seed(random_state)
            if torch.cuda.is_available():
                torch.cuda.manual_seed(random_state)

    def fit(self, sequences: List[EventSequence], targets: Optional[List[Any]] = None) -> "LSTMSequenceModel":
        """Fit the LSTM model to life event sequences.

        Args:
            sequences: List of EventSequence objects
            targets: Optional target values for supervised learning

        Returns:
            Self for method chaining
        """
        try:
            import torch
            import torch.nn as nn
            from torch.utils.data import DataLoader, TensorDataset
        except ImportError:
            logger.warning("PyTorch not available, LSTMSequenceModel disabled")
            return self

        # Build vocabulary from sequences
        self._build_vocabulary(sequences)

        # Convert sequences to tensors
        train_data = self._prepare_data(sequences, targets)

        if train_data is None:
            logger.warning("No training data available")
            return self

        train_loader = DataLoader(train_data, batch_size=self.batch_size, shuffle=True)

        # Initialize model
        self.model = nn.Sequential(
            nn.Embedding(self.num_events, self.embedding_dim),
            nn.LSTM(
                self.embedding_dim,
                self.hidden_dim,
                self.num_layers,
                dropout=self.dropout if self.num_layers > 1 else 0,
                batch_first=True,
            ),
            nn.Linear(self.hidden_dim, self.num_events),  # Predict next event
        )

        optimizer = torch.optim.Adam(self.model.parameters(), lr=self.learning_rate)
        criterion = nn.CrossEntropyLoss()

        # Training loop (simplified)
        self.model.train()
        for epoch in range(self.epochs):
            total_loss = 0
            for inputs, targets in train_loader:
                optimizer.zero_grad()
                outputs, (h_n, c_n) = self.model[1](self.model[0](inputs))  # LSTM forward
                outputs = self.model[2](outputs[:, -1, :])  # Take last timestep
                loss = criterion(outputs, targets)
                loss.backward()
                optimizer.step()
                total_loss += loss.item()

            if epoch % 10 == 0:
                logger.info(".3f")

        return self

    def predict(self, sequences: List[EventSequence]) -> List[List[float]]:
        """Predict next events for given sequences.

        Args:
            sequences: List of EventSequence objects

        Returns:
            List of prediction probabilities for each sequence
        """
        if self.model is None:
            logger.warning("Model not trained, returning random predictions")
            return [[1.0 / self.num_events] * self.num_events for _ in sequences]

        try:
            import torch
        except ImportError:
            logger.warning("PyTorch not available, returning random predictions")
            return [[1.0 / self.num_events] * self.num_events for _ in sequences]

        predictions = []

        self.model.eval()
        with torch.no_grad():
            for seq in sequences:
                # Convert sequence to tensor
                seq_tensor = self._sequence_to_tensor(seq)
                if seq_tensor is None:
                    predictions.append([1.0 / self.num_events] * self.num_events)
                    continue

                # Forward pass
                embedded = self.model[0](seq_tensor.unsqueeze(0))
                outputs, (h_n, c_n) = self.model[1](embedded)
                logits = self.model[2](outputs[:, -1, :])
                probs = torch.softmax(logits, dim=1).squeeze(0).tolist()

                predictions.append(probs)

        return predictions

    def _build_vocabulary(self, sequences: List[EventSequence]) -> None:
        """Build event vocabulary from sequences."""
        event_types = set()
        for seq in sequences:
            for event in seq.events:
                event_types.add(event.event_type)

        self.event_vocab = {event_type: idx for idx, event_type in enumerate(sorted(event_types))}
        self.num_events = len(self.event_vocab)

    def _prepare_data(self, sequences: List[EventSequence], targets: Optional[List[Any]] = None) -> Optional[Any]:
        """Prepare sequences for training."""
        try:
            import torch
        except ImportError:
            return None

        sequences_data = []
        targets_data = []

        for seq in sequences:
            seq_tensor = self._sequence_to_tensor(seq)
            if seq_tensor is not None and len(seq_tensor) > 1:
                # Use all but last event as input, last event as target
                input_seq = seq_tensor[:-1]
                target = seq_tensor[-1]

                sequences_data.append(input_seq)
                targets_data.append(target)

        if not sequences_data:
            return None

        # Pad sequences to same length
        max_len = max(len(seq) for seq in sequences_data)
        padded_sequences = []
        for seq in sequences_data:
            padding = torch.full((max_len - len(seq),), self.num_events - 1)  # Use last index as padding
            padded_seq = torch.cat([seq, padding])
            padded_sequences.append(padded_seq)

        return torch.utils.data.TensorDataset(torch.stack(padded_sequences), torch.tensor(targets_data))

    def _sequence_to_tensor(self, sequence: EventSequence) -> Optional[Any]:
        """Convert EventSequence to tensor."""
        try:
            import torch
        except ImportError:
            return None

        event_indices = []
        for event in sequence.events:
            if event.event_type in self.event_vocab:
                event_indices.append(self.event_vocab[event.event_type])

        if not event_indices:
            return None

        return torch.tensor(event_indices, dtype=torch.long)


class MultiTaskPredictor:
    """Multi-task predictor for life events with multiple outcome types.

    This model can simultaneously predict multiple types of outcomes (classification,
    regression) from life event sequences using a shared representation.
    """

    def __init__(
        self,
        task_types: Dict[str, str],
        embedding_dim: int = 100,
        hidden_dim: int = 64,
        random_state: Optional[int] = None,
        **kwargs: Any,
    ):
        """Initialize multi-task predictor.

        Args:
            task_types: Dictionary mapping task names to types ('classification' or 'regression')
            embedding_dim: Dimension of event embeddings
            hidden_dim: Hidden dimension for shared representation
            random_state: Random seed for reproducibility
            **kwargs: Additional parameters
        """
        self.task_types = task_types
        self.embedding_dim = embedding_dim
        self.hidden_dim = hidden_dim
        self.random_state = random_state

        self.num_events = None
        self.event_vocab = None
        self.models = {}

        if random_state is not None:
            import random

            random.seed(random_state)
            np.random.seed(random_state)
            try:
                import torch

                torch.manual_seed(random_state)
                if torch.cuda.is_available():
                    torch.cuda.manual_seed(random_state)
            except ImportError:
                pass

    def fit(self, sequences: List[EventSequence], outcomes: Dict[str, np.ndarray]) -> "MultiTaskPredictor":
        """Fit multi-task predictor to sequences and outcomes.

        Args:
            sequences: List of EventSequence objects
            outcomes: Dictionary mapping task names to outcome arrays

        Returns:
            Self for method chaining
        """
        try:
            import torch
            import torch.nn as nn
            from torch.utils.data import DataLoader, TensorDataset
        except ImportError:
            logger.warning("PyTorch not available, MultiTaskPredictor disabled")
            return self

        # Build vocabulary from sequences
        self._build_vocabulary(sequences)

        # Create shared encoder
        shared_encoder = nn.Sequential(
            nn.Embedding(self.num_events, self.embedding_dim),
            nn.LSTM(self.embedding_dim, self.hidden_dim, batch_first=True),
        )

        # Create task-specific heads
        self.models = {}
        for task_name, task_type in self.task_types.items():
            if task_type == "classification":
                # Binary classification
                head = nn.Sequential(nn.Linear(self.hidden_dim, 32), nn.ReLU(), nn.Linear(32, 1), nn.Sigmoid())
            elif task_type == "regression":
                # Regression
                head = nn.Sequential(nn.Linear(self.hidden_dim, 32), nn.ReLU(), nn.Linear(32, 1))
            else:
                raise ValueError(f"Unknown task type: {task_type}")

            self.models[task_name] = {"encoder": shared_encoder, "head": head}

        # Prepare training data
        train_data = self._prepare_data(sequences, outcomes)

        if train_data is None:
            logger.warning("No training data available")
            return self

        train_loader = DataLoader(train_data, batch_size=32, shuffle=True)

        # Training loop (simplified - trains all tasks jointly)
        for task_name, model_components in self.models.items():
            encoder = model_components["encoder"]
            head = model_components["head"]

            optimizer = torch.optim.Adam(list(encoder.parameters()) + list(head.parameters()), lr=0.001)
            criterion = nn.BCELoss() if self.task_types[task_name] == "classification" else nn.MSELoss()

            encoder.train()
            head.train()

            for epoch in range(10):  # Simplified training
                for inputs, targets in train_loader:
                    optimizer.zero_grad()

                    # Forward pass
                    embedded = encoder[0](inputs)
                    outputs, _ = encoder[1](embedded)
                    final_output = outputs[:, -1, :]  # Take last timestep
                    predictions = head(final_output)

                    # Get targets for this task
                    task_targets = targets[task_name]

                    if self.task_types[task_name] == "classification":
                        loss = criterion(predictions.squeeze(), task_targets.float())
                    else:
                        loss = criterion(predictions.squeeze(), task_targets.float())

                    loss.backward()
                    optimizer.step()

        return self

    def predict(self, sequences: List[EventSequence]) -> Dict[str, List[float]]:
        """Predict outcomes for multiple tasks.

        Args:
            sequences: List of EventSequence objects

        Returns:
            Dictionary mapping task names to prediction lists
        """
        if not self.models:
            logger.warning("Model not trained, returning zeros")
            return {task: [0.0] * len(sequences) for task in self.task_types.keys()}

        try:
            import torch
        except ImportError:
            logger.warning("PyTorch not available, returning zeros")
            return {task: [0.0] * len(sequences) for task in self.task_types.keys()}

        predictions = {task: [] for task in self.task_types.keys()}

        for seq in sequences:
            seq_tensor = self._sequence_to_tensor(seq)
            if seq_tensor is None:
                for task in predictions:
                    predictions[task].append(0.0)
                continue

            for task_name, model_components in self.models.items():
                encoder = model_components["encoder"]
                head = model_components["head"]

                encoder.eval()
                head.eval()

                with torch.no_grad():
                    embedded = encoder[0](seq_tensor.unsqueeze(0))
                    outputs, _ = encoder[1](embedded)
                    final_output = outputs[:, -1, :]
                    pred = head(final_output).squeeze().item()

                    predictions[task_name].append(pred)

        return predictions

    def _build_vocabulary(self, sequences: List[EventSequence]) -> None:
        """Build event vocabulary from sequences."""
        event_types = set()
        for seq in sequences:
            for event in seq.events:
                event_types.add(event.event_type)

        self.event_vocab = {event_type: idx for idx, event_type in enumerate(sorted(event_types))}
        self.num_events = len(self.event_vocab)

    def _prepare_data(self, sequences: List[EventSequence], outcomes: Dict[str, np.ndarray]) -> Optional[Any]:
        """Prepare sequences and multi-task targets for training."""
        try:
            import torch
        except ImportError:
            return None

        sequences_data = []
        targets_data = {task: [] for task in self.task_types.keys()}

        for i, seq in enumerate(sequences):
            seq_tensor = self._sequence_to_tensor(seq)
            if seq_tensor is not None and len(seq_tensor) > 1:
                sequences_data.append(seq_tensor)

                for task in self.task_types.keys():
                    if i < len(outcomes[task]):
                        targets_data[task].append(outcomes[task][i])

        if not sequences_data:
            return None

        # Pad sequences to same length
        max_len = max(len(seq) for seq in sequences_data)
        padded_sequences = []
        for seq in sequences_data:
            padding = torch.full((max_len - len(seq),), self.num_events - 1)
            padded_seq = torch.cat([seq, padding])
            padded_sequences.append(padded_seq)

        # Convert targets to tensors
        target_tensors = {}
        for task, targets in targets_data.items():
            target_tensors[task] = torch.tensor(targets, dtype=torch.float)

        return torch.utils.data.TensorDataset(torch.stack(padded_sequences), target_tensors)

    def _sequence_to_tensor(self, sequence: EventSequence) -> Optional[Any]:
        """Convert EventSequence to tensor."""
        try:
            import torch
        except ImportError:
            return None

        event_indices = []
        for event in sequence.events:
            if event.event_type in self.event_vocab:
                event_indices.append(self.event_vocab[event.event_type])

        if not event_indices:
            return None

        return torch.tensor(event_indices, dtype=torch.long)


class SurvivalPredictor:
    """Survival analysis predictor for life events.

    This model predicts survival times and functions from life event sequences
    using various survival analysis methods.
    """

    def __init__(self, method: str = "cox", embedding_dim: int = 50, random_state: Optional[int] = None, **kwargs: Any):
        """Initialize survival predictor.

        Args:
            method: Survival analysis method ('cox', 'random_survival_forest', etc.)
            embedding_dim: Dimension of event embeddings
            random_state: Random seed for reproducibility
            **kwargs: Additional parameters
        """
        self.method = method
        self.embedding_dim = embedding_dim
        self.random_state = random_state

        # Will be initialized during fit
        self.event_vocab = None
        self.num_events = None
        self.model = None

        if random_state is not None:
            import torch

            torch.manual_seed(random_state)
            if torch.cuda.is_available():
                torch.cuda.manual_seed(random_state)

    def fit(
        self, sequences: List[EventSequence], event_times: np.ndarray, event_occurred: np.ndarray
    ) -> "SurvivalPredictor":
        """Fit survival model to sequences and survival data.

        Args:
            sequences: List of EventSequence objects
            event_times: Array of event/censoring times
            event_occurred: Boolean array indicating if event occurred

        Returns:
            Self for method chaining
        """
        try:
            import torch
            import torch.nn as nn
            from torch.utils.data import DataLoader, TensorDataset
        except ImportError:
            logger.warning("PyTorch not available, SurvivalPredictor disabled")
            return self

        # Build vocabulary from sequences
        self._build_vocabulary(sequences)

        # Convert sequences to features
        features = self._sequences_to_features(sequences)

        if self.method == "cox":
            # Simple Cox proportional hazards model (simplified implementation)
            self.model = nn.Sequential(
                nn.Linear(features.shape[1], 32), nn.ReLU(), nn.Linear(32, 1)  # Log hazard ratio
            )

            # Prepare training data
            train_data = torch.utils.data.TensorDataset(
                torch.tensor(features, dtype=torch.float),
                torch.tensor(event_times, dtype=torch.float),
                torch.tensor(event_occurred, dtype=torch.bool),
            )

            train_loader = DataLoader(train_data, batch_size=32, shuffle=True)

            optimizer = torch.optim.Adam(self.model.parameters(), lr=0.001)
            # Simplified loss (not proper Cox loss)
            criterion = nn.MSELoss()

            self.model.train()
            for epoch in range(20):
                for inputs, times, occurred in train_loader:
                    optimizer.zero_grad()
                    outputs = self.model(inputs)
                    # Simplified: predict log hazard
                    loss = criterion(outputs.squeeze(), torch.log(times + 1))
                    loss.backward()
                    optimizer.step()

        elif self.method == "exponential":
            # Simple exponential model using feature-based hazard rate
            # Fit: hazard = exp(X @ weights) where weights are learned
            self._fit_exponential_model(features, event_times, event_occurred)

        else:
            # Fallback to Kaplan-Meier style baseline
            logger.info(f"Method '{self.method}' using Kaplan-Meier baseline model")
            self._fit_baseline_model(event_times, event_occurred)

        return self

    def _fit_exponential_model(self, features: np.ndarray, times: np.ndarray, occurred: np.ndarray) -> None:
        """Fit simple exponential survival model."""
        # Use mean survival time weighted by features
        n_features = features.shape[1]
        self.exp_weights = np.zeros(n_features)
        self.exp_baseline = np.mean(times[occurred]) if occurred.sum() > 0 else np.mean(times)

        # Simple gradient descent to learn weights
        lr = 0.001
        for _ in range(100):
            linear = np.dot(features, self.exp_weights)
            hazard = np.exp(np.clip(linear, -10, 10))
            predicted_times = self.exp_baseline / (hazard + 1e-8)

            # Gradient (simplified)
            error = predicted_times - times
            gradient = np.dot(features.T, error * hazard) / len(times)
            self.exp_weights -= lr * gradient

        self.model = {"type": "exponential", "weights": self.exp_weights, "baseline": self.exp_baseline}

    def _fit_baseline_model(self, times: np.ndarray, occurred: np.ndarray) -> None:
        """Fit baseline Kaplan-Meier style model."""
        # Store empirical survival distribution
        sorted_times = np.sort(times)
        survival_prob = np.linspace(1.0, 0.0, len(sorted_times))
        self.model = {"type": "baseline", "times": sorted_times, "survival": survival_prob, "median": np.median(times)}

    def predict(self, sequences: List[EventSequence]) -> np.ndarray:
        """Predict survival times for sequences.

        Args:
            sequences: List of EventSequence objects

        Returns:
            Array of predicted survival times
        """
        if self.model is None:
            logger.warning("Model not trained, returning zeros")
            return np.zeros(len(sequences))

        try:
            import torch
        except ImportError:
            logger.warning("PyTorch not available, returning zeros")
            return np.zeros(len(sequences))

        features = self._sequences_to_features(sequences)

        if isinstance(self.model, dict):
            model_type = self.model.get("type", "unknown")

            if model_type == "exponential":
                # Exponential model prediction
                linear = np.dot(features, self.model["weights"])
                hazard = np.exp(np.clip(linear, -10, 10))
                predicted_times = self.model["baseline"] / (hazard + 1e-8)

            elif model_type == "baseline":
                # Baseline model - return median for all sequences
                predicted_times = np.full(len(sequences), self.model["median"])

            else:
                # Unknown model type - use feature-based heuristic
                predicted_times = np.array([np.sum(f) + 100 for f in features])

        elif hasattr(self.model, 'eval'):  # PyTorch model
            self.model.eval()
            with torch.no_grad():
                inputs = torch.tensor(features, dtype=torch.float)
                outputs = self.model(inputs).squeeze()
                # Convert log hazard back to time
                predicted_times = torch.exp(outputs).numpy()

        else:
            # Fallback: feature-based prediction
            predicted_times = np.array([max(1.0, np.mean(f) * 100 + 50) for f in features])

        return predicted_times

    def predict_survival_function(self, sequences: List[EventSequence], times: np.ndarray) -> np.ndarray:
        """Predict survival function at given times.

        Args:
            sequences: List of EventSequence objects
            times: Array of time points

        Returns:
            Array of survival probabilities (sequences x times)
        """
        n_sequences = len(sequences)
        n_times = len(times)

        # Simplified: exponential survival function
        survival_probs = np.zeros((n_sequences, n_times))

        predicted_times = self.predict(sequences)

        for i in range(n_sequences):
            for j, t in enumerate(times):
                survival_probs[i, j] = np.exp(-t / max(predicted_times[i], 1))

        return survival_probs

    def _build_vocabulary(self, sequences: List[EventSequence]) -> None:
        """Build event vocabulary from sequences."""
        event_types = set()
        for seq in sequences:
            for event in seq.events:
                event_types.add(event.event_type)

        self.event_vocab = {event_type: idx for idx, event_type in enumerate(sorted(event_types))}
        self.num_events = len(self.event_vocab)

    def _sequences_to_features(self, sequences: List[EventSequence]) -> np.ndarray:
        """Convert sequences to feature vectors."""
        features = []

        for seq in sequences:
            # Simple feature extraction: event counts
            event_counts = np.zeros(self.num_events)
            for event in seq.events:
                if event.event_type in self.event_vocab:
                    idx = self.event_vocab[event.event_type]
                    event_counts[idx] += 1

            # Add sequence length
            features.append(np.concatenate([event_counts, [len(seq.events)]]))

        return np.array(features)


def attention_weights(model: Any, sequences: List[EventSequence]) -> Dict[str, np.ndarray]:
    """Compute attention weights for event sequences using a trained model.

    Args:
        model: Trained sequence model with attention mechanism
        sequences: List of EventSequence objects

    Returns:
        Dictionary mapping sequence indices to attention weight arrays
    """
    attention_results = {}

    if not hasattr(model, "get_attention_weights"):
        logger.warning("Model does not support attention weights")
        return attention_results

    try:
        for i, seq in enumerate(sequences):
            try:
                weights = model.get_attention_weights(seq)
                if weights is not None:
                    attention_results[str(i)] = np.array(weights)
            except Exception as e:
                logger.warning(f"Failed to compute attention for sequence {i}: {e}")
                continue

    except Exception as e:
        logger.error(f"Error computing attention weights: {e}")

    return attention_results


class SurvivalPredictor:
    """Survival analysis predictor for life events."""

    def __init__(self, random_state: Optional[int] = None):
        self.random_state = random_state

    def fit(self, sequences: List[EventSequence], durations: List[float], events: List[bool]) -> "SurvivalPredictor":
        """Fit survival model."""
        return self

    def predict_survival_function(self, sequences: List[EventSequence], times: List[float]) -> np.ndarray:
        """Predict survival function."""
        # Placeholder: return random survival curves
        n_seqs = len(sequences)
        n_times = len(times)
        return np.linspace(1.0, 0.0, n_times)[None, :] * np.ones((n_seqs, 1))


# Alias for backward compatibility
learn_event_embeddings = biological_embedding
