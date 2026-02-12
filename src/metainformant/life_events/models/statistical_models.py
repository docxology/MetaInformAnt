"""Statistical and multi-task prediction models for life events.

This module provides advanced prediction models:
- MultiTaskPredictor: Simultaneous prediction of multiple outcome types
- SurvivalPredictor: Survival analysis using Cox, exponential, or baseline methods
- attention_weights: Extract attention weights from trained models
"""

from __future__ import annotations

from typing import Any, Dict, List, Optional, Tuple, Union

import numpy as np

from metainformant.core.utils import logging

# Forward reference for type hints
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from metainformant.life_events.core.events import EventSequence

logger = logging.get_logger(__name__)


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

        elif hasattr(self.model, "eval"):  # PyTorch model
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
