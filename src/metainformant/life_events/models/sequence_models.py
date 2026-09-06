"""Neural sequence prediction models for life events.

This module provides recurrent neural network models for predicting future
life events from sequences:
- GRUSequenceModel: GRU-based sequence prediction
- LSTMSequenceModel: LSTM-based sequence prediction

Both models share their implementation via ``_RNNSequenceModel`` and require
torch as an optional dependency.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Any, Dict, List, Optional, Self

from metainformant.core.utils import logging
from metainformant.life_events.core.utils import _sequence_to_tokens

if TYPE_CHECKING:
    from metainformant.life_events.core.events import EventSequence

logger = logging.get_logger(__name__)


class _RNNSequenceModel:
    """Shared implementation for RNN-based next-event sequence models.

    Subclasses supply the recurrent cell (GRU or LSTM) via ``_make_rnn``.
    """

    rnn_name = "RNN"

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
        """Initialize the sequence model.

        Args:
            embedding_dim: Dimension of event embeddings
            hidden_dim: Hidden state dimension for the recurrent layer
            num_layers: Number of recurrent layers
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
        self.event_vocab: Dict[str, int] | None = None
        self.num_events: int | None = None
        self.model: Any = None
        self.embeddings: Any = None
        self.is_fitted = False

        if random_state is not None:
            try:
                import torch

                torch.manual_seed(random_state)
                if torch.cuda.is_available():
                    torch.cuda.manual_seed(random_state)
            except ImportError:
                pass

    def _make_rnn(self, nn_module: Any) -> Any:
        """Create the recurrent layer for this model type."""
        raise NotImplementedError

    def fit(self, sequences: List[EventSequence], targets: Optional[List[Any]] = None) -> Self:
        """Fit the model to life event sequences.

        Args:
            sequences: List of EventSequence objects
            targets: Optional target values for supervised learning

        Returns:
            Self for method chaining
        """
        try:
            import torch
            import torch.nn as nn
            from torch.utils.data import DataLoader
        except ImportError:
            logger.warning("PyTorch not available, %s disabled", type(self).__name__)
            self.is_fitted = True
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
        if self.num_events is None:
            logger.warning("Vocabulary not built; cannot initialize model")
            return self
        num_events = self.num_events
        self.model = nn.Sequential(
            nn.Embedding(num_events, self.embedding_dim),
            self._make_rnn(nn),
            nn.Linear(self.hidden_dim, self.num_events),  # Predict next event
        )

        optimizer = torch.optim.Adam(self.model.parameters(), lr=self.learning_rate)
        criterion = nn.CrossEntropyLoss()

        # Training loop (simplified)
        self.model.train()
        for epoch in range(self.epochs):
            total_loss = 0
            for inputs, batch_targets in train_loader:
                optimizer.zero_grad()
                outputs, _ = self.model[1](self.model[0](inputs))  # Recurrent forward
                outputs = self.model[2](outputs[:, -1, :])  # Take last timestep
                loss = criterion(outputs, batch_targets)
                loss.backward()
                optimizer.step()
                total_loss += loss.item()

            if epoch % 10 == 0:
                logger.info("Epoch %d, loss: %.3f", epoch, total_loss)

        self.is_fitted = True
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
            n = self.num_events or 1
            return [[1.0 / n] * n for _ in sequences]

        try:
            import torch
        except ImportError:
            logger.warning("PyTorch not available, returning random predictions")
            n = self.num_events or 1
            return [[1.0 / n] * n for _ in sequences]

        predictions = []

        self.model.eval()
        with torch.no_grad():
            for seq in sequences:
                # Convert sequence to tensor
                seq_tensor = self._sequence_to_tensor(seq)
                if seq_tensor is None:
                    n_events = self.num_events or 0
                    predictions.append([1.0 / n_events] * n_events)
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
            event_types.update(_sequence_to_tokens(seq))

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
            padding = torch.full((max_len - len(seq),), (self.num_events or 0) - 1)  # Use last index as padding
            padded_seq = torch.cat([seq, padding])
            padded_sequences.append(padded_seq)

        return torch.utils.data.TensorDataset(torch.stack(padded_sequences), torch.tensor(targets_data))

    def _sequence_to_tensor(self, sequence: EventSequence) -> Optional[Any]:
        """Convert EventSequence to tensor."""
        try:
            import torch
        except ImportError:
            return None

        vocab = self.event_vocab or {}
        event_indices = []
        for event in _sequence_to_tokens(sequence):
            if event in vocab:
                event_indices.append(vocab[event])

        if not event_indices:
            return None

        return torch.tensor(event_indices, dtype=torch.long)


class GRUSequenceModel(_RNNSequenceModel):
    """GRU-based sequence prediction model for life events.

    This model uses Gated Recurrent Units (GRUs) to predict future life events
    based on sequences of past events. It incorporates event embeddings and
    temporal patterns for improved prediction accuracy.
    """

    rnn_name = "GRU"

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
        super().__init__(
            embedding_dim=embedding_dim,
            hidden_dim=hidden_dim,
            num_layers=num_layers,
            dropout=dropout,
            epochs=epochs,
            batch_size=batch_size,
            learning_rate=learning_rate,
            random_state=random_state,
            **kwargs,
        )

    def _make_rnn(self, nn_module: Any) -> Any:
        """Create the GRU recurrent layer."""
        return nn_module.GRU(
            self.embedding_dim,
            self.hidden_dim,
            self.num_layers,
            dropout=self.dropout if self.num_layers > 1 else 0,
            batch_first=True,
        )


class LSTMSequenceModel(_RNNSequenceModel):
    """LSTM-based sequence prediction model for life events.

    This model uses Long Short-Term Memory (LSTM) networks to predict future life events
    based on sequences of past events. It incorporates event embeddings and
    temporal patterns for improved prediction accuracy.
    """

    rnn_name = "LSTM"

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
        super().__init__(
            embedding_dim=embedding_dim,
            hidden_dim=hidden_dim,
            num_layers=num_layers,
            dropout=dropout,
            epochs=epochs,
            batch_size=batch_size,
            learning_rate=learning_rate,
            random_state=random_state,
            **kwargs,
        )

    def _make_rnn(self, nn_module: Any) -> Any:
        """Create the LSTM recurrent layer."""
        return nn_module.LSTM(
            self.embedding_dim,
            self.hidden_dim,
            self.num_layers,
            dropout=self.dropout if self.num_layers > 1 else 0,
            batch_first=True,
        )
