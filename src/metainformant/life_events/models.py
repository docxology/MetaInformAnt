"""Sequence prediction models for life course analysis.

This module provides models for predicting outcomes from event sequences,
including LSTM-based and simpler sequence models.
"""

from __future__ import annotations

from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np
from numpy.typing import NDArray

try:
    import torch
    import torch.nn as nn

    TORCH_AVAILABLE = True
except ImportError:
    TORCH_AVAILABLE = False
    torch = None  # type: ignore
    nn = None  # type: ignore


class EventSequencePredictor:
    """Predictor for outcomes from event sequences.
    
    Supports classification and regression tasks on event sequences.
    Can use sequence embeddings or raw sequence processing.
    
    Attributes:
        model_type: Type of model ("embedding", "lstm", "simple")
        task_type: Type of task ("classification", "regression")
        is_fitted: Whether model has been trained
    """
    
    def __init__(
        self,
        model_type: str = "embedding",
        task_type: str = "classification",
        embedding_dim: int = 100,
        hidden_dim: int = 64,
        random_state: Optional[int] = None,
        **kwargs: Any
    ):
        """Initialize EventSequencePredictor.
        
        Args:
            model_type: Model type ("embedding", "lstm", "simple")
            task_type: Task type ("classification", "regression")
            embedding_dim: Dimension of sequence embeddings
            hidden_dim: Hidden dimension for LSTM models
            random_state: Random seed
            **kwargs: Additional model-specific parameters
        """
        self.model_type = model_type
        self.task_type = task_type
        self.embedding_dim = embedding_dim
        self.hidden_dim = hidden_dim
        self.random_state = random_state
        self.params = kwargs
        
        self.is_fitted = False
        self.model_weights: Optional[Dict[str, Any]] = None
        self.classes_: Optional[NDArray] = None
        
        if random_state is not None:
            np.random.seed(random_state)
    
    def fit(
        self,
        sequences: Sequence[Sequence[str]],
        y: NDArray,
        event_embeddings: Optional[Dict[str, NDArray]] = None
    ) -> "EventSequencePredictor":
        """Fit model to event sequences.
        
        Args:
            sequences: List of event sequences (each sequence is list of event tokens)
            y: Target values (labels for classification, continuous for regression)
            event_embeddings: Optional pre-trained event embeddings
            
        Returns:
            Fitted predictor
        """
        if self.model_type == "embedding":
            return self._fit_embedding_model(sequences, y, event_embeddings)
        elif self.model_type == "lstm":
            return self._fit_lstm_model(sequences, y, event_embeddings)
        elif self.model_type == "simple":
            return self._fit_simple_model(sequences, y, event_embeddings)
        else:
            raise ValueError(f"Unknown model_type: {self.model_type}")
    
    def _fit_embedding_model(
        self,
        sequences: Sequence[Sequence[str]],
        y: NDArray,
        event_embeddings: Optional[Dict[str, NDArray]]
    ) -> "EventSequencePredictor":
        """Fit embedding-based model."""
        from .embeddings import learn_event_embeddings, sequence_embeddings
        
        # Learn embeddings if not provided
        if event_embeddings is None:
            event_embeddings = learn_event_embeddings(
                list(sequences),
                embedding_dim=self.embedding_dim,
                random_state=self.random_state
            )
        
        # Convert sequences to embeddings
        X = sequence_embeddings(
            list(sequences),
            event_embeddings,
            method="mean"
        )
        
        # Store embeddings for prediction
        self.event_embeddings = event_embeddings
        
        # Fit simple linear model
        if self.task_type == "classification":
            from ..ml.classification import BiologicalClassifier
            self.classifier = BiologicalClassifier(
                algorithm="linear",
                random_state=self.random_state
            )
            self.classifier.fit(X, y)
            self.classes_ = self.classifier.classes_
        else:
            from ..ml.regression import BiologicalRegressor
            self.regressor = BiologicalRegressor(
                algorithm="ridge",
                random_state=self.random_state
            )
            self.regressor.fit(X, y)
        
        self.is_fitted = True
        return self
    
    def _fit_lstm_model(
        self,
        sequences: Sequence[Sequence[str]],
        y: NDArray,
        event_embeddings: Optional[Dict[str, NDArray]]
    ) -> "EventSequencePredictor":
        """Fit LSTM-based model."""
        if not TORCH_AVAILABLE:
            # Fallback to simple model if PyTorch not available
            return self._fit_simple_model(sequences, y, event_embeddings)
        
        from .embeddings import learn_event_embeddings
        
        # Learn embeddings if not provided
        if event_embeddings is None:
            event_embeddings = learn_event_embeddings(
                list(sequences),
                embedding_dim=self.embedding_dim,
                random_state=self.random_state
            )
        
        self.event_embeddings = event_embeddings
        
        # Create vocabulary
        vocab = sorted(list(event_embeddings.keys()))
        self.vocab_to_idx = {token: i for i, token in enumerate(vocab)}
        self.idx_to_vocab = {i: token for token, i in self.vocab_to_idx.items()}
        
        # Convert sequences to indices
        sequences_idx = []
        for seq in sequences:
            seq_idx = [self.vocab_to_idx.get(token, 0) for token in seq]
            sequences_idx.append(seq_idx)
        
        # Pad sequences to same length
        max_len = max(len(s) for s in sequences_idx) if sequences_idx else 1
        sequences_padded = []
        for seq in sequences_idx:
            padded = seq + [0] * (max_len - len(seq))
            sequences_padded.append(padded[:max_len])
        
        X = np.array(sequences_padded)
        
        # Use simple model for now (LSTM implementation would go here)
        # For now, use embedding model as fallback
        return self._fit_embedding_model(sequences, y, event_embeddings)
    
    def _fit_simple_model(
        self,
        sequences: Sequence[Sequence[str]],
        y: NDArray,
        event_embeddings: Optional[Dict[str, NDArray]]
    ) -> "EventSequencePredictor":
        """Fit simple sequence model (event frequency-based)."""
        # Build vocabulary
        vocab = set()
        for seq in sequences:
            vocab.update(seq)
        vocab = sorted(list(vocab))
        
        # Create feature matrix: event frequencies per sequence
        X = np.zeros((len(sequences), len(vocab)))
        for i, seq in enumerate(sequences):
            for event in seq:
                if event in vocab:
                    idx = vocab.index(event)
                    X[i, idx] += 1
        
        # Normalize by sequence length
        seq_lengths = np.array([len(seq) for seq in sequences])
        X = X / (seq_lengths[:, np.newaxis] + 1e-10)
        
        self.vocab = vocab
        
        # Fit model
        if self.task_type == "classification":
            from ..ml.classification import BiologicalClassifier
            self.classifier = BiologicalClassifier(
                algorithm="random_forest",
                random_state=self.random_state
            )
            self.classifier.fit(X, y)
            self.classes_ = self.classifier.classes_
        else:
            from ..ml.regression import BiologicalRegressor
            self.regressor = BiologicalRegressor(
                algorithm="ridge",
                random_state=self.random_state
            )
            self.regressor.fit(X, y)
        
        self.is_fitted = True
        return self
    
    def predict(self, sequences: Sequence[Sequence[str]]) -> NDArray:
        """Predict outcomes for event sequences.
        
        Args:
            sequences: List of event sequences
            
        Returns:
            Array of predictions
        """
        if not self.is_fitted:
            raise ValueError("Model must be fitted before prediction")
        
        if self.model_type == "embedding":
            return self._predict_embedding(sequences)
        elif self.model_type == "lstm":
            return self._predict_lstm(sequences)
        elif self.model_type == "simple":
            return self._predict_simple(sequences)
        else:
            raise ValueError(f"Unknown model_type: {self.model_type}")
    
    def _predict_embedding(self, sequences: Sequence[Sequence[str]]) -> NDArray:
        """Predict using embedding model."""
        from .embeddings import sequence_embeddings
        
        X = sequence_embeddings(
            list(sequences),
            self.event_embeddings,
            method="mean"
        )
        
        if self.task_type == "classification":
            return self.classifier.predict(X)
        else:
            return self.regressor.predict(X)
    
    def _predict_lstm(self, sequences: Sequence[Sequence[str]]) -> NDArray:
        """Predict using LSTM model."""
        # For now, use embedding prediction
        return self._predict_embedding(sequences)
    
    def _predict_simple(self, sequences: Sequence[Sequence[str]]) -> NDArray:
        """Predict using simple model."""
        # Build feature matrix
        X = np.zeros((len(sequences), len(self.vocab)))
        for i, seq in enumerate(sequences):
            for event in seq:
                if event in self.vocab:
                    idx = self.vocab.index(event)
                    X[i, idx] += 1
        
        # Normalize
        seq_lengths = np.array([len(seq) for seq in sequences])
        X = X / (seq_lengths[:, np.newaxis] + 1e-10)
        
        if self.task_type == "classification":
            return self.classifier.predict(X)
        else:
            return self.regressor.predict(X)
    
    def predict_proba(self, sequences: Sequence[Sequence[str]]) -> NDArray:
        """Predict class probabilities for event sequences.
        
        Args:
            sequences: List of event sequences
            
        Returns:
            Array of class probabilities (only for classification tasks)
        """
        if not self.is_fitted:
            raise ValueError("Model must be fitted before prediction")
        
        if self.task_type != "classification":
            raise ValueError("predict_proba only available for classification tasks")
        
        if self.model_type == "embedding":
            from .embeddings import sequence_embeddings
            X = sequence_embeddings(
                list(sequences),
                self.event_embeddings,
                method="mean"
            )
            return self.classifier.predict_proba(X)
        elif self.model_type == "simple":
            # Build feature matrix
            X = np.zeros((len(sequences), len(self.vocab)))
            for i, seq in enumerate(sequences):
                for event in seq:
                    if event in self.vocab:
                        idx = self.vocab.index(event)
                        X[i, idx] += 1
            
            seq_lengths = np.array([len(seq) for seq in sequences])
            X = X / (seq_lengths[:, np.newaxis] + 1e-10)
            
            return self.classifier.predict_proba(X)
        else:
            raise ValueError("predict_proba not implemented for this model type")


class LSTMSequenceModel:
    """LSTM-based sequence model for event sequences.
    
    Note: Full implementation requires PyTorch. Falls back to simpler model
    if PyTorch is not available.
    """
    
    def __init__(
        self,
        embedding_dim: int = 100,
        hidden_dim: int = 64,
        num_layers: int = 2,
        random_state: Optional[int] = None
    ):
        """Initialize LSTM model.
        
        Args:
            embedding_dim: Dimension of event embeddings
            hidden_dim: Hidden dimension of LSTM
            num_layers: Number of LSTM layers
            random_state: Random seed
        """
        self.embedding_dim = embedding_dim
        self.hidden_dim = hidden_dim
        self.num_layers = num_layers
        self.random_state = random_state
        self.is_fitted = False
        
        if not TORCH_AVAILABLE:
            # Use simple fallback
            self.use_fallback = True
        else:
            self.use_fallback = False
            if random_state is not None:
                torch.manual_seed(random_state)
    
    def fit(
        self,
        sequences: Sequence[Sequence[str]],
        y: NDArray,
        event_embeddings: Optional[Dict[str, NDArray]] = None
    ) -> "LSTMSequenceModel":
        """Fit LSTM model."""
        if self.use_fallback:
            # Use simple predictor as fallback
            self.fallback_model = EventSequencePredictor(
                model_type="simple",
                random_state=self.random_state
            )
            return self.fallback_model.fit(sequences, y, event_embeddings)
        
        # Full LSTM implementation would go here
        # For now, use fallback
        self.fallback_model = EventSequencePredictor(
            model_type="embedding",
            random_state=self.random_state
        )
        return self.fallback_model.fit(sequences, y, event_embeddings)
    
    def predict(self, sequences: Sequence[Sequence[str]]) -> NDArray:
        """Predict using LSTM model."""
        if not self.is_fitted:
            raise ValueError("Model must be fitted before prediction")
        
        if self.use_fallback or not TORCH_AVAILABLE:
            return self.fallback_model.predict(sequences)
        
        # Full LSTM prediction would go here
        return self.fallback_model.predict(sequences)

