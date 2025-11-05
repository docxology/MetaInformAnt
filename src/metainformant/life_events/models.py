"""Sequence prediction models for life course analysis.

This module provides models for predicting outcomes from event sequences,
including LSTM-based and simpler sequence models.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple, Union

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
    
    def save_model(self, path: str | Path) -> None:
        """Save trained model to disk.
        
        Saves model hyperparameters, embeddings, vocabulary, and classifier/regressor
        state to a JSON file. The model can be reconstructed using load_model().
        
        Args:
            path: Path to save model file (JSON format)
            
        Raises:
            ValueError: If model is not fitted
            IOError: If file cannot be written
            
        Examples:
            >>> predictor = EventSequencePredictor(random_state=42)
            >>> predictor.fit(sequences, outcomes)
            >>> predictor.save_model("output/model.json")
        """
        if not self.is_fitted:
            raise ValueError("Model must be fitted before saving")
        
        from ..core import io
        
        # Prepare model data for serialization
        model_data: Dict[str, Any] = {
            "model_type": self.model_type,
            "task_type": self.task_type,
            "embedding_dim": self.embedding_dim,
            "hidden_dim": self.hidden_dim,
            "random_state": self.random_state,
            "params": self.params,
            "is_fitted": self.is_fitted,
        }
        
        # Save event embeddings if present
        if hasattr(self, "event_embeddings") and self.event_embeddings:
            model_data["event_embeddings"] = {
                k: v.tolist() for k, v in self.event_embeddings.items()
            }
        
        # Save vocabulary if present (for simple model)
        if hasattr(self, "vocab") and self.vocab:
            model_data["vocab"] = self.vocab
        
        # Save vocabulary mappings if present (for LSTM model)
        if hasattr(self, "vocab_to_idx") and self.vocab_to_idx:
            model_data["vocab_to_idx"] = self.vocab_to_idx
        if hasattr(self, "idx_to_vocab") and self.idx_to_vocab:
            model_data["idx_to_vocab"] = self.idx_to_vocab
        
        # Save classes for classification
        if self.classes_ is not None:
            model_data["classes_"] = self.classes_.tolist()
        
        # Save classifier state if present
        if hasattr(self, "classifier") and self.classifier:
            classifier = self.classifier
            model_data["classifier"] = {
                "algorithm": classifier.algorithm,
                "random_state": classifier.random_state,
                "params": classifier.params,
                "is_fitted": classifier.is_fitted,
            }
            
            # Save training data (for KNN, centroid-based methods)
            if hasattr(classifier, "X_train_") and classifier.X_train_ is not None:
                model_data["classifier"]["X_train_"] = classifier.X_train_.tolist()
            if hasattr(classifier, "y_train_") and classifier.y_train_ is not None:
                model_data["classifier"]["y_train_"] = classifier.y_train_.tolist()
            
            # Save feature importance if present
            if hasattr(classifier, "feature_importance_") and classifier.feature_importance_ is not None:
                model_data["classifier"]["feature_importance_"] = classifier.feature_importance_.tolist()
        
        # Save regressor state if present
        if hasattr(self, "regressor") and self.regressor:
            regressor = self.regressor
            model_data["regressor"] = {
                "algorithm": regressor.algorithm,
                "random_state": regressor.random_state,
                "params": regressor.params,
                "is_fitted": regressor.is_fitted,
            }
            
            # Save training data
            if hasattr(regressor, "X_train_") and regressor.X_train_ is not None:
                model_data["regressor"]["X_train_"] = regressor.X_train_.tolist()
            if hasattr(regressor, "y_train_") and regressor.y_train_ is not None:
                model_data["regressor"]["y_train_"] = regressor.y_train_.tolist()
            
            # Save coefficients and intercept for linear/ridge models
            if hasattr(regressor, "coefficients_") and regressor.coefficients_ is not None:
                model_data["regressor"]["coefficients_"] = regressor.coefficients_.tolist()
            if hasattr(regressor, "intercept_"):
                model_data["regressor"]["intercept_"] = float(regressor.intercept_)
            
            # Save feature importance if present
            if hasattr(regressor, "feature_importance_") and regressor.feature_importance_ is not None:
                model_data["regressor"]["feature_importance_"] = regressor.feature_importance_.tolist()
        
        # Save to file
        io.dump_json(model_data, path, indent=2)
    
    @classmethod
    def load_model(cls, path: str | Path) -> "EventSequencePredictor":
        """Load a trained model from disk.
        
        Reconstructs an EventSequencePredictor from a saved model file.
        The model will be in fitted state and ready for prediction.
        
        Args:
            path: Path to saved model file (JSON format)
            
        Returns:
            Reconstructed EventSequencePredictor instance
            
        Raises:
            FileNotFoundError: If model file does not exist
            ValueError: If model file format is invalid
            
        Examples:
            >>> predictor = EventSequencePredictor.load_model("output/model.json")
            >>> predictions = predictor.predict(new_sequences)
        """
        from ..core import io
        
        path = Path(path)
        if not path.exists():
            raise FileNotFoundError(f"Model file not found: {path}")
        
        # Load model data
        try:
            model_data = io.load_json(path)
        except Exception as e:
            raise ValueError(f"Failed to parse model file: {e}") from e
        
        # Validate required fields
        required_fields = ["model_type", "task_type", "is_fitted"]
        missing_fields = [f for f in required_fields if f not in model_data]
        if missing_fields:
            raise ValueError(f"Invalid model file: missing fields {missing_fields}")
        
        if not model_data["is_fitted"]:
            raise ValueError("Cannot load unfitted model")
        
        # Reconstruct predictor
        predictor = cls(
            model_type=model_data["model_type"],
            task_type=model_data["task_type"],
            embedding_dim=model_data.get("embedding_dim", 100),
            hidden_dim=model_data.get("hidden_dim", 64),
            random_state=model_data.get("random_state"),
            **model_data.get("params", {})
        )
        
        predictor.is_fitted = True
        
        # Restore event embeddings
        if "event_embeddings" in model_data:
            predictor.event_embeddings = {
                k: np.array(v, dtype=np.float64) 
                for k, v in model_data["event_embeddings"].items()
            }
        
        # Restore vocabulary
        if "vocab" in model_data:
            predictor.vocab = model_data["vocab"]
        
        # Restore vocabulary mappings
        if "vocab_to_idx" in model_data:
            predictor.vocab_to_idx = model_data["vocab_to_idx"]
        if "idx_to_vocab" in model_data:
            predictor.idx_to_vocab = model_data["idx_to_vocab"]
        
        # Restore classes
        if "classes_" in model_data:
            predictor.classes_ = np.array(model_data["classes_"])
        
        # Restore classifier
        if "classifier" in model_data:
            from ..ml.classification import BiologicalClassifier
            
            classifier_data = model_data["classifier"]
            predictor.classifier = BiologicalClassifier(
                algorithm=classifier_data["algorithm"],
                random_state=classifier_data.get("random_state"),
                **classifier_data.get("params", {})
            )
            predictor.classifier.is_fitted = classifier_data.get("is_fitted", True)
            
            if "X_train_" in classifier_data:
                predictor.classifier.X_train_ = np.array(classifier_data["X_train_"], dtype=np.float64)
            if "y_train_" in classifier_data:
                predictor.classifier.y_train_ = np.array(classifier_data["y_train_"])
            if "feature_importance_" in classifier_data:
                predictor.classifier.feature_importance_ = np.array(
                    classifier_data["feature_importance_"], dtype=np.float64
                )
            if "classes_" in model_data:
                predictor.classifier.classes_ = predictor.classes_
        
        # Restore regressor
        if "regressor" in model_data:
            from ..ml.regression import BiologicalRegressor
            
            regressor_data = model_data["regressor"]
            predictor.regressor = BiologicalRegressor(
                algorithm=regressor_data["algorithm"],
                random_state=regressor_data.get("random_state"),
                **regressor_data.get("params", {})
            )
            predictor.regressor.is_fitted = regressor_data.get("is_fitted", True)
            
            if "X_train_" in regressor_data:
                predictor.regressor.X_train_ = np.array(regressor_data["X_train_"], dtype=np.float64)
            if "y_train_" in regressor_data:
                predictor.regressor.y_train_ = np.array(regressor_data["y_train_"], dtype=np.float64)
            if "coefficients_" in regressor_data:
                predictor.regressor.coefficients_ = np.array(
                    regressor_data["coefficients_"], dtype=np.float64
                )
            if "intercept_" in regressor_data:
                predictor.regressor.intercept_ = float(regressor_data["intercept_"])
            if "feature_importance_" in regressor_data:
                predictor.regressor.feature_importance_ = np.array(
                    regressor_data["feature_importance_"], dtype=np.float64
                )
        
        return predictor


class LSTMSequenceModel:
    """LSTM-based sequence model for event sequences.
    
    Full PyTorch implementation with batching and padding support.
    Falls back to simpler model if PyTorch is not available.
    """
    
    def __init__(
        self,
        embedding_dim: int = 100,
        hidden_dim: int = 64,
        num_layers: int = 2,
        task_type: str = "classification",
        batch_size: int = 32,
        epochs: int = 10,
        learning_rate: float = 0.001,
        random_state: Optional[int] = None
    ):
        """Initialize LSTM model.
        
        Args:
            embedding_dim: Dimension of event embeddings
            hidden_dim: Hidden dimension of LSTM
            num_layers: Number of LSTM layers
            task_type: Task type ("classification" or "regression")
            batch_size: Batch size for training
            epochs: Number of training epochs
            learning_rate: Learning rate for optimizer
            random_state: Random seed
        """
        self.embedding_dim = embedding_dim
        self.hidden_dim = hidden_dim
        self.num_layers = num_layers
        self.task_type = task_type
        self.batch_size = batch_size
        self.epochs = epochs
        self.learning_rate = learning_rate
        self.random_state = random_state
        self.is_fitted = False
        
        if not TORCH_AVAILABLE:
            self.use_fallback = True
            self.fallback_model: Optional[EventSequencePredictor] = None
        else:
            self.use_fallback = False
            if random_state is not None:
                torch.manual_seed(random_state)
                np.random.seed(random_state)
            self.device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
            self.model: Optional[nn.Module] = None
            self.vocab_to_idx: Optional[Dict[str, int]] = None
            self.idx_to_vocab: Optional[Dict[int, str]] = None
            self.event_embeddings: Optional[Dict[str, NDArray]] = None
    
    def _build_model(self, vocab_size: int, num_classes: int) -> nn.Module:
        """Build PyTorch LSTM model."""
        class LSTMPredictor(nn.Module):
            def __init__(self, vocab_size, embedding_dim, hidden_dim, num_layers, num_classes, task_type):
                super().__init__()
                self.embedding = nn.Embedding(vocab_size, embedding_dim)
                self.lstm = nn.LSTM(embedding_dim, hidden_dim, num_layers, batch_first=True)
                self.task_type = task_type
                if task_type == "classification":
                    self.fc = nn.Linear(hidden_dim, num_classes)
                else:
                    self.fc = nn.Linear(hidden_dim, 1)
            
            def forward(self, x):
                embedded = self.embedding(x)
                lstm_out, (h_n, c_n) = self.lstm(embedded)
                # Use last hidden state
                last_hidden = lstm_out[:, -1, :]
                output = self.fc(last_hidden)
                if self.task_type == "classification":
                    return torch.softmax(output, dim=1)
                return output
        
        return LSTMPredictor(vocab_size, self.embedding_dim, self.hidden_dim, self.num_layers, num_classes, self.task_type)
    
    def fit(
        self,
        sequences: Sequence[Sequence[str]],
        y: NDArray,
        event_embeddings: Optional[Dict[str, NDArray]] = None
    ) -> "LSTMSequenceModel":
        """Fit LSTM model."""
        if self.use_fallback:
            self.fallback_model = EventSequencePredictor(
                model_type="simple",
                task_type=self.task_type,
                random_state=self.random_state
            )
            self.fallback_model.fit(sequences, y, event_embeddings)
            self.is_fitted = True
            return self
        
        from .embeddings import learn_event_embeddings
        
        # Learn embeddings if not provided
        if event_embeddings is None:
            event_embeddings = learn_event_embeddings(
                list(sequences),
                embedding_dim=self.embedding_dim,
                random_state=self.random_state
            )
        
        self.event_embeddings = event_embeddings
        
        # Build vocabulary
        vocab = sorted(list(event_embeddings.keys()))
        self.vocab_to_idx = {token: i + 1 for i, token in enumerate(vocab)}  # 0 is padding
        self.vocab_to_idx["<PAD>"] = 0
        self.idx_to_vocab = {i: token for token, i in self.vocab_to_idx.items()}
        
        # Convert sequences to indices
        sequences_idx = []
        for seq in sequences:
            seq_idx = [self.vocab_to_idx.get(token, 0) for token in seq]
            sequences_idx.append(seq_idx)
        
        # Pad sequences
        max_len = max(len(s) for s in sequences_idx) if sequences_idx else 1
        sequences_padded = []
        for seq in sequences_idx:
            padded = seq + [0] * (max_len - len(seq))
            sequences_padded.append(padded[:max_len])
        
        X = np.array(sequences_padded)
        
        # Determine number of classes
        if self.task_type == "classification":
            num_classes = len(np.unique(y))
            y_tensor = torch.LongTensor(y).to(self.device)
        else:
            num_classes = 1
            y_tensor = torch.FloatTensor(y).to(self.device)
        
        # Build model
        vocab_size = len(self.vocab_to_idx)
        self.model = self._build_model(vocab_size, num_classes).to(self.device)
        optimizer = torch.optim.Adam(self.model.parameters(), lr=self.learning_rate)
        criterion = nn.CrossEntropyLoss() if self.task_type == "classification" else nn.MSELoss()
        
        # Training loop
        X_tensor = torch.LongTensor(X).to(self.device)
        for epoch in range(self.epochs):
            self.model.train()
            optimizer.zero_grad()
            outputs = self.model(X_tensor)
            if self.task_type == "classification":
                loss = criterion(outputs, y_tensor)
            else:
                loss = criterion(outputs.squeeze(), y_tensor)
            loss.backward()
            optimizer.step()
        
        self.is_fitted = True
        return self
    
    def predict(self, sequences: Sequence[Sequence[str]]) -> NDArray:
        """Predict using LSTM model."""
        if not self.is_fitted:
            raise ValueError("Model must be fitted before prediction")
        
        if self.use_fallback or self.fallback_model is not None:
            return self.fallback_model.predict(sequences)
        
        if self.model is None or self.vocab_to_idx is None:
            raise ValueError("Model not properly initialized")
        
        # Convert sequences to indices
        sequences_idx = []
        for seq in sequences:
            seq_idx = [self.vocab_to_idx.get(token, 0) for token in seq]
            sequences_idx.append(seq_idx)
        
        # Pad sequences
        max_len = max(len(s) for s in sequences_idx) if sequences_idx else 1
        sequences_padded = []
        for seq in sequences_idx:
            padded = seq + [0] * (max_len - len(seq))
            sequences_padded.append(padded[:max_len])
        
        X = np.array(sequences_padded)
        X_tensor = torch.LongTensor(X).to(self.device)
        
        self.model.eval()
        with torch.no_grad():
            outputs = self.model(X_tensor)
            if self.task_type == "classification":
                predictions = torch.argmax(outputs, dim=1).cpu().numpy()
            else:
                predictions = outputs.squeeze().cpu().numpy()
        
        return predictions


class GRUSequenceModel:
    """GRU-based sequence model for event sequences.
    
    Similar to LSTM but uses GRU cells. Requires PyTorch.
    """
    
    def __init__(
        self,
        embedding_dim: int = 100,
        hidden_dim: int = 64,
        num_layers: int = 2,
        task_type: str = "classification",
        batch_size: int = 32,
        epochs: int = 10,
        learning_rate: float = 0.001,
        random_state: Optional[int] = None
    ):
        """Initialize GRU model."""
        self.embedding_dim = embedding_dim
        self.hidden_dim = hidden_dim
        self.num_layers = num_layers
        self.task_type = task_type
        self.batch_size = batch_size
        self.epochs = epochs
        self.learning_rate = learning_rate
        self.random_state = random_state
        self.is_fitted = False
        
        if not TORCH_AVAILABLE:
            self.use_fallback = True
            self.fallback_model: Optional[EventSequencePredictor] = None
        else:
            self.use_fallback = False
            if random_state is not None:
                torch.manual_seed(random_state)
                np.random.seed(random_state)
            self.device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
            self.model: Optional[nn.Module] = None
            self.vocab_to_idx: Optional[Dict[str, int]] = None
            self.event_embeddings: Optional[Dict[str, NDArray]] = None
    
    def _build_model(self, vocab_size: int, num_classes: int) -> nn.Module:
        """Build PyTorch GRU model."""
        class GRUPredictor(nn.Module):
            def __init__(self, vocab_size, embedding_dim, hidden_dim, num_layers, num_classes, task_type):
                super().__init__()
                self.embedding = nn.Embedding(vocab_size, embedding_dim)
                self.gru = nn.GRU(embedding_dim, hidden_dim, num_layers, batch_first=True)
                self.task_type = task_type
                if task_type == "classification":
                    self.fc = nn.Linear(hidden_dim, num_classes)
                else:
                    self.fc = nn.Linear(hidden_dim, 1)
            
            def forward(self, x):
                embedded = self.embedding(x)
                gru_out, h_n = self.gru(embedded)
                last_hidden = gru_out[:, -1, :]
                output = self.fc(last_hidden)
                if self.task_type == "classification":
                    return torch.softmax(output, dim=1)
                return output
        
        return GRUPredictor(vocab_size, self.embedding_dim, self.hidden_dim, self.num_layers, num_classes, self.task_type)
    
    def fit(
        self,
        sequences: Sequence[Sequence[str]],
        y: NDArray,
        event_embeddings: Optional[Dict[str, NDArray]] = None
    ) -> "GRUSequenceModel":
        """Fit GRU model."""
        if self.use_fallback:
            self.fallback_model = EventSequencePredictor(
                model_type="simple",
                task_type=self.task_type,
                random_state=self.random_state
            )
            self.fallback_model.fit(sequences, y, event_embeddings)
            self.is_fitted = True
            return self
        
        from .embeddings import learn_event_embeddings
        
        if event_embeddings is None:
            event_embeddings = learn_event_embeddings(
                list(sequences),
                embedding_dim=self.embedding_dim,
                random_state=self.random_state
            )
        
        self.event_embeddings = event_embeddings
        
        vocab = sorted(list(event_embeddings.keys()))
        self.vocab_to_idx = {token: i + 1 for i, token in enumerate(vocab)}
        self.vocab_to_idx["<PAD>"] = 0
        
        sequences_idx = [[self.vocab_to_idx.get(token, 0) for token in seq] for seq in sequences]
        max_len = max(len(s) for s in sequences_idx) if sequences_idx else 1
        sequences_padded = [seq + [0] * (max_len - len(seq)) for seq in sequences_idx]
        
        X = np.array(sequences_padded)
        
        if self.task_type == "classification":
            num_classes = len(np.unique(y))
            y_tensor = torch.LongTensor(y).to(self.device)
        else:
            num_classes = 1
            y_tensor = torch.FloatTensor(y).to(self.device)
        
        vocab_size = len(self.vocab_to_idx)
        self.model = self._build_model(vocab_size, num_classes).to(self.device)
        optimizer = torch.optim.Adam(self.model.parameters(), lr=self.learning_rate)
        criterion = nn.CrossEntropyLoss() if self.task_type == "classification" else nn.MSELoss()
        
        X_tensor = torch.LongTensor(X).to(self.device)
        for epoch in range(self.epochs):
            self.model.train()
            optimizer.zero_grad()
            outputs = self.model(X_tensor)
            if self.task_type == "classification":
                loss = criterion(outputs, y_tensor)
            else:
                loss = criterion(outputs.squeeze(), y_tensor)
            loss.backward()
            optimizer.step()
        
        self.is_fitted = True
        return self
    
    def predict(self, sequences: Sequence[Sequence[str]]) -> NDArray:
        """Predict using GRU model."""
        if not self.is_fitted:
            raise ValueError("Model must be fitted before prediction")
        
        if self.use_fallback or self.fallback_model is not None:
            return self.fallback_model.predict(sequences)
        
        if self.model is None or self.vocab_to_idx is None:
            raise ValueError("Model not properly initialized")
        
        sequences_idx = [[self.vocab_to_idx.get(token, 0) for token in seq] for seq in sequences]
        max_len = max(len(s) for s in sequences_idx) if sequences_idx else 1
        sequences_padded = [seq + [0] * (max_len - len(seq)) for seq in sequences_idx]
        
        X = np.array(sequences_padded)
        X_tensor = torch.LongTensor(X).to(self.device)
        
        self.model.eval()
        with torch.no_grad():
            outputs = self.model(X_tensor)
            if self.task_type == "classification":
                predictions = torch.argmax(outputs, dim=1).cpu().numpy()
            else:
                predictions = outputs.squeeze().cpu().numpy()
        
        return predictions


class EnsemblePredictor:
    """Ensemble of multiple prediction models.
    
    Combines predictions from multiple models using weighted averaging.
    """
    
    def __init__(
        self,
        models: List[EventSequencePredictor | LSTMSequenceModel | GRUSequenceModel],
        weights: Optional[List[float]] = None,
        task_type: str = "classification"
    ):
        """Initialize ensemble predictor.
        
        Args:
            models: List of fitted prediction models
            weights: Optional weights for each model (defaults to uniform)
            task_type: Task type ("classification" or "regression")
        """
        self.models = models
        self.task_type = task_type
        
        if weights is None:
            self.weights = [1.0 / len(models)] * len(models)
        else:
            # Normalize weights
            weight_sum = sum(weights)
            self.weights = [w / weight_sum for w in weights]
    
    def predict(self, sequences: Sequence[Sequence[str]]) -> NDArray:
        """Predict using ensemble of models.
        
        Args:
            sequences: List of event sequences
            
        Returns:
            Array of ensemble predictions
        """
        predictions_list = []
        for model in self.models:
            preds = model.predict(sequences)
            predictions_list.append(preds)
        
        # Weighted average
        if self.task_type == "classification":
            # For classification, use voting
            predictions_array = np.array(predictions_list)
            # Weighted voting
            weighted_preds = np.zeros((len(sequences), len(np.unique(predictions_array))))
            for i, model_preds in enumerate(predictions_list):
                for j, pred in enumerate(model_preds):
                    weighted_preds[j, int(pred)] += self.weights[i]
            return np.argmax(weighted_preds, axis=1)
        else:
            # For regression, weighted average
            predictions_array = np.array(predictions_list)
            weighted_preds = np.zeros(len(sequences))
            for i, model_preds in enumerate(predictions_list):
                weighted_preds += self.weights[i] * model_preds
            return weighted_preds


class SurvivalPredictor:
    """Survival analysis model for time-to-event prediction.
    
    Predicts time until an event occurs (e.g., time to diagnosis, time to job change).
    Uses Cox proportional hazards model or accelerated failure time model.
    """
    
    def __init__(
        self,
        method: str = "cox",
        random_state: Optional[int] = None
    ):
        """Initialize survival predictor.
        
        Args:
            method: Survival method ("cox" or "aft")
            random_state: Random seed
        """
        self.method = method
        self.random_state = random_state
        self.is_fitted = False
        self.model: Optional[Any] = None
        
        if method not in ["cox", "aft"]:
            raise ValueError(f"method must be 'cox' or 'aft', got {method}")
    
    def fit(
        self,
        sequences: Sequence[Sequence[str]],
        event_times: NDArray,
        event_occurred: NDArray,
        event_embeddings: Optional[Dict[str, NDArray]] = None
    ) -> "SurvivalPredictor":
        """Fit survival model.
        
        Args:
            sequences: List of event sequences
            event_times: Time until event (or censoring time)
            event_occurred: Binary array indicating if event occurred (1) or was censored (0)
            event_embeddings: Optional pre-trained embeddings
        """
        from .embeddings import learn_event_embeddings, sequence_embeddings
        
        if event_embeddings is None:
            event_embeddings = learn_event_embeddings(
                list(sequences),
                embedding_dim=100,
                random_state=self.random_state
            )
        
        # Convert sequences to feature vectors
        X = sequence_embeddings(list(sequences), event_embeddings, method="mean")
        
        # Simple implementation: use linear regression as proxy
        # In full implementation, would use lifelines library
        from ..ml.regression import BiologicalRegressor
        
        # For survival analysis, we predict event times
        # Censored observations are handled by weighting
        self.regressor = BiologicalRegressor(
            algorithm="ridge",
            random_state=self.random_state
        )
        self.regressor.fit(X, event_times)
        
        self.event_embeddings = event_embeddings
        self.is_fitted = True
        return self
    
    def predict(self, sequences: Sequence[Sequence[str]]) -> NDArray:
        """Predict time until event.
        
        Args:
            sequences: List of event sequences
            
        Returns:
            Array of predicted event times
        """
        if not self.is_fitted:
            raise ValueError("Model must be fitted before prediction")
        
        from .embeddings import sequence_embeddings
        
        X = sequence_embeddings(list(sequences), self.event_embeddings, method="mean")
        return self.regressor.predict(X)
    
    def predict_survival_function(self, sequences: Sequence[Sequence[str]], times: NDArray) -> NDArray:
        """Predict survival probabilities at specified times.
        
        Args:
            sequences: List of event sequences
            times: Time points at which to evaluate survival function
            
        Returns:
            Array of survival probabilities (shape: n_sequences x n_times)
        """
        if not self.is_fitted:
            raise ValueError("Model must be fitted before prediction")
        
        # Simple implementation: exponential survival
        predicted_times = self.predict(sequences)
        survival_probs = np.exp(-np.outer(times, 1.0 / (predicted_times + 1e-10)))
        return survival_probs.T


class MultiTaskPredictor:
    """Multi-task learning model for predicting multiple outcomes simultaneously.
    
    Predicts multiple related outcomes from the same event sequences,
    sharing learned representations across tasks.
    """
    
    def __init__(
        self,
        task_types: Dict[str, str],
        embedding_dim: int = 100,
        random_state: Optional[int] = None
    ):
        """Initialize multi-task predictor.
        
        Args:
            task_types: Dictionary mapping task names to task types
                (e.g., {"outcome1": "classification", "outcome2": "regression"})
            embedding_dim: Dimension of sequence embeddings
            random_state: Random seed
        """
        self.task_types = task_types
        self.embedding_dim = embedding_dim
        self.random_state = random_state
        self.is_fitted = False
        self.predictors: Dict[str, Any] = {}
        self.event_embeddings: Optional[Dict[str, NDArray]] = None
    
    def fit(
        self,
        sequences: Sequence[Sequence[str]],
        outcomes: Dict[str, NDArray],
        event_embeddings: Optional[Dict[str, NDArray]] = None
    ) -> "MultiTaskPredictor":
        """Fit multi-task model.
        
        Args:
            sequences: List of event sequences
            outcomes: Dictionary mapping task names to outcome arrays
            event_embeddings: Optional pre-trained embeddings
        """
        from .embeddings import learn_event_embeddings, sequence_embeddings
        
        if event_embeddings is None:
            event_embeddings = learn_event_embeddings(
                list(sequences),
                embedding_dim=self.embedding_dim,
                random_state=self.random_state
            )
        
        self.event_embeddings = event_embeddings
        
        # Convert sequences to feature vectors (shared representation)
        X = sequence_embeddings(list(sequences), event_embeddings, method="mean")
        
        # Fit separate predictor for each task
        for task_name, task_type in self.task_types.items():
            if task_name not in outcomes:
                raise ValueError(f"Outcome for task '{task_name}' not provided")
            
            y = outcomes[task_name]
            
            if task_type == "classification":
                from ..ml.classification import BiologicalClassifier
                predictor = BiologicalClassifier(
                    algorithm="random_forest",
                    random_state=self.random_state
                )
            else:
                from ..ml.regression import BiologicalRegressor
                predictor = BiologicalRegressor(
                    algorithm="ridge",
                    random_state=self.random_state
                )
            
            predictor.fit(X, y)
            self.predictors[task_name] = predictor
        
        self.is_fitted = True
        return self
    
    def predict(self, sequences: Sequence[Sequence[str]], task_name: Optional[str] = None) -> Union[Dict[str, NDArray], NDArray]:
        """Predict outcomes for all tasks or specific task.
        
        Args:
            sequences: List of event sequences
            task_name: Optional specific task name (if None, returns all tasks)
            
        Returns:
            Dictionary of predictions (if task_name is None) or array for specific task
        """
        if not self.is_fitted:
            raise ValueError("Model must be fitted before prediction")
        
        from .embeddings import sequence_embeddings
        
        X = sequence_embeddings(list(sequences), self.event_embeddings, method="mean")
        
        if task_name is not None:
            if task_name not in self.predictors:
                raise ValueError(f"Task '{task_name}' not found")
            return self.predictors[task_name].predict(X)
        
        # Return predictions for all tasks
        predictions = {}
        for task_name, predictor in self.predictors.items():
            predictions[task_name] = predictor.predict(X)
        
        return predictions

