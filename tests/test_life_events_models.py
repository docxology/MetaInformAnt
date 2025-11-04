"""Tests for life_events prediction models."""

from __future__ import annotations

import numpy as np
import pytest

from metainformant.life_events import EventSequencePredictor, LSTMSequenceModel


def test_event_sequence_predictor_embedding_classification():
    """Test EventSequencePredictor with embedding model for classification."""
    sequences = [
        ["health:diagnosis", "occupation:job_change"],
        ["education:degree", "occupation:job_change"],
        ["health:diagnosis", "income:raise"],
    ]
    
    y = np.array([0, 1, 0])  # Binary classification
    
    predictor = EventSequencePredictor(
        model_type="embedding",
        task_type="classification",
        embedding_dim=50,
        random_state=42
    )
    
    predictor.fit(sequences, y)
    
    assert predictor.is_fitted
    assert predictor.classes_ is not None
    
    # Make predictions
    predictions = predictor.predict(sequences)
    assert len(predictions) == len(sequences)
    assert all(p in predictor.classes_ for p in predictions)
    
    # Predict probabilities
    probas = predictor.predict_proba(sequences)
    assert probas.shape == (len(sequences), len(predictor.classes_))


def test_event_sequence_predictor_simple_regression():
    """Test EventSequencePredictor with simple model for regression."""
    sequences = [
        ["health:diagnosis", "occupation:job_change"],
        ["education:degree", "occupation:job_change"],
    ]
    
    y = np.array([0.5, 0.8])  # Regression targets
    
    predictor = EventSequencePredictor(
        model_type="simple",
        task_type="regression",
        random_state=42
    )
    
    predictor.fit(sequences, y)
    
    assert predictor.is_fitted
    
    predictions = predictor.predict(sequences)
    assert len(predictions) == len(sequences)
    assert all(isinstance(p, (int, float)) for p in predictions)


def test_event_sequence_predictor_with_embeddings():
    """Test EventSequencePredictor with pre-computed embeddings."""
    sequences = [
        ["health:diagnosis", "occupation:job_change"],
        ["education:degree", "occupation:job_change"],
    ]
    
    y = np.array([0, 1])
    
    # Pre-compute embeddings
    event_embeddings = {
        "health:diagnosis": np.random.randn(50),
        "occupation:job_change": np.random.randn(50),
        "education:degree": np.random.randn(50),
    }
    
    predictor = EventSequencePredictor(
        model_type="embedding",
        task_type="classification",
        embedding_dim=50,
        random_state=42
    )
    
    predictor.fit(sequences, y, event_embeddings=event_embeddings)
    
    assert predictor.is_fitted
    predictions = predictor.predict(sequences)
    assert len(predictions) == len(sequences)


def test_event_sequence_predictor_not_fitted_error():
    """Test that prediction fails if model not fitted."""
    predictor = EventSequencePredictor()
    
    sequences = [["health:diagnosis"]]
    
    with pytest.raises(ValueError, match="must be fitted"):
        predictor.predict(sequences)


def test_event_sequence_predictor_predict_proba_classification_only():
    """Test that predict_proba only works for classification."""
    sequences = [["health:diagnosis"]]
    y = np.array([0.5])
    
    predictor = EventSequencePredictor(
        task_type="regression",
        random_state=42
    )
    
    predictor.fit(sequences, y)
    
    with pytest.raises(ValueError, match="only available for classification"):
        predictor.predict_proba(sequences)


def test_lstm_sequence_model():
    """Test LSTMSequenceModel (may fallback to simple model)."""
    sequences = [
        ["health:diagnosis", "occupation:job_change"],
        ["education:degree", "occupation:job_change"],
    ]
    
    y = np.array([0, 1])
    
    model = LSTMSequenceModel(
        embedding_dim=50,
        hidden_dim=32,
        random_state=42
    )
    
    model.fit(sequences, y)
    
    assert model.is_fitted
    
    predictions = model.predict(sequences)
    assert len(predictions) == len(sequences)

