"""Tests for advanced prediction models in life_events module."""

from __future__ import annotations

import numpy as np
import pytest

from metainformant.life_events import (
    EnsemblePredictor,
    EventSequencePredictor,
    GRUSequenceModel,
    LSTMSequenceModel,
    MultiTaskPredictor,
    SurvivalPredictor,
    learn_event_embeddings,
)


def test_lstm_sequence_model(tmp_path):
    """Test LSTM sequence model."""
    sequences = [
        ["health:diagnosis", "occupation:job_change"],
        ["education:degree", "occupation:job_change"],
        ["health:diagnosis", "income:raise"],
    ]
    y = np.array([0, 1, 0])
    
    model = LSTMSequenceModel(
        embedding_dim=50,
        hidden_dim=32,
        num_layers=1,
        task_type="classification",
        epochs=2,
        random_state=42
    )
    
    model.fit(sequences, y)
    predictions = model.predict(sequences)
    
    assert len(predictions) == len(sequences)
    assert all(isinstance(p, (int, np.integer)) for p in predictions)


def test_gru_sequence_model(tmp_path):
    """Test GRU sequence model."""
    sequences = [
        ["health:diagnosis", "occupation:job_change"],
        ["education:degree", "occupation:job_change"],
    ]
    y = np.array([0, 1])
    
    model = GRUSequenceModel(
        embedding_dim=50,
        hidden_dim=32,
        epochs=2,
        random_state=42
    )
    
    model.fit(sequences, y)
    predictions = model.predict(sequences)
    
    assert len(predictions) == len(sequences)


def test_ensemble_predictor(tmp_path):
    """Test ensemble predictor."""
    sequences = [
        ["health:diagnosis", "occupation:job_change"],
        ["education:degree", "occupation:job_change"],
        ["health:diagnosis", "income:raise"],
    ]
    y = np.array([0, 1, 0])
    
    # Train multiple models
    model1 = EventSequencePredictor(model_type="embedding", random_state=42)
    model1.fit(sequences, y)
    
    model2 = EventSequencePredictor(model_type="simple", random_state=43)
    model2.fit(sequences, y)
    
    ensemble = EnsemblePredictor(
        models=[model1, model2],
        task_type="classification"
    )
    
    predictions = ensemble.predict(sequences)
    assert len(predictions) == len(sequences)


def test_survival_predictor(tmp_path):
    """Test survival predictor."""
    sequences = [
        ["health:diagnosis", "occupation:job_change"],
        ["education:degree", "occupation:job_change"],
        ["health:diagnosis", "income:raise"],
    ]
    event_times = np.array([365, 730, 1095])
    event_occurred = np.array([1, 1, 0])
    
    model = SurvivalPredictor(method="cox", random_state=42)
    model.fit(sequences, event_times, event_occurred)
    
    predicted_times = model.predict(sequences)
    assert len(predicted_times) == len(sequences)
    assert all(predicted_times > 0)
    
    # Test survival function
    times = np.array([365, 730])
    survival_probs = model.predict_survival_function(sequences, times)
    assert survival_probs.shape == (len(sequences), len(times))


def test_multi_task_predictor(tmp_path):
    """Test multi-task predictor."""
    sequences = [
        ["health:diagnosis", "occupation:job_change"],
        ["education:degree", "occupation:job_change"],
        ["health:diagnosis", "income:raise"],
    ]
    
    outcomes = {
        "health_outcome": np.array([0, 1, 0]),
        "income_level": np.array([50000, 75000, 60000]),
    }
    
    model = MultiTaskPredictor(
        task_types={
            "health_outcome": "classification",
            "income_level": "regression"
        },
        random_state=42
    )
    
    model.fit(sequences, outcomes)
    
    # Predict all tasks
    all_predictions = model.predict(sequences)
    assert isinstance(all_predictions, dict)
    assert "health_outcome" in all_predictions
    assert "income_level" in all_predictions
    
    # Predict specific task
    health_predictions = model.predict(sequences, task_name="health_outcome")
    assert len(health_predictions) == len(sequences)

