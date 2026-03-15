"""Tests for life_events prediction models."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from metainformant.life_events.models.predictor import EnsemblePredictor, EventSequencePredictor
from metainformant.life_events.models.sequence_models import GRUSequenceModel, LSTMSequenceModel
from metainformant.life_events.models.statistical_models import MultiTaskPredictor, SurvivalPredictor


def test_event_sequence_predictor_embedding_classification():
    """Test EventSequencePredictor with embedding model for classification."""
    sequences = [
        ["health:diagnosis", "occupation:job_change"],
        ["education:degree", "occupation:job_change"],
        ["health:diagnosis", "income:raise"],
    ]

    y = np.array([0, 1, 0])  # Binary classification

    predictor = EventSequencePredictor(
        model_type="embedding", task_type="classification", embedding_dim=50, random_state=42
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

    predictor = EventSequencePredictor(model_type="simple", task_type="regression", random_state=42)

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
        model_type="embedding", task_type="classification", embedding_dim=50, random_state=42
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

    predictor = EventSequencePredictor(task_type="regression", random_state=42)

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

    model = LSTMSequenceModel(embedding_dim=50, hidden_dim=32, random_state=42)

    model.fit(sequences, y)

    assert model.is_fitted

    predictions = model.predict(sequences)
    assert len(predictions) == len(sequences)


# Model Persistence Tests


def test_save_model_embedding_classification(tmp_path: Path):
    """Test saving embedding model for classification."""
    sequences = [
        ["health:diagnosis", "occupation:job_change"],
        ["education:degree", "occupation:job_change"],
        ["health:diagnosis", "income:raise"],
    ]

    y = np.array([0, 1, 0])

    predictor = EventSequencePredictor(
        model_type="embedding", task_type="classification", embedding_dim=50, random_state=42
    )

    predictor.fit(sequences, y)

    model_file = tmp_path / "model.json"
    predictor.save_model(model_file)

    assert model_file.exists()

    # Verify file is valid JSON
    from metainformant.core import io

    model_data = io.load_json(model_file)
    assert model_data["model_type"] == "embedding"
    assert model_data["task_type"] == "classification"
    assert model_data["is_fitted"] is True
    assert "event_embeddings" in model_data
    assert "classifier" in model_data


def test_save_model_embedding_regression(tmp_path: Path):
    """Test saving embedding model for regression."""
    sequences = [
        ["health:diagnosis", "occupation:job_change"],
        ["education:degree", "occupation:job_change"],
    ]

    y = np.array([0.5, 0.8])

    predictor = EventSequencePredictor(
        model_type="embedding", task_type="regression", embedding_dim=50, random_state=42
    )

    predictor.fit(sequences, y)

    model_file = tmp_path / "model.json"
    predictor.save_model(model_file)

    assert model_file.exists()

    from metainformant.core import io

    model_data = io.load_json(model_file)
    assert model_data["task_type"] == "regression"
    assert "regressor" in model_data


def test_save_model_simple_classification(tmp_path: Path):
    """Test saving simple model for classification."""
    sequences = [
        ["health:diagnosis", "occupation:job_change"],
        ["education:degree", "occupation:job_change"],
    ]

    y = np.array([0, 1])

    predictor = EventSequencePredictor(model_type="simple", task_type="classification", random_state=42)

    predictor.fit(sequences, y)

    model_file = tmp_path / "model.json"
    predictor.save_model(model_file)

    assert model_file.exists()

    from metainformant.core import io

    model_data = io.load_json(model_file)
    assert model_data["model_type"] == "simple"
    assert "vocab" in model_data
    assert "classifier" in model_data


def test_load_model_roundtrip(tmp_path: Path):
    """Test save then load, verify predictions match."""
    sequences = [
        ["health:diagnosis", "occupation:job_change"],
        ["education:degree", "occupation:job_change"],
        ["health:diagnosis", "income:raise"],
    ]

    y = np.array([0, 1, 0])

    # Train and save original model
    original = EventSequencePredictor(
        model_type="embedding", task_type="classification", embedding_dim=50, random_state=42
    )
    original.fit(sequences, y)

    original_predictions = original.predict(sequences)
    original_probas = original.predict_proba(sequences)

    model_file = tmp_path / "model.json"
    original.save_model(model_file)

    # Load model
    loaded = EventSequencePredictor.load_model(model_file)

    assert loaded.is_fitted
    assert loaded.model_type == original.model_type
    assert loaded.task_type == original.task_type

    # Verify predictions match
    loaded_predictions = loaded.predict(sequences)
    np.testing.assert_array_equal(loaded_predictions, original_predictions)

    # Verify probabilities match
    loaded_probas = loaded.predict_proba(sequences)
    np.testing.assert_allclose(loaded_probas, original_probas, rtol=1e-5)


def test_load_model_invalid_file(tmp_path: Path):
    """Test handling of corrupted/missing model files."""
    # Test missing file
    missing_file = tmp_path / "nonexistent.json"
    with pytest.raises(FileNotFoundError):
        EventSequencePredictor.load_model(missing_file)

    # Test corrupted JSON
    corrupted_file = tmp_path / "corrupted.json"
    corrupted_file.write_text("{invalid json")
    with pytest.raises(ValueError, match="Failed to parse"):
        EventSequencePredictor.load_model(corrupted_file)


def test_load_model_missing_fields(tmp_path: Path):
    """Test handling of model files with missing required fields."""
    from metainformant.core import io

    # Missing model_type
    invalid_file = tmp_path / "invalid.json"
    io.dump_json({"task_type": "classification", "is_fitted": True}, invalid_file)
    with pytest.raises(ValueError, match="missing fields"):
        EventSequencePredictor.load_model(invalid_file)

    # Missing is_fitted
    invalid_file2 = tmp_path / "invalid2.json"
    io.dump_json({"model_type": "embedding", "task_type": "classification"}, invalid_file2)
    with pytest.raises(ValueError, match="missing fields"):
        EventSequencePredictor.load_model(invalid_file2)


def test_load_model_unfitted_error(tmp_path: Path):
    """Test rejection of loading unfitted models."""
    from metainformant.core import io

    invalid_file = tmp_path / "unfitted.json"
    io.dump_json({"model_type": "embedding", "task_type": "classification", "is_fitted": False}, invalid_file)

    with pytest.raises(ValueError, match="Cannot load unfitted model"):
        EventSequencePredictor.load_model(invalid_file)


def test_save_unfitted_model_error(tmp_path: Path):
    """Test prevention of saving unfitted models."""
    predictor = EventSequencePredictor()
    model_file = tmp_path / "model.json"

    with pytest.raises(ValueError, match="must be fitted before saving"):
        predictor.save_model(model_file)


def test_model_persistence_predictions_match(tmp_path: Path):
    """Verify loaded model makes same predictions as original for new sequences."""
    train_sequences = [
        ["health:diagnosis", "occupation:job_change"],
        ["education:degree", "occupation:job_change"],
        ["health:diagnosis", "income:raise"],
    ]

    y = np.array([0, 1, 0])

    # Train and save
    predictor = EventSequencePredictor(
        model_type="embedding", task_type="classification", embedding_dim=50, random_state=42
    )
    predictor.fit(train_sequences, y)

    model_file = tmp_path / "model.json"
    predictor.save_model(model_file)

    # Test on new sequences
    test_sequences = [
        ["health:diagnosis", "occupation:job_change"],
        ["education:degree", "income:raise"],
    ]

    original_preds = predictor.predict(test_sequences)
    original_probas = predictor.predict_proba(test_sequences)

    # Load and test
    loaded = EventSequencePredictor.load_model(model_file)
    loaded_preds = loaded.predict(test_sequences)
    loaded_probas = loaded.predict_proba(test_sequences)

    np.testing.assert_array_equal(loaded_preds, original_preds)
    np.testing.assert_allclose(loaded_probas, original_probas, rtol=1e-5)


def test_save_load_regression_model(tmp_path: Path):
    """Test save/load for regression model."""
    sequences = [
        ["health:diagnosis", "occupation:job_change"],
        ["education:degree", "occupation:job_change"],
    ]

    y = np.array([0.5, 0.8])

    predictor = EventSequencePredictor(
        model_type="embedding", task_type="regression", embedding_dim=50, random_state=42
    )
    predictor.fit(sequences, y)

    original_preds = predictor.predict(sequences)

    model_file = tmp_path / "model.json"
    predictor.save_model(model_file)

    loaded = EventSequencePredictor.load_model(model_file)
    loaded_preds = loaded.predict(sequences)

    np.testing.assert_allclose(loaded_preds, original_preds, rtol=1e-5)


def test_gru_sequence_model(tmp_path):
    """Test GRU sequence model."""
    sequences = [
        ["health:diagnosis", "occupation:job_change"],
        ["education:degree", "occupation:job_change"],
    ]
    y = np.array([0, 1])

    model = GRUSequenceModel(embedding_dim=50, hidden_dim=32, epochs=2, random_state=42)

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

    ensemble = EnsemblePredictor(models=[model1, model2], task_type="classification")

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
        task_types={"health_outcome": "classification", "income_level": "regression"}, random_state=42
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
