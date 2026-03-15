"""Tests for life_events interpretability functions."""

from __future__ import annotations

import numpy as np
import pytest

from metainformant.life_events.analysis.interpretability import (
    event_importance,
    feature_attribution,
    temporal_patterns,
)
from metainformant.life_events.models.embeddings import learn_event_embeddings
from metainformant.life_events.models.statistical_models import attention_weights
from metainformant.life_events.models.predictor import EventSequencePredictor


def test_event_importance_permutation():
    """Test event importance computation with permutation method."""
    sequences = [
        ["health:diagnosis", "occupation:job_change"],
        ["education:degree", "occupation:job_change"],
    ]
    y = np.array([0, 1])

    predictor = EventSequencePredictor(random_state=42)
    embeddings = learn_event_embeddings(sequences, embedding_dim=50, random_state=42)
    predictor.fit(sequences, y, event_embeddings=embeddings)

    importance = event_importance(predictor, sequences, embeddings, method="permutation")

    assert isinstance(importance, dict)
    assert len(importance) > 0
    # All scores should be normalized to [0, 1]
    assert all(0 <= score <= 1 for score in importance.values())


def test_event_importance_not_fitted():
    """Test event importance with unfitted predictor."""
    sequences = [["health:diagnosis"]]
    embeddings = {"health:diagnosis": np.random.randn(50)}

    predictor = EventSequencePredictor()

    with pytest.raises(ValueError, match="must be fitted"):
        event_importance(predictor, sequences, embeddings)


def test_event_importance_empty_sequences():
    """Test event importance with empty sequences."""
    predictor = EventSequencePredictor(random_state=42)
    sequences = [["health:diagnosis"], ["education:degree"]]
    y = np.array([0, 1])
    embeddings = learn_event_embeddings(sequences, random_state=42)
    predictor.fit(sequences, y, event_embeddings=embeddings)

    with pytest.raises(ValueError, match="cannot be empty"):
        event_importance(predictor, [], embeddings)


def test_event_importance_invalid_method():
    """Test event importance with invalid method."""
    sequences = [["health:diagnosis"], ["education:degree"]]
    y = np.array([0, 1])
    predictor = EventSequencePredictor(random_state=42)
    embeddings = learn_event_embeddings(sequences, random_state=42)
    predictor.fit(sequences, y, event_embeddings=embeddings)

    with pytest.raises(ValueError, match="Invalid method"):
        event_importance(predictor, sequences, embeddings, method="invalid")


def test_temporal_patterns():
    """Test temporal pattern analysis."""
    sequences = [
        ["health:diagnosis", "occupation:job_change"],
        ["education:degree", "occupation:job_change"],
    ]
    predictions = np.array([0.8, 0.3])

    patterns = temporal_patterns(sequences, predictions)

    assert "position_importance" in patterns
    assert "max_sequence_length" in patterns
    assert patterns["max_sequence_length"] == 2


def test_temporal_patterns_empty_sequences():
    """Test temporal patterns with empty sequences."""
    predictions = np.array([])

    with pytest.raises(ValueError, match="cannot be empty"):
        temporal_patterns([], predictions)


def test_temporal_patterns_length_mismatch():
    """Test temporal patterns with length mismatch."""
    sequences = [["health:diagnosis"]]
    predictions = np.array([0.8, 0.3])  # Wrong length

    with pytest.raises(ValueError, match="must match"):
        temporal_patterns(sequences, predictions)


def test_attention_weights():
    """Test attention weights extraction (placeholder)."""
    sequences = [
        ["health:diagnosis", "occupation:job_change"],
        ["education:degree"],
    ]
    embeddings = {"health:diagnosis": np.random.randn(50)}

    # Placeholder implementation
    attention = attention_weights(None, sequences, embeddings)

    assert attention.shape[0] == len(sequences)


def test_feature_attribution():
    """Test feature attribution computation."""
    sequences = [
        ["health:diagnosis", "occupation:job_change"],
        ["education:degree", "income:raise"],
    ]
    y = np.array([0, 1])

    predictor = EventSequencePredictor(random_state=42)
    embeddings = learn_event_embeddings(sequences, embedding_dim=50, random_state=42)
    predictor.fit(sequences, y, event_embeddings=embeddings)

    attribution = feature_attribution(predictor, sequences, embeddings, use_shap=False)

    assert isinstance(attribution, dict)
    assert len(attribution) > 0


def test_feature_attribution_not_fitted():
    """Test feature attribution with unfitted predictor."""
    sequences = [["health:diagnosis"]]
    embeddings = {"health:diagnosis": np.random.randn(50)}

    predictor = EventSequencePredictor()

    with pytest.raises(ValueError, match="must be fitted"):
        feature_attribution(predictor, sequences, embeddings)


def test_feature_attribution_empty_sequences():
    """Test feature attribution with empty sequences."""
    predictor = EventSequencePredictor(random_state=42)
    # Need at least 2 samples with 2 different classes for classification
    sequences = [["health:diagnosis"], ["education:degree"]]
    y = np.array([0, 1])
    embeddings = learn_event_embeddings(sequences, random_state=42)
    predictor.fit(sequences, y, event_embeddings=embeddings)

    with pytest.raises(ValueError, match="cannot be empty"):
        feature_attribution(predictor, [], embeddings)
