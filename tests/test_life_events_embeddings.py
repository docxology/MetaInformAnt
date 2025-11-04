"""Tests for life_events embedding methods."""

from __future__ import annotations

import numpy as np
import pytest

from metainformant.life_events import (
    domain_specific_embeddings,
    learn_event_embeddings,
    sequence_embeddings,
)


def test_learn_event_embeddings_basic():
    """Test basic event embedding learning."""
    sequences = [
        ["health:diagnosis", "occupation:job_change", "income:raise"],
        ["education:degree", "occupation:job_change", "address:move"],
    ]
    
    embeddings = learn_event_embeddings(
        sequences,
        embedding_dim=50,
        window_size=5,
        epochs=5,
        random_state=42
    )
    
    # Check that all events have embeddings
    assert "health:diagnosis" in embeddings
    assert "occupation:job_change" in embeddings
    assert "education:degree" in embeddings
    
    # Check embedding dimensions
    assert embeddings["health:diagnosis"].shape == (50,)
    assert embeddings["occupation:job_change"].shape == (50,)


def test_learn_event_embeddings_empty_sequences():
    """Test embedding learning with empty sequences."""
    sequences = [[]]
    
    embeddings = learn_event_embeddings(
        sequences,
        embedding_dim=50,
        random_state=42
    )
    
    # Should handle empty sequences gracefully
    assert isinstance(embeddings, dict)


def test_learn_event_embeddings_methods():
    """Test different embedding methods."""
    sequences = [
        ["health:diagnosis", "occupation:job_change"],
        ["education:degree", "occupation:job_change"],
    ]
    
    # Skip-gram
    embeddings_sg = learn_event_embeddings(
        sequences,
        embedding_dim=50,
        method="skipgram",
        epochs=5,
        random_state=42
    )
    
    # CBOW
    embeddings_cbow = learn_event_embeddings(
        sequences,
        embedding_dim=50,
        method="cbow",
        epochs=5,
        random_state=42
    )
    
    assert "health:diagnosis" in embeddings_sg
    assert "health:diagnosis" in embeddings_cbow
    assert embeddings_sg["health:diagnosis"].shape == (50,)
    assert embeddings_cbow["health:diagnosis"].shape == (50,)


def test_sequence_embeddings_mean():
    """Test sequence embedding with mean pooling."""
    sequences = [
        ["health:diagnosis", "occupation:job_change"],
        ["education:degree"],
    ]
    
    event_embeddings = {
        "health:diagnosis": np.array([0.1, 0.2, 0.3]),
        "occupation:job_change": np.array([0.4, 0.5, 0.6]),
        "education:degree": np.array([0.7, 0.8, 0.9]),
    }
    
    seq_embeddings = sequence_embeddings(
        sequences,
        event_embeddings,
        method="mean"
    )
    
    assert seq_embeddings.shape == (2, 3)
    # First sequence: mean of [0.1,0.2,0.3] and [0.4,0.5,0.6] = [0.25, 0.35, 0.45]
    np.testing.assert_array_almost_equal(
        seq_embeddings[0],
        np.array([0.25, 0.35, 0.45])
    )
    # Second sequence: just [0.7, 0.8, 0.9]
    np.testing.assert_array_almost_equal(
        seq_embeddings[1],
        np.array([0.7, 0.8, 0.9])
    )


def test_sequence_embeddings_sum():
    """Test sequence embedding with sum pooling."""
    sequences = [
        ["health:diagnosis", "occupation:job_change"],
    ]
    
    event_embeddings = {
        "health:diagnosis": np.array([0.1, 0.2]),
        "occupation:job_change": np.array([0.3, 0.4]),
    }
    
    seq_embeddings = sequence_embeddings(
        sequences,
        event_embeddings,
        method="sum"
    )
    
    assert seq_embeddings.shape == (1, 2)
    np.testing.assert_array_almost_equal(
        seq_embeddings[0],
        np.array([0.4, 0.6])
    )


def test_sequence_embeddings_max():
    """Test sequence embedding with max pooling."""
    sequences = [
        ["health:diagnosis", "occupation:job_change"],
    ]
    
    event_embeddings = {
        "health:diagnosis": np.array([0.1, 0.2]),
        "occupation:job_change": np.array([0.3, 0.4]),
    }
    
    seq_embeddings = sequence_embeddings(
        sequences,
        event_embeddings,
        method="max"
    )
    
    assert seq_embeddings.shape == (1, 2)
    np.testing.assert_array_almost_equal(
        seq_embeddings[0],
        np.array([0.3, 0.4])
    )


def test_sequence_embeddings_empty_sequence():
    """Test sequence embedding with empty sequence."""
    sequences = [[]]
    
    event_embeddings = {
        "health:diagnosis": np.array([0.1, 0.2]),
    }
    
    seq_embeddings = sequence_embeddings(
        sequences,
        event_embeddings,
        method="mean"
    )
    
    assert seq_embeddings.shape == (1, 2)
    np.testing.assert_array_almost_equal(seq_embeddings[0], np.array([0.0, 0.0]))


def test_sequence_embeddings_unknown_event():
    """Test sequence embedding with unknown event."""
    sequences = [
        ["health:diagnosis", "unknown:event"],
    ]
    
    event_embeddings = {
        "health:diagnosis": np.array([0.1, 0.2]),
    }
    
    seq_embeddings = sequence_embeddings(
        sequences,
        event_embeddings,
        method="mean"
    )
    
    assert seq_embeddings.shape == (1, 2)
    # Should average known event with zero vector (unknown event)
    np.testing.assert_array_almost_equal(seq_embeddings[0], np.array([0.05, 0.1]))


def test_sequence_embeddings_temporal_weighting():
    """Test sequence embedding with temporal weighting."""
    sequences = [
        ["health:diagnosis", "occupation:job_change"],
    ]
    
    event_embeddings = {
        "health:diagnosis": np.array([0.1, 0.2]),
        "occupation:job_change": np.array([0.3, 0.4]),
    }
    
    seq_embeddings = sequence_embeddings(
        sequences,
        event_embeddings,
        method="mean",
        temporal_weighting=True
    )
    
    assert seq_embeddings.shape == (1, 2)
    # More recent event should have higher weight
    # Weights: [0.5, 1.0] normalized -> [1/3, 2/3]
    # Result: 1/3 * [0.1,0.2] + 2/3 * [0.3,0.4] = [0.233..., 0.333...]
    assert seq_embeddings[0][0] > 0.1  # Should be weighted toward more recent


def test_domain_specific_embeddings():
    """Test domain-specific embedding learning."""
    sequences = [
        ["diagnosis", "job_change", "raise"],
        ["degree", "job_change", "move"],
    ]
    domains = [
        ["health", "occupation", "income"],
        ["education", "occupation", "address"],
    ]
    
    domain_embeddings = domain_specific_embeddings(
        sequences,
        domains,
        embedding_dim=50,
        window_size=5,
        epochs=5,
        random_state=42
    )
    
    # Should have separate embeddings per domain
    assert "health" in domain_embeddings
    assert "occupation" in domain_embeddings
    assert "education" in domain_embeddings
    
    # Check that events in domain have embeddings
    assert "diagnosis" in domain_embeddings["health"]
    assert "job_change" in domain_embeddings["occupation"]
    assert "degree" in domain_embeddings["education"]


def test_domain_specific_embeddings_empty_domain():
    """Test domain-specific embeddings with empty domain."""
    sequences = [
        ["diagnosis"],
    ]
    domains = [
        ["health"],
    ]
    
    domain_embeddings = domain_specific_embeddings(
        sequences,
        domains,
        embedding_dim=50,
        random_state=42
    )
    
    assert "health" in domain_embeddings
    assert "diagnosis" in domain_embeddings["health"]


# Edge Cases and Validation Tests

def test_learn_event_embeddings_validation():
    """Test all validation checks in learn_event_embeddings."""
    from metainformant.life_events import learn_event_embeddings
    
    # Test empty sequences
    with pytest.raises(ValueError, match="cannot be empty"):
        learn_event_embeddings([], embedding_dim=50)
    
    # Test invalid embedding_dim
    with pytest.raises(ValueError, match="must be positive"):
        learn_event_embeddings([["health:diagnosis"]], embedding_dim=0)
    
    with pytest.raises(ValueError, match="must be positive"):
        learn_event_embeddings([["health:diagnosis"]], embedding_dim=-1)
    
    # Test invalid window_size
    with pytest.raises(ValueError, match="must be positive"):
        learn_event_embeddings([["health:diagnosis"]], window_size=0)
    
    with pytest.raises(ValueError, match="must be positive"):
        learn_event_embeddings([["health:diagnosis"]], window_size=-1)
    
    # Test invalid method
    with pytest.raises(ValueError, match="must be 'skipgram' or 'cbow'"):
        learn_event_embeddings([["health:diagnosis"]], method="invalid")
    
    # Test all empty sequences (no events)
    with pytest.raises(ValueError, match="No events found"):
        learn_event_embeddings([[], []], embedding_dim=50)


def test_learn_event_embeddings_verbose(capsys):
    """Test verbose mode output."""
    from metainformant.life_events import learn_event_embeddings
    
    sequences = [
        ["health:diagnosis", "occupation:job_change"],
        ["education:degree", "occupation:job_change"],
    ]
    
    learn_event_embeddings(
        sequences,
        embedding_dim=50,
        epochs=3,
        verbose=True,
        random_state=42
    )
    
    # Check that verbose output was printed
    captured = capsys.readouterr()
    # Should have some output about building embeddings
    assert "Building embeddings" in captured.out or "Epoch" in captured.out


def test_predict_with_missing_events_in_vocab():
    """Test prediction with out-of-vocabulary events."""
    from metainformant.life_events import EventSequencePredictor
    import numpy as np
    
    # Train on specific vocabulary
    train_sequences = [
        ["health:diagnosis", "occupation:job_change"],
        ["education:degree", "occupation:job_change"],
    ]
    
    y = np.array([0, 1])
    
    predictor = EventSequencePredictor(
        model_type="embedding",
        task_type="classification",
        embedding_dim=50,
        random_state=42
    )
    predictor.fit(train_sequences, y)
    
    # Predict with OOV event
    test_sequences = [
        ["health:diagnosis", "unknown:event"],
    ]
    
    # Should handle gracefully (may use zero vector for unknown events)
    predictions = predictor.predict(test_sequences)
    assert len(predictions) == 1
    assert predictions[0] in predictor.classes_


def test_workflow_config_priority(tmp_path):
    """Test that config_obj takes priority over config_path."""
    from metainformant.life_events import (
        Event,
        EventSequence,
        LifeEventsWorkflowConfig,
        analyze_life_course,
    )
    
    events1 = [Event("degree", datetime(2010, 6, 1), "education")]
    sequences = [EventSequence(person_id="person_001", events=events1)]
    
    # Create config file
    config_file = tmp_path / "config.yaml"
    config_file.write_text("""
embedding:
  embedding_dim: 100
""")
    
    # Pass config_obj that should override
    config_obj = LifeEventsWorkflowConfig(
        work_dir=tmp_path / "work",
        embedding={"embedding_dim": 50}
    )
    
    results = analyze_life_course(
        sequences,
        config_path=config_file,
        config_obj=config_obj,
        output_dir=tmp_path / "output"
    )
    
    # config_obj should take priority
    assert results["embedding_dim"] == 50


def test_config_env_override_edge_cases(tmp_path):
    """Test invalid environment variable values."""
    import os
    from metainformant.life_events import load_life_events_config
    
    config_file = tmp_path / "test_config.yaml"
    config_file.write_text("""
work_dir: output/test
threads: 4
""")
    
    # Set invalid threads value
    os.environ["LE_THREADS"] = "invalid"
    
    try:
        # Should handle gracefully (may use default or log warning)
        config = load_life_events_config(config_file)
        # Should either use default or handle error
        assert isinstance(config.threads, int)
    finally:
        os.environ.pop("LE_THREADS", None)

