"""Integration tests for life_events module with other modules."""

from __future__ import annotations

from datetime import datetime

import numpy as np
import pytest

from metainformant.life_events import (
    Event,
    EventSequence,
    EventSequencePredictor,
    analyze_life_course,
    learn_event_embeddings,
    sequence_embeddings,
)


def test_integration_with_ml_module():
    """Test integration with ML module."""
    from metainformant.ml import BiologicalClassifier

    sequences = [
        ["health:diagnosis", "occupation:job_change"],
        ["education:degree", "occupation:job_change"],
        ["health:diagnosis", "income:raise"],
    ]

    # Learn embeddings
    embeddings = learn_event_embeddings(sequences, embedding_dim=50, random_state=42)

    # Convert to sequence embeddings
    X = sequence_embeddings(sequences, embeddings, method="mean")
    y = np.array([0, 1, 0])

    # Train classifier
    classifier = BiologicalClassifier(algorithm="random_forest", random_state=42)
    classifier.fit(X, y)

    # Make predictions
    predictions = classifier.predict(X)
    assert len(predictions) == len(sequences)


def test_integration_with_phenotype_module():
    """Test integration with phenotype module."""
    try:
        from metainformant.phenotype import extract_phenotypes_from_events
    except ImportError:
        pytest.skip("Phenotype integration not available")

    sequence = EventSequence(
        person_id="person_001",
        events=[
            Event("diabetes", datetime(2020, 1, 1), "health", {"severity": "moderate"}),
            Event("bachelors", datetime(2010, 6, 1), "education", {"degree": "BS"}),
        ]
    )

    phenotypes = extract_phenotypes_from_events(sequence)

    assert "health_events" in phenotypes
    assert "education_events" in phenotypes
    assert phenotypes["health_events"] == 1
    assert phenotypes["education_events"] == 1


def test_integration_with_visualization_module():
    """Test integration with visualization module."""
    try:
        from metainformant.life_events import learn_event_embeddings, plot_event_embeddings
        from metainformant.ml.dimensionality import biological_embedding
    except ImportError:
        pytest.skip("Visualization integration not available")

    sequences = [
        ["health:diagnosis", "occupation:job_change"],
        ["education:degree", "occupation:job_change"],
    ]

    embeddings = learn_event_embeddings(sequences, embedding_dim=50, random_state=42)

    # Convert to matrix for dimensionality reduction
    embedding_matrix = np.array(list(embeddings.values()))

    # Reduce dimensions
    reduced = biological_embedding(embedding_matrix, method="pca", n_components=2)

    assert "embedding" in reduced
    assert reduced["embedding"].shape[1] == 2


def test_integration_workflow_with_ml(tmp_path):
    """Test complete workflow with ML integration."""
    sequences = [
        EventSequence(
            person_id="person_001",
            events=[Event("degree", datetime(2010, 6, 1), "education")]
        ),
        EventSequence(
            person_id="person_002",
            events=[Event("diagnosis", datetime(2020, 1, 15), "health")]
        ),
    ]

    outcomes = np.array([0, 1])

    results = analyze_life_course(
        sequences,
        outcomes=outcomes,
        output_dir=tmp_path / "output"
    )

    assert "model_type" in results
    assert "predictions" in results
    assert len(results["predictions"]) == len(sequences)


def test_integration_embedding_to_classifier():
    """Test embedding to classifier pipeline."""
    from metainformant.ml import BiologicalClassifier

    sequences = [
        ["health:diagnosis", "occupation:job_change"],
        ["education:degree", "occupation:job_change"],
    ]

    # Learn embeddings
    embeddings = learn_event_embeddings(sequences, embedding_dim=50, random_state=42)

    # Get sequence embeddings
    X = sequence_embeddings(sequences, embeddings)

    # Train classifier
    y = np.array([0, 1])
    classifier = BiologicalClassifier(algorithm="random_forest", random_state=42)
    classifier.fit(X, y)

    # Verify predictions
    predictions = classifier.predict(X)
    assert len(predictions) == len(sequences)
    assert all(p in classifier.classes_ for p in predictions)

