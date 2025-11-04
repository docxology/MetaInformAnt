"""Integration tests for life_events module with other modules."""

from __future__ import annotations

from datetime import datetime
from pathlib import Path

import numpy as np
import pytest

from metainformant.life_events import (
    Event,
    EventSequence,
    EventSequencePredictor,
    LifeEventsWorkflowConfig,
    analyze_life_course,
    learn_event_embeddings,
    load_life_events_config,
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


# End-to-End Integration Tests

def test_end_to_end_save_load_predict(tmp_path: Path):
    """Test complete workflow: train -> save -> load -> predict."""
    # Train model
    train_sequences = [
        ["health:diagnosis", "occupation:job_change"],
        ["education:degree", "occupation:job_change"],
        ["health:diagnosis", "income:raise"],
    ]
    
    y_train = np.array([0, 1, 0])
    
    predictor = EventSequencePredictor(
        model_type="embedding",
        task_type="classification",
        embedding_dim=50,
        random_state=42
    )
    predictor.fit(train_sequences, y_train)
    
    # Save model
    model_file = tmp_path / "model.json"
    predictor.save_model(model_file)
    assert model_file.exists()
    
    # Load model
    loaded = EventSequencePredictor.load_model(model_file)
    assert loaded.is_fitted
    
    # Predict on new sequences
    test_sequences = [
        ["health:diagnosis", "occupation:job_change"],
        ["education:degree", "income:raise"],
    ]
    
    predictions = loaded.predict(test_sequences)
    assert len(predictions) == len(test_sequences)
    assert all(p in loaded.classes_ for p in predictions)
    
    # Verify predictions match original model
    original_preds = predictor.predict(test_sequences)
    np.testing.assert_array_equal(predictions, original_preds)


def test_config_workflow_integration(tmp_path: Path):
    """Test using config file in complete workflow."""
    # Create config file
    config_file = tmp_path / "config.yaml"
    config_content = """
work_dir: {work_dir}
threads: 4
embedding:
  embedding_dim: 50
  window_size: 3
  epochs: 5
model:
  model_type: embedding
  task_type: classification
  random_state: 42
""".format(work_dir=str(tmp_path / "work"))
    config_file.write_text(config_content)
    
    # Create sequences
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
    
    # Run workflow with config
    results = analyze_life_course(
        sequences,
        outcomes=outcomes,
        config_path=config_file,
        output_dir=tmp_path / "output"
    )
    
    assert results["n_sequences"] == 2
    assert "model" in results
    assert results["embedding_dim"] == 50
    
    # Verify model can be loaded
    model_file = Path(results["model"])
    assert model_file.exists()
    
    loaded = EventSequencePredictor.load_model(model_file)
    assert loaded.is_fitted
    assert loaded.model_type == "embedding"


def test_cli_end_to_end_workflow(tmp_path: Path):
    """Test embed -> predict -> interpret CLI workflow."""
    import subprocess
    import sys
    
    # Create sequences file
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
    
    sequences_file = tmp_path / "sequences.json"
    sequences_data = [seq.to_dict() for seq in sequences]
    from metainformant.core.io import dump_json
    dump_json(sequences_data, sequences_file)
    
    # Step 1: Train model using analyze_life_course
    outcomes = np.array([0, 1])
    results = analyze_life_course(
        sequences,
        outcomes=outcomes,
        output_dir=tmp_path / "train_output"
    )
    
    model_file = Path(results["model"])
    assert model_file.exists()
    
    # Step 2: Predict using CLI
    predict_output = tmp_path / "predict_output"
    cmd = [
        sys.executable,
        "-m",
        "metainformant",
        "life-events",
        "predict",
        f"--events={sequences_file}",
        f"--model={model_file}",
        f"--output={predict_output}",
    ]
    
    result = subprocess.run(cmd, capture_output=True, text=True, timeout=30)
    assert result.returncode == 0
    
    predictions_file = predict_output / "predictions.json"
    assert predictions_file.exists()
    
    # Step 3: Interpret using CLI
    interpret_output = tmp_path / "interpret_output"
    cmd = [
        sys.executable,
        "-m",
        "metainformant",
        "life-events",
        "interpret",
        f"--model={model_file}",
        f"--sequences={sequences_file}",
        f"--output={interpret_output}",
    ]
    
    result = subprocess.run(cmd, capture_output=True, text=True, timeout=60)
    assert result.returncode == 0
    
    report_file = interpret_output / "interpretation_report.json"
    assert report_file.exists()


def test_model_persistence_across_sessions(tmp_path: Path):
    """Verify model works after program restart (simulated)."""
    # Train and save model
    sequences = [
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
    predictor.fit(sequences, y)
    
    original_preds = predictor.predict(sequences)
    original_probas = predictor.predict_proba(sequences)
    
    model_file = tmp_path / "persistent_model.json"
    predictor.save_model(model_file)
    
    # Simulate "restart" by creating new predictor instance
    # Load model in new "session"
    loaded = EventSequencePredictor.load_model(model_file)
    
    # Verify it works identically
    loaded_preds = loaded.predict(sequences)
    loaded_probas = loaded.predict_proba(sequences)
    
    np.testing.assert_array_equal(loaded_preds, original_preds)
    np.testing.assert_allclose(loaded_probas, original_probas, rtol=1e-5)
    
    # Verify all attributes are restored
    assert loaded.model_type == predictor.model_type
    assert loaded.task_type == predictor.task_type
    assert loaded.is_fitted == predictor.is_fitted
    assert np.array_equal(loaded.classes_, predictor.classes_)

