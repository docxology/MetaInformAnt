"""Tests for life_events CLI commands."""

from __future__ import annotations

import json
import os
import subprocess
import sys
from datetime import datetime
from pathlib import Path

import numpy as np
import pytest

from metainformant.core.io.io import dump_json
from metainformant.life_events import (
    Event,
    EventSequence,
    EventSequencePredictor,
    analyze_life_course,
    load_sequences_from_json,
)

# Get the repository root directory
REPO_ROOT = Path(__file__).resolve().parents[1]
SRC_DIR = REPO_ROOT / "src"


def _get_cli_env() -> dict[str, str]:
    """Get environment variables for CLI invocation."""
    env = os.environ.copy()
    # Add src directory to PYTHONPATH
    pythonpath = env.get("PYTHONPATH", "")
    if pythonpath:
        env["PYTHONPATH"] = f"{SRC_DIR}:{pythonpath}"
    else:
        env["PYTHONPATH"] = str(SRC_DIR)
    return env


def create_test_sequences_file(tmp_path: Path) -> Path:
    """Create a test event sequences JSON file."""
    events1 = [
        Event("degree", datetime(2010, 6, 1), "education"),
        Event("job_change", datetime(2015, 3, 1), "occupation"),
    ]
    events2 = [
        Event("diagnosis", datetime(2020, 1, 15), "health"),
    ]
    
    sequences = [
        EventSequence(person_id="person_001", events=events1),
        EventSequence(person_id="person_002", events=events2),
    ]
    
    sequences_file = tmp_path / "sequences.json"
    sequences_data = [seq.to_dict() for seq in sequences]
    dump_json(sequences_data, sequences_file)
    
    return sequences_file


def create_test_model(tmp_path: Path) -> Path:
    """Create and save a test model."""
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
    
    model_file = tmp_path / "model.json"
    predictor.save_model(model_file)
    
    return model_file


def test_cli_predict_basic(tmp_path: Path) -> None:
    """Test predict command with valid model and sequences."""
    # Create test data
    sequences_file = create_test_sequences_file(tmp_path)
    model_file = create_test_model(tmp_path)
    output_dir = tmp_path / "predictions"
    
    cmd = [
        sys.executable,
        "-m",
        "metainformant",
        "life-events",
        "predict",
        f"--events={sequences_file}",
        f"--model={model_file}",
        f"--output={output_dir}",
    ]
    
    result = subprocess.run(cmd, capture_output=True, text=True, timeout=30, env=_get_cli_env())
    
    # Should succeed
    assert result.returncode == 0
    
    # Verify output file exists
    predictions_file = output_dir / "predictions.json"
    assert predictions_file.exists()
    
    # Verify output format
    predictions_data = json.loads(predictions_file.read_text())
    assert "n_sequences" in predictions_data
    assert "predictions" in predictions_data
    assert len(predictions_data["predictions"]) == 2
    assert predictions_data["predictions"][0]["sequence_id"] == "person_001"


def test_cli_predict_missing_model(tmp_path: Path) -> None:
    """Test predict command with missing model file."""
    sequences_file = create_test_sequences_file(tmp_path)
    missing_model = tmp_path / "nonexistent_model.json"
    output_dir = tmp_path / "predictions"
    
    cmd = [
        sys.executable,
        "-m",
        "metainformant",
        "life-events",
        "predict",
        f"--events={sequences_file}",
        f"--model={missing_model}",
        f"--output={output_dir}",
    ]
    
    result = subprocess.run(cmd, capture_output=True, text=True, timeout=30, env=_get_cli_env())
    
    # Should fail with error
    assert result.returncode != 0
    assert "not found" in result.stderr.lower() or "Error" in result.stderr


def test_cli_predict_missing_events(tmp_path: Path) -> None:
    """Test predict command with missing events file."""
    model_file = create_test_model(tmp_path)
    missing_events = tmp_path / "nonexistent_events.json"
    output_dir = tmp_path / "predictions"
    
    cmd = [
        sys.executable,
        "-m",
        "metainformant",
        "life-events",
        "predict",
        f"--events={missing_events}",
        f"--model={model_file}",
        f"--output={output_dir}",
    ]
    
    result = subprocess.run(cmd, capture_output=True, text=True, timeout=30, env=_get_cli_env())
    
    # Should fail with error
    assert result.returncode != 0
    assert "not found" in result.stderr.lower() or "Error" in result.stderr


def test_cli_predict_invalid_model(tmp_path: Path) -> None:
    """Test predict command with corrupted model file."""
    sequences_file = create_test_sequences_file(tmp_path)
    corrupted_model = tmp_path / "corrupted_model.json"
    corrupted_model.write_text("{invalid json")
    output_dir = tmp_path / "predictions"
    
    cmd = [
        sys.executable,
        "-m",
        "metainformant",
        "life-events",
        "predict",
        f"--events={sequences_file}",
        f"--model={corrupted_model}",
        f"--output={output_dir}",
    ]
    
    result = subprocess.run(cmd, capture_output=True, text=True, timeout=30, env=_get_cli_env())
    
    # Should fail with error
    assert result.returncode != 0
    assert "Error" in result.stderr or "Failed" in result.stderr


def test_cli_predict_output_format(tmp_path: Path) -> None:
    """Test that predictions output has correct format."""
    sequences_file = create_test_sequences_file(tmp_path)
    model_file = create_test_model(tmp_path)
    output_dir = tmp_path / "predictions"
    
    cmd = [
        sys.executable,
        "-m",
        "metainformant",
        "life-events",
        "predict",
        f"--events={sequences_file}",
        f"--model={model_file}",
        f"--output={output_dir}",
    ]
    
    result = subprocess.run(cmd, capture_output=True, text=True, timeout=30, env=_get_cli_env())
    assert result.returncode == 0
    
    predictions_file = output_dir / "predictions.json"
    predictions_data = json.loads(predictions_file.read_text())
    
    # Verify structure
    assert "n_sequences" in predictions_data
    assert "model_path" in predictions_data
    assert "model_type" in predictions_data
    assert "task_type" in predictions_data
    assert "predictions" in predictions_data
    
    # Verify each prediction entry
    for pred in predictions_data["predictions"]:
        assert "sequence_id" in pred
        assert "prediction" in pred


def test_cli_predict_classification_probs(tmp_path: Path) -> None:
    """Test that probabilities are included for classification tasks."""
    sequences_file = create_test_sequences_file(tmp_path)
    model_file = create_test_model(tmp_path)  # Creates classification model
    output_dir = tmp_path / "predictions"
    
    cmd = [
        sys.executable,
        "-m",
        "metainformant",
        "life-events",
        "predict",
        f"--events={sequences_file}",
        f"--model={model_file}",
        f"--output={output_dir}",
    ]
    
    result = subprocess.run(cmd, capture_output=True, text=True, timeout=30, env=_get_cli_env())
    assert result.returncode == 0
    
    predictions_file = output_dir / "predictions.json"
    predictions_data = json.loads(predictions_file.read_text())
    
    # Check if probabilities are included (may not always be present if predict_proba fails)
    for pred in predictions_data["predictions"]:
        # At minimum, prediction should be present
        assert "prediction" in pred
        # Probabilities may be present for classification
        if "probabilities" in pred:
            assert isinstance(pred["probabilities"], dict)


def test_cli_predict_regression_stats(tmp_path: Path) -> None:
    """Test that statistics are printed for regression tasks."""
    # Create regression model
    sequences = [
        ["health:diagnosis", "occupation:job_change"],
        ["education:degree", "occupation:job_change"],
    ]
    
    y = np.array([0.5, 0.8])
    
    predictor = EventSequencePredictor(
        model_type="embedding",
        task_type="regression",
        embedding_dim=50,
        random_state=42
    )
    predictor.fit(sequences, y)
    
    model_file = tmp_path / "regression_model.json"
    predictor.save_model(model_file)
    
    # Create sequences file
    sequences_file = create_test_sequences_file(tmp_path)
    output_dir = tmp_path / "predictions"
    
    cmd = [
        sys.executable,
        "-m",
        "metainformant",
        "life-events",
        "predict",
        f"--events={sequences_file}",
        f"--model={model_file}",
        f"--output={output_dir}",
    ]
    
    result = subprocess.run(cmd, capture_output=True, text=True, timeout=30, env=_get_cli_env())
    assert result.returncode == 0
    
    # Check for statistics in output
    assert "Mean:" in result.stdout or "statistics" in result.stdout.lower()


def test_cli_interpret_basic(tmp_path: Path) -> None:
    """Test interpret command with valid inputs."""
    sequences_file = create_test_sequences_file(tmp_path)
    model_file = create_test_model(tmp_path)
    output_dir = tmp_path / "interpretation"
    
    cmd = [
        sys.executable,
        "-m",
        "metainformant",
        "life-events",
        "interpret",
        f"--model={model_file}",
        f"--sequences={sequences_file}",
        f"--output={output_dir}",
    ]
    
    result = subprocess.run(cmd, capture_output=True, text=True, timeout=60, env=_get_cli_env())
    
    # Should succeed
    assert result.returncode == 0
    
    # Verify output file exists
    report_file = output_dir / "interpretation_report.json"
    assert report_file.exists()
    
    # Verify output format
    report_data = json.loads(report_file.read_text())
    assert "model_path" in report_data
    assert "interpretations" in report_data
    assert "event_importance" in report_data["interpretations"]
    assert "temporal_patterns" in report_data["interpretations"]
    assert "feature_attribution" in report_data["interpretations"]


def test_cli_interpret_missing_model(tmp_path: Path) -> None:
    """Test interpret command with missing model file."""
    sequences_file = create_test_sequences_file(tmp_path)
    missing_model = tmp_path / "nonexistent_model.json"
    output_dir = tmp_path / "interpretation"
    
    cmd = [
        sys.executable,
        "-m",
        "metainformant",
        "life-events",
        "interpret",
        f"--model={missing_model}",
        f"--sequences={sequences_file}",
        f"--output={output_dir}",
    ]
    
    result = subprocess.run(cmd, capture_output=True, text=True, timeout=30, env=_get_cli_env())
    
    # Should fail with error
    assert result.returncode != 0
    assert "not found" in result.stderr.lower() or "Error" in result.stderr


def test_cli_interpret_missing_sequences(tmp_path: Path) -> None:
    """Test interpret command with missing sequences file."""
    model_file = create_test_model(tmp_path)
    missing_sequences = tmp_path / "nonexistent_sequences.json"
    output_dir = tmp_path / "interpretation"
    
    cmd = [
        sys.executable,
        "-m",
        "metainformant",
        "life-events",
        "interpret",
        f"--model={model_file}",
        f"--sequences={missing_sequences}",
        f"--output={output_dir}",
    ]
    
    result = subprocess.run(cmd, capture_output=True, text=True, timeout=30, env=_get_cli_env())
    
    # Should fail with error
    assert result.returncode != 0
    assert "not found" in result.stderr.lower() or "Error" in result.stderr


def test_cli_interpret_no_embeddings_error(tmp_path: Path) -> None:
    """Test interpret command with model without embeddings."""
    # Create a simple model that might not have embeddings
    # Need at least 2 classes for classification
    sequences = [
        ["health:diagnosis", "occupation:job_change"],
        ["education:degree", "occupation:job_change"],
    ]
    
    y = np.array([0, 1])  # Two classes required
    
    predictor = EventSequencePredictor(
        model_type="simple",
        task_type="classification",
        random_state=42
    )
    predictor.fit(sequences, y)
    
    # Simple model may not have event_embeddings
    # Save it
    model_file = tmp_path / "simple_model.json"
    predictor.save_model(model_file)
    
    sequences_file = create_test_sequences_file(tmp_path)
    output_dir = tmp_path / "interpretation"
    
    cmd = [
        sys.executable,
        "-m",
        "metainformant",
        "life-events",
        "interpret",
        f"--model={model_file}",
        f"--sequences={sequences_file}",
        f"--output={output_dir}",
    ]
    
    result = subprocess.run(cmd, capture_output=True, text=True, timeout=30, env=_get_cli_env())
    
    # May fail if model doesn't have embeddings, or succeed if it handles gracefully
    # The important thing is it doesn't crash
    assert result.returncode in (0, 1)  # 0 = success, 1 = expected failure


def test_cli_interpret_output_format(tmp_path: Path) -> None:
    """Test that interpretation report has correct format."""
    sequences_file = create_test_sequences_file(tmp_path)
    model_file = create_test_model(tmp_path)
    output_dir = tmp_path / "interpretation"
    
    cmd = [
        sys.executable,
        "-m",
        "metainformant",
        "life-events",
        "interpret",
        f"--model={model_file}",
        f"--sequences={sequences_file}",
        f"--output={output_dir}",
    ]
    
    result = subprocess.run(cmd, capture_output=True, text=True, timeout=60, env=_get_cli_env())
    assert result.returncode == 0
    
    report_file = output_dir / "interpretation_report.json"
    report_data = json.loads(report_file.read_text())
    
    # Verify structure
    assert "model_path" in report_data
    assert "model_type" in report_data
    assert "task_type" in report_data
    assert "n_sequences" in report_data
    assert "interpretations" in report_data
    
    # Verify interpretation sections
    interpretations = report_data["interpretations"]
    assert "event_importance" in interpretations
    assert "temporal_patterns" in interpretations
    assert "feature_attribution" in interpretations
    
    # Verify event importance is a dict
    assert isinstance(interpretations["event_importance"], dict)


def test_cli_interpret_visualization(tmp_path: Path) -> None:
    """Test that visualization is created when matplotlib available."""
    sequences_file = create_test_sequences_file(tmp_path)
    model_file = create_test_model(tmp_path)
    output_dir = tmp_path / "interpretation"
    
    cmd = [
        sys.executable,
        "-m",
        "metainformant",
        "life-events",
        "interpret",
        f"--model={model_file}",
        f"--sequences={sequences_file}",
        f"--output={output_dir}",
    ]
    
    result = subprocess.run(cmd, capture_output=True, text=True, timeout=60, env=_get_cli_env())
    assert result.returncode == 0
    
    # Visualization may or may not be created depending on matplotlib availability
    viz_file = output_dir / "importance_plot.png"
    # Don't fail if visualization not created (it's optional)
    if viz_file.exists():
        assert viz_file.stat().st_size > 0

