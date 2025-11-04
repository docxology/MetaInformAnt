"""Tests for life_events visualization functions."""

from __future__ import annotations

from datetime import datetime
from pathlib import Path

import numpy as np
import pytest

from metainformant.life_events import Event, EventSequence


def test_plot_event_timeline(tmp_path):
    """Test event timeline plotting."""
    try:
        from metainformant.life_events import plot_event_timeline
    except ImportError:
        pytest.skip("Visualization functions require matplotlib")

    events = [
        Event("degree", datetime(2010, 6, 1), "education"),
        Event("job_change", datetime(2015, 3, 1), "occupation"),
        Event("diagnosis", datetime(2020, 1, 15), "health"),
    ]

    sequence = EventSequence(person_id="person_001", events=events)

    output_path = tmp_path / "timeline.png"
    fig = plot_event_timeline(sequence, output_path=output_path)

    assert fig is not None
    assert output_path.exists()


def test_plot_event_embeddings(tmp_path):
    """Test event embeddings plotting."""
    try:
        from metainformant.life_events import learn_event_embeddings, plot_event_embeddings
    except ImportError:
        pytest.skip("Visualization functions require matplotlib")

    sequences = [
        ["health:diagnosis", "occupation:job_change"],
        ["education:degree", "occupation:job_change"],
    ]

    embeddings = learn_event_embeddings(sequences, embedding_dim=50, random_state=42)

    output_path = tmp_path / "embeddings.png"
    fig = plot_event_embeddings(
        embeddings,
        method="pca",
        n_components=2,
        output_path=output_path
    )

    assert fig is not None
    assert output_path.exists()


def test_plot_event_embeddings_3d(tmp_path):
    """Test 3D event embeddings plotting."""
    try:
        from metainformant.life_events import learn_event_embeddings, plot_event_embeddings
    except ImportError:
        pytest.skip("Visualization functions require matplotlib")

    sequences = [
        ["health:diagnosis", "occupation:job_change"],
        ["education:degree", "occupation:job_change"],
    ]

    embeddings = learn_event_embeddings(sequences, embedding_dim=50, random_state=42)

    output_path = tmp_path / "embeddings_3d.png"
    fig = plot_event_embeddings(
        embeddings,
        method="pca",
        n_components=3,
        output_path=output_path
    )

    assert fig is not None


def test_plot_attention_heatmap(tmp_path):
    """Test attention heatmap plotting."""
    try:
        from metainformant.life_events import plot_attention_heatmap
    except ImportError:
        pytest.skip("Visualization functions require matplotlib")

    attention_weights = np.array([
        [0.5, 0.3, 0.2],
        [0.3, 0.4, 0.3],
        [0.2, 0.3, 0.5],
    ])

    event_tokens = ["health:diagnosis", "occupation:job_change", "education:degree"]

    output_path = tmp_path / "attention.png"
    fig = plot_attention_heatmap(
        attention_weights,
        event_tokens,
        output_path=output_path
    )

    assert fig is not None
    assert output_path.exists()


def test_plot_prediction_importance(tmp_path):
    """Test prediction importance plotting."""
    try:
        from metainformant.life_events import plot_prediction_importance
    except ImportError:
        pytest.skip("Visualization functions require matplotlib")

    importance = {
        "health:diagnosis": 0.8,
        "occupation:job_change": 0.6,
        "education:degree": 0.4,
        "income:raise": 0.3,
    }

    output_path = tmp_path / "importance.png"
    fig = plot_prediction_importance(
        importance,
        top_n=3,
        output_path=output_path
    )

    assert fig is not None
    assert output_path.exists()


def test_plot_prediction_importance_empty():
    """Test prediction importance plotting with empty dict."""
    try:
        from metainformant.life_events import plot_prediction_importance
    except ImportError:
        pytest.skip("Visualization functions require matplotlib")

    importance = {}

    fig = plot_prediction_importance(importance, top_n=10)
    assert fig is not None

