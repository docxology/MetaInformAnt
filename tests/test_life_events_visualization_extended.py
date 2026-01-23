"""Tests for extended visualization functions in life_events module."""

from __future__ import annotations

from datetime import datetime
from pathlib import Path

import numpy as np
import pytest

from metainformant.life_events import (
    Event,
    EventSequence,
    plot_domain_distribution,
    plot_domain_timeline,
    plot_embedding_clusters,
    plot_event_cooccurrence,
    plot_event_frequency_heatmap,
    plot_intervention_effects,
    plot_outcome_distribution,
    plot_population_comparison,
    plot_prediction_accuracy,
    plot_sequence_length_distribution,
    plot_sequence_similarity,
    plot_temporal_density,
    plot_temporal_patterns,
    plot_transition_network,
)


@pytest.fixture
def sample_sequences():
    """Create sample event sequences for testing."""
    sequences = []
    for i in range(5):
        seq = EventSequence(
            f"person_{i+1:03d}",
            [
                Event("degree", datetime(2010 + i, 1, 1), "education"),
                Event("job_change", datetime(2015 + i, 1, 1), "occupation"),
                Event("diagnosis", datetime(2020 + i, 1, 1), "health"),
            ],
        )
        sequences.append(seq)
    return sequences


def test_plot_domain_distribution(tmp_path, sample_sequences):
    """Test domain distribution plotting."""
    output_path = tmp_path / "domain_dist.png"
    fig = plot_domain_distribution(sample_sequences, output_path=output_path, plot_type="bar")
    assert output_path.exists()
    assert fig is not None


def test_plot_temporal_density(tmp_path, sample_sequences):
    """Test temporal density plotting."""
    output_path = tmp_path / "temporal_density.png"
    fig = plot_temporal_density(sample_sequences, output_path=output_path, bins=20)
    assert output_path.exists()
    assert fig is not None


def test_plot_event_cooccurrence(tmp_path, sample_sequences):
    """Test event co-occurrence plotting."""
    output_path = tmp_path / "cooccurrence.png"
    fig = plot_event_cooccurrence(sample_sequences, output_path=output_path, top_n=10)
    assert output_path.exists()
    assert fig is not None


def test_plot_outcome_distribution(tmp_path):
    """Test outcome distribution plotting."""
    outcomes = np.array([0, 1, 0, 1, 0, 1, 0])
    output_path = tmp_path / "outcome_dist.png"
    fig = plot_outcome_distribution(outcomes, output_path=output_path, plot_type="histogram")
    assert output_path.exists()
    assert fig is not None


def test_plot_sequence_similarity(tmp_path, sample_sequences):
    """Test sequence similarity plotting."""
    output_path = tmp_path / "similarity.png"
    try:
        fig = plot_sequence_similarity(sample_sequences, output_path=output_path)
        assert output_path.exists()
        assert fig is not None
    except ImportError:
        pytest.skip("sklearn not available for similarity computation")


def test_plot_transition_network(tmp_path, sample_sequences):
    """Test transition network plotting."""
    output_path = tmp_path / "transition_network.png"
    try:
        fig = plot_transition_network(sample_sequences, output_path=output_path, top_n=5)
        assert output_path.exists()
        assert fig is not None
    except ImportError:
        pytest.skip("networkx not available for network visualization")


def test_plot_domain_timeline(tmp_path, sample_sequences):
    """Test domain timeline plotting."""
    output_path = tmp_path / "domain_timeline.png"
    fig = plot_domain_timeline(sample_sequences, output_path=output_path, max_sequences=5)
    assert output_path.exists()
    assert fig is not None


def test_plot_prediction_accuracy(tmp_path):
    """Test prediction accuracy plotting."""
    y_true = np.array([0, 1, 0, 1, 0])
    y_pred = np.array([0, 1, 0, 0, 1])

    output_path = tmp_path / "prediction_accuracy.png"
    fig = plot_prediction_accuracy(y_true, y_pred, task_type="classification", output_path=output_path)
    assert output_path.exists()
    assert fig is not None


def test_plot_temporal_patterns(tmp_path, sample_sequences):
    """Test temporal patterns plotting."""
    importance_scores = {0: 0.8, 1: 0.9, 2: 0.7}
    output_path = tmp_path / "temporal_patterns.png"
    fig = plot_temporal_patterns(sample_sequences, importance_scores=importance_scores, output_path=output_path)
    assert output_path.exists()
    assert fig is not None


def test_plot_population_comparison(tmp_path, sample_sequences):
    """Test population comparison plotting."""
    group1 = sample_sequences[:2]
    group2 = sample_sequences[2:]

    output_path = tmp_path / "population_comparison.png"
    fig = plot_population_comparison(
        group1, group2, group1_label="Group 1", group2_label="Group 2", output_path=output_path
    )
    assert output_path.exists()
    assert fig is not None


def test_plot_intervention_effects(tmp_path, sample_sequences):
    """Test intervention effects plotting."""
    pre_sequences = sample_sequences[:2]
    post_sequences = sample_sequences[2:]
    pre_outcomes = np.array([0.5, 0.6])
    post_outcomes = np.array([0.8, 0.9])

    output_path = tmp_path / "intervention_effects.png"
    fig = plot_intervention_effects(
        pre_sequences, post_sequences, pre_outcomes=pre_outcomes, post_outcomes=post_outcomes, output_path=output_path
    )
    assert output_path.exists()
    assert fig is not None


def test_plot_embedding_clusters(tmp_path):
    """Test embedding clusters plotting."""
    embeddings = {
        "health:diagnosis": np.random.rand(50),
        "health:hospitalization": np.random.rand(50),
        "education:degree": np.random.rand(50),
    }

    clusters = {
        "health:diagnosis": 0,
        "health:hospitalization": 0,
        "education:degree": 1,
    }

    output_path = tmp_path / "embedding_clusters.png"
    fig = plot_embedding_clusters(embeddings, clusters=clusters, output_path=output_path)
    assert output_path.exists()
    assert fig is not None


def test_plot_sequence_length_distribution(tmp_path, sample_sequences):
    """Test sequence length distribution plotting."""
    output_path = tmp_path / "sequence_lengths.png"
    fig = plot_sequence_length_distribution(sample_sequences, output_path=output_path)
    assert output_path.exists()
    assert fig is not None


def test_plot_event_frequency_heatmap(tmp_path, sample_sequences):
    """Test event frequency heatmap plotting."""
    output_path = tmp_path / "frequency_heatmap.png"
    fig = plot_event_frequency_heatmap(sample_sequences, output_path=output_path, time_bins=5)
    assert output_path.exists()
    assert fig is not None
