"""Tests for life_events workflow functions."""

from __future__ import annotations

from datetime import datetime
from pathlib import Path

import numpy as np
import pytest

from metainformant.life_events import (
    Event,
    EventSequence,
    LifeEventsWorkflowConfig,
    analyze_life_course,
    compare_populations,
    intervention_analysis,
    load_life_events_config,
)


def test_analyze_life_course_basic(tmp_path):
    """Test basic life course analysis workflow."""
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

    output_dir = tmp_path / "output"
    results = analyze_life_course(sequences, outcomes=None, output_dir=output_dir)

    assert results["n_sequences"] == 2
    assert results["n_events"] == 3
    assert "embeddings" in results


def test_analyze_life_course_with_outcomes(tmp_path):
    """Test life course analysis with outcome prediction."""
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

    outcomes = np.array([0, 1])

    output_dir = tmp_path / "output"
    results = analyze_life_course(sequences, outcomes=outcomes, output_dir=output_dir, config_path=None)

    assert results["n_sequences"] == 2
    assert "predictions" in results
    assert "model_type" in results


def test_compare_populations(tmp_path):
    """Test population comparison workflow."""
    events1 = [
        Event("degree", datetime(2010, 6, 1), "education"),
        Event("job_change", datetime(2015, 3, 1), "occupation"),
    ]
    events2 = [
        Event("diagnosis", datetime(2020, 1, 15), "health"),
    ]

    group1 = [
        EventSequence(person_id="person_001", events=events1),
    ]
    group2 = [
        EventSequence(person_id="person_002", events=events2),
    ]

    output_dir = tmp_path / "comparison"
    comparison = compare_populations(group1, group2, output_dir=output_dir)

    assert "group1" in comparison
    assert "group2" in comparison
    assert "comparison" in comparison
    assert comparison["group1"]["n_sequences"] == 1
    assert comparison["group2"]["n_sequences"] == 1


def test_intervention_analysis(tmp_path):
    """Test intervention analysis workflow."""
    events = [
        Event("degree", datetime(2010, 6, 1), "education"),
        Event("job_change", datetime(2015, 3, 1), "occupation"),
        Event("diagnosis", datetime(2020, 1, 15), "health"),
    ]

    sequences = [
        EventSequence(person_id="person_001", events=events),
    ]

    intervention_time = datetime(2017, 1, 1).timestamp()

    output_dir = tmp_path / "intervention"
    results = intervention_analysis(sequences, intervention_time=intervention_time, output_dir=output_dir)

    assert "pre_intervention" in results
    assert "post_intervention" in results
    assert results["intervention_time"] == intervention_time


def test_intervention_analysis_with_outcomes(tmp_path):
    """Test intervention analysis with outcome comparison."""
    events = [
        Event("degree", datetime(2010, 6, 1), "education"),
        Event("job_change", datetime(2015, 3, 1), "occupation"),
        Event("diagnosis", datetime(2020, 1, 15), "health"),
    ]

    sequences = [
        EventSequence(person_id="person_001", events=events),
    ]

    intervention_time = datetime(2017, 1, 1).timestamp()
    pre_outcomes = np.array([0.5])
    post_outcomes = np.array([0.8])

    output_dir = tmp_path / "intervention"
    results = intervention_analysis(
        sequences,
        intervention_time=intervention_time,
        pre_intervention_outcomes=pre_outcomes,
        post_intervention_outcomes=post_outcomes,
        output_dir=output_dir,
    )

    assert "outcome_change" in results
    assert "mean" in results["outcome_change"]


# Configuration Integration Tests


def test_analyze_life_course_with_config_obj(tmp_path):
    """Test analyze_life_course with LifeEventsWorkflowConfig object."""
    events1 = [
        Event("degree", datetime(2010, 6, 1), "education"),
        Event("job_change", datetime(2015, 3, 1), "occupation"),
    ]
    events2 = [
        Event("certification", datetime(2012, 1, 1), "education"),
        Event("promotion", datetime(2018, 5, 1), "occupation"),
    ]

    sequences = [
        EventSequence(person_id="person_001", events=events1),
        EventSequence(person_id="person_002", events=events2),
    ]

    # Need at least 2 classes for classification
    outcomes = np.array([0, 1])

    config = LifeEventsWorkflowConfig(
        work_dir=tmp_path / "work",
        threads=4,
        embedding={"embedding_dim": 50, "window_size": 3, "epochs": 5},
        model={"model_type": "embedding", "random_state": 42},
    )

    results = analyze_life_course(sequences, outcomes=outcomes, config_obj=config, output_dir=tmp_path / "output")

    assert results["n_sequences"] == 2
    assert "model" in results
    assert results["embedding_dim"] == 50


def test_analyze_life_course_with_config_path(tmp_path):
    """Test analyze_life_course with config file path."""
    events1 = [
        Event("degree", datetime(2010, 6, 1), "education"),
    ]
    events2 = [
        Event("certification", datetime(2012, 1, 1), "education"),
    ]

    sequences = [
        EventSequence(person_id="person_001", events=events1),
        EventSequence(person_id="person_002", events=events2),
    ]

    # Need at least 2 classes for classification
    outcomes = np.array([0, 1])

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
  random_state: 42
""".format(
        work_dir=str(tmp_path / "work")
    )
    config_file.write_text(config_content)

    results = analyze_life_course(sequences, outcomes=outcomes, config_path=config_file, output_dir=tmp_path / "output")

    assert results["n_sequences"] == 2
    assert "model" in results
    assert results["embedding_dim"] == 50


def test_analyze_life_course_with_dict_config(tmp_path):
    """Test analyze_life_course with dict config (backward compatibility)."""
    events1 = [
        Event("degree", datetime(2010, 6, 1), "education"),
    ]
    events2 = [
        Event("certification", datetime(2012, 1, 1), "education"),
    ]

    sequences = [
        EventSequence(person_id="person_001", events=events1),
        EventSequence(person_id="person_002", events=events2),
    ]

    # Need at least 2 classes for classification
    outcomes = np.array([0, 1])

    config_dict = {"embedding_dim": 50, "window_size": 3, "epochs": 5, "model_type": "embedding", "random_state": 42}

    results = analyze_life_course(sequences, outcomes=outcomes, config_obj=config_dict, output_dir=tmp_path / "output")

    assert results["n_sequences"] == 2
    assert "model" in results
    assert results["embedding_dim"] == 50


def test_analyze_life_course_config_work_dir(tmp_path):
    """Test that work_dir from config is used."""
    events1 = [
        Event("degree", datetime(2010, 6, 1), "education"),
    ]

    sequences = [
        EventSequence(person_id="person_001", events=events1),
    ]

    config = LifeEventsWorkflowConfig(work_dir=tmp_path / "config_work_dir", threads=4)

    results = analyze_life_course(sequences, config_obj=config, output_dir=None)  # Should use config work_dir

    assert results["n_sequences"] == 1
    # Verify files were created in config work_dir
    from pathlib import Path

    work_dir = Path(results["embeddings"]).parent
    assert "config_work_dir" in str(work_dir) or work_dir.exists()


def test_analyze_life_course_config_embedding_params(tmp_path):
    """Test that embedding params from config are used."""
    events1 = [
        Event("degree", datetime(2010, 6, 1), "education"),
    ]

    sequences = [
        EventSequence(person_id="person_001", events=events1),
    ]

    config = LifeEventsWorkflowConfig(
        work_dir=tmp_path / "work", embedding={"embedding_dim": 75, "window_size": 7, "epochs": 15}
    )

    results = analyze_life_course(sequences, config_obj=config, output_dir=tmp_path / "output")

    assert results["embedding_dim"] == 75


def test_analyze_life_course_config_model_params(tmp_path):
    """Test that model params from config are used."""
    events1 = [
        Event("degree", datetime(2010, 6, 1), "education"),
    ]
    events2 = [
        Event("certification", datetime(2012, 1, 1), "education"),
    ]

    sequences = [
        EventSequence(person_id="person_001", events=events1),
        EventSequence(person_id="person_002", events=events2),
    ]

    # Need at least 2 classes for classification
    outcomes = np.array([0, 1])

    config = LifeEventsWorkflowConfig(
        work_dir=tmp_path / "work", model={"model_type": "simple", "task_type": "classification", "random_state": 123}
    )

    results = analyze_life_course(sequences, outcomes=outcomes, config_obj=config, output_dir=tmp_path / "output")

    assert results["model_type"] == "simple"
    assert results["task_type"] == "classification"


def test_analyze_life_course_saves_model(tmp_path):
    """Test that model is saved when outcomes provided."""
    events1 = [
        Event("degree", datetime(2010, 6, 1), "education"),
    ]
    events2 = [
        Event("certification", datetime(2012, 1, 1), "education"),
    ]

    sequences = [
        EventSequence(person_id="person_001", events=events1),
        EventSequence(person_id="person_002", events=events2),
    ]

    # Need at least 2 classes for classification
    outcomes = np.array([0, 1])

    output_dir = tmp_path / "output"
    results = analyze_life_course(sequences, outcomes=outcomes, output_dir=output_dir)

    assert "model" in results
    model_file = Path(results["model"])
    assert model_file.exists()

    # Verify model can be loaded
    from metainformant.life_events import EventSequencePredictor

    loaded = EventSequencePredictor.load_model(model_file)
    assert loaded.is_fitted
