"""Tests for life_events workflow functions."""

from __future__ import annotations

from datetime import datetime

import numpy as np
import pytest

from metainformant.life_events import (
    Event,
    EventSequence,
    analyze_life_course,
    compare_populations,
    intervention_analysis,
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
    results = analyze_life_course(
        sequences,
        outcomes=None,
        output_dir=output_dir
    )
    
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
    results = analyze_life_course(
        sequences,
        outcomes=outcomes,
        output_dir=output_dir,
        config_path=None
    )
    
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
    comparison = compare_populations(
        group1,
        group2,
        output_dir=output_dir
    )
    
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
    results = intervention_analysis(
        sequences,
        intervention_time=intervention_time,
        output_dir=output_dir
    )
    
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
        output_dir=output_dir
    )
    
    assert "outcome_change" in results
    assert "mean" in results["outcome_change"]

