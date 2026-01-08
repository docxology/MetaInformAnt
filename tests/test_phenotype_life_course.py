"""Comprehensive tests for life_course phenotype extraction functions."""

from __future__ import annotations

from datetime import datetime
from pathlib import Path

import pytest

try:
    from metainformant.life_events import Event, EventSequence
    from metainformant.phenotype import (
        aggregate_temporal_phenotypes,
        extract_phenotypes_from_events,
        map_events_to_traits,
    )
    from metainformant.core.utils.errors import ValidationError
    
    LIFE_EVENTS_AVAILABLE = True
except ImportError:
    LIFE_EVENTS_AVAILABLE = False
    pytestmark = pytest.mark.skip("life_events module not available")


@pytest.mark.skipif(not LIFE_EVENTS_AVAILABLE, reason="life_events module not available")
class TestExtractPhenotypesFromEvents:
    """Tests for extract_phenotypes_from_events function."""

    def test_basic_extraction(self):
        """Test basic phenotype extraction from event sequence."""
        events = [
            Event("diabetes", datetime(2020, 1, 1), "health"),
            Event("bachelors", datetime(2010, 6, 1), "education"),
            Event("job_change", datetime(2015, 3, 15), "occupation"),
        ]
        sequence = EventSequence(person_id="p1", events=events)
        
        phenotypes = extract_phenotypes_from_events(sequence)
        
        assert phenotypes["person_id"] == "p1"
        assert phenotypes["total_events"] == 3
        assert "health" in phenotypes["domains"]
        assert "education" in phenotypes["domains"]
        assert "occupation" in phenotypes["domains"]
        assert phenotypes["health_events"] == 1
        assert phenotypes["education_events"] == 1
        assert phenotypes["occupation_events"] == 1
        assert "diabetes" in phenotypes["health_conditions"]
        assert "bachelors" in phenotypes["education_achievements"]
        assert "job_change" in phenotypes["occupation_changes"]
        
    def test_empty_sequence(self):
        """Test extraction from empty sequence."""
        sequence = EventSequence(person_id="p1", events=[])
        
        phenotypes = extract_phenotypes_from_events(sequence)
        
        assert phenotypes["person_id"] == "p1"
        assert phenotypes["total_events"] == 0
        assert phenotypes["domains"] == []
        assert phenotypes["event_types"] == []
        assert "domain_counts" in phenotypes
        
    def test_temporal_aggregation(self):
        """Test temporal information extraction."""
        events = [
            Event("event1", datetime(2010, 1, 1), "health"),
            Event("event2", datetime(2020, 1, 1), "education"),
        ]
        sequence = EventSequence(person_id="p1", events=events)
        
        phenotypes = extract_phenotypes_from_events(sequence)
        
        assert "first_event_time" in phenotypes
        assert "last_event_time" in phenotypes
        assert "event_span_years" in phenotypes
        assert phenotypes["event_span_years"] > 0
        assert phenotypes["last_event_time"] > phenotypes["first_event_time"]
        
    def test_custom_phenotype_categories(self):
        """Test custom phenotype categories."""
        events = [
            Event("diagnosis", datetime(2020, 1, 1), "health"),
            Event("treatment", datetime(2020, 2, 1), "health"),
        ]
        sequence = EventSequence(person_id="p1", events=events)
        
        custom_categories = {
            "medical_events": ["diagnosis", "treatment", "hospitalization"]
        }
        
        phenotypes = extract_phenotypes_from_events(sequence, phenotype_categories=custom_categories)
        
        assert "medical_events_count" in phenotypes
        assert phenotypes["medical_events_count"] == 2
        
    def test_none_sequence(self):
        """Test that None sequence raises ValidationError."""
        with pytest.raises(ValidationError):
            extract_phenotypes_from_events(None)
            
    def test_invalid_custom_categories(self):
        """Test invalid custom categories."""
        events = [Event("event1", datetime(2020, 1, 1), "health")]
        sequence = EventSequence(person_id="p1", events=events)
        
        invalid_categories = {
            "test": "not a list"
        }
        
        with pytest.raises(ValidationError):
            extract_phenotypes_from_events(sequence, phenotype_categories=invalid_categories)


@pytest.mark.skipif(not LIFE_EVENTS_AVAILABLE, reason="life_events module not available")
class TestAggregateTemporalPhenotypes:
    """Tests for aggregate_temporal_phenotypes function."""

    def test_basic_aggregation(self):
        """Test basic temporal aggregation."""
        seq1 = EventSequence(
            person_id="p1",
            events=[Event("event1", datetime(2010, 1, 1), "health")]
        )
        seq2 = EventSequence(
            person_id="p2",
            events=[Event("event2", datetime(2015, 1, 1), "education")]
        )
        
        result = aggregate_temporal_phenotypes([seq1, seq2], time_window_years=5.0)
        
        assert "time_windows" in result
        assert "aggregates" in result
        assert result["aggregates"]["total_events"] == 2
        assert result["aggregates"]["total_people"] == 2
        assert result["aggregates"]["time_span_years"] > 0
        
    def test_empty_sequences(self):
        """Test aggregation with empty sequences."""
        seq1 = EventSequence(person_id="p1", events=[])
        seq2 = EventSequence(person_id="p2", events=[])
        
        result = aggregate_temporal_phenotypes([seq1, seq2], time_window_years=5.0)
        
        assert result["time_windows"] == []
        assert result["aggregates"]["total_events"] == 0
        assert result["aggregates"]["total_people"] == 0
        assert result["aggregates"]["time_span_years"] == 0.0
        
    def test_time_windows(self):
        """Test time window creation."""
        events = [
            Event("event1", datetime(2010, 1, 1), "health"),
            Event("event2", datetime(2015, 1, 1), "education"),
            Event("event3", datetime(2020, 1, 1), "occupation"),
        ]
        sequence = EventSequence(person_id="p1", events=events)
        
        result = aggregate_temporal_phenotypes([sequence], time_window_years=5.0)
        
        assert len(result["time_windows"]) > 0
        for window in result["time_windows"]:
            assert "start_time" in window
            assert "end_time" in window
            assert "n_events" in window
            assert "n_people" in window
            assert "domain_counts" in window
            
    def test_invalid_sequences_type(self):
        """Test that non-list sequences raises ValidationError."""
        with pytest.raises(ValidationError):
            aggregate_temporal_phenotypes("not a list", time_window_years=5.0)
            
    def test_invalid_time_window(self):
        """Test that invalid time window raises ValidationError."""
        seq = EventSequence(person_id="p1", events=[Event("event1", datetime(2020, 1, 1), "health")])
        
        with pytest.raises(ValidationError):
            aggregate_temporal_phenotypes([seq], time_window_years=0)
            
        with pytest.raises(ValidationError):
            aggregate_temporal_phenotypes([seq], time_window_years=-1)


@pytest.mark.skipif(not LIFE_EVENTS_AVAILABLE, reason="life_events module not available")
class TestMapEventsToTraits:
    """Tests for map_events_to_traits function."""

    def test_default_mapping(self):
        """Test default trait mapping."""
        events = [
            Event("diagnosis", datetime(2020, 1, 1), "health"),
            Event("degree", datetime(2010, 1, 1), "education"),
            Event("job_change", datetime(2015, 1, 1), "occupation"),
        ]
        sequence = EventSequence(person_id="p1", events=events)
        
        traits = map_events_to_traits(sequence)
        
        assert "health_issues" in traits
        assert "education_level" in traits
        assert "career_progression" in traits
        assert traits["health_issues"]["count"] == 1
        assert traits["education_level"]["count"] == 1
        assert traits["career_progression"]["count"] == 1
        assert "diagnosis" in traits["health_issues"]["events"]
        assert "degree" in traits["education_level"]["events"]
        assert "job_change" in traits["career_progression"]["events"]
        
    def test_custom_mapping(self):
        """Test custom trait mapping."""
        events = [
            Event("custom_event", datetime(2020, 1, 1), "health"),
        ]
        sequence = EventSequence(person_id="p1", events=events)
        
        custom_mapping = {
            "custom_trait": ["custom_event", "other_event"]
        }
        
        traits = map_events_to_traits(sequence, trait_mapping=custom_mapping)
        
        assert "custom_trait" in traits
        assert traits["custom_trait"]["count"] == 1
        assert "custom_event" in traits["custom_trait"]["events"]
        
    def test_timestamps(self):
        """Test that timestamps are included."""
        events = [
            Event("diagnosis", datetime(2020, 1, 1), "health"),
        ]
        sequence = EventSequence(person_id="p1", events=events)
        
        traits = map_events_to_traits(sequence)
        
        assert "timestamps" in traits["health_issues"]
        assert len(traits["health_issues"]["timestamps"]) == 1
        
    def test_none_sequence(self):
        """Test that None sequence raises ValidationError."""
        with pytest.raises(ValidationError):
            map_events_to_traits(None)
            
    def test_invalid_trait_mapping(self):
        """Test invalid trait mapping."""
        events = [Event("event1", datetime(2020, 1, 1), "health")]
        sequence = EventSequence(person_id="p1", events=events)
        
        invalid_mapping = {
            "trait": "not a list"
        }
        
        with pytest.raises(ValidationError):
            map_events_to_traits(sequence, trait_mapping=invalid_mapping)
            
    def test_no_matching_events(self):
        """Test trait mapping with no matching events."""
        events = [
            Event("unrelated_event", datetime(2020, 1, 1), "health"),
        ]
        sequence = EventSequence(person_id="p1", events=events)
        
        traits = map_events_to_traits(sequence)
        
        # Default mapping should still create entries with count 0
        assert "health_issues" in traits
        assert traits["health_issues"]["count"] == 0
        assert traits["health_issues"]["events"] == []


@pytest.mark.skipif(not LIFE_EVENTS_AVAILABLE, reason="life_events module not available")
class TestIntegration:
    """Integration tests for life course functions."""

    def test_full_workflow(self):
        """Test complete workflow from sequences to phenotypes."""
        sequences = [
            EventSequence(
                person_id="p1",
                events=[
                    Event("diagnosis", datetime(2020, 1, 1), "health"),
                    Event("degree", datetime(2010, 1, 1), "education"),
                ]
            ),
            EventSequence(
                person_id="p2",
                events=[
                    Event("treatment", datetime(2021, 1, 1), "health"),
                    Event("job_change", datetime(2015, 1, 1), "occupation"),
                ]
            ),
        ]
        
        # Extract phenotypes from each
        phenotypes_list = []
        for seq in sequences:
            phenotypes = extract_phenotypes_from_events(seq)
            traits = map_events_to_traits(seq)
            phenotypes_list.append(phenotypes)
            
        # Aggregate temporal patterns
        aggregated = aggregate_temporal_phenotypes(sequences, time_window_years=5.0)
        
        assert len(phenotypes_list) == 2
        assert aggregated["aggregates"]["total_people"] == 2
        assert aggregated["aggregates"]["total_events"] == 4

