"""Tests for phenotype behavior: Ethogram, BehaviorSequence, time budgets, diversity.

NO MOCKING POLICY: All tests use real implementations.
"""
from __future__ import annotations

import math

import pytest

from metainformant.core.utils.errors import ValidationError
from metainformant.phenotype.behavior.ethogram import BehaviorDefinition, Ethogram
from metainformant.phenotype.behavior.sequence import BehaviorEvent, BehaviorSequence


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _make_ethogram():
    return Ethogram(
        {
            "F": "Foraging",
            "G": "Grooming",
            "R": "Resting",
            "W": "Walking",
            "A": "Aggression",
        }
    )


def _make_events():
    return [
        BehaviorEvent(timestamp=0.0, code="F", duration=10.0),
        BehaviorEvent(timestamp=10.0, code="F", duration=5.0),
        BehaviorEvent(timestamp=15.0, code="G", duration=8.0),
        BehaviorEvent(timestamp=23.0, code="R", duration=12.0),
        BehaviorEvent(timestamp=35.0, code="W", duration=5.0),
        BehaviorEvent(timestamp=40.0, code="F", duration=10.0),
        BehaviorEvent(timestamp=50.0, code="G", duration=5.0),
        BehaviorEvent(timestamp=55.0, code="R", duration=15.0),
    ]


def _make_sequence():
    return BehaviorSequence(events=_make_events(), ethogram=_make_ethogram())


# ---------------------------------------------------------------------------
# Ethogram
# ---------------------------------------------------------------------------


class TestEthogram:
    def test_create_from_strings(self):
        eth = _make_ethogram()
        assert len(eth) == 5

    def test_create_from_definitions(self):
        behaviors = {
            "F": BehaviorDefinition(code="F", name="Foraging", description="Searching for food", category="feeding"),
            "G": BehaviorDefinition(code="G", name="Grooming", description="Self-cleaning"),
        }
        eth = Ethogram(behaviors)
        assert len(eth) == 2

    def test_validate_existing_code(self):
        eth = _make_ethogram()
        assert eth.validate("F") is True

    def test_validate_missing_code(self):
        eth = _make_ethogram()
        assert eth.validate("X") is False

    def test_get_existing(self):
        eth = _make_ethogram()
        defn = eth.get("F")
        assert defn is not None
        assert defn.code == "F"
        assert defn.description == "Foraging"

    def test_get_missing_returns_none(self):
        eth = _make_ethogram()
        assert eth.get("X") is None

    def test_invalid_definition_type_raises(self):
        with pytest.raises(ValidationError):
            Ethogram({"F": 123})


# ---------------------------------------------------------------------------
# BehaviorEvent
# ---------------------------------------------------------------------------


class TestBehaviorEvent:
    def test_basic_event(self):
        event = BehaviorEvent(timestamp=0.0, code="F", duration=10.0)
        assert event.timestamp == 0.0
        assert event.code == "F"
        assert event.duration == 10.0

    def test_modifiers(self):
        event = BehaviorEvent(timestamp=5.0, code="A", duration=2.0, modifiers={"target": "individual_2"})
        assert event.modifiers["target"] == "individual_2"


# ---------------------------------------------------------------------------
# BehaviorSequence - time_budget
# ---------------------------------------------------------------------------


class TestTimeBudget:
    def test_basic_budget(self):
        seq = _make_sequence()
        budget = seq.calculate_time_budget()
        assert "F" in budget
        assert "G" in budget
        total = sum(budget.values())
        assert total == pytest.approx(1.0, abs=0.001)

    def test_foraging_dominance(self):
        seq = _make_sequence()
        budget = seq.calculate_time_budget()
        # 10+5+10=25 out of 70 total
        assert budget["F"] == pytest.approx(25.0 / 70.0, abs=0.01)

    def test_empty_events(self):
        eth = _make_ethogram()
        events = [BehaviorEvent(timestamp=0.0, code="F", duration=0.0)]
        seq = BehaviorSequence(events=events, ethogram=eth)
        budget = seq.calculate_time_budget()
        assert budget == {}


# ---------------------------------------------------------------------------
# BehaviorSequence - transition_matrix
# ---------------------------------------------------------------------------


class TestTransitionMatrix:
    def test_transitions_exist(self):
        seq = _make_sequence()
        matrix = seq.transition_matrix()
        assert isinstance(matrix, dict)
        assert len(matrix) > 0

    def test_probabilities_sum_to_one(self):
        seq = _make_sequence()
        matrix = seq.transition_matrix()
        for source, targets in matrix.items():
            total = sum(targets.values())
            assert total == pytest.approx(1.0, abs=0.001)

    def test_known_transition(self):
        seq = _make_sequence()
        matrix = seq.transition_matrix()
        # F -> F exists (events 0->1), F -> G exists (events 1->2, 5->6)
        assert "F" in matrix


# ---------------------------------------------------------------------------
# BehaviorSequence - diversity
# ---------------------------------------------------------------------------


class TestDiversity:
    def test_shannon_positive(self):
        seq = _make_sequence()
        h = seq.shannon_diversity()
        assert h > 0

    def test_shannon_single_behavior(self):
        eth = Ethogram({"F": "Foraging"})
        events = [BehaviorEvent(timestamp=i, code="F", duration=1.0) for i in range(10)]
        seq = BehaviorSequence(events=events, ethogram=eth)
        h = seq.shannon_diversity()
        assert h == 0.0

    def test_simpson_positive(self):
        seq = _make_sequence()
        d = seq.simpson_diversity()
        assert 0.0 < d <= 1.0

    def test_simpson_single_behavior(self):
        eth = Ethogram({"F": "Foraging"})
        events = [BehaviorEvent(timestamp=i, code="F", duration=1.0) for i in range(10)]
        seq = BehaviorSequence(events=events, ethogram=eth)
        d = seq.simpson_diversity()
        assert d == 0.0

    def test_simpson_empty(self):
        eth = Ethogram({"F": "Foraging"})
        events = [BehaviorEvent(timestamp=0, code="F", duration=1.0)]
        seq = BehaviorSequence(events=events, ethogram=eth)
        d = seq.simpson_diversity()
        assert d == 0.0


# ---------------------------------------------------------------------------
# BehaviorSequence - bout_analysis
# ---------------------------------------------------------------------------


class TestBoutAnalysis:
    def test_basic_bouts(self):
        seq = _make_sequence()
        bouts = seq.bout_analysis()
        assert "F" in bouts
        assert bouts["F"]["bout_count"] >= 1

    def test_bout_with_min_duration(self):
        seq = _make_sequence()
        bouts = seq.bout_analysis(min_bout_duration=10.0)
        for code, info in bouts.items():
            assert info["mean_duration"] >= 10.0

    def test_consecutive_same_code(self):
        # F, F should form a single bout
        seq = _make_sequence()
        bouts = seq.bout_analysis()
        # First two events are F -> F, so bout of 15s
        assert bouts["F"]["bout_count"] >= 1


# ---------------------------------------------------------------------------
# BehaviorSequence - event_rate
# ---------------------------------------------------------------------------


class TestEventRate:
    def test_basic_rate(self):
        seq = _make_sequence()
        rate = seq.event_rate()
        assert rate > 0

    def test_custom_window(self):
        seq = _make_sequence()
        rate = seq.event_rate(time_window=100.0)
        assert rate == pytest.approx(8 / 100.0)

    def test_empty_returns_zero(self):
        eth = _make_ethogram()
        # Cannot create empty sequence because validation runs
        # Use a sequence with single point
        events = [BehaviorEvent(timestamp=0.0, code="F", duration=1.0)]
        seq = BehaviorSequence(events=events, ethogram=eth)
        rate = seq.event_rate()
        assert rate >= 0


# ---------------------------------------------------------------------------
# BehaviorSequence - markov_stationarity_chi2
# ---------------------------------------------------------------------------


class TestMarkovStationarity:
    def test_basic_chi2(self):
        seq = _make_sequence()
        result = seq.markov_stationarity_chi2()
        assert "chi2" in result
        assert "df" in result
        assert "sufficient_data" in result

    def test_insufficient_data(self):
        eth = Ethogram({"F": "Foraging", "G": "Grooming"})
        events = [
            BehaviorEvent(timestamp=0, code="F", duration=1.0),
            BehaviorEvent(timestamp=1, code="G", duration=1.0),
        ]
        seq = BehaviorSequence(events=events, ethogram=eth)
        result = seq.markov_stationarity_chi2()
        assert result["sufficient_data"] is False


# ---------------------------------------------------------------------------
# BehaviorSequence - latency_to_first
# ---------------------------------------------------------------------------


class TestLatencyToFirst:
    def test_first_behavior_zero_latency(self):
        seq = _make_sequence()
        latency = seq.latency_to_first("F")
        assert latency == 0.0

    def test_later_behavior_positive_latency(self):
        seq = _make_sequence()
        latency = seq.latency_to_first("G")
        assert latency > 0

    def test_absent_behavior_returns_none(self):
        seq = _make_sequence()
        latency = seq.latency_to_first("A")
        assert latency is None


# ---------------------------------------------------------------------------
# BehaviorSequence - behavior_counts
# ---------------------------------------------------------------------------


class TestBehaviorCounts:
    def test_basic_counts(self):
        seq = _make_sequence()
        counts = seq.behavior_counts()
        assert counts["F"] == 3
        assert counts["G"] == 2
        assert counts["R"] == 2
        assert counts["W"] == 1

    def test_total_matches_events(self):
        seq = _make_sequence()
        counts = seq.behavior_counts()
        assert sum(counts.values()) == 8


# ---------------------------------------------------------------------------
# BehaviorSequence - validation
# ---------------------------------------------------------------------------


class TestSequenceValidation:
    def test_invalid_code_raises(self):
        eth = _make_ethogram()
        events = [BehaviorEvent(timestamp=0.0, code="INVALID", duration=1.0)]
        with pytest.raises(ValidationError, match="not found in ethogram"):
            BehaviorSequence(events=events, ethogram=eth)
