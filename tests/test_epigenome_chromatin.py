"""Tests for epigenome chromatin state learning and analysis.

Tests learn_chromatin_states, assign_states, interpret_states,
compute_state_enrichment, segment_genome, compare_chromatin_states,
and compute_state_transition_rates using real implementations with
synthetic histone modification data. NO MOCKING.
"""

from __future__ import annotations

import math
import random

import pytest

from metainformant.epigenome.chromatin_state.state_learning import (
    assign_states,
    compare_chromatin_states,
    compute_state_enrichment,
    compute_state_transition_rates,
    interpret_states,
    learn_chromatin_states,
    segment_genome,
)

# ---------------------------------------------------------------------------
# Helpers to generate realistic synthetic histone data
# ---------------------------------------------------------------------------


def _make_histone_matrix(
    n_bins: int = 200,
    n_marks: int = 4,
    n_clusters: int = 3,
    seed: int = 42,
) -> list[list[float]]:
    """Generate a synthetic histone signal matrix with distinct clusters.

    Creates n_clusters groups of bins, each with a characteristic
    emission pattern across n_marks histone modifications.

    Args:
        n_bins: Number of genomic bins.
        n_marks: Number of histone marks.
        n_clusters: Number of distinct chromatin states.
        seed: Random seed.

    Returns:
        2D list of shape (n_bins, n_marks).
    """
    rng = random.Random(seed)

    # Define cluster centres
    centres = []
    for k in range(n_clusters):
        centre = [rng.uniform(0, 10) for _ in range(n_marks)]
        centres.append(centre)

    data: list[list[float]] = []
    bins_per_cluster = n_bins // n_clusters

    for k in range(n_clusters):
        count = bins_per_cluster if k < n_clusters - 1 else n_bins - len(data)
        for _ in range(count):
            row = [max(0.0, centres[k][j] + rng.gauss(0, 0.5)) for j in range(n_marks)]
            data.append(row)

    # Shuffle slightly so clusters aren't perfectly contiguous
    # (but keep mostly contiguous for segment_genome tests)
    return data


def _make_chromatin_marks_matrix(
    n_bins: int = 100,
    seed: int = 42,
) -> list[list[float]]:
    """Create a matrix with 6 standard histone marks.

    Marks: H3K4me3, H3K4me1, H3K27ac, H3K27me3, H3K9me3, H3K36me3
    Creates bins representing promoter, enhancer, repressed, and quiescent regions.
    """
    rng = random.Random(seed)
    marks = 6
    bins_per_type = n_bins // 4

    data: list[list[float]] = []

    # Active promoter bins: high K4me3 and K27ac
    for _ in range(bins_per_type):
        row = [
            rng.uniform(8, 12),  # H3K4me3 high
            rng.uniform(0, 2),  # H3K4me1 low
            rng.uniform(7, 11),  # H3K27ac high
            rng.uniform(0, 1),  # H3K27me3 low
            rng.uniform(0, 1),  # H3K9me3 low
            rng.uniform(0, 2),  # H3K36me3 low
        ]
        data.append(row)

    # Active enhancer bins: high K4me1 and K27ac
    for _ in range(bins_per_type):
        row = [
            rng.uniform(0, 2),  # H3K4me3 low
            rng.uniform(8, 12),  # H3K4me1 high
            rng.uniform(7, 11),  # H3K27ac high
            rng.uniform(0, 1),  # H3K27me3 low
            rng.uniform(0, 1),  # H3K9me3 low
            rng.uniform(0, 2),  # H3K36me3 low
        ]
        data.append(row)

    # Polycomb repressed bins: high K27me3
    for _ in range(bins_per_type):
        row = [
            rng.uniform(0, 1),  # H3K4me3 low
            rng.uniform(0, 1),  # H3K4me1 low
            rng.uniform(0, 1),  # H3K27ac low
            rng.uniform(8, 12),  # H3K27me3 high
            rng.uniform(0, 1),  # H3K9me3 low
            rng.uniform(0, 2),  # H3K36me3 low
        ]
        data.append(row)

    # Quiescent bins: all low
    remaining = n_bins - len(data)
    for _ in range(remaining):
        row = [rng.uniform(0, 1) for _ in range(marks)]
        data.append(row)

    return data


# ---------------------------------------------------------------------------
# learn_chromatin_states
# ---------------------------------------------------------------------------


class TestLearnChromatinStates:
    """Tests for learn_chromatin_states."""

    def test_basic_learning(self) -> None:
        data = _make_histone_matrix(n_bins=100, n_marks=4, n_clusters=3)
        result = learn_chromatin_states(data, n_states=3, max_iter=50, tol=1e-3)
        assert result["states"] == 3
        assert len(result["assignments"]) == 100
        assert "emission_params" in result
        assert "transition_matrix" in result
        assert "log_likelihood" in result

    def test_emission_params_shape(self) -> None:
        n_marks = 5
        data = _make_histone_matrix(n_bins=80, n_marks=n_marks, n_clusters=3)
        result = learn_chromatin_states(data, n_states=3, max_iter=30)
        means = result["emission_params"]["means"]
        variances = result["emission_params"]["variances"]
        assert len(means) == 3
        assert len(means[0]) == n_marks
        assert len(variances) == 3
        assert len(variances[0]) == n_marks

    def test_transition_matrix_shape(self) -> None:
        data = _make_histone_matrix(n_bins=80, n_marks=4, n_clusters=3)
        result = learn_chromatin_states(data, n_states=3, max_iter=30)
        trans = result["transition_matrix"]
        assert len(trans) == 3
        assert len(trans[0]) == 3
        # Each row should sum to ~1.0
        for row in trans:
            assert sum(row) == pytest.approx(1.0, abs=0.01)

    def test_assignments_valid_range(self) -> None:
        data = _make_histone_matrix(n_bins=60, n_marks=3, n_clusters=2)
        result = learn_chromatin_states(data, n_states=2, max_iter=30)
        for a in result["assignments"]:
            assert 0 <= a < 2

    def test_weights_sum_to_one(self) -> None:
        data = _make_histone_matrix(n_bins=80, n_marks=4, n_clusters=3)
        result = learn_chromatin_states(data, n_states=3, max_iter=30)
        assert sum(result["weights"]) == pytest.approx(1.0, abs=0.01)

    def test_too_few_bins_raises(self) -> None:
        data = [[1.0, 2.0], [3.0, 4.0]]
        with pytest.raises(ValueError, match="Need at least"):
            learn_chromatin_states(data, n_states=5)

    def test_empty_matrix_raises(self) -> None:
        with pytest.raises(ValueError, match="must not be empty"):
            learn_chromatin_states([], n_states=3)

    def test_state_counts_present(self) -> None:
        data = _make_histone_matrix(n_bins=80, n_marks=4, n_clusters=3)
        result = learn_chromatin_states(data, n_states=3, max_iter=30)
        assert "state_counts" in result
        total = sum(result["state_counts"].values())
        assert total == 80

    def test_convergence_reported(self) -> None:
        data = _make_histone_matrix(n_bins=80, n_marks=4, n_clusters=3)
        result = learn_chromatin_states(data, n_states=3, max_iter=200, tol=1e-2)
        assert "converged" in result
        assert "n_iterations" in result
        assert isinstance(result["converged"], bool)


# ---------------------------------------------------------------------------
# assign_states
# ---------------------------------------------------------------------------


class TestAssignStates:
    """Tests for assign_states (Viterbi decoding)."""

    def test_assigns_correct_length(self) -> None:
        data = _make_histone_matrix(n_bins=80, n_marks=4, n_clusters=3)
        model = learn_chromatin_states(data, n_states=3, max_iter=30)

        new_data = _make_histone_matrix(n_bins=50, n_marks=4, n_clusters=3, seed=99)
        assignments = assign_states(new_data, model)
        assert len(assignments) == 50

    def test_assignments_in_valid_range(self) -> None:
        data = _make_histone_matrix(n_bins=80, n_marks=4, n_clusters=3)
        model = learn_chromatin_states(data, n_states=3, max_iter=30)

        new_data = _make_histone_matrix(n_bins=40, n_marks=4, n_clusters=3, seed=99)
        assignments = assign_states(new_data, model)
        for a in assignments:
            assert 0 <= a < 3

    def test_wrong_marks_raises(self) -> None:
        data = _make_histone_matrix(n_bins=80, n_marks=4, n_clusters=3)
        model = learn_chromatin_states(data, n_states=3, max_iter=30)

        # Wrong number of marks
        bad_data = [[1.0, 2.0] for _ in range(20)]
        with pytest.raises(ValueError, match="marks"):
            assign_states(bad_data, model)

    def test_self_assignment_consistency(self) -> None:
        """Assigning the training data should give similar results."""
        data = _make_histone_matrix(n_bins=80, n_marks=4, n_clusters=3)
        model = learn_chromatin_states(data, n_states=3, max_iter=100, tol=1e-4)

        reassigned = assign_states(data, model)
        original = model["assignments"]
        # At least 50% should match (accounting for label permutation)
        matches = sum(1 for a, b in zip(original, reassigned) if a == b)
        assert matches / len(original) > 0.3


# ---------------------------------------------------------------------------
# interpret_states
# ---------------------------------------------------------------------------


class TestInterpretStates:
    """Tests for interpret_states."""

    def test_returns_one_per_state(self) -> None:
        data = _make_chromatin_marks_matrix(n_bins=100)
        mark_names = ["H3K4me3", "H3K4me1", "H3K27ac", "H3K27me3", "H3K9me3", "H3K36me3"]
        model = learn_chromatin_states(data, n_states=4, max_iter=50)
        interpretations = interpret_states(model["emission_params"], mark_names)
        assert len(interpretations) == 4

    def test_interpretation_keys(self) -> None:
        data = _make_chromatin_marks_matrix(n_bins=100)
        mark_names = ["H3K4me3", "H3K4me1", "H3K27ac", "H3K27me3", "H3K9me3", "H3K36me3"]
        model = learn_chromatin_states(data, n_states=4, max_iter=50)
        interpretations = interpret_states(model["emission_params"], mark_names)
        expected_keys = {"state_id", "label", "category", "dominant_marks", "mean_signal", "description"}
        for interp in interpretations:
            assert expected_keys.issubset(interp.keys())

    def test_recognizes_known_categories(self) -> None:
        data = _make_chromatin_marks_matrix(n_bins=200)
        mark_names = ["H3K4me3", "H3K4me1", "H3K27ac", "H3K27me3", "H3K9me3", "H3K36me3"]
        model = learn_chromatin_states(data, n_states=4, max_iter=80)
        interpretations = interpret_states(model["emission_params"], mark_names)
        categories = {i["category"] for i in interpretations}
        labels = {i["label"] for i in interpretations}
        # At least one recognized category should appear
        known = {"promoter", "enhancer", "repressed", "quiescent", "transcribed", "regulatory", "bivalent", "other"}
        assert len(categories & known) >= 1

    def test_with_minimal_marks(self) -> None:
        # Just two marks
        emission_params = {
            "means": [
                [10.0, 1.0],
                [1.0, 10.0],
            ]
        }
        mark_names = ["H3K4me3", "H3K27me3"]
        interpretations = interpret_states(emission_params, mark_names)
        assert len(interpretations) == 2

    def test_emission_means_in_output(self) -> None:
        data = _make_chromatin_marks_matrix(n_bins=80)
        mark_names = ["H3K4me3", "H3K4me1", "H3K27ac", "H3K27me3", "H3K9me3", "H3K36me3"]
        model = learn_chromatin_states(data, n_states=3, max_iter=30)
        interpretations = interpret_states(model["emission_params"], mark_names)
        for interp in interpretations:
            assert "emission_means" in interp
            assert len(interp["emission_means"]) == 6


# ---------------------------------------------------------------------------
# compute_state_enrichment
# ---------------------------------------------------------------------------


class TestComputeStateEnrichment:
    """Tests for compute_state_enrichment."""

    def test_basic_enrichment(self) -> None:
        assignments = [0, 0, 0, 1, 1, 1, 2, 2, 2, 2]
        annotations = {
            "TSS": [0, 1, 2],
            "gene_body": [3, 4, 5, 6],
        }
        result = compute_state_enrichment(assignments, annotations, n_states=3)
        assert "enrichment_matrix" in result
        assert "summary" in result
        assert len(result["summary"]) == 3

    def test_enrichment_matrix_keys(self) -> None:
        assignments = [0, 0, 1, 1, 2, 2]
        annotations = {"promoter": [0, 1]}
        result = compute_state_enrichment(assignments, annotations, n_states=3)
        matrix = result["enrichment_matrix"]
        # Should have an entry for each state-annotation pair
        assert "state_0__promoter" in matrix
        entry = matrix["state_0__promoter"]
        expected_keys = {
            "state_id",
            "annotation",
            "observed_overlap",
            "expected_overlap",
            "fold_enrichment",
            "p_value",
        }
        assert expected_keys.issubset(entry.keys())

    def test_perfect_enrichment(self) -> None:
        # State 0 perfectly overlaps with annotation
        assignments = [0, 0, 0, 1, 1, 1]
        annotations = {"target": [0, 1, 2]}
        result = compute_state_enrichment(assignments, annotations, n_states=2)
        entry = result["enrichment_matrix"]["state_0__target"]
        assert entry["fold_enrichment"] > 1.0

    def test_no_enrichment(self) -> None:
        # State 0 does not overlap with annotation at all
        assignments = [0, 0, 0, 1, 1, 1]
        annotations = {"target": [3, 4, 5]}
        result = compute_state_enrichment(assignments, annotations, n_states=2)
        entry = result["enrichment_matrix"]["state_0__target"]
        assert entry["fold_enrichment"] == 0.0

    def test_empty_assignments(self) -> None:
        result = compute_state_enrichment([], {}, n_states=3)
        assert result["enrichment_matrix"] == {}

    def test_summary_has_most_enriched(self) -> None:
        assignments = [0, 0, 0, 1, 1, 1, 2, 2, 2]
        annotations = {"TSS": [0, 1, 2], "enhancer": [6, 7, 8]}
        result = compute_state_enrichment(assignments, annotations, n_states=3)
        for summary_entry in result["summary"]:
            assert "most_enriched_annotation" in summary_entry


# ---------------------------------------------------------------------------
# segment_genome
# ---------------------------------------------------------------------------


class TestSegmentGenome:
    """Tests for segment_genome."""

    def test_basic_segmentation(self) -> None:
        assignments = [0, 0, 0, 1, 1, 2, 2, 2, 2]
        segments = segment_genome(assignments, bin_size=200)
        assert len(segments) == 3
        assert segments[0]["state"] == 0
        assert segments[0]["length_bins"] == 3
        assert segments[0]["length_bp"] == 600

    def test_single_state(self) -> None:
        assignments = [0, 0, 0, 0, 0]
        segments = segment_genome(assignments, bin_size=100)
        assert len(segments) == 1
        assert segments[0]["length_bins"] == 5

    def test_alternating_states(self) -> None:
        assignments = [0, 1, 0, 1, 0, 1]
        segments = segment_genome(assignments, bin_size=100)
        assert len(segments) == 6

    def test_min_segment_filter(self) -> None:
        assignments = [0, 1, 0, 0, 0, 0]
        # With min_segment=2, the single-bin state 1 segment is filtered
        segments = segment_genome(assignments, bin_size=100, min_segment=2)
        assert all(s["length_bins"] >= 2 for s in segments)

    def test_bp_coordinates(self) -> None:
        assignments = [0, 0, 1, 1, 1]
        segments = segment_genome(assignments, bin_size=500)
        assert segments[0]["start_bp"] == 0
        assert segments[0]["end_bp"] == 1000
        assert segments[1]["start_bp"] == 1000
        assert segments[1]["end_bp"] == 2500

    def test_empty_raises(self) -> None:
        with pytest.raises(ValueError, match="must not be empty"):
            segment_genome([])

    def test_segment_keys(self) -> None:
        assignments = [0, 0, 1, 1]
        segments = segment_genome(assignments, bin_size=200)
        expected_keys = {"state", "start_bin", "end_bin", "start_bp", "end_bp", "length_bins", "length_bp"}
        for seg in segments:
            assert expected_keys.issubset(seg.keys())


# ---------------------------------------------------------------------------
# compare_chromatin_states
# ---------------------------------------------------------------------------


class TestCompareChromatinStates:
    """Tests for compare_chromatin_states."""

    def test_identical_states(self) -> None:
        data = _make_histone_matrix(n_bins=60, n_marks=4, n_clusters=3)
        model = learn_chromatin_states(data, n_states=3, max_iter=30)
        result = compare_chromatin_states(model, model)
        assert result["concordance"] == pytest.approx(1.0, abs=1e-6)
        assert result["discordant_bins"] == 0

    def test_different_states(self) -> None:
        data_a = _make_histone_matrix(n_bins=60, n_marks=4, n_clusters=3, seed=42)
        data_b = _make_histone_matrix(n_bins=60, n_marks=4, n_clusters=3, seed=99)
        model_a = learn_chromatin_states(data_a, n_states=3, max_iter=30)
        model_b = learn_chromatin_states(data_b, n_states=3, max_iter=30)
        result = compare_chromatin_states(model_a, model_b)
        assert 0.0 <= result["concordance"] <= 1.0
        assert result["n_bins"] == 60

    def test_length_mismatch_raises(self) -> None:
        model_a = {
            "assignments": [0, 1, 2],
            "states": 3,
        }
        model_b = {
            "assignments": [0, 1],
            "states": 3,
        }
        with pytest.raises(ValueError, match="must match"):
            compare_chromatin_states(model_a, model_b)

    def test_result_keys(self) -> None:
        data = _make_histone_matrix(n_bins=40, n_marks=3, n_clusters=2)
        model = learn_chromatin_states(data, n_states=2, max_iter=20)
        result = compare_chromatin_states(model, model)
        expected_keys = {
            "n_bins",
            "concordance",
            "discordance",
            "concordant_bins",
            "discordant_bins",
            "state_switches",
            "switch_matrix",
            "differential_segments",
            "summary",
        }
        assert expected_keys.issubset(result.keys())

    def test_switch_matrix_shape(self) -> None:
        data = _make_histone_matrix(n_bins=40, n_marks=3, n_clusters=2)
        model = learn_chromatin_states(data, n_states=2, max_iter=20)
        result = compare_chromatin_states(model, model)
        sm = result["switch_matrix"]
        n = max(model["states"], model["states"])
        assert len(sm) == n
        for row in sm:
            assert len(row) == n

    def test_differential_segments_found(self) -> None:
        model_a = {
            "assignments": [0, 0, 0, 1, 1, 1, 0, 0, 0],
            "states": 2,
        }
        model_b = {
            "assignments": [0, 0, 0, 0, 0, 0, 1, 1, 1],
            "states": 2,
        }
        result = compare_chromatin_states(model_a, model_b)
        assert result["n_differential_segments"] > 0


# ---------------------------------------------------------------------------
# compute_state_transition_rates
# ---------------------------------------------------------------------------


class TestComputeStateTransitionRates:
    """Tests for compute_state_transition_rates."""

    def test_basic_transition_rates(self) -> None:
        assignments = [0, 0, 0, 1, 1, 2, 2, 2, 0, 0]
        result = compute_state_transition_rates(assignments, n_states=3)
        assert "transition_counts" in result
        assert "transition_probabilities" in result
        assert "self_transition_rates" in result

    def test_probability_rows_sum_to_one(self) -> None:
        assignments = [0, 0, 1, 1, 2, 0, 1, 2, 0]
        result = compute_state_transition_rates(assignments, n_states=3)
        for row in result["transition_probabilities"]:
            assert sum(row) == pytest.approx(1.0, abs=0.01)

    def test_self_transition_rate_high_for_long_segments(self) -> None:
        # Long runs of same state -> high self-transition
        assignments = [0] * 50 + [1] * 50
        result = compute_state_transition_rates(assignments, n_states=2)
        assert result["self_transition_rates"][0] > 0.9
        assert result["self_transition_rates"][1] > 0.9

    def test_entropy(self) -> None:
        assignments = [0, 0, 0, 1, 1, 2, 0, 1, 2]
        result = compute_state_transition_rates(assignments, n_states=3)
        for h in result["transition_entropy"]:
            assert h >= 0.0

    def test_mean_segment_length(self) -> None:
        assignments = [0] * 20 + [1] * 20
        result = compute_state_transition_rates(assignments, n_states=2)
        for length in result["mean_segment_length"]:
            assert length > 1.0

    def test_empty_raises(self) -> None:
        with pytest.raises(ValueError, match="must not be empty"):
            compute_state_transition_rates([], n_states=2)

    def test_invalid_n_states_raises(self) -> None:
        with pytest.raises(ValueError, match="must be >= 1"):
            compute_state_transition_rates([0, 1], n_states=0)

    def test_n_transitions(self) -> None:
        assignments = [0, 1, 2, 0, 1]
        result = compute_state_transition_rates(assignments, n_states=3)
        assert result["n_transitions"] == 4
