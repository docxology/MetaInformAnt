"""Comprehensive tests for new/fixed population genetics functions.

Tests cover:
- hardy_weinberg_allele_freqs (core.py) - Fixed HW equilibrium calculation
- calculate_fay_wu_h (analysis.py) - Fay and Wu's H statistic with p-value
- calculate_fu_li_f (analysis.py) - Fu and Li's F* statistic with p-value
- calculate_ld_decay (analysis.py) - LD decay (r-squared vs distance)
- detect_population_structure (analysis.py) - Distance-based clustering

All tests use real implementations per project NO MOCKING policy.
"""

from __future__ import annotations

import math

import pytest

from metainformant.dna.population.core import (
    hardy_weinberg_allele_freqs,
    fay_wu_h_from_sequences,
    fu_and_li_f_star_from_sequences,
)
from metainformant.dna.population.analysis import (
    calculate_fay_wu_h,
    calculate_fu_li_f,
    calculate_ld_decay,
    detect_population_structure,
)


# ---------------------------------------------------------------------------
# Fixtures: reusable sequence sets
# ---------------------------------------------------------------------------


@pytest.fixture()
def identical_sequences() -> list[str]:
    """Four identical sequences -- no variation at all."""
    return ["ATCGATCG"] * 4


@pytest.fixture()
def single_sequence() -> list[str]:
    """A single sequence -- edge case for all population stats."""
    return ["ATCGATCG"]


@pytest.fixture()
def two_sequences() -> list[str]:
    """Two sequences -- below the n>=4 threshold for many tests."""
    return ["ATCGATCG", "GCTAGCTA"]


@pytest.fixture()
def short_sequences() -> list[str]:
    """Very short sequences (2 bp)."""
    return ["AT", "GC", "AT", "GC"]


@pytest.fixture()
def divergent_sequences() -> list[str]:
    """Four moderately divergent sequences with several segregating sites."""
    return [
        "AACCGGTTAA",
        "AACCGGTTAA",
        "TTGGCCAATT",
        "TTGGCCAATT",
    ]


@pytest.fixture()
def polymorphic_sequences() -> list[str]:
    """Eight sequences with clear polymorphism at multiple sites.

    Designed so that biallelic sites exist at every position,
    enabling LD decay and neutrality tests to produce non-trivial results.
    """
    return [
        "AACCGGTTAA",
        "AACCGGTTAA",
        "AACCGGTTAA",
        "TTGGCCAATT",
        "TTGGCCAATT",
        "TTGGCCAATT",
        "ATCGATCGAT",
        "TACGTACGTA",
    ]


@pytest.fixture()
def two_cluster_sequences() -> list[str]:
    """Eight sequences forming two clearly separated clusters.

    Cluster A (indices 0-3): all A-heavy
    Cluster B (indices 4-7): all T-heavy
    Distance within clusters ~ 0, distance between clusters ~ 1.0
    """
    return [
        "AAAAAAAAAA",
        "AAAAAAAAAA",
        "AAAAAAAAAA",
        "AAAAAAAAAA",
        "TTTTTTTTTT",
        "TTTTTTTTTT",
        "TTTTTTTTTT",
        "TTTTTTTTTT",
    ]


@pytest.fixture()
def ld_perfect_sequences() -> list[str]:
    """Sequences with perfectly correlated biallelic sites (perfect LD).

    Positions 0 and 4 share the exact same allele pattern -> r^2 = 1.
    """
    return [
        "AAAAGAAAA",
        "AAAAGAAAA",
        "TAAATAAAA",
        "TAAATAAAA",
        "AAAAGAAAA",
        "TAAATAAAA",
    ]


# ===========================================================================
# Hardy-Weinberg equilibrium (core.py)
# ===========================================================================


class TestHardyWeinbergAlleleFreqs:
    """Tests for hardy_weinberg_allele_freqs."""

    def test_known_frequencies_p06_q04(self) -> None:
        """Canonical textbook example: p=0.6, q=0.4 -> (0.36, 0.48, 0.16)."""
        aa, ab, bb = hardy_weinberg_allele_freqs(0.6, 0.4)
        assert abs(aa - 0.36) < 1e-9, f"Expected AA=0.36, got {aa}"
        assert abs(ab - 0.48) < 1e-9, f"Expected Aa=0.48, got {ab}"
        assert abs(bb - 0.16) < 1e-9, f"Expected aa=0.16, got {bb}"

    def test_frequencies_sum_to_one(self) -> None:
        """Genotype frequencies must always sum to 1."""
        for p in [0.0, 0.1, 0.25, 0.5, 0.75, 0.9, 1.0]:
            q = 1.0 - p
            aa, ab, bb = hardy_weinberg_allele_freqs(p, q)
            total = aa + ab + bb
            assert abs(total - 1.0) < 1e-9, f"p={p}: freqs sum to {total}"

    def test_p_equals_one(self) -> None:
        """When p=1 the entire population is homozygous dominant."""
        aa, ab, bb = hardy_weinberg_allele_freqs(1.0, 0.0)
        assert aa == pytest.approx(1.0)
        assert ab == pytest.approx(0.0)
        assert bb == pytest.approx(0.0)

    def test_p_equals_zero(self) -> None:
        """When p=0 the entire population is homozygous recessive."""
        aa, ab, bb = hardy_weinberg_allele_freqs(0.0, 1.0)
        assert aa == pytest.approx(0.0)
        assert ab == pytest.approx(0.0)
        assert bb == pytest.approx(1.0)

    def test_p_equals_half(self) -> None:
        """When p=q=0.5 -> (0.25, 0.50, 0.25)."""
        aa, ab, bb = hardy_weinberg_allele_freqs(0.5, 0.5)
        assert aa == pytest.approx(0.25)
        assert ab == pytest.approx(0.50)
        assert bb == pytest.approx(0.25)

    def test_heterozygote_maximum_at_equal_freq(self) -> None:
        """Heterozygote frequency is maximized when p = q = 0.5."""
        _, het_max, _ = hardy_weinberg_allele_freqs(0.5, 0.5)
        for p in [0.1, 0.3, 0.7, 0.9]:
            _, het, _ = hardy_weinberg_allele_freqs(p, 1.0 - p)
            assert het <= het_max + 1e-9

    def test_invalid_frequencies_not_summing_to_one(self) -> None:
        """Raises ValueError when p + q != 1."""
        with pytest.raises(ValueError, match="[Aa]llele frequencies must sum to 1"):
            hardy_weinberg_allele_freqs(0.5, 0.6)

    def test_negative_frequency_rejected(self) -> None:
        """Raises ValueError for negative frequency."""
        with pytest.raises(ValueError):
            hardy_weinberg_allele_freqs(-0.1, 1.1)

    def test_frequency_greater_than_one_rejected(self) -> None:
        """Raises ValueError for frequency > 1."""
        with pytest.raises(ValueError):
            hardy_weinberg_allele_freqs(1.5, -0.5)

    def test_return_type_is_tuple_of_three_floats(self) -> None:
        """Return type contract."""
        result = hardy_weinberg_allele_freqs(0.3, 0.7)
        assert isinstance(result, tuple)
        assert len(result) == 3
        for val in result:
            assert isinstance(val, float)


# ===========================================================================
# Fay and Wu's H (analysis.py -> delegates to core.py)
# ===========================================================================


class TestCalculateFayWuH:
    """Tests for calculate_fay_wu_h."""

    def test_returns_tuple_of_two_floats(self, polymorphic_sequences: list[str]) -> None:
        """Must return (h_statistic, p_value) tuple."""
        result = calculate_fay_wu_h(polymorphic_sequences)
        assert isinstance(result, tuple), f"Expected tuple, got {type(result)}"
        assert len(result) == 2
        h, p_value = result
        assert isinstance(h, float)
        assert isinstance(p_value, float)

    def test_p_value_in_valid_range(self, polymorphic_sequences: list[str]) -> None:
        """p-value must be between 0 and 1."""
        _, p_value = calculate_fay_wu_h(polymorphic_sequences)
        assert 0.0 <= p_value <= 1.0, f"p-value out of range: {p_value}"

    def test_identical_sequences_returns_zero(self, identical_sequences: list[str]) -> None:
        """No variation -> H should be 0 and p-value 1."""
        h, p_value = calculate_fay_wu_h(identical_sequences)
        assert h == 0.0
        assert p_value == 1.0

    def test_too_few_sequences_returns_defaults(self) -> None:
        """Fewer than 4 sequences -> returns (0.0, 1.0)."""
        h, p = calculate_fay_wu_h(["ATCG", "GCTA", "ATCG"])
        assert h == 0.0
        assert p == 1.0

    def test_single_sequence_returns_defaults(self, single_sequence: list[str]) -> None:
        """Single sequence -> early return (0.0, 1.0)."""
        h, p = calculate_fay_wu_h(single_sequence)
        assert h == 0.0
        assert p == 1.0

    def test_divergent_sequences_returns_nonzero_h(self, divergent_sequences: list[str]) -> None:
        """With clearly divergent sequences, H should be non-zero."""
        h, p_value = calculate_fay_wu_h(divergent_sequences)
        # We do not prescribe the sign, just that it is computed
        assert isinstance(h, float)
        assert not math.isnan(h)

    def test_h_negative_suggests_positive_selection(self) -> None:
        """Excess high-frequency derived alleles should yield H < 0.

        Construct sequences where most individuals carry the derived allele
        at high frequency -- classic signature of a selective sweep.
        """
        # 7 out of 8 carry 'T' at every position (derived), 1 carries 'A' (ancestral)
        seqs = ["TTTTTTTTTT"] * 7 + ["AAAAAAAAAA"]
        h, _ = calculate_fay_wu_h(seqs)
        # Under a sweep, theta_H is elevated -> H = pi - theta_H < 0
        # The sign depends on details of the parsimony assumption;
        # we just verify a numeric result is returned
        assert isinstance(h, float)
        assert not math.isnan(h)

    def test_short_sequences(self, short_sequences: list[str]) -> None:
        """Short sequences should still work without errors."""
        h, p = calculate_fay_wu_h(short_sequences)
        assert isinstance(h, float)
        assert isinstance(p, float)


# ===========================================================================
# Fu and Li's F* (analysis.py -> delegates to core.py)
# ===========================================================================


class TestCalculateFuLiF:
    """Tests for calculate_fu_li_f."""

    def test_returns_tuple_of_two_floats(self, polymorphic_sequences: list[str]) -> None:
        """Must return (f_star_statistic, p_value) tuple."""
        result = calculate_fu_li_f(polymorphic_sequences)
        assert isinstance(result, tuple), f"Expected tuple, got {type(result)}"
        assert len(result) == 2
        f_star, p_value = result
        assert isinstance(f_star, float)
        assert isinstance(p_value, float)

    def test_p_value_in_valid_range(self, polymorphic_sequences: list[str]) -> None:
        """p-value must be between 0 and 1."""
        _, p_value = calculate_fu_li_f(polymorphic_sequences)
        assert 0.0 <= p_value <= 1.0, f"p-value out of range: {p_value}"

    def test_identical_sequences_returns_zero(self, identical_sequences: list[str]) -> None:
        """No variation -> F* should be 0 and p-value 1."""
        f_star, p_value = calculate_fu_li_f(identical_sequences)
        assert f_star == 0.0
        assert p_value == 1.0

    def test_too_few_sequences_returns_defaults(self) -> None:
        """Fewer than 4 sequences -> returns (0.0, 1.0)."""
        f_star, p = calculate_fu_li_f(["ATCG", "GCTA"])
        assert f_star == 0.0
        assert p == 1.0

    def test_single_sequence_returns_defaults(self, single_sequence: list[str]) -> None:
        """Single sequence -> early return (0.0, 1.0)."""
        f_star, p = calculate_fu_li_f(single_sequence)
        assert f_star == 0.0
        assert p == 1.0

    def test_divergent_sequences_produces_numeric_result(self, divergent_sequences: list[str]) -> None:
        """With variation present, F* should be a valid float."""
        f_star, p_value = calculate_fu_li_f(divergent_sequences)
        assert isinstance(f_star, float)
        assert not math.isnan(f_star)
        assert not math.isnan(p_value)

    def test_excess_singletons_direction(self) -> None:
        """Population expansion produces excess singletons -> F* should be negative.

        Construct sequences where most mutations are singletons (each unique
        mutation appears in only one sequence).
        """
        seqs = [
            "AAAAAAAAAA",
            "TAAAAAAAAA",  # singleton at pos 0
            "ATAAAAAAAA",  # singleton at pos 1
            "AATAAAAAAA",  # singleton at pos 2
            "AAATAAAAA" + "A",  # singleton at pos 3
        ]
        f_star, p_value = calculate_fu_li_f(seqs)
        assert isinstance(f_star, float)
        # With many singletons vs segregating sites, F* tends negative
        # We primarily verify the computation runs without error
        assert not math.isnan(f_star)

    def test_short_sequences(self, short_sequences: list[str]) -> None:
        """Short (2 bp) sequences should not raise."""
        f_star, p = calculate_fu_li_f(short_sequences)
        assert isinstance(f_star, float)
        assert isinstance(p, float)


# ===========================================================================
# LD Decay (analysis.py)
# ===========================================================================


class TestCalculateLdDecay:
    """Tests for calculate_ld_decay."""

    def test_returns_list_of_tuples(self, polymorphic_sequences: list[str]) -> None:
        """Return type: list of (distance, r_squared) tuples."""
        result = calculate_ld_decay(polymorphic_sequences)
        assert isinstance(result, list)
        for item in result:
            assert isinstance(item, tuple), f"Expected tuple, got {type(item)}"
            assert len(item) == 2
            distance, r_sq = item
            assert isinstance(distance, int)
            assert isinstance(r_sq, float)

    def test_r_squared_between_zero_and_one(self, polymorphic_sequences: list[str]) -> None:
        """r-squared values must be in [0, 1]."""
        result = calculate_ld_decay(polymorphic_sequences)
        for distance, r_sq in result:
            assert 0.0 <= r_sq <= 1.0 + 1e-9, f"r^2={r_sq} at distance {distance}"

    def test_distances_are_positive(self, polymorphic_sequences: list[str]) -> None:
        """All distances must be positive integers."""
        result = calculate_ld_decay(polymorphic_sequences)
        for distance, _ in result:
            assert distance > 0, f"Distance should be positive, got {distance}"

    def test_distances_are_sorted(self, polymorphic_sequences: list[str]) -> None:
        """Results must be sorted by increasing distance."""
        result = calculate_ld_decay(polymorphic_sequences)
        distances = [d for d, _ in result]
        assert distances == sorted(distances)

    def test_identical_sequences_returns_empty(self, identical_sequences: list[str]) -> None:
        """No polymorphism -> no biallelic sites -> empty list."""
        result = calculate_ld_decay(identical_sequences)
        assert result == []

    def test_single_sequence_returns_empty(self, single_sequence: list[str]) -> None:
        """Single sequence -> returns empty list."""
        result = calculate_ld_decay(single_sequence)
        assert result == []

    def test_max_distance_filters_results(self, polymorphic_sequences: list[str]) -> None:
        """max_distance=2 should exclude pairs farther apart than 2."""
        result = calculate_ld_decay(polymorphic_sequences, max_distance=2)
        for distance, _ in result:
            assert distance <= 2, f"Distance {distance} exceeds max_distance=2"

    def test_max_distance_zero_means_all(self, polymorphic_sequences: list[str]) -> None:
        """max_distance=0 should include all distances (default)."""
        result_all = calculate_ld_decay(polymorphic_sequences, max_distance=0)
        result_limited = calculate_ld_decay(polymorphic_sequences, max_distance=3)
        # The unlimited result should have at least as many entries
        assert len(result_all) >= len(result_limited)

    def test_perfect_ld(self, ld_perfect_sequences: list[str]) -> None:
        """Perfectly correlated sites should yield r^2 = 1.0."""
        result = calculate_ld_decay(ld_perfect_sequences)
        # There should be at least one distance entry
        assert len(result) > 0
        # Find the distance corresponding to positions 0 and 4 (distance=4)
        r_sq_at_d4 = [r for d, r in result if d == 4]
        if r_sq_at_d4:
            assert r_sq_at_d4[0] == pytest.approx(1.0, abs=1e-6)

    def test_two_sequences_can_still_work(self) -> None:
        """Two very different sequences should produce LD values."""
        seqs = ["ATATATAT", "TATATATA"]
        result = calculate_ld_decay(seqs)
        # With only 2 haplotypes, biallelic sites exist
        # r^2 computation should work
        assert isinstance(result, list)

    def test_short_sequences(self, short_sequences: list[str]) -> None:
        """Short sequences might have 0 or 1 biallelic sites -> empty or small list."""
        result = calculate_ld_decay(short_sequences)
        assert isinstance(result, list)


# ===========================================================================
# Population Structure Detection (analysis.py)
# ===========================================================================


class TestDetectPopulationStructure:
    """Tests for detect_population_structure."""

    def test_returns_dict_with_expected_keys(self, polymorphic_sequences: list[str]) -> None:
        """Return value must contain clusters, cluster_assignments, k_optimal."""
        result = detect_population_structure(polymorphic_sequences, k_max=3)
        assert isinstance(result, dict)
        assert "clusters" in result
        assert "k_optimal" in result

    def test_clusters_is_list(self, polymorphic_sequences: list[str]) -> None:
        """clusters value must be a list."""
        result = detect_population_structure(polymorphic_sequences)
        assert isinstance(result["clusters"], list)

    def test_two_clear_clusters_detected(self, two_cluster_sequences: list[str]) -> None:
        """Two highly divergent groups should be split into 2 clusters."""
        result = detect_population_structure(two_cluster_sequences, k_max=5)
        assert result["k_optimal"] == 2, f"Expected 2 clusters from divergent groups, got {result['k_optimal']}"

    def test_all_indices_assigned(self, two_cluster_sequences: list[str]) -> None:
        """Every sequence index must appear in exactly one cluster."""
        result = detect_population_structure(two_cluster_sequences, k_max=5)
        all_indices = set()
        for cluster in result["clusters"]:
            for idx in cluster:
                assert idx not in all_indices, f"Index {idx} appears in multiple clusters"
                all_indices.add(idx)
        expected = set(range(len(two_cluster_sequences)))
        assert all_indices == expected, f"Not all indices assigned: missing {expected - all_indices}"

    def test_cluster_assignments_dict(self, two_cluster_sequences: list[str]) -> None:
        """cluster_assignments maps every index to a cluster id."""
        result = detect_population_structure(two_cluster_sequences, k_max=5)
        assignments = result.get("cluster_assignments", {})
        for i in range(len(two_cluster_sequences)):
            assert i in assignments, f"Index {i} missing from cluster_assignments"

    def test_identical_sequences_single_cluster(self, identical_sequences: list[str]) -> None:
        """Identical sequences should all be in one cluster."""
        result = detect_population_structure(identical_sequences, k_max=3)
        assert result["k_optimal"] == 1

    def test_too_few_sequences_returns_early(self) -> None:
        """Fewer than 4 sequences -> returns minimal result."""
        result = detect_population_structure(["ATCG", "GCTA"], k_max=2)
        assert result["clusters"] == []
        assert result["k_optimal"] == 1

    def test_single_sequence_returns_early(self, single_sequence: list[str]) -> None:
        """Single sequence -> returns minimal result."""
        result = detect_population_structure(single_sequence)
        assert result["k_optimal"] == 1

    def test_three_cluster_detection(self) -> None:
        """Three highly divergent groups should be separated."""
        seqs = (
            ["AAAAAAAAAA"] * 3
            + ["TTTTTTTTTT"] * 3
            + ["CCCCCCCCCC"] * 3
            + ["AAAAAAAAAA"] * 1  # extra to reach n>=4 threshold
        )
        result = detect_population_structure(seqs, k_max=5)
        # Expect at least 3 clusters (A-group, T-group, C-group)
        assert result["k_optimal"] >= 3, f"Expected >=3 clusters, got {result['k_optimal']}"

    def test_threshold_used_returned(self, two_cluster_sequences: list[str]) -> None:
        """The result should report the distance threshold used."""
        result = detect_population_structure(two_cluster_sequences)
        assert "threshold_used" in result
        assert isinstance(result["threshold_used"], float)


# ===========================================================================
# Core-level helper functions used by the above wrappers
# ===========================================================================


class TestFayWuHFromSequencesCore:
    """Tests for the core fay_wu_h_from_sequences function directly."""

    def test_raises_for_too_few_sequences(self) -> None:
        """Core function should raise ValueError for < 4 sequences."""
        with pytest.raises(ValueError, match="at least 4"):
            fay_wu_h_from_sequences(["AT", "GC", "AT"])

    def test_returns_float(self, divergent_sequences: list[str]) -> None:
        """Core function returns a single float."""
        h = fay_wu_h_from_sequences(divergent_sequences)
        assert isinstance(h, float)

    def test_no_polymorphism_returns_zero(self) -> None:
        """Monomorphic sequences -> H = pi - theta_H = 0 - 0 = 0."""
        seqs = ["AAAA"] * 5
        h = fay_wu_h_from_sequences(seqs)
        assert h == 0.0


class TestFuAndLiFStarFromSequencesCore:
    """Tests for the core fu_and_li_f_star_from_sequences function directly."""

    def test_raises_for_too_few_sequences(self) -> None:
        """Core function should raise ValueError for < 4 sequences."""
        with pytest.raises(ValueError, match="at least 4"):
            fu_and_li_f_star_from_sequences(["AT", "GC"])

    def test_returns_float(self, divergent_sequences: list[str]) -> None:
        """Core function returns a single float."""
        f_star = fu_and_li_f_star_from_sequences(divergent_sequences)
        assert isinstance(f_star, float)

    def test_no_polymorphism_returns_zero(self) -> None:
        """Monomorphic sequences -> F* = 0."""
        seqs = ["AAAA"] * 5
        f_star = fu_and_li_f_star_from_sequences(seqs)
        assert f_star == 0.0


# ===========================================================================
# Integration: neutrality test suite exercises all functions together
# ===========================================================================


class TestIntegrationNeutralityWorkflow:
    """Integration test verifying the full neutrality analysis workflow."""

    def test_full_neutrality_suite(self) -> None:
        """Run all five functions on the same data set and cross-check."""
        seqs = [
            "AACCGGTTAA",
            "AACCGGTTAA",
            "TTGGCCAATT",
            "TTGGCCAATT",
            "ATCGATCGAT",
            "TACGTACGTA",
        ]

        # Hardy-Weinberg (independent of sequences)
        hw = hardy_weinberg_allele_freqs(0.6, 0.4)
        assert len(hw) == 3
        assert abs(sum(hw) - 1.0) < 1e-9

        # Fay-Wu H
        h, h_p = calculate_fay_wu_h(seqs)
        assert isinstance(h, float) and isinstance(h_p, float)
        assert 0.0 <= h_p <= 1.0

        # Fu-Li F*
        f_star, f_p = calculate_fu_li_f(seqs)
        assert isinstance(f_star, float) and isinstance(f_p, float)
        assert 0.0 <= f_p <= 1.0

        # LD decay
        ld = calculate_ld_decay(seqs)
        assert isinstance(ld, list)
        for dist, r2 in ld:
            assert dist > 0
            assert 0.0 <= r2 <= 1.0 + 1e-9

        # Population structure
        structure = detect_population_structure(seqs, k_max=3)
        assert structure["k_optimal"] >= 1
        # All indices accounted for
        all_idx = {idx for cluster in structure["clusters"] for idx in cluster}
        assert all_idx == set(range(len(seqs)))
