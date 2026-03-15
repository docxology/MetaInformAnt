"""Comprehensive tests for the k-mer analysis module.

Tests cover all 12 public functions in ``metainformant.dna.sequence.kmer``
with normal cases, edge cases (empty sequence, k > length, single base),
and mathematically verified expected values.  NO mocking -- real
implementations only.
"""

from __future__ import annotations

import math

import pytest

from metainformant.dna.sequence.kmer import (
    compare_kmer_profiles,
    count_kmers,
    find_homopolymers,
    find_microsatellites,
    find_overrepresented_kmers,
    find_unique_kmers,
    kmer_diversity,
    kmer_frequencies,
    kmer_spectrum,
    mask_low_complexity,
    sequence_complexity_profile,
    sliding_window_kmer_complexity,
)

# ===================================================================
# 1. count_kmers
# ===================================================================


class TestCountKmers:
    """Tests for count_kmers."""

    def test_basic_dinucleotides(self) -> None:
        """Docstring example: ATCGATCG with k=2."""
        result = count_kmers("ATCGATCG", 2)
        assert result == {"AT": 2, "TC": 2, "CG": 2, "GA": 1}

    def test_mononucleotides(self) -> None:
        """k=1 counts individual bases."""
        result = count_kmers("AATTCCGG", 1)
        assert result == {"A": 2, "T": 2, "C": 2, "G": 2}

    def test_k_equals_sequence_length(self) -> None:
        """When k == len(seq), only one k-mer is possible."""
        result = count_kmers("ATCG", 4)
        assert result == {"ATCG": 1}

    def test_homopolymer_counting(self) -> None:
        """A homopolymer yields a single k-mer with count n-k+1."""
        result = count_kmers("AAAAAA", 3)
        assert result == {"AAA": 4}

    def test_case_insensitivity(self) -> None:
        """Lower-case input is upper-cased before counting."""
        result = count_kmers("atcg", 2)
        assert result == {"AT": 1, "TC": 1, "CG": 1}

    def test_mixed_case(self) -> None:
        """Mixed-case input merges into upper-case k-mers."""
        result = count_kmers("AtCgAtCg", 2)
        assert result == {"AT": 2, "TC": 2, "CG": 2, "GA": 1}

    def test_empty_sequence(self) -> None:
        """Empty input returns empty dict."""
        assert count_kmers("", 2) == {}

    def test_k_zero(self) -> None:
        """k=0 returns empty dict."""
        assert count_kmers("ATCG", 0) == {}

    def test_k_negative(self) -> None:
        """Negative k returns empty dict."""
        assert count_kmers("ATCG", -1) == {}

    def test_k_greater_than_length(self) -> None:
        """k > len(seq) returns empty dict."""
        assert count_kmers("AT", 5) == {}

    def test_single_base(self) -> None:
        """Single-base sequence with k=1."""
        result = count_kmers("A", 1)
        assert result == {"A": 1}

    def test_total_count_invariant(self) -> None:
        """Total count of k-mers must equal len(seq) - k + 1."""
        seq = "ATCGATCGATCG"
        k = 3
        result = count_kmers(seq, k)
        assert sum(result.values()) == len(seq) - k + 1


# ===================================================================
# 2. kmer_frequencies
# ===================================================================


class TestKmerFrequencies:
    """Tests for kmer_frequencies."""

    def test_homopolymer_all_same(self) -> None:
        """Homopolymer has frequency 1.0 for the single k-mer."""
        result = kmer_frequencies("AAAA", 2)
        assert result == {"AA": 1.0}

    def test_frequencies_sum_to_one(self) -> None:
        """All frequencies must sum to 1.0."""
        result = kmer_frequencies("ATCGATCG", 2)
        total = sum(result.values())
        assert abs(total - 1.0) < 1e-10

    def test_known_values(self) -> None:
        """ATCGATCG k=2: AT=2/7, TC=2/7, CG=2/7, GA=1/7."""
        result = kmer_frequencies("ATCGATCG", 2)
        assert abs(result["AT"] - 2 / 7) < 1e-10
        assert abs(result["TC"] - 2 / 7) < 1e-10
        assert abs(result["CG"] - 2 / 7) < 1e-10
        assert abs(result["GA"] - 1 / 7) < 1e-10

    def test_empty_sequence(self) -> None:
        assert kmer_frequencies("", 3) == {}

    def test_k_exceeds_length(self) -> None:
        assert kmer_frequencies("AT", 10) == {}

    def test_single_kmer_window(self) -> None:
        """k == len(seq) gives a single k-mer with frequency 1.0."""
        result = kmer_frequencies("ATCG", 4)
        assert result == {"ATCG": 1.0}


# ===================================================================
# 3. kmer_spectrum
# ===================================================================


class TestKmerSpectrum:
    """Tests for kmer_spectrum."""

    def test_docstring_example(self) -> None:
        """ATCGATCG k=2: AT(2), TC(2), CG(2), GA(1) => {2: 3, 1: 1}."""
        result = kmer_spectrum("ATCGATCG", 2)
        assert result == {2: 3, 1: 1}

    def test_homopolymer_single_frequency(self) -> None:
        """All k-mers have same count => spectrum has one entry."""
        result = kmer_spectrum("AAAA", 2)
        # Only k-mer "AA" with count 3
        assert result == {3: 1}

    def test_all_unique_kmers(self) -> None:
        """When every k-mer is unique, spectrum is {1: <number of distinct k-mers>}."""
        result = kmer_spectrum("ATCG", 2)
        # AT, TC, CG each appear once
        assert result == {1: 3}

    def test_empty_sequence(self) -> None:
        assert kmer_spectrum("", 2) == {}

    def test_k_exceeds_length(self) -> None:
        assert kmer_spectrum("AT", 10) == {}

    def test_spectrum_values_sum_to_distinct_kmers(self) -> None:
        """Sum of spectrum values equals number of distinct k-mers."""
        seq = "ATCGATCGATCGAAAA"
        k = 3
        spectrum = kmer_spectrum(seq, k)
        counts = count_kmers(seq, k)
        assert sum(spectrum.values()) == len(counts)


# ===================================================================
# 4. find_unique_kmers
# ===================================================================


class TestFindUniqueKmers:
    """Tests for find_unique_kmers."""

    def test_docstring_example(self) -> None:
        """ATCGATCG k=2: only GA appears once."""
        result = find_unique_kmers("ATCGATCG", 2)
        assert result == ["GA"]

    def test_all_unique(self) -> None:
        """Short sequence where every k-mer is unique."""
        result = find_unique_kmers("ATCG", 2)
        assert sorted(result) == ["AT", "CG", "TC"]

    def test_no_unique_kmers(self) -> None:
        """Homopolymer has no unique k-mers (all repeated)."""
        result = find_unique_kmers("AAAA", 2)
        # Only "AA" appears 3 times
        assert result == []

    def test_sorted_output(self) -> None:
        """Output is sorted alphabetically."""
        result = find_unique_kmers("GCATATGC", 2)
        assert result == sorted(result)

    def test_empty_sequence(self) -> None:
        assert find_unique_kmers("", 2) == []

    def test_k_exceeds_length(self) -> None:
        assert find_unique_kmers("AT", 5) == []

    def test_single_base_k1(self) -> None:
        """Single base with k=1 is unique."""
        result = find_unique_kmers("A", 1)
        assert result == ["A"]


# ===================================================================
# 5. find_overrepresented_kmers
# ===================================================================


class TestFindOverrepresentedKmers:
    """Tests for find_overrepresented_kmers."""

    def test_docstring_example(self) -> None:
        """AAAAAA k=2: freq(AA)=1.0, expected=1/16, fold=16.0."""
        result = find_overrepresented_kmers("AAAAAA", 2, threshold=1.5)
        assert "AA" in result
        assert abs(result["AA"] - 16.0) < 1e-10

    def test_threshold_filtering(self) -> None:
        """Only k-mers at or above the threshold are returned."""
        result = find_overrepresented_kmers("AAAAAA", 2, threshold=20.0)
        # fold is 16.0, so nothing passes threshold=20
        assert result == {}

    def test_uniform_sequence_below_threshold(self) -> None:
        """A sequence with roughly uniform k-mer usage and low threshold."""
        # ATCGATCG has 4 distinct 2-mers; expected for k=2 is 1/16=0.0625
        # AT: freq=2/7=0.2857, fold=0.2857/0.0625=4.57
        result = find_overrepresented_kmers("ATCGATCG", 2, threshold=5.0)
        # All folds are ~4.57 or ~2.28 so nothing passes 5.0
        assert "GA" not in result  # GA fold = (1/7)/0.0625 = 2.29

    def test_sorted_descending_by_fold(self) -> None:
        """Results are sorted descending by fold-enrichment."""
        result = find_overrepresented_kmers("AAACCC", 2, threshold=1.0)
        folds = list(result.values())
        assert folds == sorted(folds, reverse=True)

    def test_empty_sequence(self) -> None:
        assert find_overrepresented_kmers("", 2) == {}

    def test_k_exceeds_length(self) -> None:
        assert find_overrepresented_kmers("AT", 10) == {}

    def test_known_fold_computation(self) -> None:
        """Verify fold = observed_freq / expected for known values."""
        seq = "ATCGATCG"
        k = 2
        result = find_overrepresented_kmers(seq, k, threshold=0.0)
        expected_freq = 1.0 / (4**k)  # 1/16 = 0.0625
        for kmer, fold in result.items():
            freqs = kmer_frequencies(seq, k)
            assert abs(fold - freqs[kmer] / expected_freq) < 1e-10


# ===================================================================
# 6. kmer_diversity
# ===================================================================


class TestKmerDiversity:
    """Tests for kmer_diversity."""

    def test_docstring_example(self) -> None:
        """ATCGATCG k=2: entropy approx 1.9502 bits."""
        result = kmer_diversity("ATCGATCG", 2)
        assert abs(result - 1.9502) < 0.001

    def test_exact_entropy_calculation(self) -> None:
        """Manually compute Shannon entropy for ATCGATCG k=2.

        freqs: AT=2/7, TC=2/7, CG=2/7, GA=1/7
        H = -3*(2/7)*log2(2/7) - (1/7)*log2(1/7)
        """
        p1 = 2.0 / 7.0
        p2 = 1.0 / 7.0
        expected = -(3 * p1 * math.log2(p1) + p2 * math.log2(p2))
        result = kmer_diversity("ATCGATCG", 2)
        assert abs(result - expected) < 1e-10

    def test_homopolymer_zero_entropy(self) -> None:
        """A homopolymer k-mer distribution has entropy 0 (only one k-mer type)."""
        # "AAAA" k=2: only "AA" => freq=1.0 => H = -1.0*log2(1.0) = 0.0
        result = kmer_diversity("AAAA", 2)
        assert result == 0.0

    def test_maximum_entropy_for_uniform(self) -> None:
        """When all k-mers are equally frequent, entropy is log2(n_distinct)."""
        # "ATCG" k=1: A,T,C,G each freq=0.25 => H = log2(4) = 2.0
        result = kmer_diversity("ATCG", 1)
        assert abs(result - 2.0) < 1e-10

    def test_empty_sequence(self) -> None:
        assert kmer_diversity("", 2) == 0.0

    def test_k_exceeds_length(self) -> None:
        assert kmer_diversity("AT", 10) == 0.0

    def test_entropy_non_negative(self) -> None:
        """Entropy must always be >= 0."""
        for seq in ["A", "ATCG", "AATTCCGG", "AAACCCTTTTGGG"]:
            for k in [1, 2, 3]:
                assert kmer_diversity(seq, k) >= 0.0


# ===================================================================
# 7. compare_kmer_profiles
# ===================================================================


class TestCompareKmerProfiles:
    """Tests for compare_kmer_profiles."""

    def test_identical_sequences(self) -> None:
        """Identical sequences: cosine=1, jaccard=1, bray_curtis=0."""
        result = compare_kmer_profiles("ATCGATCG", "ATCGATCG", 2)
        assert abs(result["cosine"] - 1.0) < 1e-10
        assert abs(result["jaccard"] - 1.0) < 1e-10
        assert abs(result["bray_curtis"] - 0.0) < 1e-10

    def test_completely_different_kmers(self) -> None:
        """No shared k-mers => jaccard=0, bray_curtis=1."""
        # "AAAA" k=2 => {"AA"}, "CCCC" k=2 => {"CC"} -- disjoint sets
        result = compare_kmer_profiles("AAAA", "CCCC", 2)
        assert abs(result["jaccard"] - 0.0) < 1e-10
        assert abs(result["bray_curtis"] - 1.0) < 1e-10
        # Cosine of orthogonal vectors is 0
        assert abs(result["cosine"] - 0.0) < 1e-10

    def test_both_empty(self) -> None:
        """Both empty sequences: treated as identical."""
        result = compare_kmer_profiles("", "", 2)
        assert result == {"cosine": 1.0, "jaccard": 1.0, "bray_curtis": 0.0}

    def test_one_empty(self) -> None:
        """One empty, one non-empty: maximally dissimilar."""
        result = compare_kmer_profiles("ATCG", "", 2)
        assert result == {"cosine": 0.0, "jaccard": 0.0, "bray_curtis": 1.0}

    def test_cosine_symmetry(self) -> None:
        """cosine(a, b) == cosine(b, a)."""
        r1 = compare_kmer_profiles("ATCGATCG", "AAACCC", 2)
        r2 = compare_kmer_profiles("AAACCC", "ATCGATCG", 2)
        assert abs(r1["cosine"] - r2["cosine"]) < 1e-10

    def test_jaccard_symmetry(self) -> None:
        """jaccard(a, b) == jaccard(b, a)."""
        r1 = compare_kmer_profiles("ATCGATCG", "AAACCC", 2)
        r2 = compare_kmer_profiles("AAACCC", "ATCGATCG", 2)
        assert abs(r1["jaccard"] - r2["jaccard"]) < 1e-10

    def test_bray_curtis_symmetry(self) -> None:
        """bray_curtis(a, b) == bray_curtis(b, a)."""
        r1 = compare_kmer_profiles("ATCGATCG", "AAACCC", 2)
        r2 = compare_kmer_profiles("AAACCC", "ATCGATCG", 2)
        assert abs(r1["bray_curtis"] - r2["bray_curtis"]) < 1e-10

    def test_cosine_range(self) -> None:
        """Cosine similarity should be in [0, 1] for non-negative vectors."""
        result = compare_kmer_profiles("ATCGATCG", "AATTCCGG", 2)
        assert 0.0 <= result["cosine"] <= 1.0

    def test_jaccard_range(self) -> None:
        """Jaccard index should be in [0, 1]."""
        result = compare_kmer_profiles("ATCGATCG", "AATTCCGG", 2)
        assert 0.0 <= result["jaccard"] <= 1.0

    def test_bray_curtis_range(self) -> None:
        """Bray-Curtis should be in [0, 1]."""
        result = compare_kmer_profiles("ATCGATCG", "AATTCCGG", 2)
        assert 0.0 <= result["bray_curtis"] <= 1.0

    def test_known_jaccard_value(self) -> None:
        """Manually verify Jaccard for known k-mer sets.

        AACC k=2 => {AA, AC, CC} (3 distinct)
        CCGG k=2 => {CC, CG, GG} (3 distinct)
        intersection = {CC} => 1
        union = {AA, AC, CC, CG, GG} => 5
        jaccard = 1/5 = 0.2
        """
        result = compare_kmer_profiles("AACC", "CCGG", 2)
        assert abs(result["jaccard"] - 1.0 / 5.0) < 1e-10


# ===================================================================
# 8. sliding_window_kmer_complexity
# ===================================================================


class TestSlidingWindowKmerComplexity:
    """Tests for sliding_window_kmer_complexity."""

    def test_basic_output_structure(self) -> None:
        """Returns list of (position, entropy) tuples."""
        results = sliding_window_kmer_complexity("ATCGATCGATCG", 2, 6, 3)
        assert len(results) >= 1
        for pos, entropy in results:
            assert isinstance(pos, int)
            assert isinstance(entropy, float)
            assert entropy >= 0.0

    def test_window_positions(self) -> None:
        """Verify start positions match step increments."""
        seq = "ATCGATCGATCG"  # length 12
        results = sliding_window_kmer_complexity(seq, 2, 6, 3)
        positions = [pos for pos, _ in results]
        # Expect windows at 0, 3, 6 (9 would start a 6-base window ending at 15 > 12)
        # n=12, window_size=6, step=3: range(0, 12-6+1, 3) = range(0, 7, 3) = [0, 3, 6]
        assert positions == [0, 3, 6]

    def test_step_one_gives_max_windows(self) -> None:
        """Step=1 produces the most windows: n - window_size + 1."""
        seq = "ATCGATCG"  # length 8
        results = sliding_window_kmer_complexity(seq, 2, 4, 1)
        assert len(results) == 8 - 4 + 1  # 5

    def test_window_larger_than_sequence(self) -> None:
        """Window larger than sequence returns empty."""
        assert sliding_window_kmer_complexity("ATCG", 2, 100) == []

    def test_empty_sequence(self) -> None:
        assert sliding_window_kmer_complexity("", 2, 6) == []

    def test_k_zero(self) -> None:
        assert sliding_window_kmer_complexity("ATCG", 0, 4) == []

    def test_step_zero(self) -> None:
        assert sliding_window_kmer_complexity("ATCG", 2, 4, 0) == []

    def test_window_less_than_k(self) -> None:
        """window_size < k should return empty."""
        assert sliding_window_kmer_complexity("ATCGATCG", 5, 3) == []

    def test_homopolymer_windows_zero_entropy(self) -> None:
        """Windows of a homopolymer should have entropy 0."""
        results = sliding_window_kmer_complexity("AAAAAAAAAA", 2, 5, 2)
        for _, entropy in results:
            assert entropy == 0.0


# ===================================================================
# 9. find_homopolymers
# ===================================================================


class TestFindHomopolymers:
    """Tests for find_homopolymers."""

    def test_docstring_example(self) -> None:
        """AAATTTCCG with min_length=3."""
        result = find_homopolymers("AAATTTCCG", 3)
        assert result == [("A", 0, 3), ("T", 3, 3)]

    def test_long_homopolymer(self) -> None:
        """Single long run."""
        result = find_homopolymers("AAAAAAA", 3)
        assert result == [("A", 0, 7)]

    def test_no_homopolymers(self) -> None:
        """No runs meeting the minimum length."""
        result = find_homopolymers("ATCGATCG", 3)
        assert result == []

    def test_min_length_equals_run(self) -> None:
        """Run exactly equal to min_length is included."""
        result = find_homopolymers("AAA", 3)
        assert result == [("A", 0, 3)]

    def test_min_length_exceeds_run(self) -> None:
        """Run shorter than min_length is excluded."""
        result = find_homopolymers("AA", 3)
        assert result == []

    def test_mixed_homopolymers(self) -> None:
        """Multiple homopolymers of different bases."""
        result = find_homopolymers("AAACCCGGG", 3)
        assert result == [("A", 0, 3), ("C", 3, 3), ("G", 6, 3)]

    def test_case_insensitivity(self) -> None:
        """Lower-case input is handled correctly."""
        result = find_homopolymers("aaattt", 3)
        assert result == [("A", 0, 3), ("T", 3, 3)]

    def test_empty_sequence(self) -> None:
        assert find_homopolymers("", 3) == []

    def test_min_length_zero(self) -> None:
        """min_length=0 is treated as invalid."""
        assert find_homopolymers("AAAA", 0) == []

    def test_single_base_min_one(self) -> None:
        """Single base with min_length=1 is a run."""
        result = find_homopolymers("A", 1)
        assert result == [("A", 0, 1)]

    def test_sorted_by_position(self) -> None:
        """Results should be sorted by start position."""
        result = find_homopolymers("TTTTAAACCC", 3)
        positions = [start for _, start, _ in result]
        assert positions == sorted(positions)

    def test_adjacent_different_bases(self) -> None:
        """Adjacent runs of different bases are separate entries."""
        result = find_homopolymers("AAAAACCCC", 3)
        assert len(result) == 2
        assert result[0] == ("A", 0, 5)
        assert result[1] == ("C", 5, 4)


# ===================================================================
# 10. find_microsatellites
# ===================================================================


class TestFindMicrosatellites:
    """Tests for find_microsatellites."""

    def test_docstring_example(self) -> None:
        """ACACACACTTTT with min_repeats=3, max_unit_length=2."""
        result = find_microsatellites("ACACACACTTTT", min_repeats=3, max_unit_length=2)
        assert ("AC", 0, 8, 4) in result
        assert ("T", 8, 12, 4) in result

    def test_homopolymer_as_microsatellite(self) -> None:
        """Homopolymers are captured as unit-length-1 repeats."""
        result = find_microsatellites("AAAAAA", min_repeats=3, max_unit_length=2)
        assert ("A", 0, 6, 6) in result

    def test_reducibility_check(self) -> None:
        """Reducible units (e.g. AA from AAAAAA) should not be reported as length-2 repeats.

        AAAAAA with unit 'AA' is really 'A' x 6. The module should detect
        that 'AA' is reducible to 'A' and only report the unit-1 repeat.
        """
        result = find_microsatellites("AAAAAA", min_repeats=3, max_unit_length=2)
        units = [unit for unit, _, _, _ in result]
        assert "AA" not in units
        assert "A" in units

    def test_dinucleotide_repeat(self) -> None:
        """AT repeated 4 times."""
        result = find_microsatellites("ATATATAT", min_repeats=3, max_unit_length=2)
        assert ("AT", 0, 8, 4) in result

    def test_trinucleotide_repeat(self) -> None:
        """CAG repeated 5 times."""
        result = find_microsatellites("CAGCAGCAGCAGCAG", min_repeats=3, max_unit_length=6)
        assert ("CAG", 0, 15, 5) in result

    def test_no_microsatellites(self) -> None:
        """Sequence with no tandem repeats of sufficient count."""
        result = find_microsatellites("ATCG", min_repeats=3, max_unit_length=2)
        assert result == []

    def test_min_repeats_filtering(self) -> None:
        """Repeat count below min_repeats is excluded."""
        # AC appears only twice
        result = find_microsatellites("ACAC", min_repeats=3, max_unit_length=2)
        # AC x 2 does not meet min_repeats=3
        for item in result:
            assert item[3] >= 3  # all reported must meet threshold

    def test_sorted_by_start_position(self) -> None:
        """Results are sorted by start position."""
        result = find_microsatellites("ACACACACTTTTGGGGG", min_repeats=3, max_unit_length=2)
        positions = [start for _, start, _, _ in result]
        assert positions == sorted(positions)

    def test_empty_sequence(self) -> None:
        assert find_microsatellites("", min_repeats=3) == []

    def test_min_repeats_zero(self) -> None:
        assert find_microsatellites("ATCG", min_repeats=0) == []

    def test_max_unit_length_zero(self) -> None:
        assert find_microsatellites("ATCG", max_unit_length=0) == []

    def test_case_insensitivity(self) -> None:
        """Lower-case input produces same results."""
        result = find_microsatellites("acacacactttt", min_repeats=3, max_unit_length=2)
        assert ("AC", 0, 8, 4) in result

    def test_end_is_exclusive(self) -> None:
        """Verify that end position follows Python slice convention."""
        result = find_microsatellites("ACACACACTTTT", min_repeats=3, max_unit_length=2)
        for unit, start, end, repeat_count in result:
            assert end == start + len(unit) * repeat_count


# ===================================================================
# 11. sequence_complexity_profile
# ===================================================================


class TestSequenceComplexityProfile:
    """Tests for sequence_complexity_profile."""

    def test_basic_output(self) -> None:
        """Returns non-empty list of (position, complexity) tuples."""
        seq = "ATCG" * 20  # 80 bases
        profile = sequence_complexity_profile(seq, window_size=20, step=10)
        assert len(profile) >= 1
        for pos, complexity in profile:
            assert isinstance(pos, int)
            assert 0.0 <= complexity <= 1.0

    def test_homopolymer_low_complexity(self) -> None:
        """A homopolymer window has very low complexity."""
        seq = "A" * 100
        profile = sequence_complexity_profile(seq, window_size=50, step=10)
        for _, complexity in profile:
            # Only one distinct trinucleotide "AAA" out of 48 possible slots
            # complexity = 1/min(64, 48) = 1/48 ~ 0.021
            assert complexity < 0.05

    def test_window_positions_with_step(self) -> None:
        """Verify window positions honor step parameter."""
        seq = "ATCG" * 25  # 100 bases
        profile = sequence_complexity_profile(seq, window_size=20, step=20)
        positions = [pos for pos, _ in profile]
        # range(0, 100-20+1, 20) = [0, 20, 40, 60, 80]
        assert positions == [0, 20, 40, 60, 80]

    def test_empty_sequence(self) -> None:
        assert sequence_complexity_profile("") == []

    def test_window_larger_than_sequence(self) -> None:
        assert sequence_complexity_profile("ATCG", window_size=100) == []

    def test_window_too_small_for_trinucleotides(self) -> None:
        """Window smaller than k=3 returns empty."""
        assert sequence_complexity_profile("ATCGATCG", window_size=2) == []

    def test_known_complexity_value(self) -> None:
        """Manually compute complexity for a small window.

        Window "ATCGATCG" (8 bases), k=3: trinucleotides are
        ATC, TCG, CGA, GAT, ATC, CG -- wait, 8-3+1=6 positions:
        ATC, TCG, CGA, GAT, ATC, TCG => 4 distinct out of min(64, 6)=6
        complexity = 4/6 = 0.6667
        """
        profile = sequence_complexity_profile("ATCGATCG", window_size=8, step=1)
        assert len(profile) == 1
        assert abs(profile[0][1] - 4.0 / 6.0) < 1e-10


# ===================================================================
# 12. mask_low_complexity
# ===================================================================


class TestMaskLowComplexity:
    """Tests for mask_low_complexity."""

    def test_homopolymer_masked(self) -> None:
        """A long homopolymer is entirely masked."""
        seq = "A" * 100
        result = mask_low_complexity(seq, window_size=20, threshold=0.5)
        assert result == "N" * 100

    def test_high_complexity_preserved(self) -> None:
        """A complex sequence is not masked."""
        # Build a diverse sequence that has high trinucleotide complexity
        seq = "ATCGATCGATCGATCGATCG" * 5  # 100 bases, repeating but fairly diverse
        result = mask_low_complexity(seq, window_size=20, threshold=0.01)
        # With such a low threshold almost nothing should be masked
        assert "N" not in result

    def test_returns_uppercase(self) -> None:
        """Output is always upper-cased."""
        result = mask_low_complexity("atcgatcg", window_size=4, threshold=0.01)
        assert result == result.upper()

    def test_empty_sequence(self) -> None:
        assert mask_low_complexity("") == ""

    def test_window_larger_than_sequence_low_complexity(self) -> None:
        """When window > len(seq), evaluates entire sequence as one window."""
        seq = "AAAA"
        result = mask_low_complexity(seq, window_size=100, threshold=0.5)
        # Single window over 4 bases: only 2 trinucleotides (AAA, AAA) => 1 distinct
        # complexity = 1/min(64, 2) = 0.5; 0.5 < 0.5 is False => NOT masked
        # Actually: complexity = 1/2 = 0.5, threshold is 0.5, and condition is
        # "< threshold" so 0.5 < 0.5 is False => sequence preserved
        assert result == "AAAA"

    def test_window_larger_than_sequence_passes_threshold(self) -> None:
        """When the entire-sequence complexity is below threshold, mask everything."""
        seq = "AAAA"
        result = mask_low_complexity(seq, window_size=100, threshold=0.6)
        assert result == "NNNN"

    def test_mixed_regions(self) -> None:
        """A sequence with both low and high complexity regions."""
        low = "A" * 60
        high = "ATCGATCG" * 10  # 80 bases, diverse
        seq = low + high  # 140 bases total
        result = mask_low_complexity(seq, window_size=20, threshold=0.3)
        # The homopolymer region should be predominantly masked
        low_region = result[:60]
        assert low_region.count("N") > 30  # most of the low region is masked

    def test_output_length_preserved(self) -> None:
        """Output has the same length as input."""
        seq = "ATCGATCG" * 15  # 120 bases
        result = mask_low_complexity(seq, window_size=20, threshold=0.5)
        assert len(result) == len(seq)

    def test_window_less_than_k(self) -> None:
        """Window smaller than trinucleotide k=3 returns uppercased input."""
        seq = "atcg"
        result = mask_low_complexity(seq, window_size=2)
        assert result == "ATCG"


# ===================================================================
# Cross-function integration tests
# ===================================================================


class TestKmerIntegration:
    """Integration tests combining multiple k-mer functions."""

    def test_count_and_frequency_consistency(self) -> None:
        """count_kmers and kmer_frequencies agree on the k-mer set."""
        seq = "ATCGATCGATCG"
        k = 3
        counts = count_kmers(seq, k)
        freqs = kmer_frequencies(seq, k)
        assert set(counts.keys()) == set(freqs.keys())

    def test_spectrum_matches_counts(self) -> None:
        """kmer_spectrum is consistent with count_kmers."""
        seq = "ATCGATCGATCG"
        k = 2
        counts = count_kmers(seq, k)
        spectrum = kmer_spectrum(seq, k)

        # Rebuild spectrum from counts
        from collections import Counter

        expected_spectrum = dict(Counter(counts.values()))
        assert spectrum == expected_spectrum

    def test_unique_kmers_subset_of_counts(self) -> None:
        """find_unique_kmers returns only k-mers with count 1."""
        seq = "ATCGATCG"
        k = 2
        counts = count_kmers(seq, k)
        unique = find_unique_kmers(seq, k)
        for kmer in unique:
            assert counts[kmer] == 1

    def test_overrepresented_uses_frequencies(self) -> None:
        """find_overrepresented_kmers uses the same frequencies as kmer_frequencies."""
        seq = "AAACCCTTTTGG"
        k = 2
        freqs = kmer_frequencies(seq, k)
        overrep = find_overrepresented_kmers(seq, k, threshold=0.0)
        expected_rate = 1.0 / (4**k)
        for kmer, fold in overrep.items():
            assert abs(fold - freqs[kmer] / expected_rate) < 1e-10

    def test_diversity_zero_implies_single_kmer(self) -> None:
        """Zero diversity implies only one distinct k-mer."""
        seq = "AAAAAAA"
        k = 3
        entropy = kmer_diversity(seq, k)
        assert entropy == 0.0
        counts = count_kmers(seq, k)
        assert len(counts) == 1

    def test_compare_self_perfect_similarity(self) -> None:
        """Comparing a sequence against itself is perfect similarity."""
        seq = "ATCGATCGATCGATCG"
        result = compare_kmer_profiles(seq, seq, 3)
        assert abs(result["cosine"] - 1.0) < 1e-10
        assert abs(result["jaccard"] - 1.0) < 1e-10
        assert abs(result["bray_curtis"] - 0.0) < 1e-10

    def test_mask_preserves_high_complexity_windows(self) -> None:
        """sliding_window_kmer_complexity values above threshold survive masking."""
        seq = "ATCGATCGATCGATCGATCG"  # 20 bases
        # Compute complexity profile
        profile = sequence_complexity_profile(seq, window_size=10, step=5)
        # If all windows are above threshold 0.01, nothing should be masked
        min_complexity = min(c for _, c in profile) if profile else 0.0
        if min_complexity > 0.1:
            result = mask_low_complexity(seq, window_size=10, threshold=0.05)
            assert "N" not in result

    def test_full_pipeline_on_realistic_sequence(self) -> None:
        """Run all functions on a moderately realistic 200-base sequence."""
        seq = "ATCGATCGATCG" * 10 + "AAAAAAAAAA" * 4 + "GCTAGCTAGCTA" * 5
        k = 3

        counts = count_kmers(seq, k)
        assert len(counts) > 0

        freqs = kmer_frequencies(seq, k)
        assert abs(sum(freqs.values()) - 1.0) < 1e-10

        spectrum = kmer_spectrum(seq, k)
        assert sum(spectrum.values()) == len(counts)

        unique = find_unique_kmers(seq, k)
        assert all(isinstance(u, str) for u in unique)

        overrep = find_overrepresented_kmers(seq, k, threshold=1.0)
        assert all(fold >= 1.0 for fold in overrep.values())

        diversity = kmer_diversity(seq, k)
        assert diversity > 0.0

        homopolymers = find_homopolymers(seq, min_length=5)
        # The "AAAAAAAAAA" x 4 = 40 A's is one long run
        assert any(base == "A" and length >= 10 for base, _, length in homopolymers)

        microsats = find_microsatellites(seq, min_repeats=3, max_unit_length=4)
        assert isinstance(microsats, list)

        profile = sequence_complexity_profile(seq, window_size=20, step=10)
        assert len(profile) > 0

        masked = mask_low_complexity(seq, window_size=20, threshold=0.3)
        assert len(masked) == len(seq)
