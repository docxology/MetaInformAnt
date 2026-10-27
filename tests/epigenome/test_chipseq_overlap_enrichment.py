"""ChIP-seq peak overlap, peak enrichment and overlapping motif occurrences.

* ``find_overlapping_peaks`` used a sweep that only works when the peaks inside each set do not overlap one another,
  and missed pairs otherwise. It is checked here against a brute-force comparison of every pair.
* ``calculate_peak_enrichment`` summed Poisson terms with ``expected**k / k!`` and raised ``OverflowError`` for
  realistic peak counts (and returned a p-value of exactly 0.0 where the tail was merely tiny). It is checked here
  against an independent exact computation with ``decimal``.
* ``find_motifs_in_peaks`` skipped occurrences that start inside the previous occurrence of the same motif.
* ``merge_overlapping_peaks`` took a p-value or q-value of exactly 0.0 for a missing one. Peak files store ``-log10``
  values, so every column above about 323 is read as a probability of 0.0.

Real objects and real files only. REAL IMPLEMENTATION.
"""

from __future__ import annotations

import math
import random
from decimal import Decimal, getcontext
from pathlib import Path
from typing import List, Set, Tuple

import pytest

from metainformant.epigenome.assays.chipseq import (
    ChIPPeak,
    calculate_peak_enrichment,
    find_motifs_in_peaks,
    find_overlapping_peaks,
    merge_overlapping_peaks,
)

# ---------------------------------------------------------------------------
# find_overlapping_peaks
# ---------------------------------------------------------------------------


def _brute_force(peaks1: List[ChIPPeak], peaks2: List[ChIPPeak], min_overlap: int) -> Set[Tuple[int, int]]:
    """Every (peak1, peak2) pair on the same chromosome that ``overlaps_with`` accepts, by object identity."""
    return {(id(a), id(b)) for a in peaks1 for b in peaks2 if a.overlaps_with(b, min_overlap)}


def _random_peaks(
    rng: random.Random, count: int, chromosomes: Tuple[str, ...], span: int, max_length: int
) -> List[ChIPPeak]:
    peaks = []
    for _ in range(count):
        start = rng.randrange(0, span)
        peaks.append(ChIPPeak(rng.choice(chromosomes), start, start + rng.randrange(1, max_length)))
    return peaks


class TestFindOverlappingPeaks:
    def test_nested_peaks_inside_one_set_are_all_paired(self) -> None:
        # A (0-100) holds A' (10-20). B (15-30) and C (50-60) are in the other set. The old sweep paired A with B and
        # with C, then moved on without ever comparing A' with B.
        a = ChIPPeak("chr1", 0, 100)
        a_inner = ChIPPeak("chr1", 10, 20)
        b = ChIPPeak("chr1", 15, 30)
        c = ChIPPeak("chr1", 50, 60)

        pairs = find_overlapping_peaks([a, a_inner], [b, c])

        assert {(id(x), id(y)) for x, y in pairs} == {(id(a), id(b)), (id(a), id(c)), (id(a_inner), id(b))}
        assert len(pairs) == 3

    def test_a_long_peak_overlapping_many_short_peaks(self) -> None:
        long_peak = ChIPPeak("chr1", 0, 1000)
        shorts = [ChIPPeak("chr1", 10 * i, 10 * i + 5) for i in range(40)]

        pairs = find_overlapping_peaks([long_peak], shorts)
        reverse_pairs = find_overlapping_peaks(shorts, [long_peak])

        assert len(pairs) == 40
        assert len(reverse_pairs) == 40

    @pytest.mark.parametrize("min_overlap", [1, 3, 10, 0])
    def test_matches_brute_force_on_random_peak_sets(self, min_overlap: int) -> None:
        rng = random.Random(20261027 + min_overlap)
        chromosomes = ("chr1", "chr2", "chrX")

        for _ in range(200):
            peaks1 = _random_peaks(rng, rng.randrange(0, 13), chromosomes, span=400, max_length=120)
            peaks2 = _random_peaks(rng, rng.randrange(0, 13), chromosomes, span=400, max_length=120)

            pairs = find_overlapping_peaks(peaks1, peaks2, min_overlap)

            found = {(id(x), id(y)) for x, y in pairs}
            assert found == _brute_force(peaks1, peaks2, min_overlap)
            assert len(pairs) == len(found)  # no pair is reported twice

    def test_min_overlap_counts_bases(self) -> None:
        left = ChIPPeak("chr1", 100, 200)
        right = ChIPPeak("chr1", 195, 300)  # shares 5 bases

        assert len(find_overlapping_peaks([left], [right], min_overlap=5)) == 1
        assert find_overlapping_peaks([left], [right], min_overlap=6) == []

    def test_peaks_on_other_chromosomes_never_pair(self) -> None:
        assert find_overlapping_peaks([ChIPPeak("chr1", 0, 100)], [ChIPPeak("chr2", 0, 100)]) == []

    def test_inputs_are_left_in_order(self) -> None:
        peaks1 = [ChIPPeak("chr1", 500, 600), ChIPPeak("chr1", 100, 200)]
        peaks2 = [ChIPPeak("chr1", 550, 650), ChIPPeak("chr1", 150, 250)]
        before1, before2 = list(peaks1), list(peaks2)

        find_overlapping_peaks(peaks1, peaks2)

        assert peaks1 == before1
        assert peaks2 == before2

    def test_order_is_by_chromosome_then_start(self) -> None:
        peaks1 = [ChIPPeak("chr2", 50, 90), ChIPPeak("chr1", 300, 400), ChIPPeak("chr1", 10, 60)]
        peaks2 = [ChIPPeak("chr1", 0, 500), ChIPPeak("chr2", 0, 100)]

        pairs = find_overlapping_peaks(peaks1, peaks2)

        # chromosomes in order of first appearance in peaks1, then by the start of peak1
        assert [(a.chromosome, a.start) for a, _ in pairs] == [("chr2", 50), ("chr1", 10), ("chr1", 300)]


# ---------------------------------------------------------------------------
# calculate_peak_enrichment
# ---------------------------------------------------------------------------


def _exact_poisson_upper_tail(observed: int, expected: float, precision: int = 1200) -> Decimal:
    """P(X >= observed) for X ~ Poisson(expected), summed exactly with ``decimal`` (no log tricks, no cancellation)."""
    getcontext().prec = precision
    mean = Decimal(expected)
    term = (-mean).exp()
    below = Decimal(0)
    for k in range(observed):
        below += term
        term = term * mean / (k + 1)
    return Decimal(1) - below


def _peaks(count: int) -> List[ChIPPeak]:
    return [ChIPPeak("chr1", 10 * i, 10 * i + 5) for i in range(count)]


class TestPeakEnrichment:
    @pytest.mark.parametrize(
        ("observed", "expected"),
        [(200, 100.0), (1000, 800), (3000, 1000.0), (20000, 15000.5), (200, 100)],
    )
    def test_realistic_peak_counts_do_not_overflow(self, observed: int, expected: float) -> None:
        # (200, 100.0) raised OverflowError (float ** int); (1000, 800) raised it from int / int
        result = calculate_peak_enrichment(_peaks(observed), genome_size=10**9, expected_peaks=expected)

        assert 0.0 <= result["enrichment_p_value"] <= 1.0
        assert math.isfinite(result["enrichment_log10_p_value"])
        assert result["enrichment_ratio"] == pytest.approx(observed / expected)
        assert result["fold_enrichment"] == result["enrichment_ratio"]

    @pytest.mark.parametrize(
        ("observed", "expected"),
        [(1, 2.5), (5, 3.0), (50, 40), (100, 100.0), (101, 100.0), (200, 100), (1000, 800), (1000, 1000.0)],
    )
    def test_p_value_matches_the_exact_poisson_tail(self, observed: int, expected: float) -> None:
        result = calculate_peak_enrichment(_peaks(observed), genome_size=10**9, expected_peaks=expected)

        exact = _exact_poisson_upper_tail(observed, expected)
        assert result["enrichment_p_value"] == pytest.approx(float(exact), rel=1e-9)
        assert result["enrichment_log10_p_value"] == pytest.approx(float(exact.log10()), abs=1e-9)

    def test_a_tiny_tail_is_not_reported_as_zero(self) -> None:
        # 1 - P(X < 200) at mean 100 is about 9e-19: the old subtraction returned exactly 0.0
        result = calculate_peak_enrichment(_peaks(200), genome_size=10**9, expected_peaks=100)

        assert result["enrichment_p_value"] > 0.0
        assert result["enrichment_p_value"] == pytest.approx(9.343150072988e-19, rel=1e-9)

    def test_log10_p_value_stays_finite_where_the_p_value_underflows(self) -> None:
        # The true tail is about 10 ** -564.7, far below the smallest float
        result = calculate_peak_enrichment(_peaks(3000), genome_size=10**9, expected_peaks=1000)

        assert result["enrichment_p_value"] == 0.0
        assert result["enrichment_log10_p_value"] == pytest.approx(
            float(_exact_poisson_upper_tail(3000, 1000).log10()), abs=1e-6
        )

    def test_p_value_falls_as_more_peaks_are_observed(self) -> None:
        p_values = [
            calculate_peak_enrichment(_peaks(n), genome_size=10**9, expected_peaks=100)["enrichment_p_value"]
            for n in (60, 100, 120, 140, 160)
        ]

        assert p_values == sorted(p_values, reverse=True)
        assert len(set(p_values)) == len(p_values)
        assert p_values[0] > 0.99  # far fewer peaks than expected: seeing at least 60 is almost certain

    def test_no_significance_without_a_positive_expectation(self) -> None:
        for expected in (None, 0):
            result = calculate_peak_enrichment(_peaks(10), genome_size=1000, expected_peaks=expected)
            assert "enrichment_p_value" not in result
            assert "enrichment_ratio" not in result

    def test_genome_level_figures_are_unchanged(self) -> None:
        result = calculate_peak_enrichment([ChIPPeak("chr1", 100, 300), ChIPPeak("chr1", 500, 700)], genome_size=10000)

        assert result["total_peaks"] == 2
        assert result["total_peak_bases"] == 400
        assert result["genome_coverage"] == pytest.approx(0.04)
        assert result["average_peak_length"] == 200


# ---------------------------------------------------------------------------
# find_motifs_in_peaks: overlapping occurrences
# ---------------------------------------------------------------------------


def _genome(tmp_path: Path, sequence: str, chromosome: str = "chr1") -> str:
    fasta = tmp_path / "genome.fa"
    fasta.write_text(f">{chromosome}\n{sequence}\n")
    return str(fasta)


class TestMotifOccurrencesThatOverlap:
    def test_homopolymer_motif_reports_every_start(self, tmp_path: Path) -> None:
        # AAA occurs at 0, 1 and 2 in AAAAA; finditer resumed after each match and reported only position 0
        fasta = _genome(tmp_path, "AAAAA")
        results = find_motifs_in_peaks([ChIPPeak("chr1", 0, 5)], fasta, ["AAA"], window_size=200)

        assert results["match_counts"]["AAA"] == 3
        assert [(e["position"], e["strand"]) for e in results["motif_positions"]["AAA"]] == [
            (0, "+"),
            (1, "+"),
            (2, "+"),
        ]

    def test_self_overlapping_palindrome_is_found_on_both_strands(self, tmp_path: Path) -> None:
        # ATAT is its own reverse complement and occurs at 0 and 2 in ATATAT
        fasta = _genome(tmp_path, "ATATAT")
        results = find_motifs_in_peaks([ChIPPeak("chr1", 0, 6)], fasta, ["ATAT"], window_size=200)

        assert results["match_counts"]["ATAT"] == 4
        found = sorted((e["position"], e["strand"]) for e in results["motif_positions"]["ATAT"])
        assert found == [(0, "+"), (0, "-"), (2, "+"), (2, "-")]

    def test_empty_motif_is_rejected(self, tmp_path: Path) -> None:
        # An empty pattern matches between every pair of bases, which would report a site everywhere
        fasta = _genome(tmp_path, "ACGTACGT")
        with pytest.raises(ValueError, match="empty"):
            find_motifs_in_peaks([ChIPPeak("chr1", 0, 8)], fasta, [""])


# ---------------------------------------------------------------------------
# merge_overlapping_peaks: probabilities of exactly zero
# ---------------------------------------------------------------------------


class TestMergingKeepsAProbabilityOfZero:
    """Peak files store -log10 values, so the strongest peaks (a column above about 323) are probability 0.0. Merging
    tested the probabilities for truth, took 0.0 for a missing value, and dropped both of them."""

    def test_a_zero_is_a_value_not_a_gap(self) -> None:
        peaks = [
            ChIPPeak("chr1", 100, 300, p_value=0.0, q_value=0.0),
            ChIPPeak("chr1", 200, 400, p_value=1e-5, q_value=1e-3),
        ]

        (merged,) = merge_overlapping_peaks(peaks)

        assert (merged.p_value, merged.q_value) == (0.0, 0.0)

    def test_the_lower_probability_wins_whatever_the_order(self) -> None:
        peaks = [
            ChIPPeak("chr1", 100, 300, p_value=1e-3, q_value=0.05),
            ChIPPeak("chr1", 200, 400, p_value=1e-5, q_value=0.5),
        ]

        (merged,) = merge_overlapping_peaks(peaks)

        assert (merged.p_value, merged.q_value) == (1e-5, 0.05)

    def test_a_peak_without_a_probability_still_gives_none(self) -> None:
        peaks = [ChIPPeak("chr1", 100, 300, p_value=0.0, q_value=0.0), ChIPPeak("chr1", 200, 400)]

        (merged,) = merge_overlapping_peaks(peaks)

        assert (merged.p_value, merged.q_value) == (None, None)
