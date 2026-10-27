"""Regression tests for ``IntervalIndex`` (gene and regulatory overlap, nearest gene).

``query_overlap`` scanned every interval that starts before the query end, and ``query_nearest`` looked
only at the three intervals next to the position in start order. A gene spanning 0-899,000 is far from a
variant at 900,000 in that order whenever other genes start in between, so it was never found, and even a
position *inside* such a gene was reported as far from the nearest gene.

The index is checked against a brute-force scan over adversarial data (giant intervals, zero-length ones,
duplicates, several chromosomes), against the cases from the defect report, and for how its cost grows.
"""

from __future__ import annotations

import random
import time

import pytest

from metainformant.structural_variants.annotation.overlap import (
    GenomicInterval,
    IntervalIndex,
    annotate_gene_overlap,
    find_nearest_gene,
)

CHROMS = ["chr1", "chr2", "chrX"]


def _random_intervals(rng: random.Random, n: int, span: int = 1_000_000) -> list[GenomicInterval]:
    intervals = []
    for i in range(n):
        start = rng.randint(0, span)
        kind = rng.random()
        if kind < 0.05:
            length = rng.randint(span // 2, span)  # a gene spanning most of the chromosome
        elif kind < 0.10:
            length = 0  # a point feature
        elif kind < 0.50:
            length = rng.randint(1, 50)
        else:
            length = rng.randint(50, 5_000)
        intervals.append(GenomicInterval(chrom=rng.choice(CHROMS), start=start, end=start + length, name=f"i{i}"))
    # Exact duplicates and shared starts.
    intervals += [GenomicInterval(chrom="chr1", start=500, end=900, name=f"dup{k}") for k in range(3)]
    intervals += [
        GenomicInterval(chrom="chr1", start=500, end=500 + 100 * k, name=f"same_start{k}") for k in range(1, 4)
    ]
    return intervals


def _brute_overlap(intervals: list[GenomicInterval], chrom: str, start: int, end: int) -> list[GenomicInterval]:
    on_chrom = sorted((iv for iv in intervals if iv.chrom == chrom), key=lambda iv: (iv.start, iv.end))
    return [iv for iv in on_chrom if iv.start < end and iv.end > start]


def _distance(iv: GenomicInterval, position: int) -> int:
    if iv.start <= position < iv.end:
        return 0
    return iv.start - position if position < iv.start else position - iv.end


class TestAgainstBruteForce:
    def test_query_overlap_matches_a_scan_of_every_interval(self) -> None:
        rng = random.Random(20261027)
        intervals = _random_intervals(rng, 1_500)
        index = IntervalIndex(intervals)

        non_empty = 0
        for _ in range(1_500):
            chrom = rng.choice(CHROMS + ["chrZ"])
            start = rng.randint(-100, 1_100_000)
            end = start + rng.choice([0, 1, 10, 500, 20_000, 400_000])
            expected = _brute_overlap(intervals, chrom, start, end)
            assert [iv.name for iv in index.query_overlap(chrom, start, end)] == [iv.name for iv in expected]
            non_empty += bool(expected)
        assert non_empty > 500  # the comparison is not mostly empty answers

    @pytest.mark.parametrize("max_distance", [0, 10, 1_000, 10**9])
    def test_query_nearest_matches_a_scan_of_every_interval(self, max_distance: int) -> None:
        rng = random.Random(max_distance + 3)
        intervals = _random_intervals(rng, 1_200)
        index = IntervalIndex(intervals)

        for _ in range(1_200):
            chrom = rng.choice(CHROMS + ["chrZ"])
            position = rng.randint(-100, 1_100_000)
            distances = [_distance(iv, position) for iv in intervals if iv.chrom == chrom]
            found, distance = index.query_nearest(chrom, position, max_distance)

            if not distances or min(distances) > max_distance:
                assert (found, distance) == (None, -1)
                continue
            assert found is not None and found.chrom == chrom
            assert distance == min(distances) == _distance(found, position)
            if distance == 0 and any(iv.chrom == chrom and iv.start <= position < iv.end for iv in intervals):
                assert found.start <= position < found.end  # a spanning interval is preferred to an adjacent one

    def test_degenerate_indexes(self) -> None:
        empty = IntervalIndex([])
        assert empty.query_overlap("chr1", 0, 10) == []
        assert empty.query_nearest("chr1", 5) == (None, -1)

        only = GenomicInterval("chr1", 100, 200, "only")
        single = IntervalIndex([only])
        assert single.query_overlap("chr1", 150, 151) == [only]
        assert single.query_overlap("chr1", 200, 300) == []  # end is exclusive
        assert single.query_overlap("chr1", 0, 100) == []  # and a query ending at the start does not touch it
        assert single.query_nearest("chr1", 50) == (only, 50)
        assert single.query_nearest("chr1", 250) == (only, 50)
        assert single.query_nearest("chr1", 250, max_distance=49) == (None, -1)


# ---------------------------------------------------------------------------
# The reported cases
# ---------------------------------------------------------------------------


def _genes_around_a_giant() -> list[GenomicInterval]:
    giant = GenomicInterval("chr1", 0, 899_000, "GIANT")
    # Other genes start inside the giant: in start order they sit between it and a variant at 900,000.
    inside = [GenomicInterval("chr1", 100_000 * k, 100_000 * k + 100, f"inner{k}") for k in range(1, 9)]
    after = GenomicInterval("chr1", 950_000, 960_000, "after")
    return [giant, *inside, after]


class TestReportedCases:
    def test_gene_spanning_zero_to_899000_is_the_nearest_to_900000(self) -> None:
        index = IntervalIndex(_genes_around_a_giant())
        found, distance = index.query_nearest("chr1", 900_000)
        assert found is not None and found.name == "GIANT"
        assert distance == 1_000

    def test_position_inside_the_giant_is_at_distance_zero(self) -> None:
        index = IntervalIndex(_genes_around_a_giant())
        found, distance = index.query_nearest("chr1", 450_000)
        assert found is not None and found.name == "GIANT"
        assert distance == 0

    def test_find_nearest_gene_finds_the_giant(self) -> None:
        genes = _genes_around_a_giant()
        result = find_nearest_gene({"chrom": "chr1", "start": 900_000, "end": 900_100}, genes)
        assert result["nearest_gene"] == "GIANT"
        assert 0 < result["distance"] <= 1_100
        assert result["direction"] == "downstream"

    def test_a_variant_inside_the_giant_overlaps_it(self) -> None:
        genes = _genes_around_a_giant()
        [annotated] = annotate_gene_overlap([{"chrom": "chr1", "start": 450_000, "end": 450_500}], genes)
        assert annotated["overlapping_genes"] == ["GIANT"]
        result = find_nearest_gene({"chrom": "chr1", "start": 450_000, "end": 450_500}, genes)
        assert (result["nearest_gene"], result["distance"], result["direction"]) == ("GIANT", 0, "overlapping")

    def test_a_long_interval_ending_before_the_position_is_not_hidden_by_shorter_ones(self) -> None:
        # LONG ends at 1,000; ten short intervals start after it in start order but all end before it too.
        long_interval = GenomicInterval("chr1", 0, 1_000, "LONG")
        shorts = [GenomicInterval("chr1", 10 + 90 * k, 20 + 90 * k, f"short{k}") for k in range(10)]
        index = IntervalIndex([long_interval, *shorts])
        found, distance = index.query_nearest("chr1", 1_500)
        assert found is long_interval
        assert distance == 500

    def test_overlap_reports_every_nested_interval(self) -> None:
        outer = GenomicInterval("chr1", 0, 10_000, "outer")
        inner = [GenomicInterval("chr1", 100 * k, 100 * k + 10, f"in{k}") for k in range(1, 40)]
        index = IntervalIndex([*inner, outer])
        names = [iv.name for iv in index.query_overlap("chr1", 2_000, 2_005)]
        assert names == ["outer", "in20"]


# ---------------------------------------------------------------------------
# Cost
# ---------------------------------------------------------------------------


def _timed_queries(n_intervals: int, n_queries: int = 1_000) -> float:
    """Best-of-two time for ``n_queries`` overlap and nearest queries on an index of ``n_intervals``."""
    rng = random.Random(n_intervals)
    span = n_intervals * 1_000  # constant density, so each query meets about as many intervals at any size
    intervals = []
    for i in range(n_intervals):
        start = rng.randrange(span)
        intervals.append(GenomicInterval("chr1", start, start + rng.randint(100, 5_000), f"g{i}"))
    index = IntervalIndex(intervals)
    queries = [rng.randrange(span) for _ in range(n_queries)]

    best = float("inf")
    for _ in range(2):
        started = time.perf_counter()
        for position in queries:
            index.query_overlap("chr1", position, position + 100)
            index.query_nearest("chr1", position)
        best = min(best, time.perf_counter() - started)
    return best


def test_query_cost_does_not_grow_with_the_number_of_intervals() -> None:
    """Ten times as many intervals must not make queries an order of magnitude slower.

    Scanning every interval before the query end made the cost linear in the index size, and worse than
    linear once the intervals no longer fit in the CPU cache: the ratio below was about 70 then. With the
    tree a query visits O(log n + hits) nodes; the ratio stays between 1 and about 5 (it is not 1 because a
    larger index is also slower to reach in memory), so the bound is loose on purpose. Both sizes are timed
    here, so the comparison does not depend on the machine's speed.
    """
    small = _timed_queries(10_000)
    large = _timed_queries(100_000)
    assert large / small < 25.0, f"10k intervals: {small:.3f}s, 100k intervals: {large:.3f}s"
