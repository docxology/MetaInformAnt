"""Regression tests: ``find_nearest_gene`` measures from the variant interval, and indexes the genes once.

The distance used to run from the variant's *midpoint* to the gene, so a 100 bp variant 1,000 bp past a gene was
reported 1,050 bp away, a 6 kb variant 100 bp from a gene was reported 3,100 bp away, and a gene that abutted the
variant was not "at distance 0". It is now the gap between the two intervals, as ``bedtools closest`` measures
distance (0 for overlapping features), in half-open coordinates. Every expected value below is worked out by hand
from the coordinates. The function also rebuilt its interval index on every call; it now accepts a prebuilt
``IntervalIndex``.
"""

from __future__ import annotations

import random
import time
from typing import Any

import pytest

from metainformant.structural_variants.annotation.overlap import (
    GenomicInterval,
    IntervalIndex,
    annotate_gene_overlap,
    annotate_regulatory_overlap,
    find_nearest_gene,
)

GENES = [
    {"chrom": "chr1", "start": 1_000, "end": 5_000, "name": "GENE_A"},
    {"chrom": "chr1", "start": 6_000, "end": 9_000, "name": "GENE_B"},
    {"chrom": "chr1", "start": 20_000, "end": 21_000, "name": "GENE_C"},
    {"chrom": "chr2", "start": 100, "end": 200, "name": "GENE_D"},
]


def _nearest(start: int, end: int, chrom: str = "chr1", **kwargs: Any) -> tuple[str, int, str]:
    result = find_nearest_gene({"chrom": chrom, "start": start, "end": end}, GENES, **kwargs)
    return result["nearest_gene"], result["distance"], result["direction"]


class TestDistanceFromTheVariantInterval:
    """(variant interval) -> (gene, gap in bp, direction), by hand: gap = gene.start - end or start - gene.end."""

    def test_past_the_last_base_of_a_gene(self) -> None:
        # [10000, 10100) is 10000 - 9000 = 1000 bp past GENE_B. (From the midpoint 10050 the old distance was 1050.)
        assert _nearest(10_000, 10_100) == ("GENE_B", 1_000, "downstream")

    def test_before_the_first_base_of_a_gene(self) -> None:
        # [100, 200) is 1000 - 200 = 800 bp before GENE_A. (From the midpoint 150 the old distance was 850.)
        assert _nearest(100, 200) == ("GENE_A", 800, "upstream")

    def test_between_two_genes_the_nearer_edge_decides(self) -> None:
        # [5100, 5200): GENE_A ends at 5000 (gap 100), GENE_B starts at 6000 (gap 800).
        assert _nearest(5_100, 5_200) == ("GENE_A", 100, "downstream")
        # [5400, 5900): GENE_A gap 400, GENE_B gap 6000 - 5900 = 100.
        assert _nearest(5_400, 5_900) == ("GENE_B", 100, "upstream")

    def test_a_wide_variant_is_measured_from_its_edge_not_its_middle(self) -> None:
        # [9100, 15900) is 6,800 bp wide: GENE_B ends 100 bp before it, GENE_C starts 20000 - 15900 = 4100 bp after.
        # From the midpoint 12,500 the old distance to GENE_B was 3,500.
        assert _nearest(9_100, 15_900) == ("GENE_B", 100, "downstream")
        # [12000, 19000): GENE_B is 3000 bp behind it, GENE_C 20000 - 19000 = 1000 bp ahead (from the midpoint
        # 15,500 the old distance to GENE_C was 4,500).
        assert _nearest(12_000, 19_000) == ("GENE_C", 1_000, "upstream")

    def test_a_gene_that_abuts_the_variant_is_at_distance_zero_but_not_overlapping(self) -> None:
        # GENE_A is [1000, 5000): it ends exactly where [5000, 5100) begins. (The old midpoint distance was 50.)
        assert _nearest(5_000, 5_100) == ("GENE_A", 0, "downstream")
        # ...and starts exactly where [900, 1000) ends.
        assert _nearest(900, 1_000) == ("GENE_A", 0, "upstream")

    def test_an_overlapping_variant_reports_the_gene_with_the_largest_overlap(self) -> None:
        # [4000, 6500) overlaps GENE_A by 1000 bp (4000-5000) and GENE_B by 500 bp (6000-6500).
        assert _nearest(4_000, 6_500) == ("GENE_A", 0, "overlapping")
        # [4500, 8000) overlaps GENE_A by 500 and GENE_B by 2000.
        assert _nearest(4_500, 8_000) == ("GENE_B", 0, "overlapping")

    def test_a_variant_containing_a_gene_overlaps_it(self) -> None:
        assert _nearest(19_000, 22_000) == ("GENE_C", 0, "overlapping")

    @pytest.mark.parametrize(
        ("position", "expected"),
        [
            (3_000, ("GENE_A", 0, "overlapping")),  # strictly inside
            (1_001, ("GENE_A", 0, "overlapping")),
            (4_999, ("GENE_A", 0, "overlapping")),
            (1_000, ("GENE_A", 0, "upstream")),  # on the gene's first base: the junction is before the gene
            (5_000, ("GENE_A", 0, "downstream")),  # one past the gene's last base
            (5_400, ("GENE_A", 400, "downstream")),  # 5400 - 5000 against 6000 - 5400
            (5_700, ("GENE_B", 300, "upstream")),  # 5700 - 5000 against 6000 - 5700
        ],
    )
    def test_a_point_variant(self, position: int, expected: tuple[str, int, str]) -> None:
        assert _nearest(position, position) == expected

    def test_a_point_variant_halfway_between_two_genes_is_500_bp_from_either(self) -> None:
        name, distance, direction = _nearest(5_500, 5_500)  # 5500 - 5000 == 6000 - 5500
        assert distance == 500
        assert (name, direction) in {("GENE_A", "downstream"), ("GENE_B", "upstream")}

    def test_the_search_limit_is_inclusive(self) -> None:
        assert _nearest(10_000, 10_100, max_distance=1_000) == ("GENE_B", 1_000, "downstream")
        assert _nearest(10_000, 10_100, max_distance=999) == ("", -1, "none")

    def test_the_limit_applies_to_the_gap_not_to_the_distance_from_the_midpoint(self) -> None:
        # The variant is 200,000 bp wide: its edge is 60,000 bp past the gene, its midpoint 160,000 bp.
        genes = [{"chrom": "chr1", "start": 0, "end": 1_000, "name": "NEAR_THE_EDGE"}]
        wide = {"chrom": "chr1", "start": 61_000, "end": 261_000}
        assert find_nearest_gene(wide, genes, max_distance=60_000)["distance"] == 60_000
        assert find_nearest_gene(wide, genes, max_distance=59_999)["nearest_gene"] == ""

    def test_other_chromosomes_and_empty_databases(self) -> None:
        assert _nearest(150, 160, chrom="chr2") == ("GENE_D", 0, "overlapping")
        assert _nearest(150, 160, chrom="chrX") == ("", -1, "none")
        assert find_nearest_gene({"chrom": "chr1", "start": 1, "end": 2}, []) == {
            "nearest_gene": "",
            "distance": -1,
            "direction": "none",
            "gene_interval": None,
        }

    def test_the_interval_of_the_gene_is_returned(self) -> None:
        result = find_nearest_gene({"chrom": "chr1", "start": 10_000, "end": 10_100}, GENES)
        interval = result["gene_interval"]
        assert isinstance(interval, GenomicInterval) and (interval.start, interval.end) == (6_000, 9_000)


# ---------------------------------------------------------------------------
# IntervalIndex.query_nearest_region against a scan of every interval
# ---------------------------------------------------------------------------

CHROMS = ["chr1", "chr2", "chrX"]


def _random_intervals(rng: random.Random, n: int, span: int = 1_000_000) -> list[GenomicInterval]:
    intervals = []
    for i in range(n):
        start = rng.randint(0, span)
        kind = rng.random()
        if kind < 0.004:
            length = rng.randint(span // 20, span // 10)  # a gene hiding the shorter ones that start inside it
        elif kind < 0.15:
            length = 0  # a point feature
        elif kind < 0.5:
            length = rng.randint(1, 30)
        else:
            length = rng.randint(30, 3_000)
        intervals.append(GenomicInterval(rng.choice(CHROMS), start, start + length, f"i{i}"))
    intervals += [GenomicInterval("chr1", 500, 900, f"dup{k}") for k in range(3)]  # exact duplicates
    return intervals


def _gap(interval: GenomicInterval, start: int, end: int) -> int:
    return max(interval.start - end, start - interval.end, 0)


def _overlaps(interval: GenomicInterval, start: int, end: int) -> bool:
    return interval.start < end and interval.end > start


def _overlap_bp(interval: GenomicInterval, start: int, end: int) -> int:
    return max(0, min(end, interval.end) - max(start, interval.start))


class TestQueryNearestRegion:
    @pytest.mark.parametrize("max_distance", [0, 7, 400, 10**9])
    def test_matches_a_scan_of_every_interval(self, max_distance: int) -> None:
        rng = random.Random(1_000 + max_distance)
        intervals = _random_intervals(rng, 1_200)
        index = IntervalIndex(intervals)

        n_overlapping = n_gapped = n_abutting = 0
        for trial in range(2_500):
            chrom = rng.choice(CHROMS + ["chrZ"])
            start = rng.randint(-50, 1_000_050)
            end = start + rng.choice([0, 0, 1, 20, 400, 8_000])
            if trial % 4 == 0 and chrom != "chrZ":  # regions that begin where an interval ends, or end where one begins
                anchor = rng.choice([iv for iv in intervals if iv.chrom == chrom])
                if rng.random() < 0.5:
                    start, end = anchor.end, anchor.end + (end - start)
                else:
                    start, end = anchor.start - (end - start), anchor.start
            on_chrom = [iv for iv in intervals if iv.chrom == chrom]
            found, distance = index.query_nearest_region(chrom, start, end, max_distance)

            if not on_chrom or min(_gap(iv, start, end) for iv in on_chrom) > max_distance:
                assert (found, distance) == (None, -1)
                continue

            best_gap = min(_gap(iv, start, end) for iv in on_chrom)
            assert found is not None and found.chrom == chrom
            assert distance == best_gap == _gap(found, start, end)
            overlapping = [iv for iv in on_chrom if _overlaps(iv, start, end)]
            if overlapping:
                n_overlapping += 1
                assert _overlaps(found, start, end)
                assert _overlap_bp(found, start, end) == max(_overlap_bp(iv, start, end) for iv in overlapping)
            else:
                n_gapped += distance > 0
                n_abutting += distance == 0
        if max_distance == 10**9:  # the comparison covers overlaps, gaps and abutting features
            assert n_overlapping > 300 and n_gapped > 300 and n_abutting > 100

    def test_degenerate_indexes(self) -> None:
        assert IntervalIndex([]).query_nearest_region("chr1", 0, 10) == (None, -1)
        only = GenomicInterval("chr1", 100, 200, "only")
        single = IntervalIndex([only])
        assert single.query_nearest_region("chr1", 0, 50) == (only, 50)
        assert single.query_nearest_region("chr1", 250, 300) == (only, 50)
        assert single.query_nearest_region("chr1", 150, 160) == (only, 0)
        assert single.query_nearest_region("chr1", 250, 300, max_distance=49) == (None, -1)

    def test_a_long_gene_ending_before_the_variant_is_not_hidden_by_short_ones_inside_it(self) -> None:
        long_gene = GenomicInterval("chr1", 0, 1_000, "LONG")
        shorts = [GenomicInterval("chr1", 10 + 90 * k, 20 + 90 * k, f"short{k}") for k in range(10)]
        found, distance = IntervalIndex([long_gene, *shorts]).query_nearest_region("chr1", 1_500, 1_600)
        assert found is long_gene and distance == 500


# ---------------------------------------------------------------------------
# Building the index once
# ---------------------------------------------------------------------------


class TestPrebuiltIndex:
    VARIANTS = [
        (10_000, 10_100),
        (100, 200),
        (5_100, 5_200),
        (5_400, 5_900),
        (9_100, 15_900),
        (5_000, 5_100),
        (4_000, 6_500),
        (3_000, 3_000),
        (30_000, 30_500),
    ]

    def test_a_prebuilt_index_gives_the_answers_of_a_gene_list(self) -> None:
        index = IntervalIndex.from_features(GENES, feature_type="gene")
        for start, end in self.VARIANTS:
            variant = {"chrom": "chr1", "start": start, "end": end}
            from_list = find_nearest_gene(variant, GENES)
            from_index = find_nearest_gene(variant, index)
            assert from_index["nearest_gene"] == from_list["nearest_gene"]
            assert from_index["distance"] == from_list["distance"]
            assert from_index["direction"] == from_list["direction"]

    def test_an_index_of_genomic_intervals_is_accepted(self) -> None:
        intervals = [GenomicInterval("chr1", 1_000, 5_000, "GENE_A"), GenomicInterval("chr1", 6_000, 9_000, "GENE_B")]
        result = find_nearest_gene({"chrom": "chr1", "start": 10_000, "end": 10_100}, IntervalIndex(intervals))
        assert (result["nearest_gene"], result["distance"]) == ("GENE_B", 1_000)

    def test_the_annotators_accept_a_prebuilt_index_too(self) -> None:
        variants = [{"chrom": "chr1", "start": 4_000, "end": 6_500}, {"chrom": "chr1", "start": 7_000, "end": 7_100}]
        index = IntervalIndex.from_features(GENES, feature_type="gene")
        assert [v["overlapping_genes"] for v in annotate_gene_overlap(variants, index)] == [
            ["GENE_A", "GENE_B"],
            ["GENE_B"],
        ]
        elements = IntervalIndex.from_features(
            [{"chrom": "chr1", "start": 7_050, "end": 7_060, "name": "enhancer", "feature_type": "enhancer"}],
            "regulatory",
        )
        assert [v["overlapping_regulatory"] for v in annotate_regulatory_overlap(variants, elements)] == [
            [],
            ["enhancer"],
        ]

    def test_the_feature_type_defaults_apply_to_the_features_of_an_index(self) -> None:
        index = IntervalIndex.from_features([{"chrom": "chr1", "start": 1, "end": 9, "name": "x"}], feature_type="gene")
        [interval] = index.query_overlap("chr1", 0, 10)
        assert interval.feature_type == "gene"

    def test_looking_up_many_variants_does_not_reindex_the_genes_each_time(self) -> None:
        rng = random.Random(5)
        genes = [
            {"chrom": "chr1", "start": s, "end": s + rng.randint(100, 3_000), "name": f"g{i}"}
            for i, s in enumerate(rng.sample(range(0, 5_000_000, 50), 4_000))
        ]
        variants = [{"chrom": "chr1", "start": s, "end": s + 200} for s in rng.sample(range(0, 5_000_000, 50), 60)]
        index = IntervalIndex.from_features(genes, feature_type="gene")

        def seconds(gene_db: Any) -> float:
            best = float("inf")
            for _ in range(2):
                started = time.perf_counter()
                answers = [find_nearest_gene(v, gene_db)["nearest_gene"] for v in variants]
                best = min(best, time.perf_counter() - started)
            assert len(answers) == len(variants)
            return best

        # A list is indexed on every call (4,000 genes sorted and tree-built 60 times); an index is built once.
        # Both are timed here, so the comparison does not depend on the machine's speed; the ratio is about 100.
        assert seconds(genes) / seconds(index) > 5.0

    def test_the_genes_are_not_modified(self) -> None:
        genes = [dict(g) for g in GENES]
        find_nearest_gene({"chrom": "chr1", "start": 10_000, "end": 10_100}, genes)
        assert genes == GENES
