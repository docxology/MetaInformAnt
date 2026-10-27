"""Regression tests for ``filter_by_frequency`` and ``apply_blacklist``.

* A population database in a layout the filter cannot read used to be dropped without a word: the lookup came
  out empty, no variant was ever matched, and every variant passed as if the database had said "rare".
  A database none of whose entries can be read now raises; a mostly usable one is read, and the entries that
  cannot be are counted in a logged warning.
* A zero-length variant (an insertion point, ``start == end``) has no bases to overlap, so a blacklist region
  that contained it never counted it as overlapping and the variant was kept.
"""

from __future__ import annotations

import logging

import pytest

from metainformant.structural_variants.filtering.quality_filter import (
    apply_blacklist,
    filter_by_frequency,
    filter_by_size,
)

LOGGER = "metainformant.structural_variants.filtering.quality_filter"


def _variant(chrom: str = "chr1", start: int = 1_000, end: int = 3_000, sv_type: str = "DEL") -> dict[str, object]:
    return {"chrom": chrom, "start": start, "end": end, "sv_type": sv_type}


class TestUnsupportedPopulationDatabases:
    @pytest.mark.parametrize(
        "database",
        [
            pytest.param({"chr1_1000_3000_DEL": 0.9}, id="dict-key-not-chrom:start-end"),
            pytest.param({"chr1:1000": 0.9}, id="dict-key-without-end"),
            pytest.param({"chr1:1000-3000-DEL": "common"}, id="dict-value-not-a-frequency"),
            pytest.param({"chr1:1000-3000-DEL": {"af": 0.9}}, id="dict-value-is-a-record"),
            pytest.param({"chr1:1000-3000-DEL": 90.0}, id="dict-frequency-as-a-percentage"),
            pytest.param([{"chrom": "chr1", "start": 1_000, "end": 3_000, "AF": 0.9}], id="list-af-key-uppercase"),
            pytest.param([{"chrom": "chr1", "start": 1_000, "end": 3_000}], id="list-without-a-frequency"),
            pytest.param([{"chr": "chr1", "start": 1_000, "end": 3_000, "af": 0.9}], id="list-chrom-key-named-chr"),
            pytest.param([{"chrom": "chr1", "start": "1000", "end": 3_000, "af": 0.9}], id="list-text-position"),
            pytest.param([("chr1", 1_000, 3_000, 0.9)], id="list-of-tuples"),
            pytest.param("chr1:1000-3000-DEL", id="a-string"),
        ],
    )
    def test_a_database_with_no_readable_entry_raises(self, database: object) -> None:
        with pytest.raises(ValueError, match="layout is not supported"):
            filter_by_frequency([_variant()], population_db=database, max_af=0.01)  # type: ignore[arg-type]

    def test_the_message_says_what_is_supported(self) -> None:
        with pytest.raises(ValueError) as error:
            filter_by_frequency([_variant()], population_db={"chr1_1000_3000": 0.9})
        assert "chrom:start-end" in str(error.value)
        assert "'af'" in str(error.value)
        assert "chr1_1000_3000" in str(error.value)  # the first entry that could not be read

    def test_a_database_that_is_not_a_collection_raises(self) -> None:
        with pytest.raises(ValueError, match="dict or a list of dicts"):
            filter_by_frequency([_variant()], population_db=42)  # type: ignore[arg-type]

    def test_the_error_comes_even_when_every_variant_carries_its_own_frequency(self) -> None:
        variants = [dict(_variant(), af=0.5)]
        with pytest.raises(ValueError, match="layout is not supported"):
            filter_by_frequency(variants, population_db={"unreadable": 0.9})

    def test_supported_layouts_still_filter(self) -> None:
        as_dict = {"chr1:1000-3000-DEL": 0.9}
        as_list = [{"chrom": "chr1", "start": 1_000, "end": 3_000, "allele_frequency": 0.9, "sv_type": "DEL"}]
        for database in (as_dict, as_list):
            kept, stats = filter_by_frequency([_variant(), _variant("chr2")], population_db=database, max_af=0.01)
            assert [v["chrom"] for v in kept] == ["chr2"]
            assert stats.filtered_count == 1

    def test_an_empty_database_is_valid_and_removes_nothing(self) -> None:
        for database in ({}, []):
            kept, stats = filter_by_frequency([_variant()], population_db=database, max_af=0.01)
            assert len(kept) == 1 and stats.filtered_count == 0

    def test_integral_floats_and_integer_frequencies_are_readable(self) -> None:
        database = [{"chrom": "chr1", "start": 1_000.0, "end": 3_000.0, "af": 1, "sv_type": "DEL"}]
        kept, _ = filter_by_frequency([_variant()], population_db=database, max_af=0.01)
        assert kept == []


class TestPartiallyReadableDatabases:
    def test_unreadable_entries_are_skipped_and_counted_in_a_warning(self, caplog: pytest.LogCaptureFixture) -> None:
        database = [
            {"chrom": "chr1", "start": 1_000, "end": 3_000, "af": 0.9, "sv_type": "DEL"},
            {"chrom": "chr1", "start": 1_000, "end": 3_000, "AF": 0.9, "sv_type": "DEL"},
            {"chrom": "chr1", "start": "x", "end": 3_000, "af": 0.9},
        ]
        with caplog.at_level(logging.WARNING, logger=LOGGER):
            kept, _ = filter_by_frequency([_variant(), _variant("chr2")], population_db=database, max_af=0.01)

        assert [v["chrom"] for v in kept] == ["chr2"]  # the readable entry still filters
        warnings = [r.getMessage() for r in caplog.records if r.levelno == logging.WARNING]
        assert any("skipped 2 of 3 population_db entries" in message for message in warnings)

    def test_a_fully_readable_database_logs_no_warning(self, caplog: pytest.LogCaptureFixture) -> None:
        database = [{"chrom": "chr1", "start": 1_000, "end": 3_000, "af": 0.9, "sv_type": "DEL"}]
        with caplog.at_level(logging.WARNING, logger=LOGGER):
            filter_by_frequency([_variant()], population_db=database, max_af=0.01)
        assert not [r for r in caplog.records if r.levelno == logging.WARNING]


class TestZeroLengthVariantsAndBlacklists:
    REGION = [("chr1", 1_000, 3_000)]

    @staticmethod
    def _point(position: int, chrom: str = "chr1") -> dict[str, object]:
        return {"chrom": chrom, "start": position, "end": position, "sv_type": "INS"}

    def test_a_point_variant_inside_a_region_is_removed(self) -> None:
        kept, stats = apply_blacklist([self._point(2_000)], self.REGION)
        assert kept == []
        assert stats.filtered_count == 1

    def test_the_same_holds_for_a_dict_region(self) -> None:
        region = [{"chrom": "chr1", "start": 1_000, "end": 3_000}]
        kept, _ = apply_blacklist([self._point(1_001), self._point(2_999)], region)
        assert kept == []

    @pytest.mark.parametrize("position", [999, 1_000, 3_000, 3_001])
    def test_a_point_variant_outside_or_on_the_edge_of_a_region_is_kept(self, position: int) -> None:
        # On the edge it only abuts the region, as an interval that touches it without a shared base does.
        variant = self._point(position)
        kept, _ = apply_blacklist([variant], self.REGION)
        assert kept == [variant]

    def test_a_point_variant_on_another_chromosome_is_kept(self) -> None:
        variant = self._point(2_000, chrom="chr2")
        assert apply_blacklist([variant], self.REGION)[0] == [variant]

    def test_a_point_inside_the_region_is_fully_overlapping_whatever_the_minimum_fraction(self) -> None:
        for fraction in (0.0, 0.25, 0.99):
            kept, _ = apply_blacklist([self._point(2_000)], self.REGION, min_overlap_fraction=fraction)
            assert kept == []

    def test_only_the_point_inside_a_region_among_several_is_removed(self) -> None:
        regions = [("chr1", 100, 200), ("chr1", 1_000, 3_000), ("chr1", 5_000, 6_000)]
        variants = [self._point(150), self._point(500), self._point(2_500), self._point(5_500), self._point(9_000)]
        kept, stats = apply_blacklist(variants, regions)
        assert [v["start"] for v in kept] == [500, 9_000]
        assert stats.filtered_count == 3

    def test_intervals_are_judged_as_before(self) -> None:
        inside = {"chrom": "chr1", "start": 1_500, "end": 2_500, "sv_type": "DEL"}
        abutting = {"chrom": "chr1", "start": 3_000, "end": 3_500, "sv_type": "DEL"}
        kept, _ = apply_blacklist([inside, abutting], self.REGION)
        assert kept == [abutting]


class TestSizeFilterKeepsInsertionsWithoutAMeasuredLength:
    """``call_structural_variants`` reports an insertion with size 0 (not measured) and a start-to-end span that is
    only its junction window, so a size filter has nothing to compare with the minimum."""

    @staticmethod
    def _insertion(**fields: object) -> dict[str, object]:
        return {"chrom": "chr1", "start": 1_000, "end": 1_000, "sv_type": "INS", **fields}

    def test_an_insertion_with_no_measured_size_is_kept(self) -> None:
        unmeasured = [self._insertion(size=0), self._insertion(), self._insertion(end=1_020, size=0)]
        kept, stats = filter_by_size(unmeasured, min_size=50, max_size=1_000)
        assert kept == unmeasured
        assert stats.filtered_count == 0

    def test_an_insertion_with_a_measured_size_is_filtered_on_it(self) -> None:
        small, fitting, large = (self._insertion(size=size) for size in (30, 100, 5_000))
        kept, _ = filter_by_size([small, fitting, large], min_size=50, max_size=1_000)
        assert kept == [fitting]

    def test_other_types_without_a_size_still_use_their_span(self) -> None:
        short = {"chrom": "chr1", "start": 1_000, "end": 1_020, "sv_type": "DEL"}
        long = {"chrom": "chr1", "start": 1_000, "end": 9_000, "sv_type": "DEL", "size": 0}
        kept, _ = filter_by_size([short, long], min_size=50)
        assert kept == [long]

    def test_the_type_may_be_the_enum(self) -> None:
        from metainformant.structural_variants.detection.sv_calling import SVType

        variant = self._insertion(sv_type=SVType.INS, size=0)
        assert filter_by_size([variant], min_size=50)[0] == [variant]
