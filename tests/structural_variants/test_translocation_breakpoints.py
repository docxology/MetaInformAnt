"""Regression tests: translocations are matched on both of their breakpoints.

A translocation (``TRA``, or a ``BND`` breakend pair) is a zero-length record: ``chrom``/``start`` is its first
breakpoint, ``chrom2``/``pos2`` its second, and ``end`` repeats ``start``. The merge, deduplication, frequency,
blacklist, annotation and TAD functions used to read ``chrom``/``start``/``end`` alone:

* reciprocal overlap of two zero-length records is always 0, so ``merge_callsets``, ``deduplicate_variants`` and
  ``filter_by_frequency`` never matched a translocation to anything, and ``survivor_merge`` joined two of them on
  their first breakpoints whatever their partners were, and never joined a junction written from its other side;
* ``predict_tad_disruption`` looked at the first breakpoint only, and only at boundaries on its chromosome.

A junction is the same event when the breakpoints agree on both chromosomes, each pair of positions within the
merge distance, in either orientation. Every expected value below is worked out by hand from the coordinates.
"""

from __future__ import annotations

import importlib
import random
from pathlib import Path
from typing import Any

import pytest

from metainformant.structural_variants.annotation.functional_impact import (
    predict_functional_impact,
    predict_tad_disruption,
)
from metainformant.structural_variants.annotation.overlap import (
    annotate_gene_overlap,
    annotate_regulatory_overlap,
    find_nearest_gene,
)
from metainformant.structural_variants.filtering.merge import (
    MergedVariant,
    _load_vcf_variants,
    _parse_vcf_basic,
    deduplicate_variants,
    merge_callsets,
    survivor_merge,
)
from metainformant.structural_variants.filtering.quality_filter import apply_blacklist, filter_by_frequency
from metainformant.structural_variants.population import merge_sv_callsets


def _helpers() -> Any:
    """The breakpoint helpers, imported when a test needs them."""
    return importlib.import_module("metainformant.structural_variants._breakends")


def tra(chrom: str, start: int, chrom2: str, pos2: int, sv_type: str = "TRA", **extra: Any) -> dict[str, Any]:
    """A translocation record as ``dataclasses.asdict`` of a called ``StructuralVariant`` gives it."""
    record: dict[str, Any] = {
        "chrom": chrom,
        "start": start,
        "end": start,
        "sv_type": sv_type,
        "chrom2": chrom2,
        "pos2": pos2,
    }
    record.update(extra)
    return record


# ---------------------------------------------------------------------------
# The breakpoint helpers, against a brute-force statement of "the same junction"
# ---------------------------------------------------------------------------


def _near(x: tuple[str, int], y: tuple[str, int], distance: int) -> bool:
    return x[0] == y[0] and abs(x[1] - y[1]) <= distance


def _same_junction(a: dict[str, Any], b: dict[str, Any], distance: int) -> bool:
    """Brute force: the breakpoints agree as written, or with one record's two breakpoints swapped."""
    helpers = _helpers()
    (a1, a2), (b1, b2) = helpers.breakends(a), helpers.breakends(b)
    assert a2 is not None and b2 is not None
    return (_near(a1, b1, distance) and _near(a2, b2, distance)) or (
        _near(a1, b2, distance) and _near(a2, b1, distance)
    )


def _random_junction(rng: random.Random, span: int) -> dict[str, Any]:
    chrom = rng.choice(["chr1", "chr2"])
    chrom2 = rng.choice(["chr1", "chr2"])
    return tra(chrom, rng.randint(0, span), chrom2, rng.randint(0, span))


class TestBreakendHelpers:
    @pytest.mark.parametrize("distance", [0, 1, 5, 20])
    def test_canonical_comparison_finds_exactly_the_matches_of_either_orientation(self, distance: int) -> None:
        # Positions in a range a few distances wide, on one or two chromosomes: many pairs sit exactly at, just
        # inside and just outside the distance, including both breakpoints on one chromosome in opposite orders.
        helpers = _helpers()
        rng = random.Random(20261028 + distance)
        agreed = 0
        for _ in range(8_000):
            a, b = _random_junction(rng, 3 + 3 * distance), _random_junction(rng, 3 + 3 * distance)
            expected = _same_junction(a, b, distance)
            canonical = helpers.canonical_junction(a), helpers.canonical_junction(b)
            assert helpers.junctions_agree(*canonical, distance) is expected
            agreed += expected
        assert 150 < agreed < 7_850  # neither trivially true nor trivially false

    def test_junction_links_are_the_matching_pairs_and_only_those(self) -> None:
        helpers = _helpers()
        rng = random.Random(7)
        records = [_random_junction(rng, 40) for _ in range(120)]
        links = set(helpers.junction_links([helpers.canonical_junction(r) for r in records], 4))
        expected = {
            (i, j)
            for i in range(len(records))
            for j in range(i + 1, len(records))
            if _same_junction(records[i], records[j], 4)
        }
        assert links == expected
        assert len(expected) > 20

    def test_the_second_breakpoint_defaults_to_the_first_chromosome(self) -> None:
        breakends = _helpers().breakends
        assert breakends({"chrom": "chr1", "start": 10, "pos2": 99, "chrom2": ""}) == (("chr1", 10), ("chr1", 99))
        assert breakends({"chrom": "chr1", "start": 10, "pos2": 99, "chrom2": None}) == (("chr1", 10), ("chr1", 99))

    def test_end_is_the_position_on_chrom2_for_a_record_without_pos2(self) -> None:
        # The VCF style: CHR2 on another chromosome, END the position there.
        record = {"chrom": "chr1", "start": 10, "end": 500, "chrom2": "chr2", "sv_type": "TRA"}
        assert _helpers().breakends(record) == (("chr1", 10), ("chr2", 500))

    def test_a_zero_length_record_without_a_partner_is_a_single_breakend(self) -> None:
        helpers = _helpers()
        single = {"chrom": "chr1", "start": 10, "end": 10, "sv_type": "TRA"}
        assert helpers.breakends(single) == (("chr1", 10), None)
        assert helpers.variant_regions(single) == [("chr1", 10, 10)]

    def test_regions_of_a_translocation_are_its_two_breakpoints_and_of_an_interval_variant_its_interval(self) -> None:
        variant_regions = _helpers().variant_regions
        assert variant_regions(tra("chr1", 10, "chr2", 99)) == [("chr1", 10, 10), ("chr2", 99, 99)]
        assert variant_regions({"chrom": "chr1", "start": 10, "end": 99, "sv_type": "DEL"}) == [("chr1", 10, 99)]

    def test_a_coordinate_that_is_not_a_number_is_an_error(self) -> None:
        with pytest.raises(ValueError, match="pos2"):
            _helpers().breakends(tra("chr1", 10, "chr2", "99"))  # type: ignore[arg-type]


# ---------------------------------------------------------------------------
# merge_callsets
# ---------------------------------------------------------------------------


class TestMergeCallsets:
    def test_the_same_junction_from_two_callers_is_one_variant_with_median_breakpoints(self) -> None:
        callsets = {
            "delly": [tra("chr1", 1_000, "chr2", 5_000)],
            "manta": [tra("chr1", 1_020, "chr2", 5_030)],
        }
        merged, stats = merge_callsets(callsets)

        assert stats.n_output_variants == 1 and stats.n_multi_caller == 1
        [variant] = merged
        assert (variant.n_callers, variant.caller_names) == (2, ["delly", "manta"])
        # medians: (1000, 1020) -> 1010 and (5000, 5030) -> 5015
        assert (variant.chrom, variant.start, variant.end) == ("chr1", 1_010, 1_010)
        assert (variant.chrom2, variant.pos2) == ("chr2", 5_015)
        assert variant.size == 0

    def test_a_junction_written_from_its_other_side_is_the_same_event(self) -> None:
        callsets = {
            "delly": [tra("chr1", 1_000, "chr2", 5_000)],
            "manta": [tra("chr2", 5_030, "chr1", 1_020)],
        }
        [variant], _ = merge_callsets(callsets)
        assert variant.n_callers == 2
        # The merged junction lists its breakpoints in (chromosome, position) order whichever way the callers wrote it.
        assert (variant.chrom, variant.start, variant.chrom2, variant.pos2) == ("chr1", 1_010, "chr2", 5_015)

    def test_a_different_partner_is_a_different_event(self) -> None:
        callsets = {
            "delly": [tra("chr1", 1_000, "chr2", 5_000)],
            "manta": [tra("chr1", 1_000, "chr3", 5_000), tra("chr1", 1_000, "chr2", 90_000)],
        }
        merged, _ = merge_callsets(callsets)
        assert len(merged) == 3 and all(v.n_callers == 1 for v in merged)

    @pytest.mark.parametrize(("max_distance", "n_variants"), [(500, 1), (499, 2)])
    def test_the_distance_limit_is_inclusive(self, max_distance: int, n_variants: int) -> None:
        # Both breakpoints are exactly 500 bp apart.
        callsets = {"a": [tra("chr1", 1_000, "chr2", 5_000)], "b": [tra("chr1", 1_500, "chr2", 5_500)]}
        merged, _ = merge_callsets(callsets, max_distance=max_distance)
        assert len(merged) == n_variants

    def test_one_breakpoint_within_the_distance_is_not_enough(self) -> None:
        callsets = {"a": [tra("chr1", 1_000, "chr2", 5_000)], "b": [tra("chr1", 1_000, "chr2", 5_000 + 1_001)]}
        merged, _ = merge_callsets(callsets, max_distance=1_000)
        assert len(merged) == 2

    def test_calls_of_one_caller_never_merge(self) -> None:
        callsets = {"delly": [tra("chr1", 1_000, "chr2", 5_000), tra("chr1", 1_010, "chr2", 5_010)]}
        merged, _ = merge_callsets(callsets)
        assert len(merged) == 2

    def test_a_translocation_never_merges_with_an_interval_variant(self) -> None:
        interval = {"chrom": "chr1", "start": 900, "end": 1_100, "sv_type": "DEL"}
        callsets = {"a": [tra("chr1", 1_000, "chr2", 5_000)], "b": [interval]}
        merged, _ = merge_callsets(callsets, type_match=False, min_overlap=0.0)
        assert sorted(v.sv_type for v in merged) == ["DEL", "TRA"]

    def test_tra_and_bnd_are_one_type(self) -> None:
        callsets = {
            "a": [tra("chr1", 1_000, "chr2", 5_000, sv_type="TRA")],
            "b": [tra("chr1", 1_010, "chr2", 5_010, sv_type="BND")],
        }
        merged, _ = merge_callsets(callsets, type_match=True)
        assert len(merged) == 1 and merged[0].n_callers == 2

    def test_single_breakends_merge_only_with_single_breakends(self) -> None:
        single = {"chrom": "chr1", "start": 1_000, "end": 1_000, "sv_type": "BND"}
        near_single = dict(single, start=1_010, end=1_010)
        callsets = {
            "a": [single, tra("chr1", 7_000, "chr2", 5_000)],
            "b": [near_single, tra("chr1", 7_000, "chr2", 5_000)],
        }
        merged, _ = merge_callsets(callsets)
        by_start = {v.start: v for v in merged}
        assert by_start[1_005].n_callers == 2 and by_start[1_005].pos2 is None and by_start[1_005].chrom2 == ""
        assert by_start[7_000].n_callers == 2 and by_start[7_000].pos2 == 5_000
        assert len(merged) == 2
        # A single breakend does not stand in for a full junction at the same locus.
        mixed, _ = merge_callsets({"a": [single], "b": [tra("chr1", 1_000, "chr2", 5_000)]})
        assert len(mixed) == 2

    def test_end_is_read_as_the_position_on_chrom2_when_pos2_is_absent(self) -> None:
        vcf_style = {"chrom": "chr1", "start": 1_000, "end": 5_000, "sv_type": "TRA", "chrom2": "chr2"}
        merged, _ = merge_callsets({"a": [vcf_style], "b": [tra("chr1", 1_005, "chr2", 5_005)]})
        assert len(merged) == 1

    def test_interval_variants_are_merged_as_before_and_carry_no_partner(self) -> None:
        callsets = {
            "a": [{"chrom": "chr1", "start": 1_000, "end": 5_000, "sv_type": "DEL"}],
            "b": [{"chrom": "chr1", "start": 1_050, "end": 5_050, "sv_type": "DEL"}],
        }
        [variant], _ = merge_callsets(callsets)
        assert (variant.start, variant.end, variant.n_callers) == (1_025, 5_025, 2)
        assert (variant.chrom2, variant.pos2) == ("", None)

    def test_an_enum_sv_type_is_a_translocation_too(self) -> None:
        from metainformant.structural_variants.detection.sv_calling import SVType

        callsets = {
            "a": [tra("chr1", 1_000, "chr2", 5_000, sv_type=SVType.TRA)],
            "b": [tra("chr1", 1_010, "chr2", 5_010, sv_type=SVType.TRA)],
        }
        merged, _ = merge_callsets(callsets)
        assert len(merged) == 1


# ---------------------------------------------------------------------------
# survivor_merge
# ---------------------------------------------------------------------------


class TestSurvivorMerge:
    def test_translocations_with_different_partners_are_not_joined_on_their_first_breakpoint(self) -> None:
        # chr1:1000 -> chr2:5000 and chr1:1000 -> chr3:5000 share a first breakpoint only.
        merged, _ = survivor_merge(
            [[tra("chr1", 1_000, "chr2", 5_000)], [tra("chr1", 1_000, "chr3", 5_000)]], max_distance=500, min_callers=1
        )
        assert len(merged) == 2 and all(v.n_callers == 1 for v in merged)

    def test_a_junction_written_from_its_other_side_is_joined(self) -> None:
        merged, _ = survivor_merge(
            [[tra("chr1", 1_000, "chr2", 5_000)], [tra("chr2", 5_020, "chr1", 1_030)]], max_distance=500, min_callers=2
        )
        [variant] = merged
        assert variant.n_callers == 2
        assert (variant.chrom, variant.start, variant.chrom2, variant.pos2) == ("chr1", 1_015, "chr2", 5_010)

    @pytest.mark.parametrize(("max_distance", "expected"), [(500, 1), (499, 0)])
    def test_both_breakpoints_must_be_within_the_distance_inclusive(self, max_distance: int, expected: int) -> None:
        # The first breakpoints coincide; only the second ones are 500 bp apart.
        merged, _ = survivor_merge(
            [[tra("chr1", 1_000, "chr2", 5_000)], [tra("chr1", 1_000, "chr2", 5_500)]],
            max_distance=max_distance,
            min_callers=2,
        )
        assert len(merged) == expected

    @pytest.mark.parametrize(("max_distance", "expected"), [(500, 1), (499, 0)])
    def test_the_first_breakpoints_are_held_to_the_same_limit(self, max_distance: int, expected: int) -> None:
        merged, _ = survivor_merge(
            [[tra("chr1", 1_000, "chr2", 5_000)], [tra("chr1", 1_500, "chr2", 5_000)]],
            max_distance=max_distance,
            min_callers=2,
        )
        assert len(merged) == expected

    def test_a_translocation_does_not_join_a_nearby_short_interval_variant(self) -> None:
        # The deletion's two ends are each within 100 bp of the junction's first breakpoint.
        interval = {"chrom": "chr1", "start": 900, "end": 1_100, "sv_type": "DEL"}
        merged, _ = survivor_merge(
            [[tra("chr1", 1_000, "chr2", 5_000)], [interval]], max_distance=500, min_callers=1, type_match=False
        )
        assert len(merged) == 2

    def test_strand_match_applies_to_translocations(self) -> None:
        plus, minus = tra("chr1", 1_000, "chr2", 5_000, strand="+"), tra("chr1", 1_010, "chr2", 5_010, strand="-")
        assert len(survivor_merge([[plus], [minus]], strand_match=True, min_callers=1)[0]) == 2
        assert len(survivor_merge([[plus], [minus]], strand_match=False, min_callers=2)[0]) == 1

    def test_calls_of_one_caller_never_join(self) -> None:
        merged, _ = survivor_merge(
            [[tra("chr1", 1_000, "chr2", 5_000), tra("chr1", 1_010, "chr2", 5_010)]], min_callers=1
        )
        assert len(merged) == 2

    def test_dicts_tagged_with_their_caller_are_grouped_by_caller(self) -> None:
        calls = [
            tra("chr1", 1_000, "chr2", 5_000, _caller="delly"),
            tra("chr2", 5_010, "chr1", 1_010, _caller="manta"),
        ]
        merged, stats = survivor_merge(calls, max_distance=100, min_callers=2)
        assert stats.n_input_callsets == 2 and len(merged) == 1


class TestSurvivorMergeInputs:
    def test_a_dict_without_a_caller_is_refused_rather_than_counted_as_a_caller_of_its_own(self) -> None:
        # Two calls of one caller, untagged: each used to become "caller_<index>", so a lone caller's duplicate
        # pair passed min_callers=2 as if two callers agreed.
        calls = [
            {"chrom": "chr1", "start": 1_000, "end": 5_000, "sv_type": "DEL"},
            {"chrom": "chr1", "start": 1_010, "end": 5_010, "sv_type": "DEL"},
        ]
        with pytest.raises(ValueError, match="_caller"):
            survivor_merge(calls, max_distance=500, min_callers=2)

    def test_an_item_that_is_not_a_path_dict_or_list_is_refused_not_skipped(self) -> None:
        call = {"chrom": "chr1", "start": 1_000, "end": 5_000, "sv_type": "DEL"}
        with pytest.raises(TypeError, match="tuple"):
            survivor_merge([(call,), (call,)], min_callers=1)  # type: ignore[list-item]

    def test_tagged_dicts_and_lists_are_still_accepted(self) -> None:
        call = {"chrom": "chr1", "start": 1_000, "end": 5_000, "sv_type": "DEL"}
        merged, stats = survivor_merge([dict(call, _caller="a"), dict(call, _caller="b", start=1_020)], min_callers=2)
        assert stats.n_input_callsets == 2 and len(merged) == 1


# ---------------------------------------------------------------------------
# deduplicate_variants
# ---------------------------------------------------------------------------


class TestDeduplicateVariants:
    def test_duplicates_of_a_junction_collapse_to_the_best_call_in_either_orientation(self) -> None:
        variants = [
            tra("chr1", 1_000, "chr2", 5_000, quality=10.0),
            tra("chr1", 1_030, "chr2", 5_040, quality=30.0),
            tra("chr2", 5_020, "chr1", 1_010, quality=20.0),
        ]
        deduped, removed = deduplicate_variants(variants, max_distance=100)
        assert removed == 2
        assert [v["quality"] for v in deduped] == [30.0]

    def test_a_different_partner_is_kept(self) -> None:
        variants = [
            tra("chr1", 1_000, "chr2", 5_000),
            tra("chr1", 1_000, "chr3", 5_000),
            tra("chr1", 1_000, "chr2", 9_000),
        ]
        deduped, removed = deduplicate_variants(variants, max_distance=100)
        assert removed == 0 and len(deduped) == 3

    @pytest.mark.parametrize(("max_distance", "removed"), [(100, 1), (99, 0)])
    def test_the_distance_limit_is_inclusive(self, max_distance: int, removed: int) -> None:
        variants = [tra("chr1", 1_000, "chr2", 5_000), tra("chr1", 1_100, "chr2", 5_100)]
        assert deduplicate_variants(variants, max_distance=max_distance)[1] == removed

    def test_reciprocal_overlap_does_not_apply_to_a_translocation(self) -> None:
        variants = [tra("chr1", 1_000, "chr2", 5_000), tra("chr1", 1_010, "chr2", 5_010)]
        assert deduplicate_variants(variants, max_distance=100, min_reciprocal_overlap=1.0)[1] == 1

    def test_a_translocation_is_not_a_duplicate_of_an_interval_variant(self) -> None:
        variants = [tra("chr1", 1_000, "chr2", 5_000), {"chrom": "chr1", "start": 990, "end": 1_010, "sv_type": "DEL"}]
        deduped, removed = deduplicate_variants(
            variants, max_distance=100, type_match=False, min_reciprocal_overlap=0.0
        )
        assert removed == 0 and len(deduped) == 2

    def test_calls_with_a_missing_quality_can_be_compared(self) -> None:
        variants = [
            {"chrom": "chr1", "start": 1_000, "end": 5_000, "sv_type": "DEL", "quality": None},
            {"chrom": "chr1", "start": 1_010, "end": 5_010, "sv_type": "DEL", "quality": 12.0},
        ]
        deduped, removed = deduplicate_variants(variants)
        assert removed == 1 and deduped[0]["quality"] == 12.0


class TestMergeCallsetsQuality:
    def test_a_call_with_a_none_quality_counts_as_zero(self) -> None:
        # survivor_merge and deduplicate_variants read a None quality as 0; merge_callsets raised a TypeError.
        callsets = {
            "a": [{"chrom": "chr1", "start": 1_000, "end": 5_000, "sv_type": "DEL", "quality": None}],
            "b": [{"chrom": "chr1", "start": 1_010, "end": 5_010, "sv_type": "DEL", "quality": 40.0}],
        }
        [variant], _ = merge_callsets(callsets)
        # mean(0, 40) = 20, plus 100 * log10(2 + 1)
        assert variant.consensus_quality == pytest.approx(20.0 + 100.0 * 0.47712125471966244)


# ---------------------------------------------------------------------------
# filter_by_frequency
# ---------------------------------------------------------------------------


class TestFilterByFrequency:
    DB = [
        {"chrom": "chr1", "start": 1_000, "chrom2": "chr2", "pos2": 5_000, "af": 0.5, "sv_type": "BND"},
        {"chrom": "chr1", "start": 100, "end": 300, "af": 0.2, "sv_type": "DEL"},
    ]

    def test_a_translocation_matches_a_database_translocation_on_both_breakpoints(self) -> None:
        variants = [
            tra("chr1", 1_010, "chr2", 5_050),  # near both breakpoints
            tra("chr2", 5_010, "chr1", 990),  # the junction written from its other side
            tra("chr1", 1_010, "chr3", 5_050),  # same first breakpoint, other partner
            tra("chr1", 1_010, "chr2", 9_000),  # same partner chromosome, far second breakpoint
        ]
        kept, stats = filter_by_frequency(variants, population_db=self.DB, max_af=0.01, match_window=200)
        assert [(v["chrom"], v["chrom2"], v["pos2"]) for v in kept] == [
            ("chr1", "chr3", 5_050),
            ("chr1", "chr2", 9_000),
        ]
        assert stats.filtered_count == 2

    def test_the_match_window_is_inclusive(self) -> None:
        variant = tra("chr1", 1_200, "chr2", 5_200)
        assert filter_by_frequency([variant], population_db=self.DB, max_af=0.01, match_window=200)[0] == []
        assert len(filter_by_frequency([variant], population_db=self.DB, max_af=0.01, match_window=199)[0]) == 1

    def test_the_matched_frequency_is_recorded_and_the_highest_of_several_entries_wins(self) -> None:
        database = self.DB + [
            {"chrom": "chr1", "start": 1_050, "chrom2": "chr2", "pos2": 5_050, "af": 0.7, "sv_type": "TRA"}
        ]
        kept, _ = filter_by_frequency([tra("chr1", 1_000, "chr2", 5_000)], population_db=database, max_af=0.9)
        assert kept[0]["population_af"] == 0.7

    def test_interval_and_translocation_entries_do_not_match_each_other(self) -> None:
        # The interval entry chr1:100-300 is not a translocation's breakpoint pair, and the translocation entry is
        # not an interval's.
        deletion = {"chrom": "chr1", "start": 990, "end": 1_010, "sv_type": "DEL"}
        passed_tra, _ = filter_by_frequency([tra("chr1", 200, "chr1", 250)], population_db=self.DB, max_af=0.01)
        passed_del, _ = filter_by_frequency([deletion], population_db=self.DB, max_af=0.01)
        assert len(passed_tra) == 1 and len(passed_del) == 1

    def test_an_interval_variant_is_still_filtered_by_an_interval_entry(self) -> None:
        deletion = {"chrom": "chr1", "start": 100, "end": 300, "sv_type": "DEL"}
        assert filter_by_frequency([deletion], population_db=self.DB, max_af=0.01)[0] == []

    def test_a_database_of_only_translocations_is_readable(self) -> None:
        # Its entries have no 'end': they used to count as unreadable and the whole database was refused.
        database = [self.DB[0]]
        kept, _ = filter_by_frequency([tra("chr1", 1_000, "chr2", 5_000)], population_db=database, max_af=0.01)
        assert kept == []

    def test_a_translocation_entry_without_a_start_is_unreadable(self) -> None:
        with pytest.raises(ValueError, match="layout is not supported"):
            filter_by_frequency(
                [tra("chr1", 1_000, "chr2", 5_000)],
                population_db=[{"chrom": "chr1", "chrom2": "chr2", "pos2": 5_000, "af": 0.5, "sv_type": "TRA"}],
            )


# ---------------------------------------------------------------------------
# predict_tad_disruption
# ---------------------------------------------------------------------------


class TestTadDisruption:
    CHR1_BOUNDARY = {"chrom": "chr1", "start": 99_000, "end": 101_000, "genes": ["G1"]}  # midpoint 100,000
    CHR2_BOUNDARY = {"chrom": "chr2", "start": 499_000, "end": 501_000, "genes": ["G2"]}  # midpoint 500,000

    def test_the_second_breakpoint_is_checked_against_boundaries_of_its_own_chromosome(self) -> None:
        # The first breakpoint (chr1:2,000,000) is far from every boundary; the second sits on the chr2 boundary.
        variant = tra("chr1", 2_000_000, "chr2", 505_000)
        prediction = predict_tad_disruption(variant, [self.CHR1_BOUNDARY, self.CHR2_BOUNDARY])
        assert prediction.n_boundaries_disrupted == 1
        assert prediction.disrupted_boundaries == [500_000]
        assert prediction.disrupted_chroms == ["chr2"]
        assert prediction.genes_in_affected_tads == ["G2"]
        assert (prediction.variant_chrom2, prediction.variant_pos2) == ("chr2", 505_000)

    def test_both_breakpoints_can_disrupt_a_boundary_each(self) -> None:
        variant = tra("chr1", 95_000, "chr2", 505_000)
        prediction = predict_tad_disruption(variant, [self.CHR1_BOUNDARY, self.CHR2_BOUNDARY])
        assert prediction.n_boundaries_disrupted == 2
        assert prediction.disrupted_chroms == ["chr1", "chr2"]
        # 2 boundaries x 0.4 x the TRA weight 0.95
        assert prediction.disruption_score == pytest.approx(0.76)

    def test_a_boundary_near_both_breakpoints_counts_once(self) -> None:
        variant = tra("chr1", 95_000, "chr1", 104_000)  # both within 10 kb of the midpoint 100,000
        prediction = predict_tad_disruption(variant, [self.CHR1_BOUNDARY])
        assert prediction.n_boundaries_disrupted == 1
        assert prediction.disrupted_boundaries == [100_000]

    def test_a_boundary_on_the_first_chromosome_is_not_disrupted_by_the_second_breakpoint(self) -> None:
        # chr2:100,000 is nowhere near chr2's boundary, and the chr1 boundary at 100,000 is on another chromosome.
        prediction = predict_tad_disruption(
            tra("chr1", 2_000_000, "chr2", 100_000), [self.CHR1_BOUNDARY, self.CHR2_BOUNDARY]
        )
        assert prediction.n_boundaries_disrupted == 0
        assert prediction.disruption_score == 0.0

    def test_the_window_is_inclusive_at_each_breakpoint(self) -> None:
        boundary = self.CHR2_BOUNDARY  # midpoint 500,000
        assert predict_tad_disruption(tra("chr1", 7, "chr2", 510_000), [boundary]).n_boundaries_disrupted == 1
        assert predict_tad_disruption(tra("chr1", 7, "chr2", 510_001), [boundary]).n_boundaries_disrupted == 0
        assert predict_tad_disruption(tra("chr1", 7, "chr2", 490_000), [boundary]).n_boundaries_disrupted == 1
        assert predict_tad_disruption(tra("chr1", 7, "chr2", 489_999), [boundary]).n_boundaries_disrupted == 0

    def test_a_single_breakend_is_checked_at_its_one_breakpoint(self) -> None:
        single = {"chrom": "chr1", "start": 95_000, "end": 95_000, "sv_type": "BND"}
        prediction = predict_tad_disruption(single, [self.CHR1_BOUNDARY, self.CHR2_BOUNDARY])
        assert prediction.n_boundaries_disrupted == 1 and prediction.variant_pos2 is None

    def test_interval_variants_are_unchanged(self) -> None:
        deletion = {"chrom": "chr1", "start": 50_000, "end": 150_000, "sv_type": "DEL"}
        prediction = predict_tad_disruption(deletion, [self.CHR1_BOUNDARY, self.CHR2_BOUNDARY])
        assert prediction.n_boundaries_disrupted == 1 and prediction.disrupted_chroms == ["chr1"]
        assert (prediction.variant_chrom2, prediction.variant_pos2) == ("", None)

    def test_the_functional_impact_of_a_translocation_includes_the_tad_at_its_second_breakpoint(self) -> None:
        impact = predict_functional_impact(
            tra("chr1", 2_000_000, "chr2", 505_000), [], tad_boundaries=[self.CHR2_BOUNDARY]
        )
        assert impact.tad_disrupted is True and impact.impact_type == "tad_disruption"


# ---------------------------------------------------------------------------
# The other places a translocation's second breakpoint was ignored
# ---------------------------------------------------------------------------


GENES = [
    {"chrom": "chr1", "start": 500, "end": 900, "name": "GENE_ON_CHR1", "is_coding": True},
    {"chrom": "chr2", "start": 4_000, "end": 6_000, "name": "GENE_ON_CHR2", "is_coding": True},
]


class TestMergeSvCallsets:
    def test_a_junction_from_two_samples_merges_in_either_orientation(self) -> None:
        callsets = [
            [tra("chr1", 1_000, "chr2", 5_000, sample="s1")],
            [tra("chr2", 5_010, "chr1", 995, sample="s2")],
        ]
        [merged] = merge_sv_callsets(callsets, max_breakpoint_distance=100)
        assert merged["n_samples"] == 2 and merged["samples"] == ["s1", "s2"]
        # medians: (1000, 995) -> 997 and (5000, 5010) -> 5005
        assert (merged["chrom"], merged["start"], merged["end"]) == ("chr1", 997, 997)
        assert (merged["chrom2"], merged["pos2"]) == ("chr2", 5_005)

    def test_other_partners_and_distant_breakpoints_stay_apart(self) -> None:
        callsets = [
            [tra("chr1", 1_000, "chr2", 5_000, sample="s1")],
            [tra("chr1", 1_000, "chr3", 5_000, sample="s2"), tra("chr1", 1_000, "chr2", 6_000, sample="s3")],
        ]
        assert len(merge_sv_callsets(callsets, max_breakpoint_distance=100)) == 3

    def test_a_translocation_does_not_merge_with_an_interval_variant(self) -> None:
        interval = {"chrom": "chr1", "start": 990, "end": 1_010, "sv_type": "TRA", "sample": "s2"}
        callsets = [[tra("chr1", 1_000, "chr2", 5_000, sample="s1")], [dict(interval, sv_type="DEL")]]
        assert len(merge_sv_callsets(callsets, min_overlap=0.0, max_breakpoint_distance=100)) == 2

    def test_interval_consensus_has_no_partner_keys(self) -> None:
        callsets = [[{"chrom": "chr1", "start": 10_000, "end": 15_000, "sv_type": "DEL", "sample": "s1"}]]
        [merged] = merge_sv_callsets(callsets)
        assert "chrom2" not in merged and "pos2" not in merged


class TestBlacklist:
    REGIONS = [("chr2", 4_000, 6_000)]

    def test_a_translocation_with_its_second_breakpoint_in_a_region_is_removed(self) -> None:
        # The first breakpoint (chr1:50) is in no region; the second (chr2:5,000) is inside chr2:4,000-6,000.
        kept, stats = apply_blacklist([tra("chr1", 50, "chr2", 5_000)], self.REGIONS)
        assert kept == [] and stats.filtered_count == 1

    def test_a_translocation_with_both_breakpoints_clear_is_kept(self) -> None:
        variants = [tra("chr1", 50, "chr3", 5_000), tra("chr1", 50, "chr2", 7_000)]
        kept, _ = apply_blacklist(variants, self.REGIONS)
        assert kept == variants

    def test_a_breakpoint_on_the_edge_of_a_region_is_not_inside_it(self) -> None:
        variants = [tra("chr1", 50, "chr2", 4_000), tra("chr1", 50, "chr2", 6_000)]
        kept, _ = apply_blacklist(variants, self.REGIONS)
        assert kept == variants

    def test_the_first_breakpoint_still_counts(self) -> None:
        kept, _ = apply_blacklist([tra("chr2", 5_000, "chr1", 50)], self.REGIONS)
        assert kept == []

    def test_a_breakpoint_is_only_tested_against_regions_of_its_own_chromosome(self) -> None:
        # chr1:5,000 is inside chr2's region by coordinate, not by chromosome.
        variants = [tra("chr1", 5_000, "chr3", 5_000)]
        assert apply_blacklist(variants, self.REGIONS)[0] == variants


class TestGeneAnnotation:
    def test_genes_at_both_breakpoints_are_reported(self) -> None:
        [annotated] = annotate_gene_overlap([tra("chr1", 700, "chr2", 5_000)], GENES)
        assert annotated["overlapping_genes"] == ["GENE_ON_CHR1", "GENE_ON_CHR2"]
        assert annotated["n_genes_affected"] == 2
        first, second = annotated["gene_overlaps"]
        assert (first.variant_chrom, first.variant_start, first.variant_end) == ("chr1", 700, 700)
        assert (second.variant_chrom, second.variant_start, second.variant_end) == ("chr2", 5_000, 5_000)
        assert second.relationship == "contained_in"

    def test_a_gene_holding_both_breakpoints_is_listed_once(self) -> None:
        [annotated] = annotate_gene_overlap([tra("chr2", 4_500, "chr2", 5_500)], GENES)
        assert annotated["overlapping_genes"] == ["GENE_ON_CHR2"]

    def test_regulatory_elements_at_the_second_breakpoint_are_reported(self) -> None:
        elements = [{"chrom": "chr2", "start": 4_900, "end": 5_100, "name": "enhancer", "feature_type": "enhancer"}]
        [annotated] = annotate_regulatory_overlap([tra("chr1", 700, "chr2", 5_000)], elements)
        assert annotated["overlapping_regulatory"] == ["enhancer"]
        assert annotated["regulatory_types"] == {"enhancer"}

    def test_interval_variants_are_annotated_as_before(self) -> None:
        [annotated] = annotate_gene_overlap([{"chrom": "chr1", "start": 600, "end": 800, "sv_type": "DEL"}], GENES)
        assert annotated["overlapping_genes"] == ["GENE_ON_CHR1"]
        assert annotated["gene_overlaps"][0].overlap_bp == 200


class TestNearestGeneOfATranslocation:
    def test_the_nearer_of_the_two_breakpoints_decides(self) -> None:
        # chr1:2,000 is 1,100 bp from GENE_ON_CHR1 (ends at 900); chr2:6,150 is 150 bp from GENE_ON_CHR2 (ends 6,000).
        result = find_nearest_gene(tra("chr1", 2_000, "chr2", 6_150), GENES)
        assert (result["nearest_gene"], result["distance"], result["direction"]) == ("GENE_ON_CHR2", 150, "downstream")

    def test_a_breakpoint_inside_a_gene_is_an_overlap(self) -> None:
        result = find_nearest_gene(tra("chr1", 2_000, "chr2", 5_000), GENES)
        assert (result["nearest_gene"], result["distance"], result["direction"]) == ("GENE_ON_CHR2", 0, "overlapping")

    def test_the_first_breakpoint_wins_a_tie(self) -> None:
        result = find_nearest_gene(tra("chr1", 1_000, "chr2", 6_100), GENES)  # 100 bp from each gene
        assert result["nearest_gene"] == "GENE_ON_CHR1"


class TestFunctionalImpact:
    def test_a_gene_at_the_second_breakpoint_makes_a_translocation_a_fusion_candidate(self) -> None:
        impact = predict_functional_impact(tra("chr1", 2_000, "chr2", 5_000), GENES)
        assert impact.affected_genes == ["GENE_ON_CHR2"]
        assert (impact.impact_level, impact.impact_type) == ("HIGH", "gene_fusion_candidate")

    def test_genes_at_both_breakpoints_are_both_affected(self) -> None:
        impact = predict_functional_impact(tra("chr1", 700, "chr2", 5_000), GENES)
        assert impact.affected_genes == ["GENE_ON_CHR1", "GENE_ON_CHR2"]

    def test_a_translocation_between_gene_free_regions_is_benign(self) -> None:
        impact = predict_functional_impact(tra("chr1", 2_000, "chr2", 9_000), GENES)
        assert impact.affected_genes == [] and impact.impact_type == "benign"


# ---------------------------------------------------------------------------
# VCF input: the second breakpoint of a BND or TRA record
# ---------------------------------------------------------------------------

HEADER = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"


def _vcf(path: Path, *records: str) -> str:
    path.write_text(HEADER + "".join(record + "\n" for record in records))
    return str(path)


class TestVcfTranslocations:
    @pytest.mark.parametrize("alt", ["N[chr2:5001[", "N]chr2:5001]", "]chr2:5001]N", "[chr2:5001[N"])
    def test_the_mate_of_a_bnd_alt_is_the_second_breakpoint(self, tmp_path: Path, alt: str) -> None:
        path = _vcf(tmp_path / "a.vcf", f"chr1\t1001\tbnd_1\tN\t{alt}\t50\tPASS\tSVTYPE=BND")
        [record] = _parse_vcf_basic(path)
        # VCF positions are 1-based: POS 1001 is 0-based 1000 and the mate 5001 is 5000.
        assert (record["chrom"], record["start"], record["end"]) == ("chr1", 1_000, 1_000)
        assert (record["chrom2"], record["pos2"]) == ("chr2", 5_000)

    def test_chr2_with_end_is_the_second_breakpoint_of_a_tra_record(self, tmp_path: Path) -> None:
        path = _vcf(tmp_path / "a.vcf", "chr1\t1001\t.\tN\t<TRA>\t50\tPASS\tSVTYPE=TRA;CHR2=chr2;END=5001")
        [record] = _parse_vcf_basic(path)
        assert (record["start"], record["end"], record["chrom2"], record["pos2"]) == (1_000, 1_000, "chr2", 5_000)

    def test_pos2_takes_precedence_over_end(self, tmp_path: Path) -> None:
        path = _vcf(tmp_path / "a.vcf", "chr1\t1001\t.\tN\t<TRA>\t50\tPASS\tSVTYPE=BND;CHR2=chr2;POS2=7001;END=9999")
        [record] = _parse_vcf_basic(path)
        assert (record["chrom2"], record["pos2"]) == ("chr2", 7_000)

    def test_a_single_breakend_has_no_partner(self, tmp_path: Path) -> None:
        path = _vcf(tmp_path / "a.vcf", "chr1\t1001\t.\tN\tN.\t50\tPASS\tSVTYPE=BND")
        [record] = _parse_vcf_basic(path)
        assert (record["start"], record["end"]) == (1_000, 1_000)
        assert "chrom2" not in record and "pos2" not in record

    def test_a_contig_name_with_a_colon_is_kept(self, tmp_path: Path) -> None:
        path = _vcf(tmp_path / "a.vcf", "chr1\t1001\t.\tN\tN[HLA-A*01:01:200[\t50\tPASS\tSVTYPE=BND")
        [record] = _parse_vcf_basic(path)
        assert (record["chrom2"], record["pos2"]) == ("HLA-A*01:01", 199)

    def test_interval_records_are_unchanged(self, tmp_path: Path) -> None:
        path = _vcf(tmp_path / "a.vcf", "chr1\t1001\t.\tN\t<DEL>\t50\tPASS\tSVTYPE=DEL;END=5001")
        [record] = _parse_vcf_basic(path)
        assert (record["start"], record["end"]) == (1_000, 5_001) and "chrom2" not in record

    def test_pysam_and_the_text_parser_read_breakends_alike(self, tmp_path: Path) -> None:
        pytest.importorskip("pysam")  # not installed in every environment; the text parser above always runs
        header = (
            "##fileformat=VCFv4.2\n##contig=<ID=chr1,length=100000>\n##contig=<ID=chr2,length=100000>\n"
            '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="SV type">\n'
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        )
        records = (
            "chr1\t1001\tb1\tN\tN[chr2:5001[\t50\tPASS\tSVTYPE=BND",
            "chr2\t5011\tb2\tN\t]chr1:1011]N\t40\tPASS\tSVTYPE=BND",
            "chr1\t3001\tb3\tN\tN.\t30\tPASS\tSVTYPE=BND",
        )
        path = tmp_path / "breakends.vcf"
        path.write_text(header + "".join(record + "\n" for record in records))

        keys = ("chrom", "start", "end", "sv_type", "chrom2", "pos2")
        via_pysam = [{key: v.get(key) for key in keys} for v in _load_vcf_variants(str(path))]
        via_text = [{key: v.get(key) for key in keys} for v in _parse_vcf_basic(str(path))]
        assert (
            via_pysam
            == via_text
            == [
                {"chrom": "chr1", "start": 1_000, "end": 1_000, "sv_type": "BND", "chrom2": "chr2", "pos2": 5_000},
                {"chrom": "chr2", "start": 5_010, "end": 5_010, "sv_type": "BND", "chrom2": "chr1", "pos2": 1_010},
                {"chrom": "chr1", "start": 3_000, "end": 3_000, "sv_type": "BND", "chrom2": None, "pos2": None},
            ]
        )

    def test_survivor_merge_joins_the_two_mates_of_a_junction_written_from_opposite_sides(self, tmp_path: Path) -> None:
        first = _vcf(tmp_path / "delly.vcf", "chr1\t1001\t.\tN\tN[chr2:5001[\t80\tPASS\tSVTYPE=BND")
        second = _vcf(tmp_path / "manta.vcf", "chr2\t5011\t.\tN\t]chr1:1011]N\t40\tPASS\tSVTYPE=BND")
        other = _vcf(tmp_path / "gridss.vcf", "chr1\t1001\t.\tN\tN[chr3:5001[\t40\tPASS\tSVTYPE=BND")

        merged, stats = survivor_merge([first, second, other], max_distance=100, min_callers=2)

        assert stats.n_input_variants == 3
        [variant] = merged  # the chr3 partner is a different event and, alone, below min_callers
        assert variant.n_callers == 2 and isinstance(variant, MergedVariant)
        # medians: (1000, 1010) -> 1005 and (5000, 5010) -> 5005
        assert (variant.chrom, variant.start, variant.chrom2, variant.pos2) == ("chr1", 1_005, "chr2", 5_005)
