"""``merge_tracks``: a real intersection, and the other operations merged correctly.

* ``operation="intersection"`` was ``merged_track.data[chromosome] = all_features[:10]  # Placeholder``: the first ten
  features of all the tracks together, whether or not they overlapped anything.
* ``"mean"`` and ``"max"`` compared a feature with the end of the *last* feature of the group instead of the end of
  the group, so a feature inside a long earlier one started a new group.
* ``"union"`` merged the features of the input tracks in place, rewriting the caller's tracks.
* An unknown operation silently returned an empty track.

The intersection is checked against a position-by-position model. REAL IMPLEMENTATION.
"""

from __future__ import annotations

import random
from typing import Dict, List, Optional, Tuple

import pytest

from metainformant.epigenome.analysis.tracks import GenomicTrack, merge_tracks

Feature = Tuple[str, int, int, float]


def _track(features: List[Feature], name: str = "track") -> GenomicTrack:
    track = GenomicTrack(name=name)
    for chromosome, start, end, value in features:
        track.add_feature(chromosome, start, end, value=value)
    return track


def _spans(track: GenomicTrack, chromosome: str) -> List[Tuple[int, int, float]]:
    return [(f["start"], f["end"], f["value"]) for f in track.get_features(chromosome)]


class TestIntersection:
    def test_disjoint_tracks_have_an_empty_intersection(self) -> None:
        # The placeholder returned the first ten features of the two tracks, which share no base
        first = _track([("chr1", i * 100, i * 100 + 50, 1.0) for i in range(15)])
        second = _track([("chr1", 5000 + i * 100, 5000 + i * 100 + 50, 2.0) for i in range(15)])

        merged = merge_tracks([first, second], operation="intersection")

        assert merged.get_total_features() == 0
        assert merged.get_chromosomes() == []

    def test_only_the_shared_regions_are_kept(self) -> None:
        first = _track([("chr1", 0, 200, 1.0), ("chr1", 1000, 1200, 3.0)])
        second = _track([("chr1", 100, 300, 5.0), ("chr1", 1100, 1150, 7.0)])

        merged = merge_tracks([first, second], operation="intersection")

        # [100, 200) in both; [1100, 1150) in both. The value is the lower of the two tracks' values there.
        assert _spans(merged, "chr1") == [(100, 200, 1.0), (1100, 1150, 3.0)]

    def test_three_tracks_must_all_cover_a_region(self) -> None:
        first = _track([("chr1", 0, 100, 4.0)])
        second = _track([("chr1", 50, 150, 2.0)])
        third = _track([("chr1", 80, 200, 3.0)])

        merged = merge_tracks([first, second, third], operation="intersection")

        assert _spans(merged, "chr1") == [(80, 100, 2.0)]

    def test_a_chromosome_missing_from_one_track_has_no_intersection(self) -> None:
        first = _track([("chr1", 0, 100, 1.0), ("chr2", 0, 100, 1.0)])
        second = _track([("chr1", 50, 150, 1.0)])

        merged = merge_tracks([first, second], operation="intersection")

        assert merged.get_chromosomes() == ["chr1"]
        assert _spans(merged, "chr1") == [(50, 100, 1.0)]

    def test_touching_regions_with_the_same_value_are_joined(self) -> None:
        first = _track([("chr1", 0, 100, 2.0), ("chr1", 100, 200, 2.0)])
        second = _track([("chr1", 50, 150, 2.0)])

        merged = merge_tracks([first, second], operation="intersection")

        assert _spans(merged, "chr1") == [(50, 150, 2.0)]

    def test_overlapping_features_inside_a_track_use_their_highest_value(self) -> None:
        # first: 1 over [0, 100) with a 9 over [20, 60); second: 5 over [0, 100)
        first = _track([("chr1", 0, 100, 1.0), ("chr1", 20, 60, 9.0)])
        second = _track([("chr1", 0, 100, 5.0)])

        merged = merge_tracks([first, second], operation="intersection")

        assert _spans(merged, "chr1") == [(0, 20, 1.0), (20, 60, 5.0), (60, 100, 1.0)]

    def test_a_single_track_intersects_with_itself(self) -> None:
        only = _track([("chr1", 0, 100, 1.0), ("chr1", 50, 150, 3.0)])

        merged = merge_tracks([only], operation="intersection")

        assert _spans(merged, "chr1") == [(0, 50, 1.0), (50, 150, 3.0)]

    def test_features_without_extent_are_ignored(self) -> None:
        first = _track([("chr1", 10, 10, 9.0), ("chr1", 0, 100, 1.0)])
        second = _track([("chr1", 10, 10, 9.0), ("chr1", 0, 100, 2.0)])

        assert _spans(merge_tracks([first, second], operation="intersection"), "chr1") == [(0, 100, 1.0)]

    def test_the_input_tracks_are_not_changed(self) -> None:
        first = _track([("chr1", 0, 200, 1.0)])
        second = _track([("chr1", 100, 300, 5.0)])

        merge_tracks([first, second], operation="intersection")

        assert _spans(first, "chr1") == [(0, 200, 1.0)]
        assert _spans(second, "chr1") == [(100, 300, 5.0)]

    def test_matches_a_position_by_position_model(self) -> None:
        rng = random.Random(20261027)

        for _ in range(150):
            tracks: List[GenomicTrack] = []
            covering: List[Dict[Tuple[str, int], float]] = (
                []
            )  # per track: highest value covering each (chromosome, base)
            for _track_index in range(rng.randrange(1, 4)):
                features: List[Feature] = []
                for _feature_index in range(rng.randrange(0, 6)):
                    start = rng.randrange(0, 40)
                    features.append(
                        (rng.choice(["chr1", "chr2"]), start, start + rng.randrange(1, 15), float(rng.randrange(1, 6)))
                    )
                tracks.append(_track(features))
                best: Dict[Tuple[str, int], float] = {}
                for chromosome, start, end, value in features:
                    for base in range(start, end):
                        best[(chromosome, base)] = max(best.get((chromosome, base), value), value)
                covering.append(best)

            expected: Dict[str, List[Tuple[int, int, float]]] = {}
            for chromosome in ("chr1", "chr2"):
                previous: Optional[Tuple[int, float]] = None  # (base, value) of the run being extended
                for base in range(0, 60):
                    key = (chromosome, base)
                    if all(key in best for best in covering):
                        value = min(best[key] for best in covering)
                        runs = expected.setdefault(chromosome, [])
                        if previous is not None and previous == (base - 1, value):
                            runs[-1] = (runs[-1][0], base + 1, value)
                        else:
                            runs.append((base, base + 1, value))
                        previous = (base, value)
                    else:
                        previous = None

            merged = merge_tracks(tracks, operation="intersection")

            assert {chromosome: _spans(merged, chromosome) for chromosome in merged.get_chromosomes()} == expected


class TestMeanAndMaxGroupsFollowTheGroupEnd:
    def test_a_feature_inside_an_earlier_long_one_stays_in_its_group(self) -> None:
        # [50, 60) lies inside [0, 100); comparing it with the last feature [10, 20) started a new group
        track = _track([("chr1", 0, 100, 1.0), ("chr1", 10, 20, 2.0), ("chr1", 50, 60, 4.0)])

        assert _spans(merge_tracks([track], operation="max"), "chr1") == [(0, 100, 4.0)]
        assert _spans(merge_tracks([track], operation="mean"), "chr1") == [(0, 100, pytest.approx(7 / 3))]

    def test_a_chain_of_overlaps_forms_one_group_and_a_gap_ends_it(self) -> None:
        # [28, 40) starts after [20, 25) ends but before the group [0, 30) does
        track = _track(
            [
                ("chr1", 0, 10, 1.0),
                ("chr1", 5, 30, 2.0),
                ("chr1", 20, 25, 3.0),
                ("chr1", 28, 40, 4.0),
                ("chr1", 100, 110, 5.0),
            ]
        )

        assert _spans(merge_tracks([track], operation="max"), "chr1") == [(0, 40, 4.0), (100, 110, 5.0)]

    def test_features_of_several_tracks_are_grouped_together(self) -> None:
        first = _track([("chr1", 0, 100, 10.0)])
        second = _track([("chr1", 50, 150, 20.0)])

        assert _spans(merge_tracks([first, second], operation="mean"), "chr1") == [(0, 150, 15.0)]


class TestInputsAreNotModified:
    def test_union_leaves_the_input_tracks_alone(self) -> None:
        first = _track([("chr1", 0, 100, 1.0)])
        second = _track([("chr1", 50, 300, 2.0)])
        first.data["chr1"][0]["name"] = "kept"

        merged = merge_tracks([first, second], operation="union")

        assert _spans(merged, "chr1") == [(0, 300, 2.0)]
        # the old merge extended the first feature of the first input track in place
        assert first.data["chr1"] == [{"start": 0, "end": 100, "value": 1.0, "name": "kept"}]
        assert second.data["chr1"] == [{"start": 50, "end": 300, "value": 2.0}]

    @pytest.mark.parametrize("operation", ["mean", "max"])
    def test_grouping_leaves_the_input_tracks_alone(self, operation: str) -> None:
        first = _track([("chr1", 0, 100, 1.0)])
        second = _track([("chr1", 50, 300, 2.0)])

        merge_tracks([first, second], operation=operation)

        assert _spans(first, "chr1") == [(0, 100, 1.0)]
        assert _spans(second, "chr1") == [(50, 300, 2.0)]


class TestOperations:
    def test_an_unknown_operation_is_an_error(self) -> None:
        with pytest.raises(ValueError, match="Unsupported merge operation"):
            merge_tracks([_track([("chr1", 0, 10, 1.0)])], operation="sum")

    def test_no_tracks_give_an_empty_track(self) -> None:
        assert merge_tracks([], operation="intersection").get_total_features() == 0
