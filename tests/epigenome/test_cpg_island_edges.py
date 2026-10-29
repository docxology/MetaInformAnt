"""CpG island edges are resolved to the base, not to the scanning window.

Defect: ``find_cpg_islands`` reported the union of the qualifying windows, so an island's boundaries were only as
sharp as the window grid. A window that was one-quarter flank still qualified, which pushed the reported start and end
up to a quarter of a window outside the island, and the grid could also stop short of an edge by up to one step. Sites
in the padding were counted as inside the island by ``identify_cpg_islands``. A short island could vanish altogether:
the padded union of its windows failed the GC criterion as a whole even though the island itself, and several of the
windows, met every criterion.

The merged stretches are now refined to the outermost CpGs: an island starts at the C of a CpG and ends at the G of a
CpG, and all three Gardiner-Garden & Frommer criteria hold for the island as reported.

The sequences are planted: a CpG-rich body that starts and ends with ``CG`` between flanks of A and T only (no CpG, so
the edges are sharp and the right answer is known exactly).

REAL IMPLEMENTATION.
"""

from __future__ import annotations

import random
from typing import Dict, List, Tuple

import pytest

from metainformant.epigenome.assays.methylation import (
    MethylationSite,
    _cpg_statistics,
    find_cpg_islands,
    identify_cpg_islands,
)

# Six bases, two CpGs, G + C = 2/3: every window that is at least three quarters island qualifies
_RICH_UNIT = "CGCGAT"

# Twenty-five bases with G + C = 13/25 = 0.52, a thin margin over 0.5: a window with eight or more flank bases fails
# the GC criterion, so only windows wholly inside the island qualify and the grid stops short of the edges
_THIN_UNIT = "CGCGCGCGCGCGC" + "AATTAATTAATT"


def _flank(rng: random.Random, length: int) -> str:
    """A and T only, so the flank holds no CpG and has a GC fraction of 0."""
    return "".join(rng.choice("AT") for _ in range(length))


def _rich_body(length: int) -> str:
    """A CpG-rich body of ``length`` bases (``length`` = 2 mod 6) that starts and ends with CG."""
    assert (length - 2) % len(_RICH_UNIT) == 0
    return _RICH_UNIT * ((length - 2) // len(_RICH_UNIT)) + "CG"


def _planted(left: int, body: str, right: int, seed: int = 0) -> Tuple[str, Tuple[int, int]]:
    rng = random.Random(seed)
    sequence = _flank(rng, left) + body + _flank(rng, right)
    return sequence, (left, left + len(body))


def _edges(islands: List[Dict[str, object]]) -> List[Tuple[object, object]]:
    return [(island["start"], island["end"]) for island in islands]


class TestPlantedIslandsKeepTheirEdges:
    @pytest.mark.parametrize("left", [333, 487, 251, 640])
    @pytest.mark.parametrize("length", [308, 410, 602, 998])
    @pytest.mark.parametrize("window,step", [(200, 50), (200, 25), (200, 100), (100, 20), (200, 1)])
    def test_the_reported_edges_are_the_planted_edges(self, left: int, length: int, window: int, step: int) -> None:
        sequence, planted = _planted(left, _rich_body(length), 300 + left % 97)

        islands = find_cpg_islands({"chr1": sequence}, window_size=window, step_size=step)

        assert _edges(islands) == [planted]

    def test_the_planted_island_is_not_inflated_to_the_window_grid(self) -> None:
        # Windows that are up to a quarter flank still qualify, so the merged windows run past the island on both
        # sides; the island must end where the CpGs end
        sequence, (start, end) = _planted(487, _rich_body(410), 333)

        (island,) = find_cpg_islands({"chr1": sequence}, window_size=200, step_size=50)

        assert island["start"] == start == 487
        assert island["end"] == end == 897
        assert island["length"] == 410
        assert sequence[start : start + 2] == "CG" and sequence[end - 2 : end] == "CG"

    def test_statistics_describe_the_island_not_the_padded_windows(self) -> None:
        sequence, (start, end) = _planted(487, _rich_body(410), 333)

        (island,) = find_cpg_islands({"chr1": sequence}, window_size=200, step_size=50)

        gc_content, cpg_ratio, cpg_count = _cpg_statistics(sequence[start:end])
        assert island["gc_content"] == gc_content
        assert island["cpg_ratio"] == cpg_ratio
        assert island["cpg_count"] == cpg_count == 2 * 68 + 1
        # The body is two thirds G + C; the padded windows would have reported a lower fraction
        assert island["gc_content"] == pytest.approx((2 * 68 * 2 + 2) / 410)

    def test_two_islands_keep_their_own_edges(self) -> None:
        rng = random.Random(3)
        first, second = _rich_body(410), _rich_body(308)
        sequence = _flank(rng, 400) + first + _flank(rng, 500) + second + _flank(rng, 400)

        islands = find_cpg_islands({"chr1": sequence}, window_size=200, step_size=50)

        assert _edges(islands) == [(400, 810), (1310, 1618)]

    def test_islands_on_both_strands_are_the_same_islands(self) -> None:
        # CpG is its own reverse complement, so the reverse complement holds the mirror image of every island
        sequence, (start, end) = _planted(487, _rich_body(410), 333)
        complement = {"A": "T", "T": "A", "C": "G", "G": "C"}
        reverse_complement = "".join(complement[base] for base in reversed(sequence))

        (island,) = find_cpg_islands({"chr1": reverse_complement}, window_size=200, step_size=50)

        assert (island["start"], island["end"]) == (len(sequence) - end, len(sequence) - start)

    def test_case_and_ambiguous_bases_in_the_flank_do_not_move_the_edges(self) -> None:
        sequence, planted = _planted(487, _rich_body(410), 333)
        masked = sequence[:300].lower() + "N" * 100 + sequence[400:].lower()

        assert _edges(find_cpg_islands({"chr1": masked}, window_size=200, step_size=50)) == [planted]


class TestIslandsAWindowGridMisses:
    def test_the_grid_stopping_short_of_an_edge_does_not_shorten_the_island(self) -> None:
        # G + C is 0.52 in the body, so any window with eight flank bases fails: only windows wholly inside the island
        # qualify, and the 50 bp grid leaves 37 bp uncovered at the left edge and 15 at the right
        body = _THIN_UNIT * 16 + "CG"
        sequence, planted = _planted(1013, body, 900)

        islands = find_cpg_islands({"chr1": sequence}, window_size=200, step_size=50)

        assert _edges(islands) == [planted] == [(1013, 1415)]

    @pytest.mark.parametrize("left", [1000, 1021, 1037, 1049, 1050])
    def test_the_island_is_complete_wherever_it_sits_on_the_grid(self, left: int) -> None:
        body = _THIN_UNIT * 16 + "CG"
        sequence, planted = _planted(left, body, 900)

        islands = find_cpg_islands({"chr1": sequence}, window_size=200, step_size=50)

        assert _edges(islands) == [planted]

    def test_an_island_as_short_as_the_minimum_length_is_found(self) -> None:
        # A 200 bp island: the windows that overlap it are padded with flank, and their union is 275 bp with a GC
        # fraction below 0.5, so the padded stretch failed the criteria as a whole and the island was lost
        sequence, planted = _planted(365, _rich_body(200), 435)

        islands = find_cpg_islands({"chr1": sequence}, window_size=200, step_size=25)

        assert _edges(islands) == [planted] == [(365, 565)]

    def test_a_stretch_whose_cpg_bounded_island_is_too_short_keeps_its_window_edges(self) -> None:
        # The CpG-rich body is 206 bp and sits between 60 bp blocks of G + C without any CpG. With a minimum length
        # of 250 the CpG-bounded island is too short, but the merged windows (400 bp, G + C 0.545, observed/expected
        # CpG 2.3) meet every criterion as a whole, so they are reported with their window edges, not dropped
        block = "GGCCAT" * 10
        sequence = _flank(random.Random(0), 400) + block + _rich_body(206) + block + _flank(random.Random(1), 400)

        (island,) = find_cpg_islands({"chr1": sequence}, window_size=200, step_size=50, min_length=250)

        assert (island["start"], island["end"]) == (350, 750)
        assert island["gc_content"] >= 0.5 and island["cpg_ratio"] >= 0.6

    def test_a_stretch_without_cpg_is_reported_when_it_qualifies_with_no_cpg_criterion(self) -> None:
        # G + C is 8/9 but no C is followed by G; with the CpG criterion switched off the window edges stand
        sequence = "GGGGCCCCA" * 150

        (island,) = find_cpg_islands({"chr1": sequence}, min_cpg_ratio=0.0)

        assert (island["start"], island["end"]) == (0, len(sequence))
        assert island["cpg_count"] == 0


class TestSitesAreCountedInsideTheIslandOnly:
    def test_sites_in_the_padding_of_the_windows_are_not_island_sites(self) -> None:
        sequence, (start, end) = _planted(487, _rich_body(410), 333)
        sites = {
            "chr1": [
                MethylationSite("chr1", start - 40, 1, 10),  # in the flank, inside the padded windows
                MethylationSite("chr1", start - 1, 1, 10),
                MethylationSite("chr1", start, 9, 10),  # first base of the island
                MethylationSite("chr1", end - 1, 7, 10),  # last base of the island
                MethylationSite("chr1", end, 1, 10),  # end is exclusive
                MethylationSite("chr1", end + 40, 1, 10),
            ]
        }

        (island,) = identify_cpg_islands(sites, sequences={"chr1": sequence})

        assert (island["start"], island["end"]) == (start, end)
        assert island["sites"] == 2
        assert island["avg_methylation"] == pytest.approx(0.8)


def _random_sequence(rng: random.Random, length: int) -> str:
    """Blocks of different composition: AT-rich, GC-rich without CpG, mixed and CpG-rich."""
    parts: List[str] = []
    while sum(len(part) for part in parts) < length:
        kind = rng.choice(["at", "gc", "mix", "cpg"])
        size = rng.randint(20, 400)
        if kind == "cpg":
            parts.append("".join(rng.choice(["CG", "CG", "CGC", "AT", "G"]) for _ in range(size // 2)))
        else:
            alphabet = {"at": "AAATTTTC", "gc": "GGCCAT", "mix": "ACGT"}[kind]
            parts.append("".join(rng.choice(alphabet) for _ in range(size)))
    return "".join(parts)[:length]


def _naive_refinement(
    sequence: str, start: int, end: int, low: int, high: int, min_gc: float, min_ratio: float, min_length: int
) -> object:
    """The same refinement rule, scored by slicing and counting each candidate island instead of by prefix counts."""
    cpgs: List[int] = []
    found = sequence.find("CG", low, high)
    while found != -1:
        cpgs.append(found)
        found = sequence.find("CG", found + 2, high)

    def qualifies(first: int, last: int) -> bool:
        begin, finish = cpgs[first], cpgs[last] + 2
        if finish - begin < min_length:
            return False
        gc_content, cpg_ratio, _ = _cpg_statistics(sequence[begin:finish])
        return gc_content >= min_gc and cpg_ratio >= min_ratio

    inside = [k for k, position in enumerate(cpgs) if position >= start and position + 2 <= end]
    if not inside:
        return None
    first, last = inside[0], inside[-1]
    while True:
        wider_first = next((k for k in range(first) if qualifies(k, last)), first)
        wider_last = next((k for k in range(len(cpgs) - 1, last, -1) if qualifies(wider_first, k)), last)
        if (wider_first, wider_last) == (first, last):
            break
        first, last = wider_first, wider_last
    return (cpgs[first], cpgs[last] + 2) if qualifies(first, last) else None


class TestOnRandomSequence:
    PARAMETERS = [
        # min_gc, min_cpg_ratio, min_length, window, step
        (0.4, 0.4, 100, 100, 20),
        (0.5, 0.6, 200, 200, 50),
        (0.55, 0.65, 300, 200, 25),
        (0.5, 0.6, 200, 150, 150),
        (0.4, 0.6, 300, 200, 1),
    ]

    @pytest.mark.parametrize("seed", range(40))
    def test_every_island_meets_every_criterion_and_islands_do_not_overlap(self, seed: int) -> None:
        rng = random.Random(seed)
        sequence = _random_sequence(rng, rng.randint(300, 3000))
        min_gc, min_ratio, min_length, window, step = self.PARAMETERS[seed % len(self.PARAMETERS)]

        islands = find_cpg_islands(
            {"chr1": sequence},
            window_size=window,
            step_size=step,
            min_gc=min_gc,
            min_cpg_ratio=min_ratio,
            min_length=min_length,
        )

        previous_end = 0
        for island in islands:
            start, end = island["start"], island["end"]
            assert previous_end <= start < end <= len(sequence)
            previous_end = end
            gc_content, cpg_ratio, cpg_count = _cpg_statistics(sequence[start:end])
            assert end - start >= min_length
            assert gc_content >= min_gc and cpg_ratio >= min_ratio
            assert (island["length"], island["gc_content"], island["cpg_ratio"], island["cpg_count"]) == (
                end - start,
                gc_content,
                cpg_ratio,
                cpg_count,
            )

    @pytest.mark.parametrize("seed", range(40))
    def test_refinement_agrees_with_scoring_every_candidate_by_slicing(self, seed: int) -> None:
        from metainformant.epigenome.assays.methylation import _refine_island

        rng = random.Random(1000 + seed)
        sequence = _random_sequence(rng, rng.randint(300, 3000))
        min_gc, min_ratio, min_length, window, step = self.PARAMETERS[seed % len(self.PARAMETERS)]

        starts = list(range(0, len(sequence) - window + 1, step))
        if starts[-1] != len(sequence) - window:
            starts.append(len(sequence) - window)
        stretches: List[List[int]] = []
        for window_start in starts:
            gc_content, cpg_ratio, _ = _cpg_statistics(sequence[window_start : window_start + window])
            if gc_content < min_gc or cpg_ratio < min_ratio:
                continue
            if stretches and window_start <= stretches[-1][1]:
                stretches[-1][1] = window_start + window
            else:
                stretches.append([window_start, window_start + window])

        reach = min(window, step)
        for number, (start, end) in enumerate(stretches):
            following = stretches[number + 1][0] if number + 1 < len(stretches) else len(sequence)
            # Without the neighbouring island's end the room is only bounded by the reach and the next stretch
            low, high = max(0, start - reach), min(following, end + reach)
            assert _refine_island(sequence, start, end, low, high, min_gc, min_ratio, min_length) == _naive_refinement(
                sequence, start, end, low, high, min_gc, min_ratio, min_length
            )
