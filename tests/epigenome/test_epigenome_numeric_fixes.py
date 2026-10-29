"""Numerical defects found in the epigenome package while fixing its open items.

* ``peak_summit_refinement`` moved the summit the wrong way. The offset of the parabola's vertex had the opposite
  sign, so a summit with the higher neighbour on its right was refined to the left of the highest bin. Its ``window``
  argument was documented as the window "to use for fitting" while it only bounds where the highest bin is looked for.
* ``_poisson_pvalue`` computed ``1 - cdf``, which cannot represent a p-value below about 1e-16: a count of 20 over a
  background of 1 (p = 1.6e-19) was exactly 0.0, so every strong peak tied at p = 0, a score of 300 and a q-value of 0.
* Without SciPy, ``_fisher_or_chi2`` used a chi-squared statistic that squares the deviation from the expected
  overlap, so a state *depleted* in an annotation got as small a p-value as an enriched one. With no annotation to test,
  ``compute_state_enrichment`` reported a best fold of 0.0 and a best p-value of 1.0 for a test that never happened.
* ``compute_state_transition_rates`` gave a state nothing was seen to leave a uniform row, and with it a self-transition
  rate, an entropy and a segment length that no observation supported; an assignment outside ``[0, n_states)`` was left
  out of the counts but counted in ``n_transitions``.
* ``compare_tracks`` correlated the first features of the two tracks by list position, whatever their coordinates.
* ``identify_tss_enrichment`` assumed a 3 Gb genome for every species, and left the peak length out of the chance
  expectation.

Oracles: ``numpy.polyfit`` for the parabola, ``scipy.stats.poisson`` and ``scipy.stats.fisher_exact``, exact integer
arithmetic for the hypergeometric tail, and per-base ``numpy`` arrays for the track correlation.

REAL IMPLEMENTATION.
"""

from __future__ import annotations

import json
import math
import os
import random
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Any, List, Tuple

import numpy as np
import pytest

import metainformant
from metainformant.epigenome.analysis.tracks import GenomicTrack, compare_tracks
from metainformant.epigenome.assays.atacseq import ATACPeak, identify_tss_enrichment
from metainformant.epigenome.chromatin_state.state_learning import (
    compute_state_enrichment,
    compute_state_transition_rates,
)
from metainformant.epigenome.peak_calling.peak_detection import (
    _poisson_pvalue,
    call_peaks_simple,
    peak_summit_refinement,
)

# ---------------------------------------------------------------------------
# peak_summit_refinement
# ---------------------------------------------------------------------------


def _vertex(left: float, centre: float, right: float, index: int) -> float:
    """Position of the vertex of the parabola through three consecutive points, by a least-squares polynomial fit."""
    a, b, _ = np.polyfit([index - 1, index, index + 1], [left, centre, right], 2)
    return float(-b / (2 * a))


class TestSummitRefinement:
    def test_a_higher_right_neighbour_pulls_the_summit_right(self) -> None:
        signal = [0.0] * 10
        signal[3], signal[4], signal[5] = 1.0, 3.0, 2.0

        refined = peak_summit_refinement(signal, {"start": 0, "end": 10, "summit": 4}, window=5)

        assert refined["refined_summit_position"] == pytest.approx(_vertex(1.0, 3.0, 2.0, 4), abs=5e-5)
        assert refined["refined_summit_position"] > 4  # it used to be 3.8333
        assert refined["summit"] == 4

    def test_a_higher_left_neighbour_pulls_the_summit_left(self) -> None:
        signal = [0.0] * 10
        signal[3], signal[4], signal[5] = 2.0, 3.0, 1.0

        refined = peak_summit_refinement(signal, {"start": 0, "end": 10, "summit": 4}, window=5)

        assert refined["refined_summit_position"] == pytest.approx(_vertex(2.0, 3.0, 1.0, 4), abs=5e-5)
        assert refined["refined_summit_position"] < 4

    @pytest.mark.parametrize("seed", range(20))
    def test_the_refined_position_is_the_vertex_of_the_parabola_through_the_highest_bin(self, seed: int) -> None:
        rng = random.Random(seed)
        index = rng.randint(5, 25)
        centre = rng.uniform(5.0, 50.0)
        left, right = rng.uniform(0.1, centre - 0.1), rng.uniform(0.1, centre - 0.1)
        signal = [rng.uniform(0.0, 0.05) for _ in range(31)]
        signal[index - 1], signal[index], signal[index + 1] = left, centre, right

        refined = peak_summit_refinement(signal, {"start": 0, "end": 31, "summit": index}, window=10)

        assert refined["refined_summit_position"] == pytest.approx(_vertex(left, centre, right, index), abs=5e-5)
        assert abs(refined["refined_summit_position"] - index) <= 0.5
        assert refined["summit"] == int(round(refined["refined_summit_position"]))

    def test_the_window_bounds_the_search_for_the_highest_bin_only(self) -> None:
        signal = [0.0] * 40
        signal[10], signal[11], signal[12] = 4.0, 6.0, 5.0
        signal[30] = 100.0  # far higher, but outside the window around the summit at 11

        refined = peak_summit_refinement(signal, {"start": 0, "end": 40, "summit": 11}, window=3)

        assert refined["refined_summit_position"] == pytest.approx(_vertex(4.0, 6.0, 5.0, 11), abs=5e-5)

    def test_a_highest_bin_at_the_edge_of_the_searched_region_is_not_interpolated(self) -> None:
        edge = [5.0, 3.0, 1.0, 0.0, 0.0]

        refined = peak_summit_refinement(edge, {"start": 0, "end": 5, "summit": 0}, window=2)

        assert refined["refined_summit_position"] == 0.0 and refined["summit"] == 0

    def test_a_plateau_is_refined_to_the_vertex_through_its_first_bin(self) -> None:
        # The highest bin is the first of the plateau; the parabola through (0, 0), (1, 2), (2, 2) peaks at 1.5
        refined = peak_summit_refinement([0.0, 2.0, 2.0, 2.0, 0.0], {"start": 0, "end": 5, "summit": 2}, window=2)

        assert refined["refined_summit_position"] == pytest.approx(_vertex(0.0, 2.0, 2.0, 1), abs=5e-5) == 1.5


# ---------------------------------------------------------------------------
# Poisson tail
# ---------------------------------------------------------------------------


class TestPoissonPValue:
    @pytest.mark.parametrize(
        "observed,expected,want",
        [
            (20, 1.0, 1.5875276010732616e-19),
            (60, 1.0, 4.494768303810778e-83),
            (45, 2.5, 1.0e-45),  # order of magnitude only; the exact value comes from SciPy below
        ],
    )
    def test_a_tiny_p_value_is_not_rounded_to_zero(self, observed: int, expected: float, want: float) -> None:
        got = _poisson_pvalue(observed, expected)

        assert got > 0.0
        if observed != 45:
            assert got == pytest.approx(want, rel=1e-9)

    def test_it_agrees_with_scipy_where_one_minus_the_cdf_still_works(self) -> None:
        stats = pytest.importorskip("scipy.stats")
        for observed in (0, 1, 2, 5, 12, 40, 150):
            for expected in (0.05, 1.0, 7.5, 60.0):
                want = float(stats.poisson.sf(observed - 1, expected)) if observed > 0 else 1.0
                assert _poisson_pvalue(observed, expected) == pytest.approx(want, rel=1e-8, abs=1e-300)

    def test_no_background_means_certain_enrichment_for_any_count(self) -> None:
        assert _poisson_pvalue(3, 0.0) == 0.0
        assert _poisson_pvalue(0, 0.0) == 1.0

    def test_strong_peaks_keep_distinct_scores_and_q_values(self) -> None:
        signal = [1.0] * 1000
        for position in range(300, 340):
            signal[position] = 40.0
        for position in range(700, 740):
            signal[position] = 60.0

        peaks = call_peaks_simple(signal, threshold=5.0, min_length=20, merge_distance=5, chrom="chr1")

        assert [(p["start"], p["end"]) for p in peaks] == [(300, 340), (700, 740)]
        for peak in peaks:
            assert 0.0 < peak["p_value"] < 1e-17  # below what 1 - cdf can represent (about 1e-16)
            assert peak["score"] == pytest.approx(-math.log10(peak["p_value"]), abs=1e-3)
            assert peak["score"] < 300.0  # it used to be 300 for both: p was 0.0
            assert 0.0 < peak["q_value"] < 1e-16
        assert peaks[1]["p_value"] < peaks[0]["p_value"]
        assert peaks[1]["score"] > peaks[0]["score"] + 5


# ---------------------------------------------------------------------------
# Enrichment p-values without SciPy
# ---------------------------------------------------------------------------


def _exact_enrichment_p(a: int, b: int, c: int, d: int) -> float:
    """``P(X >= a)`` for the hypergeometric table ``[[a, b], [c, d]]``, in exact integer arithmetic."""
    row, column, total = a + b, a + c, a + b + c + d
    ways = math.comb(total, row)
    tail = sum(math.comb(column, k) * math.comb(total - column, row - k) for k in range(a, min(row, column) + 1))
    return tail / ways


class TestFisherFallback:
    @pytest.mark.parametrize("seed", range(10))
    def test_the_pure_python_tail_is_the_exact_hypergeometric_tail(self, seed: int) -> None:
        from metainformant.epigenome.chromatin_state.state_learning import _fisher_exact_greater

        rng = random.Random(seed)
        for _ in range(40):
            total = rng.choice([2, 7, 30, 200, 1500])
            row, column = rng.randint(0, total), rng.randint(0, total)
            a = rng.randint(max(0, row + column - total), min(row, column))
            table = (a, row - a, column - a, total - row - column + a)

            assert _fisher_exact_greater(*table) == pytest.approx(_exact_enrichment_p(*table), rel=1e-9, abs=1e-15)

    def test_it_agrees_with_scipy_on_a_depleted_and_an_enriched_table(self) -> None:
        from metainformant.epigenome.chromatin_state.state_learning import _fisher_exact_greater

        stats = pytest.importorskip("scipy.stats")
        for table in ((2, 48, 48, 2), (40, 10, 10, 40), (0, 30, 30, 0), (25, 25, 25, 25)):
            want = float(stats.fisher_exact([[table[0], table[1]], [table[2], table[3]]], alternative="greater").pvalue)
            assert _fisher_exact_greater(*table) == pytest.approx(want, rel=1e-9)

    def test_an_empty_table_is_not_enriched(self) -> None:
        from metainformant.epigenome.chromatin_state.state_learning import _fisher_exact_greater

        assert _fisher_exact_greater(0, 0, 0, 0) == 1.0

    def test_enrichment_without_scipy_distinguishes_depletion_from_enrichment(self) -> None:
        # 100 bins; state 0 holds bins 0-49, state 1 bins 50-99. Annotation A covers 48 bins of state 1 and 2 of
        # state 0 (state 0 is depleted, state 1 enriched)
        script = r"""
import json, math, sys
sys.modules["scipy"] = None  # makes ``import scipy`` fail, as if it were not installed
from metainformant.epigenome.chromatin_state import state_learning

assert state_learning.HAS_SCIPY is False
assignments = [0] * 50 + [1] * 50
annotation = list(range(0, 2)) + list(range(52, 100))
result = state_learning.compute_state_enrichment(assignments, {"A": annotation}, 2)
entries = result["enrichment_matrix"]
depleted, enriched = entries["state_0__A"]["p_value"], entries["state_1__A"]["p_value"]
print("RESULT" + json.dumps({"depleted": depleted, "enriched": enriched}))
"""
        env = {
            "PATH": os.environ.get("PATH", ""),
            "HOME": os.environ.get("HOME", ""),
            "TMPDIR": os.environ.get("TMPDIR", tempfile.gettempdir()),
            "PYTHONPATH": str(Path(metainformant.__file__).resolve().parents[1]),
            "PYTHONDONTWRITEBYTECODE": "1",
        }
        completed = subprocess.run(
            [sys.executable, "-c", script], capture_output=True, text=True, check=True, env=env, timeout=120
        )
        line = next(line for line in completed.stdout.splitlines() if line.startswith("RESULT"))
        got = json.loads(line[len("RESULT") :])

        assert got["depleted"] == pytest.approx(_exact_enrichment_p(2, 48, 48, 2), rel=1e-9)
        assert got["depleted"] > 0.99  # it used to be about 2e-6: "enriched"
        assert got["enriched"] == pytest.approx(_exact_enrichment_p(48, 2, 2, 48), rel=1e-9)
        assert got["enriched"] < 1e-20

    def test_no_annotation_means_no_test(self) -> None:
        summary = compute_state_enrichment([0, 0, 1, 1], {}, 2)["summary"]

        assert [(s["most_enriched_annotation"], s["best_fold_enrichment"], s["best_p_value"]) for s in summary] == [
            (None, None, None),
            (None, None, None),
        ]


# ---------------------------------------------------------------------------
# Transition rates
# ---------------------------------------------------------------------------


class TestTransitionRates:
    def test_a_state_that_never_occurs_has_no_transition_statistics(self) -> None:
        result = compute_state_transition_rates([0, 0, 1, 1, 0], n_states=3)

        assert all(math.isnan(p) for p in result["transition_probabilities"][2])
        assert math.isnan(result["self_transition_rates"][2])
        assert math.isnan(result["transition_entropy"][2])
        assert math.isnan(result["mean_segment_length"][2])
        # the observed states are unchanged: 0 -> 0 once, 0 -> 1 once, 1 -> 1 once, 1 -> 0 once
        assert result["transition_probabilities"][0] == [0.5, 0.5, 0.0]
        assert result["self_transition_rates"][:2] == [0.5, 0.5]
        assert result["transition_entropy"][:2] == [1.0, 1.0]
        assert result["mean_segment_length"][:2] == [2.0, 2.0]

    def test_a_state_seen_only_in_the_last_bin_has_nothing_to_report_either(self) -> None:
        result = compute_state_transition_rates([0, 0, 0, 1], n_states=2)

        assert result["transition_counts"] == [[2, 1], [0, 0]]
        assert all(math.isnan(p) for p in result["transition_probabilities"][1])
        assert math.isnan(result["mean_segment_length"][1])
        assert result["n_transitions"] == 3

    def test_a_state_that_only_stays_has_an_infinite_segment_length(self) -> None:
        result = compute_state_transition_rates([0, 0, 0, 0], n_states=1)

        assert result["self_transition_rates"] == [1.0]
        assert result["mean_segment_length"] == [math.inf]

    def test_an_assignment_outside_the_states_is_an_error(self) -> None:
        with pytest.raises(ValueError, match=r"\[0, 2\)"):
            compute_state_transition_rates([0, 1, 2, 1], n_states=2)
        with pytest.raises(ValueError, match=r"\[0, 2\)"):
            compute_state_transition_rates([0, -1, 1], n_states=2)

    def test_the_counts_cover_every_transition(self) -> None:
        assignments = [0, 1, 2, 0, 0, 2, 1, 1, 0]
        result = compute_state_transition_rates(assignments, n_states=3)

        assert sum(sum(row) for row in result["transition_counts"]) == result["n_transitions"] == 8


# ---------------------------------------------------------------------------
# compare_tracks
# ---------------------------------------------------------------------------


def _track(features: List[Tuple[str, int, int, float]]) -> GenomicTrack:
    track = GenomicTrack(name="t")
    for chromosome, start, end, value in features:
        track.add_feature(chromosome, start, end, value=value)
    return track


def _per_base(features: List[Tuple[str, int, int, float]], chromosome: str, length: int) -> np.ndarray:
    """The track's value at each base of a chromosome (highest value where features overlap), NaN where uncovered."""
    values = np.full(length, np.nan)
    for name, start, end, value in features:
        if name == chromosome:
            values[start:end] = np.fmax(values[start:end], value)
    return values


class TestCompareTracks:
    FIRST = [("chr1", 0, 100, 1.0), ("chr1", 100, 200, 2.0), ("chr1", 200, 400, 3.0), ("chr2", 0, 300, 5.0)]
    SECOND = [("chr1", 0, 150, 10.0), ("chr1", 150, 300, 20.0), ("chr1", 300, 400, 40.0), ("chr2", 100, 300, 2.0)]

    def _expected(self, first: List[Tuple[str, int, int, float]], second: List[Tuple[str, int, int, float]]) -> Any:
        x: List[float] = []
        y: List[float] = []
        for chromosome in ("chr1", "chr2"):
            a, b = _per_base(first, chromosome, 400), _per_base(second, chromosome, 400)
            both = ~np.isnan(a) & ~np.isnan(b)
            x.extend(a[both])
            y.extend(b[both])
        return float(np.corrcoef(x, y)[0, 1]), len(x)

    def test_the_correlation_is_taken_base_by_base_over_the_shared_bases(self) -> None:
        want, shared = self._expected(self.FIRST, self.SECOND)

        result = compare_tracks(_track(self.FIRST), _track(self.SECOND))

        assert result["correlation_coefficient"] == pytest.approx(want, abs=1e-12)
        assert result["shared_bases"] == shared == 400 + 200
        assert result["shared_chromosomes"] == 2

    def test_it_does_not_depend_on_the_order_the_features_were_added_in(self) -> None:
        forward = compare_tracks(_track(self.FIRST), _track(self.SECOND))
        scrambled = compare_tracks(_track(self.FIRST[::-1]), _track(self.SECOND[2:] + self.SECOND[:2]))

        assert scrambled["correlation_coefficient"] == pytest.approx(forward["correlation_coefficient"], abs=1e-12)

    def test_where_features_of_one_track_overlap_the_highest_value_counts(self) -> None:
        first = [("chr1", 0, 100, 5.0), ("chr1", 50, 150, 9.0), ("chr1", 150, 400, 1.0)]
        second = [("chr1", 0, 200, 3.0), ("chr1", 200, 400, 8.0)]
        a, b = _per_base(first, "chr1", 400), _per_base(second, "chr1", 400)

        result = compare_tracks(_track(first), _track(second))

        assert result["correlation_coefficient"] == pytest.approx(float(np.corrcoef(a, b)[0, 1]), abs=1e-12)

    def test_a_perfectly_proportional_pair_correlates_at_one_however_it_is_cut(self) -> None:
        first = [("chr1", 0, 50, 1.0), ("chr1", 50, 250, 4.0)]
        second = [("chr1", 0, 125, 2.0), ("chr1", 125, 250, 8.0)]
        # base values: first 1 on [0,50), 4 on [50,250); second 2 on [0,125), 8 on [125,250): not proportional
        a, b = _per_base(first, "chr1", 250), _per_base(second, "chr1", 250)

        result = compare_tracks(_track(first), _track(second))

        assert result["correlation_coefficient"] == pytest.approx(float(np.corrcoef(a, b)[0, 1]), abs=1e-12)
        assert result["correlation_coefficient"] < 1.0

        cut_differently = [("chr1", 0, 125, 2.0), ("chr1", 125, 250, 8.0)]
        proportional = [("chr1", 0, 125, 1.0), ("chr1", 125, 250, 4.0)]
        assert compare_tracks(_track(proportional), _track(cut_differently))[
            "correlation_coefficient"
        ] == pytest.approx(1.0, abs=1e-12)

    def test_no_shared_base_or_a_single_value_leaves_the_correlation_undefined(self) -> None:
        apart = compare_tracks(_track([("chr1", 0, 10, 1.0)]), _track([("chr1", 10, 20, 2.0)]))
        other_chromosome = compare_tracks(_track([("chr1", 0, 10, 1.0)]), _track([("chr2", 0, 10, 2.0)]))
        constant = compare_tracks(_track([("chr1", 0, 10, 1.0), ("chr1", 10, 20, 3.0)]), _track([("chr1", 0, 20, 7.0)]))

        assert apart["correlation_coefficient"] is None and apart["shared_bases"] == 0
        assert other_chromosome["correlation_coefficient"] is None and other_chromosome["shared_chromosomes"] == 0
        assert constant["correlation_coefficient"] is None and constant["shared_bases"] == 20

    def test_chromosome_counts_are_reported(self) -> None:
        result = compare_tracks(_track([("chr1", 0, 10, 1.0), ("chr3", 0, 10, 1.0)]), _track([("chr1", 0, 10, 1.0)]))

        assert (result["shared_chromosomes"], result["unique_to_track1"], result["unique_to_track2"]) == (1, 1, 0)
        assert (result["track1_features"], result["track2_features"]) == (2, 1)


# ---------------------------------------------------------------------------
# identify_tss_enrichment
# ---------------------------------------------------------------------------


class TestTssEnrichment:
    PEAKS = [ATACPeak("chr1", 900, 1100), ATACPeak("chr1", 5000, 5400), ATACPeak("chr2", 200, 300)]
    TSS = {"chr1": [1000, 3000, 5200], "chr2": [250, 9000]}

    def test_without_a_genome_size_there_is_no_chance_expectation(self) -> None:
        result = identify_tss_enrichment(self.PEAKS, self.TSS, window_size=500)

        assert (result["total_tss"], result["enriched_tss"]) == (5, 3)
        assert result["enrichment_ratio"] == pytest.approx(3 / 5)
        assert result["expected_ratio"] is None
        assert result["fold_enrichment"] is None

    def test_the_expectation_is_the_chance_that_a_random_window_meets_a_peak(self) -> None:
        genome_size = 1_000_000
        mean_length = (200 + 400 + 100) / 3

        result = identify_tss_enrichment(self.PEAKS, self.TSS, window_size=500, genome_size=genome_size)

        want = 1.0 - math.exp(-3 * (500 + mean_length) / genome_size)
        assert result["expected_ratio"] == pytest.approx(want, rel=1e-12)
        assert result["fold_enrichment"] == pytest.approx((3 / 5) / want, rel=1e-12)

    def test_a_smaller_genome_means_a_higher_chance_expectation_and_less_enrichment(self) -> None:
        small = identify_tss_enrichment(self.PEAKS, self.TSS, window_size=500, genome_size=100_000)
        large = identify_tss_enrichment(self.PEAKS, self.TSS, window_size=500, genome_size=100_000_000)

        assert small["expected_ratio"] > large["expected_ratio"]
        assert small["fold_enrichment"] < large["fold_enrichment"]

    def test_no_peaks_leave_the_fold_undefined(self) -> None:
        result = identify_tss_enrichment([], self.TSS, window_size=500, genome_size=1_000_000)

        assert result["enriched_tss"] == 0 and result["expected_ratio"] == 0.0
        assert result["fold_enrichment"] is None

    @pytest.mark.parametrize("genome_size", [0, -5])
    def test_the_genome_size_must_be_positive(self, genome_size: int) -> None:
        with pytest.raises(ValueError, match="genome_size must be positive"):
            identify_tss_enrichment(self.PEAKS, self.TSS, genome_size=genome_size)
