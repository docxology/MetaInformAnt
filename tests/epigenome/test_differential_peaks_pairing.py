"""``differential_peaks`` pairs each peak with the peak it overlaps most and measures every row.

Defects:

* An A peak was paired with the *first* overlapping B peak in list order, so the result depended on the order of
  ``peaks_b``, a B peak that overlapped more was passed over, and the other B peaks that overlapped the A peak were
  reported as ``unique_B`` although they overlap an A peak.
* Peaks found in one condition only were given ``p_value`` 0.0, a ``log2_fold_change`` of infinity and a direction of
  "up" in their own condition whatever the signal said: a peak whose signal is higher in the *other* condition was
  labelled up in its own. Their means were already computed from the tracks; the statistics now use them.
* The p-value was 1.0 for a region with fewer than two bases, a placeholder where the test could not be run. Without
  SciPy it was also 1.0 for two constant tracks that differ (a standard error under 1e-12 was taken as "no difference"
  although the t statistic is infinite), and with SciPy those tracks drew its "nearly identical data" RuntimeWarning.
* One signal track was read for every chromosome, so peaks on several chromosomes were measured on the same bases.

Oracles: ``numpy`` means, ``numpy.log2`` and ``scipy.stats.ttest_ind(equal_var=False)`` for the statistics, and brute
force over all pairs for the pairing.

REAL IMPLEMENTATION.
"""

from __future__ import annotations

import itertools
import json
import math
import os
import random
import subprocess
import sys
import tempfile
import warnings
from pathlib import Path
from typing import Any, Dict, List, Sequence, Tuple

import numpy as np
import pytest

import metainformant
from metainformant.epigenome.peak_calling.peak_detection import differential_peaks

scipy_stats = pytest.importorskip("scipy.stats")

Record = Dict[str, Any]


def _peak(start: int, end: int, chrom: str = "chr_unknown") -> Dict[str, Any]:
    return {"chrom": chrom, "start": start, "end": end}


def _track(length: int, bumps: Sequence[Tuple[int, int, float]], seed: int = 0, noise: float = 0.2) -> List[float]:
    """Background of about 1 with seeded noise, plus ``(start, end, height)`` bumps."""
    rng = random.Random(seed)
    values = [max(0.0, 1.0 + rng.gauss(0.0, noise)) for _ in range(length)]
    for start, end, height in bumps:
        for position in range(start, end):
            values[position] += height
    return values


def _expected(track_a: Sequence[float], track_b: Sequence[float], start: int, end: int) -> Dict[str, float]:
    """Mean signals, log2 fold change and Welch p-value of the two tracks over ``[start, end)``, by numpy and scipy."""
    values_a = np.asarray(track_a[start:end], dtype=float)
    values_b = np.asarray(track_b[start:end], dtype=float)
    mean_a, mean_b = float(values_a.mean()), float(values_b.mean())
    with warnings.catch_warnings():
        # SciPy warns about "nearly identical data" when one of the samples is constant, which is a legitimate input
        warnings.simplefilter("ignore", RuntimeWarning)
        p_value = float(scipy_stats.ttest_ind(values_a, values_b, equal_var=False).pvalue)
    return {
        "mean_a": mean_a,
        "mean_b": mean_b,
        "log2fc": float(np.log2((mean_a + 1e-10) / (mean_b + 1e-10))),
        "p_value": p_value,
    }


def _signature(records: Sequence[Record]) -> List[Tuple[Any, ...]]:
    """The records as sorted tuples, so two results can be compared without regard to their order."""
    return sorted(
        (r["status"], r["chrom"], r["start"], r["end"], r["log2_fold_change"], r["direction"]) for r in records
    )


class TestPairingByMaximumOverlap:
    def setup_method(self) -> None:
        # A covers [100, 300). B1 shares only 10 bases with it and comes first in the list; B2 shares 190
        self.track_a = _track(600, [(100, 300, 10.0)], seed=1)
        self.track_b = _track(600, [(100, 290, 4.0), (290, 400, 8.0)], seed=2)
        self.peak_a = _peak(100, 300)
        self.b_small = _peak(290, 400)
        self.b_large = _peak(100, 290)

    def test_a_peak_is_compared_with_the_peak_it_shares_the_most_bases_with(self) -> None:
        records = differential_peaks([self.peak_a], [self.b_small, self.b_large], self.track_a, self.track_b)

        first = records[0]
        assert first["status"] == "shared"
        assert (first["start"], first["end"]) == (100, 300)  # the union with B2, not the union with B1 (100-400)
        want = _expected(self.track_a, self.track_b, 100, 300)
        assert first["mean_signal_a"] == pytest.approx(want["mean_a"], abs=5e-5)
        assert first["mean_signal_b"] == pytest.approx(want["mean_b"], abs=5e-5)
        assert first["log2_fold_change"] == pytest.approx(want["log2fc"], abs=5e-5)
        assert first["p_value"] == pytest.approx(want["p_value"], rel=1e-9)

    def test_a_peak_that_overlaps_something_is_never_called_unique(self) -> None:
        records = differential_peaks([self.peak_a], [self.b_small, self.b_large], self.track_a, self.track_b)

        assert [r["status"] for r in records] == ["shared", "shared"]
        assert {(r["start"], r["end"]) for r in records} == {(100, 300), (100, 400)}

    def test_the_result_does_not_depend_on_the_order_of_the_peaks(self) -> None:
        peaks_a = [_peak(100, 300), _peak(500, 560)]
        peaks_b = [self.b_small, self.b_large, _peak(520, 540), _peak(520, 600)]
        want = _signature(differential_peaks(peaks_a, peaks_b, self.track_a, self.track_b))

        for order_a in itertools.permutations(peaks_a):
            for order_b in itertools.permutations(peaks_b):
                got = _signature(differential_peaks(list(order_a), list(order_b), self.track_a, self.track_b))
                assert got == want

    def test_equal_overlaps_go_to_the_peak_that_starts_first(self) -> None:
        track_a = _track(400, [(100, 200, 5.0)], seed=3)
        track_b = _track(400, [], seed=4)
        peak_a = _peak(100, 200)
        left, right = _peak(50, 150), _peak(150, 250)  # each shares 50 bases with A

        for peaks_b in ([left, right], [right, left]):
            records = differential_peaks([peak_a], peaks_b, track_a, track_b)

            assert [r["status"] for r in records] == ["shared", "shared"]
            # A's own record is the union with the earlier peak, whatever the list order; the other peak's record
            # follows, paired with the A peak it overlaps most
            assert [(r["start"], r["end"]) for r in records] == [(50, 200), (100, 250)]

    def test_every_peak_appears_and_unique_means_no_overlap_at_all(self) -> None:
        rng = random.Random(11)
        track = _track(2000, [], seed=5)

        def random_peaks(count: int) -> List[Dict[str, Any]]:
            peaks = []
            for _ in range(count):
                start = rng.randint(0, 1800)
                peaks.append(_peak(start, start + rng.randint(1, 200)))
            return peaks

        for _ in range(30):
            peaks_a, peaks_b = random_peaks(rng.randint(0, 12)), random_peaks(rng.randint(0, 12))

            records = differential_peaks(peaks_a, peaks_b, track, track)

            def overlaps(first: Dict[str, Any], second: Dict[str, Any]) -> bool:
                return first["start"] < second["end"] and second["start"] < first["end"]

            unique_a = [r for r in records if r["status"] == "unique_A"]
            unique_b = [r for r in records if r["status"] == "unique_B"]
            assert len(unique_a) == sum(not any(overlaps(a, b) for b in peaks_b) for a in peaks_a)
            assert len(unique_b) == sum(not any(overlaps(a, b) for a in peaks_a) for b in peaks_b)
            # Every peak sits in at least one record: the A peaks in their own, the B peaks in the pair they chose
            assert len(records) >= max(len(peaks_a), len(peaks_b))

    def test_swapping_the_conditions_mirrors_the_result(self) -> None:
        peaks_a = [_peak(100, 300), _peak(500, 560), _peak(40, 60)]
        peaks_b = [self.b_small, self.b_large, _peak(520, 540), _peak(580, 599)]

        forward = differential_peaks(peaks_a, peaks_b, self.track_a, self.track_b)
        backward = differential_peaks(peaks_b, peaks_a, self.track_b, self.track_a)

        swap = {"unique_A": "unique_B", "unique_B": "unique_A", "shared": "shared"}
        mirrored = sorted((swap[r["status"]], r["start"], r["end"], round(-r["log2_fold_change"], 4)) for r in backward)
        original = sorted((r["status"], r["start"], r["end"], r["log2_fold_change"]) for r in forward)
        assert mirrored == original

    def test_the_records_of_the_a_peaks_come_first_in_list_order(self) -> None:
        peaks_a = [_peak(500, 560), _peak(100, 300), _peak(10, 20)]

        records = differential_peaks(peaks_a, [self.b_small, self.b_large], self.track_a, self.track_b)

        assert [(r["start"], r["end"]) for r in records[:3]] == [(500, 560), (100, 300), (10, 20)]
        assert [r["status"] for r in records[:3]] == ["unique_A", "shared", "unique_A"]


def _brute_best(queries: Sequence[Tuple[str, int, int]], targets: Sequence[Tuple[str, int, int]]) -> List[Any]:
    """For each query the target with the most shared bases (ties: lowest start, then end, then position)."""
    best: List[Any] = []
    for chrom, start, end in queries:
        choice = None
        key = None
        for position, (t_chrom, t_start, t_end) in enumerate(targets):
            shared = min(end, t_end) - max(start, t_start)
            if t_chrom != chrom or shared <= 0:
                continue
            candidate = (-shared, t_start, t_end, position)
            if key is None or candidate < key:
                key, choice = candidate, position
        best.append(choice)
    return best


class TestBestOverlapsAgainstBruteForce:
    @pytest.mark.parametrize("seed", range(25))
    def test_the_indexed_search_picks_what_a_scan_of_all_pairs_picks(self, seed: int) -> None:
        from metainformant.epigenome.peak_calling.peak_detection import _best_overlaps

        rng = random.Random(seed)

        def intervals(count: int) -> List[Tuple[str, int, int]]:
            rows = []
            for _ in range(count):
                start = rng.randint(0, 400)
                # long peaks nested around short ones exercise the running maximum of the ends
                rows.append((rng.choice(["chr1", "chr2", "chr3"]), start, start + rng.choice([1, 3, 10, 40, 300])))
            return rows

        queries, targets = intervals(rng.randint(0, 60)), intervals(rng.randint(0, 60))

        assert _best_overlaps(queries, targets) == _brute_best(queries, targets)


class TestUniquePeaksAreMeasured:
    def test_a_unique_peak_whose_other_condition_signal_is_higher_is_not_up_in_its_own(self) -> None:
        # A found a peak at [200, 260); B did not call one there, but its track is far higher
        track_a = _track(600, [(200, 260, 1.0)], seed=6)
        track_b = _track(600, [(200, 260, 30.0)], seed=7)

        (record,) = differential_peaks([_peak(200, 260)], [], track_a, track_b)

        want = _expected(track_a, track_b, 200, 260)
        assert record["status"] == "unique_A"
        assert record["direction"] == "up_in_B"
        assert math.isfinite(record["log2_fold_change"])
        assert record["log2_fold_change"] == pytest.approx(want["log2fc"], abs=5e-5) and want["log2fc"] < -3
        assert 0.0 < record["p_value"] < 1e-10
        assert record["p_value"] == pytest.approx(want["p_value"], rel=1e-9)

    def test_a_unique_peak_with_the_same_signal_in_both_conditions_is_unchanged(self) -> None:
        track_a = _track(600, [(200, 260, 5.0)], seed=8)
        track_b = _track(600, [(200, 260, 5.0)], seed=9)

        (record,) = differential_peaks([], [_peak(200, 260)], track_a, track_b)

        assert record["status"] == "unique_B"
        assert record["direction"] == "unchanged"
        assert abs(record["log2_fold_change"]) < 0.2 and math.isfinite(record["log2_fold_change"])
        assert record["p_value"] == pytest.approx(_expected(track_a, track_b, 200, 260)["p_value"], rel=1e-9)
        assert record["p_value"] > 0.05

    def test_a_unique_peak_with_no_signal_in_the_other_condition_has_a_finite_fold_change(self) -> None:
        track_a = _track(600, [(200, 260, 20.0)], seed=10)
        track_b = [0.0] * 600

        (record,) = differential_peaks([_peak(200, 260)], [], track_a, track_b)

        assert record["direction"] == "up_in_A"
        assert math.isfinite(record["log2_fold_change"]) and record["log2_fold_change"] > 20
        assert record["mean_signal_b"] == 0.0


class TestDegenerateTests:
    def test_constant_tracks_that_differ_are_perfectly_separated_not_p_equals_one(self) -> None:
        track_a, track_b = [5.0] * 400, [1.0] * 400

        (record,) = differential_peaks([_peak(100, 160)], [_peak(100, 160)], track_a, track_b)

        assert record["status"] == "shared"
        assert record["p_value"] == 0.0
        assert record["direction"] == "up_in_A"
        assert record["log2_fold_change"] == pytest.approx(math.log2(5.0), abs=5e-5)

    def test_identical_constant_tracks_have_a_p_value_of_one(self) -> None:
        track = [3.0] * 400

        (record,) = differential_peaks([_peak(100, 160)], [_peak(100, 160)], track, list(track))

        assert record["p_value"] == 1.0
        assert record["direction"] == "unchanged"
        assert record["log2_fold_change"] == 0.0

    def test_a_one_base_region_cannot_be_tested(self) -> None:
        track_a, track_b = [9.0] * 100, [1.0] * 100

        (record,) = differential_peaks([_peak(40, 41)], [_peak(40, 41)], track_a, track_b)

        assert math.isnan(record["p_value"])
        assert record["direction"] == "unchanged"  # an untested difference is not called
        assert record["log2_fold_change"] == pytest.approx(math.log2(9.0), abs=5e-5)

    def test_one_constant_track_against_a_noisy_one_is_still_tested(self) -> None:
        track_a = [5.0] * 400
        track_b = _track(400, [], seed=12, noise=0.5)

        (record,) = differential_peaks([_peak(100, 200)], [_peak(100, 200)], track_a, track_b)

        want = _expected(track_a, track_b, 100, 200)
        assert record["p_value"] == pytest.approx(want["p_value"], rel=1e-9)
        assert record["direction"] == "up_in_A"

    @pytest.mark.parametrize("seed", range(10))
    def test_the_p_value_is_welchs_t_test(self, seed: int) -> None:
        rng = random.Random(seed)
        track_a = [max(0.0, rng.gauss(3.0, 1.0 + seed / 5)) for _ in range(500)]
        track_b = [max(0.0, rng.gauss(3.0 + seed / 10, 0.5 + seed / 3)) for _ in range(500)]

        (record,) = differential_peaks([_peak(50, 400)], [_peak(60, 410)], track_a, track_b)

        want = _expected(track_a, track_b, 50, 410)
        assert record["p_value"] == pytest.approx(want["p_value"], rel=1e-9, abs=1e-300)


_WITHOUT_SCIPY = r"""
import json, sys
sys.modules["scipy"] = None  # makes ``import scipy`` fail, as if it were not installed
from metainformant.epigenome.peak_calling import peak_detection

assert peak_detection.HAS_SCIPY is False
peak = {"chrom": "c", "start": 100, "end": 200}
cases = {
    "separated": ([5.0] * 400, [1.0] * 400),
    "identical": ([3.0] * 400, [3.0] * 400),
    "noisy": ([float(x) for x in json.loads(sys.argv[1])], [float(x) for x in json.loads(sys.argv[2])]),
}
out = {}
for name, (a, b) in cases.items():
    (record,) = peak_detection.differential_peaks([peak], [peak], a, b)
    out[name] = [record["p_value"], record["direction"]]
print("RESULT" + json.dumps(out))
"""


class TestWithoutScipy:
    def test_the_fallback_does_not_call_separated_constant_tracks_unchanged(self) -> None:
        rng = random.Random(4)
        noisy_a = [max(0.0, rng.gauss(3.0, 1.0)) for _ in range(400)]
        noisy_b = [max(0.0, rng.gauss(3.6, 1.0)) for _ in range(400)]
        env = {
            "PATH": os.environ.get("PATH", ""),
            "HOME": os.environ.get("HOME", ""),
            "TMPDIR": os.environ.get("TMPDIR", tempfile.gettempdir()),
            "PYTHONPATH": str(Path(metainformant.__file__).resolve().parents[1]),  # the package under test
            "PYTHONDONTWRITEBYTECODE": "1",
        }

        completed = subprocess.run(
            [sys.executable, "-c", _WITHOUT_SCIPY, json.dumps(noisy_a), json.dumps(noisy_b)],
            capture_output=True,
            text=True,
            check=True,
            env=env,
            timeout=120,
        )
        line = next(line for line in completed.stdout.splitlines() if line.startswith("RESULT"))
        got = json.loads(line[len("RESULT") :])

        assert got["separated"] == [0.0, "up_in_A"]  # it was [1.0, "unchanged"]
        assert got["identical"] == [1.0, "unchanged"]
        # The normal approximation of Welch's t-test over the peak [100, 200)
        a, b = np.asarray(noisy_a[100:200]), np.asarray(noisy_b[100:200])
        t = abs(a.mean() - b.mean()) / math.sqrt(a.var(ddof=1) / a.size + b.var(ddof=1) / b.size)
        assert got["noisy"][0] == pytest.approx(math.erfc(t / math.sqrt(2.0)), rel=1e-9)


class TestScaleFreeP:
    def test_welchs_test_does_not_depend_on_the_unit_of_the_signal(self) -> None:
        # Variances of 1e-340 underflow to zero: the p-value must not turn into a placeholder or an error
        rng = random.Random(21)
        base_a = [max(0.0, rng.gauss(3.0, 1.0)) for _ in range(300)]
        base_b = [max(0.0, rng.gauss(4.0, 1.0)) for _ in range(300)]
        scale = 1e-170

        (unit,) = differential_peaks([_peak(10, 200)], [_peak(10, 200)], base_a, base_b)
        (tiny,) = differential_peaks(
            [_peak(10, 200)], [_peak(10, 200)], [v * scale for v in base_a], [v * scale for v in base_b]
        )

        assert 0.0 < unit["p_value"] < 1.0
        assert tiny["p_value"] == pytest.approx(unit["p_value"], rel=1e-6)


class TestSeveralChromosomes:
    def test_each_chromosome_is_measured_on_its_own_track(self) -> None:
        tracks_a = {"chr1": _track(300, [(100, 160, 20.0)], seed=1), "chr2": _track(300, [], seed=2)}
        tracks_b = {"chr1": _track(300, [], seed=3), "chr2": _track(300, [(100, 160, 20.0)], seed=4)}
        peaks = [_peak(100, 160, "chr1"), _peak(100, 160, "chr2")]

        records = differential_peaks(peaks, peaks, tracks_a, tracks_b)

        by_chrom = {r["chrom"]: r for r in records}
        assert by_chrom["chr1"]["direction"] == "up_in_A"
        assert by_chrom["chr2"]["direction"] == "up_in_B"
        assert by_chrom["chr1"]["mean_signal_a"] == pytest.approx(float(np.mean(tracks_a["chr1"][100:160])), abs=5e-5)
        assert by_chrom["chr2"]["mean_signal_b"] == pytest.approx(float(np.mean(tracks_b["chr2"][100:160])), abs=5e-5)

    def test_peaks_do_not_pair_across_chromosomes(self) -> None:
        track = _track(300, [], seed=5)

        records = differential_peaks(
            [_peak(100, 160, "chr1")],
            [_peak(100, 160, "chr2")],
            {"chr1": track, "chr2": track},
            {"chr1": track, "chr2": track},
        )

        assert sorted((r["status"], r["chrom"]) for r in records) == [("unique_A", "chr1"), ("unique_B", "chr2")]

    def test_one_plain_track_cannot_describe_two_chromosomes(self) -> None:
        track = _track(300, [], seed=6)

        with pytest.raises(ValueError, match="several chromosomes"):
            differential_peaks([_peak(10, 20, "chr1")], [_peak(10, 20, "chr2")], track, track)

    def test_a_chromosome_with_peaks_needs_a_track(self) -> None:
        track = _track(300, [], seed=7)

        with pytest.raises(ValueError, match="No signal track"):
            differential_peaks([_peak(10, 20, "chr3")], [], {"chr1": track}, {"chr1": track})

    def test_track_dicts_must_agree(self) -> None:
        track = _track(300, [], seed=8)

        with pytest.raises(ValueError, match="must match"):
            differential_peaks([], [], {"chr1": track}, {"chr1": track[:-1]})
        with pytest.raises(ValueError, match="must match"):
            differential_peaks([], [], {"chr1": track}, {"chr2": track})
        with pytest.raises(ValueError, match="both be lists or both be dicts"):
            differential_peaks([], [], {"chr1": track}, track)


class TestInputChecks:
    def test_signal_lengths_must_match_even_without_peaks(self) -> None:
        with pytest.raises(ValueError, match="must match"):
            differential_peaks([], [], [1.0, 2.0], [1.0])

    @pytest.mark.parametrize(
        "peak",
        [{"chrom": "c", "start": 10}, {"chrom": "c", "start": 10, "end": 10}, {"chrom": "c", "start": 20, "end": 10}],
    )
    def test_a_peak_without_extent_is_an_error_not_a_row_of_zeros(self, peak: Dict[str, Any]) -> None:
        track = [1.0] * 100

        with pytest.raises(ValueError, match="peak|Peak"):
            differential_peaks([peak], [], track, track)

    def test_a_negative_start_is_an_error(self) -> None:
        track = [1.0] * 100

        with pytest.raises(ValueError, match="0 <= start < end"):
            differential_peaks([{"chrom": "c", "start": -5, "end": 10}], [], track, track)

    def test_a_peak_wholly_beyond_the_track_is_an_error_not_a_zero_mean(self) -> None:
        track = [1.0] * 100

        with pytest.raises(ValueError, match="outside the signal tracks"):
            differential_peaks([_peak(150, 200)], [], track, track)

    def test_a_peak_running_past_the_end_of_the_track_is_measured_on_what_exists(self) -> None:
        track_a = [2.0] * 100
        track_b = [4.0] * 100

        (record,) = differential_peaks([_peak(90, 140)], [], track_a, track_b)

        assert (record["start"], record["end"]) == (90, 140)
        assert record["mean_signal_a"] == 2.0 and record["mean_signal_b"] == 4.0
        assert record["log2_fold_change"] == pytest.approx(-1.0)

    def test_negative_mean_signal_has_no_fold_change(self) -> None:
        with pytest.raises(ValueError, match="negative"):
            differential_peaks([_peak(10, 20)], [], [-1.0] * 50, [1.0] * 50)

    def test_no_peaks_gives_no_records(self) -> None:
        assert differential_peaks([], [], [1.0] * 10, [1.0] * 10) == []
