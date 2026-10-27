"""Regression tests for numpy integer read lengths in ``quality.metrics`` and the length histogram.

``read_length_stats`` and ``calculate_throughput`` accepted ``int`` only, so every element of
``np.array([...])`` (``np.int64``, which is not an ``int``) was skipped and the statistics came back
empty (``count == 0``); the length extraction behind ``plot_read_length_histogram`` had the same check.
All functions run for real on real arrays; there are no mocks.
"""

from __future__ import annotations

import dataclasses
import json
from pathlib import Path

import numpy as np
import pytest

from metainformant.longread.quality.metrics import calculate_n50, calculate_throughput, read_length_stats
from metainformant.longread.utils.batch import batch_compute_metrics
from metainformant.longread.visualization import plots

LENGTHS = [500, 1000, 2000, 3000, 5000, 12000]


@pytest.mark.parametrize("dtype", [np.int16, np.int32, np.int64, np.uint16, np.uint32, np.uint64])
def test_read_length_stats_counts_numpy_integer_arrays(dtype: type) -> None:
    stats = read_length_stats(np.array(LENGTHS, dtype=dtype))
    assert stats == read_length_stats(LENGTHS)
    assert stats.count == len(LENGTHS)


def test_the_headline_case_numpy_array_count_is_not_zero() -> None:
    stats = read_length_stats(np.array([1000, 2000, 3000]))
    assert stats.count == 3
    assert stats.total_bases == 6000
    assert stats.n50 == 3000


def test_statistics_are_plain_python_and_serialisable() -> None:
    stats = read_length_stats(np.array(LENGTHS))
    assert stats.count == len(LENGTHS)
    for name in ("count", "total_bases", "min_length", "max_length", "n50", "n90", "l50"):
        assert type(getattr(stats, name)) is int, name
    assert json.dumps(dataclasses.asdict(stats))  # np.int64 values would raise TypeError here


def test_numpy_scalars_mixed_with_python_ints_and_dicts() -> None:
    reads = [np.int64(1000), 2000, np.uint16(3000), {"length": np.int32(4000)}, {"sequence": "ACGT" * 500}]
    stats = read_length_stats(reads)
    assert (stats.count, stats.total_bases) == (5, 12000)


def test_non_positive_numpy_lengths_are_dropped_like_python_ones() -> None:
    assert read_length_stats(np.array([0, -5, 1500])) == read_length_stats([0, -5, 1500])
    assert read_length_stats(np.array([0, -5, 1500])).count == 1


def test_empty_numpy_array_gives_empty_statistics() -> None:
    assert read_length_stats(np.array([], dtype=np.int64)).count == 0


def test_calculate_throughput_counts_numpy_integer_arrays() -> None:
    from_numpy = calculate_throughput(np.array(LENGTHS), run_duration=2.0)
    assert from_numpy == calculate_throughput(LENGTHS, run_duration=2.0)
    assert from_numpy["total_reads"] == float(len(LENGTHS))
    assert from_numpy["total_bases"] == float(sum(LENGTHS))


def test_calculate_n50_of_a_numpy_array_matches_the_list() -> None:
    assert calculate_n50(np.array(LENGTHS)) == calculate_n50(LENGTHS)


def test_batch_metrics_for_numpy_inputs() -> None:
    result = batch_compute_metrics({"a": np.array(LENGTHS), "b": np.array([100, 200])})
    assert result.successful == 2
    assert result.results["a"].count == len(LENGTHS)
    assert result.results["b"].total_bases == 300


class TestLengthHistogram:
    def test_extract_lengths_accepts_numpy_integers(self) -> None:
        assert plots._extract_lengths(np.array([1000, 2000, 0, 3000])) == [1000, 2000, 3000]
        assert plots._extract_lengths([np.int32(7), 8, np.uint8(0)]) == [7, 8]

    def test_histogram_of_a_numpy_array_equals_the_histogram_of_the_list(self, tmp_path: Path) -> None:
        # with the old check the numpy input gave the empty placeholder plot, a different image
        from_list = plots.plot_read_length_histogram(LENGTHS, tmp_path / "list.png", bins=10)
        from_numpy = plots.plot_read_length_histogram(np.array(LENGTHS), tmp_path / "numpy.png", bins=10)
        empty = plots.plot_read_length_histogram([], tmp_path / "empty.png", bins=10)
        assert from_numpy.read_bytes() == from_list.read_bytes()
        assert from_numpy.read_bytes() != empty.read_bytes()
