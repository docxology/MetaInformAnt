"""Regression tests for numpy and float read lengths in ``quality.metrics``, the length histogram and the QC summary.

``read_length_stats`` and ``calculate_throughput`` accepted ``int`` only, so every element of
``np.array([...])`` (``np.int64``, which is not an ``int``) was skipped and the statistics came back
empty (``count == 0``); the length extraction behind ``plot_read_length_histogram`` had the same check.
The same check skipped float lengths (``1000.0``, a float array, a pandas column with gaps) without a
word, while a dict ``{"length": 1500.7}`` was cut to 1500 without one. Whole-number floats are counted now and
anything else raises. ``generate_qc_summary(np.array(...))`` raised on the truthiness of the array.
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
from metainformant.longread.utils.summary import generate_assembly_summary, generate_qc_summary
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


class TestFloatLengths:
    """A float length is a whole number of bases or an error; it is never skipped or truncated."""

    @pytest.mark.parametrize("dtype", [np.float32, np.float64])
    def test_whole_number_float_arrays_give_the_statistics_of_the_integers(self, dtype: type) -> None:
        stats = read_length_stats(np.array(LENGTHS, dtype=dtype))
        assert stats == read_length_stats(LENGTHS)
        assert stats.count == len(LENGTHS)
        assert type(stats.total_bases) is int and type(stats.n50) is int

    def test_whole_number_python_floats_and_numeric_text_in_a_dict(self) -> None:
        assert read_length_stats([float(n) for n in LENGTHS]) == read_length_stats(LENGTHS)
        assert read_length_stats([{"length": str(n)} for n in LENGTHS]) == read_length_stats(LENGTHS)
        assert read_length_stats([{"length": "1e3"}, {"length": "2000.0"}, {"length": " 3000 "}]).total_bases == 6000

    def test_floats_inside_dicts_and_objects(self) -> None:
        class Alignment:
            query_length = 4000.0

        reads = [{"length": 1000.0}, {"length": np.float32(2000.0)}, Alignment()]
        stats = read_length_stats(reads)
        assert (stats.count, stats.total_bases) == (3, 7000)

    def test_the_mean_of_lengths_is_a_usable_length(self) -> None:
        # a float that came out of arithmetic: 2 * 1500.0 / 2
        assert read_length_stats([np.mean([1000, 2000]) * 2 / 2 - 0.0, 500]).total_bases == 2000

    @pytest.mark.parametrize("bad", [1500.5, 0.1, float("nan"), float("inf"), -float("inf"), np.float64(2.25)])
    def test_a_fractional_or_non_finite_length_raises(self, bad: float) -> None:
        with pytest.raises(ValueError, match="whole number"):
            read_length_stats([1000, bad, 2000])
        with pytest.raises(ValueError, match="whole number"):
            calculate_throughput([1000, bad])
        with pytest.raises(ValueError, match="whole number"):
            plots._extract_lengths([1000, bad])

    def test_a_float_array_with_a_gap_raises_instead_of_dropping_the_gap(self) -> None:
        with pytest.raises(ValueError, match="whole number"):
            read_length_stats(np.array([1000.0, np.nan, 3000.0]))

    def test_a_fractional_dict_length_is_not_truncated(self) -> None:
        # int(1500.7) used to count this read as 1500 bp
        with pytest.raises(ValueError, match="whole number"):
            read_length_stats([{"length": 1500.7}])

    def test_text_that_is_not_a_number_raises(self) -> None:
        with pytest.raises(ValueError):
            read_length_stats([{"length": "abc"}])
        with pytest.raises(ValueError, match="whole number"):
            read_length_stats([{"length": "1500.5"}])

    @pytest.mark.parametrize("bad", [True, False])
    def test_a_bool_is_not_a_length(self, bad: bool) -> None:
        with pytest.raises(TypeError, match="bool"):
            read_length_stats([1000, bad])
        with pytest.raises(TypeError, match="bool"):
            read_length_stats([{"length": bad}])

    def test_a_missing_length_still_raises_a_type_error(self) -> None:
        with pytest.raises(TypeError):
            read_length_stats([{"length": None}])

    def test_throughput_counts_whole_number_floats(self) -> None:
        assert calculate_throughput(np.array(LENGTHS, dtype=float), run_duration=2.0) == calculate_throughput(
            LENGTHS, run_duration=2.0
        )

    def test_reads_without_a_length_are_still_left_out(self) -> None:
        class NoLength:
            sequence = ""

        assert read_length_stats([1000, NoLength(), {"sequence": ""}, {"read_id": "x"}]).count == 1

    def test_histogram_of_a_float_array_equals_the_histogram_of_the_integers(self, tmp_path: Path) -> None:
        assert plots._extract_lengths(np.array([1000.0, 2000.0, 0.0, 3000.0])) == [1000, 2000, 3000]
        from_ints = plots.plot_read_length_histogram(LENGTHS, tmp_path / "ints.png", bins=10)
        from_floats = plots.plot_read_length_histogram(np.array(LENGTHS, dtype=float), tmp_path / "floats.png", bins=10)
        empty = plots.plot_read_length_histogram([], tmp_path / "empty.png", bins=10)
        assert from_floats.read_bytes() == from_ints.read_bytes() != empty.read_bytes()

    def test_batch_metrics_for_float_arrays(self) -> None:
        result = batch_compute_metrics({"a": np.array(LENGTHS, dtype=float)})
        assert result.successful == 1 and result.results["a"].count == len(LENGTHS)


class TestQcSummaryOfArrays:
    """``generate_qc_summary`` tested ``not read_lengths``, which raises for a numpy array of two or more values."""

    QUALITY = [6.5, 9.0, 12.5, 18.0, 21.0, 24.5]

    def test_an_integer_array_gives_the_summary_of_the_list(self) -> None:
        assert generate_qc_summary(np.array(LENGTHS)) == generate_qc_summary(LENGTHS)
        assert generate_qc_summary(np.array(LENGTHS))["total_reads"] == len(LENGTHS)

    def test_a_float_array_gives_the_summary_of_the_list(self) -> None:
        assert generate_qc_summary(np.array(LENGTHS, dtype=float)) == generate_qc_summary(LENGTHS)

    def test_quality_scores_as_an_array(self) -> None:
        from_array = generate_qc_summary(np.array(LENGTHS), quality_scores=np.array(self.QUALITY))
        from_list = generate_qc_summary(LENGTHS, quality_scores=self.QUALITY)
        assert from_array == from_list
        assert from_array["mean_quality"] == round(sum(self.QUALITY) / len(self.QUALITY), 2)
        assert (from_array["q7_fraction"], from_array["q10_fraction"], from_array["q20_fraction"]) == (
            round(5 / 6, 4),
            round(4 / 6, 4),
            round(2 / 6, 4),
        )

    def test_the_summary_serialises_to_json_whatever_the_input_type(self) -> None:
        summary = generate_qc_summary(
            np.array(LENGTHS), quality_scores=np.array(self.QUALITY, dtype=np.float32), filtered_count=4
        )
        decoded = json.loads(json.dumps(summary))
        assert decoded["total_bases"] == sum(LENGTHS)
        assert decoded["reads_before_filter"] == len(LENGTHS) + 4

    @pytest.mark.parametrize("empty", [[], np.array([]), np.array([], dtype=np.int64), ()])
    def test_empty_inputs_report_no_data(self, empty: object) -> None:
        assert generate_qc_summary(empty)["status"] == "no_data"  # type: ignore[arg-type]

    @pytest.mark.parametrize("no_scores", [None, [], np.array([])])
    def test_missing_or_empty_quality_scores_leave_the_quality_fields_out(self, no_scores: object) -> None:
        summary = generate_qc_summary(np.array(LENGTHS), quality_scores=no_scores)  # type: ignore[arg-type]
        assert "mean_quality" not in summary and summary["status"] == "complete"

    def test_a_generator_of_lengths_works(self) -> None:
        assert generate_qc_summary(n for n in LENGTHS) == generate_qc_summary(LENGTHS)

    def test_a_fractional_length_raises(self) -> None:
        with pytest.raises(ValueError, match="whole number"):
            generate_qc_summary(np.array([1000.0, 2000.5]))

    def test_assembly_summary_of_arrays(self) -> None:
        from_array = generate_assembly_summary(np.array([9000, 2500, 400]), num_reads_used=12)
        assert from_array == generate_assembly_summary([9000, 2500, 400], num_reads_used=12)
        assert from_array["total_bases"] == 11900 and from_array["largest_contig"] == 9000
        assert all(type(from_array[key]) is int for key in ("total_contigs", "total_bases", "n50", "n90"))
        assert json.loads(json.dumps(from_array))["smallest_contig"] == 400
        assert generate_assembly_summary(np.array([9000.0, 2500.0])) == generate_assembly_summary([9000, 2500])
        assert generate_assembly_summary(np.array([], dtype=int))["status"] == "no_assembly"
        with pytest.raises(ValueError, match="whole number"):
            generate_assembly_summary(np.array([9000.0, 2500.25]))
