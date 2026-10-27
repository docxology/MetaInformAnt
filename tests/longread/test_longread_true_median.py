"""Regression tests: long-read medians must be true medians.

Four places computed ``sorted(values)[len(values) // 2]`` and called the result a median. For an even
number of values that is the *upper* of the two middle values, not their mean:

* ``detect_sv_from_long_reads`` merged clustered calls at the upper-median start, end and size;
* ``generate_sv_summary`` reported the upper median as ``median_size``;
* the SV pipeline's ``sv_summary`` step reported it as ``size_stats[type]["median"]``;
* ``plot_read_length_histogram`` drew its "Median" line, and printed its label, at the upper median.

Breakpoints and sizes of a ``StructuralVariant`` are integer base-pair values, so the half-integer median
of an even-sized cluster is truncated to an integer there (as ``int(numpy.median(...))`` would); the
summaries and the plot label keep the exact value. The oracle is ``numpy.median``. All alignments and
reads are real data run through the real functions; no test doubles.
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import Any

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pytest

from metainformant.longread.analysis.structural import detect_sv_from_long_reads
from metainformant.longread.utils.summary import generate_sv_summary
from metainformant.longread.visualization.plots import plot_read_length_histogram
from metainformant.longread.workflow.orchestrator import LongReadOrchestrator
from metainformant.longread.workflow.pipelines import get_sv_pipeline_config

matplotlib.use("Agg")


def _deletion_read(name: str, start: int, deletion: int, chrom: str = "chr1") -> dict[str, Any]:
    """A primary alignment with one ``deletion``-bp CIGAR deletion starting 500 bp after ``start``."""
    return {
        "read_name": name,
        "reference_name": chrom,
        "reference_start": start,
        "mapping_quality": 60,
        "is_unmapped": False,
        "cigar_string": f"500M{deletion}D500M",
    }


def _clustered_calls(reads: list[tuple[int, int]]) -> list[Any]:
    """Merge one read per ``(start, deletion)`` pair into calls and require them to form a single call."""
    alignments = [_deletion_read(f"read{i}", start, size) for i, (start, size) in enumerate(reads)]
    calls = detect_sv_from_long_reads(alignments, min_size=50, min_support=len(reads))
    assert len(calls) == 1
    assert calls[0].supporting_reads == len(reads)
    return calls


@pytest.mark.parametrize(
    ("reads", "expected"),
    [
        # Deletions start at 1500, 1510 and span 100, 120 bp: starts 1500/1510, ends 1600/1630, sizes 100/120.
        pytest.param([(1000, 100), (1010, 120)], (1505, 1615, 110), id="two-calls"),
        # Starts 1500, 1510, 1520, 1540; ends 1600, 1620, 1640, 1740; sizes 100, 110, 120, 200.
        pytest.param(
            [(1000, 100), (1010, 110), (1020, 120), (1040, 200)],
            (1515, 1630, 115),
            id="four-calls",
        ),
        # Identical starts; ends 1600/1601 and sizes 100/101 have the half-integer median 1600.5 / 100.5.
        pytest.param([(1000, 100), (1000, 101)], (1500, 1600, 100), id="half-integer-medians-truncate"),
        # An odd count already has a middle value, so nothing changes there.
        pytest.param([(1000, 100), (1010, 110), (1030, 160)], (1510, 1620, 110), id="three-calls"),
    ],
)
def test_merged_sv_calls_use_the_true_median_breakpoints_and_size(
    reads: list[tuple[int, int]], expected: tuple[int, int, int]
) -> None:
    (call,) = _clustered_calls(reads)

    starts = [start + 500 for start, _ in reads]
    ends = [start + 500 + size for start, size in reads]
    sizes = [size for _, size in reads]
    assert (call.start, call.end, call.size) == expected
    assert (call.start, call.end, call.size) == tuple(int(np.median(v)) for v in (starts, ends, sizes))
    assert all(isinstance(value, int) for value in (call.start, call.end, call.size))


def test_sv_summary_reports_the_true_median_size() -> None:
    sizes = [100, 200, 300, 400]
    variants = [{"sv_type": "DEL", "size": size} for size in sizes]

    summary = generate_sv_summary(variants)

    assert summary["size_distribution"]["median_size"] == pytest.approx(np.median(sizes))  # 250, not 300
    assert summary["size_distribution"]["mean_size"] == pytest.approx(250.0)


def test_sv_summary_median_of_an_odd_count_is_the_middle_size() -> None:
    summary = generate_sv_summary([{"sv_type": "DEL", "size": size} for size in (100, 200, 700)])

    assert summary["size_distribution"]["median_size"] == 200


def test_sv_pipeline_summary_reports_the_true_median_size_per_type(tmp_path: Path) -> None:
    sizes = [100, 200, 300, 400]
    alignments = []
    for index, size in enumerate(sizes):  # one deletion per locus, each supported by two reads
        for replicate in range(2):
            alignments.append(_deletion_read(f"locus{index}_read{replicate}", 10_000 * (index + 1), size))
    config = get_sv_pipeline_config(min_support=2, sv_types=["DEL"])

    result = LongReadOrchestrator(config, tmp_path / "sv").run_sv_pipeline(alignments)

    assert result.success
    stats = result.summary["sv_summary"]["size_stats"]["DEL"]
    assert stats["count"] == len(sizes)
    assert stats["median"] == pytest.approx(np.median(sizes))  # 250, not 300
    assert stats["mean"] == pytest.approx(np.mean(sizes))


def _median_label(path: Path) -> float:
    """Value printed in the histogram legend (text is kept as text in the SVG, not drawn as paths)."""
    match = re.search(r"Median = ([0-9,.]+) bp", path.read_text(encoding="utf-8"))
    assert match is not None, "the histogram has no median label"
    return float(match.group(1).replace(",", ""))


@pytest.mark.parametrize(
    "lengths",
    [
        pytest.param([1000, 2000, 3000, 4000], id="even-count"),
        pytest.param([1000, 2001, 3000, 4000], id="even-count-half-bp"),
        pytest.param([1000, 2000, 3000], id="odd-count"),
    ],
)
def test_read_length_histogram_labels_the_true_median(tmp_path: Path, lengths: list[int]) -> None:
    with plt.rc_context({"svg.fonttype": "none"}):
        saved = plot_read_length_histogram(lengths, tmp_path / "lengths.svg", log_scale=False, bins=10)

    assert _median_label(saved) == pytest.approx(np.median(lengths))
    plt.close("all")
