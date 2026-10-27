"""Regression tests: consensus breakpoints of merged SV calls must be true medians.

``merge_sv_callsets`` documents its consensus start and end as the median breakpoints of the merged calls,
but ``_build_consensus`` took ``sorted(values)[n // 2]``. For an even number of calls that is the upper of
the two middle values, not their mean, so two overlapping calls collapsed onto the later caller's
breakpoints instead of landing between them. Breakpoints are integer coordinates, so the half-integer
median of an even-sized cluster is truncated, as in ``structural_variants.filtering.merge`` (and as
``int(numpy.median(...))`` would).

The oracle is ``numpy.median``. All calls are real data run through the real merge; no test doubles.
"""

from __future__ import annotations

import random
from typing import Any

import numpy as np
import pytest

from metainformant.structural_variants.population import merge_sv_callsets


def _call(start: int, end: int, chrom: str = "chr1", sample: str | None = None) -> dict[str, Any]:
    call: dict[str, Any] = {"chrom": chrom, "start": start, "end": end, "sv_type": "DEL"}
    if sample is not None:
        call["sample"] = sample
    return call


def _merge_one_cluster(spans: list[tuple[int, int]]) -> dict[str, Any]:
    """Merge one call per caller and require that they form a single consensus call."""
    callsets = [[_call(start, end)] for start, end in spans]
    merged = merge_sv_callsets(callsets, min_overlap=0.5, max_breakpoint_distance=500)
    assert len(merged) == 1
    assert merged[0]["support_count"] == len(spans)
    return merged[0]


@pytest.mark.parametrize(
    ("spans", "expected"),
    [
        pytest.param([(10000, 15000), (10100, 15200)], (10050, 15100), id="two-calls"),
        # Starts 10000, 10020, 10040, 10100 and ends 15000, 15010, 15030, 15200.
        pytest.param(
            [(10000, 15000), (10100, 15200), (10020, 15010), (10040, 15030)],
            (10030, 15020),
            id="four-calls-in-any-order",
        ),
        # Starts 10000/10001 and ends 15000/15001 have the half-integer medians 10000.5 and 15000.5.
        pytest.param([(10000, 15000), (10001, 15001)], (10000, 15000), id="half-integer-medians-truncate"),
        # An odd count already has a middle value, so nothing changes there.
        pytest.param([(10000, 14900), (10050, 14980), (10100, 15000)], (10050, 14980), id="three-calls"),
    ],
)
def test_consensus_breakpoints_are_the_true_median(spans: list[tuple[int, int]], expected: tuple[int, int]) -> None:
    consensus = _merge_one_cluster(spans)

    assert (consensus["start"], consensus["end"]) == expected
    assert (consensus["start"], consensus["end"]) == (
        int(np.median([start for start, _ in spans])),
        int(np.median([end for _, end in spans])),
    )
    assert isinstance(consensus["start"], int) and isinstance(consensus["end"], int)


def test_consensus_matches_numpy_on_random_clusters() -> None:
    rng = random.Random(20261026)
    callsets: list[list[dict[str, Any]]] = []
    expected: dict[str, tuple[int, int, int]] = {}
    even_clusters = 0
    for index in range(30):
        chrom = f"chr{index}"  # one chromosome per cluster keeps the clusters apart
        size = rng.randint(2, 9)
        starts = [100_000 + rng.randint(-40, 40) for _ in range(size)]
        ends = [105_000 + rng.randint(-40, 40) for _ in range(size)]
        callsets.extend(
            [_call(start, end, chrom, sample=f"s{index}_{i}")] for i, (start, end) in enumerate(zip(starts, ends))
        )
        expected[chrom] = (int(np.median(starts)), int(np.median(ends)), size)
        even_clusters += size % 2 == 0

    merged = {
        entry["chrom"]: entry for entry in merge_sv_callsets(callsets, min_overlap=0.5, max_breakpoint_distance=500)
    }

    assert set(merged) == set(expected)
    assert even_clusters >= 8
    for chrom, (start, end, size) in expected.items():
        assert (merged[chrom]["start"], merged[chrom]["end"], merged[chrom]["support_count"]) == (
            start,
            end,
            size,
        ), chrom
