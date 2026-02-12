"""Tests for longread batch analysis operations not covered by test_longread.py.

Tests BatchResult, process_batch, batch_filter_reads, and batch_compute_metrics.
All tests use real implementations -- NO MOCKING.
"""

from __future__ import annotations

import random
from typing import Any

import pytest

from metainformant.longread.utils.batch import (
    BatchResult,
    batch_compute_metrics,
    batch_filter_reads,
    process_batch,
)


def _random_dna(length: int, seed: int = 42) -> str:
    """Generate a deterministic random DNA sequence."""
    rng = random.Random(seed)
    return "".join(rng.choices("ACGT", k=length))


def _phred_string(scores: list[int]) -> str:
    """Convert a list of Phred scores to an ASCII quality string (Phred+33)."""
    return "".join(chr(s + 33) for s in scores)


# ---------------------------------------------------------------------------
# BatchResult dataclass
# ---------------------------------------------------------------------------


class TestBatchResult:
    """Tests for BatchResult dataclass."""

    def test_batch_result_default_values(self) -> None:
        """BatchResult has correct defaults."""
        result = BatchResult()
        assert result.total_items == 0
        assert result.successful == 0
        assert result.failed == 0
        assert result.results == {}
        assert result.errors == {}
        assert result.duration_seconds == 0.0
        assert result.items_per_second == 0.0

    def test_batch_result_success_rate_zero_items(self) -> None:
        """success_rate returns 0.0 when total_items is 0."""
        result = BatchResult(total_items=0, successful=0)
        assert result.success_rate == 0.0

    def test_batch_result_success_rate_all_successful(self) -> None:
        """success_rate returns 1.0 when all items succeed."""
        result = BatchResult(total_items=5, successful=5)
        assert result.success_rate == pytest.approx(1.0)

    def test_batch_result_success_rate_partial(self) -> None:
        """success_rate returns correct fraction for partial success."""
        result = BatchResult(total_items=10, successful=7, failed=3)
        assert result.success_rate == pytest.approx(0.7)

    def test_batch_result_custom_results(self) -> None:
        """BatchResult stores per-item results and errors."""
        result = BatchResult(
            total_items=3,
            successful=2,
            failed=1,
            results={"sample_a": {"n50": 5000}, "sample_b": {"n50": 8000}},
            errors={"sample_c": "Processing failed"},
            duration_seconds=10.5,
            items_per_second=0.285,
        )
        assert result.results["sample_a"]["n50"] == 5000
        assert result.errors["sample_c"] == "Processing failed"
        assert result.duration_seconds == 10.5


# ---------------------------------------------------------------------------
# process_batch
# ---------------------------------------------------------------------------


class TestProcessBatch:
    """Tests for process_batch function."""

    def test_process_batch_sequential_all_success(self) -> None:
        """process_batch handles all items successfully in sequential mode."""
        items = {
            "s1": [1000, 2000, 3000],
            "s2": [5000, 8000],
            "s3": [500],
        }

        def processor(item_id: str, data: Any) -> dict[str, Any]:
            return {"count": len(data), "total": sum(data)}

        result = process_batch(items, processor, max_workers=None)
        assert result.total_items == 3
        assert result.successful == 3
        assert result.failed == 0
        assert result.success_rate == pytest.approx(1.0)
        assert "s1" in result.results
        assert result.results["s1"]["count"] == 3
        assert result.results["s2"]["total"] == 13000
        assert result.duration_seconds > 0

    def test_process_batch_with_failures_continue(self) -> None:
        """process_batch handles failures gracefully with continue_on_error=True."""
        items = {
            "good1": 100,
            "bad": -1,
            "good2": 200,
        }

        def processor(item_id: str, data: Any) -> int:
            if data < 0:
                raise ValueError(f"Negative value: {data}")
            return data * 2

        result = process_batch(items, processor, continue_on_error=True)
        assert result.total_items == 3
        assert result.successful == 2
        assert result.failed == 1
        assert "bad" in result.errors
        assert "Negative value" in result.errors["bad"]

    def test_process_batch_with_failures_stop(self) -> None:
        """process_batch stops on error when continue_on_error=False."""
        # We use an ordered dict to control processing order
        items = {
            "a": 10,
            "b": -1,
            "c": 20,
        }

        def processor(item_id: str, data: Any) -> int:
            if data < 0:
                raise ValueError("bad")
            return data

        result = process_batch(items, processor, continue_on_error=False)
        # Should stop after the failure, so not all items are processed
        assert result.failed >= 1
        assert result.successful + result.failed <= result.total_items

    def test_process_batch_empty_items(self) -> None:
        """process_batch with no items returns empty result."""

        def processor(item_id: str, data: Any) -> int:
            return 0

        result = process_batch({}, processor)
        assert result.total_items == 0
        assert result.successful == 0
        assert result.failed == 0

    def test_process_batch_parallel(self) -> None:
        """process_batch works with parallel workers."""
        items = {f"s{i}": i * 100 for i in range(5)}

        def processor(item_id: str, data: Any) -> int:
            return data * 2

        result = process_batch(items, processor, max_workers=2)
        assert result.total_items == 5
        assert result.successful == 5
        assert result.failed == 0
        for i in range(5):
            assert result.results[f"s{i}"] == i * 200


# ---------------------------------------------------------------------------
# batch_filter_reads
# ---------------------------------------------------------------------------


class TestBatchFilterReads:
    """Tests for batch_filter_reads function."""

    def test_batch_filter_reads_basic(self) -> None:
        """batch_filter_reads filters across multiple samples."""
        read_sets = {
            "sample_a": [
                {"read_id": "r1", "sequence": _random_dna(5000, seed=1)},
                {"read_id": "r2", "sequence": _random_dna(500, seed=2)},
                {"read_id": "r3", "sequence": _random_dna(3000, seed=3)},
            ],
            "sample_b": [
                {"read_id": "r4", "sequence": _random_dna(800, seed=4)},
                {"read_id": "r5", "sequence": _random_dna(10000, seed=5)},
            ],
        }

        result = batch_filter_reads(read_sets, min_length=1000, min_quality=0)
        assert result.total_items == 2
        assert result.successful == 2
        assert result.failed == 0

        # sample_a: r1 (5000) and r3 (3000) pass, r2 (500) fails length filter
        sample_a_result = result.results["sample_a"]
        assert sample_a_result["input_count"] == 3
        assert sample_a_result["output_count"] == 2
        assert sample_a_result["pass_rate"] == pytest.approx(2 / 3)

        # sample_b: r5 (10000) passes, r4 (800) fails length filter
        sample_b_result = result.results["sample_b"]
        assert sample_b_result["input_count"] == 2
        assert sample_b_result["output_count"] == 1
        assert sample_b_result["pass_rate"] == pytest.approx(0.5)

    def test_batch_filter_reads_with_quality(self) -> None:
        """batch_filter_reads filters by both length and quality."""
        high_q = _phred_string([20] * 2000)
        low_q = _phred_string([3] * 2000)

        read_sets = {
            "sample_1": [
                {"read_id": "good", "sequence": _random_dna(2000, seed=10), "quality_string": high_q},
                {"read_id": "bad_q", "sequence": _random_dna(2000, seed=11), "quality_string": low_q},
            ],
        }

        result = batch_filter_reads(read_sets, min_length=500, min_quality=10.0)
        assert result.successful == 1
        sample_result = result.results["sample_1"]
        assert sample_result["output_count"] == 1  # Only the high-quality read passes

    def test_batch_filter_reads_empty_sample(self) -> None:
        """batch_filter_reads handles empty sample read lists."""
        read_sets: dict[str, list[dict[str, Any]]] = {
            "empty_sample": [],
        }
        result = batch_filter_reads(read_sets, min_length=1000)
        assert result.successful == 1
        assert result.results["empty_sample"]["input_count"] == 0
        assert result.results["empty_sample"]["output_count"] == 0

    def test_batch_filter_reads_all_pass(self) -> None:
        """batch_filter_reads when all reads pass filters."""
        read_sets = {
            "good_sample": [{"read_id": f"r{i}", "sequence": _random_dna(5000 + i * 1000, seed=i)} for i in range(3)],
        }
        result = batch_filter_reads(read_sets, min_length=100, min_quality=0)
        assert result.results["good_sample"]["output_count"] == 3
        assert result.results["good_sample"]["pass_rate"] == pytest.approx(1.0)


# ---------------------------------------------------------------------------
# batch_compute_metrics
# ---------------------------------------------------------------------------


class TestBatchComputeMetrics:
    """Tests for batch_compute_metrics function."""

    def test_batch_compute_metrics_basic(self) -> None:
        """batch_compute_metrics computes stats for multiple samples."""
        read_sets = {
            "sample_a": [1000, 2000, 3000, 5000, 8000],
            "sample_b": [500, 1500, 2500],
        }

        result = batch_compute_metrics(read_sets)
        assert result.total_items == 2
        assert result.successful == 2
        assert result.failed == 0

        stats_a = result.results["sample_a"]
        assert stats_a.count == 5
        assert stats_a.total_bases == 19000
        assert stats_a.n50 > 0

        stats_b = result.results["sample_b"]
        assert stats_b.count == 3
        assert stats_b.total_bases == 4500

    def test_batch_compute_metrics_single_sample(self) -> None:
        """batch_compute_metrics works with a single sample."""
        read_sets = {
            "only_sample": [10000, 20000, 15000],
        }
        result = batch_compute_metrics(read_sets)
        assert result.successful == 1
        stats = result.results["only_sample"]
        assert stats.count == 3
        assert stats.total_bases == 45000
        assert stats.mean_length == pytest.approx(15000.0)

    def test_batch_compute_metrics_empty_sample(self) -> None:
        """batch_compute_metrics handles samples with no reads."""
        read_sets: dict[str, list[int]] = {
            "empty": [],
        }
        result = batch_compute_metrics(read_sets)
        assert result.successful == 1
        stats = result.results["empty"]
        assert stats.count == 0
        assert stats.total_bases == 0

    def test_batch_compute_metrics_parallel(self) -> None:
        """batch_compute_metrics works with parallel workers."""
        read_sets = {f"sample_{i}": [j * 1000 for j in range(1, 6)] for i in range(4)}
        result = batch_compute_metrics(read_sets, max_workers=2)
        assert result.total_items == 4
        assert result.successful == 4
        for key in read_sets:
            assert result.results[key].count == 5

    def test_batch_compute_metrics_realistic_longread_lengths(self) -> None:
        """batch_compute_metrics with realistic ONT read lengths."""
        rng = random.Random(42)
        # Simulate lognormal-distributed read lengths typical of ONT
        ont_lengths = [int(max(200, rng.lognormvariate(9.0, 1.0))) for _ in range(100)]
        pacbio_lengths = [int(max(500, rng.lognormvariate(9.5, 0.8))) for _ in range(100)]

        read_sets = {
            "ont_run": ont_lengths,
            "pacbio_run": pacbio_lengths,
        }

        result = batch_compute_metrics(read_sets)
        assert result.successful == 2

        ont_stats = result.results["ont_run"]
        assert ont_stats.count == 100
        assert ont_stats.n50 > 0
        assert ont_stats.mean_length > 0
        assert ont_stats.min_length >= 200

        pacbio_stats = result.results["pacbio_run"]
        assert pacbio_stats.count == 100
        assert pacbio_stats.n50 > 0
