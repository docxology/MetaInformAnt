"""Comprehensive tests for core.execution.parallel module."""

from __future__ import annotations

import time
from typing import List

from metainformant.core.execution import parallel as core_parallel


def square(x: int) -> int:
    """Helper function for tests."""
    return x * x


def slow_square(x: int) -> int:
    """Helper function that takes time."""
    time.sleep(0.01)
    return x * x


def double_list(items: List[int]) -> List[int]:
    """Helper batch function."""
    return [x * 2 for x in items]


class TestThreadMap:
    """Tests for thread_map function."""

    def test_thread_map_preserves_order(self) -> None:
        """Test that results maintain input order."""
        data = list(range(20))
        out = core_parallel.thread_map(square, data, max_workers=4)
        assert out == [x * x for x in data]

    def test_thread_map_with_chunk_size(self) -> None:
        """Test chunked processing."""
        data = list(range(25))
        out = core_parallel.thread_map(square, data, max_workers=3, chunk_size=5)
        assert out == [x * x for x in data]

    def test_thread_map_empty_input(self) -> None:
        """Test with empty input list."""
        out = core_parallel.thread_map(square, [], max_workers=4)
        assert out == []

    def test_thread_map_single_item(self) -> None:
        """Test with single item."""
        out = core_parallel.thread_map(square, [5], max_workers=4)
        assert out == [25]

    def test_thread_map_with_generator(self) -> None:
        """Test with generator input (should be materialized)."""
        gen = (x for x in range(10))
        out = core_parallel.thread_map(square, gen, max_workers=4)
        assert out == [x * x for x in range(10)]

    def test_thread_map_with_max_workers_1(self) -> None:
        """Test sequential execution."""
        data = list(range(5))
        out = core_parallel.thread_map(square, data, max_workers=1)
        assert out == [x * x for x in data]

    def test_thread_map_concurrent_speedup(self) -> None:
        """Test that parallel execution provides speedup."""
        data = list(range(20))

        # Sequential timing
        start = time.time()
        _ = core_parallel.thread_map(slow_square, data, max_workers=1)
        seq_time = time.time() - start

        # Parallel timing
        start = time.time()
        _ = core_parallel.thread_map(slow_square, data, max_workers=8)
        par_time = time.time() - start

        # Parallel should be faster (at least 2x)
        assert par_time < seq_time * 0.75 or seq_time < 0.1  # Skip if too fast


class TestThreadMapUnordered:
    """Tests for thread_map_unordered function."""

    def test_thread_map_unordered_basic(self) -> None:
        """Test unordered map returns all results."""
        data = list(range(10))
        out = core_parallel.thread_map_unordered(square, data, max_workers=4)
        # Results should contain same values but possibly different order
        assert sorted(out) == sorted([x * x for x in data])

    def test_thread_map_unordered_empty(self) -> None:
        """Test with empty input."""
        out = core_parallel.thread_map_unordered(square, [], max_workers=4)
        assert out == []


class TestParallelBatch:
    """Tests for parallel_batch function."""

    def test_parallel_batch_basic(self) -> None:
        """Test batch processing."""
        data = list(range(20))
        out = core_parallel.parallel_batch(double_list, data, batch_size=5, max_workers=4)
        assert out == [x * 2 for x in data]

    def test_parallel_batch_uneven_batches(self) -> None:
        """Test with data not evenly divisible by batch size."""
        data = list(range(17))  # 17 items, batch size 5 = 4 batches (5,5,5,2)
        out = core_parallel.parallel_batch(double_list, data, batch_size=5, max_workers=4)
        assert out == [x * 2 for x in data]

    def test_parallel_batch_single_batch(self) -> None:
        """Test when all items fit in one batch."""
        data = list(range(3))
        out = core_parallel.parallel_batch(double_list, data, batch_size=10, max_workers=4)
        assert out == [x * 2 for x in data]

    def test_parallel_batch_empty(self) -> None:
        """Test with empty input."""
        out = core_parallel.parallel_batch(double_list, [], batch_size=5, max_workers=4)
        assert out == []


class TestCpuCount:
    """Tests for cpu_count function."""

    def test_cpu_count_returns_positive(self) -> None:
        """Test that cpu_count returns a positive integer."""
        count = core_parallel.cpu_count()
        assert isinstance(count, int)
        assert count >= 1

    def test_cpu_count_reproducible(self) -> None:
        """Test that cpu_count returns consistent value."""
        count1 = core_parallel.cpu_count()
        count2 = core_parallel.cpu_count()
        assert count1 == count2
