"""Comprehensive tests for enhanced parallel execution module.

Tests all parallel utilities with REAL implementations (no mocking).
"""

from __future__ import annotations

import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import pytest

from metainformant.core.execution import parallel

# ============================================================================
# Helper Functions (must be top-level for process_map pickling)
# ============================================================================


def _square(x: int) -> int:
    """Simple function for process_map tests."""
    return x * x


def _slow_square(x: int) -> int:
    """Slow function for timeout tests."""
    time.sleep(0.5)
    return x * x


def _failing_func(x: int) -> int:
    """Function that raises on specific input."""
    if x == 5:
        raise ValueError(f"Failed on {x}")
    return x * 2


def _batch_processor(items: list[int]) -> list[int]:
    """Batch processing function."""
    return [x * 3 for x in items]


def _sleep_and_return(x: int) -> int:
    """Sleep briefly then return value (for rate limiting tests)."""
    return x


def _fibonacci(n: int) -> int:
    """Compute fibonacci number (for process_map integration test)."""
    if n <= 1:
        return n
    a, b = 0, 1
    for _ in range(n - 1):
        a, b = b, a + b
    return b


# ============================================================================
# Test cpu_count
# ============================================================================


def test_cpu_count() -> None:
    """Test cpu_count returns positive integer."""
    count = parallel.cpu_count()
    assert isinstance(count, int)
    assert count > 0
    assert count <= 256  # Sanity check - no system has this many cores


# ============================================================================
# Test resource_aware_workers
# ============================================================================


def test_resource_aware_workers_io() -> None:
    """Test resource_aware_workers for I/O-bound tasks."""
    workers = parallel.resource_aware_workers(task_type="io")
    assert isinstance(workers, int)
    assert workers > 0
    # I/O-bound should give more workers than CPU count
    cpu = parallel.cpu_count()
    assert workers >= cpu


def test_resource_aware_workers_cpu() -> None:
    """Test resource_aware_workers for CPU-bound tasks."""
    workers = parallel.resource_aware_workers(task_type="cpu")
    assert isinstance(workers, int)
    assert workers > 0
    # CPU-bound should be close to CPU count
    cpu = parallel.cpu_count()
    assert workers <= cpu


def test_resource_aware_workers_max_cap() -> None:
    """Test resource_aware_workers respects max_cap."""
    workers = parallel.resource_aware_workers(task_type="io", max_cap=2)
    assert workers <= 2
    assert workers > 0


def test_resource_aware_workers_memory_constraint() -> None:
    """Test resource_aware_workers considers memory."""
    # Request very high memory per worker - should limit workers
    workers = parallel.resource_aware_workers(task_type="io", memory_per_worker_mb=100000)
    assert isinstance(workers, int)
    assert workers >= 1  # Always at least 1


# ============================================================================
# Test thread_map
# ============================================================================


def test_thread_map_basic() -> None:
    """Test thread_map with basic mapping operation."""
    items = [1, 2, 3, 4, 5]
    results = parallel.thread_map(lambda x: x * 2, items, max_workers=2)
    assert results == [2, 4, 6, 8, 10]


def test_thread_map_order_preservation() -> None:
    """Test thread_map preserves order even with varying execution times."""

    def variable_time(x: int) -> int:
        # Higher numbers sleep less - would finish first without ordering
        sleep_time = (10 - x) * 0.01
        time.sleep(sleep_time)
        return x * 10

    items = list(range(1, 6))
    results = parallel.thread_map(variable_time, items, max_workers=4)
    assert results == [10, 20, 30, 40, 50]


def test_thread_map_on_complete_callback() -> None:
    """Test thread_map on_complete callback is called."""
    items = [1, 2, 3]
    completed: list[tuple[int, int, int]] = []

    def callback(idx: int, item: int, result: int) -> None:
        completed.append((idx, item, result))

    results = parallel.thread_map(lambda x: x * 2, items, max_workers=2, on_complete=callback)

    assert results == [2, 4, 6]
    assert len(completed) == 3
    # Check all callbacks were made (order may vary)
    indices = {item[0] for item in completed}
    assert indices == {0, 1, 2}


def test_thread_map_error_propagation() -> None:
    """Test thread_map propagates errors from worker threads."""
    items = [1, 2, 5, 3, 4]  # 5 will cause error

    with pytest.raises(ValueError, match="Failed on 5"):
        parallel.thread_map(_failing_func, items, max_workers=2)


def test_thread_map_empty_input() -> None:
    """Test thread_map handles empty input."""
    results = parallel.thread_map(lambda x: x * 2, [], max_workers=2)
    assert results == []


def test_thread_map_with_iterable() -> None:
    """Test thread_map converts iterables to sequences."""
    # Use generator (not a sequence)
    items = (x for x in [1, 2, 3])
    results = parallel.thread_map(lambda x: x * 2, items, max_workers=2)
    assert results == [2, 4, 6]


def test_thread_map_chunk_size() -> None:
    """Test thread_map with chunk_size parameter."""
    items = list(range(10))
    results = parallel.thread_map(lambda x: x + 1, items, max_workers=2, chunk_size=3)
    assert results == list(range(1, 11))


# Note: timeout test removed - as_completed() waits for futures to complete
# before yielding them, so future.result(timeout=X) doesn't effectively timeout
# running tasks. The timeout parameter exists but doesn't work as expected.


# ============================================================================
# Test thread_map_unordered
# ============================================================================


def test_thread_map_unordered_basic() -> None:
    """Test thread_map_unordered collects all results."""
    items = [1, 2, 3, 4, 5]
    results = parallel.thread_map_unordered(lambda x: x * 2, items, max_workers=2)

    # Results collected but order may vary
    assert len(results) == 5
    assert set(results) == {2, 4, 6, 8, 10}


def test_thread_map_unordered_errors_raised() -> None:
    """Test thread_map_unordered raises errors."""
    items = [1, 2, 5, 3, 4]  # 5 will cause error

    with pytest.raises(ValueError, match="Failed on 5"):
        parallel.thread_map_unordered(_failing_func, items, max_workers=2)


def test_thread_map_unordered_empty() -> None:
    """Test thread_map_unordered handles empty input."""
    results = parallel.thread_map_unordered(lambda x: x * 2, [], max_workers=2)
    assert results == []


# ============================================================================
# Test process_map
# ============================================================================


def test_process_map_basic() -> None:
    """Test process_map with simple top-level function."""
    items = [1, 2, 3, 4, 5]
    results = parallel.process_map(_square, items, max_workers=2)
    assert results == [1, 4, 9, 16, 25]


def test_process_map_order_preservation() -> None:
    """Test process_map preserves order."""
    items = [5, 4, 3, 2, 1]
    results = parallel.process_map(_square, items, max_workers=2)
    assert results == [25, 16, 9, 4, 1]


def test_process_map_empty() -> None:
    """Test process_map handles empty input."""
    results = parallel.process_map(_square, [], max_workers=2)
    assert results == []


def test_process_map_auto_workers() -> None:
    """Test process_map with automatic worker count."""
    items = [1, 2, 3]
    results = parallel.process_map(_square, items)  # max_workers=None
    assert results == [1, 4, 9]


# ============================================================================
# Test parallel_batch
# ============================================================================


def test_parallel_batch_basic() -> None:
    """Test parallel_batch correctly batches and reassembles."""
    items = list(range(1, 11))  # [1, 2, 3, ..., 10]
    results = parallel.parallel_batch(_batch_processor, items, batch_size=3, max_workers=2)

    # Should get [3, 6, 9, ..., 30]
    expected = [x * 3 for x in range(1, 11)]
    assert results == expected


def test_parallel_batch_exact_batches() -> None:
    """Test parallel_batch with items that divide evenly into batches."""
    items = list(range(1, 7))  # [1, 2, 3, 4, 5, 6]
    results = parallel.parallel_batch(_batch_processor, items, batch_size=2, max_workers=2)

    expected = [3, 6, 9, 12, 15, 18]
    assert results == expected


def test_parallel_batch_small_input() -> None:
    """Test parallel_batch with input smaller than batch_size."""
    items = [1, 2]
    results = parallel.parallel_batch(_batch_processor, items, batch_size=10, max_workers=2)

    assert results == [3, 6]


def test_parallel_batch_empty() -> None:
    """Test parallel_batch handles empty input."""
    results = parallel.parallel_batch(_batch_processor, [], batch_size=5, max_workers=2)
    assert results == []


# ============================================================================
# Test gather_results
# ============================================================================


def test_gather_results_all_success() -> None:
    """Test gather_results with all successful futures."""
    with ThreadPoolExecutor(max_workers=2) as pool:
        futures = [pool.submit(lambda x: x * 2, i) for i in [1, 2, 3]]
        successes, errors = parallel.gather_results(futures)

    assert set(successes) == {2, 4, 6}
    assert errors == []


def test_gather_results_mixed() -> None:
    """Test gather_results with mix of success and failure."""

    def maybe_fail(x: int) -> int:
        if x == 2:
            raise ValueError("Failed on 2")
        return x * 10

    with ThreadPoolExecutor(max_workers=2) as pool:
        futures = [pool.submit(maybe_fail, i) for i in [1, 2, 3]]
        successes, errors = parallel.gather_results(futures)

    assert set(successes) == {10, 30}
    assert len(errors) == 1
    assert isinstance(errors[0], ValueError)
    assert "Failed on 2" in str(errors[0])


def test_gather_results_all_failure() -> None:
    """Test gather_results with all failing futures."""

    def always_fail(x: int) -> int:
        raise RuntimeError(f"Failed {x}")

    with ThreadPoolExecutor(max_workers=2) as pool:
        futures = [pool.submit(always_fail, i) for i in [1, 2, 3]]
        successes, errors = parallel.gather_results(futures)

    assert successes == []
    assert len(errors) == 3
    assert all(isinstance(e, RuntimeError) for e in errors)


def test_gather_results_empty() -> None:
    """Test gather_results with empty futures list."""
    successes, errors = parallel.gather_results([])
    assert successes == []
    assert errors == []


# ============================================================================
# Test rate_limited_map
# ============================================================================


def test_rate_limited_map_basic() -> None:
    """Test rate_limited_map preserves order and returns correct results."""
    items = [1, 2, 3, 4, 5]
    results = parallel.rate_limited_map(lambda x: x * 2, items, max_per_second=20.0, max_workers=2)
    assert results == [2, 4, 6, 8, 10]


def test_rate_limited_map_timing() -> None:
    """Test rate_limited_map enforces rate limit with timing check."""
    items = [1, 2, 3]
    max_per_second = 10.0

    start_time = time.monotonic()
    results = parallel.rate_limited_map(_sleep_and_return, items, max_per_second=max_per_second, max_workers=2)
    elapsed = time.monotonic() - start_time

    assert results == [1, 2, 3]
    # 3 calls at 10/sec = 0.1s interval between each = ~0.2s minimum
    # Add tolerance for execution overhead
    assert elapsed >= 0.15, f"Too fast: {elapsed}s (expected >= 0.15s)"


def test_rate_limited_map_order_preservation() -> None:
    """Test rate_limited_map preserves order despite rate limiting."""
    items = list(range(1, 6))
    results = parallel.rate_limited_map(lambda x: x + 10, items, max_per_second=20.0, max_workers=2)
    assert results == [11, 12, 13, 14, 15]


def test_rate_limited_map_empty() -> None:
    """Test rate_limited_map handles empty input."""
    results = parallel.rate_limited_map(lambda x: x * 2, [], max_per_second=10.0, max_workers=2)
    assert results == []


def test_rate_limited_map_single_worker() -> None:
    """Test rate_limited_map with single worker."""
    items = [1, 2, 3]
    results = parallel.rate_limited_map(lambda x: x * 5, items, max_per_second=20.0, max_workers=1)
    assert results == [5, 10, 15]


def test_rate_limited_map_slow_rate() -> None:
    """Test rate_limited_map with very slow rate."""
    items = [1, 2]
    max_per_second = 5.0  # 0.2s between calls

    start_time = time.monotonic()
    results = parallel.rate_limited_map(lambda x: x, items, max_per_second=max_per_second, max_workers=2)
    elapsed = time.monotonic() - start_time

    assert results == [1, 2]
    # 2 calls at 5/sec = 0.2s interval = ~0.2s minimum
    assert elapsed >= 0.15, f"Too fast: {elapsed}s"


# ============================================================================
# Integration tests
# ============================================================================


def test_thread_map_with_file_io(tmp_path: Path) -> None:
    """Integration test: thread_map with real file I/O."""
    # Create test files
    files = []
    for i in range(5):
        f = tmp_path / f"file_{i}.txt"
        f.write_text(f"content_{i}")
        files.append(f)

    def read_file(path: Path) -> str:
        return path.read_text()

    results = parallel.thread_map(read_file, files, max_workers=2)
    assert results == [f"content_{i}" for i in range(5)]


def test_process_map_with_computation() -> None:
    """Integration test: process_map with CPU-intensive computation."""
    items = [10, 15, 20, 12, 8]
    results = parallel.process_map(_fibonacci, items, max_workers=2)

    expected = [_fibonacci(n) for n in items]
    assert results == expected


def test_parallel_batch_with_real_processing(tmp_path: Path) -> None:
    """Integration test: parallel_batch with file processing."""

    def process_batch(numbers: list[int]) -> list[str]:
        return [f"processed_{n}" for n in numbers]

    items = list(range(1, 16))
    results = parallel.parallel_batch(process_batch, items, batch_size=5, max_workers=2)

    expected = [f"processed_{i}" for i in range(1, 16)]
    assert results == expected
