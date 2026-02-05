from __future__ import annotations

import multiprocessing
import os
import time
import threading
from collections.abc import Callable, Iterable, Sequence
from concurrent.futures import (
    Future,
    ProcessPoolExecutor,
    ThreadPoolExecutor,
    as_completed,
)
from typing import TypeVar

T = TypeVar("T")
U = TypeVar("U")


def cpu_count() -> int:
    """Get the number of CPU cores available."""
    return os.cpu_count() or multiprocessing.cpu_count()


def resource_aware_workers(
    *,
    task_type: str = "io",
    max_cap: int | None = None,
    memory_per_worker_mb: int = 256,
) -> int:
    """Return recommended worker count based on CPU cores and available memory.

    Args:
        task_type: "io" for I/O-bound (more workers) or "cpu" for CPU-bound (match cores)
        max_cap: Hard upper limit on workers (None for no cap)
        memory_per_worker_mb: Estimated memory per worker in MB
    """
    cores = cpu_count()

    if task_type == "cpu":
        # CPU-bound: use core count minus 1, leave headroom for OS
        workers = max(1, cores - 1)
    else:
        # I/O-bound: 2-4x cores is reasonable
        workers = min(cores * 4, 32)

    # Memory constraint: check available system memory
    try:
        import psutil

        available_mb = psutil.virtual_memory().available / (1024 * 1024)
        mem_limited = max(1, int(available_mb * 0.7 / memory_per_worker_mb))
        workers = min(workers, mem_limited)
    except ImportError:
        pass  # psutil not available, skip memory check

    if max_cap is not None:
        workers = min(workers, max_cap)

    return max(1, workers)


def thread_map(
    func: Callable[[T], U],
    items: Sequence[T] | Iterable[T],
    *,
    max_workers: int = 8,
    chunk_size: int | None = None,
    timeout: float | None = None,
    ordered: bool = True,
    on_complete: Callable[[int, T, U], None] | None = None,
) -> list[U]:
    """Map a function across items using threads, preserving order.

    Args:
        func: Function to apply to each item
        items: Items to process (will be materialized if not a sequence)
        max_workers: Maximum number of worker threads
        chunk_size: Size of chunks for batch processing (None for auto)
        timeout: Timeout in seconds for each task (None for no timeout)
        ordered: Whether to preserve input order in results
        on_complete: Optional callback(index, input_item, result) called after each task completes
    """
    if not isinstance(items, Sequence):
        items = list(items)
    if not items:
        return []

    results: list[U] = [None] * len(items)  # type: ignore[assignment]
    errors: list[tuple[int, Exception]] = []

    with ThreadPoolExecutor(max_workers=max_workers) as pool:
        # Submit all futures, tracking index -> future mapping
        future_to_idx: dict[Future[U], int] = {}
        if chunk_size and chunk_size > 1:
            for start in range(0, len(items), chunk_size):
                end = min(start + chunk_size, len(items))
                for i in range(start, end):
                    future = pool.submit(func, items[i])
                    future_to_idx[future] = i
        else:
            for i, x in enumerate(items):
                future = pool.submit(func, x)
                future_to_idx[future] = i

        # Collect results using Future.result(timeout=) -- thread-safe, cross-platform
        for future in as_completed(future_to_idx):
            idx = future_to_idx[future]
            try:
                result = future.result(timeout=timeout)
                results[idx] = result
                if on_complete is not None:
                    on_complete(idx, items[idx], result)
            except TimeoutError:
                err = TimeoutError(f"Task {idx} timed out after {timeout}s")
                errors.append((idx, err))
                results[idx] = err  # type: ignore[assignment]
            except Exception as e:
                errors.append((idx, e))
                results[idx] = e  # type: ignore[assignment]

    if errors:
        # Raise the first error to preserve existing behavior
        raise errors[0][1]

    return results


def thread_map_unordered(
    func: Callable[[T], U],
    items: Iterable[T],
    *,
    max_workers: int = 8,
    timeout: float | None = None,
) -> list[U]:
    """Map a function across items using threads, without preserving order.

    Uses as_completed for correct result collection (no race conditions).

    Args:
        func: Function to apply to each item
        items: Items to process
        max_workers: Maximum number of worker threads
        timeout: Timeout in seconds for each task (None for no timeout)
    """
    results: list[U] = []

    with ThreadPoolExecutor(max_workers=max_workers) as pool:
        futures = [pool.submit(func, item) for item in items]

        for future in as_completed(futures):
            try:
                result = future.result(timeout=timeout)
                results.append(result)
            except Exception as e:
                results.append(e)  # type: ignore[arg-type]
                raise

    return results


def process_map(
    func: Callable[[T], U],
    items: Sequence[T] | Iterable[T],
    *,
    max_workers: int | None = None,
    timeout: float | None = None,
    ordered: bool = True,
) -> list[U]:
    """Map a function across items using processes for CPU-bound work.

    Args:
        func: Function to apply (must be picklable -- top-level or classmethod)
        items: Items to process
        max_workers: Number of worker processes (None = cpu_count - 1)
        timeout: Timeout in seconds for each task
        ordered: Whether to preserve input order
    """
    if not isinstance(items, Sequence):
        items = list(items)
    if not items:
        return []

    if max_workers is None:
        max_workers = max(1, cpu_count() - 1)

    results: list[U] = [None] * len(items)  # type: ignore[assignment]

    with ProcessPoolExecutor(max_workers=max_workers) as pool:
        future_to_idx: dict[Future[U], int] = {}
        for i, x in enumerate(items):
            future = pool.submit(func, x)
            future_to_idx[future] = i

        for future in as_completed(future_to_idx):
            idx = future_to_idx[future]
            result = future.result(timeout=timeout)
            results[idx] = result

    return results


def parallel_batch(
    func: Callable[[list[T]], list[U]],
    items: list[T],
    *,
    batch_size: int = 10,
    max_workers: int = 4,
) -> list[U]:
    """Process items in batches using parallel execution.

    Args:
        func: Function that takes a batch of items and returns processed results
        items: Items to process
        batch_size: Number of items per batch
        max_workers: Maximum number of worker threads
    """
    batches = [items[i : i + batch_size] for i in range(0, len(items), batch_size)]
    batch_results = thread_map(func, batches, max_workers=max_workers)
    return [result for batch_result in batch_results for result in batch_result]


def gather_results(
    futures: Sequence[Future[U]],
    *,
    timeout: float | None = None,
) -> tuple[list[U], list[Exception]]:
    """Collect results from multiple futures, separating successes from errors.

    Args:
        futures: Futures to collect results from
        timeout: Timeout per future in seconds

    Returns:
        Tuple of (successes, errors) -- successes in completion order
    """
    successes: list[U] = []
    errors: list[Exception] = []

    for future in as_completed(futures):
        try:
            result = future.result(timeout=timeout)
            successes.append(result)
        except Exception as e:
            errors.append(e)

    return successes, errors


def rate_limited_map(
    func: Callable[[T], U],
    items: Sequence[T] | Iterable[T],
    *,
    max_per_second: float = 10.0,
    max_workers: int = 4,
    timeout: float | None = None,
) -> list[U]:
    """Map with rate limiting -- max N calls per second. Useful for API rate limits.

    Preserves input order. Rate limiting is enforced at submission time via a token bucket.

    Args:
        func: Function to apply
        items: Items to process
        max_per_second: Maximum function invocations per second
        max_workers: Maximum concurrent workers
        timeout: Timeout per task in seconds
    """
    if not isinstance(items, Sequence):
        items = list(items)
    if not items:
        return []

    interval = 1.0 / max_per_second
    results: list[U] = [None] * len(items)  # type: ignore[assignment]

    # Token bucket for rate limiting submission
    lock = threading.Lock()
    last_submit_time = [0.0]

    def rate_limited_func(idx_and_item: tuple[int, T]) -> tuple[int, U]:
        idx, item = idx_and_item
        # Enforce rate limit before executing
        with lock:
            now = time.monotonic()
            wait = interval - (now - last_submit_time[0])
            if wait > 0:
                time.sleep(wait)
            last_submit_time[0] = time.monotonic()
        return idx, func(item)

    indexed_items: list[tuple[int, T]] = list(enumerate(items))

    with ThreadPoolExecutor(max_workers=max_workers) as pool:
        future_to_idx: dict[Future[tuple[int, U]], int] = {}
        for pair in indexed_items:
            future = pool.submit(rate_limited_func, pair)
            future_to_idx[future] = pair[0]

        for future in as_completed(future_to_idx):
            idx, result = future.result(timeout=timeout)
            results[idx] = result

    return results
