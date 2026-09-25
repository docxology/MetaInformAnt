from __future__ import annotations

import multiprocessing
import os
import threading
import time
from collections.abc import Callable, Iterable, Sequence
from concurrent.futures import (
    Future,
    ProcessPoolExecutor,
    ThreadPoolExecutor,
    as_completed,
)
from typing import TypeVar, cast

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
        if task_type != "cpu" and memory_per_worker_mb <= 256:
            # Default I/O workloads should keep at least one worker per core.
            # Callers can still request strict memory limiting by raising the
            # per-worker estimate.
            workers = max(workers, cores)
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
    """Map a function across items using threads.

    Args:
        func: Function to apply to each item
        items: Items to process (will be materialized if not a sequence)
        max_workers: Maximum number of worker threads
        chunk_size: Size of chunks for batch processing (None for auto)
        timeout: Maximum wall-clock seconds to wait for ALL tasks; when the
            deadline passes, queued tasks are cancelled and ``TimeoutError`` is
            raised (None disables the deadline)
        ordered: ``True`` (default) returns results in input order; ``False``
            returns results in completion order
        on_complete: Optional callback(index, input_item, result) called after
            each successful task completes; index is the input position

    Raises:
        TimeoutError: If ``timeout`` elapses before every task finishes.
    """
    if not isinstance(items, Sequence):
        items = list(items)
    if not items:
        return []

    ordered_results: list[U] = cast("list[U]", [None] * len(items))
    completion_results: list[U] = []
    errors: list[tuple[int, Exception]] = []
    finished = 0

    pool = ThreadPoolExecutor(max_workers=max_workers)
    try:
        # Submit all futures, tracking future -> input index mapping
        future_to_idx: dict[Future[U], int] = {}
        if chunk_size and chunk_size > 1:
            for start in range(0, len(items), chunk_size):
                end = min(start + chunk_size, len(items))
                for i in range(start, end):
                    future_to_idx[pool.submit(func, items[i])] = i
        else:
            for i, item in enumerate(items):
                future_to_idx[pool.submit(func, item)] = i

        try:
            # as_completed yields futures in completion order and enforces the
            # overall wall-clock deadline (raises TimeoutError when exceeded).
            for future in as_completed(future_to_idx, timeout=timeout):
                idx = future_to_idx[future]
                try:
                    result = future.result()
                except Exception as e:
                    errors.append((idx, e))
                    ordered_results[idx] = cast("U", e)
                else:
                    if ordered:
                        ordered_results[idx] = result
                    else:
                        completion_results.append(result)
                    if on_complete is not None:
                        on_complete(idx, items[idx], result)
                finished += 1
        except TimeoutError:
            # Deadline exceeded: drop queued work instead of blocking on it.
            pool.shutdown(wait=False, cancel_futures=True)
            raise TimeoutError(
                f"thread_map timed out after {timeout}s: "
                f"{finished} of {len(items)} task(s) finished"
            ) from None
        pool.shutdown(wait=True)
    except BaseException:
        # Loop aborted early (timeout above or callback failure): avoid
        # blocking on queued work that will never be consumed.
        pool.shutdown(wait=False, cancel_futures=True)
        raise

    if errors:
        # Raise the first error to preserve existing behavior
        raise errors[0][1]

    return ordered_results if ordered else completion_results


def thread_map_unordered(
    func: Callable[[T], U],
    items: Iterable[T],
    *,
    max_workers: int = 8,
    timeout: float | None = None,
) -> list[U]:
    """Map a function across items using threads, without preserving order.

    Results are returned in completion order.

    Args:
        func: Function to apply to each item
        items: Items to process
        max_workers: Maximum number of worker threads
        timeout: Maximum wall-clock seconds to wait for ALL tasks; when the
            deadline passes, queued tasks are cancelled and ``TimeoutError`` is
            raised (None disables the deadline)

    Raises:
        TimeoutError: If ``timeout`` elapses before every task finishes.
    """
    results: list[U] = []

    pool = ThreadPoolExecutor(max_workers=max_workers)
    try:
        futures = [pool.submit(func, item) for item in items]
        try:
            for future in as_completed(futures, timeout=timeout):
                results.append(future.result())
        except TimeoutError:
            pool.shutdown(wait=False, cancel_futures=True)
            raise TimeoutError(
                f"thread_map_unordered timed out after {timeout}s: "
                f"{len(results)} of {len(futures)} task(s) finished"
            ) from None
        pool.shutdown(wait=True)
    except BaseException:
        pool.shutdown(wait=False, cancel_futures=True)
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

    # None-padded placeholder slots; every slot is overwritten before return.
    results: list[U] = cast("list[U]", [None] * len(items))

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
    # None-padded placeholder slots; every slot is overwritten before return.
    results: list[U] = cast("list[U]", [None] * len(items))

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
