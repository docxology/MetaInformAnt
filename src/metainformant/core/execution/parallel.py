from __future__ import annotations

from collections.abc import Callable, Iterable, Sequence
from concurrent.futures import ThreadPoolExecutor
from typing import TypeVar

T = TypeVar("T")
U = TypeVar("U")


def thread_map(
    func: Callable[[T], U],
    items: Sequence[T] | Iterable[T],
    *,
    max_workers: int = 8,
    chunk_size: int | None = None,
    timeout: float | None = None,
    ordered: bool = True,
) -> list[U]:
    """Map a function across items using threads, preserving order.

    Args:
        func: Function to apply to each item
        items: Items to process (will be materialized if not a sequence)
        max_workers: Maximum number of worker threads
        chunk_size: Size of chunks for batch processing (None for auto)
        timeout: Timeout in seconds for each task (None for no timeout)
        ordered: Whether to preserve input order in results

    Returns:
        List of results in the same order as input items

    If items is not a sequence, it will be materialized to preserve order.
    """
    if not isinstance(items, Sequence):
        items = list(items)
    results: list[U] = [None] * len(items)  # type: ignore[assignment]

    def _wrap(idx: int, x: T) -> None:
        try:
            if timeout is not None:
                import signal
                def timeout_handler(signum, frame):
                    raise TimeoutError(f"Task timed out after {timeout} seconds")
                old_handler = signal.signal(signal.SIGALRM, timeout_handler)
                signal.alarm(int(timeout))
                try:
                    results[idx] = func(x)
                finally:
                    signal.alarm(0)
                    signal.signal(signal.SIGALRM, old_handler)
            else:
                results[idx] = func(x)
        except Exception as e:
            results[idx] = e  # Store exception for later
            raise

    with ThreadPoolExecutor(max_workers=max_workers) as pool:
        if chunk_size and chunk_size > 1:
            for start in range(0, len(items), chunk_size):
                end = min(start + chunk_size, len(items))
                for i in range(start, end):
                    pool.submit(_wrap, i, items[i])
            pool.shutdown(wait=True)
        else:
            for i, x in enumerate(items):
                pool.submit(_wrap, i, x)
            pool.shutdown(wait=True)
    return results


def thread_map_unordered(
    func: Callable[[T], U],
    items: Iterable[T],
    *,
    max_workers: int = 8,
    timeout: float | None = None,
) -> list[U]:
    """Map a function across items using threads, without preserving order.

    Args:
        func: Function to apply to each item
        items: Items to process
        max_workers: Maximum number of worker threads
        timeout: Timeout in seconds for each task (None for no timeout)

    Returns:
        List of results (order not preserved for better performance)
    """
    results = []

    def collect_result(future):
        try:
            if timeout is not None:
                import signal
                def timeout_handler(signum, frame):
                    raise TimeoutError(f"Task timed out after {timeout} seconds")
                old_handler = signal.signal(signal.SIGALRM, timeout_handler)
                signal.alarm(int(timeout))
                try:
                    result = future.result()
                    results.append(result)
                finally:
                    signal.alarm(0)
                    signal.signal(signal.SIGALRM, old_handler)
            else:
                result = future.result()
                results.append(result)
        except Exception as e:
            results.append(e)
            raise

    with ThreadPoolExecutor(max_workers=max_workers) as pool:
        futures = [pool.submit(func, item) for item in items]
        for future in futures:
            future.add_done_callback(lambda f: collect_result(f))
        pool.shutdown(wait=True)

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

    Returns:
        List of all results
    """
    batches = [items[i:i + batch_size] for i in range(0, len(items), batch_size)]
    batch_results = thread_map(func, batches, max_workers=max_workers)
    return [result for batch_result in batch_results for result in batch_result]


def cpu_count() -> int:
    """Get the number of CPU cores available."""
    try:
        import multiprocessing
        return multiprocessing.cpu_count()
    except ImportError:
        # Fallback for systems without multiprocessing
        import os
        return os.cpu_count() or 1
