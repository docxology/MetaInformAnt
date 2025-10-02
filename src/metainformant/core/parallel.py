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
