from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor
from collections.abc import Callable, Iterable, Sequence
from typing import TypeVar


T = TypeVar("T")
U = TypeVar("U")


def thread_map(
    func: Callable[[T], U],
    items: Sequence[T] | Iterable[T],
    *,
    max_workers: int = 8,
    chunk_size: int | None = None,
) -> list[U]:
    """Map a function across items using threads, preserving order.

    If items is not a sequence, it will be materialized to preserve order.
    """
    if not isinstance(items, Sequence):
        items = list(items)
    results: list[U] = [None] * len(items)  # type: ignore[assignment]

    def _wrap(idx: int, x: T) -> None:
        results[idx] = func(x)

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
