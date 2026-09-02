"""Memory-hygiene helpers for campaign-scale frames and iterables.

Campaign-scale outputs (per-sample expression tables of 8k+ genes x 100+
samples, ortholog bridge tables with hundreds of thousands of rows) should be
processed in bounded chunks rather than materialized into intermediate copies
whole. These helpers provide API-stable, dependency-light chunking:

- ``chunked``: split any iterable into fixed-size lists.
- ``iter_frame_chunks``: lazily iterate a pandas DataFrame by row blocks
  without copying the frame (zero-copy ``iloc`` views per chunk).
- ``frame_memory_mb`` / ``frame_memory_report``: measure and report the
  memory footprint of a frame so callers can pick a chunk size deliberately.

All functions are side-effect free and use real data structures only (no
mocking); ``pandas`` is imported lazily so the iterable helpers work in
environments without it.
"""

from __future__ import annotations

import itertools
from typing import TYPE_CHECKING, Any, Iterator, Sequence, TypeVar

if TYPE_CHECKING:  # pragma: no cover - typing-only import
    import pandas as pd

T = TypeVar("T")

DEFAULT_CHUNK_ROWS = 50_000


def chunked(iterable: Sequence[T] | Iterator[T], size: int) -> Iterator[list[T]]:
    """Yield lists of at most ``size`` items from ``iterable``.

    Works with sequences and lazy iterators; the input iterator is consumed
    exactly once. The final chunk may be shorter than ``size``.

    Args:
        iterable: Items to chunk (list, generator, file handle, ...).
        size: Maximum items per chunk; must be >= 1.

    Yields:
        Lists of up to ``size`` items.

    Raises:
        ValueError: If ``size`` < 1.

    Example:
        >>> list(chunked([1, 2, 3, 4, 5], 2))
        [[1, 2], [3, 4], [5]]
    """
    if size < 1:
        raise ValueError(f"chunk size must be >= 1, got {size}")
    iterator = iter(iterable)
    while True:
        block = list(itertools.islice(iterator, size))
        if not block:
            return
        yield block


def iter_frame_chunks(
    frame: pd.DataFrame,
    rows_per_chunk: int = DEFAULT_CHUNK_ROWS,
) -> Iterator[pd.DataFrame]:
    """Lazily yield row-block views of ``frame`` without copying it.

    Each yielded chunk is an ``iloc`` slice (a view-backed frame), so peak
    memory stays close to O(chunk) rather than O(frame) beyond the original
    allocation. Mutating a yielded chunk may or may not propagate to the
    parent frame depending on pandas copy-on-write semantics; callers must
    treat chunks as read-only inputs.

    Args:
        frame: DataFrame to iterate in bounded row blocks.
        rows_per_chunk: Rows per chunk; must be >= 1.

    Yields:
        DataFrame slices of at most ``rows_per_chunk`` rows, in order.

    Raises:
        ValueError: If ``rows_per_chunk`` < 1.
        TypeError: If ``frame`` is not a pandas DataFrame.

    Example:
        for chunk in iter_frame_chunks(expression_df, rows_per_chunk=10_000):
            process(chunk)
    """
    import pandas as pd

    if not isinstance(frame, pd.DataFrame):
        raise TypeError(f"iter_frame_chunks expects a pandas DataFrame, got {type(frame).__name__}")
    if rows_per_chunk < 1:
        raise ValueError(f"rows_per_chunk must be >= 1, got {rows_per_chunk}")
    n_rows = len(frame)
    for start in range(0, n_rows, rows_per_chunk):
        yield frame.iloc[start : start + rows_per_chunk]


def frame_memory_mb(frame: pd.DataFrame) -> float:
    """Return the deep memory footprint of ``frame`` in MiB.

    Uses ``memory_usage(deep=True)`` so object-dtype columns (typical for
    transcript-ID and orthogroup tables) are measured by their real payload.

    Args:
        frame: pandas DataFrame to measure.

    Returns:
        Deep memory footprint in MiB (megabytes of 2**20 bytes).
    """
    return float(frame.memory_usage(deep=True).sum()) / (1024.0 * 1024.0)


def frame_memory_report(frame: pd.DataFrame) -> dict[str, Any]:
    """Build a descriptive memory report for a frame.

    Args:
        frame: pandas DataFrame to describe.

    Returns:
        Dict with ``n_rows``, ``n_cols``, ``dtypes`` (column -> dtype string),
        ``memory_mb`` (deep footprint), and ``memory_per_row_bytes``.
    """
    total_bytes = float(frame.memory_usage(deep=True).sum())
    n_rows = int(len(frame))
    return {
        "n_rows": n_rows,
        "n_cols": int(frame.shape[1]),
        "dtypes": {str(col): str(dtype) for col, dtype in frame.dtypes.items()},
        "memory_mb": total_bytes / (1024.0 * 1024.0),
        "memory_per_row_bytes": (total_bytes / n_rows) if n_rows else 0.0,
    }
