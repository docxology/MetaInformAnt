from __future__ import annotations

from time import sleep
from metainformant.core import parallel as core_parallel


def square(x: int) -> int:
    return x * x


def test_thread_map_preserves_order() -> None:
    data = list(range(10))
    out = core_parallel.thread_map(square, data, max_workers=4)
    assert out == [x * x for x in data]


def test_thread_map_chunked() -> None:
    # Ensure chunking works without errors
    data = list(range(25))
    out = core_parallel.thread_map(square, data, max_workers=3, chunk_size=5)
    assert out == [x * x for x in data]


