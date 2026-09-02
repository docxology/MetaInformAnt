"""Zero-mock tests for metainformant.core.utils.batches."""

from __future__ import annotations

import pandas as pd
import pytest

from metainformant.core.utils import batches


def test_chunked_splits_evenly_and_handles_remainder() -> None:
    assert list(batches.chunked([1, 2, 3, 4, 5], 2)) == [[1, 2], [3, 4], [5]]
    assert list(batches.chunked([], 3)) == []
    assert list(batches.chunked("abc", 2)) == [["a", "b"], ["c"]]


def test_chunked_consumes_generator_once() -> None:
    gen = iter(range(7))
    out = list(batches.chunked(gen, 3))
    assert out == [[0, 1, 2], [3, 4, 5], [6]]


def test_chunked_rejects_nonpositive_size() -> None:
    with pytest.raises(ValueError, match="chunk size"):
        list(batches.chunked([1], 0))


def test_iter_frame_chunks_yields_ordered_row_blocks() -> None:
    frame = pd.DataFrame({"gene": [f"g{i}" for i in range(250)], "value": range(250)})
    chunks = list(batches.iter_frame_chunks(frame, rows_per_chunk=100))
    assert len(chunks) == 3
    assert [len(c) for c in chunks] == [100, 100, 50]
    combined = pd.concat(chunks, ignore_index=True)
    pd.testing.assert_frame_equal(combined, frame)


def test_iter_frame_chunks_empty_frame_yields_nothing() -> None:
    frame = pd.DataFrame({"a": []})
    assert list(batches.iter_frame_chunks(frame, rows_per_chunk=10)) == []


def test_iter_frame_chunks_validates_inputs() -> None:
    with pytest.raises(ValueError, match="rows_per_chunk"):
        list(batches.iter_frame_chunks(pd.DataFrame({"a": [1]}), rows_per_chunk=0))
    with pytest.raises(TypeError, match="DataFrame"):
        list(batches.iter_frame_chunks([[1, 2]], rows_per_chunk=1))  # type: ignore[arg-type]


def test_frame_memory_mb_and_report_are_self_consistent() -> None:
    frame = pd.DataFrame(
        {
            "gene": [f"gene_{i}" for i in range(1_000)],
            "count": range(1_000),
        }
    )
    mb = batches.frame_memory_mb(frame)
    assert mb > 0.0
    report = batches.frame_memory_report(frame)
    assert report["n_rows"] == 1_000
    assert report["n_cols"] == 2
    assert report["memory_mb"] == pytest.approx(mb)
    assert report["memory_per_row_bytes"] > 0
    assert report["dtypes"]["count"] == "int64"
