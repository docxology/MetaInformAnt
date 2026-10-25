from __future__ import annotations

import gzip
import time
from pathlib import Path

import pytest

from metainformant.rna.retrieval.ena_staging import (
    StagingDeadline,
    validate_fastq_records,
)


@pytest.mark.parametrize(
    "payload,expected",
    [
        (b"@r\nACGT\n+\n!!!!\n", True),
        (b"@r\r\nACGT\r\n+\r\n!!!!\r\n", True),
        (b"", False),
        (b"@\nACGT\n+\n!!!!\n", False),
        (b"@/1\nACGT\n+\n!!!!\n", False),
        (b"r\nACGT\n+\n!!!!\n", False),
        (b"@r\nACGT\nwrong\n!!!!\n", False),
        (b"@r\nACGT\n+\n!!!\n", False),
        (b"@r\nACGT\n+\n!!!!", False),
        (b"@r\nACGT\n+\n", False),
    ],
)
def test_record_shapes_with_real_gzip_and_plain_files(
    tmp_path: Path, payload: bytes, expected: bool
) -> None:
    compressed, plain = tmp_path / "reads.fastq.gz", tmp_path / "reads.fastq"
    compressed.write_bytes(gzip.compress(payload, mtime=0))
    plain.write_bytes(payload)
    assert validate_fastq_records([compressed]) is expected
    assert validate_fastq_records([plain]) is expected


@pytest.mark.parametrize(
    "header_length", [1024 * 1024 - 30, 1024 * 1024 - 2, 1024 * 1024 + 2]
)
def test_records_cross_chunk_boundaries_and_enforce_line_bound(
    tmp_path: Path, header_length: int
) -> None:
    first, second = tmp_path / "one.fastq.gz", tmp_path / "two.fastq.gz"
    stem = b"@" + b"x" * header_length
    for path, mate in ((first, b"1"), (second, b"2")):
        path.write_bytes(
            gzip.compress(stem + b"/" + mate + b"\nACGT\n+\n!!!!\n", mtime=0)
        )
    assert validate_fastq_records([first, second]) is (header_length + 3 <= 1024 * 1024)


def test_invalid_gzip_missing_files_deadline_and_path_cardinality(
    tmp_path: Path,
) -> None:
    malformed = tmp_path / "bad.fastq.gz"
    malformed.write_bytes(b"not gzip")
    assert not validate_fastq_records([malformed])
    assert not validate_fastq_records([tmp_path / "missing.fastq.gz"])
    assert not validate_fastq_records([])
    assert not validate_fastq_records([malformed] * 3)
    valid = tmp_path / "valid.fastq"
    valid.write_bytes(b"@r\nACGT\n+\n!!!!\n")
    with pytest.raises(StagingDeadline):
        validate_fastq_records([valid], deadline_at=time.monotonic() - 1)


def test_casava_matching_ids_and_labels(tmp_path: Path) -> None:
    first, second = tmp_path / "one.fastq.gz", tmp_path / "two.fastq.gz"
    first.write_bytes(gzip.compress(b"@r 1:N:0:1\nACGT\n+\n!!!!\n", mtime=0))
    second.write_bytes(gzip.compress(b"@r 2:N:0:1\nACGT\n+\n!!!!\n", mtime=0))
    assert validate_fastq_records([first, second])
    second.write_bytes(gzip.compress(b"@r 1:N:0:1\nACGT\n+\n!!!!\n", mtime=0))
    assert not validate_fastq_records([first, second])
