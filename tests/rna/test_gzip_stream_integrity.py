from __future__ import annotations

import gzip
from pathlib import Path

import pytest

from metainformant.rna.retrieval.ena_downloader import verify_gzip_integrity


@pytest.mark.parametrize("suffix", [".gz", ".gz.part"])
@pytest.mark.parametrize("damage", ["none", "crc", "truncated", "second-member"])
def test_late_corruption_and_multiple_members(tmp_path: Path, suffix: str, damage: str) -> None:
    payload = gzip.compress(b"@read\nACGT\n+\n!!!!\n" * 200_000, mtime=0)
    if damage == "crc":
        payload = payload[:-8] + bytes([payload[-8] ^ 1]) + payload[-7:]
    elif damage == "truncated":
        payload = payload[:-5]
    elif damage == "second-member":
        payload += gzip.compress(b"second member", mtime=0)[:-4]
    else:
        payload += gzip.compress(b"second member", mtime=0)
    target = tmp_path / f"reads.fastq{suffix}"
    target.write_bytes(payload)
    assert verify_gzip_integrity(target) is (damage == "none")
    assert target.read_bytes() == payload
