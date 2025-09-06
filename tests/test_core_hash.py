from __future__ import annotations

from pathlib import Path

from metainformant.core import hash as core_hash


def test_content_hash_and_file_hash(tmp_path: Path) -> None:
    content = b"hello world\n"
    h1 = core_hash.sha256_bytes(content)
    path = tmp_path / "file.txt"
    path.write_bytes(content)
    h2 = core_hash.sha256_file(path)
    assert h1 == h2 and len(h1) == 64
