"""Index digests: aliases share one hash only while the resolved file is unchanged."""

from __future__ import annotations

import os
from pathlib import Path

import pytest

from metainformant.rna.engine.acquisition_manifest import sha256_file
from metainformant.rna.engine.reference_digest import IndexDigestCache


def test_aliases_share_only_unchanged_resolved_index(tmp_path: Path) -> None:
    original = tmp_path / "reference.idx"
    original.write_bytes(b"original reference")
    alias = tmp_path / "alias.idx"
    alias.symlink_to(original.name)
    digests = IndexDigestCache()
    expected = sha256_file(original)
    assert digests.digest(original) == digests.digest(alias) == expected
    assert len(digests.values) == 1

    before = original.stat()
    original.write_bytes(b"modified reference")
    os.utime(original, ns=(before.st_atime_ns, before.st_mtime_ns))
    assert original.stat().st_size == before.st_size
    assert digests.digest(alias) == sha256_file(original) != expected

    replacement = tmp_path / "replacement.idx"
    replacement.write_bytes(b"different reference")
    alias.unlink()
    alias.symlink_to(replacement.name)
    assert digests.digest(alias) == sha256_file(replacement)


def test_separate_files_and_calls_are_hashed_independently(tmp_path: Path) -> None:
    paths = [tmp_path / name for name in ("first.idx", "second.idx")]
    for path in paths:
        path.write_bytes(b"same content")
    digests = IndexDigestCache()
    assert digests.digest(paths[0]) == digests.digest(paths[1])
    assert len(digests.values) == 2
    fresh = IndexDigestCache()
    assert fresh.values == {}
    assert fresh.digest(paths[0]) == sha256_file(paths[0])


def test_disappearing_index_is_never_accepted_from_cache(tmp_path: Path) -> None:
    index = tmp_path / "removed.idx"
    index.write_bytes(b"reference")
    digests = IndexDigestCache()
    digests.digest(index)
    index.unlink()
    with pytest.raises(FileNotFoundError):
        digests.digest(index)
