"""Real-file custody checks for operation-local hash verification."""

from __future__ import annotations

import hashlib
from contextlib import chdir
from pathlib import Path

import pytest

from metainformant.core.io.verification import BoundReads, VerificationChangedError


@pytest.mark.parametrize("boundary", ["repeat", "finish"])
@pytest.mark.parametrize("mutation", ["append", "replace", "symlink", "mode"])
def test_hash_snapshot_refuses_physical_file_drift(tmp_path: Path, boundary: str, mutation: str) -> None:
    path = tmp_path / "source"
    path.write_bytes(b"preserved source")
    reader = BoundReads()
    expected = reader(path)
    assert reader(path) == expected
    if mutation == "append":
        with path.open("ab") as handle:
            handle.write(b"changed")
    elif mutation == "replace":
        path.rename(tmp_path / "retained")
        path.write_bytes(b"replacement")
    elif mutation == "symlink":
        path.rename(tmp_path / "retained")
        path.symlink_to(tmp_path / "retained")
    else:
        path.chmod(0o400)
    with pytest.raises(VerificationChangedError):
        reader(path) if boundary == "repeat" else reader.finish()


@pytest.mark.parametrize("kind", ["relative", "missing", "directory", "symlink"])
def test_hash_snapshot_requires_physical_absolute_regular_file(tmp_path: Path, kind: str) -> None:
    path = tmp_path / "source"
    if kind == "directory":
        path.mkdir()
    elif kind in {"relative", "symlink"}:
        path.write_bytes(b"source")
    if kind == "symlink":
        link = tmp_path / "link"
        link.symlink_to(path)
        path = link
    with chdir(tmp_path), pytest.raises(VerificationChangedError):
        BoundReads()(Path("source") if kind == "relative" else path)


def test_hash_cache_is_limited_to_one_verification_instance(tmp_path: Path) -> None:
    path = tmp_path / "source"
    path.write_bytes(b"first")
    with BoundReads() as first:
        a = first(path)
        assert first(path) == a
        assert first.bytes_hashed == 5
    path.write_bytes(b"second")
    with BoundReads() as second:
        b = second(path)
    assert a != b and b == hashlib.sha256(b"second").hexdigest()


def test_context_exit_rechecks_paths_read_only_once(tmp_path: Path) -> None:
    path = tmp_path / "source"
    path.write_bytes(b"first")
    with pytest.raises(VerificationChangedError):
        with BoundReads() as reader:
            reader(path)
            path.write_bytes(b"later")


def test_existing_caller_failure_is_not_replaced_by_exit_check(tmp_path: Path) -> None:
    path = tmp_path / "source"
    path.write_bytes(b"first")
    with pytest.raises(ValueError, match="caller failed"):
        with BoundReads() as reader:
            reader(path)
            path.rename(tmp_path / "retained")
            raise ValueError("caller failed")
