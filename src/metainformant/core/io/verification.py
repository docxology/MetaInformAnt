"""Physical-file verification snapshots with a cache limited to one operation."""

from __future__ import annotations

import hashlib
import os
import stat
from pathlib import Path
from types import TracebackType
from typing import Literal, Self

FileIdentity = tuple[int, int, int, int, int, int, int, int]


class VerificationChangedError(RuntimeError):
    """A required physical file disappeared, redirected, or changed."""


def _identity(path: Path) -> FileIdentity:
    try:
        value = path.lstat()
        if not path.is_absolute() or not stat.S_ISREG(value.st_mode) or path.resolve(strict=True) != path:
            raise VerificationChangedError(f"Expected an absolute physical regular file: {path}")
        return (
            value.st_dev,
            value.st_ino,
            value.st_mode,
            value.st_uid,
            value.st_gid,
            value.st_size,
            value.st_mtime_ns,
            value.st_ctime_ns,
        )
    except OSError as error:
        raise VerificationChangedError(f"Required file is unavailable: {path}") from error


class BoundReads:
    """Hash each path once while checking its physical identity throughout.

    Reusing a path checks its identity again. Successful context exit rechecks
    every path, including paths read only once. A new instance always rereads
    bytes. This observation does not replace filesystem locks or an independent
    preservation/activation decision.
    """

    def __init__(self) -> None:
        self._files: dict[Path, tuple[FileIdentity, str]] = {}
        self._bytes_hashed = 0

    @property
    def bytes_hashed(self) -> int:
        """Return bytes actually hashed by this instance."""
        return self._bytes_hashed

    def __enter__(self) -> Self:
        return self

    def __call__(self, path: Path) -> str:
        """Return SHA-256, refusing a changed, redirected, or nonregular path.

        Args:
            path: Absolute path with no symlink components.

        Returns:
            The file's full SHA-256 digest.

        Raises:
            VerificationChangedError: File identity or byte-read custody changed.
            OSError: Opening or reading the verified descriptor failed.
        """
        path = Path(path)
        before = _identity(path)
        if path in self._files:
            expected, digest = self._files[path]
            if before != expected:
                raise VerificationChangedError(f"File changed within verification: {path}")
            return digest
        descriptor = os.open(path, os.O_RDONLY | os.O_NOFOLLOW | os.O_CLOEXEC)
        with os.fdopen(descriptor, "rb") as handle:
            opened = os.fstat(handle.fileno())
            stamp = (
                opened.st_dev,
                opened.st_ino,
                opened.st_mode,
                opened.st_uid,
                opened.st_gid,
                opened.st_size,
                opened.st_mtime_ns,
                opened.st_ctime_ns,
            )
            if stamp != before:
                raise VerificationChangedError(f"File changed before read: {path}")
            digest = hashlib.file_digest(handle, "sha256").hexdigest()
            after = os.fstat(handle.fileno())
        finished = (
            after.st_dev,
            after.st_ino,
            after.st_mode,
            after.st_uid,
            after.st_gid,
            after.st_size,
            after.st_mtime_ns,
            after.st_ctime_ns,
        )
        if before != finished or _identity(path) != before:
            raise VerificationChangedError(f"File changed during read: {path}")
        self._files[path] = before, digest
        self._bytes_hashed += before[5]
        return digest

    def finish(self) -> None:
        """Refuse any file change since its last verified byte read."""
        for path, (expected, _) in self._files.items():
            if _identity(path) != expected:
                raise VerificationChangedError(f"File changed after read: {path}")

    def __exit__(
        self,
        kind: type[BaseException] | None,
        error: BaseException | None,
        trace: TracebackType | None,
    ) -> Literal[False]:
        if kind is None:
            self.finish()
        return False
