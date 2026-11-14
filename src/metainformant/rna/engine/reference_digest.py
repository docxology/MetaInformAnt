"""Content digests of reference files, computed once per file identity.

Kallisto indexes are tens to hundreds of megabytes on a USB data volume. Reuse
classification binds outputs to index *bytes*, not to index paths, so the same
file is hashed many times unless digests are cached. The cache key is the file
identity (device, inode, size, mtime, ctime), never the path: a symlink or a
hard link to an already-hashed file costs one ``stat``, and any in-place edit
changes ``ctime`` and therefore forces a fresh digest.
"""

from __future__ import annotations

import threading
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from metainformant.core.utils.hash import sha256_file

FileSignature = tuple[int, int, int, int, int]


def file_signature(path: Path) -> FileSignature:
    """Return the identity of the file a path currently resolves to."""

    stat = path.stat()
    return (stat.st_dev, stat.st_ino, stat.st_size, stat.st_mtime_ns, stat.st_ctime_ns)


@dataclass(slots=True)
class IndexDigestCache:
    """Thread-safe SHA-256 cache keyed by file identity rather than path.

    ``values`` maps a file signature to its digest. ``memo`` holds derived,
    per-pass results (for example "every index named by this reference manifest
    hashes to X"); callers that need fresh evidence build a new cache.
    """

    values: dict[FileSignature, str] = field(default_factory=dict)
    memo: dict[Any, Any] = field(default_factory=dict)
    _lock: threading.Lock = field(default_factory=threading.Lock, repr=False, compare=False)

    def digest(self, path: Path) -> str:
        """Return the SHA-256 of the file's bytes; raise ``OSError`` if unreadable.

        Raises ``ValueError`` when the file changes while it is being hashed,
        so a digest is never attributed to bytes that no longer exist.
        """

        resolved = Path(path).resolve(strict=True)
        signature = file_signature(resolved)
        with self._lock:
            cached = self.values.get(signature)
        if cached is not None:
            return cached
        digest = sha256_file(resolved)
        if file_signature(resolved) != signature:
            raise ValueError(f"reference file changed while hashing: {resolved}")
        with self._lock:
            self.values[signature] = digest
        return digest

    def digest_or_none(self, path: Path) -> str | None:
        """Return the digest, or ``None`` when the file is missing, unreadable or changing."""

        try:
            return self.digest(path)
        except OSError, ValueError, RuntimeError:
            return None
