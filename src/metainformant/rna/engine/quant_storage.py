"""Append-only content-addressed storage for quantification recovery."""

from __future__ import annotations

import hashlib
import json
import os
import re
import tempfile
from pathlib import Path, PurePosixPath
from typing import Any, Protocol

ACCESSION = re.compile(r"^(SRR|ERR|DRR)\d+$")
IDENTIFIER = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.-]*$")


def _digest(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def _encode(value: Any) -> bytes:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()


def safe_key(key: str) -> str:
    path = PurePosixPath(key)
    if not key or path.is_absolute() or "\\" in key or any(p in (".", "..", "") for p in key.split("/")):
        raise ValueError(f"unsafe object key: {key!r}")
    return key


def receipt_key(cohort: str, species: str, accession: str) -> str:
    if not IDENTIFIER.fullmatch(cohort) or not IDENTIFIER.fullmatch(species) or not ACCESSION.fullmatch(accession):
        raise ValueError("invalid cohort, species or run accession")
    return f"{cohort}/receipts/{species}/{accession}.json"


class ObjectStore(Protocol):
    def get(self, key: str) -> bytes: ...
    def put(self, key: str, data: bytes) -> None: ...


class DirectoryStore:
    """Append-only local store; publication uses a same-filesystem hard link."""

    def __init__(self, root: Path):
        self.root = root.resolve()
        self.root.mkdir(parents=True, exist_ok=True)

    def _path(self, key: str) -> Path:
        path = self.root / safe_key(key)
        if not path.resolve().is_relative_to(self.root):
            raise ValueError("object path escapes store")
        return path

    def get(self, key: str) -> bytes:
        return self._path(key).read_bytes()

    def put(self, key: str, data: bytes) -> None:
        target = self._path(key)
        target.parent.mkdir(parents=True, exist_ok=True)
        fd, name = tempfile.mkstemp(prefix=".publish-", dir=target.parent)
        temporary = Path(name)
        try:
            with os.fdopen(fd, "wb") as handle:
                handle.write(data)
                handle.flush()
                os.fsync(handle.fileno())
            temporary.chmod(0o444)
            try:
                os.link(temporary, target)
            except FileExistsError:
                if target.read_bytes() != data:
                    raise FileExistsError(f"immutable object conflict: {key}") from None
            directory_fd = os.open(target.parent, os.O_RDONLY)
            try:
                os.fsync(directory_fd)
            finally:
                os.close(directory_fd)
        finally:
            temporary.unlink(missing_ok=True)
