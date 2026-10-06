"""Append-only content-addressed storage for quantification recovery."""

from __future__ import annotations
import base64
import hashlib
import json
import os
import re
import tempfile
import time
from pathlib import Path, PurePosixPath
from typing import Any, Protocol

ACCESSION = re.compile(r"^(SRR|ERR|DRR)\d+$")
IDENTIFIER = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.-]*$")


def _digest(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def _encode(value: Any) -> bytes:
    return json.dumps(
        value, sort_keys=True, separators=(",", ":"), allow_nan=False
    ).encode()


def safe_key(key: str) -> str:
    path = PurePosixPath(key)
    if (
        not key
        or path.is_absolute()
        or "\\" in key
        or any(p in (".", "..", "") for p in key.split("/"))
    ):
        raise ValueError(f"unsafe object key: {key!r}")
    return key


def receipt_key(cohort: str, species: str, accession: str) -> str:
    if (
        not IDENTIFIER.fullmatch(cohort)
        or not IDENTIFIER.fullmatch(species)
        or not ACCESSION.fullmatch(accession)
    ):
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


class S3Store:
    """S3 conditional writes with server checksum and explicit readback."""

    def __init__(
        self,
        bucket: str,
        prefix: str,
        *,
        profile: str | None = None,
        region: str | None = None,
    ):
        import boto3
        from botocore.config import Config

        self.bucket = bucket
        self.prefix = safe_key(prefix.strip("/"))
        self._published: set[str] = set()
        self.client = boto3.Session(profile_name=profile, region_name=region).client(
            "s3",
            config=Config(
                connect_timeout=10,
                read_timeout=60,
                retries={"mode": "standard", "max_attempts": 5},
            ),
        )

    def _key(self, key: str) -> str:
        return f"{self.prefix}/{safe_key(key)}"

    def get(self, key: str) -> bytes:
        for attempt in range(6):
            try:
                response = self.client.get_object(
                    Bucket=self.bucket, Key=self._key(key)
                )
            except self.client.exceptions.NoSuchKey as exc:
                if key not in self._published or attempt == 5:
                    raise FileNotFoundError(key) from exc
                time.sleep(min(4.0, 0.25 * 2**attempt))
            else:
                break
        with response["Body"] as body:
            return body.read()

    def put(self, key: str, data: bytes) -> None:
        try:
            self.client.put_object(
                Bucket=self.bucket,
                Key=self._key(key),
                Body=data,
                IfNoneMatch="*",
                ChecksumSHA256=base64.b64encode(hashlib.sha256(data).digest()).decode(),
                ServerSideEncryption="AES256",
            )
        except self.client.exceptions.ClientError as exc:
            # S3 does not model conditional-write failures as typed exceptions.
            # Only a failed precondition is an idempotence check; all other errors propagate.
            if exc.response["ResponseMetadata"]["HTTPStatusCode"] != 412:
                raise
            self._published.add(key)
            if self.get(key) != data:
                raise FileExistsError(f"immutable object conflict: {key}") from exc
        self._published.add(key)
        for attempt in range(6):
            try:
                written = self.get(key)
            except FileNotFoundError:
                # Production readback occasionally returned NoSuchKey immediately
                # after a successful PUT. Retry boundedly; never publish on absence.
                if attempt == 5:
                    raise
                time.sleep(min(4.0, 0.25 * 2**attempt))
                continue
            if written != data:
                raise ValueError(f"S3 object readback mismatch: {key}")
            break
