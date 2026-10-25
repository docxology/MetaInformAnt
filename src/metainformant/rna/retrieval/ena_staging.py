from __future__ import annotations

import fcntl
import gzip
import hashlib
import json
import os
import shutil
import stat
import tempfile
import time
import zlib
from collections.abc import Callable, Iterator
from contextlib import ExitStack
from io import BufferedIOBase
from itertools import zip_longest
from pathlib import Path
from typing import Protocol

from metainformant.rna.retrieval.ena_downloader import (
    ENADownloader,
    _remote_content_length,
    record_invalid_transfer,
)


class StagingError(RuntimeError):
    def __init__(self, reason: str) -> None:
        self.reason = reason
        super().__init__(reason)


class FirstPassDeferred(RuntimeError):
    def __init__(self, accession: str, reason: str) -> None:
        self.accession = accession
        self.reason = reason
        super().__init__(f"{accession}: first-pass acquisition deferred: {reason}")


class StagingDeadline(StagingError):
    pass


class PromotionLock(Protocol):
    def acquire(self, blocking: bool = True, timeout: float = -1) -> bool: ...
    def release(self) -> None: ...


def _check_deadline(deadline: float | None) -> None:
    if deadline is not None and time.monotonic() >= deadline:
        raise StagingDeadline("SSD acquisition time budget exhausted; sources retained")


def _fastq_record_keys(reader: BufferedIOBase, mate: int | None, deadline: float | None) -> Iterator[bytes | None]:
    carry = b""
    pending: list[bytes] = []
    limit = 1024 * 1024
    while True:
        _check_deadline(deadline)
        chunk = reader.read(limit)
        if not chunk:
            break
        lines = (carry + chunk).split(b"\n")
        carry = lines.pop()
        if len(carry) > limit:
            yield None
            return
        pending.extend(lines)
        end = len(pending) // 4 * 4
        for index in range(0, end, 4):
            if index % 4096 == 0:
                _check_deadline(deadline)
            header, sequence, plus, quality = pending[index : index + 4]
            sequence, quality = (
                sequence.removesuffix(b"\r"),
                quality.removesuffix(b"\r"),
            )
            if (
                max(len(header), len(sequence), len(plus), len(quality)) > limit
                or not header.startswith(b"@")
                or not plus.startswith(b"+")
                or not sequence
                or len(sequence) != len(quality)
            ):
                yield None
                return
            tokens = header[1:].split(None, 1)
            if not tokens:
                yield None
                return
            key = tokens[0]
            if key.endswith((b"/1", b"/2")):
                if mate is not None and not key.endswith(f"/{mate}".encode()):
                    yield None
                    return
                key = key[:-2]
            if not key or (
                mate is not None
                and len(tokens) > 1
                and tokens[1].startswith((b"1:", b"2:"))
                and not tokens[1].startswith(f"{mate}:".encode())
            ):
                yield None
                return
            yield key
        pending = pending[end:]
    if carry or pending:
        yield None


def validate_fastq_records(paths: list[Path], *, deadline_at: float | None = None) -> bool:
    if len(paths) not in (1, 2):
        return False
    count = 0
    try:
        with ExitStack() as stack:
            readers: list[BufferedIOBase] = []
            for path in paths:
                _signature(_plain_path(path))
                if path.name.endswith(".gz"):
                    readers.append(stack.enter_context(gzip.open(path, "rb")))
                else:
                    readers.append(stack.enter_context(path.open("rb")))
            streams = [
                _fastq_record_keys(reader, index if len(paths) == 2 else None, deadline_at)
                for index, reader in enumerate(readers, start=1)
            ]
            for records in zip_longest(*streams):
                if any(key is None for key in records) or len(set(records)) != 1:
                    return False
                count += 1
    except (OSError, EOFError, zlib.error):
        return False
    return count > 0


def _plain_path(path: Path) -> Path:
    absolute = path.absolute()
    if any(part.is_symlink() for part in (absolute, *absolute.parents)):
        raise StagingError(f"staging refuses a symlink path: {path}")
    return absolute


def prepare_sample_directory(path: Path, root: Path, identity: tuple[int, int]) -> None:
    parts = _plain_path(path).relative_to(_plain_path(root)).parts
    if any(part in {".", ".."} for part in parts):
        raise StagingError("sample directory escapes the selected root")
    flags = os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW
    descriptor = os.open(root, flags)
    try:
        root_stat = os.fstat(descriptor)
        if (root_stat.st_dev, root_stat.st_ino) != identity:
            raise StagingError("selected acquisition root changed before directory creation")
        for part in parts:
            try:
                os.mkdir(part, dir_fd=descriptor)
            except FileExistsError:
                pass
            child = os.open(part, flags, dir_fd=descriptor)
            os.close(descriptor)
            descriptor = child
    finally:
        os.close(descriptor)


def _signature(path: Path) -> tuple[int, int, int, int, int]:
    value = path.stat()
    if not stat.S_ISREG(value.st_mode):
        raise StagingError(f"staging requires a regular file: {path}")
    return (
        value.st_dev,
        value.st_ino,
        value.st_size,
        value.st_mtime_ns,
        value.st_ctime_ns,
    )


def _digest(path: Path, *, algorithm: str = "sha256", deadline: float | None = None) -> str:
    digest = hashlib.new(algorithm, usedforsecurity=False)
    with path.open("rb") as handle:
        while chunk := handle.read(1024 * 1024):
            _check_deadline(deadline)
            digest.update(chunk)
    return digest.hexdigest()


def publish_verified_file(
    source: Path,
    destination: Path,
    *,
    deadline: float | None = None,
    expected_signature: tuple[int, int, int, int, int] | None = None,
) -> None:
    source, destination = _plain_path(source), _plain_path(destination)
    before = _signature(source)
    if expected_signature is not None and before != expected_signature:
        raise StagingError("staged source changed after record validation")
    if before[2] == 0:
        raise StagingError(f"staged file is empty: {source}")
    if destination.exists():
        if _signature(destination)[2] != before[2] or _digest(source, deadline=deadline) != _digest(
            destination, deadline=deadline
        ):
            raise StagingError(f"canonical FASTQ conflicts with staged content: {destination}")
        if _signature(source) != before:
            raise StagingError(f"staged source changed: {source}")
        return
    if not destination.parent.is_dir():
        raise StagingError(f"canonical destination directory is unavailable: {destination.parent}")
    with tempfile.NamedTemporaryFile(
        prefix=f".{destination.name}.", suffix=".copying", dir=destination.parent
    ) as temporary:
        digest = hashlib.sha256()
        with source.open("rb") as handle:
            while chunk := handle.read(1024 * 1024):
                _check_deadline(deadline)
                temporary.write(chunk)
                digest.update(chunk)
        temporary.flush()
        os.fsync(temporary.fileno())
        copied = Path(temporary.name)
        if _signature(source) != before or copied.stat().st_size != before[2]:
            raise StagingError(f"staged source changed during copy: {source}")
        if _digest(copied, deadline=deadline) != digest.hexdigest():
            raise StagingError(f"canonical staging-copy checksum mismatch: {destination}")
        os.utime(copied, ns=(source.stat().st_atime_ns, before[3]))
        _check_deadline(deadline)
        os.link(copied, destination)
        directory_fd = os.open(destination.parent, os.O_RDONLY)
        try:
            os.fsync(directory_fd)
        finally:
            os.close(directory_fd)


def download_staged_ena(
    downloader: ENADownloader,
    urls: list[str],
    canonical_dir: Path,
    scratch_dir: Path,
    *,
    reserve: Callable[[int], bool],
    release: Callable[[int], None],
    promotion_lock: PromotionLock,
    validate: Callable[[Path], list[Path]],
    expected_paired: bool | None = None,
    deadline_at: float | None = None,
    minimum_free_bytes: int = 0,
    validation_lock: PromotionLock | None = None,
) -> tuple[bool, str, list[Path]] | None:
    """Stage bounded fresh transfers; None retains the existing bulk resume path.

    Retained SSD parts are keyed by destination, exact URLs and advertised sizes.
    The shared scratch ledger owns the reservation through serialized verified
    promotion. A failed transfer or publication never removes its SSD source.
    """
    canonical = _plain_path(canonical_dir)
    scratch = _plain_path(scratch_dir)
    if not canonical.is_dir():
        raise StagingError("canonical sample directory must already exist")
    destination_identity = canonical.stat().st_dev, canonical.stat().st_ino
    if scratch.is_relative_to(canonical) or canonical.is_relative_to(scratch):
        raise StagingError("canonical and staging roots must be isolated")
    names = [url.rsplit("/", 1)[-1] for url in urls]
    allowed = {f"{canonical.name}{suffix}.fastq.gz" for suffix in ("", "_1", "_2")}
    if not urls or len(set(names)) != len(names) or any(name not in allowed for name in names):
        raise StagingError("staged ENA files must be distinct accession FASTQ names")
    if expected_paired is not False and len(names) == 1 and names[0] == f"{canonical.name}.fastq.gz":
        return None
    deadline = time.monotonic() + downloader.run_timeout_seconds if downloader.run_timeout_seconds else None
    if deadline_at is not None:
        deadline = min(deadline or float("inf"), deadline_at)
    sizes: dict[str, int] = {}
    for name, url in zip(names, urls, strict=True):
        remaining = deadline - time.monotonic() if deadline is not None else 30
        if remaining <= 0:
            return False, "ENA staging discovery exhausted the run time budget", []
        size = downloader.file_sizes.get(name)
        if size is None:
            size = _remote_content_length(url, timeout=min(30, remaining))
        if size is None or size <= 0:
            return None
        sizes[name] = size
    checksums = {name: downloader.file_md5[name] for name in names if name in downloader.file_md5}
    identity = {
        "destination": str(canonical),
        "urls": urls,
        "sizes": sizes,
        "md5": checksums,
        "paired": expected_paired,
    }
    raw_identity = json.dumps(identity, sort_keys=True).encode()
    key = hashlib.sha256(raw_identity).hexdigest()
    workspace = scratch / key
    sample = workspace / "getfastq" / canonical.name
    if (
        canonical.is_dir()
        and any((canonical / name).exists() or (canonical / f"{name}.part").exists() for name in names)
        and not workspace.exists()
    ):
        return None
    if canonical.is_dir() and any((canonical / f"{name}.part").exists() for name in names):
        return None
    reservation = sum(sizes.values()) + len(names) * (2 * downloader.invalid_witness_max_bytes + 1024**2)
    if not reserve(reservation):
        return None
    try:
        _plain_path(sample)
        scratch.mkdir(exist_ok=True)
        workspace.mkdir(exist_ok=True)
        (workspace / "getfastq").mkdir(exist_ok=True)
        sample.mkdir(exist_ok=True)
        lock_path = _plain_path(workspace / ".lock")
        if lock_path.exists():
            _signature(lock_path)
        with lock_path.open("a") as owner:
            fcntl.flock(owner, fcntl.LOCK_EX | fcntl.LOCK_NB)
            binding = _plain_path(workspace / "identity.json")
            if binding.exists() and (
                _signature(binding)[2] != len(raw_identity) or binding.read_bytes() != raw_identity
            ):
                raise StagingError("staging identity changed")
            if not binding.exists():
                with tempfile.NamedTemporaryFile(prefix=".identity-", dir=workspace) as temporary:
                    temporary.write(raw_identity)
                    temporary.flush()
                    os.fsync(temporary.fileno())
                    os.link(temporary.name, binding)
            for name in names:
                for candidate in (sample / name, sample / f"{name}.part"):
                    _plain_path(candidate)
                    if candidate.exists():
                        if _signature(candidate)[2] > sizes[name]:
                            raise StagingError(f"retained SSD input exceeds its bound: {candidate}")
                        if candidate.stat().st_size and name not in checksums:
                            return None
            success, message, _ = downloader._download_files(
                urls,
                sample,
                max_file_bytes=sizes,
                deadline_at=deadline,
                require_existing_output=True,
            )
            if not success:
                return False, message, []
            for name, expected in checksums.items():
                candidate = sample / name
                if _digest(candidate, algorithm="md5", deadline=deadline) != expected:
                    record_invalid_transfer(
                        candidate,
                        reason="ena_md5_mismatch",
                        max_witness_bytes=downloader.invalid_witness_max_bytes,
                        deadline_at=deadline,
                    )
                    return False, "SSD transfer failed ENA content checksum", []
            signatures = {name: _signature(sample / name) for name in names}
            acquired_validation = validation_lock is None or (
                validation_lock.acquire(timeout=max(0, deadline - time.monotonic()))
                if deadline is not None
                else validation_lock.acquire()
            )
            if not acquired_validation:
                raise StagingDeadline("SSD record validation slot time budget exhausted; sources retained")
            try:
                if not validate_fastq_records([sample / name for name in sorted(names)], deadline_at=deadline):
                    return (
                        False,
                        "SSD FASTQs failed full record or paired-read validation",
                        [],
                    )
            finally:
                if validation_lock is not None:
                    validation_lock.release()
            inputs = validate(sample)
            if len(inputs) != len(names) or {path.name for path in inputs} != set(names):
                return False, "SSD FASTQs failed the declared layout contract", []
            if any(path.parent != sample or _signature(_plain_path(path))[2] != sizes[path.name] for path in inputs):
                raise StagingError("validated SSD inputs differ from the bounded transfer set")
            _check_deadline(deadline)
            acquired = (
                promotion_lock.acquire(timeout=max(0, deadline - time.monotonic()))
                if deadline is not None
                else promotion_lock.acquire()
            )
            if not acquired:
                raise StagingDeadline("SSD promotion slot time budget exhausted; sources retained")
            try:
                for path in inputs:
                    if (
                        canonical.stat().st_dev,
                        canonical.stat().st_ino,
                    ) != destination_identity:
                        raise StagingError("canonical directory identity changed during acquisition")
                    missing_bytes = sum(sizes[name] for name in names if not (canonical / name).exists())
                    if shutil.disk_usage(canonical).free < minimum_free_bytes + missing_bytes:
                        raise StagingError("bulk capacity floor prevents SSD publication; sources retained")
                    publish_verified_file(
                        path,
                        canonical / path.name,
                        deadline=deadline,
                        expected_signature=signatures[path.name],
                    )
                for path in inputs:
                    path.unlink()
            finally:
                promotion_lock.release()
            return (
                True,
                "Downloaded and verified SSD-staged FASTQs",
                [canonical / path.name for path in inputs],
            )
    finally:
        release(reservation)
