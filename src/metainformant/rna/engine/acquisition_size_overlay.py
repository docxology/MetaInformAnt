"""Hash-bound size supplements for a frozen local acquisition task set."""
from __future__ import annotations

import hashlib
import json
import os
import re
import stat
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path, PurePosixPath
from types import MappingProxyType
from typing import Any
from urllib.parse import unquote, urlsplit

@dataclass(frozen=True)
class SizeOverlay:
    path: Path | None
    sha256: str | None
    sizes: Mapping[str, int]

    def receipt(self) -> dict[str, Any]:
        if self.path is None:
            return {}
        return {"size_overlay_path": str(self.path), "size_overlay_sha256": self.sha256,
                "size_overlay_records": len(self.sizes)}


def effective_fastq_bytes(task: Mapping[str, Any], sizes: Mapping[str, int] | None = None) -> int:
    original = int(task.get("fastq_bytes", 0))
    identity = str(task["task_id"])
    if sizes is None or identity not in sizes:
        return original
    size = sizes[identity]
    if original != 0 or type(size) is not int or not 0 < size < 2**63:
        raise ValueError("Size supplements may only fill a zero weight with positive integer bytes")
    return size


def _object(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for key, value in pairs:
        if key in result:
            raise ValueError("Duplicate size-overlay JSON key")
        result[key] = value
    return result


def _read_bound(path: Path, expected_sha256: str) -> bytes:
    if (not path.is_absolute() or ".." in path.parts
            or any(p.is_symlink() for p in (path, *path.parents))
            or not path.is_file() or path.resolve(strict=True) != path):
        raise ValueError("Size overlay must be an existing absolute nonsymlink regular file")
    before = path.stat()
    if not 0 < before.st_size <= 1024 * 1024:
        raise ValueError("Size overlay exceeds its bounded metadata size")
    descriptor = os.open(path, os.O_RDONLY | os.O_NOFOLLOW)
    with os.fdopen(descriptor, "rb") as handle:
        opened = os.fstat(handle.fileno())
        payload = handle.read(1024 * 1024 + 1)
        after = os.fstat(handle.fileno())
    signature = lambda s: (s.st_dev, s.st_ino, s.st_size, s.st_mtime_ns, s.st_ctime_ns)
    if (not stat.S_ISREG(opened.st_mode) or signature(before) != signature(opened)
            or signature(opened) != signature(after) or signature(after) != signature(path.stat())
            or any(p.is_symlink() for p in (path, *path.parents))
            or hashlib.sha256(payload).hexdigest() != expected_sha256):
        raise ValueError("Size overlay changed or differs from its pinned hash")
    return payload


def load_size_overlay(
    manifest_path: Path, tasks: Sequence[Mapping[str, Any]], *,
    path: str | None, expected_sha256: str | None,
) -> SizeOverlay:
    if not path and not expected_sha256:
        return SizeOverlay(None, None, MappingProxyType({}))
    if not path or not expected_sha256 or re.fullmatch("[0-9a-f]{64}", expected_sha256) is None:
        raise ValueError("Size overlay path and SHA-256 must both be configured")
    location = Path(path)
    document = json.loads(_read_bound(location, expected_sha256), object_pairs_hook=_object)
    manifest_bytes = manifest_path.read_bytes()
    if (not isinstance(document, dict)
            or set(document) != {"schema", "frozen_manifest_sha256", "records"}
            or document["schema"] != "metainformant.rna.size_overlay.v1"
            or document["frozen_manifest_sha256"] != hashlib.sha256(manifest_bytes).hexdigest()):
        raise ValueError("Size overlay does not bind the selected frozen manifest")
    records = document["records"]
    if not isinstance(records, list) or not records or len(records) > len(tasks):
        raise ValueError("Size overlay must contain a bounded nonempty record list")
    by_id = {str(task["task_id"]): task for task in tasks}
    if len(by_id) != len(tasks):
        raise ValueError("Duplicate frozen task identity")
    frozen_tasks = [json.loads(line) for line in manifest_bytes.splitlines() if line.strip()]
    frozen_by_id = {str(task["task_id"]): task for task in frozen_tasks}
    if len(frozen_by_id) != len(frozen_tasks) or frozen_by_id != by_id:
        raise ValueError("Provided tasks differ from the hash-bound frozen task records")
    sizes: dict[str, int] = {}
    for record in records:
        fields = {"task_id", "accession", "expected_paired", "reference_index_sha256", "fastq_bytes", "files"}
        if not isinstance(record, dict) or set(record) != fields:
            raise ValueError("Invalid size-overlay record fields")
        identity = record["task_id"]
        if not isinstance(identity, str) or identity not in by_id or identity in sizes:
            raise ValueError("Unknown or duplicate size-overlay task")
        task = by_id[identity]
        if (record["accession"] != task["accession"]
                or type(record["expected_paired"]) is not bool
                or record["expected_paired"] != task["expected_paired"]
                or record["reference_index_sha256"] != task["reference_index_sha256"]):
            raise ValueError("Size-overlay run, layout or reference identity differs")
        size = effective_fastq_bytes(task, {identity: record["fastq_bytes"]})
        files = record["files"]
        if not isinstance(files, list) or not 1 <= len(files) <= 2:
            raise ValueError("Unsupported size-overlay FASTQ file set")
        names: set[str] = set()
        urls: set[str] = set()
        total = 0
        for entry in files:
            if not isinstance(entry, dict) or set(entry) != {"name", "url", "bytes", "md5"}:
                raise ValueError("Invalid size-overlay file fields")
            if (not isinstance(entry["name"], str) or not isinstance(entry["url"], str)
                    or not isinstance(entry["md5"], str) or type(entry["bytes"]) is not int
                    or not 0 < entry["bytes"] < 2**63
                    or re.fullmatch("[0-9a-f]{32}", entry["md5"]) is None):
                raise ValueError("Invalid size-overlay file bytes or checksum")
            parsed = urlsplit(entry["url"])
            decoded = unquote(parsed.path)
            if (parsed.scheme != "https" or parsed.netloc != "ftp.sra.ebi.ac.uk"
                    or parsed.query or parsed.fragment or not decoded.startswith("/vol1/fastq/")
                    or ".." in PurePosixPath(decoded).parts or PurePosixPath(decoded).name != entry["name"]
                    or entry["name"] in names or entry["url"] in urls):
                raise ValueError("Invalid or duplicate size-overlay file URL")
            names.add(entry["name"]); urls.add(entry["url"]); total += entry["bytes"]
        accession = str(task["accession"])
        valid_names = ({accession + "_1.fastq.gz", accession + "_2.fastq.gz"} if task["expected_paired"] else None)
        if ((valid_names is not None and names != valid_names)
                or (valid_names is None and names not in ({accession + ".fastq.gz"}, {accession + "_1.fastq.gz"}))
                or total != size):
            raise ValueError("Size-overlay mate names or aggregate bytes differ")
        sizes[identity] = size
    return SizeOverlay(location, expected_sha256, MappingProxyType(sizes))


def configured_size_overlay(manifest_path: Path, tasks: Sequence[Mapping[str, Any]]) -> SizeOverlay:
    return load_size_overlay(
        manifest_path, tasks, path=os.environ.get("AMALGKIT_LOCAL_SIZE_OVERLAY_PATH"),
        expected_sha256=os.environ.get("AMALGKIT_LOCAL_SIZE_OVERLAY_SHA256"),
    )
