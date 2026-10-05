"""Immutable, content-addressed sample outputs with portable validation.

A receipt is published only after every blob is read back and verified. This
is a quantification recovery boundary, not a manuscript promotion certificate.
"""

from __future__ import annotations

import base64
import csv
import hashlib
import json
import math
import os
import re
import shutil
import tempfile
import time
from pathlib import Path, PurePosixPath
from typing import Any, Protocol

from metainformant.rna.amalgkit import (
    AMALGKIT_RELEASE_TAG,
    AMALGKIT_SOURCE_REVISION,
    REQUIRED_AMALGKIT_VERSION,
)
from metainformant.rna.engine.provenance import (
    QUANT_PROVENANCE_FILENAME,
    quantification_contract_id,
)

SCHEMA = "metainformant.rna.locked_quant.v1"
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


def validate_quantification(
    sample_dir: Path,
    species: str,
    accession: str,
    *,
    expected_config_sha256: str | None = None,
    expected_reference_sha256: str | None = None,
) -> tuple[dict[str, Any], list[Path], int]:
    """Validate portable provenance, unique features and real numeric outputs."""
    receipt_key("validation", species, accession)
    provenance_path = sample_dir / QUANT_PROVENANCE_FILENAME
    provenance = json.loads(provenance_path.read_text())
    if (
        provenance.get("species") != species
        or provenance.get("run_accession") != accession
    ):
        raise ValueError("sample provenance identity mismatch")
    if any(
        provenance.get(k) != v
        for k, v in {
            "amalgkit_version": REQUIRED_AMALGKIT_VERSION,
            "amalgkit_release_tag": AMALGKIT_RELEASE_TAG,
            "amalgkit_source_revision": AMALGKIT_SOURCE_REVISION,
        }.items()
    ):
        raise ValueError("sample runtime differs from current contract")
    contract = quantification_contract_id(provenance)
    if not contract or provenance.get("quant_contract_id") != contract:
        raise ValueError("quantification contract checksum mismatch")
    name = provenance.get("quantification_file")
    if not isinstance(name, str) or Path(name).name != name or name in ("", ".", ".."):
        raise ValueError("unsafe quantification filename")
    abundance = sample_dir / name
    if _digest(abundance.read_bytes()) != provenance.get("quantification_file_sha256"):
        raise ValueError("abundance checksum mismatch")
    for expected, key in [
        (expected_config_sha256, "config_sha256"),
        (expected_reference_sha256, "reference_manifest_sha256"),
    ]:
        if expected is not None and provenance.get(key) != expected:
            raise ValueError(f"{key} differs from expected input")
    features: set[str] = set()
    total_counts = 0.0
    with abundance.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"target_id", "length", "eff_length", "est_counts", "tpm"}
        if not required.issubset(reader.fieldnames or []):
            raise ValueError("abundance table lacks required columns")
        for row in reader:
            target = row["target_id"]
            if not target or target in features:
                raise ValueError("empty or duplicated feature identifier")
            features.add(target)
            for column in required - {"target_id"}:
                value = float(row[column])
                if not math.isfinite(value) or value < 0:
                    raise ValueError(f"invalid abundance {column}")
            total_counts += float(row["est_counts"])
    if not features or total_counts <= 0:
        raise ValueError("empty or zero-count quantification")
    info_path = sample_dir / f"{accession}_run_info.json"
    if not info_path.is_file():
        info_path = sample_dir / "run_info.json"
    info = json.loads(info_path.read_text())
    if int(info.get("n_processed", 0)) <= 0 or int(info.get("n_pseudoaligned", 0)) <= 0:
        raise ValueError("Kallisto run has no processed or aligned reads")
    files = [provenance_path, abundance, info_path]
    files.extend(sorted(p for p in sample_dir.glob("*.h5") if p.is_file()))
    if any(p.is_symlink() for p in files):
        raise ValueError("quantification files must be regular, owned files")
    return provenance, files, len(features)


def lock_quantification(
    store: ObjectStore,
    cohort: str,
    species: str,
    accession: str,
    sample_dir: Path,
    *,
    expected_config_sha256: str | None = None,
    expected_reference_sha256: str | None = None,
) -> dict[str, Any]:
    key = receipt_key(cohort, species, accession)
    provenance, paths, rows = validate_quantification(
        sample_dir,
        species,
        accession,
        expected_config_sha256=expected_config_sha256,
        expected_reference_sha256=expected_reference_sha256,
    )
    files = []
    for path in paths:
        payload = path.read_bytes()
        digest = _digest(payload)
        blob_key = f"blobs/sha256/{digest}"
        store.put(blob_key, payload)
        if _digest(store.get(blob_key)) != digest:
            raise ValueError("stored blob checksum mismatch")
        files.append(
            {"name": path.name, "sha256": digest, "size": len(payload), "key": blob_key}
        )
    receipt = {
        "schema": SCHEMA,
        "cohort": cohort,
        "species": species,
        "accession": accession,
        "contract_id": provenance["quant_contract_id"],
        "config_sha256": provenance["config_sha256"],
        "reference_manifest_sha256": provenance.get("reference_manifest_sha256"),
        "feature_count": rows,
        "files": files,
    }
    store.put(key, _encode(receipt))
    if store.get(key) != _encode(receipt):
        raise ValueError("receipt readback mismatch")
    return receipt


def restore_quantification(
    store: ObjectStore,
    cohort: str,
    species: str,
    accession: str,
    destination: Path,
) -> dict[str, Any]:
    """Restore all blobs to fresh staging; reject corruption before publication."""
    receipt = json.loads(store.get(receipt_key(cohort, species, accession)))
    if receipt.get("schema") != SCHEMA or (
        receipt.get("cohort"),
        receipt.get("species"),
        receipt.get("accession"),
    ) != (cohort, species, accession):
        raise ValueError("stored receipt identity mismatch")
    destination.parent.mkdir(parents=True, exist_ok=True)
    staging = Path(tempfile.mkdtemp(prefix=".restore-", dir=destination.parent))
    try:
        if not isinstance(receipt.get("files"), list) or not receipt["files"]:
            raise ValueError("receipt has no output files")
        for record in receipt["files"]:
            name = record["name"]
            if (
                not isinstance(name, str)
                or Path(name).name != name
                or name in ("", ".", "..")
            ):
                raise ValueError("unsafe stored filename")
            payload = store.get(record["key"])
            if _digest(payload) != record["sha256"] or len(payload) != record["size"]:
                raise ValueError("stored blob checksum mismatch")
            (staging / name).write_bytes(payload)
        validate_quantification(
            staging,
            species,
            accession,
            expected_config_sha256=receipt["config_sha256"],
            expected_reference_sha256=receipt.get("reference_manifest_sha256"),
        )
        if destination.exists():
            for record in receipt["files"]:
                if (
                    not (destination / record["name"]).is_file()
                    or _digest((destination / record["name"]).read_bytes())
                    != record["sha256"]
                ):
                    raise FileExistsError(
                        "restore would overwrite different sample output"
                    )
        else:
            os.rename(staging, destination)
    finally:
        if staging.exists():
            shutil.rmtree(staging)
    return receipt
