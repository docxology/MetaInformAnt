"""Immutable quant receipt publication and verified restoration.

Storage and numerical validation are separate; legacy imports remain public.
"""

from __future__ import annotations
import json
import os
import shutil
import tempfile
from pathlib import Path
from typing import Any
from metainformant.rna.engine.quant_storage import (
    ACCESSION as ACCESSION,
    IDENTIFIER as IDENTIFIER,
    DirectoryStore as DirectoryStore,
    S3Store as S3Store,
    ObjectStore as ObjectStore,
    _digest as _digest,
    _encode as _encode,
    safe_key as safe_key,
    receipt_key as receipt_key,
)
from metainformant.rna.engine.quant_validation import (
    validate_quantification as validate_quantification,
)

SCHEMA = "metainformant.rna.locked_quant.v1"


def bound_receipt_key(cohort: str, species: str, accession: str) -> str:
    """Keep verified index bindings separate from immutable recovery receipts."""
    return receipt_key(cohort, species, accession).replace(
        "/receipts/", "/reference-bound-receipts/"
    )


def lock_quantification(
    store: ObjectStore,
    cohort: str,
    species: str,
    accession: str,
    sample_dir: Path,
    *,
    expected_config_sha256: str | None = None,
    expected_reference_sha256: str | None = None,
    expected_reference_index_sha256: str | None = None,
) -> dict[str, Any]:
    key = (
        bound_receipt_key
        if expected_reference_index_sha256 is not None
        else receipt_key
    )(cohort, species, accession)
    provenance, paths, rows = validate_quantification(
        sample_dir,
        species,
        accession,
        expected_config_sha256=expected_config_sha256,
        expected_reference_sha256=expected_reference_sha256,
    )
    reference_binding = None
    if expected_reference_index_sha256 is not None:
        manifest_path = provenance.get("reference_manifest_path")
        if not manifest_path:
            raise ValueError("quantification lacks a reference index binding")
        manifest_bytes = Path(manifest_path).read_bytes()
        if _digest(manifest_bytes) != provenance.get("reference_manifest_sha256"):
            raise ValueError("reference manifest checksum mismatch")
        manifest = json.loads(manifest_bytes)
        index_path = manifest.get("kallisto_index")
        if (
            not index_path
            or manifest.get("species") != species
            or manifest.get("status") != "complete"
        ):
            raise ValueError("reference manifest lacks a complete species index")
        index_hash = _digest(Path(index_path).read_bytes())
        if index_hash != expected_reference_index_sha256:
            raise ValueError("reference index differs from frozen inventory")
        manifest_key = f"blobs/sha256/{_digest(manifest_bytes)}"
        store.put(manifest_key, manifest_bytes)
        reference_binding = {
            "reference_index_sha256": index_hash,
            "reference_manifest_key": manifest_key,
        }
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
    if reference_binding is not None:
        receipt.update(reference_binding)
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
    *,
    expected_config_sha256: str | None = None,
    expected_reference_index_sha256: str | None = None,
) -> dict[str, Any]:
    """Restore all blobs to fresh staging; reject corruption before publication."""
    key = (
        bound_receipt_key
        if expected_reference_index_sha256 is not None
        else receipt_key
    )(cohort, species, accession)
    receipt = json.loads(store.get(key))
    if receipt.get("schema") != SCHEMA or (
        receipt.get("cohort"),
        receipt.get("species"),
        receipt.get("accession"),
    ) != (cohort, species, accession):
        raise ValueError("stored receipt identity mismatch")
    if (
        expected_config_sha256 is not None
        and receipt.get("config_sha256") != expected_config_sha256
    ):
        raise ValueError("stored configuration differs from frozen inventory")
    if expected_reference_index_sha256 is not None:
        if receipt.get("reference_index_sha256") != expected_reference_index_sha256:
            raise ValueError("stored reference index differs from frozen inventory")
        manifest = store.get(receipt["reference_manifest_key"])
        if _digest(manifest) != receipt.get("reference_manifest_sha256"):
            raise ValueError("stored reference manifest checksum mismatch")
        reference = json.loads(manifest)
        if reference.get("species") != species or reference.get("status") != "complete":
            raise ValueError("stored reference manifest identity mismatch")
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
