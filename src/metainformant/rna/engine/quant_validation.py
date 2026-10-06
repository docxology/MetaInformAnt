"""Portable quantification input and output validation."""

from __future__ import annotations
import csv
import json
import math
from pathlib import Path
from typing import Any
from metainformant.rna.amalgkit import (
    AMALGKIT_RELEASE_TAG,
    AMALGKIT_SOURCE_REVISION,
    REQUIRED_AMALGKIT_VERSION,
)
from metainformant.rna.engine.provenance import (
    QUANT_PROVENANCE_FILENAME,
    quantification_contract_id,
)
from metainformant.rna.engine.quant_storage import _digest, receipt_key


class QuantValidationError(ValueError):
    """A quantification fact violates its numeric or structural contract."""

    def __init__(self, field: str, reason: str) -> None:
        self.field = field
        self.reason = reason
        super().__init__(f"{field}: {reason}")


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
    total_tpm = 0.0
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
                    raise QuantValidationError(column, "must be finite and nonnegative")
                if column in {"length", "eff_length"} and value == 0:
                    raise QuantValidationError(column, "must be positive")
            total_counts += float(row["est_counts"])
            total_tpm += float(row["tpm"])
    if not features:
        raise QuantValidationError("target_id", "empty quantification")
    for column, total in [("est_counts", total_counts), ("tpm", total_tpm)]:
        if not math.isfinite(total) or total <= 0:
            raise QuantValidationError(column, "aggregate must be finite and positive")
    info_path = sample_dir / f"{accession}_run_info.json"
    if not info_path.is_file():
        info_path = sample_dir / "run_info.json"
    info = json.loads(info_path.read_text())
    if not isinstance(info, dict):
        raise QuantValidationError("run_info", "must be a JSON object")
    for counter in ("n_processed", "n_pseudoaligned"):
        value = info.get(counter)
        if type(value) is not int or value <= 0:
            raise QuantValidationError(counter, "must be a positive integer")
    if info["n_pseudoaligned"] > info["n_processed"]:
        raise QuantValidationError("n_pseudoaligned", "cannot exceed processed reads")
    if "n_targets" in info:
        value = info["n_targets"]
        if type(value) is not int or value != len(features):
            raise QuantValidationError(
                "n_targets", "must match the abundance feature count"
            )
    files = [provenance_path, abundance, info_path]
    files.extend(sorted(p for p in sample_dir.glob("*.h5") if p.is_file()))
    if any(p.is_symlink() for p in files):
        raise ValueError("quantification files must be regular, owned files")
    return provenance, files, len(features)
