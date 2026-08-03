"""Current-method provenance for per-sample RNA quantification outputs.

Earlier Amalgkit runs can leave readable ``abundance.tsv`` files in the same
external tree as a current campaign. File readability alone cannot establish
that a sample was processed under the pinned project contract, so the
streaming runner writes a small sidecar after successful current quantification
and all release-facing audits can require that sidecar.
"""

from __future__ import annotations

import hashlib
import json
import os
import tempfile
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Mapping, Sequence

from metainformant.rna.amalgkit import AMALGKIT_RELEASE_TAG, AMALGKIT_SOURCE_REVISION, REQUIRED_AMALGKIT_VERSION
from metainformant.rna.core.sample_utils import find_quantification_file

QUANT_PROVENANCE_FILENAME = ".metainformant_quant_provenance.json"
QUANT_PROVENANCE_SCHEMA = "metainformant.rna.quantification.v1"
METADATA_PROVENANCE_FILENAME = ".metainformant_metadata_provenance.json"
METADATA_PROVENANCE_SCHEMA = "metainformant.rna.metadata.v1"
DOWNSTREAM_PROVENANCE_FILENAME = ".metainformant_downstream_provenance.json"
DOWNSTREAM_PROVENANCE_SCHEMA = "metainformant.rna.downstream.v2"
DOWNSTREAM_STEPS = ("merge", "wsfilter", "finalize", "sanity")
DOWNSTREAM_OUTPUT_STAGES = DOWNSTREAM_STEPS


def quant_provenance_path(sample_dir: str | Path) -> Path:
    """Return the current-method provenance sidecar path for a sample."""

    return Path(sample_dir) / QUANT_PROVENANCE_FILENAME


def digest_file(path: Path) -> str | None:
    """Return a SHA-256 digest for a readable file, or ``None``."""

    try:
        digest = hashlib.sha256()
        with path.open("rb") as handle:
            for block in iter(lambda: handle.read(1024 * 1024), b""):
                digest.update(block)
        return digest.hexdigest()
    except OSError:
        return None


def _canonical_digest(value: Any) -> str:
    """Hash a JSON-compatible value with stable serialization."""

    encoded = json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=True)
    return hashlib.sha256(encoded.encode("utf-8")).hexdigest()


def _normalise_steps(steps: Sequence[str]) -> list[str] | None:
    """Return a valid current step list, or ``None``."""

    normalised = [str(step).strip().lower() for step in steps if str(step).strip()]
    if len(normalised) > len(DOWNSTREAM_STEPS):
        return None
    if any(step not in DOWNSTREAM_STEPS for step in normalised):
        return None
    if len(set(normalised)) != len(normalised):
        return None
    return normalised


def _normalise_prefix(steps: Sequence[str]) -> list[str] | None:
    """Return a valid current workflow prefix, or ``None``."""

    normalised = _normalise_steps(steps)
    if normalised is None or normalised != list(DOWNSTREAM_STEPS[: len(normalised)]):
        return None
    return normalised


def _file_state(path: Path, root: Path, *, include_hash: bool) -> dict[str, Any] | None:
    """Return a portable file state record, optionally including its digest."""

    try:
        stat = path.stat()
    except OSError:
        return None
    record: dict[str, Any] = {
        "path": path.relative_to(root).as_posix(),
        "size": stat.st_size,
        "mtime_ns": stat.st_mtime_ns,
    }
    if include_hash:
        digest = digest_file(path)
        if digest is None:
            return None
        record["sha256"] = digest
    return record


def _stage_files(work_path: Path, stage: str) -> list[Path]:
    """List regular files owned by one downstream stage in stable order."""

    stage_path = work_path / stage
    if not stage_path.is_dir():
        return []
    files: list[Path] = []
    try:
        for path in stage_path.rglob("*"):
            if not path.is_file() or path.is_symlink():
                continue
            if path.name.startswith(".") or ".tmp-" in path.name:
                continue
            files.append(path)
    except OSError:
        return []
    return sorted(files)


def downstream_output_manifest(
    work_dir: str | Path,
    steps: Sequence[str],
    *,
    include_hash: bool = True,
) -> dict[str, list[dict[str, Any]]]:
    """Return hashed output records for the requested downstream prefix."""

    work_path = Path(work_dir)
    normalised = _normalise_prefix(steps)
    if normalised is None:
        raise ValueError(f"Downstream steps must be a current prefix: {steps!r}")
    manifest: dict[str, list[dict[str, Any]]] = {}
    for stage in normalised:
        records: list[dict[str, Any]] = []
        for path in _stage_files(work_path, stage):
            record = _file_state(path, work_path, include_hash=include_hash)
            if record is None:
                raise OSError(f"Unable to record downstream output: {path}")
            records.append(record)
        manifest[stage] = records
    return manifest


def _quantification_snapshot(work_dir: Path, *, verify_hashes: bool) -> list[dict[str, Any]] | None:
    """Return the current quantification table snapshot for one work tree."""

    quant_dir = work_dir / "quant"
    if not quant_dir.is_dir():
        return []
    snapshot: list[dict[str, Any]] = []
    try:
        sample_dirs = sorted(path for path in quant_dir.iterdir() if path.is_dir())
    except OSError:
        return None
    for sample_dir in sample_dirs:
        payload = read_quant_provenance(sample_dir)
        if payload is None or not is_current_quantification(
            sample_dir,
            sample_dir.name,
            verify_content=verify_hashes,
        ):
            continue
        quantification_file = payload.get("quantification_file")
        recorded_hash = payload.get("quantification_file_sha256")
        if not isinstance(quantification_file, str) or not isinstance(recorded_hash, str):
            return None
        path = sample_dir / quantification_file
        try:
            stat = path.stat()
        except OSError:
            return None
        if verify_hashes and digest_file(path) != recorded_hash:
            return None
        snapshot.append(
            {
                "run_accession": sample_dir.name,
                "path": path.relative_to(work_dir).as_posix(),
                "size": stat.st_size,
                "mtime_ns": stat.st_mtime_ns,
                "sha256": recorded_hash,
            }
        )
    return snapshot


def _metadata_snapshot(work_path: Path, *, include_hash: bool) -> list[dict[str, Any]]:
    """Return metadata input states, omitting absent optional inputs."""

    records: list[dict[str, Any]] = []
    for name in ("metadata.tsv", "metadata_selected.tsv"):
        path = work_path / "metadata" / name
        if not path.is_file():
            continue
        record = _file_state(path, work_path, include_hash=include_hash)
        if record is None:
            return []
        records.append(record)
    return records


def _snapshot_matches(
    work_path: Path,
    recorded: list[dict[str, Any]],
    *,
    kind: str,
) -> bool:
    """Compare a recorded snapshot cheaply, then verify changed contents."""

    if kind == "quant":
        current = _quantification_snapshot(work_path, verify_hashes=False)
    else:
        current = _metadata_snapshot(work_path, include_hash=False)
    if current is None:
        return False
    if current == recorded:
        return True
    if kind == "quant":
        verified = _quantification_snapshot(work_path, verify_hashes=True)
    else:
        verified = _metadata_snapshot(work_path, include_hash=True)
    return verified == recorded


def _output_manifest_matches(
    work_path: Path,
    recorded: Mapping[str, Sequence[Mapping[str, Any]]],
    steps: Sequence[str],
) -> bool:
    """Compare output metadata first and rehash only changed files."""

    current = downstream_output_manifest(work_path, steps, include_hash=False)
    recorded_plain = {
        stage: [
            {
                "path": item.get("path"),
                "size": item.get("size"),
                "mtime_ns": item.get("mtime_ns"),
            }
            for item in items
        ]
        for stage, items in recorded.items()
    }
    if current == recorded_plain:
        return True
    verified = downstream_output_manifest(work_path, steps, include_hash=True)
    return verified == {stage: list(items) for stage, items in recorded.items()}


def read_quant_provenance(sample_dir: str | Path) -> dict[str, Any] | None:
    """Read a provenance sidecar, returning ``None`` for invalid/missing data."""

    path = quant_provenance_path(sample_dir)
    try:
        payload = json.loads(path.read_text())
    except (OSError, json.JSONDecodeError, TypeError):
        return None
    return payload if isinstance(payload, dict) else None


def is_current_quantification(
    sample_dir: str | Path,
    run_accession: str | None = None,
    *,
    verify_content: bool = True,
) -> bool:
    """Return whether a sample sidecar records the exact current contract."""

    payload = read_quant_provenance(sample_dir)
    if payload is None:
        return False
    required = {
        "schema": QUANT_PROVENANCE_SCHEMA,
        "amalgkit_version": REQUIRED_AMALGKIT_VERSION,
        "amalgkit_release_tag": AMALGKIT_RELEASE_TAG,
        "amalgkit_source_revision": AMALGKIT_SOURCE_REVISION,
    }
    if any(payload.get(key) != value for key, value in required.items()):
        return False
    if run_accession is not None and payload.get("run_accession") != run_accession:
        return False

    recorded_hash = payload.get("quantification_file_sha256")
    recorded_name = payload.get("quantification_file")
    if not isinstance(recorded_hash, str) or not isinstance(recorded_name, str):
        return False
    relative_path = Path(recorded_name)
    if relative_path.is_absolute() or relative_path.parent != Path("."):
        return False
    quantification_file = Path(sample_dir) / relative_path
    if not quantification_file.is_file():
        return False
    return not verify_content or digest_file(quantification_file) == recorded_hash


def metadata_provenance_path(work_dir: str | Path) -> Path:
    """Return the metadata/selection provenance sidecar path."""

    return Path(work_dir) / "metadata" / METADATA_PROVENANCE_FILENAME


def is_current_metadata(work_dir: str | Path) -> bool:
    """Return whether metadata and selection record the exact current contract."""

    work_path = Path(work_dir)
    if not (work_path / "metadata" / "metadata.tsv").is_file():
        return False
    if not (work_path / "metadata" / "metadata_selected.tsv").is_file():
        return False
    try:
        payload = json.loads(metadata_provenance_path(work_path).read_text())
    except (OSError, json.JSONDecodeError, TypeError):
        return False
    if not isinstance(payload, dict) or not all(
        payload.get(key) == value
        for key, value in {
            "schema": METADATA_PROVENANCE_SCHEMA,
            "amalgkit_version": REQUIRED_AMALGKIT_VERSION,
            "amalgkit_release_tag": AMALGKIT_RELEASE_TAG,
            "amalgkit_source_revision": AMALGKIT_SOURCE_REVISION,
        }.items()
    ):
        return False

    recorded_hashes = {
        "metadata_sha256": payload.get("metadata_sha256"),
        "selected_metadata_sha256": payload.get("selected_metadata_sha256"),
    }
    if not all(isinstance(value, str) for value in recorded_hashes.values()):
        return False
    expected_hashes = {
        "metadata_sha256": digest_file(work_path / "metadata" / "metadata.tsv"),
        "selected_metadata_sha256": digest_file(work_path / "metadata" / "metadata_selected.tsv"),
    }
    return recorded_hashes == expected_hashes


def write_metadata_provenance(
    work_dir: str | Path,
    *,
    species: str,
    config_path: str | Path,
    selection_rules_path: str | Path,
) -> Path:
    """Atomically write metadata and selection provenance."""

    work_path = Path(work_dir)
    metadata_dir = work_path / "metadata"
    metadata_dir.mkdir(parents=True, exist_ok=True)

    payload: Mapping[str, Any] = {
        "schema": METADATA_PROVENANCE_SCHEMA,
        "species": species,
        "amalgkit_version": REQUIRED_AMALGKIT_VERSION,
        "amalgkit_release_tag": AMALGKIT_RELEASE_TAG,
        "amalgkit_source_revision": AMALGKIT_SOURCE_REVISION,
        "config_path": str(config_path),
        "config_sha256": digest_file(Path(config_path)),
        "selection_rules_path": str(selection_rules_path),
        "selection_rules_sha256": digest_file(Path(selection_rules_path)),
        # Bind the metadata contract to the bytes that were actually consumed.
        # This catches partial writes or out-of-band edits without forcing old
        # version-bound sidecars to invalidate an otherwise resumable campaign.
        "metadata_sha256": digest_file(metadata_dir / "metadata.tsv"),
        "selected_metadata_sha256": digest_file(metadata_dir / "metadata_selected.tsv"),
        "generated_at": datetime.now(timezone.utc).isoformat(),
    }
    destination = metadata_provenance_path(work_path)
    fd, temporary_name = tempfile.mkstemp(prefix=f".{destination.name}.", dir=str(metadata_dir))
    try:
        with os.fdopen(fd, "w", encoding="utf-8") as handle:
            json.dump(payload, handle, indent=2, sort_keys=True)
            handle.write("\n")
        os.replace(temporary_name, destination)
    except Exception:
        try:
            os.unlink(temporary_name)
        except OSError:
            pass
        raise
    return destination


def downstream_provenance_path(work_dir: str | Path) -> Path:
    """Return the downstream-stage provenance sidecar path."""

    return Path(work_dir) / DOWNSTREAM_PROVENANCE_FILENAME


def read_downstream_provenance(work_dir: str | Path) -> dict[str, Any] | None:
    """Read a downstream sidecar, returning ``None`` if it is invalid."""

    try:
        payload = json.loads(downstream_provenance_path(work_dir).read_text())
    except (OSError, json.JSONDecodeError, TypeError):
        return None
    return payload if isinstance(payload, dict) else None


def _downstream_payload_is_current_base(
    work_path: Path,
    payload: Mapping[str, Any],
    *,
    required_steps: set[str] | None = None,
    expected_species: str | None = None,
    expected_config_path: Path | None = None,
) -> bool:
    """Validate the common downstream sidecar contract."""

    if any(
        payload.get(key) != value
        for key, value in {
            "schema": DOWNSTREAM_PROVENANCE_SCHEMA,
            "amalgkit_version": REQUIRED_AMALGKIT_VERSION,
            "amalgkit_release_tag": AMALGKIT_RELEASE_TAG,
            "amalgkit_source_revision": AMALGKIT_SOURCE_REVISION,
        }.items()
    ):
        return False
    if expected_species is not None and payload.get("species") != expected_species:
        return False
    declared = payload.get("steps")
    if not isinstance(declared, list):
        return False
    normalised = _normalise_steps(declared)
    if normalised is None or not normalised:
        return False
    required = {str(step).strip().lower() for step in (required_steps or set()) if str(step).strip()}
    if not required <= set(normalised):
        return False

    recorded_config_path = payload.get("config_path")
    recorded_config_hash = payload.get("config_sha256")
    if not isinstance(recorded_config_path, str):
        return False
    config_path = expected_config_path or Path(recorded_config_path)
    if config_path.is_file():
        if not isinstance(recorded_config_hash, str):
            return False
        if digest_file(config_path) != recorded_config_hash:
            return False
    elif expected_config_path is not None:
        return False

    if payload.get("quantified_samples", 0) < 0:
        return False
    return True


def current_downstream_steps(work_dir: str | Path) -> set[str]:
    """Return the declared current downstream steps, or an empty set.

    The step list is trusted only after the v2 base contract is valid. Strict
    sidecars additionally bind metadata, quantification, and stage-output
    snapshots; declared test fixtures remain useful for structural report tests
    but are never produced by the campaign entry points.
    """

    work_path = Path(work_dir)
    payload = read_downstream_provenance(work_path)
    if payload is None or not _downstream_payload_is_current_base(work_path, payload):
        return set()
    return set(_normalise_steps(payload["steps"]) or [])


def is_current_downstream(
    work_dir: str | Path,
    required_steps: set[str] | None = None,
    *,
    species: str | None = None,
    config_path: str | Path | None = None,
) -> bool:
    """Return whether downstream evidence matches the current contract.

    The sidecar must include an explicit, usable ``steps`` list.  When
    ``required_steps`` is supplied, every named stage must also be present in
    that list.  This prevents a sidecar for one downstream stage from being
    used as evidence for a different stage, and prevents an incomplete sidecar
    from being treated as generic downstream evidence.
    """

    work_path = Path(work_dir)
    payload = read_downstream_provenance(work_path)
    if payload is None or not _downstream_payload_is_current_base(
        work_path,
        payload,
        required_steps=required_steps,
        expected_species=species,
        expected_config_path=Path(config_path).expanduser().resolve() if config_path else None,
    ):
        return False

    if payload.get("output_contract") != "strict":
        return True

    declared_steps = _normalise_prefix(payload["steps"])
    assert declared_steps is not None
    metadata_inputs = payload.get("metadata_inputs")
    quant_inputs = payload.get("quantification_inputs")
    output_manifest = payload.get("output_manifest")
    if not isinstance(metadata_inputs, list) or not isinstance(quant_inputs, list):
        return False
    if not isinstance(output_manifest, dict):
        return False
    if not _snapshot_matches(work_path, metadata_inputs, kind="metadata"):
        return False
    current_quant = _quantification_snapshot(work_path, verify_hashes=False)
    if current_quant is None:
        return False
    if current_quant != quant_inputs:
        current_quant = _quantification_snapshot(work_path, verify_hashes=True)
        if current_quant != quant_inputs:
            return False
    if int(payload.get("quantified_samples", -1)) != len(quant_inputs):
        return False
    try:
        if not _output_manifest_matches(work_path, output_manifest, declared_steps):
            return False
    except (OSError, ValueError, TypeError):
        return False
    return True


def write_downstream_provenance(
    work_dir: str | Path,
    *,
    species: str,
    config_path: str | Path,
    quantified_samples: int,
    steps: list[str],
    strict: bool = False,
) -> Path:
    """Atomically write a resumable downstream checkpoint.

    ``strict=True`` is the production contract. It records exact metadata,
    current quantification, and stage-output snapshots. ``strict=False`` is
    reserved for lightweight synthetic fixtures that exercise report parsing
    without constructing a full work tree.
    """

    work_path = Path(work_dir)
    work_path.mkdir(parents=True, exist_ok=True)
    normalised_steps = _normalise_steps(steps)
    if normalised_steps is None or not normalised_steps:
        raise ValueError(f"Downstream steps must be a non-empty current step list: {steps!r}")
    if strict and _normalise_prefix(normalised_steps) is None:
        raise ValueError(f"Strict downstream steps must be a current prefix: {steps!r}")
    metadata_inputs = _metadata_snapshot(work_path, include_hash=True) if strict else []
    quant_inputs = _quantification_snapshot(work_path, verify_hashes=True) if strict else []
    if quant_inputs is None:
        raise OSError(f"Unable to record current quantification inputs under {work_path}")
    if strict and len(metadata_inputs) != 2:
        raise OSError(f"Downstream checkpoint requires metadata.tsv and metadata_selected.tsv under {work_path}")
    if strict and not quant_inputs:
        raise OSError(f"Downstream checkpoint requires current quantification inputs under {work_path}")
    if strict and quantified_samples != len(quant_inputs):
        raise ValueError(
            "quantified_samples does not match the current quantification inputs: "
            f"{quantified_samples} != {len(quant_inputs)}"
        )
    output_manifest = downstream_output_manifest(work_path, normalised_steps, include_hash=True) if strict else {}
    if strict and any(not output_manifest.get(step) for step in normalised_steps):
        missing = [step for step in normalised_steps if not output_manifest.get(step)]
        raise OSError("Downstream checkpoint has no readable output for stage(s): " + ", ".join(missing))
    config_file = Path(config_path).expanduser().resolve()
    if strict and digest_file(config_file) is None:
        raise OSError(f"Downstream checkpoint config is missing or unreadable: {config_file}")
    destination = downstream_provenance_path(work_path)
    payload: Mapping[str, Any] = {
        "schema": DOWNSTREAM_PROVENANCE_SCHEMA,
        "species": species,
        "amalgkit_version": REQUIRED_AMALGKIT_VERSION,
        "amalgkit_release_tag": AMALGKIT_RELEASE_TAG,
        "amalgkit_source_revision": AMALGKIT_SOURCE_REVISION,
        "config_path": str(config_file),
        "config_sha256": digest_file(config_file),
        "quantified_samples": quantified_samples,
        "steps": normalised_steps,
        "output_contract": "strict" if strict else "declared",
        "metadata_inputs": metadata_inputs,
        "quantification_inputs": quant_inputs,
        "output_manifest": output_manifest,
        "generated_at": datetime.now(timezone.utc).isoformat(),
    }
    fd, temporary_name = tempfile.mkstemp(prefix=f".{destination.name}.", dir=str(work_path))
    try:
        with os.fdopen(fd, "w", encoding="utf-8") as handle:
            json.dump(payload, handle, indent=2, sort_keys=True)
            handle.write("\n")
        os.replace(temporary_name, destination)
    except Exception:
        try:
            os.unlink(temporary_name)
        except OSError:
            pass
        raise
    return destination


def write_quant_provenance(
    sample_dir: str | Path,
    *,
    species: str,
    run_accession: str,
    config_path: str | Path,
    command: list[str],
    reference_manifest_path: str | Path | None = None,
    quantification_file: str | Path | None = None,
) -> Path:
    """Atomically write the current-method quantification sidecar.

    When the output is omitted, the recognized abundance table in the sample
    directory is resolved deterministically. A current sidecar can never be
    written without binding an actual quantification file.
    """

    sample_path = Path(sample_dir)
    sample_path.mkdir(parents=True, exist_ok=True)
    config_file = Path(config_path)
    try:
        config_sha256 = hashlib.sha256(config_file.read_bytes()).hexdigest()
    except OSError:
        config_sha256 = None
    quant_file = (
        Path(quantification_file)
        if quantification_file is not None
        else find_quantification_file(sample_path, run_accession)
    )
    if quant_file is None:
        raise OSError(f"No recognized quantification output found in {sample_path}")
    if quant_file is not None and not quant_file.is_absolute():
        quant_file = sample_path / quant_file
    quant_file_name: str | None = None
    quant_file_sha256: str | None = None
    if quant_file is not None:
        try:
            relative_quant_file = quant_file.resolve().relative_to(sample_path.resolve())
        except (OSError, ValueError):
            relative_quant_file = Path(quant_file.name)
        if relative_quant_file.parent != Path(".") or not quant_file.is_file():
            raise OSError(f"Quantification file must be a file directly in {sample_path}")
        quant_file_name = relative_quant_file.name
        quant_file_sha256 = digest_file(quant_file)
        if quant_file_sha256 is None:
            raise OSError(f"Unable to hash quantification file: {quant_file}")
    payload: Mapping[str, Any] = {
        "schema": QUANT_PROVENANCE_SCHEMA,
        "species": species,
        "run_accession": run_accession,
        "amalgkit_version": REQUIRED_AMALGKIT_VERSION,
        "amalgkit_release_tag": AMALGKIT_RELEASE_TAG,
        "amalgkit_source_revision": AMALGKIT_SOURCE_REVISION,
        "config_path": str(config_file),
        "config_sha256": config_sha256,
        "reference_manifest_path": str(reference_manifest_path) if reference_manifest_path else None,
        "reference_manifest_sha256": digest_file(Path(reference_manifest_path)) if reference_manifest_path else None,
        "quantification_file": quant_file_name,
        "quantification_file_sha256": quant_file_sha256,
        "command": command,
        "generated_at": datetime.now(timezone.utc).isoformat(),
    }
    destination = quant_provenance_path(sample_path)
    fd, temporary_name = tempfile.mkstemp(prefix=f".{destination.name}.", dir=str(sample_path))
    try:
        with os.fdopen(fd, "w", encoding="utf-8") as handle:
            json.dump(payload, handle, indent=2, sort_keys=True)
            handle.write("\n")
        os.replace(temporary_name, destination)
    except Exception:
        try:
            os.unlink(temporary_name)
        except OSError:
            pass
        raise
    return destination
