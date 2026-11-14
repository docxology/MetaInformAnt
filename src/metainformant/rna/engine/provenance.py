"""Current-method provenance for per-sample RNA quantification outputs.

Earlier Amalgkit runs can leave readable ``abundance.tsv`` files in the same
external tree as a current campaign. File readability alone cannot establish
that a sample was processed under the pinned project contract, so the
streaming runner writes a small sidecar after successful current quantification
and all release-facing audits can require that sidecar.

Reuse is content-addressed. A sidecar binds its output to the sha256 of the
configuration bytes, the sha256 of the reference manifest bytes and, through a
reference-binding witness, restored-input witness or reference manifest, to
the sha256 of the kallisto index bytes. Recorded absolute paths are
informational: moving the repository or the data root, or deleting the
original checkout, never turns a valid output into an invalid one when the
bytes still match. A different Amalgkit runtime is recorded as observed drift
(``version_drift_compatible``), not as a contract failure.
"""

from __future__ import annotations

import errno
import functools
import hashlib
import importlib.metadata
import json
import os
import stat as stat_module
import sys
import tempfile
import threading
from collections.abc import Iterator
from contextlib import contextmanager
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Mapping, Sequence

from metainformant.core.utils.hash import sha256_bytes
from metainformant.rna.amalgkit import AMALGKIT_RELEASE_TAG, AMALGKIT_SOURCE_REVISION, REQUIRED_AMALGKIT_VERSION
from metainformant.rna.core.sample_utils import find_quantification_file
from metainformant.rna.engine.reference_digest import IndexDigestCache, file_signature

QUANT_PROVENANCE_FILENAME = ".metainformant_quant_provenance.json"
QUANT_PROVENANCE_SCHEMA = "metainformant.rna.quantification.v1"
QUANT_CONTRACT_SCHEMA = "metainformant.rna.quantification.contract.v1"
QUANT_STATUS_CURRENT = "current"
QUANT_STATUS_VERSION_DRIFT = "version_drift_compatible"
QUANT_STATUS_LEGACY = "legacy_unverified"
QUANT_STATUS_INVALID = "invalid"
REUSABLE_QUANT_STATUSES = frozenset({QUANT_STATUS_CURRENT, QUANT_STATUS_VERSION_DRIFT})
METADATA_PROVENANCE_FILENAME = ".metainformant_metadata_provenance.json"
METADATA_PROVENANCE_SCHEMA = "metainformant.rna.metadata.v1"
DOWNSTREAM_PROVENANCE_FILENAME = ".metainformant_downstream_provenance.json"
DOWNSTREAM_PROVENANCE_SCHEMA = "metainformant.rna.downstream.v2"
DOWNSTREAM_STEPS = ("merge", "wsfilter", "finalize", "sanity")
DOWNSTREAM_OUTPUT_STAGES = DOWNSTREAM_STEPS
RESTORED_INPUTS_DIR = ".metainformant_restored_inputs"
RESTORED_INPUTS_SCHEMA = "metainformant.rna.restored_inputs.v1"
REFERENCE_BINDING_FILENAME = ".metainformant_reference_binding.json"
REFERENCE_BINDING_SCHEMA = "metainformant.rna.reference_binding.v1"
# Written by the October 2026 local campaign after an index-bound restore.
LOCAL_RESTORE_PROOF_FILENAME = ".local_restore_verified.json"
REFERENCE_MANIFEST_NAME = "reference_aliases.json"
# A real kallisto abundance table is never this small; smaller files are partial writes.
MIN_QUANTIFICATION_BYTES = 101
RUNTIME_FIELDS = ("amalgkit_version", "amalgkit_release_tag", "amalgkit_source_revision")
EVIDENCE_BINDING_WITNESS = "binding_witness"
EVIDENCE_RESTORED_WITNESS = "restored_witness"
EVIDENCE_LOCAL_RESTORE_PROOF = "local_restore_proof"
EVIDENCE_REFERENCE_MANIFEST = "reference_manifest"
# Not a verdict: a file the classification needed could not be read (I/O error, permission).
QUANT_STATUS_UNAVAILABLE = "unavailable"
# errnos that describe the content (nothing there, or not a regular file), never the medium.
_ABSENT_ERRNOS = frozenset({errno.ENOENT, errno.ENOTDIR, errno.ELOOP, errno.EISDIR, errno.ENAMETOOLONG})
_STRICT_READS = threading.local()


class QuantificationReadError(OSError):
    """A file a classification needs exists but could not be read: no verdict about the output is possible."""


def is_unreadable_error(error: BaseException) -> bool:
    """True for an ``OSError`` about the medium or permissions (EIO, EACCES, ...), not about absent content."""

    if isinstance(error, QuantificationReadError):
        return True
    return isinstance(error, OSError) and error.errno is not None and error.errno not in _ABSENT_ERRNOS


@contextmanager
def strict_reads() -> Iterator[None]:
    """In this thread, make classification raise :class:`QuantificationReadError` instead of judging unreadable files.

    By default an unreadable sidecar, witness, manifest or index is folded
    into a verdict (``legacy_unverified``, ``invalid``), which is right for a
    one-off check but wrong for a recorded verdict: a USB hiccup during a
    reconcile pass would quarantine a valid output. Reconcile classifies
    inside this context and records nothing for an unreadable output.
    """

    previous = getattr(_STRICT_READS, "active", False)
    _STRICT_READS.active = True
    try:
        yield
    finally:
        _STRICT_READS.active = previous


def _raise_if_unreadable(error: BaseException, path: Path | str) -> None:
    """Inside :func:`strict_reads`, raise :class:`QuantificationReadError` for an unreadable ``path``.

    Absence is not an error here; I/O and permission errors are.
    """

    if not getattr(_STRICT_READS, "active", False) or not is_unreadable_error(error):
        return
    if isinstance(error, QuantificationReadError):
        raise error
    assert isinstance(error, OSError)
    raise QuantificationReadError(
        error.errno, f"cannot read {Path(path).name}: {error.strerror or error}", str(path)
    ) from error


# ---------------------------------------------------------------------------
# Runtime identity: the pin versus what is actually installed
# ---------------------------------------------------------------------------


def pinned_amalgkit_runtime() -> dict[str, str]:
    """Return the Amalgkit runtime this code is validated against."""

    return {
        "amalgkit_version": REQUIRED_AMALGKIT_VERSION,
        "amalgkit_release_tag": AMALGKIT_RELEASE_TAG,
        "amalgkit_source_revision": AMALGKIT_SOURCE_REVISION,
    }


@functools.lru_cache(maxsize=1)
def _installed_amalgkit() -> tuple[str | None, str | None, str | None]:
    """Return ``(version, commit, requested_revision)`` of the installed distribution."""

    try:
        distribution = importlib.metadata.distribution("amalgkit")
    except importlib.metadata.PackageNotFoundError:
        return None, None, None
    commit: str | None = None
    requested: str | None = None
    try:
        direct_url = json.loads(distribution.read_text("direct_url.json") or "{}")
    except OSError, json.JSONDecodeError, TypeError:
        direct_url = {}
    vcs_info = direct_url.get("vcs_info") if isinstance(direct_url, dict) else None
    if isinstance(vcs_info, dict):
        commit = vcs_info.get("commit_id") if isinstance(vcs_info.get("commit_id"), str) else None
        requested = vcs_info.get("requested_revision") if isinstance(vcs_info.get("requested_revision"), str) else None
    return distribution.version, commit, requested


def observed_amalgkit_runtime() -> dict[str, str | None]:
    """Return the Amalgkit runtime actually installed, as recorded in new sidecars.

    When the installed distribution is exactly the pinned one, the pinned
    fields are returned unchanged, so sidecars stay byte-identical to those
    written before observation was introduced. A different installation is
    recorded as it is; classification then reports ``version_drift_compatible``.
    Only when no distribution metadata is available does this fall back to the
    pin, and ``amalgkit_runtime_source`` says so.
    """

    pinned = pinned_amalgkit_runtime()
    version, commit, requested = _installed_amalgkit()
    if version is None:
        return {**pinned, "amalgkit_runtime_source": "pinned"}
    if version == pinned["amalgkit_version"] and commit == pinned["amalgkit_source_revision"]:
        return {**pinned, "amalgkit_runtime_source": "installed"}
    tag = requested if requested and requested != commit else None
    return {
        "amalgkit_version": version,
        "amalgkit_release_tag": tag,
        "amalgkit_source_revision": commit,
        "amalgkit_runtime_source": "installed",
    }


def quantification_runtime_status(payload: Mapping[str, Any]) -> str:
    """Return ``current`` for the exact pinned runtime, otherwise ``version_drift_compatible``."""

    pinned = pinned_amalgkit_runtime()
    exact = all(payload.get(key) == pinned[key] for key in RUNTIME_FIELDS)
    return QUANT_STATUS_CURRENT if exact else QUANT_STATUS_VERSION_DRIFT


def _observed_kallisto(sample_path: Path, run_accession: str) -> dict[str, Any]:
    """Return the kallisto version and index format recorded by the run that produced the output."""

    for name in (f"{run_accession}_run_info.json", "run_info.json"):
        try:
            info = json.loads((sample_path / name).read_text())
        except OSError, json.JSONDecodeError, TypeError, UnicodeDecodeError:
            continue
        if not isinstance(info, dict):
            continue
        version = info.get("kallisto_version")
        index_version = info.get("index_version")
        return {
            "kallisto_version": version if isinstance(version, str) else None,
            "kallisto_index_version": index_version if type(index_version) is int else None,
        }
    return {"kallisto_version": None, "kallisto_index_version": None}


# ---------------------------------------------------------------------------
# Durable writes (the data volume is a USB disk that has dropped off the bus)
# ---------------------------------------------------------------------------


def _durable_fsync(fd: int) -> None:
    """Flush a file descriptor to stable storage, including the drive cache on macOS."""

    os.fsync(fd)
    if sys.platform == "darwin":
        import fcntl

        full_fsync = getattr(fcntl, "F_FULLFSYNC", None)
        if full_fsync is not None:
            try:
                fcntl.fcntl(fd, full_fsync)
            except OSError:
                pass


def _fsync_directory(directory: Path) -> None:
    """Persist a directory entry change (rename, link) where the platform allows it."""

    try:
        fd = os.open(directory, os.O_RDONLY)
    except OSError:
        return
    try:
        _durable_fsync(fd)
    except OSError:
        pass
    finally:
        os.close(fd)


def _atomic_write_bytes(destination: Path, data: bytes) -> None:
    """Write ``data`` to a temp file, fsync it, replace ``destination`` and fsync the directory."""

    fd, temporary_name = tempfile.mkstemp(prefix=f".{destination.name}.", dir=str(destination.parent))
    try:
        with os.fdopen(fd, "wb") as handle:
            handle.write(data)
            handle.flush()
            _durable_fsync(handle.fileno())
        os.replace(temporary_name, destination)
    except BaseException:
        try:
            os.unlink(temporary_name)
        except OSError:
            pass
        raise
    _fsync_directory(destination.parent)


def _json_document_bytes(payload: Mapping[str, Any]) -> bytes:
    """Serialize a sidecar as indented, key-sorted UTF-8 JSON with a trailing newline."""

    return (json.dumps(payload, indent=2, sort_keys=True) + "\n").encode("utf-8")


def _canonical_json_bytes(payload: Mapping[str, Any]) -> bytes:
    """Serialize a witness canonically: sorted keys, no whitespace, ASCII."""

    return (json.dumps(payload, sort_keys=True, separators=(",", ":"), ensure_ascii=True) + "\n").encode("ascii")


def quantification_input_paths(sample_path: Path, payload: Mapping[str, Any]) -> tuple[Path, Path | None]:
    """Resolve original inputs or explicitly restored, contract-bound copies."""
    restored = sample_path / RESTORED_INPUTS_DIR
    if not restored.exists() and not restored.is_symlink():
        reference = payload.get("reference_manifest_path")
        return Path(payload["config_path"]), Path(reference) if reference else None
    if not restored.resolve().is_relative_to(sample_path.resolve()):
        raise ValueError("restored inputs escape sample directory")
    manifest = json.loads((restored / "manifest.json").read_text())
    expected = {
        "schema": "metainformant.rna.restored_inputs.v1",
        "contract_id": quantification_contract_id(payload),
        "config_sha256": payload.get("config_sha256"),
        "reference_manifest_sha256": payload.get("reference_manifest_sha256"),
    }
    if manifest != expected:
        raise ValueError("restored input witness differs from quantification contract")
    config = restored / "config.yaml"
    reference = restored / "reference_manifest.json" if payload.get("reference_manifest_sha256") else None
    for path in (config, reference, restored / "manifest.json"):
        if path is not None and not path.resolve().is_relative_to(restored.resolve()):
            raise ValueError("restored input file escapes witness directory")
    return config, reference


def quant_provenance_path(sample_dir: str | Path) -> Path:
    """Return the current-method provenance sidecar path for a sample."""

    return Path(sample_dir) / QUANT_PROVENANCE_FILENAME


def digest_file(path: Path) -> str | None:
    """Return a SHA-256 digest for a readable file, or ``None``.

    Delegates to the canonical :func:`metainformant.core.utils.hash.sha256_file`;
    this wrapper returns ``None`` for an unreadable file.
    """

    from metainformant.core.utils.hash import sha256_file

    try:
        return sha256_file(path)
    except OSError as error:
        _raise_if_unreadable(error, path)
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
        classification = classify_quantification(
            sample_dir,
            sample_dir.name,
            verify_content=verify_hashes,
        )
        if payload is None or classification["status"] not in {
            QUANT_STATUS_CURRENT,
            QUANT_STATUS_VERSION_DRIFT,
        }:
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
                "quantification_status": classification["status"],
                "amalgkit_version": classification["observed_amalgkit_version"],
                "amalgkit_release_tag": classification["observed_release_tag"],
                "amalgkit_source_revision": classification["observed_source_revision"],
                "quant_contract_id": classification["contract_id"],
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
    except OSError, json.JSONDecodeError, TypeError:
        return None
    return payload if isinstance(payload, dict) else None


def quantification_contract_payload(payload: Mapping[str, Any]) -> dict[str, Any] | None:
    """Return runtime-independent quantification contract fields."""

    required = (
        "schema",
        "species",
        "run_accession",
        "config_sha256",
        "quantification_file",
        "quantification_file_sha256",
        "command",
    )
    if any(key not in payload for key in required):
        return None
    if payload.get("schema") != QUANT_PROVENANCE_SCHEMA:
        return None
    if not isinstance(payload.get("command"), list) or not payload["command"]:
        return None
    return {
        "schema": QUANT_CONTRACT_SCHEMA,
        "species": payload["species"],
        "run_accession": payload["run_accession"],
        "config_sha256": payload["config_sha256"],
        "reference_manifest_sha256": payload.get("reference_manifest_sha256"),
        "quantification_file": payload["quantification_file"],
        "quantification_file_sha256": payload["quantification_file_sha256"],
        "command": [str(item) for item in payload["command"]],
        "declared_library_layout": payload.get("declared_library_layout"),
        "effective_library_layout": payload.get("effective_library_layout"),
        "layout_reconciliation": payload.get("layout_reconciliation"),
    }


def quantification_contract_id(payload: Mapping[str, Any]) -> str | None:
    """Return a stable identifier for the non-runtime quantification contract."""

    contract = quantification_contract_payload(payload)
    return _canonical_digest(contract) if contract is not None else None


def _read_sidecar(sample_path: Path) -> tuple[bytes | None, dict[str, Any] | None]:
    """Return the raw sidecar bytes and the parsed object, or ``(None, None)``."""

    try:
        data = quant_provenance_path(sample_path).read_bytes()
        payload = json.loads(data)
    except OSError as error:
        _raise_if_unreadable(error, quant_provenance_path(sample_path))
        return None, None
    except json.JSONDecodeError, TypeError, UnicodeDecodeError:
        return None, None
    return (data, payload) if isinstance(payload, dict) else (None, None)


def _read_json_object(path: Path) -> dict[str, Any] | None:
    """Read a small JSON object from a regular (non-symlink) file, or ``None``."""

    try:
        if not stat_module.S_ISREG(path.lstat().st_mode):
            return None
        payload = json.loads(path.read_bytes())
    except OSError as error:
        _raise_if_unreadable(error, path)
        return None
    except json.JSONDecodeError, TypeError, UnicodeDecodeError:
        return None
    return payload if isinstance(payload, dict) else None


def _present(path: Path) -> bool:
    return path.exists() or path.is_symlink()


def _base_result() -> dict[str, Any]:
    return {
        "status": QUANT_STATUS_LEGACY,
        "reason": "missing or unreadable provenance sidecar",
        "observed_amalgkit_version": None,
        "observed_release_tag": None,
        "observed_source_revision": None,
        "contract_id": None,
        "evidence": None,
    }


def _record_observed(result: dict[str, Any], payload: Mapping[str, Any]) -> None:
    result.update(
        observed_amalgkit_version=payload.get("amalgkit_version"),
        observed_release_tag=payload.get("amalgkit_release_tag"),
        observed_source_revision=payload.get("amalgkit_source_revision"),
        observed_kallisto_version=payload.get("kallisto_version"),
    )


def _reference_manifest_candidates(sample_path: Path, payload: Mapping[str, Any]) -> list[tuple[Path, str]]:
    """Return manifest locations to try, the sample's own data root first.

    ``<species>/work/quant/<run>`` sits next to ``<species>/work/reference``, so
    the relative location finds the manifest after the data root moves. The
    recorded absolute location is the fallback. A candidate is evidence only if
    its bytes hash to the recorded ``reference_manifest_sha256``.
    """

    relative = sample_path.absolute().parent.parent / "reference" / REFERENCE_MANIFEST_NAME
    candidates = [(relative, "relative")]
    recorded = payload.get("reference_manifest_path")
    if isinstance(recorded, str) and recorded and Path(recorded).is_absolute() and Path(recorded) != relative:
        candidates.append((Path(recorded), "recorded"))
    return candidates


def _find_reference_manifest(
    sample_path: Path,
    payload: Mapping[str, Any],
    digests: IndexDigestCache,
) -> tuple[Path, str] | None:
    """Return the first manifest candidate whose bytes match the sidecar, or ``None``."""

    expected = payload.get("reference_manifest_sha256")
    if not isinstance(expected, str):
        return None
    for candidate, mode in _reference_manifest_candidates(sample_path, payload):
        try:
            digest: str | None = digests.digest(candidate)
        except OSError as error:
            _raise_if_unreadable(error, candidate)
            digest = None
        except ValueError, RuntimeError:
            digest = None
        if digest == expected:
            return candidate, mode
    return None


def _manifest_index_paths(manifest: Mapping[str, Any]) -> list[str]:
    """Every kallisto index a reference manifest can hand to ``amalgkit quant``."""

    paths: list[str] = []
    main = manifest.get("kallisto_index")
    if isinstance(main, str) and main:
        paths.append(main)
    aliases = manifest.get("aliases")
    for alias in aliases if isinstance(aliases, list) else []:
        if not isinstance(alias, dict):
            continue
        values = [alias.get("index")] if "index" in alias else alias.get("index_aliases") or []
        paths.extend(value for value in values if isinstance(value, str) and value)
    return list(dict.fromkeys(paths))


def _root_mapping(found: Path, mode: str, recorded: Any) -> tuple[Path, Path] | None:
    """Return ``(old_root, new_root)`` when a manifest was found under a moved data root."""

    if mode != "relative" or not isinstance(recorded, str) or not Path(recorded).is_absolute():
        return None
    old_parts, new_parts = Path(recorded).parts, found.absolute().parts
    common = 0
    while common < min(len(old_parts), len(new_parts)) and old_parts[-1 - common] == new_parts[-1 - common]:
        common += 1
    # ``work/reference/reference_aliases.json`` at least; otherwise the roots are unrelated.
    if common < 3:
        return None
    old_root, new_root = Path(*old_parts[:-common]), Path(*new_parts[:-common])
    return None if old_root == new_root else (old_root, new_root)


def _manifest_index_digests(
    manifest_path: Path,
    mapping: tuple[Path, Path] | None,
    digests: IndexDigestCache,
) -> tuple[frozenset[str], str | None]:
    """Hash every index a manifest names (memoised per manifest and index identity).

    Returns ``(digests, error)``. Index paths under the recorded data root are
    re-rooted to the root the manifest was found in, so a copied or moved data
    root is checked against its own index bytes, never another root's. A
    memoised result is reused only while every index still has the identity
    (device, inode, size, mtime, ctime) it had when it was hashed, so a sample
    costs one ``stat`` per index instead of a re-hash, and a rebuilt index is
    noticed.
    """

    try:
        key = ("manifest_index_digests", file_signature(manifest_path), mapping)
    except OSError as unreadable:
        _raise_if_unreadable(unreadable, manifest_path)
        return frozenset(), f"reference manifest unreadable: {unreadable}"
    outcome: tuple[frozenset[str], str | None]
    cached = digests.memo.get(key)
    if cached is not None:
        cached_outcome, signatures = cached
        try:
            if all(file_signature(path) == signature for path, signature in signatures):
                return (frozenset(cached_outcome[0]), cached_outcome[1])
        except OSError:
            pass
    manifest = _read_json_object(manifest_path)
    signatures_seen: list[tuple[Path, Any]] = []
    if manifest is None:
        outcome = (frozenset(), "reference manifest is not a JSON object")
    else:
        found: set[str] = set()
        error: str | None = None
        index_paths = _manifest_index_paths(manifest)
        if not index_paths:
            error = "reference manifest names no kallisto index"
        for value in index_paths:
            index = Path(value)
            if mapping is not None and index.is_relative_to(mapping[0]):
                index = mapping[1] / index.relative_to(mapping[0])
            try:
                signature = file_signature(index)
                digest = digests.digest(index)
                if file_signature(index) != signature:
                    raise ValueError("index changed while hashing")
            except (OSError, ValueError, RuntimeError) as failure:
                _raise_if_unreadable(failure, index)
                error = f"reference index unavailable: {index.name}"
                break
            signatures_seen.append((index, signature))
            found.add(digest)
        outcome = (frozenset(found), error)
    if outcome[1] is None:
        digests.memo[key] = (outcome, tuple(signatures_seen))
    return outcome


def _binding_witness_error(
    binding: Mapping[str, Any] | None,
    sidecar_bytes: bytes,
    payload: Mapping[str, Any],
    contract_id: str,
) -> str | None:
    """Return why a binding witness is not self-consistent with its sidecar, or ``None``."""

    if binding is None:
        return "unreadable"
    if binding.get("schema") != REFERENCE_BINDING_SCHEMA:
        return "unknown schema"
    checks = {
        "species": payload.get("species"),
        "run_accession": payload.get("run_accession"),
        "quant_contract_id": contract_id,
        "sidecar_sha256": sha256_bytes(sidecar_bytes),
        "config_sha256": payload.get("config_sha256"),
    }
    for key, expected in checks.items():
        if binding.get(key) != expected:
            return f"{key} differs from the quantification sidecar"
    index_hash = binding.get("reference_index_sha256")
    if not isinstance(index_hash, str) or len(index_hash) != 64:
        return "reference_index_sha256 is missing"
    return None


def classify_quantification(
    sample_dir: str | Path,
    run_accession: str | None = None,
    *,
    verify_content: bool = True,
    expected_config_path: str | Path | None = None,
    expected_reference_manifest_path: str | Path | None = None,
    expected_command: Sequence[str] | None = None,
    expected_config_sha256: str | None = None,
    index_digests: IndexDigestCache | None = None,
) -> dict[str, Any]:
    """Classify quantification without treating runtime version drift as failure.

    The configuration is checked by content. With ``expected_config_sha256``
    (the sha256 of the current runtime configuration bytes, computed once per
    species) the recorded ``config_sha256`` must equal it and the recorded
    ``config_path`` is never read. Without it, the recorded path is hashed as
    before, and a missing path is still accepted when a restored-input witness
    or a self-consistent reference-binding witness carries the content binding.

    The reference manifest is located next to the sample (``../../reference``)
    or at its recorded path; either is accepted only when its bytes hash to the
    recorded ``reference_manifest_sha256``.
    """

    sample_path = Path(sample_dir)
    result = _base_result()
    sidecar_bytes, payload = _read_sidecar(sample_path)
    if payload is None or sidecar_bytes is None:
        return result
    _record_observed(result, payload)
    contract = quantification_contract_payload(payload)
    contract_id = quantification_contract_id(payload)
    result["contract_id"] = contract_id
    if contract is None or contract_id is None:
        result["reason"] = "provenance sidecar lacks complete quantification contract"
        return result
    if run_accession is not None and payload.get("run_accession") != run_accession:
        result.update(status=QUANT_STATUS_INVALID, reason="provenance accession mismatch")
        return result

    quant_name = str(contract["quantification_file"])
    quant_path = sample_path / quant_name
    if Path(quant_name).is_absolute() or Path(quant_name).parent != Path("."):
        result.update(status=QUANT_STATUS_INVALID, reason="quantification file escapes sample directory")
        return result
    result["quantification_file"] = quant_name
    if not quant_path.is_file():
        result.update(status=QUANT_STATUS_INVALID, reason="quantification file is missing")
        return result
    if verify_content and digest_file(quant_path) != contract["quantification_file_sha256"]:
        result.update(status=QUANT_STATUS_INVALID, reason="quantification checksum mismatch")
        return result

    config_path = payload.get("config_path")
    config_hash = payload.get("config_sha256")
    if not isinstance(config_path, str) or not isinstance(config_hash, str):
        result["reason"] = "configuration provenance is incomplete"
        return result
    reference_path = payload.get("reference_manifest_path")
    reference_hash = payload.get("reference_manifest_sha256")
    if reference_path and (not isinstance(reference_path, str) or not isinstance(reference_hash, str)):
        result["reason"] = "reference provenance is incomplete"
        return result

    digests = index_digests if index_digests is not None else IndexDigestCache()
    binding_path = sample_path / REFERENCE_BINDING_FILENAME
    binding: dict[str, Any] | None = None
    if _present(binding_path):
        binding = _read_json_object(binding_path)
        error = _binding_witness_error(binding, sidecar_bytes, payload, contract_id)
        if error is not None:
            result.update(status=QUANT_STATUS_INVALID, reason=f"invalid reference binding witness: {error}")
            return result

    restored = sample_path / RESTORED_INPUTS_DIR
    has_witness = _present(restored)
    local_reference: Path | None = None
    if has_witness:
        try:
            local_config, local_reference = quantification_input_paths(sample_path, payload)
        except (OSError, ValueError, TypeError, KeyError, json.JSONDecodeError) as error:
            _raise_if_unreadable(error, sample_path / RESTORED_INPUTS_DIR)
            result.update(status=QUANT_STATUS_INVALID, reason=f"invalid restored input witness: {error}")
            return result
        if digest_file(local_config.expanduser().resolve()) != config_hash:
            result.update(status=QUANT_STATUS_INVALID, reason="configuration checksum mismatch")
            return result
        config_evidence = EVIDENCE_RESTORED_WITNESS
    elif expected_config_sha256 is not None:
        config_evidence = "content"
    elif digest_file(Path(config_path).expanduser()) == config_hash:
        config_evidence = "recorded_path"
    elif binding is not None:
        config_evidence = EVIDENCE_BINDING_WITNESS
    else:
        result.update(status=QUANT_STATUS_INVALID, reason="configuration checksum mismatch")
        return result
    if expected_config_sha256 is not None and config_hash != expected_config_sha256:
        result.update(status=QUANT_STATUS_INVALID, reason="configuration differs from current runtime configuration")
        return result

    reference_evidence: str | None = None
    if reference_hash:
        if has_witness:
            if local_reference is None or digest_file(local_reference.expanduser().resolve()) != reference_hash:
                result.update(status=QUANT_STATUS_INVALID, reason="reference manifest checksum mismatch")
                return result
            reference_evidence = EVIDENCE_RESTORED_WITNESS
        else:
            found = _find_reference_manifest(sample_path, payload, digests)
            if found is not None:
                reference_evidence = f"{EVIDENCE_REFERENCE_MANIFEST}:{found[1]}"
            elif binding is not None:
                reference_evidence = EVIDENCE_BINDING_WITNESS
            else:
                result.update(status=QUANT_STATUS_INVALID, reason="reference manifest checksum mismatch")
                return result
    if expected_config_path is not None and digest_file(Path(expected_config_path)) != config_hash:
        result.update(status=QUANT_STATUS_INVALID, reason="current configuration differs")
        return result
    if expected_reference_manifest_path is not None:
        if not reference_hash or digest_file(Path(expected_reference_manifest_path)) != reference_hash:
            result.update(status=QUANT_STATUS_INVALID, reason="current reference manifest differs")
            return result
    if expected_command is not None and contract["command"] != [str(item) for item in expected_command]:
        result.update(status=QUANT_STATUS_INVALID, reason="quantification command differs")
        return result

    status = quantification_runtime_status(payload)
    result.update(
        status=status,
        reason=(
            "exact current runtime"
            if status == QUANT_STATUS_CURRENT
            else "contract verified; runtime version drift recorded"
        ),
        evidence={"config": config_evidence, "reference": reference_evidence},
    )
    return result


def is_reusable_quantification(
    sample_dir: str | Path,
    run_accession: str | None = None,
    *,
    verify_content: bool = True,
    expected_config_path: str | Path | None = None,
    expected_reference_manifest_path: str | Path | None = None,
    expected_command: Sequence[str] | None = None,
    expected_config_sha256: str | None = None,
    index_digests: IndexDigestCache | None = None,
) -> bool:
    """Return whether quantification is safe to reuse without re-quantifying."""

    status = classify_quantification(
        sample_dir,
        run_accession,
        verify_content=verify_content,
        expected_config_path=expected_config_path,
        expected_reference_manifest_path=expected_reference_manifest_path,
        expected_command=expected_command,
        expected_config_sha256=expected_config_sha256,
        index_digests=index_digests,
    )["status"]
    return status in REUSABLE_QUANT_STATUSES


def is_current_quantification(
    sample_dir: str | Path,
    run_accession: str | None = None,
    *,
    verify_content: bool = True,
) -> bool:
    """Return whether a sample sidecar records the exact current runtime."""

    status: str = classify_quantification(sample_dir, run_accession, verify_content=verify_content)["status"]
    return status == QUANT_STATUS_CURRENT


# ---------------------------------------------------------------------------
# Cohort (frozen-envelope) reuse contract v2
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class CohortExpectation:
    """What one frozen cohort task requires of a reusable output.

    Every value comes from the frozen cohort definition, never from the sample
    directory being judged.
    """

    cohort_id: str
    species: str
    run_accession: str
    config_sha256: str
    reference_index_sha256: str
    reference_index_name: str


def _run_info_error(sample_path: Path, run_accession: str) -> str | None:
    """Return why kallisto run counters are unusable, or ``None`` when they are consistent."""

    for name in (f"{run_accession}_run_info.json", "run_info.json"):
        info = _read_json_object(sample_path / name)
        if info is None:
            continue
        processed, aligned = info.get("n_processed"), info.get("n_pseudoaligned")
        if type(processed) is not int or processed <= 0:
            return "run_info n_processed is not a positive integer"
        if type(aligned) is not int or aligned <= 0:
            return "run_info n_pseudoaligned is not a positive integer"
        if aligned > processed:
            return "run_info n_pseudoaligned exceeds n_processed"
        return None
    return "run_info is missing or unreadable"


def _cohort_reference_evidence(
    sample_path: Path,
    payload: Mapping[str, Any],
    sidecar_bytes: bytes,
    contract_id: str,
    expectation: CohortExpectation,
    digests: IndexDigestCache,
) -> tuple[dict[str, Any] | None, str | None]:
    """Return ``(evidence, None)``, ``(None, error)`` or ``(None, None)`` when nothing binds the index.

    Evidence is checked in a fixed order. A witness that is present but does
    not match is an error, never a fall-through, so tampering cannot be hidden
    behind weaker evidence.
    """

    binding_path = sample_path / REFERENCE_BINDING_FILENAME
    if _present(binding_path):
        binding = _read_json_object(binding_path)
        error = _binding_witness_error(binding, sidecar_bytes, payload, contract_id)
        if error is not None or binding is None:
            return None, f"invalid reference binding witness: {error}"
        for key, expected in {
            "cohort_id": expectation.cohort_id,
            "config_sha256": expectation.config_sha256,
            "reference_index_sha256": expectation.reference_index_sha256,
        }.items():
            if binding.get(key) != expected:
                return None, f"reference binding witness {key} differs from the cohort"
        return {"kind": EVIDENCE_BINDING_WITNESS, "binding": binding.get("evidence")}, None

    if _present(sample_path / RESTORED_INPUTS_DIR):
        try:
            local_config, local_reference = quantification_input_paths(sample_path, payload)
        except (OSError, ValueError, TypeError, KeyError, json.JSONDecodeError) as error:
            _raise_if_unreadable(error, sample_path / RESTORED_INPUTS_DIR)
            return None, f"invalid restored input witness: {error}"
        if digest_file(local_config.resolve()) != expectation.config_sha256:
            return None, "restored configuration differs from the cohort configuration"
        reference_hash = payload.get("reference_manifest_sha256")
        if not isinstance(reference_hash, str) or local_reference is None:
            return None, "restored input witness lacks a reference manifest"
        manifest_bytes_ok = digest_file(local_reference.resolve()) == reference_hash
        manifest = _read_json_object(local_reference) if manifest_bytes_ok else None
        if manifest is None:
            return None, "restored reference manifest checksum mismatch"
        index = manifest.get("kallisto_index")
        if manifest.get("species") != expectation.species or manifest.get("status") != "complete":
            return None, "restored reference manifest is not a complete manifest for this species"
        if not isinstance(index, str) or Path(index).name != expectation.reference_index_name:
            return None, "restored reference manifest names a different index"
        return {"kind": EVIDENCE_RESTORED_WITNESS}, None

    proof_path = sample_path / LOCAL_RESTORE_PROOF_FILENAME
    if _present(proof_path):
        proof = _read_json_object(proof_path)
        if proof is None or proof.get("contract_id") != contract_id:
            return None, "local restore proof names a different quantification contract"
        if proof.get("reference_index_sha256") != expectation.reference_index_sha256:
            return None, "local restore proof names a different reference index"
        return {"kind": EVIDENCE_LOCAL_RESTORE_PROOF}, None

    found = _find_reference_manifest(sample_path, payload, digests)
    if found is None:
        return None, None
    manifest_path, mode = found
    mapping = _root_mapping(manifest_path, mode, payload.get("reference_manifest_path"))
    manifest = _read_json_object(manifest_path)
    if manifest is None or manifest.get("species") != expectation.species or manifest.get("status") != "complete":
        return None, "reference manifest is not a complete manifest for this species"
    index_digests, index_error = _manifest_index_digests(manifest_path, mapping, digests)
    if index_error is not None:
        return None, index_error
    if index_digests != {expectation.reference_index_sha256}:
        return None, "reference index bytes differ from the cohort index"
    return {
        "kind": EVIDENCE_REFERENCE_MANIFEST,
        "manifest": mode,
        "index_paths": "rerooted" if mapping is not None else "recorded",
    }, None


def classify_cohort_quantification(
    sample_dir: str | Path,
    expectation: CohortExpectation,
    *,
    verify_content: bool = False,
    index_digests: IndexDigestCache | None = None,
) -> dict[str, Any]:
    """Classify one output against a frozen cohort task (quant reuse contract v2).

    The decision uses only content: the sidecar bytes and fields, the
    configuration sha256, and a reference-index binding (binding witness,
    restored-input witness, local restore proof, or a reference manifest whose
    index bytes hash to the cohort index). Recorded paths and the command are
    informational. Resume uses ``verify_content=False`` (stat only); import and
    deep audits pass ``True`` to re-hash the abundance table and check the
    kallisto run counters. Pass one ``index_digests`` cache per species so each
    index and manifest is hashed once per pass.
    """

    sample_path = Path(sample_dir)
    result = _base_result()
    result["sidecar_sha256"] = None
    if sample_path.name != expectation.run_accession:
        result.update(status=QUANT_STATUS_INVALID, reason="sample directory name differs from run accession")
        return result
    sidecar_bytes, payload = _read_sidecar(sample_path)
    if payload is None or sidecar_bytes is None:
        return result
    _record_observed(result, payload)
    result["sidecar_sha256"] = sha256_bytes(sidecar_bytes)
    contract = quantification_contract_payload(payload)
    contract_id = quantification_contract_id(payload)
    result["contract_id"] = contract_id
    if contract is None or contract_id is None:
        result["reason"] = "provenance sidecar lacks complete quantification contract"
        return result
    recorded_contract = payload.get("quant_contract_id")
    if recorded_contract is not None and recorded_contract != contract_id:
        result.update(status=QUANT_STATUS_INVALID, reason="quantification contract checksum mismatch")
        return result
    if payload.get("species") != expectation.species or payload.get("run_accession") != expectation.run_accession:
        result.update(status=QUANT_STATUS_INVALID, reason="provenance identity differs from the cohort task")
        return result

    quant_name = contract["quantification_file"]
    if not isinstance(quant_name, str) or Path(quant_name).name != quant_name or quant_name in ("", ".", ".."):
        result.update(status=QUANT_STATUS_INVALID, reason="quantification file escapes sample directory")
        return result
    quant_path = sample_path / quant_name
    result["quantification_file"] = quant_name
    try:
        quant_stat = quant_path.lstat()
    except OSError as error:
        _raise_if_unreadable(error, quant_path)
        result.update(status=QUANT_STATUS_INVALID, reason="quantification file is missing")
        return result
    if not stat_module.S_ISREG(quant_stat.st_mode):
        result.update(status=QUANT_STATUS_INVALID, reason="quantification file is not a regular file")
        return result
    if quant_stat.st_size < MIN_QUANTIFICATION_BYTES:
        result.update(status=QUANT_STATUS_INVALID, reason="quantification file is truncated")
        return result
    if verify_content:
        if digest_file(quant_path) != contract["quantification_file_sha256"]:
            result.update(status=QUANT_STATUS_INVALID, reason="quantification checksum mismatch")
            return result
        run_info_error = _run_info_error(sample_path, expectation.run_accession)
        if run_info_error is not None:
            result.update(status=QUANT_STATUS_INVALID, reason=run_info_error)
            return result

    if payload.get("config_sha256") != expectation.config_sha256:
        result.update(status=QUANT_STATUS_INVALID, reason="configuration differs from the cohort configuration")
        return result

    digests = index_digests if index_digests is not None else IndexDigestCache()
    evidence, evidence_error = _cohort_reference_evidence(
        sample_path, payload, sidecar_bytes, contract_id, expectation, digests
    )
    if evidence_error is not None:
        result.update(status=QUANT_STATUS_INVALID, reason=evidence_error)
        return result
    if evidence is None:
        result.update(status=QUANT_STATUS_INVALID, reason="no reference-index binding evidence")
        return result

    status = quantification_runtime_status(payload)
    result.update(
        status=status,
        reason=(
            "exact current runtime"
            if status == QUANT_STATUS_CURRENT
            else "contract verified; runtime version drift recorded"
        ),
        evidence=evidence,
    )
    return result


def _reject_absolute_paths(value: Any, where: str = "evidence") -> None:
    """Refuse absolute paths inside witness evidence; witnesses must survive a root move."""

    if isinstance(value, str) and (value.startswith("/") or (len(value) > 2 and value[1:3] in (":\\", ":/"))):
        raise ValueError(f"{where} must not contain absolute paths: {value!r}")
    if isinstance(value, Mapping):
        for key, item in value.items():
            _reject_absolute_paths(item, f"{where}.{key}")
    elif isinstance(value, (list, tuple)):
        for item in value:
            _reject_absolute_paths(item, where)


def reference_binding_payload(
    sidecar_bytes: bytes,
    sidecar: Mapping[str, Any],
    expectation: CohortExpectation,
    evidence: Mapping[str, Any],
) -> dict[str, Any]:
    """Build a reference-binding witness for a sidecar that matches its cohort task.

    ``sidecar`` must be the object parsed from ``sidecar_bytes``. ``evidence``
    says how the index binding was established (for example ``{"kind":
    "produced"}`` or ``{"kind": "receipt", "receipt_cohort": ...}``); it may
    carry labels and relative paths only.
    """

    if json.loads(sidecar_bytes) != dict(sidecar):
        raise ValueError("sidecar object does not match sidecar bytes")
    if sidecar.get("species") != expectation.species or sidecar.get("run_accession") != expectation.run_accession:
        raise ValueError("sidecar identity differs from the cohort task")
    contract_id = quantification_contract_id(sidecar)
    if contract_id is None:
        raise ValueError("sidecar lacks a complete quantification contract")
    if sidecar.get("quant_contract_id") not in (None, contract_id):
        raise ValueError("sidecar quantification contract checksum mismatch")
    if sidecar.get("config_sha256") != expectation.config_sha256:
        raise ValueError("sidecar configuration differs from the cohort configuration")
    if not isinstance(evidence.get("kind"), str) or not evidence["kind"]:
        raise ValueError("binding evidence needs a kind")
    _reject_absolute_paths(dict(evidence))
    return {
        "schema": REFERENCE_BINDING_SCHEMA,
        "cohort_id": expectation.cohort_id,
        "species": expectation.species,
        "run_accession": expectation.run_accession,
        "quant_contract_id": contract_id,
        "sidecar_sha256": sha256_bytes(sidecar_bytes),
        "config_sha256": expectation.config_sha256,
        "reference_index_sha256": expectation.reference_index_sha256,
        "reference_index_name": expectation.reference_index_name,
        "evidence": dict(evidence),
    }


def write_reference_binding(sample_dir: str | Path, payload: Mapping[str, Any]) -> Path:
    """Create the binding witness once; identical rewrites are no-ops, different ones raise.

    The witness is written to a temp file, fsynced, then hard-linked into
    place, so an existing witness is never replaced and a crash leaves either
    no witness or a complete one.
    """

    if payload.get("schema") != REFERENCE_BINDING_SCHEMA:
        raise ValueError("not a reference binding payload")
    sample_path = Path(sample_dir)
    destination = sample_path / REFERENCE_BINDING_FILENAME
    data = _canonical_json_bytes(payload)

    def existing_matches() -> bool:
        try:
            return stat_module.S_ISREG(destination.lstat().st_mode) and destination.read_bytes() == data
        except OSError:
            return False

    if _present(destination):
        if existing_matches():
            return destination
        raise FileExistsError(f"a different reference binding witness exists: {destination}")
    fd, temporary_name = tempfile.mkstemp(prefix=f".{destination.name}.", dir=str(sample_path))
    try:
        with os.fdopen(fd, "wb") as handle:
            handle.write(data)
            handle.flush()
            _durable_fsync(handle.fileno())
        try:
            os.link(temporary_name, destination)
        except FileExistsError:
            if not existing_matches():
                raise FileExistsError(f"a different reference binding witness exists: {destination}") from None
    finally:
        try:
            os.unlink(temporary_name)
        except OSError:
            pass
    _fsync_directory(sample_path)
    return destination


def read_reference_binding(sample_dir: str | Path) -> dict[str, Any] | None:
    """Read a binding witness, returning ``None`` if it is missing or unreadable."""

    binding = _read_json_object(Path(sample_dir) / REFERENCE_BINDING_FILENAME)
    return binding if binding is not None and binding.get("schema") == REFERENCE_BINDING_SCHEMA else None


@dataclass(frozen=True)
class SpeciesReuseContext:
    """Per-species expectations, computed once and shared by every sample classification.

    Hashing the configuration and kallisto index once per species keeps a
    reconcile pass O(samples): each sample then costs a sidecar read, a few
    ``stat`` calls and, for witnesses, a few small reads. With ``cohort_id``,
    ``reference_index_sha256`` and ``reference_index_name`` set, samples are
    judged by :func:`classify_cohort_quantification`; without them, by
    :func:`classify_quantification` with the configuration bound by content.
    """

    species: str
    config_sha256: str
    cohort_id: str | None = None
    reference_index_sha256: str | None = None
    reference_index_name: str | None = None
    index_digests: IndexDigestCache = field(default_factory=IndexDigestCache, compare=False, repr=False)

    @property
    def is_cohort(self) -> bool:
        return bool(self.cohort_id and self.reference_index_sha256 and self.reference_index_name)

    def expectation(self, run_accession: str) -> CohortExpectation:
        if not self.is_cohort:
            raise ValueError("cohort expectations need cohort_id, reference_index_sha256 and reference_index_name")
        return CohortExpectation(
            cohort_id=str(self.cohort_id),
            species=self.species,
            run_accession=run_accession,
            config_sha256=self.config_sha256,
            reference_index_sha256=str(self.reference_index_sha256),
            reference_index_name=str(self.reference_index_name),
        )

    def classify(
        self,
        sample_dir: str | Path,
        run_accession: str | None = None,
        *,
        verify_content: bool = False,
    ) -> dict[str, Any]:
        """Classify one sample directory against this species' expectations."""

        sample_path = Path(sample_dir)
        accession = run_accession or sample_path.name
        if self.is_cohort:
            return classify_cohort_quantification(
                sample_path,
                self.expectation(accession),
                verify_content=verify_content,
                index_digests=self.index_digests,
            )
        return classify_quantification(
            sample_path,
            accession,
            verify_content=verify_content,
            expected_config_sha256=self.config_sha256,
            index_digests=self.index_digests,
        )


def build_species_reuse_context(
    species: str,
    config_path: str | Path,
    *,
    cohort_id: str | None = None,
    reference_index_sha256: str | None = None,
    reference_index_name: str | None = None,
    reference_index_path: str | Path | None = None,
    index_digests: IndexDigestCache | None = None,
) -> SpeciesReuseContext:
    """Hash the runtime configuration (and optionally one index) once for a species.

    ``reference_index_path`` is hashed and must equal ``reference_index_sha256``
    when both are given; given alone, it supplies the expected index digest and
    name. Raises ``OSError`` for an unreadable configuration and ``ValueError``
    for an index that differs from the expectation.
    """

    digests = index_digests if index_digests is not None else IndexDigestCache()
    config_sha256 = sha256_bytes(Path(config_path).read_bytes())
    if reference_index_path is not None:
        index_path = Path(reference_index_path)
        actual = digests.digest(index_path)
        if reference_index_sha256 is not None and actual != reference_index_sha256:
            raise ValueError(f"reference index differs from the expected digest: {index_path.name}")
        reference_index_sha256 = actual
        reference_index_name = reference_index_name or index_path.name
    return SpeciesReuseContext(
        species=species,
        config_sha256=config_sha256,
        cohort_id=cohort_id,
        reference_index_sha256=reference_index_sha256,
        reference_index_name=reference_index_name,
        index_digests=digests,
    )


def metadata_provenance_path(work_dir: str | Path) -> Path:
    """Return the metadata/selection provenance sidecar path."""

    return Path(work_dir) / "metadata" / METADATA_PROVENANCE_FILENAME


def read_metadata_provenance(work_dir: str | Path) -> dict[str, Any] | None:
    """Read the metadata/selection sidecar, returning ``None`` for invalid or missing data."""

    try:
        payload = json.loads(metadata_provenance_path(work_dir).read_text())
    except OSError, json.JSONDecodeError, TypeError, UnicodeDecodeError:
        return None
    return payload if isinstance(payload, dict) else None


def is_current_metadata(
    work_dir: str | Path,
    *,
    expected_config_sha256: str | None = None,
    expected_selection_rules_sha256: str | None = None,
) -> bool:
    """Return whether metadata and selection are bound to the bytes on disk.

    Currency is a content contract: the sidecar schema, and the sha256 of
    ``metadata.tsv`` and ``metadata_selected.tsv`` as recorded when they were
    produced. The Amalgkit runtime that produced them is recorded, not
    compared, so bumping the Amalgkit pin never forces a fresh NCBI
    ``metadata``/``select`` query (which could change the cohort); see
    :func:`metadata_runtime_status`. Callers that know the current
    configuration or selection-rule bytes can also require those hashes.
    """

    work_path = Path(work_dir)
    if not (work_path / "metadata" / "metadata.tsv").is_file():
        return False
    if not (work_path / "metadata" / "metadata_selected.tsv").is_file():
        return False
    payload = read_metadata_provenance(work_path)
    if payload is None or payload.get("schema") != METADATA_PROVENANCE_SCHEMA:
        return False
    for expected, key in (
        (expected_config_sha256, "config_sha256"),
        (expected_selection_rules_sha256, "selection_rules_sha256"),
    ):
        if expected is not None and payload.get(key) != expected:
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


def metadata_runtime_status(work_dir: str | Path) -> str | None:
    """Return ``current`` or ``version_drift_compatible`` for current metadata, else ``None``."""

    if not is_current_metadata(work_dir):
        return None
    payload = read_metadata_provenance(work_dir) or {}
    return quantification_runtime_status(payload)


def write_metadata_provenance(
    work_dir: str | Path,
    *,
    species: str,
    config_path: str | Path,
    selection_rules_path: str | Path,
) -> Path:
    """Atomically and durably write metadata and selection provenance.

    The runtime fields name the Amalgkit that produced the tables. When the
    existing sidecar already binds the same table bytes, its runtime is kept,
    so a later runtime does not claim tables it did not produce. Otherwise the
    observed (installed) runtime is recorded.
    """

    work_path = Path(work_dir)
    metadata_dir = work_path / "metadata"
    metadata_dir.mkdir(parents=True, exist_ok=True)
    metadata_sha256 = digest_file(metadata_dir / "metadata.tsv")
    selected_sha256 = digest_file(metadata_dir / "metadata_selected.tsv")
    previous = read_metadata_provenance(work_path)
    if (
        previous is not None
        and previous.get("schema") == METADATA_PROVENANCE_SCHEMA
        and metadata_sha256 is not None
        and previous.get("metadata_sha256") == metadata_sha256
        and previous.get("selected_metadata_sha256") == selected_sha256
        and all(isinstance(previous.get(key), (str, type(None))) for key in RUNTIME_FIELDS)
        and isinstance(previous.get("amalgkit_version"), str)
    ):
        runtime = {key: previous.get(key) for key in RUNTIME_FIELDS}
    else:
        observed = observed_amalgkit_runtime()
        runtime = {key: observed.get(key) for key in RUNTIME_FIELDS}

    # Content-deterministic payload: no wall-clock fields.  Recency lives in
    # the orchestrator log and the progress DB, not in hashed/compared sidecars.
    payload: Mapping[str, Any] = {
        "schema": METADATA_PROVENANCE_SCHEMA,
        "species": species,
        **runtime,
        "config_path": str(config_path),
        "config_sha256": digest_file(Path(config_path)),
        "selection_rules_path": str(selection_rules_path),
        "selection_rules_sha256": digest_file(Path(selection_rules_path)),
        # Bind the metadata contract to the bytes that were actually consumed.
        # This catches partial writes or out-of-band edits without forcing old
        # version annotations to invalidate an otherwise resumable campaign.
        "metadata_sha256": metadata_sha256,
        "selected_metadata_sha256": selected_sha256,
    }
    destination = metadata_provenance_path(work_path)
    _atomic_write_bytes(destination, _json_document_bytes(payload))
    return destination


def downstream_provenance_path(work_dir: str | Path) -> Path:
    """Return the downstream-stage provenance sidecar path."""

    return Path(work_dir) / DOWNSTREAM_PROVENANCE_FILENAME


def read_downstream_provenance(work_dir: str | Path) -> dict[str, Any] | None:
    """Read a downstream sidecar, returning ``None`` if it is invalid."""

    try:
        payload = json.loads(downstream_provenance_path(work_dir).read_text())
    except OSError, json.JSONDecodeError, TypeError:
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
    except OSError, ValueError, TypeError:
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
    config_hash = digest_file(config_file)
    if strict and config_hash is None:
        raise OSError(f"Downstream checkpoint config is missing or unreadable: {config_file}")
    if strict and any(
        (read_quant_provenance(work_path / "quant" / sample["run_accession"]) or {}).get("config_sha256") != config_hash
        for sample in quant_inputs
    ):
        raise ValueError("Downstream configuration differs from quantified sample inputs")
    destination = downstream_provenance_path(work_path)
    # Content-deterministic payload: no wall-clock fields.  Any restart-varying
    # byte would make the sidecar look rewritten on every resume.  Recency
    # lives in the orchestrator log and the progress DB, not here.
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
    }
    _atomic_write_bytes(destination, _json_document_bytes(payload))
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
    """Atomically and durably write the current-method quantification sidecar.

    When the output is omitted, the recognized abundance table in the sample
    directory is resolved deterministically. A current sidecar can never be
    written without binding an actual quantification file.

    The sidecar records the Amalgkit runtime actually installed and the
    kallisto version from the run's ``run_info.json``, plus portable fields
    (``config_name`` and ``reference_manifest_relpath``) next to the
    absolute paths. None of these enter the quantification
    contract, so contract ids are unchanged. The quantification outputs in
    the sample directory (abundance table, ``.h5``, run info) are flushed to
    stable storage first, so a crash can never leave a sidecar that vouches
    for bytes still in the page cache. The sidecar is fsynced before it
    replaces the previous one, and the directory entry is fsynced after.
    """

    sample_path = Path(sample_dir)
    sample_path.mkdir(parents=True, exist_ok=True)
    config_file = Path(config_path).expanduser().resolve()
    reference_file = Path(reference_manifest_path).expanduser().resolve() if reference_manifest_path else None
    try:
        config_sha256 = sha256_bytes(config_file.read_bytes())
    except OSError:
        config_sha256 = None
    quant_file = (
        Path(quantification_file)
        if quantification_file is not None
        else find_quantification_file(sample_path, run_accession)
    )
    if quant_file is None:
        raise OSError(f"No recognized quantification output found in {sample_path}")
    if not quant_file.is_absolute():
        quant_file = sample_path / quant_file
    quant_file_name: str | None = None
    quant_file_sha256: str | None = None
    try:
        relative_quant_file = quant_file.resolve().relative_to(sample_path.resolve())
    except OSError, ValueError:
        relative_quant_file = Path(quant_file.name)
    if relative_quant_file.parent != Path(".") or not quant_file.is_file():
        raise OSError(f"Quantification file must be a file directly in {sample_path}")
    quant_file_name = relative_quant_file.name
    quant_file_sha256 = digest_file(quant_file)
    if quant_file_sha256 is None:
        raise OSError(f"Unable to hash quantification file: {quant_file}")
    reference_relpath: str | None = None
    if reference_file is not None:
        resolved_sample = sample_path.resolve()
        if reference_file.is_relative_to(resolved_sample.parent.parent):
            reference_relpath = Path(os.path.relpath(reference_file, resolved_sample)).as_posix()
    # Content-deterministic payload: no wall-clock fields.  This sidecar is
    # re-verified fail-closed on every resume; any restart-varying byte would
    # invalidate reusable quantification work.  Recency lives in the
    # orchestrator log and the progress DB, not here.
    payload: dict[str, Any] = {
        "schema": QUANT_PROVENANCE_SCHEMA,
        "species": species,
        "run_accession": run_accession,
        **observed_amalgkit_runtime(),
        **_observed_kallisto(sample_path, run_accession),
        "config_path": str(config_file),
        "config_name": config_file.name,
        "config_sha256": config_sha256,
        "reference_manifest_path": str(reference_file) if reference_file else None,
        "reference_manifest_relpath": reference_relpath,
        "reference_manifest_sha256": digest_file(reference_file) if reference_file else None,
        "quantification_file": quant_file_name,
        "quantification_file_sha256": quant_file_sha256,
        "quant_contract_schema": QUANT_CONTRACT_SCHEMA,
        "command": command,
    }
    payload["quant_contract_id"] = quantification_contract_id(payload)
    sync_sample_outputs(sample_path)
    destination = quant_provenance_path(sample_path)
    _atomic_write_bytes(destination, _json_document_bytes(payload))
    return destination


def sync_sample_outputs(sample_dir: str | Path) -> int:
    """Flush every regular file directly in ``sample_dir`` (not dot files) and the directory; return the count.

    Raises ``OSError`` when a file cannot be flushed: a sidecar must not be
    published over outputs that may not be on the disk.
    """

    sample_path = Path(sample_dir)
    synced = 0
    with os.scandir(sample_path) as entries:
        names = sorted(
            entry.name for entry in entries if not entry.name.startswith(".") and entry.is_file(follow_symlinks=False)
        )
    for name in names:
        fd = os.open(sample_path / name, os.O_RDONLY)
        try:
            _durable_fsync(fd)
        finally:
            os.close(fd)
        synced += 1
    _fsync_directory(sample_path)
    return synced


def bind_produced_output(
    sample_dir: str | Path,
    expectation: CohortExpectation,
    *,
    index_digests: IndexDigestCache | None = None,
) -> Path:
    """Write the reference-binding witness for an output this producer just quantified under a cohort task.

    The output is first judged with content verification against the
    cohort task; the witness then records the sidecar's sha256 and the cohort
    index sha256, so later reuse does not depend on the mutable
    ``reference_aliases.json`` staying byte-identical. Raises ``ValueError``
    when the output does not verify (no witness is written).
    """

    sample_path = Path(sample_dir)
    verdict = classify_cohort_quantification(sample_path, expectation, verify_content=True, index_digests=index_digests)
    if verdict["status"] not in REUSABLE_QUANT_STATUSES:
        raise ValueError(f"produced output does not verify against cohort {expectation.cohort_id}: {verdict['reason']}")
    sidecar_bytes = quant_provenance_path(sample_path).read_bytes()
    evidence = verdict.get("evidence") or {}
    payload = reference_binding_payload(
        sidecar_bytes,
        json.loads(sidecar_bytes),
        expectation,
        {"kind": "produced", "verified_by": str(evidence.get("kind") or "unknown")},
    )
    return write_reference_binding(sample_path, payload)


def discard_reference_binding(sample_dir: str | Path) -> bool:
    """Remove the binding witness of an output that is about to be replaced; True when one was removed.

    A witness names one sidecar's sha256: once quantification rewrites the
    sidecar, the old witness would make the new output invalid.
    """

    destination = Path(sample_dir) / REFERENCE_BINDING_FILENAME
    try:
        destination.unlink()
    except FileNotFoundError:
        return False
    _fsync_directory(Path(sample_dir))
    return True
