"""First-class evidence bundles for release-grade RNA pipeline artifacts.

An evidence bundle is an immutable, single-file ledger (JSON) that binds
together every piece of evidence a campaign produced: immutable input and
configuration hashes, software/tool versions, lock and heartbeat history,
terminal and unresolved task counts, per-step receipts, artifact manifests,
result matrices, the durable error taxonomy, and command transcripts.

Design invariants carried over from the existing evidence vocabulary:

- Hashes are SHA-256; content digests use stable canonical JSON serialization
  (same convention as :mod:`metainformant.rna.engine.provenance`).
- Bundles are written atomically: staged in a temporary file inside the
  destination directory, then moved into place with ``os.replace``.
- All evidence artifacts live under exactly one data root. Evidence recorded
  from more than one root is rejected (fail-closed) rather than merged.
- Roles are distinct states, never conflated: acquisition, recovery, and
  quantification are operational evidence; descriptive analysis and
  biological inference are separate artifact classes. The public API offers
  one registration method per role, so a caller cannot register one artifact
  under two roles or invent a sixth role.

The schema is ``metainformant.rna.evidence_bundle.v1``.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import tempfile
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from metainformant.rna.engine.progress_db import classify_sample_error

EVIDENCE_BUNDLE_SCHEMA = "metainformant.rna.evidence_bundle.v1"
EVIDENCE_BUNDLE_FILENAME = "evidence_bundle.json"

# ---------- Artifact roles (distinct states, never conflated) ----------

ROLE_ACQUISITION = "acquisition"
ROLE_RECOVERY = "recovery"
ROLE_QUANTIFICATION = "quantification"
ROLE_DESCRIPTIVE_ANALYSIS = "descriptive_analysis"
ROLE_BIOLOGICAL_INFERENCE = "biological_inference"

ALL_ROLES = frozenset(
    {
        ROLE_ACQUISITION,
        ROLE_RECOVERY,
        ROLE_QUANTIFICATION,
        ROLE_DESCRIPTIVE_ANALYSIS,
        ROLE_BIOLOGICAL_INFERENCE,
    }
)
"""Operational evidence records how data was obtained and quantified."""

OPERATIONAL_ROLES = frozenset({ROLE_ACQUISITION, ROLE_RECOVERY, ROLE_QUANTIFICATION})
DESCRIPTIVE_ROLES = frozenset({ROLE_DESCRIPTIVE_ANALYSIS})
INFERENTIAL_ROLES = frozenset({ROLE_BIOLOGICAL_INFERENCE})

ARTIFACT_KINDS = frozenset({"receipt", "manifest", "matrix"})

# ---------- Task-count state partition (mirrors ProgressDB states) ----------

TERMINAL_TASK_STATES = frozenset({"quantified", "failed"})
UNRESOLVED_TASK_STATES = frozenset({"pending", "downloading", "downloaded", "quantifying", "quarantined"})

# ---------- Exceptions ----------


class EvidenceBundleError(Exception):
    """Base class for evidence-bundle failures."""


class BundleSchemaError(EvidenceBundleError):
    """The payload is not a recognisable evidence bundle."""


class MixedRootError(EvidenceBundleError):
    """Evidence was recorded from more than one data root."""


class RoleConflictError(EvidenceBundleError):
    """One artifact was bound to two distinct roles."""


class StaleEvidenceError(EvidenceBundleError):
    """A recorded digest no longer matches the bytes on disk."""


# ---------- Hashing helpers (same conventions as provenance.py) ----------


def digest_file(path: str | Path) -> str:
    """Return the SHA-256 digest of a file's bytes."""

    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _canonical_digest(value: Any) -> str:
    """Hash a JSON-compatible value with stable serialization."""

    encoded = json.dumps(value, sort_keys=True, separators=(",", ":"), default=str)
    return hashlib.sha256(encoded.encode("utf-8")).hexdigest()


# ---------- Builder ----------


class EvidenceBundleBuilder:
    """Accumulate evidence and emit an immutable bundle payload.

    Every artifact-bearing registration binds the artifact's resolved path to
    exactly one role. Re-binding the same path under a different role raises
    :class:`RoleConflictError`; registering a path outside the bundle's data
    root raises :class:`MixedRootError`.
    """

    def __init__(self, *, data_root: str | Path) -> None:
        self._root = Path(data_root).expanduser().resolve()
        self._inputs: list[dict[str, Any]] = []
        self._software: list[dict[str, Any]] = []
        self._locks: list[dict[str, Any]] = []
        self._task_counts: dict[str, dict[str, int]] = {}
        self._artifacts: list[dict[str, Any]] = []
        self._errors: list[dict[str, Any]] = []
        self._commands: list[dict[str, Any]] = []
        self._role_index: dict[Path, str] = {}

    # -- internal plumbing --

    def _resolve_in_root(self, path: str | Path) -> Path:
        """Resolve an artifact path and refuse anything outside the data root."""

        resolved = Path(path).expanduser().resolve()
        if not resolved.is_relative_to(self._root):
            raise MixedRootError(
                f"artifact {resolved} is outside the bundle data root {self._root}; "
                "evidence from a second root is never merged"
            )
        return resolved

    def _bind(
        self,
        *,
        role: str,
        kind: str,
        path: str | Path,
        extra: Mapping[str, Any] | None = None,
    ) -> dict[str, Any]:
        """Hash one artifact and bind it to a single role."""

        if role not in ALL_ROLES:  # defensive; public methods fix the role
            raise ValueError(f"unknown artifact role: {role!r}")
        if kind not in ARTIFACT_KINDS:
            raise ValueError(f"unknown artifact kind: {kind!r}")
        resolved = self._resolve_in_root(path)
        existing_role = self._role_index.get(resolved)
        if existing_role is not None and existing_role != role:
            raise RoleConflictError(
                f"artifact {resolved} is already bound to role {existing_role!r}; "
                f"it cannot also be bound to {role!r}"
            )
        if existing_role is None:
            self._role_index[resolved] = role
        record: dict[str, Any] = {
            "path": resolved.as_posix(),
            "role": role,
            "kind": kind,
            "bytes": resolved.stat().st_size,
            "sha256": digest_file(resolved),
        }
        if extra:
            record.update(extra)
        self._artifacts.append(record)
        return record

    # -- immutable inputs and software contract --

    def add_input(self, *, name: str, path: str | Path) -> dict[str, Any]:
        """Hash one immutable input (config, metadata, rules) by reference.

        Inputs are external sources consumed by the campaign, so unlike
        evidence artifacts they may live outside the data root; they are
        recorded by hash, not merged as evidence.
        """

        resolved = Path(path).expanduser().resolve()
        record = {
            "name": name,
            "path": resolved.as_posix(),
            "bytes": resolved.stat().st_size,
            "sha256": digest_file(resolved),
        }
        self._inputs.append(record)
        return record

    def add_software(self, *, name: str, version: str, source_revision: str | None = None) -> dict[str, Any]:
        """Record one software/tool version in the bundle's software contract."""

        record: dict[str, Any] = {"name": name, "version": version}
        if source_revision is not None:
            record["source_revision"] = source_revision
        self._software.append(record)
        return record

    # -- lock and heartbeat history --

    def record_lock(
        self,
        *,
        name: str,
        holder: str,
        acquired_at_utc: str,
        released_at_utc: str | None = None,
        heartbeats: Sequence[Mapping[str, Any]] = (),
    ) -> dict[str, Any]:
        """Record one lock acquisition with its heartbeat history."""

        record: dict[str, Any] = {
            "name": name,
            "holder": holder,
            "acquired_at_utc": acquired_at_utc,
            "released_at_utc": released_at_utc,
            "heartbeats": [dict(heartbeat) for heartbeat in heartbeats],
        }
        self._locks.append(record)
        return record

    # -- terminal and unresolved task counts --

    def record_task_counts(
        self,
        *,
        terminal: Mapping[str, int] | None = None,
        unresolved: Mapping[str, int] | None = None,
    ) -> dict[str, dict[str, int]]:
        """Record terminal and unresolved task counts.

        A state may appear in exactly one partition; counting a terminal
        state as unresolved (or vice versa) is rejected at the API boundary.
        """

        terminal = dict(terminal or {})
        unresolved = dict(unresolved or {})
        overlap = sorted(set(terminal) & set(unresolved))
        if overlap:
            raise EvidenceBundleError(
                "terminal and unresolved task states must be disjoint; " f"both counted: {', '.join(overlap)}"
            )
        unknown = sorted((set(terminal) | set(unresolved)) - TERMINAL_TASK_STATES - UNRESOLVED_TASK_STATES)
        if unknown:
            raise EvidenceBundleError(f"unknown task states: {', '.join(unknown)}")
        self._task_counts = {"terminal": terminal, "unresolved": unresolved}
        return self._task_counts

    # -- per-role artifact registration (one method per role) --

    def add_acquisition_receipt(self, *, step: str, path: str | Path) -> dict[str, Any]:
        """Bind a transfer/download receipt as acquisition evidence."""

        return self._bind(role=ROLE_ACQUISITION, kind="receipt", path=path, extra={"step": step})

    def add_recovery_receipt(self, *, step: str, path: str | Path) -> dict[str, Any]:
        """Bind a retry/quarantine audit as recovery evidence."""

        return self._bind(role=ROLE_RECOVERY, kind="receipt", path=path, extra={"step": step})

    def add_quantification_receipt(self, *, step: str, path: str | Path) -> dict[str, Any]:
        """Bind a quantification-sidecar receipt as quantification evidence."""

        return self._bind(role=ROLE_QUANTIFICATION, kind="receipt", path=path, extra={"step": step})

    def add_descriptive_manifest(self, *, path: str | Path) -> dict[str, Any]:
        """Bind a manifest of generated inputs as descriptive-analysis evidence."""

        return self._bind(role=ROLE_DESCRIPTIVE_ANALYSIS, kind="manifest", path=path)

    def add_descriptive_matrix(self, *, path: str | Path) -> dict[str, Any]:
        """Bind a descriptive result table as descriptive-analysis evidence."""

        return self._bind(role=ROLE_DESCRIPTIVE_ANALYSIS, kind="matrix", path=path)

    def add_inference_manifest(self, *, path: str | Path) -> dict[str, Any]:
        """Bind an inference-input manifest as biological-inference evidence."""

        return self._bind(role=ROLE_BIOLOGICAL_INFERENCE, kind="manifest", path=path)

    def add_inference_matrix(self, *, path: str | Path) -> dict[str, Any]:
        """Bind an inference result table as biological-inference evidence."""

        return self._bind(role=ROLE_BIOLOGICAL_INFERENCE, kind="matrix", path=path)

    # -- error taxonomy and command transcripts --

    def record_error(
        self,
        message: str,
        *,
        species: str | None = None,
        sample: str | None = None,
    ) -> dict[str, Any]:
        """Record one failure message classified by the durable error taxonomy."""

        record: dict[str, Any] = {
            "message": message,
            "error_class": classify_sample_error(message),
        }
        if species is not None:
            record["species"] = species
        if sample is not None:
            record["sample"] = sample
        self._errors.append(record)
        return record

    def record_command(
        self,
        argv: Sequence[str],
        *,
        exit_code: int = 0,
        started_at_utc: str | None = None,
        finished_at_utc: str | None = None,
    ) -> dict[str, Any]:
        """Record one executed command transcript entry."""

        record: dict[str, Any] = {
            "argv": [str(part) for part in argv],
            "exit_code": int(exit_code),
        }
        if started_at_utc is not None:
            record["started_at_utc"] = started_at_utc
        if finished_at_utc is not None:
            record["finished_at_utc"] = finished_at_utc
        self._commands.append(record)
        return record

    # -- emission --

    def payload(self) -> dict[str, Any]:
        """Return the complete bundle payload including its content digest."""

        body: dict[str, Any] = {
            "schema": EVIDENCE_BUNDLE_SCHEMA,
            "data_root": self._root.as_posix(),
            "inputs": list(self._inputs),
            "software": list(self._software),
            "locks": list(self._locks),
            "task_counts": {
                "terminal": dict(self._task_counts.get("terminal", {})),
                "unresolved": dict(self._task_counts.get("unresolved", {})),
            },
            "artifacts": list(self._artifacts),
            "errors": list(self._errors),
            "commands": list(self._commands),
        }
        payload = dict(body)
        payload["bundle_id"] = _canonical_digest(body)
        return payload

    def build(self, destination: str | Path) -> Path:
        """Atomically write the bundle: stage a temporary file, then rename."""

        destination = Path(destination).expanduser().resolve()
        destination.parent.mkdir(parents=True, exist_ok=True)
        rendered = json.dumps(self.payload(), indent=2, sort_keys=True) + "\n"
        fd, temporary_name = tempfile.mkstemp(
            prefix=f".{destination.name}.", suffix=".tmp", dir=str(destination.parent)
        )
        try:
            with os.fdopen(fd, "w", encoding="utf-8") as handle:
                handle.write(rendered)
                handle.flush()
                os.fsync(handle.fileno())
            os.replace(temporary_name, destination)
        except Exception:
            try:
                os.unlink(temporary_name)
            except OSError:
                pass
            raise
        return destination


# ---------- Validation ----------


@dataclass
class ValidationReport:
    """Outcome of validating one evidence bundle payload."""

    schema_ok: bool = False
    ok: bool = False
    failures: list[str] = field(default_factory=list)
    sections: dict[str, list[dict[str, Any]]] = field(default_factory=dict)

    def fail(self, category: str, detail: str) -> None:
        self.failures.append(f"{category}: {detail}")
        self.ok = False


def _verify_recorded_digest(
    report: ValidationReport,
    label: str,
    record: Mapping[str, Any],
) -> None:
    """Recompute one recorded digest against the bytes on disk."""

    path = Path(str(record.get("path", "")))
    recorded = record.get("sha256")
    if not path.is_file():
        report.fail("missing_artifact", f"{label} {path} is not present on disk")
        return
    try:
        actual = digest_file(path)
    except OSError as exc:
        report.fail("stale", f"{label} {path} could not be re-hashed: {exc}")
        return
    if actual != recorded:
        report.fail(
            "stale",
            f"{label} {path} changed since the bundle was built " f"(recorded {recorded}, actual {actual})",
        )


def validate_bundle(
    payload: Mapping[str, Any],
    *,
    verify_hashes: bool = True,
) -> ValidationReport:
    """Validate one evidence bundle payload.

    Operational, descriptive, and inferential artifacts are separated into
    distinct report sections and are never evaluated as one class. Any
    inconsistency (schema drift, digest drift, mixed roots, role conflation,
    state conflation) fails validation closed.
    """

    report = ValidationReport()
    if not isinstance(payload, Mapping) or payload.get("schema") != EVIDENCE_BUNDLE_SCHEMA:
        report.fail(
            "schema",
            f"payload is not {EVIDENCE_BUNDLE_SCHEMA}; "
            f"got {payload.get('schema') if isinstance(payload, Mapping) else type(payload)!r}",
        )
        return report
    report.schema_ok = True

    body = {key: value for key, value in payload.items() if key != "bundle_id"}
    recorded_id = payload.get("bundle_id")
    if not isinstance(recorded_id, str) or _canonical_digest(body) != recorded_id:
        report.fail("bundle_id", "content digest does not match the bundle body")

    # -- task counts: terminal and unresolved states never overlap --
    task_counts = payload.get("task_counts")
    if not isinstance(task_counts, Mapping):
        report.fail("task_counts", "task_counts section is missing or not an object")
    else:
        terminal = task_counts.get("terminal")
        unresolved = task_counts.get("unresolved")
        if not isinstance(terminal, Mapping) or not isinstance(unresolved, Mapping):
            report.fail("task_counts", "terminal and unresolved count tables are required")
        else:
            overlap = sorted(set(terminal) & set(unresolved))
            if overlap:
                report.fail(
                    "state_conflation",
                    "states counted as both terminal and unresolved: " + ", ".join(overlap),
                )
            unknown = sorted((set(terminal) | set(unresolved)) - TERMINAL_TASK_STATES - UNRESOLVED_TASK_STATES)
            if unknown:
                report.fail("task_counts", "unknown task states: " + ", ".join(unknown))

    # -- artifacts: roles, root confinement, digest freshness --
    artifacts = payload.get("artifacts")
    if not isinstance(artifacts, list):
        report.fail("artifacts", "artifacts section is missing or not a list")
        artifacts = []

    data_root = payload.get("data_root")
    if not isinstance(data_root, str) or not data_root:
        report.fail("data_root", "data_root is missing")
        data_root = ""

    role_index: dict[str, str] = {}
    for record in artifacts:
        if not isinstance(record, Mapping):
            report.fail("artifacts", "artifact record is not an object")
            continue
        path = str(record.get("path", ""))
        role = record.get("role")
        if role not in ALL_ROLES:
            report.fail("unknown_role", f"{path} declares role {role!r}")
            continue
        previous = role_index.get(path)
        if previous is not None and previous != role:
            report.fail(
                "role_conflation",
                f"{path} is bound to both {previous!r} and {role!r}",
            )
        role_index[path] = role
        if data_root and not Path(path).is_relative_to(Path(data_root)):
            report.fail(
                "mixed_root",
                f"{path} is outside the bundle data root {data_root}",
            )

    if verify_hashes:
        for record in artifacts:
            if isinstance(record, Mapping) and record.get("role") in ALL_ROLES:
                _verify_recorded_digest(report, "artifact", record)
        inputs = payload.get("inputs")
        if isinstance(inputs, list):
            for record in inputs:
                if isinstance(record, Mapping):
                    _verify_recorded_digest(report, "input", record)

    # -- role separation into report sections --
    report.sections = {
        "operational": [
            dict(record)
            for record in artifacts
            if isinstance(record, Mapping) and record.get("role") in OPERATIONAL_ROLES
        ],
        "descriptive": [
            dict(record)
            for record in artifacts
            if isinstance(record, Mapping) and record.get("role") in DESCRIPTIVE_ROLES
        ],
        "inferential": [
            dict(record)
            for record in artifacts
            if isinstance(record, Mapping) and record.get("role") in INFERENTIAL_ROLES
        ],
    }

    report.ok = not report.failures
    return report


def load_bundle(path: str | Path) -> dict[str, Any]:
    """Load and schema-check a bundle file, raising :class:`BundleSchemaError`."""

    try:
        payload = json.loads(Path(path).read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise BundleSchemaError(f"bundle {path} is unreadable: {exc}") from exc
    if not isinstance(payload, dict) or payload.get("schema") != EVIDENCE_BUNDLE_SCHEMA:
        raise BundleSchemaError(f"bundle {path} is not {EVIDENCE_BUNDLE_SCHEMA}")
    return payload


def require_valid_bundle(
    payload: Mapping[str, Any],
    *,
    verify_hashes: bool = True,
) -> ValidationReport:
    """Validate a bundle and raise the mapped exception on any failure."""

    report = validate_bundle(payload, verify_hashes=verify_hashes)
    if not report.ok:
        first = report.failures[0]
        category, _, detail = first.partition(": ")
        if category == "stale":
            raise StaleEvidenceError(detail)
        if category == "mixed_root":
            raise MixedRootError(detail)
        if category == "role_conflation":
            raise RoleConflictError(detail)
        raise EvidenceBundleError(first)
    return report


def main(argv: Sequence[str] | None = None) -> int:
    """CLI entry point: validate one evidence bundle file."""

    parser = argparse.ArgumentParser(description="Validate a first-class RNA evidence bundle.")
    parser.add_argument("bundle", help="path to evidence_bundle.json")
    parser.add_argument(
        "--no-verify-hashes",
        action="store_true",
        help="skip re-hashing artifacts and inputs against the bundle digests",
    )
    args = parser.parse_args(argv)

    try:
        payload = load_bundle(args.bundle)
    except BundleSchemaError as exc:
        print(f"INVALID: {exc}")
        return 1
    report = validate_bundle(payload, verify_hashes=not args.no_verify_hashes)
    print(f"bundle_id: {payload.get('bundle_id')}")
    print(f"data_root: {payload.get('data_root')}")
    for section in ("operational", "descriptive", "inferential"):
        print(f"{section} artifacts: {len(report.sections.get(section, []))}")
    if report.ok:
        print("VALID")
        return 0
    for failure in report.failures:
        print(f"FAIL {failure}")
    print("INVALID")
    return 1


__all__ = [
    "ALL_ROLES",
    "ARTIFACT_KINDS",
    "DESCRIPTIVE_ROLES",
    "EVIDENCE_BUNDLE_FILENAME",
    "EVIDENCE_BUNDLE_SCHEMA",
    "INFERENTIAL_ROLES",
    "OPERATIONAL_ROLES",
    "TERMINAL_TASK_STATES",
    "UNRESOLVED_TASK_STATES",
    "BundleSchemaError",
    "EvidenceBundleBuilder",
    "EvidenceBundleError",
    "MixedRootError",
    "RoleConflictError",
    "StaleEvidenceError",
    "ValidationReport",
    "digest_file",
    "load_bundle",
    "main",
    "require_valid_bundle",
    "validate_bundle",
]
