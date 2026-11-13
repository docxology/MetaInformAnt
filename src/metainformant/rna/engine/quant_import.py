"""Import verified quantifications from other data roots into a campaign data root.

A source is another data root (``<root>/<species>/work/quant/<run>``, layout
``work``) or a restore target (``<root>/<species>/<run>``, layout ``flat``).
Sources are only read. For every sample directory the importer decides, with
reasons, one of:

``import``            a cohort task, valid in its source context, bound to the
                      cohort configuration and index bytes, absent from the
                      destination (``imported`` once executed);
``exists_identical``  the destination already holds the same files;
``exists_different``  the destination holds different files (never
                      overwritten; both versions' bindings are reported);
``skipped``           with ``reference_mismatch``, ``config_mismatch``,
                      ``not_in_cohort``, ``permanent_drop``,
                      ``no_reference_evidence``, ``receipt_mismatch`` or
                      ``invalid_source:<reason>``.

Validation in the source context is the reuse contract v2
(:func:`classify_cohort_quantification` with content verification) plus the
structural and numeric checks of :func:`validate_quantification`, plus, when a
receipt store is configured, a reference-bound receipt whose sidecar, files,
configuration and index digest must all agree. The destination must quantify
the species against the same index bytes (else every sample of that species is
``reference_mismatch``).

Execution copies only the contract files (sidecar, abundance table, run info,
HDF5 and an existing restored-input witness) into a staging directory under the
destination ``quant`` directory, re-hashes every copy, writes a portable
reference-binding witness (labels and relative paths only), checks the staged
sample with the v2 classifier, and publishes it with a single ``rename``. A
destination sample directory is never replaced. Every decision is appended to
a JSON-lines manifest (``<dest>/results/imports/`` by default) and fsynced.
Progress-database adoption is optional; by default the producer's next
reconcile adopts imported outputs.
"""

from __future__ import annotations

import argparse
import concurrent.futures
import contextlib
import json
import os
import re
import secrets
import shutil
import stat as stat_module
import sys
import time
from collections.abc import Callable, Iterable, Iterator, Mapping, Sequence
from dataclasses import asdict, dataclass, field
from pathlib import Path, PurePosixPath
from typing import Any

from metainformant.core.utils.hash import sha256_bytes
from metainformant.core.utils.logging import get_logger
from metainformant.rna.engine.cohort import (
    CohortError,
    DataRootBusy,
    FrozenCohort,
    ReferenceCheck,
    data_root_maintenance_lock,
    load_frozen_cohort,
    select_cohort_species,
    verify_cohort_configs,
    verify_species_reference,
)
from metainformant.rna.engine.provenance import (
    EVIDENCE_BINDING_WITNESS,
    QUANT_PROVENANCE_FILENAME,
    REFERENCE_BINDING_FILENAME,
    RESTORED_INPUTS_DIR,
    REUSABLE_QUANT_STATUSES,
    CohortExpectation,
    _durable_fsync,
    _fsync_directory,
    classify_cohort_quantification,
    quantification_input_paths,
    reference_binding_payload,
    write_reference_binding,
)
from metainformant.rna.engine.quant_validation import validate_quantification
from metainformant.rna.engine.reference_digest import IndexDigestCache

logger = get_logger(__name__)

IMPORT_SCHEMA = "metainformant.rna.quant_import.v1"
STAGING_PREFIX = ".import-"
STAGING_MARKER = ".metainformant_import_staging.json"
STAGING_SCHEMA = "metainformant.rna.quant_import_staging.v1"
RECEIPT_SCHEMA = "metainformant.rna.locked_quant.v1"
RESTORED_INPUT_FILES = ("config.yaml", "manifest.json", "reference_manifest.json")

ACTION_IMPORT = "import"
ACTION_IMPORTED = "imported"
ACTION_EXISTS_IDENTICAL = "exists_identical"
ACTION_EXISTS_DIFFERENT = "exists_different"
ACTION_SKIPPED = "skipped"

REASON_REFERENCE_MISMATCH = "reference_mismatch"
REASON_CONFIG_MISMATCH = "config_mismatch"
REASON_NOT_IN_COHORT = "not_in_cohort"
REASON_PERMANENT_DROP = "permanent_drop"
REASON_NO_REFERENCE_EVIDENCE = "no_reference_evidence"
REASON_RECEIPT_MISMATCH = "receipt_mismatch"

_ACCESSION = re.compile(r"^(SRR|ERR|DRR)\d+$")
_LABEL = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.-]*$")
# Classifier reasons that positively contradict the cohort (as opposed to "no evidence found").
_CONFIG_MISMATCH_MARKERS = ("differs from the cohort configuration", "witness config_sha256 differs")
_REFERENCE_MISMATCH_MARKERS = (
    "different index",
    "different reference index",
    "index bytes differ",
    "reference_index_sha256 differs",
)
# Evidence that is absent or unreachable here (rather than contradicting the cohort) may come from a receipt.
_RESCUABLE_MARKERS = (
    "no reference-index binding evidence",
    "reference index unavailable",
    "reference manifest unreadable",
    "reference manifest is not a complete manifest",
    "reference manifest names no kallisto index",
    "witness cohort_id differs",
)


@dataclass(frozen=True)
class ReceiptStore:
    """A read-only reference-bound receipt store (``<root>/<cohort>/reference-bound-receipts/...``)."""

    label: str
    root: Path
    cohort: str


@dataclass(frozen=True)
class ImportSource:
    """A read-only data root to import from."""

    label: str
    root: Path
    layout: str = "work"

    def species_dir(self, species: str) -> Path:
        return self.root / species / "work" / "quant" if self.layout == "work" else self.root / species

    def relpath(self, species: str, run: str) -> str:
        return f"{species}/work/quant/{run}" if self.layout == "work" else f"{species}/{run}"


def parse_source(spec: str) -> ImportSource:
    """Parse ``label=path[:layout]`` (layout ``work`` or ``flat``)."""

    label, separator, rest = spec.partition("=")
    if not separator or not _LABEL.fullmatch(label) or not rest:
        raise ValueError(f"source must be label=path[:layout]: {spec!r}")
    layout = "work"
    path, colon, suffix = rest.rpartition(":")
    if colon and suffix in ("work", "flat"):
        rest, layout = path, suffix
    return ImportSource(label, Path(rest), layout)


def parse_receipt_store(spec: str) -> ReceiptStore:
    """Parse ``label=path:cohort``."""

    label, separator, rest = spec.partition("=")
    path, colon, cohort = rest.rpartition(":")
    if not separator or not colon or not _LABEL.fullmatch(label) or not _LABEL.fullmatch(cohort) or not path:
        raise ValueError(f"receipt store must be label=path:cohort: {spec!r}")
    return ReceiptStore(label, Path(path), cohort)


@dataclass
class ImportDecision:
    """What the importer decided for one source sample directory, and why."""

    species: str
    accession: str
    source: str
    source_relpath: str
    action: str
    reason: str
    destination_relpath: str
    evidence: dict[str, Any] | None = None
    files: list[dict[str, Any]] = field(default_factory=list)
    details: dict[str, Any] = field(default_factory=dict)

    def as_record(self) -> dict[str, Any]:
        record = asdict(self)
        record["schema"] = IMPORT_SCHEMA
        record["kind"] = "decision"
        return record


# ---------------------------------------------------------------------------
# Reading evidence (never writes)
# ---------------------------------------------------------------------------


def _file_sha256(path: Path) -> str:
    from metainformant.core.utils.hash import sha256_file

    return sha256_file(path)


def receipt_path(store: ReceiptStore, species: str, accession: str) -> Path:
    return store.root / store.cohort / "reference-bound-receipts" / species / f"{accession}.json"


def verify_receipt(
    store: ReceiptStore,
    sample_dir: Path,
    expectation: CohortExpectation,
    *,
    sidecar_bytes: bytes,
    verify_content: bool,
) -> tuple[str, dict[str, Any]]:
    """Check one reference-bound receipt against a sample directory, reading only.

    Returns ``("absent" | "verified" | "reference_mismatch" | "config_mismatch" |
    "mismatch", details)``. A verified receipt binds the sample's sidecar (and,
    with ``verify_content``, every listed file) to the cohort index digest.
    """

    path = receipt_path(store, expectation.species, expectation.run_accession)
    try:
        data = path.read_bytes()
    except FileNotFoundError:
        return "absent", {}
    except OSError as error:
        return "mismatch", {"receipt_error": f"unreadable: {error.strerror}"}
    details: dict[str, Any] = {"receipt_store": store.label, "receipt_cohort": store.cohort}
    details["receipt_sha256"] = sha256_bytes(data)
    try:
        receipt = json.loads(data)
    except ValueError:
        return "mismatch", {**details, "receipt_error": "not JSON"}
    if not isinstance(receipt, dict) or receipt.get("schema") != RECEIPT_SCHEMA:
        return "mismatch", {**details, "receipt_error": "unknown schema"}
    if receipt.get("species") != expectation.species or receipt.get("accession") != expectation.run_accession:
        return "mismatch", {**details, "receipt_error": "identity differs"}
    if receipt.get("reference_index_sha256") != expectation.reference_index_sha256:
        return "reference_mismatch", {
            **details,
            "receipt_reference_index_sha256": receipt.get("reference_index_sha256"),
        }
    if receipt.get("config_sha256") != expectation.config_sha256:
        return "config_mismatch", {**details, "receipt_config_sha256": receipt.get("config_sha256")}
    sidecar = json.loads(sidecar_bytes)
    if receipt.get("contract_id") != sidecar.get("quant_contract_id"):
        return "mismatch", {**details, "receipt_error": "contract differs from the sidecar"}
    files = receipt.get("files")
    if not isinstance(files, list) or not files:
        return "mismatch", {**details, "receipt_error": "lists no files"}
    names = set()
    for entry in files:
        name = entry.get("name") if isinstance(entry, dict) else None
        if not isinstance(name, str) or PurePosixPath(name).name != name or name in ("", ".", ".."):
            return "mismatch", {**details, "receipt_error": "unsafe file name"}
        names.add(name)
        target = sample_dir / name
        if name == QUANT_PROVENANCE_FILENAME:
            if sha256_bytes(sidecar_bytes) != entry.get("sha256"):
                return "mismatch", {**details, "receipt_error": "sidecar bytes differ"}
            continue
        try:
            info = target.lstat()
        except OSError:
            return "mismatch", {**details, "receipt_error": f"listed file missing: {name}"}
        if not stat_module.S_ISREG(info.st_mode) or info.st_size != entry.get("size"):
            return "mismatch", {**details, "receipt_error": f"listed file size differs: {name}"}
        if verify_content and _file_sha256(target) != entry.get("sha256"):
            return "mismatch", {**details, "receipt_error": f"listed file bytes differ: {name}"}
    if QUANT_PROVENANCE_FILENAME not in names:
        return "mismatch", {**details, "receipt_error": "sidecar not listed"}
    return "verified", details


def _first_receipt(
    stores: Sequence[ReceiptStore],
    sample_dir: Path,
    expectation: CohortExpectation,
    sidecar_bytes: bytes,
    verify_content: bool,
) -> tuple[str, dict[str, Any]]:
    """Return the verdict of the first store holding a receipt for this task (``absent`` if none)."""

    for store in stores:
        verdict, details = verify_receipt(
            store, sample_dir, expectation, sidecar_bytes=sidecar_bytes, verify_content=verify_content
        )
        if verdict != "absent":
            return verdict, details
    return "absent", {}


@dataclass(frozen=True)
class SourceVerdict:
    """Whether a sample directory is a reusable output of its cohort task, in its own context."""

    accepted: bool
    reason: str
    evidence: dict[str, Any] | None
    classification: dict[str, Any]
    files: tuple[Path, ...] = ()
    restored_inputs: bool = False
    sidecar_sha256: str | None = None


def judge_sample(
    sample_dir: Path,
    expectation: CohortExpectation,
    stores: Sequence[ReceiptStore],
    digests: IndexDigestCache,
    *,
    source_label: str,
    verify_content: bool = True,
) -> SourceVerdict:
    """Validate one output in its own data-root context and pick the binding evidence (reads only)."""

    classification = classify_cohort_quantification(
        sample_dir, expectation, verify_content=verify_content, index_digests=digests
    )
    status, why = classification.get("status"), str(classification.get("reason"))
    sidecar_path = sample_dir / QUANT_PROVENANCE_FILENAME
    try:
        sidecar_bytes = sidecar_path.read_bytes()
    except OSError:
        return SourceVerdict(False, f"invalid_source:{why}", None, classification)
    try:
        provenance, files, _ = validate_quantification(
            sample_dir,
            expectation.species,
            expectation.run_accession,
            expected_config_sha256=expectation.config_sha256,
        )
    except (OSError, ValueError, KeyError, TypeError) as error:
        text = str(error)
        reason = REASON_CONFIG_MISMATCH if "config_sha256 differs" in text else f"invalid_source:{text}"
        return SourceVerdict(False, reason, None, classification)
    if any(marker in why for marker in _CONFIG_MISMATCH_MARKERS):
        return SourceVerdict(False, REASON_CONFIG_MISMATCH, None, classification)
    if any(marker in why for marker in _REFERENCE_MISMATCH_MARKERS):
        return SourceVerdict(False, REASON_REFERENCE_MISMATCH, None, classification)

    receipt_verdict, receipt_details = _first_receipt(stores, sample_dir, expectation, sidecar_bytes, verify_content)
    if receipt_verdict == "reference_mismatch":
        return SourceVerdict(False, REASON_REFERENCE_MISMATCH, None, {**classification, **receipt_details})
    if receipt_verdict == "config_mismatch":
        return SourceVerdict(False, REASON_CONFIG_MISMATCH, None, {**classification, **receipt_details})
    if receipt_verdict == "mismatch":
        return SourceVerdict(False, REASON_RECEIPT_MISMATCH, None, {**classification, **receipt_details})

    reusable = status in REUSABLE_QUANT_STATUSES
    if not reusable:
        # Only absent or unreachable evidence may be replaced by a receipt; a witness that is present
        # but inconsistent (possible tampering) or a content problem is never rescued.
        if not any(marker in why for marker in _RESCUABLE_MARKERS):
            return SourceVerdict(False, f"invalid_source:{why}", None, classification)
        if receipt_verdict != "verified":
            return SourceVerdict(False, REASON_NO_REFERENCE_EVIDENCE, None, classification)
    source_kind = (classification.get("evidence") or {}).get("kind") if reusable else None
    evidence: dict[str, Any] = {"imported_from": source_label, "source_evidence": source_kind}
    if receipt_verdict == "verified":
        evidence.update(kind="receipt", index_bytes_verified=True, **receipt_details)
    else:
        # A restored witness binds the index by construction: the restore compared the receipt's index
        # digest with the frozen one, which cannot be repeated here without that receipt.
        evidence.update(kind=str(source_kind), index_bytes_verified=source_kind != "restored_witness")

    restored = False
    if (sample_dir / RESTORED_INPUTS_DIR).is_dir() and not (sample_dir / RESTORED_INPUTS_DIR).is_symlink():
        try:
            quantification_input_paths(sample_dir, provenance)
            restored = all(
                (sample_dir / RESTORED_INPUTS_DIR / name).is_file()
                and not (sample_dir / RESTORED_INPUTS_DIR / name).is_symlink()
                for name in RESTORED_INPUT_FILES
            )
        except (OSError, ValueError, KeyError, TypeError):
            restored = False
    return SourceVerdict(
        True,
        "verified",
        evidence,
        classification,
        tuple(files),
        restored,
        sha256_bytes(sidecar_bytes),
    )


# ---------------------------------------------------------------------------
# Planning (reads only)
# ---------------------------------------------------------------------------


def _contract_file_names(sample_dir: Path) -> list[str]:
    """Names of the files that make up one output (sidecar, abundance, run info, HDF5)."""

    try:
        sidecar = json.loads((sample_dir / QUANT_PROVENANCE_FILENAME).read_bytes())
    except (OSError, ValueError):
        return []
    names = [QUANT_PROVENANCE_FILENAME]
    quant_name = sidecar.get("quantification_file") if isinstance(sidecar, dict) else None
    if isinstance(quant_name, str) and PurePosixPath(quant_name).name == quant_name:
        names.append(quant_name)
    with contextlib.suppress(OSError):
        for entry in os.scandir(sample_dir):
            if entry.name.endswith("run_info.json") or entry.name.endswith(".h5"):
                names.append(entry.name)
    return sorted(set(names))


def _compare_outputs(source_dir: Path, destination_dir: Path, *, deep: bool) -> tuple[bool, dict[str, Any]]:
    """Return ``(identical, details)``; identical means same contract files with the same bytes.

    Without ``deep``, identity is the sidecar bytes (which pin the abundance
    sha256) plus equal sizes of every file; ``deep`` hashes every file.
    """

    names = _contract_file_names(source_dir)
    details: dict[str, Any] = {"identity_basis": "sha256" if deep else "sidecar_bytes+sizes"}
    if not names:
        return False, {**details, "difference": "source has no readable sidecar"}
    try:
        source_sidecar = (source_dir / QUANT_PROVENANCE_FILENAME).read_bytes()
        destination_sidecar = (destination_dir / QUANT_PROVENANCE_FILENAME).read_bytes()
    except OSError:
        return False, {**details, "difference": "destination has no readable sidecar"}
    details["source_sidecar_sha256"] = sha256_bytes(source_sidecar)
    details["destination_sidecar_sha256"] = sha256_bytes(destination_sidecar)
    if source_sidecar != destination_sidecar:
        return False, {**details, "difference": "sidecar bytes differ"}
    for name in names:
        source_file, destination_file = source_dir / name, destination_dir / name
        try:
            source_info, destination_info = source_file.lstat(), destination_file.lstat()
        except OSError:
            return False, {**details, "difference": f"file missing on one side: {name}"}
        if not stat_module.S_ISREG(destination_info.st_mode) or source_info.st_size != destination_info.st_size:
            return False, {**details, "difference": f"file size differs: {name}"}
        if deep and _file_sha256(source_file) != _file_sha256(destination_file):
            return False, {**details, "difference": f"file bytes differ: {name}"}
    return True, details


def _bound_summary(source_bound: bool, destination_bound: bool) -> str:
    if source_bound and destination_bound:
        return "both"
    if source_bound:
        return "source_only"
    if destination_bound:
        return "destination_only"
    return "neither"


def _quantification_relation(source_sha: str | None, destination_sha: str | None) -> str:
    """``identical``/``different`` abundance tables (by their sidecar-recorded sha256), else ``unknown``."""

    if not source_sha or not destination_sha:
        return "unknown"
    return "identical" if source_sha == destination_sha else "different"


def _index_bytes_proven(verdict: SourceVerdict) -> bool:
    """True when the accepted evidence checked the reference index bytes (not only the index name)."""

    return bool(verdict.accepted and verdict.evidence and verdict.evidence.get("index_bytes_verified") is True)


def _abundance_sha(sample_dir: Path) -> str | None:
    try:
        sidecar = json.loads((sample_dir / QUANT_PROVENANCE_FILENAME).read_bytes())
    except (OSError, ValueError):
        return None
    return sidecar.get("quantification_file_sha256") if isinstance(sidecar, dict) else None


@dataclass
class ImportPlan:
    """Decisions for every source sample plus per-species destination reference checks."""

    cohort_id: str
    definition_sha256: str
    destination: Path
    decisions: list[ImportDecision]
    references: dict[str, ReferenceCheck]
    ignored: dict[str, int]

    def counts(self) -> dict[str, dict[str, int]]:
        table: dict[str, dict[str, int]] = {}
        for decision in self.decisions:
            key = decision.action if decision.action != ACTION_SKIPPED else f"skipped:{decision.reason}"
            row = table.setdefault(f"{decision.source}/{decision.species}", {})
            row[key] = row.get(key, 0) + 1
        return table


def _iter_source_samples(source: ImportSource, species: str, ignored: dict[str, int]) -> Iterator[tuple[str, Path]]:
    directory = source.species_dir(species)
    if not directory.is_dir():
        return
    with os.scandir(directory) as entries:
        names = sorted((entry.name, Path(entry.path)) for entry in entries if entry.is_dir(follow_symlinks=False))
    for name, path in names:
        if name.startswith("."):
            ignored["dot_entries"] = ignored.get("dot_entries", 0) + 1
            continue
        if not _ACCESSION.fullmatch(name):
            ignored["non_accession_entries"] = ignored.get("non_accession_entries", 0) + 1
            continue
        yield name, path


def plan_imports(
    cohort: FrozenCohort,
    destination: str | Path,
    sources: Sequence[ImportSource],
    *,
    receipt_stores: Sequence[ReceiptStore] = (),
    destination_receipt_stores: Sequence[ReceiptStore] = (),
    species: Iterable[str] | None = None,
    verify_content: bool = True,
    deep_compare: bool = False,
    workers: int = 4,
    index_digests: IndexDigestCache | None = None,
) -> ImportPlan:
    """Decide, reading only, what importing ``sources`` into ``destination`` would do.

    Sources are processed in order; a run already planned for import from an
    earlier source is reported as ``exists_identical``/``exists_different``
    against that planned copy's source only when the destination already holds
    it, otherwise as ``skipped`` with ``duplicate_source:<label>``.
    ``destination_receipt_stores`` are consulted when judging a destination
    version that conflicts with a source.
    """

    root = Path(destination)
    digests = index_digests if index_digests is not None else IndexDigestCache()
    selected = select_cohort_species(cohort, species)
    references = {name: verify_species_reference(cohort, root, name, index_digests=digests) for name in selected}
    ignored: dict[str, int] = {}
    work: list[tuple[ImportSource, str, str, Path]] = []
    for source in sources:
        for name in selected:
            for accession, path in _iter_source_samples(source, name, ignored):
                work.append((source, name, accession, path))

    def decide(item: tuple[ImportSource, str, str, Path]) -> ImportDecision:
        source, name, accession, path = item
        destination_rel = f"{name}/work/quant/{accession}"
        decision = ImportDecision(
            species=name,
            accession=accession,
            source=source.label,
            source_relpath=source.relpath(name, accession),
            action=ACTION_SKIPPED,
            reason="",
            destination_relpath=destination_rel,
        )
        if cohort.is_excluded(name, accession):
            decision.reason = REASON_PERMANENT_DROP
            return decision
        if cohort.task(name, accession) is None:
            decision.reason = REASON_NOT_IN_COHORT
            return decision
        expectation = cohort.expectation(name, accession)
        destination_dir = root / destination_rel
        exists = destination_dir.exists() or destination_dir.is_symlink()
        if exists:
            identical, comparison = _compare_outputs(path, destination_dir, deep=deep_compare)
            decision.details["comparison"] = comparison
            if identical:
                decision.action, decision.reason = ACTION_EXISTS_IDENTICAL, "destination holds the same files"
                return decision
            source_abundance, destination_abundance = _abundance_sha(path), _abundance_sha(destination_dir)
            quantification = _quantification_relation(source_abundance, destination_abundance)
            decision.action = ACTION_EXISTS_DIFFERENT
            decision.reason = (
                "destination holds the same quantification with different runtime files"
                if quantification == "identical"
                else "destination holds a different quantification"
            )
            verdict = judge_sample(
                path, expectation, receipt_stores, digests, source_label=source.label, verify_content=verify_content
            )
            destination_verdict = judge_sample(
                destination_dir,
                expectation,
                destination_receipt_stores,
                digests,
                source_label="destination",
                verify_content=verify_content,
            )
            decision.details.update(
                source_bound=verdict.accepted,
                source_verdict=verdict.reason,
                source_evidence=verdict.evidence,
                source_contract_id=verdict.classification.get("contract_id"),
                source_abundance_sha256=source_abundance,
                destination_bound=destination_verdict.accepted,
                destination_verdict=destination_verdict.reason,
                destination_evidence=destination_verdict.evidence,
                destination_contract_id=destination_verdict.classification.get("contract_id"),
                destination_abundance_sha256=destination_abundance,
                quantification=quantification,
                bound=_bound_summary(verdict.accepted, destination_verdict.accepted),
                index_bytes_proven=_bound_summary(
                    _index_bytes_proven(verdict), _index_bytes_proven(destination_verdict)
                ),
                destination_reference=references[name].reason,
            )
            return decision
        if not references[name].ok:
            decision.reason = REASON_REFERENCE_MISMATCH
            decision.details["destination_reference"] = references[name].reason
            return decision
        verdict = judge_sample(
            path, expectation, receipt_stores, digests, source_label=source.label, verify_content=verify_content
        )
        decision.details["source_status"] = verdict.classification.get("status")
        if not verdict.accepted:
            decision.reason = verdict.reason
            decision.details["source_reason"] = verdict.classification.get("reason")
            return decision
        decision.action, decision.reason = ACTION_IMPORT, "verified in source context"
        decision.evidence = verdict.evidence
        decision.files = [
            {"name": file.name, "size": file.stat().st_size, "sha256": _file_sha256(file) if verify_content else None}
            for file in verdict.files
        ]
        decision.details["restored_inputs"] = verdict.restored_inputs
        decision.details["sidecar_sha256"] = verdict.sidecar_sha256
        return decision

    if workers > 1 and len(work) > 1:
        with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as pool:
            decisions = list(pool.map(decide, work))
    else:
        decisions = [decide(item) for item in work]

    planned: dict[tuple[str, str], str] = {}
    for decision in decisions:
        if decision.action != ACTION_IMPORT:
            continue
        key = (decision.species, decision.accession)
        if key in planned:
            decision.action, decision.reason = ACTION_SKIPPED, f"duplicate_source:{planned[key]}"
        else:
            planned[key] = decision.source
    return ImportPlan(cohort.cohort_id, cohort.definition_sha256, root, decisions, references, ignored)


# ---------------------------------------------------------------------------
# Manifest
# ---------------------------------------------------------------------------


class _Manifest:
    """Append-only JSON-lines writer; every line is flushed and fsynced."""

    def __init__(self, path: Path) -> None:
        self.path = path
        path.parent.mkdir(parents=True, exist_ok=True)
        self._handle = open(path, "x", encoding="utf-8")

    def write(self, record: Mapping[str, Any]) -> None:
        self._handle.write(json.dumps(record, sort_keys=True, ensure_ascii=True) + "\n")
        self._handle.flush()
        _durable_fsync(self._handle.fileno())

    def close(self) -> None:
        self._handle.close()
        _fsync_directory(self.path.parent)


def _publishing_record(manifest: "_Manifest", decision: ImportDecision) -> Callable[[], None]:
    """The intent record written (and fsynced) just before a staged copy is published."""

    def announce() -> None:
        manifest.write(
            {
                "kind": "publishing",
                "species": decision.species,
                "accession": decision.accession,
                "source": decision.source,
                "source_relpath": decision.source_relpath,
                "destination_relpath": decision.destination_relpath,
            }
        )

    return announce


def _utc_stamp() -> str:
    return time.strftime("%Y%m%dT%H%M%SZ", time.gmtime())


def default_manifest_path(destination: Path, label: str) -> Path:
    return destination / "results" / "imports" / f"{_utc_stamp()}-{label}-{secrets.token_hex(3)}.jsonl"


def _source_key(root: str | Path) -> str:
    return os.path.normpath(os.path.abspath(os.path.expanduser(str(root))))


def executed_imports(destination: str | Path, *, definition_sha256: str | None = None) -> dict[str, Any]:
    """The ``--execute`` imports recorded in ``<destination>/results/imports/*.jsonl``; reads only.

    Returns ``{"finished": {source_root: record}, "unfinished": [manifest, ...]}``.
    ``finished`` maps each source root (absolute, normalised) to the latest
    execute manifest that listed it and ran to its summary record: ``manifest``,
    ``label``, ``cohort_definition_sha256``, ``finished_at`` and that run's
    ``imported`` and ``failed`` counts. An execute manifest without a summary
    (a killed import) is listed under ``unfinished``. With
    ``definition_sha256`` only manifests of that cohort definition count. A
    launcher uses this to refuse a first producer run that would recompute
    outputs an import should have brought in.
    """

    directory = Path(destination) / "results" / "imports"
    finished: dict[str, dict[str, Any]] = {}
    unfinished: list[str] = []
    try:
        manifests = sorted(directory.glob("*.jsonl"))
    except OSError:
        manifests = []
    for manifest in manifests:
        header: dict[str, Any] | None = None
        summary: dict[str, Any] | None = None
        try:
            with manifest.open("r", encoding="utf-8") as handle:
                for line in handle:
                    try:
                        record = json.loads(line)
                    except ValueError:
                        continue  # a killed import can leave a partial last line
                    if not isinstance(record, dict):
                        continue
                    if header is None:
                        if record.get("kind") != "header" or record.get("mode") != "execute":
                            break
                        header = record
                    elif record.get("kind") == "summary" and isinstance(record.get("execution"), dict):
                        summary = record
        except OSError:
            continue
        if header is None or header.get("schema") != IMPORT_SCHEMA:
            continue
        if definition_sha256 is not None and header.get("cohort_definition_sha256") != definition_sha256:
            continue
        if summary is None:
            unfinished.append(str(manifest))
            continue
        execution = summary["execution"]
        for source in header.get("sources") or ():
            if not isinstance(source, dict) or not source.get("root"):
                continue
            key = _source_key(source["root"])
            record = {
                "manifest": str(manifest),
                "label": source.get("label"),
                "cohort_definition_sha256": header.get("cohort_definition_sha256"),
                "finished_at": summary.get("finished_at"),
                "imported": int(execution.get("imported", 0) or 0),
                "failed": int(execution.get("failed", 0) or 0),
            }
            previous = finished.get(key)
            if previous is None or str(record["finished_at"] or "") >= str(previous["finished_at"] or ""):
                finished[key] = record
    return {"finished": finished, "unfinished": unfinished}


def _refuse_inside(path: Path, roots: Iterable[Path]) -> None:
    resolved = path.absolute().resolve()
    for root in roots:
        if resolved.is_relative_to(Path(root).absolute().resolve()):
            raise ValueError(f"manifest path must not be inside a source root: {path}")


def write_plan_manifest(plan: ImportPlan, path: Path, sources: Sequence[ImportSource], *, mode: str = "plan") -> Path:
    """Write a dry-run manifest to ``path`` (never inside a source root)."""

    _refuse_inside(path, [source.root for source in sources])
    manifest = _Manifest(path)
    try:
        manifest.write(_header(plan, sources, mode))
        for decision in plan.decisions:
            manifest.write(decision.as_record())
        manifest.write(_footer(plan))
    finally:
        manifest.close()
    return path


def _header(plan: ImportPlan, sources: Sequence[ImportSource], mode: str) -> dict[str, Any]:
    return {
        "schema": IMPORT_SCHEMA,
        "kind": "header",
        "mode": mode,
        "cohort_id": plan.cohort_id,
        "cohort_definition_sha256": plan.definition_sha256,
        "destination": str(plan.destination),
        "sources": [{"label": s.label, "root": str(s.root), "layout": s.layout} for s in sources],
        "destination_references": {name: asdict(check) for name, check in sorted(plan.references.items())},
        "started_at": _utc_stamp(),
    }


def _footer(plan: ImportPlan) -> dict[str, Any]:
    return {
        "schema": IMPORT_SCHEMA,
        "kind": "summary",
        "counts": plan.counts(),
        "ignored": plan.ignored,
        "finished_at": _utc_stamp(),
    }


# ---------------------------------------------------------------------------
# Execution (writes only under the destination root)
# ---------------------------------------------------------------------------


def _clean_stale_staging(quant_dir: Path) -> list[str]:
    """Remove ``.import-*`` staging directories that carry this module's marker; return their names."""

    removed: list[str] = []
    if not quant_dir.is_dir():
        return removed
    with os.scandir(quant_dir) as entries:
        candidates = [
            Path(e.path) for e in entries if e.name.startswith(STAGING_PREFIX) and e.is_dir(follow_symlinks=False)
        ]
    for candidate in candidates:
        try:
            marker = json.loads((candidate / STAGING_MARKER).read_bytes())
        except (OSError, ValueError):
            continue
        if isinstance(marker, dict) and marker.get("schema") == STAGING_SCHEMA:
            shutil.rmtree(candidate)
            removed.append(candidate.name)
    if removed:
        _fsync_directory(quant_dir)
    return removed


def _copy_verified(source: Path, destination: Path, expected_sha256: str | None) -> str:
    """Copy one regular file into a new destination, fsync it and return its sha256."""

    if source.is_symlink() or not source.is_file():
        raise ValueError(f"not a regular file: {source.name}")
    with source.open("rb") as reader, open(destination, "xb") as writer:
        shutil.copyfileobj(reader, writer, 4 * 1024 * 1024)
        writer.flush()
        _durable_fsync(writer.fileno())
    digest = _file_sha256(destination)
    if expected_sha256 is not None and digest != expected_sha256:
        raise ValueError(f"copied file differs from its source digest: {source.name}")
    if digest != _file_sha256(source):
        raise ValueError(f"source file changed while copying: {source.name}")
    return digest


def _stage_and_publish(
    cohort: FrozenCohort,
    destination_root: Path,
    source: ImportSource,
    decision: ImportDecision,
    sample_dir: Path,
    announce: Callable[[], None] | None = None,
) -> dict[str, Any]:
    """Copy one verified sample into the destination without ever replacing an existing directory.

    ``announce`` runs (and must have made its record durable) just before the
    staged copy is renamed into place, so a kill between the rename and the
    import's own record still leaves a record of the publication.
    """

    quant_dir = destination_root / decision.species / "work" / "quant"
    quant_dir.mkdir(parents=True, exist_ok=True)
    final = quant_dir / decision.accession
    if final.exists() or final.is_symlink():
        raise FileExistsError(f"destination appeared before publishing: {decision.destination_relpath}")
    staging = quant_dir / f"{STAGING_PREFIX}{decision.accession}-{secrets.token_hex(4)}"
    staging.mkdir()
    marker = {
        "schema": STAGING_SCHEMA,
        "accession": decision.accession,
        "source": source.label,
        "pid": os.getpid(),
        "created_at": _utc_stamp(),
    }
    (staging / STAGING_MARKER).write_text(json.dumps(marker, sort_keys=True) + "\n", encoding="utf-8")
    _fsync_directory(quant_dir)
    staged = staging / decision.accession
    staged.mkdir()
    try:
        copied: list[dict[str, Any]] = []
        for entry in decision.files:
            digest = _copy_verified(sample_dir / entry["name"], staged / entry["name"], entry.get("sha256"))
            copied.append({"name": entry["name"], "sha256": digest, "size": (staged / entry["name"]).stat().st_size})
        if decision.details.get("restored_inputs"):
            (staged / RESTORED_INPUTS_DIR).mkdir()
            for name in RESTORED_INPUT_FILES:
                _copy_verified(sample_dir / RESTORED_INPUTS_DIR / name, staged / RESTORED_INPUTS_DIR / name, None)
            _fsync_directory(staged / RESTORED_INPUTS_DIR)
        sidecar_bytes = (staged / QUANT_PROVENANCE_FILENAME).read_bytes()
        if sha256_bytes(sidecar_bytes) != decision.details.get("sidecar_sha256"):
            raise ValueError("staged sidecar differs from the verified sidecar")
        expectation = cohort.expectation(decision.species, decision.accession)
        evidence = {key: value for key, value in (decision.evidence or {}).items() if value is not None}
        evidence["source_relpath"] = decision.source_relpath
        payload = reference_binding_payload(sidecar_bytes, json.loads(sidecar_bytes), expectation, evidence)
        write_reference_binding(staged, payload)
        _fsync_directory(staged)
        check = classify_cohort_quantification(staged, expectation, verify_content=True)
        if check.get("status") not in REUSABLE_QUANT_STATUSES or (check.get("evidence") or {}).get("kind") != (
            EVIDENCE_BINDING_WITNESS
        ):
            raise ValueError(f"staged copy does not classify as reusable: {check.get('reason')}")
        if final.exists() or final.is_symlink():
            raise FileExistsError(f"destination appeared before publishing: {decision.destination_relpath}")
        if announce is not None:
            announce()
        # rename(2) refuses a non-empty destination directory, so existing outputs can never be replaced.
        os.rename(staged, final)
        _fsync_directory(quant_dir)
    except BaseException:
        shutil.rmtree(staging, ignore_errors=True)
        raise
    shutil.rmtree(staging)
    _fsync_directory(quant_dir)
    return {"copied": copied, "binding_sha256": sha256_bytes((final / REFERENCE_BINDING_FILENAME).read_bytes())}


def execute_imports(
    cohort: FrozenCohort,
    destination: str | Path,
    sources: Sequence[ImportSource],
    *,
    receipt_stores: Sequence[ReceiptStore] = (),
    destination_receipt_stores: Sequence[ReceiptStore] = (),
    species: Iterable[str] | None = None,
    manifest_path: Path | None = None,
    adopt_db: Path | None = None,
    workers: int = 4,
) -> tuple[ImportPlan, Path, dict[str, Any]]:
    """Plan and import under the data-root maintenance lock; return the plan, manifest and summary.

    Each ``import`` decision is re-checked at copy time; a destination that
    appeared in the meantime is recorded as ``exists_different`` and left
    alone. With ``adopt_db``, imported runs are marked ``quantified`` in that
    progress database in one transaction per species (no producer can be
    running: the maintenance lock is the producer's writer lock).
    """

    root = Path(destination)
    species = None if species is None else tuple(species)
    summary: dict[str, Any] = {"imported": 0, "failed": 0, "stale_staging_removed": [], "adopted": 0}
    with data_root_maintenance_lock(root, purpose="import"):
        # A killed import leaves its staging directory (a full copy of the sample) behind: remove every one in the
        # selected species now, whether or not this run imports into that quantification directory again.
        cleaned: set[Path] = set()
        for name in select_cohort_species(cohort, species):
            quant_dir = root / name / "work" / "quant"
            summary["stale_staging_removed"].extend(_clean_stale_staging(quant_dir))
            cleaned.add(quant_dir)
        plan = plan_imports(
            cohort,
            root,
            sources,
            receipt_stores=receipt_stores,
            destination_receipt_stores=destination_receipt_stores,
            species=species,
            verify_content=True,
            workers=workers,
        )
        label = "+".join(source.label for source in sources) or "none"
        path = manifest_path or default_manifest_path(root, label)
        _refuse_inside(path, [source.root for source in sources])
        manifest = _Manifest(path)
        by_label = {source.label: source for source in sources}
        adopted: dict[str, list[dict[str, Any]]] = {}
        try:
            manifest.write(_header(plan, sources, "execute"))
            for decision in plan.decisions:
                if decision.action == ACTION_IMPORT:
                    source = by_label[decision.source]
                    sample_dir = source.root / decision.source_relpath
                    quant_dir = root / decision.species / "work" / "quant"
                    if quant_dir not in cleaned:
                        summary["stale_staging_removed"].extend(_clean_stale_staging(quant_dir))
                        cleaned.add(quant_dir)
                    try:
                        decision.details["published"] = _stage_and_publish(
                            cohort, root, source, decision, sample_dir, announce=_publishing_record(manifest, decision)
                        )
                        decision.action = ACTION_IMPORTED
                        summary["imported"] += 1
                        adopted.setdefault(decision.species, []).append({"srr_id": decision.accession})
                    except FileExistsError as error:
                        decision.action, decision.reason = ACTION_EXISTS_DIFFERENT, str(error)
                    except (OSError, ValueError) as error:
                        decision.action, decision.reason = ACTION_SKIPPED, f"import_failed:{error}"
                        summary["failed"] += 1
                manifest.write(decision.as_record())
            if adopt_db is not None and adopted:
                summary["adopted"] = _adopt(cohort, root, adopt_db, adopted)
                summary["adoption"] = "adopted"
            manifest.write({**_footer(plan), "execution": {k: v for k, v in summary.items()}})
        finally:
            manifest.close()
    return plan, path, summary


def _adopt(cohort: FrozenCohort, root: Path, db_path: Path, adopted: Mapping[str, list[dict[str, Any]]]) -> int:
    from metainformant.rna.engine.progress_db import ProgressDB

    count = 0
    with ProgressDB(db_path) as db:
        for name, items in adopted.items():
            entries = []
            for item in items:
                expectation = cohort.expectation(name, item["srr_id"])
                audit = classify_cohort_quantification(root / name / "work" / "quant" / item["srr_id"], expectation)
                entries.append({"srr_id": item["srr_id"], **audit})
            result = db.adopt_quantified(
                name, [entry for entry in entries if entry["status"] in REUSABLE_QUANT_STATUSES]
            )
            count += len(result["adopted"])
    return count


# ---------------------------------------------------------------------------
# Command line
# ---------------------------------------------------------------------------


def format_counts(plan: ImportPlan) -> str:
    """Render per-source/species action counts as aligned text."""

    lines = []
    for key, row in sorted(plan.counts().items()):
        lines.append(f"{key:45s} " + ", ".join(f"{name}={count}" for name, count in sorted(row.items())))
    return "\n".join(lines)


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Import verified quantifications into a campaign data root.")
    parser.add_argument("--data-root", type=Path, required=True, help="Destination data root")
    parser.add_argument("--cohort", type=Path, required=True, help="Frozen cohort directory")
    parser.add_argument("--cohort-definition", type=Path, help="Definition file (default: <cohort>/COHORT.json)")
    parser.add_argument("--cohort-definition-sha256", help="Refuse unless the definition has this sha256")
    parser.add_argument("--config-dir", type=Path, required=True, help="Directory of the cohort's species YAMLs")
    parser.add_argument("--source", action="append", default=[], required=True, help="label=path[:work|flat]")
    parser.add_argument("--receipt-store", action="append", default=[], help="label=path:cohort (source evidence)")
    parser.add_argument(
        "--destination-receipt-store", action="append", default=[], help="label=path:cohort (conflict reporting)"
    )
    parser.add_argument("--species", action="append", help="Restrict to these cohort species")
    parser.add_argument("--manifest", type=Path, help="Manifest path (required for a dry run)")
    parser.add_argument("--execute", action="store_true", help="Copy verified samples (default: dry run)")
    parser.add_argument("--adopt-db", type=Path, help="Also mark imported runs quantified in this progress DB")
    parser.add_argument("--workers", type=int, default=4)
    args = parser.parse_args(argv)
    try:
        sources = [parse_source(spec) for spec in args.source]
        stores = [parse_receipt_store(spec) for spec in args.receipt_store]
        destination_stores = [parse_receipt_store(spec) for spec in args.destination_receipt_store]
        cohort = load_frozen_cohort(
            args.cohort,
            definition_path=args.cohort_definition,
            expected_definition_sha256=args.cohort_definition_sha256,
        )
        selected = select_cohort_species(cohort, args.species)
        verify_cohort_configs(cohort, args.config_dir, selected)
        if args.execute:
            plan, path, summary = execute_imports(
                cohort,
                args.data_root,
                sources,
                receipt_stores=stores,
                destination_receipt_stores=destination_stores,
                species=selected,
                manifest_path=args.manifest,
                adopt_db=args.adopt_db,
                workers=args.workers,
            )
            print(format_counts(plan))
            print(json.dumps({"manifest": str(path), **summary}, sort_keys=True))
            return 1 if summary["failed"] else 0
        if args.manifest is None:
            raise ValueError("a dry run needs --manifest (it never writes into the destination or a source)")
        plan = plan_imports(
            cohort,
            args.data_root,
            sources,
            receipt_stores=stores,
            destination_receipt_stores=destination_stores,
            species=selected,
            workers=args.workers,
        )
        write_plan_manifest(plan, args.manifest, sources)
        print(format_counts(plan))
        print(json.dumps({"manifest": str(args.manifest)}))
        return 0
    except (CohortError, DataRootBusy, ValueError, OSError) as error:
        print(f"error: {error}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
