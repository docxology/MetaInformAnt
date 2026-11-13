"""Frozen RNA-seq cohorts: a hash-pinned, read-only task universe.

A frozen cohort fixes which runs a campaign quantifies, against which
configuration bytes and which kallisto index bytes. It is the answer to "what
is the work?" that never silently changes: no NCBI metadata re-query, no
rewritten selection table, no pruned progress rows.

On disk a cohort directory holds the frozen acquisition snapshot exactly as it
was frozen, plus one definition file::

    <cohort_dir>/
        COHORT.json                                   # definition (this module)
        snapshot.json  manifest.jsonl  inventory.json # frozen envelope, verbatim
        inputs/<species>/work/metadata/metadata_selected.tsv
        inputs/<species>/work/index/<Stem>.idx

``COHORT.json`` (schema ``metainformant.rna.frozen_cohort.v1``) names the
cohort, lists every file with its sha256 and size, records per-species
configuration/metadata/index digests and carries the permanent exclusions. It
contains no timestamps of its own, so building it twice from the same frozen
source gives identical bytes and its sha256 can be pinned in a repository.

:func:`load_frozen_cohort` verifies every listed file by content, the manifest
and inventory against the snapshot, and that each task's ``batch_index`` is its
1-based row position in the frozen selection table (which ``amalgkit quant
--batch`` relies on). Any mismatch raises :class:`CohortError`; nothing is
repaired.
"""

from __future__ import annotations

import argparse
import contextlib
import csv
import json
import os
import re
import shutil
import stat as stat_module
import sys
import tempfile
from collections.abc import Iterable, Iterator, Mapping, Sequence
from dataclasses import dataclass, field
from pathlib import Path, PurePosixPath
from typing import Any

from metainformant.core.utils.hash import sha256_bytes, sha256_file
from metainformant.rna.engine.acquisition_manifest import load_snapshot
from metainformant.rna.engine.campaign_lock import (
    LEGACY_DIRECTORY_LOCKS,
    LOCK_FILENAME as CAMPAIGN_LOCK_FILENAME,
    CampaignLockBusy,
    acquire_writer,
)
from metainformant.rna.engine.provenance import (
    REFERENCE_MANIFEST_NAME,
    CohortExpectation,
    _durable_fsync,
    _fsync_directory,
    _manifest_index_digests,
)
from metainformant.rna.engine.reference_digest import IndexDigestCache

COHORT_SCHEMA = "metainformant.rna.frozen_cohort.v1"
COHORT_DEFINITION_FILENAME = "COHORT.json"
SNAPSHOT_FILENAME = "snapshot.json"
MANIFEST_FILENAME = "manifest.jsonl"
INVENTORY_FILENAME = "inventory.json"
PERMANENT_DROP = "permanent_drop"
# Maintenance commands that write into a data root (import, seal, DB migration) hold the campaign writer
# flock, the same one the producer holds (see metainformant.rna.engine.campaign_lock).
DATA_ROOT_LOCK_FILENAME = CAMPAIGN_LOCK_FILENAME
# Directory locks created by the shell launchers; their presence means a producer may be writing.
LEGACY_LAUNCHER_LOCKS = LEGACY_DIRECTORY_LOCKS

_IDENTIFIER = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.-]*$")
_SPECIES = re.compile(r"^[a-z0-9][a-z0-9_]*$")
# The values ``amalgkit`` (strtobool) treats as true in ``is_sampled``.
_TRUE_FLAGS = frozenset({"yes", "y", "true", "t", "on", "1"})


class CohortError(ValueError):
    """The frozen cohort definition, or a file it names, cannot be trusted."""


@dataclass(frozen=True)
class CohortTask:
    """One frozen (species, run) task with the facts the producer needs."""

    species: str
    accession: str
    batch_index: int
    expected_paired: bool | None
    total_bases: float
    fastq_bytes: int
    config_name: str
    config_sha256: str
    reference_index_sha256: str


@dataclass(frozen=True)
class CohortSpecies:
    """Per-species digests pinned by the cohort definition."""

    species: str
    config_name: str
    config_sha256: str
    metadata_relpath: str
    metadata_sha256: str
    index_relpath: str
    index_name: str
    index_sha256: str
    task_count: int
    exclusion_count: int


@dataclass(frozen=True, eq=False)
class FrozenCohort:
    """A verified frozen cohort. Construct it with :func:`load_frozen_cohort`."""

    cohort_id: str
    root: Path
    definition_sha256: str
    definition: Mapping[str, Any]
    species: Mapping[str, CohortSpecies]
    tasks: tuple[CohortTask, ...]
    exclusions: Mapping[str, Mapping[str, str]]
    _index: Mapping[tuple[str, str], CohortTask] = field(repr=False)

    @property
    def species_names(self) -> tuple[str, ...]:
        return tuple(self.species)

    @property
    def config_names(self) -> tuple[str, ...]:
        return tuple(entry.config_name for entry in self.species.values())

    def entry(self, species: str) -> CohortSpecies:
        try:
            return self.species[species]
        except KeyError:
            raise CohortError(f"species is not part of cohort {self.cohort_id}: {species}") from None

    def species_for_config(self, config_name: str) -> str:
        for entry in self.species.values():
            if entry.config_name == config_name:
                return entry.species
        raise CohortError(f"configuration is not part of cohort {self.cohort_id}: {config_name}")

    def metadata_path(self, species: str) -> Path:
        return self.root / self.entry(species).metadata_relpath

    def index_path(self, species: str) -> Path:
        return self.root / self.entry(species).index_relpath

    def tasks_for(self, species: str) -> tuple[CohortTask, ...]:
        self.entry(species)
        return tuple(sorted((t for t in self.tasks if t.species == species), key=lambda t: t.batch_index))

    def task(self, species: str, accession: str) -> CohortTask | None:
        return self._index.get((species, accession))

    def exclusions_for(self, species: str) -> Mapping[str, str]:
        return self.exclusions.get(species, {})

    def is_excluded(self, species: str, accession: str) -> bool:
        return accession in self.exclusions.get(species, {})

    def expectation(self, species: str, accession: str) -> CohortExpectation:
        """The reuse expectation for one task (raises ``KeyError`` for a non-task)."""

        task = self._index.get((species, accession))
        if task is None:
            raise KeyError(f"{species}/{accession} is not a task of cohort {self.cohort_id}")
        entry = self.species[species]
        return CohortExpectation(
            cohort_id=self.cohort_id,
            species=species,
            run_accession=accession,
            config_sha256=entry.config_sha256,
            reference_index_sha256=entry.index_sha256,
            reference_index_name=entry.index_name,
        )


# ---------------------------------------------------------------------------
# Small helpers
# ---------------------------------------------------------------------------


def _check_identifier(value: str, label: str) -> str:
    if not isinstance(value, str) or not _IDENTIFIER.fullmatch(value):
        raise CohortError(f"unsafe {label}: {value!r}")
    return value


def _safe_relpath(value: Any) -> str:
    """Return a normalized relative POSIX path, refusing absolute or parent references."""

    if not isinstance(value, str) or not value:
        raise CohortError(f"cohort file path must be a non-empty string: {value!r}")
    path = PurePosixPath(value)
    if path.is_absolute() or "\\" in value or any(part in ("", ".", "..") for part in value.split("/")):
        raise CohortError(f"cohort file path must be relative and normalized: {value!r}")
    return str(path)


def _definition_bytes(definition: Mapping[str, Any]) -> bytes:
    """Serialize a definition deterministically (sorted keys, two-space indent, ASCII)."""

    return (json.dumps(definition, indent=2, sort_keys=True, ensure_ascii=True) + "\n").encode("ascii")


def _regular_file(root: Path, relpath: str) -> Path:
    """Resolve ``relpath`` under ``root`` to a regular, non-symlink file inside the root."""

    candidate = root / relpath
    try:
        mode = candidate.lstat().st_mode
    except OSError as error:
        raise CohortError(f"cohort file is missing: {relpath} ({error.strerror})") from None
    if not stat_module.S_ISREG(mode):
        raise CohortError(f"cohort file is not a regular file: {relpath}")
    if not candidate.resolve().is_relative_to(root.resolve()):
        raise CohortError(f"cohort file escapes the cohort directory: {relpath}")
    return candidate


def _is_true_flag(value: str) -> bool:
    return value.strip().lower() in _TRUE_FLAGS


def _table_lines(data: bytes, relpath: str) -> list[bytes]:
    """Split a TSV into physical lines (line endings kept), refusing quoted fields."""

    if b'"' in data:
        # A quoted field may hide a tab or newline; byte-exact row filtering would then be unsafe.
        raise CohortError(f"frozen selection table contains quoted fields: {relpath}")
    lines = data.splitlines(keepends=True)
    if not lines:
        raise CohortError(f"frozen selection table is empty: {relpath}")
    return lines


def _split_row(line: bytes) -> list[str]:
    return line.rstrip(b"\r\n").decode("utf-8").split("\t")


@dataclass(frozen=True)
class _TableScan:
    columns: list[str]
    runs: list[str]
    sampled: list[bool]
    excluded_flag: list[str]


def _scan_table(data: bytes, relpath: str) -> _TableScan:
    lines = _table_lines(data, relpath)
    columns = _split_row(lines[0])
    if "run" not in columns:
        raise CohortError(f"frozen selection table has no run column: {relpath}")
    if "is_sampled" not in columns:
        raise CohortError(f"frozen selection table has no is_sampled column: {relpath}")
    run_at, sampled_at = columns.index("run"), columns.index("is_sampled")
    exclusion_at = columns.index("exclusion") if "exclusion" in columns else None
    runs: list[str] = []
    sampled: list[bool] = []
    flags: list[str] = []
    for number, line in enumerate(lines[1:], start=2):
        fields = _split_row(line)
        if len(fields) != len(columns):
            raise CohortError(f"frozen selection table row {number} has {len(fields)} fields: {relpath}")
        runs.append(fields[run_at].strip())
        sampled.append(_is_true_flag(fields[sampled_at]))
        flags.append(fields[exclusion_at].strip().lower() if exclusion_at is not None else "")
    return _TableScan(columns, runs, sampled, flags)


# ---------------------------------------------------------------------------
# Building a definition from a frozen source (read-only)
# ---------------------------------------------------------------------------


def _inventory_species(inventory: Mapping[str, Any]) -> list[Mapping[str, Any]]:
    rows = inventory.get("species")
    if not isinstance(rows, list) or not rows:
        raise CohortError("frozen inventory has no species list")
    return rows


def build_cohort_definition(
    source_root: str | Path, *, cohort_id: str, index_digests: IndexDigestCache | None = None
) -> dict[str, Any]:
    """Describe a frozen acquisition snapshot as a cohort definition, reading only.

    ``source_root`` holds ``snapshot.json``, ``manifest.jsonl``, ``inventory.json``
    and the ``inputs/`` tree written by the completion-inventory freeze. Every
    input file is hashed and the result is checked exactly as
    :func:`load_frozen_cohort` checks an installed cohort, so a definition is
    never produced for inconsistent inputs.
    """

    root = Path(source_root)
    _check_identifier(cohort_id, "cohort id")
    digests = index_digests if index_digests is not None else IndexDigestCache()
    try:
        snapshot, _ = load_snapshot(root / MANIFEST_FILENAME)
    except (OSError, ValueError, json.JSONDecodeError) as error:
        raise CohortError(f"frozen snapshot is not usable: {error}") from None
    inventory_path = _regular_file(root, INVENTORY_FILENAME)
    inventory = json.loads(inventory_path.read_bytes())
    if not isinstance(inventory, dict):
        raise CohortError("frozen inventory must be a JSON object")
    inputs = {_safe_relpath(record.get("path")): record.get("sha256") for record in snapshot.get("input_files", [])}

    species_entries: list[dict[str, Any]] = []
    exclusions: list[dict[str, str]] = []
    for row in _inventory_species(inventory):
        species = str(row.get("species"))
        if not _SPECIES.fullmatch(species):
            raise CohortError(f"unsafe species name in inventory: {species!r}")
        index_name = str(row.get("index_name"))
        if PurePosixPath(index_name).name != index_name or not index_name.endswith(".idx"):
            raise CohortError(f"unsafe index name for {species}: {index_name!r}")
        metadata_rel = f"inputs/{species}/work/metadata/metadata_selected.tsv"
        index_rel = f"inputs/{species}/work/index/{index_name}"
        if inputs.get(metadata_rel) != row.get("metadata_sha256"):
            raise CohortError(f"snapshot and inventory disagree on the metadata digest of {species}")
        if inputs.get(index_rel) != row.get("index_sha256"):
            raise CohortError(f"snapshot and inventory disagree on the index digest of {species}")
        excluded = row.get("excluded") or []
        species_entries.append(
            {
                "species": species,
                "config_name": row.get("config_name"),
                "config_sha256": row.get("config_sha256"),
                "metadata_path": metadata_rel,
                "metadata_sha256": row.get("metadata_sha256"),
                "index_path": index_rel,
                "index_name": index_name,
                "index_sha256": row.get("index_sha256"),
                "tasks": len(row.get("tasks") or []),
                "excluded": len(excluded),
            }
        )
        for item in excluded:
            exclusions.append(
                {"species": species, "accession": str(item.get("accession")), "reason_code": str(item.get("reason"))}
            )

    files: list[dict[str, Any]] = []
    for relpath in [SNAPSHOT_FILENAME, MANIFEST_FILENAME, INVENTORY_FILENAME, *sorted(inputs)]:
        path = _regular_file(root, relpath)
        digest = digests.digest(path) if relpath.endswith(".idx") else sha256_file(path)
        expected = inputs.get(relpath)
        if expected is not None and digest != expected:
            raise CohortError(f"frozen input differs from its snapshot digest: {relpath}")
        files.append({"path": relpath, "sha256": digest, "size": path.stat().st_size})

    by_path = {record["path"]: record["sha256"] for record in files}
    if snapshot.get("inventory_sha256") is not None and snapshot["inventory_sha256"] != by_path[INVENTORY_FILENAME]:
        raise CohortError("frozen inventory differs from the snapshot inventory digest")
    definition: dict[str, Any] = {
        "schema": COHORT_SCHEMA,
        "cohort_id": cohort_id,
        "frozen_at": inventory.get("frozen_at"),
        "source_schemas": {"snapshot": snapshot.get("schema"), "inventory": inventory.get("schema")},
        "snapshot_sha256": by_path[SNAPSHOT_FILENAME],
        "manifest_sha256": by_path[MANIFEST_FILENAME],
        "inventory_sha256": by_path[INVENTORY_FILENAME],
        "task_count": int(snapshot.get("task_count", -1)),
        "exclusion_count": len(exclusions),
        "amalgkit": {
            "version": snapshot.get("amalgkit_version"),
            "release_tag": snapshot.get("amalgkit_release_tag"),
            "source_revision": snapshot.get("amalgkit_source_revision"),
        },
        "species": species_entries,
        "exclusions": sorted(exclusions, key=lambda item: (item["species"], item["accession"])),
        "files": files,
    }
    # Validate exactly as a load would, so a bad source never yields a definition.
    _validated_cohort(root, definition, sha256_bytes(_definition_bytes(definition)), digests, files_verified=True)
    return definition


def write_cohort_definition(definition: Mapping[str, Any], destination: str | Path) -> str:
    """Write a definition with its canonical serialization; return its sha256."""

    data = _definition_bytes(definition)
    path = Path(destination)
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, temporary = tempfile.mkstemp(prefix=f".{path.name}.", dir=str(path.parent))
    try:
        with os.fdopen(fd, "wb") as handle:
            handle.write(data)
            handle.flush()
            _durable_fsync(handle.fileno())
        os.chmod(temporary, 0o644)
        os.replace(temporary, path)
    except BaseException:
        with contextlib.suppress(OSError):
            os.unlink(temporary)
        raise
    _fsync_directory(path.parent)
    return sha256_bytes(data)


# ---------------------------------------------------------------------------
# Loading and verifying
# ---------------------------------------------------------------------------


def _verify_files(root: Path, files: Any, digests: IndexDigestCache) -> None:
    if not isinstance(files, list) or not files:
        raise CohortError("cohort definition lists no files")
    seen: set[str] = set()
    for record in files:
        if not isinstance(record, dict):
            raise CohortError("cohort file record must be an object")
        relpath = _safe_relpath(record.get("path"))
        if relpath in seen or relpath == COHORT_DEFINITION_FILENAME:
            raise CohortError(f"cohort file listed twice or self-referential: {relpath}")
        seen.add(relpath)
        path = _regular_file(root, relpath)
        if path.stat().st_size != record.get("size"):
            raise CohortError(f"cohort file size differs from its definition: {relpath}")
        try:
            digest = digests.digest(path)
        except (OSError, ValueError) as error:
            raise CohortError(f"cohort file could not be hashed: {relpath} ({error})") from None
        if digest != record.get("sha256"):
            raise CohortError(f"cohort file content differs from its definition: {relpath}")


def _validated_cohort(
    root: Path,
    definition: Mapping[str, Any],
    definition_sha256: str,
    digests: IndexDigestCache,
    *,
    files_verified: bool = False,
) -> FrozenCohort:
    """Check a definition against the files under ``root`` and build the cohort."""

    if definition.get("schema") != COHORT_SCHEMA:
        raise CohortError(f"unsupported cohort definition schema: {definition.get('schema')!r}")
    cohort_id = _check_identifier(str(definition.get("cohort_id")), "cohort id")
    if not files_verified:
        _verify_files(root, definition.get("files"), digests)
    by_path = {record["path"]: record["sha256"] for record in definition["files"]}
    for name, key in (
        (SNAPSHOT_FILENAME, "snapshot_sha256"),
        (MANIFEST_FILENAME, "manifest_sha256"),
        (INVENTORY_FILENAME, "inventory_sha256"),
    ):
        if by_path.get(name) is None or by_path[name] != definition.get(key):
            raise CohortError(f"cohort definition {key} does not match its file list")
    try:
        snapshot, raw_tasks = load_snapshot(root / MANIFEST_FILENAME)
    except (OSError, ValueError, json.JSONDecodeError) as error:
        raise CohortError(f"frozen manifest is not usable: {error}") from None
    if snapshot.get("manifest_sha256") != definition["manifest_sha256"]:
        raise CohortError("snapshot manifest digest differs from the cohort definition")
    if snapshot.get("inventory_sha256") not in (None, definition["inventory_sha256"]):
        raise CohortError("snapshot inventory digest differs from the cohort definition")
    if len(raw_tasks) != definition.get("task_count"):
        raise CohortError("frozen manifest task count differs from the cohort definition")

    species_rows = definition.get("species")
    if not isinstance(species_rows, list) or not species_rows:
        raise CohortError("cohort definition has no species")
    species: dict[str, CohortSpecies] = {}
    for row in species_rows:
        name = str(row.get("species"))
        if not _SPECIES.fullmatch(name) or name in species:
            raise CohortError(f"unsafe or duplicated cohort species: {name!r}")
        config_name = str(row.get("config_name"))
        if PurePosixPath(config_name).name != config_name or not config_name.endswith(".yaml"):
            raise CohortError(f"unsafe configuration name for {name}: {config_name!r}")
        metadata_rel = _safe_relpath(row.get("metadata_path"))
        index_rel = _safe_relpath(row.get("index_path"))
        if by_path.get(metadata_rel) != row.get("metadata_sha256"):
            raise CohortError(f"cohort metadata digest is not bound to a listed file: {name}")
        if by_path.get(index_rel) != row.get("index_sha256") or PurePosixPath(index_rel).name != row.get("index_name"):
            raise CohortError(f"cohort index digest is not bound to a listed file: {name}")
        for key in ("config_sha256", "metadata_sha256", "index_sha256"):
            if not isinstance(row.get(key), str) or len(row[key]) != 64:
                raise CohortError(f"cohort {key} is missing for {name}")
        species[name] = CohortSpecies(
            species=name,
            config_name=config_name,
            config_sha256=row["config_sha256"],
            metadata_relpath=metadata_rel,
            metadata_sha256=row["metadata_sha256"],
            index_relpath=index_rel,
            index_name=str(row["index_name"]),
            index_sha256=row["index_sha256"],
            task_count=int(row.get("tasks", -1)),
            exclusion_count=int(row.get("excluded", -1)),
        )
    if list(species) != [str(name) for name in snapshot.get("species", list(species))]:
        raise CohortError("cohort species differ from the frozen snapshot species")

    exclusions: dict[str, dict[str, str]] = {name: {} for name in species}
    for item in definition.get("exclusions") or []:
        name, accession, reason = str(item.get("species")), str(item.get("accession")), str(item.get("reason_code"))
        if name not in species or not _IDENTIFIER.fullmatch(accession) or reason != PERMANENT_DROP:
            raise CohortError(f"invalid cohort exclusion: {item!r}")
        if accession in exclusions[name]:
            raise CohortError(f"duplicated cohort exclusion: {name}/{accession}")
        exclusions[name][accession] = reason
    if sum(len(values) for values in exclusions.values()) != definition.get("exclusion_count"):
        raise CohortError("cohort exclusion count differs from its exclusion list")

    inventory = json.loads(_regular_file(root, INVENTORY_FILENAME).read_bytes())
    for row in _inventory_species(inventory):
        name = str(row.get("species"))
        entry = species.get(name)
        if entry is None:
            raise CohortError(f"inventory species is not part of the cohort definition: {name}")
        for key, value in (
            ("config_name", entry.config_name),
            ("config_sha256", entry.config_sha256),
            ("metadata_sha256", entry.metadata_sha256),
            ("index_name", entry.index_name),
            ("index_sha256", entry.index_sha256),
        ):
            if row.get(key) != value:
                raise CohortError(f"inventory {key} differs from the cohort definition for {name}")
        inventory_excluded = {str(item.get("accession")) for item in row.get("excluded") or []}
        if inventory_excluded != set(exclusions[name]):
            raise CohortError(f"inventory exclusions differ from the cohort definition for {name}")

    tasks: list[CohortTask] = []
    index: dict[tuple[str, str], CohortTask] = {}
    accessions_seen: dict[str, str] = {}
    for raw in raw_tasks:
        name = str(raw["species"])
        entry = species.get(name)
        if entry is None:
            raise CohortError(f"manifest task species is not part of the cohort: {name}")
        accession = str(raw["accession"])
        if not _IDENTIFIER.fullmatch(accession):
            raise CohortError(f"unsafe accession in manifest: {accession!r}")
        if accession in accessions_seen:
            raise CohortError(f"accession assigned to two species: {accession}")
        accessions_seen[accession] = name
        if raw.get("config_name") != entry.config_name:
            raise CohortError(f"task configuration name differs from its species entry: {name}/{accession}")
        if raw.get("config_sha256", entry.config_sha256) != entry.config_sha256:
            raise CohortError(f"task configuration digest differs from its species entry: {name}/{accession}")
        if raw.get("reference_index_sha256") != entry.index_sha256:
            raise CohortError(f"task reference index differs from its species entry: {name}/{accession}")
        if accession in exclusions[name]:
            raise CohortError(f"excluded accession is also a task: {name}/{accession}")
        batch_index = raw.get("batch_index")
        if type(batch_index) is not int or batch_index < 1:
            raise CohortError(f"task batch index must be a positive integer: {name}/{accession}")
        paired = raw.get("expected_paired")
        task = CohortTask(
            species=name,
            accession=accession,
            batch_index=batch_index,
            expected_paired=paired if isinstance(paired, bool) else None,
            total_bases=float(raw.get("total_bases") or 0.0),
            fastq_bytes=int(raw.get("fastq_bytes") or 0),
            config_name=entry.config_name,
            config_sha256=entry.config_sha256,
            reference_index_sha256=entry.index_sha256,
        )
        tasks.append(task)
        index[(name, accession)] = task

    for name, entry in species.items():
        own = [task for task in tasks if task.species == name]
        if entry.task_count != len(own) or entry.exclusion_count != len(exclusions[name]):
            raise CohortError(f"cohort task or exclusion count differs from the manifest for {name}")
        relpath = entry.metadata_relpath
        scan = _scan_table(_regular_file(root, relpath).read_bytes(), relpath)
        positions = {run: number for number, run in enumerate(scan.runs, start=1)}
        if len(positions) != len(scan.runs):
            raise CohortError(f"frozen selection table repeats a run: {relpath}")
        expected_runs = {task.accession for task in own} | set(exclusions[name])
        if set(positions) != expected_runs:
            raise CohortError(f"frozen selection table rows differ from tasks plus exclusions: {relpath}")
        for task in own:
            position = positions[task.accession]
            if position != task.batch_index:
                raise CohortError(
                    f"task batch index {task.batch_index} is not the table position {position}: {name}/{task.accession}"
                )
            if not scan.sampled[position - 1]:
                raise CohortError(f"frozen task row is not sampled: {name}/{task.accession}")
            if scan.excluded_flag[position - 1] not in ("", "no"):
                raise CohortError(f"frozen task row carries an Amalgkit exclusion: {name}/{task.accession}")

    return FrozenCohort(
        cohort_id=cohort_id,
        root=root,
        definition_sha256=definition_sha256,
        definition=definition,
        species=species,
        tasks=tuple(tasks),
        exclusions={name: dict(values) for name, values in exclusions.items()},
        _index=index,
    )


def load_frozen_cohort(
    cohort_dir: str | Path,
    *,
    definition_path: str | Path | None = None,
    expected_definition_sha256: str | None = None,
    index_digests: IndexDigestCache | None = None,
) -> FrozenCohort:
    """Load and fully verify a frozen cohort, refusing on any mismatch.

    ``definition_path`` defaults to ``<cohort_dir>/COHORT.json``; passing a
    separate (for example repository-pinned) definition lets an unmodified
    frozen source directory be used read-only as the cohort directory. With
    ``expected_definition_sha256`` the definition bytes themselves are pinned.
    Index files are hashed through ``index_digests`` (keyed by file identity),
    so a second load in the same process costs one ``stat`` per index.
    """

    root = Path(cohort_dir)
    path = Path(definition_path) if definition_path is not None else root / COHORT_DEFINITION_FILENAME
    try:
        data = path.read_bytes()
    except OSError as error:
        raise CohortError(f"cohort definition is unreadable: {path} ({error.strerror})") from None
    digest = sha256_bytes(data)
    if expected_definition_sha256 is not None and digest != expected_definition_sha256:
        raise CohortError(f"cohort definition digest {digest} differs from the pinned {expected_definition_sha256}")
    try:
        definition = json.loads(data)
    except (json.JSONDecodeError, UnicodeDecodeError) as error:
        raise CohortError(f"cohort definition is not JSON: {error}") from None
    if not isinstance(definition, dict):
        raise CohortError("cohort definition must be a JSON object")
    digests = index_digests if index_digests is not None else IndexDigestCache()
    return _validated_cohort(root, definition, digest, digests)


def verify_cohort_configs(
    cohort: FrozenCohort, config_dir: str | Path, species: Iterable[str] | None = None
) -> dict[str, str]:
    """Require every (selected) species configuration to hash to its frozen digest.

    The configuration bytes are part of the quantification contract; a changed
    YAML would invalidate every output of that species. Returns
    ``{species: sha256}``; raises :class:`CohortError` listing every mismatch.
    """

    directory = Path(config_dir)
    selected = list(species) if species is not None else list(cohort.species_names)
    problems: list[str] = []
    digests: dict[str, str] = {}
    for name in selected:
        entry = cohort.entry(name)
        path = directory / entry.config_name
        try:
            if path.is_symlink() or not path.is_file():
                problems.append(f"{entry.config_name}: missing or not a regular file")
                continue
            digest = sha256_file(path)
        except OSError as error:
            problems.append(f"{entry.config_name}: unreadable ({error.strerror})")
            continue
        if digest != entry.config_sha256:
            problems.append(f"{entry.config_name}: sha256 {digest} != frozen {entry.config_sha256}")
            continue
        digests[name] = digest
    if problems:
        raise CohortError("configuration differs from the frozen cohort: " + "; ".join(problems))
    return digests


def select_cohort_species(cohort: FrozenCohort, requested: Iterable[str] | None) -> tuple[str, ...]:
    """Narrow the cohort to requested species (cohort order); unknown names are refused."""

    if not requested:
        return cohort.species_names
    wanted = list(dict.fromkeys(str(name) for name in requested))
    unknown = [name for name in wanted if name not in cohort.species]
    if unknown:
        raise CohortError(f"species not in cohort {cohort.cohort_id}: {', '.join(unknown)}")
    return tuple(name for name in cohort.species_names if name in wanted)


@dataclass(frozen=True)
class ReferenceCheck:
    """Whether a data root's reference for one species is the cohort's index bytes."""

    species: str
    ok: bool
    reason: str
    manifest_relpath: str


def verify_species_reference(
    cohort: FrozenCohort,
    data_root: str | Path,
    species: str,
    *,
    index_digests: IndexDigestCache | None = None,
    index_dir: str | Path | None = None,
) -> ReferenceCheck:
    """Check, reading only, that a data root quantifies ``species`` against the cohort index bytes.

    ``<data_root>/<species>/work/reference/reference_aliases.json`` must be a
    complete manifest for the species whose main index has the cohort index
    name and whose every named index file hashes to the cohort index sha256.
    With ``index_dir`` (the directory ``amalgkit quant --index_dir`` will use),
    the cohort index name inside it must hash to the same digest.
    """

    entry = cohort.entry(species)
    digests = index_digests if index_digests is not None else IndexDigestCache()
    relpath = f"{species}/work/reference/{REFERENCE_MANIFEST_NAME}"
    manifest_path = Path(data_root) / relpath

    def fail(reason: str) -> ReferenceCheck:
        return ReferenceCheck(species, False, reason, relpath)

    try:
        manifest = json.loads(manifest_path.read_bytes())
    except (OSError, ValueError) as error:
        return fail(f"reference manifest unreadable: {error}")
    if not isinstance(manifest, dict) or manifest.get("species") != species or manifest.get("status") != "complete":
        return fail("reference manifest is not a complete manifest for this species")
    main_index = manifest.get("kallisto_index")
    if not isinstance(main_index, str) or PurePosixPath(main_index).name != entry.index_name:
        return fail(f"reference manifest main index is not {entry.index_name}")
    found, index_error = _manifest_index_digests(manifest_path, None, digests)
    if index_error is not None:
        return fail(index_error)
    if found != {entry.index_sha256}:
        return fail("reference index bytes differ from the cohort index")
    if index_dir is not None:
        candidate = Path(index_dir) / entry.index_name
        try:
            digest = digests.digest(candidate)
        except (OSError, ValueError) as error:
            return fail(f"quantification index unavailable: {candidate.name} ({error})")
        if digest != entry.index_sha256:
            return fail("quantification index bytes differ from the cohort index")
    return ReferenceCheck(species, True, "reference index bytes equal the cohort index", relpath)


# ---------------------------------------------------------------------------
# Producer input: the per-species quantification task table
# ---------------------------------------------------------------------------


def cohort_task_table_name(cohort_id: str) -> str:
    """File name of the derived ``amalgkit quant --metadata`` table for a cohort."""

    return f"metadata_cohort_{_check_identifier(cohort_id, 'cohort id')}.tsv"


def cohort_task_table_bytes(
    cohort: FrozenCohort, species: str, *, excluded: Iterable[str] = ()
) -> tuple[bytes, dict[str, int]]:
    """Return the frozen selection rows of a species' schedulable tasks and their batch indexes.

    Rows are byte-identical to the frozen table, in frozen order; rows of
    cohort exclusions and of ``excluded`` accessions are left out. The returned
    mapping gives each kept accession's 1-based position, which is the
    ``amalgkit quant --batch`` value for that table.
    """

    entry = cohort.entry(species)
    drop = set(cohort.exclusions_for(species)) | {str(accession) for accession in excluded}
    data = _regular_file(cohort.root, entry.metadata_relpath).read_bytes()
    lines = _table_lines(data, entry.metadata_relpath)
    columns = _split_row(lines[0])
    run_at = columns.index("run")
    kept = [lines[0]]
    positions: dict[str, int] = {}
    for line in lines[1:]:
        accession = _split_row(line)[run_at].strip()
        if accession in drop:
            continue
        if cohort.task(species, accession) is None:
            raise CohortError(f"frozen selection row is not a cohort task: {species}/{accession}")
        positions[accession] = len(kept)
        kept.append(line if line.endswith(b"\n") else line + b"\n")
    return b"".join(kept), positions


def quant_batch_run(metadata_path: str | Path, batch_index: int) -> str | None:
    """Return the run ``amalgkit quant --batch batch_index --metadata metadata_path`` selects.

    Mirrors Amalgkit 0.16.60: rows whose ``is_sampled`` is true (strtobool:
    yes/y/true/t/on/1) are kept in file order and the ``batch_index``-th one is
    used (all rows when the column is absent). ``None`` when the table is
    unreadable or the batch is out of range.
    """

    try:
        with open(metadata_path, encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            sampled = 0
            for row in reader:
                flag = row.get("is_sampled")
                if flag is not None and not _is_true_flag(str(flag)):
                    continue
                sampled += 1
                if sampled == batch_index:
                    return str(row.get("run") or "").strip()
    except (OSError, UnicodeDecodeError, csv.Error):
        return None
    return None


def write_cohort_task_table(
    cohort: FrozenCohort, species: str, destination: str | Path, *, excluded: Iterable[str] = ()
) -> dict[str, int]:
    """Write the species task table atomically (only when its bytes change); return batch indexes.

    This is a separate file next to, never instead of, the species'
    ``metadata_selected.tsv``. Writing identical bytes is skipped, so a resume
    leaves the file (and its mtime) untouched.
    """

    data, positions = cohort_task_table_bytes(cohort, species, excluded=excluded)
    path = Path(destination)
    try:
        if path.is_file() and not path.is_symlink() and path.read_bytes() == data:
            return positions
    except OSError:
        pass
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, temporary = tempfile.mkstemp(prefix=f".{path.name}.", dir=str(path.parent))
    try:
        with os.fdopen(fd, "wb") as handle:
            handle.write(data)
            handle.flush()
            _durable_fsync(handle.fileno())
        os.chmod(temporary, 0o644)
        os.replace(temporary, path)
    except BaseException:
        with contextlib.suppress(OSError):
            os.unlink(temporary)
        raise
    _fsync_directory(path.parent)
    return positions


# ---------------------------------------------------------------------------
# Installing a cohort directory (copy, never move; never overwrite)
# ---------------------------------------------------------------------------


INSTALL_TEMP_SUFFIX = ".cohort-install.tmp"
# mkstemp's random part (8 characters) after ".<name>." is how temps written before INSTALL_TEMP_SUFFIX look.
_LEGACY_INSTALL_TEMP = re.compile(r"^\.(?P<name>.+)\.[A-Za-z0-9_]{8}$")


def _remove_install_temporaries(target: Path, expected: Iterable[str]) -> list[str]:
    """Remove temp files a killed install left in ``target``; return their relative paths.

    Only files this module names are touched: ``*.cohort-install.tmp`` temps,
    and ``.<name>.<8 random characters>`` temps of a file the definition
    lists in the same directory (the shape before the suffix existed).
    """

    if not target.is_dir():
        return []
    expected_by_dir: dict[Path, set[str]] = {}
    for relpath in [*expected, COHORT_DEFINITION_FILENAME]:
        path = target / relpath
        expected_by_dir.setdefault(path.parent, set()).add(path.name)
    removed: list[str] = []
    for directory in sorted({target, *expected_by_dir}):
        if not directory.is_dir():
            continue
        names = expected_by_dir.get(directory, set())
        with os.scandir(directory) as entries:
            for entry in entries:
                if not entry.name.startswith(".") or not entry.is_file(follow_symlinks=False):
                    continue
                legacy = _LEGACY_INSTALL_TEMP.match(entry.name)
                if entry.name.endswith(INSTALL_TEMP_SUFFIX) or (legacy is not None and legacy["name"] in names):
                    os.unlink(entry.path)
                    removed.append(str(Path(entry.path).relative_to(target)))
        if removed:
            _fsync_directory(directory)
    return sorted(removed)


def _link_new_file(source: Path, destination: Path, expected_sha256: str) -> str:
    """Copy ``source`` to a new ``destination`` (temp + fsync + link); return ``copied`` or ``identical``."""

    if destination.exists() or destination.is_symlink():
        if destination.is_file() and not destination.is_symlink() and sha256_file(destination) == expected_sha256:
            return "identical"
        raise FileExistsError(f"a different file already exists: {destination}")
    destination.parent.mkdir(parents=True, exist_ok=True)
    fd, temporary = tempfile.mkstemp(
        prefix=f".{destination.name}.", suffix=INSTALL_TEMP_SUFFIX, dir=str(destination.parent)
    )
    try:
        with os.fdopen(fd, "wb") as handle, source.open("rb") as reader:
            shutil.copyfileobj(reader, handle, 4 * 1024 * 1024)
            handle.flush()
            _durable_fsync(handle.fileno())
        if sha256_file(Path(temporary)) != expected_sha256:
            raise CohortError(f"copied file differs from its digest: {destination}")
        os.chmod(temporary, 0o444)
        try:
            os.link(temporary, destination)
        except FileExistsError:
            if sha256_file(destination) != expected_sha256:
                raise FileExistsError(f"a different file already exists: {destination}") from None
            return "identical"
    finally:
        with contextlib.suppress(OSError):
            os.unlink(temporary)
    _fsync_directory(destination.parent)
    return "copied"


def install_cohort(source_root: str | Path, cohort_dir: str | Path, *, cohort_id: str) -> dict[str, Any]:
    """Copy a frozen source into a cohort directory and write its definition.

    The source is only read. Every file is created with temp + fsync + hard
    link, so an existing file is never replaced: identical files are skipped
    and a different existing file raises ``FileExistsError``. Files are made
    read-only. The installed directory is then loaded and fully verified.
    Re-running the install is a no-op, except that temp files a killed
    install left behind are removed first (``stale_temporaries_removed``).
    """

    source = Path(source_root)
    target = Path(cohort_dir)
    if target.resolve() == source.resolve():
        raise CohortError("cohort directory must differ from the frozen source")
    digests = IndexDigestCache()
    definition = build_cohort_definition(source, cohort_id=cohort_id, index_digests=digests)
    stale = _remove_install_temporaries(target, [record["path"] for record in definition["files"]])
    counts = {"copied": 0, "identical": 0}
    for record in definition["files"]:
        outcome = _link_new_file(_regular_file(source, record["path"]), target / record["path"], record["sha256"])
        counts[outcome] += 1
    data = _definition_bytes(definition)
    definition_path = target / COHORT_DEFINITION_FILENAME
    if definition_path.exists() or definition_path.is_symlink():
        if definition_path.is_symlink() or definition_path.read_bytes() != data:
            raise FileExistsError(f"a different cohort definition already exists: {definition_path}")
        counts["identical"] += 1
    else:
        fd, temporary = tempfile.mkstemp(
            prefix=f".{COHORT_DEFINITION_FILENAME}.", suffix=INSTALL_TEMP_SUFFIX, dir=str(target)
        )
        try:
            with os.fdopen(fd, "wb") as handle:
                handle.write(data)
                handle.flush()
                _durable_fsync(handle.fileno())
            os.chmod(temporary, 0o444)
            os.link(temporary, definition_path)
        finally:
            with contextlib.suppress(OSError):
                os.unlink(temporary)
        _fsync_directory(target)
        counts["copied"] += 1
    cohort = load_frozen_cohort(target, index_digests=IndexDigestCache())
    return {
        "cohort_id": cohort.cohort_id,
        "definition_sha256": cohort.definition_sha256,
        "files": len(definition["files"]) + 1,
        **counts,
        "stale_temporaries_removed": stale,
    }


# ---------------------------------------------------------------------------
# Data-root maintenance lock
# ---------------------------------------------------------------------------


DataRootBusy = CampaignLockBusy
"""Another process (producer, finalization or maintenance) may be writing the data root."""


@contextlib.contextmanager
def data_root_maintenance_lock(data_root: str | Path, *, purpose: str) -> Iterator[Path]:
    """Hold the campaign writer lock of the data root for a maintenance write.

    This is :func:`metainformant.rna.engine.campaign_lock.acquire_writer` with
    ``kind=purpose``: it refuses while another writer (a producer, a live shell
    launcher lock, the acquisition worker's flock, or another maintenance
    command) holds the root. The kernel releases the flock when the process
    dies, so a crash never leaves a stale lock behind. A command started under
    ``campaign_lock run`` is inside its parent's lock.
    """

    with acquire_writer(data_root, purpose) as lock:
        yield lock.path


# ---------------------------------------------------------------------------
# Producer command-line support (scripts/rna/run_all_species.py, process_species.py)
# ---------------------------------------------------------------------------

DEFAULT_MAX_SAMPLE_GB = 1000.0


def default_max_sample_gb(environ: Mapping[str, str] | None = None) -> float:
    """Per-sample size budget: ``AMALGKIT_MAX_SAMPLE_GB`` (the launcher's), then ``PIPELINE_MAX_GB``, then 1000."""

    values = os.environ if environ is None else environ
    for name in ("AMALGKIT_MAX_SAMPLE_GB", "PIPELINE_MAX_GB"):
        raw = values.get(name, "").strip()
        if raw:
            try:
                value = float(raw)
            except ValueError:
                raise ValueError(f"{name} must be a number of gigabytes: {raw!r}") from None
            if not value > 0:
                raise ValueError(f"{name} must be positive: {raw!r}")
            return value
    return DEFAULT_MAX_SAMPLE_GB


def add_producer_cohort_arguments(parser: argparse.ArgumentParser) -> None:
    """Add ``--cohort``, ``--cohort-definition``, ``--cohort-definition-sha256`` and ``--refresh-metadata``."""

    parser.add_argument(
        "--cohort",
        type=Path,
        help="Frozen cohort directory: build tasks only from it (no NCBI metadata query, no selection rewrite)",
    )
    parser.add_argument("--cohort-definition", type=Path, help="Cohort definition file (default: <cohort>/COHORT.json)")
    parser.add_argument("--cohort-definition-sha256", help="Refuse unless the cohort definition has this sha256")
    parser.add_argument(
        "--refresh-metadata",
        action="store_true",
        help="Allow Amalgkit metadata/select (an NCBI query) when metadata is missing or stale; not with --cohort",
    )


def producer_cohort_from_args(args: argparse.Namespace) -> FrozenCohort | None:
    """Load the verified cohort named on the command line (``None`` without ``--cohort``)."""

    if args.cohort is None:
        if args.cohort_definition is not None or args.cohort_definition_sha256 is not None:
            raise CohortError("--cohort-definition and --cohort-definition-sha256 need --cohort")
        return None
    if args.refresh_metadata:
        raise CohortError("--refresh-metadata cannot be combined with --cohort (a frozen cohort never re-queries NCBI)")
    return load_frozen_cohort(
        args.cohort.expanduser(),
        definition_path=args.cohort_definition,
        expected_definition_sha256=args.cohort_definition_sha256,
    )


# ---------------------------------------------------------------------------
# Command line
# ---------------------------------------------------------------------------


def _load_from_args(args: argparse.Namespace) -> FrozenCohort:
    return load_frozen_cohort(
        args.cohort,
        definition_path=args.cohort_definition,
        expected_definition_sha256=args.cohort_definition_sha256,
    )


def _add_cohort_arguments(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--cohort", type=Path, required=True, help="Cohort directory (frozen inputs)")
    parser.add_argument("--cohort-definition", type=Path, help="Definition file (default: <cohort>/COHORT.json)")
    parser.add_argument("--cohort-definition-sha256", help="Refuse unless the definition has this sha256")


def _read_exclusions_tsv(path: Path) -> list[dict[str, str | None]]:
    rows: list[dict[str, str | None]] = []
    lines = path.read_text(encoding="utf-8").splitlines()
    header = lines[0].split("\t") if lines else []
    run_column = "accession" if "accession" in header else "srr_id"
    if "species" not in header or run_column not in header:
        raise CohortError(f"exclusions table needs species and accession (or srr_id) columns: {path}")
    for line in lines[1:]:
        if not line.strip():
            continue
        values = dict(zip(header, line.split("\t"), strict=False))
        rows.append(
            {
                "species": values["species"].strip(),
                "srr_id": values[run_column].strip(),
                "reason_code": (values.get("reason_code") or PERMANENT_DROP).strip(),
                "reason_detail": (values.get("reason_detail") or "").strip() or None,
            }
        )
    return rows


def main(argv: Sequence[str] | None = None) -> int:
    """``python -m metainformant.rna.engine.cohort {define,install,verify,migrate-db}``."""

    parser = argparse.ArgumentParser(
        prog="python -m metainformant.rna.engine.cohort", description=__doc__.split("\n")[0]
    )
    commands = parser.add_subparsers(dest="command", required=True)

    define = commands.add_parser("define", help="Describe a frozen source as COHORT.json (reads only)")
    define.add_argument("--source", type=Path, required=True)
    define.add_argument("--cohort-id", required=True)
    define.add_argument("--output", type=Path, required=True, help="Where to write the definition")

    install = commands.add_parser("install", help="Copy a frozen source into a cohort directory")
    install.add_argument("--source", type=Path, required=True)
    install.add_argument("--cohort-dir", type=Path, required=True)
    install.add_argument("--cohort-id", required=True)

    verify = commands.add_parser("verify", help="Verify a cohort (and optionally configurations)")
    _add_cohort_arguments(verify)
    verify.add_argument("--config-dir", type=Path, help="Also require these YAMLs to match the cohort")

    migrate = commands.add_parser("migrate-db", help="Bind a progress DB to the cohort (adds rows, never deletes)")
    _add_cohort_arguments(migrate)
    migrate.add_argument("--db", type=Path, required=True)
    migrate.add_argument("--backup", type=Path, help="Backup path (default: <db>.pre-<cohort>.bak)")
    migrate.add_argument("--extra-exclusions", type=Path, help="TSV (species, accession[, reason_code]) to add")
    migrate.add_argument(
        "--reset-transfer-failures",
        action="store_true",
        help="Also reset 'Download Failed (all sources)' rows to pending (default: only environmental failures)",
    )
    args = parser.parse_args(argv)

    try:
        if args.command == "define":
            definition = build_cohort_definition(args.source, cohort_id=args.cohort_id)
            digest = write_cohort_definition(definition, args.output)
            print(json.dumps({"definition": str(args.output), "sha256": digest, "tasks": definition["task_count"]}))
            return 0
        if args.command == "install":
            print(json.dumps(install_cohort(args.source, args.cohort_dir, cohort_id=args.cohort_id), sort_keys=True))
            return 0
        cohort = _load_from_args(args)
        if args.command == "verify":
            if args.config_dir is not None:
                verify_cohort_configs(cohort, args.config_dir)
            print(
                json.dumps(
                    {
                        "cohort_id": cohort.cohort_id,
                        "definition_sha256": cohort.definition_sha256,
                        "species": len(cohort.species),
                        "tasks": len(cohort.tasks),
                        "exclusions": sum(len(v) for v in cohort.exclusions.values()),
                    },
                    sort_keys=True,
                )
            )
            return 0
        from metainformant.rna.engine.progress_db import (
            RETRYABLE_ENVIRONMENT_CLASSES,
            TRANSFER_RESET_CLASSES,
            ProgressDB,
        )

        db_path = args.db.expanduser().absolute()
        if not db_path.is_file():
            raise CohortError(f"progress database does not exist: {db_path}")
        backup = args.backup or db_path.with_name(f"{db_path.name}.pre-{cohort.cohort_id}.bak")
        reset_classes = RETRYABLE_ENVIRONMENT_CLASSES | (
            TRANSFER_RESET_CLASSES if args.reset_transfer_failures else frozenset()
        )
        extra = _read_exclusions_tsv(args.extra_exclusions) if args.extra_exclusions else []
        # The maintenance lock is the producer's own writer lock: a running producer (or any other writer) makes
        # it raise DataRootBusy before anything is read or written.
        with data_root_maintenance_lock(db_path.parent, purpose="migrate-db"), ProgressDB(db_path) as db:
            summary = db.migrate_to_cohort(
                cohort, backup_path=backup, reset_classes=reset_classes, extra_exclusions=extra
            )
        print(json.dumps(summary, indent=2, sort_keys=True))
        return 0
    except (CohortError, DataRootBusy, FileExistsError, OSError) as error:
        print(f"error: {error}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
