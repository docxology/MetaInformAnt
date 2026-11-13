"""Export verified, reusable quantification outputs to a portable, versionable snapshot.

A snapshot is ``<dest>/<species>/<run>/`` holding a sample's quantification
payload and its witnesses, plus a deterministic inventory under
``<inventory-dir>``. Each sample is classified with the content-addressed reuse
contract (:func:`metainformant.rna.engine.provenance.build_species_reuse_context`,
``verify_content=True``) against the frozen cohort's configuration and kallisto
index digests; progress-table state is never evidence. Only ``current`` and
``version_drift_compatible`` outputs of non-excluded cohort tasks are exported,
and only an allowlist of files (abundance table and HDF5, run info, provenance
sidecar, reference-binding and restored-input witnesses) is copied. Kallisto
indexes, raw reads and logs are never part of a snapshot.

Idempotent and conservative:

* A sample whose destination already holds identical bytes is ``present`` and is
  not copied again. A destination that differs is a ``conflict`` and is never
  overwritten (exit 65). A destination that matches but lacks witness files the
  source now has (a reference binding, restored inputs) is ``extend``: only the
  missing witness files are added.
* ``--git-tree REPO:TREEISH`` takes presence from a Git tree as well as the
  working tree (blob ids; files matching ``--lfs-glob`` compare as Git LFS
  pointers), so a sparse checkout that does not materialize the snapshot can
  still be compared and extended. Copied but uncommitted samples are ``staged``.
* New samples are copied to a hidden staging directory, re-hashed, and
  published with one rename; added witness files are published one rename each.
  Source files are only read. Leftover staging from an interrupted run is
  removed at the start of the next ``--execute``.
* The inventory is per species, so hosts that export different species never
  touch the same files: ``<species>.tsv`` lists every file (size, sha256, Git
  blob id) verified at a data root, and ``<species>.tasks.tsv`` each cohort
  task's state (``present`` when exported here or already in the Git tree).
  Rows of runs the Git tree holds but this data root lacks are carried forward.
  No timestamps: re-running on unchanged inputs rewrites identical bytes.

The default is a read-only plan printed as JSON; ``--execute`` copies and
writes the inventory. Exit codes: 0 consistent, 65 conflicts need an owner
decision, 2 usage or cohort errors.
"""

from __future__ import annotations

import argparse
import fnmatch
import hashlib
import json
import os
import shutil
import subprocess
from collections import Counter
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Iterable, Mapping

from metainformant.rna.engine.cohort import CohortError, FrozenCohort, load_frozen_cohort
from metainformant.rna.engine.provenance import (
    QUANT_PROVENANCE_FILENAME,
    REFERENCE_BINDING_FILENAME,
    RESTORED_INPUTS_DIR,
    REUSABLE_QUANT_STATUSES,
    SpeciesReuseContext,
    build_species_reuse_context,
)
from metainformant.rna.engine.species import configured_data_root
from metainformant.rna.engine.streaming_orchestrator import _species_work_dir

SCHEMA = "metainformant.rna.quant_snapshot.v1"
EXIT_OK = 0
EXIT_CONFLICT = 65
EXIT_USAGE = 2
LFS_POINTER_TEMPLATE = "version https://git-lfs.github.com/spec/v1\noid sha256:{sha256}\nsize {size}\n"
RESTORED_INPUT_FILES = ("config.yaml", "manifest.json", "reference_manifest.json")
WITNESS_RELPATHS = frozenset(
    {REFERENCE_BINDING_FILENAME, *(f"{RESTORED_INPUTS_DIR}/{name}" for name in RESTORED_INPUT_FILES)}
)
EXPORTED_STATES = frozenset({"present", "staged", "new", "extend"})
INVENTORY_COLUMNS = ("run", "status", "file", "size", "sha256", "git_oid")
TASK_COLUMNS = ("run", "cohort", "snapshot")
STAGING_SUFFIX = ".snapshot-staging"
TEMPORARY_SUFFIX = ".snapshot-tmp"
_CHUNK = 1 << 20


def payload_relpaths(run: str) -> tuple[str, ...]:
    """The files a sample snapshot may carry, relative to its directory (an allowlist)."""

    return (
        f"{run}_abundance.tsv",
        f"{run}_abundance.h5",
        f"{run}_run_info.json",
        QUANT_PROVENANCE_FILENAME,
        *sorted(WITNESS_RELPATHS),
    )


def required_relpaths(run: str) -> tuple[str, ...]:
    return (f"{run}_abundance.tsv", f"{run}_run_info.json", QUANT_PROVENANCE_FILENAME)


@dataclass(frozen=True)
class FileDigest:
    relpath: str
    size: int
    sha256: str
    git_oid: str


def lfs_pointer(sha256: str, size: int) -> bytes:
    """The Git LFS pointer file that stands for content with this digest and size."""

    return LFS_POINTER_TEMPLATE.format(sha256=sha256, size=size).encode()


def git_blob_oid(data: bytes) -> str:
    """The Git blob id of ``data``."""

    return hashlib.sha1(f"blob {len(data)}\0".encode() + data, usedforsecurity=False).hexdigest()


def _is_lfs(relpath: str, lfs_globs: tuple[str, ...]) -> bool:
    name = relpath.rsplit("/", 1)[-1]
    return any(fnmatch.fnmatch(name, pattern) for pattern in lfs_globs)


def digest_file(path: Path, relpath: str, lfs_globs: tuple[str, ...]) -> FileDigest:
    """SHA-256 of a file and the blob id it has in a Git tree (its LFS pointer's id for LFS files)."""

    size = path.stat().st_size
    sha256 = hashlib.sha256()
    blob = hashlib.sha1(f"blob {size}\0".encode(), usedforsecurity=False)
    read = 0
    with path.open("rb") as handle:
        while chunk := handle.read(_CHUNK):
            sha256.update(chunk)
            blob.update(chunk)
            read += len(chunk)
    if read != size:
        raise OSError(f"{path} changed size while it was hashed ({size} -> {read} bytes)")
    content_sha256 = sha256.hexdigest()
    oid = git_blob_oid(lfs_pointer(content_sha256, size)) if _is_lfs(relpath, lfs_globs) else blob.hexdigest()
    return FileDigest(relpath=relpath, size=size, sha256=content_sha256, git_oid=oid)


@dataclass
class SampleExport:
    species: str
    run: str
    source: Path
    status: str = "missing"  # reuse classification of the source output
    files: list[FileDigest] = field(default_factory=list)
    snapshot: str = ""  # present | staged | new | extend | conflict | not_reusable | missing_payload
    reason: str = ""
    to_copy: list[FileDigest] = field(default_factory=list)


def read_git_tree(repo: Path, treeish: str, prefix: str) -> dict[tuple[str, str], dict[str, str]]:
    """Blob ids under ``prefix`` of a Git tree, grouped by ``(species, run)``.

    Only tree metadata is read (``git ls-tree``), so a partial or sparse clone
    fetches no file contents.
    """

    prefix = prefix.strip("/")
    result = subprocess.run(
        ["git", "-C", str(repo), "ls-tree", "-r", "-z", "--full-tree", treeish, "--", prefix],
        check=False,
        capture_output=True,
    )
    if result.returncode != 0:
        raise SystemExit(f"git ls-tree {treeish} in {repo} failed: {result.stderr.decode(errors='replace').strip()}")
    grouped: dict[tuple[str, str], dict[str, str]] = {}
    base = f"{prefix}/"
    for record in result.stdout.split(b"\0"):
        meta, _, raw_path = record.partition(b"\t")
        fields = meta.split()
        if len(fields) != 3 or fields[1] != b"blob":
            continue
        path = raw_path.decode()
        if not path.startswith(base):
            continue
        parts = path[len(base) :].split("/", 2)
        if len(parts) == 3:
            grouped.setdefault((parts[0], parts[1]), {})[parts[2]] = fields[2].decode()
    return grouped


def working_tree_oids(sample_dir: Path, run: str, lfs_globs: tuple[str, ...]) -> dict[str, str]:
    """Blob ids of the allowlisted files present in a destination sample directory."""

    found: dict[str, str] = {}
    for relpath in payload_relpaths(run):
        path = sample_dir / relpath
        if path.is_file() and not path.is_symlink():
            found[relpath] = digest_file(path, relpath, lfs_globs).git_oid
    return found


def compare_destination(
    expected: list[FileDigest],
    committed: Mapping[str, str],
    on_disk: Mapping[str, str],
    *,
    tree_mode: bool,
) -> tuple[str, str, list[FileDigest]]:
    """``(state, reason, files to copy)`` of one sample against its destination.

    ``committed`` and ``on_disk`` map relative paths to Git blob ids, from the
    tree and from the working tree. Existing destination bytes are never
    replaced: any difference is a conflict, and only witness files may be added
    next to an otherwise identical payload.
    """

    for relpath, oid in on_disk.items():
        if relpath in committed and committed[relpath] != oid:
            return "conflict", f"working tree differs from the tree: {relpath}", []
    combined = {**on_disk, **committed}
    if not combined:
        return "new", "", list(expected)
    changed = [d.relpath for d in expected if d.relpath in combined and combined[d.relpath] != d.git_oid]
    if changed:
        return "conflict", f"different bytes: {', '.join(changed)}", []
    missing = [d for d in expected if d.relpath not in combined]
    if missing:
        names = [d.relpath for d in missing]
        if all(name in WITNESS_RELPATHS for name in names):
            return "extend", f"adds witnesses: {', '.join(names)}", missing
        return "conflict", f"destination lacks payload files: {', '.join(names)}", []
    if tree_mode and any(d.relpath not in committed for d in expected):
        return "staged", "", []
    return "present", "", []


def classify_sample(
    context: SpeciesReuseContext,
    species: str,
    run: str,
    sample_dir: Path,
    lfs_globs: tuple[str, ...],
) -> SampleExport:
    """Classify a source output and, when reusable, hash its allowlisted files."""

    item = SampleExport(species=species, run=run, source=sample_dir)
    try:
        verdict = context.classify(sample_dir, run, verify_content=True)
    except OSError as error:
        item.status, item.snapshot, item.reason = "unreadable", "not_reusable", str(error)
        return item
    item.status = str(verdict.get("status"))
    if item.status not in REUSABLE_QUANT_STATUSES:
        item.snapshot, item.reason = "not_reusable", str(verdict.get("reason") or "")
        return item
    missing = [relpath for relpath in required_relpaths(run) if not (sample_dir / relpath).is_file()]
    if missing:
        item.snapshot, item.reason = "missing_payload", f"missing: {', '.join(missing)}"
        return item
    try:
        item.files = [
            digest_file(sample_dir / relpath, relpath, lfs_globs)
            for relpath in payload_relpaths(run)
            if (sample_dir / relpath).is_file() and not (sample_dir / relpath).is_symlink()
        ]
    except OSError as error:
        item.status, item.snapshot, item.reason = "unreadable", "not_reusable", str(error)
    return item


def _fsync_path(path: Path) -> None:
    descriptor = os.open(path, os.O_RDONLY)
    try:
        os.fsync(descriptor)
    finally:
        os.close(descriptor)


def _copy_verified(source: Path, destination: Path, digest: FileDigest, lfs_globs: tuple[str, ...]) -> None:
    shutil.copyfile(source, destination)
    if digest_file(destination, digest.relpath, lfs_globs) != digest:
        raise OSError(f"{source} changed while it was exported")
    _fsync_path(destination)


def remove_stale_staging(dest: Path) -> list[str]:
    """Remove staging left by an interrupted ``--execute`` (only this module's own names)."""

    removed: list[str] = []
    if not dest.is_dir():
        return removed
    for path in sorted(dest.glob(f"*/.*{STAGING_SUFFIX}")):
        shutil.rmtree(path)
        removed.append(str(path.relative_to(dest)))
    for pattern in (f"*/*/.*{TEMPORARY_SUFFIX}", f"*/*/{RESTORED_INPUTS_DIR}/.*{TEMPORARY_SUFFIX}"):
        for path in sorted(dest.glob(pattern)):
            path.unlink()
            removed.append(str(path.relative_to(dest)))
    return removed


def publish_sample(item: SampleExport, dest: Path, lfs_globs: tuple[str, ...]) -> None:
    """Publish a ``new`` sample with one rename, or add an ``extend`` sample's witness files."""

    species_dir = dest / item.species
    species_dir.mkdir(parents=True, exist_ok=True)
    target = species_dir / item.run
    if item.snapshot == "extend":
        for digest in item.to_copy:
            final = target / digest.relpath
            final.parent.mkdir(parents=True, exist_ok=True)
            temporary = final.with_name(f".{final.name}{TEMPORARY_SUFFIX}")
            try:
                _copy_verified(item.source / digest.relpath, temporary, digest, lfs_globs)
                if final.exists():
                    raise FileExistsError(f"{final} appeared while it was staged")
                os.rename(temporary, final)
            finally:
                temporary.unlink(missing_ok=True)
            _fsync_path(final.parent)
        return
    if item.snapshot != "new":
        raise ValueError(f"nothing to publish for a {item.snapshot!r} sample")
    if target.exists():
        if any(target.iterdir()):
            raise FileExistsError(f"{target} holds unrecognized files; refusing to publish into it")
        target.rmdir()
    staging = species_dir / f".{item.run}{STAGING_SUFFIX}"
    if staging.exists():
        shutil.rmtree(staging)
    try:
        staging.mkdir()
        for digest in item.to_copy:
            copied = staging / digest.relpath
            copied.parent.mkdir(parents=True, exist_ok=True)
            _copy_verified(item.source / digest.relpath, copied, digest, lfs_globs)
        for directory in sorted({p.parent for p in staging.rglob("*") if p.is_file()}, reverse=True):
            _fsync_path(directory)
        os.rename(staging, target)
        _fsync_path(species_dir)
    except BaseException:
        shutil.rmtree(staging, ignore_errors=True)
        raise


def _write_text_if_changed(path: Path, text: str) -> bool:
    path.parent.mkdir(parents=True, exist_ok=True)
    data = text.encode()
    if path.is_file() and path.read_bytes() == data:
        return False
    temporary = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    try:
        temporary.write_bytes(data)
        os.replace(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)
    return True


def _digest_rows(item: SampleExport) -> list[str]:
    return [
        "\t".join((item.run, item.status, digest.relpath, str(digest.size), digest.sha256, digest.git_oid))
        for digest in sorted(item.files, key=lambda entry: entry.relpath)
    ]


def _previous_rows(path: Path) -> dict[str, list[str]]:
    """Rows of an existing species inventory, grouped by run (empty when there is none)."""

    grouped: dict[str, list[str]] = {}
    if not path.is_file():
        return grouped
    lines = path.read_text(encoding="utf-8").splitlines()
    if not lines or lines[0] != "\t".join(INVENTORY_COLUMNS):
        raise ValueError(f"{path} is not a {SCHEMA} species inventory")
    for line in lines[1:]:
        grouped.setdefault(line.split("\t", 1)[0], []).append(line)
    return grouped


def task_state(
    cohort: FrozenCohort,
    items: Mapping[tuple[str, str], SampleExport],
    tree: Mapping[tuple[str, str], Mapping[str, str]],
    species: str,
    run: str,
) -> str:
    """The snapshot state of one cohort task, as recorded in ``<species>.tasks.tsv``.

    ``present`` when this data root exports it or the Git tree already holds
    it (another host may have produced it); otherwise the reason it is not in
    the snapshot, or ``pending`` when this data root has no output for it.
    """

    if cohort.is_excluded(species, run):
        return "excluded"
    item = items.get((species, run))
    if item is not None and (item.snapshot in EXPORTED_STATES or item.snapshot == "conflict"):
        return "present" if item.snapshot != "conflict" else "conflict"
    if (species, run) in tree:
        return "present"
    return "pending" if item is None else item.snapshot


def write_inventory(
    inventory_dir: Path,
    cohort: FrozenCohort,
    items: Mapping[tuple[str, str], SampleExport],
    tree: Mapping[tuple[str, str], Mapping[str, str]],
    species_names: Iterable[str],
) -> list[str]:
    """Write ``<species>.tsv`` (files) and ``<species>.tasks.tsv`` (task states) per species.

    Rows of runs the Git tree holds but this data root does not export are
    carried forward from the previous inventory, so a host that only produced
    part of a species never drops another host's records.
    """

    changed: list[str] = []
    for species in sorted(species_names):
        inventory = inventory_dir / f"{species}.tsv"
        previous = _previous_rows(inventory)
        exported = {
            run: item for (owner, run), item in items.items() if owner == species and item.snapshot in EXPORTED_STATES
        }
        committed = {run for owner, run in tree if owner == species}
        rows = ["\t".join(INVENTORY_COLUMNS)]
        for run in sorted(set(exported) | (committed & set(previous))):
            rows.extend(_digest_rows(exported[run]) if run in exported else previous[run])
        if _write_text_if_changed(inventory, "\n".join(rows) + "\n"):
            changed.append(inventory.name)
        task_rows = ["\t".join(TASK_COLUMNS)]
        for task in sorted(cohort.tasks_for(species), key=lambda entry: entry.accession):
            cohort_state = "excluded" if cohort.is_excluded(species, task.accession) else "task"
            task_rows.append(
                "\t".join((task.accession, cohort_state, task_state(cohort, items, tree, species, task.accession)))
            )
        tasks = inventory_dir / f"{species}.tasks.tsv"
        if _write_text_if_changed(tasks, "\n".join(task_rows) + "\n"):
            changed.append(tasks.name)
    return changed


def build_summary(
    cohort: FrozenCohort,
    items: Mapping[tuple[str, str], SampleExport],
    tree: Mapping[tuple[str, str], Mapping[str, str]],
    species_names: Iterable[str],
) -> dict[str, Any]:
    """Counts of the snapshot after this export (task states as ``<species>.tasks.tsv`` records them)."""

    by_species: dict[str, dict[str, int]] = {}
    totals: Counter[str] = Counter()
    for species in species_names:
        counts: Counter[str] = Counter(
            task_state(cohort, items, tree, species, task.accession) for task in cohort.tasks_for(species)
        )
        totals.update(counts)
        by_species[species] = dict(sorted(counts.items()))
    exported = [item for item in items.values() if item.snapshot in EXPORTED_STATES]
    return {
        "schema": SCHEMA,
        "cohort_id": cohort.cohort_id,
        "cohort_definition_sha256": cohort.definition_sha256,
        "task_states": dict(sorted(totals.items())),
        "by_species": by_species,
        "exported_samples": len(exported),
        "exported_files": sum(len(item.files) for item in exported),
        "exported_bytes": sum(digest.size for item in exported for digest in item.files),
        "reuse_statuses": dict(sorted(Counter(item.status for item in exported).items())),
    }


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--data-root", type=Path, help="Amalgkit data root (default: AMALGKIT_DATA_ROOT)")
    parser.add_argument("--config-dir", type=Path, required=True, help="Species YAML directory bound by the cohort")
    parser.add_argument("--cohort", type=Path, required=True, help="Installed frozen cohort directory")
    parser.add_argument("--cohort-definition-sha256", help="Refuse unless COHORT.json has this sha256")
    parser.add_argument("--dest", type=Path, required=True, help="Snapshot directory: <dest>/<species>/<run>/")
    parser.add_argument("--inventory-dir", type=Path, help="Write per-species inventories here on --execute")
    parser.add_argument(
        "--git-tree",
        metavar="REPO:TREEISH",
        help="Also take presence from this Git tree; --dest must lie inside REPO's working tree",
    )
    parser.add_argument(
        "--lfs-glob",
        action="append",
        help="File-name glob stored as Git LFS pointers in --git-tree (repeatable; default *.h5)",
    )
    parser.add_argument("--species", action="append", help="Restrict to these species")
    parser.add_argument("--workers", type=int, default=8, help="Parallel sample classifications (default 8)")
    parser.add_argument("--execute", action="store_true", help="Copy new samples and write the inventory")
    return parser


def main(argv: list[str] | None = None) -> int:
    """Run the snapshot CLI; ``--data-root`` is exported only for the duration of the call."""

    args = _parser().parse_args(argv)
    previous = os.environ.get("AMALGKIT_DATA_ROOT")
    try:
        return _run(args)
    except CohortError as error:
        print(f"cohort error: {error}")
        return EXIT_USAGE
    finally:
        if previous is None:
            os.environ.pop("AMALGKIT_DATA_ROOT", None)
        else:
            os.environ["AMALGKIT_DATA_ROOT"] = previous


def _run(args: argparse.Namespace) -> int:
    if args.data_root:
        os.environ["AMALGKIT_DATA_ROOT"] = str(args.data_root.expanduser().resolve())
    data_root = configured_data_root().expanduser().resolve()
    config_dir = args.config_dir.expanduser().resolve()
    dest = args.dest.expanduser().resolve()
    if dest == data_root or data_root in dest.parents:
        raise SystemExit(f"--dest {dest} lies inside the data root; refusing to write into it")
    cohort = load_frozen_cohort(args.cohort.expanduser(), expected_definition_sha256=args.cohort_definition_sha256)
    lfs_globs = tuple(args.lfs_glob or ["*.h5"])
    requested = set(args.species or [])
    unknown = requested - set(cohort.species_names)
    if unknown:
        raise SystemExit(f"unknown species: {', '.join(sorted(unknown))}")
    species_names = [name for name in cohort.species_names if not requested or name in requested]

    tree: dict[tuple[str, str], dict[str, str]] = {}
    tree_mode = bool(args.git_tree)
    if args.git_tree:
        repo_text, _, treeish = args.git_tree.rpartition(":")
        if not repo_text or not treeish:
            raise SystemExit("--git-tree needs REPO:TREEISH")
        repo = Path(repo_text).expanduser().resolve()
        try:
            prefix = dest.relative_to(repo).as_posix()
        except ValueError:
            raise SystemExit(f"--dest {dest} is not inside {repo}") from None
        tree = read_git_tree(repo, treeish, prefix)

    def examine(species: str, context: SpeciesReuseContext, run: str, quant_root: Path) -> SampleExport:
        item = classify_sample(context, species, run, quant_root / run, lfs_globs)
        if item.snapshot:
            return item
        on_disk = working_tree_oids(dest / species / run, run, lfs_globs)
        item.snapshot, item.reason, item.to_copy = compare_destination(
            item.files, tree.get((species, run), {}), on_disk, tree_mode=tree_mode
        )
        return item

    items: dict[tuple[str, str], SampleExport] = {}
    for species in species_names:
        entry = cohort.entry(species)
        context = build_species_reuse_context(
            species,
            config_dir / entry.config_name,
            cohort_id=cohort.cohort_id,
            reference_index_sha256=entry.index_sha256,
            reference_index_name=entry.index_name,
        )
        if context.config_sha256 != entry.config_sha256:
            raise SystemExit(f"{config_dir / entry.config_name} does not hash to the cohort digest for {species}")
        quant_root = _species_work_dir(species) / "quant"
        runs = [
            task.accession
            for task in cohort.tasks_for(species)
            if not cohort.is_excluded(species, task.accession) and (quant_root / task.accession).is_dir()
        ]
        with ThreadPoolExecutor(max_workers=max(1, args.workers)) as pool:
            for item in pool.map(lambda run: examine(species, context, run, quant_root), runs):
                items[(species, item.run)] = item

    to_publish = sorted(
        (item for item in items.values() if item.snapshot in {"new", "extend"}),
        key=lambda entry: (entry.species, entry.run),
    )
    summary = build_summary(cohort, items, tree, species_names)
    report: dict[str, Any] = {
        **summary,
        "mode": "execute" if args.execute else "plan",
        "dest": str(dest),
        "snapshot_states": dict(sorted(Counter(item.snapshot for item in items.values()).items())),
        "to_publish": {
            "new": sum(1 for item in to_publish if item.snapshot == "new"),
            "extend": sum(1 for item in to_publish if item.snapshot == "extend"),
            "bytes": sum(digest.size for item in to_publish for digest in item.to_copy),
        },
        "conflicts": [
            {"species": item.species, "run": item.run, "reason": item.reason}
            for item in sorted(items.values(), key=lambda entry: (entry.species, entry.run))
            if item.snapshot == "conflict"
        ],
        "not_exported": [
            {"species": item.species, "run": item.run, "status": item.status, "reason": item.reason}
            for item in sorted(items.values(), key=lambda entry: (entry.species, entry.run))
            if item.snapshot in {"not_reusable", "missing_payload"}
        ],
    }
    if args.execute:
        report["removed_stale_staging"] = remove_stale_staging(dest)
        for item in to_publish:
            publish_sample(item, dest, lfs_globs)
        report["published"] = len(to_publish)
        if args.inventory_dir is not None:
            report["inventory_changed"] = write_inventory(
                args.inventory_dir.expanduser().resolve(), cohort, items, tree, species_names
            )
    print(json.dumps(report, indent=2, sort_keys=True))
    return EXIT_CONFLICT if report["conflicts"] else EXIT_OK
