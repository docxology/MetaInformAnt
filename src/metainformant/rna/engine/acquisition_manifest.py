"""Hash-bound Amalgkit acquisition snapshots and disjoint task selections."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any

def sha256_file(path: Path) -> str:
    """Hash one file in bounded chunks."""

    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()

def load_snapshot(manifest_path: Path) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    """Load and validate a static manifest and its snapshot metadata."""

    snapshot_path = manifest_path.with_name("snapshot.json")
    snapshot = json.loads(snapshot_path.read_text(encoding="utf-8"))
    if snapshot.get("schema") not in {"metainformant.hymenoptera.gcp_snapshot.v1", "metainformant.rna.acquisition_snapshot.v1"}:
        raise ValueError(f"unsupported snapshot schema: {snapshot.get('schema')!r}")
    if snapshot.get("cloud_launch_policy") not in {"canary_only", "checkpointed"}:
        raise ValueError("snapshot does not declare an allowed cloud launch policy")
    expected_hash = snapshot.get("manifest_sha256")
    if (
        not isinstance(expected_hash, str)
        or sha256_file(manifest_path) != expected_hash
    ):
        raise ValueError("task manifest hash does not match snapshot metadata")
    tasks: list[dict[str, Any]] = []
    task_ids: set[str] = set()
    with manifest_path.open(encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip():
                continue
            task = json.loads(line)
            if task.get("schema") not in {"metainformant.hymenoptera.gcp_task_manifest.v1", "metainformant.rna.acquisition_task.v1"}:
                raise ValueError(f"unsupported task schema on line {line_number}")
            for key in ("species", "config_name", "accession", "batch_index"):
                if key not in task:
                    raise ValueError(f"task line {line_number} lacks {key}")
            species = str(task["species"])
            config_name = str(task["config_name"])
            accession = str(task["accession"])
            task_id = str(task.get("task_id", f"{species}/{accession}"))
            if not species or any(part in species for part in ("/", "\\", "..")):
                raise ValueError(f"unsafe species on line {line_number}: {species!r}")
            if Path(config_name).name != config_name or not config_name.endswith(
                ".yaml"
            ):
                raise ValueError(
                    f"unsafe config name on line {line_number}: {config_name!r}"
                )
            if not accession or any(part in accession for part in ("/", "\\", "..")):
                raise ValueError(
                    f"unsafe accession on line {line_number}: {accession!r}"
                )
            expected_task_id = f"{species}/{accession}"
            if task_id != expected_task_id:
                raise ValueError(
                    f"task_id must equal species/accession on line {line_number}: "
                    f"{task_id!r} != {expected_task_id!r}"
                )
            if task_id in task_ids:
                raise ValueError(f"duplicate task_id on line {line_number}: {task_id}")
            task_ids.add(task_id)
            tasks.append(task)
    if len(tasks) != int(snapshot.get("task_count", -1)):
        raise ValueError("snapshot task count does not match manifest")
    if not tasks:
        raise ValueError("acquisition manifest must not be empty")
    return snapshot, tasks

def load_task_selection(
    selection_path: Path | None,
    snapshot: dict[str, Any],
    tasks: list[dict[str, Any]],
    *,
    snapshot_sha256: str | None = None,
) -> tuple[dict[str, Any] | None, list[dict[str, Any]]]:
    """Load a hash-bound partition sidecar and select its manifest tasks.

    The full manifest remains the immutable input envelope.  A partition is a
    small explicit task-id set, rather than an offset into a mutable queue, so
    retries remain disjoint even when the planner uses non-contiguous,
    base-balanced assignments.
    """

    all_ids = [
        str(task.get("task_id") or f"{task['species']}/{task['accession']}")
        for task in tasks
    ]
    if selection_path is None:
        return None, tasks
    payload = json.loads(
        selection_path.expanduser().resolve().read_text(encoding="utf-8")
    )
    if (
        not isinstance(payload, dict)
        or payload.get("schema") not in {"metainformant.hymenoptera.gcp_partition.v1", "metainformant.rna.acquisition_partition.v1"}
    ):
        raise ValueError("unsupported cloud partition selection schema")
    if payload.get("manifest_sha256") != snapshot.get("manifest_sha256"):
        raise ValueError("partition manifest hash does not match the worker snapshot")
    expected_snapshot = payload.get("snapshot_sha256")
    if expected_snapshot is not None and (
        snapshot_sha256 is None or expected_snapshot != snapshot_sha256
    ):
        raise ValueError("partition snapshot hash does not match the worker snapshot")
    task_ids = payload.get("task_ids")
    if (
        not isinstance(task_ids, list)
        or not task_ids
        or not all(isinstance(item, str) and item for item in task_ids)
    ):
        raise ValueError("partition selection must contain a non-empty task_ids list")
    if len(set(task_ids)) != len(task_ids):
        raise ValueError("partition selection contains duplicate task ids")
    by_id = dict(zip(all_ids, tasks, strict=True))
    unknown = sorted(set(task_ids) - set(by_id))
    if unknown:
        raise ValueError(
            f"partition selection contains tasks absent from manifest: {', '.join(unknown)}"
        )
    if int(payload.get("task_count", -1)) != len(task_ids):
        raise ValueError("partition task_count does not match task_ids")
    selected = [by_id[task_id] for task_id in task_ids]
    return payload, selected

def verify_input_files(snapshot: dict[str, Any], snapshot_dir: Path) -> None:
    """Verify the staged metadata/index files before starting downloads."""

    for record in snapshot.get("input_files", []):
        relative = Path(str(record["path"]))
        if relative.is_absolute() or ".." in relative.parts:
            raise ValueError(f"snapshot input path must be relative: {relative}")
        candidate = snapshot_dir / relative
        if candidate.is_symlink() or not candidate.resolve().is_relative_to(snapshot_dir.resolve()):
            raise ValueError("snapshot input escapes its owned directory")
        if not candidate.is_file():
            raise FileNotFoundError(
                f"staged snapshot input is missing: {record['path']}"
            )
        expected = str(record.get("sha256", ""))
        if expected and sha256_file(candidate) != expected:
            raise ValueError(f"staged input hash mismatch: {candidate}")


def verify_worker_configs(snapshot: dict[str, Any], tasks: list[dict[str, Any]], config_dir: Path) -> None:
    """Bind actual worker YAML bytes to the generic envelope before acquisition."""
    config_hashes = {Path(record["path"]).name: record["sha256"] for record in snapshot.get("input_files", [])
                     if record["path"].startswith("config/amalgkit/")}
    for task in tasks:
        expected = task.get("config_sha256") or config_hashes.get(task["config_name"])
        if snapshot["schema"] == "metainformant.rna.acquisition_snapshot.v1" and not expected:
            raise ValueError("generic acquisition task lacks a frozen configuration checksum")
        path = config_dir / task["config_name"]
        if expected and (path.is_symlink() or sha256_file(path) != expected):
            raise ValueError("worker configuration differs from frozen acquisition envelope")
