"""Hash-bound frozen acquisition snapshots: ``snapshot.json`` plus its ``manifest.jsonl`` task list."""

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
    if snapshot.get("schema") != "metainformant.rna.acquisition_snapshot.v1":
        raise ValueError(f"unsupported snapshot schema: {snapshot.get('schema')!r}")
    # ``cloud_launch_policy`` is a field of the frozen snapshot format; every
    # published snapshot, the installed cohort's included, records "checkpointed".
    if snapshot.get("cloud_launch_policy") != "checkpointed":
        raise ValueError("snapshot does not declare the checkpointed launch policy")
    expected_hash = snapshot.get("manifest_sha256")
    if not isinstance(expected_hash, str) or sha256_file(manifest_path) != expected_hash:
        raise ValueError("task manifest hash does not match snapshot metadata")
    tasks: list[dict[str, Any]] = []
    task_ids: set[str] = set()
    with manifest_path.open(encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip():
                continue
            task = json.loads(line)
            # The frozen Hymenoptera cohort's tasks carry the
            # ``metainformant.hymenoptera.gcp_task_manifest.v1`` schema identifier.
            if task.get("schema") not in {
                "metainformant.hymenoptera.gcp_task_manifest.v1",
                "metainformant.rna.acquisition_task.v1",
            }:
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
            if Path(config_name).name != config_name or not config_name.endswith(".yaml"):
                raise ValueError(f"unsafe config name on line {line_number}: {config_name!r}")
            if not accession or any(part in accession for part in ("/", "\\", "..")):
                raise ValueError(f"unsafe accession on line {line_number}: {accession!r}")
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
