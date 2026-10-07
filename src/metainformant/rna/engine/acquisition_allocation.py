"""Immutable, disjoint local/AWS allocations for one frozen acquisition envelope."""

from __future__ import annotations

import hashlib
import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Literal

from metainformant.rna.engine.acquisition_manifest import load_snapshot, sha256_file
from metainformant.rna.engine.quant_storage import DirectoryStore

Backend = Literal["local", "aws", "hybrid"]


class AcquisitionAllocationError(ValueError):
    """A plan cannot safely assign the frozen acquisition tasks."""

    def __init__(self, reason: str) -> None:
        self.reason = reason
        super().__init__(reason)


@dataclass(frozen=True, slots=True)
class AcquisitionAllocation:
    schema: str
    manifest_sha256: str
    snapshot_sha256: str
    inventory_sha256: str | None
    backend: Backend
    eligible: int
    completed_task_ids: tuple[str, ...]
    reserved_task_ids: tuple[str, ...]
    local_task_ids: tuple[str, ...]
    aws_task_ids: tuple[str, ...]


def allocate_tasks(
    manifest: Path,
    *,
    backend: Backend,
    completed: frozenset[str] = frozenset(),
    reserved: frozenset[str] = frozenset(),
    local_fraction: float = 0.5,
) -> AcquisitionAllocation:
    """Exclude receipts and active owners before assigning deterministic task IDs."""
    snapshot, tasks = load_snapshot(manifest)
    ids = frozenset(task["task_id"] for task in tasks)
    if not ids or not (completed | reserved) <= ids:
        raise AcquisitionAllocationError("Empty manifest or observations outside the frozen task universe")
    if isinstance(local_fraction, bool) or not math.isfinite(local_fraction) or not 0 <= local_fraction <= 1:
        raise AcquisitionAllocationError("local_fraction must be finite and between zero and one")
    remaining = sorted(ids - completed - reserved)
    fractions = {"local": 1.0, "aws": 0.0, "hybrid": local_fraction}
    if backend not in fractions:
        raise AcquisitionAllocationError("Unsupported acquisition backend")
    count = int(len(remaining) * fractions[backend])
    # Spread local assignments throughout species/accession order rather than
    # putting every early species on one lane. This is deterministic count
    # balancing, not an assertion that sample sizes or processing times match.
    local_ids = (
        tuple(
            task for i, task in enumerate(remaining) if (i + 1) * count // len(remaining) > i * count // len(remaining)
        )
        if remaining
        else ()
    )
    local_set = frozenset(local_ids)
    aws_ids = tuple(task for task in remaining if task not in local_set)
    return AcquisitionAllocation(
        "metainformant.rna.acquisition_allocation.v1",
        snapshot["manifest_sha256"],
        sha256_file(manifest.with_name("snapshot.json")),
        snapshot.get("inventory_sha256"),
        backend,
        len(ids),
        tuple(sorted(completed)),
        tuple(sorted(reserved - completed)),
        local_ids,
        aws_ids,
    )


def write_allocation(allocation: AcquisitionAllocation, directory: Path) -> Path:
    """Publish immutable plan and hash-bound lane selections; reject changed reuse."""
    store = DirectoryStore(directory)
    payload = json.dumps(asdict(allocation), sort_keys=True, indent=2, allow_nan=False).encode() + b"\n"
    store.put("allocation.json", payload)
    for lane, ids in (("local", allocation.local_task_ids), ("aws", allocation.aws_task_ids)):
        if not ids:
            continue
        partition = {
            "schema": "metainformant.rna.acquisition_partition.v1",
            "manifest_sha256": allocation.manifest_sha256,
            "snapshot_sha256": allocation.snapshot_sha256,
            "task_count": len(ids),
            "task_ids": ids,
            "allocation_sha256": hashlib.sha256(payload).hexdigest(),
            "backend": lane,
        }
        store.put(f"{lane}_partition.json", json.dumps(partition, sort_keys=True, indent=2).encode() + b"\n")
    return directory / "allocation.json"


def aws_allocation_ids(path: Path, inventory_ids: frozenset[str], inventory_sha256: str | None) -> frozenset[str]:
    """Validate a controller allocation before any launch or source publication."""
    payload = json.loads(path.read_text())
    if not isinstance(payload, dict) or payload.get("schema") != "metainformant.rna.acquisition_allocation.v1":
        raise AcquisitionAllocationError("Unsupported acquisition allocation")
    if payload.get("inventory_sha256") != inventory_sha256:
        raise AcquisitionAllocationError("Allocation is not bound to the controller's frozen inventory bytes")
    groups = []
    for field in ("completed_task_ids", "reserved_task_ids", "local_task_ids", "aws_task_ids"):
        ids = payload.get(field)
        if not isinstance(ids, list) or not all(isinstance(x, str) for x in ids) or len(set(ids)) != len(ids):
            raise AcquisitionAllocationError(f"Malformed allocation field: {field}")
        groups.append(frozenset(ids))
    union: set[str] = set()
    for group in groups:
        if union & group:
            raise AcquisitionAllocationError("Overlapping acquisition allocation")
        union.update(group)
    if union != inventory_ids or payload.get("eligible") != len(inventory_ids):
        raise AcquisitionAllocationError("Allocation does not cover the complete frozen inventory")
    # Reserved tasks came from the existing AWS owners. Preserve their retry
    # responsibility after termination; occupied-task guards prevent overlap.
    return groups[-1] | groups[1]


def verify_local_allocation(
    selection_path: Path | None,
    manifest: Path,
    all_ids: frozenset[str],
    selected_ids: frozenset[str],
    inventory_sha256: str | None,
) -> None:
    """Require the AWS coordinator's allocation acknowledgment before a hybrid local lane.

    An unactivated proposal must not start local work while an old unrestricted
    AWS controller can still claim the same tasks. Legacy GCP selections retain
    their existing isolated-worker contract.
    """
    ledger_path = manifest.parent / "aws_controller.json"
    if selection_path is None:
        if ledger_path.exists():
            raise AcquisitionAllocationError("a cloud campaign requires an explicit disjoint local allocation")
        return
    selection = json.loads(selection_path.read_text())
    if selection.get("schema") != "metainformant.rna.acquisition_partition.v1":
        return
    allocation_path = selection_path.parent / "allocation.json"
    encoded = allocation_path.read_bytes()
    digest = hashlib.sha256(encoded).hexdigest()
    if selection.get("allocation_sha256") != digest or selection.get("backend") != "local":
        raise AcquisitionAllocationError("local selection is not bound to its immutable local allocation")
    aws_allocation_ids(allocation_path, all_ids, inventory_sha256)
    allocation = json.loads(encoded)
    if frozenset(allocation["local_task_ids"]) != selected_ids:
        raise AcquisitionAllocationError("local selection differs from the coordinated task assignment")
    if allocation["aws_task_ids"] or allocation["reserved_task_ids"] or ledger_path.exists():
        if not ledger_path.exists():
            raise AcquisitionAllocationError("start the allocated AWS coordinator before the hybrid local worker")
        ledger = json.loads(ledger_path.read_text())
        if ledger.get("allocation_sha256") != digest:
            raise AcquisitionAllocationError("AWS coordinator has not acknowledged this disjoint local allocation")
