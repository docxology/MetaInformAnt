"""Real manifest, immutable-plan, staging and controller-binding controls."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pytest

from metainformant.rna.engine.acquisition_allocation import (
    AcquisitionAllocationError,
    allocate_tasks,
    aws_allocation_ids,
    write_allocation,
)
from metainformant.rna.engine.acquisition_manifest import (
    load_snapshot,
    load_task_selection,
    sha256_file,
)
from metainformant.rna.engine.acquisition_snapshot import (
    create_campaign_manifest,
    stage_worker_inputs,
)


def envelope(tmp_path: Path, count: int = 4) -> Path:
    tmp_path.mkdir(exist_ok=True)
    manifest = tmp_path / "manifest.jsonl"
    tasks = [
        {
            "schema": "metainformant.rna.acquisition_task.v1",
            "task_id": f"ant_a/SRR{i}",
            "species": "ant_a",
            "accession": f"SRR{i}",
            "config_name": "amalgkit_ant_a.yaml",
            "batch_index": i,
        }
        for i in range(1, count + 1)
    ]
    manifest.write_text("".join(json.dumps(t) + "\n" for t in tasks))
    (tmp_path / "snapshot.json").write_text(
        json.dumps(
            {
                "schema": "metainformant.rna.acquisition_snapshot.v1",
                "manifest_sha256": sha256_file(manifest),
                "inventory_sha256": "a" * 64,
                "task_count": count,
                "cloud_launch_policy": "checkpointed",
                "input_files": [],
            }
        )
    )
    return manifest


def test_completed_reserved_and_lane_marginals_are_disjoint(tmp_path: Path) -> None:
    manifest = envelope(tmp_path)
    plan = allocate_tasks(
        manifest,
        backend="hybrid",
        completed=frozenset({"ant_a/SRR1"}),
        reserved=frozenset({"ant_a/SRR2"}),
    )
    assert plan.local_task_ids == ("ant_a/SRR4",)
    assert plan.aws_task_ids == ("ant_a/SRR3",)
    assert (
        sum(
            map(
                len,
                (
                    plan.completed_task_ids,
                    plan.reserved_task_ids,
                    plan.local_task_ids,
                    plan.aws_task_ids,
                ),
            )
        )
        == 4
    )
    path = write_allocation(plan, tmp_path / "plan")
    assert write_allocation(plan, path.parent) == path
    snapshot, tasks = load_snapshot(manifest)
    _, selected = load_task_selection(
        path.parent / "local_partition.json",
        snapshot,
        tasks,
        snapshot_sha256=sha256_file(manifest.with_name("snapshot.json")),
    )
    assert [t["task_id"] for t in selected] == ["ant_a/SRR4"]
    assert aws_allocation_ids(path, frozenset(t["task_id"] for t in tasks), "a" * 64) == frozenset(
        {"ant_a/SRR2", "ant_a/SRR3"}
    )


def test_changed_assignment_cannot_overwrite_existing_plan(tmp_path: Path) -> None:
    manifest = envelope(tmp_path)
    write_allocation(allocate_tasks(manifest, backend="local"), tmp_path / "plan")
    with pytest.raises(FileExistsError):
        write_allocation(allocate_tasks(manifest, backend="aws"), tmp_path / "plan")


def test_foreign_and_empty_task_universes_refused(tmp_path: Path) -> None:
    manifest = envelope(tmp_path)
    with pytest.raises(AcquisitionAllocationError):
        allocate_tasks(manifest, backend="aws", completed=frozenset({"foreign/SRR9"}))
    with pytest.raises(ValueError, match="empty"):
        allocate_tasks(envelope(tmp_path / "empty", 0), backend="local")


def test_controller_refuses_cross_inventory_binding_and_overlap(tmp_path: Path) -> None:
    manifest = envelope(tmp_path)
    path = write_allocation(allocate_tasks(manifest, backend="hybrid"), tmp_path / "plan")
    ids = frozenset(f"ant_a/SRR{i}" for i in range(1, 5))
    with pytest.raises(AcquisitionAllocationError, match="bound"):
        aws_allocation_ids(path, ids, "b" * 64)
    payload = json.loads(path.read_text())
    payload["aws_task_ids"].append(payload["local_task_ids"][0])
    corrupted = tmp_path / "corrupted.json"
    corrupted.write_text(json.dumps(payload))
    with pytest.raises(AcquisitionAllocationError, match="Overlapping"):
        aws_allocation_ids(corrupted, ids, "a" * 64)


def test_generic_inventory_conversion_staging_and_conflict(tmp_path: Path) -> None:
    root = tmp_path / "campaign"
    root.mkdir()
    index = root / "inputs/ant_a/work/index/Ant_a.idx"
    index.parent.mkdir(parents=True)
    index.write_bytes(b"index checksum fixture")
    metadata = root / "inputs/ant_a/work/metadata/metadata_selected.tsv"
    metadata.parent.mkdir()
    metadata.write_text("run\tscientific_name\nSRR1\tAnt a\n")
    config = "1" * 64
    inventory = {
        "species_count": 1,
        "task_count": 1,
        "species": [
            {
                "species": "ant_a",
                "config_name": "amalgkit_ant_a.yaml",
                "config_sha256": config,
                "index_name": "Ant_a.idx",
                "index_sha256": sha256_file(index),
                "metadata_sha256": sha256_file(metadata),
                "tasks": [
                    {
                        "schema": "metainformant.rna.acquisition_task.v1",
                        "task_id": "ant_a/SRR1",
                        "species": "ant_a",
                        "accession": "SRR1",
                        "config_name": "amalgkit_ant_a.yaml",
                        "batch_index": 1,
                    }
                ],
            }
        ],
    }
    (root / "inventory.json").write_text(json.dumps(inventory))
    manifest = create_campaign_manifest(root)
    snapshot, tasks = load_snapshot(manifest)
    assert tasks[0]["config_sha256"] == config
    assert snapshot["inventory_sha256"] == hashlib.sha256((root / "inventory.json").read_bytes()).hexdigest()
    data = tmp_path / "local"
    stage_worker_inputs(manifest, data, frozenset({"ant_a"}))
    stage_worker_inputs(manifest, data, frozenset({"ant_a"}))
    assert (data / "ant_a/work/index/Ant_a.idx").read_bytes() == index.read_bytes()
    (data / "ant_a/work/index/Ant_a.idx").write_bytes(b"changed")
    with pytest.raises(FileExistsError):
        stage_worker_inputs(manifest, data, frozenset({"ant_a"}))


def test_hybrid_local_requires_ack_and_cannot_substitute_aws_selection(tmp_path: Path) -> None:
    from metainformant.rna.engine.acquisition_allocation import verify_local_allocation

    manifest = envelope(tmp_path)
    plan = allocate_tasks(manifest, backend="hybrid")
    allocation = write_allocation(plan, tmp_path / "plan")
    selection = allocation.parent / "local_partition.json"
    ids = frozenset(f"ant_a/SRR{i}" for i in range(1, 5))
    selected = frozenset(plan.local_task_ids)
    with pytest.raises(AcquisitionAllocationError, match="coordinator"):
        verify_local_allocation(selection, manifest, ids, selected, "a" * 64)
    ledger = manifest.parent / "aws_controller.json"
    ledger.write_text(json.dumps({"allocation_sha256": "wrong"}))
    with pytest.raises(AcquisitionAllocationError, match="acknowledged"):
        verify_local_allocation(selection, manifest, ids, selected, "a" * 64)
    ledger.write_text(json.dumps({"allocation_sha256": sha256_file(allocation)}))
    verify_local_allocation(selection, manifest, ids, selected, "a" * 64)
    with pytest.raises(AcquisitionAllocationError, match="coordinated"):
        verify_local_allocation(selection, manifest, ids, frozenset(plan.aws_task_ids), "a" * 64)
    with pytest.raises(AcquisitionAllocationError, match="local allocation"):
        verify_local_allocation(
            allocation.parent / "aws_partition.json", manifest, ids, frozenset(plan.aws_task_ids), "a" * 64
        )
    with pytest.raises(AcquisitionAllocationError, match="explicit"):
        verify_local_allocation(None, manifest, ids, ids, "a" * 64)
