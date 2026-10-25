from __future__ import annotations

import csv
import gzip
import json
import os
import subprocess
import sys
from pathlib import Path
from typing import TypedDict

import pytest

from metainformant.rna.engine.acquisition_allocation import (
    allocate_tasks,
    write_allocation,
)
from metainformant.rna.engine.acquisition_manifest import sha256_file
from metainformant.rna.engine.acquisition_order import plan_acquisition_order
from metainformant.rna.engine.acquisition_worker import run_manifest

MIB = 1024**2
GIB = 1024**3


class SizeTask(TypedDict):
    task_id: str
    species: str
    fastq_bytes: int


def sample(species: str, name: str, size: int) -> SizeTask:
    return {"task_id": f"{species}/{name}", "species": species, "fastq_bytes": size}


def test_global_bands_clear_small_work_across_species_before_large_or_unknown() -> None:
    pending = [
        sample("ant", "small", 90 * MIB),
        sample("ant", "medium", 200 * MIB),
        sample("bee", "small", 70 * MIB),
        sample("bee", "large", 25 * GIB),
        sample("wasp", "unknown", 0),
    ]
    before = json.dumps(pending, sort_keys=True)
    result = plan_acquisition_order(iter(pending), iter(["sealed/b", "sealed/a"]), global_size_bands=True)
    assert result == [
        "sealed/a",
        "ant/small",
        "sealed/b",
        "bee/small",
        "ant/medium",
        "bee/large",
        "wasp/unknown",
    ]
    assert json.dumps(pending, sort_keys=True) == before
    assert plan_acquisition_order(pending, []) == [
        "ant/small",
        "bee/small",
        "wasp/unknown",
        "ant/medium",
        "bee/large",
    ]


def test_size_boundaries_are_inclusive_and_species_rotate_inside_each_band() -> None:
    pending = [
        sample("ant", "upper", 100 * MIB),
        sample("ant", "next", 101 * MIB),
        sample("bee", "lower", MIB),
        sample("bee", "next", 500 * MIB),
        sample("ant", "one_gib", GIB),
        sample("bee", "two_gib", 2 * GIB),
        sample("ant", "five_gib", 5 * GIB),
        sample("bee", "twenty_gib", 20 * GIB),
        sample("ant", "huge", 21 * GIB),
    ]
    assert plan_acquisition_order(pending, [], global_size_bands=True) == [
        "ant/upper",
        "bee/lower",
        "ant/next",
        "bee/next",
        "ant/one_gib",
        "bee/two_gib",
        "ant/five_gib",
        "bee/twenty_gib",
        "ant/huge",
    ]


@pytest.mark.parametrize(
    "pending,checkpoints",
    [
        ([], ["duplicate", "duplicate"]),
        ([sample("ant", "x", 1)], ["ant/x"]),
        ([sample("ant", "x", 1), sample("ant", "x", 10 * GIB)], []),
    ],
)
def test_global_bands_refuse_duplicate_identities(pending: list[SizeTask], checkpoints: list[str]) -> None:
    with pytest.raises(ValueError, match="Duplicate acquisition task"):
        plan_acquisition_order(pending, checkpoints, global_size_bands=True)


def test_empty_lanes_and_checkpoint_tail_preserve_every_identity() -> None:
    assert plan_acquisition_order([], [], global_size_bands=True) == []
    assert plan_acquisition_order([], ["c", "a", "b"], global_size_bands=True) == [
        "a",
        "b",
        "c",
    ]
    assert plan_acquisition_order([sample("ant", "x", 1)], ["c", "a", "b"], global_size_bands=True) == [
        "a",
        "ant/x",
        "b",
        "c",
    ]


def test_full_cohort_sized_plan_is_deterministic_and_preserves_all_tasks() -> None:
    pending = [sample(f"species_{number % 27:02d}", str(number), (number % 31) * 900 * MIB) for number in range(18200)]
    first = plan_acquisition_order(pending, [], global_size_bands=True)
    assert len(first) == len(set(first)) == 18200
    assert set(first) == {task["task_id"] for task in pending}
    assert first == plan_acquisition_order(reversed(pending), [], global_size_bands=True)
    sizes = {task["task_id"]: task["fastq_bytes"] for task in pending}
    unknown = [task_id for task_id in first if sizes[task_id] == 0]
    assert first[-len(unknown) :] == unknown


def test_malformed_size_is_refused_instead_of_silently_treated_as_small() -> None:
    malformed = json.loads('[{"task_id":"ant/x","species":"ant","fastq_bytes":"broken"}]')
    with pytest.raises(ValueError):
        plan_acquisition_order(malformed, [], global_size_bands=True)


def test_continuous_cli_refuses_unknown_order_before_control_writes(
    tmp_path: Path,
) -> None:
    repo = Path(__file__).resolve().parents[2]
    operations = repo / "output/hymenoptera_local_20261019"
    if not (operations / "continuous.py").is_file():
        pytest.skip("operator-local continuous CLI is absent; output/ is not versioned")
    paths = [operations / name for name in (".run.lock", "active_runtime.json", "active_scope.json")]
    before = {str(path): path.read_bytes() if path.exists() else None for path in paths}
    completed = subprocess.run(
        [
            sys.executable,
            str(operations / "continuous.py"),
            "--data-root",
            str(tmp_path / "data"),
            "--config-dir",
            str(tmp_path / "config"),
        ],
        env={**os.environ, "AMALGKIT_LOCAL_ACQUISITION_ORDER": "broken"},
        capture_output=True,
        text=True,
    )
    assert completed.returncode != 0
    assert "AMALGKIT_LOCAL_ACQUISITION_ORDER must be species_fair or size_bands" in completed.stderr
    assert {str(path): path.read_bytes() if path.exists() else None for path in paths} == before


@pytest.mark.external_tool
def test_size_ordered_native_partition_executes_in_order_and_replays_without_rewrite(
    tmp_path: Path, native_acquisition_factory
) -> None:
    fixture = native_acquisition_factory(accessions=("SRR125", "SRR123", "SRR124"))
    counts = {"SRR123": 20, "SRR124": 40, "SRR125": 90}
    raw = {}
    for accession, count in counts.items():
        path = fixture.work / "getfastq" / accession / f"{accession}.fastq.gz"
        with gzip.open(path, "rt") as handle:
            lines = handle.readlines()[: 4 * count]
        with gzip.open(path, "wt") as handle:
            handle.writelines(lines)
        raw[accession] = path
    with fixture.metadata.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        columns, rows = reader.fieldnames, list(reader)
    for row in rows:
        row.update(
            total_spots=str(counts[row["run"]]),
            total_bases=str(75 * counts[row["run"]]),
            size=str(raw[row["run"]].stat().st_size),
        )
    with fixture.metadata.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    tasks = [json.loads(line) for line in fixture.manifest.read_text().splitlines()]
    for task in tasks:
        task["fastq_bytes"] = raw[task["accession"]].stat().st_size
    fixture.manifest.write_text("".join(json.dumps(task) + "\n" for task in tasks))
    snapshot_path = fixture.manifest.with_name("snapshot.json")
    snapshot = json.loads(snapshot_path.read_text())
    snapshot["manifest_sha256"] = sha256_file(fixture.manifest)
    for record in snapshot["input_files"]:
        record["sha256"] = sha256_file(tmp_path / record["path"])
    snapshot_path.write_text(json.dumps(snapshot))
    ordered = plan_acquisition_order(tasks, [], global_size_bands=True)
    assert ordered == [
        "test_species/SRR123",
        "test_species/SRR124",
        "test_species/SRR125",
    ]
    allocation = write_allocation(allocate_tasks(fixture.manifest, backend="local"), tmp_path / "allocation")
    selection = allocation.parent / "small-first.json"
    selection.write_text(
        json.dumps(
            {
                "schema": "metainformant.rna.acquisition_partition.v1",
                "manifest_sha256": sha256_file(fixture.manifest),
                "snapshot_sha256": sha256_file(snapshot_path),
                "task_ids": ordered,
                "task_count": len(ordered),
                "backend": "local",
                "allocation_sha256": sha256_file(allocation),
            }
        )
    )
    source_bytes = fixture.manifest.read_bytes(), selection.read_bytes()
    kwargs = {
        "manifest_path": fixture.manifest,
        "data_root": fixture.data,
        "config_dir": fixture.config_dir,
        "task_selection": selection,
        "workers": 1,
        "threads": 2,
        "max_in_flight": 1,
        "fastq_threads": 1,
        "compression_threads": 1,
        "validation_slots": 1,
        "quant_slots": 1,
        "fasterq_slots": 1,
    }
    first = run_manifest(**kwargs)
    assert all(row["quantified"] and not row.get("skipped") for row in first["results"])
    journal = next((fixture.data / "acquisition_runs").glob("*/task_results.jsonl"))
    assert [json.loads(line)["task_id"] for line in journal.read_text().splitlines()] == ordered
    protected = [*raw.values(), *fixture.work.glob("quant/*/*abundance.tsv")]
    assert len(protected) == 6
    before = {str(path): (sha256_file(path), path.stat().st_mtime_ns) for path in protected}
    replay = run_manifest(**kwargs)
    assert all(row["quantified"] and row["skipped"] for row in replay["results"])
    assert {str(path): (sha256_file(path), path.stat().st_mtime_ns) for path in protected} == before
    assert (fixture.manifest.read_bytes(), selection.read_bytes()) == source_bytes
