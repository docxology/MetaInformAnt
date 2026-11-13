from __future__ import annotations

import csv
import gzip
import json
from pathlib import Path
from typing import TypedDict

import pytest

from metainformant.rna.engine.acquisition_manifest import sha256_file
from metainformant.rna.engine.acquisition_order import plan_acquisition_order

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


@pytest.mark.external_tool
@pytest.mark.parametrize("strict", [False, True], ids=["size_bands", "smallest_first"])
def test_native_samples_quantify_in_planned_order_and_replay_without_rewrite(
    tmp_path: Path, native_acquisition_factory, strict: bool
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
    ordered = plan_acquisition_order(tasks, [], global_size_bands=not strict, global_size_order=strict)
    assert ordered == [
        "test_species/SRR123",
        "test_species/SRR124",
        "test_species/SRR125",
    ]
    accessions = [task_id.split("/", 1)[1] for task_id in ordered]
    source_bytes = fixture.manifest.read_bytes()
    first = [fixture.produce(accession) for accession in accessions]
    assert [row["srr"] for row in first] == accessions
    assert all(row["quantified"] and not row["skipped"] and row["error"] is None for row in first), first
    protected = [*raw.values(), *fixture.work.glob("quant/*/*abundance.tsv")]
    assert len(protected) == 6
    before = {str(path): (sha256_file(path), path.stat().st_mtime_ns) for path in protected}
    replay = [fixture.produce(accession) for accession in accessions]
    assert all(row["quantified"] and row["skipped"] for row in replay), replay
    assert {str(path): (sha256_file(path), path.stat().st_mtime_ns) for path in protected} == before
    assert fixture.manifest.read_bytes() == source_bytes
