"""Deterministic admission arithmetic with real sizes and deadline counterexamples."""

from __future__ import annotations

import math
from dataclasses import replace

import pytest

from metainformant.rna.engine.acquisition_estimates import AcquisitionEstimateError
from metainformant.rna.engine.acquisition_scheduling import (
    PlanningAssumptions,
    Workload,
    fits_worker_deadline,
    positive_size,
    workload_seconds,
)
from metainformant.rna.engine.aws_completion import choose_deadline_partition


def assumptions() -> PlanningAssumptions:
    return PlanningAssumptions(
        100,
        1000,
        100,
        "deterministic arithmetic scenario",
        setup_seconds=10,
        drain_seconds=10,
        task_overhead_seconds=10,
        safety_factor=1,
    )


def test_batch_accounts_for_total_transfer_and_stage_capacity() -> None:
    tasks = [Workload(str(i), 1000, 1000, requires_extraction=True) for i in range(4)]
    # 40 transfer +4 extraction +17.5 order-independent quant bound +40 overhead +20 reserve.
    assert workload_seconds(tasks, assumptions()) == 122
    assert workload_seconds(tasks, replace(assumptions(), quant_slots=1)) == 144
    assert workload_seconds(tasks, replace(assumptions(), extraction_slots=2)) == 120


@pytest.mark.parametrize("value", [None, "", "nan", "inf", math.nan, math.inf, -1, 0, True, 1.2])
def test_unknown_or_nonfinite_sizes_cannot_enter_admission(
    value: float | str | None,
) -> None:
    with pytest.raises(AcquisitionEstimateError):
        positive_size(value, "size")


@pytest.mark.parametrize("rate", [0, -1, math.nan, math.inf, True])
def test_invalid_rates_cannot_create_planning_envelope(rate: float) -> None:
    with pytest.raises(AcquisitionEstimateError):
        replace(assumptions(), quant_bases_per_second=rate)


def task(i: int, raw: int | float = 1000, bases: int | None = 1000) -> dict:
    return dict(task_id=f"s/SRR{i}", accession=f"SRR{i}", fastq_bytes=raw, total_bases=bases)


def test_batch_time_cap_splits_work_without_mutating_input() -> None:
    tasks = [task(i) for i in range(4)]
    selected, unresolved = choose_deadline_partition(
        tasks,
        set(),
        assumptions=assumptions(),
        max_bytes=10000,
        max_tasks=10,
        target_seconds=100,
        maximum_seconds=200,
    )
    assert len(selected) == 3
    assert unresolved == {}
    assert all("planning_seconds" not in row for row in tasks)


def test_unsupported_tasks_do_not_discard_supported_work() -> None:
    selected, unresolved = choose_deadline_partition(
        [task(0, bases=None), task(1), task(2, raw=math.inf), task(3, raw=100000)],
        set(),
        assumptions=assumptions(),
        max_bytes=10000,
        max_tasks=10,
        target_seconds=100,
        maximum_seconds=200,
    )
    assert [t["task_id"] for t in selected] == ["s/SRR1"]
    assert set(unresolved) == {"s/SRR0", "s/SRR2", "s/SRR3"}


def test_reviewable_oversized_singleton_has_full_duration() -> None:
    selected, unresolved = choose_deadline_partition(
        [task(1, raw=10000)],
        set(),
        assumptions=assumptions(),
        max_bytes=1000,
        max_tasks=10,
        target_seconds=100,
        maximum_seconds=200,
    )
    assert len(selected) == 1 and unresolved == {}
    assert selected[0]["planning_seconds"] == 120


def test_deadline_includes_unfinished_tasks_and_drain() -> None:
    assert fits_worker_deadline(now=100, deadline=160, task_seconds=30, reserved_seconds=20, drain_seconds=10)
    assert not fits_worker_deadline(now=101, deadline=160, task_seconds=30, reserved_seconds=20, drain_seconds=10)
    assert not fits_worker_deadline(now=150, deadline=160, task_seconds=30, reserved_seconds=0, drain_seconds=10)


@pytest.mark.parametrize("value", [0, -1, math.nan, math.inf, True])
def test_invalid_task_estimate_never_starts(value: float) -> None:
    with pytest.raises(AcquisitionEstimateError):
        fits_worker_deadline(
            now=100,
            deadline=200,
            task_seconds=value,
            reserved_seconds=0,
            drain_seconds=10,
        )


def test_source_transfer_is_separate_from_raw_scratch_reservation() -> None:
    work = Workload("s/SRR1", raw_bytes=1000000, total_bases=1000, transfer_bytes=1000)
    assert workload_seconds([work], assumptions()) == 50
    assert work.raw_bytes == 1000000


def test_existing_inflight_ownership_survives_new_selection() -> None:
    from metainformant.rna.engine.aws_fleet import in_flight_tasks

    jobs = [
        {"status": state, "task_ids": [f"s/SRR{i}"]} for i, state in enumerate(("admitting", "running", "terminating"))
    ]
    selected, _ = choose_deadline_partition(
        [task(i) for i in range(4)],
        in_flight_tasks(jobs),
        assumptions=assumptions(),
        max_bytes=10000,
        max_tasks=10,
        target_seconds=100,
        maximum_seconds=200,
    )
    assert [row["task_id"] for row in selected] == ["s/SRR3"]
    assert [job["status"] for job in jobs] == ["admitting", "running", "terminating"]


def test_missing_extraction_calibration_blocks_only_sra() -> None:
    sra = {**task(2), "source_evidence_sha256": "a" * 64, "sra_bytes": 100}
    selected, unresolved = choose_deadline_partition(
        [task(1), sra],
        set(),
        assumptions=replace(assumptions(), extraction_bases_per_second=0),
        max_bytes=10000,
        max_tasks=10,
        target_seconds=100,
        maximum_seconds=200,
    )
    assert [t["task_id"] for t in selected] == ["s/SRR1"]
    assert set(unresolved) == {"s/SRR2"}


def test_sra_cannot_claim_ena_duration_without_extraction_evidence() -> None:
    with pytest.raises(AcquisitionEstimateError):
        workload_seconds(
            [Workload("s/SRR1", 1000, 1000, requires_extraction=True)],
            replace(assumptions(), extraction_bases_per_second=0),
        )


def test_legacy_partition_unknown_sizes_do_not_hide_valid_tasks() -> None:
    from metainformant.rna.engine.aws_completion import choose_partition

    rows = [task(1, raw=math.nan), task(2), task(3, raw=math.inf)]
    assert [row["task_id"] for row in choose_partition(rows, set())] == ["s/SRR2"]


def test_integer_workload_counts_are_not_rounded_through_float() -> None:
    assert positive_size("9007199254740993", "bytes") == 9007199254740993


def test_parallel_worker_reserves_shared_stage_capacity_not_serial_quant() -> None:
    from metainformant.rna.engine.acquisition_scheduling import task_workload

    tasks = [task_workload(task(i)) for i in range(4)]
    aggregate = workload_seconds(tasks, assumptions(), include_setup=False)
    serialized = sum(workload_seconds([t], assumptions(), include_setup=False) for t in tasks)
    assert aggregate == 98 and serialized == 120
    assert fits_worker_deadline(now=100, deadline=210, task_seconds=aggregate, reserved_seconds=0, drain_seconds=10)
    assert not fits_worker_deadline(
        now=100, deadline=210, task_seconds=serialized, reserved_seconds=0, drain_seconds=10
    )


def test_archive_binds_planning_profile_without_changing_frozen_inputs(tmp_path) -> None:
    import hashlib
    import json
    import tarfile
    from dataclasses import asdict

    from metainformant.rna.engine.aws_inputs import _inputs_bundle

    work = tmp_path / "inputs/s/work"
    (work / "metadata").mkdir(parents=True)
    (work / "index").mkdir()
    metadata = work / "metadata/metadata_selected.tsv"
    metadata.write_text("run\ttotal_bases\nSRR1\t1000\n")
    index = work / "index/s.idx"
    index.write_bytes(b"opaque reference binding for archive test")
    config = tmp_path / "s.yaml"
    config.write_text("species_list: [s]\n")
    hashes = {p: hashlib.sha256(p.read_bytes()).hexdigest() for p in (metadata, index, config)}
    species = dict(
        species="s",
        index_name="s.idx",
        metadata_sha256=hashes[metadata],
        index_sha256=hashes[index],
        config_sha256=hashes[config],
    )
    bundle, digest = _inputs_bundle(
        tmp_path, species, [task(1)], tmp_path / "job", config_path=config, planning=assumptions()
    )
    with tarfile.open(bundle) as archive:
        snapshot = json.load(archive.extractfile("snapshot.json"))
    assert snapshot["planning_assumptions"] == asdict(assumptions())
    assert snapshot["task_count"] == 1
    assert digest == hashlib.sha256(bundle.read_bytes()).hexdigest()
    assert all(hashlib.sha256(p.read_bytes()).hexdigest() == expected for p, expected in hashes.items())


def test_workload_parses_read_only_scalar_manifest_mapping() -> None:
    from types import MappingProxyType

    from metainformant.rna.engine.acquisition_scheduling import task_workload

    manifest = MappingProxyType({"task_id": "s/SRR1", "fastq_bytes": "1000", "total_bases": 1000.0})
    assert task_workload(manifest) == Workload("s/SRR1", 1000, 1000, transfer_bytes=1000)


def test_workload_rejects_nonstring_manifest_identity() -> None:
    from metainformant.rna.engine.acquisition_scheduling import task_workload

    with pytest.raises(AcquisitionEstimateError, match="task_id"):
        task_workload({"task_id": 123, "fastq_bytes": 1000, "total_bases": 1000})
