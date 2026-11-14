"""Acquisition weights: frozen ``fastq_bytes`` and the supplements that may fill a zero weight."""

from __future__ import annotations

import copy

import pytest

from metainformant.rna.engine.acquisition_order import plan_acquisition_order
from metainformant.rna.engine.acquisition_size_overlay import effective_fastq_bytes

TASKS = [
    {"task_id": "species_b/SRR1", "species": "species_b", "fastq_bytes": 0},
    {"task_id": "species_a/SRR2", "species": "species_a", "fastq_bytes": 90},
]


def test_frozen_weight_is_used_without_or_outside_supplements() -> None:
    assert effective_fastq_bytes(TASKS[0]) == 0
    assert effective_fastq_bytes(TASKS[1]) == 90
    assert effective_fastq_bytes(TASKS[1], {"species_b/SRR1": 30}) == 90
    assert effective_fastq_bytes({"task_id": "x/SRR3", "species": "x"}) == 0


def test_supplement_fills_a_zero_weight_without_rewriting_the_task() -> None:
    task = copy.deepcopy(TASKS[0])
    assert effective_fastq_bytes(task, {"species_b/SRR1": 30}) == 30
    assert task == TASKS[0]


@pytest.mark.parametrize(
    ("task", "size"),
    [
        (TASKS[1], 30),  # a supplement never overwrites a positive frozen weight
        (TASKS[0], 0),
        (TASKS[0], -1),
        (TASKS[0], True),
        (TASKS[0], 30.0),
        (TASKS[0], "30"),
        (TASKS[0], 2**63),
    ],
)
def test_bad_supplement_refuses(task: dict, size: object) -> None:
    with pytest.raises(ValueError, match="only fill a zero weight"):
        effective_fastq_bytes(task, {task["task_id"]: size})


def test_order_uses_supplements_without_rewriting_tasks() -> None:
    tasks = copy.deepcopy(TASKS)
    assert plan_acquisition_order(tasks, ["checkpoint/old"], global_size_order=True) == [
        "checkpoint/old",
        "species_a/SRR2",
        "species_b/SRR1",
    ]
    assert plan_acquisition_order(
        tasks, ["checkpoint/old"], global_size_order=True, size_overrides={"species_b/SRR1": 30}
    ) == ["checkpoint/old", "species_b/SRR1", "species_a/SRR2"]
    assert plan_acquisition_order(tasks, [], global_size_bands=True, size_overrides={"species_b/SRR1": 30}) == [
        "species_a/SRR2",
        "species_b/SRR1",
    ]
    assert tasks == TASKS
    with pytest.raises(ValueError, match="only fill a zero weight"):
        plan_acquisition_order(tasks, [], global_size_order=True, size_overrides={"species_a/SRR2": 30})
