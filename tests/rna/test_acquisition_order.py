from __future__ import annotations

import pytest

from metainformant.rna.engine.acquisition_order import plan_acquisition_order


def test_checkpoint_work_is_interleaved_with_size_ordered_species_fair_acquisition() -> None:
    pending = [
        {"task_id": "ant/large", "species": "ant", "fastq_bytes": 900},
        {"task_id": "bee/unknown", "species": "bee"},
        {"task_id": "ant/small", "species": "ant", "fastq_bytes": 10},
        {"task_id": "bee/small", "species": "bee", "fastq_bytes": 20},
    ]
    checkpoints = ["sealed/c", "sealed/a", "sealed/b"]
    assert plan_acquisition_order(pending, checkpoints) == [
        "sealed/a",
        "ant/small",
        "sealed/b",
        "bee/small",
        "sealed/c",
        "ant/large",
        "bee/unknown",
    ]
    assert pending[0]["task_id"] == "ant/large"
    assert checkpoints == ["sealed/c", "sealed/a", "sealed/b"]


def test_empty_lanes_and_remaining_checkpoint_tail_are_preserved() -> None:
    assert plan_acquisition_order([], []) == []
    assert plan_acquisition_order([], ["b", "a"]) == ["a", "b"]
    assert plan_acquisition_order([{"task_id": "a", "species": "ant"}], []) == ["a"]


@pytest.mark.parametrize(
    "pending,checkpoints",
    [
        ([], ["a", "a"]),
        ([{"task_id": "a", "species": "ant"}], ["a"]),
        ([{"task_id": "a", "species": "ant"}, {"task_id": "a", "species": "ant"}], []),
    ],
)
def test_duplicate_tasks_fail_instead_of_silently_changing_scope(pending: list[dict], checkpoints: list[str]) -> None:
    with pytest.raises(ValueError, match="Duplicate acquisition task"):
        plan_acquisition_order(pending, checkpoints)


def test_global_size_order_is_monotone_across_species_with_interleaved_checkpoints() -> None:
    pending = [
        {"task_id": "ant/large", "species": "ant", "fastq_bytes": 80},
        {"task_id": "bee/tiny", "species": "bee", "fastq_bytes": 10},
        {"task_id": "ant/larger", "species": "ant", "fastq_bytes": 90},
        {"task_id": "bee/small", "species": "bee", "fastq_bytes": 20},
    ]
    original = [dict(task) for task in pending]
    checkpoints = ["sealed/b", "sealed/a"]
    assert plan_acquisition_order(pending, checkpoints, global_size_order=True) == [
        "sealed/a",
        "bee/tiny",
        "sealed/b",
        "bee/small",
        "ant/large",
        "ant/larger",
    ]
    assert pending == original and checkpoints == ["sealed/b", "sealed/a"]


def test_global_size_order_ties_unknowns_and_checkpoint_tail() -> None:
    pending = [
        {"task_id": "z/tie", "species": "z", "fastq_bytes": 10},
        {"task_id": "a/unknown", "species": "a"},
        {"task_id": "a/tie", "species": "a", "fastq_bytes": 10},
        {"task_id": "z/zero", "species": "z", "fastq_bytes": 0},
        {"task_id": "z/nonpositive", "species": "z", "fastq_bytes": -1},
    ]
    expected = ["a/tie", "z/tie", "z/nonpositive", "a/unknown", "z/zero"]
    assert plan_acquisition_order(pending, [], global_size_order=True) == expected
    assert plan_acquisition_order(reversed(pending), [], global_size_order=True) == expected
    assert plan_acquisition_order([], ["b", "a"], global_size_order=True) == ["a", "b"]
    assert plan_acquisition_order([], [], global_size_order=True) == []
    assert plan_acquisition_order(pending[:1], ["c", "b", "a"], global_size_order=True) == ["a", "z/tie", "b", "c"]


@pytest.mark.parametrize(
    "pending,checkpoints",
    [
        ([], ["a", "a"]),
        ([{"task_id": "a", "species": "ant"}], ["a"]),
        ([{"task_id": "a", "species": "ant"}, {"task_id": "a", "species": "ant"}], []),
    ],
)
def test_global_size_order_refuses_duplicate_identity(pending: list[dict], checkpoints: list[str]) -> None:
    with pytest.raises(ValueError, match="Duplicate acquisition task"):
        plan_acquisition_order(pending, checkpoints, global_size_order=True)


@pytest.mark.parametrize("size", ["broken", None])
def test_global_size_order_refuses_malformed_weights(size: object) -> None:
    with pytest.raises((TypeError, ValueError)):
        plan_acquisition_order([{"task_id": "a/x", "species": "a", "fastq_bytes": size}], [], global_size_order=True)


def test_global_size_policies_are_mutually_exclusive() -> None:
    with pytest.raises(ValueError, match="not both"):
        plan_acquisition_order([], [], global_size_order=True, global_size_bands=True)


def test_global_size_order_preserves_full_cohort_sized_scope_under_permutation() -> None:
    tasks = [
        {
            "task_id": f"species{i % 27}/run{i:05d}",
            "species": f"species{i % 27}",
            "fastq_bytes": (i * 7919) % 100003 if i % 89 else 0,
        }
        for i in range(18200)
    ]
    expected = [
        task["task_id"]
        for task in sorted(tasks, key=lambda task: (task["fastq_bytes"] <= 0, task["fastq_bytes"], task["task_id"]))
    ]
    ordered = plan_acquisition_order(iter(tasks), [], global_size_order=True)
    assert ordered == expected
    assert ordered == plan_acquisition_order(reversed(tasks), [], global_size_order=True)
    assert len(set(ordered)) == 18200
