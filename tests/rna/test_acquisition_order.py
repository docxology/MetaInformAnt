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
        "sealed/a", "ant/small", "sealed/b", "bee/small", "sealed/c", "ant/large", "bee/unknown"
    ]
    assert pending[0]["task_id"] == "ant/large"
    assert checkpoints == ["sealed/c", "sealed/a", "sealed/b"]


def test_empty_lanes_and_remaining_checkpoint_tail_are_preserved() -> None:
    assert plan_acquisition_order([], []) == []
    assert plan_acquisition_order([], ["b", "a"]) == ["a", "b"]
    assert plan_acquisition_order([{"task_id": "a", "species": "ant"}], []) == ["a"]


@pytest.mark.parametrize("pending,checkpoints", [
    ([], ["a", "a"]),
    ([{"task_id": "a", "species": "ant"}], ["a"]),
    ([{"task_id": "a", "species": "ant"}, {"task_id": "a", "species": "ant"}], []),
])
def test_duplicate_tasks_fail_instead_of_silently_changing_scope(pending: list[dict], checkpoints: list[str]) -> None:
    with pytest.raises(ValueError, match="Duplicate acquisition task"):
        plan_acquisition_order(pending, checkpoints)
