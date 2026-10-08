"""Disjoint fleet admission and shared deadline reservations with real ledgers."""

from __future__ import annotations

import pytest

from metainformant.rna.engine.aws_completion import budget_allows, choose_partition
from metainformant.rna.engine.aws_fleet import in_flight_tasks, reserved_future_charge


def test_active_task_is_not_reissued_before_receipt_arrives() -> None:
    jobs = [{"status": "running", "task_ids": ["species/SRR1"]}]
    tasks = [{"task_id": f"species/SRR{i}", "accession": f"SRR{i}", "fastq_bytes": 10} for i in [1, 2]]
    partition = choose_partition(tasks, in_flight_tasks(jobs))
    assert [task["task_id"] for task in partition] == ["species/SRR2"]


@pytest.mark.parametrize("status", ["admitting", "running", "terminating"])
def test_every_live_phase_retains_task_ownership(status: str) -> None:
    assert in_flight_tasks([{"status": status, "task_ids": ["species/SRR1"]}]) == {"species/SRR1"}


def test_terminal_jobs_release_tasks_for_bounded_retry() -> None:
    assert in_flight_tasks([{"status": "terminated", "task_ids": ["species/SRR1"]}]) == set()


def test_duplicate_live_task_ownership_is_rejected() -> None:
    jobs = [{"status": "running", "task_ids": ["species/SRR1"]}] * 2
    with pytest.raises(ValueError, match="overlapping"):
        in_flight_tasks(jobs)


def test_existing_deadlines_cannot_be_spent_again_on_new_worker() -> None:
    jobs = [
        {"status": "running", "deadline": 14400, "hourly_upper_bound": 2.0},
        {"status": "running", "deadline": 14400, "hourly_upper_bound": 3.0},
        {"status": "terminated", "deadline": 14400, "hourly_upper_bound": 9.0},
    ]
    reserved = reserved_future_charge(jobs, now=3600, historical_hourly_bound=0.55)
    assert reserved == pytest.approx(15.0)
    assert budget_allows(720, 750, 14400, 2.0)
    assert not budget_allows(720 + reserved, 750, 14400, 2.0)


def test_legacy_job_preserves_rate_and_overdue_job_retains_shutdown_allowance() -> None:
    jobs = [
        {"status": "running", "deadline": 7200},
        {"status": "terminating", "deadline": 1000, "hourly_upper_bound": 3.0},
    ]
    assert reserved_future_charge(jobs, 3600, 0.55) == pytest.approx(0.65)


@pytest.mark.parametrize("now", [float("nan"), float("inf"), -1])
def test_invalid_clock_cannot_hide_reservations(now: float) -> None:
    with pytest.raises(ValueError):
        reserved_future_charge([], now, 0.55)


@pytest.mark.parametrize("deadline", [float("nan"), float("inf"), -1, True])
def test_invalid_deadline_cannot_drop_existing_reservation(deadline: float) -> None:
    with pytest.raises(ValueError):
        reserved_future_charge([{"status": "running", "deadline": deadline}], 3600, 0.55)


def test_once_wait_persists_live_workers_without_sleeping(tmp_path) -> None:
    import argparse
    import json

    from metainformant.rna.engine.aws_completion import _wait_for_fleet

    state = {
        "jobs": [{"status": "running"}],
        "locked_count": 3,
        "eligible_count": 4,
        "spent_upper_bound": 1.0,
        "reserved_future_upper_bound": 2.0,
    }
    path = tmp_path / "ledger.json"
    args = argparse.Namespace(once=True, poll_seconds=3600, max_workers=4)
    _wait_for_fleet(state, path, args, "waiting_for_budget_reservations")
    saved = json.loads(path.read_text())
    assert saved["status"] == "waiting_for_budget_reservations"
    assert saved["jobs"] == [{"status": "running"}]
