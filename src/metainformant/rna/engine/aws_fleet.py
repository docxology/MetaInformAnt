"""Disjoint task ownership and outstanding reservations for a bounded EC2 fleet."""

from __future__ import annotations

import math
from typing import NotRequired, Sequence, TypedDict

from metainformant.rna.engine.aws_resources import runtime_charge

LIVE_STATUSES = frozenset({"admitting", "running", "terminating"})
SHUTDOWN_RESERVE_SECONDS = 120


class TaskOwner(TypedDict):
    status: str
    task_ids: list[str]


class WorkerReservation(TypedDict):
    status: str
    deadline: float
    hourly_upper_bound: NotRequired[float]


def in_flight_tasks(jobs: Sequence[TaskOwner]) -> set[str]:
    """Keep tasks owned until a worker's termination is observed."""
    occupied: set[str] = set()
    for job in jobs:
        if job["status"] in LIVE_STATUSES:
            task_ids = set(job["task_ids"])
            if occupied.intersection(task_ids) or len(task_ids) != len(job["task_ids"]):
                raise ValueError("overlapping live task assignments")
            occupied.update(task_ids)
    return occupied


def reserved_future_charge(jobs: Sequence[WorkerReservation], now: float, legacy_hourly_bound: float) -> float:
    """Hold unspent worker deadlines plus an overdue shutdown allowance."""
    runtime_charge(0, now, legacy_hourly_bound)
    charges = []
    for job in jobs:
        if job["status"] in {"running", "terminating"}:
            deadline = job["deadline"]
            if isinstance(deadline, bool) or not math.isfinite(deadline) or deadline < 0:
                raise ValueError("invalid worker deadline")
            remaining = max(SHUTDOWN_RESERVE_SECONDS, deadline - now)
            charges.append(runtime_charge(0, remaining, job.get("hourly_upper_bound", legacy_hourly_bound)))
    return math.fsum(charges)
