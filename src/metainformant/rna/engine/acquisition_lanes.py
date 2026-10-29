"""Bounded manifest execution with sample ownership across recovery transitions."""

from __future__ import annotations

import concurrent.futures
from collections import deque
from collections.abc import Callable, Sequence
from dataclasses import dataclass
from typing import Any

Task = dict[str, Any]
Result = dict[str, Any]


@dataclass(frozen=True)
class AdmissionDecision:
    """Admit, wait for existing reservations, or reject with a terminal reason."""

    wait: bool = False
    reason: str | None = None


def run_acquisition_lanes(
    tasks: Sequence[Task],
    *,
    workers: int,
    max_in_flight: int,
    recovery_slots: int,
    execute: Callable[[Task, bool], Result],
    recover: Callable[[Task], bool],
    uses_local_sra: Callable[[Task], bool],
    admission: Callable[[Task, Sequence[Task]], AdmissionDecision],
    finished: Callable[[Task, Result], None],
    should_stop: Callable[[], bool] | None = None,
    raw_bytes_limit: int = 0,
    dedicated_recovery: bool = True,
) -> dict[str, int]:
    """Drain admitted samples; STOP affects new admissions only.

    A primary future reserves its potential recovery position at admission:
    the entire admitted ledger, including queued transitions, is capped by
    ``max_in_flight``. A local archive is admitted only into a free recovery
    slot, allowing the scheduler to find ENA work further down the manifest.
    ``finished`` must persist/checkpoint before ownership is released. If it
    raises, executor shutdown drains submitted work without admitting or
    dispatching anything else and the original failure propagates.
    """

    if min(workers, max_in_flight, recovery_slots) < 1 or raw_bytes_limit < 0:
        raise ValueError("lane capacities must be positive and raw limit nonnegative")
    identities = [str(task["task_id"]) for task in tasks]
    if len(set(identities)) != len(identities):
        raise ValueError("duplicate manifest task identities")
    if raw_bytes_limit and any(
        not 0 < int(t.get("fastq_bytes", 0)) <= raw_bytes_limit for t in tasks
    ):
        raise ValueError("cloud task lacks a bounded raw-byte reservation")
    primary_pending: deque[Task] = deque()
    local_pending: deque[Task] = deque()
    for task in tasks:
        target = (
            local_pending
            if dedicated_recovery and uses_local_sra(task)
            else primary_pending
        )
        target.append(task)
    live: dict[str, Task] = {}
    ready: deque[Task] = deque()
    recovery: deque[Task] = deque()
    futures: dict[concurrent.futures.Future[Any], tuple[Task, str]] = {}
    reserved = 0
    peak_live = peak_raw = peak_recovery = 0
    deferred = recovered = 0

    def terminal(task: Task, result: Result) -> None:
        nonlocal reserved
        finished(task, result)
        del live[str(task["task_id"])]
        reserved -= int(task.get("fastq_bytes", 0))

    with (
        concurrent.futures.ThreadPoolExecutor(max_workers=workers) as primary_pool,
        concurrent.futures.ThreadPoolExecutor(
            max_workers=recovery_slots
        ) as recovery_pool,
    ):
        while primary_pending or local_pending or live:
            primary_count = sum(
                stage in ("primary", "ready") for _, stage in futures.values()
            )
            recovery_count = len(futures) - primary_count
            while ready and primary_count < workers:
                task = ready.popleft()
                futures[primary_pool.submit(execute, task, dedicated_recovery)] = (
                    task,
                    "ready",
                )
                primary_count += 1
            while recovery and recovery_count < recovery_slots:
                task = recovery.popleft()
                futures[recovery_pool.submit(recover, task)] = (task, "recover")
                recovery_count += 1

            for local_sra, pending in ((False, primary_pending), (True, local_pending)):
                while pending and len(live) < max_in_flight:
                    if should_stop is not None and should_stop():
                        break
                    if (local_sra and recovery_count >= recovery_slots) or (
                        not local_sra and primary_count >= workers
                    ):
                        break
                    task = pending[0]
                    needed = int(task.get("fastq_bytes", 0))
                    if raw_bytes_limit and reserved + needed > raw_bytes_limit:
                        break
                    decision = admission(task, list(live.values()))
                    if decision.wait:
                        break
                    pending.popleft()
                    if decision.reason:
                        finished(
                            task,
                            {
                                "quantified": False,
                                "status": "unresolved_deadline",
                                "error": decision.reason,
                            },
                        )
                        continue
                    live[str(task["task_id"])] = task
                    reserved += needed
                    if local_sra:
                        futures[recovery_pool.submit(execute, task, False)] = (
                            task,
                            "local",
                        )
                        recovery_count += 1
                    else:
                        futures[
                            primary_pool.submit(execute, task, dedicated_recovery)
                        ] = (task, "primary")
                        primary_count += 1
                    peak_live = max(peak_live, len(live))
                    peak_raw = max(peak_raw, reserved)
            peak_recovery = max(peak_recovery, len(recovery) + recovery_count)
            if not futures:
                if live:
                    raise RuntimeError("admitted samples have no executable stage")
                break
            done, _ = concurrent.futures.wait(
                tuple(futures), return_when=concurrent.futures.FIRST_COMPLETED
            )
            for future in done:
                task, stage = futures.pop(future)
                try:
                    outcome = future.result()
                except Exception as exc:  # noqa: BLE001 - external sample execution boundary
                    terminal(
                        task, {"quantified": False, "error": f"worker exception: {exc}"}
                    )
                    continue
                if stage == "recover":
                    if outcome:
                        ready.append(task)
                        recovered += 1
                    else:
                        terminal(
                            task,
                            {
                                "quantified": False,
                                "error": "dedicated NCBI recovery failed",
                            },
                        )
                elif outcome.get("fallback_deferred"):
                    if stage == "primary":
                        recovery.append(task)
                        deferred += 1
                    else:
                        terminal(
                            task,
                            {
                                "quantified": False,
                                "error": "recovery inputs unavailable after successful acquisition",
                            },
                        )
                else:
                    terminal(task, outcome)
            peak_recovery = max(
                peak_recovery,
                len(recovery)
                + sum(stage in ("recover", "local") for _, stage in futures.values()),
            )
    return {
        "unsubmitted": len(primary_pending) + len(local_pending),
        "peak_admitted": peak_live,
        "peak_reserved_raw_bytes": peak_raw,
        "peak_recovery_backlog": peak_recovery,
        "deferred": deferred,
        "recovered": recovered,
    }
