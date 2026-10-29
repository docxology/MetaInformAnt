from __future__ import annotations

import concurrent.futures
from collections import deque
from collections.abc import Callable, Mapping, Sequence
from enum import StrEnum
from typing import TypedDict, assert_never

from metainformant.rna.engine.acquisition_lanes import AdmissionDecision, Result, Task
from metainformant.rna.engine.acquisition_size_overlay import effective_fastq_bytes


class Stage(StrEnum):
    ACQUIRING = "acquiring"
    RECOVERING = "recovering"
    REVALIDATING = "revalidating"
    QUANTIFYING = "quantifying"


class FrontierSummary(TypedDict):
    unsubmitted: int
    peak_admitted: int
    peak_reserved_raw_bytes: int
    peak_recovery_backlog: int
    peak_frontier: int
    deferred: int
    recovered: int
    reused_without_acquisition: int


def run_acquisition_frontier(
    tasks: Sequence[Task],
    *,
    acquisition_slots: int,
    quant_slots: int,
    recovery_slots: int,
    max_in_flight: int,
    frontier_size: int,
    raw_bytes_limit: int,
    acquire: Callable[[Task, bool], Result],
    quantify: Callable[[Task], Result],
    recover: Callable[[Task], bool],
    uses_local_sra: Callable[[Task], bool],
    admission: Callable[[Task, Sequence[Task]], AdmissionDecision],
    finished: Callable[[Task, Result], None],
    transition: Callable[[Task, str], None],
    should_stop: Callable[[], bool] | None = None,
    reuse: Callable[[Task], Result | None] | None = None,
    size_overrides: Mapping[str, int] | None = None,
) -> FrontierSummary:
    if min(acquisition_slots, quant_slots, recovery_slots, max_in_flight, frontier_size, raw_bytes_limit) <= 0:
        raise ValueError("frontier capacities and raw reservation limit must be positive")
    identities = [str(task["task_id"]) for task in tasks]
    if len(set(identities)) != len(identities):
        raise ValueError("duplicate manifest task identities")
    pending = deque(tasks)
    planned: deque[Task] = deque()
    planned_local: dict[str, bool] = {}
    ready: deque[Task] = deque()
    recovery: deque[Task] = deque()
    revalidate: deque[Task] = deque()
    live: dict[str, Task] = {}
    reserved_weights: dict[str, int] = {}
    futures: dict[concurrent.futures.Future[Result], tuple[Task, Stage]] = {}
    reuse_probed: set[str] = set()
    reused_without_acquisition = 0
    reserved = peak_live = peak_raw = peak_recovery = peak_frontier = deferred = recovered = 0

    def stopped() -> bool:
        return should_stop is not None and should_stop()

    def terminal(task: Task, result: Result) -> None:
        nonlocal reserved
        finished(task, result)
        transition(task, "checkpointed" if result.get("quantified") else "deferred_or_failed")
        del live[str(task["task_id"])]
        reserved -= reserved_weights.pop(str(task["task_id"]))

    def recover_result(task: Task) -> Result:
        return {"recovered": recover(task)}

    with (
        concurrent.futures.ThreadPoolExecutor(max_workers=acquisition_slots) as downloads,
        concurrent.futures.ThreadPoolExecutor(max_workers=quant_slots) as quantifiers,
        concurrent.futures.ThreadPoolExecutor(max_workers=recovery_slots) as extractors,
    ):
        while pending or planned or live:
            acquiring = sum(stage in (Stage.ACQUIRING, Stage.REVALIDATING) for _, stage in futures.values())
            quantifying = sum(stage is Stage.QUANTIFYING for _, stage in futures.values())
            recovering = sum(stage is Stage.RECOVERING for _, stage in futures.values())
            while ready and quantifying < quant_slots:
                task = ready.popleft()
                transition(task, Stage.QUANTIFYING)
                futures[quantifiers.submit(quantify, task)] = task, Stage.QUANTIFYING
                quantifying += 1
            while recovery and recovering < recovery_slots:
                task = recovery.popleft()
                transition(task, Stage.RECOVERING)
                futures[extractors.submit(recover_result, task)] = task, Stage.RECOVERING
                recovering += 1
            while revalidate and acquiring < acquisition_slots:
                task = revalidate.popleft()
                transition(task, Stage.REVALIDATING)
                futures[downloads.submit(acquire, task, False)] = task, Stage.REVALIDATING
                acquiring += 1
            while pending and len(planned) + len(live) < frontier_size and not stopped():
                task = pending.popleft()
                transition(task, "planned")
                planned.append(task)
                planned_local[str(task["task_id"])] = uses_local_sra(task)
                peak_frontier = max(peak_frontier, len(planned))
            while planned and len(live) < max_in_flight and not stopped():
                available = next(
                    (
                        index
                        for index, entry in enumerate(planned)
                        if (planned_local[str(entry["task_id"])] and recovering < recovery_slots)
                        or (not planned_local[str(entry["task_id"])] and acquiring < acquisition_slots)
                    ),
                    None,
                )
                if available is None:
                    break
                task = planned[available]
                weight = effective_fastq_bytes(task, size_overrides)
                if weight <= 0 or weight > raw_bytes_limit:
                    del planned[available]
                    del planned_local[str(task["task_id"])]
                    finished(
                        task,
                        {
                            "quantified": False,
                            "status": "unresolved_deadline",
                            "error": "raw reservation missing or exceeds the bounded frontier",
                        },
                    )
                    transition(task, "reservation_refused")
                    continue
                local = planned_local[str(task["task_id"])]
                if reserved + weight > raw_bytes_limit:
                    break
                decision = admission(task, list(live.values()))
                if decision.wait:
                    break
                del planned[available]
                del planned_local[str(task["task_id"])]
                if decision.reason:
                    finished(task, {"quantified": False, "status": "unresolved_deadline", "error": decision.reason})
                    transition(task, "admission_refused")
                    continue
                live[str(task["task_id"])] = task
                reserved_weights[str(task["task_id"])] = weight
                reserved += weight
                peak_live = max(peak_live, len(live))
                peak_raw = max(peak_raw, reserved)
                if local:
                    transition(task, Stage.RECOVERING)
                    futures[extractors.submit(recover_result, task)] = task, Stage.RECOVERING
                    recovering += 1
                else:
                    transition(task, Stage.ACQUIRING)
                    futures[downloads.submit(acquire, task, True)] = task, Stage.ACQUIRING
                    acquiring += 1
            probe_progress = False
            if reuse is not None:
                for task in tuple(planned):
                    identity = str(task["task_id"])
                    weight = effective_fastq_bytes(task, size_overrides)
                    if stopped() or len(live) >= max_in_flight:
                        break
                    if (
                        identity in reuse_probed
                        or not 0 < weight <= raw_bytes_limit
                        or reserved + weight > raw_bytes_limit
                    ):
                        continue
                    decision = admission(task, list(live.values()))
                    if decision.wait or decision.reason:
                        continue
                    reuse_probed.add(identity)
                    probe_progress = True
                    result = reuse(task)
                    if result is None:
                        break
                    if result.get("quantified") is not True or result.get("skipped") is not True:
                        raise ValueError(f"reuse callback must return a verified reused success: {identity}")
                    if stopped():
                        break
                    planned.remove(task)
                    del planned_local[identity]
                    live[identity] = task
                    reserved_weights[identity] = weight
                    reserved += weight
                    peak_live = max(peak_live, len(live))
                    peak_raw = max(peak_raw, reserved)
                    transition(task, "reusing")
                    terminal(task, result)
                    reused_without_acquisition += 1
                    break
            peak_recovery = max(peak_recovery, len(recovery) + recovering)
            if not futures:
                if probe_progress:
                    continue
                if live:
                    raise RuntimeError("admitted frontier has no executable stage")
                if not planned and pending and not stopped():
                    continue
                if planned and not stopped():
                    raise RuntimeError("unadmitted frontier blocked without active reservations")
                break
            done, _ = concurrent.futures.wait(
                tuple(futures),
                timeout=0 if probe_progress else None,
                return_when=concurrent.futures.FIRST_COMPLETED,
            )
            for future in done:
                task, stage = futures.pop(future)
                try:
                    result = future.result()
                except (OSError, RuntimeError, ValueError) as exc:
                    terminal(task, {"quantified": False, "error": f"stage execution failed: {exc}"})
                    continue
                match stage:
                    case Stage.RECOVERING:
                        if result["recovered"]:
                            recovered += 1
                            revalidate.append(task)
                            transition(task, "recovered_ready_for_validation")
                        else:
                            terminal(task, {"quantified": False, "error": "dedicated NCBI recovery failed"})
                    case Stage.ACQUIRING | Stage.REVALIDATING:
                        if result.get("ready_for_quant"):
                            ready.append(task)
                            transition(task, "ready_for_quant")
                        elif result.get("fallback_deferred") and stage is Stage.ACQUIRING:
                            deferred += 1
                            recovery.append(task)
                            transition(task, "recovery_queued")
                        elif result.get("fallback_deferred"):
                            terminal(task, {"quantified": False, "error": "recovered inputs failed revalidation"})
                        else:
                            terminal(task, result)
                    case Stage.QUANTIFYING:
                        terminal(task, result)
                    case unreachable:
                        assert_never(unreachable)
    return {
        "unsubmitted": len(pending) + len(planned),
        "peak_admitted": peak_live,
        "peak_reserved_raw_bytes": peak_raw,
        "peak_recovery_backlog": peak_recovery,
        "peak_frontier": peak_frontier,
        "deferred": deferred,
        "recovered": recovered,
        "reused_without_acquisition": reused_without_acquisition,
    }
