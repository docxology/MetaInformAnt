from __future__ import annotations

import json
import threading
from pathlib import Path

import pytest

from metainformant.rna.engine.acquisition_lanes import (
    AdmissionDecision,
    run_acquisition_lanes,
)


def tasks(count: int) -> list[dict]:
    return [{"task_id": str(i), "fastq_bytes": 10} for i in range(count)]


@pytest.mark.parametrize("dedicated", [False, True])
def test_local_recovery_progress_has_serial_negative_control(
    tmp_path: Path, dedicated: bool
) -> None:
    entered, release, checkpoint = (
        threading.Event(),
        threading.Event(),
        threading.Event(),
    )
    errors: list[BaseException] = []
    results: list[str] = []

    def execute(task: dict, defer: bool) -> dict:
        if task["task_id"] == "0":
            entered.set()
            if not release.wait(10):
                raise TimeoutError("test recovery was never released")
        return {"quantified": True}

    def finish(task: dict, result: dict) -> None:
        results.append(task["task_id"])
        (tmp_path / task["task_id"]).write_text(json.dumps(result))
        if task["task_id"] == "1":
            checkpoint.set()

    def run() -> None:
        try:
            run_acquisition_lanes(
                tasks(2),
                workers=1,
                max_in_flight=2,
                recovery_slots=1,
                execute=execute,
                recover=lambda task: True,
                uses_local_sra=lambda task: task["task_id"] == "0",
                admission=lambda task, live: AdmissionDecision(),
                finished=finish,
                dedicated_recovery=dedicated,
            )
        except (AssertionError, OSError, RuntimeError, ValueError) as exc:
            errors.append(exc)

    worker = threading.Thread(target=run)
    worker.start()
    try:
        assert entered.wait(2)
        assert checkpoint.wait(0.5) is dedicated
    finally:
        release.set()
        worker.join(10)
    assert not worker.is_alive()
    assert errors == []
    assert sorted(results) == ["0", "1"]
    assert all(
        json.loads((tmp_path / str(i)).read_text())["quantified"] for i in range(2)
    )


def test_deferrals_preserve_raw_and_deadline_reservations_through_checkpoint() -> None:
    admitted: list[tuple[str, tuple[str, ...]]] = []
    attempts: dict[str, int] = {}
    terminals: list[str] = []

    def execute(task: dict, defer: bool) -> dict:
        identity = task["task_id"]
        attempts[identity] = attempts.get(identity, 0) + 1
        return (
            {"fallback_deferred": True}
            if attempts[identity] == 1
            else {"quantified": True}
        )

    def admission(task: dict, live: list[dict]) -> AdmissionDecision:
        admitted.append((task["task_id"], tuple(t["task_id"] for t in live)))
        assert sum(t["fastq_bytes"] for t in live) + task["fastq_bytes"] <= 20
        return AdmissionDecision()

    summary = run_acquisition_lanes(
        tasks(12),
        workers=2,
        max_in_flight=3,
        recovery_slots=1,
        execute=execute,
        recover=lambda task: True,
        uses_local_sra=lambda task: False,
        admission=admission,
        finished=lambda task, result: terminals.append(task["task_id"]),
        raw_bytes_limit=20,
    )
    assert len(terminals) == len(set(terminals)) == 12
    assert len(admitted) == 12
    assert set(attempts.values()) == {2}
    assert summary["peak_admitted"] == 2
    assert summary["peak_reserved_raw_bytes"] == 20
    assert summary["peak_recovery_backlog"] <= 2
    assert summary["deferred"] == summary["recovered"] == 12


def test_stop_drains_admitted_recovery_without_new_admission() -> None:
    stop = threading.Event()
    attempts: dict[str, int] = {}
    terminal: list[str] = []

    def execute(task: dict, defer: bool) -> dict:
        identity = task["task_id"]
        attempts[identity] = attempts.get(identity, 0) + 1
        if attempts[identity] == 1:
            stop.set()
            return {"fallback_deferred": True}
        return {"quantified": True}

    summary = run_acquisition_lanes(
        tasks(5),
        workers=1,
        max_in_flight=2,
        recovery_slots=1,
        execute=execute,
        recover=lambda task: True,
        uses_local_sra=lambda task: False,
        admission=lambda task, live: AdmissionDecision(),
        finished=lambda task, result: terminal.append(task["task_id"]),
        should_stop=stop.is_set,
    )
    assert terminal == ["0"]
    assert attempts == {"0": 2}
    assert summary["unsubmitted"] == 4


@pytest.mark.parametrize("recover_result", [False, True])
def test_failed_or_repeated_recovery_is_terminal(recover_result: bool) -> None:
    results: list[dict] = []
    summary = run_acquisition_lanes(
        tasks(1),
        workers=1,
        max_in_flight=1,
        recovery_slots=1,
        execute=lambda task, defer: {"fallback_deferred": True},
        recover=lambda task: recover_result,
        uses_local_sra=lambda task: False,
        admission=lambda task, live: AdmissionDecision(),
        finished=lambda task, result: results.append(result),
    )
    assert summary["deferred"] == 1
    assert len(results) == 1
    assert results[0]["quantified"] is False
    assert "recovery" in results[0]["error"]


def test_checkpoint_failure_prevents_replacement_admission(tmp_path: Path) -> None:
    executed: list[str] = []

    def execute(task: dict, defer: bool) -> dict:
        executed.append(task["task_id"])
        (tmp_path / task["task_id"]).write_bytes(b"retained raw input")
        return {"quantified": True}

    def checkpoint(task: dict, result: dict) -> None:
        raise OSError("retention unavailable")

    with pytest.raises(OSError, match="retention unavailable"):
        run_acquisition_lanes(
            tasks(3),
            workers=1,
            max_in_flight=1,
            recovery_slots=1,
            execute=execute,
            recover=lambda task: True,
            uses_local_sra=lambda task: False,
            admission=lambda task, live: AdmissionDecision(),
            finished=checkpoint,
        )
    assert executed == ["0"]
    assert (tmp_path / "0").read_bytes() == b"retained raw input"


def test_rejected_deadline_and_worker_exception_are_recorded_once() -> None:
    results: list[tuple[str, dict]] = []

    def execute(task: dict, defer: bool) -> dict:
        raise ValueError("real execution failure")

    run_acquisition_lanes(
        tasks(2),
        workers=1,
        max_in_flight=2,
        recovery_slots=1,
        execute=execute,
        recover=lambda task: True,
        uses_local_sra=lambda task: False,
        admission=lambda task, live: (
            AdmissionDecision(reason="deadline")
            if task["task_id"] == "0"
            else AdmissionDecision()
        ),
        finished=lambda task, result: results.append((task["task_id"], result)),
    )
    assert [identity for identity, _ in results] == ["0", "1"]
    assert results[0][1]["status"] == "unresolved_deadline"
    assert "real execution failure" in results[1][1]["error"]


@pytest.mark.parametrize("count,limit", [(0, -1), (1, 9)])
def test_invalid_reservation_contract_rejects_before_execution(
    count: int, limit: int
) -> None:
    with pytest.raises(ValueError):
        run_acquisition_lanes(
            tasks(count),
            workers=1,
            max_in_flight=1,
            recovery_slots=1,
            execute=lambda task, defer: {},
            recover=lambda task: True,
            uses_local_sra=lambda task: False,
            admission=lambda task, live: AdmissionDecision(),
            finished=lambda task, result: None,
            raw_bytes_limit=limit,
        )


def test_routing_classifies_each_pending_sample_once() -> None:
    classified: list[str] = []

    def route(task: dict) -> bool:
        classified.append(task["task_id"])
        return False

    summary = run_acquisition_lanes(
        tasks(300),
        workers=1,
        max_in_flight=2,
        recovery_slots=1,
        execute=lambda task, defer: {"quantified": True},
        recover=lambda task: True,
        uses_local_sra=route,
        admission=lambda task, live: AdmissionDecision(),
        finished=lambda task, result: None,
    )
    assert classified == [str(i) for i in range(300)]
    assert summary["unsubmitted"] == 0
