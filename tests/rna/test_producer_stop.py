"""A stopping producer never frees the writer lock while its child process groups can still write.

Each test runs a real producer process (``StreamingPipelineOrchestrator.run_all``
with the real lock, scheduler, stop handling and process-group registry). Its
per-sample worker starts a real child through the production
``_run_command_in_process_group`` (a new session, like curl, fasterq-dump and
kallisto): the child heartbeats into a file until it is killed, and in
``stubborn`` mode it ignores SIGTERM like a process stuck in uninterruptible
I/O. The tests signal the producer exactly as a launcher, an operator's
Ctrl-C or a force quit would, and probe the lock from this process.
"""

from __future__ import annotations

import json
import os
import signal
import subprocess
import sys
import time
from pathlib import Path

import pytest

from metainformant.rna.engine.campaign_lock import (
    EXIT_RETRYABLE,
    TOKEN_ENV,
    CampaignLockBusy,
    acquire,
    is_quiescent,
)
from metainformant.rna.engine.producer_stop import ProducerStop
from metainformant.rna.engine.progress_db import ProgressDB
from metainformant.rna.retrieval import process_groups

REPO_ROOT = Path(__file__).resolve().parents[2]
RUNS = ["SRR1", "SRR2", "SRR3", "SRR4"]

CHILD = """
import os, signal, sys, time
beat, mode = sys.argv[1], sys.argv[2]
if mode == "stubborn":
    signal.signal(signal.SIGTERM, signal.SIG_IGN)
with open(beat + ".tmp", "w") as handle:
    handle.write(str(os.getpid()))
os.replace(beat + ".tmp", beat + ".pid")
while True:
    with open(beat, "a") as handle:
        handle.write("x")
    time.sleep(0.05)
"""

PRODUCER = """
import json, os, sys
from pathlib import Path
root, mode = Path(sys.argv[1]), sys.argv[2]
from metainformant.rna.engine import streaming_orchestrator as so
from metainformant.rna.retrieval.ena_downloader import _run_command_in_process_group

def preflight(data_root, **_kwargs):  # the preflight has its own tests; this host is taken as ready
    return {"data_root": str(data_root)}

so.run_campaign_preflight = preflight
orchestrator = so.StreamingPipelineOrchestrator(root / "config", root / "logs", db_path=root / "data" / "progress.db")
runs = %(runs)r

def discover(_config_name, *_args):
    orchestrator.db.init_species("apis", runs)
    return [
        {"srr": run, "batch_idx": i + 1, "fastq_dir": root / "fastq", "config_path": root / "apis.yaml",
         "species_name": "apis", "expected_paired": False}
        for i, run in enumerate(runs)
    ]

def worker(srr, *_args, **_kwargs):
    orchestrator.db.set_state("apis", srr, "downloading")
    try:
        _run_command_in_process_group([sys.executable, "-c", %(child)r, str(root / "beats" / srr), mode], 600)
    except Exception as error:  # a refused spawn after the stop began, or a timeout
        print("worker", srr, type(error).__name__, flush=True)
    orchestrator.db.set_state("apis", srr, "failed", error="Download Failed (all sources: ENA FTP/HTTP, NCBI)")
    return {"quantified": False, "skipped": False, "error": "download interrupted"}

orchestrator.discover_species_tasks = discover
orchestrator.process_single_sample = worker
summary = orchestrator.run_all(["amalgkit_apis.yaml"], max_gb=1000.0, workers=2, threads=1, discovery_workers=1)
print("SUMMARY " + json.dumps({key: summary[key] for key in ("status", "exit_code", "stopped_by")}), flush=True)
sys.exit(summary["exit_code"])
""" % {
    "runs": RUNS,
    "child": CHILD,
}


class Producer:
    """A producer process over a fresh data root, plus probes of its children and the lock."""

    def __init__(self, base: Path, mode: str, *, grace_seconds: float) -> None:
        self.base = base
        self.data_root = base / "data"
        self.data_root.mkdir(parents=True)
        (base / "beats").mkdir()
        environment = os.environ.copy()
        environment.pop(TOKEN_ENV, None)
        environment["AMALGKIT_DATA_ROOT"] = str(self.data_root)
        environment["AMALGKIT_STOP_GRACE_SECONDS"] = str(grace_seconds)
        environment["PYTHONPATH"] = os.pathsep.join(filter(None, [str(REPO_ROOT), environment.get("PYTHONPATH")]))
        self.process = subprocess.Popen(
            [sys.executable, "-c", PRODUCER, str(base), mode],
            env=environment,
            stdout=subprocess.PIPE,
            stderr=subprocess.DEVNULL,
            text=True,
            cwd=REPO_ROOT,
        )

    def child_pids(self, count: int, timeout: float = 60.0) -> list[int]:
        deadline = time.monotonic() + timeout
        while time.monotonic() < deadline:
            pids = sorted(int(path.read_text()) for path in (self.base / "beats").glob("*.pid"))
            if len(pids) >= count:
                return pids
            assert self.process.poll() is None, "the producer exited before its children started"
            time.sleep(0.05)
        raise AssertionError(f"fewer than {count} children started")

    def lock_free(self) -> bool:
        try:
            lock = acquire(self.data_root, "producer")
        except CampaignLockBusy:
            return False
        lock.release()
        return True

    def finish(self, timeout: float = 120.0) -> tuple[int, dict[str, object]]:
        output, _ = self.process.communicate(timeout=timeout)
        lines = [line for line in output.splitlines() if line.startswith("SUMMARY ")]
        return self.process.returncode, json.loads(lines[-1][len("SUMMARY ") :]) if lines else {}

    def kill(self) -> None:
        if self.process.poll() is None:
            self.process.kill()
            self.process.wait(timeout=30)
        for pid_file in (self.base / "beats").glob("*.pid"):
            try:
                os.killpg(int(pid_file.read_text()), signal.SIGKILL)
            except (ProcessLookupError, PermissionError, ValueError):
                pass


def _alive(pid: int) -> bool:
    try:
        os.kill(pid, 0)
    except ProcessLookupError:
        return False
    return True


def _wait_dead(pids: list[int], timeout: float = 30.0) -> None:
    deadline = time.monotonic() + timeout
    while time.monotonic() < deadline:
        if not any(_alive(pid) for pid in pids):
            return
        time.sleep(0.05)
    raise AssertionError(f"children still alive: {[pid for pid in pids if _alive(pid)]}")


@pytest.fixture
def producers(tmp_path: Path, monkeypatch: pytest.MonkeyPatch):  # type: ignore[no-untyped-def]
    monkeypatch.delenv(TOKEN_ENV, raising=False)
    started: list[Producer] = []

    def start(mode: str, *, grace_seconds: float = 30.0) -> Producer:
        producer = Producer(tmp_path / f"run{len(started)}", mode, grace_seconds=grace_seconds)
        started.append(producer)
        return producer

    yield start
    for producer in started:
        producer.kill()


def _rows(producer: Producer) -> dict[str, tuple[str, str | None]]:
    with ProgressDB(producer.data_root / "progress.db") as db:
        rows = db._conn.execute("SELECT srr_id, state, error FROM samples WHERE species = 'apis'").fetchall()
    return {run: (state, error) for run, state, error in rows}


@pytest.mark.timeout(240)
def test_sigterm_stops_children_before_the_lock_is_released(producers) -> None:  # type: ignore[no-untyped-def]
    producer = producers("stubborn", grace_seconds=3.0)
    children = producer.child_pids(2)
    producer.process.send_signal(signal.SIGTERM)

    # The children ignore SIGTERM: until the grace ends and they are killed, the root stays locked.
    time.sleep(1.0)
    assert all(_alive(pid) for pid in children)
    assert not producer.lock_free()

    code, summary = producer.finish()

    assert code == EXIT_RETRYABLE
    assert summary == {"status": "stopped", "exit_code": EXIT_RETRYABLE, "stopped_by": "SIGTERM"}
    assert not any(_alive(pid) for pid in children)  # gone before the producer (and its lock) went
    assert is_quiescent(producer.data_root) and producer.lock_free()
    rows = _rows(producer)
    interrupted = [run for run, (state, error) in rows.items() if error == "interrupted: producer stopped by SIGTERM"]
    assert len(interrupted) >= 2  # their attempt never counts as a download failure
    assert all(state == "pending" for state, _ in rows.values())


@pytest.mark.timeout(240)
def test_a_second_sigint_kills_the_children_without_waiting_for_the_grace(  # type: ignore[no-untyped-def]
    producers,
) -> None:
    producer = producers("stubborn", grace_seconds=600.0)
    children = producer.child_pids(2)
    started = time.monotonic()
    producer.process.send_signal(signal.SIGINT)
    time.sleep(0.5)
    assert not producer.lock_free()
    producer.process.send_signal(signal.SIGINT)

    code, summary = producer.finish()

    assert time.monotonic() - started < 60  # not the 600 s grace
    assert code == EXIT_RETRYABLE and summary["stopped_by"] == "SIGINT"
    assert not any(_alive(pid) for pid in children)
    assert producer.lock_free()
    assert all(state == "pending" for state, _ in _rows(producer).values())


@pytest.mark.timeout(240)
def test_children_of_a_killed_producer_keep_the_root_locked_until_they_exit(  # type: ignore[no-untyped-def]
    producers,
) -> None:
    producer = producers("obedient")
    children = producer.child_pids(2)
    producer.process.kill()  # kill -9: no handler runs, nothing is stopped
    producer.process.wait(timeout=30)

    assert all(_alive(pid) for pid in children)
    # The children inherited the writer lock's descriptor: no second producer can start beside them.
    assert not producer.lock_free()

    for pid in children:
        os.killpg(pid, signal.SIGKILL)
    _wait_dead(children)
    deadline = time.monotonic() + 30
    while not producer.lock_free():
        assert time.monotonic() < deadline, "the lock outlived every process that held it"
        time.sleep(0.05)


def test_no_group_starts_after_shutdown_began(tmp_path: Path) -> None:
    process_groups.reset()
    try:
        assert process_groups.begin_shutdown("test") is True
        assert process_groups.begin_shutdown("again") is False
        with pytest.raises(process_groups.ProcessGroupsStopping):
            process_groups.spawn([sys.executable, "-c", "pass"])
        assert process_groups.registered() == ()
    finally:
        process_groups.reset()


def test_spawned_groups_inherit_the_named_descriptors_and_are_signalled(tmp_path: Path) -> None:
    marker = tmp_path / "inherited"
    read_end, write_end = os.pipe()
    process_groups.reset()
    process_groups.set_inherited_fds([write_end, None])
    try:
        child = process_groups.spawn(
            [
                sys.executable,
                "-c",
                f"import os, time; os.write({write_end}, b'x'); open({str(marker)!r}, 'w'); time.sleep(600)",
            ]
        )
        assert os.read(read_end, 1) == b"x"  # the descriptor number is valid in the child
        assert process_groups.registered() == (child.pid,)
        assert os.getpgid(child.pid) == child.pid  # its own session and group
        assert process_groups.signal_all(signal.SIGKILL) == [child.pid]
        child.wait(timeout=30)
        process_groups.unregister(child.pid)
        assert process_groups.wait_until_empty(0.1) is True
        assert marker.exists()
    finally:
        process_groups.set_inherited_fds(())
        process_groups.reset()
        os.close(read_end)
        os.close(write_end)


def test_stop_escalates_to_sigkill_after_the_grace(tmp_path: Path) -> None:
    process_groups.reset()
    try:
        child = process_groups.spawn(
            [
                sys.executable,
                "-c",
                "import signal, time; signal.signal(signal.SIGTERM, signal.SIG_IGN); time.sleep(600)",
            ]
        )
        time.sleep(0.5)  # let the child install its handler
        stop = ProducerStop(grace_seconds=0.5)
        stop.request(signal.SIGTERM)
        stop.stop_children()  # SIGTERM, ignored
        with pytest.raises(subprocess.TimeoutExpired):
            child.wait(timeout=0.3)
        time.sleep(0.5)
        stop.stop_children()  # the grace is over: SIGKILL
        assert child.wait(timeout=30) == -signal.SIGKILL
        process_groups.unregister(child.pid)
        assert stop.signal_name == "SIGTERM"
    finally:
        process_groups.reset()
