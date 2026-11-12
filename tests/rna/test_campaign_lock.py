"""The campaign writer lock: one writer per data root, released by the kernel on any death.

Every test uses real processes, real ``flock`` locks and real directories
under ``tmp_path``; holders are separate Python processes that are killed with
SIGKILL exactly as a crashed or force-quit producer would be.
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

import metainformant.rna.engine.streaming_orchestrator as orchestrator_module
from metainformant.rna.engine.campaign_lock import (
    EXIT_ENVIRONMENT,
    EXIT_LOCK_BUSY,
    LOCK_FILENAME,
    LOCK_SCHEMA,
    TOKEN_ENV,
    CampaignLockBusy,
    acquire,
    acquire_writer,
    is_quiescent,
    main,
    read_owner,
    writer_status,
)
from metainformant.rna.engine.completion_inventory import freeze_inventory
from metainformant.rna.engine.streaming_orchestrator import (
    StreamingPipelineOrchestrator,
    producer_exit_code_for_exception,
)

REPO_ROOT = Path(__file__).resolve().parents[2]

HOLDER = """
import sys, time
from metainformant.rna.engine.campaign_lock import acquire
lock = acquire(sys.argv[1], sys.argv[2])
print("held", lock.token, flush=True)
time.sleep(600)
"""

INSIDE_CHILD = """
import os, sys, time
from metainformant.rna.engine.campaign_lock import acquire, writer_status
lock = acquire(sys.argv[1], "producer")
status = writer_status(sys.argv[1])
print("inside", os.getpid(), lock.inherited, status["inside_lock"], status["quiescent"], flush=True)
time.sleep(600)
"""

ACQUISITION_FLOCK = """
import fcntl, sys, time
handle = open(sys.argv[1], "a")
fcntl.flock(handle, fcntl.LOCK_EX | fcntl.LOCK_NB)
print("held", flush=True)
time.sleep(600)
"""


def _environment() -> dict[str, str]:
    environment = os.environ.copy()
    environment.pop(TOKEN_ENV, None)
    source = str(REPO_ROOT / "src")
    environment["PYTHONPATH"] = os.pathsep.join(filter(None, [source, environment.get("PYTHONPATH")]))
    return environment


def _python(code: str, *args: str) -> subprocess.Popen[str]:
    return subprocess.Popen(
        [sys.executable, "-c", code, *args], env=_environment(), stdout=subprocess.PIPE, text=True, cwd=REPO_ROOT
    )


def _first_line(process: subprocess.Popen[str]) -> str:
    assert process.stdout is not None
    line = process.stdout.readline().strip()
    assert line, f"child exited before reporting (exit {process.poll()})"
    return line


def _lock_cli(*args: str) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [sys.executable, "-m", "metainformant.rna.engine.campaign_lock", *args],
        env=_environment(),
        capture_output=True,
        text=True,
        timeout=120,
        check=False,
        cwd=REPO_ROOT,
    )


def _kill(process: subprocess.Popen[str]) -> None:
    if process.poll() is None:
        process.kill()
    process.wait(timeout=30)
    if process.stdout is not None:
        process.stdout.close()


def _wait_gone(pid: int, timeout: float = 30.0) -> None:
    deadline = time.monotonic() + timeout
    while time.monotonic() < deadline:
        try:
            os.kill(pid, 0)
        except ProcessLookupError:
            return
        time.sleep(0.05)
    raise AssertionError(f"process {pid} did not exit")


@pytest.fixture
def data_root(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> Path:
    root = tmp_path / "data"
    root.mkdir()
    monkeypatch.delenv(TOKEN_ENV, raising=False)
    return root


def test_owner_record_nesting_and_release(data_root: Path) -> None:
    with acquire_writer(data_root, "producer") as outer:
        owner = read_owner(data_root)
        assert owner is not None
        assert owner["schema"] == LOCK_SCHEMA and owner["kind"] == "producer"
        assert owner["pid"] == os.getpid() and owner["host"] and owner["started_at"].endswith("Z")
        assert owner["token"] == outer.token and len(outer.token) == 32
        with acquire_writer(data_root, "producer") as inner:  # same kind nests
            assert inner.inherited and inner.token == outer.token
        with pytest.raises(CampaignLockBusy, match="refusing a nested 'import' writer"):
            acquire(data_root, "import")
        status = writer_status(data_root)
        assert status["held"] and status["inside_lock"] and status["quiescent"]
        # The flock survived the inner release: another process is still refused.
        busy = _lock_cli("status", "--data-root", str(data_root))
        assert busy.returncode == EXIT_LOCK_BUSY
        assert json.loads(busy.stdout)["held"] is True
    free = _lock_cli("status", "--data-root", str(data_root))
    assert free.returncode == 0, free.stderr
    assert json.loads(free.stdout)["quiescent"] is True
    assert (data_root / LOCK_FILENAME).is_file()  # the last owner record stays for humans


def test_missing_root_is_never_created(tmp_path: Path) -> None:
    absent = tmp_path / "unmounted"
    with pytest.raises(FileNotFoundError):
        acquire(absent, "producer")
    assert main(["run", "--data-root", str(absent), "--kind", "producer", "--", "true"]) == EXIT_ENVIRONMENT
    assert not absent.exists()
    with pytest.raises(ValueError, match="short identifier"):
        acquire(tmp_path, "two words")


def test_second_run_all_refuses_while_a_holder_lives_and_runs_after_its_sigkill(
    data_root: Path, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setenv("AMALGKIT_DATA_ROOT", str(data_root))
    preflights: list[Path] = []

    def record_preflight(root: Path, **_kwargs: object) -> dict[str, str]:
        preflights.append(Path(root))
        return {"data_root": str(root), "amalgkit_cli": "recorded"}

    monkeypatch.setattr(orchestrator_module, "run_campaign_preflight", record_preflight)
    orchestrator = StreamingPipelineOrchestrator(
        tmp_path / "config", tmp_path / "logs", db_path=data_root / "pipeline_progress.db"
    )
    orchestrator.db.init_species("test_species", ["SRR1"])
    orchestrator.db.set_state("test_species", "SRR1", "downloading")

    holder = _python(HOLDER, str(data_root), "producer")
    try:
        assert _first_line(holder).startswith("held ")
        with pytest.raises(CampaignLockBusy, match="locked by another writer") as excinfo:
            orchestrator.run_all(["amalgkit_definitely_absent.yaml"], 1.0, 1, 1)
        assert excinfo.value.holder is not None and excinfo.value.holder["pid"] == holder.pid
        assert producer_exit_code_for_exception(excinfo.value) == EXIT_LOCK_BUSY
        # Refused before the preflight and before any state change.
        assert preflights == []
        assert orchestrator.db.get_state("test_species", "SRR1") == "downloading"
        assert not is_quiescent(data_root)
    finally:
        holder.send_signal(signal.SIGKILL)
        _kill(holder)

    summary = orchestrator.run_all(["amalgkit_definitely_absent.yaml"], 1.0, 1, 1)
    assert summary["exit_code"] == 0 and summary["status"] == "drained"
    assert summary["startup"]["inflight_reset"] == {"downloading": 1}
    assert orchestrator.db.get_state("test_species", "SRR1") == "pending"
    assert [path.resolve() for path in preflights] == [data_root.resolve()]
    owner = read_owner(data_root)
    assert owner is not None and owner["pid"] == os.getpid() and owner["kind"] == "producer"
    assert is_quiescent(data_root)
    orchestrator.db.close()


def test_run_wrapper_hands_the_lock_to_its_child(data_root: Path) -> None:
    wrapper = subprocess.Popen(
        [
            sys.executable,
            "-m",
            "metainformant.rna.engine.campaign_lock",
            "run",
            "--data-root",
            str(data_root),
            "--kind",
            "producer",
            "--",
            sys.executable,
            "-c",
            INSIDE_CHILD,
            str(data_root),
        ],
        env=_environment(),
        stdout=subprocess.PIPE,
        text=True,
        cwd=REPO_ROOT,
    )
    child_pid = None
    try:
        words = _first_line(wrapper).split()
        assert words[0] == "inside"
        child_pid = int(words[1])
        # The child carries the holder token: it is inside the lock and sees the root as quiescent.
        assert words[2:] == ["True", "True", "True"]
        assert read_owner(data_root)["argv"][0] == sys.executable  # type: ignore[index]
        with pytest.raises(CampaignLockBusy):
            acquire(data_root, "import")
        assert _lock_cli("run", "--data-root", str(data_root), "--kind", "import", "--", "true").returncode == (
            EXIT_LOCK_BUSY
        )
        # SIGKILL of the wrapper alone: the child still holds the passed descriptor and keeps writing safely.
        wrapper.send_signal(signal.SIGKILL)
        wrapper.wait(timeout=30)
        with pytest.raises(CampaignLockBusy):
            acquire(data_root, "import")
        os.kill(child_pid, signal.SIGKILL)
        _wait_gone(child_pid)
        child_pid = None
        with acquire_writer(data_root, "import") as lock:
            assert not lock.inherited
    finally:
        _kill(wrapper)
        if child_pid is not None:
            os.kill(child_pid, signal.SIGKILL)


def test_run_wrapper_propagates_the_child_exit_code(data_root: Path) -> None:
    result = _lock_cli(
        "run", "--data-root", str(data_root), "--kind", "seal", "--", sys.executable, "-c", "raise SystemExit(5)"
    )
    assert result.returncode == 5, result.stderr
    assert is_quiescent(data_root)


def test_live_launcher_directory_lock_blocks_and_stale_or_ancestor_ones_do_not(data_root: Path) -> None:
    lock_dir = data_root / "results" / ".full_campaign.lock"
    lock_dir.mkdir(parents=True)
    sleeper = _python("import time; print('up', flush=True); time.sleep(600)")
    try:
        _first_line(sleeper)
        (lock_dir / "pid").write_text(f"{sleeper.pid}\n")
        with pytest.raises(CampaignLockBusy, match=r"full_campaign\.lock \(live launcher process, pid"):
            acquire(data_root, "producer")
        assert not (data_root / LOCK_FILENAME).exists()  # refused before anything was created
        status = writer_status(data_root)
        assert not status["quiescent"] and status["other_writers"][0]["pid"] == sleeper.pid
        assert _lock_cli("status", "--data-root", str(data_root)).returncode == EXIT_LOCK_BUSY
    finally:
        _kill(sleeper)
    # A dead launcher's leftover directory is not a writer.
    with acquire_writer(data_root, "producer"):
        pass
    # The launcher that started this producer (an ancestor) is not another writer.
    (lock_dir / "pid").write_text(f"{os.getppid()}\n")
    with acquire_writer(data_root, "producer"):
        pass
    # Between the launcher's mkdir and its pid write: busy for a grace period, then ignored.
    (lock_dir / "pid").unlink()
    os.utime(lock_dir, None)
    with pytest.raises(CampaignLockBusy, match="without a pid file"):
        acquire(data_root, "producer")
    old = time.time() - 3600
    os.utime(lock_dir, (old, old))
    with acquire_writer(data_root, "producer"):
        pass


def test_held_acquisition_worker_flock_blocks_the_root(data_root: Path) -> None:
    worker = _python(ACQUISITION_FLOCK, str(data_root / ".acquisition.lock"))
    try:
        assert _first_line(worker) == "held"
        with pytest.raises(CampaignLockBusy, match=r"\.acquisition\.lock"):
            acquire(data_root, "producer")
        assert not is_quiescent(data_root)
    finally:
        _kill(worker)
    assert (data_root / ".acquisition.lock").exists()  # an unheld lock file is not a writer
    assert is_quiescent(data_root)
    with acquire_writer(data_root, "producer"):
        pass


def test_freeze_inventory_refuses_while_another_writer_holds_the_root(data_root: Path, tmp_path: Path) -> None:
    holder = _python(HOLDER, str(data_root), "producer")
    try:
        _first_line(holder)
        with pytest.raises(RuntimeError, match="another writer holds the data root"):
            freeze_inventory(data_root, tmp_path / "config", tmp_path / "snapshot", expected_species_count=None)
        assert not (tmp_path / "snapshot").exists()
    finally:
        _kill(holder)
    # The former check looked for ``<root>/.full_campaign.lock``; the launcher's lock lives under results/.
    (data_root / "results" / ".full_campaign.lock").mkdir(parents=True)
    with pytest.raises(RuntimeError, match="another writer holds the data root"):
        freeze_inventory(data_root, tmp_path / "config", tmp_path / "snapshot", expected_species_count=None)
