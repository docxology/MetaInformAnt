"""Producer startup, environment circuit breaker, run summary exit codes and run logging.

The campaigns are real frozen cohorts on real data roots under ``tmp_path``
(``Campaign`` from the cohort discovery tests). Discovery, the progress DB, the
lock and the scheduler are the production code; only the per-sample worker is
replaced by a small function that records what it saw and writes the outcome
a real worker would write. The environment failures are real ``OSError``s from
a read-only directory under the data root.
"""

from __future__ import annotations

import errno
import logging
import os
import stat
import subprocess
import sys
from pathlib import Path
from typing import Any

import pytest

import metainformant.rna.engine.streaming_orchestrator as orchestrator_module
from metainformant.rna.engine.campaign_lock import (
    EXIT_ENVIRONMENT,
    EXIT_INTERNAL,
    EXIT_LOCK_BUSY,
    EXIT_OK,
    EXIT_RETRYABLE,
    TOKEN_ENV,
    CampaignLockBusy,
    is_quiescent,
)
from metainformant.rna.engine.cohort import CohortError
from metainformant.rna.engine.preflight import PreflightError
from metainformant.rna.engine.progress_db import ProgressDB
from metainformant.rna.engine.streaming_orchestrator import (
    DiscoveryFailed,
    EnvironmentCircuitBreaker,
    StreamingPipelineOrchestrator,
    _write_text_durably,
    configure_producer_logging,
    environment_failure_class,
    producer_exit_code_for_exception,
)
from tests.rna.frozen_cohort_support import CONFIG_NAME, SPECIES
from tests.rna.test_cohort_discovery import Campaign

REPO_ROOT = Path(__file__).resolve().parents[2]
EPERM = "Unexpected sample task error: [Errno 1] Operation not permitted: '/data/test_species/work/x'"
KALLISTO = "Quantification Failed: kallisto exited 1"
TASKS = ["SRR101", "SRR102", "SRR103", "SRR105"]


def _stub_preflight(data_root: Path, **_kwargs: object) -> dict[str, str]:
    """The preflight has its own tests; here the host is taken as ready."""

    return {"data_root": str(data_root), "amalgkit_cli": "preflight-tested-elsewhere"}


@pytest.fixture
def campaign(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> Campaign:
    monkeypatch.delenv(TOKEN_ENV, raising=False)
    monkeypatch.delenv("AMALGKIT_ENV_FAILURE_TRIP", raising=False)
    monkeypatch.setattr(orchestrator_module, "run_campaign_preflight", _stub_preflight)
    return Campaign(tmp_path, monkeypatch)


def _errors(orchestrator: StreamingPipelineOrchestrator) -> dict[str, str | None]:
    connection = orchestrator.db._conn
    return dict(connection.execute("SELECT srr_id, error FROM samples WHERE species = ?", (SPECIES,)).fetchall())


def test_startup_resets_orphans_and_environment_failures_then_drains_with_exit_0(campaign: Campaign) -> None:
    orchestrator = campaign.orchestrator()
    db = orchestrator.db
    db.init_species(SPECIES, TASKS)
    db.set_state(SPECIES, "SRR101", "downloading")  # orphaned by a SIGKILLed producer seconds ago
    db.set_state(SPECIES, "SRR102", "quantifying")
    db.set_state(SPECIES, "SRR103", "failed", error=EPERM)  # the host's fault, not the sample's
    db.set_state(SPECIES, "SRR105", "failed", error=KALLISTO)  # the sample's own failure
    seen: dict[str, str | None] = {}

    def worker(srr: str, *_args: object, **_kwargs: object) -> dict[str, Any]:
        seen[srr] = db.get_state(SPECIES, srr)
        db.set_state(SPECIES, srr, "quantified")
        return {"quantified": True, "skipped": False, "error": None}

    orchestrator.process_single_sample = worker  # type: ignore[method-assign]
    summary = orchestrator.run_all([CONFIG_NAME], 10_000.0, workers=1, threads=1, discovery_workers=1)

    assert seen == {"SRR101": "pending", "SRR102": "pending", "SRR103": "pending", "SRR105": "failed"}
    assert summary["startup"]["inflight_reset"] == {"downloading": 1, "quantifying": 1}
    assert summary["startup"]["environment_failures_reset"] == {"environment_write_denied": 1}
    assert summary["exit_code"] == EXIT_OK and summary["status"] == "drained"
    assert summary["tasks_scheduled"] == 4 and summary["quantified_this_run"] == 4
    assert summary["remaining_by_state"] == {} and summary["cohort_id"] == campaign.cohort.cohort_id
    assert summary["breaker"]["tripped"] is False
    assert is_quiescent(campaign.data.base)
    db.close()


def test_sample_failures_leave_exit_75_until_the_operator_excludes_them(campaign: Campaign) -> None:
    orchestrator = campaign.orchestrator()
    db = orchestrator.db
    calls: list[str] = []

    def worker(srr: str, *_args: object, **_kwargs: object) -> dict[str, Any]:
        calls.append(srr)
        db.set_state(SPECIES, srr, "failed", error=KALLISTO)
        return {"quantified": False, "skipped": False, "error": KALLISTO}

    orchestrator.process_single_sample = worker  # type: ignore[method-assign]
    summary = orchestrator.run_all([CONFIG_NAME], 10_000.0, workers=1, threads=1, discovery_workers=1)
    assert sorted(calls) == TASKS  # deterministic sample failures never trip the environment breaker
    assert summary["exit_code"] == EXIT_RETRYABLE and summary["status"] == "retryable_work_remains"
    assert summary["remaining_by_state"] == {"failed": 4}
    assert summary["failed_by_class"] == {"quantification_failed": 4}
    assert summary["breaker"]["tripped"] is False and summary["breaker"]["environment_failures"] == 0

    db.record_exclusions(SPECIES, [{"srr_id": run, "reason_code": "permanent_drop"} for run in TASKS])
    calls.clear()
    drained = orchestrator.run_all([CONFIG_NAME], 10_000.0, workers=1, threads=1, discovery_workers=1)
    assert calls == [] and drained["tasks_scheduled"] == 0
    assert drained["exit_code"] == EXIT_OK
    db.close()


def test_environment_breaker_stops_admission_and_leaves_rows_pending(
    campaign: Campaign, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setenv("AMALGKIT_ENV_FAILURE_TRIP", "3")
    denied = campaign.data.base / "test_species" / "work" / "getfastq"
    denied.mkdir(parents=True, exist_ok=True)
    os.chmod(denied, stat.S_IRUSR | stat.S_IXUSR)
    orchestrator = campaign.orchestrator()
    calls: list[str] = []

    def worker(srr: str, *_args: object, **_kwargs: object) -> dict[str, Any]:
        calls.append(srr)
        with open(denied / f"{srr}.fastq.gz.part", "wb"):  # a real EACCES under the data root
            pass
        raise AssertionError("the directory is read-only")

    orchestrator.process_single_sample = worker  # type: ignore[method-assign]
    try:
        summary = orchestrator.run_all(
            [CONFIG_NAME], 10_000.0, workers=1, threads=1, max_in_flight=1, discovery_workers=1
        )
    finally:
        os.chmod(denied, stat.S_IRWXU)

    assert len(calls) == 3  # the fourth task was never admitted
    assert summary["exit_code"] == EXIT_ENVIRONMENT and summary["status"] == "environment_failure"
    assert summary["breaker"]["tripped"] is True and summary["breaker"]["consecutive"] == 3
    assert summary["environment_deferred"] == 3
    assert {orchestrator.db.get_state(SPECIES, run) for run in TASKS} == {"pending"}
    errors = _errors(orchestrator)
    for run in calls:
        assert errors[run].startswith("environment-deferred[environment_access_denied]: ")  # type: ignore[union-attr]
    assert summary["remaining_by_state"] == {"pending": 4}
    assert is_quiescent(campaign.data.base)
    orchestrator.db.close()


def test_environment_errors_reported_by_the_worker_are_deferred_too(
    campaign: Campaign, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setenv("AMALGKIT_ENV_FAILURE_TRIP", "2")
    orchestrator = campaign.orchestrator()
    db = orchestrator.db
    enxio = f"Quant exception batch 1: [Errno 6] Device not configured: '{campaign.data.base}/test_species/work'"

    def worker(srr: str, *_args: object, **_kwargs: object) -> dict[str, Any]:
        db.set_state(SPECIES, srr, "failed", error=enxio)  # what the real worker records before returning
        return {"quantified": False, "skipped": False, "error": enxio}

    orchestrator.process_single_sample = worker  # type: ignore[method-assign]
    summary = orchestrator.run_all([CONFIG_NAME], 10_000.0, workers=1, threads=1, max_in_flight=1)
    assert summary["exit_code"] == EXIT_ENVIRONMENT
    assert summary["environment_deferred"] == 2
    assert {db.get_state(SPECIES, run) for run in TASKS} == {"pending"}
    assert summary["failed_by_class"] == {}
    db.close()


def test_environment_failure_classes() -> None:
    root = Path("/data/root")
    assert environment_failure_class(PermissionError(errno.EPERM, "denied")) == "environment_write_denied"
    assert environment_failure_class(OSError(errno.ENXIO, "gone", "/data/root/x"), root) == (
        "environment_device_unavailable"
    )
    assert environment_failure_class(OSError(errno.EROFS, "ro", "/data/root/x"), root) == "environment_read_only"
    # A permission problem outside the data root is the sample's (or the tool's) problem, not the volume's.
    assert environment_failure_class(PermissionError(errno.EACCES, "denied", "/elsewhere/x"), root) is None
    assert environment_failure_class(OSError(errno.ENOENT, "missing", "/data/root/x"), root) is None
    assert environment_failure_class(EPERM) == "environment_write_denied"
    assert environment_failure_class("[Errno 2] No such file or directory: 'amalgkit'") == "environment_missing_tool"
    assert environment_failure_class("[Errno 2] No such file or directory: '/data/root/SRR1.fastq'") is None
    assert environment_failure_class(KALLISTO) is None and environment_failure_class(None) is None


def test_breaker_counts_only_consecutive_environment_failures() -> None:
    breaker = EnvironmentCircuitBreaker(3)
    assert not breaker.record_environment_failure("environment_write_denied", "a")
    assert not breaker.record_environment_failure("environment_write_denied", "b")
    breaker.record_other_outcome()
    assert not breaker.record_environment_failure("environment_write_denied", "c")
    assert not breaker.record_environment_failure("environment_write_denied", "d")
    assert breaker.record_environment_failure("environment_write_denied", "e")  # trips exactly once
    assert not breaker.record_environment_failure("environment_write_denied", "f")
    state = breaker.state()
    assert state["tripped"] and state["environment_failures"] == 6 and len(state["recent"]) == 3
    with pytest.raises(ValueError):
        EnvironmentCircuitBreaker(0)


def test_exit_codes_for_exceptions() -> None:
    assert producer_exit_code_for_exception(CampaignLockBusy("busy")) == EXIT_LOCK_BUSY == 73
    assert producer_exit_code_for_exception(PreflightError("no tools")) == EXIT_ENVIRONMENT == 69
    assert producer_exit_code_for_exception(CohortError("bytes differ")) == EXIT_ENVIRONMENT
    assert producer_exit_code_for_exception(DiscoveryFailed("one species")) == EXIT_ENVIRONMENT
    assert producer_exit_code_for_exception(OSError(errno.ENXIO, "gone")) == EXIT_ENVIRONMENT
    assert producer_exit_code_for_exception(KeyboardInterrupt()) == EXIT_RETRYABLE == 75
    assert producer_exit_code_for_exception(ValueError("bug")) == EXIT_INTERNAL == 70


def _script(*args: str, path: str | None = None) -> subprocess.CompletedProcess[str]:
    environment = os.environ.copy()
    environment.pop(TOKEN_ENV, None)
    environment["PYTHONPATH"] = str(REPO_ROOT / "src")
    environment["AMALGKIT_MIN_EXTERNAL_FREE_GB"] = "0"
    environment["AMALGKIT_MIN_SYSTEM_FREE_GB"] = "0"
    if path is not None:
        environment["PATH"] = path
    return subprocess.run(
        [sys.executable, str(REPO_ROOT / "scripts" / "rna" / "run_all_species.py"), *args],
        cwd=REPO_ROOT,
        env=environment,
        capture_output=True,
        text=True,
        timeout=300,
        check=False,
    )


def test_run_all_species_exit_codes_for_busy_root_and_broken_host(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    campaign = Campaign(tmp_path, monkeypatch, forbid_subprocess=False)
    common = [
        "--config-dir",
        str(campaign.config_dir),
        "--data-root",
        str(campaign.data.base),
        "--cohort",
        str(campaign.base / "cohort"),
        "--species",
        SPECIES,
    ]
    holder_code = (
        "import sys, time\n"
        "from metainformant.rna.engine.campaign_lock import acquire\n"
        "lock = acquire(sys.argv[1], 'import')\n"
        "print('held', flush=True)\n"
        "time.sleep(600)\n"
    )
    environment = os.environ.copy()
    environment.pop(TOKEN_ENV, None)
    environment["PYTHONPATH"] = str(REPO_ROOT / "src")
    holder = subprocess.Popen(
        [sys.executable, "-c", holder_code, str(campaign.data.base)], env=environment, stdout=subprocess.PIPE, text=True
    )
    try:
        assert holder.stdout is not None and holder.stdout.readline().strip() == "held"
        busy = _script(*common)
        assert busy.returncode == EXIT_LOCK_BUSY, busy.stderr
        assert "Another writer holds the data root" in busy.stderr and "kind=import" in busy.stderr
        assert not (campaign.data.base / "logs").exists()  # refused before anything was written
    finally:
        holder.kill()
        holder.wait(timeout=30)
        assert holder.stdout is not None
        holder.stdout.close()

    empty_bin = tmp_path / "empty_bin"
    empty_bin.mkdir()
    broken = _script(*common, path=str(empty_bin))
    assert broken.returncode == EXIT_ENVIRONMENT, broken.stdout + broken.stderr
    assert "PreflightError" in broken.stderr
    logs = sorted((campaign.data.base / "logs").glob("run_all_species_*.log"))
    assert len(logs) == 1
    text = logs[0].read_text()
    assert "Holding the campaign writer lock" in text and "Campaign preflight failed" in text
    assert is_quiescent(campaign.data.base)  # the lock went with the process


def test_run_all_species_exits_0_when_every_task_is_already_quantified(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    campaign = Campaign(tmp_path, monkeypatch, forbid_subprocess=False)
    for run in TASKS:
        campaign.data.write_sample(run, config=campaign.yaml)  # reusable under the cohort contract
    tools = tmp_path / "tools"
    tools.mkdir()
    for tool in ("amalgkit", "kallisto", "fasterq-dump", "prefetch", "vdb-validate", "pigz"):
        script = tools / tool
        body = 'if [ "$1" = version ]; then echo "kallisto, version 0.52.0"; fi' if tool == "kallisto" else "exit 0"
        script.write_text(f"#!/bin/sh\n{body}\n")
        script.chmod(0o755)
    result = _script(
        "--config-dir",
        str(campaign.config_dir),
        "--data-root",
        str(campaign.data.base),
        "--cohort",
        str(campaign.base / "cohort"),
        "--species",
        SPECIES,
        path=f"{tools}{os.pathsep}/usr/bin{os.pathsep}/bin",
    )
    assert result.returncode == EXIT_OK, result.stdout + result.stderr
    assert "Producer finished: drained (exit 0)" in result.stdout
    log = next((campaign.data.base / "logs").glob("run_all_species_*.log")).read_text()
    assert "Campaign preflight passed" in log and "kallisto_version=0.52.0" in log
    assert "Producer summary" in log
    with ProgressDB(campaign.data.base / "pipeline_progress.db") as db:
        assert {db.get_state(SPECIES, run) for run in TASKS} == {"quantified"}
    assert is_quiescent(campaign.data.base)


def test_producer_logging_writes_info_to_a_per_run_file(tmp_path: Path) -> None:
    root = logging.getLogger()
    level, handlers = root.level, list(root.handlers)
    try:
        path = configure_producer_logging(tmp_path / "logs", prefix="run_all_species")
        assert path is not None and path.parent == tmp_path / "logs" and path.name.startswith("run_all_species_")
        assert root.level == logging.INFO
        orchestrator_module.logger.info("progress line for the run log")
        for handler in root.handlers:
            handler.flush()
        assert "progress line for the run log" in path.read_text()
    finally:
        for handler in list(root.handlers):
            if handler not in handlers:
                root.removeHandler(handler)
                handler.close()
        root.setLevel(level)


def test_default_stall_timeout_is_valid_and_silent(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.delenv("AMALGKIT_PIPELINE_QUANT_STALL_TIMEOUT_SECONDS", raising=False)
    records: list[logging.LogRecord] = []

    class Collect(logging.Handler):
        def emit(self, record: logging.LogRecord) -> None:
            records.append(record)

    handler = Collect(level=logging.WARNING)
    orchestrator_module.logger.addHandler(handler)
    try:
        orchestrator = StreamingPipelineOrchestrator(tmp_path / "config", tmp_path / "logs", db_path=tmp_path / "p.db")
    finally:
        orchestrator_module.logger.removeHandler(handler)
    assert orchestrator.quant_stall_timeout_seconds == 0
    assert not [record for record in records if "STALL" in record.getMessage()]
    orchestrator.db.close()


def test_durable_text_write_replaces_atomically_and_cleans_up(tmp_path: Path) -> None:
    destination = tmp_path / "marker.json"
    _write_text_durably(destination, '{"a": 1}\n', tmp_path / ".marker.json.tmp")
    _write_text_durably(destination, '{"a": 2}\n', tmp_path / ".marker.json.tmp")
    assert destination.read_text() == '{"a": 2}\n'
    assert sorted(path.name for path in tmp_path.iterdir()) == ["marker.json"]
    with pytest.raises(IsADirectoryError):
        _write_text_durably(tmp_path, "x", tmp_path / ".dir.tmp")
    assert sorted(path.name for path in tmp_path.iterdir()) == ["marker.json"]
