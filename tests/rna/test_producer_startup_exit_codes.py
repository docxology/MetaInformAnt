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
import shutil
import sqlite3
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
    EXIT_NEEDS_OPERATOR,
    EXIT_OK,
    EXIT_RETRYABLE,
    TOKEN_ENV,
    CampaignLockBusy,
    is_quiescent,
)
from metainformant.rna.engine.cohort import CohortError
from metainformant.rna.engine.preflight import PreflightError
from metainformant.rna.engine.progress_db import NON_ENVIRONMENT_OS_ERROR, ProgressDB, classify_sample_error
from metainformant.rna.engine.provenance import QUANT_PROVENANCE_FILENAME, REFERENCE_BINDING_FILENAME
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
    assert "needs_operator={} unverified_outputs=0" in result.stdout  # the launcher copies these counts
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
    # A failed rename (here onto a directory: EISDIR on Linux and macOS) leaves no temporary behind.
    occupied = tmp_path / "occupied"
    occupied.mkdir()
    with pytest.raises(OSError):
        _write_text_durably(occupied, "x", tmp_path / ".dir.tmp")
    assert sorted(path.name for path in tmp_path.iterdir()) == ["marker.json", "occupied"]
    assert not any(occupied.iterdir())


# ---------------------------------------------------------------------------
# Cohort completeness, requantification, unreadable outputs, error mapping
# ---------------------------------------------------------------------------


def _quantify(campaign: Campaign, db: ProgressDB, calls: list[str]) -> Any:
    """A worker that writes what a real quantification leaves: a cohort-valid output and the quantified row."""

    def worker(srr: str, *_args: object, **_kwargs: object) -> dict[str, Any]:
        calls.append(srr)
        if not campaign.data.sample(srr).exists():
            campaign.data.write_sample(srr, config=campaign.yaml)
        db.set_state(SPECIES, srr, "quantified")
        return {"quantified": True, "skipped": False, "error": None}

    return worker


def test_rows_only_the_operator_can_move_exit_65_with_counts_instead_of_drained(campaign: Campaign) -> None:
    campaign.data.write_sample("SRR102", config=campaign.data.config)  # other configuration bytes: quarantined
    orchestrator = campaign.orchestrator()
    db = orchestrator.db
    db.init_species(SPECIES, ["SRR102"])
    db.set_state(SPECIES, "SRR102", "quantified")
    calls: list[str] = []
    orchestrator.process_single_sample = _quantify(campaign, db, calls)  # type: ignore[method-assign]

    # In cohort mode the reuse dry-run is on by default: the only recorded output would be quarantined, so the
    # producer refuses before discovery changes anything.
    with pytest.raises(PreflightError, match="1 of 1 sampled quantified outputs would be quarantined"):
        orchestrator.run_all([CONFIG_NAME], 1000.0, workers=1, threads=1, discovery_workers=1)
    assert calls == [] and db.get_state(SPECIES, "SRR102") == "quantified"

    # The operator has looked at it and opts out (as AMALGKIT_PREFLIGHT_REUSE_DRY_RUN=0 does).
    # SRR105 has 5,000 Gbases: above --max-gb 1000 it is never scheduled.
    summary = orchestrator.run_all(
        [CONFIG_NAME], 1000.0, workers=1, threads=1, discovery_workers=1, reuse_dry_run_fraction=0
    )
    assert summary["startup"]["reuse_dry_run_fraction"] == 0

    assert sorted(calls) == ["SRR101", "SRR103"]
    assert summary["exit_code"] == EXIT_NEEDS_OPERATOR and summary["status"] == "needs_operator"
    assert summary["needs_operator"] == {"deferred_max_gb": 1, "quarantined": 1}
    assert summary["remaining_by_state"] == {}

    # The operator raises the budget and drops the quarantined run: the cohort is complete.
    db.record_exclusions(SPECIES, [{"srr_id": "SRR102", "reason_code": "permanent_drop"}])
    drained = orchestrator.run_all([CONFIG_NAME], 10_000.0, workers=1, threads=1, discovery_workers=1)
    assert sorted(calls) == ["SRR101", "SRR103", "SRR105"]
    assert drained["exit_code"] == EXIT_OK and drained["needs_operator"] == {}
    db.close()


@pytest.mark.parametrize(
    ("explicit", "env", "cohort_mode", "expected"),
    [
        (None, None, True, 0.01),  # cohort producer: on by default
        (None, None, False, 0.0),  # legacy lane: off by default
        (None, "0", True, 0.0),  # explicit opt-out
        (None, "0.05", False, 0.05),
        (None, "not-a-number", True, 0.01),
        (None, "-1", True, 0.01),
        (0.0, "0.5", True, 0.0),  # the flag wins over the environment
        (0.2, None, False, 0.2),
    ],
)
def test_reuse_dry_run_fraction_defaults_on_in_cohort_mode(
    monkeypatch: pytest.MonkeyPatch, explicit: float | None, env: str | None, cohort_mode: bool, expected: float
) -> None:
    if env is None:
        monkeypatch.delenv("AMALGKIT_PREFLIGHT_REUSE_DRY_RUN", raising=False)
    else:
        monkeypatch.setenv("AMALGKIT_PREFLIGHT_REUSE_DRY_RUN", env)
    assert orchestrator_module.reuse_dry_run_fraction_setting(explicit, cohort_mode=cohort_mode) == expected


def test_a_cohort_producer_runs_the_reuse_dry_run_by_default(campaign: Campaign) -> None:
    for run in TASKS:
        campaign.data.write_sample(run, config=campaign.yaml)
    orchestrator = campaign.orchestrator()
    db = orchestrator.db
    db.init_species(SPECIES, TASKS)
    for run in TASKS:
        db.set_state(SPECIES, run, "quantified")
    calls: list[str] = []
    orchestrator.process_single_sample = _quantify(campaign, db, calls)  # type: ignore[method-assign]

    summary = orchestrator.run_all([CONFIG_NAME], 10_000.0, workers=1, threads=1, discovery_workers=1)

    assert summary["startup"]["reuse_dry_run_fraction"] == 0.01
    assert summary["startup"]["reuse_dry_run"][SPECIES]["sampled"] == 1  # ceil(1% of 4)
    assert summary["startup"]["reuse_dry_run"][SPECIES]["would_quarantine"] == 0
    assert summary["exit_code"] == EXIT_OK and calls == []
    db.close()


def _rewrite_output(campaign: Campaign, srr: str) -> None:
    sample = campaign.data.sample(srr)
    shutil.rmtree(sample)
    campaign.data.write_sample(srr, config=campaign.yaml, counts_offset=1)


def test_requantification_policy_all_requantifies_then_drains(campaign: Campaign) -> None:
    for run in TASKS:
        campaign.data.write_sample(run, config=campaign.yaml)
    orchestrator = campaign.orchestrator()
    db = orchestrator.db
    quantified: list[str] = []

    def download(srr: str, *_args: object, **_kwargs: object) -> bool:
        return True

    def quant(_config: Path, _batch: int, _species: str, _threads: int, srr: str, *_a: object, **_k: object) -> Any:
        quantified.append(srr)
        _rewrite_output(campaign, srr)
        return True, None

    # The real per-sample worker runs; only the network transfer and the kallisto run are replaced.
    orchestrator.download_fastq = download  # type: ignore[method-assign]
    orchestrator.quant_sample = quant  # type: ignore[method-assign]

    first = orchestrator.run_all(
        [CONFIG_NAME], 10_000.0, workers=1, threads=1, discovery_workers=1, requantification_policy="all"
    )
    assert sorted(quantified) == TASKS  # every verified output was quantified again, none skipped
    assert first["exit_code"] == EXIT_OK and first["status"] == "drained"
    assert all(db.get_state(SPECIES, run) == "quantified" for run in TASKS)

    quantified.clear()
    preserved = orchestrator.run_all(
        [CONFIG_NAME], 10_000.0, workers=1, threads=1, discovery_workers=1, requantification_policy="preserve"
    )
    assert quantified == [] and preserved["exit_code"] == EXIT_OK
    db.close()


def test_a_reusable_output_found_by_the_worker_is_recorded_quantified(campaign: Campaign) -> None:
    orchestrator = campaign.orchestrator()
    db = orchestrator.db
    orchestrator.discover_species_tasks(CONFIG_NAME, max_gb=10_000.0, threads=1)  # registers the reuse contract
    assert db.get_state(SPECIES, "SRR101") == "pending"
    campaign.data.write_sample("SRR101", config=campaign.yaml)  # appeared after discovery (an import, a peer)

    result = orchestrator.process_single_sample("SRR101", 2, campaign.data.work / "getfastq", campaign.yaml, SPECIES, 1)

    assert result["skipped"] is True and db.get_state(SPECIES, "SRR101") == "quantified"
    db.close()


@pytest.mark.skipif(os.geteuid() == 0, reason="root reads files whatever their mode")
def test_an_unreadable_output_leaves_the_run_retryable_and_is_never_quarantined(campaign: Campaign) -> None:
    for run in TASKS:
        campaign.data.write_sample(run, config=campaign.yaml)
    orchestrator = campaign.orchestrator()
    db = orchestrator.db
    db.init_species(SPECIES, TASKS)
    for run in TASKS:
        db.set_state(SPECIES, run, "quantified")  # recorded before any fingerprinted audit, like the canonical DB
    calls: list[str] = []
    orchestrator.process_single_sample = _quantify(campaign, db, calls)  # type: ignore[method-assign]
    sidecar = campaign.data.sample("SRR102") / QUANT_PROVENANCE_FILENAME
    os.chmod(sidecar, 0)
    try:
        flaky = orchestrator.run_all([CONFIG_NAME], 10_000.0, workers=1, threads=1, discovery_workers=1)
    finally:
        os.chmod(sidecar, 0o644)

    assert flaky["exit_code"] == EXIT_RETRYABLE and flaky["unverified_outputs"] == 1
    assert db.get_state(SPECIES, "SRR102") == "quantified" and calls == []
    readable = orchestrator.run_all([CONFIG_NAME], 10_000.0, workers=1, threads=1, discovery_workers=1)
    assert readable["exit_code"] == EXIT_OK and readable["unverified_outputs"] == 0
    assert db.get_state(SPECIES, "SRR102") == "quantified"
    db.close()


def test_an_os_error_outside_the_data_root_is_counted_and_never_reset_at_startup(campaign: Campaign) -> None:
    orchestrator = campaign.orchestrator()
    db = orchestrator.db
    calls: list[str] = []
    outside = "/nonexistent-outside-root/SRR103_1.fastq.gz"

    quantify = _quantify(campaign, db, calls)

    def worker(srr: str, *args: object, **kwargs: object) -> dict[str, Any]:
        if srr == "SRR103":
            calls.append(srr)
            raise PermissionError(errno.EPERM, "Operation not permitted", outside)
        return quantify(srr, *args, **kwargs)

    orchestrator.process_single_sample = worker  # type: ignore[method-assign]
    first = orchestrator.run_all([CONFIG_NAME], 10_000.0, workers=1, threads=1, discovery_workers=1)
    assert first["breaker"]["environment_failures"] == 0  # the producer said: not the host
    record = db.retry_records(SPECIES)[(SPECIES, "SRR103")]
    assert record.attempts == 1 and record.last_error_class == NON_ENVIRONMENT_OS_ERROR

    calls.clear()
    second = orchestrator.run_all([CONFIG_NAME], 10_000.0, workers=1, threads=1, discovery_workers=1)
    # The DB's own text classifier would call the message environmental; the recorded verdict wins.
    assert second["startup"]["environment_failures_reset"] == {}
    assert calls == [] and second["deferred_backoff"] == 1 and second["exit_code"] == EXIT_RETRYABLE
    db.close()


def test_a_missing_file_path_is_not_a_missing_tool_in_either_classifier(tmp_path: Path) -> None:
    path_missing = f"Unexpected sample task error: [Errno 2] No such file or directory: '{tmp_path}/SRR1_1.fastq.gz'"
    tool_missing = "Unexpected sample task error: [Errno 2] No such file or directory: 'kallisto'"
    assert environment_failure_class(path_missing, tmp_path) is None
    assert classify_sample_error(path_missing) != "environment_missing_tool"
    assert environment_failure_class(tool_missing) == classify_sample_error(tool_missing) == "environment_missing_tool"


@pytest.mark.parametrize(
    ("message", "expected"),
    [
        ("disk I/O error", EXIT_ENVIRONMENT),
        ("unable to open database file", EXIT_ENVIRONMENT),
        ("attempt to write a readonly database", EXIT_ENVIRONMENT),
        ("database or disk is full", EXIT_ENVIRONMENT),
        ("no such table: samples", EXIT_INTERNAL),
    ],
)
def test_a_lost_database_volume_maps_to_69_not_70(message: str, expected: int) -> None:
    assert producer_exit_code_for_exception(sqlite3.OperationalError(message)) == expected


def test_a_produced_cohort_output_gets_a_binding_witness(campaign: Campaign) -> None:
    from metainformant.rna.engine.provenance import bind_produced_output, read_reference_binding

    orchestrator = campaign.orchestrator()
    orchestrator.discover_species_tasks(CONFIG_NAME, max_gb=10_000.0, threads=1)
    context = orchestrator._reuse_context(SPECIES)
    assert context is not None and context.is_cohort
    sample = campaign.data.write_sample("SRR101", config=campaign.yaml)

    bind_produced_output(sample, context.expectation("SRR101"), index_digests=context.index_digests)

    binding = read_reference_binding(sample)
    assert binding is not None and binding["cohort_id"] == campaign.cohort.cohort_id
    assert binding["evidence"] == {"kind": "produced", "verified_by": "reference_manifest"}
    # Reuse no longer depends on the reference manifest's bytes.
    campaign.data.manifest.write_text(campaign.data.manifest.read_text() + "\n")
    assert context.classify(sample, "SRR101")["status"] == "current"
    assert (sample / REFERENCE_BINDING_FILENAME).is_file()
    orchestrator.db.close()


def test_an_output_that_does_not_verify_gets_no_witness(campaign: Campaign) -> None:
    from metainformant.rna.engine.provenance import bind_produced_output

    orchestrator = campaign.orchestrator()
    orchestrator.discover_species_tasks(CONFIG_NAME, max_gb=10_000.0, threads=1)
    context = orchestrator._reuse_context(SPECIES)
    assert context is not None
    sample = campaign.data.write_sample("SRR101", config=campaign.data.config)  # not the cohort configuration
    with pytest.raises(ValueError, match="does not verify"):
        bind_produced_output(sample, context.expectation("SRR101"), index_digests=context.index_digests)
    assert not (sample / REFERENCE_BINDING_FILENAME).exists()
    orchestrator.db.close()
