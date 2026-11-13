"""Campaign scheduling and throughput: global order, streaming discovery, retry bookkeeping, host profile.

Real ProgressDB files, real FASTQ bytes and real orchestrator runs; only the
per-species discovery and the per-sample worker are replaced by small
deterministic functions (they are the scheduler's inputs, not the code under
test).
"""

from __future__ import annotations

import gzip
import json
import logging
import os
import subprocess
import sys
import threading
import time
from pathlib import Path
from typing import Any

import pytest

from metainformant.rna.engine import streaming_orchestrator as orchestrator_module
from metainformant.rna.engine.progress_db import (
    GATE_BACKOFF,
    GATE_HELD,
    GATE_READY,
    ProgressDB,
    RetryPolicy,
    classify_sample_error,
    retry_class_for,
    retry_gate,
    retry_status,
)
from metainformant.rna.engine.streaming_orchestrator import (
    FAST_FAIL_ENA_ENVIRONMENT,
    HostFacts,
    StreamingPipelineOrchestrator,
    _layout_contradiction,
    _local_sra_path,
    _prefetch_no_archive_reason,
    _raw_rejection_marker,
    _validated_local_fastq_inputs,
    apply_fast_fail_ena_environment,
    build_host_resource_profile,
    build_pipeline_resource_profile,
    order_campaign_tasks,
    resource_profile_environment,
)

REPO_ROOT = Path(__file__).resolve().parents[2]
GiB = 1024**3
MiB = 1024**2
DOWNLOAD_FAILED = "Download Failed (all sources: ENA FTP/HTTP, NCBI)"
KALLISTO_EXIT_1 = "ERROR: kallisto quant failed with exit code 1 for SRR_K."

PROFILE_VARIABLES = (
    "AMALGKIT_PIPELINE_ACQUISITION_SLOTS",
    "AMALGKIT_PIPELINE_ENA_FILE_WORKERS",
    "AMALGKIT_ENA_MAX_CONNECTIONS",
    "AMALGKIT_ENA_API_CONNECTIONS",
    "AMALGKIT_PIPELINE_QUANT_SLOTS",
    "AMALGKIT_PIPELINE_FASTQ_THREADS",
    "AMALGKIT_PIPELINE_FASTQ_SLOTS",
    "AMALGKIT_PIPELINE_COMPRESSION_THREADS",
    "AMALGKIT_PIPELINE_VALIDATION_SLOTS",
    "AMALGKIT_PIPELINE_MAX_IN_FLIGHT",
    "AMALGKIT_PIPELINE_SCHEDULING_ORDER",
    "AMALGKIT_PIPELINE_RAW_RESERVATION_GB",
    "AMALGKIT_RETRY_BACKOFF_BASE_SECONDS",
    "AMALGKIT_RETRY_BACKOFF_CAP_SECONDS",
    "AMALGKIT_MAX_SAMPLE_ATTEMPTS",
    "AMALGKIT_RETRY_WAIT_MAX_SECONDS",
    "PIPELINE_WORKERS",
    "PIPELINE_THREADS",
)


@pytest.fixture(autouse=True)
def _clean_profile_environment(monkeypatch: pytest.MonkeyPatch) -> None:
    for variable in PROFILE_VARIABLES:
        monkeypatch.delenv(variable, raising=False)


def _task(species: str, run: str, *, total_bases: int = 0, fastq_bytes: int = 0, batch: int = 1) -> dict[str, Any]:
    task: dict[str, Any] = {
        "srr": run,
        "batch_idx": batch,
        "species_name": species,
        "total_bases": total_bases,
        "_raw_input_priority": 2,
        "_raw_input_bytes": 0,
        "_uses_local_sra": False,
    }
    if fastq_bytes:
        task["fastq_bytes"] = fastq_bytes
    return task


def _runs(tasks: list[dict[str, Any]]) -> list[str]:
    return [f"{task['species_name']}/{task['srr']}" for task in tasks]


# ---------------------------------------------------------------------------
# Global ordering
# ---------------------------------------------------------------------------


def test_size_bands_order_two_species_smallest_band_first_with_round_robin() -> None:
    """Size bands run smallest first across species; inside a band the species alternate."""

    tasks = [
        _task("apis", "A_BIG", fastq_bytes=3 * GiB),
        _task("apis", "A_S2", fastq_bytes=60 * MiB),
        _task("apis", "A_S1", fastq_bytes=10 * MiB),
        _task("apis", "A_S3", fastq_bytes=90 * MiB),
        _task("bombus", "B_S1", fastq_bytes=20 * MiB),
        _task("bombus", "B_MID", fastq_bytes=300 * MiB),
        # No ENA size: total_bases x 0.5 (2e8 bases ~ 95 MiB, band 0).
        _task("bombus", "B_BASES", total_bases=200_000_000),
        _task("bombus", "B_UNKNOWN"),
    ]

    ordered = order_campaign_tasks(tasks, order="size_bands")

    assert _runs(ordered) == [
        # band <100 MiB: apis and bombus alternate, each smallest first
        "apis/A_S1",
        "bombus/B_S1",
        "apis/A_S2",
        "bombus/B_BASES",
        "apis/A_S3",
        # band 100-500 MiB, then 2-5 GiB, then unknown sizes last
        "bombus/B_MID",
        "apis/A_BIG",
        "bombus/B_UNKNOWN",
    ]


def test_smallest_first_and_species_order_are_selectable() -> None:
    tasks = [
        _task("apis", "A2", fastq_bytes=50),
        _task("apis", "A1", fastq_bytes=40),
        _task("bombus", "B1", fastq_bytes=10),
    ]

    assert _runs(order_campaign_tasks(tasks, order="smallest_first")) == ["bombus/B1", "apis/A1", "apis/A2"]
    assert _runs(order_campaign_tasks(tasks, order="species_order")) == ["apis/A2", "apis/A1", "bombus/B1"]
    with pytest.raises(ValueError, match="scheduling order"):
        order_campaign_tasks(tasks, order="random")


def test_order_is_independent_of_which_species_discovery_finished_first() -> None:
    """Ties break on (species, run), so the queue never depends on discovery completion order."""

    apis = [_task("apis", f"A{i}", fastq_bytes=100) for i in range(3)]
    bombus = [_task("bombus", f"B{i}", fastq_bytes=100) for i in range(3)]

    first = order_campaign_tasks([*apis, *bombus])
    second = order_campaign_tasks([*bombus, *apis])

    assert _runs(first) == _runs(second)
    assert _runs(first) == ["apis/A0", "bombus/B0", "apis/A1", "bombus/B1", "apis/A2", "bombus/B2"]


def test_local_inputs_and_partials_stay_ahead_of_fresh_acquisitions() -> None:
    fresh_small = _task("apis", "FRESH_SMALL", fastq_bytes=1)
    reusable = _task("bombus", "REUSE", fastq_bytes=5 * GiB) | {"_raw_input_priority": 0, "_raw_input_bytes": 9}
    partial = _task("apis", "PARTIAL", fastq_bytes=5 * GiB) | {"_raw_input_priority": 1, "_raw_input_bytes": 7}

    assert _runs(order_campaign_tasks([fresh_small, partial, reusable])) == [
        "bombus/REUSE",
        "apis/PARTIAL",
        "apis/FRESH_SMALL",
    ]


# ---------------------------------------------------------------------------
# Streaming discovery and run_all scheduling
# ---------------------------------------------------------------------------


def _stub_preflight(_data_root: object, **_kwargs: object) -> dict[str, str]:
    return {"data_root": "stub", "amalgkit_cli": "stub"}


def _orchestrator(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> StreamingPipelineOrchestrator:
    data_root = tmp_path / "data_root"
    data_root.mkdir(exist_ok=True)
    monkeypatch.setenv("AMALGKIT_DATA_ROOT", str(data_root))
    monkeypatch.setattr(orchestrator_module, "run_campaign_preflight", _stub_preflight)
    return StreamingPipelineOrchestrator(tmp_path / "config", tmp_path / "logs", db_path=data_root / "progress.db")


def _species_task(tmp_path: Path, species: str, run: str, **extra: Any) -> dict[str, Any]:
    return {
        "srr": run,
        "batch_idx": 1,
        "fastq_dir": tmp_path / "fastq" / species,
        "config_path": tmp_path / f"{species}.yaml",
        "species_name": species,
        "expected_paired": False,
        **extra,
    }


def _discoverer(
    orchestrator: StreamingPipelineOrchestrator,
    tasks_by_config: dict[str, list[dict[str, Any]]],
    *,
    exclude_done: bool = True,
) -> Any:
    def discover(config_name: str, *_args: object) -> list[dict[str, Any]]:
        # Like the cohort discovery: register the rows, then leave out quantified and quarantined ones.
        tasks = [dict(task) for task in tasks_by_config[config_name]]
        for species in {task["species_name"] for task in tasks}:
            orchestrator.db.init_species(species, [task["srr"] for task in tasks if task["species_name"] == species])
        if not exclude_done:
            return tasks
        return [
            task
            for task in tasks
            if orchestrator.db.get_state(task["species_name"], task["srr"]) not in ("quantified", "quarantined")
        ]

    return discover


def test_streaming_dispatch_starts_work_before_a_slow_discovery_finishes(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """The fast species' sample runs while the slow species is still being discovered."""

    orchestrator = _orchestrator(tmp_path, monkeypatch)
    fast_sample_started = threading.Event()
    slow_discovery_saw_execution: list[bool] = []
    tasks = {
        "amalgkit_fast.yaml": [_species_task(tmp_path, "fast", "SRR_FAST", total_bases=2_000_000_000)],
        "amalgkit_slow.yaml": [_species_task(tmp_path, "slow", "SRR_SLOW", total_bases=1_000)],
    }
    discover = _discoverer(orchestrator, tasks)

    def fake_discover(config_name: str, *args: object) -> list[dict[str, Any]]:
        if config_name == "amalgkit_slow.yaml":
            # A barrier-style scheduler would never start SRR_FAST before this returns.
            slow_discovery_saw_execution.append(fast_sample_started.wait(timeout=10))
        return discover(config_name, *args)

    executed: list[str] = []

    def fake_process(srr: str, _batch: int, fastq_dir: Path, _config: Path, species: str, *_a: Any, **_k: Any) -> Any:
        executed.append(srr)
        if srr == "SRR_FAST":
            fast_sample_started.set()
        orchestrator.db.set_state(species, srr, "quantified")
        return {"quantified": True, "skipped": False, "error": None, "input_bytes": 1_000_000}

    monkeypatch.setattr(orchestrator, "discover_species_tasks", fake_discover)
    monkeypatch.setattr(orchestrator, "process_single_sample", fake_process)

    summary = orchestrator.run_all(
        ["amalgkit_fast.yaml", "amalgkit_slow.yaml"], max_gb=1.0, workers=2, threads=2, discovery_workers=2
    )

    assert slow_discovery_saw_execution == [True]
    assert executed == ["SRR_FAST", "SRR_SLOW"]
    assert summary["exit_code"] == 0
    assert summary["throughput"]["quantified_new"] == 2
    assert summary["throughput"]["fastq_gb_per_hour"] > 0
    assert summary["resource_profile"]["scheduling_order"] == "size_bands"


def test_raw_reservation_admits_only_what_fits_in_flight(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """With a 150-byte raw reservation, 100-byte acquisitions run one at a time despite four workers."""

    orchestrator = _orchestrator(tmp_path, monkeypatch)
    tasks = {
        "amalgkit_apis.yaml": [_species_task(tmp_path, "apis", f"SRR{i}", fastq_bytes=100) for i in range(4)],
    }
    monkeypatch.setattr(orchestrator, "discover_species_tasks", _discoverer(orchestrator, tasks))
    active = 0
    peak = 0
    lock = threading.Lock()

    def fake_process(srr: str, _batch: int, _fastq: Path, _config: Path, species: str, *_a: Any, **_k: Any) -> Any:
        nonlocal active, peak
        with lock:
            active += 1
            peak = max(peak, active)
        time.sleep(0.05)
        with lock:
            active -= 1
        orchestrator.db.set_state(species, srr, "quantified")
        return {"quantified": True, "skipped": False, "error": None}

    monkeypatch.setattr(orchestrator, "process_single_sample", fake_process)
    profile = build_pipeline_resource_profile(4, 4, cpu_count=8)
    from dataclasses import replace

    summary = orchestrator.run_all(
        ["amalgkit_apis.yaml"], max_gb=1.0, resource_profile=replace(profile, raw_reservation_bytes=150)
    )

    assert summary["exit_code"] == 0
    assert summary["quantified_this_run"] == 4
    assert peak == 1


def test_transient_failure_backs_off_and_is_retried_in_run_when_due(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """A transient failure waits out its backoff, then the same run retries it."""

    monkeypatch.setenv("AMALGKIT_RETRY_BACKOFF_BASE_SECONDS", "1")
    monkeypatch.setenv("AMALGKIT_RETRY_BACKOFF_CAP_SECONDS", "1")
    monkeypatch.setenv("AMALGKIT_RETRY_WAIT_MAX_SECONDS", "5")
    orchestrator = _orchestrator(tmp_path, monkeypatch)
    tasks = {"amalgkit_apis.yaml": [_species_task(tmp_path, "apis", "SRR_FLAKY")]}
    monkeypatch.setattr(orchestrator, "discover_species_tasks", _discoverer(orchestrator, tasks))
    attempts: list[float] = []

    def fake_process(srr: str, _batch: int, _fastq: Path, _config: Path, species: str, *_a: Any, **_k: Any) -> Any:
        attempts.append(time.time())
        if len(attempts) == 1:
            orchestrator.db.set_state(species, srr, "failed", error=DOWNLOAD_FAILED)
            return {"quantified": False, "skipped": False, "error": DOWNLOAD_FAILED}
        orchestrator.db.set_state(species, srr, "quantified")
        return {"quantified": True, "skipped": False, "error": None}

    monkeypatch.setattr(orchestrator, "process_single_sample", fake_process)

    summary = orchestrator.run_all(["amalgkit_apis.yaml"], max_gb=1.0, workers=1, threads=1)

    assert len(attempts) == 2
    assert attempts[1] - attempts[0] >= 0.5
    assert summary["exit_code"] == 0
    assert summary["deferred_backoff"] == 0
    # Success clears the bookkeeping.
    assert orchestrator.db.retry_records("apis") == {}


def test_backoff_defers_the_task_at_the_next_start(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """Without an in-process wait the run exits 75 with the due time; the next start skips the task until then."""

    orchestrator = _orchestrator(tmp_path, monkeypatch)
    tasks = {
        "amalgkit_apis.yaml": [
            _species_task(tmp_path, "apis", "SRR_FLAKY"),
            _species_task(tmp_path, "apis", "SRR_OK"),
        ]
    }
    monkeypatch.setattr(orchestrator, "discover_species_tasks", _discoverer(orchestrator, tasks))
    calls: list[str] = []

    def fake_process(srr: str, _batch: int, _fastq: Path, _config: Path, species: str, *_a: Any, **_k: Any) -> Any:
        calls.append(srr)
        if srr == "SRR_FLAKY":
            orchestrator.db.set_state(species, srr, "failed", error=DOWNLOAD_FAILED)
            return {"quantified": False, "skipped": False, "error": DOWNLOAD_FAILED}
        orchestrator.db.set_state(species, srr, "quantified")
        return {"quantified": True, "skipped": False, "error": None}

    monkeypatch.setattr(orchestrator, "process_single_sample", fake_process)
    started = time.time()

    first = orchestrator.run_all(["amalgkit_apis.yaml"], max_gb=1.0, workers=1, threads=1)

    assert sorted(calls) == ["SRR_FLAKY", "SRR_OK"]
    assert first["exit_code"] == 75
    assert first["deferred_backoff"] == 1
    assert first["remaining_by_state"] == {"failed": 1}
    due = time.mktime(time.strptime(first["next_attempt_due_at"], "%Y-%m-%dT%H:%M:%SZ")) - time.timezone
    assert started + 3600 - 2 <= due <= time.time() + 3600 + 1
    record = orchestrator.db.retry_records("apis")[("apis", "SRR_FLAKY")]
    assert (record.attempts, record.last_error_class) == (1, "transfer_all_sources_failed")

    calls.clear()
    second = orchestrator.run_all(["amalgkit_apis.yaml"], max_gb=1.0, workers=1, threads=1)

    assert calls == []
    assert second["exit_code"] == 75
    assert second["deferred_backoff"] == 1
    assert orchestrator.db.retry_records("apis")[("apis", "SRR_FLAKY")].attempts == 1


def test_layout_mismatch_is_quarantined_and_not_retried_at_every_start(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    orchestrator = _orchestrator(tmp_path, monkeypatch)
    tasks = {"amalgkit_apis.yaml": [_species_task(tmp_path, "apis", "SRR_LAYOUT")]}
    monkeypatch.setattr(orchestrator, "discover_species_tasks", _discoverer(orchestrator, tasks))
    calls: list[str] = []
    error = f"{DOWNLOAD_FAILED}; layout_mismatch: ENA files: declared single, source delivered mate files X_1, X_2"

    def fake_process(srr: str, _batch: int, _fastq: Path, _config: Path, species: str, *_a: Any, **_k: Any) -> Any:
        calls.append(srr)
        orchestrator.db.set_state(species, srr, "failed", error=error)
        return {"quantified": False, "skipped": False, "error": error}

    monkeypatch.setattr(orchestrator, "process_single_sample", fake_process)

    first = orchestrator.run_all(["amalgkit_apis.yaml"], max_gb=1.0, workers=1, threads=1)

    assert calls == ["SRR_LAYOUT"]
    assert orchestrator.db.get_state("apis", "SRR_LAYOUT") == "quarantined"
    assert first["held_by_reason"] == {"layout_mismatch": 1}
    assert first["exit_code"] == 0

    second = orchestrator.run_all(["amalgkit_apis.yaml"], max_gb=1.0, workers=1, threads=1)
    assert calls == ["SRR_LAYOUT"]
    assert second["exit_code"] == 0
    # A discovery that still lists the quarantined row (the legacy metadata path) meets the retry gate instead.
    monkeypatch.setattr(orchestrator, "discover_species_tasks", _discoverer(orchestrator, tasks, exclude_done=False))
    third = orchestrator.run_all(["amalgkit_apis.yaml"], max_gb=1.0, workers=1, threads=1)
    assert calls == ["SRR_LAYOUT"]
    assert third["held_by_reason"] == {"layout_mismatch": 1}
    assert third["exit_code"] == 0

    # The operator reconciles the inputs and releases the hold: the next start schedules the run again.
    released = orchestrator.db.release_retry_holds(policy=RetryPolicy())
    assert released == {"layout_mismatch": 1}
    assert orchestrator.db.get_state("apis", "SRR_LAYOUT") == "pending"
    orchestrator.run_all(["amalgkit_apis.yaml"], max_gb=1.0, workers=1, threads=1)
    assert calls == ["SRR_LAYOUT", "SRR_LAYOUT"]


def test_kallisto_exit_1_is_permanent_and_never_retried(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    orchestrator = _orchestrator(tmp_path, monkeypatch)
    tasks = {"amalgkit_apis.yaml": [_species_task(tmp_path, "apis", "SRR_K")]}
    monkeypatch.setattr(orchestrator, "discover_species_tasks", _discoverer(orchestrator, tasks))
    calls: list[str] = []

    def fake_process(srr: str, _batch: int, _fastq: Path, _config: Path, species: str, *_a: Any, **_k: Any) -> Any:
        calls.append(srr)
        orchestrator.db.set_state(species, srr, "failed", error=KALLISTO_EXIT_1)
        return {"quantified": False, "skipped": False, "error": KALLISTO_EXIT_1}

    monkeypatch.setattr(orchestrator, "process_single_sample", fake_process)

    first = orchestrator.run_all(["amalgkit_apis.yaml"], max_gb=1.0, workers=1, threads=1)
    second = orchestrator.run_all(["amalgkit_apis.yaml"], max_gb=1.0, workers=1, threads=1)

    assert calls == ["SRR_K"]
    assert orchestrator.db.get_state("apis", "SRR_K") == "failed"
    assert first["held_by_reason"] == second["held_by_reason"] == {"permanent": 1}
    assert first["exit_code"] == second["exit_code"] == 0


# ---------------------------------------------------------------------------
# Retry bookkeeping in the progress database
# ---------------------------------------------------------------------------


def test_backoff_policy_is_capped_exponential() -> None:
    policy = RetryPolicy(base_seconds=3600, cap_seconds=86400, max_attempts=6)

    assert [policy.backoff_seconds(n) for n in range(0, 7)] == [0, 3600, 7200, 14400, 28800, 57600, 86400]
    assert policy.backoff_seconds(60) == 86400
    with pytest.raises(ValueError):
        RetryPolicy(base_seconds=10, cap_seconds=5)


def test_error_classes_map_to_retry_classes() -> None:
    assert retry_class_for(classify_sample_error(KALLISTO_EXIT_1)) == "permanent"
    assert retry_class_for(classify_sample_error(f"{DOWNLOAD_FAILED}; layout_mismatch: x")) == "layout_mismatch"
    assert retry_class_for(classify_sample_error(DOWNLOAD_FAILED)) == "transient"
    assert retry_class_for(classify_sample_error("[Errno 1] Operation not permitted: '/Volumes/x'")) == "environment"
    # Messages that only mention kallisto keep their previous classes.
    assert classify_sample_error("Quantification Failed: kallisto exited 1") == "quantification_failed"
    assert classify_sample_error("kallisto exited with exit code 1") == "unclassified"


def test_record_attempt_outcome_schedules_backoff_then_holds_after_the_budget(tmp_path: Path) -> None:
    policy = RetryPolicy(base_seconds=60, cap_seconds=120, max_attempts=3)
    now = 1_800_000_000.0
    with ProgressDB(tmp_path / "progress.db") as db:
        db.init_species("apis", ["SRR1"])
        assert db.retry_records("apis") == {}  # migrates the empty database first
        decisions = []
        for attempt in range(3):
            db.set_state("apis", "SRR1", "failed", error=DOWNLOAD_FAILED)
            decisions.append(db.record_attempt_outcome("apis", "SRR1", DOWNLOAD_FAILED, policy=policy, now=now))
        record = db.retry_records("apis")[("apis", "SRR1")]

    assert [(d["decision"], d["attempts"]) for d in decisions] == [("backoff", 1), ("backoff", 2), ("held", 3)]
    assert decisions[0]["next_attempt_at"] == now + 60
    assert decisions[1]["next_attempt_at"] == now + 120
    assert retry_gate(record, policy, now) == (GATE_HELD, "attempts_exhausted:transfer_all_sources_failed")
    two = record.__class__(**{**record.__dict__, "attempts": 2, "next_attempt_at": now + 120})
    assert retry_gate(two, policy, now) == (GATE_BACKOFF, "backoff:transfer_all_sources_failed")
    assert retry_gate(two, policy, now + 121) == (GATE_READY, None)


def test_environment_failures_do_not_spend_attempts(tmp_path: Path) -> None:
    with ProgressDB(tmp_path / "progress.db") as db:
        db.init_species("apis", ["SRR1"])
        outcome = db.record_attempt_outcome(
            "apis", "SRR1", "[Errno 30] Read-only file system: '/Volumes/x'", policy=RetryPolicy()
        )
    assert outcome["decision"] == "retry_next_start"
    assert outcome["attempts"] == 0


def test_migration_backfills_failed_rows_and_status_reads_them(tmp_path: Path) -> None:
    """The additive columns classify existing failures: kallisto exit 1 is held, transfers are due at once."""

    db_path = tmp_path / "progress.db"
    with ProgressDB(db_path) as db:
        db.init_species("apis", ["SRR_K", "SRR_T", "SRR_Q"])
        db.set_state("apis", "SRR_K", "failed", error=KALLISTO_EXIT_1)
        db.set_state("apis", "SRR_T", "failed", error=DOWNLOAD_FAILED)
        db.set_state("apis", "SRR_Q", "quantified")
    assert retry_status(db_path)["held_by_reason"] == {}  # read-only: no migration yet

    with ProgressDB(db_path) as db:
        records = db.retry_records("apis")
    assert sorted(records) == [("apis", "SRR_K"), ("apis", "SRR_T")]
    assert records[("apis", "SRR_K")].attempts == 1
    assert retry_gate(records[("apis", "SRR_K")], RetryPolicy(), time.time())[0] == GATE_HELD
    assert retry_gate(records[("apis", "SRR_T")], RetryPolicy(), time.time()) == (GATE_READY, None)
    status = retry_status(db_path)
    assert status["held_by_reason"] == {"permanent": 1}
    assert status["backoff"] == 0


def test_release_holds_cli_returns_held_rows_to_pending(tmp_path: Path) -> None:
    data_root = tmp_path / "root"
    data_root.mkdir()
    db_path = data_root / "pipeline_progress.db"
    with ProgressDB(db_path) as db:
        db.init_species("apis", ["SRR_K"])
        db.set_state("apis", "SRR_K", "failed", error=KALLISTO_EXIT_1)
        db.record_attempt_outcome("apis", "SRR_K", KALLISTO_EXIT_1, policy=RetryPolicy())

    result = subprocess.run(
        [sys.executable, "-m", "metainformant.rna.engine.progress_db", "release-holds", "--db", str(db_path)],
        capture_output=True,
        text=True,
        cwd=REPO_ROOT,
        env={**os.environ, "PYTHONPATH": str(REPO_ROOT / "src")},
        timeout=120,
    )

    assert result.returncode == 0, result.stderr
    assert json.loads(result.stdout) == {"released": {"permanent": 1}}
    with ProgressDB(db_path) as db:
        assert db.get_state("apis", "SRR_K") == "pending"
        assert db.retry_records("apis") == {}


# ---------------------------------------------------------------------------
# Negative raw-validation witness
# ---------------------------------------------------------------------------


def test_rejected_one_file_paired_input_is_not_rescanned_until_it_changes(
    tmp_path: Path, caplog: pytest.LogCaptureFixture
) -> None:
    sample_dir = tmp_path / "apis" / "work" / "getfastq" / "SRR9"
    sample_dir.mkdir(parents=True)
    source = sample_dir / "SRR9.fastq.gz"
    # Declared paired, one file, and the second record has an invalid separator: deterministic content failure.
    with gzip.open(source, "wt") as handle:
        handle.write("@r1/1\nACGT\n+\nIIII\n@r1/2\nACGT\nX\nIIII\n")
    caplog.set_level(logging.INFO, logger=orchestrator_module.logger.name)

    assert _validated_local_fastq_inputs(sample_dir, "SRR9", True, require_records=False) == []
    marker = _raw_rejection_marker(sample_dir, "SRR9")
    payload = json.loads(marker.read_text())
    assert payload["files"][0]["name"] == "SRR9.fastq.gz"
    assert set(payload["files"][0]) == {"name", "inode", "size", "mtime_ns"}
    assert "Rejecting one-file paired input" in caplog.text

    caplog.clear()
    assert _validated_local_fastq_inputs(sample_dir, "SRR9", True, require_records=False) == []
    assert "Skipping the rescan of SRR9" in caplog.text
    assert "Rejecting one-file paired input" not in caplog.text

    # Rewriting the file (a new download) invalidates the witness: the bytes are scanned again.
    os.utime(source, ns=(time.time_ns(), time.time_ns() + 5_000_000_000))
    caplog.clear()
    assert _validated_local_fastq_inputs(sample_dir, "SRR9", True, require_records=False) == []
    assert "Rejecting one-file paired input" in caplog.text


def test_layout_contradiction_names_the_mismatch() -> None:
    mates = [Path("SRR1_1.fastq.gz"), Path("SRR1_2.fastq.gz")]
    assert "declared single" in str(_layout_contradiction(mates, "SRR1", False))
    assert "one mate" in str(_layout_contradiction(mates[:1], "SRR1", True))
    assert _layout_contradiction(mates, "SRR1", True) is None
    assert _layout_contradiction([Path("SRR1.fastq.gz")], "SRR1", False) is None


# ---------------------------------------------------------------------------
# prefetch .sralite
# ---------------------------------------------------------------------------


def test_sralite_archive_is_a_local_sra_source(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setenv("AMALGKIT_DATA_ROOT", str(tmp_path / "root"))
    sample_dir = tmp_path / "getfastq" / "SRR5"
    archive = sample_dir / ".prefetch" / "SRR5" / "SRR5" / "SRR5.sralite"
    archive.parent.mkdir(parents=True)
    archive.write_bytes(b"lite")

    assert _local_sra_path(sample_dir, sample_dir.parent, "SRR5") == archive


def test_prefetch_sralite_output_is_accepted(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """sra-tools 3.4.1 names SRA Lite downloads <acc>.sralite; a zero-exit prefetch that wrote one succeeded."""

    data_root = tmp_path / "root"
    data_root.mkdir()
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    prefetch = bin_dir / "prefetch"
    # A stand-in for the external tool: writes <output-directory>/<acc>/<acc>.sralite like prefetch does.
    prefetch.write_text(
        "#!/bin/sh\n"
        'while [ "$#" -gt 1 ]; do\n'
        '  if [ "$1" = "--output-directory" ]; then out="$2"; fi\n'
        "  shift\n"
        "done\n"
        'mkdir -p "$out/$1" && printf lite > "$out/$1/$1.sralite"\n'
    )
    prefetch.chmod(0o755)
    monkeypatch.setenv("PATH", f"{bin_dir}{os.pathsep}{os.environ['PATH']}")
    orchestrator = _orchestrator(tmp_path, monkeypatch)
    orchestrator.min_external_free_gb = 0.0
    sample_dir = tmp_path / "getfastq" / "SRR6"
    sample_dir.mkdir(parents=True)

    archive = orchestrator._prefetch_sra_accession("SRR6", sample_dir)

    assert archive == sample_dir / ".prefetch" / "SRR6" / "SRR6" / "SRR6.sralite"


def test_prefetch_without_archive_names_the_reason(tmp_path: Path) -> None:
    cached = tmp_path / "SRR7.sra"
    assert str(cached) in _prefetch_no_archive_reason(
        "2026-11-12T00:00:00 prefetch.3.4.1: 'SRR7' is found locally", cached
    )
    assert "max-size" in _prefetch_no_archive_reason("'SRR7' (9 GB) is larger than maximum allowed: skipped", None)
    assert _prefetch_no_archive_reason("", None) == ""


# ---------------------------------------------------------------------------
# Host-derived resource profile
# ---------------------------------------------------------------------------

M3_HOST = HostFacts(cpu_count=8, memory_bytes=16 * GiB, performance_cores=4, efficiency_cores=4)
NYLANDERI_INDEX_BYTES = 80_800_000


def test_host_profile_for_an_8_core_16_gib_host_matches_the_measured_c19_settings() -> None:
    profile = build_host_resource_profile(
        host=M3_HOST, kallisto_index_bytes=NYLANDERI_INDEX_BYTES, data_root_free_bytes=4 * 1024 * GiB
    )

    assert profile.source == "host"
    assert (profile.acquisition_slots, profile.ena_file_workers) == (8, 2)
    assert (profile.ena_max_connections, profile.ena_api_connections) == (16, 2)
    assert (profile.quant_slots, profile.quant_threads_per_worker) == (8, 1)
    assert profile.effective_quant_threads <= M3_HOST.cpu_count
    assert profile.validation_slots == 2
    assert (profile.fasterq_slots, profile.fasterq_threads, profile.compression_threads) == (1, 4, 2)
    assert (profile.workers, profile.max_in_flight) == (16, 16)
    assert profile.raw_reservation_bytes == 256 * GiB
    assert profile.scheduling_order == "size_bands"
    # 50% of 16 GiB over (256 MiB + 4 x 80.8 MB) per kallisto run: 14 runs would fit, so 8 stands.
    assert profile.memory_quant_slot_cap == 14
    assert profile.as_dict()["host_performance_cores"] == 4


def test_host_profile_memory_guard_and_free_space_bound_a_small_host() -> None:
    small = HostFacts(cpu_count=4, memory_bytes=2 * GiB)

    profile = build_host_resource_profile(
        host=small, kallisto_index_bytes=NYLANDERI_INDEX_BYTES, data_root_free_bytes=100 * GiB, quant_slots=4
    )

    assert profile.memory_quant_slot_cap == 1
    assert profile.quant_slots == 1
    assert profile.raw_reservation_bytes == 25 * GiB
    assert profile.fasterq_threads == 2


def test_host_profile_precedence_is_argument_then_environment_then_default(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setenv("AMALGKIT_PIPELINE_ACQUISITION_SLOTS", "4")
    monkeypatch.setenv("AMALGKIT_PIPELINE_SCHEDULING_ORDER", "smallest_first")
    monkeypatch.setenv("AMALGKIT_PIPELINE_RAW_RESERVATION_GB", "10")

    from_environment = build_host_resource_profile(host=M3_HOST)
    explicit = build_host_resource_profile(host=M3_HOST, acquisition_slots=3, scheduling_order="species_order")

    assert (from_environment.acquisition_slots, from_environment.ena_max_connections) == (4, 8)
    assert from_environment.workers == 4 + 8
    assert from_environment.scheduling_order == "smallest_first"
    assert from_environment.raw_reservation_bytes == 10 * GiB
    assert (explicit.acquisition_slots, explicit.scheduling_order) == (3, "species_order")
    with pytest.raises(ValueError, match="scheduling order"):
        build_host_resource_profile(host=M3_HOST, scheduling_order="random")


def test_profile_environment_is_scoped_to_the_run(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setenv("AMALGKIT_ENA_MAX_CONNECTIONS", "3")
    profile = build_host_resource_profile(host=M3_HOST, ena_max_connections=16)

    with resource_profile_environment(profile):
        assert os.environ["AMALGKIT_ENA_MAX_CONNECTIONS"] == "16"
        assert os.environ["AMALGKIT_PIPELINE_ACQUISITION_SLOTS"] == "8"
    assert os.environ["AMALGKIT_ENA_MAX_CONNECTIONS"] == "3"
    assert "AMALGKIT_PIPELINE_ACQUISITION_SLOTS" not in os.environ


def test_fast_fail_ena_is_opt_in_and_keeps_explicit_settings() -> None:
    environment = {"AMALGKIT_PIPELINE_ENA_RETRIES": "2"}

    applied = apply_fast_fail_ena_environment(environment)

    assert environment["AMALGKIT_PIPELINE_ENA_RETRIES"] == "2"
    assert applied == {
        key: value for key, value in FAST_FAIL_ENA_ENVIRONMENT.items() if key != "AMALGKIT_PIPELINE_ENA_RETRIES"
    }
    assert environment["AMALGKIT_PIPELINE_DOWNLOAD_SPEED_LIMIT_BYTES"] == "65536"


def test_print_profile_resolves_without_the_lock_or_database(tmp_path: Path) -> None:
    data_root = tmp_path / "root"
    data_root.mkdir()
    environment = {key: value for key, value in os.environ.items() if key not in PROFILE_VARIABLES}
    environment["AMALGKIT_DATA_ROOT"] = str(data_root)

    result = subprocess.run(
        [
            sys.executable,
            str(REPO_ROOT / "scripts" / "rna" / "run_all_species.py"),
            "--print-profile",
            "--acquisition-slots",
            "3",
            "--ena-file-workers",
            "1",
            "--scheduling-order",
            "smallest_first",
        ],
        capture_output=True,
        text=True,
        env=environment,
        timeout=120,
    )

    assert result.returncode == 0, result.stderr
    profile = json.loads(result.stdout)
    assert profile["source"] == "host"
    assert (profile["acquisition_slots"], profile["ena_file_workers"], profile["ena_max_connections"]) == (3, 1, 3)
    assert profile["scheduling_order"] == "smallest_first"
    assert profile["effective_quant_threads"] <= profile["host_cpu_count"]
    assert sorted(path.name for path in data_root.iterdir()) == []


def test_process_single_sample_reports_layout_mismatch_and_input_bytes(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """A noted layout contradiction reaches the stored failure; a quantified sample reports its FASTQ bytes."""

    orchestrator = _orchestrator(tmp_path, monkeypatch)
    orchestrator.db.init_species("apis", ["SRR_BAD", "SRR_GOOD"])
    fastq_dir = tmp_path / "apis" / "work" / "getfastq"

    def acquire(srr_id: str, out_dir: Path, *_args: Any, **_kwargs: Any) -> bool:
        if srr_id == "SRR_BAD":
            orchestrator._note_layout_mismatch(srr_id, "ENA files: declared single, source delivered mate files")
            return False
        (out_dir / srr_id).mkdir(parents=True, exist_ok=True)
        (out_dir / srr_id / f"{srr_id}.fastq.gz").write_bytes(b"x" * 1234)
        return True

    monkeypatch.setattr(orchestrator, "download_fastq", acquire)
    monkeypatch.setattr(orchestrator, "quant_sample", lambda *_args, **_kwargs: (True, None))
    monkeypatch.setattr(orchestrator, "_reclaim_after_current_quantification", lambda *_args: None)

    bad = orchestrator.process_single_sample("SRR_BAD", 1, fastq_dir, tmp_path / "c.yaml", "apis", 1, False)
    good = orchestrator.process_single_sample("SRR_GOOD", 2, fastq_dir, tmp_path / "c.yaml", "apis", 1, False)

    assert bad["error"].startswith(DOWNLOAD_FAILED + "; layout_mismatch: ENA files")
    assert classify_sample_error(bad["error"]) == "layout_mismatch"
    assert orchestrator.db.get_failed("apis")[0]["error"] == bad["error"]
    assert good["quantified"] is True
    assert good["input_bytes"] == 1234
