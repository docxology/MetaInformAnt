"""Bounded, idempotent acquisition and quantification on local or AWS hosts."""
from __future__ import annotations

import concurrent.futures
import fcntl
import json
import os
import threading
import time
from collections import Counter, deque
from dataclasses import asdict
from pathlib import Path
from typing import Any
from uuid import uuid4
from metainformant.rna.engine.acquisition_manifest import sha256_file,load_snapshot,load_task_selection,verify_input_files
from metainformant.rna.engine.acquisition_references import prepare_reference_inputs
from metainformant.rna.engine.acquisition_manifest import verify_worker_configs
from metainformant.rna.engine.acquisition_allocation import verify_local_allocation
from metainformant.rna.engine.acquisition_sample import execute_manifest_task
from metainformant.rna.engine.streaming_orchestrator import StreamingPipelineOrchestrator,build_pipeline_resource_profile

def _run_manifest_owned(
    *,
    manifest_path: Path,
    data_root: Path,
    config_dir: Path,
    workers: int,
    threads: int,
    fastq_threads: int,
    compression_threads: int,
    validation_slots: int,
    quant_slots: int | None = None,
    fasterq_slots: int | None = None,
    max_in_flight: int | None = None,
    task_selection: Path | None = None,
) -> dict[str, Any]:
    """Run all manifest tasks with bounded idempotent concurrency."""

    snapshot, manifest_tasks = load_snapshot(manifest_path)
    selection, tasks = load_task_selection(
        task_selection,
        snapshot,
        manifest_tasks,
        snapshot_sha256=sha256_file(manifest_path.with_name("snapshot.json")),
    )
    verify_local_allocation(task_selection, manifest_path,
                            frozenset(task["task_id"] for task in manifest_tasks),
                            frozenset(task["task_id"] for task in tasks), snapshot.get("inventory_sha256"))
    profile = build_pipeline_resource_profile(
        workers,
        threads,
        quant_slots=quant_slots,
        fasterq_slots=fasterq_slots,
        fasterq_threads=fastq_threads,
        compression_threads=compression_threads,
        validation_slots=validation_slots,
        max_in_flight=max_in_flight,
    )
    verify_input_files(snapshot, manifest_path.parent)
    verify_worker_configs(snapshot, tasks, config_dir)
    os.environ["AMALGKIT_DATA_ROOT"] = str(data_root.resolve())
    os.environ.setdefault("AMALGKIT_RECLAIM_RAW_AFTER_QUANT", "no")
    os.environ.setdefault("AMALGKIT_MIN_EXTERNAL_FREE_GB", "8")
    os.environ.setdefault("AMALGKIT_MIN_SYSTEM_FREE_GB", "4")
    os.environ["AMALGKIT_PIPELINE_FASTQ_THREADS"] = str(fastq_threads)
    os.environ["AMALGKIT_PIPELINE_COMPRESSION_THREADS"] = str(compression_threads)
    os.environ["AMALGKIT_PIPELINE_VALIDATION_SLOTS"] = str(validation_slots)

    data_root.mkdir(parents=True, exist_ok=True)
    reference_preflight = prepare_reference_inputs(
        data_root=data_root,
        config_dir=config_dir,
        tasks=tasks,
    )
    durable_store = None
    durable_cohort = os.environ.get("AMALGKIT_DURABLE_COHORT", "")
    if os.environ.get("AMALGKIT_DURABLE_BUCKET"):
        if not durable_cohort:
            raise ValueError(
                "AMALGKIT_DURABLE_COHORT is required with a durable bucket"
            )
        from metainformant.rna.engine.durable_quant import S3Store

        durable_store = S3Store(
            os.environ["AMALGKIT_DURABLE_BUCKET"],
            os.environ.get("AMALGKIT_DURABLE_PREFIX", "locked-quant-v1"),
            region=os.environ.get("AWS_DEFAULT_REGION"),
        )
    orchestrator = StreamingPipelineOrchestrator(
        config_dir=config_dir,
        log_dir=data_root / "logs",
        db_path=data_root / "pipeline_progress.db",
    )
    try:
        orchestrator._resource_profile = profile
        orchestrator._quant_semaphore = threading.BoundedSemaphore(profile.quant_slots)
        orchestrator._fasterq_semaphore = threading.BoundedSemaphore(profile.fasterq_slots)
        orchestrator._raw_validation_semaphore = threading.BoundedSemaphore(
            profile.validation_slots
        )
        by_species: dict[str, list[dict[str, Any]]] = {}
        for task in tasks:
            by_species.setdefault(str(task["species"]), []).append(task)
        for species, species_tasks in by_species.items():
            orchestrator.db.init_species(
                species, [str(task["accession"]) for task in species_tasks]
            )

        started = time.time()
        monotonic_started = time.monotonic()
        counts: Counter[str] = Counter()
        results: list[dict[str, Any]] = []
        run_directory = data_root / "acquisition_runs" / uuid4().hex
        run_directory.mkdir(parents=True, exist_ok=False)
        journal_path = run_directory / "task_results.jsonl"
        journal_path.touch(exist_ok=False)

        def append_journal(result: dict[str, Any]) -> None:
            """Persist one completed task before scheduling its replacement."""

            with journal_path.open("a", encoding="utf-8") as handle:
                handle.write(json.dumps(result, sort_keys=True) + "\n")
                handle.flush()
                os.fsync(handle.fileno())

        def execute(task: dict[str, Any]) -> dict[str, Any]:
            return execute_manifest_task(task, orchestrator, data_root, config_dir, durable_store, durable_cohort, profile.quant_threads_per_worker)

        with concurrent.futures.ThreadPoolExecutor(max_workers=profile.workers) as executor:
            pending_tasks = deque(tasks)
            future_to_task: dict[concurrent.futures.Future[Any], dict[str, Any]] = {}
            raw_bytes_limit = int(os.environ.get("AMALGKIT_CLOUD_MAX_RAW_BYTES", "0"))
            reserved_bytes = 0
            if raw_bytes_limit:
                for task in tasks:
                    if (
                        int(task.get("fastq_bytes", 0)) <= 0
                        or int(task["fastq_bytes"]) > raw_bytes_limit
                    ):
                        raise ValueError("cloud task lacks a bounded raw-byte reservation")

            def submit_next() -> bool:
                nonlocal reserved_bytes
                if not pending_tasks:
                    return False
                task = pending_tasks[0]
                needed = int(task.get("fastq_bytes", 0))
                if raw_bytes_limit and reserved_bytes + needed > raw_bytes_limit:
                    return False
                pending_tasks.popleft()
                reserved_bytes += needed
                future_to_task[executor.submit(execute, task)] = task
                return True

            for _ in range(min(profile.max_in_flight, len(tasks))):
                submit_next()
            while future_to_task:
                done, _ = concurrent.futures.wait(
                    tuple(future_to_task),
                    return_when=concurrent.futures.FIRST_COMPLETED,
                )
                for future in done:
                    task = future_to_task.pop(future)
                    reserved_bytes -= int(task.get("fastq_bytes", 0))
                    try:
                        result = future.result()
                    except Exception as exc:  # noqa: BLE001 - record failure evidence at the external task boundary
                        result = {
                            "srr": task["accession"],
                            "species": task["species"],
                            "quantified": False,
                            "error": f"worker exception: {exc}",
                        }
                    result["species"] = task["species"]
                    result["task_id"] = task.get("task_id")
                    results.append(result)
                    append_journal(result)
                    counts["quantified" if result.get("quantified") else "failed"] += 1
                    if result.get("quantified"):
                        counts["reused" if result.get("skipped") else "newly_quantified"] += 1
                while len(future_to_task) < profile.max_in_flight and submit_next():
                    pass

        results.sort(key=lambda result: str(result.get("task_id", "")))

        summary = {
            "schema": {"metainformant.hymenoptera.gcp_snapshot.v1": "metainformant.hymenoptera.gcp_worker_result.v1",
                       "metainformant.rna.acquisition_snapshot.v1": "metainformant.rna.acquisition_result.v1"}[snapshot["schema"]],
            "started_at_epoch": started,
            "finished_at_epoch": time.time(),
            "elapsed_seconds": time.monotonic() - monotonic_started,
            "snapshot_sha256": sha256_file(manifest_path.with_name("snapshot.json")),
            "manifest_sha256": snapshot["manifest_sha256"],
            "snapshot": snapshot,
            "task_selection": selection,
            "reference_preflight": reference_preflight,
            "data_root": str(data_root.resolve()),
            "workers": profile.workers,
            "quant_slots": profile.quant_slots,
            "quant_threads_per_sample": profile.quant_threads_per_worker,
            "resource_profile": {
                **asdict(profile),
                "effective_quant_threads": profile.effective_quant_threads,
                "peak_stage_threads": profile.peak_stage_threads,
            },
            "task_results_journal": str(journal_path),
            "task_results_journal_sha256": sha256_file(journal_path),
            "task_count": len(tasks),
            "manifest_task_count": len(manifest_tasks),
            "counts": dict(counts),
            "results": results,
        }
        result_path = data_root / "cloud_worker_result.json"
        (run_directory / "result.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        result_path.write_text(
            json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )
        return summary

    finally:
        orchestrator.db.close()


_PROCESS_WORKER_LOCK = threading.Lock()

def run_manifest(*, manifest_path: Path, data_root: Path, config_dir: Path, workers: int, threads: int, fastq_threads: int, compression_threads: int, validation_slots: int, quant_slots: int | None=None, fasterq_slots: int | None=None, max_in_flight: int | None=None, task_selection: Path | None=None) -> dict[str, Any]:
    """Own the target root exclusively for one idempotent worker invocation."""
    if not _PROCESS_WORKER_LOCK.acquire(blocking=False):
        raise RuntimeError("manifest workers require separate processes for concurrent data roots")
    try:
        return _run_with_root_lock(manifest_path=manifest_path, data_root=data_root, config_dir=config_dir, workers=workers, threads=threads, fastq_threads=fastq_threads, compression_threads=compression_threads, validation_slots=validation_slots, quant_slots=quant_slots, fasterq_slots=fasterq_slots, max_in_flight=max_in_flight, task_selection=task_selection)
    finally:
        _PROCESS_WORKER_LOCK.release()

def _run_with_root_lock(*, manifest_path: Path, data_root: Path, config_dir: Path, workers: int, threads: int, fastq_threads: int, compression_threads: int, validation_slots: int, quant_slots: int | None = None, fasterq_slots: int | None = None, max_in_flight: int | None = None, task_selection: Path | None = None) -> dict[str, Any]:
    data_root.mkdir(parents=True, exist_ok=True)
    with (data_root / ".acquisition.lock").open("a") as owned_lock:
        fcntl.flock(owned_lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        return _run_manifest_owned(
            manifest_path=manifest_path,
            data_root=data_root,
            config_dir=config_dir,
            workers=workers,
            threads=threads,
            fastq_threads=fastq_threads,
            compression_threads=compression_threads,
            validation_slots=validation_slots,
            quant_slots=quant_slots,
            fasterq_slots=fasterq_slots,
            max_in_flight=max_in_flight,
            task_selection=task_selection,
        )
