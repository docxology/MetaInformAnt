"""Bounded, idempotent acquisition and quantification on local or AWS hosts."""

from __future__ import annotations

import fcntl
import json
import os
import threading
import time
from collections import Counter
from collections.abc import Callable, Sequence
from dataclasses import asdict, replace
from pathlib import Path
from typing import Any
from uuid import uuid4

from metainformant.rna.engine.acquisition_allocation import verify_local_allocation
from metainformant.rna.engine.acquisition_lanes import (
    AdmissionDecision,
    run_acquisition_lanes,
)
from metainformant.rna.engine.acquisition_manifest import (
    load_snapshot,
    load_task_selection,
    sha256_file,
    verify_input_files,
    verify_worker_configs,
)
from metainformant.rna.engine.acquisition_prerequisites import (
    classify_task_prerequisites,
    require_worker_prerequisites,
)
from metainformant.rna.engine.acquisition_references import prepare_reference_inputs
from metainformant.rna.engine.acquisition_sample import execute_manifest_task, reuse_manifest_quantification
from metainformant.rna.engine.fastq_compression import compression_level
from metainformant.rna.engine.progress_db import configured_progress_db_path
from metainformant.rna.engine.streaming_orchestrator import (
    StreamingPipelineOrchestrator,
    _cached_sra_validation_result,
    _local_sra_path,
    build_pipeline_resource_profile,
)


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
    on_quantified: Callable[[str, str], None] | None = None,
    should_stop: Callable[[], bool] | None = None,
) -> dict[str, Any]:
    """Run all manifest tasks with bounded idempotent concurrency."""

    snapshot, manifest_tasks = load_snapshot(manifest_path)
    from metainformant.rna.engine.acquisition_size_overlay import configured_size_overlay
    size_overlay = configured_size_overlay(manifest_path, manifest_tasks)
    frontier_size = int(os.environ.get("AMALGKIT_LOCAL_FRONTIER_SIZE", "0"))
    frontier_raw_bytes = int(os.environ.get("AMALGKIT_LOCAL_RAW_RESERVATION_BYTES", "0"))
    if not 0 <= frontier_size <= 256 or frontier_raw_bytes < 0:
        raise ValueError("invalid local frontier bounds")
    if size_overlay.path is not None and not frontier_size:
        raise ValueError("Size overlays require the bounded local acquisition frontier")
    if frontier_size and (
        not frontier_raw_bytes
        or os.environ.get("AMALGKIT_DURABLE_BUCKET")
        or os.environ.get("AMALGKIT_WORKER_DEADLINE_EPOCH")
    ):
        raise ValueError("local frontier requires a raw-byte bound and cannot run cloud/deadline allocations")
    selection, tasks = load_task_selection(
        task_selection,
        snapshot,
        manifest_tasks,
        snapshot_sha256=sha256_file(manifest_path.with_name("snapshot.json")),
    )
    verify_local_allocation(
        task_selection,
        manifest_path,
        frozenset(task["task_id"] for task in manifest_tasks),
        frozenset(task["task_id"] for task in tasks),
        snapshot.get("inventory_sha256"),
    )
    profile = build_pipeline_resource_profile(
        min(workers, max(1, len(tasks))),
        threads,
        quant_slots=quant_slots,
        fasterq_slots=fasterq_slots,
        fasterq_threads=fastq_threads,
        compression_threads=compression_threads,
        validation_slots=validation_slots,
        max_in_flight=max_in_flight,
    )
    verify_input_files(snapshot, manifest_path.parent)
    selected_compression_level = compression_level()
    verify_worker_configs(snapshot, tasks, config_dir)
    os.environ["AMALGKIT_DATA_ROOT"] = str(data_root.resolve())
    os.environ.setdefault("AMALGKIT_RECLAIM_RAW_AFTER_QUANT", "no")
    os.environ.setdefault("AMALGKIT_MIN_EXTERNAL_FREE_GB", "8")
    os.environ.setdefault("AMALGKIT_MIN_SYSTEM_FREE_GB", "4")
    os.environ["AMALGKIT_PIPELINE_FASTQ_THREADS"] = str(fastq_threads)
    os.environ["AMALGKIT_PIPELINE_COMPRESSION_THREADS"] = str(compression_threads)
    os.environ["AMALGKIT_PIPELINE_VALIDATION_SLOTS"] = str(validation_slots)

    data_root.mkdir(parents=True, exist_ok=True)
    prerequisite_results = []
    for species in sorted({str(task["species"]) for task in tasks}):
        prerequisite_results.extend(
            classify_task_prerequisites(
                metadata_path=data_root
                / species
                / "work/metadata/metadata_selected.tsv",
                tasks=[task for task in tasks if task["species"] == species],
            )
        )
    require_worker_prerequisites(prerequisite_results)
    reference_preflight = prepare_reference_inputs(
        data_root=data_root,
        config_dir=config_dir,
        tasks=tasks,
        snapshot=snapshot,
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
        db_path=configured_progress_db_path(data_root / "pipeline_progress.db", os.environ.get("AMALGKIT_PROGRESS_DB_PATH")),
    )
    try:
        orchestrator._resource_profile = profile
        orchestrator._quant_semaphore = threading.BoundedSemaphore(profile.quant_slots)
        orchestrator._fasterq_semaphore = threading.BoundedSemaphore(
            profile.fasterq_slots
        )
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

        def execute(task: dict[str, Any], defer: bool, *, acquisition_only: bool = False,
                    quantification_only: bool = False) -> dict[str, Any]:
            return execute_manifest_task(
                task,
                orchestrator,
                data_root,
                config_dir,
                durable_store,
                durable_cohort,
                profile.quant_threads_per_worker,
                defer_ncbi_fallback=defer,
                acquisition_only=acquisition_only,
                quantification_only=quantification_only,
            )

        raw_bytes_limit = int(os.environ.get("AMALGKIT_CLOUD_MAX_RAW_BYTES", "0"))
        deadline_text = os.environ.get("AMALGKIT_WORKER_DEADLINE_EPOCH", "")
        worker_deadline = float(deadline_text) if deadline_text else None
        from metainformant.rna.engine.acquisition_scheduling import PlanningAssumptions

        worker_planning = None
        if snapshot.get("planning_assumptions") is not None:
            worker_planning = PlanningAssumptions(**snapshot["planning_assumptions"])
            worker_planning = replace(
                worker_planning,
                quant_slots=min(
                    worker_planning.quant_slots,
                    profile.quant_slots,
                    profile.workers,
                    profile.max_in_flight,
                ),
                extraction_slots=min(
                    worker_planning.extraction_slots,
                    profile.fasterq_slots,
                    profile.workers,
                    profile.max_in_flight,
                ),
            )

        def admission(
            task: dict[str, Any], live: Sequence[dict[str, Any]]
        ) -> AdmissionDecision:
            from metainformant.rna.engine.acquisition_estimates import (
                AcquisitionEstimateError,
            )
            from metainformant.rna.engine.acquisition_scheduling import (
                fits_worker_deadline,
                task_workload,
                workload_seconds,
            )

            if raw_bytes_limit and worker_deadline is None:
                return AdmissionDecision(reason="cloud worker has no hard deadline")
            if worker_deadline is None:
                return AdmissionDecision()
            try:
                duration = task.get("planning_seconds", 0)
                reserved_duration = sum(t.get("planning_seconds", 0) for t in live)
                drain = task.get("planning_drain_seconds", 300)
                if worker_planning is not None:
                    duration = workload_seconds(
                        [task_workload(t) for t in [*live, task]],
                        worker_planning,
                        include_setup=False,
                    )
                    reserved_duration = 0
                    drain = worker_planning.drain_seconds
                if fits_worker_deadline(
                    now=time.time(),
                    deadline=worker_deadline,
                    task_seconds=duration,
                    reserved_seconds=reserved_duration,
                    drain_seconds=drain,
                ):
                    return AdmissionDecision()
                if live:
                    return AdmissionDecision(wait=True)
                return AdmissionDecision(
                    reason="insufficient remaining worker deadline including unfinished reservations"
                )
            except AcquisitionEstimateError as exc:
                return AdmissionDecision(reason=str(exc))

        def local_sra(task: dict[str, Any]) -> bool:
            out_dir = data_root / str(task["species"]) / "work/getfastq"
            accession = str(task["accession"])
            source = _local_sra_path(out_dir / accession, out_dir, accession)
            return (
                source is not None
                and _cached_sra_validation_result(source) is not False
            )

        def recover(task: dict[str, Any]) -> bool:
            return orchestrator._download_fastq_ncbi_only(
                str(task["accession"]),
                data_root / str(task["species"]) / "work/getfastq",
                task.get("expected_paired"),
            )

        def finished(task: dict[str, Any], result: dict[str, Any]) -> None:
            result.update(
                species=task["species"], task_id=task["task_id"], srr=task["accession"]
            )
            if (
                not result.get("quantified")
                and result.get("status") not in {"unresolved_deadline", "deferred_download"}
            ):
                orchestrator.db.set_state(
                    str(task["species"]),
                    str(task["accession"]),
                    "failed",
                    error=str(result.get("error", "sample execution failed")),
                )
            append_journal(result)
            results.append(result)
            if result.get("status") == "unresolved_deadline":
                counts["unresolved"] += 1
            elif result.get("status") == "deferred_download":
                counts["deferred_download"] += 1
            else:
                counts["quantified" if result.get("quantified") else "failed"] += 1
                if result.get("quantified"):
                    counts[
                        "reused" if result.get("skipped") else "newly_quantified"
                    ] += 1
                    if on_quantified is not None:
                        on_quantified(str(task["species"]), str(task["accession"]))

        dedicated_recovery = os.environ.get(
            "AMALGKIT_PIPELINE_DEDICATED_RECOVERY", "no"
        ).lower() in ("yes", "true", "1")
        if frontier_size:
            from metainformant.rna.engine.acquisition_frontier import run_acquisition_frontier

            if orchestrator.acquisition_slots < 1 or not dedicated_recovery:
                raise ValueError("local frontier requires bounded acquisition slots and dedicated recovery")
            frontier_manifest_hash = sha256_file(manifest_path)
            frontier_selection_hash = sha256_file(task_selection) if task_selection else None
            source_manifest = os.environ.get("AMALGKIT_RUNTIME_SOURCE_MANIFEST")
            frontier_source_hash = sha256_file(Path(source_manifest)) if source_manifest else None

            def transition(task: dict[str, Any], stage: str) -> None:
                with (run_directory / "frontier_events.jsonl").open("a", encoding="utf-8") as handle:
                    handle.write(json.dumps({
                        **size_overlay.receipt(),
                        "task_id": task["task_id"], "stage": stage, "at_epoch": time.time(),
                        "manifest_sha256": frontier_manifest_hash,
                        "selection_sha256": frontier_selection_hash,
                        "source_manifest_sha256": frontier_source_hash,
                        "frontier_size": frontier_size, "raw_byte_limit": frontier_raw_bytes,
                    }, sort_keys=True) + "\n")
                    handle.flush()
                    os.fsync(handle.fileno())

            lane_summary = run_acquisition_frontier(
                tasks, acquisition_slots=orchestrator.acquisition_slots,
                quant_slots=profile.quant_slots, recovery_slots=profile.fasterq_slots,
                max_in_flight=profile.max_in_flight, frontier_size=frontier_size,
                raw_bytes_limit=frontier_raw_bytes,
                acquire=lambda task, defer: execute(task, defer, acquisition_only=True),
                quantify=lambda task: execute(task, True, quantification_only=True),
                recover=recover, uses_local_sra=local_sra, admission=admission,
                finished=finished, should_stop=should_stop, transition=transition,
                reuse=lambda task: reuse_manifest_quantification(task, data_root, config_dir),
                size_overrides=size_overlay.sizes,
            )
        else:
            lane_summary = run_acquisition_lanes(
                tasks, workers=profile.workers, max_in_flight=profile.max_in_flight,
                recovery_slots=profile.fasterq_slots, execute=execute,
                recover=recover, uses_local_sra=local_sra, admission=admission,
                finished=finished, should_stop=should_stop, raw_bytes_limit=raw_bytes_limit,
                dedicated_recovery=dedicated_recovery,
            )

        results.sort(key=lambda result: str(result.get("task_id", "")))

        summary = {
            **size_overlay.receipt(),
            "schema": {
                "metainformant.hymenoptera.gcp_snapshot.v1": "metainformant.hymenoptera.gcp_worker_result.v1",
                "metainformant.rna.acquisition_snapshot.v1": "metainformant.rna.acquisition_result.v1",
            }[snapshot["schema"]],
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
                "compression_level": selected_compression_level,
                "ena_file_workers": orchestrator.ena_file_workers,
                "download_first_pass_seconds": orchestrator.download_first_pass_seconds,
                "local_ena_scratch_dir": str(orchestrator.local_ena_scratch_dir) if orchestrator.local_ena_scratch_dir else None,
                "ena_promotion_slots": 1,
                "strict_fastq_record_validation": orchestrator.record_validation_required,
                "effective_quant_threads": profile.effective_quant_threads,
                "peak_stage_threads": profile.peak_stage_threads,
            },
            "task_results_journal": str(journal_path),
            "task_results_journal_sha256": sha256_file(journal_path),
            "task_count": len(tasks),
            "unsubmitted_task_count": lane_summary["unsubmitted"],
            "stopped_early": bool(lane_summary["unsubmitted"]),
            "acquisition_lanes": {
                "dedicated_recovery": dedicated_recovery,
                **lane_summary,
            },
            "manifest_task_count": len(manifest_tasks),
            "counts": dict(counts),
            "results": results,
        }
        result_path = data_root / "cloud_worker_result.json"
        (run_directory / "result.json").write_text(
            json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )
        result_path.write_text(
            json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )
        return summary

    finally:
        orchestrator.db.close()


_PROCESS_WORKER_LOCK = threading.Lock()


def run_manifest(
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
    on_quantified: Callable[[str, str], None] | None = None,
    should_stop: Callable[[], bool] | None = None,
) -> dict[str, Any]:
    """Own the target root exclusively for one idempotent worker invocation."""
    if not _PROCESS_WORKER_LOCK.acquire(blocking=False):
        raise RuntimeError(
            "manifest workers require separate processes for concurrent data roots"
        )
    try:
        return _run_with_root_lock(
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
            on_quantified=on_quantified,
            should_stop=should_stop,
        )
    finally:
        _PROCESS_WORKER_LOCK.release()


def _run_with_root_lock(
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
    on_quantified: Callable[[str, str], None] | None = None,
    should_stop: Callable[[], bool] | None = None,
) -> dict[str, Any]:
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
            on_quantified=on_quantified,
            should_stop=should_stop,
        )
