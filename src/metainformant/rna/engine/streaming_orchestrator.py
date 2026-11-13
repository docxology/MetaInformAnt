"""Streaming RNA-seq Orchestrator implementation.

This module provides robust orchestration for running the amalgkit pipeline
with ENA-first download streaming, concurrent processing, and automatic cleanup.
"""

from __future__ import annotations

import concurrent.futures
import csv
import errno
import gzip
import heapq
import itertools
import json
import logging
import math
import os
import re
import shutil
import sqlite3
import stat
import subprocess
import tempfile
import threading
import time
import urllib.request
from collections.abc import Iterator, MutableMapping
from contextlib import contextmanager
from dataclasses import dataclass, fields, replace
from functools import lru_cache
from pathlib import Path
from typing import Any, Dict, List, Mapping, Optional, Sequence, cast

import yaml

from metainformant.core.utils import logging as log_utils
from metainformant.core.utils.hash import sha256_file
from metainformant.rna.amalgkit.sra_environment import (
    build_sra_environment,
    ensure_ncbi_settings,
    resolve_data_root,
)
from metainformant.rna.amalgkit.tissue_normalizer import apply_tissue_normalization
from metainformant.rna.core.sample_utils import find_quantification_file
from metainformant.rna.engine.fastq_compression import compression_level, pigz_command
from metainformant.rna.engine.campaign_lock import (
    EXIT_ENVIRONMENT,
    EXIT_INTERNAL,
    EXIT_LOCK_BUSY,
    EXIT_NEEDS_OPERATOR,
    EXIT_OK,
    EXIT_RETRYABLE,
    CampaignLockBusy,
)
from metainformant.rna.engine.campaign_lock import acquire as acquire_campaign_lock
from metainformant.rna.engine.acquisition_order import plan_acquisition_order
from metainformant.rna.engine.preflight import PreflightError, ReuseProbe, reuse_dry_run, run_campaign_preflight
from metainformant.rna.engine.producer_stop import ProducerStop
from metainformant.rna.engine.cohort import (
    CohortError,
    FrozenCohort,
    cohort_task_table_name,
    quant_batch_run,
    verify_cohort_configs,
    verify_species_reference,
    write_cohort_task_table,
)
from metainformant.rna.engine.progress_db import (
    GATE_BACKOFF,
    GATE_HELD,
    NON_ENVIRONMENT_OS_ERROR,
    RETRYABLE_ENVIRONMENT_CLASSES,
    ProgressDB,
    ReconcileScope,
    RetryPolicy,
    RetryRecord,
    classify_sample_error,
    retry_gate,
)
from metainformant.rna.engine.provenance import (
    QUANT_STATUS_CURRENT,
    QUANT_STATUS_VERSION_DRIFT,
    REFERENCE_MANIFEST_NAME,
    REUSABLE_QUANT_STATUSES,
    SpeciesReuseContext,
    _durable_fsync,
    _fsync_directory,
    bind_produced_output,
    build_species_reuse_context,
    classify_quantification,
    discard_reference_binding,
    is_current_metadata,
    is_reusable_quantification,
    pinned_amalgkit_runtime,
    write_metadata_provenance,
    write_quant_provenance,
)
from metainformant.rna.engine.raw_cleanup import reclaim_sample_raw_inputs
from metainformant.rna.engine.species import resolve_species_work_dir
from metainformant.rna.engine.reference_digest import IndexDigestCache, file_signature
from metainformant.rna.retrieval import process_groups
from metainformant.rna.retrieval.ena_downloader import (
    ENADownloader,
    _run_command_in_process_group,
    verify_gzip_integrity,
)
from metainformant.rna.retrieval.ena_staging import (
    FirstPassDeferred,
    StagingDeadline,
    StagingError,
    _check_deadline,
    _plain_path,
    download_staged_ena,
    prepare_sample_directory,
    validate_fastq_records,
)

# Determine project root if possible, or use relative paths
PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent.parent.parent
NESTED_PROJECT_ROOT = PROJECT_ROOT / "projects/hymenoptera_amalgkit"
CONFIG_DIR = (
    NESTED_PROJECT_ROOT / "config/amalgkit"
    if (NESTED_PROJECT_ROOT / "config/amalgkit").is_dir()
    else PROJECT_ROOT / "config/amalgkit"
)
# Keep the repository-local default for tests and fresh checkouts, while
# allowing production runs to route the large tree to the mounted volume.
LOG_DIR = Path(os.environ.get("AMALGKIT_LOG_DIR", str(PROJECT_ROOT / "output/amalgkit")))


logger = log_utils.get_logger(__name__)

REQUANTIFICATION_POLICIES = frozenset({"preserve", "version-drift", "all"})
DEFAULT_DISCOVERY_WORKERS = 4
# How often one run queues a task again after its acquisition admission or staging deadline expired.
DEFAULT_DEFERRAL_REQUEUES = 20

_INDEX_CACHE_LOCK = threading.Lock()
_INDEX_STAGE_LOCK = threading.Lock()
INDEX_STAGE_MANIFEST = ".metainformant_index_stage.json"
_SRA_VALIDATION_SCHEMA = 1

# Revision of the resume-reconcile verdict rules; part of every ReconcileScope so a rule change re-judges outputs.
RECONCILE_SCOPE_REVISION = "reconcile.v1"
# Consecutive environment-class failures (EPERM, EACCES, EROFS, ENXIO, ENOSPC) that stop admission of new work.
DEFAULT_ENVIRONMENT_FAILURE_TRIP = 5
_ENVIRONMENT_ERRNO_CLASSES = {
    errno.EPERM: "environment_write_denied",
    errno.EACCES: "environment_access_denied",
    errno.EROFS: "environment_read_only",
    errno.ENXIO: "environment_device_unavailable",
    errno.ENOSPC: "environment_no_space",
}
_MISSING_TOOL = re.compile(r"No such file or directory: '([^'/]+)'")
# sqlite3.OperationalError texts that mean the database's volume, not the producer, failed.
_SQLITE_ENVIRONMENT_MARKERS = (
    "disk i/o error",
    "unable to open database file",
    "attempt to write a readonly database",
    "database or disk is full",
)
# Seconds between PASSIVE WAL checkpoints during a run (AMALGKIT_WAL_CHECKPOINT_SECONDS): a killed producer then
# leaves at most this much WAL, and immutable read-only observers of the database file lag by at most this much.
WAL_CHECKPOINT_INTERVAL_SECONDS = 600.0
_TRUTHY = frozenset({"1", "true", "yes", "on"})


def _write_text_durably(destination: Path, text: str, temporary: Path) -> None:
    """Write ``text`` to ``temporary``, fsync it, rename it over ``destination`` and fsync the directory."""

    try:
        with open(temporary, "w", encoding="utf-8") as handle:
            handle.write(text)
            handle.flush()
            _durable_fsync(handle.fileno())
        os.replace(temporary, destination)
    except BaseException:
        temporary.unlink(missing_ok=True)
        raise
    _fsync_directory(destination.parent)


def _path_is_under(path: str | os.PathLike[str], root: Path) -> bool:
    try:
        candidate = os.path.abspath(os.fspath(path))
    except (TypeError, ValueError):
        return False
    base = os.path.abspath(os.fspath(root))
    return candidate == base or candidate.startswith(base.rstrip(os.sep) + os.sep)


def environment_failure_class(failure: BaseException | str | None, data_root: Path | None = None) -> str | None:
    """The environment failure class of an exception or stored error message, else ``None``.

    Environment failures are caused by the host, not the sample: EPERM (macOS
    TCC denial), EACCES, EROFS, ENXIO (the volume dropped off the bus) and
    ENOSPC, plus a missing executable. An ``OSError`` carrying a path counts
    only when the path is under ``data_root``; one without a path always
    counts (the 2026-09-03 TCC denials carried none).
    """

    if failure is None:
        return None
    if isinstance(failure, OSError) and failure.errno in _ENVIRONMENT_ERRNO_CLASSES:
        paths = [path for path in (failure.filename, failure.filename2) if path not in (None, "")]
        if data_root is None or not paths or any(_path_is_under(path, data_root) for path in paths):
            return _ENVIRONMENT_ERRNO_CLASSES[failure.errno]
        return None
    text = str(failure)
    failure_class = classify_sample_error(text)
    if failure_class == "environment_missing_tool":
        return failure_class if _MISSING_TOOL.search(text) else None
    return failure_class if failure_class in RETRYABLE_ENVIRONMENT_CLASSES else None


def data_root_unavailable(data_root: Path) -> str | None:
    """A reason when the data root cannot be stat-ed as a directory (unmounted volume, TCC denial), else ``None``."""

    try:
        info = os.stat(data_root)
    except OSError as exc:
        return environment_failure_class(exc) or "environment_device_unavailable"
    if not stat.S_ISDIR(info.st_mode):
        return "environment_device_unavailable"
    return None


class EnvironmentCircuitBreaker:
    """Stops admission after ``threshold`` consecutive environment-class failures.

    A lost volume permission or a dropped USB volume fails every remaining
    sample within minutes (2026-09-03: 6,173 rows). Rows hit by an environment
    failure go back to ``pending``; after ``threshold`` in a row no new work is
    admitted, in-flight work drains, and the producer exits with
    ``EXIT_ENVIRONMENT`` so the operator fixes the host instead of the cohort.
    """

    def __init__(self, threshold: int) -> None:
        if threshold < 1:
            raise ValueError("environment failure threshold must be at least 1")
        self.threshold = threshold
        self.consecutive = 0
        self.total = 0
        self.tripped = False
        self.events: list[dict[str, str]] = []
        self._lock = threading.Lock()

    def record_environment_failure(self, failure_class: str, detail: str) -> bool:
        """Count one environment failure; return True when this one trips the breaker."""

        with self._lock:
            self.consecutive += 1
            self.total += 1
            self.events.append({"class": failure_class, "detail": detail[:300]})
            del self.events[: -self.threshold * 2]
            if not self.tripped and self.consecutive >= self.threshold:
                self.tripped = True
                return True
            return False

    def record_other_outcome(self) -> None:
        with self._lock:
            self.consecutive = 0

    def state(self) -> dict[str, Any]:
        with self._lock:
            return {
                "threshold": self.threshold,
                "tripped": self.tripped,
                "consecutive": self.consecutive,
                "environment_failures": self.total,
                "recent": list(self.events[-self.threshold :]),
            }


def configure_producer_logging(log_dir: Path, *, prefix: str, level: str | None = None) -> Path | None:
    """Log INFO (``CORE_LOG_LEVEL`` overrides) from every module to ``<log_dir>/<prefix>_<utc>_<pid>.log``.

    Module loggers keep their own console handlers; this sets the root level
    (which they inherit) and adds one root file handler, so progress lines
    (Phase timings, ``[n/N] ... Done``) are never dropped again. Returns the
    log path, or ``None`` (with a warning) when the file cannot be created.
    """

    log_utils.configure_logging_from_env(level or "INFO")
    try:
        log_dir.mkdir(parents=True, exist_ok=True)
        path = log_dir / f"{prefix}_{time.strftime('%Y%m%dT%H%M%SZ', time.gmtime())}_{os.getpid()}.log"
        handler = logging.FileHandler(path)
    except OSError as exc:
        logger.warning("Cannot open a run log under %s: %s", log_dir, exc)
        return None
    handler.setFormatter(logging.Formatter("%(asctime)s | %(levelname)s | %(name)s | %(message)s"))
    handler.set_name(f"metainformant-run-log:{path}")
    logging.getLogger().addHandler(handler)
    return path


def _resume_verify_hashes() -> bool:
    """``AMALGKIT_RESUME_VERIFY_HASHES``: re-judge and re-hash every output at startup (a deep audit)."""

    return os.environ.get("AMALGKIT_RESUME_VERIFY_HASHES", "").strip().lower() in _TRUTHY


def _reconcile_scope_label(mode: str) -> str:
    """The non-content part of a reconcile scope: rule revision, mode and the pinned Amalgkit runtime."""

    runtime = pinned_amalgkit_runtime()
    return (
        f"{RECONCILE_SCOPE_REVISION}|{mode}|amalgkit={runtime['amalgkit_version']}"
        f"@{runtime['amalgkit_source_revision']}"
    )


class DiscoveryFailed(RuntimeError):
    """One or more species could not be discovered; no sample task was submitted."""


def producer_exit_code_for_exception(error: BaseException) -> int:
    """Map an exception that ended a producer run to its exit code (73 busy, 69 environment, 75 interrupted, 70)."""

    if isinstance(error, CampaignLockBusy):
        return EXIT_LOCK_BUSY
    if isinstance(error, (PreflightError, CohortError, DiscoveryFailed)):
        return EXIT_ENVIRONMENT
    if isinstance(error, KeyboardInterrupt):
        return EXIT_RETRYABLE
    if isinstance(error, OSError):
        return EXIT_ENVIRONMENT
    if isinstance(error, sqlite3.OperationalError) and any(
        marker in str(error).lower() for marker in _SQLITE_ENVIRONMENT_MARKERS
    ):
        # The progress database lives on the data volume: a dropped or full volume is an environment failure.
        return EXIT_ENVIRONMENT
    return EXIT_INTERNAL


@dataclass(frozen=True)
class PipelineResourceProfile:
    """Bounded resources assigned to one streaming campaign.

    ``workers`` is the number of sample tasks that may be active at once.
    Quantification has its own semaphore: ``quant_slots`` is the maximum
    number of Kallisto/Amalgkit subprocesses, and
    ``quant_threads_per_worker`` is the thread count passed to each one.  The
    download fallback and FASTQ compression are separate child processes, so
    their per-process and concurrent-slot budgets are recorded explicitly
    instead of being hidden hard-coded multipliers. ``validation_slots``
    limits full local gzip scans because they compete with Kallisto for the
    mounted volume. Keeping task workers separate from quant/fallback slots
    lets network I/O feed the CPU-bound stage without silently increasing the
    quantification or external-disk contention budget.
    """

    requested_workers: int
    requested_quant_threads: int
    workers: int
    quant_slots: int
    quant_threads_per_worker: int
    fasterq_threads: int
    fasterq_slots: int
    compression_threads: int
    validation_slots: int
    max_in_flight: int
    host_cpu_count: int
    # The fields below are resolved by :func:`build_host_resource_profile`; an explicit
    # ``build_pipeline_resource_profile`` keeps the historical behavior (0 = not managed by the profile).
    acquisition_slots: int = 0
    ena_file_workers: int = 1
    ena_max_connections: int = 0
    ena_api_connections: int = 0
    scheduling_order: str = "size_bands"
    raw_reservation_bytes: int = 0
    host_memory_bytes: int | None = None
    host_performance_cores: int | None = None
    host_efficiency_cores: int | None = None
    kallisto_index_bytes: int | None = None
    quant_memory_bytes_per_slot: int | None = None
    memory_quant_slot_cap: int | None = None
    source: str = "explicit"

    def as_dict(self) -> Dict[str, Any]:
        """The resolved profile, including derived totals, for logs and ``--print-profile``."""

        payload: Dict[str, Any] = {field.name: getattr(self, field.name) for field in fields(self)}
        payload["effective_quant_threads"] = self.effective_quant_threads
        payload["peak_stage_threads"] = self.peak_stage_threads
        return payload

    @property
    def effective_quant_threads(self) -> int:
        """Return the maximum simultaneous Amalgkit quant thread count."""

        return self.quant_slots * self.quant_threads_per_worker

    @property
    def peak_stage_threads(self) -> int:
        """Return the conservative process-thread ceiling across stages."""

        fallback_threads = self.fasterq_slots * max(
            self.fasterq_threads,
            self.compression_threads,
        )
        return max(self.effective_quant_threads, fallback_threads)


def _resource_int(value: Any, default: int, label: str) -> int:
    """Parse a positive resource setting, falling back with a warning."""

    if value is None or value == "":
        return default
    try:
        parsed = int(value)
    except (TypeError, ValueError):
        logger.warning("Invalid %s=%r; using %d", label, value, default)
        return default
    if parsed <= 0:
        logger.warning("Invalid %s=%r; using %d", label, value, default)
        return default
    return parsed


def _format_duration(seconds: int) -> str:
    """Render a timeout compactly for logs and persisted failure messages."""

    if seconds % 3600 == 0:
        return f"{seconds // 3600}h"
    if seconds % 60 == 0:
        return f"{seconds // 60}m"
    return f"{seconds}s"


def build_pipeline_resource_profile(
    workers: int,
    threads: int,
    *,
    quant_slots: Optional[int] = None,
    fasterq_threads: Optional[int] = None,
    fasterq_slots: Optional[int] = None,
    compression_threads: Optional[int] = None,
    validation_slots: Optional[int] = None,
    max_in_flight: Optional[int] = None,
    cpu_count: Optional[int] = None,
) -> PipelineResourceProfile:
    """Build a reproducible, bounded resource profile.

    ``threads`` is a total quantification budget, not a per-worker value.
    ``workers`` controls all active sample tasks (including local validation
    and acquisition), while a quantification semaphore limits Kallisto to at
    most ``min(workers, threads)`` slots.  Stage-specific settings can be
    supplied explicitly or through
    ``AMALGKIT_PIPELINE_QUANT_SLOTS``,
    ``AMALGKIT_PIPELINE_FASTQ_THREADS``,
    ``AMALGKIT_PIPELINE_FASTQ_SLOTS``,
    ``AMALGKIT_PIPELINE_COMPRESSION_THREADS``, and
    ``AMALGKIT_PIPELINE_VALIDATION_SLOTS``, and
    ``AMALGKIT_PIPELINE_MAX_IN_FLIGHT``.
    """

    if workers <= 0:
        raise ValueError("workers must be positive")
    if threads <= 0:
        raise ValueError("threads must be positive")

    host_cpu_count = max(1, cpu_count or (os.cpu_count() or 1))
    # Keep sample-task concurrency independent from the total quantification
    # budget.  This allows more I/O/validation tasks to feed a bounded set of
    # Kallisto processes without turning a low thread budget into hidden
    # quantification oversubscription.
    configured_quant_slots: Any = quant_slots
    if configured_quant_slots is None:
        configured_quant_slots = os.environ.get("AMALGKIT_PIPELINE_QUANT_SLOTS")
    quant_slots = min(
        workers,
        threads,
        _resource_int(configured_quant_slots, min(workers, threads), "AMALGKIT_PIPELINE_QUANT_SLOTS"),
    )
    quant_threads_per_worker = max(1, threads // quant_slots)

    default_fasterq_threads = max(1, min(4, quant_threads_per_worker))
    default_compression_threads = max(1, min(2, quant_threads_per_worker // 2 or 1))
    configured_fasterq: Any = fasterq_threads
    if configured_fasterq is None:
        configured_fasterq = os.environ.get("AMALGKIT_PIPELINE_FASTQ_THREADS")
    configured_fasterq_slots: Any = fasterq_slots
    if configured_fasterq_slots is None:
        configured_fasterq_slots = os.environ.get("AMALGKIT_PIPELINE_FASTQ_SLOTS")
    configured_compression: Any = compression_threads
    if configured_compression is None:
        configured_compression = os.environ.get("AMALGKIT_PIPELINE_COMPRESSION_THREADS")
    configured_validation: Any = validation_slots
    if configured_validation is None:
        configured_validation = os.environ.get("AMALGKIT_PIPELINE_VALIDATION_SLOTS")
    configured_in_flight: Any = max_in_flight
    if configured_in_flight is None:
        configured_in_flight = os.environ.get("AMALGKIT_PIPELINE_MAX_IN_FLIGHT")

    profile = PipelineResourceProfile(
        requested_workers=workers,
        requested_quant_threads=threads,
        workers=workers,
        quant_slots=quant_slots,
        quant_threads_per_worker=quant_threads_per_worker,
        fasterq_threads=_resource_int(
            configured_fasterq,
            default_fasterq_threads,
            "AMALGKIT_PIPELINE_FASTQ_THREADS",
        ),
        fasterq_slots=min(
            workers,
            _resource_int(
                configured_fasterq_slots,
                min(4, workers),
                "AMALGKIT_PIPELINE_FASTQ_SLOTS",
            ),
        ),
        compression_threads=_resource_int(
            configured_compression,
            default_compression_threads,
            "AMALGKIT_PIPELINE_COMPRESSION_THREADS",
        ),
        validation_slots=max(
            1,
            min(
                workers,
                _resource_int(
                    configured_validation,
                    min(4, workers),
                    "AMALGKIT_PIPELINE_VALIDATION_SLOTS",
                ),
            ),
        ),
        max_in_flight=max(
            1,
            _resource_int(
                configured_in_flight,
                workers * 2,
                "AMALGKIT_PIPELINE_MAX_IN_FLIGHT",
            ),
        ),
        host_cpu_count=host_cpu_count,
    )
    if profile.workers > profile.quant_slots:
        logger.info(
            "Using %d sample-task workers with %d quantification slots; local validation/acquisition "
            "can overlap with the bounded Kallisto budget",
            profile.workers,
            profile.quant_slots,
        )
    if profile.effective_quant_threads > profile.host_cpu_count:
        logger.warning(
            "Quantification budget is %d threads on a host reporting %d CPUs; "
            "expect contention unless this host has an intentionally larger scheduler allocation",
            profile.effective_quant_threads,
            profile.host_cpu_count,
        )
    if profile.peak_stage_threads > profile.host_cpu_count:
        logger.warning(
            "Peak stage process budget is %d threads on a host reporting %d CPUs; "
            "fallback extraction/compression may contend with quantification",
            profile.peak_stage_threads,
            profile.host_cpu_count,
        )
    return profile


# ---------------------------------------------------------------------------
# Host-derived campaign profile
# ---------------------------------------------------------------------------

SCHEDULING_ORDERS = ("size_bands", "smallest_first", "species_order")
"""Fresh-acquisition order: global size bands with species round-robin (the measured default), one global
smallest-first queue, or the historical species-by-species order."""

# Measured host and workload facts behind the host profile (the campaign host: 8 cores as 4P+4E, 16 GiB RAM,
# USB APFS data volume, about 270 Mbit/s downlink):
#
# * c19 run cfba5dce (2026-10-24 22:06Z to 10-25 23:04Z, 24.97 h): 808 new + 843 reused quantifications from
#   504.6 GiB of validated FASTQ, 32.4 samples/h and 6.03 MB/s on average, best 6 h 46.5 samples/h at 9.23 MB/s,
#   with 8 whole-sample acquisitions, 16 ENA file connections, 2 ENA API connections, 8 quant slots x 1 kallisto
#   thread, 2 record validators, 1 recovery (NCBI extraction) slot, a 48-task frontier, a 256 GiB raw
#   reservation, size-band ordering and scratch on the data root.
# * Profiles with 4 acquisitions (control-plan 76e1dea0, 20c86570, 21364178) gave 6.3-10 samples/h; 12 workers
#   with max_in_flight 24 and no frontier (e64e2c1a) gave 3.4 samples/h with 16 admission-queue exhaustions;
#   internal-SSD scratch with under 20 GiB internal free (6acba5f4) collapsed to 1.5 samples/h.
# * Quantification is not the bottleneck: on 2026-10-25, 433 single-thread kallisto runs had a median of 51 s
#   (mean 103.5 s), 12.45 thread-hours in 24 h, at most 5 running at once. The 2026-10-10 shadow benchmark
#   (SRR25343810) took 121.3 s at 2 threads, 50.3 s at 4 and 36.6 s at 8: extra threads per run buy little,
#   so slots get one thread each and the slot count, not the thread count, uses the cores.
# * kallisto 0.52.0 peak memory (``/usr/bin/time -l``, index load plus a small read set, 2026-11-12): 322 MB
#   for the 31.8 MB Apis index, 484-485 MB for the 71.5 MB Linepithema and 80.8 MB Temnothorax nylanderi
#   indexes (the largest of the 27 frozen indexes), i.e. about 216 MB + 3.3 x index bytes.
HOST_ACQUISITION_SLOTS = 8
HOST_ENA_FILE_WORKERS = 2
HOST_ENA_API_CONNECTIONS = 2
HOST_VALIDATION_SLOTS = 2
HOST_RECOVERY_SLOTS = 1
HOST_MAX_QUANT_SLOTS = 8
HOST_RAW_RESERVATION_BYTES = 256 * 1024**3
HOST_RAW_RESERVATION_FREE_FRACTION = 0.25
# Memory guard: concurrent kallisto runs may use at most this share of physical memory, each estimated at
# 256 MiB + 4 x index bytes (above the measured 216 MB + 3.3 x, for the EM over real read sets).
HOST_QUANT_MEMORY_FRACTION = 0.5
QUANT_MEMORY_BASE_BYTES = 256 * 1024**2
QUANT_MEMORY_PER_INDEX_BYTE = 4
DEFAULT_KALLISTO_INDEX_BYTES = 96 * 1024**2


@dataclass(frozen=True)
class HostFacts:
    """CPU and memory of the host the producer runs on (``None`` when the platform does not say)."""

    cpu_count: int
    memory_bytes: int | None = None
    performance_cores: int | None = None
    efficiency_cores: int | None = None


def _darwin_sysctl_int(name: str) -> int | None:
    """One integer ``sysctlbyname`` value on macOS, else ``None`` (no subprocess)."""

    if os.uname().sysname != "Darwin":
        return None
    try:
        import ctypes

        libc = ctypes.CDLL(None)
        value = ctypes.c_int64(0)
        size = ctypes.c_size_t(ctypes.sizeof(value))
        if libc.sysctlbyname(name.encode(), ctypes.byref(value), ctypes.byref(size), None, ctypes.c_size_t(0)):
            return None
        return int(value.value)
    except (OSError, AttributeError, ValueError):
        return None


def detect_host_facts() -> HostFacts:
    """Logical CPUs (``os.process_cpu_count``), physical memory and, on Apple silicon, P/E core counts."""

    counter = getattr(os, "process_cpu_count", None) or os.cpu_count
    cpu_count = max(1, counter() or 1)
    memory: int | None = None
    try:
        memory = int(os.sysconf("SC_PAGE_SIZE")) * int(os.sysconf("SC_PHYS_PAGES"))
    except (ValueError, OSError, AttributeError):
        memory = _darwin_sysctl_int("hw.memsize")
    return HostFacts(
        cpu_count=cpu_count,
        memory_bytes=memory if memory and memory > 0 else None,
        performance_cores=_darwin_sysctl_int("hw.perflevel0.physicalcpu"),
        efficiency_cores=_darwin_sysctl_int("hw.perflevel1.physicalcpu"),
    )


def quant_memory_bytes_per_slot(kallisto_index_bytes: int | None) -> int:
    """The memory one kallisto run is budgeted for (see :data:`QUANT_MEMORY_BASE_BYTES`)."""

    index_bytes = kallisto_index_bytes if kallisto_index_bytes and kallisto_index_bytes > 0 else None
    return QUANT_MEMORY_BASE_BYTES + QUANT_MEMORY_PER_INDEX_BYTE * (index_bytes or DEFAULT_KALLISTO_INDEX_BYTES)


def _explicit_env_or(explicit: Any, variable: str, default: Any) -> Any:
    """Precedence: an explicit argument, then a set environment variable, then the host ``default``."""

    if explicit is not None:
        return explicit
    configured = os.environ.get(variable, "").strip()
    return configured if configured else default


def _nonnegative_profile_int(value: Any, default: int, label: str, maximum: int | None = None) -> int:
    """Parse a setting where 0 means unbounded; invalid values fall back to ``default`` with a warning."""

    try:
        parsed = int(value)
    except (TypeError, ValueError):
        logger.warning("Invalid %s=%r; using %d", label, value, default)
        return default
    if parsed < 0 or (maximum is not None and parsed > maximum):
        logger.warning("Invalid %s=%r; using %d", label, value, default)
        return default
    return parsed


def build_host_resource_profile(
    *,
    host: HostFacts | None = None,
    kallisto_index_bytes: int | None = None,
    data_root_free_bytes: int | None = None,
    workers: int | None = None,
    threads: int | None = None,
    quant_slots: int | None = None,
    fasterq_threads: int | None = None,
    fasterq_slots: int | None = None,
    compression_threads: int | None = None,
    validation_slots: int | None = None,
    max_in_flight: int | None = None,
    acquisition_slots: int | None = None,
    ena_file_workers: int | None = None,
    ena_max_connections: int | None = None,
    ena_api_connections: int | None = None,
    scheduling_order: str | None = None,
    raw_reservation_bytes: int | None = None,
) -> PipelineResourceProfile:
    """The campaign profile for this host: measured defaults, overridable per value.

    Each value comes from its argument, else its environment variable, else
    the host default below (measurement sources in the comment above
    :data:`HOST_ACQUISITION_SLOTS`):

    * ``quant_slots`` = min(logical CPUs, 8, memory cap), one kallisto thread
      each, so total kallisto threads never exceed the cores. The memory cap
      is ``50% of RAM / (256 MiB + 4 x largest index bytes)`` (16 GiB and the
      80.8 MB nylanderi index: 13 slots, so 8 stands).
    * ``acquisition_slots`` = 8 whole-sample acquisitions
      (``AMALGKIT_PIPELINE_ACQUISITION_SLOTS``); ``ena_file_workers`` = 2, so
      both mates of a paired run transfer at once
      (``AMALGKIT_PIPELINE_ENA_FILE_WORKERS``); ``ena_max_connections`` = 16
      file and ``ena_api_connections`` = 2 metadata requests
      (``AMALGKIT_ENA_MAX_CONNECTIONS``/``AMALGKIT_ENA_API_CONNECTIONS``).
    * ``workers`` = acquisition slots + quant slots: every acquisition slot
      can be busy while each quant slot holds a sample whose inputs are
      ready; more threads would only queue on the acquisition semaphore.
      ``max_in_flight`` = ``workers``, so nothing waits pre-submitted and a
      newly discovered species re-sorts the whole unsubmitted queue.
    * ``validation_slots`` = 2 full FASTQ scans, ``fasterq_slots`` = 1
      recovery (NCBI extraction) with ``fasterq_threads`` = min(4, CPUs/2):
      both are USB-disk bound and compete with downloads.
    * ``raw_reservation_bytes`` = min(256 GiB, 25% of the data root's free
      space): the compressed FASTQ bytes all admitted samples may hold at once
      (``AMALGKIT_PIPELINE_RAW_RESERVATION_GB``).
    * ``scheduling_order`` = ``size_bands`` (``AMALGKIT_PIPELINE_SCHEDULING_ORDER``).
    """

    facts = host or detect_host_facts()
    cpu_count = max(1, facts.cpu_count)
    per_slot = quant_memory_bytes_per_slot(kallisto_index_bytes)
    memory_cap: int | None = None
    if facts.memory_bytes:
        memory_cap = max(1, int(facts.memory_bytes * HOST_QUANT_MEMORY_FRACTION) // per_slot)
    default_quant_slots = max(1, min(cpu_count, HOST_MAX_QUANT_SLOTS, memory_cap or HOST_MAX_QUANT_SLOTS))

    # 0 keeps the orchestrator's historical "unbounded" meaning for acquisitions and ENA connections.
    resolved_acquisition = _nonnegative_profile_int(
        _explicit_env_or(acquisition_slots, "AMALGKIT_PIPELINE_ACQUISITION_SLOTS", HOST_ACQUISITION_SLOTS),
        HOST_ACQUISITION_SLOTS,
        "AMALGKIT_PIPELINE_ACQUISITION_SLOTS",
        256,
    )
    resolved_file_workers = _resource_int(
        _explicit_env_or(ena_file_workers, "AMALGKIT_PIPELINE_ENA_FILE_WORKERS", HOST_ENA_FILE_WORKERS),
        HOST_ENA_FILE_WORKERS,
        "AMALGKIT_PIPELINE_ENA_FILE_WORKERS",
    )
    if resolved_file_workers not in (1, 2):
        raise ValueError("ENA file workers must be 1 or 2")
    resolved_quant_slots = _resource_int(
        _explicit_env_or(quant_slots, "AMALGKIT_PIPELINE_QUANT_SLOTS", default_quant_slots),
        default_quant_slots,
        "AMALGKIT_PIPELINE_QUANT_SLOTS",
    )
    if memory_cap is not None and resolved_quant_slots > memory_cap:
        logger.warning(
            "Capping quant slots from %d to %d: concurrent kallisto runs budgeted at %.0f MiB each may use at "
            "most %.0f%% of %.1f GiB RAM",
            resolved_quant_slots,
            memory_cap,
            per_slot / 1024**2,
            HOST_QUANT_MEMORY_FRACTION * 100,
            (facts.memory_bytes or 0) / 1024**3,
        )
        resolved_quant_slots = memory_cap
    # One kallisto thread per slot unless a total thread budget was asked for; the total never exceeds the cores
    # through the default, and build_pipeline_resource_profile warns when an explicit budget does.
    resolved_threads = threads if threads is not None else resolved_quant_slots
    # Every acquisition slot busy plus one sample per quant slot whose inputs are ready.
    resolved_workers = workers if workers is not None else max(1, resolved_acquisition) + resolved_quant_slots
    profile = build_pipeline_resource_profile(
        resolved_workers,
        resolved_threads,
        quant_slots=resolved_quant_slots,
        fasterq_threads=_explicit_env_or(
            fasterq_threads, "AMALGKIT_PIPELINE_FASTQ_THREADS", max(1, min(4, cpu_count // 2))
        ),
        fasterq_slots=_explicit_env_or(fasterq_slots, "AMALGKIT_PIPELINE_FASTQ_SLOTS", HOST_RECOVERY_SLOTS),
        compression_threads=_explicit_env_or(
            compression_threads, "AMALGKIT_PIPELINE_COMPRESSION_THREADS", max(1, min(2, cpu_count // 4))
        ),
        validation_slots=_explicit_env_or(
            validation_slots, "AMALGKIT_PIPELINE_VALIDATION_SLOTS", HOST_VALIDATION_SLOTS
        ),
        max_in_flight=_explicit_env_or(max_in_flight, "AMALGKIT_PIPELINE_MAX_IN_FLIGHT", resolved_workers),
        cpu_count=cpu_count,
    )
    order = str(_explicit_env_or(scheduling_order, "AMALGKIT_PIPELINE_SCHEDULING_ORDER", SCHEDULING_ORDERS[0]))
    if order not in SCHEDULING_ORDERS:
        raise ValueError(f"scheduling order must be one of {SCHEDULING_ORDERS}, not {order!r}")
    if raw_reservation_bytes is None:
        configured_gb = os.environ.get("AMALGKIT_PIPELINE_RAW_RESERVATION_GB", "").strip()
        if configured_gb:
            raw_reservation_bytes = int(float(configured_gb) * 1024**3)
        else:
            raw_reservation_bytes = HOST_RAW_RESERVATION_BYTES
            if data_root_free_bytes is not None:
                raw_reservation_bytes = min(
                    raw_reservation_bytes, int(data_root_free_bytes * HOST_RAW_RESERVATION_FREE_FRACTION)
                )
    if raw_reservation_bytes < 0:
        raise ValueError("raw reservation must be nonnegative")
    # Both mates of every admitted acquisition may hold an ENA file connection (8 x 2 = 16 in c19).
    default_connections = resolved_acquisition * resolved_file_workers
    connections = _nonnegative_profile_int(
        _explicit_env_or(ena_max_connections, "AMALGKIT_ENA_MAX_CONNECTIONS", default_connections),
        default_connections,
        "AMALGKIT_ENA_MAX_CONNECTIONS",
        256,
    )
    api_connections = _nonnegative_profile_int(
        _explicit_env_or(ena_api_connections, "AMALGKIT_ENA_API_CONNECTIONS", HOST_ENA_API_CONNECTIONS),
        HOST_ENA_API_CONNECTIONS,
        "AMALGKIT_ENA_API_CONNECTIONS",
        256,
    )
    return replace(
        profile,
        acquisition_slots=resolved_acquisition,
        ena_file_workers=resolved_file_workers,
        ena_max_connections=connections,
        ena_api_connections=api_connections,
        scheduling_order=order,
        raw_reservation_bytes=raw_reservation_bytes,
        host_memory_bytes=facts.memory_bytes,
        host_performance_cores=facts.performance_cores,
        host_efficiency_cores=facts.efficiency_cores,
        kallisto_index_bytes=kallisto_index_bytes,
        quant_memory_bytes_per_slot=per_slot,
        memory_quant_slot_cap=memory_cap,
        source="host",
    )


PROFILE_ENVIRONMENT = (
    ("acquisition_slots", "AMALGKIT_PIPELINE_ACQUISITION_SLOTS"),
    ("ena_file_workers", "AMALGKIT_PIPELINE_ENA_FILE_WORKERS"),
    ("ena_max_connections", "AMALGKIT_ENA_MAX_CONNECTIONS"),
    ("ena_api_connections", "AMALGKIT_ENA_API_CONNECTIONS"),
)
"""Profile fields the acquisition code reads from the environment (the ENA transport gates are process-wide)."""


@contextmanager
def resource_profile_environment(profile: PipelineResourceProfile) -> Iterator[None]:
    """Export a host profile's acquisition budgets for the duration of a run, restoring the previous values."""

    if profile.source != "host":
        yield
        return
    previous = {variable: os.environ.get(variable) for _, variable in PROFILE_ENVIRONMENT}
    try:
        for name, variable in PROFILE_ENVIRONMENT:
            os.environ[variable] = str(getattr(profile, name))
        yield
    finally:
        for variable, value in previous.items():
            if value is None:
                os.environ.pop(variable, None)
            else:
                os.environ[variable] = value


FAST_FAIL_ENA_ENVIRONMENT: Mapping[str, str] = {
    "AMALGKIT_PIPELINE_ENA_RETRIES": "0",
    "AMALGKIT_NCBI_PREFETCH_FIRST": "yes",
    "AMALGKIT_PIPELINE_DOWNLOAD_SPEED_LIMIT_BYTES": str(64 * 1024),
    "AMALGKIT_PIPELINE_DOWNLOAD_SPEED_TIME_SECONDS": "20",
}
"""The opt-in fast-fail ENA profile: no ENA retries, NCBI prefetch first, and a 64 KiB/s floor over 20 s.

ENA connections already give up after 30 s (curl ``--connect-timeout 30`` and
the 30 s ENA API timeout). It is not the default: every ENA hiccup then goes
straight to the single recovery (NCBI prefetch + fasterq-dump) slot, whose
extraction on the USB volume is slower than an ENA resume, and a 20 s stall
window ends transfers that the default (1 KiB/s over 600 s, 5 resumable
retries) completes. The c19 measurements used the default transfer settings.
"""


def apply_fast_fail_ena_environment(environ: MutableMapping[str, str] | None = None) -> Dict[str, str]:
    """Set the :data:`FAST_FAIL_ENA_ENVIRONMENT` values that are not already set; return what was applied."""

    target = os.environ if environ is None else environ
    applied: Dict[str, str] = {}
    for variable, value in FAST_FAIL_ENA_ENVIRONMENT.items():
        if not target.get(variable, "").strip():
            target[variable] = value
            applied[variable] = value
    return applied


def format_resource_profile(profile: PipelineResourceProfile) -> str:
    """One log line describing every resolved budget."""

    memory = f"{profile.host_memory_bytes / 1024**3:.1f} GiB" if profile.host_memory_bytes else "unknown"
    cores = (
        f"{profile.host_cpu_count} CPUs ({profile.host_performance_cores}P+{profile.host_efficiency_cores}E)"
        if profile.host_performance_cores is not None
        else f"{profile.host_cpu_count} CPUs"
    )
    per_slot = profile.quant_memory_bytes_per_slot
    quant_memory = f"{per_slot / 1024**2:.0f} MiB/slot" if per_slot else "unguarded"
    reservation = profile.raw_reservation_bytes
    raw_reservation = f"{reservation / 1024**3:.0f} GiB" if reservation else "off"
    return (
        f"resource profile [{profile.source}]: host {cores}, {memory} RAM; workers={profile.workers}, "
        f"max_in_flight={profile.max_in_flight}, acquisition_slots={profile.acquisition_slots or 'unbounded'}, "
        f"ena_file_workers={profile.ena_file_workers}, ena_connections={profile.ena_max_connections or 'unbounded'} "
        f"file/{profile.ena_api_connections or 'unbounded'} api, quant={profile.quant_slots}x"
        f"{profile.quant_threads_per_worker} threads ({profile.effective_quant_threads} total), "
        f"quant_memory={quant_memory} cap={profile.memory_quant_slot_cap or '-'}, "
        f"validation_slots={profile.validation_slots}, recovery={profile.fasterq_slots}x{profile.fasterq_threads} "
        f"threads, compression_threads={profile.compression_threads}, raw_reservation={raw_reservation}, "
        f"order={profile.scheduling_order}"
    )


def _species_name_from_config(config_name: str) -> str:
    """Derive the output species key from an amalgkit config filename."""
    return config_name.replace("amalgkit_", "").replace(".yaml", "")


def _data_root() -> Path:
    """Return the streaming data root, honoring the external volume override."""

    return resolve_data_root()


def _ensure_ncbi_settings_for_data_root(data_root: Path) -> Path:
    """Compatibility facade for the centralized campaign SRA environment."""

    return ensure_ncbi_settings(data_root)


def _species_work_dir(species_name: str) -> Path:
    """Return the active streaming work directory for a species.

    The external data root contains two historical aliases (Apis and
    Pogonomyrmex).  Prefer a directory that already contains current workflow
    evidence over an empty canonical directory created by an interrupted
    preparation attempt.  This keeps resume runs attached to the existing
    quantifications without copying or deleting data.
    """

    # One rule for producers, finalization and reports (rna.engine.species.resolve_species_work_dir); the
    # evidence markers are cheap existence checks so discovery stays responsive over thousands of rows.
    return resolve_species_work_dir(_data_root(), species_name)


def _remap_species_workspace_value(value: Any, species_name: str, active_dir_name: str) -> Any:
    """Remap repository-relative paths to an active external alias.

    Configurations remain canonical and reviewable in Git.  When a mounted
    data root uses an established alias, the runtime YAML must point every
    species-specific path (work, logs, FASTQ, merge, filters, and genome) at
    that alias.  Non-path values and shared/reference paths are left intact.
    """

    if isinstance(value, dict):
        return {key: _remap_species_workspace_value(item, species_name, active_dir_name) for key, item in value.items()}
    if isinstance(value, list):
        return [_remap_species_workspace_value(item, species_name, active_dir_name) for item in value]
    if not isinstance(value, str) or active_dir_name == species_name:
        return value

    prefixes = (f"output/amalgkit/{species_name}", f"./output/amalgkit/{species_name}")
    for prefix in prefixes:
        if value == prefix or value.startswith(prefix + "/"):
            suffix = value[len(prefix) :].lstrip("/")
            replacement = f"output/amalgkit/{active_dir_name}"
            return f"{replacement}/{suffix}" if suffix else replacement
    return value


def _runtime_config_path(config_path: Path, species_name: str) -> Path:
    """Return a config path whose workspace paths target the active alias.

    For canonical data directories the checked-in YAML is used directly.  An
    alias-specific YAML is written only under the active external work tree;
    it is deterministic, ignored by the repository, and allows all current
    workflow stages to share the same paths as streaming quantification.
    """

    active_work_dir = _species_work_dir(species_name)
    active_dir_name = active_work_dir.parent.name
    if active_dir_name == species_name:
        return config_path

    try:
        source = yaml.safe_load(config_path.read_text()) or {}
        remapped = _remap_species_workspace_value(source, species_name, active_dir_name)
        payload = yaml.safe_dump(remapped, sort_keys=False)
    except (OSError, yaml.YAMLError) as exc:
        logger.warning("Unable to construct alias runtime config for %s: %s", species_name, exc)
        return config_path

    runtime_dir = active_work_dir / ".metainformant_runtime"
    runtime_path = runtime_dir / config_path.name
    try:
        runtime_dir.mkdir(parents=True, exist_ok=True)
        if not runtime_path.exists() or runtime_path.read_text() != payload:
            runtime_path.write_text(payload)
        logger.info("Using alias-aware runtime config for %s: %s", species_name, runtime_path)
        return runtime_path
    except OSError as exc:
        logger.warning("Unable to persist alias runtime config for %s: %s", species_name, exc)
        return config_path


def _resolve_metadata_path(work_dir: Path) -> Path:
    """Prefer the selected cohort, then fall back to the discovery table.

    ``metadata.tsv`` is the discovery universe; ``metadata_selected.tsv`` is
    the Amalgkit selection output consumed by acquisition and quantification.
    """
    selected_path = work_dir / "metadata/metadata_selected.tsv"
    if selected_path.exists():
        return selected_path
    return work_dir / "metadata/metadata.tsv"


def _filter_metadata_by_size(df: Any, max_gb: float) -> Any:
    """Return metadata rows whose total_bases fit the per-sample size budget."""
    import pandas as pd

    filtered_df = df.copy()
    filtered_df["total_bases"] = pd.to_numeric(filtered_df["total_bases"], errors="coerce").fillna(0)
    max_bases = max_gb * 1e9
    return filtered_df[filtered_df["total_bases"] <= max_bases].copy().sort_values("total_bases")


def _sample_run_column(df: Any) -> str:
    """Resolve the sample accession column used by amalgkit metadata."""
    return "run" if "run" in df.columns else "run_accession"


_RUN_ACCESSION_PATTERN = re.compile(r"^(?:SRR|ERR|DRR)\d+$")

# These are the final raw-input names that Amalgkit can consume directly.  A
# ``.part`` file is intentionally a separate scheduler tier: it is not a
# valid quantification input, but resuming it is cheaper and safer than
# starting a new transfer.
# ``prefetch`` (sra-tools 3.4.1) writes ``<accession>.sralite`` instead of ``<accession>.sra`` when NCBI serves the
# SRA Lite normalized-quality form (its output-name format ``%s.sralite``); fasterq-dump reads both, and kallisto
# ignores base qualities, so a lite archive is a complete local source.
_SRA_ARCHIVE_SUFFIXES = (".sra", ".sra.cache", ".sralite")
_RAW_INPUT_SUFFIXES = (".fastq.gz", ".fastq", ".fq.gz", ".fq", *_SRA_ARCHIVE_SUFFIXES)
_RAW_PARTIAL_SUFFIXES = tuple(f"{suffix}.part" for suffix in _RAW_INPUT_SUFFIXES)
_FASTQ_INPUT_SUFFIXES = (".fastq.gz", ".fastq", ".fq.gz", ".fq")
_RAW_VALIDATION_SCHEMA = "metainformant.rna.raw_validation.v2"
_RAW_REJECTION_SCHEMA = "metainformant.rna.raw_validation_rejected.v1"


def _library_layout_is_paired(value: Any) -> bool | None:
    """Resolve an Amalgkit library-layout value without guessing silently."""

    if not isinstance(value, str):
        return None
    normalized = value.strip().lower().replace("_", "-")
    if normalized in {"paired", "paired-end", "pair", "pe"}:
        return True
    if normalized in {"single", "single-end", "se"}:
        return False
    return None


def _raw_fastq_path(sample_dir: Path, accession: str, mate: int | None = None) -> Path | None:
    """Return one non-empty local FASTQ path for an accession and mate."""

    stem = f"{accession}_{mate}" if mate is not None else accession
    for suffix in _FASTQ_INPUT_SUFFIXES:
        candidate = sample_dir / f"{stem}{suffix}"
        try:
            if candidate.is_file() and candidate.stat().st_size > 0:
                return candidate
        except OSError:
            continue
    return None


def _is_accession_fastq_file(path: Path, accession: str) -> bool:
    """Return whether a final FASTQ name belongs to the exact accession.

    Completed mates are useful independently resumable transfer checkpoints
    even though a paired library is not quantifiable until both are present.
    Diagnostic suffixes and ``.part`` files intentionally do not match.
    """

    stems = (accession, f"{accession}_1", f"{accession}_2")
    return any(path.name == f"{stem}{suffix}" for stem in stems for suffix in _FASTQ_INPUT_SUFFIXES)


def _raw_validation_marker(sample_dir: Path, accession: str) -> Path:
    """Return the durable validation marker for one sample's raw inputs.

    The marker lives beside ``getfastq`` under the species work directory, not
    inside the sample directory.  This keeps the sample directory removable
    after quantification while retaining a small audit trail for any future
    reacquisition.  The path is deterministic for both supported
    ``work/getfastq`` and ``work/fastq/getfastq`` layouts.
    """

    return sample_dir.parent.parent / "raw_validation" / f"{accession}.json"


def _raw_validation_identity(path: Path) -> dict[str, int]:
    identity = path.lstat()
    if not stat.S_ISREG(identity.st_mode):
        raise OSError(f"Refusing non-regular FASTQ validation input: {path}")
    return {
        "size": identity.st_size,
        "mtime_ns": identity.st_mtime_ns,
        "ctime_ns": identity.st_ctime_ns,
        "device": identity.st_dev,
        "inode": identity.st_ino,
    }


def _raw_validation_marker_matches(
    marker: Path,
    candidates: list[Path],
    expected_paired: bool | None,
    *,
    require_records: bool = False,
) -> bool:
    """Return whether a cached full validation still describes the files."""

    try:
        marker_stat = marker.stat()
        if marker.is_symlink() or not marker.is_file() or marker_stat.st_size > 65536:
            return False
        payload = json.loads(marker.read_text())
    except (OSError, ValueError, TypeError):
        return False
    if not isinstance(payload, dict) or payload.get("schema") != _RAW_VALIDATION_SCHEMA:
        return False
    if require_records and not str(payload.get("layout_validation", "")).startswith("full_fastq_record_scan"):
        return False
    declared_layout = "paired" if expected_paired is True else "single" if expected_paired is False else "unknown"
    recorded_layout = payload.get("declared_library_layout")
    if recorded_layout is not None and recorded_layout != declared_layout:
        return False
    effective_layout = "paired" if len(candidates) == 2 else "single"
    recorded_effective_layout = payload.get("effective_library_layout")
    if recorded_effective_layout is not None and recorded_effective_layout != effective_layout:
        return False
    expected_reconciliation = (
        # Historical witness name retained for compatibility. This records
        # the observed layout, not an invalidating runtime version.
        "amalgkit_0.16.38_paired_metadata_single_fastq"
        if declared_layout == "paired" and effective_layout == "single"
        else "none"
    )
    recorded_reconciliation = payload.get("layout_reconciliation")
    if recorded_reconciliation is not None and recorded_reconciliation != expected_reconciliation:
        return False
    # Older witnesses proved gzip integrity but did not prove that a lone
    # FASTQ under paired metadata was a complete read-1-only stream rather
    # than malformed or non-adjacent paired data. Force one bounded full
    # record scan before such a witness can be reused.
    if (
        expected_reconciliation == "amalgkit_0.16.38_paired_metadata_single_fastq"
        and payload.get("layout_validation") != "full_fastq_record_scan_read1_only"
    ):
        return False
    recorded = payload.get("files")
    if not isinstance(recorded, list) or len(recorded) != len(candidates):
        return False
    by_name = {entry.get("name"): entry for entry in recorded if isinstance(entry, dict)}
    if len(by_name) != len(candidates):
        return False
    for path in candidates:
        entry = by_name.get(path.name)
        if entry is None:
            return False
        try:
            identity = _raw_validation_identity(path)
        except OSError:
            return False
        if any(entry.get(key) != value for key, value in identity.items()):
            return False
    return True


def _write_raw_validation_marker(
    marker: Path,
    candidates: list[Path],
    expected_paired: bool | None,
    *,
    layout_validation: str = "payload_stream_and_filename_cardinality",
) -> None:
    """Atomically record a successful local FASTQ validation."""

    files = []
    for path in candidates:
        files.append({"name": path.name, **_raw_validation_identity(path)})
    declared_layout = "paired" if expected_paired is True else "single" if expected_paired is False else "unknown"
    effective_layout = "paired" if len(candidates) == 2 else "single"
    payload = {
        "schema": _RAW_VALIDATION_SCHEMA,
        "validation": "full_gzip_stream" if all(path.name.endswith(".gz") for path in candidates) else "nonempty_fastq",
        "declared_library_layout": declared_layout,
        "effective_library_layout": effective_layout,
        "layout_reconciliation": (
            "amalgkit_0.16.38_paired_metadata_single_fastq"
            if declared_layout == "paired" and effective_layout == "single"
            else "none"
        ),
        "layout_validation": layout_validation,
        "files": files,
    }
    # Content-deterministic payload: no wall-clock fields.  Any restart-varying
    # byte in this marker would make a still-valid cached validation look
    # rewritten and force one bounded full record rescan per sample on every
    # resume.  Recency lives in the orchestrator log and the progress DB.
    marker.parent.mkdir(exist_ok=True)
    if marker.is_symlink() or (marker.exists() and not marker.is_file()):
        raise OSError(f"Refusing unsafe raw validation witness: {marker}")
    with tempfile.NamedTemporaryFile(mode="w", prefix=f".{marker.name}.", dir=marker.parent, delete=False) as handle:
        temporary = Path(handle.name)
        try:
            handle.write(json.dumps(payload, indent=2, sort_keys=True) + "\n")
            handle.flush()
            _durable_fsync(handle.fileno())
            os.replace(temporary, marker)
        finally:
            temporary.unlink(missing_ok=True)
    _fsync_directory(marker.parent)


def _local_fastq_payload_is_valid(path: Path) -> bool:
    """Validate a local FASTQ file without treating plain FASTQ as gzip."""

    if path.name.endswith(".gz"):
        return verify_gzip_integrity(path)
    try:
        with path.open("rb") as handle:
            return bool(handle.read(1))
    except OSError:
        return False


def _paths_share_filesystem(first: Path, second: Path) -> bool:
    """Return whether two existing paths reside on the same filesystem."""

    return first.stat().st_dev == second.stat().st_dev


def _sra_source_fingerprint(source: Path) -> dict[str, int]:
    """Return the source fields that make an SRA validation witness current."""

    source_stat = source.stat()
    return {
        "size": source_stat.st_size,
        "mtime_ns": source_stat.st_mtime_ns,
    }


def _sra_validation_marker_path(source: Path) -> Path:
    """Return the hidden validation-witness path for one cached SRA archive."""

    return source.parent / ".metainformant_sra_validation" / f"{source.name}.json"


def _cached_sra_validation_result(source: Path) -> bool | None:
    """Read a current valid/invalid SRA witness, if one exists."""

    marker = _sra_validation_marker_path(source)
    try:
        payload = json.loads(marker.read_text(encoding="utf-8"))
        if payload.get("schema") != _SRA_VALIDATION_SCHEMA:
            return None
        if payload.get("source") != source.name:
            return None
        if payload.get("fingerprint") != _sra_source_fingerprint(source):
            return None
        result = payload.get("valid")
        return result if isinstance(result, bool) else None
    except (OSError, ValueError, TypeError):
        return None


def _write_sra_validation_result(
    source: Path,
    *,
    fingerprint: dict[str, int],
    valid: bool,
    returncode: int,
    detail: str,
) -> None:
    """Persist an atomic source-bound SRA validation witness."""

    marker = _sra_validation_marker_path(source)
    payload = {
        "schema": _SRA_VALIDATION_SCHEMA,
        "source": source.name,
        "fingerprint": fingerprint,
        "valid": valid,
        "validator": "vdb-validate",
        "returncode": returncode,
        "detail": detail[:2000],
    }
    # Content-deterministic payload: no wall-clock fields.  Any restart-varying
    # byte would make a still-valid vdb-validate witness look rewritten and
    # force an expensive re-validation of the cached SRA archive on every
    # resume.  Recency lives in the orchestrator log and the progress DB.
    marker.parent.mkdir(parents=True, exist_ok=True)
    temporary = marker.with_name(f".{marker.name}.{os.getpid()}.{threading.get_ident()}.tmp")
    _write_text_durably(marker, json.dumps(payload, indent=2, sort_keys=True) + "\n", temporary)


def _campaign_cached_sra_path(accession: str) -> Path | None:
    """Return a complete-looking accession file in the campaign VDB cache."""

    cache_dir = _data_root() / ".sra-cache" / "sra"
    for suffix in _SRA_ARCHIVE_SUFFIXES:
        candidate = cache_dir / f"{accession}{suffix}"
        try:
            if candidate.is_file() and candidate.stat().st_size > 0:
                return candidate
        except OSError:
            continue
    return None


def _local_sra_path(sample_dir: Path, out_dir: Path, accession: str) -> Path | None:
    """Find an exact local SRA archive in nested, flat, or campaign-cache layout.

    ``prefetch`` may materialize a complete accession as either ``.sra`` or
    ``.sra.cache`` under the campaign-local VDB cache.  Those files are valid
    local databases even though they are outside Amalgkit's per-species
    ``getfastq`` tree.  Prefer them before ENA so an interrupted or previous
    campaign does not pay for the same transfer twice.  Presence is only a
    fast-path hint; ``fasterq-dump`` remains the authoritative extraction and
    layout check, and a failed extraction falls through to normal recovery.
    """

    roots = (sample_dir / ".prefetch" / accession, _data_root() / ".sra-cache/prefetch" / accession)
    candidates = [
        root / nested / f"{accession}{suffix}"
        for root in roots
        for nested in ("", accession)
        for suffix in _SRA_ARCHIVE_SUFFIXES
    ]
    candidates.extend(
        directory / f"{accession}{suffix}" for directory in (sample_dir, out_dir) for suffix in (".sra", ".sralite")
    )
    cached_sra = _campaign_cached_sra_path(accession)
    if cached_sra is not None:
        candidates.append(cached_sra)
    seen: set[Path] = set()
    for candidate in candidates:
        if candidate in seen:
            continue
        seen.add(candidate)
        try:
            if candidate.is_file() and candidate.stat().st_size > 0:
                return candidate
        except OSError:
            continue
    return None


def _prefetch_no_archive_reason(output: str, cached: Path | None) -> str:
    """Explain a zero-exit ``prefetch`` that wrote nothing, from sra-tools 3.4.1's own messages.

    prefetch exits 0 without writing when it resolves the accession to a copy
    it already has (``'<acc>' is found locally``: the user-repository cache
    under ``<data root>/.sra-cache``) or when the run exceeds ``--max-size``
    (``... is larger than maximum allowed: skipped``).
    """

    if "is found locally" in output:
        if cached is not None:
            return (
                f": prefetch resolved its local copy {cached}, which this campaign's vdb-validate witness rejects; "
                "remove or rename that cache file to let prefetch fetch it again"
            )
        return ": prefetch resolved a local copy outside the campaign output directory"
    if "larger than maximum allowed" in output:
        return ": the run is larger than the free-space bound passed as --max-size"
    return ""


def _validated_local_fastq_inputs(
    sample_dir: Path,
    accession: str,
    expected_paired: bool | None,
    *,
    assume_payload_valid: bool = False,
    deadline_at: float | None = None,
    require_records: bool | None = None,
    fresh_record_scan: bool = False,
) -> list[Path]:
    """Find validated local FASTQs without making an ENA metadata request.

    Existing files are only adopted when the metadata layout is compatible:
    paired libraries use both numbered mates when present, while one
    unnumbered FASTQ is retained for Amalgkit's explicit
    paired-metadata/single-file reconciliation. A lone numbered mate remains
    incomplete. Single libraries require the unnumbered FASTQ, and unknown
    layouts are accepted only when the on-disk naming is unambiguous. Gzip
    streams are fully checked before bypassing the downloader; incomplete
    ``.part`` files therefore remain resumable rather than being mistaken for
    usable input.
    """

    _check_deadline(deadline_at)
    if not sample_dir.is_dir():
        return []

    mate_1 = _raw_fastq_path(sample_dir, accession, 1)
    mate_2 = _raw_fastq_path(sample_dir, accession, 2)
    unnumbered = _raw_fastq_path(sample_dir, accession)

    if expected_paired is True:
        if mate_1 is not None or mate_2 is not None:
            if mate_1 is None or mate_2 is None:
                return []
            candidates = [mate_1, mate_2]
        elif unnumbered is not None:
            candidates = [unnumbered]
        else:
            return []
    elif expected_paired is False:
        if unnumbered is None:
            return []
        candidates = [unnumbered]
    elif mate_1 is not None or mate_2 is not None:
        if mate_1 is None or mate_2 is None:
            return []
        candidates = [mate_1, mate_2]
    elif unnumbered is not None:
        candidates = [unnumbered]
    else:
        return []

    if require_records is None:
        require_records = os.environ.get("AMALGKIT_PIPELINE_VALIDATE_FASTQ_RECORDS", "no").lower() in {
            "1",
            "yes",
            "true",
        } or bool(os.environ.get("AMALGKIT_LOCAL_ENA_SCRATCH_DIR"))
    record_proof = fresh_record_scan or require_records or (deadline_at is not None and not assume_payload_valid)
    marker = _raw_validation_marker(sample_dir, accession)
    if not fresh_record_scan and _raw_validation_marker_matches(
        marker, candidates, expected_paired, require_records=record_proof
    ):
        return candidates
    if not fresh_record_scan:
        rejected = _raw_rejection_reason(sample_dir, accession, candidates)
        if rejected is not None:
            # The same bytes failed layout validation before: rescanning them
            # (a full multi-gigabyte read on the USB volume) cannot change the
            # answer. Any rewrite of the files invalidates this witness.
            logger.info("Skipping the rescan of %s: unchanged inputs were rejected before (%s)", accession, rejected)
            return []
    try:
        before_scan = [_raw_validation_identity(path) for path in candidates]
        if record_proof and not validate_fastq_records(candidates, deadline_at=deadline_at):
            return []
        if before_scan != [_raw_validation_identity(path) for path in candidates]:
            return []
    except OSError:
        return []

    layout_validation = "full_fastq_record_scan" if record_proof else "payload_stream_and_filename_cardinality"
    layout_proved_payload = record_proof
    if expected_paired is True and len(candidates) == 1 and candidates[0].name == f"{accession}.fastq.gz":
        # A completed one-file checkpoint may have been left by an interrupted
        # older run before the raw witness was written. Classify it here,
        # before adoption: true interleaved reads are split atomically, while
        # Amalgkit's read-1-only case is accepted only after a complete
        # FASTQ record scan. This path also centralizes classification for new
        # ENA and fasterq outputs, avoiding a second multi-gigabyte scan.
        try:
            normalized = _normalize_interleaved_paired_fastq_layout(
                candidates,
                accession,
                sample_dir,
                expected_paired,
                deadline_at=deadline_at,
            )
        except OSError as exc:
            # I/O errors (a dropped USB volume, a truncated gzip member) may not repeat: no negative witness.
            logger.warning(
                "Rejecting one-file paired input for %s after layout validation: %s",
                accession,
                exc,
            )
            return []
        except (ValueError, EOFError) as exc:
            logger.warning(
                "Rejecting one-file paired input for %s after layout validation: %s",
                accession,
                exc,
            )
            try:
                _write_raw_rejection_marker(
                    sample_dir,
                    accession,
                    candidates,
                    expected_paired,
                    f"one-file paired input failed layout validation: {exc}",
                )
            except OSError as marker_error:
                logger.debug("Could not persist the raw rejection witness for %s: %s", accession, marker_error)
            return []
        candidates = normalized
        layout_proved_payload = True
        layout_validation = (
            "full_fastq_record_scan_interleaved_split" if len(candidates) == 2 else "full_fastq_record_scan_read1_only"
        )

    if (
        not assume_payload_valid
        and not layout_proved_payload
        and not all(_local_fastq_payload_is_valid(path) for path in candidates)
    ):
        return []
    try:
        _write_raw_validation_marker(
            marker,
            candidates,
            expected_paired,
            layout_validation=layout_validation,
        )
    except OSError as exc:
        # A validation cache is an optimization, never a reason to reject a
        # usable read.  The next run will conservatively revalidate it.
        logger.debug("Could not persist raw validation marker %s: %s", marker, exc)
    try:
        _raw_rejection_marker(sample_dir, accession).unlink(missing_ok=True)
    except OSError:
        pass
    return candidates


def _acquired_input_bytes(sample_dir: Path, accession: str) -> int:
    """Bytes of the canonical FASTQ inputs present for ``accession`` (0 when none or unreadable)."""

    total = 0
    for mate in (None, 1, 2):
        path = _raw_fastq_path(sample_dir, accession, mate)
        if path is None:
            continue
        try:
            total += path.stat().st_size
        except OSError:
            continue
    return total


def _raw_rejection_marker(sample_dir: Path, accession: str) -> Path:
    """The negative witness beside :func:`_raw_validation_marker`: these exact files failed layout validation."""

    return sample_dir.parent.parent / "raw_validation" / f"{accession}.rejected.json"


def _raw_rejection_identity(path: Path) -> dict[str, int]:
    identity = path.lstat()
    if not stat.S_ISREG(identity.st_mode):
        raise OSError(f"Refusing non-regular FASTQ validation input: {path}")
    return {"inode": identity.st_ino, "size": identity.st_size, "mtime_ns": identity.st_mtime_ns}


def _raw_rejection_reason(sample_dir: Path, accession: str, candidates: list[Path] | None = None) -> str | None:
    """The recorded reason when a negative witness still describes the files on disk, else ``None``.

    The witness is keyed by each file's ``(inode, size, mtime_ns)``: replacing
    or rewriting a file (a new download, an extraction, an operator repair)
    invalidates it, so only byte-identical inputs skip the rescan.
    """

    marker = _raw_rejection_marker(sample_dir, accession)
    try:
        if marker.is_symlink() or not marker.is_file() or marker.stat().st_size > 65536:
            return None
        payload = json.loads(marker.read_text(encoding="utf-8"))
    except (OSError, ValueError, TypeError):
        return None
    if not isinstance(payload, dict) or payload.get("schema") != _RAW_REJECTION_SCHEMA:
        return None
    recorded = payload.get("files")
    if not isinstance(recorded, list) or not recorded:
        return None
    names = [entry.get("name") for entry in recorded if isinstance(entry, dict)]
    if len(names) != len(recorded) or not all(isinstance(name, str) and "/" not in name for name in names):
        return None
    if candidates is not None and sorted(str(name) for name in names) != sorted(path.name for path in candidates):
        return None
    for entry in recorded:
        name = entry["name"]
        try:
            identity = _raw_rejection_identity(sample_dir / name)
        except OSError:
            return None
        if any(entry.get(key) != value for key, value in identity.items()):
            return None
    reason = payload.get("reason")
    return reason if isinstance(reason, str) and reason else "layout validation rejected these files"


def _write_raw_rejection_marker(
    sample_dir: Path, accession: str, candidates: list[Path], expected_paired: bool | None, reason: str
) -> None:
    """Atomically record that ``candidates`` failed layout validation (a deterministic, content-only payload)."""

    marker = _raw_rejection_marker(sample_dir, accession)
    payload = {
        "schema": _RAW_REJECTION_SCHEMA,
        "accession": accession,
        "declared_library_layout": (
            "paired" if expected_paired is True else "single" if expected_paired is False else "unknown"
        ),
        "files": [{"name": path.name, **_raw_rejection_identity(path)} for path in candidates],
        "reason": reason[:500],
    }
    marker.parent.mkdir(exist_ok=True)
    if marker.is_symlink() or (marker.exists() and not marker.is_file()):
        raise OSError(f"Refusing unsafe raw validation witness: {marker}")
    temporary = marker.with_name(f".{marker.name}.{os.getpid()}.{threading.get_ident()}.tmp")
    _write_text_durably(marker, json.dumps(payload, indent=2, sort_keys=True) + "\n", temporary)


def _layout_contradiction(files: list[Path], accession: str, expected_paired: bool | None) -> str | None:
    """Describe downloaded FASTQ names that cannot satisfy the declared layout, else ``None``."""

    names = sorted(path.name for path in files)
    numbered = [name for name in names if re.match(rf"^{re.escape(accession)}_[12]\.", name)]
    unnumbered = [name for name in names if name.startswith(f"{accession}.")]
    if expected_paired is False and numbered and not unnumbered:
        return f"declared single, source delivered mate files {', '.join(numbered)}"
    if expected_paired is True and len(numbered) == 1 and not unnumbered:
        return f"declared paired, source delivered one mate {numbered[0]}"
    return None


def _promote_staged_fastq_file(staged: Path, destination: Path) -> None:
    """Promote one validated FASTQ atomically, including across filesystems.

    ``os.replace`` is atomic on one filesystem. Local extraction scratch can
    live on another volume, so an ``EXDEV`` result falls back to an exact-size
    copy into a hidden destination-side temporary file followed by an atomic
    rename. The validated source is removed only after promotion succeeds.
    """

    try:
        os.replace(staged, destination)
        return
    except OSError as exc:
        if exc.errno != errno.EXDEV:
            raise

    temporary = destination.with_name(f".{destination.name}.{os.getpid()}.{threading.get_ident()}.copying")
    try:
        temporary.unlink(missing_ok=True)
        source_size = staged.stat().st_size
        shutil.copy2(staged, temporary)
        if temporary.stat().st_size != source_size:
            raise OSError(f"staged FASTQ size mismatch for {staged.name}")
        os.replace(temporary, destination)
        staged.unlink()
    finally:
        temporary.unlink(missing_ok=True)


def _promote_validated_fastq_inputs(
    staged_inputs: list[Path], sample_dir: Path, accession: str, expected_paired: bool | None
) -> bool:
    """Promote layout-validated extraction outputs and validate the actual destination."""

    # Prove the exact retained/extracted union before moving or discarding
    # anything. Individually intact gzip streams can still be unrelated mates.
    prospective_by_stem: dict[str, Path] = {}
    retained_destinations: set[Path] = set()
    for mate in (None, 1, 2):
        existing = _raw_fastq_path(sample_dir, accession, mate)
        if existing is None or not _local_fastq_payload_is_valid(existing):
            continue
        if expected_paired is True and mate is None:
            try:
                _validate_single_fastq_stream(existing)
            except (OSError, EOFError, ValueError):
                continue
        stem = f"{accession}_{mate}" if mate is not None else accession
        prospective_by_stem[stem] = existing
        retained_destinations.add(existing)
    for staged in staged_inputs:
        staged.stat()
        _plain_path(staged)
        if not _is_accession_fastq_file(staged, accession):
            return False
        stem = next(
            staged.name.removesuffix(suffix) for suffix in _FASTQ_INPUT_SUFFIXES if staged.name.endswith(suffix)
        )
        destination = sample_dir / staged.name
        _plain_path(destination)
        if destination not in retained_destinations:
            prospective_by_stem[stem] = staged
    first, second = (prospective_by_stem.get(f"{accession}_{mate}") for mate in (1, 2))
    single = prospective_by_stem.get(accession)
    if expected_paired is False:
        prospective = [single] if single is not None else []
    elif first is not None or second is not None:
        prospective = [first, second] if first is not None and second is not None else []
    else:
        prospective = [single] if single is not None else []
    if not prospective or not validate_fastq_records(prospective):
        logger.error("Conflicting extracted/retained FASTQs for %s; all inputs preserved", accession)
        return False

    for staged in staged_inputs:
        destination = sample_dir / staged.name
        if destination.exists():
            existing_is_valid = destination in retained_destinations
            if existing_is_valid:
                logger.info(
                    "Retaining existing valid FASTQ while promoting complementary " "fallback output for %s: %s",
                    accession,
                    destination,
                )
                staged.unlink(missing_ok=True)
                continue
            invalid = destination.with_name(f"{destination.name}.invalid")
            counter = 1
            while invalid.exists():
                invalid = destination.with_name(f"{destination.name}.invalid.{counter}")
                counter += 1
            destination.replace(invalid)
        _promote_staged_fastq_file(staged, destination)

    # The actual destination may include retained files with different bytes.
    # Preserve its full layout check; a staging witness cannot authorize them.
    return bool(
        _validated_local_fastq_inputs(
            sample_dir,
            accession,
            expected_paired,
            assume_payload_valid=True,
        )
    )


def _fastq_pair_key(header: str) -> str | None:
    """Return a mate-insensitive identifier from a FASTQ header."""

    if not header.startswith("@"):
        return None
    fields = header[1:].strip().split()
    if not fields:
        return None
    # ENA's interleaved paired FASTQs use either ``read/1 read/2`` or the
    # Casava ``read 1:N:... read 2:N:...`` convention. The first token is
    # stable for both forms. Do not strip ``.1``/``.2``: SRA uses those as
    # consecutive spot identifiers, not mate suffixes.
    return re.sub(r"/[12]$", "", fields[0])


class _PairedFastqNotInterleaved(ValueError):
    """The stream is structurally FASTQ but adjacent records are not mates."""


def _fastq_header_explicit_mate(header: str) -> int | None:
    """Return an explicit mate number from common FASTQ header conventions."""

    fields = header[1:].strip().split()
    if not fields:
        return None
    slash_suffix = re.search(r"/([12])$", fields[0])
    if slash_suffix:
        return int(slash_suffix.group(1))
    if len(fields) > 1:
        casava = re.match(r"([12]):[YN]:", fields[1])
        if casava:
            return int(casava.group(1))
        sra_suffix = re.search(r"/([12])$", fields[1])
        if sra_suffix:
            return int(sra_suffix.group(1))
    return None


def _validate_single_fastq_stream(source: Path, *, deadline_at: float | None = None) -> int:
    """Validate a complete one-read stream before layout reconciliation."""

    record_count = 0
    with gzip.open(source, "rt", encoding="utf-8", newline="") as reader:
        while True:
            _check_deadline(deadline_at)
            record = [reader.readline() for _ in range(4)]
            if not record[0]:
                break
            if any(not line for line in record):
                raise ValueError(f"truncated FASTQ record in {source.name}")
            if not record[0].startswith("@"):
                raise ValueError(f"invalid FASTQ header in {source.name}")
            if not record[2].startswith("+"):
                raise ValueError(f"invalid FASTQ separator in {source.name}")
            sequence = record[1].rstrip("\r\n")
            quality = record[3].rstrip("\r\n")
            if not sequence or len(sequence) != len(quality):
                raise ValueError(f"sequence/quality length mismatch in {source.name}")
            if _fastq_header_explicit_mate(record[0]) == 2:
                raise ValueError(f"non-adjacent mate-2 record in one-file paired input {source.name}")
            record_count += 1
    if record_count == 0:
        raise ValueError(f"FASTQ contains no complete records: {source.name}")
    return record_count


def _split_interleaved_fastq(
    source: Path,
    accession: str,
    sample_dir: Path,
    *,
    deadline_at: float | None = None,
) -> list[Path]:
    """Split one interleaved paired FASTQ into numbered mates.

    ENA and some valid SRA archives can expose a paired run as one interleaved
    ``<accession>.fastq.gz`` object even though the portal reports
    ``library_layout=PAIRED``.  Amalgkit/Kallisto require two mate files, so
    convert that representation only after checking every adjacent record
    pair has matching read IDs.  Temporary outputs are never adopted as
    inputs, which keeps interruption resume-safe.
    """

    if not source.is_file() or source.name != f"{accession}.fastq.gz":
        raise ValueError(f"not an unnumbered ENA FASTQ for {accession}: {source.name}")

    output_1 = sample_dir / f"{accession}_1.fastq.gz"
    output_2 = sample_dir / f"{accession}_2.fastq.gz"
    temporary_1 = output_1.with_name(f".{output_1.name}.part")
    temporary_2 = output_2.with_name(f".{output_2.name}.part")
    if output_1.exists() or output_2.exists():
        raise FileExistsError(f"cannot split interleaved FASTQ over existing mate output(s) for {accession}")

    pair_count = 0
    try:
        with gzip.open(source, "rt", encoding="utf-8", newline="") as reader:
            with gzip.open(temporary_1, "wt", compresslevel=1, encoding="utf-8", newline="") as writer_1:
                with gzip.open(temporary_2, "wt", compresslevel=1, encoding="utf-8", newline="") as writer_2:
                    while True:
                        _check_deadline(deadline_at)
                        first = [reader.readline() for _ in range(4)]
                        if not first[0]:
                            break
                        second = [reader.readline() for _ in range(4)]
                        if any(not line for line in (*first, *second)):
                            raise ValueError(f"truncated FASTQ record in {source.name}")
                        if not first[0].startswith("@") or not second[0].startswith("@"):
                            raise ValueError(f"invalid FASTQ header in {source.name}")
                        if not first[2].startswith("+") or not second[2].startswith("+"):
                            raise ValueError(f"invalid FASTQ separator in {source.name}")
                        for record in (first, second):
                            sequence = record[1].rstrip("\r\n")
                            quality = record[3].rstrip("\r\n")
                            if not sequence or len(sequence) != len(quality):
                                raise ValueError(f"sequence/quality length mismatch in {source.name}")
                        if _fastq_pair_key(first[0]) != _fastq_pair_key(second[0]):
                            raise _PairedFastqNotInterleaved(
                                f"adjacent FASTQ records are not a mate pair in {source.name}"
                            )
                        first_mate = _fastq_header_explicit_mate(first[0])
                        second_mate = _fastq_header_explicit_mate(second[0])
                        if first_mate == 2 or second_mate == 1:
                            raise ValueError(f"explicit FASTQ mate order is reversed in {source.name}")
                        if first_mate is not None and second_mate is not None and (first_mate, second_mate) != (1, 2):
                            raise ValueError(f"explicit FASTQ mate labels are inconsistent in {source.name}")
                        writer_1.writelines(first)
                        writer_2.writelines(second)
                        pair_count += 1
        if pair_count == 0:
            raise ValueError(f"interleaved FASTQ contains no complete pairs: {source.name}")
        temporary_1.replace(output_1)
        temporary_2.replace(output_2)
        source.unlink()
    except Exception:
        temporary_1.unlink(missing_ok=True)
        temporary_2.unlink(missing_ok=True)
        raise
    return [output_1, output_2]


def _normalize_interleaved_paired_fastq_layout(
    fastq_files: list[Path],
    accession: str,
    sample_dir: Path,
    expected_paired: bool | None,
    *,
    deadline_at: float | None = None,
) -> list[Path]:
    """Normalize a one-file interleaved paired representation."""

    if expected_paired is not True or len(fastq_files) != 1:
        return fastq_files
    source = fastq_files[0]
    if source.name != f"{accession}.fastq.gz":
        return fastq_files
    try:
        return _split_interleaved_fastq(source, accession, sample_dir, deadline_at=deadline_at)
    except _PairedFastqNotInterleaved:
        record_count = _validate_single_fastq_stream(source, deadline_at=deadline_at)
        logger.warning(
            "Run %s is declared paired but contains one validated read-1-only "
            "FASTQ (%d records); retaining the unnumbered file for Amalgkit "
            "0.16.60 layout reconciliation.",
            accession,
            record_count,
        )
        return fastq_files


def _metadata_requires_fastq_input_stats(metadata_path: Path, accession: str) -> bool:
    """Return whether Amalgkit needs read-derived statistics for one accession.

    Amalgkit 0.16.60 infers ``spot_length`` from ``total_bases / total_spots``
    when the public SRA metadata omits the explicit field.  Re-counting a
    validated multi-gigabyte FASTQ merely because ``spot_length`` is missing
    therefore adds substantial I/O/CPU without changing the quantification
    inputs.  Direct FASTQ statistics remain necessary when either count is
    absent or unusable, because ``get_sra_stat`` cannot construct its
    quantification contract without them.
    """

    try:
        with metadata_path.open(newline="", encoding="utf-8") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                if str(row.get("run", "")).strip() != accession:
                    continue
                for column in ("total_spots", "total_bases"):
                    try:
                        value = float(row.get(column, ""))
                    except (TypeError, ValueError):
                        return True
                    if not math.isfinite(value) or value <= 0:
                        return True
                return False
    except (OSError, csv.Error):
        return False
    return False


def _count_fastq_records_and_bases(path: Path) -> tuple[int, int]:
    """Count records and sequence bases in one validated FASTQ payload."""

    opener = gzip.open if path.name.endswith(".gz") else open
    records = 0
    bases = 0
    with opener(path, "rb") as handle:
        while True:
            header = handle.readline()
            if not header:
                break
            sequence = handle.readline()
            plus = handle.readline()
            quality = handle.readline()
            if not sequence or not plus or not quality:
                raise ValueError(f"truncated FASTQ record: {path}")
            records += 1
            bases += len(sequence.rstrip(b"\r\n"))
    return records, bases


def _write_fastq_input_stats(sample_dir: Path, accession: str) -> bool:
    """Write a minimal Amalgkit quant-only stats row for direct FASTQ inputs.

    Direct ENA acquisition does not run Amalgkit's ``getfastq`` stage, so a
    metadata row with zero or missing SRA statistics otherwise fails before
    Kallisto can inspect the already validated FASTQ.  The row is derived from
    the immutable local payload, written atomically, and consumed only as a
    quantification fallback; source metadata is never modified.
    """

    stats_path = sample_dir / "getfastq_stats.tsv"
    try:
        if stats_path.is_file():
            with stats_path.open(newline="", encoding="utf-8") as handle:
                for row in csv.DictReader(handle, delimiter="\t"):
                    if str(row.get("run", "")).strip() != accession:
                        continue
                    try:
                        if float(row.get("num_written", "")) > 0:
                            return True
                    except (TypeError, ValueError):
                        break
    except (OSError, csv.Error):
        logger.debug("Unable to inspect existing FASTQ statistics for %s", accession)

    try:
        fastq_files = sorted(
            path
            for path in sample_dir.iterdir()
            if path.is_file() and path.name.startswith(accession) and path.name.endswith(_FASTQ_INPUT_SUFFIXES)
        )
    except OSError:
        return False
    if not fastq_files:
        return False

    try:
        measured = [_count_fastq_records_and_bases(path) for path in fastq_files]
    except (OSError, EOFError, gzip.BadGzipFile, ValueError) as exc:
        logger.warning("Unable to derive FASTQ statistics for %s: %s", accession, exc)
        return False

    paired = [
        result
        for path, result in zip(fastq_files, measured)
        if re.search(r"_[12](?:\.fastq|\.fq)(?:\.gz)?$", path.name)
    ]
    num_written = max((records for records, _ in paired), default=measured[0][0])
    bp_written = sum(bases for _, bases in measured)
    payload = "run\tnum_written\tbp_written\tbp_fastp_in\n" f"{accession}\t{num_written}\t{bp_written}\t{bp_written}\n"
    temporary = stats_path.with_name(f".{stats_path.name}.{os.getpid()}.tmp")
    try:
        sample_dir.mkdir(parents=True, exist_ok=True)
        _write_text_durably(stats_path, payload, temporary)
    except OSError as exc:
        logger.warning("Unable to persist FASTQ statistics for %s: %s", accession, exc)
        return False
    logger.info(
        "Derived quant-only FASTQ statistics for %s: %d spots, %d bases",
        accession,
        num_written,
        bp_written,
    )
    return True


def _is_valid_run_accession(value: Any) -> bool:
    """Return whether *value* is a concrete SRA/ENA run accession.

    Amalgkit metadata can contain experiment or SRA-primary records without a
    public read-run accession.  Those rows cannot be downloaded or quantified
    by the streaming runner and must not be allowed to become filesystem path
    components (pandas represents blank TSV fields as floating-point NaN).
    """

    return isinstance(value, str) and bool(_RUN_ACCESSION_PATTERN.fullmatch(value.strip()))


def _raw_file_priority(name: str) -> int | None:
    """Return the priority represented by a raw-input filename."""

    if name.endswith(_RAW_INPUT_SUFFIXES):
        return 0
    if name.endswith(_RAW_PARTIAL_SUFFIXES):
        return 1
    return None


def _is_diagnostic_transfer_artifact(path: Path) -> bool:
    """Return whether a failed-transfer artifact should be retained."""

    name = path.name
    return name.endswith(".part") or ".invalid" in name


def _raw_input_state_index(fastq_dir: Path, accessions: set[str]) -> dict[str, tuple[int, int]]:
    """Index raw state and payload size with bounded directory scans.

    The returned tuple is ``(priority, payload_bytes)``: valid FASTQs use
    priority 0, resumable partials priority 1, and fresh acquisitions priority
    2. The size is used only to improve queue order; it never changes the
    provenance or completion contract.
    """

    states = {accession: (2, 0) for accession in accessions}

    def record(path: Path, accession: str) -> None:
        try:
            file_stat = path.stat()
            if not stat.S_ISREG(file_stat.st_mode) or file_stat.st_size <= 0:
                return
            size = file_stat.st_size
        except OSError:
            return
        priority = _raw_file_priority(path.name)
        if priority is None:
            return
        current_priority, current_bytes = states[accession]
        if priority < current_priority:
            states[accession] = (priority, size)
        elif priority == current_priority:
            states[accession] = (current_priority, current_bytes + size)

    try:
        for accession in accessions:
            sample_dir = fastq_dir / accession
            if sample_dir.is_dir():
                for path in sample_dir.iterdir():
                    record(path, accession)

        # Some historical Amalgkit trees keep files directly under
        # ``getfastq``. Read this directory once and recover the accession
        # from the ENA filename instead of rescanning it for every task.
        if fastq_dir.is_dir():
            for path in fastq_dir.iterdir():
                match = re.match(r"^(SRR|ERR|DRR)\d+", path.name)
                if match and match.group(0) in states:
                    record(path, match.group(0))

        # A campaign-local NCBI cache is shared across species and may hold
        # complete SRA databases after an earlier fallback or a separate
        # acquisition pass.  Index exact accession paths without recursively
        # scanning the potentially very large cache tree.
        for accession in accessions:
            cached_sra = _campaign_cached_sra_path(accession)
            if cached_sra is not None:
                record(cached_sra, accession)
    except OSError:
        # Discovery must not abort a campaign because one pre-contract mount
        # entry became unreadable. The normal downloader will report failure.
        return states
    return states


def _task_has_local_sra(task: Dict[str, Any]) -> bool:
    """Return whether a task has an SRA source that uses the fallback lane."""

    accession = str(task["srr"]).strip()
    fastq_dir = Path(task["fastq_dir"])
    sample_dir = fastq_dir / accession
    expected_paired = task.get("expected_paired")
    unnumbered = _raw_fastq_path(sample_dir, accession)
    mate_1 = _raw_fastq_path(sample_dir, accession, 1)
    mate_2 = _raw_fastq_path(sample_dir, accession, 2)
    if expected_paired is True:
        has_fastq = (mate_1 is not None and mate_2 is not None) or (
            mate_1 is None and mate_2 is None and unnumbered is not None
        )
    elif expected_paired is False:
        has_fastq = unnumbered is not None
    else:
        has_fastq = unnumbered is not None or (mate_1 is not None and mate_2 is not None)
    return not has_fastq and _local_sra_path(sample_dir, fastq_dir, accession) is not None


# Median gzip-FASTQ bytes per sequenced base across the 18,167 frozen Hymenoptera tasks with both values (0.498):
# the acquisition weight of a task whose metadata has total_bases but no ENA fastq_bytes.
FASTQ_BYTES_PER_BASE_NUMERATOR = 1
FASTQ_BYTES_PER_BASE_DENOMINATOR = 2


def task_identity(task: Mapping[str, Any]) -> str:
    """The scheduling identity ``species/run`` (also the deterministic tie-break)."""

    return f"{task.get('species_name') or task.get('species') or ''}/{str(task['srr']).strip()}"


def task_acquisition_weight(task: Mapping[str, Any]) -> int:
    """Estimated compressed FASTQ bytes a fresh acquisition writes (0 when unknown).

    Uses the frozen cohort's ENA ``fastq_bytes`` when present, else
    ``total_bases`` x 0.5 (see :data:`FASTQ_BYTES_PER_BASE_NUMERATOR`).
    """

    for key, scale in (("fastq_bytes", None), ("total_bases", True)):
        try:
            value = int(task.get(key) or 0)
        except (TypeError, ValueError):
            value = 0
        if value > 0:
            if scale:
                return value * FASTQ_BYTES_PER_BASE_NUMERATOR // FASTQ_BYTES_PER_BASE_DENOMINATOR
            return value
    return 0


def _annotate_raw_input_state(tasks: List[Dict[str, Any]]) -> None:
    """Record each task's raw-input tier and payload from bounded directory scans."""

    accessions_by_dir: dict[Path, set[str]] = {}
    for task in tasks:
        fastq_dir = Path(task["fastq_dir"])
        accessions_by_dir.setdefault(fastq_dir, set()).add(str(task["srr"]).strip())
    states_by_dir = {
        fastq_dir: _raw_input_state_index(fastq_dir, accessions) for fastq_dir, accessions in accessions_by_dir.items()
    }
    for task in tasks:
        # Cache the filesystem observation on the transient task record.  A
        # campaign can contain tens of thousands of rows; the index above
        # keeps this lookup O(1) per task after bounded directory scans.
        fastq_dir = Path(task["fastq_dir"])
        accession = str(task["srr"]).strip()
        priority, payload_bytes = states_by_dir[fastq_dir].get(accession, (2, 0))
        task["_raw_input_priority"] = priority
        task["_raw_input_bytes"] = payload_bytes
        task["_uses_local_sra"] = _task_has_local_sra(task)


def _fresh_acquisition_order(tasks: List[Dict[str, Any]], order: str) -> List[Dict[str, Any]]:
    """Order tier-2 (fresh) tasks with the validated acquisition planner.

    ``size_bands`` (the c19 default) visits global size bands smallest first
    and round-robins species inside a band, so every species advances;
    ``smallest_first`` is one global size queue; ``species_order`` keeps the
    discovery order. Ties break on ``species/run``; unknown sizes go last.
    """

    if order == "species_order" or len(tasks) < 2:
        return list(tasks)
    by_identity: dict[str, Dict[str, Any]] = {}
    duplicates: List[Dict[str, Any]] = []
    entries = []
    for task in tasks:
        identity = task_identity(task)
        if identity in by_identity:
            duplicates.append(task)
            continue
        by_identity[identity] = task
        entries.append(
            {
                "task_id": identity,
                "species": task.get("species_name") or task.get("species") or "",
                "fastq_bytes": task_acquisition_weight(task),
            }
        )
    planned = plan_acquisition_order(
        entries,
        (),
        global_size_bands=order == "size_bands",
        global_size_order=order == "smallest_first",
    )
    return [by_identity[identity] for identity in planned] + duplicates


def order_campaign_tasks(
    tasks: List[Dict[str, Any]],
    *,
    workers: int | None = None,
    fasterq_slots: int | None = None,
    order: str = "size_bands",
) -> List[Dict[str, Any]]:
    """Order annotated tasks: reusable inputs, then partials, then fresh acquisitions.

    Tier 0 (valid local FASTQs) runs smallest first; tier 1 (resumable
    partials) runs most complete first; tier 2 (fresh) follows ``order``.
    Local SRA databases are admitted in bounded fallback-sized chunks,
    interleaved with partial and new acquisitions, so the extraction
    semaphore does not leave every sample worker waiting on a mounted-disk
    conversion while ENA work could continue independently. Pure: reads only
    the ``_raw_input_*`` annotations, so a campaign can re-order its pending
    queue whenever a species' discovery completes.
    """

    if order not in SCHEDULING_ORDERS:
        raise ValueError(f"scheduling order must be one of {SCHEDULING_ORDERS}, not {order!r}")
    reusable = [task for task in tasks if task.get("_raw_input_priority", 2) == 0]
    partial = [task for task in tasks if task.get("_raw_input_priority", 2) == 1]
    fresh = [task for task in tasks if task.get("_raw_input_priority", 2) not in (0, 1)]

    def tie_break(task: Mapping[str, Any]) -> tuple[str, int, str]:
        # Independent of which species' discovery finished first: species, its
        # discovery position, then the run. ``species_order`` keeps the
        # historical stable input order (the sort is stable).
        if order == "species_order":
            return "", 0, ""
        species = str(task.get("species_name") or task.get("species") or "")
        try:
            position = int(task.get("batch_idx") or 0)
        except (TypeError, ValueError):
            position = 0
        return species, position, str(task["srr"]).strip()

    # Shortest reusable reads first improves time-to-first-results and
    # prevents a few multi-gigabyte libraries from monopolizing every
    # quant slot while smaller validated samples wait.
    reusable.sort(key=lambda task: (task.get("_raw_input_bytes", 0), tie_break(task)))
    # A larger partial has more acquired data and is usually closer to
    # becoming reusable; finish the most complete resumptions first.
    partial.sort(key=lambda task: (-task.get("_raw_input_bytes", 0), tie_break(task)))
    ordered = reusable + partial + _fresh_acquisition_order(fresh, order)
    if not ordered:
        return ordered

    # A small initial FASTQ-reuse tier should drain immediately. After that,
    # interleave at most one fallback-sized cache chunk with the non-cache
    # queue. This preserves the preference for downloaded material while
    # keeping network/partial-transfer workers useful during long extraction.
    reusable_fastq = [task for task in reusable if not task.get("_uses_local_sra")]
    cached_sra = [task for task in reusable if task.get("_uses_local_sra")]
    non_cached = ordered[len(reusable) :]
    if not cached_sra or not non_cached:
        return ordered

    worker_budget = max(1, workers or len(ordered))
    fallback_budget = max(1, min(fasterq_slots or 1, worker_budget))
    non_cached_chunk = max(1, worker_budget - fallback_budget)
    result = list(reusable_fastq)
    cache_index = 0
    non_cached_index = 0
    while cache_index < len(cached_sra) or non_cached_index < len(non_cached):
        result.extend(cached_sra[cache_index : cache_index + fallback_budget])
        cache_index += fallback_budget
        result.extend(non_cached[non_cached_index : non_cached_index + non_cached_chunk])
        non_cached_index += non_cached_chunk
    return result


def _prioritize_tasks_by_raw_state(
    tasks: List[Dict[str, Any]],
    *,
    workers: int | None = None,
    fasterq_slots: int | None = None,
    order: str = "size_bands",
) -> List[Dict[str, Any]]:
    """Annotate raw-input state, then order with :func:`order_campaign_tasks`."""

    _annotate_raw_input_state(tasks)
    return order_campaign_tasks(tasks, workers=workers, fasterq_slots=fasterq_slots, order=order)


def _filter_valid_run_rows(df: Any, species_name: str, metadata_path: Path) -> Any:
    """Keep downloadable run rows and write an auditable exclusion ledger.

    Rows without a concrete ``run``/``run_accession`` are retained in the
    source metadata provenance only through the ledger; they are excluded from
    the streaming metadata used for batch quantification.  This makes the
    unavailable-data boundary explicit and prevents a missing accession from
    aborting discovery for an entire species.
    """

    run_column = _sample_run_column(df)
    valid_mask = df[run_column].map(_is_valid_run_accession)
    invalid = df.loc[~valid_mask].copy()
    if invalid.empty:
        return df.loc[valid_mask].copy()

    ledger = metadata_path.parent / "validation" / "invalid_run_accessions.tsv"
    ledger.parent.mkdir(parents=True, exist_ok=True)
    invalid.insert(0, "invalid_run_reason", "missing_or_invalid_run_accession")
    invalid.to_csv(ledger, sep="\t", index=False)
    logger.warning(
        "%s: excluding %d metadata rows without a concrete %s accession; " "see %s",
        species_name,
        len(invalid),
        run_column,
        ledger,
    )
    return df.loc[valid_mask].copy()


def _resolve_quant_metadata_path(work_dir: Path) -> Optional[str]:
    """Find the metadata file consumed by amalgkit quant."""
    for metadata_path in [work_dir / "metadata/metadata_selected.tsv", work_dir / "metadata/metadata.tsv"]:
        if metadata_path.exists():
            return str(metadata_path)
    return None


def _resolve_data_config_path(value: Any) -> str:
    """Resolve a config path under the selected external data root."""

    if not isinstance(value, str) or not value:
        return ""
    expanded = os.path.expanduser(os.path.expandvars(value))
    relative = expanded[2:] if expanded.startswith("./") else expanded
    prefix = "output/amalgkit/"
    if relative == "output/amalgkit":
        return str(_data_root())
    if relative.startswith(prefix):
        return str(_data_root() / relative[len(prefix) :])
    return expanded


def _quant_index_candidates(cfg: Dict[str, Any], species_name: str) -> List[str]:
    """Build candidate Kallisto index directories from the current config."""
    quant_index_dir = _resolve_data_config_path(cfg.get("steps", {}).get("quant", {}).get("index_dir", ""))
    genome_dest = _resolve_data_config_path(cfg.get("genome", {}).get("dest_dir", ""))
    genome_index_dir = f"{genome_dest}/index" if genome_dest else ""
    explicit_index_dir = _resolve_data_config_path(cfg.get("genome", {}).get("index_dir", ""))

    candidates = [
        quant_index_dir,
        genome_index_dir,
        explicit_index_dir,
        str(_species_work_dir(species_name) / "index"),
        str(_data_root() / "shared" / "genome" / species_name / "index"),
    ]
    for species in cfg.get("species_list", []):
        candidates.append(str(_data_root() / "shared" / "genome" / species / "index"))
    return [str(Path(candidate).expanduser()) for candidate in candidates if candidate]


def _resolve_index_dir(cfg: Dict[str, Any], species_name: str) -> str:
    """Resolve the first existing index directory containing index files."""
    explicit_index_dir = _resolve_data_config_path(cfg.get("genome", {}).get("index_dir", ""))
    explicit_path = Path(explicit_index_dir).expanduser() if explicit_index_dir else Path()
    if explicit_index_dir and explicit_path.exists():
        return str(explicit_path)

    for candidate in _quant_index_candidates(cfg, species_name):
        if candidate and Path(candidate).exists():
            usable_indexes = [
                index_path
                for index_path in Path(candidate).glob("*.idx")
                if index_path.is_file() and index_path.stat().st_size > 0
            ]
            if usable_indexes:
                return candidate
    return ""


def _normalized_reference_stem(value: Any) -> str:
    """Normalize a metadata target name to Amalgkit's index-file stem."""

    return re.sub(r"[^A-Za-z0-9]+", "_", str(value).strip()).strip("_")


def _find_reference_index(index_dir: Path, stem: str) -> Optional[Path]:
    """Find a non-empty current or compatibility index for one stem."""

    for candidate in (index_dir / f"{stem}.idx", index_dir / f"{stem}_transcripts.idx"):
        try:
            if candidate.is_file() and candidate.stat().st_size > 0:
                return candidate
        except OSError:
            continue
    return None


KALLISTO_DEFAULT_KMER_SIZE = 31


@lru_cache(maxsize=128)
def _kallisto_index_kmer(index_path: str) -> int:
    """Read the k-mer size from a Kallisto index, failing closed to 31.

    Kallisto cannot pseudoalign reads shorter than the index k-mer.  The
    executable project indexes currently use 31-mers, but reading the index
    keeps the eligibility gate correct if a reference is rebuilt with another
    supported k-mer size.  A missing inspection tool must not make discovery
    non-deterministic, so the documented Kallisto default is used as a strict
    fallback and recorded in the eligibility manifest.
    """

    try:
        result = subprocess.run(
            ["kallisto", "inspect", str(index_path)],
            capture_output=True,
            text=True,
            timeout=120,
        )
        combined = f"{result.stdout}\n{result.stderr}"
        match = re.search(r"k-mer length:\s*(\d+)", combined, flags=re.IGNORECASE)
        if match:
            return int(match.group(1))
    except (OSError, subprocess.SubprocessError, ValueError) as exc:
        logger.warning("Unable to inspect Kallisto index %s: %s", index_path, exc)
    logger.warning(
        "Using the documented Kallisto k-mer fallback (%d) for %s",
        KALLISTO_DEFAULT_KMER_SIZE,
        index_path,
    )
    return KALLISTO_DEFAULT_KMER_SIZE


def _reference_index_for_kmer(index_dir: Optional[Path], target_names: list[str]) -> Optional[Path]:
    """Choose one readable index representing the active species reference."""

    if index_dir is None:
        return None
    for target in target_names:
        candidate = _find_reference_index(index_dir, _normalized_reference_stem(target))
        if candidate is not None:
            return candidate
    try:
        for candidate in sorted(index_dir.glob("*.idx")):
            if candidate.is_file() and candidate.stat().st_size > 0:
                return candidate
    except OSError:
        return None
    return None


def _filter_kallisto_ineligible_reads(dataframe: Any, kmer_size: int) -> tuple[Any, Any, Any]:
    """Exclude only known reads shorter than the active Kallisto k-mer.

    ``spot_length`` is preferred, with ``total_bases / total_spots`` as the
    documented Amalgkit fallback when the public metadata row omits it.  Rows
    with no usable length estimate are retained for Amalgkit to resolve from
    the downloaded input; silently treating unknown lengths as eligible would
    be less auditable than recording the unresolved value.
    """

    import pandas as pd

    read_length = pd.Series(float("nan"), index=dataframe.index, dtype="float64")
    if "spot_length" in dataframe.columns:
        read_length = pd.to_numeric(dataframe["spot_length"], errors="coerce")
    if {"total_bases", "total_spots"}.issubset(dataframe.columns):
        total_bases = pd.to_numeric(dataframe["total_bases"], errors="coerce")
        total_spots = pd.to_numeric(dataframe["total_spots"], errors="coerce")
        ratio = total_bases.div(total_spots.replace(0, pd.NA))
        read_length = read_length.where(read_length.gt(0), ratio)

    ineligible = read_length.notna() & read_length.lt(int(kmer_size))
    eligible = dataframe.loc[~ineligible].copy()
    excluded = dataframe.loc[ineligible].copy()
    excluded["quant_eligibility_mean_read_length"] = read_length.loc[ineligible].astype(float)
    excluded["quant_eligibility_kmer_size"] = int(kmer_size)
    excluded["quant_eligibility_reason"] = "mean_read_length_below_kallisto_kmer"
    return eligible, excluded, read_length


def _write_quant_eligibility_audit(
    work_dir: Path,
    excluded: Any,
    *,
    kmer_size: int,
    index_path: Optional[Path],
) -> Optional[Path]:
    """Write the fail-closed short-read exclusion table atomically."""

    if excluded is None or excluded.empty:
        return None
    destination = work_dir / "metadata" / "quant_eligibility_exclusions.tsv"
    destination.parent.mkdir(parents=True, exist_ok=True)
    output = excluded.copy()
    output["quant_eligibility_reference_index"] = str(index_path) if index_path else ""
    output["quant_eligibility_kmer_size"] = int(kmer_size)
    temporary = destination.with_name(f".{destination.name}.tmp")
    try:
        output.to_csv(temporary, sep="\t", index=False)
        os.replace(temporary, destination)
    except OSError:
        temporary.unlink(missing_ok=True)
        raise
    return destination


def _write_reference_alias_manifest(
    species_name: str,
    work_dir: Path,
    index_dir: Optional[Path],
    targets: list[str],
    aliases: list[dict[str, Any]],
    missing: list[str],
    *,
    kallisto_index: Optional[Path] = None,
    kallisto_kmer_size: Optional[int] = None,
) -> Path:
    """Persist the reference-target audit used by quantification provenance."""

    destination = work_dir / "reference" / "reference_aliases.json"
    destination.parent.mkdir(parents=True, exist_ok=True)
    # Content-deterministic payload: no wall-clock fields.  This manifest's
    # SHA-256 is recorded in per-sample quantification provenance sidecars and
    # re-verified on every resume, so any restart-varying byte (a timestamp)
    # would invalidate every previously quantified sample of the species and
    # quarantine it for redundant re-quantification (observed 2026-08-30:
    # 1,620 apis_mellifera samples re-queued by one restart).  Recency lives
    # in the orchestrator log and the progress DB, not here.
    payload = {
        "schema": "metainformant.rna.reference_aliases.v1",
        "species": species_name,
        "index_dir": str(index_dir) if index_dir else None,
        "targets": targets,
        "aliases": aliases,
        "missing": missing,
        "status": "complete" if not missing else "incomplete",
        "kallisto_index": str(kallisto_index) if kallisto_index else None,
        "kallisto_kmer_size": kallisto_kmer_size,
    }
    text = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    try:
        if destination.read_bytes() == text.encode("utf-8"):
            # Identical bytes: never touch the file. Outputs bound to it through their sidecar's
            # reference_manifest_sha256 stay valid, and the file's identity (inode, mtime) is unchanged.
            return destination
    except FileNotFoundError:
        pass
    else:
        logger.warning(
            "%s: reference manifest %s changes; outputs bound only to its old bytes (no binding witness) will "
            "no longer verify through it",
            species_name,
            destination,
        )
    temporary = destination.with_name(f".{destination.name}.{os.getpid()}.{threading.get_ident()}.tmp")
    _write_text_durably(destination, text, temporary)
    return destination


def _ensure_reference_alias_indexes(
    cfg: Dict[str, Any], species_name: str, target_names: list[str]
) -> tuple[bool, list[str], Optional[Path]]:
    """Ensure every metadata target resolves to a current index.

    Amalgkit resolves an index by metadata target name.  Public metadata can
    contain subspecies or breed labels even when a project intentionally uses
    a species-level reference.  Such mappings must be explicit: this helper
    creates lightweight symlink aliases only for mappings declared in the
    configuration and records them under the external work tree.  Unmapped
    targets fail closed before any reads are downloaded.
    """

    work_dir = _species_work_dir(species_name)
    index_value = _resolve_index_dir(cfg, species_name)
    index_dir = Path(index_value) if index_value else None
    configured = cfg.get("reference_aliases", {})
    aliases_config = configured if isinstance(configured, dict) else {}
    aliases_by_target = {
        str(key).strip().casefold(): str(value).strip()
        for key, value in aliases_config.items()
        if str(key).strip() and str(value).strip()
    }
    targets = sorted({str(target).strip() for target in target_names if str(target).strip()})
    aliases: list[dict[str, Any]] = []
    missing: list[str] = []

    if index_dir is None:
        missing = [f"{target} (no usable index directory)" for target in targets]
    else:
        index_dir.mkdir(parents=True, exist_ok=True)
        for target in targets:
            target_stem = _normalized_reference_stem(target)
            target_index = _find_reference_index(index_dir, target_stem)
            if target_index is not None and not (target_index.is_symlink() and target.casefold() in aliases_by_target):
                aliases.append(
                    {
                        "metadata_target": target,
                        "target_stem": target_stem,
                        "reference_stem": target_stem,
                        "index": str(target_index),
                        "mode": "native_index",
                    }
                )
                continue

            source_value = aliases_by_target.get(target.casefold())
            if source_value is None:
                missing.append(target)
                continue
            source_stem = _normalized_reference_stem(source_value)
            source_index = _find_reference_index(index_dir, source_stem)
            if source_index is None:
                missing.append(f"{target} (alias source unavailable: {source_value})")
                continue

            alias_paths = [
                index_dir / f"{target_stem}.idx",
                index_dir / f"{target_stem}_transcripts.idx",
            ]
            conflict = False
            for alias_path in alias_paths:
                try:
                    if alias_path.is_symlink() and alias_path.resolve(strict=False) == source_index.resolve():
                        continue
                    if alias_path.exists():
                        # Never overwrite an existing index with a potentially
                        # different reference.  Require explicit cleanup or a
                        # corrected config instead of silently changing reads.
                        conflict = True
                        continue
                    alias_path.symlink_to(os.path.relpath(source_index, alias_path.parent))
                except OSError as exc:
                    missing.append(f"{target} (unable to create alias: {exc})")
                    conflict = True
            if conflict:
                continue
            aliases.append(
                {
                    "metadata_target": target,
                    "target_stem": target_stem,
                    "reference_stem": source_stem,
                    "source_index": str(source_index),
                    "index_aliases": [str(path) for path in alias_paths],
                    "mode": "explicit_species_level_alias",
                    "note": (
                        "Metadata target lacks an independent configured index; "
                        "reads use the declared reference source."
                    ),
                }
            )

    reference_index = _reference_index_for_kmer(index_dir, targets)
    reference_kmer = _kallisto_index_kmer(str(reference_index)) if reference_index else None
    manifest = _write_reference_alias_manifest(
        species_name,
        work_dir,
        index_dir,
        targets,
        aliases,
        missing,
        kallisto_index=reference_index,
        kallisto_kmer_size=reference_kmer,
    )
    if missing:
        logger.error(
            "%s: reference preflight failed for metadata targets: %s; see %s",
            species_name,
            "; ".join(missing),
            manifest,
        )
        return False, missing, manifest
    if aliases:
        explicit = [entry for entry in aliases if entry.get("mode") == "explicit_species_level_alias"]
        if explicit:
            logger.warning(
                "%s: using %d explicitly declared species-level reference aliases; see %s",
                species_name,
                len(explicit),
                manifest,
            )
    return True, missing, manifest


def _build_quant_command(
    cfg: Dict[str, Any],
    species_name: str,
    batch_index: int,
    threads: int,
    metadata_path: str,
    reference_target: str | None = None,
) -> List[str]:
    """Build the amalgkit quant command for one streaming batch."""
    cmd = [
        "amalgkit",
        "quant",
        "--out_dir",
        str(_species_work_dir(species_name)),
        "--metadata",
        metadata_path,
        "--threads",
        str(threads),
        "--batch",
        str(batch_index),
        # A task is scheduled only when no current-method provenance sidecar
        # exists. Replace only explicitly selected re-quantification candidates
        # before writing the new sidecar.
        "--redo",
        "yes",
    ]

    # Keep raw-input deletion behind the provenance gate in
    # ``process_single_sample``.  Amalgkit's own clean-up can run before the
    # current-method sidecar is written, which would make a failed sidecar
    # write indistinguishable from a safely reclaimed sample.  The explicit
    # cleanup helper remains idempotent and removes both reads and redundant
    # ``.safely_removed`` markers after current quantification is proven.
    cmd.extend(["--clean_fastq", "no"])

    # The index_dir handed to amalgkit must contain only non-colliding
    # per-species indexes, so the resolved directory (optionally cached) is
    # re-staged with the metadata target, including explicit subspecies aliases.
    index_dir = _staged_exact_index_dir(
        _cache_index_directory(_resolve_index_dir(cfg, species_name), species_name), reference_target or species_name
    )
    if index_dir:
        cmd.extend(["--index_dir", index_dir])

    return cmd


def _copy_index_verified(source_file: Path, destination: Path) -> dict[str, Any]:
    """Copy an index through a unique, fsynced temp file, verify its sha256 against the source, then promote it.

    Returns the copy's record for the stage manifest: the source's identity
    (device, inode, size, mtime, ctime), the sha256 of both, and the copy's
    own ``(inode, size, mtime_ns, ctime_ns)`` so a later rewrite, replacement
    or truncation of the copy is noticed without re-hashing it (``copy2``
    gives the copy the source's mtime, and ``utime`` can restore an mtime,
    but no write leaves the ctime unchanged).
    """

    source_signature = file_signature(source_file)
    descriptor, temporary_name = tempfile.mkstemp(prefix=f".{destination.name}.", suffix=".tmp", dir=destination.parent)
    os.close(descriptor)
    temporary = Path(temporary_name)
    try:
        shutil.copy2(source_file, temporary)
        with open(temporary, "rb") as handle:
            _durable_fsync(handle.fileno())
        source_sha256 = sha256_file(source_file)
        copy_sha256 = sha256_file(temporary)
        if copy_sha256 != source_sha256 or file_signature(source_file) != source_signature:
            raise OSError(f"index copy of {source_file} does not match its source bytes")
        os.replace(temporary, destination)
    except BaseException:
        temporary.unlink(missing_ok=True)
        raise
    _fsync_directory(destination.parent)
    info = destination.lstat()
    return {
        "source": str(source_file),
        "source_signature": list(source_signature),
        "sha256": copy_sha256,
        "copy": [info.st_ino, info.st_size, info.st_mtime_ns, info.st_ctime_ns],
    }


def _index_copies_current(target: Path, sources: Sequence[Path], records: Any) -> bool:
    """True when every copy in ``target`` is the one recorded for its unchanged source."""

    if not isinstance(records, dict):
        return False
    try:
        for source_file in sources:
            record = records.get(source_file.name)
            if not isinstance(record, dict) or record.get("source_signature") != list(file_signature(source_file)):
                return False
            info = (target / source_file.name).lstat()
            if record.get("copy") != [info.st_ino, info.st_size, info.st_mtime_ns, info.st_ctime_ns]:
                return False
    except OSError:
        return False
    return True


def _write_index_records(manifest_path: Path, payload: Mapping[str, Any]) -> None:
    temporary = manifest_path.with_name(f".{manifest_path.name}.{os.getpid()}.{threading.get_ident()}.tmp")
    _write_text_durably(manifest_path, json.dumps(payload, indent=2, sort_keys=True) + "\n", temporary)


def _cache_index_directory(index_dir: str, species_name: str) -> str:
    """Return a validated local copy of a reference index when configured.

    Kallisto repeatedly memory-maps the reference index while the FASTQ files
    remain on the external data volume.  ``AMALGKIT_LOCAL_INDEX_CACHE_DIR`` is
    an opt-in performance lane for hosts where that volume has poor random-read
    latency.  Only ``*.idx`` files are copied, each copy is promoted atomically,
    and a size/mtime manifest invalidates the cache when the source changes.
    With no cache variable, the original path is returned unchanged.
    """

    if not index_dir:
        return index_dir
    cache_value = os.environ.get("AMALGKIT_LOCAL_INDEX_CACHE_DIR", "").strip()
    if not cache_value:
        return index_dir

    source = Path(index_dir).expanduser()
    if not source.is_dir():
        return index_dir
    try:
        source_files = sorted(path for path in source.glob("*.idx") if path.is_file() and path.stat().st_size > 0)
    except OSError as exc:
        logger.warning("Unable to inspect reference index directory %s: %s", source, exc)
        return index_dir
    if not source_files:
        return index_dir

    target = Path(cache_value).expanduser() / species_name / "index"
    manifest_path = target / ".metainformant_index_cache.json"
    expected = {path.name: {"size": path.stat().st_size, "mtime_ns": path.stat().st_mtime_ns} for path in source_files}

    with _INDEX_CACHE_LOCK:
        try:
            manifest = json.loads(manifest_path.read_text()) if manifest_path.is_file() else {}
        except (OSError, json.JSONDecodeError, TypeError):
            manifest = {}
        cached_files = manifest.get("files") if isinstance(manifest, dict) else None
        if (
            isinstance(cached_files, dict)
            and manifest.get("source") == str(source)
            and cached_files == expected
            and _index_copies_current(target, source_files, manifest.get("copies"))
        ):
            return str(target)

        try:
            target.mkdir(parents=True, exist_ok=True)
            copies = {
                source_file.name: _copy_index_verified(source_file, target / source_file.name)
                for source_file in source_files
            }
            payload = {
                "schema": "metainformant.rna.index_cache.v2",
                "source": str(source),
                "files": expected,
                "copies": copies,
            }
            _write_index_records(manifest_path, payload)
            logger.info("Using local reference-index cache for %s: %s", species_name, target)
            return str(target)
        except OSError as exc:
            logger.warning("Unable to cache reference indexes for %s; using %s: %s", species_name, source, exc)
            return index_dir


def _staged_exact_index_dir(index_dir: str, species_name: str) -> str:
    """Return a directory holding only the non-colliding index for one species.

    A species such as ``bombus_terrestris`` may share a reference volume with
    its subspecies, so one directory can contain both ``Bombus_terrestris.idx``
    and ``Bombus_terrestris_audax.idx``.  A file ``X.idx`` prefix-collides with
    ``stem`` when ``X != stem`` and the two names share an underscore-delimited
    prefix (``X.startswith(stem + "_")`` or ``stem.startswith(X + "_")``), so
    either a parent or a subspecies index makes the directory ambiguous.  When
    the exact ``{stem}.idx`` file exists, it (plus a matching
    ``{stem}_transcripts.idx``) is staged alone into
    ``<index_dir>_exact_stem/<stem>/index/`` so Amalgkit and Kallisto can never
    pick the wrong index by prefix.  Each file is copied to a unique temp
    name, fsynced, verified by sha256 against its source and promoted with
    ``os.replace``; a stage manifest records the source identity and the
    copy's identity, and a copy is reused only while both are unchanged (a
    rebuilt source or a damaged copy is staged again).  When no
    usable index can be isolated, the original directory is returned unchanged
    and the failure is left to worker triage instead of raising here.
    """

    if not index_dir:
        return index_dir
    stem = _normalized_reference_stem(species_name)
    if not stem:
        return index_dir
    source = Path(index_dir).expanduser()
    if not source.is_dir():
        return index_dir
    try:
        index_files = sorted(path for path in source.glob("*.idx") if path.is_file() and path.stat().st_size > 0)
    except OSError as exc:
        logger.warning("Unable to inspect reference index directory %s: %s", source, exc)
        return index_dir

    # Species keys arrive lowercase (e.g. ``nasonia_vitripennis``) while staged
    # index files keep the reference case (``Nasonia_vitripennis.idx``); match
    # stems case-insensitively or every species silently degrades here.
    folded_names = {path.stem.casefold(): path for path in index_files}
    folded_stem = stem.casefold()
    exact_file = folded_names.get(folded_stem)
    colliding = sorted(
        name
        for name in folded_names
        if name != folded_stem and (name.startswith(folded_stem + "_") or folded_stem.startswith(name + "_"))
    )
    if exact_file is not None and not colliding:
        return index_dir

    selected = [exact_file] if exact_file is not None else None
    if selected is None:
        forward_matches = sorted(
            name for name in folded_names if name != folded_stem and name.startswith(folded_stem + "_")
        )
        if len(forward_matches) == 1:
            # No exact-stem file: a single forward-prefix match is unambiguous.
            selected = [folded_names[forward_matches[0]]]
        else:
            logger.warning(
                "index_unavailable: found %d candidate Kallisto index files for species %s in %s",
                len(forward_matches),
                species_name,
                source,
            )
            return index_dir

    transcripts_file = folded_names.get(f"{folded_stem}_transcripts")
    if transcripts_file is not None:
        selected.append(transcripts_file)

    target = Path(f"{source}_exact_stem") / stem / "index"
    manifest_path = target / INDEX_STAGE_MANIFEST
    with _INDEX_STAGE_LOCK:
        try:
            try:
                records = json.loads(manifest_path.read_bytes()).get("copies")
            except (OSError, ValueError, AttributeError):
                records = None
            if _index_copies_current(target, selected, records):
                return str(target)
            # Missing, stale (the source changed) or tampered (the copy changed): copy again, verified.
            target.mkdir(parents=True, exist_ok=True)
            copies = {
                source_file.name: _copy_index_verified(source_file, target / source_file.name)
                for source_file in selected
            }
            _write_index_records(
                manifest_path, {"schema": "metainformant.rna.index_stage.v1", "source": str(source), "copies": copies}
            )
        except OSError as exc:
            logger.warning("Unable to stage exact-stem indexes for %s; using %s: %s", species_name, source, exc)
            return index_dir
    logger.info("Staged exact-stem Kallisto index for %s: %s", species_name, target)
    return str(target)


def _build_sample_tasks(
    filtered: Any, srr_col: str, fastq_dir: Path, config_path: Path, species_name: str
) -> List[Dict[str, Any]]:
    """Build task dictionaries from filtered metadata rows."""
    tasks = []
    for i, (_, row) in enumerate(filtered.iterrows()):
        library_layout = row.get("lib_layout")
        if not isinstance(library_layout, str):
            library_layout = row.get("library_layout")
        try:
            total_bases = max(0, int(float(row.get("total_bases", 0) or 0)))
        except (TypeError, ValueError):
            total_bases = 0
        tasks.append(
            {
                "srr": row[srr_col],
                "batch_idx": i + 1,
                "fastq_dir": fastq_dir,
                "config_path": config_path,
                "species_name": species_name,
                "expected_paired": _library_layout_is_paired(library_layout),
                # Carried for global smallest-first scheduling across species.
                "total_bases": total_bases,
            }
        )
    return tasks


class StreamingPipelineOrchestrator:
    """Orchestrator for streaming ENA download -> Quant processing."""

    def __init__(
        self,
        config_dir: Path = CONFIG_DIR,
        log_dir: Path = LOG_DIR,
        db_path: Optional[Path] = None,
        *,
        cohort: FrozenCohort | None = None,
        allow_metadata_refresh: bool = True,
    ):
        """Create the orchestrator.

        With ``cohort``, discovery builds tasks only from the frozen cohort (no
        Amalgkit metadata/select, no rewritten selection table, no pruning).
        Without it, ``allow_metadata_refresh=False`` refuses to query NCBI
        (Amalgkit metadata/select) when a species' metadata is missing or stale.
        """
        from metainformant.rna.engine.ehd_storage import capture_ehd_write_policy

        self.cohort = cohort
        self.allow_metadata_refresh = allow_metadata_refresh
        self._cohort_digests = IndexDigestCache()

        # Validate every managed write root before logs or SQLite are opened.
        self._ehd_write_policy = capture_ehd_write_policy(Path(config_dir), Path(log_dir), db_path)
        self.acquisition_slots = int(os.environ.get("AMALGKIT_PIPELINE_ACQUISITION_SLOTS", "0"))
        self.acquisition_queue_timeout_seconds = float(
            os.environ.get("AMALGKIT_PIPELINE_ACQUISITION_QUEUE_TIMEOUT_SECONDS", "7200")
        )
        if not 0 <= self.acquisition_slots <= 256:
            raise ValueError("AMALGKIT_PIPELINE_ACQUISITION_SLOTS must be between 0 and 256")
        if not math.isfinite(self.acquisition_queue_timeout_seconds) or self.acquisition_queue_timeout_seconds <= 0:
            raise ValueError("AMALGKIT_PIPELINE_ACQUISITION_QUEUE_TIMEOUT_SECONDS must be finite and positive")
        self._acquisition_semaphore = (
            threading.BoundedSemaphore(self.acquisition_slots) if self.acquisition_slots else None
        )
        # Why an accession's acquired files contradict its declared layout (see _note_layout_mismatch).
        self._acquisition_notes: Dict[str, str] = {}
        self._acquisition_notes_lock = threading.Lock()
        self.record_validation_required = os.environ.get(
            "AMALGKIT_PIPELINE_VALIDATE_FASTQ_RECORDS", "no"
        ).strip().lower() in {"1", "yes", "true"}
        local_ena_cache = os.environ.get("AMALGKIT_LOCAL_ENA_SCRATCH_DIR", "").strip()
        self.local_ena_scratch_dir = Path(local_ena_cache).expanduser().absolute() if local_ena_cache else None
        self._ena_data_root: Path | None = None
        self._ena_data_identity: tuple[int, int] | None = None
        self._ena_scratch_identity: tuple[int, int] | None = None
        self._ena_promotion_lock = threading.Lock()
        if self.local_ena_scratch_dir is not None or self.record_validation_required:
            root_text = os.environ.get("AMALGKIT_DATA_ROOT", "").strip()
            if not root_text:
                raise ValueError("ENA scratch requires an explicit AMALGKIT_DATA_ROOT")
            self._ena_data_root = Path(root_text).resolve(strict=True)
            if Path(root_text).absolute() != self._ena_data_root or (
                self.local_ena_scratch_dir is not None
                and self.local_ena_scratch_dir != self.local_ena_scratch_dir.resolve()
            ):
                raise StagingError("ENA roots must not contain symlink components")
            data_stat = self._ena_data_root.stat()
            self._ena_data_identity = (data_stat.st_dev, data_stat.st_ino)
            if self.local_ena_scratch_dir is not None:
                self.local_ena_scratch_dir.parent.resolve(strict=True)
                self.local_ena_scratch_dir.mkdir(exist_ok=True)
                scratch_stat = self.local_ena_scratch_dir.stat()
                self._ena_scratch_identity = (scratch_stat.st_dev, scratch_stat.st_ino)
        self.config_dir = Path(config_dir)
        self.log_dir = Path(log_dir)
        self.log_dir.mkdir(parents=True, exist_ok=True)

        # SQLite progress database
        self.db = ProgressDB(db_path) if db_path else ProgressDB()
        if cohort is not None:
            bound = self.db.cohort_binding().get("cohort_id")
            if bound not in (None, cohort.cohort_id):
                self.db.close()
                raise CohortError(f"progress database is bound to cohort {bound}, not {cohort.cohort_id}")

        # Raw reads are transient once the exact current-contract
        # quantification sidecar has been written.  Keep an explicit opt-out
        # for users who need raw retention, while making the default match the
        # external-volume storage contract.
        self.reclaim_raw_after_quant = os.environ.get(
            "AMALGKIT_RECLAIM_RAW_AFTER_QUANT", "yes"
        ).strip().lower() not in {"0", "false", "no", "off"}
        self.requantification_policy = os.environ.get("AMALGKIT_REQUANTIFICATION_POLICY", "preserve").strip().lower()
        if self.requantification_policy not in REQUANTIFICATION_POLICIES:
            raise ValueError(
                f"Invalid AMALGKIT_REQUANTIFICATION_POLICY={self.requantification_policy!r}; "
                f"expected one of {sorted(REQUANTIFICATION_POLICIES)}"
            )

        # The launcher checks the host volume before starting, but a long
        # campaign can still consume host-side cache/temp space after launch.
        # Keep the same floor in the worker loop so new downloads pause before
        # macOS reaches an emergency-full state.  On Linux/CI, a nonexistent
        # macOS path disables this host-specific guard without affecting the
        # external data-root throttle below.
        self.system_root = (
            None
            if self._ehd_write_policy is not None
            else Path(os.environ.get("AMALGKIT_SYSTEM_ROOT", "/System/Volumes/Data")).expanduser()
        )
        try:
            self.min_system_free_gb = max(
                0.0,
                float(os.environ.get("AMALGKIT_MIN_SYSTEM_FREE_GB", "4")),
            )
        except ValueError:
            logger.warning("Invalid AMALGKIT_MIN_SYSTEM_FREE_GB; using 4 GiB")
            self.min_system_free_gb = 4.0

        # Keep the mounted-data safety floor explicit and configurable. A
        # temporary test directory or a deliberately local scratch root may
        # live on the system volume rather than the production data volume;
        # callers can set this to zero for such isolated runs without
        # weakening the production default.
        try:
            self.min_external_free_gb = max(
                0.0,
                float(os.environ.get("AMALGKIT_MIN_EXTERNAL_FREE_GB", "8")),
            )
        except ValueError:
            logger.warning("Invalid AMALGKIT_MIN_EXTERNAL_FREE_GB; using 8 GiB")
            self.min_external_free_gb = 8.0

        def _duration_setting(name: str, default: int) -> int:
            raw = os.environ.get(name, str(default)).strip()
            try:
                value = int(raw)
            except ValueError:
                logger.warning("Invalid %s=%r; using %d seconds", name, raw, default)
                return default
            if value <= 0:
                logger.warning("Invalid %s=%r; using %d seconds", name, raw, default)
                return default
            return value

        def _positive_int_setting(name: str, default: int) -> int:
            raw = os.environ.get(name, str(default)).strip()
            try:
                value = int(raw)
            except ValueError:
                logger.warning("Invalid %s=%r; using %d", name, raw, default)
                return default
            if value <= 0:
                logger.warning("Invalid %s=%r; using %d", name, raw, default)
                return default
            return value

        def _nonnegative_int_setting(name: str, default: int) -> int:
            raw = os.environ.get(name, str(default)).strip()
            try:
                value = int(raw)
            except ValueError:
                logger.warning("Invalid %s=%r; using %d", name, raw, default)
                return default
            if value < 0:
                logger.warning("Invalid %s=%r; using %d", name, raw, default)
                return default
            return value

        def _nonnegative_float_setting(name: str, default: float) -> float:
            raw = os.environ.get(name, str(default)).strip()
            try:
                value = float(raw)
            except ValueError:
                logger.warning("Invalid %s=%r; using %.2f", name, raw, default)
                return default
            if value < 0:
                logger.warning("Invalid %s=%r; using %.2f", name, raw, default)
                return default
            return value

        self.ena_file_workers = int(os.environ.get("AMALGKIT_PIPELINE_ENA_FILE_WORKERS", "1"))
        if self.ena_file_workers not in (1, 2):
            raise ValueError("AMALGKIT_PIPELINE_ENA_FILE_WORKERS must be 1 or 2")
        self.download_timeout_seconds = _duration_setting("AMALGKIT_PIPELINE_DOWNLOAD_TIMEOUT_SECONDS", 7200)
        self.download_first_pass_seconds = int(os.environ.get("AMALGKIT_PIPELINE_DOWNLOAD_FIRST_PASS_SECONDS", "0"))
        if self.download_first_pass_seconds < 0:
            raise ValueError("AMALGKIT_PIPELINE_DOWNLOAD_FIRST_PASS_SECONDS must be nonnegative")
        self.download_speed_limit_bytes = _positive_int_setting("AMALGKIT_PIPELINE_DOWNLOAD_SPEED_LIMIT_BYTES", 1024)
        self.download_speed_time_seconds = _duration_setting("AMALGKIT_PIPELINE_DOWNLOAD_SPEED_TIME_SECONDS", 600)
        self.download_retries = _nonnegative_int_setting("AMALGKIT_PIPELINE_ENA_RETRIES", 5)
        self.download_integrity_retries = _nonnegative_int_setting("AMALGKIT_PIPELINE_ENA_INTEGRITY_RETRIES", 1)
        self.invalid_witness_max_bytes = (
            _nonnegative_int_setting("AMALGKIT_PIPELINE_INVALID_WITNESS_MAX_MB", 16) * 1024 * 1024
        )
        self.download_retry_delay_seconds = _duration_setting("AMALGKIT_PIPELINE_ENA_RETRY_DELAY_SECONDS", 5)
        self.ena_api_retries = _nonnegative_int_setting("AMALGKIT_PIPELINE_ENA_API_RETRIES", 2)
        self.ena_api_retry_delay_seconds = _duration_setting("AMALGKIT_PIPELINE_ENA_API_RETRY_DELAY_SECONDS", 2)
        self.ncbi_prefetch_first = os.environ.get("AMALGKIT_NCBI_PREFETCH_FIRST", "no").strip().lower() in {
            "1",
            "true",
            "yes",
            "on",
        }
        self.ncbi_prefetch_timeout_seconds = _duration_setting(
            "AMALGKIT_PIPELINE_NCBI_PREFETCH_TIMEOUT_SECONDS", self.download_timeout_seconds
        )
        self.fasterq_timeout_seconds = _duration_setting("AMALGKIT_PIPELINE_FASTQ_TIMEOUT_SECONDS", 7200)
        self.sra_validation_timeout_seconds = _duration_setting("AMALGKIT_PIPELINE_SRA_VALIDATE_TIMEOUT_SECONDS", 600)
        self.compression_timeout_seconds = _duration_setting("AMALGKIT_PIPELINE_COMPRESSION_TIMEOUT_SECONDS", 1800)
        self.quant_timeout_seconds = _duration_setting("AMALGKIT_PIPELINE_QUANT_TIMEOUT_SECONDS", 7200)
        # 0 (the default) disables the stall watchdog; it is a valid value, not a warning.
        self.quant_stall_timeout_seconds = _nonnegative_int_setting("AMALGKIT_PIPELINE_QUANT_STALL_TIMEOUT_SECONDS", 0)

        local_quant_cache = os.environ.get("AMALGKIT_LOCAL_QUANT_SCRATCH_DIR", "").strip()
        self.local_quant_scratch_dir = Path(local_quant_cache).expanduser() if local_quant_cache else None
        local_fasterq_cache = os.environ.get("AMALGKIT_LOCAL_FASTERQ_SCRATCH_DIR", "").strip()
        self.local_fasterq_scratch_dir = Path(local_fasterq_cache).expanduser() if local_fasterq_cache else None
        # NCBI documents an 8x--10x accession-size rule of thumb for combined
        # fasterq output and temporary space. Local extraction reserves that
        # full estimate plus one source-sized ephemeral copy when the
        # authoritative SRA is on another volume.
        self.local_fasterq_scratch_multiplier = _nonnegative_float_setting(
            "AMALGKIT_LOCAL_FASTERQ_SCRATCH_MULTIPLIER", 10.0
        )
        self.local_fasterq_scratch_min_gb = _nonnegative_float_setting("AMALGKIT_LOCAL_FASTERQ_SCRATCH_MIN_GB", 2.0)
        self.local_fasterq_remote_reserve_gb = _nonnegative_float_setting(
            "AMALGKIT_LOCAL_FASTERQ_REMOTE_RESERVE_GB", 16.0
        )
        self.local_fastq_stage = os.environ.get("AMALGKIT_LOCAL_FASTQ_STAGE", "").strip().lower() in {
            "1",
            "true",
            "yes",
            "on",
        }
        try:
            self.local_quant_scratch_reserve_gb = max(
                0.0,
                float(
                    os.environ.get(
                        "AMALGKIT_LOCAL_QUANT_SCRATCH_RESERVE_GB",
                        os.environ.get("AMALGKIT_MIN_SYSTEM_FREE_GB", "4"),
                    )
                ),
            )
        except ValueError:
            logger.warning("Invalid AMALGKIT_LOCAL_QUANT_SCRATCH_RESERVE_GB; using 4 GiB")
            self.local_quant_scratch_reserve_gb = 4.0

        # Staleness window for resumable .part payloads: a live curl refreshes
        # its .part continuously, so a payload untouched for this many minutes
        # is owned by no active download and is safe to evict under disk
        # pressure (its row's retry budget is already exhausted).
        self._PART_STALENESS_MINUTES = 10

        # Lock for filesystem operations (symlink cleanup, directory creation)
        self._fs_lock = threading.Lock()
        # Local FASTQ staging and optional fasterq temporary storage can share
        # the host volume. Reserve bytes through one ledger so concurrent
        # stages cannot all observe the same free-space headroom and overfill
        # that volume during a burst of large libraries.
        self._local_scratch_reservation_lock = threading.Lock()
        self._local_scratch_reserved_bytes = 0
        # FASTQ staging occupies the local volume for the entire quantification
        # workspace lifetime, not just while the copy is in progress. Keep one
        # reservation per sample workspace so concurrent quant slots cannot
        # collectively cross the configured host-volume reserve.
        self._local_fastq_stage_reservations: Dict[Path, int] = {}
        # A campaign may use more sample-task workers than quantification
        # slots so local validation and acquisition can feed Kallisto.  This
        # semaphore is installed by ``run_all`` and remains absent for small
        # direct-method calls and unit tests.
        self._quant_semaphore: threading.BoundedSemaphore | None = None
        # SRA extraction is usually external-volume I/O, not network I/O. Keep
        # its concurrency independent from sample workers so a large cache of
        # local archives cannot starve Kallisto or saturate the mounted disk.
        self._fasterq_semaphore: threading.BoundedSemaphore | None = None
        self._raw_validation_semaphore: threading.BoundedSemaphore | None = None

        # Direct calls to ``download_fastq`` (including recovery utilities and
        # tests) still receive an explicit bounded fallback profile.  A full
        # campaign replaces this with the requested profile at the beginning
        # of ``run_all``.
        self._resource_profile = build_pipeline_resource_profile(1, 1)

        # Configure logging to file as well
        timestamp = time.strftime("%Y%m%d_%H%M%S")
        log_file = self.log_dir / f"streaming_orchestrator_{timestamp}.log"
        file_handler = logging.FileHandler(log_file)
        file_handler.setFormatter(logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s"))
        logger.addHandler(file_handler)

    def query_ena_fastq_urls(self, srr_id: str) -> List[str]:
        """Query ENA API for direct FASTQ download URLs."""
        url = (
            f"https://www.ebi.ac.uk/ena/portal/api/filereport"
            f"?accession={srr_id}&result=read_run"
            f"&fields=run_accession,fastq_ftp&format=tsv"
        )
        try:
            with urllib.request.urlopen(url, timeout=30) as resp:  # nosec B310
                lines = resp.read().decode().strip().split("\n")
                if len(lines) < 2:
                    return []
                ftp_field = lines[1].split("\t")[1] if "\t" in lines[1] else ""
                if not ftp_field:
                    return []
                # Ensure protocol is https
                return [f"https://{p}" if not p.startswith("http") else p for p in ftp_field.split(";") if p]
        except Exception as e:
            logger.warning(f"ENA query failed for {srr_id}: {e}")
            return []

    def _validate_local_sra_archive(
        self,
        source: Path,
        srr_id: str,
    ) -> bool | None:
        """Validate one unchanged local SRA once and cache the result.

        ``False`` is authoritative evidence that the current source fingerprint
        is invalid. ``None`` means validation could not run or timed out, so
        fasterq-dump remains the fallback authority rather than rejecting a
        potentially usable archive.
        """

        cached_result = _cached_sra_validation_result(source)
        if cached_result is not None:
            logger.info(
                "Reusing current SRA validation witness for %s: %s",
                srr_id,
                "valid" if cached_result else "invalid",
            )
            return cached_result

        validator = shutil.which("vdb-validate")
        if validator is None:
            logger.warning(
                "vdb-validate is unavailable for local SRA preflight %s; " "continuing with fasterq-dump validation",
                srr_id,
            )
            return None

        try:
            before = _sra_source_fingerprint(source)
            result = _run_command_in_process_group(
                [validator, str(source)],
                timeout=self.sra_validation_timeout_seconds,
                env=build_sra_environment(_data_root()),
            )
            after = _sra_source_fingerprint(source)
        except subprocess.TimeoutExpired:
            logger.warning(
                "SRA preflight validation timed out for %s (>%s); " "continuing with fasterq-dump validation",
                srr_id,
                _format_duration(self.sra_validation_timeout_seconds),
            )
            return None
        except OSError as exc:
            logger.warning(
                "SRA preflight validation could not inspect %s; " "continuing with fasterq-dump validation: %s",
                srr_id,
                exc,
            )
            return None

        if before != after:
            logger.warning(
                "SRA source changed during preflight validation for %s; "
                "discarding the witness and continuing with fasterq-dump",
                srr_id,
            )
            return None

        detail = (result.stderr or result.stdout or "").strip()
        valid = result.returncode == 0
        try:
            _write_sra_validation_result(
                source,
                fingerprint=after,
                valid=valid,
                returncode=result.returncode,
                detail=detail,
            )
        except OSError as exc:
            logger.warning(
                "Could not persist SRA validation witness for %s: %s",
                srr_id,
                exc,
            )
        if valid:
            logger.info("Validated local SRA archive for %s before extraction.", srr_id)
        else:
            logger.warning(
                "Skipping invalid local SRA archive for %s after vdb-validate " "exit %d: %s",
                srr_id,
                result.returncode,
                detail[:300] or "no validator detail",
            )
        return valid

    def _extract_sra_to_fastq_bounded(
        self,
        source: Path | str,
        sample_dir: Path,
        out_dir: Path,
        srr_id: str,
        expected_paired: bool | None,
    ) -> bool:
        """Validate and extract one SRA run under the fallback-slot budget."""

        semaphore = self._fasterq_semaphore
        if semaphore is not None and isinstance(source, Path) and _cached_sra_validation_result(source) is False:
            logger.info("Skipping unchanged invalid SRA for %s without waiting for an extraction slot.", srr_id)
            return False
        if semaphore is None:
            if isinstance(source, Path) and self._validate_local_sra_archive(source, srr_id) is False:
                return False
            return self._extract_sra_to_fastq(
                source,
                sample_dir,
                out_dir,
                srr_id,
                expected_paired,
            )
        with semaphore:
            if isinstance(source, Path) and self._validate_local_sra_archive(source, srr_id) is False:
                return False
            return self._extract_sra_to_fastq(
                source,
                sample_dir,
                out_dir,
                srr_id,
                expected_paired,
            )

    def _extract_sra_to_fastq(
        self,
        source: Path | str,
        sample_dir: Path,
        out_dir: Path,
        srr_id: str,
        expected_paired: bool | None,
    ) -> bool:
        """Extract one local or remote SRA source into validated FASTQs.

        Extraction happens in a per-run staging directory.  Existing sample
        outputs are therefore never overwritten by a partially completed
        ``fasterq-dump`` invocation.  A local ``.sra`` archive remains in
        place until quantification succeeds, so interruption is resumable and
        the current-contract cleanup gate can reclaim it later.
        """

        temp_root = out_dir / ".fasterq_tmp" / srr_id
        staging_dir = temp_root / "output"
        fasterq_temp_dir = temp_root / "temp"
        fasterq_source: Path | str = source
        local_temp_root: Path | None = None
        local_temp_reservation = 0
        source_bytes: int | None = None
        # An interrupted producer can leave only this accession's ephemeral
        # workspace behind. The extraction semaphore is already held here, so
        # no current invocation can own the same staging path. Remove that
        # exact workspace before retrying; never touch the authoritative SRA
        # cache, resumable ENA partials, or promoted sample directory.
        if temp_root.is_symlink() or (temp_root.exists() and not temp_root.is_dir()):
            temp_root.unlink()
        elif temp_root.exists():
            shutil.rmtree(temp_root)
            logger.info("Removed stale fasterq staging workspace for %s: %s", srr_id, temp_root)
        staging_dir.mkdir(parents=True, exist_ok=True)
        if self.local_fasterq_scratch_dir is not None:
            if isinstance(source, Path):
                try:
                    source_bytes = source.stat().st_size
                except OSError as exc:
                    logger.warning(
                        "Unable to size local SRA source for SSD fasterq scratch %s: %s",
                        srr_id,
                        exc,
                    )
                else:
                    # Reserve the documented output/temp workspace estimate
                    # plus one ephemeral source copy. Copying a cross-volume
                    # SRA archive onto the same SSD prevents fasterq-dump's
                    # repeated/random source reads from contending with ENA,
                    # validation, promotion, and quantification on the mounted
                    # data volume.
                    local_temp_reservation = max(
                        int(self.local_fasterq_scratch_min_gb * (1024**3)),
                        math.ceil(source_bytes * (self.local_fasterq_scratch_multiplier + 1.0)),
                    )
            else:
                # A bare accession has no local source size before the SRA
                # Toolkit resolves it. Use a fixed conservative reservation;
                # fasterq-dump retains its own enabled output/temp size check.
                local_temp_reservation = int(self.local_fasterq_remote_reserve_gb * (1024**3))

            if local_temp_reservation and self._reserve_local_scratch(
                self.local_fasterq_scratch_dir,
                local_temp_reservation,
                purpose=f"Local fasterq scratch {srr_id}",
            ):
                try:
                    local_temp_root = self.local_fasterq_scratch_dir / srr_id
                    if local_temp_root.is_symlink() or (local_temp_root.exists() and not local_temp_root.is_dir()):
                        local_temp_root.unlink()
                    elif local_temp_root.exists():
                        shutil.rmtree(local_temp_root)
                    staging_dir = local_temp_root / "output"
                    fasterq_temp_dir = local_temp_root / "temp"
                    staging_dir.mkdir(parents=True, exist_ok=True)
                    fasterq_temp_dir.mkdir(parents=True, exist_ok=True)
                    if isinstance(source, Path) and source_bytes is not None:
                        copying_source: Path | None = None
                        try:
                            if not _paths_share_filesystem(source, local_temp_root):
                                source_dir = local_temp_root / "source"
                                source_dir.mkdir(parents=True, exist_ok=True)
                                staged_source = source_dir / source.name
                                copying_source = source_dir / (
                                    f".{source.name}.{os.getpid()}." f"{threading.get_ident()}.copying"
                                )
                                shutil.copyfile(source, copying_source)
                                copied_bytes = copying_source.stat().st_size
                                if copied_bytes != source_bytes:
                                    raise OSError(
                                        "local SRA source copy size mismatch: "
                                        f"expected {source_bytes}, observed {copied_bytes}"
                                    )
                                os.replace(copying_source, staged_source)
                                fasterq_source = staged_source
                                logger.info(
                                    "Staged SRA source for %s onto local fasterq "
                                    "scratch: %s (%.2f GiB); authoritative source "
                                    "retained at %s",
                                    srr_id,
                                    staged_source,
                                    source_bytes / (1024**3),
                                    source,
                                )
                        except OSError as exc:
                            if copying_source is not None:
                                copying_source.unlink(missing_ok=True)
                            fasterq_source = source
                            logger.warning(
                                "Could not stage SRA source for %s onto local "
                                "fasterq scratch; extracting from the authoritative "
                                "source instead: %s",
                                srr_id,
                                exc,
                            )
                except OSError as exc:
                    logger.warning(
                        "Unable to prepare local fasterq scratch for %s; " "using mounted-volume temporary storage: %s",
                        srr_id,
                        exc,
                    )
                    if local_temp_root is not None:
                        shutil.rmtree(local_temp_root, ignore_errors=True)
                    self._release_local_scratch(
                        local_temp_reservation,
                        purpose=f"Local fasterq scratch {srr_id}",
                    )
                    local_temp_root = None
                    local_temp_reservation = 0
                    fasterq_temp_dir = temp_root / "temp"
                else:
                    logger.info(
                        "Using local fasterq extraction scratch for %s: %s "
                        "(reserved %.2f GiB; ephemeral output and temp are "
                        "local; authoritative SRA/cache and promoted FASTQs "
                        "remain on %s)",
                        srr_id,
                        fasterq_temp_dir,
                        local_temp_reservation / (1024**3),
                        out_dir,
                    )
            else:
                local_temp_reservation = 0
        try:
            fasterq_temp_dir.mkdir(parents=True, exist_ok=True)
            fqd_cmd = [
                "fasterq-dump",
                str(fasterq_source),
                "--outdir",
                str(staging_dir),
                "--temp",
                str(fasterq_temp_dir),
                "--split-3",
                "--threads",
                str(self._resource_profile.fasterq_threads),
                "--skip-technical",
            ]
            fqd_result = _run_command_in_process_group(
                fqd_cmd,
                timeout=self.fasterq_timeout_seconds,
                env=build_sra_environment(_data_root()),
            )
            if fqd_result.returncode != 0:
                logger.error(
                    "fasterq-dump failed for %s (exit %d): %s",
                    srr_id,
                    fqd_result.returncode,
                    fqd_result.stderr[:200],
                )
                return False

            uncompressed = sorted(staging_dir.glob("*.fastq"))
            if not uncompressed:
                logger.error("fasterq-dump produced no FASTQ files for %s", srr_id)
                return False
            compressed_outputs: list[Path] = []
            for fastq in uncompressed:
                compression = _run_command_in_process_group(
                    pigz_command(fastq, threads=self._resource_profile.compression_threads, level=compression_level()),
                    timeout=self.compression_timeout_seconds,
                )
                gz_path = fastq.with_suffix(".fastq.gz")
                if compression.returncode != 0 or not gz_path.is_file() or gz_path.stat().st_size == 0:
                    logger.error("FASTQ compression failed for %s (%s)", srr_id, fastq.name)
                    return False
                compressed_outputs.append(gz_path)

            staged_inputs = _validated_local_fastq_inputs(staging_dir, srr_id, expected_paired)
            if not staged_inputs:
                logger.error(
                    "fasterq-dump output for %s did not match the declared library layout",
                    srr_id,
                )
                return False

            promoted = _promote_validated_fastq_inputs(staged_inputs, sample_dir, srr_id, expected_paired)
            if not promoted:
                # Keep useful extraction evidence outside the retry workspace,
                # whose owned temporary files are cleaned below and on retry.
                retained = Path(tempfile.mkdtemp(prefix=f"{srr_id}-conflict-", dir=staging_dir.parent.parent))
                staging_dir.rename(retained / "output")
                logger.error("Retained extracted FASTQs for %s at %s", srr_id, retained)
            return promoted
        except subprocess.TimeoutExpired:
            logger.error(
                "fasterq-dump timeout for %s (>%s)",
                srr_id,
                _format_duration(self.fasterq_timeout_seconds),
            )
            return False
        except OSError as exc:
            logger.error("SRA extraction filesystem error for %s: %s", srr_id, exc)
            return False
        finally:
            shutil.rmtree(temp_root, ignore_errors=True)
            if local_temp_root is not None:
                shutil.rmtree(local_temp_root, ignore_errors=True)
            if local_temp_reservation:
                self._release_local_scratch(
                    local_temp_reservation,
                    purpose=f"Local fasterq scratch {srr_id}",
                )

    def _check_ena_roots(self) -> None:
        if self._ena_data_root is None:
            raise StagingError("Acquisition data root is not configured")
        data_stat = self._ena_data_root.stat()
        if (data_stat.st_dev, data_stat.st_ino) != self._ena_data_identity:
            raise StagingError("Acquisition data root identity changed")
        if self.local_ena_scratch_dir is not None:
            scratch_stat = self.local_ena_scratch_dir.stat()
            if (scratch_stat.st_dev, scratch_stat.st_ino) != self._ena_scratch_identity:
                raise StagingError("ENA scratch root identity changed")

    def _reserve_ena_scratch(self, size: int, accession: str) -> bool:
        self._check_ena_roots()
        # Only the staged ENA download reserves, and it runs only with a scratch root.
        scratch_dir = cast(Path, self.local_ena_scratch_dir)
        return self._reserve_local_scratch(scratch_dir, size, purpose=f"ENA staging {accession}")

    def _note_layout_mismatch(self, srr_id: str, reason: str) -> None:
        """Remember that ``srr_id``'s acquired files contradict its declared layout, for its failure message."""

        logger.warning("Layout mismatch for %s: %s", srr_id, reason)
        with self._acquisition_notes_lock:
            self._acquisition_notes[srr_id] = reason[:300]

    def _layout_failure_suffix(self, srr_id: str) -> str:
        """``"; layout_mismatch: <reason>"`` when this run saw a layout contradiction for ``srr_id``, else ``""``.

        Appended to the acquisition failure message, it makes the progress DB
        classify the row ``layout_mismatch`` (quarantined for reconciliation,
        never re-downloaded at every start) instead of a transient transfer
        failure.
        """

        with self._acquisition_notes_lock:
            reason = self._acquisition_notes.get(srr_id)
        return f"; layout_mismatch: {reason}" if reason else ""

    def _forget_acquisition_note(self, srr_id: str) -> None:
        with self._acquisition_notes_lock:
            self._acquisition_notes.pop(srr_id, None)

    def _validate_acquisition_fastqs(
        self,
        directory: Path,
        accession: str,
        expected_paired: bool | None,
        deadline: float | None,
        *,
        assume_payload_valid: bool = False,
        require_records: bool | None = None,
        fresh_record_scan: bool = False,
    ) -> list[Path]:
        _check_deadline(deadline)
        semaphore = self._raw_validation_semaphore
        acquired = semaphore is None or (
            semaphore.acquire(timeout=max(0, deadline - time.monotonic()))
            if deadline is not None
            else semaphore.acquire()
        )
        if not acquired:
            raise StagingDeadline("acquisition validation slot time budget exhausted; inputs retained")
        try:
            _check_deadline(deadline)
            return _validated_local_fastq_inputs(
                directory,
                accession,
                expected_paired,
                assume_payload_valid=assume_payload_valid,
                deadline_at=deadline,
                require_records=require_records,
                fresh_record_scan=fresh_record_scan,
            )
        finally:
            if semaphore is not None:
                semaphore.release()

    def download_fastq(
        self,
        srr_id: str,
        out_dir: Path,
        expected_paired: bool | None = None,
        *,
        defer_ncbi_fallback: bool = False,
    ) -> bool | None:
        semaphore = self._acquisition_semaphore
        if semaphore is None:
            return self._download_fastq_admitted(
                srr_id, out_dir, expected_paired, defer_ncbi_fallback=defer_ncbi_fallback
            )
        queued_at = time.monotonic()
        if not semaphore.acquire(timeout=self.acquisition_queue_timeout_seconds):
            raise StagingDeadline("sample acquisition admission queue exhausted; inputs retained")
        try:
            admitted_at = time.monotonic()
            logger.info(
                "Acquisition admitted for %s after %.3f queue seconds; slots=%d",
                srr_id,
                admitted_at - queued_at,
                self.acquisition_slots,
            )
            return self._download_fastq_admitted(
                srr_id,
                out_dir,
                expected_paired,
                defer_ncbi_fallback=defer_ncbi_fallback,
                admitted_at=admitted_at,
            )
        finally:
            semaphore.release()

    def _download_fastq_admitted(
        self,
        srr_id: str,
        out_dir: Path,
        expected_paired: bool | None = None,
        *,
        defer_ncbi_fallback: bool = False,
        admitted_at: float | None = None,
    ) -> bool | None:
        """Download FASTQ files using ENADownloader.

        Uses the ENA Portal API to discover FASTQ URLs and downloads
        via curl with retries. When ``defer_ncbi_fallback`` is true, a run
        without a local SRA archive returns ``None`` after ENA exhaustion so a
        dedicated fallback executor can process it without occupying an ENA
        worker.

        Args:
            srr_id: SRA accession ID (e.g., SRR12345)
            out_dir: Base directory for downloads (sample subdir created automatically)

        Returns:
            True when reusable FASTQs are ready, False after all enabled
            acquisition routes fail, or None when NCBI acquisition was
            deliberately deferred to the scheduler's dedicated fallback lane.
        """
        # Admission logging consumes service time too; only the semaphore queue is excluded.
        deadline = (
            (admitted_at if admitted_at is not None else time.monotonic()) + self.download_first_pass_seconds
            if self.download_first_pass_seconds
            else None
        )
        sample_dir = out_dir / srr_id
        guarded_acquisition = (
            self.local_ena_scratch_dir is not None or self.record_validation_required or deadline is not None
        )
        if guarded_acquisition:
            if self._ena_data_root is not None:
                self._check_ena_roots()
            destination_root = self._ena_data_root or _data_root().resolve(strict=True)
            destination_stat = destination_root.stat()
            destination_identity = (destination_stat.st_dev, destination_stat.st_ino)
            _plain_path(sample_dir)
            if not sample_dir.absolute().is_relative_to(destination_root):
                raise StagingError("ENA destination escapes the selected data root")
            if sample_dir.exists() and not sample_dir.is_dir():
                raise StagingError("ENA destination conflicts with a non-directory")
            _plain_path(_raw_validation_marker(sample_dir, srr_id))
            marker = _raw_validation_marker(sample_dir, srr_id)
            legacy_temporary = marker.with_name(f".{marker.name}.tmp")
            _plain_path(legacy_temporary)
            if any(
                path.exists() and (not path.is_file() or path.stat().st_size > 65536)
                for path in (marker, legacy_temporary)
            ):
                raise StagingError("Unsafe raw validation witness object; existing inputs preserved")
            if sample_dir.is_dir() and any(path.is_symlink() for path in sample_dir.iterdir()):
                raise StagingError("ENA sample contains a symlink; existing inputs preserved")
        # Remove broken symlinks in the path chain (amalgkit preprocessing creates
        # symlinks like work/getfastq -> /local/path that are invalid inside Docker).
        # Use a lock because 28+ threads may race on the same parent directory.
        acquired_fs = (
            self._fs_lock.acquire(timeout=max(0, deadline - time.monotonic()))
            if deadline is not None
            else self._fs_lock.acquire()
        )
        if not acquired_fs:
            raise StagingDeadline("acquisition filesystem slot time budget exhausted; inputs retained")
        try:
            for check_path in [out_dir, sample_dir]:
                if check_path.is_symlink() and not check_path.exists():
                    logger.warning(f"Removing broken symlink: {check_path} -> {os.readlink(check_path)}")
                    check_path.unlink()
            if guarded_acquisition:
                prepare_sample_directory(sample_dir, destination_root, destination_identity)
            else:
                sample_dir.mkdir(parents=True, exist_ok=True)
        finally:
            self._fs_lock.release()

        # A valid local sample is the cheapest and most reproducible input.
        # Do this before disk throttling and ENA URL discovery: large active
        # campaigns otherwise spend every worker slot waiting on the ENA API
        # even when the FASTQs are already mounted locally.
        local_fastqs = self._validate_acquisition_fastqs(
            sample_dir,
            srr_id,
            expected_paired,
            deadline,
            require_records=self.record_validation_required or self.local_ena_scratch_dir is not None,
        )
        if local_fastqs:
            logger.info(
                "Reusing %d validated local FASTQ file(s) for %s; skipping ENA acquisition.",
                len(local_fastqs),
                srr_id,
            )
            return True
        rejected_inputs = _raw_rejection_reason(sample_dir, srr_id)
        if rejected_inputs is not None:
            self._note_layout_mismatch(srr_id, f"local inputs: {rejected_inputs}")

        local_sra = _local_sra_path(sample_dir, out_dir, srr_id)
        local_sra_attempted = False
        if local_sra is not None and _cached_sra_validation_result(local_sra) is not False:
            if defer_ncbi_fallback:
                return None
            logger.info(
                "Extracting existing local SRA archive for %s before any network acquisition: %s",
                srr_id,
                local_sra,
            )
            local_sra_attempted = True
            if self._extract_sra_to_fastq_bounded(
                local_sra,
                sample_dir,
                out_dir,
                srr_id,
                expected_paired,
            ):
                return True
            logger.warning(
                "Local SRA validation/extraction did not produce usable FASTQ files "
                "for %s; continuing with ENA/NCBI recovery.",
                srr_id,
            )

        # Disk space throttle: Prevent starting new enormous downloads if
        # either the external data volume or the host system volume is nearing
        # capacity.  Quantification workers can reclaim external FASTQ space
        # while a paused downloader waits, so this remains resume-safe.
        throttle_attempts = 0
        while True:
            _check_deadline(deadline)
            try:
                external_free_gb = shutil.disk_usage(out_dir).free / (1024**3)
                system_free_gb = None
                if self._ehd_write_policy is not None:
                    self._ehd_write_policy.revalidate(out_dir)
                if self.system_root is not None and self.system_root.is_dir():
                    system_free_gb = shutil.disk_usage(self.system_root).free / (1024**3)
                system_low = system_free_gb is not None and system_free_gb < self.min_system_free_gb
                if external_free_gb < self.min_external_free_gb or system_low:
                    self._maybe_reclaim_part_files(out_dir)
                    if throttle_attempts % 6 == 0:  # Log every ~5 mins
                        system_text = f"; system {system_free_gb:.1f} GiB free" if system_free_gb is not None else ""
                        logger.warning(
                            "[Disk Throttle] External volume %.1f GiB free "
                            "(threshold %.1f GiB)%s; "
                            "worker pausing to allow quantification/reclamation to proceed...",
                            external_free_gb,
                            self.min_external_free_gb,
                            system_text,
                        )
                    time.sleep(min(60, max(0, deadline - time.monotonic())) if deadline is not None else 60)
                    throttle_attempts += 1
                else:
                    if throttle_attempts > 0:
                        logger.info(
                            "[Disk Throttle] Space available again (external %.1f GiB%s). " "Resuming download for %s.",
                            external_free_gb,
                            (f"; system {system_free_gb:.1f} GiB" if system_free_gb is not None else ""),
                            srr_id,
                        )
                    break
            except StagingDeadline:
                raise
            except Exception as e:
                logger.error(f"Failed to check disk space: {e}")
                if (
                    self._ehd_write_policy is not None
                    or self.local_ena_scratch_dir is not None
                    or self.record_validation_required
                    or deadline is not None
                ):
                    raise StagingError("acquisition capacity check failed; inputs retained") from e
                break  # Proceed anyway if check fails

        downloader = ENADownloader(
            timeout=self.download_timeout_seconds,
            retries=self.download_retries,
            integrity_retries=self.download_integrity_retries,
            invalid_witness_max_bytes=self.invalid_witness_max_bytes,
            retry_delay_seconds=self.download_retry_delay_seconds,
            speed_limit_bytes=self.download_speed_limit_bytes,
            speed_time_seconds=self.download_speed_time_seconds,
            api_retries=self.ena_api_retries,
            api_retry_delay_seconds=self.ena_api_retry_delay_seconds,
            file_workers=self.ena_file_workers,
            run_timeout_seconds=self.download_first_pass_seconds,
            api_base_url=os.environ.get("AMALGKIT_ENA_API_BASE_URL", ENADownloader.ENA_API_BASE),
        )
        staged_success = False
        if self.local_ena_scratch_dir is None:
            urls = (
                downloader.get_fastq_urls(srr_id, deadline_at=deadline)
                if deadline is not None
                else downloader.get_fastq_urls(srr_id)
            )
            if self._ena_data_root is not None:
                self._check_ena_roots()
            success, message, downloaded_files = (
                downloader._download_files(
                    urls, sample_dir, deadline_at=deadline, require_existing_output=guarded_acquisition
                )
                if urls
                else (False, downloader.last_metadata_error or "Not found on ENA", [])
            )
        else:
            urls = (
                downloader.get_fastq_urls(srr_id, deadline_at=deadline)
                if deadline is not None
                else downloader.get_fastq_urls(srr_id)
            )
            outcome = None
            if urls:

                def validate_staged(directory: Path) -> list[Path]:
                    self._check_ena_roots()
                    return self._validate_acquisition_fastqs(
                        directory, srr_id, expected_paired, deadline, assume_payload_valid=True
                    )

                try:
                    outcome = download_staged_ena(
                        downloader,
                        urls,
                        sample_dir,
                        self.local_ena_scratch_dir,
                        reserve=lambda size: self._reserve_ena_scratch(size, srr_id),
                        release=lambda size: self._release_local_scratch(size, purpose=f"ENA staging {srr_id}"),
                        promotion_lock=self._ena_promotion_lock,
                        validate=validate_staged,
                        expected_paired=expected_paired,
                        deadline_at=deadline,
                        minimum_free_bytes=math.ceil(self.min_external_free_gb * 1024**3),
                        validation_lock=self._raw_validation_semaphore,
                    )
                except StagingDeadline as exc:
                    outcome = False, str(exc), []
                staged_success = outcome is not None and outcome[0]
            if outcome is None:
                self._check_ena_roots()
                outcome = (
                    downloader._download_files(urls, sample_dir, deadline_at=deadline, require_existing_output=True)
                    if urls
                    else (False, downloader.last_metadata_error or "Not found on ENA", [])
                )
            success, message, downloaded_files = outcome
        if success:
            for fq in downloaded_files:
                try:
                    file_size = fq.stat().st_size
                except OSError:
                    file_size = 0
                sz_gb = file_size / (1024**3)
                logger.info(f"Downloaded {fq.name}: {sz_gb:.2f} GB")
                if sz_gb == 0:
                    logger.error(f"Downloaded file {fq.name} is empty (0.00 GB). Marking as failed.")
                    success = False
            if success:
                # ENA's transfer validator proves the downloaded payload is a
                # readable gzip stream, but the current-contract witness also
                # needs to bind the accession to the expected single/paired
                # FASTQ layout.  Direct ENA success used to return before this
                # marker was written, leaving an otherwise valid cloud result
                # non-promotable after raw-input reclamation.
                validated_fastqs = self._validate_acquisition_fastqs(
                    sample_dir,
                    srr_id,
                    expected_paired,
                    deadline,
                    assume_payload_valid=True,
                    require_records=not staged_success
                    and (
                        self.local_ena_scratch_dir is not None
                        or self.record_validation_required
                        or deadline is not None
                    ),
                )
                if not validated_fastqs:
                    logger.error(
                        "Downloaded FASTQ layout for %s does not match the declared library layout",
                        srr_id,
                    )
                    contradiction = _raw_rejection_reason(sample_dir, srr_id) or _layout_contradiction(
                        downloaded_files, srr_id, expected_paired
                    )
                    if contradiction is not None:
                        self._note_layout_mismatch(srr_id, f"ENA files: {contradiction}")
                    success = False
                else:
                    downloaded_files = validated_fastqs
        if not success and self.download_first_pass_seconds:
            raise FirstPassDeferred(srr_id, message)
        if not success:
            logger.warning(f"ENA download failed for {srr_id}: {message}")
            logger.info(f"Falling back to NCBI fasterq-dump for {srr_id}...")
            # Retain completed accession FASTQs as per-mate checkpoints. A
            # paired layout remains non-quantifiable until both mates pass the
            # full layout validator, but discarding mate 1 when mate 2 fails
            # needlessly repeats a potentially multi-gigabyte transfer.
            # Retain ENA ``.part`` and diagnostic ``.invalid`` files as well.
            for f in sample_dir.iterdir():
                if (
                    f.is_file()
                    and not _is_diagnostic_transfer_artifact(f)
                    and not f.name.endswith((".sra", ".sralite"))
                    and not _is_accession_fastq_file(f, srr_id)
                ):
                    f.unlink(missing_ok=True)

            if defer_ncbi_fallback:
                logger.info(
                    "Deferring NCBI fallback for %s to the dedicated extraction queue.",
                    srr_id,
                )
                return None
            if local_sra_attempted:
                logger.warning(
                    "Skipping duplicate NCBI fallback extraction for %s: the campaign-local "
                    "SRA archive already failed in this invocation and remains available "
                    "for an idempotent retry.",
                    srr_id,
                )
            else:
                try:
                    success = self._download_fastq_ncbi_only(srr_id, out_dir, expected_paired)
                    if success:
                        logger.info(f"NCBI fallback succeeded for {srr_id}")
                except Exception as e:
                    logger.error(f"NCBI fasterq-dump exception for {srr_id}: {e}")

        if not success:
            # Final cleanup on total failure. Retain canonical completed mates,
            # resumable ENA ``.part`` files, and diagnostic ``.invalid``
            # payloads. A later ENA attempt can skip each validated complete
            # mate, and a fallback promotion can keep it while supplying the
            # missing mate.
            for f in sample_dir.iterdir():
                if (
                    f.is_file()
                    and not _is_diagnostic_transfer_artifact(f)
                    and not f.name.endswith((".sra", ".sralite"))
                    and not _is_accession_fastq_file(f, srr_id)
                ):
                    f.unlink(missing_ok=True)
            # Under disk pressure the retention above becomes the lane's
            # dominant ENOSPC term (the 2026-09-28/29 production tail lanes
            # drowned in ~250 terminal rows' partial payloads). When free
            # space is already low, this row's resumable partial is worth
            # less than the scratch it holds: reclaim the space-pressure
            # share of the .part debt, this row's own payload included.
            self._maybe_reclaim_part_files(sample_dir.parent)

        return success

    def _evict_stale_part_files(self, out_dir: Path) -> float:
        """Delete resumable ``.part`` payloads that no active writer owns.

        ``.part`` files exist to resume interrupted multi-gigabyte transfers
        across invocations that share the disk. On ephemeral cloud scratch
        (and on any disk under space pressure) a terminal-failed row's
        partial payload is pure disk debt: its retry budget is exhausted, so
        a recovery pass restarts the transfer from byte 0 regardless.
        Accumulated ``.part`` debt is what exhausted the 256 GB scratch on
        the 2026-09-28/29 production tail lanes (~250 failed rows x 1-3 GB).

        Safety rule: only evict parts whose mtime is older than
        ``_PART_STALENESS_MINUTES``. A live curl refreshes its ``.part``
        continuously, so a fresh mtime means an active download owns the
        payload and it must survive.

        Returns the number of bytes freed.
        """

        freed = 0
        now = time.time()
        cutoff = now - self._PART_STALENESS_MINUTES * 60
        try:
            candidates = list(out_dir.rglob("*.part"))
        except OSError as exc:
            logger.warning("Unable to scan %s for stale .part payloads: %s", out_dir, exc)
            return 0.0
        for part in candidates:
            try:
                stat = part.stat()
                if stat.st_mtime > cutoff:
                    continue
                size = stat.st_size
            except OSError:
                continue
            try:
                part.unlink(missing_ok=True)
            except OSError:
                continue
            freed += size
        if freed:
            logger.warning(
                "[Disk Reclaim] Evicted %.2f GB of stale resumable .part payloads "
                "under %s (terminal-failed rows re-download from byte 0 on the "
                "next pass; active downloads are untouched).",
                freed / (1024**3),
                out_dir,
            )
        return freed

    def _maybe_reclaim_part_files(self, out_dir: Path) -> None:
        """Evict stale ``.part`` payloads only when free space is low."""

        if (
            self.local_ena_scratch_dir is not None
            or self.download_first_pass_seconds
            or os.environ.get("AMALGKIT_PRESERVE_RESUMABLE_PARTS", "no").lower() in {"1", "yes", "true"}
        ):
            return
        try:
            if shutil.disk_usage(out_dir).free / (1024**3) >= self.min_external_free_gb:
                return
        except OSError as exc:
            logger.warning("Unable to check %s capacity for .part reclaim: %s", out_dir, exc)
            return
        self._evict_stale_part_files(out_dir)

    def _prefetch_sra_accession(self, srr_id: str, sample_dir: Path | None = None) -> Path | None:
        """Resume and verify one NCBI accession before FASTQ extraction."""

        cached = (
            _local_sra_path(sample_dir, sample_dir.parent, srr_id)
            if sample_dir is not None
            else _campaign_cached_sra_path(srr_id)
        )
        if cached is not None and _cached_sra_validation_result(cached) is not False:
            return cached
        prefetch = shutil.which("prefetch")
        if prefetch is None:
            logger.error("NCBI prefetch is unavailable for %s", srr_id)
            return None

        root = (
            sample_dir / ".prefetch" / srr_id
            if sample_dir is not None
            else _data_root() / ".sra-cache" / "prefetch" / srr_id
        )
        root.mkdir(parents=True, exist_ok=True)
        available = shutil.disk_usage(root).free - int(self.min_external_free_gb * 1024**3)
        if available < 1024:
            logger.error("Insufficient bulk capacity to prefetch %s", srr_id)
            return None
        for rejected in (
            root / nested / f"{srr_id}{suffix}" for nested in ("", srr_id) for suffix in (".sra", ".sralite")
        ):
            if rejected.is_file() and _cached_sra_validation_result(rejected) is False:
                rejected.rename(rejected.with_name(f"{rejected.name}.invalid-{time.time_ns()}"))
        command = [
            prefetch,
            "--transport",
            "http",
            "--resume",
            "yes",
            "--max-size",
            str(available // 1024),
            "--verify",
            "yes",
            "--output-directory",
            str(root),
            srr_id,
        ]
        logger.info("Prefetching NCBI accession %s with resumable verified transport", srr_id)
        try:
            result = _run_command_in_process_group(
                command,
                timeout=self.ncbi_prefetch_timeout_seconds,
                env=build_sra_environment(_data_root()),
            )
        except subprocess.TimeoutExpired:
            logger.error("NCBI prefetch timed out for %s", srr_id)
            return None
        if result.returncode != 0:
            detail = (result.stderr or result.stdout or "").strip()
            logger.error("NCBI prefetch failed for %s (exit %d): %s", srr_id, result.returncode, detail[:300])
            return None

        candidates = [
            root / nested / f"{srr_id}{suffix}" for nested in ("", srr_id) for suffix in _SRA_ARCHIVE_SUFFIXES
        ]
        for candidate in candidates:
            if candidate.is_file() and candidate.stat().st_size > 0:
                logger.info("Prefetch completed for %s: %s", srr_id, candidate)
                return candidate
        cached = _campaign_cached_sra_path(srr_id)
        if cached is not None and _cached_sra_validation_result(cached) is not False:
            return cached
        logger.error(
            "NCBI prefetch reported success but produced no SRA archive for %s%s",
            srr_id,
            _prefetch_no_archive_reason(f"{result.stdout or ''}\n{result.stderr or ''}", cached),
        )
        return None

    def _download_fastq_ncbi_only(
        self,
        srr_id: str,
        out_dir: Path,
        expected_paired: bool | None,
    ) -> bool:
        """Run only the bounded NCBI extraction lane for a deferred sample."""

        sample_dir = out_dir / srr_id
        with self._fs_lock:
            if self._ena_data_root is not None:
                self._check_ena_roots()
                # _check_ena_roots() has just verified the recorded root identity.
                data_identity = cast("tuple[int, int]", self._ena_data_identity)
                prepare_sample_directory(sample_dir, self._ena_data_root, data_identity)
                _plain_path(_raw_validation_marker(sample_dir, srr_id))
                if any(path.is_symlink() for path in sample_dir.iterdir()):
                    raise StagingError("Recovery sample contains a symlink; existing inputs preserved")
            else:
                sample_dir.mkdir(parents=True, exist_ok=True)
        if self._raw_validation_semaphore is None:
            local_fastqs = _validated_local_fastq_inputs(
                sample_dir,
                srr_id,
                expected_paired,
                require_records=self.record_validation_required or self.local_ena_scratch_dir is not None,
            )
        else:
            with self._raw_validation_semaphore:
                local_fastqs = _validated_local_fastq_inputs(
                    sample_dir,
                    srr_id,
                    expected_paired,
                    require_records=self.record_validation_required or self.local_ena_scratch_dir is not None,
                )
        if local_fastqs:
            logger.info(
                "Deferred fallback found %d validated FASTQ file(s) for %s; skipping extraction.",
                len(local_fastqs),
                srr_id,
            )
            return True

        # A complete pair rejected above is a canonical conflict, not missing
        # input. Extraction cannot silently combine it with another source.
        if (self.record_validation_required or self.local_ena_scratch_dir is not None) and (
            (_raw_fastq_path(sample_dir, srr_id, 1) is not None and _raw_fastq_path(sample_dir, srr_id, 2) is not None)
            or _raw_fastq_path(sample_dir, srr_id) is not None
        ):
            logger.error("Recovery refuses invalid complete FASTQs for %s; inputs preserved", srr_id)
            rejected_inputs = _raw_rejection_reason(sample_dir, srr_id)
            if rejected_inputs is not None:
                # Deterministic: the same bytes failed layout validation. Other rejections (an I/O error during
                # the scan) stay transient and go through the backoff budget.
                self._note_layout_mismatch(srr_id, f"recovery refuses rejected complete FASTQs: {rejected_inputs}")
            return False

        local_sra = _local_sra_path(sample_dir, out_dir, srr_id)
        if self.ncbi_prefetch_first and (local_sra is None or _cached_sra_validation_result(local_sra) is False):
            local_sra = self._prefetch_sra_accession(srr_id, sample_dir)
        source: Path | str = local_sra if local_sra is not None else srr_id
        if self.ncbi_prefetch_first and local_sra is None:
            return False
        logger.info("Running dedicated NCBI fallback for %s from %s", srr_id, source)
        success = self._extract_sra_to_fastq_bounded(
            source,
            sample_dir,
            out_dir,
            srr_id,
            expected_paired,
        )
        if success:
            logger.info("Dedicated NCBI fallback succeeded for %s", srr_id)
        return success

    def quant_sample(
        self,
        config_path: Path,
        batch_index: int,
        species_name: str,
        threads: int,
        srr_id: str,
        expected_paired: bool | None = None,
        *,
        metadata_path: Path | None = None,
    ) -> tuple:
        """Run amalgkit quant for a single sample using --batch.

        ``metadata_path`` (a frozen-cohort task table in the species'
        ``work/metadata``) replaces the selected metadata table. Either way the
        ``batch_index``-th sampled row must be ``srr_id``; otherwise the run is
        refused instead of quantifying another sample.

        Returns:
            Tuple of (success: bool, error_msg: str or None).
        """
        with open(config_path) as f:
            cfg = yaml.safe_load(f)

        work_dir = _species_work_dir(species_name)
        if metadata_path is not None:
            table = Path(metadata_path)
            if table.parent != work_dir / "metadata" or table.is_symlink() or not table.is_file():
                return False, f"Quant task table is not a regular file in {work_dir / 'metadata'}: {table}"
            meta_path: Optional[str] = str(table)
        else:
            meta_path = _resolve_quant_metadata_path(work_dir)

        if not meta_path:
            logger.error(f"No metadata found for {species_name} batch {batch_index}")
            return False, f"No metadata found for {species_name}"
        batch_run = quant_batch_run(Path(meta_path), batch_index)
        if batch_run != srr_id:
            return False, (
                f"Quant batch mismatch: row {batch_index} of {Path(meta_path).name} is {batch_run!r}, not {srr_id}"
            )

        # Direct ENA acquisition bypasses Amalgkit's getfastq stage.  For
        # public metadata rows with missing/zero spot statistics, derive the
        # quant-only fallback from the validated FASTQ before invoking
        # Amalgkit 0.16.60; otherwise it rejects usable reads before Kallisto.
        metadata_file = Path(meta_path)
        sample_fastq_dir = work_dir / "getfastq" / srr_id
        if _metadata_requires_fastq_input_stats(metadata_file, srr_id):
            _write_fastq_input_stats(sample_fastq_dir, srr_id)

        with metadata_file.open(encoding="utf-8", newline="") as handle:
            reference_targets = [
                row["scientific_name"] for row in csv.DictReader(handle, delimiter="\t") if row["run"] == srr_id
            ]
        if len(reference_targets) != 1 or not reference_targets[0].strip():
            return False, f"No unique metadata reference target for {srr_id}"
        cmd = _build_quant_command(
            cfg, species_name, batch_index, threads, meta_path, reference_target=reference_targets[0]
        )
        scratch_context = self._prepare_local_quant_workspace(
            work_dir,
            species_name,
            srr_id,
            meta_path,
            cmd,
            expected_paired=expected_paired,
        )
        run_cmd = scratch_context[0] if scratch_context else cmd
        scratch_root = scratch_context[1] if scratch_context else None

        log_path = self.log_dir / f"{species_name}_quant.log"

        try:
            with open(log_path, "a") as log_f:
                log_f.write(
                    f"\n[{time.strftime('%Y-%m-%d %H:%M:%S')}] Quant batch {batch_index} START: {' '.join(run_cmd)}\n"
                )

            quant_output_dirs = [work_dir / "quant" / srr_id]
            if scratch_root is not None:
                quant_output_dirs.append(Path(scratch_root))
            result = _run_command_in_process_group(
                run_cmd,
                self.quant_timeout_seconds,
                output_dirs=quant_output_dirs,
                stall_timeout=self.quant_stall_timeout_seconds,
            )

            with open(log_path, "a") as log_f:
                log_f.write(
                    f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] Quant batch {batch_index} END (Exit {result.returncode})\n"
                )
                if result.stdout:
                    log_f.write(f"-- STDOUT --\n{result.stdout}\n")
                if result.stderr:
                    log_f.write(f"-- STDERR --\n{result.stderr}\n")

            if result.returncode == 0:
                if scratch_root is not None and not self._promote_local_quant_output(work_dir, srr_id, scratch_root):
                    return False, f"Unable to promote local quantification output for {srr_id}"
                sample_dir = work_dir / "quant" / srr_id
                quantification_file = find_quantification_file(sample_dir, srr_id)
                if quantification_file is None:
                    return False, f"No abundance output produced for {srr_id}"
                context = self._reuse_context(species_name)
                cohort_context = context if context is not None and context.is_cohort else None
                try:
                    reference_manifest = work_dir / "reference" / "reference_aliases.json"
                    if cohort_context is not None:
                        # The previous output's witness names the sidecar about to be replaced.
                        discard_reference_binding(sample_dir)
                    write_quant_provenance(
                        sample_dir,
                        species=species_name,
                        run_accession=srr_id,
                        config_path=config_path,
                        command=run_cmd,
                        reference_manifest_path=reference_manifest if reference_manifest.exists() else None,
                        quantification_file=quantification_file,
                    )
                except OSError as exc:
                    return False, f"Unable to write current-method provenance: {exc}"
                if cohort_context is not None:
                    # Bind the new output to the cohort index bytes now, while the reference is known to match,
                    # so its reuse never depends on reference_aliases.json staying byte-identical.
                    try:
                        bind_produced_output(
                            sample_dir, cohort_context.expectation(srr_id), index_digests=cohort_context.index_digests
                        )
                    except (OSError, ValueError) as exc:
                        return False, f"Unable to bind the output to cohort {cohort_context.cohort_id}: {exc}"
                return True, None

            # Extract meaningful error from amalgkit output
            error_msg = "Quantification Failed"
            for line in (result.stdout + result.stderr).splitlines():
                line_lower = line.strip().lower()
                if any(
                    kw in line_lower for kw in ["error", "exiting", "no sample", "not found", "failed", "exception"]
                ):
                    error_msg = line.strip()[:120]
                    break
            return False, error_msg
        except subprocess.TimeoutExpired as exc:
            stalled = exc.output == b"stalled"
            label = "Quant stalled" if stalled else "Quant timeout"
            logger.error(
                "%s batch %s after %s",
                label,
                batch_index,
                _format_duration(self.quant_stall_timeout_seconds if stalled else self.quant_timeout_seconds),
            )
            with open(log_path, "a") as log_f:
                log_f.write(f"{'STALL' if stalled else 'TIMEOUT'} detected for batch {batch_index}\n")
            if stalled:
                stall_desc = _format_duration(self.quant_stall_timeout_seconds)
                return False, f"Quant stalled (no output growth for {stall_desc})"
            return False, f"Quant timeout (>{_format_duration(self.quant_timeout_seconds)})"
        except Exception as e:
            logger.error(f"Quant exception batch {batch_index}: {e}")
            with open(log_path, "a") as log_f:
                log_f.write(f"EXCEPTION: {e}\n")
            return False, f"Exception: {e}"
        finally:
            if scratch_root is not None:
                self._cleanup_local_quant_workspace(scratch_root)

    def _reserve_local_scratch(
        self,
        scratch_dir: Path,
        payload_bytes: int,
        *,
        purpose: str,
    ) -> bool:
        """Reserve host scratch capacity across all local-I/O optimizations."""

        if payload_bytes <= 0:
            return False
        reserve_bytes = int(self.local_quant_scratch_reserve_gb * (1024**3))
        with self._local_scratch_reservation_lock:
            try:
                scratch_dir.mkdir(parents=True, exist_ok=True)
                free_bytes = shutil.disk_usage(scratch_dir).free
            except OSError as exc:
                logger.warning(
                    "Unable to inspect %s capacity at %s: %s",
                    purpose,
                    scratch_dir,
                    exc,
                )
                return False
            available_after_reservation = free_bytes - self._local_scratch_reserved_bytes - payload_bytes
            if available_after_reservation <= reserve_bytes:
                logger.info(
                    "%s would leave %.2f GiB free, below the %.2f GiB reserve; "
                    "using mounted-volume storage for this sample.",
                    purpose,
                    available_after_reservation / (1024**3),
                    self.local_quant_scratch_reserve_gb,
                )
                return False
            self._local_scratch_reserved_bytes += payload_bytes
            logger.info(
                "[%s] reserved %.2f GiB for this sample; "
                "in-flight reservation %.2f GiB, free before copy %.2f GiB, reserve %.2f GiB",
                purpose,
                payload_bytes / (1024**3),
                self._local_scratch_reserved_bytes / (1024**3),
                free_bytes / (1024**3),
                self.local_quant_scratch_reserve_gb,
            )
            return True

    def _release_local_scratch(self, payload_bytes: int, *, purpose: str) -> None:
        """Release capacity from the shared host-scratch reservation ledger."""

        with self._local_scratch_reservation_lock:
            self._local_scratch_reserved_bytes = max(
                0,
                self._local_scratch_reserved_bytes - max(payload_bytes, 0),
            )
            logger.debug(
                "[%s] released %.2f GiB; in-flight reservation %.2f GiB",
                purpose,
                max(payload_bytes, 0) / (1024**3),
                self._local_scratch_reserved_bytes / (1024**3),
            )

    def _reserve_local_fastq_stage(self, payload_bytes: int) -> bool:
        """Reserve local scratch capacity for one FASTQ staging copy."""

        if self.local_quant_scratch_dir is None:
            return False
        return self._reserve_local_scratch(
            self.local_quant_scratch_dir,
            payload_bytes,
            purpose="Local FASTQ stage",
        )

    def _release_local_fastq_stage(self, payload_bytes: int) -> None:
        """Release a temporary local FASTQ staging reservation."""

        self._release_local_scratch(
            payload_bytes,
            purpose="Local FASTQ stage",
        )

    def _retain_local_fastq_stage_reservation(self, scratch_root: Path, payload_bytes: int) -> None:
        """Keep a successful FASTQ-stage reservation until workspace cleanup."""

        with self._local_scratch_reservation_lock:
            self._local_fastq_stage_reservations[scratch_root.absolute()] = payload_bytes

    def _release_local_fastq_stage_for_workspace(self, scratch_root: Path) -> None:
        """Release the FASTQ bytes retained by one completed/aborted workspace."""

        with self._local_scratch_reservation_lock:
            payload_bytes = self._local_fastq_stage_reservations.pop(scratch_root.absolute(), 0)
        if payload_bytes:
            self._release_local_fastq_stage(payload_bytes)

    def _stage_local_fastq_inputs(
        self,
        work_dir: Path,
        scratch_root: Path,
        srr_id: str,
        expected_paired: bool | None = None,
    ) -> bool:
        """Copy validated FASTQs to local scratch for one quantification.

        The external raw inputs remain authoritative and are never moved. A
        staged copy is promoted only after all files have been copied with
        matching sizes; any failure falls back to the normal mounted-volume
        symlink and leaves the external inputs intact.
        """

        if not self.local_fastq_stage or self.local_quant_scratch_dir is None:
            return False
        source_dir = work_dir / "getfastq" / srr_id
        candidates = _validated_local_fastq_inputs(
            source_dir,
            srr_id,
            expected_paired,
        )
        if not candidates:
            logger.debug("No validated FASTQs available to stage for %s", srr_id)
            return False
        try:
            payload_bytes = sum(path.stat().st_size for path in candidates)
        except OSError as exc:
            logger.warning("Unable to size FASTQs for local staging %s: %s", srr_id, exc)
            return False
        if not self._reserve_local_fastq_stage(payload_bytes):
            return False

        staged_root = scratch_root / "getfastq"
        staged_sample_dir = staged_root / srr_id
        staged_successfully = False
        try:
            staged_sample_dir.mkdir(parents=True, exist_ok=True)
            reused = 0
            for source in candidates:
                destination = staged_sample_dir / source.name
                source_stat = source.stat()
                try:
                    destination_stat = destination.stat()
                except OSError:
                    destination_stat = None
                if (
                    destination_stat is not None
                    and destination_stat.st_size == source_stat.st_size
                    and destination_stat.st_mtime_ns == source_stat.st_mtime_ns
                ):
                    reused += 1
                    continue
                temporary = destination.with_name(f".{destination.name}.{os.getpid()}.{threading.get_ident()}.copying")
                shutil.copy2(source, temporary)
                if temporary.stat().st_size != source_stat.st_size:
                    raise OSError(f"staged FASTQ size mismatch for {source.name}")
                os.replace(temporary, destination)
            logger.info(
                "Staged %d validated FASTQ file(s) locally for %s (%.2f GB; reused %d)",
                len(candidates),
                srr_id,
                payload_bytes / (1024**3),
                reused,
            )
            self._retain_local_fastq_stage_reservation(scratch_root, payload_bytes)
            staged_successfully = True
            return True
        except OSError as exc:
            logger.warning(
                "Local FASTQ staging failed for %s; using mounted inputs: %s",
                srr_id,
                exc,
            )
            shutil.rmtree(staged_root, ignore_errors=True)
            return False
        finally:
            if not staged_successfully:
                self._release_local_fastq_stage(payload_bytes)

    def _prepare_local_quant_workspace(
        self,
        work_dir: Path,
        species_name: str,
        srr_id: str,
        metadata_path: str,
        command: List[str],
        *,
        expected_paired: bool | None = None,
    ) -> tuple[List[str], Path] | None:
        """Prepare a local output workspace and optionally stage its inputs.

        Amalgkit resolves FASTQ inputs relative to ``--out_dir``.  A temporary
        local workspace therefore uses a directory symlink for ``metadata``
        and either a validated local FASTQ copy or a ``getfastq`` symlink.
        The completed sample is copied atomically back into the external work
        tree; no authoritative raw input is moved or deleted by this
        optimization.
        """

        if self.local_quant_scratch_dir is None:
            return None
        scratch_root = self.local_quant_scratch_dir / species_name / srr_id
        try:
            self.local_quant_scratch_dir.mkdir(parents=True, exist_ok=True)
            available_gb = shutil.disk_usage(self.local_quant_scratch_dir).free / (1024**3)
            if available_gb <= self.local_quant_scratch_reserve_gb:
                logger.warning(
                    "Local quant scratch has %.2f GiB free; reserve is %.2f GiB. "
                    "Using external quant output for %s/%s.",
                    available_gb,
                    self.local_quant_scratch_reserve_gb,
                    species_name,
                    srr_id,
                )
                return None
            if scratch_root.is_symlink():
                scratch_root.unlink()
            elif scratch_root.exists() and not scratch_root.is_dir():
                scratch_root.unlink()
            scratch_root.mkdir(parents=True, exist_ok=True)
            # Keep a completed staged FASTQ directory across an interrupted
            # quantification, but reset only transient metadata/output trees.
            # This makes a stopped run resume without copying the same large
            # external FASTQ back to the local SSD.
            for transient_name in ("metadata", "quant"):
                transient = scratch_root / transient_name
                if transient.is_symlink() or (transient.exists() and not transient.is_dir()):
                    transient.unlink()
                elif transient.is_dir():
                    shutil.rmtree(transient)
            staged_getfastq = scratch_root / "getfastq"
            if staged_getfastq.is_symlink() or (staged_getfastq.exists() and not staged_getfastq.is_dir()):
                staged_getfastq.unlink()
            (scratch_root / "metadata").symlink_to(work_dir / "metadata", target_is_directory=True)
            if not self._stage_local_fastq_inputs(
                work_dir,
                scratch_root,
                srr_id,
                expected_paired,
            ):
                # Amalgkit v0.16.60 rejects a symlinked input/output root. If
                # the local FASTQ copy cannot fit under the host reserve, use
                # the canonical external workspace instead of constructing a
                # local workspace that is guaranteed to fail validation.
                self._cleanup_local_quant_workspace(scratch_root)
                logger.info(
                    "Local quant scratch unavailable for %s/%s; using canonical mounted workspace",
                    species_name,
                    srr_id,
                )
                return None
            (scratch_root / "quant").mkdir()
            local_command = list(command)
            local_command[local_command.index("--out_dir") + 1] = str(scratch_root)
            local_command[local_command.index("--metadata") + 1] = str(
                scratch_root / "metadata" / Path(metadata_path).name
            )
            logger.debug("Using local quant scratch for %s/%s: %s", species_name, srr_id, scratch_root)
            return local_command, scratch_root
        except OSError as exc:
            logger.warning(
                "Unable to prepare local quant scratch for %s/%s; using external output: %s",
                species_name,
                srr_id,
                exc,
            )
            if scratch_root.is_dir() and not scratch_root.is_symlink():
                shutil.rmtree(scratch_root, ignore_errors=True)
            self._release_local_fastq_stage_for_workspace(scratch_root)
            return None

    def _promote_local_quant_output(self, work_dir: Path, srr_id: str, scratch_root: Path) -> bool:
        """Copy one validated local quant result into the canonical data-root path."""

        local_sample_dir = scratch_root / "quant" / srr_id
        if find_quantification_file(local_sample_dir, srr_id) is None:
            return False
        external_quant_root = work_dir / "quant"
        external_quant_root.mkdir(parents=True, exist_ok=True)
        external_sample_dir = external_quant_root / srr_id
        temporary = external_quant_root / f".{srr_id}.promote.{os.getpid()}.{threading.get_ident()}"
        try:
            if temporary.exists() or temporary.is_symlink():
                if temporary.is_dir() and not temporary.is_symlink():
                    shutil.rmtree(temporary)
                else:
                    temporary.unlink()
            shutil.copytree(local_sample_dir, temporary, symlinks=False)
            if find_quantification_file(temporary, srr_id) is None:
                shutil.rmtree(temporary, ignore_errors=True)
                return False
            if external_sample_dir.exists() or external_sample_dir.is_symlink():
                archive_dir = work_dir / "archive" / "quant_replaced" / srr_id / str(time.time_ns())
                archive_dir.parent.mkdir(parents=True, exist_ok=True)
                external_sample_dir.replace(archive_dir)
            temporary.replace(external_sample_dir)
            return True
        except OSError as exc:
            logger.error("Unable to promote local quantification for %s: %s", srr_id, exc)
            if temporary.exists() and temporary.is_dir() and not temporary.is_symlink():
                shutil.rmtree(temporary, ignore_errors=True)
            return False

    def _cleanup_local_quant_workspace(self, scratch_root: Path) -> None:
        """Remove only the exact per-sample local scratch directory."""

        cleanup_succeeded = False
        try:
            if scratch_root.is_symlink():
                scratch_root.unlink()
            elif scratch_root.exists():
                shutil.rmtree(scratch_root)
            cleanup_succeeded = True
        except OSError as exc:
            logger.warning("Unable to clean local quant scratch %s: %s", scratch_root, exc)
        if cleanup_succeeded:
            self._release_local_fastq_stage_for_workspace(scratch_root)

    def _register_reuse_context(
        self,
        species_name: str,
        source_config_path: Path,
        runtime_config_path: Path,
    ) -> SpeciesReuseContext | None:
        """Bind quantification reuse to the configuration bytes, hashed once per species.

        Only the checked-in YAML is content-bound. An alias runtime YAML has
        different bytes from the checked-in one, and outputs recorded under
        either keep the historical path-based check.
        """

        if runtime_config_path != source_config_path:
            return None
        try:
            context = build_species_reuse_context(species_name, runtime_config_path)
        except OSError as exc:
            logger.warning("Unable to hash runtime configuration for %s: %s", species_name, exc)
            return None
        contexts: dict[str, SpeciesReuseContext] | None = getattr(self, "_reuse_contexts", None)
        if contexts is None:
            contexts = {}
            self._reuse_contexts = contexts
        contexts[species_name] = context
        return context

    def _reuse_context(self, species_name: str) -> SpeciesReuseContext | None:
        return getattr(self, "_reuse_contexts", {}).get(species_name)

    def _reuse_kwargs(self, species_name: str) -> dict[str, Any]:
        context = self._reuse_context(species_name)
        if context is None:
            return {}
        return {"expected_config_sha256": context.config_sha256, "index_digests": context.index_digests}

    def is_quantified(self, species_name: str, srr_id: str) -> bool:
        """Check for a reusable quantification, including compatible version drift."""
        quant_dir = _species_work_dir(species_name) / "quant" / srr_id
        quant_file = find_quantification_file(quant_dir, srr_id)
        if not quant_dir.exists() or quant_file is None or not quant_file.is_file() or quant_file.stat().st_size <= 100:
            return False
        context = self._reuse_context(species_name)
        if context is not None and context.is_cohort:
            return context.classify(quant_dir, srr_id).get("status") in REUSABLE_QUANT_STATUSES
        return is_reusable_quantification(quant_dir, srr_id, **self._reuse_kwargs(species_name))

    def _reclaim_after_current_quantification(self, species_name: str, srr_id: str) -> None:
        """Reclaim one sample's raw inputs only after current evidence exists."""

        if not self.reclaim_raw_after_quant or not self.is_quantified(species_name, srr_id):
            return
        cleanup = reclaim_sample_raw_inputs(_species_work_dir(species_name), srr_id)
        if cleanup["errors"]:
            logger.warning(
                "Raw-input reclamation for %s/%s completed with %d errors: %s",
                species_name,
                srr_id,
                len(cleanup["errors"]),
                "; ".join(cleanup["errors"][:3]),
            )
        elif cleanup["files_deleted"]:
            logger.info(
                "Reclaimed %.2f GB of raw inputs for current quantification %s/%s",
                cleanup["bytes_freed"] / (1024**3),
                species_name,
                srr_id,
            )

    def _reconcile_legacy_species(self, species_name: str, config_path: Path, quant_dir: Path) -> Dict[str, Any]:
        """One reconcile pass over every recorded run of a species (non-cohort discovery).

        Replaces the former ``reconcile`` + ``_reset_non_current_quantified_states``
        pair (two classifications of every quantified output and one commit per
        row): the quant directory is listed once, unchanged outputs whose
        sidecar fingerprint and scope match their audit are trusted with one
        ``lstat``, new or changed ones are classified, and every change lands in
        one transaction. Reusable outputs (current and version drift) become
        ``quantified`` unless the requantification policy requeues them;
        recorded ``quantified`` rows whose output is not reusable become
        ``quarantined``.
        """

        reuse_kwargs = self._reuse_kwargs(species_name)
        verify_content = _resume_verify_hashes()
        context = self._reuse_context(species_name)
        try:
            config_sha256 = context.config_sha256 if context is not None else sha256_file(config_path)
        except OSError as exc:
            logger.warning("Cannot hash %s for %s; re-judging every output: %s", config_path, species_name, exc)
            config_sha256 = None
        manifest = quant_dir.parent / "reference" / REFERENCE_MANIFEST_NAME
        try:
            reference_sha256 = sha256_file(manifest)
        except OSError:
            reference_sha256 = "absent"
        scope = (
            ReconcileScope(
                config_sha256=config_sha256,
                reference_sha256=reference_sha256,
                scope=_reconcile_scope_label("legacy:" + ("content" if context is not None else "path")),
            )
            if config_sha256 is not None
            else None
        )
        policy = getattr(self, "requantification_policy", "preserve")
        requeue = (
            REUSABLE_QUANT_STATUSES
            if policy == "all"
            else frozenset({QUANT_STATUS_VERSION_DRIFT}) if policy == "version-drift" else frozenset()
        )
        recorded = [srr_id for srr_id, _ in self.db.states_for_species(species_name)]
        profile = getattr(self, "_resource_profile", None)
        return self.db.reconcile_outputs(
            species_name,
            quant_dir,
            recorded,
            lambda path, run: classify_quantification(path, run, verify_content=verify_content, **reuse_kwargs),
            scope=scope,
            requantify_statuses=requeue,
            deep=verify_content,
            workers=max(1, getattr(profile, "validation_slots", 1)),
        )

    def _reset_non_current_quantified_states(self, species_name: str) -> int:
        """Audit quantified outputs without invalidating runtime-version drift."""

        sample_ids = self.db.get_samples(species_name, "quantified")
        if not sample_ids:
            return 0

        profile = getattr(self, "_resource_profile", None)
        validation_slots = max(
            1,
            min(
                len(sample_ids),
                getattr(profile, "validation_slots", 1),
            ),
        )

        # Re-digesting every quantified payload on every resume made the
        # Phase-1 audit O(campaign bytes): a large species (apis_mellifera,
        # ~6.7 GB across 3190 sample dirs) stalled discovery for 15+ minutes
        # on cold external-volume reads, so no task was ever submitted.
        # The contract audit (sidecar, accession, config, reference manifest)
        # remains fail-closed; full content re-verification is opt-in.
        verify_resume_hashes = os.environ.get("AMALGKIT_RESUME_VERIFY_HASHES", "").strip().lower() in {
            "1",
            "true",
            "yes",
            "on",
        }

        reuse_kwargs = self._reuse_kwargs(species_name)

        def classify(srr_id: str) -> tuple[str, dict[str, Any]]:
            quant_dir = _species_work_dir(species_name) / "quant" / srr_id
            return (
                srr_id,
                classify_quantification(quant_dir, srr_id, verify_content=verify_resume_hashes, **reuse_kwargs),
            )

        if validation_slots == 1:
            classifications = [classify(srr_id) for srr_id in sample_ids]
        else:
            with concurrent.futures.ThreadPoolExecutor(max_workers=validation_slots) as validation_executor:
                classifications = list(validation_executor.map(classify, sample_ids))

        reset = 0
        requantification_candidates = 0
        quarantined = 0
        for srr_id, audit in classifications:
            status = audit["status"]
            if status in {QUANT_STATUS_CURRENT, QUANT_STATUS_VERSION_DRIFT}:
                should_requantify = self.requantification_policy == "all" or (
                    self.requantification_policy == "version-drift" and status == QUANT_STATUS_VERSION_DRIFT
                )
                self.db.record_quantification_audit(
                    species_name,
                    srr_id,
                    status=status,
                    reason=audit["reason"],
                    contract_id=audit.get("contract_id"),
                    observed_amalgkit_version=audit.get("observed_amalgkit_version"),
                    observed_release_tag=audit.get("observed_release_tag"),
                    observed_source_revision=audit.get("observed_source_revision"),
                    state="pending" if should_requantify else "quantified",
                )
                if should_requantify:
                    reset += 1
                    requantification_candidates += 1
            else:
                self.db.record_quantification_audit(
                    species_name,
                    srr_id,
                    status=status,
                    reason=audit["reason"],
                    contract_id=audit.get("contract_id"),
                    observed_amalgkit_version=audit.get("observed_amalgkit_version"),
                    observed_release_tag=audit.get("observed_release_tag"),
                    observed_source_revision=audit.get("observed_source_revision"),
                    state="quarantined",
                )
                reset += 1
                quarantined += 1
        if len(sample_ids) > 1:
            logger.info(
                "Validated %d recorded quantification states for %s with %d provenance-read slot(s); " "policy=%s",
                len(sample_ids),
                species_name,
                validation_slots,
                self.requantification_policy,
            )
        if requantification_candidates:
            logger.info(
                "Queued %d explicit re-quantification candidates for %s",
                requantification_candidates,
                species_name,
            )
        if quarantined:
            logger.info(
                "Quarantined %d incomplete or unverifiable quantification outputs for %s",
                quarantined,
                species_name,
            )
        return reset

    def process_single_sample(
        self,
        srr_id: str,
        batch_idx: int,
        fastq_dir: Path,
        config_path: Path,
        species_name: str,
        threads: int,
        expected_paired: bool | None = None,
        *,
        defer_ncbi_fallback: bool = False,
        acquisition_only: bool = False,
        quantification_only: bool = False,
        metadata_path: Path | None = None,
        requantify: bool = False,
    ) -> Dict[str, Any]:
        """Processing unit for a single sample (``metadata_path``: frozen-cohort task table).

        ``requantify`` (a row the requantification policy requeued) bypasses
        the reuse shortcut so the sample is quantified again; any other
        sample with a reusable output on disk is recorded ``quantified`` and
        skipped.
        """
        if acquisition_only and quantification_only:
            raise ValueError("acquisition_only and quantification_only are mutually exclusive")
        result = {
            "srr": srr_id,
            "batch": batch_idx,
            "downloaded": False,
            "quantified": False,
            "skipped": False,
            "fallback_deferred": False,
            "error": None,
        }

        if not requantify and self.is_quantified(species_name, srr_id):
            # The row may be pending (a reconcile that has not seen this output yet): record what is on disk,
            # or the run summary counts it as remaining work and the producer never drains.
            if self.db.get_state(species_name, srr_id) != "quantified":
                self.db.set_state(species_name, srr_id, "quantified")
            self._reclaim_after_current_quantification(species_name, srr_id)
            result.update({"downloaded": True, "quantified": True, "skipped": True})
            return result

        # Download
        if not quantification_only:
            self.db.set_state(species_name, srr_id, "downloading")
        download_outcome: bool | None  # None: NCBI fallback queued after ENA failure
        try:
            if quantification_only:
                deadline = time.monotonic() + self.quant_timeout_seconds if self.quant_timeout_seconds else None
                download_outcome = bool(
                    self._validate_acquisition_fastqs(
                        fastq_dir / srr_id,
                        srr_id,
                        expected_paired,
                        deadline,
                        require_records=True,
                        fresh_record_scan=True,
                    )
                )
            else:
                download_outcome = self.download_fastq(
                    srr_id,
                    fastq_dir,
                    expected_paired,
                    defer_ncbi_fallback=defer_ncbi_fallback,
                )
        except (FirstPassDeferred, StagingDeadline) as exc:
            result.update(
                status="deferred_download",
                download_deferred=True,
                # A first-pass deferral waits for a later pass by design; a staging deadline is only a full queue.
                deferral="first_pass" if isinstance(exc, FirstPassDeferred) else "staging_deadline",
                error=str(exc),
                acquisition_budget_seconds=self.download_first_pass_seconds,
                deferred_at_epoch=time.time(),
            )
            self.db.set_state(species_name, srr_id, "pending", error=str(exc))
            return result
        if download_outcome is None:
            self.db.set_state(
                species_name,
                srr_id,
                "pending",
                error="NCBI fallback queued after ENA acquisition failure",
            )
            result["fallback_deferred"] = True
            return result
        if not download_outcome:
            error_msg = (
                "Acquired FASTQs missing or invalid; no network attempted in quantification stage"
                if quantification_only
                else "Download Failed (all sources: ENA FTP/HTTP, NCBI)" + self._layout_failure_suffix(srr_id)
            )
            self.db.set_state(species_name, srr_id, "failed", error=error_msg)
            result["error"] = error_msg
            return result
        self.db.set_state(species_name, srr_id, "downloaded")
        result["downloaded"] = True
        # Measured before quantification reclaims the reads: the run's FASTQ GB/hour.
        result["input_bytes"] = _acquired_input_bytes(fastq_dir / srr_id, srr_id)
        if acquisition_only:
            result["ready_for_quant"] = True
            return result

        # Quantify
        self.db.set_state(species_name, srr_id, "quantifying")
        quant_kwargs: Dict[str, Any] = {"metadata_path": metadata_path} if metadata_path is not None else {}
        if self._quant_semaphore is None:
            success, quant_error = self.quant_sample(
                config_path,
                batch_idx,
                species_name,
                threads,
                srr_id,
                expected_paired,
                **quant_kwargs,
            )
        else:
            with self._quant_semaphore:
                success, quant_error = self.quant_sample(
                    config_path,
                    batch_idx,
                    species_name,
                    threads,
                    srr_id,
                    expected_paired,
                    **quant_kwargs,
                )
        if success:
            self.db.set_state(species_name, srr_id, "quantified")
            self._reclaim_after_current_quantification(species_name, srr_id)
        else:
            error_msg = quant_error or "Quantification Failed"
            self.db.set_state(species_name, srr_id, "failed", error=error_msg)
        result["quantified"] = success
        if not success:
            result["error"] = quant_error or "Quantification Failed"

        return result

    def verify_genome_index(self, config_path: Path, species_name: str) -> bool:
        """Verify Kallisto index exists."""
        try:
            with open(config_path) as f:
                cfg = yaml.safe_load(f)
        except Exception as e:
            logger.error(f"Failed to load config {config_path}: {e}")
            return False

        logger.info(f"Verifying index for {species_name}...")
        for d in _quant_index_candidates(cfg, species_name):
            if not d:
                continue

            path_obj = Path(d)
            abs_path = path_obj.resolve() if path_obj.exists() else path_obj.absolute()

            logger.info(f"  Checking {d} -> {abs_path}")

            if not path_obj.exists():
                logger.warning(f"  Path does not exist: {d}")
                continue

            # Check permissions by attempting listdir
            try:
                files = [
                    candidate
                    for candidate in path_obj.glob("*.idx")
                    if candidate.is_file() and candidate.stat().st_size > 0
                ]
                if files:
                    logger.info(f"  Genome index found in {d}: {[f.name for f in files]}")
                    return True
                else:
                    logger.warning(f"  Path exists but no .idx files found in {d}")
                    try:
                        contents = os.listdir(d)
                        logger.info(f"  Directory contents: {contents[:5]}...")
                    except Exception as e:
                        logger.error(f"  Failed to list directory {d}: {e}")
            except Exception as e:
                logger.error(f"  Permission/Error accessing {d}: {e}")

        logger.error(f"No genome index found for {species_name}")
        return False

    def run_tissue_normalization(self, metadata_path: Path) -> bool:
        """Run tissue normalization atomically and preserve source labels.

        The orchestrator may revisit the same selected table on every resume.
        ``tissue_original`` therefore remains the first observed value instead
        of being overwritten by an already-normalized value. A failed write or
        mapping operation is reported to the caller so discovery can fail
        closed rather than proceeding with an ambiguous metadata state.
        """
        mapping_path = self.config_dir / "tissue_mapping.yaml"
        patches_path = self.config_dir / "tissue_patches.yaml"

        if not mapping_path.exists():
            logger.warning(f"Tissue mapping not found at {mapping_path}, skipping normalization.")
            return True

        logger.info(f"Normalizing tissues in {metadata_path}")
        temporary: Path | None = None
        try:
            import pandas as pd

            df = pd.read_csv(metadata_path, sep="\t", low_memory=False)
            if "tissue" not in df.columns:
                logger.warning("Metadata has no tissue column; skipping tissue normalization.")
                return True
            original_tissue = df["tissue_original"].copy() if "tissue_original" in df.columns else df["tissue"].copy()

            # Apply normalization
            df_norm = apply_tissue_normalization(
                df,
                mapping_path=mapping_path,
                patches_path=patches_path if patches_path.exists() else None,
                output_column="tissue_normalized",
            )

            # Keep the normalized tissue in the canonical field while retaining
            # the source value for provenance.
            if "tissue_normalized" in df_norm.columns:
                df_norm["tissue_original"] = original_tissue
                df_norm["tissue"] = df_norm["tissue_normalized"]
                # Drop temporary column if desired, or keep it

            temporary = metadata_path.with_name(f".{metadata_path.name}.{os.getpid()}.{threading.get_ident()}.tmp")
            df_norm.to_csv(temporary, sep="\t", index=False)
            # Avoid changing the source table when normalization is already
            # idempotent; otherwise promote the complete file atomically.
            if not metadata_path.exists() or temporary.read_bytes() != metadata_path.read_bytes():
                os.replace(temporary, metadata_path)
            else:
                temporary.unlink(missing_ok=True)
            logger.info("Tissue normalization complete and saved.")
            return True

        except Exception as e:
            logger.error(f"Tissue normalization failed: {e}")
            if temporary is not None:
                temporary.unlink(missing_ok=True)
            return False

    def _discover_cohort_species_tasks(self, config_name: str, max_gb: float) -> List[Dict[str, Any]]:
        """Build one species' tasks from the frozen cohort, never from a metadata query.

        Reads the cohort (already hash-verified), checks that the checked-in
        configuration and the data root's reference index bytes are the
        cohort's, reconciles outputs on disk against the cohort reuse contract,
        and returns the schedulable tasks: cohort tasks minus ``permanent_drop``
        exclusions (cohort and progress DB) minus quantified and quarantined
        runs. It runs no Amalgkit metadata/select or tissue normalization,
        never writes ``metadata_selected.tsv`` or the metadata witness, and
        prunes nothing. ``amalgkit quant --metadata`` receives a separate task
        table (frozen rows of the schedulable tasks, written atomically only
        when its bytes change), and each task's batch index is its row in that
        table. Any inconsistency raises, which fails the whole discovery phase.
        """

        cohort = cast(FrozenCohort, self.cohort)
        species_name = cohort.species_for_config(config_name)
        entry = cohort.entry(species_name)
        source_config_path = self.config_dir / config_name
        if source_config_path.is_symlink() or not source_config_path.is_file():
            raise CohortError(f"cohort configuration is missing: {source_config_path}")
        verify_cohort_configs(cohort, self.config_dir, [species_name])
        work_dir = _species_work_dir(species_name)
        if work_dir.parent.name != species_name or work_dir.parent.parent != _data_root():
            raise CohortError(
                f"{species_name}: cohort mode needs the canonical work directory {_data_root() / species_name}, "
                f"not the alias {work_dir}"
            )
        logger.info("=== Discovering %s from frozen cohort %s ===", species_name, cohort.cohort_id)

        try:
            runtime_cfg = yaml.safe_load(source_config_path.read_text()) or {}
        except (OSError, yaml.YAMLError) as exc:
            raise CohortError(f"{species_name}: configuration is unreadable: {exc}") from exc
        index_dir = _resolve_index_dir(runtime_cfg, species_name)
        if not index_dir:
            raise CohortError(f"{species_name}: no kallisto index directory resolves from the configuration")
        reference = verify_species_reference(
            cohort, _data_root(), species_name, index_digests=self._cohort_digests, index_dir=index_dir
        )
        if not reference.ok:
            raise CohortError(f"{species_name}: {reference.reason} ({reference.manifest_relpath})")

        context = build_species_reuse_context(
            species_name,
            source_config_path,
            cohort_id=cohort.cohort_id,
            reference_index_sha256=entry.index_sha256,
            reference_index_name=entry.index_name,
            index_digests=self._cohort_digests,
        )
        contexts: dict[str, SpeciesReuseContext] | None = getattr(self, "_reuse_contexts", None)
        if contexts is None:
            contexts = {}
            self._reuse_contexts = contexts
        contexts[species_name] = context

        tasks = cohort.tasks_for(species_name)
        self.db.init_species(species_name, [task.accession for task in tasks])
        self.db.ensure_exclusions(
            species_name,
            [
                {
                    "srr_id": accession,
                    "reason_code": reason,
                    "reason_detail": f"frozen cohort {cohort.cohort_id}",
                    "recorded_by": f"cohort:{cohort.cohort_id}",
                }
                for accession, reason in sorted(cohort.exclusions_for(species_name).items())
            ],
        )
        dropped = self.db.get_excluded_srr_ids(species_name, reason_code="permanent_drop")
        schedulable = [task for task in tasks if task.accession not in dropped]

        policy = getattr(self, "requantification_policy", "preserve")
        requeue = (
            REUSABLE_QUANT_STATUSES
            if policy == "all"
            else frozenset({QUANT_STATUS_VERSION_DRIFT}) if policy == "version-drift" else frozenset()
        )
        verify_content = _resume_verify_hashes()
        profile = getattr(self, "_resource_profile", None)
        reconciled = self.db.reconcile_outputs(
            species_name,
            work_dir / "quant",
            [task.accession for task in schedulable],
            lambda path, run: context.classify(path, run, verify_content=verify_content),
            scope=ReconcileScope(
                config_sha256=context.config_sha256,
                reference_sha256=entry.index_sha256,
                scope=_reconcile_scope_label(f"cohort:{cohort.cohort_id}"),
            ),
            requantify_statuses=requeue,
            deep=verify_content,
            workers=max(1, getattr(profile, "validation_slots", 1)),
        )

        metadata_path = work_dir / "metadata" / cohort_task_table_name(cohort.cohort_id)
        positions = write_cohort_task_table(cohort, species_name, metadata_path, excluded=dropped)
        done = set(self.db.get_samples(species_name, "quantified")) | set(
            self.db.get_samples(species_name, "quarantined")
        )
        max_bases = max_gb * 1e9
        deferred = [task for task in schedulable if task.accession not in done and task.total_bases > max_bases]
        requantify = set(reconciled.get("requeued_accessions") or ())
        self._record_discovery(
            species_name,
            accessions=[task.accession for task in schedulable],
            deferred_max_gb=[task.accession for task in deferred],
            unavailable=list(reconciled.get("unavailable") or ()),
        )
        fastq_dir = work_dir / "getfastq"
        built = [
            {
                "srr": task.accession,
                "batch_idx": positions[task.accession],
                "fastq_dir": fastq_dir,
                "config_path": source_config_path,
                "species_name": species_name,
                "expected_paired": task.expected_paired,
                "metadata_path": metadata_path,
                "cohort_id": cohort.cohort_id,
                "total_bases": task.total_bases,
                "fastq_bytes": task.fastq_bytes,
                **({"requantify": True} if task.accession in requantify else {}),
            }
            for task in sorted(schedulable, key=lambda task: (task.total_bases, task.batch_index))
            if task.accession not in done and task.total_bases <= max_bases
        ]
        logger.info(
            "%s: cohort %s has %d tasks; %d permanently excluded, %d reusable on disk, %d quarantined, "
            "%d above %.0f GB deferred, %d scheduled (task table %s)",
            species_name,
            cohort.cohort_id,
            len(tasks),
            len(tasks) - len(schedulable),
            reconciled["quantified"],
            len(reconciled["quarantined"]),
            len(deferred),
            max_gb,
            len(built),
            metadata_path.name,
        )
        return built

    def _record_discovery(self, species_name: str, **fields: Sequence[str]) -> None:
        """Remember what one species' discovery found for the run summary.

        ``accessions``: cohort tasks after exclusions; ``deferred_max_gb``:
        tasks above ``max_gb``; ``unavailable``: outputs reconcile could not
        read (left unjudged).
        """

        lock = self.__dict__.setdefault("_discovery_record_lock", threading.Lock())
        with lock:
            record = self.__dict__.setdefault("_cohort_discovered", {}).setdefault(species_name, {})
            for name, values in fields.items():
                record[name] = tuple(sorted(values))

    def discover_species_tasks(self, config_name: str, max_gb: float, threads: int) -> List[Dict[str, Any]]:
        """Main processing loop for a species.
        Returns a list of tasks instead of executing them immediately.
        """
        if getattr(self, "cohort", None) is not None:
            return self._discover_cohort_species_tasks(config_name, max_gb)
        source_config_path = self.config_dir / config_name
        if not source_config_path.exists():
            logger.error(f"Config not found: {source_config_path}")
            return []

        species_name = _species_name_from_config(config_name)
        config_path = _runtime_config_path(source_config_path, species_name)
        self._register_reuse_context(species_name, source_config_path, config_path)
        logger.info(f"=== Processing {species_name} ===")

        # Autonomous Preprocessing Detection
        work_dir = _species_work_dir(species_name)
        metadata_path = _resolve_metadata_path(work_dir)

        needs_metadata = not metadata_path.exists() or not is_current_metadata(work_dir)
        if needs_metadata and not getattr(self, "allow_metadata_refresh", True):
            raise RuntimeError(
                f"Metadata for {species_name} is missing or not current and metadata refresh is disabled; "
                "refusing to query NCBI (run with a frozen cohort, or pass --refresh-metadata explicitly)"
            )
        needs_index = not self.verify_genome_index(config_path, species_name)
        needs_prep = needs_metadata or needs_index
        if needs_metadata:
            logger.info(f"Metadata missing for {species_name}. Queuing autonomous preprocessing...")
        if needs_index:
            logger.info(f"Genome index missing for {species_name}. Queuing autonomous preprocessing...")

        if needs_prep:
            logger.info(
                "Running current Amalgkit preparation stages (reference -> metadata -> select) for %s...",
                species_name,
            )

            # Amalgkit 0.16.60 builds the quantification index through the
            # reference-preparation helper; it no longer exposes the removed
            # `config` or `index` CLI stages.
            if needs_index:
                try:
                    from metainformant.rna.engine.workflow import load_workflow_config
                    from metainformant.rna.engine.workflow_planning import (
                        prepare_reference_genome,
                    )

                    runtime_config = load_workflow_config(config_path)
                    if not prepare_reference_genome(runtime_config):
                        logger.error("Reference preparation failed for %s.", species_name)
                        return []
                except Exception as exc:
                    logger.error("Reference preparation exception for %s: %s", species_name, exc)
                    return []

            # Taxonomy acquisition is owned by Amalgkit 0.16.60.  When the
            # campaign launcher provides AMALGKIT_SHARED_DOWNLOAD_DIR, its
            # native taxdump/database cache is reused across species.  The
            # older per-species seed path is intentionally not recreated.
            # Without the shared setting, Amalgkit's own downloader remains
            # responsible for creating its inferred download directory.
            try:
                from metainformant.rna.amalgkit.amalgkit import build_amalgkit_command
                from metainformant.rna.engine.workflow import load_workflow_config
                from metainformant.rna.engine.workflow_planning import (
                    apply_step_defaults,
                    plan_workflow,
                    sanitize_params_for_cli,
                )

                runtime_config = apply_step_defaults(load_workflow_config(config_path))
                planned_steps = dict(plan_workflow(runtime_config))
                # Re-query NCBI only when the metadata itself is missing or stale: an index-only
                # preparation must not silently change the selected run universe.
                for prep_step in ("metadata", "select") if needs_metadata else ():
                    raw_params = planned_steps.get(prep_step)
                    if raw_params is None:
                        logger.error("Current workflow plan omitted %s for %s", prep_step, species_name)
                        return []
                    prep_params = sanitize_params_for_cli(prep_step, dict(raw_params))
                    prep_params["redo"] = "yes"
                    prep_cmd = build_amalgkit_command(prep_step, prep_params)
                    prep_result = subprocess.run(
                        prep_cmd,
                        capture_output=True,
                        text=True,
                        timeout=7200,
                        check=False,
                    )
                    if prep_result.returncode != 0:
                        logger.error(
                            "%s preparation failed for %s:\n%s\n%s",
                            prep_step,
                            species_name,
                            prep_result.stderr,
                            prep_result.stdout,
                        )
                        return []
                logger.info("Metadata/select preparation complete for %s.", species_name)
            except subprocess.TimeoutExpired:
                logger.error("Metadata/select preparation timed out after 2h for %s.", species_name)
                return []
            except Exception as exc:
                logger.error("Metadata/select preparation exception for %s: %s", species_name, exc)
                return []

            # Post-prep explicit verification
            if needs_index and not self.verify_genome_index(config_path, species_name):
                logger.error(f"Genome index still missing for {species_name} after preparation.")
                return []

            metadata_path = _resolve_metadata_path(work_dir)

            if needs_metadata and not metadata_path.exists():
                logger.error(f"No metadata generated for {species_name} after preparation.")
                return []
        # Apply Normalization Step
        if not self.run_tissue_normalization(metadata_path):
            logger.error("Metadata normalization did not complete for %s", species_name)
            return []

        # Some valid Amalgkit 0.16.60 runs emit only ``metadata.tsv`` when
        # selection does not create a separate table.  Materialize the
        # canonical selected-table boundary once so the current metadata
        # witness can remain strict and every later resume can distinguish
        # the source table from the exact acquisition/quantification input.
        selected_metadata_path = work_dir / "metadata" / "metadata_selected.tsv"
        if not selected_metadata_path.is_file():
            try:
                shutil.copy2(metadata_path, selected_metadata_path)
                metadata_path = selected_metadata_path
                logger.info(
                    "Materialized canonical selected metadata for %s: %s",
                    species_name,
                    selected_metadata_path,
                )
            except OSError as exc:
                logger.error(
                    "Unable to materialize selected metadata for %s: %s",
                    species_name,
                    exc,
                )
                return []

        # Load samples
        import pandas as pd

        df = pd.read_csv(metadata_path, sep="\t", low_memory=False)

        filtered = _filter_metadata_by_size(df, max_gb)

        logger.info(f"Samples: {len(df)} total -> {len(filtered)} filtered (<= {max_gb} GB)")

        if filtered.empty:
            logger.info("No samples to process.")
            return []

        before_accession_filter = len(filtered)
        filtered = _filter_valid_run_rows(filtered, species_name, metadata_path)
        if len(filtered) != before_accession_filter:
            logger.info(
                "%s: %d rows remain after concrete run-accession validation",
                species_name,
                len(filtered),
            )
        if filtered.empty:
            logger.warning("No downloadable run accessions remain for %s.", species_name)
            return []

        try:
            runtime_cfg = yaml.safe_load(config_path.read_text()) or {}
        except (OSError, yaml.YAMLError) as exc:
            logger.error("Unable to load runtime config for reference preflight %s: %s", species_name, exc)
            return []
        target_column = next(
            (column for column in ("scientific_name", "organism", "species") if column in filtered.columns),
            None,
        )
        if target_column:
            target_names = sorted(
                str(value).strip() for value in filtered[target_column].dropna().unique() if str(value).strip()
            )
        else:
            target_names = [str(value).strip() for value in runtime_cfg.get("species_list", []) if str(value).strip()]
        references_ready, _, _ = _ensure_reference_alias_indexes(runtime_cfg, species_name, target_names)
        if not references_ready:
            return []

        # Kallisto cannot pseudoalign reads shorter than the active index
        # k-mer.  Exclude only rows with a known incompatible length before
        # acquisition, preserve the complete Amalgkit selection table for
        # auditability, and leave unknown lengths for Amalgkit's own
        # quant-time recovery from downloaded-read statistics.
        index_dir_value = _resolve_index_dir(runtime_cfg, species_name)
        index_dir = Path(index_dir_value) if index_dir_value else None
        reference_index = _reference_index_for_kmer(index_dir, target_names)
        kmer_size = _kallisto_index_kmer(str(reference_index)) if reference_index else KALLISTO_DEFAULT_KMER_SIZE
        quant_eligible, quant_excluded, _ = _filter_kallisto_ineligible_reads(filtered, kmer_size)
        if not quant_excluded.empty:
            selected_snapshot = metadata_path.with_name("metadata_selected_all.tsv")
            try:
                # Refresh the complete-selection snapshot whenever metadata
                # was regenerated.  On a resume without metadata work,
                # preserve the prior full snapshot while the canonical
                # selected table contains only quant-eligible rows.
                if needs_metadata or not selected_snapshot.exists():
                    shutil.copy2(metadata_path, selected_snapshot)
                exclusion_audit = _write_quant_eligibility_audit(
                    work_dir,
                    quant_excluded,
                    kmer_size=kmer_size,
                    index_path=reference_index,
                )
            except OSError as exc:
                logger.error("Unable to preserve quantification eligibility audit for %s: %s", species_name, exc)
                return []
            logger.warning(
                "%s: excluding %d samples before quantification because their known mean read length is below "
                "the Kallisto k-mer size (%d); see %s",
                species_name,
                len(quant_excluded),
                kmer_size,
                exclusion_audit,
            )
            filtered = quant_eligible

        # Register samples in progress DB and reconcile existing results
        srr_col = _sample_run_column(filtered)
        all_srr_ids = filtered[srr_col].tolist()
        self.db.init_species(species_name, all_srr_ids)
        self.db.prune_species(species_name, all_srr_ids)

        quant_dir = _species_work_dir(species_name) / "quant"
        legacy_reconciled = self._reconcile_legacy_species(species_name, config_path, quant_dir)
        requantify = set(legacy_reconciled.get("requeued_accessions") or ())
        self._record_discovery(species_name, unavailable=list(legacy_reconciled.get("unavailable") or ()))

        # Mark all filtered samples as sampled so amalgkit quant processes them
        filtered["is_sampled"] = "yes"

        # Write sorted metadata for batch processing
        filtered.to_csv(metadata_path, sep="\t", index=False)

        # Write the metadata witness only after the final selected table has
        # been written.  The selected table is the acquisition/quantification
        # input and may be filtered for invalid accessions or incompatible
        # read lengths above; hashing it before this write made every restart
        # appear stale and needlessly repeated metadata discovery.
        try:
            rules_path = self.config_dir / "select_rules.tsv"
            write_metadata_provenance(
                work_dir,
                species=species_name,
                config_path=config_path,
                selection_rules_path=rules_path,
            )
        except OSError as exc:
            logger.error("Unable to write metadata provenance for %s: %s", species_name, exc)
            return []

        fastq_dir = work_dir / "getfastq"
        # ``_reconcile_legacy_species`` has already established the current
        # filesystem state in the progress database and quarantined stale
        # quantified states. Re-read that authoritative set
        # once instead of calling ``is_quantified`` for every accession. On
        # large external-volume cohorts, thousands of repeated directory/stat
        # calls here can delay execution before the first reusable sample is
        # even submitted.
        quantified_accessions = set(self.db.get_samples(species_name, "quantified"))
        permanently_excluded = self.db.get_excluded_srr_ids(species_name, reason_code="permanent_drop")
        if permanently_excluded:
            logger.info(
                "Skipping %d permanently excluded accessions for %s per sample_exclusions",
                len(permanently_excluded),
                species_name,
            )
        tasks = [
            task
            for task in _build_sample_tasks(filtered, srr_col, fastq_dir, config_path, species_name)
            if str(task["srr"]).strip() not in quantified_accessions
            and str(task["srr"]).strip() not in permanently_excluded
        ]
        for task in tasks:
            if str(task["srr"]).strip() in requantify:
                task["requantify"] = True

        return tasks

    def run_all(
        self,
        species_list: List[str],
        max_gb: float,
        workers: Optional[int] = None,
        threads: Optional[int] = None,
        *,
        quant_slots: Optional[int] = None,
        fasterq_threads: Optional[int] = None,
        fasterq_slots: Optional[int] = None,
        compression_threads: Optional[int] = None,
        max_in_flight: Optional[int] = None,
        validation_slots: Optional[int] = None,
        requantification_policy: Optional[str] = None,
        discovery_workers: Optional[int] = None,
        reuse_dry_run_fraction: Optional[float] = None,
        resource_profile: Optional[PipelineResourceProfile] = None,
        scheduling_order: Optional[str] = None,
    ) -> Dict[str, Any]:
        """Run the current pipeline with bounded, explicit stage resources; return a run summary.

        The whole run holds the campaign writer lock of the data root
        (:mod:`metainformant.rna.engine.campaign_lock`, kind ``producer``);
        another writer makes this raise :class:`CampaignLockBusy` before any
        state is touched. Under the lock the preflight runs first (a broken
        host changes nothing); then startup resets every in-flight row
        (``downloading``/``downloaded``/``quantifying``) and every
        environment-class failure to ``pending`` (``permanent_drop`` rows
        stay excluded) and, when ``reuse_dry_run_fraction`` (or
        ``AMALGKIT_PREFLIGHT_REUSE_DRY_RUN``) is set, runs the reuse dry-run.

        Resources: ``resource_profile`` if given; else explicit ``workers``
        and ``threads`` build the historical profile
        (:func:`build_pipeline_resource_profile`); else the profile is derived
        from this host (:func:`build_host_resource_profile`: measured c19
        defaults, which also set the acquisition slots and the ENA
        connection budgets for the run). The resolved profile is logged.

        Discovery streams: each species' tasks join the queue as soon as its
        own discovery finishes, so execution starts before the slowest
        species is discovered. The unsubmitted queue is re-ordered on every
        arrival (:func:`order_campaign_tasks`): reusable local inputs first,
        then resumable partials, then fresh acquisitions in
        ``scheduling_order`` (default global size bands with species
        round-robin), ties broken by species and run. A species whose
        discovery raises is skipped for this run and the exit code becomes
        ``EXIT_RETRYABLE`` (75); the other species proceed. The sample
        executor uses a bounded submission window, a raw-byte reservation
        when the profile has one, and its quantification budget stays
        separate from the fallback ``fasterq-dump`` and ``pigz`` budgets.

        Retry bookkeeping (:class:`~metainformant.rna.engine.progress_db.RetryPolicy`):
        each finished attempt is classified. Environment failures are
        retried at the next start; transient failures back off
        exponentially (re-queued in this run when the backoff elapses,
        skipped by later starts until then); layout mismatches are
        quarantined and kallisto exit-1 failures and spent attempt budgets
        are held for the operator
        (:meth:`~metainformant.rna.engine.progress_db.ProgressDB.release_retry_holds`).
        Held rows do not count as remaining work.

        ``AMALGKIT_ENV_FAILURE_TRIP`` (default 5) consecutive environment-class
        failures (EPERM, EACCES, EROFS, ENXIO, ENOSPC) stop admission; those
        rows stay ``pending`` and the summary's ``exit_code`` is
        ``EXIT_ENVIRONMENT`` (69). Otherwise ``exit_code`` is ``EXIT_OK`` (0)
        when no scheduled task is left retryable and ``EXIT_RETRYABLE`` (75)
        when the same command should be run again (``next_attempt_due_at`` in
        the summary says when backed-off work becomes due). The WAL is
        checkpointed at exit.
        """

        lock = acquire_campaign_lock(_data_root(), "producer")
        if not lock.inherited:
            holder = f"pid {os.getpid()}"
        elif lock.fd is not None:
            holder = f"already held by this process, pid {os.getpid()}"
        else:
            holder = f"inherited from the launching process, pid {lock.owner.get('pid')}"
        logger.info("Holding the campaign writer lock %s (%s)", lock.path, holder)
        # Every transfer/extraction/quant process group inherits the locked descriptor: the kernel keeps the
        # root locked while any of them lives, even if this process is SIGKILLed.
        process_groups.reset()
        process_groups.set_inherited_fds([lock.child_fd])
        stop = ProducerStop()
        self._producer_stop = stop
        try:
            with stop.installed():
                return self._run_all_locked(
                    species_list,
                    max_gb,
                    workers,
                    threads,
                    quant_slots=quant_slots,
                    fasterq_threads=fasterq_threads,
                    fasterq_slots=fasterq_slots,
                    compression_threads=compression_threads,
                    max_in_flight=max_in_flight,
                    validation_slots=validation_slots,
                    requantification_policy=requantification_policy,
                    discovery_workers=discovery_workers,
                    reuse_dry_run_fraction=reuse_dry_run_fraction,
                    resource_profile=resource_profile,
                    scheduling_order=scheduling_order,
                )
        finally:
            # The executors are joined inside _run_with_profile; never release the lock while a child group lives.
            if process_groups.registered():
                stop.request(0)
                while not process_groups.wait_until_empty(1.0):
                    stop.stop_children()
            process_groups.set_inherited_fds(())
            process_groups.reset()
            self._producer_stop = None
            try:
                self.db.checkpoint_wal()
            except Exception as exc:  # noqa: BLE001 - a dropped volume must not mask the run's own outcome
                logger.warning("Final WAL checkpoint failed: %s", exc)
            lock.release()

    def _preflight_targets(self, species_list: List[str]) -> tuple[Dict[str, Path], Dict[str, Path]]:
        """Per species: the quant directory (one run_info.json is compared with kallisto) and, in cohort mode,
        the frozen index whose header carries the kallisto index format."""

        quant_dirs: Dict[str, Path] = {}
        index_paths: Dict[str, Path] = {}
        cohort = getattr(self, "cohort", None)
        for config_name in species_list:
            try:
                if cohort is not None:
                    species_name = cohort.species_for_config(config_name)
                    entry = cohort.entry(species_name)
                    index_paths[species_name] = Path(cohort.root) / entry.index_relpath
                    quant_dirs[species_name] = _data_root() / species_name / "work" / "quant"
                else:
                    species_name = _species_name_from_config(config_name)
                    quant_dirs[species_name] = _species_work_dir(species_name) / "quant"
            except (CohortError, OSError, KeyError) as exc:
                logger.warning("Preflight skips %s: %s", config_name, exc)
        return quant_dirs, index_paths

    def _reuse_probes(self, species_list: List[str]) -> List[ReuseProbe]:
        """The reuse dry-run inputs: each species' recorded-quantified runs and the resume classifier."""

        probes: List[ReuseProbe] = []
        cohort = getattr(self, "cohort", None)
        for config_name in species_list:
            source = self.config_dir / config_name
            if cohort is not None:
                species_name = cohort.species_for_config(config_name)
                entry = cohort.entry(species_name)
                context = build_species_reuse_context(
                    species_name,
                    source,
                    cohort_id=cohort.cohort_id,
                    reference_index_sha256=entry.index_sha256,
                    reference_index_name=entry.index_name,
                    index_digests=self._cohort_digests,
                )
                tasks = {task.accession for task in cohort.tasks_for(species_name)}
                recorded = [run for run in self.db.get_samples(species_name, "quantified") if run in tasks]
                quant_dir = _data_root() / species_name / "work" / "quant"
                classify: Any = context.classify
            else:
                species_name = _species_name_from_config(config_name)
                work_dir = _species_work_dir(species_name)
                quant_dir = work_dir / "quant"
                recorded = self.db.get_samples(species_name, "quantified")
                if source.is_file() and work_dir.parent.name == species_name:
                    classify = build_species_reuse_context(species_name, source).classify
                else:
                    classify = classify_quantification
            probes.append(ReuseProbe(species_name, quant_dir, recorded, classify))
        return probes

    def _resolve_resource_profile(
        self,
        workers: Optional[int],
        threads: Optional[int],
        *,
        index_paths: Mapping[str, Path],
        data_root: Path,
        quant_slots: Optional[int] = None,
        fasterq_threads: Optional[int] = None,
        fasterq_slots: Optional[int] = None,
        compression_threads: Optional[int] = None,
        max_in_flight: Optional[int] = None,
        validation_slots: Optional[int] = None,
    ) -> PipelineResourceProfile:
        """Explicit ``workers`` and ``threads`` keep the historical profile; either one missing derives it
        from this host (:func:`build_host_resource_profile`)."""

        if workers is not None and threads is not None:
            profile = build_pipeline_resource_profile(
                workers,
                threads,
                quant_slots=quant_slots,
                fasterq_threads=fasterq_threads,
                fasterq_slots=fasterq_slots,
                compression_threads=compression_threads,
                validation_slots=validation_slots,
                max_in_flight=max_in_flight,
            )
            order = os.environ.get("AMALGKIT_PIPELINE_SCHEDULING_ORDER", "").strip() or SCHEDULING_ORDERS[0]
            if order not in SCHEDULING_ORDERS:
                raise ValueError(f"scheduling order must be one of {SCHEDULING_ORDERS}, not {order!r}")
            reservation_gb = os.environ.get("AMALGKIT_PIPELINE_RAW_RESERVATION_GB", "").strip()
            return replace(
                profile,
                scheduling_order=order,
                raw_reservation_bytes=int(float(reservation_gb) * 1024**3) if reservation_gb else 0,
            )
        index_sizes = []
        for path in index_paths.values():
            try:
                index_sizes.append(path.stat().st_size)
            except OSError:
                continue
        try:
            free_bytes: int | None = shutil.disk_usage(data_root).free
        except OSError:
            free_bytes = None
        return build_host_resource_profile(
            kallisto_index_bytes=max(index_sizes) if index_sizes else None,
            data_root_free_bytes=free_bytes,
            workers=workers,
            threads=threads,
            quant_slots=quant_slots,
            fasterq_threads=fasterq_threads,
            fasterq_slots=fasterq_slots,
            compression_threads=compression_threads,
            validation_slots=validation_slots,
            max_in_flight=max_in_flight,
        )

    def _apply_resource_profile(self, profile: PipelineResourceProfile) -> None:
        """Install the profile's stage semaphores; a host profile also sets the acquisition budgets."""

        self._resource_profile = profile
        self._quant_semaphore = threading.BoundedSemaphore(profile.quant_slots)
        self._fasterq_semaphore = threading.BoundedSemaphore(profile.fasterq_slots)
        self._raw_validation_semaphore = threading.BoundedSemaphore(profile.validation_slots)
        if profile.source == "host":
            self.acquisition_slots = profile.acquisition_slots
            self._acquisition_semaphore = (
                threading.BoundedSemaphore(profile.acquisition_slots) if profile.acquisition_slots else None
            )
            self.ena_file_workers = profile.ena_file_workers

    def _discover_and_annotate(self, config_name: str, max_gb: float, threads: int) -> List[Dict[str, Any]]:
        """One species' discovery plus its raw-input scan, both off the scheduling thread."""

        tasks = self.discover_species_tasks(config_name, max_gb, threads)
        _annotate_raw_input_state(tasks)
        return tasks

    def _run_all_locked(
        self,
        species_list: List[str],
        max_gb: float,
        workers: Optional[int],
        threads: Optional[int],
        *,
        quant_slots: Optional[int] = None,
        fasterq_threads: Optional[int] = None,
        fasterq_slots: Optional[int] = None,
        compression_threads: Optional[int] = None,
        max_in_flight: Optional[int] = None,
        validation_slots: Optional[int] = None,
        requantification_policy: Optional[str] = None,
        discovery_workers: Optional[int] = None,
        reuse_dry_run_fraction: Optional[float] = None,
        resource_profile: Optional[PipelineResourceProfile] = None,
        scheduling_order: Optional[str] = None,
    ) -> Dict[str, Any]:
        """:meth:`run_all` while the campaign writer lock is held: resolve the profile, then run under it."""

        run_started = time.monotonic()
        self._cohort_discovered = {}
        data_root = _data_root()

        if requantification_policy is not None:
            requantification_policy = requantification_policy.strip().lower()
            if requantification_policy not in REQUANTIFICATION_POLICIES:
                raise ValueError(
                    f"Invalid requantification policy {requantification_policy!r}; "
                    f"expected one of {sorted(REQUANTIFICATION_POLICIES)}"
                )
            self.requantification_policy = requantification_policy

        quant_dirs, index_paths = self._preflight_targets(species_list)
        profile = resource_profile or self._resolve_resource_profile(
            workers,
            threads,
            index_paths=index_paths,
            data_root=data_root,
            quant_slots=quant_slots,
            fasterq_threads=fasterq_threads,
            fasterq_slots=fasterq_slots,
            compression_threads=compression_threads,
            max_in_flight=max_in_flight,
            validation_slots=validation_slots,
        )
        if scheduling_order is not None:
            if scheduling_order not in SCHEDULING_ORDERS:
                raise ValueError(f"scheduling order must be one of {SCHEDULING_ORDERS}, not {scheduling_order!r}")
            profile = replace(profile, scheduling_order=scheduling_order)
        if profile.source != "host":
            # An explicit profile leaves acquisition to the environment the orchestrator read; record what applies.
            profile = replace(
                profile,
                acquisition_slots=self.acquisition_slots,
                ena_file_workers=self.ena_file_workers,
                ena_max_connections=_nonnegative_profile_int(
                    os.environ.get("AMALGKIT_ENA_MAX_CONNECTIONS", "0") or 0, 0, "AMALGKIT_ENA_MAX_CONNECTIONS"
                ),
                ena_api_connections=_nonnegative_profile_int(
                    os.environ.get("AMALGKIT_ENA_API_CONNECTIONS", "0") or 0, 0, "AMALGKIT_ENA_API_CONNECTIONS"
                ),
            )
        self._apply_resource_profile(profile)
        logger.info(format_resource_profile(profile))
        with resource_profile_environment(profile):
            return self._run_with_profile(
                species_list,
                max_gb,
                profile,
                data_root=data_root,
                quant_dirs=quant_dirs,
                index_paths=index_paths,
                discovery_workers=discovery_workers,
                reuse_dry_run_fraction=reuse_dry_run_fraction,
                run_started=run_started,
            )

    def _run_with_profile(
        self,
        species_list: List[str],
        max_gb: float,
        profile: PipelineResourceProfile,
        *,
        data_root: Path,
        quant_dirs: Dict[str, Path],
        index_paths: Dict[str, Path],
        discovery_workers: Optional[int],
        reuse_dry_run_fraction: Optional[float],
        run_started: float,
    ) -> Dict[str, Any]:
        """Preflight, startup resets, then streaming discovery and bounded sample execution."""
        from collections import defaultdict, deque

        # A producer that starts into a broken environment poisons the cohort
        # instead of processing it (observed 2026-09-03: 7,291 tasks failed
        # with EPERM because the process lacked external-volume write access,
        # and quant batches failed with ENOENT because 'amalgkit' was not on
        # PATH). Verify both before any discovery or scheduling work.
        try:
            preflight_facts = run_campaign_preflight(
                data_root, kallisto_quant_dirs=quant_dirs, kallisto_index_paths=index_paths
            )
        except PreflightError as exc:
            logger.error("Campaign preflight failed; refusing to start the producer.\n%s", exc)
            raise
        else:
            logger.info(
                "Campaign preflight passed: %s",
                ", ".join(f"{key}={value}" for key, value in sorted(preflight_facts.items())),
            )
        # Under the exclusive writer lock no other producer exists: every in-flight row is orphaned, and an
        # environment-class failure says nothing about its sample.
        startup: Dict[str, Any] = {
            "inflight_reset": self.db.reset_inflight_states(),
            "environment_failures_reset": self.db.reset_environment_failures(),
        }
        if os.environ.get("AMALGKIT_PIPELINE_STALE_DOWNLOADING_SECONDS", "").strip() not in ("", "0"):
            logger.info(
                "AMALGKIT_PIPELINE_STALE_DOWNLOADING_SECONDS is ignored: under the campaign writer lock every "
                "in-flight row is orphaned and was reset at startup."
            )
        if reuse_dry_run_fraction is None:
            raw_fraction = os.environ.get("AMALGKIT_PREFLIGHT_REUSE_DRY_RUN", "").strip()
            reuse_dry_run_fraction = float(raw_fraction) if raw_fraction else 0.0
        if reuse_dry_run_fraction > 0:
            try:
                startup["reuse_dry_run"] = reuse_dry_run(
                    self._reuse_probes(species_list), fraction=min(1.0, reuse_dry_run_fraction)
                )
            except PreflightError as exc:
                logger.error("Reuse dry-run refused the resume; no state was changed by discovery.\n%s", exc)
                raise
            logger.info(
                "Reuse dry-run passed: %s",
                ", ".join(
                    f"{species}={report['would_quarantine']}/{report['sampled']}"
                    for species, report in sorted(startup["reuse_dry_run"].items())
                ),
            )

        retry_policy = RetryPolicy.from_environment()
        # How long a producer whose only remaining work is in backoff waits in-process for the earliest retry
        # before exiting 75 (default 0: exit at once and report next_attempt_due_at to the supervisor).
        try:
            retry_wait_max = max(0.0, float(os.environ.get("AMALGKIT_RETRY_WAIT_MAX_SECONDS", "0") or 0))
        except ValueError:
            logger.warning("Invalid AMALGKIT_RETRY_WAIT_MAX_SECONDS; not waiting for backoffs")
            retry_wait_max = 0.0

        configured_discovery_workers: Any = discovery_workers
        if configured_discovery_workers is None:
            configured_discovery_workers = os.environ.get("AMALGKIT_PIPELINE_DISCOVERY_WORKERS")
        discovery_worker_count = _resource_int(
            configured_discovery_workers,
            DEFAULT_DISCOVERY_WORKERS,
            "AMALGKIT_PIPELINE_DISCOVERY_WORKERS",
        )
        discovery_worker_count = min(max(1, discovery_worker_count), max(1, len(species_list)))
        logger.info(
            "=== Streaming discovery and bounded sample execution: %d species, discovery workers=%d, "
            "workers=%d (requested %d), quant_slots=%d, quant_threads=%d total/%d per worker, "
            "validation_slots=%d, fasterq_threads=%d, fasterq_slots=%d in a dedicated fallback queue, "
            "compression_threads=%d, peak_stage_threads=%d, max_in_flight=%d, host_cpus=%d, order=%s, "
            "retry backoff %ds..%ds over %d attempts ===",
            len(species_list),
            discovery_worker_count,
            profile.workers,
            profile.requested_workers,
            profile.quant_slots,
            profile.effective_quant_threads,
            profile.quant_threads_per_worker,
            profile.validation_slots,
            profile.fasterq_threads,
            profile.fasterq_slots,
            profile.compression_threads,
            profile.peak_stage_threads,
            profile.max_in_flight,
            profile.host_cpu_count,
            profile.scheduling_order,
            retry_policy.base_seconds,
            retry_policy.cap_seconds,
            retry_policy.max_attempts,
        )
        quantified_by_species: Dict[str, int] = defaultdict(int)

        all_tasks: List[Dict[str, Any]] = []
        total_discovered = 0
        discovery_started = time.monotonic()
        discovery_elapsed = 0.0
        discovery_failures: Dict[str, str] = {}
        held_keys: set[tuple[str, str]] = set()
        held_by_reason: Dict[str, int] = defaultdict(int)
        # (due POSIX time, sequence, task): transient failures waiting out their backoff.
        deferred_heap: List[tuple[float, int, Dict[str, Any]]] = []
        sequence = itertools.count()
        # The unsubmitted queue, kept in order_campaign_tasks order as species arrive.
        pending: deque[Dict[str, Any]] = deque()

        breaker = EnvironmentCircuitBreaker(
            _resource_int(
                os.environ.get("AMALGKIT_ENV_FAILURE_TRIP"),
                DEFAULT_ENVIRONMENT_FAILURE_TRIP,
                "AMALGKIT_ENV_FAILURE_TRIP",
            )
        )
        execution_started = time.monotonic()
        completed = 0
        successful = 0
        reused = 0
        quantified_input_bytes = 0
        quantified_bases = 0
        environment_deferred = 0
        inflight_raw_bytes = 0

        def task_key(task: Mapping[str, Any]) -> tuple[str, str]:
            return str(task["species_name"]), str(task["srr"])

        def gate_tasks(tasks: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
            """Split freshly discovered tasks into ready, backoff (heap) and held (counted, never scheduled)."""

            records: Dict[tuple[str, str], RetryRecord] = {}
            for species in sorted({str(task["species_name"]) for task in tasks}):
                records.update(self.db.retry_records(species))
            now = time.time()
            ready: List[Dict[str, Any]] = []
            for task in tasks:
                record = records.get(task_key(task))
                gate, reason = retry_gate(record, retry_policy, now)
                if gate == GATE_HELD:
                    held_keys.add(task_key(task))
                    held_by_reason[(reason or "held").split(":", 1)[0]] += 1
                elif gate == GATE_BACKOFF and record is not None and record.next_attempt_at is not None:
                    heapq.heappush(deferred_heap, (record.next_attempt_at, next(sequence), task))
                else:
                    ready.append(task)
            return ready

        def enqueue(tasks: List[Dict[str, Any]]) -> None:
            """Merge tasks into the unsubmitted queue and re-order it (in-flight work is never preempted)."""

            nonlocal pending
            if tasks:
                pending = deque(
                    order_campaign_tasks(
                        [*pending, *tasks],
                        workers=profile.workers,
                        fasterq_slots=profile.fasterq_slots,
                        order=profile.scheduling_order,
                    )
                )

        def promote_due_retries() -> None:
            now = time.time()
            due: List[Dict[str, Any]] = []
            while deferred_heap and deferred_heap[0][0] <= now:
                due.append(heapq.heappop(deferred_heap)[2])
            if due:
                # Raw inputs may have changed while the task waited (a partial grew, a fallback finished).
                _annotate_raw_input_state(due)
                logger.info("Backoff elapsed for %d task(s); re-queued", len(due))
                enqueue(due)

        def raw_weight(task: Mapping[str, Any]) -> int:
            """Compressed FASTQ bytes a task may add to the data root (local complete inputs add none)."""

            return 0 if task.get("_raw_input_priority") == 0 else task_acquisition_weight(task)

        # A stop request (SIGINT/SIGTERM/SIGHUP via run_all, or an internal error) drains instead of exiting.
        stop: ProducerStop = getattr(self, "_producer_stop", None) or ProducerStop()
        interrupted_tasks: List[Dict[str, Any]] = []
        # Tasks whose acquisition slot or staging deadline expired are queued again in this run (bounded).
        deferral_requeues: Dict[tuple[str, str], int] = defaultdict(int)
        # AMALGKIT_PIPELINE_DEFERRAL_REQUEUES (default 20; 0 disables): how often one task whose admission or staging
        # deadline expired is queued again in the same run before it is left pending for the next start.
        deferral_requeue_limit = DEFAULT_DEFERRAL_REQUEUES
        raw_requeues = os.environ.get("AMALGKIT_PIPELINE_DEFERRAL_REQUEUES", "").strip()
        if raw_requeues:
            try:
                deferral_requeue_limit = int(raw_requeues)
            except ValueError:
                deferral_requeue_limit = -1
            if deferral_requeue_limit < 0:
                logger.warning(
                    "Invalid AMALGKIT_PIPELINE_DEFERRAL_REQUEUES=%r; using %d", raw_requeues, DEFAULT_DEFERRAL_REQUEUES
                )
                deferral_requeue_limit = DEFAULT_DEFERRAL_REQUEUES
        deferred_left_pending = 0
        last_wal_checkpoint = time.monotonic()
        wal_checkpoint_seconds = WAL_CHECKPOINT_INTERVAL_SECONDS
        raw_interval = os.environ.get("AMALGKIT_WAL_CHECKPOINT_SECONDS", "").strip()
        if raw_interval:
            try:
                wal_checkpoint_seconds = max(0.0, float(raw_interval))
            except ValueError:
                logger.warning(
                    "Invalid AMALGKIT_WAL_CHECKPOINT_SECONDS=%r; using %.0f", raw_interval, wal_checkpoint_seconds
                )

        with (
            concurrent.futures.ThreadPoolExecutor(max_workers=discovery_worker_count) as discovery_executor,
            concurrent.futures.ThreadPoolExecutor(max_workers=profile.workers) as executor,
            concurrent.futures.ThreadPoolExecutor(max_workers=profile.fasterq_slots) as fallback_executor,
        ):
            # Each species' tasks join the queue as soon as its own discovery finishes: a slow metadata or
            # reference step for one species never idles the acquisition and quantification slots.
            discovery_futures: Dict[concurrent.futures.Future[Any], str] = {}
            discovery_start_times: Dict[str, float] = {}
            discovery_threads = profile.requested_quant_threads
            for config_name in species_list:
                discovery_start_times[config_name] = time.monotonic()
                discovery_futures[
                    discovery_executor.submit(self._discover_and_annotate, config_name, max_gb, discovery_threads)
                ] = config_name

            ready_primary_tasks: deque[Dict[str, Any]] = deque()
            fallback_queue: deque[Dict[str, Any]] = deque()
            primary_futures: Dict[concurrent.futures.Future[Any], Dict[str, Any]] = {}
            fallback_futures: Dict[concurrent.futures.Future[Any], Dict[str, Any]] = {}

            def submit_next_primary() -> bool:
                nonlocal inflight_raw_bytes
                if breaker.tripped:
                    return False
                if ready_primary_tasks:
                    task = ready_primary_tasks.popleft()
                    task["_admitted_raw_bytes"] = 0
                elif pending:
                    weight = raw_weight(pending[0])
                    if (
                        profile.raw_reservation_bytes
                        and primary_futures
                        and inflight_raw_bytes + weight > profile.raw_reservation_bytes
                    ):
                        # The raw reservation is full: wait for an in-flight sample to finish (and reclaim its
                        # reads) instead of skipping ahead, so the size order holds. An oversized task still runs
                        # alone.
                        return False
                    task = pending.popleft()
                    task["_admitted_raw_bytes"] = weight
                    inflight_raw_bytes += weight
                else:
                    return False
                cohort_kwargs: Dict[str, Any] = (
                    {"metadata_path": task["metadata_path"]} if task.get("metadata_path") is not None else {}
                )
                if task.get("requantify"):
                    cohort_kwargs["requantify"] = True
                future = executor.submit(
                    self.process_single_sample,
                    task["srr"],
                    task["batch_idx"],
                    task["fastq_dir"],
                    task["config_path"],
                    task["species_name"],
                    profile.quant_threads_per_worker,
                    task.get("expected_paired"),
                    # A cached SRA is deliberately handled in this lane first:
                    # it may be reusable immediately and, if corrupt, must get
                    # one ENA recovery attempt without repeating the same
                    # expensive local extraction in this invocation. Network
                    # acquisitions defer NCBI work so they never wait behind
                    # the single mounted-disk extraction slot.
                    defer_ncbi_fallback=not task.get("_uses_local_sra", False),
                    **cohort_kwargs,
                )
                primary_futures[future] = task
                return True

            def submit_fallbacks() -> None:
                while not breaker.tripped and fallback_queue and len(fallback_futures) < profile.fasterq_slots:
                    task = fallback_queue.popleft()
                    self.db.set_state(
                        task["species_name"],
                        task["srr"],
                        "downloading",
                        error="running dedicated NCBI fallback after ENA acquisition failure",
                    )
                    future = fallback_executor.submit(
                        self._download_fastq_ncbi_only,
                        task["srr"],
                        task["fastq_dir"],
                        task.get("expected_paired"),
                    )
                    fallback_futures[future] = task

            fallback_error = "Download Failed (all sources: ENA FTP/HTTP, NCBI)"

            def defer_for_environment(task: Dict[str, Any], failure_class: str, detail: str) -> None:
                """An environment failure leaves the row pending and counts toward the circuit breaker."""

                nonlocal environment_deferred
                environment_deferred += 1
                srr = task["srr"]
                sp = task["species_name"]
                message = f"environment-deferred[{failure_class}]: {detail[:500]}"
                try:
                    self.db.set_state(sp, srr, "pending", error=message)
                except Exception as db_error:  # noqa: BLE001 - the DB may live on the volume that just failed
                    # Counted once below with the original failure: one failed sample is one breaker step.
                    logger.error("Cannot record %s %s as pending: %s", sp, srr, db_error)
                logger.error(
                    "[%d/%d] %s %s: environment failure (%s), left pending — %s",
                    completed,
                    total_discovered,
                    sp,
                    srr,
                    failure_class,
                    detail[:300],
                )
                if breaker.record_environment_failure(failure_class, f"{sp}/{srr}: {detail}"):
                    logger.critical(
                        "Environment circuit breaker tripped after %d consecutive environment failures (last: %s). "
                        "No new samples are admitted; in-flight samples drain and their rows stay pending. Fix the "
                        "host (volume mounted and writable, Full Disk Access/Removable Volumes permission, free "
                        "space, tools on PATH), then rerun the same command.",
                        breaker.threshold,
                        failure_class,
                    )

            def environment_class_for(failure: BaseException | str | None) -> str | None:
                found = environment_failure_class(failure, data_root)
                if found is None and failure is not None:
                    found = data_root_unavailable(data_root)
                return found

            def account_attempt(task: Dict[str, Any], error: str | None) -> None:
                """Retry bookkeeping for one finished attempt: clear, back off (in-run heap), or hold."""

                self._forget_acquisition_note(str(task["srr"]))
                error_class: str | None = None
                if error is not None and classify_sample_error(error) in RETRYABLE_ENVIRONMENT_CLASSES:
                    # This attempt reached account_attempt, so environment_class_for already said it is NOT an
                    # environment failure (for example ENOENT/EPERM on a path outside the data root). Record that
                    # verdict so the startup reset never re-reads the message as environmental and loops forever.
                    error_class = NON_ENVIRONMENT_OS_ERROR
                try:
                    outcome = self.db.record_attempt_outcome(
                        task["species_name"], str(task["srr"]), error, policy=retry_policy, error_class=error_class
                    )
                except Exception as exc:  # noqa: BLE001 - bookkeeping must not abort the campaign
                    logger.error("Cannot record the attempt of %s %s: %s", task["species_name"], task["srr"], exc)
                    return
                decision = outcome.get("decision")
                if decision == "backoff" and outcome.get("next_attempt_at") is not None:
                    heapq.heappush(deferred_heap, (float(outcome["next_attempt_at"]), next(sequence), task))
                    logger.info(
                        "%s %s: %s failure %d/%d; next attempt after %s UTC",
                        task["species_name"],
                        task["srr"],
                        outcome.get("error_class"),
                        outcome.get("attempts"),
                        retry_policy.max_attempts,
                        time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime(float(outcome["next_attempt_at"]))),
                    )
                elif decision == "held":
                    reason = {"layout_mismatch": "layout_mismatch", "permanent": "permanent"}.get(
                        str(outcome.get("retry_class")), "attempts_exhausted"
                    )
                    held_keys.add(task_key(task))
                    held_by_reason[reason] += 1
                    logger.warning(
                        "%s %s: held (%s, %s) after %d attempt(s); not retried automatically",
                        task["species_name"],
                        task["srr"],
                        reason,
                        outcome.get("error_class"),
                        outcome.get("attempts", 0),
                    )

            def progress_rates() -> str:
                elapsed = max(time.monotonic() - execution_started, 1e-9)
                hours = elapsed / 3600
                return (
                    f"{successful - reused} quantified + {reused} reused this run; "
                    f"{(successful - reused) / hours:.1f} samples/h, "
                    f"{quantified_input_bytes / 1e9 / hours:.2f} GB/h FASTQ, "
                    f"{quantified_bases / 1e9 / hours:.2f} Gbases/h"
                )

            def record_terminal(
                task: Dict[str, Any], res: Dict[str, Any] | None, error: Exception | None = None
            ) -> None:
                nonlocal completed, successful, reused, quantified_input_bytes, quantified_bases
                completed += 1
                srr = task["srr"]
                sp = task["species_name"]
                if error is not None:
                    failure_class = environment_class_for(error)
                    if failure_class is not None:
                        defer_for_environment(task, failure_class, f"Unexpected sample task error: {error}")
                        return
                    breaker.record_other_outcome()
                    error_message = f"Unexpected sample task error: {error}"
                    self.db.set_state(sp, srr, "failed", error=error_message)
                    account_attempt(task, error_message)
                    logger.error(
                        "[%d/%d] %s %s: Unexpected error — %s",
                        completed,
                        total_discovered,
                        sp,
                        srr,
                        error,
                    )
                    return
                if res is None:
                    logger.error(
                        "[%d/%d] %s %s: Unexpected empty task result",
                        completed,
                        total_discovered,
                        sp,
                        srr,
                    )
                    return
                if not res["quantified"] and res.get("error"):
                    failure_class = environment_class_for(str(res["error"]))
                    if failure_class is not None:
                        defer_for_environment(task, failure_class, str(res["error"]))
                        return
                breaker.record_other_outcome()
                if res["quantified"]:
                    successful += 1
                    quantified_by_species[sp] += 1
                    if res.get("skipped"):
                        reused += 1
                    else:
                        quantified_input_bytes += int(res.get("input_bytes") or 0)
                        quantified_bases += int(task.get("total_bases") or 0)
                    account_attempt(task, None)
                    status = "Skipped (Done)" if res.get("skipped") else "Done"
                    logger.info(
                        "[%d/%d] %s %s: %s (%s)",
                        completed,
                        total_discovered,
                        sp,
                        srr,
                        status,
                        progress_rates(),
                    )
                elif res["error"]:
                    if not res.get("download_deferred"):
                        account_attempt(task, str(res["error"]))
                    logger.warning(
                        "[%d/%d] %s %s: Failed (%s)",
                        completed,
                        total_discovered,
                        sp,
                        srr,
                        res["error"],
                    )

            def finish_discovery(future: concurrent.futures.Future[Any]) -> None:
                nonlocal total_discovered, discovery_elapsed
                config_name = discovery_futures.pop(future)
                try:
                    tasks = future.result()
                except Exception as exc:  # noqa: BLE001 - one species never aborts the campaign
                    discovery_failures[config_name] = str(exc)[:500]
                    logger.exception(
                        "Discovery failed for %s; skipping this species for this run (exit code 75)", config_name
                    )
                else:
                    all_tasks.extend(tasks)
                    total_discovered += len(tasks)
                    ready = gate_tasks(tasks)
                    enqueue(ready)
                    logger.info(
                        "Discovery complete for %s: %d tasks (%d ready, %d held, %d in backoff overall) in %.1fs; "
                        "total discovered %d, queued %d",
                        config_name,
                        len(tasks),
                        len(ready),
                        sum(1 for task in tasks if task_key(task) in held_keys),
                        len(deferred_heap),
                        time.monotonic() - discovery_start_times[config_name],
                        total_discovered,
                        len(pending),
                    )
                if not discovery_futures:
                    discovery_elapsed = max(time.monotonic() - discovery_started, 1e-9)
                    reusable_count = sum(task.get("_raw_input_priority") == 0 for task in all_tasks)
                    resumable_count = sum(task.get("_raw_input_priority") == 1 for task in all_tasks)
                    logger.info(
                        "=== Discovery complete: %d of %d species, %d tasks (%d with existing raw inputs, %d "
                        "resumable partials), %d held, %d in backoff, %.1fs wall time ===",
                        len(species_list) - len(discovery_failures),
                        len(species_list),
                        total_discovered,
                        reusable_count,
                        resumable_count,
                        len(held_keys),
                        len(deferred_heap),
                        discovery_elapsed,
                    )

            def requeue_deferred(task: Dict[str, Any], res: Mapping[str, Any]) -> bool:
                """Queue a task whose admission or staging deadline expired again in this run (bounded)."""

                nonlocal deferred_left_pending
                key = task_key(task)
                if res.get("deferral") == "first_pass":
                    return False
                if stop.requested or deferral_requeues[key] >= deferral_requeue_limit:
                    deferred_left_pending += 1
                    return False
                deferral_requeues[key] += 1
                _annotate_raw_input_state([task])
                enqueue([task])
                logger.info(
                    "%s %s: %s; queued again in this run (%d/%d)",
                    task["species_name"],
                    task["srr"],
                    res.get("error"),
                    deferral_requeues[key],
                    deferral_requeue_limit,
                )
                return True

            def drain_after_stop() -> None:
                """Cancel queued work, stop child process groups and wait for every running future."""

                for futures in (primary_futures, fallback_futures, discovery_futures):
                    for future in list(futures):
                        future.cancel()
                for pool in (discovery_executor, executor, fallback_executor):
                    pool.shutdown(wait=False, cancel_futures=True)
                stop.stop_children()
                running = {
                    future
                    for future in (*primary_futures, *fallback_futures, *discovery_futures)
                    if not future.cancelled()
                }
                while running:
                    finished, running = concurrent.futures.wait(running, timeout=1.0)
                    stop.stop_children()
                for future, task in [*primary_futures.items(), *fallback_futures.items()]:
                    if not future.cancelled():
                        interrupted_tasks.append(task)
                    else:
                        ready_primary_tasks.append(task)
                primary_futures.clear()
                fallback_futures.clear()
                discovery_futures.clear()
                logger.warning(
                    "Stopped (%s): %d in-flight sample(s) interrupted; every worker and child process has exited",
                    stop.signal_name,
                    len(interrupted_tasks),
                )

            try:
                while True:
                    if stop.requested:
                        break
                    if time.monotonic() - last_wal_checkpoint >= wal_checkpoint_seconds:
                        last_wal_checkpoint = time.monotonic()
                        try:
                            self.db.checkpoint_wal(mode="PASSIVE")
                        except Exception as exc:  # noqa: BLE001 - a periodic checkpoint never stops the run
                            logger.warning("Periodic WAL checkpoint failed: %s", exc)
                    promote_due_retries()
                    submit_fallbacks()
                    while len(primary_futures) < profile.max_in_flight:
                        if not submit_next_primary():
                            break

                    active_futures = set(primary_futures) | set(fallback_futures) | set(discovery_futures)
                    if not active_futures:
                        if breaker.tripped:
                            break
                        if pending or ready_primary_tasks or fallback_queue:
                            # This should be unreachable, but protects against a
                            # scheduler spin if future routing changes later.
                            logger.error("Execution queues are non-empty but no futures are active; stopping safely.")
                            break
                        if not deferred_heap:
                            break
                        wait_seconds = deferred_heap[0][0] - time.time()
                        if wait_seconds > retry_wait_max:
                            break
                        # Short sleeps: a stop request is noticed within a second.
                        time.sleep(max(0.0, min(wait_seconds, 1.0)))
                        continue
                    timeout = min(1.0, max(0.0, deferred_heap[0][0] - time.time())) if deferred_heap else 1.0
                    done, _ = concurrent.futures.wait(
                        active_futures,
                        timeout=timeout,
                        return_when=concurrent.futures.FIRST_COMPLETED,
                    )
                    for future in done:
                        if future in discovery_futures:
                            finish_discovery(future)
                            continue
                        if future in primary_futures:
                            task = primary_futures.pop(future)
                            inflight_raw_bytes -= int(task.pop("_admitted_raw_bytes", 0) or 0)
                            try:
                                res = future.result()
                            except Exception as exc:
                                record_terminal(task, None, exc)
                                continue
                            if res.get("fallback_deferred"):
                                fallback_queue.append(task)
                            elif res.get("download_deferred") and requeue_deferred(task, res):
                                continue
                            else:
                                record_terminal(task, res)
                            continue

                        task = fallback_futures.pop(future)
                        try:
                            fallback_ready = bool(future.result())
                        except Exception as exc:
                            fallback_ready = False
                            logger.error(
                                "Dedicated NCBI fallback raised for %s %s: %s",
                                task["species_name"],
                                task["srr"],
                                exc,
                            )
                            failure_class = environment_class_for(exc)
                            if failure_class is not None:
                                completed += 1
                                defer_for_environment(task, failure_class, f"NCBI fallback error: {exc}")
                                continue
                        if fallback_ready:
                            # Return the now-local sample to the ordinary lane.
                            # download_fastq() revalidates and reuses the FASTQ,
                            # then quantification competes under the shared CPU
                            # semaphore while the fallback executor immediately
                            # starts the next extraction.
                            ready_primary_tasks.appendleft(task)
                        else:
                            failure = fallback_error + self._layout_failure_suffix(str(task["srr"]))
                            self.db.set_state(task["species_name"], task["srr"], "failed", error=failure)
                            record_terminal(task, {"quantified": False, "error": failure})

            except BaseException:
                # An error in the scheduling loop must not leave queued samples running without the lock.
                if not stop.requested:
                    stop.request(0)
                drain_after_stop()
                raise
            if stop.requested:
                drain_after_stop()

            # Samples admitted to neither lane when the breaker tripped (or a stop arrived) keep no in-flight state.
            admission_note = (
                f"deferred: producer stopped by {stop.signal_name}"
                if stop.requested
                else "deferred: environment circuit breaker stopped admission"
            )
            for task in list(ready_primary_tasks) + list(fallback_queue):
                try:
                    self.db.set_state(task["species_name"], task["srr"], "pending", error=admission_note)
                except Exception as exc:  # noqa: BLE001 - reset at the next startup in any case
                    logger.error("Cannot return %s %s to pending: %s", task["species_name"], task["srr"], exc)
            # Samples interrupted by the stop: whatever their worker recorded (a killed download looks like a
            # failure), the attempt never counts and the row is pending again unless it finished.
            if interrupted_tasks:
                try:
                    interrupted_states = self.db.states_for([task_key(task) for task in interrupted_tasks])
                except Exception as exc:  # noqa: BLE001 - reset_inflight_states covers in-flight rows at startup
                    logger.error("Cannot read interrupted task states: %s", exc)
                    interrupted_states = {}
                for task in interrupted_tasks:
                    if interrupted_states.get(task_key(task)) in ("quantified", "quarantined", None):
                        continue
                    try:
                        self.db.set_state(
                            task["species_name"],
                            task["srr"],
                            "pending",
                            error=f"interrupted: producer stopped by {stop.signal_name}",
                        )
                    except Exception as exc:  # noqa: BLE001 - reset at the next startup in any case
                        logger.error("Cannot return %s %s to pending: %s", task["species_name"], task["srr"], exc)

        if not discovery_elapsed:
            discovery_elapsed = max(time.monotonic() - discovery_started, 1e-9)
        if breaker.tripped:
            logger.critical(
                "=== Stopped by the environment circuit breaker: %d environment failures, %d samples completed. "
                "Nothing was marked failed for the environment; rerun after fixing the host. ===",
                breaker.total,
                completed,
            )
        else:
            logger.info(
                "=== Acquisition and quantification complete for this invocation (%s); "
                "run the current downstream checkpoint launcher after the producer stops. ===",
                progress_rates(),
            )

        # Global progress summary mapping
        logger.info("=== Final Results Summary ===")
        try:
            all_counts = self.db.get_counts()
            for sp, c in all_counts.items():
                total = sum(c.values())
                q = c.get("quantified", 0)
                f = c.get("failed", 0)
                p = c.get("pending", 0)
                d = c.get("downloaded", 0) + c.get("downloading", 0)
                pct = (q / total * 100) if total > 0 else 0
                logger.info(f"{sp.ljust(30)} | {q}/{total} ({pct:3.0f}%) | " f"Fail: {f:3} | Pend: {p:3} | DL: {d:3}")
        except Exception as exc:
            logger.debug("Final results summary unavailable: %s", exc)

        execution_seconds = time.monotonic() - execution_started
        hours = max(execution_seconds, 1e-9) / 3600
        return self._run_summary(
            all_tasks,
            breaker=breaker,
            startup=startup,
            completed=completed,
            successful=successful,
            environment_deferred=environment_deferred,
            discovery_seconds=discovery_elapsed,
            execution_seconds=execution_seconds,
            run_seconds=time.monotonic() - run_started,
            held_keys=held_keys,
            held_by_reason=dict(held_by_reason),
            deferred_backoff=len(deferred_heap),
            next_attempt_due_at=min((entry[0] for entry in deferred_heap), default=None),
            discovery_failures=discovery_failures,
            throughput={
                "quantified_new": successful - reused,
                "reused": reused,
                "samples_per_hour": round((successful - reused) / hours, 2),
                "fastq_gb_per_hour": round(quantified_input_bytes / 1e9 / hours, 3),
                "gbases_per_hour": round(quantified_bases / 1e9 / hours, 3),
            },
            resource_profile=profile,
            stopped=stop.signal_name if stop.requested else None,
        )

    def _run_summary(
        self,
        all_tasks: List[Dict[str, Any]],
        *,
        breaker: EnvironmentCircuitBreaker,
        startup: Mapping[str, Any],
        completed: int,
        successful: int,
        environment_deferred: int,
        discovery_seconds: float,
        execution_seconds: float,
        run_seconds: float,
        held_keys: Optional[set[tuple[str, str]]] = None,
        held_by_reason: Optional[Mapping[str, int]] = None,
        deferred_backoff: int = 0,
        next_attempt_due_at: Optional[float] = None,
        discovery_failures: Optional[Mapping[str, str]] = None,
        throughput: Optional[Mapping[str, Any]] = None,
        resource_profile: Optional[PipelineResourceProfile] = None,
        stopped: Optional[str] = None,
    ) -> Dict[str, Any]:
        """The run summary and its exit code.

        * 69 ``environment_failure``: the breaker tripped or the database
          could not be read;
        * 75 ``stopped``: a signal (or an internal error) stopped the run
          before it drained; interrupted rows are pending;
        * 75 ``retryable_work_remains``: a scheduled task whose final row is
          not ``quantified`` or ``quarantined``, is not a ``permanent_drop``
          and is not held by the retry bookkeeping (a layout mismatch, a
          kallisto exit-1 failure or a spent attempt budget), a species whose
          discovery failed, or an output reconcile could not read
          (``unverified_outputs``): the same command should run again;
        * 65 ``needs_operator`` (cohort mode): nothing is left to run, but
          cohort rows are held, quarantined or above ``max_gb``
          (``needs_operator`` counts them); finalization would refuse them;
        * 0 ``drained``: every cohort row (cohort mode) or every scheduled
          task is quantified or excluded.
        """

        held = held_keys or set()
        keys = sorted({(str(task["species_name"]), str(task["srr"])) for task in all_tasks})
        remaining: Dict[str, int] = {}
        failed_by_class: Dict[str, int] = {}
        database_error: str | None = None
        try:
            final_states = self.db.states_for(keys)
            species_names = sorted({species for species, _ in keys})
            dropped = {
                (species, run)
                for species in species_names
                for run in self.db.get_excluded_srr_ids(species, reason_code="permanent_drop")
            }
            for key, state in final_states.items():
                if key in dropped or key in held or state in ("quantified", "quarantined"):
                    continue
                label = state or "missing"
                remaining[label] = remaining.get(label, 0) + 1
            wanted = set(keys)
            for species in species_names:
                for row in self.db.get_failed(species):
                    if (
                        (species, row["srr_id"]) in wanted
                        and (species, row["srr_id"]) not in dropped
                        and (species, row["srr_id"]) not in held
                    ):
                        error_class = classify_sample_error(row["error"])
                        failed_by_class[error_class] = failed_by_class.get(error_class, 0) + 1
            needs_operator, unverified = self._cohort_completeness(held)
        except Exception as exc:  # noqa: BLE001 - the DB may have gone with the volume
            database_error = str(exc)
            needs_operator, unverified = {}, 0
            logger.error("Cannot read final task states: %s", exc)
        if breaker.tripped or database_error is not None:
            status, exit_code = "environment_failure", EXIT_ENVIRONMENT
        elif stopped is not None:
            status, exit_code = "stopped", EXIT_RETRYABLE
        elif remaining or discovery_failures or unverified:
            status, exit_code = "retryable_work_remains", EXIT_RETRYABLE
        elif needs_operator:
            status, exit_code = "needs_operator", EXIT_NEEDS_OPERATOR
        else:
            status, exit_code = "drained", EXIT_OK
        summary: Dict[str, Any] = {
            "status": status,
            "exit_code": exit_code,
            "cohort_id": getattr(getattr(self, "cohort", None), "cohort_id", None),
            "tasks_scheduled": len(keys),
            "completed": completed,
            "quantified_this_run": successful,
            "environment_deferred": environment_deferred,
            "remaining_by_state": dict(sorted(remaining.items())),
            "failed_by_class": dict(sorted(failed_by_class.items())),
            "held_by_reason": dict(sorted((held_by_reason or {}).items())),
            "needs_operator": dict(sorted(needs_operator.items())),
            "unverified_outputs": unverified,
            "stopped_by": stopped,
            "deferred_backoff": deferred_backoff,
            "next_attempt_due_at": (
                time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime(next_attempt_due_at))
                if next_attempt_due_at is not None
                else None
            ),
            "discovery_failures": dict(sorted((discovery_failures or {}).items())),
            "throughput": dict(throughput or {}),
            "resource_profile": resource_profile.as_dict() if resource_profile is not None else None,
            "breaker": breaker.state(),
            "startup": dict(startup),
            "database_error": database_error,
            "discovery_seconds": round(discovery_seconds, 1),
            "execution_seconds": round(execution_seconds, 1),
            "run_seconds": round(run_seconds, 1),
        }
        logger.info("=== Producer summary: %s ===", json.dumps(summary, sort_keys=True, default=str))
        return summary

    def _cohort_completeness(self, held: set[tuple[str, str]]) -> tuple[Dict[str, int], int]:
        """Cohort rows nothing will run without the operator, by kind, and the outputs reconcile could not read.

        Kinds: ``quarantined`` (reuse contract or layout mismatch), ``held``
        (kallisto exit 1, spent attempt budget) and ``deferred_max_gb``
        (``total_bases`` above ``--max-gb``). Only cohort mode defines the
        complete set; elsewhere the first value is empty.
        """

        discovered: Mapping[str, Mapping[str, Sequence[str]]] = dict(getattr(self, "_cohort_discovered", {}) or {})
        unverified = sum(len(record.get("unavailable", ())) for record in discovered.values())
        counts: Dict[str, int] = {}
        if getattr(self, "cohort", None) is None:
            return counts, unverified
        for species, record in sorted(discovered.items()):
            accessions = record.get("accessions")
            if not accessions:
                continue
            deferred_max_gb = set(record.get("deferred_max_gb", ()))
            excluded = set(self.db.get_excluded_srr_ids(species, reason_code="permanent_drop"))
            states = self.db.states_for([(species, accession) for accession in accessions])
            for accession in accessions:
                state = states.get((species, accession))
                if accession in excluded or state == "quantified":
                    continue
                if state == "quarantined":
                    kind = "quarantined"
                elif (species, accession) in held:
                    kind = "held"
                elif accession in deferred_max_gb:
                    kind = "deferred_max_gb"
                else:
                    continue  # scheduled work: counted under remaining_by_state
                counts[kind] = counts.get(kind, 0) + 1
        return counts, unverified
