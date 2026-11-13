"""SQLite-backed progress tracking for the current RNA-seq pipeline.

The current progress database provides:
- O(1) status queries (no filesystem scanning)
- Thread-safe concurrent writes from parallel workers
- Instant dashboard: ``SELECT state, COUNT(*) GROUP BY species, state``
- Robust resume: state survives crashes, no data loss

Sample state machine::

    pending → downloading → downloaded → quantifying → quantified
                ↘ failed       ↘ failed        ↘ failed
"""

from __future__ import annotations

import calendar
import concurrent.futures
import json
import math
import os
import re
import sqlite3
import stat as stat_module
import sys
import threading
import time
from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING, Any, Callable, Dict, Iterable, List, Mapping, Optional

if TYPE_CHECKING:
    from metainformant.rna.engine.cohort import FrozenCohort

from metainformant.core.utils.logging import get_logger
from metainformant.rna.core.sample_utils import find_quantification_file
from metainformant.rna.engine.provenance import (
    QUANT_PROVENANCE_FILENAME,
    QUANT_STATUS_UNAVAILABLE,
    REUSABLE_QUANT_STATUSES,
    classify_quantification,
    is_unreadable_error,
    strict_reads,
)
from metainformant.rna.engine.reference_digest import IndexDigestCache

logger = get_logger(__name__)

# ---------- Constants ----------

VALID_STATES = frozenset(
    {
        "pending",
        "downloading",
        "downloaded",
        "quantifying",
        "quantified",
        "failed",
        "quarantined",
    }
)

EXCLUSION_REASON_CODES = frozenset({"permanent_drop", "re_download"})

DEFAULT_DB_PATH = Path("output/amalgkit/pipeline_progress.db")


def configured_progress_db_path(default: Path, configured: str | None) -> Path:
    """Resolve an explicitly migrated progress store, refusing implicit creation.

    An override must name an existing absolute regular file without symlinks or
    lexical parent traversal. Migration and storage-root attestation belong to
    the campaign operator; this helper never copies or initializes that store.
    """
    if configured is None:
        return default
    path = Path(configured)
    if (
        not configured
        or not path.is_absolute()
        or ".." in path.parts
        or any(part.is_symlink() for part in (path, *path.parents))
        or not path.is_file()
        or path.resolve(strict=True) != path
    ):
        raise ValueError("Configured progress database must be an existing absolute nonsymlink file")
    return path


SAMPLE_ERROR_CLASSES: dict[str, str] = {
    "environment_write_denied": "Operation not permitted",
    "environment_access_denied": "Permission denied",
    "environment_read_only": "Read-only file system",
    "environment_no_space": "No space left on device",
    "environment_device_unavailable": "Device not configured",
    "environment_missing_tool": "No such file or directory: '",
    # Before the transfer marker: a layout note is appended to "Download Failed (all sources ...)".
    "layout_mismatch": "layout_mismatch:",
    "transfer_all_sources_failed": "Download Failed (all sources",
    "extraction_timeout": "fasterq-dump timeout",
    "quantification_timeout": "Quant timeout",
    # Amalgkit 0.16.60 and the 0.16.108 pin (quant.py) both report a failed kallisto run as
    # "kallisto quant failed with exit code 1 for <run>.".
    "quantification_kallisto_exit_1": "kallisto quant failed with exit code 1",
    "quantification_failed": "Quantification Failed",
    "quantification_exception": "Quant exception",
}
"""Order-sensitive substring markers for durable failure classes.

The M-02 terminal-failure audit must distinguish environmental failures
(a broken producer environment that failed the whole cohort at once and is
fully retryable) from genuine per-sample transfer or quantification failures.
Markers are matched against the stored ``samples.error`` text; the first
matching class wins, so more specific markers precede generic ones.
"""


NON_ENVIRONMENT_OS_ERROR = "non_environment_os_error"
"""``last_error_class`` the producer records when an attempt's message reads like an environment failure
(``Permission denied``, ``No such file or directory: '...'``) but its own classifier decided otherwise (a
path outside the data root, a missing input file). A transient class: counted, backed off, finally held."""

# A missing *executable* is reported by a bare name ("No such file or directory: 'amalgkit'"); a missing file
# carries a path and says something about the sample, not the host.
_MISSING_TOOL_MESSAGE = re.compile(r"No such file or directory: '([^'/]+)'")


RETRYABLE_ENVIRONMENT_CLASSES = frozenset(
    {
        "environment_write_denied",
        "environment_access_denied",
        "environment_read_only",
        "environment_no_space",
        "environment_device_unavailable",
        "environment_missing_tool",
    }
)
"""Failure classes caused by the producer's host (EPERM, EACCES, EROFS, ENOSPC, ENXIO, ENOENT), not by the sample."""

TRANSFER_RESET_CLASSES = frozenset({"transfer_all_sources_failed"})
"""Transfer failures: possibly per-sample, reset only on explicit request."""

RETRYABLE_RESET_CLASSES = RETRYABLE_ENVIRONMENT_CLASSES | TRANSFER_RESET_CLASSES

COHORT_MIGRATION_ERROR_PREFIX = "reset-by-cohort-migration"

INFLIGHT_STATES = ("downloading", "downloaded", "quantifying")
"""States only a live producer may hold; under the exclusive campaign lock any such row is orphaned."""

STARTUP_INFLIGHT_RESET_PREFIX = "reset-inflight-at-startup"
STARTUP_ENVIRONMENT_RESET_PREFIX = "reset-environment-failure-at-startup"

AUDIT_FINGERPRINT_COLUMNS: tuple[tuple[str, str], ...] = (
    ("sidecar_ino", "INTEGER"),
    ("sidecar_size", "INTEGER"),
    ("sidecar_mtime_ns", "INTEGER"),
    ("config_sha256", "TEXT"),
    ("reference_sha256", "TEXT"),
    ("scope", "TEXT"),
    ("evidence", "TEXT"),
    ("output_fingerprint", "TEXT"),
)
"""Nullable ``quantification_audit`` columns that let a restart trust an unchanged classified output."""

RECONCILE_HEARTBEAT_EVERY = 1000

# ---------- Retry bookkeeping ----------

RETRY_COLUMNS: tuple[tuple[str, str], ...] = (
    ("attempts", "INTEGER NOT NULL DEFAULT 0"),
    ("last_error_class", "TEXT"),
    ("next_attempt_at", "TEXT"),
)
"""Additive ``samples`` columns: counted sample attempts, the last failure class and the UTC
(``YYYY-MM-DD HH:MM:SS``, the ``datetime('now')`` format) time before which the scheduler must not retry."""

RETRY_ENVIRONMENT = "environment"
RETRY_TRANSIENT = "transient"
RETRY_LAYOUT_MISMATCH = "layout_mismatch"
RETRY_PERMANENT = "permanent"

PERMANENT_ERROR_CLASSES = frozenset({"quantification_kallisto_exit_1"})
"""Failures that repeat byte-for-byte on every attempt: never retried automatically.

kallisto exits 1 for the same accessions on every run (canonical root: the
106 apis ``permanent_drop`` rows; October root: apis SRR21601882-SRR21601887
and 13 bombus SRR237263xx/SRR237264xx, each failing again in every launch).
"""

LAYOUT_MISMATCH_ERROR_CLASSES = frozenset({"layout_mismatch"})
"""The acquired files contradict the frozen library layout. Downloading again yields the same files
(2026-09-25..27 logs: SRR7704201, SRR18219008, SRR15496110 re-downloaded and rescanned at every start),
so the row is quarantined for reconciliation instead of being retried."""

RETRY_RELEASED_PREFIX = "released-by-operator"


def retry_class_for(error_class: str | None) -> str:
    """The retry class of a :func:`classify_sample_error` class.

    ``environment`` (the host, not the sample: retried immediately at the next
    start), ``layout_mismatch`` and ``permanent`` (held until an operator
    releases them) or ``transient`` (network, source, timeouts and anything
    unrecognized: exponential backoff, then held once the attempt budget is
    spent).
    """

    if error_class in RETRYABLE_ENVIRONMENT_CLASSES:
        return RETRY_ENVIRONMENT
    if error_class in LAYOUT_MISMATCH_ERROR_CLASSES:
        return RETRY_LAYOUT_MISMATCH
    if error_class in PERMANENT_ERROR_CLASSES:
        return RETRY_PERMANENT
    return RETRY_TRANSIENT


def _setting_int(name: str, default: int, minimum: int) -> int:
    raw = os.environ.get(name, "").strip()
    if not raw:
        return default
    try:
        value = int(raw)
    except ValueError:
        logger.warning("Invalid %s=%r; using %d", name, raw, default)
        return default
    if value < minimum:
        logger.warning("Invalid %s=%r; using %d", name, raw, default)
        return default
    return value


@dataclass(frozen=True)
class RetryPolicy:
    """Exponential backoff for transient sample failures, capped, with a bounded attempt budget.

    Defaults: the first retry waits 1 h, then 2 h, 4 h, 8 h, 16 h (cap 24 h);
    the 6th failed attempt holds the row for the operator, 31 h after the
    first. The window is wider than an overnight network or ENA outage (a
    sample failing every hour of an 8 h outage spends 4 attempts), but a
    source that is really gone (the 27 rows that failed "Download Failed (all
    sources)" from 2026-09-14 to 09-29 and were re-attempted at every start)
    stops consuming acquisition slots. Overrides:
    ``AMALGKIT_RETRY_BACKOFF_BASE_SECONDS``,
    ``AMALGKIT_RETRY_BACKOFF_CAP_SECONDS``, ``AMALGKIT_MAX_SAMPLE_ATTEMPTS``.
    """

    base_seconds: int = 3600
    cap_seconds: int = 86400
    max_attempts: int = 6

    def __post_init__(self) -> None:
        if self.base_seconds < 1 or self.cap_seconds < self.base_seconds or self.max_attempts < 1:
            raise ValueError("retry policy needs 1 <= base <= cap seconds and at least one attempt")

    @classmethod
    def from_environment(cls) -> "RetryPolicy":
        base = _setting_int("AMALGKIT_RETRY_BACKOFF_BASE_SECONDS", cls.base_seconds, 1)
        cap = max(base, _setting_int("AMALGKIT_RETRY_BACKOFF_CAP_SECONDS", cls.cap_seconds, 1))
        return cls(base, cap, _setting_int("AMALGKIT_MAX_SAMPLE_ATTEMPTS", cls.max_attempts, 1))

    def backoff_seconds(self, attempts: int) -> int:
        """Delay after ``attempts`` counted failures: ``base * 2**(attempts - 1)``, at most ``cap``."""

        if attempts <= 0:
            return 0
        return min(self.cap_seconds, self.base_seconds * (1 << min(attempts - 1, 40)))


def _utc_text(epoch: float) -> str:
    return time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime(epoch))


def _utc_epoch(text: str | None) -> float | None:
    if not text:
        return None
    try:
        return float(calendar.timegm(time.strptime(text, "%Y-%m-%d %H:%M:%S")))
    except ValueError:
        return None


@dataclass(frozen=True)
class RetryRecord:
    """One row's retry bookkeeping (``next_attempt_at`` as POSIX seconds)."""

    species: str
    srr_id: str
    state: str
    attempts: int
    last_error_class: Optional[str]
    next_attempt_at: Optional[float]

    @property
    def retry_class(self) -> str | None:
        return retry_class_for(self.last_error_class) if self.last_error_class else None


GATE_READY = "ready"
GATE_BACKOFF = "backoff"
GATE_HELD = "held"


def retry_gate(record: RetryRecord | None, policy: RetryPolicy, now: float) -> tuple[str, str | None]:
    """Whether a task may run now: ``(ready|backoff|held, reason)``.

    Held: a layout mismatch, a permanent failure, or a transient failure whose
    attempt budget is spent. Backoff: a transient failure whose
    ``next_attempt_at`` is still in the future. Environment failures and rows
    without bookkeeping are ready.
    """

    if record is None:
        return GATE_READY, None
    retry_class = record.retry_class
    if retry_class == RETRY_LAYOUT_MISMATCH:
        return GATE_HELD, RETRY_LAYOUT_MISMATCH
    if retry_class == RETRY_PERMANENT:
        return GATE_HELD, f"{RETRY_PERMANENT}:{record.last_error_class}"
    if retry_class == RETRY_TRANSIENT and record.attempts >= policy.max_attempts:
        return GATE_HELD, f"attempts_exhausted:{record.last_error_class}"
    if retry_class != RETRY_ENVIRONMENT and record.next_attempt_at is not None and record.next_attempt_at > now:
        return GATE_BACKOFF, f"backoff:{record.last_error_class}"
    return GATE_READY, None


@dataclass(frozen=True)
class ReconcileScope:
    """What a recorded classification was made against; a restart trusts it only under the same scope.

    ``config_sha256`` is the species' runtime configuration bytes,
    ``reference_sha256`` its reference identity (the cohort kallisto index
    sha256, or the reference manifest bytes outside a cohort) and ``scope``
    everything else the verdict depends on (cohort id, pinned Amalgkit
    runtime, classifier revision).
    """

    config_sha256: str
    reference_sha256: str
    scope: str


SidecarFingerprint = tuple[int, int, int]


def sidecar_fingerprint(path: Path) -> SidecarFingerprint | None:
    """``(st_ino, st_size, st_mtime_ns)`` of a regular, non-symlink file (one ``lstat``), else ``None``."""

    try:
        info = os.lstat(path)
    except OSError:
        return None
    if not stat_module.S_ISREG(info.st_mode):
        return None
    return (info.st_ino, info.st_size, info.st_mtime_ns)


def output_fingerprint(sample_dir: Path, quantification_file: Any) -> str | None:
    """``name:ino:size:mtime_ns`` of the abundance table a verdict covered (one ``lstat``), else ``None``.

    Recorded next to the sidecar fingerprint, so a trusted verdict also
    notices an abundance table that was truncated, replaced or removed after
    it was judged.
    """

    if not isinstance(quantification_file, str) or Path(quantification_file).name != quantification_file:
        return None
    if quantification_file in ("", ".", ".."):
        return None
    fingerprint = sidecar_fingerprint(sample_dir / quantification_file)
    if fingerprint is None:
        return None
    return f"{quantification_file}:{fingerprint[0]}:{fingerprint[1]}:{fingerprint[2]}"


def _recorded_output_matches(sample_dir: Path, recorded: Any) -> bool:
    if not isinstance(recorded, str) or recorded.count(":") < 3:
        return False
    name = recorded.rsplit(":", 3)[0]
    return output_fingerprint(sample_dir, name) == recorded


def classify_sample_error(error: str | None) -> str:
    """Classify a stored sample failure message into a coarse retry class."""

    if not error:
        return "unrecorded"
    for class_name, marker in SAMPLE_ERROR_CLASSES.items():
        if marker in error:
            if class_name == "environment_missing_tool" and not _MISSING_TOOL_MESSAGE.search(error):
                # A path (not a bare executable name) is missing: per-sample, not the host.
                continue
            return class_name
    return "unclassified"


_SCHEMA = """
CREATE TABLE IF NOT EXISTS samples (
    species     TEXT    NOT NULL,
    srr_id      TEXT    NOT NULL,
    state       TEXT    NOT NULL DEFAULT 'pending',
    error       TEXT,
    updated_at  TEXT    NOT NULL DEFAULT (datetime('now')),
    PRIMARY KEY (species, srr_id)
);

CREATE INDEX IF NOT EXISTS idx_species_state
    ON samples(species, state);

CREATE TABLE IF NOT EXISTS quantification_audit (
    species                    TEXT NOT NULL,
    srr_id                     TEXT NOT NULL,
    status                     TEXT NOT NULL,
    reason                     TEXT NOT NULL,
    contract_id                TEXT,
    observed_amalgkit_version  TEXT,
    observed_release_tag       TEXT,
    observed_source_revision   TEXT,
    checked_at                 TEXT NOT NULL DEFAULT (datetime('now')),
    PRIMARY KEY (species, srr_id)
);

CREATE INDEX IF NOT EXISTS idx_quantification_audit_status
    ON quantification_audit(status);

CREATE TABLE IF NOT EXISTS sample_exclusions (
    species        TEXT    NOT NULL,
    srr_id         TEXT    NOT NULL,
    reason_code    TEXT    NOT NULL,
    reason_detail  TEXT,
    recorded_by    TEXT,
    recorded_at    TEXT    NOT NULL DEFAULT (datetime('now')),
    PRIMARY KEY (species, srr_id)
);

CREATE INDEX IF NOT EXISTS idx_sample_exclusions_reason
    ON sample_exclusions(species, reason_code);
"""


# ---------- ProgressDB ----------


class ProgressDB:
    """SQLite-backed sample progress tracker.

    Thread-safe: uses a lock around all writes and ``check_same_thread=False``.

    Usage::

        db = ProgressDB()  # uses default path
        db.init_species("amellifera", ["SRR1", "SRR2", "SRR3"])
        db.set_state("amellifera", "SRR1", "downloading")
        print(db.get_counts())
    """

    def __init__(self, db_path: Optional[Path] = None, *, journal_mode: str | None = None):
        selected_mode = (
            journal_mode if journal_mode is not None else os.environ.get("AMALGKIT_PROGRESS_DB_JOURNAL_MODE", "WAL")
        )
        if selected_mode not in {"WAL", "DELETE"}:
            raise ValueError("Progress database journal mode must be WAL or DELETE")
        self.db_path = Path(db_path) if db_path else DEFAULT_DB_PATH
        self.db_path.parent.mkdir(parents=True, exist_ok=True)
        self._lock = threading.Lock()
        self._conn = sqlite3.connect(
            str(self.db_path),
            check_same_thread=False,
            isolation_level="DEFERRED",
        )
        try:
            self._conn.execute("PRAGMA busy_timeout=5000")
            actual_mode = self._conn.execute(f"PRAGMA journal_mode={selected_mode}").fetchone()
            if actual_mode != (selected_mode.lower(),):
                raise RuntimeError("Progress database journal mode was not applied")
            if selected_mode == "DELETE":
                self._conn.execute("PRAGMA synchronous=FULL")
                self._conn.execute("PRAGMA mmap_size=0")
                if self._conn.execute("PRAGMA synchronous").fetchone() != (2,) or self._conn.execute(
                    "PRAGMA mmap_size"
                ).fetchone() != (0,):
                    raise RuntimeError("Progress database durable rollback policy was not applied")
            self._conn.executescript(_SCHEMA)
            self._conn.commit()
        except BaseException:
            self._conn.close()
            raise
        logger.debug(f"ProgressDB opened at {self.db_path}")

    # ---- Write operations ----

    def init_species(self, species: str, srr_ids: List[str]) -> int:
        """Bulk-insert samples as 'pending'. Skips IDs that already exist.

        Returns:
            Number of newly inserted samples.
        """
        with self._lock:
            cur = self._conn.executemany(
                "INSERT OR IGNORE INTO samples (species, srr_id, state) VALUES (?, ?, 'pending')",
                [(species, srr) for srr in srr_ids],
            )
            self._conn.commit()
            inserted = cur.rowcount
            logger.info(f"init_species({species}): {inserted} new / {len(srr_ids)} total")
            return inserted

    def prune_species(self, species: str, srr_ids: List[str]) -> int:
        """Report progress rows outside the current task universe; never delete them.

        Rows outside ``srr_ids`` used to be deleted when a refreshed metadata
        query changed the set of public runs. That silently dropped quantified
        states, failure history and the rows behind exclusions, so this is now
        report-only: every row is kept with its state and error, and the
        scheduler simply never builds a task for it. Task building, not row
        deletion, defines the work; finalization must count against the task
        universe (for example a frozen cohort), not against raw table rows.

        Returns:
            Number of rows outside ``srr_ids`` (all kept).
        """

        outside = self.rows_outside(species, srr_ids)
        total = sum(outside.values())
        if total:
            logger.info(
                "prune_species(%s): kept %d progress rows outside the current task universe (%s)",
                species,
                total,
                ", ".join(f"{state}={count}" for state, count in sorted(outside.items())),
            )
        return total

    def rows_outside(self, species: str, srr_ids: Iterable[str]) -> Dict[str, int]:
        """Count a species' progress rows that are not in ``srr_ids``, by state."""

        keep = {str(srr_id) for srr_id in srr_ids}
        with self._lock:
            rows = self._conn.execute("SELECT srr_id, state FROM samples WHERE species = ?", (species,)).fetchall()
        counts: Dict[str, int] = {}
        for srr_id, state in rows:
            if srr_id not in keep:
                counts[state] = counts.get(state, 0) + 1
        return counts

    def set_state(self, species: str, srr_id: str, state: str, error: Optional[str] = None) -> None:
        """Atomically transition a sample to a new state."""
        if state not in VALID_STATES:
            raise ValueError(f"Invalid state '{state}'. Must be one of {VALID_STATES}")
        for attempt in range(3):
            try:
                with self._lock:
                    self._conn.execute(
                        "UPDATE samples SET state = ?, error = ?, updated_at = datetime('now') "
                        "WHERE species = ? AND srr_id = ?",
                        (state, error, species, srr_id),
                    )
                    self._conn.commit()
                return
            except sqlite3.OperationalError:
                if attempt < 2:
                    time.sleep(0.1 * (attempt + 1))
                else:
                    raise

    def record_quantification_audit(
        self,
        species: str,
        srr_id: str,
        *,
        status: str,
        reason: str,
        contract_id: Optional[str] = None,
        observed_amalgkit_version: Optional[str] = None,
        observed_release_tag: Optional[str] = None,
        observed_source_revision: Optional[str] = None,
        state: Optional[str] = None,
    ) -> None:
        """Record provenance compatibility and optionally transition sample state atomically."""

        if state is not None and state not in VALID_STATES:
            raise ValueError(f"Invalid state '{state}'. Must be one of {VALID_STATES}")
        with self._lock:
            self._conn.execute(
                """
                INSERT INTO quantification_audit (
                    species, srr_id, status, reason, contract_id,
                    observed_amalgkit_version, observed_release_tag,
                    observed_source_revision, checked_at
                ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, datetime('now'))
                ON CONFLICT(species, srr_id) DO UPDATE SET
                    status = excluded.status,
                    reason = excluded.reason,
                    contract_id = excluded.contract_id,
                    observed_amalgkit_version = excluded.observed_amalgkit_version,
                    observed_release_tag = excluded.observed_release_tag,
                    observed_source_revision = excluded.observed_source_revision,
                    checked_at = excluded.checked_at
                """,
                (
                    species,
                    srr_id,
                    status,
                    reason,
                    contract_id,
                    observed_amalgkit_version,
                    observed_release_tag,
                    observed_source_revision,
                ),
            )
            self._clear_audit_fingerprint_locked(species, srr_id)
            if state is not None:
                self._conn.execute(
                    "UPDATE samples SET state = ?, error = ?, updated_at = datetime('now') "
                    "WHERE species = ? AND srr_id = ?",
                    (state, reason if state == "quarantined" else None, species, srr_id),
                )
            self._conn.commit()

    def get_quantification_audit_counts(self, species: Optional[str] = None) -> Dict[str, Dict[str, int]]:
        """Return provenance compatibility counts grouped by species and status."""

        with self._lock:
            if species:
                rows = self._conn.execute(
                    "SELECT species, status, COUNT(*) FROM quantification_audit "
                    "WHERE species = ? GROUP BY species, status",
                    (species,),
                ).fetchall()
            else:
                rows = self._conn.execute(
                    "SELECT species, status, COUNT(*) FROM quantification_audit GROUP BY species, status"
                ).fetchall()
        result: Dict[str, Dict[str, int]] = {}
        for sp, status, count in rows:
            result.setdefault(sp, {})[status] = count
        return result

    def bulk_set_state(self, species: str, srr_ids: List[str], state: str, error: Optional[str] = None) -> None:
        """Set state for multiple samples at once."""
        if state not in VALID_STATES:
            raise ValueError(f"Invalid state '{state}'. Must be one of {VALID_STATES}")
        with self._lock:
            self._conn.executemany(
                "UPDATE samples SET state = ?, error = ?, updated_at = datetime('now') "
                "WHERE species = ? AND srr_id = ?",
                [(state, error, species, srr) for srr in srr_ids],
            )
            self._conn.commit()

    def record_exclusions(
        self,
        species: str,
        entries: List[Dict[str, str]],
    ) -> int:
        """Record or update sample exclusions (upsert by species + accession).

        Each entry is a mapping with keys ``srr_id`` (required) plus optional
        ``reason_code`` (default ``permanent_drop``), ``reason_detail``, and
        ``recorded_by``. Reason codes: ``permanent_drop`` removes the sample
        from acquisition/quantification eligibility; ``re_download`` marks it
        for a fresh ENA transfer without blocking eligibility.
        """
        rows = []
        for entry in entries:
            reason_code = entry.get("reason_code", "permanent_drop")
            if reason_code not in EXCLUSION_REASON_CODES:
                raise ValueError(
                    f"Invalid exclusion reason code '{reason_code}'. Must be one of {sorted(EXCLUSION_REASON_CODES)}"
                )
            rows.append(
                (
                    species,
                    entry["srr_id"],
                    reason_code,
                    entry.get("reason_detail"),
                    entry.get("recorded_by"),
                )
            )
        with self._lock:
            self._conn.executemany(
                """
                INSERT INTO sample_exclusions (
                    species, srr_id, reason_code, reason_detail, recorded_by, recorded_at
                ) VALUES (?, ?, ?, ?, ?, datetime('now'))
                ON CONFLICT(species, srr_id) DO UPDATE SET
                    reason_code = excluded.reason_code,
                    reason_detail = excluded.reason_detail,
                    recorded_by = excluded.recorded_by,
                    recorded_at = excluded.recorded_at
                """,
                rows,
            )
            self._conn.commit()
        logger.info("record_exclusions(%s): %d entries", species, len(rows))
        return len(rows)

    def get_exclusions(self, species: Optional[str] = None) -> List[Dict[str, str]]:
        """Return recorded exclusions, optionally restricted to one species."""
        with self._lock:
            if species:
                rows = self._conn.execute(
                    "SELECT species, srr_id, reason_code, reason_detail, recorded_by, recorded_at "
                    "FROM sample_exclusions WHERE species = ? ORDER BY srr_id",
                    (species,),
                ).fetchall()
            else:
                rows = self._conn.execute(
                    "SELECT species, srr_id, reason_code, reason_detail, recorded_by, recorded_at "
                    "FROM sample_exclusions ORDER BY species, srr_id"
                ).fetchall()
        return [
            {
                "species": r[0],
                "srr_id": r[1],
                "reason_code": r[2],
                "reason_detail": r[3],
                "recorded_by": r[4],
                "recorded_at": r[5],
            }
            for r in rows
        ]

    def get_excluded_srr_ids(self, species: str, reason_code: Optional[str] = None) -> set[str]:
        """Return the set of excluded SRR IDs for a species.

        With ``reason_code`` given, only exclusions carrying that code are
        returned; the orchestrator passes ``permanent_drop`` so re-download
        markers never block eligibility.
        """
        with self._lock:
            if reason_code is not None:
                rows = self._conn.execute(
                    "SELECT srr_id FROM sample_exclusions WHERE species = ? AND reason_code = ?",
                    (species, reason_code),
                ).fetchall()
            else:
                rows = self._conn.execute(
                    "SELECT srr_id FROM sample_exclusions WHERE species = ?",
                    (species,),
                ).fetchall()
        return {r[0] for r in rows}

    def remove_exclusion(self, species: str, srr_id: str) -> bool:
        """Delete a recorded exclusion. Returns True when a row was removed."""
        with self._lock:
            cur = self._conn.execute(
                "DELETE FROM sample_exclusions WHERE species = ? AND srr_id = ?",
                (species, srr_id),
            )
            self._conn.commit()
        return cur.rowcount > 0

    # ---- Read operations ----

    def get_state(self, species: str, srr_id: str) -> Optional[str]:
        """Get the current state of a single sample, or None if not tracked."""
        with self._lock:
            row = self._conn.execute(
                "SELECT state FROM samples WHERE species = ? AND srr_id = ?",
                (species, srr_id),
            ).fetchone()
        return row[0] if row else None

    def get_counts(self, species: Optional[str] = None) -> Dict[str, Dict[str, int]]:
        """Get sample counts grouped by species and state.

        Returns::

            {
                "amellifera": {"pending": 2900, "quantified": 254, ...},
                "atta_cephalotes": {"pending": 215, "downloaded": 5, ...},
            }
        """
        with self._lock:
            if species:
                rows = self._conn.execute(
                    "SELECT species, state, COUNT(*) FROM samples " "WHERE species = ? GROUP BY species, state",
                    (species,),
                ).fetchall()
            else:
                rows = self._conn.execute(
                    "SELECT species, state, COUNT(*) FROM samples GROUP BY species, state"
                ).fetchall()

        result: Dict[str, Dict[str, int]] = {}
        for sp, state, count in rows:
            if sp not in result:
                result[sp] = {}
            result[sp][state] = count
        return result

    def get_total_counts(self) -> Dict[str, int]:
        """Get total sample count per species.

        Returns::
            {"amellifera": 3154, "atta_cephalotes": 220, ...}
        """
        with self._lock:
            rows = self._conn.execute("SELECT species, COUNT(*) FROM samples GROUP BY species").fetchall()
        return {sp: count for sp, count in rows}

    def get_samples(self, species: str, state: str) -> List[str]:
        """Get all SRR IDs in a given state for a species."""
        with self._lock:
            rows = self._conn.execute(
                "SELECT srr_id FROM samples WHERE species = ? AND state = ?",
                (species, state),
            ).fetchall()
        return [r[0] for r in rows]

    def get_failed(self, species: Optional[str] = None) -> List[Dict[str, str]]:
        """Get all failed samples with their error messages."""
        with self._lock:
            if species:
                rows = self._conn.execute(
                    "SELECT species, srr_id, error, updated_at FROM samples "
                    "WHERE state = 'failed' AND species = ? ORDER BY updated_at DESC",
                    (species,),
                ).fetchall()
            else:
                rows = self._conn.execute(
                    "SELECT species, srr_id, error, updated_at FROM samples "
                    "WHERE state = 'failed' ORDER BY species, updated_at DESC",
                ).fetchall()
        return [{"species": r[0], "srr_id": r[1], "error": r[2], "updated_at": r[3]} for r in rows]

    def reset_stale_downloading(self, stale_seconds: int, error: str | None = None) -> int:
        """Reset stale downloading states to pending for resume resilience."""

        if stale_seconds <= 0:
            return 0
        if error is None:
            error = f"reset stale downloading heartbeat > {stale_seconds}s"
        with self._lock:
            rows = self._conn.execute(
                "SELECT species, srr_id FROM samples WHERE state = 'downloading' "
                "AND (strftime('%s', 'now') - strftime('%s', updated_at)) > ?",
                (stale_seconds,),
            ).fetchall()
            if not rows:
                return 0
            self._conn.executemany(
                "UPDATE samples SET state = 'pending', error = ?, updated_at = datetime('now') "
                "WHERE species = ? AND srr_id = ?",
                [(error, species, srr_id) for species, srr_id in rows],
            )
            self._conn.commit()
        return len(rows)

    # ---- Startup (only under the exclusive campaign writer lock) ----

    def reset_inflight_states(self) -> Dict[str, int]:
        """Return every ``downloading``/``downloaded``/``quantifying`` row to ``pending``, with no age threshold.

        Call this only while holding the campaign writer lock
        (:mod:`metainformant.rna.engine.campaign_lock`): then no other producer
        exists, so every in-flight row was orphaned by a previous producer that
        died (SIGKILL, crash, USB drop, reboot). Raw inputs stay on disk and are
        revalidated and reused. The previous error is kept after a
        ``reset-inflight-at-startup[<state>]`` prefix. Returns counts by the
        state each row was reset from.
        """

        placeholders = ",".join("?" for _ in INFLIGHT_STATES)
        with self._lock:
            self._conn.commit()
            self._conn.execute("BEGIN IMMEDIATE")
            try:
                counts = {
                    state: count
                    for state, count in self._conn.execute(
                        f"SELECT state, COUNT(*) FROM samples WHERE state IN ({placeholders}) GROUP BY state",
                        INFLIGHT_STATES,
                    )
                }
                self._conn.execute(
                    "UPDATE samples SET state = 'pending', "
                    f"error = '{STARTUP_INFLIGHT_RESET_PREFIX}[' || state || ']' "
                    "|| CASE WHEN error IS NULL OR error = '' THEN '' ELSE ': ' || substr(error, 1, 500) END, "
                    f"updated_at = datetime('now') WHERE state IN ({placeholders})",
                    INFLIGHT_STATES,
                )
                self._conn.commit()
            except BaseException:
                self._conn.rollback()
                raise
        if counts:
            logger.info(
                "Reset %d orphaned in-flight rows to pending at startup (%s)",
                sum(counts.values()),
                ", ".join(f"{state}={count}" for state, count in sorted(counts.items())),
            )
        return dict(sorted(counts.items()))

    def reset_environment_failures(self, classes: Iterable[str] = RETRYABLE_ENVIRONMENT_CLASSES) -> Dict[str, int]:
        """Return ``failed`` rows whose error is environmental to ``pending``; ``permanent_drop`` rows stay failed.

        An environmental failure (EPERM, EACCES, EROFS, ENOSPC, ENXIO, a
        missing tool) says nothing about the sample, so it must never stay a
        terminal state. The previous error is kept after a
        ``reset-environment-failure-at-startup[<class>]`` prefix. Returns counts
        by failure class. Call it only under the campaign writer lock.
        """

        wanted = frozenset(classes)
        unknown = wanted - set(SAMPLE_ERROR_CLASSES)
        if unknown:
            raise ValueError(f"unknown error classes: {sorted(unknown)}")
        counts: Dict[str, int] = {}
        with self._lock:
            self._ensure_retry_columns_locked()
            self._conn.commit()
            self._conn.execute("BEGIN IMMEDIATE")
            try:
                dropped = {
                    (species, srr_id)
                    for species, srr_id in self._conn.execute(
                        "SELECT species, srr_id FROM sample_exclusions WHERE reason_code = 'permanent_drop'"
                    )
                }
                updates: list[tuple[str, str, str]] = []
                for species, srr_id, error, recorded_class in self._conn.execute(
                    "SELECT species, srr_id, error, last_error_class FROM samples WHERE state = 'failed'"
                ).fetchall():
                    error_class = classify_sample_error(error)
                    if error_class not in wanted or (species, srr_id) in dropped:
                        continue
                    if recorded_class is not None and recorded_class not in RETRYABLE_ENVIRONMENT_CLASSES:
                        # The producer classified this attempt as not environmental when it happened (for
                        # example EPERM on a path outside the data root); re-reading its text would retry it
                        # at every start without ever counting an attempt.
                        continue
                    counts[error_class] = counts.get(error_class, 0) + 1
                    message = f"{STARTUP_ENVIRONMENT_RESET_PREFIX}[{error_class}]: {(error or '')[:500]}"
                    updates.append((message, species, srr_id))
                self._conn.executemany(
                    "UPDATE samples SET state = 'pending', error = ?, updated_at = datetime('now') "
                    "WHERE species = ? AND srr_id = ? AND state = 'failed'",
                    updates,
                )
                self._conn.commit()
            except BaseException:
                self._conn.rollback()
                raise
        if counts:
            logger.info(
                "Reset %d environment-class failures to pending at startup (%s)",
                sum(counts.values()),
                ", ".join(f"{name}={count}" for name, count in sorted(counts.items())),
            )
        return dict(sorted(counts.items()))

    # ---- Retry bookkeeping (attempts, failure class, backoff) ----

    def _retry_columns_present_locked(self) -> bool:
        if getattr(self, "_retry_columns_ready", False):
            return True
        present = {row[1] for row in self._conn.execute("PRAGMA table_info(samples)")}
        ready = all(name in present for name, _ in RETRY_COLUMNS)
        self._retry_columns_ready = ready
        return ready

    def _ensure_retry_columns_locked(self) -> None:
        """Add :data:`RETRY_COLUMNS` to ``samples`` once and backfill the failed rows (an additive migration).

        Readers that select named columns are unaffected. Each ``failed`` row
        gets the class of its stored error and one counted attempt, without a
        backoff: it has waited since it failed, so a transient failure is due
        at once, while a permanent or layout failure is held at once instead
        of spending another acquisition on it.
        """

        if self._retry_columns_present_locked():
            return
        self._conn.commit()
        self._conn.execute("BEGIN IMMEDIATE")
        try:
            present = {row[1] for row in self._conn.execute("PRAGMA table_info(samples)")}
            for name, column_type in RETRY_COLUMNS:
                if name not in present:
                    self._conn.execute(f"ALTER TABLE samples ADD COLUMN {name} {column_type}")
            failed = self._conn.execute(
                "SELECT species, srr_id, error FROM samples WHERE state = 'failed' AND last_error_class IS NULL"
            ).fetchall()
            self._conn.executemany(
                "UPDATE samples SET attempts = MAX(attempts, 1), last_error_class = ? WHERE species = ? AND srr_id = ?",
                [(classify_sample_error(error), species, srr_id) for species, srr_id, error in failed],
            )
            self._conn.commit()
        except BaseException:
            self._conn.rollback()
            raise
        self._retry_columns_ready = True
        if failed:
            logger.info("Retry bookkeeping added to %s; %d failed rows classified", self.db_path, len(failed))

    def retry_records(self, species: Optional[str] = None) -> Dict[tuple[str, str], RetryRecord]:
        """Rows that carry retry bookkeeping (attempts, a failure class or a backoff), keyed by ``(species, srr_id)``."""

        query = (
            "SELECT species, srr_id, state, attempts, last_error_class, next_attempt_at FROM samples "
            "WHERE (attempts > 0 OR last_error_class IS NOT NULL OR next_attempt_at IS NOT NULL)"
        )
        with self._lock:
            self._ensure_retry_columns_locked()
            rows = (
                self._conn.execute(query + " AND species = ?", (species,)).fetchall()
                if species is not None
                else self._conn.execute(query).fetchall()
            )
        return {
            (row[0], row[1]): RetryRecord(row[0], row[1], row[2], int(row[3] or 0), row[4], _utc_epoch(row[5]))
            for row in rows
        }

    def record_attempt_outcome(
        self,
        species: str,
        srr_id: str,
        error: Optional[str],
        *,
        policy: RetryPolicy,
        now: Optional[float] = None,
        error_class: Optional[str] = None,
    ) -> Dict[str, Any]:
        """Account one finished sample attempt and return the scheduling decision.

        ``error_class`` overrides the text classification of ``error`` (the
        producer passes :data:`NON_ENVIRONMENT_OS_ERROR` when its own classifier
        decided an environment-looking message is not environmental); the
        override is stored in ``last_error_class``, which the startup reset
        then trusts.

        ``error=None`` (quantified) clears the bookkeeping. Otherwise the error
        text is classified (:func:`classify_sample_error`, :func:`retry_class_for`):

        * ``environment``: no attempt is counted (the host failed, not the
          sample); decision ``retry_next_start``;
        * ``transient``: one more attempt; ``next_attempt_at`` is ``now`` plus
          the policy backoff (decision ``backoff``), or, once
          ``policy.max_attempts`` attempts are spent, decision ``held``;
        * ``layout_mismatch``: one more attempt, the row becomes
          ``quarantined`` with the error kept (decision ``held``): the inputs
          stay on disk for reconciliation and are never re-downloaded at
          every start;
        * ``permanent`` (kallisto exit 1): one more attempt, decision ``held``.

        The state written by the worker is otherwise kept. Returns
        ``{decision, error_class, retry_class, attempts, next_attempt_at}``
        (``next_attempt_at`` in POSIX seconds or ``None``).
        """

        moment = time.time() if now is None else now
        with self._lock:
            self._ensure_retry_columns_locked()
            row = self._conn.execute(
                "SELECT state, attempts FROM samples WHERE species = ? AND srr_id = ?", (species, srr_id)
            ).fetchone()
            if row is None:
                return {"decision": "missing", "error_class": None, "retry_class": None, "attempts": 0}
            if error is None:
                self._conn.execute(
                    "UPDATE samples SET attempts = 0, last_error_class = NULL, next_attempt_at = NULL "
                    "WHERE species = ? AND srr_id = ? "
                    "AND (attempts != 0 OR last_error_class IS NOT NULL OR next_attempt_at IS NOT NULL)",
                    (species, srr_id),
                )
                self._conn.commit()
                return {
                    "decision": "succeeded",
                    "error_class": None,
                    "retry_class": None,
                    "attempts": 0,
                    "next_attempt_at": None,
                }
            error_class = error_class or classify_sample_error(error)
            retry_class = retry_class_for(error_class)
            attempts = int(row[1] or 0)
            next_attempt: Optional[float] = None
            state = row[0]
            if retry_class == RETRY_ENVIRONMENT:
                decision = "retry_next_start"
            else:
                attempts += 1
                if retry_class == RETRY_TRANSIENT and attempts < policy.max_attempts:
                    # Rounded up: the stored text has whole seconds and must never shorten the backoff.
                    next_attempt = float(math.ceil(moment + policy.backoff_seconds(attempts)))
                    decision = "backoff"
                else:
                    decision = "held"
                if retry_class == RETRY_LAYOUT_MISMATCH:
                    state = "quarantined"
            self._conn.execute(
                "UPDATE samples SET attempts = ?, last_error_class = ?, next_attempt_at = ?, state = ?, "
                "error = CASE WHEN ? THEN ? ELSE error END, updated_at = datetime('now') "
                "WHERE species = ? AND srr_id = ?",
                (
                    attempts,
                    error_class,
                    _utc_text(next_attempt) if next_attempt is not None else None,
                    state,
                    state != row[0],
                    error[:2000],
                    species,
                    srr_id,
                ),
            )
            self._conn.commit()
        return {
            "decision": decision,
            "error_class": error_class,
            "retry_class": retry_class,
            "attempts": attempts,
            # Stored at whole-second resolution; report what a later start will read.
            "next_attempt_at": _utc_epoch(_utc_text(next_attempt)) if next_attempt is not None else None,
        }

    def release_retry_holds(
        self,
        *,
        policy: RetryPolicy,
        species: Optional[str] = None,
        accessions: Optional[Iterable[str]] = None,
        reasons: Iterable[str] = (RETRY_LAYOUT_MISMATCH, RETRY_PERMANENT, "attempts_exhausted"),
    ) -> Dict[str, int]:
        """Operator action: clear the bookkeeping of held rows so the next start schedules them again.

        A held ``failed`` row, and a ``quarantined`` row held for a layout
        mismatch, become ``pending`` with the previous error kept after a
        ``released-by-operator[<reason>]`` prefix. Quarantined rows whose
        output on disk failed the reuse contract are never touched. Only rows
        whose hold reason (``layout_mismatch``, ``permanent`` or
        ``attempts_exhausted``) is in ``reasons`` are released. Call it only
        under the campaign writer lock. Returns counts by hold reason.
        """

        wanted_reasons = frozenset(reasons)
        wanted_runs = None if accessions is None else {str(run) for run in accessions}
        counts: Dict[str, int] = {}
        moment = time.time()
        with self._lock:
            self._ensure_retry_columns_locked()
            self._conn.commit()
            self._conn.execute("BEGIN IMMEDIATE")
            try:
                query = (
                    "SELECT species, srr_id, state, attempts, last_error_class, next_attempt_at, error FROM samples "
                    "WHERE last_error_class IS NOT NULL"
                )
                parameters: tuple[Any, ...] = ()
                if species is not None:
                    query += " AND species = ?"
                    parameters = (species,)
                updates: list[tuple[Any, ...]] = []
                for row in self._conn.execute(query, parameters).fetchall():
                    if wanted_runs is not None and row[1] not in wanted_runs:
                        continue
                    record = RetryRecord(row[0], row[1], row[2], int(row[3] or 0), row[4], _utc_epoch(row[5]))
                    gate, reason = retry_gate(record, policy, moment)
                    kind = (reason or "").split(":", 1)[0]
                    if gate != GATE_HELD or kind not in wanted_reasons:
                        continue
                    releasable_state = row[2] == "failed" or (row[2] == "quarantined" and kind == RETRY_LAYOUT_MISMATCH)
                    new_state = "pending" if releasable_state else row[2]
                    message = (
                        f"{RETRY_RELEASED_PREFIX}[{reason}]: {(row[6] or '')[:500]}" if releasable_state else row[6]
                    )
                    updates.append((new_state, message, row[0], row[1]))
                    counts[kind] = counts.get(kind, 0) + 1
                self._conn.executemany(
                    "UPDATE samples SET attempts = 0, last_error_class = NULL, next_attempt_at = NULL, state = ?, "
                    "error = ?, updated_at = datetime('now') WHERE species = ? AND srr_id = ?",
                    updates,
                )
                self._conn.commit()
            except BaseException:
                self._conn.rollback()
                raise
        if counts:
            logger.info(
                "Released %d held rows (%s)",
                sum(counts.values()),
                ", ".join(f"{name}={count}" for name, count in sorted(counts.items())),
            )
        return dict(sorted(counts.items()))

    def states_for_species(self, species: str) -> List[tuple[str, str]]:
        """Every ``(srr_id, state)`` recorded for ``species``, in accession order."""

        with self._lock:
            return [
                (srr_id, state)
                for srr_id, state in self._conn.execute(
                    "SELECT srr_id, state FROM samples WHERE species = ? ORDER BY srr_id", (species,)
                )
            ]

    def states_for(self, keys: Iterable[tuple[str, str]]) -> Dict[tuple[str, str], Optional[str]]:
        """The current state of each ``(species, srr_id)`` (``None`` when the row does not exist)."""

        wanted = {(str(species), str(srr_id)) for species, srr_id in keys}
        species_names = sorted({species for species, _ in wanted})
        found: Dict[tuple[str, str], Optional[str]] = dict.fromkeys(wanted)
        with self._lock:
            for species in species_names:
                for srr_id, state in self._conn.execute(
                    "SELECT srr_id, state FROM samples WHERE species = ?", (species,)
                ):
                    if (species, srr_id) in found:
                        found[(species, srr_id)] = state
        return found

    def checkpoint_wal(self, mode: str = "TRUNCATE") -> tuple[int, int, int] | None:
        """``PRAGMA wal_checkpoint(<mode>)`` so read-only observers of the file see every commit.

        ``TRUNCATE`` (the default, at exit) also empties the WAL file; the
        producer's periodic in-run checkpoint uses ``PASSIVE``, which never
        waits for readers or writers.
        """

        if mode not in ("PASSIVE", "FULL", "RESTART", "TRUNCATE"):
            raise ValueError(f"unknown WAL checkpoint mode: {mode!r}")
        with self._lock:
            self._conn.commit()
            try:
                row = self._conn.execute(f"PRAGMA wal_checkpoint({mode})").fetchone()
            except sqlite3.OperationalError as exc:
                logger.warning("WAL checkpoint of %s failed: %s", self.db_path, exc)
                return None
        return tuple(row) if row is not None else None  # type: ignore[return-value]

    # ---- Reconciliation ----

    def reconcile(
        self,
        species: str,
        quant_dir: Path,
        *,
        verify_hashes: bool = False,
        config_path: Path | None = None,
        reference_manifest_path: Path | None = None,
        expected_config_sha256: str | None = None,
        index_digests: IndexDigestCache | None = None,
    ) -> int:
        """Reconcile only provenance-qualified reusable quantification outputs.

        Reusable means ``current`` or ``version_drift_compatible``: an output
        whose contract verifies but which a different Amalgkit build produced
        (for example every 0.16.60 output after the 0.16.108 pin) is promoted,
        exactly as :meth:`reconcile_outputs` does, so an Amalgkit bump never
        turns valid work back into pending tasks. Readable abundance files
        without a reusable provenance sidecar remain visible on disk but are
        not promoted to ``quantified`` in the DB.

        Args:
            species: Species name.
            quant_dir: Path to the quant output directory
                       (e.g. ``output/amalgkit/amellifera/work/quant``).

        Returns:
            Number of samples reconciled.

        ``verify_hashes`` defaults to False: the resume-time filesystem scan
        classifies by provenance contract (sidecar completeness, accession
        binding, output presence) without re-digesting every quantification
        payload. Re-hashing the full quantified corpus (multi-GB for a large
        species) inside Phase-1 discovery stalled task submission for longer
        than the campaign stall watchdogs on every restart. Pass True only for
        an explicit deep audit.

        ``expected_config_sha256`` binds the configuration by content (the
        sha256 of the current runtime configuration, computed once per
        species) instead of re-reading the path recorded in every sidecar, so
        outputs whose original checkout moved or was deleted still reconcile.
        One ``index_digests`` cache per species avoids re-hashing the shared
        reference manifest for every sample.
        """
        if not quant_dir.exists():
            return 0
        digests = index_digests if index_digests is not None else IndexDigestCache()

        reconciled = []
        rejected = []
        try:
            for subdir in quant_dir.iterdir():
                if not subdir.is_dir():
                    continue
                if find_quantification_file(subdir, subdir.name) is None:
                    continue
                classification = classify_quantification(
                    subdir,
                    subdir.name,
                    verify_content=verify_hashes,
                    expected_config_path=config_path,
                    expected_reference_manifest_path=reference_manifest_path,
                    expected_config_sha256=expected_config_sha256,
                    index_digests=digests,
                )
                if classification.get("status") in REUSABLE_QUANT_STATUSES:
                    reconciled.append(subdir.name)
                else:
                    rejected.append((subdir.name, classification.get("status", "unknown")))
        except Exception as e:
            logger.warning(f"Reconciliation scan failed for {species}: {e}")
            return 0

        if reconciled:
            self.bulk_set_state(species, reconciled, "quantified")
            logger.info(f"Reconciled {len(reconciled)} quantified samples for {species}")
        if rejected:
            logger.info(
                "Left %d readable but non-reusable quantifications visible for %s: %s",
                len(rejected),
                species,
                ", ".join(f"{sample}={status}" for sample, status in rejected[:10]),
            )

        return len(reconciled)

    # ---- Frozen cohort binding ----

    def _ensure_meta_table(self) -> None:
        self._conn.execute(
            "CREATE TABLE IF NOT EXISTS campaign_meta ("
            " key TEXT PRIMARY KEY, value TEXT NOT NULL, updated_at TEXT NOT NULL DEFAULT (datetime('now')))"
        )

    def cohort_binding(self) -> Dict[str, str]:
        """Return the frozen-cohort binding recorded by :meth:`migrate_to_cohort` (empty when unbound)."""

        with self._lock:
            present = self._conn.execute(
                "SELECT 1 FROM sqlite_master WHERE type = 'table' AND name = 'campaign_meta'"
            ).fetchone()
            if not present:
                return {}
            rows = self._conn.execute("SELECT key, value FROM campaign_meta").fetchall()
        return {key: value for key, value in rows}

    def _state_counts_locked(self) -> Dict[str, int]:
        counts: Dict[str, int] = {}
        for state, count in self._conn.execute("SELECT state, COUNT(*) FROM samples GROUP BY state"):
            counts[state] = count
        return counts

    def migrate_to_cohort(
        self,
        cohort: "FrozenCohort",
        *,
        backup_path: Path | None = None,
        reset_classes: Iterable[str] = RETRYABLE_ENVIRONMENT_CLASSES,
        extra_exclusions: Iterable[Mapping[str, Any]] = (),
        recorded_by: str | None = None,
    ) -> Dict[str, Any]:
        """Bind this progress store to a frozen cohort without losing any recorded state.

        In one ``BEGIN IMMEDIATE`` transaction this:

        * adds a ``pending`` row for every cohort task that has none;
        * records every cohort exclusion (and each operator-supplied
          ``extra_exclusions`` entry) as ``permanent_drop`` when the accession
          has no exclusion yet, leaving existing exclusions untouched;
        * resets ``failed`` rows of cohort tasks whose error is in
          ``reset_classes`` to ``pending``, keeping the old message after a
          ``reset-by-cohort-migration[<class>]`` prefix;
        * records the cohort id and definition digest in ``campaign_meta``.

        No row is deleted or otherwise changed: quantified rows, exclusions,
        quarantines, in-flight rows and rows outside the cohort keep their state
        and are reported. A store already bound to a different cohort is
        refused. With ``backup_path`` the database is first copied there with
        SQLite's online backup API (an existing backup is never overwritten).
        Re-running the migration is idempotent.
        """

        from metainformant.rna.engine.cohort import CohortError, PERMANENT_DROP

        resets = frozenset(reset_classes)
        unknown = resets - set(SAMPLE_ERROR_CLASSES)
        if unknown:
            raise ValueError(f"unknown reset classes: {sorted(unknown)}")
        binding = self.cohort_binding()
        if binding.get("cohort_id") not in (None, cohort.cohort_id):
            raise CohortError(f"progress store is bound to cohort {binding['cohort_id']}, not {cohort.cohort_id}")
        if binding.get("cohort_definition_sha256") not in (None, cohort.definition_sha256):
            raise CohortError("progress store is bound to a different definition of this cohort")
        recorder = recorded_by or f"cohort:{cohort.cohort_id}"
        extra_rows: list[tuple[str, str, str, str | None, str]] = []
        for entry in extra_exclusions:
            species, accession = str(entry["species"]), str(entry["srr_id"])
            reason_code = str(entry.get("reason_code") or PERMANENT_DROP)
            if reason_code not in EXCLUSION_REASON_CODES:
                raise ValueError(f"Invalid exclusion reason code '{reason_code}'")
            if cohort.task(species, accession) is None and not cohort.is_excluded(species, accession):
                raise CohortError(f"extra exclusion is not part of cohort {cohort.cohort_id}: {species}/{accession}")
            detail = entry.get("reason_detail")
            extra_rows.append(
                (
                    species,
                    accession,
                    reason_code,
                    str(detail) if detail else None,
                    str(entry.get("recorded_by") or recorder),
                )
            )

        backup: str | None = None
        if backup_path is not None:
            target = Path(backup_path)
            if target.exists() or target.is_symlink():
                if binding.get("cohort_id") != cohort.cohort_id:
                    raise FileExistsError(f"backup already exists: {target}")
                backup = f"kept existing {target}"
            else:
                target.parent.mkdir(parents=True, exist_ok=True)
                with self._lock, sqlite3.connect(str(target)) as copy:
                    self._conn.backup(copy)
                copy.close()
                backup = str(target)

        task_rows = [(task.species, task.accession) for task in cohort.tasks]
        task_keys = set(task_rows)
        excluded_keys = {(species, accession) for species, values in cohort.exclusions.items() for accession in values}
        cohort_exclusion_rows = [
            (species, accession, reason, f"frozen cohort {cohort.cohort_id}", recorder)
            for species, values in sorted(cohort.exclusions.items())
            for accession, reason in sorted(values.items())
        ]
        summary: Dict[str, Any] = {"cohort_id": cohort.cohort_id, "definition_sha256": cohort.definition_sha256}
        with self._lock:
            self._conn.commit()
            self._conn.execute("BEGIN IMMEDIATE")
            try:
                self._ensure_meta_table()
                rows_before = self._conn.execute("SELECT COUNT(*) FROM samples").fetchone()[0]
                states_before = self._state_counts_locked()
                cursor = self._conn.executemany(
                    "INSERT OR IGNORE INTO samples (species, srr_id, state) VALUES (?, ?, 'pending')", task_rows
                )
                inserted = cursor.rowcount

                existing_exclusions = {
                    (species, accession): reason
                    for species, accession, reason in self._conn.execute(
                        "SELECT species, srr_id, reason_code FROM sample_exclusions"
                    )
                }
                conflicts = sorted(
                    f"{species}/{accession}={existing_exclusions[(species, accession)]}"
                    for species, accession, reason, _, _ in cohort_exclusion_rows + extra_rows
                    if (species, accession) in existing_exclusions
                    and existing_exclusions[(species, accession)] != reason
                )
                insert_exclusion = (
                    "INSERT OR IGNORE INTO sample_exclusions "
                    "(species, srr_id, reason_code, reason_detail, recorded_by, recorded_at) "
                    "VALUES (?, ?, ?, ?, ?, datetime('now'))"
                )
                cohort_inserted = self._conn.executemany(insert_exclusion, cohort_exclusion_rows).rowcount
                extra_inserted = self._conn.executemany(insert_exclusion, extra_rows).rowcount if extra_rows else 0
                dropped = {
                    (species, accession)
                    for species, accession, reason in self._conn.execute(
                        "SELECT species, srr_id, reason_code FROM sample_exclusions"
                    )
                    if reason == PERMANENT_DROP
                }

                reset_by_class: Dict[str, int] = {}
                failed_kept_by_class: Dict[str, int] = {}
                updates: list[tuple[str, str, str]] = []
                for species, accession, error in self._conn.execute(
                    "SELECT species, srr_id, error FROM samples WHERE state = 'failed'"
                ).fetchall():
                    error_class = classify_sample_error(error)
                    key = (species, accession)
                    if key in task_keys and key not in dropped and error_class in resets:
                        reset_by_class[error_class] = reset_by_class.get(error_class, 0) + 1
                        message = f"{COHORT_MIGRATION_ERROR_PREFIX}[{error_class}]: {(error or '')[:500]}"
                        updates.append((message, species, accession))
                    else:
                        failed_kept_by_class[error_class] = failed_kept_by_class.get(error_class, 0) + 1
                self._conn.executemany(
                    "UPDATE samples SET state = 'pending', error = ?, updated_at = datetime('now') "
                    "WHERE species = ? AND srr_id = ? AND state = 'failed'",
                    updates,
                )

                outside: Dict[str, int] = {}
                excluded_states: Dict[str, int] = {}
                in_flight: Dict[str, int] = {}
                for species, accession, state in self._conn.execute("SELECT species, srr_id, state FROM samples"):
                    key = (species, accession)
                    if key in excluded_keys or key in dropped:
                        excluded_states[state] = excluded_states.get(state, 0) + 1
                    elif key not in task_keys:
                        outside[state] = outside.get(state, 0) + 1
                    elif state in {"downloading", "downloaded", "quantifying"}:
                        in_flight[state] = in_flight.get(state, 0) + 1
                missing = [
                    f"{species}/{accession}"
                    for species, accession in task_rows
                    if not self._conn.execute(
                        "SELECT 1 FROM samples WHERE species = ? AND srr_id = ?", (species, accession)
                    ).fetchone()
                ]
                if missing:
                    raise RuntimeError(f"cohort tasks still lack progress rows: {missing[:5]}")
                for meta_key, meta_value in (
                    ("cohort_id", cohort.cohort_id),
                    ("cohort_definition_sha256", cohort.definition_sha256),
                ):
                    self._conn.execute(
                        "INSERT INTO campaign_meta (key, value, updated_at) VALUES (?, ?, datetime('now')) "
                        "ON CONFLICT(key) DO UPDATE SET value = excluded.value, updated_at = excluded.updated_at "
                        "WHERE campaign_meta.value <> excluded.value",
                        (meta_key, meta_value),
                    )
                rows_after = self._conn.execute("SELECT COUNT(*) FROM samples").fetchone()[0]
                states_after = self._state_counts_locked()
                if rows_after != rows_before + inserted:
                    raise RuntimeError("cohort migration changed the row count unexpectedly")
                self._conn.commit()
            except BaseException:
                self._conn.rollback()
                raise
        summary.update(
            backup=backup,
            rows_before=rows_before,
            rows_after=rows_after,
            states_before=dict(sorted(states_before.items())),
            states_after=dict(sorted(states_after.items())),
            inserted_task_rows=inserted,
            cohort_exclusions_inserted=cohort_inserted,
            extra_exclusions_inserted=extra_inserted,
            exclusion_reason_conflicts=conflicts,
            reset_by_class=dict(sorted(reset_by_class.items())),
            failed_kept_by_class=dict(sorted(failed_kept_by_class.items())),
            excluded_rows_by_state=dict(sorted(excluded_states.items())),
            in_flight_task_rows=dict(sorted(in_flight.items())),
            outside_cohort_rows=dict(sorted(outside.items())),
        )
        logger.info(
            "migrate_to_cohort(%s): +%d task rows, reset %d failed rows, %d rows outside the cohort kept",
            cohort.cohort_id,
            inserted,
            sum(reset_by_class.values()),
            sum(outside.values()),
        )
        return summary

    def ensure_exclusions(self, species: str, entries: Iterable[Mapping[str, Any]]) -> int:
        """Insert exclusions that are not recorded yet; never change an existing one. Returns rows added."""

        rows = []
        for entry in entries:
            reason_code = str(entry.get("reason_code") or "permanent_drop")
            if reason_code not in EXCLUSION_REASON_CODES:
                raise ValueError(f"Invalid exclusion reason code '{reason_code}'")
            rows.append(
                (species, str(entry["srr_id"]), reason_code, entry.get("reason_detail"), entry.get("recorded_by"))
            )
        if not rows:
            return 0
        with self._lock:
            cursor = self._conn.executemany(
                "INSERT OR IGNORE INTO sample_exclusions "
                "(species, srr_id, reason_code, reason_detail, recorded_by, recorded_at) "
                "VALUES (?, ?, ?, ?, ?, datetime('now'))",
                rows,
            )
            self._conn.commit()
        return cursor.rowcount

    def _upsert_audit_locked(self, species: str, srr_id: str, audit: Mapping[str, Any]) -> None:
        self._conn.execute(
            """
            INSERT INTO quantification_audit (
                species, srr_id, status, reason, contract_id,
                observed_amalgkit_version, observed_release_tag,
                observed_source_revision, checked_at
            ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, datetime('now'))
            ON CONFLICT(species, srr_id) DO UPDATE SET
                status = excluded.status,
                reason = excluded.reason,
                contract_id = excluded.contract_id,
                observed_amalgkit_version = excluded.observed_amalgkit_version,
                observed_release_tag = excluded.observed_release_tag,
                observed_source_revision = excluded.observed_source_revision,
                checked_at = excluded.checked_at
            """,
            (
                species,
                srr_id,
                str(audit.get("status")),
                str(audit.get("reason")),
                audit.get("contract_id"),
                audit.get("observed_amalgkit_version"),
                audit.get("observed_release_tag"),
                audit.get("observed_source_revision"),
            ),
        )
        self._clear_audit_fingerprint_locked(species, srr_id)

    def adopt_quantified(self, species: str, entries: Iterable[Mapping[str, Any]]) -> Dict[str, List[str]]:
        """Mark verified imported outputs ``quantified`` (with their audit) in one transaction.

        Each entry carries ``srr_id`` plus the classification fields (``status``,
        ``reason``, ``contract_id``, ``observed_*``); ``status`` must be
        reusable. Accessions with a ``permanent_drop`` exclusion are left
        unchanged and reported.
        """

        items = [dict(entry) for entry in entries]
        for item in items:
            if item.get("status") not in REUSABLE_QUANT_STATUSES:
                raise ValueError(f"cannot adopt a non-reusable quantification: {item.get('srr_id')}")
        adopted: List[str] = []
        excluded: List[str] = []
        with self._lock:
            self._conn.commit()
            self._conn.execute("BEGIN IMMEDIATE")
            try:
                dropped = {
                    row[0]
                    for row in self._conn.execute(
                        "SELECT srr_id FROM sample_exclusions WHERE species = ? AND reason_code = 'permanent_drop'",
                        (species,),
                    )
                }
                for item in items:
                    srr_id = str(item["srr_id"])
                    if srr_id in dropped:
                        excluded.append(srr_id)
                        continue
                    self._conn.execute(
                        "INSERT OR IGNORE INTO samples (species, srr_id, state) VALUES (?, ?, 'pending')",
                        (species, srr_id),
                    )
                    self._upsert_audit_locked(species, srr_id, item)
                    self._conn.execute(
                        "UPDATE samples SET state = 'quantified', error = NULL, updated_at = datetime('now') "
                        "WHERE species = ? AND srr_id = ?",
                        (species, srr_id),
                    )
                    adopted.append(srr_id)
                self._conn.commit()
            except BaseException:
                self._conn.rollback()
                raise
        return {"adopted": adopted, "excluded": excluded}

    def _audit_fingerprint_columns_present_locked(self) -> bool:
        if getattr(self, "_audit_fingerprint_ready", False):
            return True
        present = {row[1] for row in self._conn.execute("PRAGMA table_info(quantification_audit)")}
        ready = all(name in present for name, _ in AUDIT_FINGERPRINT_COLUMNS)
        self._audit_fingerprint_ready = ready
        return ready

    def _ensure_audit_fingerprint_columns_locked(self) -> None:
        """Add the nullable fingerprint columns to ``quantification_audit`` once (an additive migration)."""

        if self._audit_fingerprint_columns_present_locked():
            return
        present = {row[1] for row in self._conn.execute("PRAGMA table_info(quantification_audit)")}
        for name, column_type in AUDIT_FINGERPRINT_COLUMNS:
            if name not in present:
                self._conn.execute(f"ALTER TABLE quantification_audit ADD COLUMN {name} {column_type}")
        self._audit_fingerprint_ready = True

    def _clear_audit_fingerprint_locked(self, species: str, srr_id: str) -> None:
        """An audit written without a fingerprint must never be trusted by a later reconcile."""

        if self._audit_fingerprint_columns_present_locked():
            self._conn.execute(
                "UPDATE quantification_audit SET sidecar_ino = NULL, sidecar_size = NULL, sidecar_mtime_ns = NULL, "
                "config_sha256 = NULL, reference_sha256 = NULL, scope = NULL, evidence = NULL, "
                "output_fingerprint = NULL WHERE species = ? AND srr_id = ?",
                (species, srr_id),
            )

    def reconcile_outputs(
        self,
        species: str,
        quant_dir: Path,
        accessions: Iterable[str],
        classify: Callable[[Path, str], Mapping[str, Any]],
        *,
        scope: ReconcileScope | None = None,
        requantify_statuses: Iterable[str] = (),
        deep: bool = False,
        workers: int = 1,
        sidecar_name: str = QUANT_PROVENANCE_FILENAME,
        heartbeat_every: int = RECONCILE_HEARTBEAT_EVERY,
    ) -> Dict[str, Any]:
        """Bring task states in line with the outputs on disk in one pass, proportional to what changed.

        Only ``accessions`` are considered. The quantification directory is
        listed once (one ``scandir``). For each listed task directory:

        * **Unchanged, already judged.** When the row is recorded with an
          audit whose sidecar fingerprint ``(st_ino, st_size, st_mtime_ns)``
          equals the sidecar on disk (one ``lstat``), whose abundance-table
          fingerprint equals the table on disk (one more ``lstat``) and whose
          ``config_sha256``/``reference_sha256``/``scope`` equal ``scope``,
          the recorded verdict is reused without reading the sidecar.
        * **New or changed.** Otherwise ``classify(path, accession)`` judges
          it and the verdict is recorded with both fingerprints.
        * **Unreadable.** ``classify`` runs under
          :func:`~metainformant.rna.engine.provenance.strict_reads`; an output
          whose files cannot be read (EIO, EACCES, a dropped volume) gets no
          verdict: its row and audit are left as they were, nothing is cached,
          and it is listed under ``unavailable`` so the caller can report an
          incomplete pass. The next pass judges it again.

        A reusable verdict (current or version drift) makes the row
        ``quantified`` (``pending`` when its status is in
        ``requantify_statuses``); a non-reusable verdict for a row recorded
        ``quantified`` makes it ``quarantined`` (the bytes stay on disk); a
        row recorded ``quantified`` without a directory becomes ``pending``.
        Other rows with a non-reusable output are left alone (the producer
        re-quantifies them).

        All writes happen in one transaction, and only for rows whose state,
        error or audit actually changes: an unchanged row keeps its
        ``updated_at`` and its audit ``checked_at``. ``deep=True`` (the
        ``AMALGKIT_RESUME_VERIFY_HASHES`` audit) re-judges every directory.
        Without ``scope`` nothing is trusted (and no fingerprint is recorded).
        A heartbeat is logged every ``heartbeat_every`` directories.
        """

        started = time.monotonic()
        wanted = {str(accession) for accession in accessions}
        requeue = frozenset(requantify_statuses)
        directories: Dict[str, Path] = {}
        listed = 0
        if quant_dir.is_dir():
            with os.scandir(quant_dir) as entries:
                for entry in entries:
                    listed += 1
                    if entry.name in wanted and entry.is_dir(follow_symlinks=False):
                        directories[entry.name] = Path(entry.path)
        with self._lock:
            self._ensure_audit_fingerprint_columns_locked()
            self._conn.commit()
            states: Dict[str, str] = {}
            errors: Dict[str, Optional[str]] = {}
            for srr_id, state, error in self._conn.execute(
                "SELECT srr_id, state, error FROM samples WHERE species = ?", (species,)
            ):
                if srr_id in wanted:
                    states[srr_id] = state
                    errors[srr_id] = error
            audits: Dict[str, Dict[str, Any]] = {}
            if scope is not None and not deep:
                for row in self._conn.execute(
                    "SELECT srr_id, status, reason, contract_id, observed_amalgkit_version, observed_release_tag, "
                    "observed_source_revision, sidecar_ino, sidecar_size, sidecar_mtime_ns, config_sha256, "
                    "reference_sha256, scope, evidence, output_fingerprint FROM quantification_audit "
                    "WHERE species = ? AND sidecar_ino IS NOT NULL AND output_fingerprint IS NOT NULL "
                    "AND config_sha256 = ? AND reference_sha256 = ? AND scope = ?",
                    (species, scope.config_sha256, scope.reference_sha256, scope.scope),
                ):
                    if row[0] in directories:
                        audits[row[0]] = {
                            "status": row[1],
                            "reason": row[2],
                            "contract_id": row[3],
                            "observed_amalgkit_version": row[4],
                            "observed_release_tag": row[5],
                            "observed_source_revision": row[6],
                            "fingerprint": (row[7], row[8], row[9]),
                            "evidence": json.loads(row[13]) if row[13] else None,
                            "output_fingerprint": row[14],
                        }

        def judge(item: tuple[str, Path]) -> tuple[str, Mapping[str, Any], bool, SidecarFingerprint | None]:
            accession, path = item
            sidecar = path / sidecar_name
            recorded = audits.get(accession)
            if (
                recorded is not None
                and sidecar_fingerprint(sidecar) == recorded["fingerprint"]
                and _recorded_output_matches(path, recorded["output_fingerprint"])
            ):
                return accession, recorded, False, recorded["fingerprint"]
            before = sidecar_fingerprint(sidecar) if scope is not None else None
            try:
                with strict_reads():
                    verdict = dict(classify(path, accession))
            except (OSError, ValueError) as error:
                if is_unreadable_error(error):
                    # No verdict: an I/O or permission failure says nothing about the output's content.
                    return accession, {"status": QUANT_STATUS_UNAVAILABLE, "reason": f"unreadable: {error}"}, True, None
                verdict = {"status": "invalid", "reason": f"classification failed: {error}"}
            after = sidecar_fingerprint(sidecar) if scope is not None else None
            if scope is not None:
                verdict["output_fingerprint"] = output_fingerprint(path, verdict.get("quantification_file"))
            return accession, verdict, True, before if before is not None and before == after else None

        pairs = sorted(directories.items())
        judged: list[tuple[str, Mapping[str, Any], bool, SidecarFingerprint | None]] = []
        fresh_count = 0

        def heartbeat(done: int) -> None:
            if heartbeat_every > 0 and done % heartbeat_every == 0:
                logger.info(
                    "reconcile(%s): %d/%d output directories checked (%d re-judged) in %.1fs",
                    species,
                    done,
                    len(pairs),
                    fresh_count,
                    time.monotonic() - started,
                )

        if workers > 1 and len(pairs) > 1:
            with concurrent.futures.ThreadPoolExecutor(max_workers=min(workers, len(pairs))) as pool:
                for result in pool.map(judge, pairs, chunksize=32):
                    judged.append(result)
                    fresh_count += result[2]
                    heartbeat(len(judged))
        else:
            for pair in pairs:
                result = judge(pair)
                judged.append(result)
                fresh_count += result[2]
                heartbeat(len(judged))

        audit_rows: list[tuple[Any, ...]] = []
        state_rows: list[tuple[str, Optional[str], str, str, str, Optional[str]]] = []
        reusable: list[str] = []
        requeued: list[str] = []
        quarantined: list[str] = []
        newly_quarantined: list[str] = []
        unavailable: list[str] = []
        for accession, verdict, fresh, fingerprint in judged:
            status = verdict.get("status")
            if status == QUANT_STATUS_UNAVAILABLE:
                unavailable.append(accession)
                logger.warning("reconcile(%s): %s left as %s: %s", species, accession, states.get(accession), verdict)
                continue
            if fresh:
                evidence = verdict.get("evidence")
                keyed = scope is not None and fingerprint is not None and verdict.get("output_fingerprint") is not None
                audit_rows.append(
                    (
                        species,
                        accession,
                        str(status),
                        str(verdict.get("reason")),
                        verdict.get("contract_id"),
                        verdict.get("observed_amalgkit_version"),
                        verdict.get("observed_release_tag"),
                        verdict.get("observed_source_revision"),
                        fingerprint[0] if keyed else None,
                        fingerprint[1] if keyed else None,
                        fingerprint[2] if keyed else None,
                        scope.config_sha256 if keyed else None,
                        scope.reference_sha256 if keyed else None,
                        scope.scope if keyed else None,
                        json.dumps(evidence, sort_keys=True, default=str) if evidence else None,
                        verdict.get("output_fingerprint") if keyed else None,
                    )
                )
            target: Optional[str] = None
            error: Optional[str] = None
            if status in REUSABLE_QUANT_STATUSES:
                if status in requeue:
                    target, error = "pending", f"requantify: {status}"
                    requeued.append(accession)
                else:
                    target = "quantified"
                    reusable.append(accession)
            elif states.get(accession) == "quantified":
                target, error = "quarantined", str(verdict.get("reason"))
                quarantined.append(accession)
                newly_quarantined.append(accession)
            elif states.get(accession) == "quarantined":
                quarantined.append(accession)
            if target is not None and (states.get(accession), errors.get(accession)) != (target, error):
                state_rows.append((target, error, species, accession, target, error))
        vanished = sorted(
            accession for accession, state in states.items() if state == "quantified" and accession not in directories
        )
        for accession in vanished:
            error = "reconcile: recorded quantification has no output directory"
            state_rows.append(("pending", error, species, accession, "pending", error))

        if audit_rows or state_rows:
            with self._lock:
                self._conn.commit()
                self._conn.execute("BEGIN IMMEDIATE")
                try:
                    self._conn.executemany(
                        "INSERT OR IGNORE INTO samples (species, srr_id, state) VALUES (?, ?, 'pending')",
                        [(row[2], row[3]) for row in state_rows],
                    )
                    self._conn.executemany(
                        """
                        INSERT INTO quantification_audit (
                            species, srr_id, status, reason, contract_id,
                            observed_amalgkit_version, observed_release_tag, observed_source_revision,
                            sidecar_ino, sidecar_size, sidecar_mtime_ns, config_sha256, reference_sha256,
                            scope, evidence, output_fingerprint, checked_at
                        ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, datetime('now'))
                        ON CONFLICT(species, srr_id) DO UPDATE SET
                            status = excluded.status,
                            reason = excluded.reason,
                            contract_id = excluded.contract_id,
                            observed_amalgkit_version = excluded.observed_amalgkit_version,
                            observed_release_tag = excluded.observed_release_tag,
                            observed_source_revision = excluded.observed_source_revision,
                            sidecar_ino = excluded.sidecar_ino,
                            sidecar_size = excluded.sidecar_size,
                            sidecar_mtime_ns = excluded.sidecar_mtime_ns,
                            config_sha256 = excluded.config_sha256,
                            reference_sha256 = excluded.reference_sha256,
                            scope = excluded.scope,
                            evidence = excluded.evidence,
                            output_fingerprint = excluded.output_fingerprint,
                            checked_at = excluded.checked_at
                        """,
                        audit_rows,
                    )
                    self._conn.executemany(
                        "UPDATE samples SET state = ?, error = ?, updated_at = datetime('now') "
                        "WHERE species = ? AND srr_id = ? AND NOT (state = ? AND error IS ?)",
                        state_rows,
                    )
                    self._conn.commit()
                except BaseException:
                    self._conn.rollback()
                    raise
        summary: Dict[str, Any] = {
            "listed": listed,
            "scanned": len(judged),
            "trusted": len(judged) - fresh_count,
            "classified": fresh_count,
            "quantified": len(reusable),
            "requeued": len(requeued),
            "requeued_accessions": sorted(requeued),
            "unavailable": sorted(unavailable),
            "quarantined": sorted(quarantined),
            "newly_quarantined": sorted(newly_quarantined),
            "vanished": vanished,
            "rows_changed": len(state_rows),
            "audits_written": len(audit_rows),
            "seconds": round(time.monotonic() - started, 3),
        }
        logger.info(
            "reconcile(%s): %d task directories (%d trusted unchanged, %d re-judged), %d reusable, %d requeued, "
            "%d quarantined, %d vanished, %d unreadable (left unjudged), %d rows changed in %.1fs",
            species,
            summary["scanned"],
            summary["trusted"],
            summary["classified"],
            summary["quantified"],
            summary["requeued"],
            len(summary["quarantined"]),
            len(vanished),
            len(unavailable),
            summary["rows_changed"],
            summary["seconds"],
        )
        return summary

    def reconcile_cohort(
        self,
        species: str,
        quant_dir: Path,
        accessions: Iterable[str],
        classify: Callable[[Path, str], Mapping[str, Any]],
        *,
        requantify_statuses: Iterable[str] = (),
        workers: int = 1,
        scope: ReconcileScope | None = None,
        deep: bool = False,
    ) -> Dict[str, Any]:
        """Cohort reconcile (:meth:`reconcile_outputs`) reporting ``scanned``/``quantified``/``requeued``/
        ``quarantined``/``vanished``; ``quarantined`` lists every task run whose output is quarantined."""

        summary = self.reconcile_outputs(
            species,
            quant_dir,
            accessions,
            classify,
            scope=scope,
            requantify_statuses=requantify_statuses,
            deep=deep,
            workers=workers,
        )
        return {key: summary[key] for key in ("scanned", "quantified", "requeued", "quarantined", "vanished")}

    # ---- Lifecycle ----

    def close(self) -> None:
        """Close the database connection."""
        self._conn.close()

    def __enter__(self) -> "ProgressDB":
        return self

    def __exit__(self, exc_type: object, exc: object, tb: object) -> None:
        self.close()

    def __repr__(self) -> str:
        counts = self.get_total_counts()
        total = sum(counts.values())
        return f"ProgressDB({self.db_path}, {len(counts)} species, {total} samples)"


# ---------- Module-level convenience ----------

_db_instance: Optional[ProgressDB] = None
_db_lock = threading.Lock()


def get_db(db_path: Optional[Path] = None) -> ProgressDB:
    """Get or create the global ProgressDB singleton.

    Args:
        db_path: Custom path. If provided, creates a new non-singleton instance.

    Returns:
        A ProgressDB instance.
    """
    global _db_instance
    if db_path:
        return ProgressDB(db_path)

    with _db_lock:
        if _db_instance is None:
            _db_instance = ProgressDB()
        return _db_instance


# ---------- Retry bookkeeping CLI ----------


def retry_status(db_path: Path, *, policy: Optional[RetryPolicy] = None, now: Optional[float] = None) -> Dict[str, Any]:
    """Read-only summary of the retry bookkeeping: held rows by reason and rows in backoff.

    Opens the database with ``mode=ro`` (no schema creation, no migration), so
    it is safe while a producer runs. A database without the retry columns
    reports nothing held.
    """

    active = policy or RetryPolicy.from_environment()
    moment = time.time() if now is None else now
    held: Dict[str, int] = {}
    backoff = 0
    earliest: Optional[float] = None
    connection = sqlite3.connect(f"file:{Path(db_path).absolute()}?mode=ro", uri=True)
    try:
        columns = {row[1] for row in connection.execute("PRAGMA table_info(samples)")}
        rows = (
            connection.execute(
                "SELECT species, srr_id, state, attempts, last_error_class, next_attempt_at FROM samples "
                "WHERE last_error_class IS NOT NULL OR next_attempt_at IS NOT NULL"
            ).fetchall()
            if all(name in columns for name, _ in RETRY_COLUMNS)
            else []
        )
    finally:
        connection.close()
    for row in rows:
        record = RetryRecord(row[0], row[1], row[2], int(row[3] or 0), row[4], _utc_epoch(row[5]))
        gate, reason = retry_gate(record, active, moment)
        if gate == GATE_HELD:
            kind = (reason or "held").split(":", 1)[0]
            held[kind] = held.get(kind, 0) + 1
        elif gate == GATE_BACKOFF:
            backoff += 1
            if record.next_attempt_at is not None:
                earliest = record.next_attempt_at if earliest is None else min(earliest, record.next_attempt_at)
    return {
        "held_by_reason": dict(sorted(held.items())),
        "backoff": backoff,
        "next_attempt_due_at": _utc_text(earliest) if earliest is not None else None,
        "policy": {
            "base_seconds": active.base_seconds,
            "cap_seconds": active.cap_seconds,
            "max_attempts": active.max_attempts,
        },
    }


PROGRESS_SNAPSHOT_SCHEMA = "metainformant.rna.progress_snapshot.v1"
EXCLUDED_STATE = "excluded"
"""Effective state of a ``permanent_drop`` row in :func:`progress_snapshot` (whatever its recorded state)."""
EXIT_COMPLETE = 0
EXIT_INCOMPLETE = 3
EXIT_NO_DATABASE = 2


def _readonly_connection(db_path: Path) -> sqlite3.Connection:
    """Open an existing progress DB for reading only: no schema runs, nothing migrates, nothing checkpoints.

    A WAL-mode database with a ``-wal`` file (a live or crashed writer) is
    opened ``mode=ro`` so the un-checkpointed rows are read (``immutable=1``
    would ignore them). Without a ``-wal`` file every committed row is in the
    main file, and ``mode=ro`` would *create* empty ``-wal``/``-shm`` files
    next to it (a write on the data volume, and a failure on a read-only
    mount), so that case reads the main file with ``immutable=1``. A
    rollback-journal database is opened ``mode=ro``, which creates nothing.
    """

    path = Path(db_path).absolute()
    if not path.is_file():
        raise FileNotFoundError(f"progress database does not exist: {path}")
    with path.open("rb") as handle:
        header = handle.read(20)
    # SQLite header bytes 18/19 (file format write/read version) are 2 in WAL mode.
    wal_mode = len(header) == 20 and header[18] == 2 and header[19] == 2
    has_wal = path.with_name(path.name + "-wal").exists()
    mode = "immutable=1" if wal_mode and not has_wal else "mode=ro"
    connection = sqlite3.connect(f"{path.as_uri()}?{mode}", uri=True, timeout=5.0)
    connection.execute("PRAGMA busy_timeout=5000")
    return connection


def _table_columns(connection: sqlite3.Connection, table: str) -> set[str]:
    return {row[1] for row in connection.execute(f"PRAGMA table_info({table})")}


def progress_snapshot(
    db_path: Path,
    *,
    species: Optional[Iterable[str]] = None,
    policy: Optional[RetryPolicy] = None,
    now: Optional[float] = None,
) -> Dict[str, Any]:
    """Read-only per-species completion and failure summary of a progress database.

    Safe while a producer runs (``mode=ro``); tolerates databases written by
    older schemas (no ``sample_exclusions`` table, no retry columns). For
    every species (or only ``species``):

    * ``states``: recorded row states; ``effective_states`` moves every
      ``permanent_drop`` exclusion to ``excluded`` (it is not part of the
      analysable cohort, whatever its recorded state);
    * ``failure_classes``: :func:`classify_sample_error` over the failed rows
      that are not excluded;
    * ``held_by_reason`` / ``backoff``: the retry bookkeeping
      (:func:`retry_gate`) of the non-excluded, non-quantified rows;
    * ``nonterminal``: rows neither quantified nor excluded, and ``complete``
      (rows exist and none is nonterminal) -- the gate downstream merging uses.

    Raises :class:`FileNotFoundError` when the database does not exist.
    """

    active = policy or RetryPolicy.from_environment()
    moment = time.time() if now is None else now
    wanted = None if species is None else list(dict.fromkeys(species))
    path = Path(db_path).absolute()
    connection = _readonly_connection(path)
    try:
        columns = _table_columns(connection, "samples")
        retry_columns = all(name in columns for name, _ in RETRY_COLUMNS)
        select = "species, srr_id, state, error, updated_at" + (
            ", attempts, last_error_class, next_attempt_at" if retry_columns else ""
        )
        if wanted is None:
            rows = connection.execute(f"SELECT {select} FROM samples").fetchall()
        else:
            rows = []
            for name in wanted:
                rows.extend(connection.execute(f"SELECT {select} FROM samples WHERE species = ?", (name,)).fetchall())
        exclusions: Dict[tuple[str, str], str] = {}
        if _table_columns(connection, "sample_exclusions"):
            for name, srr_id, reason_code in connection.execute(
                "SELECT species, srr_id, reason_code FROM sample_exclusions"
            ):
                exclusions[(name, srr_id)] = reason_code
        binding: Dict[str, str] = {}
        if _table_columns(connection, "campaign_meta"):
            binding = {key: value for key, value in connection.execute("SELECT key, value FROM campaign_meta")}
    finally:
        connection.close()

    per_species: Dict[str, Dict[str, Any]] = {}
    for name in wanted or []:
        per_species[name] = {}
    latest: Optional[str] = None
    for row in rows:
        name, srr_id, state, error, updated_at = row[:5]
        entry = per_species.setdefault(name, {})
        entry.setdefault("rows", 0)
        entry["rows"] += 1
        states = entry.setdefault("states", {})
        states[state] = states.get(state, 0) + 1
        excluded = exclusions.get((name, srr_id)) == "permanent_drop"
        effective = EXCLUDED_STATE if excluded else state
        effective_states = entry.setdefault("effective_states", {})
        effective_states[effective] = effective_states.get(effective, 0) + 1
        if updated_at is not None and (latest is None or updated_at > latest):
            latest = updated_at
        if excluded or state == "quantified":
            continue
        if state == "failed":
            failure_class = classify_sample_error(error)
            classes = entry.setdefault("failure_classes", {})
            classes[failure_class] = classes.get(failure_class, 0) + 1
        if retry_columns:
            record = RetryRecord(name, srr_id, state, int(row[5] or 0), row[6], _utc_epoch(row[7]))
            gate, reason = retry_gate(record, active, moment)
            if gate == GATE_HELD:
                kind = (reason or "held").split(":", 1)[0]
                held = entry.setdefault("held_by_reason", {})
                held[kind] = held.get(kind, 0) + 1
            elif gate == GATE_BACKOFF:
                entry["backoff"] = entry.get("backoff", 0) + 1

    totals: Dict[str, int] = {}
    effective_totals: Dict[str, int] = {}
    failure_totals: Dict[str, int] = {}
    for name, entry in per_species.items():
        rows_count = int(entry.get("rows", 0))
        effective_states = dict(sorted(entry.get("effective_states", {}).items()))
        nonterminal = rows_count - effective_states.get("quantified", 0) - effective_states.get(EXCLUDED_STATE, 0)
        per_species[name] = {
            "rows": rows_count,
            "states": dict(sorted(entry.get("states", {}).items())),
            "effective_states": effective_states,
            "quantified": effective_states.get("quantified", 0),
            "excluded": effective_states.get(EXCLUDED_STATE, 0),
            "nonterminal": nonterminal,
            "failure_classes": dict(sorted(entry.get("failure_classes", {}).items())),
            "held_by_reason": dict(sorted(entry.get("held_by_reason", {}).items())),
            "backoff": int(entry.get("backoff", 0)),
            "complete": rows_count > 0 and nonterminal == 0,
        }
        for state, count in per_species[name]["states"].items():
            totals[state] = totals.get(state, 0) + count
        for state, count in effective_states.items():
            effective_totals[state] = effective_totals.get(state, 0) + count
        for failure_class, count in per_species[name]["failure_classes"].items():
            failure_totals[failure_class] = failure_totals.get(failure_class, 0) + count
    exclusions_by_reason: Dict[str, int] = {}
    for (name, _), reason_code in exclusions.items():
        if wanted is None or name in per_species:
            exclusions_by_reason[reason_code] = exclusions_by_reason.get(reason_code, 0) + 1
    return {
        "schema": PROGRESS_SNAPSHOT_SCHEMA,
        "db_path": str(path),
        "cohort_id": binding.get("cohort_id"),
        "rows": sum(entry["rows"] for entry in per_species.values()),
        "states": dict(sorted(totals.items())),
        "effective_states": dict(sorted(effective_totals.items())),
        "nonterminal": sum(entry["nonterminal"] for entry in per_species.values()),
        "failure_classes": dict(sorted(failure_totals.items())),
        "exclusions_by_reason": dict(sorted(exclusions_by_reason.items())),
        "latest_update": latest,
        "species": dict(sorted(per_species.items())),
    }


def completion_message(species: str, entry: Mapping[str, Any]) -> str:
    """One human line explaining a species' completion gate (see :func:`progress_snapshot`)."""

    if not entry.get("rows"):
        return f"{species}: no progress rows; refusing untracked downstream work"
    if entry.get("complete"):
        return (
            f"{species}: complete ({entry['quantified']} quantified, {entry['excluded']} permanently excluded "
            f"of {entry['rows']} rows)"
        )
    waiting = {
        state: count
        for state, count in entry.get("effective_states", {}).items()
        if state not in ("quantified", EXCLUDED_STATE)
    }
    detail = ", ".join(f"{state} {count}" for state, count in sorted(waiting.items()))
    held = entry.get("held_by_reason") or {}
    held_text = f"; held {sum(held.values())} ({', '.join(f'{k} {v}' for k, v in held.items())})" if held else ""
    return (
        f"{species}: incomplete, {entry['nonterminal']} of {entry['rows']} rows neither quantified nor "
        f"permanently excluded ({detail}){held_text}"
    )


def main(argv: Optional[Iterable[str]] = None) -> int:
    """``python -m metainformant.rna.engine.progress_db {retry-status,snapshot,completion,release-holds}``."""

    import argparse

    parser = argparse.ArgumentParser(
        prog="python -m metainformant.rna.engine.progress_db",
        description="Retry bookkeeping of a campaign progress database",
    )
    commands = parser.add_subparsers(dest="command", required=True)
    status = commands.add_parser("retry-status", help="Held and backed-off rows (read-only)")
    status.add_argument("--db", type=Path, required=True)
    snapshot = commands.add_parser("snapshot", help="Per-species states, exclusions and completion as JSON (read-only)")
    snapshot.add_argument("--db", type=Path, required=True)
    snapshot.add_argument("--species", action="append", help="Only this species (repeatable)")
    completion = commands.add_parser(
        "completion",
        help="Exit 0 when every row of the species is quantified or permanently excluded, else 3 (read-only)",
    )
    completion.add_argument("--db", type=Path, required=True)
    completion.add_argument("--species", action="append", required=True, help="Species to check (repeatable)")
    release = commands.add_parser(
        "release-holds", help="Return held rows to pending (under the campaign writer lock of the DB's data root)"
    )
    release.add_argument("--db", type=Path, required=True)
    release.add_argument("--species")
    release.add_argument("--accession", action="append", help="Only these runs (repeatable)")
    release.add_argument(
        "--reason",
        action="append",
        choices=(RETRY_LAYOUT_MISMATCH, RETRY_PERMANENT, "attempts_exhausted"),
        help="Only rows held for this reason (repeatable; default: all three)",
    )
    args = parser.parse_args(list(argv) if argv is not None else None)
    db_path = args.db.expanduser().absolute()
    if not db_path.is_file():
        print(f"error: progress database does not exist: {db_path}", file=sys.stderr)
        return EXIT_NO_DATABASE
    if args.command in ("snapshot", "completion"):
        try:
            summary = progress_snapshot(db_path, species=args.species)
        except sqlite3.Error as error:
            print(f"error: progress database is unreadable: {db_path} ({error})", file=sys.stderr)
            return EXIT_NO_DATABASE
        if args.command == "snapshot":
            print(json.dumps(summary, indent=2, sort_keys=True))
            return 0
        complete = True
        for name in args.species:
            entry = summary["species"][name]
            print(completion_message(name, entry))
            complete = complete and bool(entry["complete"])
        return EXIT_COMPLETE if complete else EXIT_INCOMPLETE
    if args.command == "retry-status":
        print(json.dumps(retry_status(db_path), indent=2, sort_keys=True))
        return 0
    from metainformant.rna.engine.cohort import DataRootBusy, data_root_maintenance_lock

    try:
        with data_root_maintenance_lock(db_path.parent, purpose="release-holds"), ProgressDB(db_path) as db:
            released = db.release_retry_holds(
                policy=RetryPolicy.from_environment(),
                species=args.species,
                accessions=args.accession,
                reasons=(
                    tuple(args.reason)
                    if args.reason
                    else (RETRY_LAYOUT_MISMATCH, RETRY_PERMANENT, "attempts_exhausted")
                ),
            )
    except (DataRootBusy, OSError) as error:
        print(f"error: {error}", file=sys.stderr)
        return 2
    print(json.dumps({"released": released}, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
