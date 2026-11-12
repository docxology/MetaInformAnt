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

import concurrent.futures
import os
import sqlite3
import threading
import time
from pathlib import Path
from typing import TYPE_CHECKING, Any, Callable, Dict, Iterable, List, Mapping, Optional

if TYPE_CHECKING:
    from metainformant.rna.engine.cohort import FrozenCohort

from metainformant.core.utils.logging import get_logger
from metainformant.rna.core.sample_utils import find_quantification_file
from metainformant.rna.engine.provenance import (
    QUANT_STATUS_CURRENT,
    REUSABLE_QUANT_STATUSES,
    classify_quantification,
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
    "environment_missing_tool": "No such file or directory: '",
    "transfer_all_sources_failed": "Download Failed (all sources",
    "extraction_timeout": "fasterq-dump timeout",
    "quantification_timeout": "Quant timeout",
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


RETRYABLE_ENVIRONMENT_CLASSES = frozenset(
    {
        "environment_write_denied",
        "environment_access_denied",
        "environment_read_only",
        "environment_no_space",
        "environment_missing_tool",
    }
)
"""Failure classes caused by the producer's host (EPERM, EACCES, EROFS, ENOSPC, ENOENT), not by the sample."""

TRANSFER_RESET_CLASSES = frozenset({"transfer_all_sources_failed"})
"""Transfer failures: possibly per-sample, reset only on explicit request."""

RETRYABLE_RESET_CLASSES = RETRYABLE_ENVIRONMENT_CLASSES | TRANSFER_RESET_CLASSES

COHORT_MIGRATION_ERROR_PREFIX = "reset-by-cohort-migration"


def classify_sample_error(error: str | None) -> str:
    """Classify a stored sample failure message into a coarse retry class."""

    if not error:
        return "unrecorded"
    for class_name, marker in SAMPLE_ERROR_CLASSES.items():
        if marker in error:
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
        """Reconcile only provenance-qualified current quantification outputs.

        Readable abundance files without a current provenance sidecar remain
        visible on disk but are not promoted to ``quantified`` in the DB.

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
                if classification.get("status") == QUANT_STATUS_CURRENT:
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
                "Left %d readable but non-current quantifications visible for %s: %s",
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

    def reconcile_cohort(
        self,
        species: str,
        quant_dir: Path,
        accessions: Iterable[str],
        classify: Callable[[Path, str], Mapping[str, Any]],
        *,
        requantify_statuses: Iterable[str] = (),
        workers: int = 1,
    ) -> Dict[str, Any]:
        """Bring task states in line with the outputs on disk, in one transaction.

        Only ``accessions`` (the species' schedulable cohort tasks) are
        considered; nothing else is read or changed. For each task:

        * a reusable output (``classify`` status current or version drift)
          becomes ``quantified`` with its audit, unless its status is in
          ``requantify_statuses`` (then ``pending``);
        * a recorded ``quantified`` state whose output is not reusable becomes
          ``quarantined`` with the reason (the bytes stay on disk);
        * a recorded ``quantified`` state without an output directory becomes
          ``pending``.

        Rows in other states with a non-reusable output are left alone (the
        producer re-quantifies them). Returns counts and the quarantined runs.
        """

        wanted = {str(accession) for accession in accessions}
        requeue = frozenset(requantify_statuses)
        directories: Dict[str, Path] = {}
        if quant_dir.is_dir():
            with os.scandir(quant_dir) as entries:
                for entry in entries:
                    if entry.name in wanted and entry.is_dir(follow_symlinks=False):
                        directories[entry.name] = Path(entry.path)
        with self._lock:
            states = {
                srr_id: state
                for srr_id, state in self._conn.execute(
                    "SELECT srr_id, state FROM samples WHERE species = ?", (species,)
                )
                if srr_id in wanted
            }

        def judge(item: tuple[str, Path]) -> tuple[str, Mapping[str, Any]]:
            accession, path = item
            try:
                return accession, classify(path, accession)
            except (OSError, ValueError) as error:
                return accession, {"status": "invalid", "reason": f"classification failed: {error}"}

        pairs = sorted(directories.items())
        if workers > 1 and len(pairs) > 1:
            with concurrent.futures.ThreadPoolExecutor(max_workers=min(workers, len(pairs))) as pool:
                judged = list(pool.map(judge, pairs))
        else:
            judged = [judge(pair) for pair in pairs]

        quantified: List[tuple[str, Mapping[str, Any]]] = []
        requeued: List[tuple[str, Mapping[str, Any]]] = []
        quarantined: List[tuple[str, Mapping[str, Any]]] = []
        for accession, audit in judged:
            status = audit.get("status")
            if status in REUSABLE_QUANT_STATUSES:
                (requeued if status in requeue else quantified).append((accession, audit))
            elif states.get(accession) == "quantified":
                quarantined.append((accession, audit))
        vanished = [
            accession for accession, state in states.items() if state == "quantified" and accession not in directories
        ]
        with self._lock:
            self._conn.commit()
            self._conn.execute("BEGIN IMMEDIATE")
            try:
                for accession, audit in quantified + requeued + quarantined:
                    self._conn.execute(
                        "INSERT OR IGNORE INTO samples (species, srr_id, state) VALUES (?, ?, 'pending')",
                        (species, accession),
                    )
                    self._upsert_audit_locked(species, accession, audit)
                for accession, _ in quantified:
                    self._conn.execute(
                        "UPDATE samples SET state = 'quantified', error = NULL, updated_at = datetime('now') "
                        "WHERE species = ? AND srr_id = ? AND state <> 'quantified'",
                        (species, accession),
                    )
                for accession, audit in requeued:
                    self._conn.execute(
                        "UPDATE samples SET state = 'pending', error = ?, updated_at = datetime('now') "
                        "WHERE species = ? AND srr_id = ?",
                        (f"requantify: {audit.get('status')}", species, accession),
                    )
                for accession, audit in quarantined:
                    self._conn.execute(
                        "UPDATE samples SET state = 'quarantined', error = ?, updated_at = datetime('now') "
                        "WHERE species = ? AND srr_id = ?",
                        (str(audit.get("reason")), species, accession),
                    )
                for accession in vanished:
                    self._conn.execute(
                        "UPDATE samples SET state = 'pending', error = ?, updated_at = datetime('now') "
                        "WHERE species = ? AND srr_id = ? AND state = 'quantified'",
                        ("reconcile: recorded quantification has no output directory", species, accession),
                    )
                self._conn.commit()
            except BaseException:
                self._conn.rollback()
                raise
        summary: Dict[str, Any] = {
            "scanned": len(judged),
            "quantified": len(quantified),
            "requeued": len(requeued),
            "quarantined": sorted(accession for accession, _ in quarantined),
            "vanished": sorted(vanished),
        }
        logger.info(
            "reconcile_cohort(%s): %d outputs scanned, %d reusable, %d requeued, %d quarantined, %d vanished",
            species,
            summary["scanned"],
            summary["quantified"],
            summary["requeued"],
            len(summary["quarantined"]),
            len(summary["vanished"]),
        )
        return summary

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
