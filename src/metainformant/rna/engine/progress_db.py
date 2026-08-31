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

import sqlite3
import threading
import time
from pathlib import Path
from typing import Dict, List, Optional

from metainformant.core.utils.logging import get_logger
from metainformant.rna.core.sample_utils import find_quantification_file
from metainformant.rna.engine.provenance import QUANT_STATUS_CURRENT, classify_quantification

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

    def __init__(self, db_path: Optional[Path] = None):
        self.db_path = Path(db_path) if db_path else DEFAULT_DB_PATH
        self.db_path.parent.mkdir(parents=True, exist_ok=True)
        self._lock = threading.Lock()
        self._conn = sqlite3.connect(
            str(self.db_path),
            check_same_thread=False,
            isolation_level="DEFERRED",
        )
        self._conn.execute("PRAGMA journal_mode=WAL")
        self._conn.execute("PRAGMA busy_timeout=5000")
        self._conn.executescript(_SCHEMA)
        self._conn.commit()
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
        """Remove progress rows outside the current metadata cohort.

        A refreshed Amalgkit metadata query can legitimately change the set of
        public runs.  Keeping rows from an older universe would make a fully
        quantified current cohort appear permanently pending and could block
        downstream finalization.  This only removes internal SQLite state; raw
        reads and quantification outputs are never touched.
        """

        keep = {str(srr_id) for srr_id in srr_ids}
        with self._lock:
            existing = [
                row[0]
                for row in self._conn.execute(
                    "SELECT srr_id FROM samples WHERE species = ?",
                    (species,),
                ).fetchall()
            ]
            obsolete = [srr_id for srr_id in existing if srr_id not in keep]
            if obsolete:
                self._conn.executemany(
                    "DELETE FROM samples WHERE species = ? AND srr_id = ?",
                    [(species, srr_id) for srr_id in obsolete],
                )
                self._conn.commit()
        if obsolete:
            logger.info("prune_species(%s): removed %d stale progress rows", species, len(obsolete))
        return len(obsolete)

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

    def get_excluded_srr_ids(self, species: str, reason_code: Optional[str] = None) -> set:
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
        """
        if not quant_dir.exists():
            return 0

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

    # ---- Lifecycle ----

    def close(self) -> None:
        """Close the database connection."""
        self._conn.close()

    def __enter__(self):
        return self

    def __exit__(self, *args):
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
