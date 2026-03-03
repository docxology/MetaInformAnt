"""SQLite-backed progress tracking for RNA-seq pipeline.

Replaces the JSON-based ProgressTracker with an SQLite database for:
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
from typing import Any, Dict, List, Optional, Set

from metainformant.core.utils.logging import get_logger

logger = get_logger(__name__)

# ---------- Constants ----------

VALID_STATES = frozenset({
    "pending",
    "downloading",
    "downloaded",
    "quantifying",
    "quantified",
    "failed",
})

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

    def set_state(
        self, species: str, srr_id: str, state: str, error: Optional[str] = None
    ) -> None:
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

    def bulk_set_state(
        self, species: str, srr_ids: List[str], state: str, error: Optional[str] = None
    ) -> None:
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

    # ---- Read operations ----

    def get_state(self, species: str, srr_id: str) -> Optional[str]:
        """Get the current state of a single sample, or None if not tracked."""
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
        if species:
            rows = self._conn.execute(
                "SELECT species, state, COUNT(*) FROM samples "
                "WHERE species = ? GROUP BY species, state",
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
        rows = self._conn.execute(
            "SELECT species, COUNT(*) FROM samples GROUP BY species"
        ).fetchall()
        return {sp: count for sp, count in rows}

    def get_samples(self, species: str, state: str) -> List[str]:
        """Get all SRR IDs in a given state for a species."""
        rows = self._conn.execute(
            "SELECT srr_id FROM samples WHERE species = ? AND state = ?",
            (species, state),
        ).fetchall()
        return [r[0] for r in rows]

    def get_failed(self, species: Optional[str] = None) -> List[Dict[str, str]]:
        """Get all failed samples with their error messages."""
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
        return [
            {"species": r[0], "srr_id": r[1], "error": r[2], "updated_at": r[3]}
            for r in rows
        ]

    # ---- Reconciliation ----

    def reconcile(self, species: str, quant_dir: Path) -> int:
        """Detect already-quantified samples from filesystem and update DB.

        Scans ``quant_dir`` for subdirectories containing ``*_abundance.tsv``
        and marks those samples as 'quantified' in the database.

        Args:
            species: Species name.
            quant_dir: Path to the quant output directory
                       (e.g. ``output/amalgkit/amellifera/work/quant``).

        Returns:
            Number of samples reconciled.
        """
        if not quant_dir.exists():
            return 0

        reconciled = []
        try:
            for subdir in quant_dir.iterdir():
                if not subdir.is_dir():
                    continue
                if any(subdir.glob("*_abundance.tsv")) or (subdir / "quant.sf").exists():
                    reconciled.append(subdir.name)
        except Exception as e:
            logger.warning(f"Reconciliation scan failed for {species}: {e}")
            return 0

        if reconciled:
            self.bulk_set_state(species, reconciled, "quantified")
            logger.info(f"Reconciled {len(reconciled)} quantified samples for {species}")

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
