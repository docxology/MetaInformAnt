"""Native on-disk progress-store migration and refusal controls."""

from __future__ import annotations

import sqlite3
import threading
from pathlib import Path

import pytest

from metainformant.rna.engine.progress_db import ProgressDB


@pytest.mark.parametrize("mode", ["WAL", "DELETE"])
def test_journal_modes_preserve_states_and_concurrent_writes(tmp_path: Path, mode: str) -> None:
    path = tmp_path / "progress.db"
    db = ProgressDB(path, journal_mode=mode)
    db.init_species("ant", [f"SRR{i}" for i in range(24)])
    errors = []

    def write(i: int) -> None:
        try:
            db.set_state("ant", f"SRR{i}", "failed" if i % 2 else "quantified", error=f"retained-{i}")
        except Exception as exc:
            errors.append(exc)

    workers = [threading.Thread(target=write, args=(i,)) for i in range(24)]
    for worker in workers:
        worker.start()
    for worker in workers:
        worker.join()
    assert not errors
    assert db._conn.execute("PRAGMA journal_mode").fetchone() == (mode.lower(),)
    if mode == "DELETE":
        assert db._conn.execute("PRAGMA synchronous").fetchone() == (2,)
        assert db._conn.execute("PRAGMA mmap_size").fetchone() == (0,)
        assert not Path(str(path) + "-shm").exists()
        assert not Path(str(path) + "-wal").exists()
    db.close()
    reopened = ProgressDB(path, journal_mode=mode)
    try:
        assert reopened.init_species("ant", [f"SRR{i}" for i in range(24)]) == 0
        assert reopened.get_counts()["ant"]["failed"] == 12
        assert reopened.get_counts()["ant"]["quantified"] == 12
        assert reopened._conn.execute("SELECT error FROM samples WHERE srr_id='SRR1'").fetchone() == ("retained-1",)
    finally:
        reopened.close()


@pytest.mark.parametrize("mode", ["OFF", "MEMORY", "TRUNCATE", "", "delete", "WAL; DROP TABLE samples"])
def test_unsupported_journal_refuses_before_creation(tmp_path: Path, mode: str) -> None:
    path = tmp_path / "must-not-create" / "progress.db"
    with pytest.raises(ValueError, match="journal"):
        ProgressDB(path, journal_mode=mode)
    assert not path.parent.exists()


def test_default_wal_is_retained(tmp_path: Path) -> None:
    with ProgressDB(tmp_path / "default.db") as db:
        assert db._conn.execute("PRAGMA journal_mode").fetchone() == ("wal",)


def test_native_sqlite_backup_keeps_all_schema_and_rows(tmp_path: Path) -> None:
    old, migrated = tmp_path / "old.db", tmp_path / "migrated.db"
    with ProgressDB(old) as db:
        db.init_species("ant", ["SRR1", "SRR2"])
        db.set_state("ant", "SRR1", "failed", error="retry history stays")
        before = db._conn.execute("SELECT * FROM samples ORDER BY srr_id").fetchall()
        out = sqlite3.connect(migrated)
        db._conn.backup(out)
        out.close()
    with ProgressDB(migrated, journal_mode="DELETE") as db:
        assert db._conn.execute("SELECT * FROM samples ORDER BY srr_id").fetchall() == before
        assert db._conn.execute("PRAGMA integrity_check").fetchone() == ("ok",)


@pytest.mark.external_tool
def test_native_producer_uses_migrated_store_and_replays(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, native_acquisition_factory
) -> None:
    fixture = native_acquisition_factory()
    migrated = tmp_path / "migrated-progress.db"
    with ProgressDB(migrated, journal_mode="DELETE"):
        pass
    monkeypatch.setenv("AMALGKIT_PROGRESS_DB_JOURNAL_MODE", "DELETE")
    first = fixture.produce(db_path=migrated)
    assert first["quantified"] and not first["skipped"] and first["error"] is None, first
    assert not (fixture.data / "pipeline_progress.db").exists()
    connection = sqlite3.connect(migrated.as_uri() + "?mode=ro", uri=True)
    try:
        assert connection.execute("PRAGMA journal_mode").fetchone() == ("delete",)
        assert connection.execute("SELECT state FROM samples WHERE srr_id='SRR123'").fetchone() == ("quantified",)
    finally:
        connection.close()
    abundance = next((fixture.work / "quant/SRR123").glob("*abundance.tsv"))
    preserved = abundance.read_bytes(), abundance.stat().st_mtime_ns
    replay = fixture.produce(db_path=migrated)
    assert replay["quantified"] and replay["skipped"], replay
    assert (abundance.read_bytes(), abundance.stat().st_mtime_ns) == preserved
    assert not Path(str(migrated) + "-shm").exists()
    assert not Path(str(migrated) + "-wal").exists()
