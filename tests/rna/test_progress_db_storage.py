"""Native on-disk progress-store migration and refusal controls."""
from __future__ import annotations

import sqlite3
import threading
from pathlib import Path

import pytest

from metainformant.rna.engine.progress_db import ProgressDB, configured_progress_db_path


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


def test_existing_override_and_default_paths(tmp_path: Path) -> None:
    default = tmp_path / "default.db"
    assert configured_progress_db_path(default, None) == default
    migrated = tmp_path / "migrated.db"
    with ProgressDB(migrated, journal_mode="DELETE"):
        pass
    assert configured_progress_db_path(default, str(migrated)) == migrated


@pytest.mark.parametrize("kind", ["relative", "empty", "missing", "directory", "symlink", "parent-symlink", "dotdot"])
def test_bad_override_paths_refuse_without_creating(tmp_path: Path, kind: str) -> None:
    good = tmp_path / "existing.db"
    good.write_text("Existing input is not opened by path resolution.")
    default = tmp_path / "unused.db"
    paths = {"relative": "relative.db", "empty": "", "missing": str(tmp_path / "missing.db"),
             "directory": str(tmp_path), "symlink": str(tmp_path / "link.db"),
             "parent-symlink": str(tmp_path / "link-dir/existing.db"),
             "dotdot": str(tmp_path / "unused/../existing.db")}
    (tmp_path / "link.db").symlink_to(good)
    (tmp_path / "link-dir").symlink_to(tmp_path, target_is_directory=True)
    with pytest.raises(ValueError, match="progress"):
        configured_progress_db_path(default, paths[kind])
    assert not default.exists()
    assert not (tmp_path / "missing.db").exists()


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
def test_native_worker_uses_migrated_store_and_replays(tmp_path: Path, monkeypatch: pytest.MonkeyPatch, native_acquisition_factory) -> None:
    from metainformant.rna.engine.acquisition_worker import run_manifest
    fixture = native_acquisition_factory()
    migrated = tmp_path / "migrated-progress.db"
    with ProgressDB(migrated, journal_mode="DELETE"):
        pass
    monkeypatch.setenv("AMALGKIT_PROGRESS_DB_PATH", str(migrated))
    monkeypatch.setenv("AMALGKIT_PROGRESS_DB_JOURNAL_MODE", "DELETE")
    args = dict(manifest_path=fixture.manifest, data_root=fixture.data, config_dir=fixture.config_dir,
                workers=1, threads=1, fastq_threads=1, compression_threads=1,
                validation_slots=1, quant_slots=1, fasterq_slots=1, max_in_flight=1)
    first = run_manifest(**args)
    assert first["counts"].get("failed", 0) == 0, first
    assert first["counts"]["newly_quantified"] == 1
    assert not (fixture.data / "pipeline_progress.db").exists()
    connection = sqlite3.connect(migrated.as_uri() + "?mode=ro", uri=True)
    try:
        assert connection.execute("PRAGMA journal_mode").fetchone() == ("delete",)
        assert connection.execute("SELECT state FROM samples WHERE srr_id='SRR123'").fetchone() == ("quantified",)
    finally:
        connection.close()
    abundance = next((fixture.work / "quant/SRR123").glob("*abundance.tsv"))
    preserved = abundance.read_bytes(), abundance.stat().st_mtime_ns
    replay = run_manifest(**args)
    assert replay["counts"]["reused"] == 1
    assert (abundance.read_bytes(), abundance.stat().st_mtime_ns) == preserved
    assert not Path(str(migrated) + "-shm").exists()
    assert not Path(str(migrated) + "-wal").exists()


