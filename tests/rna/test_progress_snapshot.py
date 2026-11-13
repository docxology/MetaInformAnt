"""Read-only progress snapshot and the downstream completion gate.

Real SQLite databases written through :class:`ProgressDB` (and one minimal
legacy schema), no mocks: the snapshot must count permanent exclusions as
terminal, report retry holds, see a live writer's WAL, and never write.
"""

from __future__ import annotations

import hashlib
import json
import sqlite3
import subprocess
import sys
from pathlib import Path

import pytest

from metainformant.rna.engine.progress_db import (
    EXIT_COMPLETE,
    EXIT_INCOMPLETE,
    EXIT_NO_DATABASE,
    ProgressDB,
    RetryPolicy,
    completion_message,
    progress_snapshot,
)

KALLISTO_EXIT_1 = "Quant exception batch 3: kallisto quant failed with exit code 1"


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _campaign_db(path: Path) -> Path:
    """Two species: apis has quantified, permanently dropped, pending and a held kallisto exit-1 row."""

    with ProgressDB(path) as db:
        db.init_species("apis_mellifera", ["SRR1", "SRR2", "SRR3", "SRR4", "SRR5"])
        db.init_species("bombus_terrestris", ["ERR1", "ERR2"])
        db.set_state("apis_mellifera", "SRR1", "quantified")
        db.set_state("apis_mellifera", "SRR2", "failed", KALLISTO_EXIT_1)
        db.record_exclusions("apis_mellifera", [{"srr_id": "SRR2", "reason_code": "permanent_drop"}])
        db.set_state(
            "apis_mellifera", "SRR3", "failed", "Unexpected sample task error: [Errno 1] Operation not permitted"
        )
        db.set_state("apis_mellifera", "SRR4", "failed", KALLISTO_EXIT_1)
        db.record_attempt_outcome("apis_mellifera", "SRR4", KALLISTO_EXIT_1, policy=RetryPolicy())
        db.record_exclusions("apis_mellifera", [{"srr_id": "SRR5", "reason_code": "re_download"}])
        db.bulk_set_state("bombus_terrestris", ["ERR1", "ERR2"], "quantified")
    return path


def test_snapshot_counts_permanent_drops_as_terminal_and_reports_holds(tmp_path: Path) -> None:
    db_path = _campaign_db(tmp_path / "pipeline_progress.db")

    snapshot = progress_snapshot(db_path)

    apis = snapshot["species"]["apis_mellifera"]
    assert apis["rows"] == 5
    assert apis["states"] == {"failed": 3, "pending": 1, "quantified": 1}
    assert apis["effective_states"] == {"excluded": 1, "failed": 2, "pending": 1, "quantified": 1}
    assert (apis["quantified"], apis["excluded"], apis["nonterminal"], apis["complete"]) == (1, 1, 3, False)
    # The excluded kallisto failure is not counted; the held one and the EPERM row are.
    assert apis["failure_classes"] == {"environment_write_denied": 1, "quantification_kallisto_exit_1": 1}
    assert apis["held_by_reason"] == {"permanent": 1}
    bombus = snapshot["species"]["bombus_terrestris"]
    assert bombus["complete"] is True and bombus["nonterminal"] == 0
    assert snapshot["rows"] == 7
    assert snapshot["nonterminal"] == 3
    assert snapshot["exclusions_by_reason"] == {"permanent_drop": 1, "re_download": 1}
    assert snapshot["latest_update"] is not None
    assert "complete (2 quantified, 0 permanently excluded of 2 rows)" in completion_message(
        "bombus_terrestris", bombus
    )
    message = completion_message("apis_mellifera", apis)
    assert "3 of 5 rows neither quantified nor permanently excluded" in message
    assert "held 1 (permanent 1)" in message


def test_snapshot_is_read_only(tmp_path: Path) -> None:
    db_path = _campaign_db(tmp_path / "pipeline_progress.db")
    before = {path.name: _sha(path) for path in tmp_path.iterdir()}

    progress_snapshot(db_path, species=["apis_mellifera"])

    assert {path.name: _sha(path) for path in tmp_path.iterdir()} == before
    with pytest.raises(FileNotFoundError):
        progress_snapshot(tmp_path / "missing.db")
    assert not (tmp_path / "missing.db").exists()


def test_snapshot_reads_a_live_writers_wal(tmp_path: Path) -> None:
    db_path = tmp_path / "pipeline_progress.db"
    writer = ProgressDB(db_path)
    try:
        writer.init_species("apis_mellifera", ["SRR1"])
        writer.set_state("apis_mellifera", "SRR1", "quantified")
        assert (tmp_path / "pipeline_progress.db-wal").stat().st_size > 0
        snapshot = progress_snapshot(db_path)
        assert snapshot["species"]["apis_mellifera"]["complete"] is True
    finally:
        writer.close()


def test_snapshot_tolerates_a_legacy_schema_and_unknown_species(tmp_path: Path) -> None:
    db_path = tmp_path / "legacy.db"
    connection = sqlite3.connect(db_path)
    connection.execute(
        "CREATE TABLE samples (species TEXT NOT NULL, srr_id TEXT NOT NULL, state TEXT NOT NULL, error TEXT, "
        "updated_at TEXT, PRIMARY KEY (species, srr_id))"
    )
    connection.execute("INSERT INTO samples VALUES ('nasonia_vitripennis', 'SRR9', 'failed', NULL, '2026-10-13')")
    connection.commit()
    connection.close()

    snapshot = progress_snapshot(db_path, species=["nasonia_vitripennis", "never_seen"])

    assert snapshot["species"]["nasonia_vitripennis"]["failure_classes"] == {"unrecorded": 1}
    assert snapshot["species"]["nasonia_vitripennis"]["held_by_reason"] == {}
    assert snapshot["species"]["never_seen"]["rows"] == 0
    assert snapshot["species"]["never_seen"]["complete"] is False
    assert "no progress rows" in completion_message("never_seen", snapshot["species"]["never_seen"])


def _cli(*args: str) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [sys.executable, "-m", "metainformant.rna.engine.progress_db", *args],
        capture_output=True,
        text=True,
        timeout=120,
        check=False,
    )


def test_completion_and_snapshot_cli_exit_codes(tmp_path: Path) -> None:
    db_path = _campaign_db(tmp_path / "pipeline_progress.db")

    complete = _cli("completion", "--db", str(db_path), "--species", "bombus_terrestris")
    assert complete.returncode == EXIT_COMPLETE, complete.stderr
    assert "bombus_terrestris: complete" in complete.stdout

    incomplete = _cli(
        "completion", "--db", str(db_path), "--species", "bombus_terrestris", "--species", "apis_mellifera"
    )
    assert incomplete.returncode == EXIT_INCOMPLETE
    assert "apis_mellifera: incomplete" in incomplete.stdout

    missing = _cli("completion", "--db", str(tmp_path / "absent.db"), "--species", "apis_mellifera")
    assert missing.returncode == EXIT_NO_DATABASE
    assert not (tmp_path / "absent.db").exists()

    snapshot = _cli("snapshot", "--db", str(db_path))
    assert snapshot.returncode == 0, snapshot.stderr
    payload = json.loads(snapshot.stdout)
    assert payload["schema"] == "metainformant.rna.progress_snapshot.v1"
    assert payload["species"]["apis_mellifera"]["excluded"] == 1


def test_snapshot_of_a_rollback_journal_database_creates_nothing(tmp_path: Path) -> None:
    db_path = tmp_path / "pipeline_progress.db"
    with ProgressDB(db_path, journal_mode="DELETE") as db:
        db.init_species("apis_mellifera", ["SRR1"])
    before = {path.name: _sha(path) for path in tmp_path.iterdir()}

    assert progress_snapshot(db_path)["species"]["apis_mellifera"]["effective_states"] == {"pending": 1}
    assert {path.name: _sha(path) for path in tmp_path.iterdir()} == before


def test_kallisto_exit_1_marker_matches_the_pinned_amalgkit_message() -> None:
    """The permanent kallisto exit-1 class must keep matching what the pinned Amalgkit actually prints."""

    import inspect

    from metainformant.rna.engine.progress_db import SAMPLE_ERROR_CLASSES, classify_sample_error

    quant = pytest.importorskip("amalgkit.quant", reason="the rna extra (pinned Amalgkit) is not installed")
    template = "kallisto quant failed with exit code {} for {}."
    assert template in inspect.getsource(quant)  # the message format of the installed (pinned) build
    message = "ERROR: " + template.format(1, "SRR21601882")
    assert SAMPLE_ERROR_CLASSES["quantification_kallisto_exit_1"] in message
    assert classify_sample_error(f"Quant exception batch 3: {message}") == "quantification_kallisto_exit_1"
    assert classify_sample_error("ERROR: " + template.format(2, "SRR21601882")) != "quantification_kallisto_exit_1"
