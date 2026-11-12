"""Progress DB startup resets and the change-proportional resume reconcile.

Real SQLite databases, real quantification outputs written by the production
``write_quant_provenance`` (via ``DataRoot``) and the real reuse classifier.
System calls are counted by wrappers that delegate to the real functions.
"""

from __future__ import annotations

import os
import sqlite3
from collections.abc import Callable, Mapping
from pathlib import Path
from typing import Any

import pytest

from metainformant.rna.engine.progress_db import (
    INFLIGHT_STATES,
    QUANT_PROVENANCE_FILENAME,
    ProgressDB,
    ReconcileScope,
    classify_sample_error,
    sidecar_fingerprint,
)
from metainformant.rna.engine.provenance import build_species_reuse_context
from tests.rna.frozen_cohort_support import SPECIES, DataRoot

EPERM = "Unexpected sample task error: [Errno 1] Operation not permitted: '/data/x'"
ENXIO = "Unexpected sample task error: [Errno 6] Device not configured: '/data/test_species/work'"
ENOSPC = "fasterq-dump failed: [Errno 28] No space left on device"
KALLISTO = "Quantification Failed: kallisto exited 1"
SENTINEL = "2000-01-01 00:00:00"


def _row(db_path: Path, srr_id: str) -> tuple[str, str | None, str]:
    with sqlite3.connect(db_path) as connection:
        row = connection.execute(
            "SELECT state, error, updated_at FROM samples WHERE species = ? AND srr_id = ?", (SPECIES, srr_id)
        ).fetchone()
    assert row is not None
    return row


def _set_updated_at(db_path: Path, value: str) -> None:
    with sqlite3.connect(db_path) as connection:
        connection.execute("UPDATE samples SET updated_at = ?", (value,))
        connection.execute("UPDATE quantification_audit SET checked_at = ?", (value,))


def test_reset_inflight_states_has_no_age_threshold(tmp_path: Path) -> None:
    db = ProgressDB(tmp_path / "p.db")
    runs = ["SRR1", "SRR2", "SRR3", "SRR4", "SRR5", "SRR6"]
    db.init_species(SPECIES, runs)
    db.set_state(SPECIES, "SRR1", "downloading", error="ENA attempt 1")
    db.set_state(SPECIES, "SRR2", "downloaded")
    db.set_state(SPECIES, "SRR3", "quantifying", error="batch 3")
    db.set_state(SPECIES, "SRR4", "quantified")
    db.set_state(SPECIES, "SRR5", "failed", error=KALLISTO)
    # Rows updated a second ago are reset too: under the writer lock nothing else is in flight.
    assert set(INFLIGHT_STATES) == {"downloading", "downloaded", "quantifying"}

    counts = db.reset_inflight_states()

    assert counts == {"downloaded": 1, "downloading": 1, "quantifying": 1}
    assert db.get_state(SPECIES, "SRR1") == "pending"
    assert _row(db.db_path, "SRR1")[1] == "reset-inflight-at-startup[downloading]: ENA attempt 1"
    assert _row(db.db_path, "SRR2")[1] == "reset-inflight-at-startup[downloaded]"
    assert _row(db.db_path, "SRR3")[:2] == ("pending", "reset-inflight-at-startup[quantifying]: batch 3")
    assert db.get_state(SPECIES, "SRR4") == "quantified"
    assert _row(db.db_path, "SRR5")[:2] == ("failed", KALLISTO)
    assert db.get_state(SPECIES, "SRR6") == "pending"
    assert db.reset_inflight_states() == {}
    db.close()


def test_reset_environment_failures_keeps_sample_failures_and_permanent_drops(tmp_path: Path) -> None:
    db = ProgressDB(tmp_path / "p.db")
    db.init_species(SPECIES, ["SRR1", "SRR2", "SRR3", "SRR4", "SRR5"])
    db.set_state(SPECIES, "SRR1", "failed", error=EPERM)
    db.set_state(SPECIES, "SRR2", "failed", error=ENXIO)
    db.set_state(SPECIES, "SRR3", "failed", error=ENOSPC)
    db.set_state(SPECIES, "SRR4", "failed", error=KALLISTO)
    db.set_state(SPECIES, "SRR5", "failed", error=EPERM)
    db.record_exclusions(SPECIES, [{"srr_id": "SRR5", "reason_code": "permanent_drop"}])
    assert classify_sample_error(ENXIO) == "environment_device_unavailable"

    counts = db.reset_environment_failures()

    assert counts == {
        "environment_device_unavailable": 1,
        "environment_no_space": 1,
        "environment_write_denied": 1,
    }
    assert _row(db.db_path, "SRR1")[:2] == (
        "pending",
        f"reset-environment-failure-at-startup[environment_write_denied]: {EPERM}",
    )
    assert db.get_state(SPECIES, "SRR2") == "pending"
    assert db.get_state(SPECIES, "SRR3") == "pending"
    assert _row(db.db_path, "SRR4")[:2] == ("failed", KALLISTO)  # a sample failure stays terminal
    assert _row(db.db_path, "SRR5")[:2] == ("failed", EPERM)  # the operator excluded it
    assert db.reset_environment_failures() == {}
    with pytest.raises(ValueError, match="unknown error classes"):
        db.reset_environment_failures(["not_a_class"])
    db.close()


class Outputs:
    """A data root with real quantification outputs and the real resume classifier, counting judgements."""

    def __init__(self, base: Path, runs: list[str]) -> None:
        self.data = DataRoot(base / "data")
        self.quant = self.data.work / "quant"
        self.runs = runs
        for run in runs:
            self.data.write_sample(run)
        context = build_species_reuse_context(SPECIES, self.data.config)
        self.scope = ReconcileScope(context.config_sha256, "reference-sha256", "reconcile.v1|test")
        self.judged: list[str] = []
        self._classify = context.classify

    def classify(self, path: Path, run: str) -> Mapping[str, Any]:
        self.judged.append(run)
        return self._classify(path, run)

    def reconcile(self, db: ProgressDB, **kwargs: Any) -> dict[str, Any]:
        self.judged.clear()
        kwargs.setdefault("scope", self.scope)
        return db.reconcile_outputs(SPECIES, self.quant, self.runs, self.classify, **kwargs)


def test_restart_trusts_unchanged_outputs_and_rejudges_only_changes(tmp_path: Path) -> None:
    runs = [f"SRR{i:03d}" for i in range(1, 9)]
    outputs = Outputs(tmp_path, runs)
    db = ProgressDB(tmp_path / "p.db")
    db.init_species(SPECIES, runs + ["SRR900"])
    db.set_state(SPECIES, "SRR900", "quantified")  # recorded, but its output directory is gone
    outputs.runs = runs + ["SRR900"]

    first = outputs.reconcile(db, workers=3)
    assert first["classified"] == 8 and first["trusted"] == 0 and first["quantified"] == 8
    assert first["vanished"] == ["SRR900"] and db.get_state(SPECIES, "SRR900") == "pending"
    assert all(db.get_state(SPECIES, run) == "quantified" for run in runs)

    _set_updated_at(db.db_path, SENTINEL)
    second = outputs.reconcile(db, workers=3)
    assert outputs.judged == []
    assert second["trusted"] == 8 and second["classified"] == 0
    assert second["rows_changed"] == 0 and second["audits_written"] == 0
    # Unchanged rows keep updated_at and checked_at: a restart rewrites nothing.
    assert all(_row(db.db_path, run)[2] == SENTINEL for run in runs + ["SRR900"])
    with sqlite3.connect(db.db_path) as connection:
        assert {row[0] for row in connection.execute("SELECT checked_at FROM quantification_audit")} == {SENTINEL}

    # One output rewritten with other configuration bytes: only it is re-judged, and it is quarantined.
    other_config = tmp_path / "other.yaml"
    other_config.write_bytes(b"species_list: [Other]\n")
    changed = outputs.quant / "SRR003"
    for child in changed.iterdir():
        child.unlink()
    changed.rmdir()
    outputs.data.write_sample("SRR003", config=other_config)
    third = outputs.reconcile(db)
    assert outputs.judged == ["SRR003"]
    assert third["newly_quarantined"] == ["SRR003"] and third["rows_changed"] == 1
    assert db.get_state(SPECIES, "SRR003") == "quarantined"
    assert all(_row(db.db_path, run)[2] == SENTINEL for run in runs if run != "SRR003")

    # Another scope (configuration, reference, runtime or rule revision changed) re-judges everything once.
    other_scope = ReconcileScope(outputs.scope.config_sha256, "another-reference", outputs.scope.scope)
    fourth = outputs.reconcile(db, scope=other_scope)
    assert fourth["classified"] == 8
    assert outputs.reconcile(db, scope=other_scope)["classified"] == 0
    # The deep audit (AMALGKIT_RESUME_VERIFY_HASHES=1) re-judges every output regardless of fingerprints.
    assert outputs.reconcile(db, scope=other_scope, deep=True)["classified"] == 8
    # Without a scope nothing is trusted.
    assert outputs.reconcile(db, scope=None)["classified"] == 8
    db.close()


def test_requantification_policy_requeues_through_the_same_pass(tmp_path: Path) -> None:
    runs = ["SRR001", "SRR002"]
    outputs = Outputs(tmp_path, runs)
    db = ProgressDB(tmp_path / "p.db")
    db.init_species(SPECIES, runs)
    summary = outputs.reconcile(db, requantify_statuses={"current"})
    assert summary["requeued"] == 2
    assert _row(db.db_path, "SRR001")[:2] == ("pending", "requantify: current")
    _set_updated_at(db.db_path, SENTINEL)
    again = outputs.reconcile(db, requantify_statuses={"current"})
    assert again["classified"] == 0 and again["rows_changed"] == 0
    assert _row(db.db_path, "SRR001")[2] == SENTINEL
    db.close()


def test_fingerprint_columns_are_added_once_and_unfingerprinted_audits_are_never_trusted(tmp_path: Path) -> None:
    outputs = Outputs(tmp_path, ["SRR001"])
    db = ProgressDB(tmp_path / "p.db")
    db.init_species(SPECIES, ["SRR001"])
    outputs.reconcile(db)
    with sqlite3.connect(db.db_path) as connection:
        columns = {row[1] for row in connection.execute("PRAGMA table_info(quantification_audit)")}
    assert {"sidecar_ino", "sidecar_size", "sidecar_mtime_ns", "config_sha256", "reference_sha256", "scope"} <= (
        columns
    )
    fingerprint = sidecar_fingerprint(outputs.quant / "SRR001" / QUANT_PROVENANCE_FILENAME)
    assert fingerprint is not None and fingerprint[1] > 0
    assert sidecar_fingerprint(outputs.quant / "SRR001") is None  # a directory is not a sidecar

    # Another writer (adoption, the producer's own audit) rewrites the audit without a fingerprint.
    db.record_quantification_audit(SPECIES, "SRR001", status="current", reason="re-audited")
    assert outputs.reconcile(db)["classified"] == 1
    assert outputs.reconcile(db)["classified"] == 0
    db.close()


def test_restart_stat_cost_is_one_listing_and_one_lstat_per_output(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    runs = [f"SRR{i:04d}" for i in range(300)]
    outputs = Outputs(tmp_path, runs)
    db = ProgressDB(tmp_path / "p.db")
    db.init_species(SPECIES, runs)
    assert outputs.reconcile(db, workers=4)["classified"] == 300

    calls = {"scandir": 0, "lstat": 0, "stat": 0}

    def counted(name: str, real: Callable[..., Any]) -> Callable[..., Any]:
        def wrapper(*args: Any, **kwargs: Any) -> Any:
            calls[name] += 1
            return real(*args, **kwargs)

        return wrapper

    monkeypatch.setattr(os, "scandir", counted("scandir", os.scandir))
    monkeypatch.setattr(os, "lstat", counted("lstat", os.lstat))
    monkeypatch.setattr(os, "stat", counted("stat", os.stat))
    restart = outputs.reconcile(db, workers=4)
    monkeypatch.undo()

    assert restart["trusted"] == 300 and outputs.judged == []
    assert calls["scandir"] == 1
    assert calls["lstat"] <= len(runs) + 5
    assert calls["stat"] <= 5
    db.close()
