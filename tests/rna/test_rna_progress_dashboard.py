"""Tests for the RNA pipeline progress dashboard (descriptive-only)."""

from __future__ import annotations

import sqlite3
from pathlib import Path

import matplotlib.pyplot as plt  # noqa: E402  (progress_dashboard sets Agg on import)
import pytest

from metainformant.rna.engine.progress_dashboard import (
    load_counts,
    load_failed_details,
    plot_overall_donut,
    plot_species_bars,
)
from metainformant.rna.engine.progress_db import ProgressDB


def _build_progress_db(db_path: Path) -> ProgressDB:
    """Create a real progress DB with one quantified, one failed, one pending sample."""
    db = ProgressDB(db_path=db_path)
    db.init_species("apis", ["SRR1", "SRR2", "SRR3"])
    db.set_state("apis", "SRR1", "quantified")
    db.set_state("apis", "SRR2", "failed", error="Quantification Failed")
    return db


def test_load_counts_from_real_progress_db(tmp_path: Path) -> None:
    """load_counts aggregates per-species state counts from a real ProgressDB file."""
    db_path = tmp_path / "progress.db"
    db = _build_progress_db(db_path)
    try:
        counts = load_counts(db_path)
    finally:
        db.close()

    assert counts == {"apis": {"quantified": 1, "failed": 1, "pending": 1}}


def test_load_failed_details_from_real_progress_db(tmp_path: Path) -> None:
    """load_failed_details returns the failed sample with species, id, and error text."""
    db_path = tmp_path / "progress.db"
    db = _build_progress_db(db_path)
    try:
        details = load_failed_details(db_path)
    finally:
        db.close()

    assert len(details) == 1
    row = details[0]
    assert row["species"] == "apis"
    assert row["srr_id"] == "SRR2"
    assert "Quantification Failed" in row["error"]


def test_load_counts_missing_db_raises_operational_error(tmp_path: Path) -> None:
    """sqlite3.connect creates an empty file for a missing DB path, then the
    samples-table query raises sqlite3.OperationalError (pre-existing semantic
    of load_counts; the module intentionally does not guard this)."""
    db_path = tmp_path / "nonexistent.db"
    with pytest.raises(sqlite3.OperationalError):
        load_counts(db_path)


def test_plot_functions_smoke_under_agg() -> None:
    """plot_species_bars and plot_overall_donut render without error under the Agg backend."""
    counts: dict[str, dict[str, int]] = {"apis": {"quantified": 1, "failed": 1}}
    fig, ax = plt.subplots()
    plot_species_bars(ax, counts)
    plot_overall_donut(ax, counts)
    plt.close(fig)
