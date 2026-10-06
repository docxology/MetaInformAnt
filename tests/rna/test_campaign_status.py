"""Real SQLite, filesystem, subprocess and reconciliation controls."""
from __future__ import annotations

import json
import sqlite3
import subprocess
import sys
from pathlib import Path

import pytest

from metainformant.rna.engine.campaign_status import Inventory, Observation, StatusError, load_inventory, markdown_tables, reconcile
from metainformant.rna.engine.campaign_status_io import REQUIRED_FILES, file_coverage, parse_probe, probe_script, read_database


def inventory() -> Inventory:
    return load_inventory(json.dumps({"species_count": 2, "task_count": 4, "species": [
        {"species": "ant_a", "tasks": [{"task_id": f"ant_a/SRR{i}", "accession": f"SRR{i}"} for i in (1, 2)]},
        {"species": "ant_b", "tasks": [{"task_id": f"ant_b/SRR{i}", "accession": f"SRR{i}"} for i in (3, 4)]}]}).encode())


def test_marginals_locked_precedence_and_transfer_gap() -> None:
    report = reconcile(inventory(), frozenset({"ant_a/SRR1", "ant_b/SRR3"}),
                       frozenset({"ant_a/SRR1", "ant_a/SRR2", "ant_b/SRR3"}),
                       (Observation("ant_a/SRR1", "quantifying", "worker"), Observation("ant_a/SRR2", "downloading", "worker")),
                       (Observation("ant_a/SRR1", "failed", "local"), Observation("outside/SRR5", "pending", "local")),
                       frozenset({"ant_a/SRR1"}), frozenset({"ant_b/SRR3"}), frozenset({"ant_b/SRR3"}))
    assert report.totals.eligible == 4
    assert report.totals.cloud["locked"] == 2
    assert report.totals.cloud["downloading"] == 1
    assert report.totals.cloud["unassigned"] == 1
    assert sum(report.totals.cloud.values()) == 4
    assert sum(report.totals.local.values()) == 4
    assert report.totals.coverage["transfer_gap"] == 1
    assert report.local_outside_inventory == 1
    assert "| TOTAL |" in markdown_tables(report)


def test_missing_worker_is_unknown_not_pending() -> None:
    report = reconcile(inventory(), frozenset(), frozenset({"ant_a/SRR1"}), (), (), frozenset(), frozenset(), frozenset())
    assert report.totals.cloud["worker_unknown"] == 1
    assert report.totals.cloud["pending"] == 0


@pytest.mark.parametrize("case", ["outside", "duplicate", "bad_state", "unassigned", "coverage_overlap"])
def test_reconciliation_refuses_incoherence(case: str) -> None:
    locked = frozenset({"outside/SRR9"}) if case == "outside" else frozenset()
    records = (Observation("ant_a/SRR1", "invented" if case == "bad_state" else "pending", "w"),)
    if case == "duplicate":
        records *= 2
    present = frozenset({"ant_a/SRR1"}) if case == "coverage_overlap" else frozenset()
    assigned = frozenset() if case == "unassigned" else frozenset({"ant_a/SRR1"})
    with pytest.raises(StatusError):
        reconcile(inventory(), locked, assigned, records, (), present, present, frozenset())


def test_invalid_inventory_refused() -> None:
    with pytest.raises(StatusError):
        Inventory(inventory().species, 2, 5).task_ids()


def test_real_readonly_sqlite_and_worker_probe(tmp_path: Path) -> None:
    root = tmp_path / "worker"
    root.mkdir()
    db = root / "pipeline_progress.db"
    with sqlite3.connect(db) as connection:
        connection.execute("CREATE TABLE samples(species TEXT,srr_id TEXT,state TEXT)")
        connection.executemany("INSERT INTO samples VALUES(?,?,?)", [("ant_a", "SRR1", "downloading"), ("ant_a", "SRR2", "pending")])
    before = db.read_bytes()
    assert len(read_database(db)) == 2
    result = subprocess.run([sys.executable, "-c", probe_script(frozenset({"ant_a/SRR1"}), str(root))], check=True, capture_output=True, text=True)
    probe = parse_probe(result.stdout)
    assert [(r.task_id, r.state) for r in probe.rows] == [("ant_a/SRR1", "downloading")]
    assert db.read_bytes() == before
    with pytest.raises(sqlite3.OperationalError):
        read_database(root / "absent.db")
    assert not (root / "absent.db").exists()


def test_real_file_presence_and_partial_output(tmp_path: Path) -> None:
    quant = tmp_path / "ant_a/work/quant"
    for run in ("SRR1", "SRR2"):
        (quant / run).mkdir(parents=True)
    for name in REQUIRED_FILES:
        (quant / "SRR1" / name).write_text(json.dumps({"fixture": "presence-only"}))
    (quant / "SRR2/abundance.tsv").write_text("target_id\tlength\n")
    complete, partial = file_coverage(inventory(), tmp_path)
    assert complete == frozenset({"ant_a/SRR1"})
    assert partial == frozenset({"ant_a/SRR2"})
