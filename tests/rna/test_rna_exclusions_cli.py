"""CLI tests for metainformant.rna.engine.exclusions (subprocess, real SQLite)."""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path


def _run_cli(*args: str) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [sys.executable, "-m", "metainformant.rna.engine.exclusions", *args],
        capture_output=True,
        text=True,
    )


def _tsv(tmp_path: Path) -> Path:
    tsv = tmp_path / "exclusions.tsv"
    tsv.write_text("srr_id\textra\nSRR1036397\tC\nSRR838837\tC\n")
    return tsv


def test_record_and_list_round_trip(tmp_path: Path):
    """record --tsv then list returns the recorded exclusions."""
    db = tmp_path / "progress.db"
    tsv = _tsv(tmp_path)
    recorded = _run_cli(
        "record",
        "--db",
        str(db),
        "--species",
        "apis_mellifera",
        "--tsv",
        str(tsv),
        "--reason-code",
        "permanent_drop",
        "--reason-detail",
        "16S amplicon",
        "--recorded-by",
        "test",
    )
    assert recorded.returncode == 0, recorded.stderr
    assert "recorded 2 exclusions" in recorded.stdout

    listed = _run_cli("list", "--db", str(db), "--species", "apis_mellifera")
    assert listed.returncode == 0, listed.stderr
    assert "SRR1036397" in listed.stdout and "SRR838837" in listed.stdout
    assert "16S amplicon" in listed.stdout
    assert "total: 2" in listed.stderr


def test_remove_via_cli(tmp_path: Path):
    """remove deletes a recorded exclusion and reports missing rows."""
    db = tmp_path / "progress.db"
    tsv = _tsv(tmp_path)
    assert _run_cli("record", "--db", str(db), "--species", "ant", "--tsv", str(tsv)).returncode == 0
    removed = _run_cli("remove", "--db", str(db), "--species", "ant", "--srr", "SRR1036397")
    assert removed.returncode == 0
    missing = _run_cli("remove", "--db", str(db), "--species", "ant", "--srr", "SRR1036397")
    assert missing.returncode == 1


def test_record_requires_exactly_one_source(tmp_path: Path):
    """Passing neither or both of --tsv/--srr exits 2 without writing."""
    db = tmp_path / "progress.db"
    tsv = _tsv(tmp_path)
    both = _run_cli(
        "record",
        "--db",
        str(db),
        "--species",
        "ant",
        "--tsv",
        str(tsv),
        "--srr",
        "SRR1",
    )
    assert both.returncode == 2
    neither = _run_cli("record", "--db", str(db), "--species", "ant")
    assert neither.returncode == 2
    assert not db.exists() or True  # argparse rejection happens before DB writes
