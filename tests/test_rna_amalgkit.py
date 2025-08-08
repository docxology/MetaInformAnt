from __future__ import annotations

import os
import shutil
from pathlib import Path

import pytest


def test_build_cli_args_transforms_types():
    from metainformant.rna.amalgkit import build_cli_args

    params = {
        "db": "sra",
        "threads": 8,
        "dry_run": True,
        "optional": None,
        "species_list": ["Homo_sapiens", "Mus_musculus"],
        "out_dir": Path("/tmp/out"),
    }

    args = build_cli_args(params)

    # string/number
    assert "--db" in args and args[args.index("--db") + 1] == "sra"
    assert "--threads" in args and args[args.index("--threads") + 1] == "8"

    # booleans → flag present
    assert "--dry-run" in args

    # None → omitted
    assert "--optional" not in args

    # list → repeated flags
    species_idx = [i for i, v in enumerate(args) if v == "--species-list"]
    assert len(species_idx) == 2
    assert args[species_idx[0] + 1] == "Homo_sapiens"
    assert args[species_idx[1] + 1] == "Mus_musculus"

    # Path → string
    assert "--out-dir" in args and args[args.index("--out-dir") + 1] == "/tmp/out"


def test_build_amalgkit_command_prefix_and_order():
    from metainformant.rna.amalgkit import build_amalgkit_command

    cmd = build_amalgkit_command(
        subcommand="metadata",
        params={"db": "sra", "threads": 4},
    )

    assert cmd[:2] == ["amalgkit", "metadata"]
    assert "--db" in cmd and cmd[cmd.index("--db") + 1] == "sra"
    assert "--threads" in cmd and cmd[cmd.index("--threads") + 1] == "4"


@pytest.mark.skipif(shutil.which("amalgkit") is None, reason="amalgkit not installed")
def test_check_cli_available_runs_help():
    from metainformant.rna.amalgkit import check_cli_available

    ok, version_text = check_cli_available()
    assert ok is True
    assert isinstance(version_text, str) and len(version_text) > 0


def test_curate_summary_counts_from_fixture(tmp_path: Path):
    from metainformant.rna.pipeline import summarize_curate_tables

    # Use repo test data under tests/data/rna/curate/Apis_mellifera/tables
    repo_root = Path(__file__).resolve().parents[1]
    tables_dir = repo_root / "tests" / "data" / "rna" / "curate" / "Apis_mellifera" / "tables"

    counts = summarize_curate_tables(tables_dir)
    # At least the two expected TSVs should be counted
    assert counts.get("Apis_mellifera.metadata.tsv", 0) >= 1
    assert counts.get("Apis_mellifera.uncorrected.tc.tsv", 0) >= 1


