"""Tests for amalgkit logging functionality.

Tests that run_amalgkit creates log files for stdout and stderr.
"""

from __future__ import annotations

from pathlib import Path


def test_run_amalgkit_writes_logs(tmp_path: Path):
    """Test that run_amalgkit writes log files for stdout and stderr."""
    from metainformant.rna.amalgkit import check_cli_available, run_amalgkit

    # Verify amalgkit is available (fixture ensures this)
    ok, _ = check_cli_available()
    assert ok, "amalgkit CLI must be available (ensured by fixture)"

    logs = tmp_path / "logs"
    res = run_amalgkit(
        "metadata",
        {"threads": 1},
        work_dir=tmp_path / "work",
        log_dir=logs,
        step_name="metadata",
        check=False,
    )
    # Accept non-zero if amalgkit returns an error; still assert logs are produced
    written = list(logs.glob("*.metadata.*.log"))
    assert any(p.name.endswith("stdout.log") for p in written)
    assert any(p.name.endswith("stderr.log") for p in written)
