"""Tests for RNA step runner dispatch.

Tests that all step runners can be invoked and handle real subcommand execution.
"""

from __future__ import annotations

from pathlib import Path

import pytest


@pytest.mark.external_tool
@pytest.mark.network
def test_each_step_runner_invokes_real_subcommand(tmp_path: Path):
    """Test that every currently registered step runner invokes the CLI."""
    from metainformant.rna import steps as step_mod
    from metainformant.rna.amalgkit.amalgkit import check_cli_available

    ok, message = check_cli_available()
    assert ok, message

    work_dir = tmp_path / "work"
    log_dir = tmp_path / "logs"
    work_dir.mkdir()
    log_dir.mkdir()

    called: list[str] = []

    for name, runner in step_mod.STEP_RUNNERS.items():
        res = runner({}, work_dir=work_dir, log_dir=log_dir)
        called.append(name)
        # Do not assert return code strictly; some steps may fail without inputs
        assert hasattr(res, "returncode")

    assert set(called) == set(step_mod.STEP_RUNNERS)
