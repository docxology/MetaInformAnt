"""Contract tests for the current ENA-first Amalgkit producer surface.

These tests exercise only command discovery and dry-run behavior. They do not
download data, mutate the campaign data root, or require external services.
"""

from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]
CONFIG_DIR = REPO_ROOT / "projects" / "hymenoptera_amalgkit" / "config" / "amalgkit"
RUN_ALL = REPO_ROOT / "scripts" / "rna" / "run_all_species.py"
PROCESS_ONE = REPO_ROOT / "scripts" / "rna" / "process_species.py"

pytestmark = pytest.mark.skipif(
    not CONFIG_DIR.is_dir(),
    reason="Hymenoptera Amalgkit submodule is unavailable in this parent-only checkout",
)


def _run(script: Path, *args: str) -> subprocess.CompletedProcess[str]:
    """Run a producer command with the repository source tree importable."""

    environment = os.environ.copy()
    source_path = str(REPO_ROOT / "src")
    environment["PYTHONPATH"] = os.pathsep.join(
        [source_path, environment["PYTHONPATH"]] if environment.get("PYTHONPATH") else [source_path]
    )
    return subprocess.run(
        [sys.executable, str(script), *args],
        cwd=REPO_ROOT,
        env=environment,
        capture_output=True,
        text=True,
        timeout=30,
        check=False,
    )


def test_current_producer_scripts_exist() -> None:
    """The cohort and single-species producers are present."""

    assert RUN_ALL.is_file()
    assert PROCESS_ONE.is_file()


def test_run_all_species_help_is_current() -> None:
    """The cohort producer exposes the current bounded execution options."""

    result = _run(RUN_ALL, "--help")
    assert result.returncode == 0, result.stderr
    assert "--config-dir" in result.stdout
    assert "--data-root" in result.stdout
    assert "--dry-run" in result.stdout
    assert "--max-in-flight" in result.stdout


def test_process_species_help_is_current() -> None:
    """The single-species producer shares the current option contract."""

    result = _run(PROCESS_ONE, "--help")
    assert result.returncode == 0, result.stderr
    assert "--species" in result.stdout
    assert "--config-dir" in result.stdout
    assert "--dry-run" in result.stdout


def test_run_all_species_dry_run_discovers_the_configured_cohort(tmp_path: Path) -> None:
    """A dry run discovers configs without contacting ENA or writing outputs."""

    result = _run(
        RUN_ALL,
        "--config-dir",
        str(CONFIG_DIR),
        "--data-root",
        str(tmp_path / "data"),
        "--dry-run",
    )
    assert result.returncode == 0, result.stderr
    assert "Config directory:" in result.stdout
    assert "Species configs:" in result.stdout
    assert "amalgkit_apis_mellifera.yaml" in result.stdout
    assert not (tmp_path / "data" / "pipeline_progress.db").exists()


def test_single_species_dry_run_resolves_a_current_config(tmp_path: Path) -> None:
    """A single-species dry run resolves one current config idempotently."""

    args = (
        "--species",
        "apis_mellifera",
        "--config-dir",
        str(CONFIG_DIR),
        "--data-root",
        str(tmp_path / "data"),
        "--dry-run",
    )
    first = _run(PROCESS_ONE, *args)
    second = _run(PROCESS_ONE, *args)
    assert first.returncode == second.returncode == 0
    assert first.stdout == second.stdout
    assert not (tmp_path / "data" / "pipeline_progress.db").exists()


def test_current_documentation_describes_ena_and_the_producers() -> None:
    """Current script documentation names the current ENA-first producers."""

    readme = (REPO_ROOT / "scripts" / "rna" / "README.md").read_text(encoding="utf-8")
    assert "run_all_species.py" in readme
    assert "process_species.py" in readme
    assert "ENA" in readme
