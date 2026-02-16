"""Tests for amalgkit curate command wrapper.

Follows NO_MOCKING_POLICY: all tests use real build_cli_args() and
AmalgkitParams to verify command construction without mocking subprocess.
"""

import pytest
from pathlib import Path

from metainformant.rna.amalgkit.amalgkit import (
    AmalgkitParams,
    build_amalgkit_command,
    build_cli_args,
)


class TestAmalgkitCurate:
    """Tests for amalgkit curate command construction using real functions."""

    @pytest.fixture
    def temp_dir(self, tmp_path: Path) -> Path:
        return tmp_path

    def test_curate_basic_command(self, temp_dir: Path) -> None:
        """Test basic curate command generation via build_amalgkit_command."""
        params = AmalgkitParams(
            work_dir=temp_dir,
            dist="correlation",
            eval="r2",
        )

        command = build_amalgkit_command("curate", params)

        assert command[0] == "amalgkit"
        assert command[1] == "curate"
        assert "--out_dir" in command
        assert str(temp_dir) in command
        assert "--dist" in command
        assert "correlation" in command
        assert "--eval" in command
        assert "r2" in command

    def test_curate_with_dict_params(self, temp_dir: Path) -> None:
        """Test curate with dictionary parameters."""
        params = {
            "work_dir": temp_dir,
            "threshold": 0.5,
            "skip_curation": True,
        }

        args = build_cli_args(params, subcommand="curate")

        assert "--out_dir" in args
        assert "--threshold" in args
        assert "0.5" in args
        assert "--skip_curation" in args
        assert "yes" in args

    def test_curate_threads_filtered(self, temp_dir: Path) -> None:
        """Verify curate does not pass --threads (R script uses internal parallelism)."""
        params = AmalgkitParams(work_dir=temp_dir, threads=4)

        args = build_cli_args(params, subcommand="curate")

        # curate should NOT contain --threads
        assert "--threads" not in args

    def test_curate_full_command_structure(self, temp_dir: Path) -> None:
        """Verify the full command list structure is subprocess-ready."""
        params = AmalgkitParams(work_dir=temp_dir)

        command = build_amalgkit_command("curate", params)

        # First two elements are always the executable and subcommand
        assert isinstance(command, list)
        assert len(command) >= 2
        assert command[0] == "amalgkit"
        assert command[1] == "curate"
        # All elements should be strings (subprocess requirement)
        assert all(isinstance(c, str) for c in command)
