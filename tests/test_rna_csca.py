"""Tests for amalgkit csca command wrapper.

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


class TestAmalgkitCsca:
    """Tests for csca command construction using real functions."""

    @pytest.fixture
    def temp_dir(self, tmp_path: Path) -> Path:
        return tmp_path

    def test_csca_basic(self, temp_dir: Path) -> None:
        """Test basic csca command generation."""
        params = AmalgkitParams(
            work_dir=temp_dir,
            threads=8,
            species_list=["Species_A"],
            batch_effect_alg="sva",
        )

        command = build_amalgkit_command("csca", params)

        assert command[0] == "amalgkit"
        assert command[1] == "csca"
        assert "--out_dir" in command
        assert str(temp_dir) in command
        assert "--threads" not in command  # csca does not support threads
        assert "--batch_effect_alg" in command
        assert "sva" in command
        assert "--species" in command
        assert "Species_A" in command

    def test_csca_dict_params(self, temp_dir: Path) -> None:
        """Test csca with dictionary parameters."""
        params = {
            "work_dir": temp_dir,
            "threads": 4,  # Should be filtered for csca
            "batch_effect_alg": "combatseq",
            "dir_busco": "/path/to/busco",
        }

        args = build_cli_args(params, subcommand="csca")

        assert "--threads" not in args
        assert "--batch_effect_alg" in args
        assert "combatseq" in args
        assert "--dir_busco" in args
        assert "/path/to/busco" in args

    def test_csca_species_list_handling(self, temp_dir: Path) -> None:
        """Verify species list is correctly passed to CLI args."""
        params = AmalgkitParams(
            work_dir=temp_dir,
            species_list=["Sp_A", "Sp_B"],
        )

        args = build_cli_args(params, subcommand="csca")

        assert "--species" in args
        assert "Sp_A" in args
        assert "Sp_B" in args

    def test_csca_all_strings(self, temp_dir: Path) -> None:
        """Verify all command elements are strings for subprocess compatibility."""
        params = AmalgkitParams(work_dir=temp_dir, species_list=["Test"])

        command = build_amalgkit_command("csca", params)

        assert all(isinstance(c, str) for c in command)
