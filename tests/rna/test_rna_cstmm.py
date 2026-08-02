"""Tests for amalgkit cstmm command wrapper.

Follows real-implementation policy: all tests use real build_cli_args() and
AmalgkitParams to verify command construction without replacing subprocess.
"""

from pathlib import Path

import pytest

from metainformant.rna.amalgkit.amalgkit import (
    AmalgkitParams,
    build_amalgkit_command,
    build_cli_args,
)


class TestAmalgkitCstmm:
    """Tests for cstmm command construction using real functions."""

    @pytest.fixture
    def temp_dir(self, tmp_path: Path) -> Path:
        return tmp_path

    def test_cstmm_basic(self, temp_dir: Path) -> None:
        """Test basic cstmm command generation."""
        params = AmalgkitParams(
            work_dir=temp_dir,
            threads=8,
            species_list=["Species_A"],
            orthogroup_table="groups.tsv",
        )

        command = build_amalgkit_command("cstmm", params)

        assert command[0] == "amalgkit"
        assert command[1] == "cstmm"
        assert "--out_dir" in command
        assert str(temp_dir) in command
        assert "--threads" not in command  # cstmm does not support threads
        assert "--orthogroup_table" in command
        assert "groups.tsv" in command
        assert "--species" not in command

    def test_cstmm_dict_params(self, temp_dir: Path) -> None:
        """Test cstmm with dictionary parameters."""
        params = {
            "work_dir": temp_dir,
            "threads": 4,  # Should be filtered for cstmm
            "dir_busco": "/path/to/busco",
        }

        args = build_cli_args(params, subcommand="cstmm")

        assert "--threads" not in args
        assert "--dir_busco" in args
        assert "/path/to/busco" in args

    def test_cstmm_omits_species_list_for_current_cli(self, temp_dir: Path) -> None:
        """Verify species lists are not leaked into the current cstmm CLI."""
        params = AmalgkitParams(
            work_dir=temp_dir,
            species_list=["Sp_A", "Sp_B", "Sp_C"],
        )

        args = build_cli_args(params, subcommand="cstmm")

        assert "--species" not in args
        assert all(species not in args for species in params.species_list)

    def test_cstmm_all_strings(self, temp_dir: Path) -> None:
        """Verify all command elements are strings for subprocess compatibility."""
        params = AmalgkitParams(work_dir=temp_dir, species_list=["Test"])

        command = build_amalgkit_command("cstmm", params)

        assert all(isinstance(c, str) for c in command)
