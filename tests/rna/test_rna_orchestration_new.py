"""Tests for the multi-species RNA-seq orchestrator.

Follows NO_MOCKING_POLICY: all tests use real YAML configurations,
real class instantiation, and real method calls.
"""

import pytest
from pathlib import Path

from metainformant.rna.engine.orchestration_multi_species import PipelineOrchestrator
from metainformant.rna.amalgkit.amalgkit import AmalgkitParams, build_amalgkit_command


class TestPipelineOrchestrator:
    """Tests for PipelineOrchestrator using real configs and methods."""

    @pytest.fixture
    def config_path(self, tmp_path: Path) -> Path:
        """Create a real YAML config file for testing."""
        config = tmp_path / "orchestration_config.yaml"
        config.write_text(
            """
work_dir: {work_dir}
threads: 4
run_integrate: false
species:
  - species_A
  - name: species_B
    threads: 2
    extra_params:
      redo: yes
""".format(work_dir=str(tmp_path / "output"))
        )
        return config

    @pytest.fixture
    def orchestrator(self, config_path: Path) -> PipelineOrchestrator:
        return PipelineOrchestrator(config_path)

    def test_initialization(self, orchestrator: PipelineOrchestrator, config_path: Path) -> None:
        """Test orchestrator initializes with real config values."""
        assert orchestrator.config_path == config_path
        assert "metadata" in orchestrator.steps
        assert "curate" in orchestrator.steps
        assert "getfastq" in orchestrator.steps
        assert "quant" in orchestrator.steps

    def test_load_config(self, orchestrator: PipelineOrchestrator) -> None:
        """Test configuration loads correctly from real YAML file."""
        config = orchestrator.config
        assert len(config["species"]) == 2
        assert config["species"][0] == "species_A"
        assert config["species"][1]["name"] == "species_B"
        assert config["threads"] == 4
        assert config.get("run_integrate") is False

    def test_config_not_found(self, tmp_path: Path) -> None:
        """Test FileNotFoundError for missing config."""
        bad_path = tmp_path / "nonexistent.yaml"
        with pytest.raises(FileNotFoundError):
            PipelineOrchestrator(bad_path)

    def test_step_list_completeness(self, orchestrator: PipelineOrchestrator) -> None:
        """Verify all expected pipeline steps are present in order."""
        expected = ["metadata", "integrate", "getfastq", "quant", "merge", "cstmm", "csca", "curate"]
        assert orchestrator.steps == expected

    def test_species_params_construction(self, orchestrator: PipelineOrchestrator) -> None:
        """Verify that AmalgkitParams can be constructed from config species entries."""
        config = orchestrator.config

        # Species A: simple string name
        sp_a = config["species"][0]
        params_a = AmalgkitParams(
            work_dir=orchestrator.work_dir,
            threads=config["threads"],
            species_list=[sp_a],
        )
        assert params_a.threads == 4
        assert params_a.species_list == ["species_A"]

        # Species B: dict with overrides
        sp_b = config["species"][1]
        params_b = AmalgkitParams(
            work_dir=orchestrator.work_dir,
            threads=sp_b.get("threads", config["threads"]),
            species_list=[sp_b["name"]],
            **sp_b.get("extra_params", {}),
        )
        assert params_b.threads == 2
        assert params_b.species_list == ["species_B"]
        # YAML loads 'yes' → True (bool), but as extra_params it passes through
        assert params_b.extra_params.get("redo") is not None

    def test_run_step_unknown_step(self, orchestrator: PipelineOrchestrator) -> None:
        """Verify run_step returns False for unknown step names."""
        params = AmalgkitParams(work_dir="test", species_list=["test"])
        result = orchestrator.run_step("nonexistent_step", params, "test_species")
        assert result is False

    def test_integrate_skipped_when_disabled(self, orchestrator: PipelineOrchestrator) -> None:
        """Verify integrate step is skipped when run_integrate is false."""
        params = AmalgkitParams(work_dir="test", species_list=["test"])
        # integrate should succeed (by skipping) when disabled in config
        result = orchestrator.run_step("integrate", params, "test_species")
        assert result is True

    def test_command_construction_for_metadata(self) -> None:
        """Verify build_amalgkit_command produces correct metadata command.

        Note: metadata subcommand does NOT include --threads per build_cli_args logic.
        """
        params = AmalgkitParams(
            work_dir="/tmp/test",
            threads=4,
            species_list=["Apis_mellifera"],
        )
        command = build_amalgkit_command("metadata", params)
        assert command[0] == "amalgkit"
        assert command[1] == "metadata"
        assert "--out_dir" in command
        # metadata subcommand intentionally excludes --threads
        assert "--threads" not in command
