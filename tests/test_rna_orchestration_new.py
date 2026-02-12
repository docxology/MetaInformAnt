"""Tests for the multi-species RNA-seq orchestrator."""

import pytest
from pathlib import Path
from unittest.mock import MagicMock, patch

from metainformant.rna.engine.orchestration_multi_species import PipelineOrchestrator
from metainformant.rna.amalgkit import amalgkit

class TestPipelineOrchestrator:
    @pytest.fixture
    def mock_config(self, tmp_path):
        config_path = tmp_path / "config.yaml"
        with open(config_path, "w") as f:
            f.write("""
work_dir: output/test_amalgkit
threads: 4
run_integrate: false
species:
  - species_A
  - name: species_B
    threads: 2
    extra_params:
      redo: yes
""")
        return config_path

    @pytest.fixture
    def orchestrator(self, mock_config):
        return PipelineOrchestrator(mock_config)

    def test_initialization(self, orchestrator, mock_config):
        assert orchestrator.config_path == mock_config
        assert orchestrator.work_dir == Path("output/test_amalgkit")
        assert orchestrator.log_dir == Path("output/test_amalgkit/logs")
        assert "metadata" in orchestrator.steps
        assert "curate" in orchestrator.steps

    def test_load_config(self, orchestrator):
        config = orchestrator.config
        assert config["work_dir"] == "output/test_amalgkit"
        assert len(config["species"]) == 2
        assert config["species"][0] == "species_A"
        assert config["species"][1]["name"] == "species_B"

    @pytest.mark.parametrize("step_name,expected_call", [
        ("metadata", True),
        ("integrate", False), # run_integrate is false in mock_config
        ("getfastq", True),
    ])
    def test_run_step_logic(self, orchestrator, step_name, expected_call):
        # Mock the actual step function in amalgkit
        mock_step_func = MagicMock()
        mock_step_func.return_value.returncode = 0
        
        with patch.object(amalgkit, step_name, mock_step_func, create=True):
            params = amalgkit.AmalgkitParams(work_dir="test", species_list=["test"])
            result = orchestrator.run_step(step_name, params, "test_species")
            
            assert result is True
            if expected_call:
                mock_step_func.assert_called_once()
            else:
                mock_step_func.assert_not_called()

    def test_run_step_failure(self, orchestrator):
        mock_step_func = MagicMock()
        mock_step_func.return_value.returncode = 1
        mock_step_func.return_value.stderr = "Error message"
        
        with patch.object(amalgkit, "metadata", mock_step_func):
            params = amalgkit.AmalgkitParams(work_dir="test", species_list=["test"])
            result = orchestrator.run_step("metadata", params, "test_species")
            
            assert result is False
            mock_step_func.assert_called_once()

    def test_run_orchestration_flow(self, orchestrator):
        """Test the full run loop with mocked steps."""
        # Mock all steps
        with patch.object(orchestrator, "run_step") as mock_run_step:
            mock_run_step.return_value = True
            
            orchestrator.run()
            
            # Should be called for each species * each enabled step
            # Species A: metadata, getfastq, quant, merge, cstmm, csca, curate (integrate skipped) check implementation
            # logic in run() calls run_step for ALL steps in self.steps
            # run_step handles skipping internally
            
            expected_steps_per_species = len(orchestrator.steps)
            assert mock_run_step.call_count == 2 * expected_steps_per_species

    def test_run_orchestration_stops_on_failure(self, orchestrator):
        """Test that orchestration stops for a species if a step fails."""
        with patch.object(orchestrator, "run_step") as mock_run_step:
            # First step fails
            mock_run_step.side_effect = [False, True, True] * 10 
            
            orchestrator.run()
            
            # Should handle failure gracefully. 
            # In run() loop:
            # Species A: Step 1 (fails) -> break
            # Species B: Step 1 (fails - based on side_effect) -> break
            
            # Actually side_effect iterator is global.
            # Call 1: Species A, Step 1 -> False. Break A.
            # Call 2: Species B, Step 1 -> True. Continue B.
            # Call 3: Species B, Step 2 -> True. Continue B.
            
            # We want to verify it breaks. 
            # Let's start FRESH with a specific side effect logic
            
            def side_effect(step, params, species):
                if species == "species_A" and step == "metadata":
                    return False
                return True
            
            mock_run_step.side_effect = side_effect
            mock_run_step.reset_mock()
            
            orchestrator.run()
            
            # Species A should call only metadata
            # Species B should call all steps
            
            calls = mock_run_step.call_args_list
            species_a_calls = [c for c in calls if c[0][2] == "species_A"]
            species_b_calls = [c for c in calls if c[0][2] == "species_B"]
            
            assert len(species_a_calls) == 1
            assert species_a_calls[0][0][0] == "metadata"
            
            assert len(species_b_calls) == len(orchestrator.steps)
