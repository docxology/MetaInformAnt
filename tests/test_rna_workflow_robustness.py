
import pytest
import csv
import tempfile
from pathlib import Path
from unittest.mock import Mock, patch
from metainformant.rna.engine.workflow_execution import _execute_streaming_mode, AmalgkitWorkflowConfig

class TestWorkflowRobustness:
    @patch("metainformant.rna.engine.workflow_execution.cleanup_fastqs")
    def test_cleanup_skipped_on_failure(self, mock_cleanup, tmp_path):
        """Verify that cleanup_fastqs is NOT called if a step fails."""
        
        # Setup mocks
        config = Mock(spec=AmalgkitWorkflowConfig)
        config.work_dir = tmp_path / "work_dir"
        config.work_dir.mkdir()
        (config.work_dir / "metadata").mkdir()
        
        # Mock step function that fails
        def fail_step(params, **kwargs):
            return Mock(returncode=1, stderr="Mock failure")
            
        # "quant" is recognized as a chunk step
        step_functions = {"quant": fail_step}
        
        # Run streaming mode
        with tempfile.NamedTemporaryFile(mode='w', delete=False) as f:
            writer = csv.DictWriter(f, fieldnames=["run"], delimiter="\t")
            writer.writeheader()
            writer.writerow({"run": "SRR_TEST"})
            meta_path = Path(f.name)
            
        try:
            _execute_streaming_mode(
                config=config,
                steps_remaining=[("quant", {})],
                metadata_file=meta_path,
                chunk_size=1,
                step_functions=step_functions,
                check=False
            )
            
            # Assert cleanup was NOT called
            mock_cleanup.assert_not_called()
            
        finally:
            if meta_path.exists():
                meta_path.unlink()

    @patch("metainformant.rna.engine.workflow_execution.cleanup_fastqs")
    def test_cleanup_called_on_success(self, mock_cleanup, tmp_path):
        """Verify that cleanup_fastqs IS called if a step succeeds."""
        
        config = Mock(spec=AmalgkitWorkflowConfig)
        config.work_dir = tmp_path / "work_dir"
        config.work_dir.mkdir()
        (config.work_dir / "metadata").mkdir()
        
        # Mock step function that succeeds
        def success_step(params, **kwargs):
            return Mock(returncode=0)
            
        step_functions = {"quant": success_step}
        
        with tempfile.NamedTemporaryFile(mode='w', delete=False) as f:
            writer = csv.DictWriter(f, fieldnames=["run"], delimiter="\t")
            writer.writeheader()
            writer.writerow({"run": "SRR_TEST"})
            meta_path = Path(f.name)
            
        try:
            _execute_streaming_mode(
                config=config,
                steps_remaining=[("quant", {})],
                metadata_file=meta_path,
                chunk_size=1,
                step_functions=step_functions,
                check=False
            )
            
            # Assert cleanup WAS called
            mock_cleanup.assert_called()
            
        finally:
            if meta_path.exists():
                meta_path.unlink()
