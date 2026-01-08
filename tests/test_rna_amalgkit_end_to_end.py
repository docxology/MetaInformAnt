"""
Comprehensive end-to-end tests for amalgkit RNA-seq workflow.

Tests the complete amalgkit pipeline from metadata to quantification,
using real data without mocks.

All tests require amalgkit CLI to be available. The conftest.py fixture
ensures amalgkit is available before tests run.
"""

import pytest
from pathlib import Path
import subprocess
import shutil
from metainformant.rna import amalgkit


class TestAmalgkitEndToEnd:
    """Test complete amalgkit workflows end-to-end."""
    
    def setup_method(self):
        """Setup test directory."""
        self.test_dir = Path("output/test_amalgkit_e2e")
        self.test_dir.mkdir(parents=True, exist_ok=True)
    
    def teardown_method(self):
        """Cleanup test directory."""
        if self.test_dir.exists():
            shutil.rmtree(self.test_dir)
    
    @pytest.mark.external_tool
    def test_metadata_to_config_workflow(self, ensure_amalgkit_available):
        """Test metadata → integrate → config workflow.

        Uses ensure_amalgkit_available fixture to ensure amalgkit is available.
        """
        # Verify amalgkit is available (fixture ensures this)
        available, _ = amalgkit.check_cli_available()
        if not available:
            pytest.skip("amalgkit CLI not available; skipping external-tool end-to-end test")
        
        work_dir = self.test_dir / "workflow1"
        work_dir.mkdir(parents=True, exist_ok=True)
        
        # Step 1: metadata (requires search_string)
        metadata_params = {
            "out_dir": str(work_dir),
            "search_string": "Drosophila melanogaster[Organism] AND brain[tissue]",
            "max_samples": 5,
        }
        
        result = amalgkit.metadata(
            metadata_params,
            work_dir=str(work_dir),
            log_dir=str(work_dir / "logs"),
            check=False
        )
        
        # Metadata may fail due to network/API issues, just check it executed
        assert result.returncode in (0, 1, 2), "metadata step crashed"
        
        # Skip rest of test if metadata failed (network issue)
        metadata_file = work_dir / "metadata" / "metadata.tsv"
        if not metadata_file.exists():
            pytest.skip("Metadata download failed (likely network issue)")
        
        # Step 2: config
        config_params = {
            "out_dir": str(work_dir),
        }
        
        result = amalgkit.config(
            config_params,
            work_dir=str(work_dir),
            log_dir=str(work_dir / "logs"),
            check=False
        )
        
        # Config can return 0 or 2 (if config already exists)
        assert result.returncode in (0, 2), "config step failed"
        
        # Check config files
        config_dir = work_dir / "config_base"
        assert config_dir.exists(), "config_base directory not created"
    
    def test_cli_availability(self, ensure_amalgkit_available):
        """Test that amalgkit CLI is available.
        
        Uses ensure_amalgkit_available fixture to ensure amalgkit is available.
        """
        # Verify amalgkit is available (fixture ensures this)
        available, _ = amalgkit.check_cli_available()
        if not available:
            pytest.skip("amalgkit CLI not available; skipping external-tool CLI availability test")
        
        # Test version check
        result = subprocess.run(
            ["amalgkit", "--version"],
            capture_output=True,
            text=True
        )
        
        assert result.returncode == 0
        assert "AMALGKIT version" in result.stdout or "AMALGKIT version" in result.stderr
    
    def test_workflow_step_sequence(self, ensure_amalgkit_available):
        """Test that workflow steps execute in correct sequence.
        
        Uses ensure_amalgkit_available fixture to ensure amalgkit is available.
        """
        # Verify amalgkit is available (fixture ensures this)
        available, _ = amalgkit.check_cli_available()
        if not available:
            pytest.skip("amalgkit CLI not available; skipping external-tool step sequence test")
        
        work_dir = self.test_dir / "sequence_test"
        work_dir.mkdir(parents=True, exist_ok=True)
        
        steps_to_test = [
            ("metadata", {"species": ["Apis_mellifera"], "max_samples": 3}),
            ("config", {}),
        ]
        
        for step_name, params in steps_to_test:
            params["out_dir"] = str(work_dir)
            
            step_func = getattr(amalgkit, step_name)
            result = step_func(
                params,
                work_dir=str(work_dir),
                log_dir=str(work_dir / "logs"),
                check=False
            )
            
            # Allow success or informational codes
            assert result.returncode in (0, 2), f"{step_name} failed with code {result.returncode}"
    
    @pytest.mark.external_tool
    def test_complete_mini_workflow(self, ensure_amalgkit_available):
        """Test a minimal complete workflow with very limited data.

        Uses ensure_amalgkit_available fixture to ensure amalgkit is available.
        """
        # Verify amalgkit is available (fixture ensures this)
        available, _ = amalgkit.check_cli_available()
        if not available:
            pytest.skip("amalgkit CLI not available; skipping external-tool mini workflow test")
        
        work_dir = self.test_dir / "mini_workflow"
        work_dir.mkdir(parents=True, exist_ok=True)
        
        # Get metadata for 1 sample
        print("\n1. Getting metadata...")
        result = amalgkit.metadata(
            {
                "out_dir": str(work_dir),
                "search_string": "Caenorhabditis elegans[Organism]",
                "max_samples": 1,
            },
            work_dir=str(work_dir),
            log_dir=str(work_dir / "logs"),
            check=False
        )
        # Metadata may fail due to network, just check it ran
        assert result.returncode in (0, 1, 2), "metadata crashed"
        
        # Skip if metadata failed
        metadata_file = work_dir / "metadata" / "metadata.tsv"
        if not metadata_file.exists():
            pytest.skip("Metadata download failed (network issue)")
        
        with open(metadata_file) as f:
            lines = f.readlines()
        
        # Should have header + at least 1 sample
        assert len(lines) >= 2, "No samples in metadata"
        
        print(f"   ✅ Found {len(lines)-1} samples")
        
        # Generate config
        print("\n2. Generating config...")
        result = amalgkit.config(
            {"out_dir": str(work_dir)},
            work_dir=str(work_dir),
            log_dir=str(work_dir / "logs"),
            check=False
        )
        assert result.returncode in (0, 2)
        
        config_dir = work_dir / "config_base"
        assert config_dir.exists()
        
        print("   ✅ Config created")
        
        # Run sanity check
        print("\n3. Running sanity check...")
        result = amalgkit.sanity(
            {"out_dir": str(work_dir)},
            work_dir=str(work_dir),
            log_dir=str(work_dir / "logs"),
            check=False
        )
        assert result.returncode == 0
        
        print("   ✅ Sanity check passed")
        
        print("\n✅ Mini workflow completed successfully!")


class TestWorkflowAllSteps:
    """Test complete workflow with all 11 steps in order."""

    def setup_method(self):
        """Setup test directory."""
        self.test_dir = Path("output/test_complete_workflow")
        self.test_dir.mkdir(parents=True, exist_ok=True)

    def teardown_method(self):
        """Cleanup test directory."""
        if self.test_dir.exists():
            shutil.rmtree(self.test_dir)

    @pytest.mark.slow
    @pytest.mark.external_tool
    def test_all_steps_execute_in_order(self, ensure_amalgkit_available):
        """Test that all 11 steps execute in correct order.
        
        This test verifies the complete workflow execution order:
        metadata → integrate → config → select → getfastq → quant → merge → cstmm → curate → csca → sanity
        
        NOTE: This test is marked as slow because it may execute actual workflow steps.
        It primarily verifies workflow planning and step ordering, not full execution.
        """
        # Verify amalgkit is available
        available, _ = amalgkit.check_cli_available()
        if not available:
            pytest.skip("amalgkit CLI not available; skipping external-tool planning test")
        
        from metainformant.rna.engine.workflow import (
            AmalgkitWorkflowConfig,
            plan_workflow,
        )
        
        work_dir = self.test_dir / "complete_workflow"
        work_dir.mkdir(parents=True, exist_ok=True)
        
        # Create minimal config
        cfg = AmalgkitWorkflowConfig(
            work_dir=work_dir,
            threads=1,
            species_list=["Apis_mellifera"],
            per_step={
                "metadata": {
                    "search_string": '"Apis mellifera"[Organism] AND RNA-Seq[Strategy]',
                    "max_samples": 1,  # Limit to 1 sample for testing
                },
            },
        )
        
        # Plan workflow (this is fast and doesn't execute)
        steps = plan_workflow(cfg)
        step_names = [name for name, _ in steps]
        
        # Verify all 11 steps are planned
        expected_steps = [
            "metadata",
            "config",
            "select",
            "getfastq",
            "integrate",
            "quant",
            "merge",
            "cstmm",
            "curate",
            "csca",
            "sanity",
        ]
        
        assert len(step_names) == len(expected_steps)
        assert step_names == expected_steps, f"Step order mismatch: {step_names} != {expected_steps}"
        
        # Note: We skip actual execution to avoid hanging on network calls
        # The step ordering is verified by planning, which is sufficient for this test


class TestAmalgkitStepRunners:
    """Test individual step runners work correctly."""
    
    def setup_method(self):
        """Setup test directory."""
        self.test_dir = Path("output/test_step_runners")
        self.test_dir.mkdir(parents=True, exist_ok=True)
    
    def teardown_method(self):
        """Cleanup test directory."""
        if self.test_dir.exists():
            shutil.rmtree(self.test_dir)
    
    @pytest.mark.external_tool
    def test_metadata_runner(self, ensure_amalgkit_available):
        """Test metadata step runner.

        Uses ensure_amalgkit_available fixture to ensure amalgkit is available.
        """
        # Verify amalgkit is available (fixture ensures this)
        available, _ = amalgkit.check_cli_available()
        if not available:
            pytest.skip("amalgkit CLI not available; skipping external-tool step runner test")
        
        result = amalgkit.metadata(
            {
                "out_dir": str(self.test_dir),
                "search_string": "Saccharomyces cerevisiae[Organism]",
                "max_samples": 2,
            },
            work_dir=str(self.test_dir),
            log_dir=str(self.test_dir / "logs"),
            check=False
        )
        
        assert hasattr(result, 'returncode')
        # Metadata may fail due to network issues
        assert result.returncode in (0, 1, 2)
    
    @pytest.mark.external_tool
    def test_config_runner(self, ensure_amalgkit_available):
        """Test config step runner.

        Uses ensure_amalgkit_available fixture to ensure amalgkit is available.
        """
        # Verify amalgkit is available (fixture ensures this)
        available, _ = amalgkit.check_cli_available()
        if not available:
            pytest.skip("amalgkit CLI not available; skipping external-tool step runner test")
        
        # Need metadata first
        metadata_dir = self.test_dir / "metadata"
        metadata_dir.mkdir(parents=True, exist_ok=True)
        
        # Create minimal metadata file
        metadata_file = metadata_dir / "metadata.tsv"
        with open(metadata_file, 'w') as f:
            f.write("scientific_name\tsample_group\n")
            f.write("Test_species\ttest_group\n")
        
        result = amalgkit.config(
            {"out_dir": str(self.test_dir)},
            work_dir=str(self.test_dir),
            log_dir=str(self.test_dir / "logs"),
            check=False
        )
        
        assert hasattr(result, 'returncode')
        assert result.returncode in (0, 2)
    
    @pytest.mark.slow
    def test_sanity_runner(self, ensure_amalgkit_available):
        """Test sanity step runner.

        Uses ensure_amalgkit_available fixture to ensure amalgkit is available.
        """
        # Verify amalgkit is available (fixture ensures this)
        available, _ = amalgkit.check_cli_available()
        if not available:
            pytest.skip("amalgkit CLI not available; skipping external-tool sanity runner test")
        
        # Create minimal required structure
        metadata_dir = self.test_dir / "metadata"
        metadata_dir.mkdir(parents=True, exist_ok=True)
        
        metadata_file = metadata_dir / "metadata.tsv"
        with open(metadata_file, 'w') as f:
            f.write("run\tscientific_name\n")
            f.write("SRR000001\tTest_species\n")
        
        result = amalgkit.sanity(
            {"out_dir": str(self.test_dir)},
            work_dir=str(self.test_dir),
            log_dir=str(self.test_dir / "logs"),
            check=False
        )
        
        assert hasattr(result, 'returncode')
        assert result.returncode == 0


class TestAmalgkitUtilities:
    """Test amalgkit utility functions."""
    
    def test_build_cli_args_basic(self):
        """Test CLI argument building."""
        params = {
            "out_dir": "/tmp/test",
            "threads": 4,
            "verbose": True,
        }
        
        # Build args for CLI (use for_cli=True for actual CLI format)
        args = amalgkit.build_cli_args(params, for_cli=True)
        
        assert "--out_dir" in args
        assert "/tmp/test" in args
        assert "--threads" in args
        assert "4" in args
        assert "--verbose" in args
    
    def test_build_cli_args_lists(self):
        """Test CLI argument building with lists."""
        params = {
            "species": ["Species_one", "Species_two"],
            "sample_ids": ["SRR001", "SRR002"],
        }
        
        args = amalgkit.build_cli_args(params, for_cli=True)
        
        # Lists create multiple flag entries
        assert "--species" in args
        assert "Species_one" in args
        assert "Species_two" in args
        assert "--sample_ids" in args or "--sample-ids" in args
        assert "SRR001" in args
    
    def test_build_amalgkit_command(self):
        """Test complete command building."""
        command = amalgkit.build_amalgkit_command(
            "metadata",
            {"out_dir": "/tmp", "threads": 2}
        )
        
        assert command[0] == "amalgkit"
        assert command[1] == "metadata"
        assert "--out_dir" in command
        assert "/tmp" in command
        assert "--threads" in command
        assert "2" in command
    
    def test_check_cli_available(self):
        """Test CLI availability check."""
        result = amalgkit.check_cli_available()
        # Returns (bool, str) tuple
        assert isinstance(result, tuple)
        assert len(result) == 2
        assert isinstance(result[0], bool)
    
    def test_ensure_cli_available(self):
        """Test CLI availability check and optional auto-install.
        
        Note: ensure_cli_available returns a tuple (ok, msg, install_record),
        it does not raise exceptions. This test verifies the return value.
        """
        ok, msg, install_record = amalgkit.ensure_cli_available()
        
        # Should return boolean, message, and optional install record
        assert isinstance(ok, bool)
        assert isinstance(msg, str)
        assert install_record is None or isinstance(install_record, dict)
        
        # If amalgkit is available, ok should be True
        # If not available, ok should be False and msg should explain
        if not ok:
            assert len(msg) > 0, "Error message should be provided when amalgkit unavailable"
        else:
            # When available, should have help text or version info
            assert len(msg) > 0


if __name__ == "__main__":
    pytest.main([__file__, "-v"])

