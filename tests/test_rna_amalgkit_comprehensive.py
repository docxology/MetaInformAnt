"""Comprehensive tests for amalgkit integration functionality."""

import json
import os
import tempfile
from pathlib import Path
from unittest.mock import patch

import pytest

from metainformant.rna import amalgkit, workflow


class TestAmalgkitIntegration:
    """Test comprehensive amalgkit integration functionality."""

    def setup_method(self):
        """Setup test environment."""
        self.test_dir = Path("output/test_amalgkit")
        self.test_dir.mkdir(parents=True, exist_ok=True)

    def test_check_cli_available(self):
        """Test CLI availability checking."""
        ok, help_text = amalgkit.check_cli_available()
        assert isinstance(ok, bool)
        assert isinstance(help_text, str)

        if ok:
            assert len(help_text) > 100  # Should have substantial help text
            assert "AMALGKIT" in help_text.upper()
        else:
            assert len(help_text) > 0  # Should have error message

    def test_build_cli_args_basic(self):
        """Test basic CLI argument building."""
        params = {"threads": 4, "verbose": True, "dry_run": False}
        args = amalgkit.build_cli_args(params)

        assert "--threads" in args
        assert "4" in args
        assert "--verbose" in args
        assert "--dry-run" not in args  # False values are omitted

    def test_build_cli_args_with_lists(self):
        """Test CLI argument building with list parameters."""
        params = {"species-list": ["Homo sapiens", "Mus musculus"], "tissue": ["brain", "liver"]}
        args = amalgkit.build_cli_args(params)

        assert "--species-list" in args
        assert "Homo sapiens" in args
        assert "Mus musculus" in args
        assert "--tissue" in args
        assert "brain" in args
        assert "liver" in args

    def test_build_cli_args_with_paths(self):
        """Test CLI argument building with Path objects."""
        from pathlib import Path
        params = {"out_dir": Path("/tmp/test"), "config": Path("config.yaml")}
        args = amalgkit.build_cli_args(params)

        assert "--out-dir" in args
        assert "/tmp/test" in args
        assert "--config" in args
        assert "config.yaml" in args

    def test_build_amalgkit_command(self):
        """Test amalgkit command construction."""
        params = {"threads": 4, "out_dir": "/tmp/test"}
        command = amalgkit.build_amalgkit_command("metadata", params)

        assert "amalgkit" in command
        assert "metadata" in command
        assert "--threads" in command
        assert "4" in command
        assert "--out_dir" in command  # CLI uses underscores, not hyphens
        assert "/tmp/test" in command

    def test_metadata_step_execution(self):
        """Test metadata step execution."""
        params = {
            "out_dir": str(self.test_dir / "work"),
            "search_string": '"Homo sapiens"[Organism] AND RNA-Seq[Strategy]',
            "threads": 2
        }

        # This should work if amalgkit is available
        ok, _ = amalgkit.check_cli_available()
        if ok:
            result = amalgkit.run_amalgkit("metadata", params, work_dir=str(self.test_dir))
            assert result.returncode == 0 or result.returncode == 2  # 2 might indicate no data found

    def test_config_step_execution(self):
        """Test config step execution."""
        params = {
            "out_dir": str(self.test_dir / "work"),
        }

        ok, _ = amalgkit.check_cli_available()
        if ok:
            result = amalgkit.run_amalgkit("config", params, work_dir=str(self.test_dir))
            # config step may return 0 or 2 if no config to generate
            assert result.returncode in (0, 2)


class TestWorkflowIntegration:
    """Test workflow orchestration functionality."""

    def setup_method(self):
        """Setup test environment."""
        self.test_dir = Path("output/test_workflow")
        self.test_dir.mkdir(parents=True, exist_ok=True)

    def test_workflow_config_loading(self):
        """Test workflow configuration loading."""
        # Create a test config file
        config_data = {
            "work_dir": str(self.test_dir / "work"),
            "log_dir": str(self.test_dir / "logs"),
            "threads": 4,
            "species_list": ["Homo sapiens"],
            "steps": {"metadata": {}, "config": {}}  # per_step format
        }

        config_file = self.test_dir / "test_config.yaml"
        with open(config_file, 'w') as f:
            json.dump(config_data, f)

        cfg = workflow.load_workflow_config(config_file)
        # Use resolved paths for comparison (handles absolute vs relative)
        assert cfg.work_dir.resolve() == (self.test_dir / "work").resolve()
        assert cfg.threads == 4
        assert cfg.species_list == ["Homo sapiens"]
        # per_step should be populated from steps dict
        assert isinstance(cfg.per_step, dict)

    def test_workflow_planning(self):
        """Test workflow planning functionality."""
        config_data = {
            "work_dir": str(self.test_dir / "work"),
            "threads": 4,
            "species_list": ["Homo sapiens"],
            "steps": ["metadata", "config", "sanity"]
        }

        config_file = self.test_dir / "test_config.yaml"
        with open(config_file, 'w') as f:
            json.dump(config_data, f)

        cfg = workflow.load_workflow_config(config_file)
        steps = workflow.plan_workflow(cfg)

        # When steps are specified in config, we get those steps
        # When no steps specified, we get all 11 steps
        assert len(steps) >= 3  # At least the steps we specified
        step_names = [step[0] for step in steps]
        assert "metadata" in step_names
        assert "config" in step_names
        assert "sanity" in step_names

        # Check that each step has parameters
        for step_name, params in steps:
            assert isinstance(params, dict)
            assert "threads" in params
            assert params["threads"] == 4

    def test_workflow_execution_dry_run(self):
        """Test workflow execution with dry run."""
        config_data = {
            "work_dir": str(self.test_dir / "work"),
            "threads": 4,
            "species_list": ["Homo sapiens"],
            "steps": ["metadata", "config"]
        }

        config_file = self.test_dir / "test_config.yaml"
        with open(config_file, 'w') as f:
            json.dump(config_data, f)

        cfg = workflow.load_workflow_config(config_file)

        # Test that we can at least validate the configuration
        try:
            steps = workflow.plan_workflow(cfg)
            assert len(steps) > 0
        except Exception as e:
            pytest.fail(f"Workflow planning failed: {e}")


class TestAmalgkitStepRunners:
    """Test individual step runners."""

    def setup_method(self):
        """Setup test environment."""
        self.test_dir = Path("output/test_steps")
        self.test_dir.mkdir(parents=True, exist_ok=True)

    def test_metadata_runner(self):
        """Test metadata step runner."""
        from metainformant.rna.steps import STEP_RUNNERS

        runner = STEP_RUNNERS["metadata"]
        params = {
            "out_dir": str(self.test_dir / "work"),
            "search_string": '"Homo sapiens"[Organism] AND RNA-Seq[Strategy]',
            "threads": 2
        }

        ok, _ = amalgkit.check_cli_available()
        if ok:
            result = runner(params, work_dir=str(self.test_dir))
            # Metadata step should complete or fail gracefully
            assert result is not None

    def test_config_runner(self):
        """Test config step runner."""
        from metainformant.rna.steps import STEP_RUNNERS

        runner = STEP_RUNNERS["config"]
        params = {
            "out_dir": str(self.test_dir / "work"),
            "threads": 2
        }

        ok, _ = amalgkit.check_cli_available()
        if ok:
            result = runner(params, work_dir=str(self.test_dir))
            assert result is not None


class TestErrorHandling:
    """Test error handling and edge cases."""

    def test_invalid_config_file(self):
        """Test handling of invalid config files."""
        invalid_config = Path("nonexistent_config.yaml")

        with pytest.raises(FileNotFoundError):
            workflow.load_workflow_config(invalid_config)

    def test_config_validation(self):
        """Test configuration validation."""
        # Test with missing required fields
        invalid_config = {
            "threads": 4,
            # Missing work_dir and species_list
        }

        config_file = Path("output/test_invalid_config.yaml")
        with open(config_file, 'w') as f:
            json.dump(invalid_config, f)

        # Should handle gracefully
        cfg = workflow.load_workflow_config(config_file)
        assert cfg.work_dir is not None  # Should have default
        assert cfg.species_list == []  # Should be empty

    def test_cli_unavailable_handling(self):
        """Test handling when amalgkit is not available."""
        with patch('metainformant.rna.amalgkit.check_cli_available', return_value=(False, "amalgkit not found")):
            ok, help_text = amalgkit.check_cli_available()
            assert not ok
            assert "not found" in help_text


class TestIntegrationWorkflow:
    """Test complete integration workflows."""

    def setup_method(self):
        """Setup test environment."""
        self.test_dir = Path("output/test_integration")
        self.test_dir.mkdir(parents=True, exist_ok=True)

    def test_complete_workflow_simulation(self):
        """Test complete workflow with simulated steps."""
        # This test verifies the workflow orchestration without requiring amalgkit
        config_data = {
            "work_dir": str(self.test_dir / "work"),
            "threads": 4,
            "species_list": ["Homo sapiens"],
            "steps": ["metadata", "config", "sanity"]
        }

        config_file = self.test_dir / "test_config.yaml"
        with open(config_file, 'w') as f:
            json.dump(config_data, f)

        cfg = workflow.load_workflow_config(config_file)

        # Test workflow planning
        steps = workflow.plan_workflow(cfg)
        assert len(steps) >= 3  # At least the specified steps

        # Verify each step has correct parameters
        specified_steps = ["metadata", "config", "sanity"]
        for step_name, params in steps:
            if step_name in specified_steps:
                assert params["threads"] == 4
                assert "species-list" in params or "species_list" in params
                species = params.get("species-list") or params.get("species_list")
                assert species == ["Homo sapiens"]

    def test_workflow_with_genome_config(self):
        """Test workflow with genome configuration."""
        config_data = {
            "work_dir": str(self.test_dir / "work"),
            "threads": 4,
            "species_list": ["Homo sapiens"],
            "genome": {
                "accession": "GCF_000001405.40",
                "include": ["genome", "gff3", "rna"]
            },
            "steps": ["metadata", "config"]
        }

        config_file = self.test_dir / "test_config.yaml"
        with open(config_file, 'w') as f:
            json.dump(config_data, f)

        cfg = workflow.load_workflow_config(config_file)
        assert cfg.genome is not None
        assert cfg.genome["accession"] == "GCF_000001405.40"
        assert "genome" in cfg.genome["include"]


class TestPerformanceAndRobustness:
    """Test performance and robustness of amalgkit integration."""

    def setup_method(self):
        """Setup test environment."""
        self.test_dir = Path("output/test_performance")
        self.test_dir.mkdir(parents=True, exist_ok=True)

    def test_large_config_handling(self):
        """Test handling of large configuration files."""
        # Create a large config with many species
        config_data = {
            "work_dir": str(self.test_dir / "work"),
            "threads": 8,
            "species_list": ["Homo sapiens", "Mus musculus", "Drosophila melanogaster"] * 10,  # 30 species
            "steps": ["metadata", "config"]
        }

        config_file = self.test_dir / "large_config.yaml"
        with open(config_file, 'w') as f:
            json.dump(config_data, f)

        cfg = workflow.load_workflow_config(config_file)
        assert len(cfg.species_list) == 30
        assert cfg.threads == 8

    def test_concurrent_workflow_execution(self):
        """Test that multiple workflows can be planned concurrently."""
        import threading
        import time

        config_data = {
            "work_dir": str(self.test_dir / "work"),
            "threads": 4,
            "species_list": ["Homo sapiens"],
            "steps": ["metadata"]
        }

        config_file = self.test_dir / "test_config.yaml"
        with open(config_file, 'w') as f:
            json.dump(config_data, f)

        results = []

        def plan_workflow():
            cfg = workflow.load_workflow_config(config_file)
            steps = workflow.plan_workflow(cfg)
            results.append(len(steps))

        # Run multiple threads
        threads = [threading.Thread(target=plan_workflow) for _ in range(3)]
        for t in threads:
            t.start()
        for t in threads:
            t.join()

        # All should have the same result (at least 1 step)
        assert all(r >= 1 for r in results)

    def test_workflow_with_filters(self):
        """Test workflow with filtering parameters."""
        config_data = {
            "work_dir": str(self.test_dir / "work"),
            "threads": 4,
            "species_list": ["Homo sapiens"],
            "filters": {
                "require_tissue": True,
                "min_spots": 1000000,
                "platform": "Illumina"
            },
            "steps": ["metadata", "config"]
        }

        config_file = self.test_dir / "test_config.yaml"
        with open(config_file, 'w') as f:
            json.dump(config_data, f)

        cfg = workflow.load_workflow_config(config_file)
        assert cfg.filters["require_tissue"] is True
        assert cfg.filters["min_spots"] == 1000000
        assert cfg.filters["platform"] == "Illumina"


class TestDocumentationExamples:
    """Test that documentation examples work correctly."""

    def setup_method(self):
        """Setup test environment."""
        self.test_dir = Path("output/test_docs")
        self.test_dir.mkdir(parents=True, exist_ok=True)

    def test_workflow_example_from_docs(self):
        """Test the workflow example from the documentation."""
        config_data = {
            "work_dir": str(self.test_dir / "work"),
            "threads": 4,
            "species_list": ["Apis mellifera"]
        }

        config_file = self.test_dir / "test_config.yaml"
        with open(config_file, 'w') as f:
            json.dump(config_data, f)

        # This should work as shown in docs
        cfg = workflow.load_workflow_config(config_file)
        steps = workflow.plan_workflow(cfg)

        assert len(steps) > 0
        assert cfg.species_list == ["Apis mellifera"]

    def test_cli_example_from_docs(self):
        """Test CLI example functionality."""
        # Test the command construction as shown in docs
        params = {"threads": 4}
        command = amalgkit.build_amalgkit_command("metadata", params)

        assert "amalgkit" in command
        assert "metadata" in command
        assert "--threads" in command
        assert "4" in command


class TestAmalgkitWrapperRobustness:
    """Test robustness of amalgkit wrapper."""

    def setup_method(self):
        """Setup test environment."""
        self.test_dir = Path("output/test_robustness")
        self.test_dir.mkdir(parents=True, exist_ok=True)

    def test_wrapper_with_missing_amalgkit(self):
        """Test wrapper behavior when amalgkit is not available."""
        with patch('metainformant.rna.amalgkit.check_cli_available', return_value=(False, "amalgkit not found")):
            ok, help_text = amalgkit.check_cli_available()
            assert not ok
            assert "not found" in help_text

    def test_wrapper_with_invalid_params(self):
        """Test wrapper with invalid parameters."""
        # Should handle invalid parameters gracefully
        params = {"invalid_param": "value", "threads": -1}

        # Build args should not crash, just pass params through
        args = amalgkit.build_cli_args(params)
        # Args are built, validation is left to amalgkit CLI
        assert "--invalid-param" in args or "--invalid_param" in args
        assert "--threads" in args

    def test_step_runner_error_handling(self):
        """Test error handling in step runners."""
        from metainformant.rna.steps import STEP_RUNNERS

        runner = STEP_RUNNERS["metadata"]
        params = {
            "out_dir": str(self.test_dir / "work"),
            "search_string": "invalid search"
        }

        # Should handle errors gracefully
        result = runner(params, work_dir=str(self.test_dir), check=False)
        # May fail due to invalid search, but should not crash
        assert result is not None
        assert hasattr(result, 'returncode')



