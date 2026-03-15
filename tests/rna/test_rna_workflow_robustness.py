"""Tests for RNA workflow robustness and edge cases.

Follows NO_MOCKING_POLICY: all tests use real filesystem operations,
real configuration objects, and real function calls.
"""

import csv
import tempfile
import pytest
from pathlib import Path

from metainformant.rna.engine.workflow_core import (
    AmalgkitWorkflowConfig,
    WorkflowExecutionResult,
    WorkflowStepResult,
    validate_workflow_config,
)


class TestWorkflowRobustness:
    """Tests for workflow configuration validation and result handling."""

    def test_config_creation_from_dict(self, tmp_path: Path) -> None:
        """Test real config creation from dictionary."""
        config = AmalgkitWorkflowConfig(
            work_dir=tmp_path,
            threads=4,
            species_list=["species_A"],
        )
        assert config.work_dir == tmp_path
        assert config.threads == 4
        assert config.species_list == ["species_A"]

    def test_config_validation_valid(self, tmp_path: Path) -> None:
        """Test validation passes for a well-formed config."""
        config = AmalgkitWorkflowConfig(
            work_dir=tmp_path,
            threads=4,
            species_list=["species_A"],
        )
        is_valid, errors = validate_workflow_config(config)
        assert is_valid is True
        assert errors == []

    def test_config_validation_empty_species(self, tmp_path: Path) -> None:
        """Test validation fails when species_list is empty."""
        config = AmalgkitWorkflowConfig(
            work_dir=tmp_path,
            threads=4,
            species_list=[],
        )
        is_valid, errors = validate_workflow_config(config)
        assert is_valid is False
        assert any("species" in e.lower() for e in errors)

    def test_config_validation_invalid_threads(self, tmp_path: Path) -> None:
        """Test validation fails when threads < 1."""
        config = AmalgkitWorkflowConfig(
            work_dir=tmp_path,
            threads=0,
            species_list=["sp"],
        )
        is_valid, errors = validate_workflow_config(config)
        assert is_valid is False
        assert any("threads" in e.lower() for e in errors)

    def test_workflow_result_success(self) -> None:
        """Test WorkflowExecutionResult correctly reports success."""
        steps = [
            WorkflowStepResult(step_name="metadata", return_code=0, success=True),
            WorkflowStepResult(step_name="quant", return_code=0, success=True),
        ]
        result = WorkflowExecutionResult(
            steps_executed=steps,
            success=True,
            total_steps=2,
            successful_steps=2,
            failed_steps=0,
        )
        assert result.success is True
        assert len(result) == 2
        assert result.return_codes == [0, 0]

    def test_workflow_result_failure(self) -> None:
        """Test WorkflowExecutionResult correctly reports failure."""
        steps = [
            WorkflowStepResult(step_name="metadata", return_code=0, success=True),
            WorkflowStepResult(
                step_name="quant", return_code=1, success=False,
                error_message="Quantification failed",
            ),
        ]
        result = WorkflowExecutionResult(
            steps_executed=steps,
            success=False,
            total_steps=2,
            successful_steps=1,
            failed_steps=1,
        )
        assert result.success is False
        assert result.failed_steps == 1
        assert result.return_codes == [0, 1]

    def test_workflow_result_get_by_name(self) -> None:
        """Test WorkflowExecutionResult.get() retrieves step by name."""
        steps = [
            WorkflowStepResult(step_name="metadata", return_code=0, success=True),
            WorkflowStepResult(step_name="quant", return_code=1, success=False),
        ]
        result = WorkflowExecutionResult(
            steps_executed=steps,
            success=False,
            total_steps=2,
            successful_steps=1,
            failed_steps=1,
        )
        assert result.get("metadata") == 0
        assert result.get("quant") == 1
        assert result.get("nonexistent") is None
        assert result.get("nonexistent", -1) == -1

    def test_config_serialization_roundtrip(self, tmp_path: Path) -> None:
        """Test config can be serialized to dict and back."""
        original = AmalgkitWorkflowConfig(
            work_dir=tmp_path,
            threads=8,
            species_list=["Apis_mellifera"],
            search_string="RNA-Seq",
        )
        d = original.to_dict()
        restored = AmalgkitWorkflowConfig.from_dict(d)

        assert str(restored.work_dir) == str(original.work_dir)
        assert restored.threads == original.threads
        assert restored.species_list == original.species_list

    def test_metadata_file_creation(self, tmp_path: Path) -> None:
        """Test real metadata TSV file creation and parsing."""
        metadata_path = tmp_path / "metadata.tsv"
        with open(metadata_path, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=["run", "tissue"], delimiter="\t")
            writer.writeheader()
            writer.writerow({"run": "SRR_TEST_001", "tissue": "brain"})
            writer.writerow({"run": "SRR_TEST_002", "tissue": "antenna"})

        # Verify the file can be read back
        with open(metadata_path) as f:
            reader = csv.DictReader(f, delimiter="\t")
            rows = list(reader)

        assert len(rows) == 2
        assert rows[0]["run"] == "SRR_TEST_001"
        assert rows[1]["tissue"] == "antenna"
