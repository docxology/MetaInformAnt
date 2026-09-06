"""Tests for RNA workflow robustness and edge cases.

Follows real-implementation policy: all tests use real filesystem operations,
real configuration objects, and real function calls.
"""

import csv
import json
import logging
from pathlib import Path

from metainformant.rna.engine.workflow_core import (
    AmalgkitWorkflowConfig,
    WorkflowExecutionResult,
    WorkflowStepResult,
    apply_config_defaults,
    create_sample_config,
    load_workflow_config,
    validate_workflow_config,
    validate_workflow_outputs,
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
                step_name="quant",
                return_code=1,
                success=False,
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


class TestWorkflowCoreValidation:
    """Tests for never-covered public functions in workflow_core."""

    def test_validate_workflow_outputs_missing_files(self, tmp_path: Path) -> None:
        """Validation fails with one error per missing expected output file."""
        config = AmalgkitWorkflowConfig(work_dir=str(tmp_path))

        is_valid, errors = validate_workflow_outputs(config)

        assert is_valid is False
        assert len(errors) == 3
        assert all("missing" in error.lower() for error in errors)
        assert any("metadata.tsv" in error for error in errors)
        assert any("expression_matrix.tsv" in error for error in errors)
        assert any("sanity_check.txt" in error for error in errors)

    def test_validate_workflow_outputs_all_present(self, tmp_path: Path) -> None:
        """Validation passes once all three expected output files exist on disk."""
        config = AmalgkitWorkflowConfig(work_dir=str(tmp_path))
        (tmp_path / "metadata.tsv").write_text("run\ttissue\nSRR_TEST_001\tbrain\n")
        (tmp_path / "expression_matrix.tsv").write_text("gene\tSRR_TEST_001\ngene1\t10\n")
        (tmp_path / "sanity_check.txt").write_text("ok\n")

        is_valid, errors = validate_workflow_outputs(config)

        assert is_valid is True
        assert errors == []

    def test_create_sample_config_basic_yaml_roundtrip(self, tmp_path: Path) -> None:
        """Basic sample config written as YAML round-trips through load_workflow_config."""
        output_path = tmp_path / "sample.yaml"

        create_sample_config(output_path, sample_type="basic")

        assert output_path.exists()
        config = load_workflow_config(output_path)
        assert config.work_dir.name == "work"
        assert config.work_dir.parts[-3:] == ("output", "amalgkit", "work")
        assert config.threads == 8
        assert config.species_list == ["Apis_mellifera"]

    def test_create_sample_config_advanced_yaml_roundtrip(self, tmp_path: Path) -> None:
        """Advanced sample config carries both species and higher thread count."""
        output_path = tmp_path / "advanced.yaml"

        create_sample_config(output_path, sample_type="advanced")

        assert output_path.exists()
        config = load_workflow_config(output_path)
        assert config.work_dir.name == "work"
        assert config.work_dir.parts[-3:] == ("output", "amalgkit", "work")
        assert config.threads == 12
        assert config.species_list == ["Apis_mellifera", "Pogonomyrmex_barbatus"]

    def test_create_sample_config_json_suffix(self, tmp_path: Path) -> None:
        """A .json suffix path writes JSON that json.load can read back."""
        output_path = tmp_path / "sample.json"

        create_sample_config(output_path, sample_type="basic")

        assert output_path.exists()
        data = json.loads(output_path.read_text())
        assert data["work_dir"] == "output/amalgkit/work"
        assert data["threads"] == 8

    def test_apply_config_defaults_env_thread_override(self, monkeypatch) -> None:
        """AMALGKIT_PIPELINE_THREADS env var overrides the default thread count."""
        monkeypatch.setenv("AMALGKIT_PIPELINE_THREADS", "2")

        result = apply_config_defaults({})

        assert result["threads"] == 2

    def test_apply_config_defaults_invalid_env_value_ignored(self, monkeypatch) -> None:
        """A non-integer AMALGKIT_PIPELINE_THREADS leaves the default of 8 intact."""
        monkeypatch.setenv("AMALGKIT_PIPELINE_THREADS", "abc")
        logging.getLogger("metainformant.rna.engine.workflow_core").setLevel(logging.CRITICAL)

        result = apply_config_defaults({})

        assert result["threads"] == 8
