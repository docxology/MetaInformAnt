"""RNA-seq workflow configuration and data classes.

This module provides the core configuration classes, data structures,
and validation functions for RNA-seq analysis workflows.
"""

from __future__ import annotations

import os
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, Union

from metainformant.core.utils import logging
from metainformant.core import io
from metainformant.core.utils.config import load_mapping_from_file

logger = logging.get_logger(__name__)


@dataclass
class WorkflowStepResult:
    """Result of executing a single workflow step."""

    step_name: str
    return_code: int
    success: bool
    error_message: Optional[str] = None
    command: Optional[str] = None


class WorkflowExecutionResult(list):
    """Result of executing the entire amalgkit workflow.

    Inherits from list to satisfy legacy tests checking isinstance(result, list).
    """

    def __init__(
        self,
        steps_executed: List[WorkflowStepResult],
        success: bool,
        total_steps: int,
        successful_steps: int,
        failed_steps: int,
    ):
        # Initialize list with return codes
        super().__init__([s.return_code for s in steps_executed])
        self.steps_executed = steps_executed
        self.success = success
        self.total_steps = total_steps
        self.successful_steps = successful_steps
        self.failed_steps = failed_steps

    @property
    def return_codes(self) -> List[int]:
        """Return codes for backward compatibility."""
        return [s.return_code for s in self.steps_executed]

    def __len__(self) -> int:
        return len(self.steps_executed)

    def __getitem__(self, index: Any) -> Any:
        if isinstance(index, (int, slice)):
            # If we want the step result object
            if isinstance(index, int):
                return self.steps_executed[index]
            return self.steps_executed[index]
        return super().__getitem__(index)

    def get(self, step_name: str, default: Any = None) -> Any:
        """Get step result by name for backward compatibility."""
        for step in self.steps_executed:
            if step.step_name == step_name:
                return step.return_code
        return default


class AmalgkitWorkflowConfig:
    """Configuration class for amalgkit workflows."""

    def __init__(
        self,
        work_dir: Union[str, Path] = ".",
        threads: int = 8,
        species_list: Optional[List[str]] = None,
        search_string: Optional[str] = None,
        max_samples: Optional[int] = None,
        genome: Optional[Dict[str, Any]] = None,
        taxon_id: Optional[str] = None,
        auto_install_amalgkit: bool = True,
        log_dir: Optional[str] = None,
        **kwargs,
    ):
        """Initialize workflow configuration.

        Args:
            work_dir: Working directory for the workflow
            threads: Number of threads to use
            species_list: List of species to process
            search_string: Search string for sample discovery
            max_samples: Maximum number of samples per species
            genome: Genome configuration dictionary
            **kwargs: Additional configuration options
        """
        self.work_dir = Path(work_dir)
        self.threads = threads
        self.species_list = species_list or []
        self.search_string = search_string
        self.max_samples = max_samples
        self.genome = genome or {}
        self.taxon_id = taxon_id
        self.auto_install_amalgkit = auto_install_amalgkit
        self.log_dir = Path(log_dir) if log_dir else self.work_dir / "logs"

        # New attributes for workflow control
        self.per_step = kwargs.get("per_step", {}) or kwargs.get("steps", {})
        self.filters = kwargs.get("filters", {})

        self.extra_config = kwargs

    @property
    def manifest_path(self) -> Path:
        """Path to the workflow manifest file."""
        return self.work_dir / "amalgkit.manifest.jsonl"

    @property
    def log_file(self) -> Path:
        """Path to the main workflow log file."""
        return self.log_dir / "workflow.log"

    def to_dict(self) -> Dict[str, Any]:
        """Convert configuration to dictionary."""
        return {
            "work_dir": str(self.work_dir),
            "threads": self.threads,
            "species_list": self.species_list,
            "search_string": self.search_string,
            "max_samples": self.max_samples,
            "genome": self.genome,
            "taxon_id": self.taxon_id,
            "auto_install_amalgkit": self.auto_install_amalgkit,
            "log_dir": str(self.log_dir),
            "per_step": self.per_step,
            "filters": self.filters,
            **self.extra_config,
        }

    @classmethod
    def from_dict(cls, config_dict: Dict[str, Any]) -> AmalgkitWorkflowConfig:
        """Create configuration from dictionary."""
        return cls(**config_dict)


def load_workflow_config(config_file: Union[str, Path]) -> AmalgkitWorkflowConfig:
    """Load workflow configuration from YAML/TOML/JSON file.

    Args:
        config_file: Path to configuration file

    Returns:
        AmalgkitWorkflowConfig instance

    Raises:
        FileNotFoundError: If config file doesn't exist
        ValueError: If config format is invalid
    """
    config_path = Path(config_file)

    if not config_path.exists():
        raise FileNotFoundError(f"Configuration file not found: {config_path}")

    try:
        # Load config using core.config which supports YAML/TOML/JSON
        config_dict = load_mapping_from_file(config_path)

        # Apply defaults
        config_dict = apply_config_defaults(config_dict)

        return AmalgkitWorkflowConfig.from_dict(config_dict)

    except (OSError, IOError) as e:
        raise ValueError(f"Error reading configuration file {config_path}: {e}") from e
    except (KeyError, TypeError) as e:
        raise ValueError(f"Invalid configuration format in {config_path}: {e}") from e


def apply_config_defaults(config: Dict[str, Any]) -> Dict[str, Any]:
    """Apply default values to configuration dictionary.

    Args:
        config: Configuration dictionary

    Returns:
        Configuration with defaults applied
    """
    defaults = {
        "threads": 8,
        "species_list": [],
        "genome": {},
        "resolve_names": "yes",  # v0.12.20+ feature
        "mark_missing_rank": "species",  # v0.12.20+ feature
    }

    # Merge defaults with provided config
    result = defaults.copy()
    result.update(config)

    # Environment variable overrides
    env_threads = os.environ.get("AK_THREADS")
    if env_threads:
        try:
            result["threads"] = int(env_threads)
        except ValueError:
            pass

    return result


def validate_workflow_config(config: AmalgkitWorkflowConfig) -> Tuple[bool, List[str]]:
    """Validate workflow configuration.

    Args:
        config: Workflow configuration to validate

    Returns:
        Tuple of (is_valid, error_messages)
    """
    errors = []

    # Check required fields
    if not config.work_dir:
        errors.append("work_dir is required")

    if not config.species_list:
        errors.append("species_list cannot be empty")

    if config.threads < 1:
        errors.append("threads must be >= 1")

    # Check work directory
    if config.work_dir.exists() and not config.work_dir.is_dir():
        errors.append(f"work_dir exists but is not a directory: {config.work_dir}")

    # Validate species names (basic check)
    for species in config.species_list:
        if not isinstance(species, str) or not species.strip():
            errors.append(f"Invalid species name: {species}")

    return len(errors) == 0, errors


def validate_workflow_outputs(config: AmalgkitWorkflowConfig) -> Tuple[bool, List[str]]:
    """Validate workflow outputs.

    Args:
        config: Workflow configuration

    Returns:
        Tuple of (is_valid, error_messages)
    """
    errors = []

    # Check for expected output files
    expected_outputs = [
        config.work_dir / "metadata.tsv",
        config.work_dir / "expression_matrix.tsv",
        config.work_dir / "sanity_check.txt",
    ]

    for output_file in expected_outputs:
        if not output_file.exists():
            errors.append(f"Expected output file missing: {output_file}")

    return len(errors) == 0, errors


def create_sample_config(output_path: Union[str, Path], sample_type: str = "basic") -> None:
    """Create a sample workflow configuration file.

    Args:
        output_path: Path to write sample configuration
        sample_type: Type of sample configuration ('basic', 'advanced')
    """
    if sample_type == "basic":
        sample_config = {
            "work_dir": "output/amalgkit/work",
            "threads": 8,
            "species_list": ["Apis_mellifera"],
            "genome": {"accession": "GCF_003254395.2", "assembly_name": "Amel_HAv3.1", "annotation_release": 104},
            "resolve_names": "yes",
            "mark_missing_rank": "species",
        }
    else:
        sample_config = {
            "work_dir": "output/amalgkit/work",
            "threads": 12,
            "species_list": ["Apis_mellifera", "Pogonomyrmex_barbatus"],
            "search_string": "RNA-seq",
            "max_samples": 100,
            "genome": {"accession": "GCF_003254395.2", "assembly_name": "Amel_HAv3.1", "annotation_release": 104},
            "resolve_names": "yes",
            "mark_missing_rank": "species",
        }

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    if output_path.suffix.lower() in [".yaml", ".yml"]:
        io.dump_yaml(sample_config, str(output_path))
    elif output_path.suffix.lower() == ".json":
        io.dump_json(sample_config, output_path)
    else:
        # Default to YAML
        io.dump_yaml(sample_config, str(output_path))

    logger.info(f"Sample configuration created: {output_path}")
