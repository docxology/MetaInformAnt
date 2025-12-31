"""RNA-seq workflow orchestration and configuration management.

This module provides high-level functions for managing RNA-seq analysis workflows,
including configuration loading, workflow planning, and execution orchestration.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, Union

from metainformant.core import io, logging

logger = logging.get_logger(__name__)


class AmalgkitWorkflowConfig:
    """Configuration class for amalgkit workflows."""

    def __init__(self,
                 work_dir: Path,
                 threads: int = 8,
                 species_list: Optional[List[str]] = None,
                 search_string: Optional[str] = None,
                 max_samples: Optional[int] = None,
                 genome: Optional[Dict[str, Any]] = None,
                 **kwargs):
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
        self.extra_config = kwargs

    def to_dict(self) -> Dict[str, Any]:
        """Convert configuration to dictionary."""
        return {
            'work_dir': str(self.work_dir),
            'threads': self.threads,
            'species_list': self.species_list,
            'search_string': self.search_string,
            'max_samples': self.max_samples,
            'genome': self.genome,
            **self.extra_config
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
        if config_path.suffix.lower() in ['.yaml', '.yml']:
            config_dict = io.load_yaml(str(config_path))
        elif config_path.suffix.lower() == '.toml':
            config_dict = io.load_toml(str(config_path))
        elif config_path.suffix.lower() == '.json':
            config_dict = io.load_json(str(config_path))
        else:
            raise ValueError(f"Unsupported config format: {config_path.suffix}")

        # Apply defaults
        config_dict = apply_config_defaults(config_dict)

        return AmalgkitWorkflowConfig.from_dict(config_dict)

    except Exception as e:
        raise ValueError(f"Error loading configuration from {config_path}: {e}")


def apply_config_defaults(config: Dict[str, Any]) -> Dict[str, Any]:
    """Apply default values to configuration dictionary.

    Args:
        config: Configuration dictionary

    Returns:
        Configuration with defaults applied
    """
    defaults = {
        'threads': 8,
        'species_list': [],
        'genome': {},
        'resolve_names': 'yes',  # v0.12.20+ feature
        'mark_missing_rank': 'species',  # v0.12.20+ feature
    }

    # Merge defaults with provided config
    result = defaults.copy()
    result.update(config)

    return result


def apply_step_defaults(config: AmalgkitWorkflowConfig) -> AmalgkitWorkflowConfig:
    """Apply default values to workflow step configurations.

    Args:
        config: Workflow configuration

    Returns:
        Configuration with step defaults applied
    """
    # This is a placeholder - in practice, step-specific defaults would be applied
    return config


def plan_workflow(config: AmalgkitWorkflowConfig) -> List[Tuple[str, Any]]:
    """Plan the workflow steps and their parameters.

    Args:
        config: Workflow configuration

    Returns:
        List of (step_name, parameters) tuples
    """
    steps = []

    # Define the standard amalgkit workflow steps
    workflow_steps = [
        'metadata',
        'integrate',
        'config',
        'select',
        'getfastq',
        'quant',
        'merge',
        'cstmm',
        'curate',
        'csca',
        'sanity'
    ]

    for step in workflow_steps:
        step_params = config.to_dict()
        # Step-specific parameter overrides would go here
        steps.append((step, step_params))

    logger.info(f"Planned workflow with {len(steps)} steps")
    return steps


def plan_workflow_with_params(config: AmalgkitWorkflowConfig, **param_overrides: Any) -> List[Tuple[str, Any]]:
    """Plan workflow with parameter overrides.

    Args:
        config: Base workflow configuration
        **param_overrides: Parameter overrides

    Returns:
        List of (step_name, parameters) tuples
    """
    # Apply overrides to config
    override_config = AmalgkitWorkflowConfig(
        **{**config.to_dict(), **param_overrides}
    )

    return plan_workflow(override_config)


def sanitize_params_for_cli(subcommand: str, params: Dict[str, Any]) -> Dict[str, Any]:
    """Sanitize parameters for CLI usage.

    Args:
        subcommand: Amalgkit subcommand
        params: Parameter dictionary

    Returns:
        Sanitized parameter dictionary
    """
    # Remove or transform parameters that shouldn't go to CLI
    sanitized = params.copy()

    # Remove internal Python objects
    to_remove = []
    for key, value in sanitized.items():
        if isinstance(value, (Path, type)):
            to_remove.append(key)

    for key in to_remove:
        del sanitized[key]

    # Convert Path objects to strings
    for key, value in sanitized.items():
        if isinstance(value, Path):
            sanitized[key] = str(value)

    return sanitized


def execute_workflow(config: AmalgkitWorkflowConfig, *,
                    check: bool = False,
                    walk: bool = False,
                    progress: bool = True,
                    show_commands: bool = False) -> List[int]:
    """Execute the complete amalgkit workflow.

    Args:
        config: Workflow configuration
        check: Stop on first failure if True
        walk: Dry run mode
        progress: Show progress indicators
        show_commands: Print commands being executed

    Returns:
        List of return codes from each step
    """
    from metainformant.rna.amalgkit import (
        AmalgkitParams, metadata, integrate, config, select,
        getfastq, quant, merge, cstmm, curate, csca, sanity
    )

    logger.info(f"Starting amalgkit workflow for species: {config.species_list}")
    logger.info(f"Working directory: {config.work_dir}")

    # Create working directory
    config.work_dir.mkdir(parents=True, exist_ok=True)

    # Plan workflow
    steps = plan_workflow(config)

    return_codes = []
    step_functions = {
        'metadata': metadata,
        'integrate': integrate,
        'config': config,
        'select': select,
        'getfastq': getfastq,
        'quant': quant,
        'merge': merge,
        'cstmm': cstmm,
        'curate': curate,
        'csca': csca,
        'sanity': sanity
    }

    for step_name, step_params in steps:
        if progress:
            logger.info(f"Executing step: {step_name}")

        if walk:
            logger.info(f"Would execute: {step_name}")
            return_codes.append(0)
            continue

        try:
            # Create step parameters
            params = AmalgkitParams(**step_params)

            # Get step function
            step_func = step_functions.get(step_name)
            if not step_func:
                logger.error(f"Unknown step: {step_name}")
                return_codes.append(1)
                if check:
                    break
                continue

            # Execute step
            if show_commands:
                from metainformant.rna.amalgkit import build_amalgkit_command
                command = build_amalgkit_command(step_name, params)
                logger.info(f"Command: {' '.join(command)}")

            result = step_func(params)

            return_codes.append(result.returncode)

            if result.returncode != 0:
                logger.error(f"Step {step_name} failed with return code {result.returncode}")
                if result.stderr:
                    logger.error(f"Error output: {result.stderr}")
                if check:
                    break
            else:
                logger.info(f"Step {step_name} completed successfully")

        except Exception as e:
            logger.error(f"Error executing step {step_name}: {e}")
            return_codes.append(1)
            if check:
                break

    success_count = sum(1 for rc in return_codes if rc == 0)
    total_count = len(return_codes)

    if progress:
        logger.info(f"Workflow completed: {success_count}/{total_count} steps successful")

    return return_codes


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
        config.work_dir / "sanity_check.txt"
    ]

    for output_file in expected_outputs:
        if not output_file.exists():
            errors.append(f"Expected output file missing: {output_file}")

    return len(errors) == 0, errors


def create_sample_config(output_path: Union[str, Path],
                        sample_type: str = "basic") -> None:
    """Create a sample workflow configuration file.

    Args:
        output_path: Path to write sample configuration
        sample_type: Type of sample configuration ('basic', 'advanced')
    """
    if sample_type == "basic":
        sample_config = {
            'work_dir': 'output/amalgkit/work',
            'threads': 8,
            'species_list': ['Apis_mellifera'],
            'genome': {
                'accession': 'GCF_003254395.2',
                'assembly_name': 'Amel_HAv3.1',
                'annotation_release': 104
            },
            'resolve_names': 'yes',
            'mark_missing_rank': 'species'
        }
    else:
        sample_config = {
            'work_dir': 'output/amalgkit/work',
            'threads': 12,
            'species_list': ['Apis_mellifera', 'Pogonomyrmex_barbatus'],
            'search_string': 'RNA-seq',
            'max_samples': 100,
            'genome': {
                'accession': 'GCF_003254395.2',
                'assembly_name': 'Amel_HAv3.1',
                'annotation_release': 104
            },
            'resolve_names': 'yes',
            'mark_missing_rank': 'species'
        }

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    if output_path.suffix.lower() in ['.yaml', '.yml']:
        io.dump_yaml(sample_config, str(output_path))
    elif output_path.suffix.lower() == '.json':
        io.dump_json(sample_config, output_path)
    else:
        # Default to YAML
        io.dump_yaml(sample_config, str(output_path))

    logger.info(f"Sample configuration created: {output_path}")


# Helper functions for backwards compatibility
def run_config_based_workflow(config_path: Union[str, Path], **kwargs: Any) -> Dict[str, Any]:
    """Run workflow from configuration file.

    Args:
        config_path: Path to configuration file
        **kwargs: Additional workflow options

    Returns:
        Workflow results dictionary
    """
    config = load_workflow_config(config_path)

    # Apply any overrides
    for key, value in kwargs.items():
        if hasattr(config, key):
            setattr(config, key, value)

    return_codes = execute_workflow(config, **kwargs)

    return {
        'config': config.to_dict(),
        'return_codes': return_codes,
        'success': all(rc == 0 for rc in return_codes)
    }


