"""RNA-seq workflow orchestration and configuration management.

This module provides high-level functions for managing RNA-seq analysis workflows,
including configuration loading, workflow planning, and execution orchestration.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, Union
from dataclasses import dataclass

from metainformant.core import io, logging
from metainformant.core.config import load_mapping_from_file
from metainformant.rna.metadata_filter import filter_selected_metadata

logger = logging.get_logger(__name__)


@dataclass
class WorkflowStepResult:
    """Result of executing a single workflow step."""
    step_name: str
    return_code: int
    success: bool
    error_message: Optional[str] = None
    command: Optional[str] = None


@dataclass
class WorkflowExecutionResult:
    """Result of executing a complete workflow."""
    steps_executed: List[WorkflowStepResult]
    success: bool
    total_steps: int
    successful_steps: int
    failed_steps: int

    @property
    def return_codes(self) -> List[int]:
        """Return list of return codes for backward compatibility."""
        return [step.return_code for step in self.steps_executed]


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
        # Load config using core.config which supports YAML/TOML/JSON
        config_dict = load_mapping_from_file(config_path)

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
        'config',
        'select',
        'getfastq',
        'integrate',  # After getfastq to integrate downloaded FASTQ files
        'quant',
        'merge',
        'cstmm',
        'curate',
        'csca',
        'sanity'
    ]

    # Get per-step overrides if provided (support both 'steps' and 'per_step' keys)
    per_step = config.extra_config.get('steps', {}) or config.extra_config.get('per_step', {})

    # Check if cstmm/csca required parameters are provided
    has_ortholog_params = (
        (per_step.get('cstmm', {}).get('orthogroup_table') or
         per_step.get('cstmm', {}).get('dir_busco')) or
        (per_step.get('csca', {}).get('orthogroup_table') or
         per_step.get('csca', {}).get('dir_busco')) or
        config.extra_config.get('orthogroup_table') or
        config.extra_config.get('dir_busco')
    )

    # Filter out cstmm and csca if required parameters not provided
    if not has_ortholog_params:
        workflow_steps = [step for step in workflow_steps if step not in ('cstmm', 'csca')]
        if 'cstmm' in per_step or 'csca' in per_step:
            logger.warning("Skipping cstmm and csca steps: required parameters (orthogroup_table or dir_busco) not provided")

    for step in workflow_steps:
        # Start with only step-specific params if provided, otherwise use common config
        if step in per_step:
            # Use step-specific parameters
            step_params = per_step[step].copy()
            # Ensure work_dir is set from common config if not in step params
            if 'out_dir' not in step_params and 'work_dir' in config.to_dict():
                step_params['out_dir'] = str(config.work_dir)
        else:
            # Use common config params
            step_params = {
                'out_dir': str(config.work_dir),
                'threads': config.threads,
            }
            # Add species for steps that need it
            if config.species_list:
                step_params['species'] = config.species_list

        # Auto-inject metadata paths for steps that need them
        # Use filtered metadata if it exists (created after select step)
        if step in ('getfastq', 'integrate', 'merge'):
            filtered_metadata = config.work_dir / 'metadata' / 'metadata_selected.tsv'
            if filtered_metadata.exists():
                step_params['metadata'] = str(filtered_metadata)
            elif 'metadata' not in step_params:
                step_params['metadata'] = str(config.work_dir / 'metadata' / 'metadata.tsv')
        elif step == 'curate' and 'metadata' not in step_params:
            # Curate uses metadata from merge directory
            step_params['metadata'] = str(config.work_dir / 'merge' / 'metadata' / 'metadata.tsv')

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

    # Remove invalid parameters per subcommand
    INVALID_PARAMS = {
        'quant': {'keep_fastq'},  # keep_fastq is not supported by amalgkit CLI
        # Add more as needed
    }

    if subcommand in INVALID_PARAMS:
        for param in INVALID_PARAMS[subcommand]:
            sanitized.pop(param, None)

    return sanitized


def execute_workflow(config: AmalgkitWorkflowConfig, *,
                    steps: Optional[List[str]] = None,
                    check: bool = False,
                    walk: bool = False,
                    progress: bool = True,
                    show_commands: bool = False) -> WorkflowExecutionResult:
    """Execute the complete amalgkit workflow.

    Args:
        config: Workflow configuration
        steps: Specific steps to run (if None, run all steps)
        check: Stop on first failure if True
        walk: Dry run mode
        progress: Show progress indicators
        show_commands: Print commands being executed

    Returns:
        WorkflowExecutionResult with detailed step results
    """
    from metainformant.rna.amalgkit import (
        AmalgkitParams, metadata, integrate, select,
        getfastq, quant, merge, cstmm, curate, csca, sanity
    )
    from metainformant.rna.amalgkit import config as amalgkit_config

    logger.info(f"Starting amalgkit workflow for species: {config.species_list}")
    logger.info(f"Working directory: {config.work_dir}")

    # Create working directory
    config.work_dir.mkdir(parents=True, exist_ok=True)

    # Plan workflow
    steps_planned = plan_workflow(config)

    # Filter steps if specific steps requested
    if steps:
        original_count = len(steps_planned)
        steps_planned = [(name, params) for name, params in steps_planned if name in steps]
        logger.info(f"Filtered workflow from {original_count} to {len(steps_planned)} steps: {[s[0] for s in steps_planned]}")

    if not steps_planned:
        logger.warning("No steps to execute after filtering")
        return WorkflowExecutionResult(
            steps_executed=[],
            success=True,
            total_steps=0,
            successful_steps=0,
            failed_steps=0
        )

    step_results = []
    step_functions = {
        'metadata': metadata,
        'integrate': integrate,
        'config': amalgkit_config,
        'select': select,
        'getfastq': getfastq,
        'quant': quant,
        'merge': merge,
        'cstmm': cstmm,
        'curate': curate,
        'csca': csca,
        'sanity': sanity
    }

    logger.info(f"Starting execution of {len(steps_planned)} steps")
    for i, (step_name, step_params) in enumerate(steps_planned, 1):
        if progress:
            logger.info(f"Step {i}/{len(steps_planned)}: {step_name}")

        if walk:
            logger.info(f"Would execute: {step_name}")
            step_results.append(WorkflowStepResult(
                step_name=step_name,
                return_code=0,
                success=True,
                command="(dry run)"
            ))
            continue

        try:
            # Get step function
            step_func = step_functions.get(step_name)
            if not step_func:
                error_msg = f"Unknown step: {step_name}"
                logger.error(error_msg)
                step_results.append(WorkflowStepResult(
                    step_name=step_name,
                    return_code=1,
                    success=False,
                    error_message=error_msg
                ))
                if check:
                    break
                continue

            # Execute step (pass params dict directly, not AmalgkitParams)
            command_str = None
            if show_commands:
                from metainformant.rna.amalgkit import build_amalgkit_command
                sanitized_params = sanitize_params_for_cli(step_name, step_params)
                command = build_amalgkit_command(step_name, sanitized_params)
                command_str = ' '.join(command)
                logger.info(f"Command: {command_str}")

            # Validate filtered metadata exists for steps that need it
            if step_name in ('getfastq', 'integrate', 'merge'):
                filtered_metadata = config.work_dir / 'metadata' / 'metadata_selected.tsv'
                if not filtered_metadata.exists():
                    error_msg = f"Filtered metadata not found: {filtered_metadata}. Run 'select' step first."
                    logger.error(error_msg)
                    step_results.append(WorkflowStepResult(
                        step_name=step_name,
                        return_code=1,
                        success=False,
                        error_message=error_msg
                    ))
                    if check:
                        break
                    continue

            # Sanitize parameters before passing to step function
            sanitized_params = sanitize_params_for_cli(step_name, step_params)
            result = step_func(sanitized_params)

            # After select step, create filtered metadata for downstream steps
            if step_name == 'select' and result.returncode == 0:
                try:
                    logger.info("Creating filtered metadata for selected samples")
                    filter_selected_metadata(
                        config.work_dir / 'metadata' / 'metadata.tsv',
                        config.work_dir / 'metadata' / 'metadata_selected.tsv'
                    )
                except Exception as e:
                    logger.error(f"Failed to create filtered metadata: {e}")
                    if check:
                        step_results.append(WorkflowStepResult(
                            step_name="metadata_filtering",
                            return_code=1,
                            success=False,
                            error_message=str(e)
                        ))
                        break

            # Create step result
            error_message = None
            if result.returncode != 0:
                error_message = f"Step failed with return code {result.returncode}"
                if result.stderr:
                    error_message += f": {result.stderr}"
                logger.error(f"Step {step_name} failed with return code {result.returncode}")
                if result.stderr:
                    logger.error(f"Error output: {result.stderr}")
            else:
                logger.info(f"Step {step_name} completed successfully")

            step_results.append(WorkflowStepResult(
                step_name=step_name,
                return_code=result.returncode,
                success=result.returncode == 0,
                error_message=error_message,
                command=command_str
            ))

            if result.returncode != 0 and check:
                break

        except Exception as e:
            error_msg = f"Exception during execution: {e}"
            logger.error(f"Error executing step {step_name}: {e}")
            step_results.append(WorkflowStepResult(
                step_name=step_name,
                return_code=1,
                success=False,
                error_message=error_msg
            ))
            if check:
                break

    successful_steps = sum(1 for sr in step_results if sr.success)
    total_steps = len(step_results)
    failed_steps = total_steps - successful_steps

    if progress:
        logger.info(f"Workflow completed: {successful_steps}/{total_steps} steps successful")

    return WorkflowExecutionResult(
        steps_executed=step_results,
        success=failed_steps == 0,
        total_steps=total_steps,
        successful_steps=successful_steps,
        failed_steps=failed_steps
    )


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







