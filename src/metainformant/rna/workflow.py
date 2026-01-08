"""RNA-seq workflow orchestration and configuration management.

This module provides high-level functions for managing RNA-seq analysis workflows,
including configuration loading, workflow planning, and execution orchestration.
"""

from __future__ import annotations

import json
import shutil
import os
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
    
    This function handles automatic path resolution for workflow steps:
    
    - **getfastq**: Creates FASTQ files in {out_dir}/getfastq/{sample_id}/ (automatically creates getfastq subdirectory)
    - **integrate**: Automatically adjusts fastq_dir to point to {fastq_dir}/getfastq/ if it exists
    - **quant**: Should use out_dir = work_dir so it can find getfastq output in {out_dir}/getfastq/
    - **merge**: Looks for abundance files in {out_dir}/quant/{sample_id}/{sample_id}_abundance.tsv
    
    See docs/rna/amalgkit/PATH_RESOLUTION.md for complete path resolution documentation.

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
            
            # Warn if redo: yes is set for getfastq (usually unnecessary)
            if step == 'getfastq':
                redo_value = step_params.get('redo', 'no')
                # Handle both string and boolean values
                if isinstance(redo_value, bool):
                    redo_is_yes = redo_value
                else:
                    redo_is_yes = str(redo_value).lower() in ('yes', 'true', '1')
                
                if redo_is_yes:
                    logger.warning(
                        "getfastq step has redo: yes - this will force re-download of all files, "
                        "even if they already exist. Consider using redo: no to skip already-downloaded files. "
                        "Use redo: yes only if files are corrupted or you need to refresh data."
                    )
        else:
            # Use common config params
            step_params = {
                'out_dir': str(config.work_dir),
                'threads': config.threads,
            }
            # Add species for steps that need it
            if config.species_list:
                step_params['species'] = config.species_list

        # Parallel download configuration for getfastq
        if step == 'getfastq' and 'jobs' not in step_params:
            # Default heuristic: Use 1 job per 2 threads, max 4 jobs
            # Adjust threads per job to keep total thread usage roughly within limit
            total_threads = config.threads
            if total_threads >= 2:
                jobs = min(4, total_threads)
                threads_per_job = max(1, total_threads // jobs)
                
                step_params['jobs'] = jobs
                step_params['threads'] = threads_per_job
                logger.debug(f"Auto-configured parallel getfastq: {jobs} jobs with {threads_per_job} threads each (total threads: {total_threads})")

        # Auto-inject metadata path for select step
        if step == 'select' and 'metadata' not in step_params:
            # Select step needs metadata.tsv as input
            metadata_file = config.work_dir / 'metadata' / 'metadata.tsv'
            if metadata_file.exists():
                step_params['metadata'] = str(metadata_file.resolve())
        
        # Auto-inject metadata paths for steps that need them
        # Use filtered metadata if it exists (created after select step)
        if step in ('getfastq', 'integrate', 'merge'):
            filtered_metadata = config.work_dir / 'metadata' / 'metadata_selected.tsv'
            if filtered_metadata.exists():
                step_params['metadata'] = str(filtered_metadata)
            elif 'metadata' not in step_params:
                step_params['metadata'] = str(config.work_dir / 'metadata' / 'metadata.tsv')
        elif step == 'quant':
            # Quant uses metadata from integrate step (amalgkit creates metadata.tsv after integrate)
            # Prefer integrated metadata if it exists, otherwise use selected metadata
            integrated_metadata = config.work_dir / 'metadata' / 'metadata.tsv'
            filtered_metadata = config.work_dir / 'metadata' / 'metadata_selected.tsv'
            if integrated_metadata.exists():
                step_params['metadata'] = str(integrated_metadata)
            elif filtered_metadata.exists():
                step_params['metadata'] = str(filtered_metadata)
            elif 'metadata' not in step_params:
                step_params['metadata'] = str(config.work_dir / 'metadata' / 'metadata.tsv')
        elif step == 'curate' and 'metadata' not in step_params:
            # Curate uses metadata from merge directory
            step_params['metadata'] = str(config.work_dir / 'merge' / 'metadata' / 'metadata.tsv')

        # PATH RESOLUTION: Auto-adjust integrate fastq_dir to include getfastq subdirectory
        # 
        # Context: amalgkit getfastq automatically creates a 'getfastq/' subdirectory within
        # the specified out_dir. For example, if getfastq uses out_dir: output/fastq,
        # then FASTQ files are actually in output/fastq/getfastq/{sample_id}/.
        #
        # The integrate step needs to point to this getfastq subdirectory to find the files.
        # This code automatically adjusts the path if the getfastq subdirectory exists.
        #
        # See docs/rna/amalgkit/PATH_RESOLUTION.md for complete path resolution guide.
        if step == 'integrate':
            if 'fastq_dir' in step_params:
                fastq_dir = Path(step_params['fastq_dir'])
                # Check if getfastq subdirectory exists (amalgkit getfastq creates this)
                getfastq_subdir = fastq_dir / "getfastq"
                if getfastq_subdir.exists():
                    step_params['fastq_dir'] = str(getfastq_subdir)
                    logger.debug(f"Adjusted integrate fastq_dir to include getfastq subdirectory: {getfastq_subdir}")
                # Also check if fastq_dir itself is the getfastq directory
                elif fastq_dir.name == "getfastq":
                    # Already pointing to getfastq, no adjustment needed
                    logger.debug(f"Integrate fastq_dir already points to getfastq directory: {fastq_dir}")
            else:
                # If fastq_dir not specified, try to infer from getfastq step
                getfastq_out_dir = per_step.get('getfastq', {}).get('out_dir')
                if getfastq_out_dir:
                    fastq_dir = Path(getfastq_out_dir)
                    getfastq_subdir = fastq_dir / "getfastq"
                    if getfastq_subdir.exists():
                        step_params['fastq_dir'] = str(getfastq_subdir)
                        logger.debug(f"Inferred integrate fastq_dir from getfastq step: {getfastq_subdir}")

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


def _log_getfastq_summary(output_text: str, logger: Any) -> None:
    """Parse amalgkit getfastq output and log summary of skipped vs downloaded files.
    
    Args:
        output_text: Combined stdout/stderr output from amalgkit getfastq command
        logger: Logger instance for logging summary
    """
    if not output_text:
        return
    
    # Count patterns in amalgkit output
    # Pattern 1: Files that were already downloaded (skipped)
    skipped_count = output_text.count("Previously-downloaded sra file was detected")
    # Pattern 2: Files that needed to be downloaded
    downloaded_count = output_text.count("Previously-downloaded sra file was not detected")
    # Pattern 3: Total files processed
    processed_count = output_text.count("Processing SRA ID:")
    
    # Try to extract total samples from metadata if available
    total_samples = None
    if "Number of SRAs to be processed:" in output_text:
        try:
            for line in output_text.split('\n'):
                if "Number of SRAs to be processed:" in line:
                    total_samples = int(line.split(":")[-1].strip())
                    break
        except (ValueError, IndexError):
            pass
    
    # Log summary only if we found relevant information
    if skipped_count > 0 or downloaded_count > 0 or processed_count > 0:
        summary_parts = []
        if skipped_count > 0:
            summary_parts.append(f"{skipped_count} file(s) skipped (already exists)")
        if downloaded_count > 0:
            summary_parts.append(f"{downloaded_count} file(s) downloaded")
        if processed_count > 0 and (skipped_count == 0 and downloaded_count == 0):
            # Only show processed count if we couldn't determine skipped/downloaded
            summary_parts.append(f"{processed_count} file(s) processed")
        
        if summary_parts:
            summary = ", ".join(summary_parts)
            if total_samples:
                logger.info(f"Step getfastq summary: {summary} (total: {total_samples} samples)")
            else:
                logger.info(f"Step getfastq summary: {summary}")


def _cleanup_incorrectly_placed_sra_files(getfastq_dir: Path) -> None:
    """Find and move SRA files from wrong locations to correct location.
    
    Args:
        getfastq_dir: Directory where amalgkit expects SRA files
    """
    # Common locations where prefetch might download SRA files
    default_locations = [
        Path.home() / "ncbi" / "public" / "sra",
        Path("/tmp") / "ncbi" / "public" / "sra",
    ]
    
    moved_count = 0
    for default_loc in default_locations:
        if not default_loc.exists():
            continue
        
        # Find SRA files in default location
        for sra_file in default_loc.rglob("*.sra"):
            try:
                # Extract sample ID from filename (e.g., SRR34065661.sra -> SRR34065661)
                sample_id = sra_file.stem
                target_dir = getfastq_dir / sample_id
                target_dir.mkdir(parents=True, exist_ok=True)
                target_file = target_dir / sra_file.name
                
                # Only move if target doesn't exist
                if not target_file.exists():
                    logger.info(f"Moving SRA file from wrong location: {sra_file} -> {target_file}")
                    shutil.move(str(sra_file), str(target_file))
                    moved_count += 1
                else:
                    # Target exists, remove duplicate
                    logger.debug(f"Target SRA file already exists, removing duplicate: {sra_file}")
                    sra_file.unlink()
            except Exception as e:
                logger.warning(f"Could not move SRA file {sra_file}: {e}")
    
    if moved_count > 0:
        logger.info(f"Moved {moved_count} SRA files from wrong locations to correct location")


def _cleanup_temp_files(tmp_dir: Path, max_size_gb: float = 50.0) -> None:
    """Clean up temporary files if directory gets too large.
    
    Args:
        tmp_dir: Temporary directory to clean
        max_size_gb: Maximum size in GB before cleanup
    """
    if not tmp_dir.exists():
        return
    
    try:
        # Calculate directory size
        total_size = sum(f.stat().st_size for f in tmp_dir.rglob('*') if f.is_file())
        size_gb = total_size / (1024 ** 3)
        
        if size_gb > max_size_gb:
            logger.warning(f"Temporary directory {tmp_dir} is {size_gb:.2f} GB (max: {max_size_gb} GB), cleaning up...")
            # Remove all files in temp directory
            for item in tmp_dir.iterdir():
                if item.is_file():
                    item.unlink()
                elif item.is_dir():
                    shutil.rmtree(item)
            logger.info(f"Cleaned up temporary directory {tmp_dir}")
    except Exception as e:
        logger.warning(f"Could not clean up temporary directory {tmp_dir}: {e}")


def _check_disk_space(path: Path, min_free_gb: float = 10.0) -> bool:
    """Check if there's sufficient disk space.
    
    Args:
        path: Path to check disk space for
        min_free_gb: Minimum free space required in GB
        
    Returns:
        True if sufficient space, False otherwise
    """
    try:
        stat = os.statvfs(path)
        free_gb = (stat.f_bavail * stat.f_frsize) / (1024 ** 3)
        
        if free_gb < min_free_gb:
            logger.warning(f"Low disk space: {free_gb:.2f} GB free (minimum: {min_free_gb} GB) at {path}")
            return False
        else:
            logger.debug(f"Disk space check: {free_gb:.2f} GB free at {path}")
            return True
    except Exception as e:
        logger.warning(f"Could not check disk space: {e}")
        return True  # Assume OK if we can't check


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

    # Convert Path objects to absolute path strings
    for key, value in sanitized.items():
        if isinstance(value, Path):
            # Use absolute path to ensure amalgkit can find files regardless of working directory
            sanitized[key] = str(value.resolve())

    # Remove invalid parameters per subcommand
    INVALID_PARAMS = {
        'quant': {'keep_fastq'},  # keep_fastq is not supported by amalgkit CLI
        # Add more as needed
    }

    if subcommand in INVALID_PARAMS:
        for param in INVALID_PARAMS[subcommand]:
            sanitized.pop(param, None)

    return sanitized


def _is_step_completed(step_name: str, step_params: dict, config: AmalgkitWorkflowConfig) -> tuple[bool, Optional[str]]:
    """Check if a workflow step has already completed.
    
    Args:
        step_name: Name of the workflow step
        step_params: Parameters for the step
        config: Workflow configuration
        
    Returns:
        Tuple of (is_completed: bool, completion_indicator: Optional[str])
        completion_indicator is a file path or description indicating why step is considered complete
    """
    work_dir = config.work_dir
    steps_config = config.extra_config.get('steps', {})
    
    if step_name == 'metadata':
        metadata_file = work_dir / 'metadata' / 'metadata.tsv'
        if metadata_file.exists():
            return True, str(metadata_file)
        return False, None
    
    elif step_name == 'config':
        config_dir = work_dir / 'config_base'
        if config_dir.exists():
            config_files = list(config_dir.glob('*.config'))
            if config_files:
                return True, f"{len(config_files)} config files in {config_dir}"
        return False, None
    
    elif step_name == 'select':
        selected_metadata = work_dir / 'metadata' / 'metadata_selected.tsv'
        if selected_metadata.exists():
            return True, str(selected_metadata)
        return False, None
    
    elif step_name == 'getfastq':
        # Check for FASTQ files in getfastq output directory
        fastq_dir = Path(step_params.get('out_dir', work_dir / 'fastq'))
        if fastq_dir.name != 'getfastq':
            fastq_dir = fastq_dir / 'getfastq'
        
        if fastq_dir.exists():
            fastq_files = list(fastq_dir.glob('**/*.fastq*'))
            if fastq_files:
                return True, f"{len(fastq_files)} FASTQ files in {fastq_dir}"
        return False, None
    
    elif step_name == 'integrate':
        # Check for integrated metadata
        integrated_meta = work_dir / 'integration' / 'integrated_metadata.json'
        if integrated_meta.exists():
            return True, str(integrated_meta)
        
        # Also check if metadata.tsv was updated (integrate updates it)
        metadata_tsv = work_dir / 'metadata' / 'metadata.tsv'
        if metadata_tsv.exists():
            # Check if there's an integration directory or if metadata has been processed
            integration_dir = work_dir / 'integration'
            if integration_dir.exists() and any(integration_dir.iterdir()):
                return True, f"Integration directory exists: {integration_dir}"
        return False, None
    
    elif step_name == 'quant':
        # Check for quantification output files
        quant_dir = Path(step_params.get('out_dir', work_dir))
        quant_files = list(quant_dir.glob('quant/**/abundance.tsv'))
        if quant_files:
            return True, f"{len(quant_files)} abundance files in {quant_dir / 'quant'}"
        return False, None
    
    elif step_name == 'merge':
        # Check for merged abundance file
        merge_out = step_params.get('out')
        if merge_out:
            merge_path = Path(merge_out)
            if merge_path.exists():
                return True, str(merge_path)
        
        merge_dir = Path(step_params.get('out_dir', work_dir / 'merged'))
        merged_file = merge_dir / 'merged_abundance.tsv'
        if merged_file.exists():
            return True, str(merged_file)
        return False, None
    
    elif step_name == 'curate':
        curate_dir = Path(step_params.get('out_dir', work_dir / 'curate'))
        if curate_dir.exists():
            # Check for expected curate output files
            expected_files = ['curated_abundance.tsv', 'curated_metadata.tsv']
            found_files = [f for f in expected_files if (curate_dir / f).exists()]
            if found_files:
                return True, f"Curate outputs found: {', '.join(found_files)} in {curate_dir}"
        return False, None
    
    elif step_name == 'sanity':
        sanity_file = work_dir / 'sanity_check.txt'
        if sanity_file.exists():
            return True, str(sanity_file)
        return False, None
    
    return False, None


def execute_workflow(config: AmalgkitWorkflowConfig, *,
                    steps: Optional[List[str]] = None,
                    check: bool = False,
                    walk: bool = False,
                    progress: bool = True,
                    show_commands: bool = False) -> WorkflowExecutionResult:
    """Execute the complete amalgkit workflow.
    
    This function automatically configures vdb-config repository path to use
    repository .tmp/vdb directory (on external drive with sufficient space) to
    prevent "disk-limit exceeded" errors during SRA extraction.
    
    The function checks if steps are already completed and skips them unless
    redo: yes is explicitly set in the configuration.

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

    # Configure vdb-config repository path for getfastq step (if getfastq is in workflow)
    # prefetch downloads SRA files to vdb-config repository root, so we need to set it to
    # the amalgkit fastq output directory where amalgkit expects the files
    if any(step == 'getfastq' for step, _ in steps_planned):
        try:
            import subprocess
            
            # Find getfastq step parameters to get out_dir
            getfastq_params = next((params for step, params in steps_planned if step == 'getfastq'), None)
            if getfastq_params:
                # Get fastq output directory from getfastq step params
                fastq_out_dir = getfastq_params.get('out_dir', str(config.work_dir / "fastq"))
                fastq_dir = Path(fastq_out_dir)
                # Amalgkit creates getfastq subdirectory, so prefetch should download to that
                getfastq_dir = fastq_dir / "getfastq" if fastq_dir.name != "getfastq" else fastq_dir
                getfastq_dir.mkdir(parents=True, exist_ok=True)
                
                # Check disk space before proceeding
                if not _check_disk_space(getfastq_dir, min_free_gb=20.0):
                    logger.warning("Low disk space detected, but continuing with getfastq step")
                
                # Clean up any incorrectly placed SRA files before configuring vdb-config
                _cleanup_incorrectly_placed_sra_files(getfastq_dir)
                
                # Clean up temp files if needed
                repo_root = Path(__file__).resolve().parent.parent.parent.parent
                tmp_dir = repo_root / ".tmp" / "fasterq-dump"
                _cleanup_temp_files(tmp_dir, max_size_gb=50.0)
                
                # Set vdb-config repository root to getfastq directory (where prefetch downloads SRA files)
                result = subprocess.run(
                    ['vdb-config', '-s', f'/repository/user/main/public/root={getfastq_dir}'],
                    capture_output=True,
                    text=True,
                    timeout=5
                )
                if result.returncode == 0:
                    logger.info(f"Configured vdb-config repository path for prefetch: {getfastq_dir}")
                else:
                    logger.warning(f"Could not set vdb-config repository path (may require interactive): {result.stderr[:100]}")
                    logger.info(f"Will rely on TMPDIR and VDB_CONFIG environment variables")
            else:
                logger.warning("Could not find getfastq step parameters, skipping vdb-config setup")
        except Exception as e:
            logger.warning(f"Could not configure vdb-config: {e}")
            logger.info("Will rely on TMPDIR and VDB_CONFIG environment variables")

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

    # Pre-execution status reporting: check which steps are already completed
    completed_steps = []
    steps_to_run = []
    steps_to_skip = []
    
    for step_name, step_params in steps_planned:
        is_completed, completion_indicator = _is_step_completed(step_name, step_params, config)
        
        # Check redo setting
        redo_value = step_params.get('redo', 'no')
        if isinstance(redo_value, bool):
            force_redo = redo_value
        else:
            force_redo = str(redo_value).lower() in ('yes', 'true', '1')
        
        if is_completed and not force_redo:
            completed_steps.append((step_name, completion_indicator))
            steps_to_skip.append(step_name)
        else:
            steps_to_run.append(step_name)
    
    # Log status summary
    logger.info(f"Workflow status summary:")
    logger.info(f"  Total steps planned: {len(steps_planned)}")
    if completed_steps:
        logger.info(f"  Steps already completed ({len(completed_steps)}): {', '.join([s[0] for s in completed_steps])}")
        for step_name, indicator in completed_steps:
            logger.info(f"    - {step_name}: {indicator}")
    if steps_to_run:
        logger.info(f"  Steps to run ({len(steps_to_run)}): {', '.join(steps_to_run)}")
    else:
        logger.info(f"  All steps already completed - nothing to run")

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

        # Check if step is already completed
        is_completed, completion_indicator = _is_step_completed(step_name, step_params, config)
        
        # Check redo setting
        redo_value = step_params.get('redo', 'no')
        if isinstance(redo_value, bool):
            force_redo = redo_value
        else:
            force_redo = str(redo_value).lower() in ('yes', 'true', '1')
        
        if is_completed and not force_redo:
            logger.info(f"Step {i}/{len(steps_planned)}: {step_name} - Already completed, skipping")
            logger.info(f"  Completion indicator: {completion_indicator}")
            
            # For config step, ensure symlink exists even if step was skipped
            if step_name == 'config':
                config_base_dir = config.work_dir / 'config_base'
                config_dir = config.work_dir / 'config'
                if config_base_dir.exists():
                    if config_dir.exists() or config_dir.is_symlink():
                        # Symlink or directory already exists - verify it points to config_base
                        if config_dir.is_symlink():
                            target = config_dir.readlink()
                            if target == config_base_dir or target.name == 'config_base':
                                logger.debug(f"Config symlink already exists and is correct: {config_dir}")
                            else:
                                logger.warning(f"Config symlink exists but points to wrong target: {target}")
                        else:
                            logger.debug(f"Config directory already exists: {config_dir}")
                    else:
                        try:
                            # Use absolute path for symlink to ensure it resolves correctly
                            config_dir.symlink_to(config_base_dir.resolve())
                            logger.info(f"Created symlink: {config_dir} -> {config_base_dir.resolve()} (for select step compatibility)")
                        except (OSError, FileExistsError) as e:
                            # Symlink already exists or file exists - check if it's correct
                            if config_dir.exists() or config_dir.is_symlink():
                                logger.debug(f"Config symlink/directory already exists: {config_dir}")
                            else:
                                logger.warning(f"Could not create config symlink: {e}")
                                # Try creating a regular directory and copying files as fallback
                                try:
                                    config_dir.mkdir(parents=True, exist_ok=True)
                                    import shutil
                                    for config_file in config_base_dir.glob('*.config'):
                                        shutil.copy2(config_file, config_dir / config_file.name)
                                    logger.info(f"Copied config files from {config_base_dir} to {config_dir} as fallback")
                                except Exception as e2:
                                    logger.warning(f"Could not copy config files as fallback: {e2}")
            
            step_results.append(WorkflowStepResult(
                step_name=step_name,
                return_code=0,
                success=True,
                command="(skipped - already completed)"
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

            # Ensure config symlink exists before select step
            if step_name == 'select':
                config_base_dir = config.work_dir / 'config_base'
                config_dir = config.work_dir / 'config'
                if config_base_dir.exists() and not (config_dir.exists() or config_dir.is_symlink()):
                    try:
                        config_dir.symlink_to(config_base_dir.resolve())
                        logger.info(f"Created config symlink before select step: {config_dir} -> {config_base_dir.resolve()}")
                    except Exception as e:
                        logger.warning(f"Could not create config symlink: {e}")
                        # Try copying as fallback
                        try:
                            config_dir.mkdir(parents=True, exist_ok=True)
                            import shutil
                            for config_file in config_base_dir.glob('*.config'):
                                shutil.copy2(config_file, config_dir / config_file.name)
                            logger.info(f"Copied config files as fallback")
                        except Exception as e2:
                            logger.warning(f"Could not copy config files: {e2}")
                # Update step_params to use absolute path to config directory
                # Always prefer config/ (symlink) over config_base/ for select step compatibility
                # Use absolute path of config_dir itself (not resolved) so amalgkit uses the symlink
                if (config_dir.exists() or config_dir.is_symlink()) and config_base_dir.exists():
                    step_params['config_dir'] = str(config_dir.absolute())
                    logger.debug(f"Using config directory (symlink): {config_dir.absolute()}")
                elif config_base_dir.exists():
                    # Fallback: use config_base directly if symlink doesn't exist
                    step_params['config_dir'] = str(config_base_dir.absolute())
                    logger.debug(f"Using config_base directory (fallback): {config_base_dir.absolute()}")
            
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
            
            # Pre-step prerequisite validation
            if step_name == 'integrate':
                # Check if FASTQ files exist before integrate
                steps_config = config.extra_config.get('steps', {})
                fastq_dir_raw = steps_config.get('getfastq', {}).get('out_dir', config.work_dir / "fastq")
                fastq_dir = Path(fastq_dir_raw)
                # Check for amalgkit getfastq/ subdirectory structure
                if fastq_dir.name != "getfastq":
                    getfastq_subdir = fastq_dir / "getfastq"
                    if getfastq_subdir.exists():
                        fastq_dir = getfastq_subdir
                
                # Check if any FASTQ files exist
                fastq_files = list(fastq_dir.glob("**/*.fastq*"))
                if not fastq_files:
                    error_msg = (
                        f"PREREQUISITE CHECK FAILED: No FASTQ files found before integrate step.\n"
                        f"  - Expected location: {fastq_dir}\n"
                        f"  - Checked for: *.fastq, *.fastq.gz, *.fq, *.fq.gz\n\n"
                        f"REMEDIATION:\n"
                        f"  1. Ensure getfastq step completed successfully\n"
                        f"  2. Check validation report: {config.work_dir / 'validation' / 'getfastq_validation.json'}\n"
                        f"  3. Re-run getfastq step if FASTQ files are missing\n"
                        f"  4. Verify amalgkit getfastq extracted files correctly"
                    )
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
            
            elif step_name == 'quant':
                # Check if integrate completed (quant needs integrated metadata)
                integrated_metadata = config.work_dir / 'integration' / 'integrated_metadata.json'
                # Also check for metadata.tsv that amalgkit integrate creates
                metadata_tsv = config.work_dir / 'metadata' / 'metadata.tsv'
                if not integrated_metadata.exists() and not metadata_tsv.exists():
                    error_msg = (
                        f"PREREQUISITE CHECK FAILED: Integrate step must complete before quant.\n"
                        f"  - Expected: {integrated_metadata} or {metadata_tsv}\n\n"
                        f"REMEDIATION:\n"
                        f"  1. Ensure integrate step completed successfully\n"
                        f"  2. Re-run integrate step if needed\n"
                    )
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
                
                # Check if FASTQ files exist (quant needs them for quantification)
                steps_config = config.extra_config.get('steps', {})
                fastq_dir_raw = steps_config.get('getfastq', {}).get('out_dir', config.work_dir / "fastq")
                fastq_dir = Path(fastq_dir_raw)
                if fastq_dir.name != "getfastq":
                    getfastq_subdir = fastq_dir / "getfastq"
                    if getfastq_subdir.exists():
                        fastq_dir = getfastq_subdir
                
                fastq_files = list(fastq_dir.glob("**/*.fastq*"))
                if not fastq_files:
                    error_msg = (
                        f"PREREQUISITE CHECK FAILED: No FASTQ files found before quant step.\n"
                        f"  - Expected location: {fastq_dir}\n\n"
                        f"REMEDIATION:\n"
                        f"  1. Ensure getfastq and integrate steps completed successfully\n"
                        f"  2. Check validation reports in: {config.work_dir / 'validation'}\n"
                        f"  3. Re-run getfastq step if FASTQ files are missing"
                    )
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
                
                # Check quantification tools availability
                try:
                    from metainformant.rna.deps import check_quantification_tools
                    quant_tools = check_quantification_tools()
                    available_tools = [tool for tool, (avail, _) in quant_tools.items() if avail]
                    if not available_tools:
                        error_msg = (
                            f"PREREQUISITE CHECK FAILED: No quantification tools available.\n"
                            f"  - Checked: kallisto, salmon\n"
                            f"  - Status: {quant_tools}\n\n"
                            f"REMEDIATION:\n"
                            f"  1. Install kallisto: conda install -c bioconda kallisto\n"
                            f"     OR: apt-get install kallisto\n"
                            f"  2. Install salmon: conda install -c bioconda salmon\n"
                            f"     OR: apt-get install salmon\n"
                            f"  3. Verify tools are in PATH: which kallisto / which salmon"
                        )
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
                except Exception as e:
                    logger.warning(f"Could not check quantification tools: {e}")
            
            elif step_name == 'merge':
                # Check if quant files exist before merge
                steps_config = config.extra_config.get('steps', {})
                quant_dir_raw = steps_config.get('quant', {}).get('out_dir', config.work_dir / "quant")
                quant_dir = Path(quant_dir_raw)
                
                # Check for quant output files (abundance.tsv or quant.sf)
                quant_files = list(quant_dir.glob("**/abundance.tsv")) + list(quant_dir.glob("**/quant.sf"))
                if not quant_files:
                    error_msg = (
                        f"PREREQUISITE CHECK FAILED: No quantification files found before merge step.\n"
                        f"  - Expected location: {quant_dir}\n"
                        f"  - Checked for: abundance.tsv, quant.sf\n\n"
                        f"REMEDIATION:\n"
                        f"  1. Ensure quant step completed successfully\n"
                        f"  2. Check validation report: {config.work_dir / 'validation' / 'quant_validation.json'}\n"
                        f"  3. Re-run quant step if quantification files are missing"
                    )
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
                
                # Check R dependencies for merge (optional but recommended)
                try:
                    from metainformant.rna.environment import check_rscript
                    r_available, r_message = check_rscript()
                    if not r_available:
                        logger.warning(
                            f"R/Rscript not available - merge step may fail if R plotting is required.\n"
                            f"  Status: {r_message}\n"
                            f"  Installation: apt-get install r-base-core or conda install r-base\n"
                            f"  Note: Merge can work without R, but plots will be skipped."
                        )
                    else:
                        # Check if ggplot2 package is available
                        try:
                            import subprocess
                            check_ggplot2 = subprocess.run(
                                ["Rscript", "-e", "library(ggplot2)"],
                                capture_output=True,
                                text=True,
                                timeout=10
                            )
                            if check_ggplot2.returncode != 0:
                                logger.warning(
                                    f"R package 'ggplot2' not available - merge step plotting may fail.\n"
                                    f"  Installation: Rscript -e \"install.packages('ggplot2', repos='https://cloud.r-project.org')\"\n"
                                    f"  Note: Merge can work without ggplot2, but plots will be skipped."
                                )
                        except Exception:
                            pass  # R check failed, but that's okay - merge can still work
                except Exception as e:
                    logger.debug(f"Could not check R dependencies: {e}")

            # Sanitize parameters before passing to step function
            sanitized_params = sanitize_params_for_cli(step_name, step_params)
            # Pass monitoring/progress preferences through to amalgkit wrappers.
            # `show_progress` controls tqdm progress bars (if available) and heartbeat cadence.
            result = step_func(
                sanitized_params,
                show_progress=progress,
                heartbeat_interval=5,
            )

            # After config step, create symlink from config/ to config_base/ for select step compatibility
            if step_name == 'config':
                config_base_dir = config.work_dir / 'config_base'
                config_dir = config.work_dir / 'config'
                if config_base_dir.exists() and not config_dir.exists():
                    try:
                        # Use absolute path for symlink to ensure it resolves correctly
                        config_dir.symlink_to(config_base_dir.resolve())
                        logger.info(f"Created symlink: {config_dir} -> {config_base_dir.resolve()} (for select step compatibility)")
                    except Exception as e:
                        logger.warning(f"Could not create config symlink: {e}")
                        # Try creating a regular directory and copying files as fallback
                        try:
                            config_dir.mkdir(parents=True, exist_ok=True)
                            import shutil
                            for config_file in config_base_dir.glob('*.config'):
                                shutil.copy2(config_file, config_dir / config_file.name)
                            logger.info(f"Copied config files from {config_base_dir} to {config_dir} as fallback")
                        except Exception as e2:
                            logger.warning(f"Could not copy config files as fallback: {e2}")

            # After select step, create filtered metadata for downstream steps
            if step_name == 'select' and result.returncode == 0:
                try:
                    logger.info("Creating filtered metadata for selected samples (excluding LITE files)")
                    try:
                        filter_selected_metadata(
                            config.work_dir / 'metadata' / 'metadata.tsv',
                            config.work_dir / 'metadata' / 'metadata_selected.tsv',
                            exclude_lite_files=True  # Automatically filter out LITE SRA files
                        )
                    except ValueError as e:
                        if "No samples meet the filtering criteria" in str(e):
                            logger.warning("All selected samples are LITE files. Creating metadata_selected.tsv without LITE filtering for this run.")
                            # Fall back to not excluding LITE files if all samples would be filtered out
                            filter_selected_metadata(
                                config.work_dir / 'metadata' / 'metadata.tsv',
                                config.work_dir / 'metadata' / 'metadata_selected.tsv',
                                exclude_lite_files=False  # Don't exclude LITE files if all would be filtered
                            )
                        else:
                            raise
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
                # Check if step actually produced expected outputs despite error
                is_completed, completion_indicator = _is_step_completed(step_name, step_params, config)
                if is_completed:
                    logger.warning(f"Step {step_name} reported error (return code {result.returncode}) but outputs exist: {completion_indicator}")
                    logger.warning(f"Continuing workflow - step appears to have completed successfully despite error")
                    # Mark as successful despite error
                    result.returncode = 0
                    logger.info(f"Step {step_name} marked as successful based on output validation")
                else:
                    error_message = f"Step failed with return code {result.returncode}"
                    if result.stderr:
                        error_message += f": {result.stderr}"
                    logger.error(f"Step {step_name} failed with return code {result.returncode}")
                    if result.stderr:
                        logger.error(f"Error output: {result.stderr}")
            else:
                logger.info(f"Step {step_name} completed successfully")
                
                # Add step summary for getfastq to show skipped vs downloaded files
                if step_name == 'getfastq':
                    # Check both stdout and stderr (amalgkit may output to either)
                    output_text = (result.stdout or "") + (result.stderr or "")
                    if output_text:
                        _log_getfastq_summary(output_text, logger)
                    
                    # Run validation after getfastq step
                    try:
                        from metainformant.rna.validation import validate_all_samples, save_validation_report
                        validation_result = validate_all_samples(config, stage='extraction')
                        validation_dir = config.work_dir / "validation"
                        validation_dir.mkdir(parents=True, exist_ok=True)
                        save_validation_report(validation_result, validation_dir / "getfastq_validation.json")
                        
                        # Log validation summary
                        total = validation_result.get('total_samples', 0)
                        validated = validation_result.get('validated', 0)
                        failed = validation_result.get('failed', 0)
                        logger.info(f"Validation after getfastq: {validated}/{total} samples have FASTQ files extracted")
                        
                        # Early exit if critical failure: no FASTQ files extracted
                        if validated == 0 and failed > 0 and total > 0:
                            # Check if selected samples are LITE files
                            lite_check_msg = ""
                            try:
                                import pandas as pd
                                selected_meta = config.work_dir / 'metadata' / 'metadata_selected.tsv'
                                if selected_meta.exists():
                                    df = pd.read_csv(selected_meta, sep='\t')
                                    if 'run' in df.columns:
                                        runs = df['run'].tolist()
                                        # Check AWS/GCP links for .lite indicators
                                        lite_runs = []
                                        for _, row in df.iterrows():
                                            aws_link = str(row.get('AWS_Link', ''))
                                            gcp_link = str(row.get('GCP_Link', ''))
                                            if '.lite' in aws_link or '.lite' in gcp_link or '.sralite' in aws_link or '.sralite' in gcp_link:
                                                lite_runs.append(row.get('run', 'unknown'))
                                        if lite_runs:
                                            lite_check_msg = (
                                                f"\n  ⚠️  LITE FILE DETECTION:\n"
                                                f"  - {len(lite_runs)}/{total} selected samples are LITE files (metadata-only, no sequence data)\n"
                                                f"  - LITE samples: {', '.join(lite_runs[:5])}{'...' if len(lite_runs) > 5 else ''}\n"
                                                f"  - LITE files cannot be extracted to FASTQ format\n"
                                                f"  - SOLUTION: Re-run select step with LITE filtering enabled or select different samples\n"
                                            )
                            except Exception:
                                pass  # Ignore errors in LITE detection
                            
                            error_msg = (
                                f"CRITICAL: getfastq step failed to extract FASTQ files for any samples.\n"
                                f"  - Total samples: {total}\n"
                                f"  - Samples with FASTQ files: {validated}\n"
                                f"  - Samples missing FASTQ files: {failed}\n"
                                f"{lite_check_msg}\n"
                                f"REMEDIATION STEPS:\n"
                                f"  1. Check if selected samples are LITE files (see above)\n"
                                f"  2. Check amalgkit getfastq logs: {config.work_dir / 'logs' / 'getfastq.log'}\n"
                                f"  3. Verify SRA files were downloaded: {config.work_dir / 'fastq' / 'getfastq'}\n"
                                f"  4. Check if fastp/fasterq-dump tools are available:\n"
                                f"     - fasterq-dump: shutil.which('fasterq-dump')\n"
                                f"     - fastp: shutil.which('fastp')\n"
                                f"  5. Check amalgkit command output for extraction errors\n"
                                f"  6. Try re-running with redo: yes if files may be corrupted\n"
                                f"  7. Verify disk space is sufficient for FASTQ extraction\n\n"
                                f"Workflow stopped to prevent cascading failures in downstream steps."
                            )
                            logger.error(error_msg)
                            step_results.append(WorkflowStepResult(
                                step_name="getfastq_validation",
                                return_code=1,
                                success=False,
                                error_message=error_msg
                            ))
                            # Stop workflow unless check=False (user wants to continue anyway)
                            if check:
                                break
                        elif failed > 0:
                            logger.warning(f"{failed} samples missing FASTQ files after getfastq step")
                    except Exception as e:
                        logger.warning(f"Validation after getfastq failed: {e}")
                
                # Run validation after quant step
                if step_name == 'quant':
                    try:
                        from metainformant.rna.validation import validate_all_samples, save_validation_report
                        validation_result = validate_all_samples(config, stage='quantification')
                        validation_dir = config.work_dir / "validation"
                        validation_dir.mkdir(parents=True, exist_ok=True)
                        save_validation_report(validation_result, validation_dir / "quant_validation.json")
                        
                        # Log validation summary
                        total = validation_result.get('total_samples', 0)
                        validated = validation_result.get('validated', 0)
                        failed = validation_result.get('failed', 0)
                        logger.info(f"Validation after quant: {validated}/{total} samples quantified")
                        if failed > 0:
                            logger.warning(f"{failed} samples missing quantification files after quant step")
                    except Exception as e:
                        logger.warning(f"Validation after quant failed: {e}")

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
            
            # If exception is tqdm-related, the subprocess might still be running
            # Wait longer for subprocess to complete and files to be written
            import time
            if "tqdm" in str(e).lower() or "refresh" in str(e).lower():
                logger.warning(f"tqdm error detected - waiting for subprocess to complete (up to 30 seconds)")
                # Wait longer for subprocess to finish (metadata step can take time)
                for i in range(30):
                    time.sleep(1)
                    is_completed, completion_indicator = _is_step_completed(step_name, step_params, config)
                    if is_completed:
                        break
                    if i % 5 == 0:
                        logger.debug(f"Waiting for {step_name} outputs... ({i+1}/30 seconds)")
            else:
                # For other exceptions, wait briefly
                time.sleep(2)
            
            # Check if step actually produced expected outputs despite exception
            is_completed, completion_indicator = _is_step_completed(step_name, step_params, config)
            if is_completed:
                logger.warning(f"Step {step_name} raised exception but outputs exist: {completion_indicator}")
                logger.warning(f"Continuing workflow - step appears to have completed successfully despite exception")
                # Mark as successful despite exception
                step_results.append(WorkflowStepResult(
                    step_name=step_name,
                    return_code=0,
                    success=True,
                    error_message=f"Exception occurred but outputs validated: {error_msg}",
                    command=command_str
                ))
                logger.info(f"Step {step_name} marked as successful based on output validation")
                # Continue to next step (don't break workflow)
                continue
            else:
                step_results.append(WorkflowStepResult(
                    step_name=step_name,
                    return_code=1,
                    success=False,
                    error_message=error_msg
                ))
                if check:
                    break
                continue

    successful_steps = sum(1 for sr in step_results if sr.success)
    total_steps = len(step_results)
    failed_steps = total_steps - successful_steps

    # Generate workflow summary with remediation steps
    if progress:
        logger.info(f"Workflow completed: {successful_steps}/{total_steps} steps successful")
        
        if failed_steps > 0:
            logger.info("=" * 80)
            logger.info("WORKFLOW SUMMARY")
            logger.info("=" * 80)
            
            failed_step_results = [sr for sr in step_results if not sr.success]
            logger.info(f"\nFailed steps: {len(failed_step_results)}/{total_steps}")
            
            for failed_step in failed_step_results:
                logger.info(f"\n  ❌ {failed_step.step_name}")
                if failed_step.error_message:
                    # Log error message with indentation for readability
                    for line in failed_step.error_message.split('\n'):
                        logger.info(f"     {line}")
            
            logger.info("\n" + "=" * 80)
            logger.info("REMEDIATION STEPS")
            logger.info("=" * 80)
            
            # Provide step-specific remediation guidance
            remediation_steps = []
            
            for failed_step in failed_step_results:
                step_name = failed_step.step_name
                
                if step_name == 'getfastq' or step_name == 'getfastq_validation':
                    remediation_steps.append(
                        f"\n  For {step_name}:\n"
                        f"    1. Check logs: {config.work_dir / 'logs' / 'getfastq.log'}\n"
                        f"    2. Check validation: {config.work_dir / 'validation' / 'getfastq_validation.json'}\n"
                        f"    3. Verify SRA files downloaded: ls -lh {config.work_dir / 'fastq' / 'getfastq'}\n"
                        f"    4. Check tool availability: which fasterq-dump which fastp\n"
                        f"    5. Re-run: python3 scripts/rna/run_workflow.py {config.work_dir.parent.name} --steps getfastq"
                    )
                elif step_name == 'integrate':
                    remediation_steps.append(
                        f"\n  For {step_name}:\n"
                        f"    1. Ensure getfastq completed successfully\n"
                        f"    2. Check FASTQ files exist: find {config.work_dir / 'fastq'} -name '*.fastq*'\n"
                        f"    3. Re-run: python3 scripts/rna/run_workflow.py {config.work_dir.parent.name} --steps getfastq integrate"
                    )
                elif step_name == 'quant':
                    remediation_steps.append(
                        f"\n  For {step_name}:\n"
                        f"    1. Ensure getfastq and integrate completed successfully\n"
                        f"    2. Check quantification tools: which kallisto which salmon\n"
                        f"    3. Check validation: {config.work_dir / 'validation' / 'quant_validation.json'}\n"
                        f"    4. Re-run: python3 scripts/rna/run_workflow.py {config.work_dir.parent.name} --steps quant"
                    )
                elif step_name == 'merge':
                    remediation_steps.append(
                        f"\n  For {step_name}:\n"
                        f"    1. Ensure quant completed successfully\n"
                        f"    2. Check quant files: find {config.work_dir / 'quant'} -name 'abundance.tsv' -o -name 'quant.sf'\n"
                        f"    3. If R plotting failed, install: Rscript -e \"install.packages('ggplot2', repos='https://cloud.r-project.org')\"\n"
                        f"    4. Re-run: python3 scripts/rna/run_workflow.py {config.work_dir.parent.name} --steps merge"
                    )
                elif step_name == 'curate':
                    remediation_steps.append(
                        f"\n  For {step_name}:\n"
                        f"    1. Ensure merge completed successfully\n"
                        f"    2. Check merge output: {config.work_dir / 'merge' / 'merge' / 'metadata.tsv'}\n"
                        f"    3. Install R if missing: apt-get install r-base-core\n"
                        f"    4. Re-run: python3 scripts/rna/run_workflow.py {config.work_dir.parent.name} --steps curate"
                    )
                else:
                    remediation_steps.append(
                        f"\n  For {step_name}:\n"
                        f"    1. Check logs: {config.work_dir / 'logs' / f'{step_name}.log'}\n"
                        f"    2. Review error message above for specific guidance\n"
                        f"    3. Re-run: python3 scripts/rna/run_workflow.py {config.work_dir.parent.name} --steps {step_name}"
                    )
            
            for remediation in remediation_steps:
                logger.info(remediation)
            
            logger.info("\n" + "=" * 80)
            logger.info("INDEPENDENT STEP RE-RUN")
            logger.info("=" * 80)
            logger.info(
                "You can re-run individual steps without re-running the entire workflow:\n"
                f"  python3 scripts/rna/run_workflow.py <config_file> --steps <step_name>\n\n"
                "Example:\n"
                f"  python3 scripts/rna/run_workflow.py config/amalgkit/amalgkit_pbarbatus_5sample.yaml --steps getfastq"
            )
            logger.info("=" * 80)

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







