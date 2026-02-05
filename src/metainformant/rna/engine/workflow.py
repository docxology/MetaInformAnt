"""RNA-seq workflow orchestration and configuration management.

This module provides high-level functions for managing RNA-seq analysis workflows,
including configuration loading, workflow planning, and execution orchestration.
"""

from __future__ import annotations

import csv
import json
import math
import os
import shutil
import time as time_mod
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, Union

from metainformant.core import io, logging
from metainformant.core.utils.config import load_mapping_from_file
from metainformant.rna.amalgkit.metadata_filter import filter_selected_metadata
from metainformant.rna.amalgkit.metadata_utils import deduplicate_metadata
from metainformant.rna.engine.sra_extraction import (
    extract_sra_directly,
    manual_integration_fallback,
)
from metainformant.rna.engine.workflow_cleanup import (
    check_disk_space,
    check_disk_space_or_fail,
    cleanup_after_quant,
    cleanup_fastqs,
    cleanup_incorrectly_placed_sra_files,
    cleanup_temp_files,
    filter_metadata_for_unquantified,
    get_quantified_samples,
)
from metainformant.rna.engine.workflow_steps import (
    check_step_completion_status,
    handle_post_step_actions,
    log_workflow_summary,
    setup_vdb_config,
    validate_step_prerequisites,
)

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
    import os

    env_threads = os.environ.get("AK_THREADS")
    if env_threads:
        try:
            result["threads"] = int(env_threads)
        except ValueError:
            pass

    return result


def apply_step_defaults(config: AmalgkitWorkflowConfig) -> AmalgkitWorkflowConfig:
    """Apply default values to workflow step configurations.

    Each amalgkit workflow step has sensible defaults that are applied
    if not explicitly configured. This function ensures consistent
    behavior across workflows.

    Args:
        config: Workflow configuration

    Returns:
        Configuration with step defaults applied

    Step defaults applied:
        - getfastq: out_dir, layout detection, retry settings
        - integrate: fastq_dir, out_dir
        - quant: out_dir, index settings, bootstrap count
        - merge: out_dir, batch size
        - cstmm: out_dir, normalization method
        - curate: out_dir, metadata path
    """
    # Define step-specific defaults
    step_defaults: Dict[str, Dict[str, Any]] = {
        "getfastq": {
            "out_dir": str(config.work_dir / "fastq"),
            "layout": "auto",  # Auto-detect paired/single-end
            "num_retries": 3,
            "retry_delay": 30,  # seconds between retries
            "validate_md5": True,
        },
        "integrate": {
            "out_dir": str(config.work_dir),
            # fastq_dir is computed dynamically in plan_workflow
        },
        "quant": {
            "out_dir": str(config.work_dir),
            "bootstrap": 100,  # Bootstrap samples for kallisto
            "fragment_length": 200,  # For single-end reads
            "fragment_sd": 20,  # Standard deviation for single-end
        },
        "merge": {
            "out_dir": str(config.work_dir / "merged"),
            "batch_size": 100,  # Samples per batch for large datasets
            "output_format": "tsv",
        },
        "cstmm": {
            "out_dir": str(config.work_dir),
            "normalization": "tpm",  # TPM normalization by default
            "min_expression": 1.0,  # Minimum expression threshold
        },
        "curate": {
            "out_dir": str(config.work_dir / "curate"),
            # metadata path computed from merge output in plan_workflow
        },
        "csca": {
            "out_dir": str(config.work_dir),
            "correlation_method": "pearson",
        },
    }

    # Get current per_step configuration
    per_step = config.per_step.copy() if config.per_step else {}

    # Apply defaults for each step that has defaults defined
    for step_name, defaults in step_defaults.items():
        if step_name not in per_step:
            per_step[step_name] = {}

        # Apply defaults only for keys not already set
        for key, default_value in defaults.items():
            if key not in per_step[step_name]:
                per_step[step_name][key] = default_value

    # Handle environment variable overrides for step settings
    env_prefix_map = {
        "getfastq": "AK_GETFASTQ_",
        "integrate": "AK_INTEGRATE_",
        "quant": "AK_QUANT_",
        "merge": "AK_MERGE_",
        "cstmm": "AK_CSTMM_",
        "curate": "AK_CURATE_",
    }

    for step_name, prefix in env_prefix_map.items():
        if step_name in per_step:
            for key in per_step[step_name]:
                env_key = f"{prefix}{key.upper()}"
                env_value = os.environ.get(env_key)
                if env_value is not None:
                    # Try to parse as int or float if applicable
                    try:
                        if "." in env_value:
                            per_step[step_name][key] = float(env_value)
                        else:
                            per_step[step_name][key] = int(env_value)
                    except ValueError:
                        # Keep as string if not numeric
                        per_step[step_name][key] = env_value

    # Update config with merged per_step settings
    config.per_step = per_step

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
        "metadata",
        "config",
        "select",
        "getfastq",
        "integrate",  # After getfastq to integrate downloaded FASTQ files
        "quant",
        "merge",
        "cstmm",
        "curate",
        "csca",
        "sanity",
    ]

    # Get per-step overrides if provided
    per_step = config.per_step or config.extra_config.get("steps", {}) or config.extra_config.get("per_step", {})

    # per_step is a dict of step_name -> params
    per_step_dict = per_step if isinstance(per_step, dict) else {}
    # per_step_list is a list of step names to run
    per_step_list = per_step if isinstance(per_step, list) else []

    # Check if cstmm/csca required parameters are provided
    has_ortholog_params = (
        (
            per_step_dict.get("cstmm", {}).get("orthogroup_table") is not None
            or per_step_dict.get("cstmm", {}).get("dir_busco") is not None
        )
        or (
            per_step_dict.get("csca", {}).get("orthogroup_table") is not None
            or per_step_dict.get("csca", {}).get("dir_busco") is not None
        )
        or config.extra_config.get("orthogroup_table") is not None
        or config.extra_config.get("dir_busco") is not None
    )

    # Filter out cstmm and csca if required parameters not provided
    if not has_ortholog_params:
        workflow_steps = [step for step in workflow_steps if step not in ("cstmm", "csca")]
        if "cstmm" in per_step_dict or "csca" in per_step_dict or "cstmm" in per_step_list or "csca" in per_step_list:
            logger.warning(
                "Skipping cstmm and csca steps: required parameters (orthogroup_table or dir_busco) not provided"
            )

    # If steps were explicitly specified as a list, only run those steps
    if per_step_list:
        workflow_steps = [step for step in workflow_steps if step in per_step_list]

    for step in workflow_steps:
        # Start with defaults from common config
        step_params = {
            "out_dir": str(config.work_dir),
        }
        # Propagate global redo flag if set
        if config.extra_config.get("redo"):
            step_params["redo"] = "yes"
        # Only add threads for steps that support it (merge, curate, metadata don't)
        if step not in ("merge", "curate", "metadata", "sanity", "select", "cstmm", "csca", "config"):
            step_params["threads"] = config.threads
        # Only add species for steps that support it (merge, curate, integrate, getfastq, quant don't)
        if config.species_list and step not in (
            "integrate",
            "merge",
            "curate",
            "metadata",
            "getfastq",
            "quant",
            "sanity",
            "select",
            "cstmm",
            "csca",
            "config",
        ):
            step_params["species"] = config.species_list

        # Apply step-specific overrides
        if step in per_step_dict:
            overrides = per_step_dict[step].copy()

            step_params.update(overrides)
            logger.debug(f"Workflow step {step} params: {step_params}")

            # Warn if redo: yes is set for getfastq (usually unnecessary)
            if step == "getfastq":
                redo_value = step_params.get("redo", "no")
                if isinstance(redo_value, bool):
                    redo_is_yes = redo_value
                else:
                    redo_is_yes = str(redo_value).lower() in ("yes", "true", "1")

                if redo_is_yes:
                    logger.warning(
                        "getfastq step has redo: yes - this will force re-download of all files, "
                        "even if they already exist."
                    )

        # Ensure work_dir is always set (crucial for steps like merge to find inputs)
        if "work_dir" not in step_params:
            step_params["work_dir"] = str(config.work_dir)
        if "out_dir" not in step_params:
            step_params["out_dir"] = str(config.work_dir)

        # INTELLIGENT REDO LOGIC:
        if step == "getfastq":
            getfastq_out_dir = Path(step_params.get("out_dir", str(config.work_dir / "fastq")))
            if getfastq_out_dir.name != "getfastq":
                getfastq_subdir = getfastq_out_dir / "getfastq"
            else:
                getfastq_subdir = getfastq_out_dir

            sra_files_found = 0
            if getfastq_subdir.exists():
                sra_files_found = len(
                    list(getfastq_subdir.glob("**/*.sra")) + list(getfastq_subdir.glob("**/*.sra.part"))
                )

            if sra_files_found > 0:
                current_redo = step_params.get("redo", "no")
                is_redo_yes = str(current_redo).lower() in ("yes", "true", "1")

                if is_redo_yes:
                    logger.warning(
                        f"Found {sra_files_found} SRA files in {getfastq_subdir}. "
                        f"Overriding 'redo: yes' to 'redo: no' to prevent deletion/redownload. "
                        f"NOTE: If extraction skips inappropriately, please delete SRA files or force redo."
                    )
                    step_params["redo"] = "no"
                else:
                    logger.info(f"Verified {sra_files_found} SRA files exist. Keeping 'redo: no' to enable extraction.")

        # Auto-inject metadata path for select step
        if step == "select" and "metadata" not in step_params:
            # Select step needs metadata.tsv as input
            metadata_file = config.work_dir / "metadata" / "metadata.tsv"
            if metadata_file.exists():
                step_params["metadata"] = str(metadata_file.resolve())

        # Auto-inject metadata paths for steps that need them
        # Use filtered metadata if it exists (created after select step)
        if step in ("getfastq", "integrate", "merge"):
            filtered_metadata = config.work_dir / "metadata" / "metadata_selected.tsv"
            if filtered_metadata.exists():
                step_params["metadata"] = str(filtered_metadata)
            elif "metadata" not in step_params:
                step_params["metadata"] = str(config.work_dir / "metadata" / "metadata.tsv")
        elif step == "quant":
            # Quant uses metadata from integrate step (amalgkit creates metadata_updated_for_private_fastq.tsv after integrate)
            # Prefer integrated metadata if it exists, otherwise use selected metadata
            integrated_metadata = config.work_dir / "metadata" / "metadata_updated_for_private_fastq.tsv"
            old_integrated_metadata = config.work_dir / "metadata" / "metadata.tsv"
            filtered_metadata = config.work_dir / "metadata" / "metadata_selected.tsv"

            if integrated_metadata.exists():
                step_params["metadata"] = str(integrated_metadata)
            elif old_integrated_metadata.exists():
                step_params["metadata"] = str(old_integrated_metadata)
            elif filtered_metadata.exists():
                step_params["metadata"] = str(filtered_metadata)
            elif "metadata" not in step_params:
                step_params["metadata"] = str(config.work_dir / "metadata" / "metadata.tsv")

            # Inject fasta_dir and index_dir from config if not present
            if "fasta_dir" not in step_params:
                # Use same logic as prepare_reference_genome
                dest_dir = config.genome.get("dest_dir") if config.genome else None
                if dest_dir:
                    step_params["fasta_dir"] = str(dest_dir)
                else:
                    # Default fallback
                    step_params["fasta_dir"] = str(config.work_dir.parent / "genome")

            if "index_dir" not in step_params:
                step_params["index_dir"] = str(config.work_dir / "index")

        elif step == "curate" and "metadata" not in step_params:
            # Curate uses metadata from merge directory
            # We need to find where merge step puts its output
            merge_params = per_step.get("merge", {})
            merge_dir = Path(merge_params.get("out_dir", config.work_dir / "merged"))
            # amalgkit merge creates a 'merge' subdirectory inside out_dir
            step_params["metadata"] = str(merge_dir / "merge" / "metadata.tsv")

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
        if step == "integrate":
            if "fastq_dir" in step_params:
                fastq_dir = Path(step_params["fastq_dir"])

                # Check if fastq_dir itself is the getfastq directory FIRST
                if fastq_dir.name == "getfastq":
                    # Already pointing to getfastq, no adjustment needed
                    logger.debug(f"Integrate fastq_dir already points to getfastq directory: {fastq_dir}")
                else:
                    # Check if getfastq subdirectory exists (amalgkit getfastq creates this)
                    getfastq_subdir = fastq_dir / "getfastq"
                    if getfastq_subdir.exists():
                        step_params["fastq_dir"] = str(getfastq_subdir)
                        logger.debug(
                            f"Adjusted integrate fastq_dir to include getfastq subdirectory: {getfastq_subdir}"
                        )
            else:
                # If fastq_dir not specified, try to infer from getfastq step
                getfastq_out_dir = per_step_dict.get("getfastq", {}).get("out_dir")
                if getfastq_out_dir:
                    fastq_dir = Path(getfastq_out_dir)
                    # If specified out_dir ends in getfastq, use it directly
                    if fastq_dir.name == "getfastq":
                        step_params["fastq_dir"] = str(fastq_dir)
                        logger.debug(f"Inferred integrate fastq_dir from getfastq out_dir (direct): {fastq_dir}")
                    else:
                        getfastq_subdir = fastq_dir / "getfastq"
                        if getfastq_subdir.exists():
                            step_params["fastq_dir"] = str(getfastq_subdir)
                            logger.debug(f"Inferred integrate fastq_dir from getfastq step: {getfastq_subdir}")

        steps.append((step, step_params))

    logger.info(f"Planned workflow with {len(steps)} steps")
    return steps


def plan_workflow_with_params(
    config: AmalgkitWorkflowConfig, param_overrides: Optional[Dict[str, Any]] = None, **kwargs: Any
) -> List[Tuple[str, Any]]:
    """Plan workflow with parameter overrides.

    Args:
        config: Base workflow configuration
        param_overrides: Parameter overrides for specific steps
        **kwargs: Additional overrides

    Returns:
        List of (step_name, parameters) tuples
    """
    if param_overrides:
        if config.per_step is None:
            config.per_step = {}
        config.per_step.update(param_overrides)

    return plan_workflow(config)


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
            for line in output_text.split("\n"):
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


# Private aliases for backward compatibility with internal callers
_cleanup_incorrectly_placed_sra_files = cleanup_incorrectly_placed_sra_files
_cleanup_temp_files = cleanup_temp_files
_check_disk_space = check_disk_space
_check_disk_space_or_fail = check_disk_space_or_fail


def _log_heartbeat(step_name: str, start_time: float, message: str = "") -> None:
    elapsed = time_mod.time() - start_time
    elapsed_str = f"{int(elapsed // 3600)}h{int((elapsed % 3600) // 60)}m{int(elapsed % 60)}s"

    # Check disk space
    ok, free_gb = _check_disk_space(Path.cwd(), min_free_gb=5.0)
    disk_status = f"{free_gb:.1f}GB free" if free_gb >= 0 else "unknown"

    msg = f"[HEARTBEAT] {step_name}: {elapsed_str} elapsed, disk: {disk_status}"
    if message:
        msg += f" | {message}"
    logger.info(msg)
    # Force flush to ensure visibility during long blocking operations
    for handler in logger.handlers:
        handler.flush()


def prepare_extraction_directories(getfastq_dir: Path) -> int:
    """Prepare directories for FASTQ extraction by cleaning empty sample dirs.

    Amalgkit with redo: no skips extraction if sample directories exist,
    even if they're empty or only contain SRA files without FASTQ output.
    This function:
    1. Finds all sample directories (SRR*, ERR*, DRR*)
    2. For each dir: if SRA present but no .amalgkit.fastq.gz -> move SRA to root sra/ dir
    3. Delete the empty sample directory so amalgkit will re-create and extract

    Args:
        getfastq_dir: Path to the getfastq output directory

    Returns:
        Number of directories cleaned
    """
    if not getfastq_dir.exists():
        logger.debug(f"getfastq directory does not exist: {getfastq_dir}")
        return 0

    # Ensure we have a central sra/ directory for preserving downloaded files
    sra_backup_dir = getfastq_dir / "sra"
    sra_backup_dir.mkdir(parents=True, exist_ok=True)

    cleaned_count = 0
    sample_dirs = list(getfastq_dir.glob("SRR*")) + list(getfastq_dir.glob("ERR*")) + list(getfastq_dir.glob("DRR*"))

    for sample_dir in sample_dirs:
        if not sample_dir.is_dir():
            continue

        # Check for FASTQ output
        fastq_files = list(sample_dir.glob("*.amalgkit.fastq.gz"))
        if fastq_files:
            # FASTQ already extracted, skip
            continue

        # Check for SRA files
        sra_files = list(sample_dir.glob("*.sra"))

        if sra_files:
            # SRA present but no FASTQ - need to force extraction
            # Move SRA to backup location before deleting dir
            for sra_file in sra_files:
                backup_path = sra_backup_dir / sra_file.name
                if not backup_path.exists():
                    logger.info(f"Moving SRA for extraction: {sra_file.name} -> sra/")
                    shutil.move(str(sra_file), str(backup_path))
                else:
                    # Already backed up, just remove
                    sra_file.unlink()

        # Remove the empty sample directory so amalgkit will re-create it during extraction
        try:
            # Remove any remaining files in the directory
            for item in sample_dir.iterdir():
                if item.is_file():
                    item.unlink()
                elif item.is_dir():
                    shutil.rmtree(item)
            sample_dir.rmdir()
            logger.info(f"Cleaned empty sample directory for re-extraction: {sample_dir.name}")
            cleaned_count += 1
        except Exception as e:
            logger.warning(f"Could not clean directory {sample_dir}: {e}")

    if cleaned_count > 0:
        logger.info(f"Prepared {cleaned_count} directories for FASTQ extraction")
        logger.info(f"SRA files preserved in: {sra_backup_dir}")

    return cleaned_count


def create_extraction_metadata(getfastq_dir: Path, source_metadata: Path, output_metadata: Path) -> int:
    """Create filtered metadata containing only samples with existing SRA files.

    This function:
    1. Reads the source metadata file
    2. Scans getfastq_dir/sra/ for available SRA files
    3. Filters metadata to only include samples with SRA files present
    4. Moves SRA files from sra/ to individual sample directories for amalgkit
    5. Writes filtered metadata for extraction-only run

    Args:
        getfastq_dir: Path to the getfastq output directory
        source_metadata: Path to the source metadata TSV file
        output_metadata: Path where filtered metadata will be written

    Returns:
        Number of samples with available SRA files
    """
    import csv

    sra_dir = getfastq_dir / "sra"
    if not sra_dir.exists():
        logger.warning(f"No SRA directory found at {sra_dir}")
        return 0

    # Find all available SRA files
    available_sra = {}
    for sra_file in sra_dir.glob("*.sra"):
        # Extract sample ID from filename (e.g., SRR12345.sra -> SRR12345)
        sample_id = sra_file.stem
        available_sra[sample_id] = sra_file

    if not available_sra:
        logger.warning("No SRA files found in sra/ directory")
        return 0

    logger.info(f"Found {len(available_sra)} SRA files available for extraction")

    # Read source metadata and filter
    try:
        with open(source_metadata, "r", newline="") as f:
            reader = csv.DictReader(f, delimiter="\t")
            fieldnames = reader.fieldnames
            rows = list(reader)
    except Exception as e:
        logger.error(f"Could not read metadata: {e}")
        return 0

    # Filter rows to only those with SRA files present
    filtered_rows = []
    moved_count = 0
    for row in rows:
        run_id = row.get("run", "")
        if run_id in available_sra:
            filtered_rows.append(row)

            # Move SRA file to sample directory where amalgkit expects it
            sra_source = available_sra[run_id]
            sample_dir = getfastq_dir / run_id
            sample_dir.mkdir(parents=True, exist_ok=True)
            sra_dest = sample_dir / f"{run_id}.sra"

            if not sra_dest.exists():
                try:
                    shutil.copy2(str(sra_source), str(sra_dest))
                    logger.info(f"Copied SRA to sample dir: {run_id}.sra")
                    moved_count += 1
                except Exception as e:
                    logger.warning(f"Could not copy SRA file {run_id}: {e}")

    if not filtered_rows:
        logger.warning("No samples in metadata have corresponding SRA files")
        return 0

    # Write filtered metadata
    try:
        with open(output_metadata, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
            writer.writeheader()
            writer.writerows(filtered_rows)
        logger.info(f"Created extraction metadata with {len(filtered_rows)} samples: {output_metadata}")
    except Exception as e:
        logger.error(f"Could not write filtered metadata: {e}")
        return 0

    return len(filtered_rows)


def prepare_reference_genome(config: AmalgkitWorkflowConfig) -> bool:
    """Download reference genome and build index if needed.

    This function checks if the kallisto index for the species currently exists.
    If not, it downloads the transcriptome FASTA from the URL specified in
    the config and builds the index using 'kallisto index'.

    Args:
        config: Workflow configuration

    Returns:
        True if index is ready (exists or created), False on failure
    """
    if not config.genome or "ftp_url" not in config.genome:
        logger.debug("No genome configuration or FTP URL found, skipping automated genome prep")
        return True

    try:
        species_name = config.species_list[0] if config.species_list else "species"

        # 1. Define paths
        # Index location: {work_dir}/index/{Species}_transcripts.idx
        index_dir = config.work_dir / "index"
        index_file = index_dir / f"{species_name}_transcripts.idx"

        # Also check if it exists in the genome dest_dir/index
        dest_dir = Path(config.genome.get("dest_dir", config.work_dir.parent / "genome"))
        shared_index_dir = dest_dir / "index"
        shared_index_file = shared_index_dir / f"{species_name}_transcripts.idx"

        if index_file.exists():
            logger.info(f"Reference genome index found at: {index_file}")
            return True

        if shared_index_file.exists():
            logger.info(f"Found shared reference genome index at: {shared_index_file}")
            # Ensure work index dir exists
            index_dir.mkdir(parents=True, exist_ok=True)
            # Link or copy it to the work index dir
            import shutil

            logger.info(f"Copying shared index to: {index_file}")
            shutil.copy2(shared_index_file, index_file)
            return True

        logger.info(f"Reference index missing at {index_file}. Preparing reference genome...")

        # 2. Get download URL
        # URL construction: ftp_url + files.transcriptome_fasta
        ftp_url = config.genome["ftp_url"].rstrip("/")
        transcriptome_file = config.genome.get("files", {}).get("transcriptome_fasta")

        if not transcriptome_file:
            logger.warning("No transcriptome_fasta file name in config genome.files")
            return False

        download_url = f"{ftp_url}/{transcriptome_file}"

        # 3. Create destination
        dest_dir = Path(config.genome.get("dest_dir", config.work_dir.parent / "genome"))
        dest_dir.mkdir(parents=True, exist_ok=True)
        dest_file = dest_dir / transcriptome_file

        # 4. Download file
        if not dest_file.exists():
            logger.info(f"Downloading transcriptome from {download_url}...")
            # Use curl for simplicity
            import subprocess

            cmd = ["curl", "-L", "-o", str(dest_file), download_url]
            result = subprocess.run(cmd, check=True)
            if result.returncode != 0:
                logger.error("Download failed")
                return False
            logger.info("Download complete")
        else:
            logger.info(f"Transcriptome file already exists at {dest_file}")

        # 5. Build Index
        index_dir.mkdir(parents=True, exist_ok=True)
        logger.info(f"Building kallisto index at {index_file}...")

        # Find kallisto executable (assume in path or env)
        kallisto_cmd = ["kallisto", "index", "-i", str(index_file), str(dest_file)]

        import subprocess

        result = subprocess.run(kallisto_cmd, capture_output=True, text=True)

        if result.returncode == 0:
            logger.info("Genome index built successfully")
            # Also copy it back to shared index dir for future use
            try:
                import shutil

                shared_index_dir.mkdir(parents=True, exist_ok=True)
                logger.info(f"Saving a copy of index to shared directory: {shared_index_file}")
                shutil.copy2(index_file, shared_index_file)
            except Exception as e:
                logger.warning(f"Failed to save copy to shared index dir: {e}")
            return True
        else:
            logger.error(f"Kallisto index failed: {result.stderr}")
            return False

    except Exception as e:
        logger.error(f"Failed to prepare reference genome: {e}")
        return False

        return False


def verify_getfastq_prerequisites(config: AmalgkitWorkflowConfig, steps_planned: List[Tuple[str, Any]]) -> None:
    """Verify prerequisites before running getfastq (e.g. SRA files exist).

    Args:
        config: Workflow configuration
        steps_planned: Planned steps
    """
    # Only run if getfastq is in the plan
    getfastq_params = next((params for step, params in steps_planned if step == "getfastq"), None)
    if not getfastq_params:
        return

    logger.info("--- Pre-flight Check: getfastq ---")

    # 1. Check directories
    out_dir = Path(getfastq_params.get("out_dir", str(config.work_dir / "fastq")))
    getfastq_dir = out_dir / "getfastq" if out_dir.name != "getfastq" else out_dir

    logger.info(f"Target Directory: {getfastq_dir}")
    if getfastq_dir.exists():
        sras = list(getfastq_dir.glob("**/*.sra"))
        parts = list(getfastq_dir.glob("**/*.sra.part"))
        fastqs = list(getfastq_dir.glob("**/*.fastq.gz*"))

        logger.info(f"Files Found:")
        logger.info(f"  - SRA files (.sra): {len(sras)}")
        logger.info(f"  - Partial downloads (.part): {len(parts)}")
        logger.info(f"  - Extracted FASTQ (.fastq.gz): {len(fastqs)}")

        if len(sras) > 0:
            logger.info("✓ SRA files detected. 'redo: no' should trigger extraction.")
        elif len(parts) > 0:
            logger.warning("! Only partial downloads found. Robust download may be incomplete.")
        elif len(fastqs) == 0:
            logger.warning("! No input SRA or output FASTQ files found. getfastq will attempt fresh download.")
    else:
        logger.info(f"Directory {getfastq_dir} does not exist yet.")

    logger.info("----------------------------------")


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

    # ALWAYS remove work_dir - amalgkit CLI does not accept --work-dir (only --out_dir)
    sanitized.pop("work_dir", None)

    # Remove invalid parameters per subcommand
    INVALID_PARAMS = {
        "quant": {"keep_fastq"},  # keep_fastq is not supported by amalgkit CLI
        "getfastq": {"num_download_workers"},  # Internal worker count alias
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
    steps_config = config.extra_config.get("steps", {})

    if step_name == "metadata":
        metadata_file = work_dir / "metadata" / "metadata.tsv"
        if metadata_file.exists():
            return True, str(metadata_file)
        return False, None

    elif step_name == "config":
        config_dir = work_dir / "config_base"
        if config_dir.exists():
            config_files = list(config_dir.glob("*.config"))
            if config_files:
                return True, f"{len(config_files)} config files in {config_dir}"
        return False, None

    elif step_name == "select":
        selected_metadata = work_dir / "metadata" / "metadata_selected.tsv"
        if selected_metadata.exists():
            return True, str(selected_metadata)
        return False, None

    elif step_name == "getfastq":
        # Defer to amalgkit's internal skipping logic (redo: no)
        # We don't want to skip the whole step if only some files are missing
        return False, None

    elif step_name == "integrate":
        # Check for integrated metadata
        integrated_meta = work_dir / "integration" / "integrated_metadata.json"
        if integrated_meta.exists():
            return True, str(integrated_meta)

        # Also check if metadata.tsv was updated (integrate updates it)
        metadata_tsv = work_dir / "metadata" / "metadata.tsv"
        if metadata_tsv.exists():
            # Check if there's an integration directory or if metadata has been processed
            integration_dir = work_dir / "integration"
            if integration_dir.exists() and any(integration_dir.iterdir()):
                return True, f"Integration directory exists: {integration_dir}"
        return False, None

    elif step_name == "quant":
        # Defer to amalgkit's internal skipping logic (redo: no)
        # We don't want to skip the whole step if only some specific samples are missing
        return False, None

    elif step_name == "merge":
        # Check for merged abundance file
        merge_out = step_params.get("out")
        if merge_out:
            merge_path = Path(merge_out)
            if merge_path.exists():
                return True, str(merge_path)

        merge_dir = Path(step_params.get("out_dir", work_dir / "merged"))
        merged_file = merge_dir / "merged_abundance.tsv"
        if merged_file.exists():
            return True, str(merged_file)
        return False, None

    elif step_name == "curate":
        curate_dir = Path(step_params.get("out_dir", work_dir / "curate"))
        if curate_dir.exists():
            # Check for expected curate output files
            expected_files = ["curated_abundance.tsv", "curated_metadata.tsv"]
            found_files = [f for f in expected_files if (curate_dir / f).exists()]
            if found_files:
                return True, f"Curate outputs found: {', '.join(found_files)} in {curate_dir}"
        return False, None

    elif step_name == "sanity":
        sanity_file = work_dir / "sanity_check.txt"
        if sanity_file.exists():
            return True, str(sanity_file)
        return False, None

    return False, None


def _execute_streaming_mode(
    config: AmalgkitWorkflowConfig,
    steps_remaining: List[Tuple[str, Any]],
    metadata_file: Path,
    chunk_size: int,
    step_functions: Dict[str, Any],
    check: bool = False,
    walk: bool = False,
    progress: bool = True,
    show_commands: bool = False,
) -> WorkflowExecutionResult:
    """Execute remaining steps in streaming mode (chunked by sample).

    This handles the getfastq -> integrate -> quant loop for small batches of samples
    to preserve disk space.
    """
    import csv

    logger.info(f"ENTERING STREAMING MODE: Processing remaining steps in chunks of {chunk_size} samples")

    # Identify which steps are per-chunk (getfastq, integrate, quant) and which are post-process (merge, etc)
    chunk_steps = []
    post_process_steps = []

    for name, params in steps_remaining:
        if name in ("getfastq", "integrate", "quant"):
            chunk_steps.append((name, params))
        else:
            post_process_steps.append((name, params))

    logger.info(f"Chunk loop steps: {[s[0] for s in chunk_steps]}")
    logger.info(f"Post-process steps: {[s[0] for s in post_process_steps]}")

    # Read samples
    with open(metadata_file, "r") as f:
        reader = csv.DictReader(f, delimiter="\t")
        fieldnames = reader.fieldnames
        samples = list(reader)

    total_samples = len(samples)
    num_chunks = math.ceil(total_samples / chunk_size)

    all_step_results = []
    success = True

    # Chunk loop
    for chunk_idx in range(num_chunks):
        start_idx = chunk_idx * chunk_size
        end_idx = min((chunk_idx + 1) * chunk_size, total_samples)
        chunk_samples = samples[start_idx:end_idx]

        logger.info(f"Processing chunk {chunk_idx+1}/{num_chunks} (Samples {start_idx+1}-{end_idx})")

        # Create chunk metadata
        chunk_meta_file = config.work_dir / "metadata" / f"metadata_chunk_{chunk_idx}.tsv"
        with open(chunk_meta_file, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
            writer.writeheader()
            writer.writerows(chunk_samples)

        # Run chunk steps
        for step_name, step_params in chunk_steps:
            logger.info(f"  > Chunk step: {step_name}")

            # Override params with chunk metadata
            chunk_params = step_params.copy()
            chunk_params["metadata"] = str(chunk_meta_file)

            # Sanitization and command building just for logging/dry run
            if show_commands or walk:
                from metainformant.rna.amalgkit.amalgkit import build_amalgkit_command

                sanitized = sanitize_params_for_cli(step_name, chunk_params)
                cmd = build_amalgkit_command(step_name, sanitized)
                if show_commands:
                    logger.info(f"  Command: {' '.join(cmd)}")

            if walk:
                all_step_results.append(
                    WorkflowStepResult(step_name=f"{step_name}_chunk{chunk_idx}", return_code=0, success=True)
                )
                continue

            # Execute
            try:
                # Add check/fail logic for disk space within the chunk loop too
                if step_name == "getfastq":
                    _check_disk_space_or_fail(
                        config.work_dir, min_free_gb=5.0, step_name=f"{step_name} (chunk {chunk_idx})"
                    )

                step_func = step_functions.get(step_name)
                # Call step function (returns CompletedProcess)
                # Note: step_func in amalgkit.py handles CLI arg building internally from dict
                result = step_func(chunk_params)

                step_res = WorkflowStepResult(
                    step_name=f"{step_name}_chunk{chunk_idx}",
                    return_code=result.returncode,
                    success=result.returncode == 0,
                    error_message=result.stderr if result.returncode != 0 else None,
                )
                all_step_results.append(step_res)

                if result.returncode != 0:
                    logger.error(f"  Step {step_name} failed for chunk {chunk_idx}")
                    success = False
                    if check:
                        return WorkflowExecutionResult(all_step_results, False, len(all_step_results), 0, 1)
                    # Break chunk loop to prevents cascading failures (e.g. running quant on empty inputs)
                    break
            except Exception as e:
                logger.error(f"  Exception in {step_name} chunk {chunk_idx}: {e}")
                all_step_results.append(WorkflowStepResult(f"{step_name}_chunk{chunk_idx}", 1, False, str(e)))
                success = False

        # Aggressive cleanup after chunk
        # Delete fastq files for this chunk to free space
        # Note: 'quant' usually handles this if keep_fastq: no is set.
        # But we force it here to be safe
        logger.info(f"  Cleaning up chunk {chunk_idx} FASTQ files...")
        chunk_ids = [row.get("run", "") for row in chunk_samples]
        cleanup_fastqs(config, chunk_ids)

    # After all chunks, run post-process steps (merge, etc)
    if success or not check:
        logger.info("Streaming chunks completed. processing post-chunk steps...")
        for step_name, step_params in post_process_steps:
            # Logic same as main loop...
            # We should probably recurse or just execute them here
            # For simplicity, let's just execute them using the main metadata (selected)
            # because merge needs ALL samples.

            logger.info(f"Step: {step_name}")
            if walk:
                all_step_results.append(WorkflowStepResult(step_name, 0, True))
                continue

            try:
                step_func = step_functions.get(step_name)
                result = step_func(step_params)  # Use original params with full metadata
                all_step_results.append(WorkflowStepResult(step_name, result.returncode, result.returncode == 0))
            except Exception as e:
                all_step_results.append(WorkflowStepResult(step_name, 1, False, str(e)))

    return WorkflowExecutionResult(
        steps_executed=all_step_results,
        success=success and all(s.return_code == 0 for s in all_step_results),
        total_steps=len(all_step_results),
        successful_steps=sum(1 for s in all_step_results if s.success),
        failed_steps=sum(1 for s in all_step_results if not s.success),
    )


def execute_workflow(
    config: AmalgkitWorkflowConfig,
    steps: Optional[List[str]] = None,
    check: bool = False,
    walk: bool = False,
    progress: bool = True,
    show_commands: bool = False,
    **kwargs: Any,
) -> WorkflowExecutionResult:
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
    from metainformant.rna.amalgkit.amalgkit import (
        AmalgkitParams,
    )
    from metainformant.rna.amalgkit.amalgkit import config as amalgkit_config
    from metainformant.rna.amalgkit.amalgkit import (
        csca,
        cstmm,
        curate,
        getfastq,
        integrate,
        merge,
        metadata,
        quant,
        sanity,
        select,
    )
    from metainformant.rna.core.deps import check_step_dependencies

    # Check if we need to prepare reference genome (for quant step)
    has_quant_or_merge = (not steps) or any(s in ("quant", "merge") for s in steps)
    if not steps or has_quant_or_merge:
        # If running steps is None (all) or includes quant/merge, ensure genome is ready
        prepare_reference_genome(config)

    logger.info(f"Starting amalgkit workflow for species: {config.species_list}")
    logger.info(f"Working directory: {config.work_dir}")

    # Create working directory
    config.work_dir.mkdir(parents=True, exist_ok=True)

    steps_planned = plan_workflow(config)
    steps_config = config.extra_config.get("steps", {}) or config.extra_config.get("per_step", {})

    # Filter steps if specific steps requested
    if steps:
        original_count = len(steps_planned)
        steps_planned = [(name, params) for name, params in steps_planned if name in steps]
        logger.info(
            f"Filtered workflow from {original_count} to {len(steps_planned)} steps: {[s[0] for s in steps_planned]}"
        )

    if not steps_planned:
        logger.warning("No steps to execute after filtering")
        return WorkflowExecutionResult(
            steps_executed=[], success=True, total_steps=0, successful_steps=0, failed_steps=0
        )

    # Run pre-flight checks
    if any(step == "getfastq" for step, _ in steps_planned):
        verify_getfastq_prerequisites(config, steps_planned)

    # Configure vdb-config and prepare getfastq step
    if any(step == "getfastq" for step, _ in steps_planned):
        steps_planned = setup_vdb_config(config, steps_planned)

    step_results = []
    step_functions = {
        "metadata": metadata,
        "integrate": integrate,
        "config": amalgkit_config,
        "select": select,
        "getfastq": getfastq,
        "quant": quant,
        "merge": merge,
        "cstmm": cstmm,
        "curate": curate,
        "csca": csca,
        "sanity": sanity,
    }

    # Pre-execution status reporting
    completed_steps, steps_to_run = check_step_completion_status(steps_planned, config, **kwargs)

    logger.info(f"Starting execution of {len(steps_planned)} steps")
    manifest_path = config.work_dir / "amalgkit.manifest.jsonl"

    for i, (step_name, step_params) in enumerate(steps_planned, 1):
        if progress:
            logger.info(f"Step {i}/{len(steps_planned)}: {step_name}")

        # Check dependencies
        deps_ok, deps_msg = check_step_dependencies(step_name, step_params, config)
        if not deps_ok:
            logger.warning(f"Skipping step {step_name}: Dependencies not satisfied: {deps_msg}")
            res = WorkflowStepResult(
                step_name=step_name,
                return_code=126,  # Magic code for dependency skip in some tests
                success=False,
                error_message=deps_msg,
            )
            step_results.append(res)
            # Record in manifest
            with open(manifest_path, "a") as f:
                f.write(
                    json.dumps(
                        {
                            "step": step_name,
                            "status": "skipped",
                            "return_code": 126,
                            "duration_seconds": 0,
                            "timestamp": time_mod.time(),
                            "note": f"Dependency skip: {deps_msg}",
                        }
                    )
                    + "\n"
                )
            continue

        if walk:
            logger.info(f"Would execute: {step_name}")
            step_results.append(
                WorkflowStepResult(step_name=step_name, return_code=0, success=True, command="(dry run)")
            )
            continue

        # Check if step is already completed
        is_completed, completion_indicator = _is_step_completed(step_name, step_params, config)

        # Check redo setting
        redo_value = step_params.get("redo", "no")
        if isinstance(redo_value, bool):
            force_redo = redo_value
        else:
            force_redo = str(redo_value).lower() in ("yes", "true", "1")

        if is_completed and not force_redo:
            logger.info(f"Step {i}/{len(steps_planned)}: {step_name} - Already completed, skipping")
            logger.info(f"  Completion indicator: {completion_indicator}")

            # Record skipped step in manifest
            with open(manifest_path, "a") as f:
                f.write(
                    json.dumps(
                        {
                            "step": step_name,
                            "status": "skipped",
                            "return_code": 0,
                            "duration_seconds": 0,
                            "timestamp": time_mod.time(),
                        }
                    )
                    + "\n"
                )

            # For config step, ensure symlink exists even if step was skipped
            if step_name == "config":
                config_base_dir = config.work_dir / "config_base"
                config_dir = config.work_dir / "config"
                if config_base_dir.exists():
                    if config_dir.exists() or config_dir.is_symlink():
                        # Symlink or directory already exists - verify it points to config_base
                        if config_dir.is_symlink():
                            target = config_dir.readlink()
                            if target == config_base_dir or target.name == "config_base":
                                logger.debug(f"Config symlink already exists and is correct: {config_dir}")
                            else:
                                logger.warning(f"Config symlink exists but points to wrong target: {target}")
                        else:
                            logger.debug(f"Config directory already exists: {config_dir}")
                    else:
                        try:
                            # Use absolute path for symlink to ensure it resolves correctly
                            config_dir.symlink_to(config_base_dir.resolve())
                            logger.info(
                                f"Created symlink: {config_dir} -> {config_base_dir.resolve()} (for select step compatibility)"
                            )
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

                                    for config_file in config_base_dir.glob("*.config"):
                                        shutil.copy2(config_file, config_dir / config_file.name)
                                    logger.info(
                                        f"Copied config files from {config_base_dir} to {config_dir} as fallback"
                                    )
                                except Exception as e2:
                                    logger.warning(f"Could not copy config files as fallback: {e2}")

            step_results.append(
                WorkflowStepResult(
                    step_name=step_name, return_code=0, success=True, command="(skipped - already completed)"
                )
            )
            continue

        # Quant/Integrate path fix: Ensure getfastq directory is available in work_dir
        if step_name in ("quant", "integrate", "merge") and any(s[0] == "getfastq" for s in steps_planned):
            # Find getfastq output dir
            getfastq_params = next((p for n, p in steps_planned if n == "getfastq"), None)
            if getfastq_params:
                raw_gf_out = Path(getfastq_params.get("out_dir", config.work_dir / "fastq"))
                real_gf_dir = raw_gf_out / "getfastq" if raw_gf_out.name != "getfastq" else raw_gf_out

                target_link = config.work_dir / "getfastq"

                if real_gf_dir.exists() and not target_link.exists() and real_gf_dir != target_link:
                    try:
                        if target_link.is_symlink():
                            target_link.unlink()
                        target_link.symlink_to(real_gf_dir)
                        logger.info(f"Symlinked getfastq dir for compatibility: {real_gf_dir} -> {target_link}")
                    except Exception as e:
                        logger.warning(f"Could not symlink getfastq dir: {e}")

        try:
            # Late-binding metadata path: if metadata_selected.tsv was created after planning, use it
            if step_name in ("getfastq", "integrate", "merge", "quant", "curate", "sanity"):
                selected_metadata = (config.work_dir / "metadata" / "metadata_selected.tsv").absolute()
                planned_metadata = step_params.get("metadata")
                if selected_metadata.exists() and (
                    not planned_metadata or "metadata_selected.tsv" not in str(planned_metadata)
                ):
                    # For quant, only override if integrated metadata doesn't exist yet
                    if step_name == "quant":
                        integrated_metadata = config.work_dir / "metadata" / "metadata_updated_for_private_fastq.tsv"
                        if not integrated_metadata.exists():
                            step_params["metadata"] = str(selected_metadata)
                            logger.info(
                                f"Dynamically updated {step_name} metadata to filtered version: {selected_metadata}"
                            )
                    else:
                        step_params["metadata"] = str(selected_metadata)
                        logger.info(
                            f"Dynamically updated {step_name} metadata to filtered version: {selected_metadata}"
                        )

            # Get step function
            step_func = step_functions.get(step_name)
            if not step_func:
                error_msg = f"Unknown step: {step_name}"
                logger.error(error_msg)
                step_results.append(
                    WorkflowStepResult(step_name=step_name, return_code=1, success=False, error_message=error_msg)
                )
                if check:
                    break
                continue

            # Execute step (pass params dict directly, not AmalgkitParams)
            command_str = None
            if show_commands:
                from metainformant.rna.amalgkit.amalgkit import build_amalgkit_command

                sanitized_params = sanitize_params_for_cli(step_name, step_params)
                command = build_amalgkit_command(step_name, sanitized_params)
                command_str = " ".join(command)
                logger.info(f"Command: {command_str}")

            # Ensure config symlink exists before select step
            if step_name == "select":
                config_base_dir = config.work_dir / "config_base"
                config_dir = config.work_dir / "config"
                if config_base_dir.exists() and not (config_dir.exists() or config_dir.is_symlink()):
                    try:
                        config_dir.symlink_to(config_base_dir.resolve())
                        logger.info(
                            f"Created config symlink before select step: {config_dir} -> {config_base_dir.resolve()}"
                        )
                    except Exception as e:
                        logger.warning(f"Could not create config symlink: {e}")
                        # Try copying as fallback
                        try:
                            config_dir.mkdir(parents=True, exist_ok=True)
                            import shutil

                            for config_file in config_base_dir.glob("*.config"):
                                shutil.copy2(config_file, config_dir / config_file.name)
                            logger.info(f"Copied config files as fallback")
                        except Exception as e2:
                            logger.warning(f"Could not copy config files: {e2}")
                # Update step_params to use absolute path to config directory
                # Always prefer config/ (symlink) over config_base/ for select step compatibility
                # Use absolute path of config_dir itself (not resolved) so amalgkit uses the symlink
                if (config_dir.exists() or config_dir.is_symlink()) and config_base_dir.exists():
                    step_params["config_dir"] = str(config_dir.absolute())
                    logger.debug(f"Using config directory (symlink): {config_dir.absolute()}")
                elif config_base_dir.exists():
                    # Fallback: use config_base directly if symlink doesn't exist
                    step_params["config_dir"] = str(config_base_dir.absolute())
                    logger.debug(f"Using config_base directory (fallback): {config_base_dir.absolute()}")

            # STREAMING MODE DETECTION
            # If we are about to run getfastq, check if we should switch to streaming mode
            # This happens if we have many samples and limited disk space
            if step_name == "getfastq" and any(s[0] == "quant" for s in steps_planned):
                # We have getfastq AND quant in the plan - candidate for streaming

                # Check remaining samples
                unquantified_metadata = config.work_dir / "metadata" / "metadata_unquantified.tsv"
                if not unquantified_metadata.exists():
                    # Try to create it if selected metadata exists
                    selected_metadata = config.work_dir / "metadata" / "metadata_selected.tsv"
                    if selected_metadata.exists():
                        filter_metadata_for_unquantified(config, selected_metadata, unquantified_metadata)

                if unquantified_metadata.exists():
                    # Count samples
                    with open(unquantified_metadata, "r") as f:
                        reader = csv.DictReader(f, delimiter="\t")
                        remaining_count = sum(1 for _ in reader)

                    # Check disk space
                    ok, free_gb = _check_disk_space(config.work_dir, min_free_gb=10.0)
                    estimated_need = remaining_count * 3.0  # 3GB per sample conservative estimate

                    logger.debug(
                        f"Streaming detection - Remaining: {remaining_count}, Free: {free_gb:.2f}GB, Need: {estimated_need:.2f}GB"
                    )

                    if remaining_count > 5 and estimated_need > (free_gb * 0.8):
                        logger.warning(
                            f"Insufficient disk space for all-at-once processing: "
                            f"Need ~{estimated_need:.1f}GB, have {free_gb:.1f}GB. "
                            f"Switching to STREAMING MODE (chunk size 5)."
                        )

                        # EXECUTE STREAMING MODE
                        # This processes the remaining steps (getfastq, integrate, quant) in chunks
                        # Returns the combined result
                        streaming_result = _execute_streaming_mode(
                            config,
                            steps_planned[i - 1 :],
                            unquantified_metadata,
                            chunk_size=5,
                            step_functions=step_functions,
                            check=check,
                            walk=walk,
                            progress=progress,
                            show_commands=show_commands,
                        )

                        # Merge results and return complete workflow result
                        # Note: we need to merge with any steps already executed (steps 0 to i-1)
                        # But execute_workflow returns WorkflowExecutionResult which aggregates.
                        # We'll just return the streaming result combined with what we have?
                        # Actually easier to just break here and return combined list
                        combined_steps = step_results + streaming_result.steps_executed

                        return WorkflowExecutionResult(
                            steps_executed=combined_steps,
                            success=streaming_result.success and all(s.return_code == 0 for s in step_results),
                            total_steps=len(combined_steps),
                            successful_steps=sum(1 for s in combined_steps if s.return_code == 0),
                            failed_steps=sum(1 for s in combined_steps if s.return_code != 0),
                        )

            # Validate filtered metadata exists for steps that need it
            if step_name in ("getfastq", "integrate", "merge"):
                filtered_metadata = (config.work_dir / "metadata" / "metadata_selected.tsv").absolute()
                if not filtered_metadata.exists():
                    error_msg = f"Filtered metadata not found: {filtered_metadata}. Run 'select' step first."
                    logger.error(error_msg)
                    step_results.append(
                        WorkflowStepResult(step_name=step_name, return_code=1, success=False, error_message=error_msg)
                    )
                    if check:
                        break
                    continue

            # Pre-step prerequisite validation
            prereq_error = validate_step_prerequisites(step_name, step_params, config, steps_planned, steps_config)
            if prereq_error:
                logger.error(prereq_error)
                step_results.append(
                    WorkflowStepResult(step_name=step_name, return_code=1, success=False, error_message=prereq_error)
                )
                if check:
                    break
                continue

            # Sanitize parameters before passing to step function
            sanitized_params = sanitize_params_for_cli(step_name, step_params)
            # Pass monitoring/progress preferences through to amalgkit wrappers.
            # `show_progress` controls tqdm progress bars (if available) and heartbeat cadence.
            result = step_func(
                sanitized_params,
                show_progress=progress,
                heartbeat_interval=5,
                check=check,
            )

            # Create step result
            error_message = None
            if result.returncode != 0:
                # Check if step actually produced expected outputs despite error
                is_completed, completion_indicator = _is_step_completed(step_name, step_params, config)
                if is_completed:
                    logger.warning(
                        f"Step {step_name} reported error (return code {result.returncode}) but outputs exist: {completion_indicator}"
                    )
                    logger.warning(f"Continuing workflow - step appears to have completed successfully despite error")
                    # Mark as successful despite error
                    result.returncode = 0
                    logger.info(f"Step {step_name} marked as successful based on output validation")
                elif step_name == "integrate":
                    # FALLBACK: If integrate fails (e.g. KeyError), try manual fallback
                    logger.warning(f"Step integrate reported error. Attempting manual integration fallback...")
                    try:
                        if manual_integration_fallback(config):
                            logger.info("Manual integration fallback successful")
                            result.returncode = 0
                        else:
                            logger.error("Manual integration fallback failed")
                    except Exception as e:
                        logger.error(f"Manual integration fallback exception: {e}")

                if result.returncode != 0:
                    error_message = f"Step failed with return code {result.returncode}"
                    if result.stderr:
                        error_message += f": {result.stderr}"
                    logger.error(f"Step {step_name} failed with return code {result.returncode}")
                    if result.stderr:
                        logger.error(f"Error output: {result.stderr}")
            else:
                logger.info(f"Step {step_name} completed successfully")

            # Run post-step actions (metadata dedup, symlinks, validation, cleanup)
            handle_post_step_actions(
                step_name,
                step_params,
                result,
                config,
                steps_planned,
                steps_config,
                check,
                step_results,
            )

            step_results.append(
                WorkflowStepResult(
                    step_name=step_name,
                    return_code=result.returncode,
                    success=result.returncode == 0,
                    error_message=error_message,
                    command=command_str,
                )
            )

            if result.returncode != 0 and check:
                # Bubble up the exception if check=True
                raise subprocess.CalledProcessError(
                    result.returncode,
                    command_str or step_name,
                    output=getattr(result, "stdout", ""),
                    stderr=getattr(result, "stderr", ""),
                )

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
                step_results.append(
                    WorkflowStepResult(
                        step_name=step_name,
                        return_code=0,
                        success=True,
                        error_message=f"Exception occurred but outputs validated: {error_msg}",
                        command=command_str,
                    )
                )
                logger.info(f"Step {step_name} marked as successful based on output validation")
                # Continue to next step (don't break workflow)
                continue
            else:
                step_results.append(
                    WorkflowStepResult(step_name=step_name, return_code=1, success=False, error_message=error_msg)
                )
                if check:
                    # Re-raise the original exception if check=True
                    raise e
                continue

    successful_steps = sum(1 for sr in step_results if sr.success)
    total_steps = len(step_results)
    failed_steps = total_steps - successful_steps

    # Generate workflow summary with remediation steps
    if progress:
        log_workflow_summary(step_results, config)

    # After loop, ensure we record results in manifest
    # This ensures consistency with tests expecting a manifest file
    try:
        manifest_path = config.manifest_path
        with open(manifest_path, "a") as f:
            for res in step_results:
                f.write(
                    json.dumps(
                        {
                            "step": res.step_name,
                            "return_code": res.return_code,
                            "success": res.success,
                            "duration_seconds": getattr(res, "duration_seconds", 0),
                            "timestamp": time_mod.time(),
                        }
                    )
                    + "\n"
                )
        logger.debug(f"Workflow manifest updated: {manifest_path}")
    except Exception as e:
        logger.warning(f"Could not write workflow manifest: {e}")

    return WorkflowExecutionResult(
        steps_executed=step_results,
        success=failed_steps == 0,
        total_steps=total_steps,
        successful_steps=successful_steps,
        failed_steps=failed_steps,
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

    return {"config": config.to_dict(), "return_codes": return_codes, "success": all(rc == 0 for rc in return_codes)}
