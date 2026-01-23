"""RNA-seq workflow orchestration and configuration management.

This module provides high-level functions for managing RNA-seq analysis workflows,
including configuration loading, workflow planning, and execution orchestration.
"""

from __future__ import annotations

import json
import shutil
import os
import math
import csv
import time as time_mod
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, Union
from dataclasses import dataclass

from metainformant.core import io, logging
from metainformant.core.utils.config import load_mapping_from_file
from metainformant.rna.amalgkit.metadata_filter import filter_selected_metadata
from metainformant.rna.amalgkit.metadata_utils import deduplicate_metadata

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
            print(f"DEBUG workflow step {step} params: {step_params}")

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
        total_size = sum(f.stat().st_size for f in tmp_dir.rglob("*") if f.is_file())
        size_gb = total_size / (1024**3)

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


def _check_disk_space(path: Path, min_free_gb: float = 10.0) -> Tuple[bool, float]:
    """Check if there's sufficient disk space.

    Args:
        path: Path to check disk space for
        min_free_gb: Minimum free space required in GB

    Returns:
        Tuple of (is_sufficient, free_gb)
    """
    try:
        stat = os.statvfs(path)
        free_gb = (stat.f_bavail * stat.f_frsize) / (1024**3)

        if free_gb < min_free_gb:
            logger.warning(f"Low disk space: {free_gb:.2f} GB free (minimum: {min_free_gb} GB) at {path}")
            return False, free_gb
        else:
            logger.debug(f"Disk space check: {free_gb:.2f} GB free at {path}")
            return True, free_gb
    except Exception as e:
        logger.warning(f"Could not check disk space: {e}")
        return True, 100.0  # Assume OK if we can't check


def _check_disk_space_or_fail(path: Path, min_free_gb: float = 5.0, step_name: str = "") -> float:
    ok, free_gb = _check_disk_space(path, min_free_gb)
    if not ok and free_gb >= 0:
        raise RuntimeError(
            f"CRITICAL: Disk space too low to continue {step_name}. "
            f"Only {free_gb:.2f} GB free (need {min_free_gb} GB minimum). "
            f"Free up disk space and re-run with --steps {step_name} to resume."
        )
    return free_gb


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


def cleanup_fastqs(config, sample_ids):
    """Delete FASTQ files for specific samples."""
    # Determine configured getfastq output directory
    getfastq_conf_dir = None
    if config.per_step and "getfastq" in config.per_step:
        gf_out = config.per_step["getfastq"].get("out_dir")
        if gf_out:
            getfastq_conf_dir = Path(gf_out)

    for sample_id in sample_ids:
        if not sample_id:
            continue
        # Possible locations
        paths = [
            config.work_dir / "fastq" / "getfastq" / sample_id,
            config.work_dir / "fastq" / sample_id,
            config.work_dir / "getfastq" / sample_id,
        ]

        # Add configured path if available (amalgkit adds 'getfastq' subdir)
        if getfastq_conf_dir:
            paths.append(getfastq_conf_dir / "getfastq" / sample_id)

        for p in paths:
            if p.exists() and p.is_dir():
                import shutil

                try:
                    shutil.rmtree(p, ignore_errors=True)
                except Exception as e:
                    pass  # Best effort cleanup


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
        metadata,
        integrate,
        select,
        getfastq,
        quant,
        merge,
        cstmm,
        curate,
        csca,
        sanity,
    )
    from metainformant.rna.amalgkit.amalgkit import config as amalgkit_config
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

    # Configure vdb-config repository path for getfastq step (if getfastq is in workflow)
    # prefetch downloads SRA files to vdb-config repository root, so we need to set it to
    # the amalgkit fastq output directory where amalgkit expects the files
    if any(step == "getfastq" for step, _ in steps_planned):
        try:
            import subprocess

            # Find getfastq step parameters to get out_dir
            getfastq_params = next((params for step, params in steps_planned if step == "getfastq"), None)
            if getfastq_params:
                # Get fastq output directory from getfastq step params
                fastq_out_dir = getfastq_params.get("out_dir", str(config.work_dir / "fastq"))
                fastq_dir = Path(fastq_out_dir)
                # Amalgkit creates getfastq subdirectory, so prefetch should download to that
                getfastq_dir = fastq_dir / "getfastq" if fastq_dir.name != "getfastq" else fastq_dir
                getfastq_dir.mkdir(parents=True, exist_ok=True)

                # CRITICAL: Check disk space before proceeding - fail if too low
                free_gb = _check_disk_space_or_fail(getfastq_dir, min_free_gb=5.0, step_name="getfastq")
                if free_gb < 20.0:
                    logger.warning(f"Low disk space ({free_gb:.1f}GB) - workflow may fail during downloads")

                # Clean up any incorrectly placed SRA files before configuring vdb-config
                _cleanup_incorrectly_placed_sra_files(getfastq_dir)

                # Clean up temp files if needed
                repo_root = Path(__file__).resolve().parent.parent.parent.parent
                tmp_dir = repo_root / ".tmp" / "fasterq-dump"
                _cleanup_temp_files(tmp_dir, max_size_gb=50.0)

                # Set vdb-config repository root to getfastq directory (where prefetch downloads SRA files)
                result = subprocess.run(
                    ["vdb-config", "-s", f"/repository/user/main/public/root={getfastq_dir}"],
                    capture_output=True,
                    text=True,
                    timeout=5,
                )
                if result.returncode == 0:
                    logger.info(f"Configured vdb-config repository path for prefetch: {getfastq_dir}")
                else:
                    logger.warning(
                        f"Could not set vdb-config repository path (may require interactive): {result.stderr[:100]}"
                    )
                    logger.info(f"Will rely on TMPDIR and VDB_CONFIG environment variables")

                # CRITICAL: Prepare directories for extraction by cleaning empty sample dirs
                # This fixes the issue where amalgkit skips extraction if sample directories exist
                cleaned = prepare_extraction_directories(getfastq_dir)
                if cleaned > 0:
                    logger.info(f"Prepared {cleaned} sample directories for extraction")

                # Create filtered metadata with only samples that have SRA files available
                # This prevents amalgkit from crashing on network failures
                source_metadata = config.work_dir / "metadata" / "metadata_selected.tsv"
                extraction_metadata = config.work_dir / "metadata" / "metadata_extraction.tsv"

                # SKIP LOGIC: Filter out samples that are already quantified
                # This prevents re-downloading samples that have already been processed
                unquantified_metadata = config.work_dir / "metadata" / "metadata_unquantified.tsv"
                if source_metadata.exists():
                    remaining = filter_metadata_for_unquantified(config, source_metadata, unquantified_metadata)
                    if remaining == 0:
                        # All samples already quantified - skip getfastq entirely
                        logger.info("All samples already quantified - skipping getfastq step")
                        # Remove getfastq from steps_planned
                        steps_planned = [(name, params) for name, params in steps_planned if name != "getfastq"]
                    elif remaining < sum(1 for _ in open(source_metadata)) - 1:
                        # Some samples already quantified - use filtered metadata
                        source_metadata = unquantified_metadata
                        logger.info(f"Using unquantified metadata ({remaining} samples remaining)")

                if source_metadata.exists() and any(s[0] == "getfastq" for s in steps_planned):
                    num_samples = create_extraction_metadata(getfastq_dir, source_metadata, extraction_metadata)
                    if num_samples > 0:
                        # Update the getfastq step to use the extraction-only metadata
                        for i, (step_name, step_params) in enumerate(steps_planned):
                            if step_name == "getfastq":
                                steps_planned[i] = (step_name, {**step_params, "metadata": str(extraction_metadata)})
                                logger.info(f"Using extraction metadata ({num_samples} samples) for getfastq step")
            else:
                logger.warning("Could not find getfastq step parameters, skipping vdb-config setup")
        except Exception as e:
            logger.warning(f"Could not configure vdb-config: {e}")
            logger.info("Will rely on TMPDIR and VDB_CONFIG environment variables")

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

    # Pre-execution status reporting: check which steps are already completed
    completed_steps = []
    steps_to_run = []
    steps_to_skip = []

    for step_name, step_params in steps_planned:
        is_completed, completion_indicator = _is_step_completed(step_name, step_params, config)

        # Check redo setting (override from kwargs if present)
        redo_value = kwargs.get("redo", step_params.get("redo", "no"))
        if isinstance(redo_value, bool):
            force_redo = redo_value
        else:
            force_redo = str(redo_value).lower() in ("yes", "true", "1")

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

                    logger.info(
                        f"DEBUG: Streaming detection - Remaining: {remaining_count}, Free: {free_gb:.2f}GB, Need: {estimated_need:.2f}GB"
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
            if step_name == "integrate":
                # Check if FASTQ files exist before integrate
                # Check multiple possible locations for FASTQ files
                fastq_locations = []

                # 1. From config (if absolute or relative to CWD)
                fastq_dir_raw = steps_config.get("getfastq", {}).get("out_dir", config.work_dir / "fastq")
                fastq_locations.append(Path(fastq_dir_raw))

                # 2. Relative to work_dir parent (common layout)
                fastq_locations.append(config.work_dir.parent / "fastq")

                # 3. Inside work_dir (legacy/default layout)
                fastq_locations.append(config.work_dir)
                fastq_locations.append(config.work_dir / "fastq")

                fastq_files = []
                checked_dirs = []

                for f_dir in fastq_locations:
                    if not f_dir.exists():
                        continue

                    # Try both root and getfastq/ subdirectory
                    search_dirs = [f_dir]
                    if f_dir.name != "getfastq":
                        getfastq_subdir = f_dir / "getfastq"
                        if getfastq_subdir.exists():
                            search_dirs.append(getfastq_subdir)

                    for s_dir in search_dirs:
                        checked_dirs.append(str(s_dir))
                        found = list(s_dir.glob("**/*.fastq*"))
                        if found:
                            fastq_files.extend(found)
                            fastq_dir = s_dir  # Select for logging

                if not fastq_files:
                    error_msg = (
                        f"PREREQUISITE CHECK FAILED: No FASTQ files found before integrate step.\n"
                        f"  - Checked locations: {', '.join(checked_dirs)}\n"
                        f"  - Checked for: *.fastq, *.fastq.gz, *.fq, *.fq.gz\n\n"
                        f"REMEDIATION:\n"
                        f"  1. Ensure getfastq step completed successfully\n"
                        f"  2. Check validation report: {config.work_dir / 'validation' / 'getfastq_validation.json'}\n"
                        f"  3. Re-run getfastq step if FASTQ files are missing\n"
                        f"  4. Verify amalgkit getfastq extracted files correctly"
                    )
                    logger.error(error_msg)
                    step_results.append(
                        WorkflowStepResult(step_name=step_name, return_code=1, success=False, error_message=error_msg)
                    )
                    if check:
                        break
                    continue

            elif step_name == "quant":
                # Check if integrate completed (quant needs integrated metadata)
                integrated_metadata = config.work_dir / "integration" / "integrated_metadata.json"
                # Also check for metadata.tsv that amalgkit integrate creates
                metadata_tsv = config.work_dir / "metadata" / "metadata.tsv"
                if not integrated_metadata.exists() and not metadata_tsv.exists():
                    error_msg = (
                        f"PREREQUISITE CHECK FAILED: Integrate step must complete before quant.\n"
                        f"  - Expected: {integrated_metadata} or {metadata_tsv}\n\n"
                        f"REMEDIATION:\n"
                        f"  1. Ensure integrate step completed successfully\n"
                        f"  2. Re-run integrate step if needed\n"
                    )
                    logger.error(error_msg)
                    step_results.append(
                        WorkflowStepResult(step_name=step_name, return_code=1, success=False, error_message=error_msg)
                    )
                    if check:
                        break
                    continue

                # Check if FASTQ files exist (quant needs them for quantification)
                # Use same robust check for quant
                steps_config = config.extra_config.get("steps", {})
                fastq_locations = []
                fastq_locations.append(Path(steps_config.get("getfastq", {}).get("out_dir", config.work_dir / "fastq")))
                fastq_locations.append(config.work_dir.parent / "fastq")
                fastq_locations.append(config.work_dir)
                fastq_locations.append(config.work_dir / "fastq")

                fastq_files = []
                checked_dirs = []
                for f_dir in fastq_locations:
                    if not f_dir.exists():
                        continue
                    search_dirs = [f_dir]
                    if f_dir.name != "getfastq":
                        getfastq_subdir = f_dir / "getfastq"
                        if getfastq_subdir.exists():
                            search_dirs.append(getfastq_subdir)
                    for s_dir in search_dirs:
                        checked_dirs.append(str(s_dir))
                        found = list(s_dir.glob("**/*.fastq*"))
                        if found:
                            fastq_files.extend(found)
                            fastq_dir = s_dir

                if not fastq_files:
                    error_msg = (
                        f"PREREQUISITE CHECK FAILED: No FASTQ files found before quant step.\n"
                        f"  - Checked locations: {', '.join(checked_dirs)}\n\n"
                        f"REMEDIATION:\n"
                        f"  1. Ensure getfastq and integrate steps completed successfully\n"
                        f"  2. Check validation reports in: {config.work_dir / 'validation'}\n"
                        f"  3. Re-run getfastq step if FASTQ files are missing"
                    )
                    logger.error(error_msg)
                    step_results.append(
                        WorkflowStepResult(step_name=step_name, return_code=1, success=False, error_message=error_msg)
                    )
                    if check:
                        break
                    continue

                # Check quantification tools availability
                try:
                    from metainformant.rna.core.deps import check_quantification_tools

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
                        step_results.append(
                            WorkflowStepResult(
                                step_name=step_name, return_code=1, success=False, error_message=error_msg
                            )
                        )
                        if check:
                            break
                        continue
                except Exception as e:
                    logger.warning(f"Could not check quantification tools: {e}")

            elif step_name == "merge":
                # Check if quant files exist before merge
                steps_config = config.extra_config.get("steps", {})
                quant_dir_raw = steps_config.get("quant", {}).get("out_dir", config.work_dir / "quant")
                quant_dir = Path(quant_dir_raw)
                # Resolve to the actual 'quant' directory where files are stored (amalgkit creates a 'quant' subdir)
                if quant_dir.name != "quant":
                    actual_quant_dir = quant_dir / "quant"
                    if actual_quant_dir.exists():
                        quant_dir = actual_quant_dir

                # Check for quant output files (abundance.tsv, *_abundance.tsv, or quant.sf)
                # amalgkit quant can name them either way
                quant_files = (
                    list(quant_dir.glob("**/abundance.tsv"))
                    + list(quant_dir.glob("**/*_abundance.tsv"))
                    + list(quant_dir.glob("**/quant.sf"))
                )
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
                    step_results.append(
                        WorkflowStepResult(step_name=step_name, return_code=1, success=False, error_message=error_msg)
                    )
                    if check:
                        break
                    continue

                # Check R dependencies for merge (mandatory per user request)
                try:
                    from metainformant.rna.core.environment import check_rscript

                    r_available, r_message = check_rscript()
                    if not r_available:
                        error_msg = (
                            f"PREREQUISITE CHECK FAILED: R/Rscript not available for merge step.\n"
                            f"  Status: {r_message}\n\n"
                            f"REMEDIATION:\n"
                            f"  1. Install R: brew install r (on Mac) or apt-get install r-base-core\n"
                        )
                        logger.error(error_msg)
                        step_results.append(
                            WorkflowStepResult(
                                step_name=step_name, return_code=1, success=False, error_message=error_msg
                            )
                        )
                        if check:
                            break
                        continue
                    else:
                        # Check if ggplot2 package is available
                        try:
                            import subprocess

                            check_ggplot2 = subprocess.run(
                                ["Rscript", "-e", "library(ggplot2)"], capture_output=True, text=True, timeout=10
                            )
                            if check_ggplot2.returncode != 0:
                                error_msg = (
                                    f"PREREQUISITE CHECK FAILED: R package 'ggplot2' not available for merge step.\n"
                                    f"  Installation: Rscript -e \"install.packages('ggplot2', repos='https://cloud.r-project.org')\"\n"
                                )
                                logger.error(error_msg)
                                step_results.append(
                                    WorkflowStepResult(
                                        step_name=step_name, return_code=1, success=False, error_message=error_msg
                                    )
                                )
                                if check:
                                    break
                                continue
                        except Exception as e:
                            logger.warning(f"Failed to check R packages: {e}")
                except Exception as e:
                    logger.debug(f"Could not check R dependencies: {e}")

                # Path mapping fix: amalgkit merge expects 'quant' directory in ITS out_dir.
                # If quant outputs were written to a different location (like work_dir/quant),
                # we link them into the merge out_dir so amalgkit can find them.
                merge_out_dir_path = Path(step_params.get("out_dir", config.work_dir / "merge"))
                merge_quant_link = merge_out_dir_path / "quant"

                if quant_dir.resolve() != merge_quant_link.resolve():
                    if not merge_quant_link.exists() or not list(merge_quant_link.glob("**/*.tsv")):
                        logger.info(f"Bridging quantification results for merge: {quant_dir} -> {merge_quant_link}")
                        try:
                            merge_quant_link.parent.mkdir(parents=True, exist_ok=True)
                            if merge_quant_link.exists() and not merge_quant_link.is_symlink():
                                # If it's an empty directory, remove it to replace with symlink
                                if merge_quant_link.is_dir() and not any(merge_quant_link.iterdir()):
                                    merge_quant_link.rmdir()

                            if not merge_quant_link.exists():
                                # Use absolute path for symlink
                                merge_quant_link.symlink_to(quant_dir.resolve())
                        except Exception as e:
                            logger.warning(f"Could not bridge quant results for merge: {e}")

            # Check R availability for R-dependent steps (mandatory per user request)
            if step_name in ["curate", "cstmm"]:
                import shutil

                if not shutil.which("Rscript"):
                    error_msg = f"PREREQUISITE CHECK FAILED: Rscript not found for {step_name} step."
                    logger.error(error_msg)
                    step_results.append(
                        WorkflowStepResult(step_name=step_name, return_code=1, success=False, error_message=error_msg)
                    )
                    if check:
                        break
                    continue

                # Bridging for curate: it expects a 'merge' subdirectory in its out_dir
                merge_params = steps_config.get("merge", {})
                merge_out_dir = Path(merge_params.get("out_dir", config.work_dir / "merged"))
                merge_results_dir = merge_out_dir / "merge"

                curate_out_dir = Path(step_params.get("out_dir", config.work_dir / "curate"))
                curate_merge_link = curate_out_dir / "merge"

                if merge_results_dir.exists() and merge_results_dir.resolve() != curate_merge_link.resolve():
                    if not curate_merge_link.exists():
                        logger.info(f"Bridging merge results for curate: {merge_results_dir} -> {curate_merge_link}")
                        try:
                            curate_merge_link.parent.mkdir(parents=True, exist_ok=True)
                            curate_merge_link.symlink_to(merge_results_dir.resolve())
                        except Exception as e:
                            logger.warning(f"Could not bridge merge results for curate: {e}")

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

            # After integrate step, deduplicate metadata to avoid quant failure
            # amalgkit integrate often duplicates rows when adding private FASTQ paths
            if step_name == "integrate" and result.returncode == 0:
                metadata_dir = config.work_dir / "metadata"
                # Check for common metadata filenames produced/used by integrate
                metadata_files = [
                    metadata_dir / "metadata_updated_for_private_fastq.tsv",
                    metadata_dir / "metadata.tsv",
                ]
                for m_file in metadata_files:
                    if m_file.exists():
                        logger.info(f"Deduplicating metadata after integrate: {m_file}")
                        deduplicate_metadata(m_file)

            # After config step, create symlink from config/ to config_base/ for select step compatibility
            if step_name == "config":
                config_base_dir = config.work_dir / "config_base"
                config_dir = config.work_dir / "config"
                if config_base_dir.exists() and not config_dir.exists():
                    try:
                        # Use absolute path for symlink to ensure it resolves correctly
                        config_dir.symlink_to(config_base_dir.resolve())
                        logger.info(
                            f"Created symlink: {config_dir} -> {config_base_dir.resolve()} (for select step compatibility)"
                        )
                    except Exception as e:
                        logger.warning(f"Could not create config symlink: {e}")
                        # Try creating a regular directory and copying files as fallback
                        try:
                            config_dir.mkdir(parents=True, exist_ok=True)
                            import shutil

                            for config_file in config_base_dir.glob("*.config"):
                                shutil.copy2(config_file, config_dir / config_file.name)
                            logger.info(f"Copied config files from {config_base_dir} to {config_dir} as fallback")
                        except Exception as e2:
                            logger.warning(f"Could not copy config files as fallback: {e2}")

            # After select step, create filtered metadata for downstream steps
            if step_name == "select" and result.returncode == 0:
                try:
                    logger.info("Creating filtered metadata for selected samples (excluding LITE files)")
                    try:
                        filter_selected_metadata(
                            config.work_dir / "metadata" / "metadata.tsv",
                            config.work_dir / "metadata" / "metadata_selected.tsv",
                            exclude_lite_files=True,  # Automatically filter out LITE SRA files
                        )
                    except ValueError as e:
                        if "No samples meet the filtering criteria" in str(e):
                            logger.warning(
                                "All selected samples are LITE files. Creating metadata_selected.tsv without LITE filtering for this run."
                            )
                            # Fall back to not excluding LITE files if all samples would be filtered out
                            filter_selected_metadata(
                                config.work_dir / "metadata" / "metadata.tsv",
                                config.work_dir / "metadata" / "metadata_selected.tsv",
                                exclude_lite_files=False,  # Don't exclude LITE files if all would be filtered
                            )
                        else:
                            raise
                except Exception as e:
                    logger.error(f"Failed to create filtered metadata: {e}")
                    if check:
                        step_results.append(
                            WorkflowStepResult(
                                step_name="metadata_filtering", return_code=1, success=False, error_message=str(e)
                            )
                        )
                        break

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

                # Add step summary for getfastq to show skipped vs downloaded files
                if step_name == "getfastq":
                    # Check both stdout and stderr (amalgkit may output to either)
                    output_text = (result.stdout or "") + (result.stderr or "")
                    if output_text:
                        _log_getfastq_summary(output_text, logger)

                    # Run validation after getfastq step
                    try:
                        from metainformant.rna.analysis.validation import validate_all_samples, save_validation_report

                        validation_result = validate_all_samples(config, stage="extraction")
                        validation_dir = config.work_dir / "validation"
                        validation_dir.mkdir(parents=True, exist_ok=True)
                        save_validation_report(validation_result, validation_dir / "getfastq_validation.json")

                        # Log validation summary
                        total = validation_result.get("total_samples", 0)
                        validated = validation_result.get("validated", 0)
                        failed = validation_result.get("failed", 0)
                        logger.info(
                            f"Validation after getfastq: {validated}/{total} samples have FASTQ files extracted"
                        )

                        # Early exit if critical failure: no FASTQ files extracted
                        if validated == 0 and failed > 0 and total > 0:
                            # FALLBACK: Check if SRA files exist and attempt direct extraction
                            # This handles the case where amalgkit skipped extraction due to redo:no
                            # but SRA files were successfully downloaded
                            sra_dir = config.work_dir / "fastq" / "getfastq"
                            # Check config for explicit output directory
                            steps_config = config.extra_config.get("steps", {})
                            if "getfastq" in steps_config and "out_dir" in steps_config["getfastq"]:
                                step_out = Path(steps_config["getfastq"]["out_dir"])
                                if step_out.name != "getfastq":
                                    sra_dir = step_out / "getfastq"
                                else:
                                    sra_dir = step_out
                            elif not sra_dir.exists():
                                # Try standard layout relative to work_dir parent
                                sra_dir = config.work_dir.parent / "fastq" / "getfastq"

                            if sra_dir.exists():
                                # Check for SRA files in root, sra/ subdir, or sample subdirs
                                sra_files_exist = list(sra_dir.glob("*.sra"))
                                sra_files_exist.extend(list(sra_dir.glob("sra/*.sra")))
                                sra_files_exist.extend(list(sra_dir.glob("*/*.sra")))

                                if sra_files_exist:
                                    logger.warning(
                                        f"Validation failed (0 FASTQ) but found {len(sra_files_exist)} SRA files."
                                    )
                                    logger.info("Triggering direct fallback extraction...")
                                    try:
                                        extracted = extract_sra_directly(config, sra_dir, sra_dir)

                                        if extracted > 0:
                                            logger.info(f"Fallback extraction recovered {extracted} samples")
                                            # Re-validate
                                            validation_result = validate_all_samples(config, stage="extraction")
                                            validated = validation_result.get("validated", 0)
                                            failed = validation_result.get("failed", 0)
                                            logger.info(f"Re-validation: {validated}/{total} samples ready")

                                            if validated > 0:
                                                logger.info("Fallback recovery successful - continuing workflow")
                                                continue
                                    except Exception as e:
                                        logger.error(f"Fallback extraction failed: {e}")

                            lite_check_msg = ""
                            try:
                                import pandas as pd

                                selected_meta = config.work_dir / "metadata" / "metadata_selected.tsv"
                                if selected_meta.exists():
                                    df = pd.read_csv(selected_meta, sep="\t")
                                    if "run" in df.columns:
                                        runs = df["run"].tolist()
                                        # Check AWS/GCP links for .lite indicators
                                        lite_runs = []
                                        for _, row in df.iterrows():
                                            aws_link = str(row.get("AWS_Link", ""))
                                            gcp_link = str(row.get("GCP_Link", ""))
                                            if (
                                                ".lite" in aws_link
                                                or ".lite" in gcp_link
                                                or ".sralite" in aws_link
                                                or ".sralite" in gcp_link
                                            ):
                                                lite_runs.append(row.get("run", "unknown"))
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
                            step_results.append(
                                WorkflowStepResult(
                                    step_name="getfastq_validation",
                                    return_code=1,
                                    success=False,
                                    error_message=error_msg,
                                )
                            )
                            # Stop workflow unless check=False (user wants to continue anyway)
                            if check:
                                break
                        elif failed > 0:
                            logger.warning(f"{failed} samples missing FASTQ files after getfastq step")
                    except Exception as e:
                        logger.warning(f"Validation after getfastq failed: {e}")

                # Run validation after quant step
                if step_name == "quant":
                    try:
                        from metainformant.rna.analysis.validation import validate_all_samples, save_validation_report

                        validation_result = validate_all_samples(config, stage="quantification")
                        validation_dir = config.work_dir / "validation"
                        validation_dir.mkdir(parents=True, exist_ok=True)
                        save_validation_report(validation_result, validation_dir / "quant_validation.json")

                        # Log validation summary
                        total = validation_result.get("total_samples", 0)
                        validated = validation_result.get("validated", 0)
                        failed = validation_result.get("failed", 0)
                        logger.info(f"Validation after quant: {validated}/{total} samples quantified")
                        if failed > 0:
                            logger.warning(f"{failed} samples missing quantification files after quant step")

                        # CLEANUP: Delete FASTQ/SRA files for successfully quantified samples
                        # This implements the workflow_resilience.md requirement:
                        # "Only remove raw SRA and intermediate FASTQ data once the quant step
                        # has successfully produced abundance files"
                        if validated > 0:
                            cleanup_config = config.extra_config.get("cleanup_after_quant", True)
                            if cleanup_config:
                                logger.info("Running post-quant cleanup to free disk space...")
                                try:
                                    cleanup_result = cleanup_after_quant(config)
                                    if cleanup_result["samples_cleaned"] > 0:
                                        freed_gb = cleanup_result["bytes_freed"] / 1024**3
                                        logger.info(f"Post-quant cleanup freed {freed_gb:.2f} GB")
                                except Exception as cleanup_err:
                                    logger.warning(f"Post-quant cleanup failed (non-fatal): {cleanup_err}")
                            else:
                                logger.info("Post-quant cleanup disabled in config")
                    except Exception as e:
                        logger.warning(f"Validation after quant failed: {e}")

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
                    for line in failed_step.error_message.split("\n"):
                        logger.info(f"     {line}")

            logger.info("\n" + "=" * 80)
            logger.info("REMEDIATION STEPS")
            logger.info("=" * 80)

            # Provide step-specific remediation guidance
            remediation_steps = []

            for failed_step in failed_step_results:
                step_name = failed_step.step_name

                if step_name == "getfastq" or step_name == "getfastq_validation":
                    remediation_steps.append(
                        f"\n  For {step_name}:\n"
                        f"    1. Check logs: {config.work_dir / 'logs' / 'getfastq.log'}\n"
                        f"    2. Check validation: {config.work_dir / 'validation' / 'getfastq_validation.json'}\n"
                        f"    3. Verify SRA files downloaded: ls -lh {config.work_dir / 'fastq' / 'getfastq'}\n"
                        f"    4. Check tool availability: which fasterq-dump which fastp\n"
                        f"    5. Re-run: python3 scripts/rna/run_workflow.py {config.work_dir.parent.name} --steps getfastq"
                    )
                elif step_name == "integrate":
                    remediation_steps.append(
                        f"\n  For {step_name}:\n"
                        f"    1. Ensure getfastq completed successfully\n"
                        f"    2. Check FASTQ files exist: find {config.work_dir / 'fastq'} -name '*.fastq*'\n"
                        f"    3. Re-run: python3 scripts/rna/run_workflow.py {config.work_dir.parent.name} --steps getfastq integrate"
                    )
                elif step_name == "quant":
                    remediation_steps.append(
                        f"\n  For {step_name}:\n"
                        f"    1. Ensure getfastq and integrate completed successfully\n"
                        f"    2. Check quantification tools: which kallisto which salmon\n"
                        f"    3. Check validation: {config.work_dir / 'validation' / 'quant_validation.json'}\n"
                        f"    4. Re-run: python3 scripts/rna/run_workflow.py {config.work_dir.parent.name} --steps quant"
                    )
                elif step_name == "merge":
                    remediation_steps.append(
                        f"\n  For {step_name}:\n"
                        f"    1. Ensure quant completed successfully\n"
                        f"    2. Check quant files: find {config.work_dir / 'quant'} -name 'abundance.tsv' -o -name 'quant.sf'\n"
                        f"    3. If R plotting failed, install: Rscript -e \"install.packages('ggplot2', repos='https://cloud.r-project.org')\"\n"
                        f"    4. Re-run: python3 scripts/rna/run_workflow.py {config.work_dir.parent.name} --steps merge"
                    )
                elif step_name == "curate":
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


def extract_sra_directly(config: AmalgkitWorkflowConfig, sra_dir: Path, output_dir: Path) -> int:
    """Manually extract SRA files using fasterq-dump when amalgkit fails.

    This acts as a fallback when amalgkit getfastq skips extraction due to
    existing SRA files (redo: no) or other internal checks.

    Args:
        config: Workflow configuration
        sra_dir: Directory containing .sra files
        output_dir: Directory to output .fastq.gz files

    Returns:
        Number of successfully extracted samples
    """
    import shutil
    import subprocess
    import concurrent.futures

    from metainformant.core.io.download import monitor_subprocess_directory_growth

    # Check if fasterq-dump is available
    fasterq_dump = shutil.which("fasterq-dump")
    if not fasterq_dump:
        logger.error("fasterq-dump not found in PATH - cannot run fallback extraction")
        return 0

    # Check for gzip/pigz
    gzip_cmd = shutil.which("pigz") or shutil.which("gzip")
    if not gzip_cmd:
        logger.error("gzip/pigz not found - cannot compress FASTQ output")
        return 0

    # Find SRA files (root, sra/, or sample subdirs)
    sra_files = list(sra_dir.glob("*.sra"))
    sra_files.extend(list(sra_dir.glob("sra/*.sra")))
    sra_files.extend(list(sra_dir.glob("*/*.sra")))
    # Deduplicate by path
    sra_files = list({str(f): f for f in sra_files}.values())
    if not sra_files:
        logger.warning(f"No SRA files found in {sra_dir} for fallback extraction")
        return 0

    logger.info(f"Attempting fallback extraction for {len(sra_files)} files using {fasterq_dump}...")

    # Configure path for fasterq-dump temp files
    # Use repo .tmp directory to avoid filling up system /tmp or /var
    repo_root = Path(__file__).resolve().parent.parent.parent.parent.parent.parent
    tmp_dir = repo_root / ".tmp" / "fasterq-dump"
    tmp_dir.mkdir(parents=True, exist_ok=True)

    extracted_count = 0

    def process_sra(sra_file):
        sample_id = sra_file.stem
        sample_out_dir = output_dir / sample_id
        sample_out_dir.mkdir(parents=True, exist_ok=True)

        # Check if already extracted
        if list(sample_out_dir.glob("*.fastq.gz")):
            return True

        try:
            # Run fasterq-dump
            # Note: fasterq-dump outputs to CWD or specified dir.
            # We output to sample dir directly.
            cmd = [
                fasterq_dump,
                "--split-3",
                "--threads",
                str(min(config.threads, 4)),  # Limit threads per job
                "--outdir",
                str(sample_out_dir),
                "--temp",
                str(tmp_dir),
                str(sra_file),
            ]

            logger.info(f"Extracting {sample_id} with fasterq-dump...")

            # Use monitor_subprocess_directory_growth for progress + heartbeat
            heartbeat_file = config.work_dir / "heartbeat" / f"fallback_{sample_id}.json"
            heartbeat_file.parent.mkdir(exist_ok=True, parents=True)

            process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)

            rc, _ = monitor_subprocess_directory_growth(
                process=process,
                watch_dir=sample_out_dir,
                heartbeat_path=heartbeat_file,
                desc=f"Extracting {sample_id}",
                heartbeat_interval=10,
            )

            if rc != 0:
                raise RuntimeError(f"fasterq-dump exited with {rc}")

            # Gzip output files
            fastqs = list(sample_out_dir.glob("*.fastq"))
            if not fastqs:
                logger.warning(f"No FASTQ exported for {sample_id}")
                return False

            for fq in fastqs:
                subprocess.run([gzip_cmd, "-f", str(fq)], check=True)

            return True

        except Exception as e:
            logger.error(f"Failed to extract {sample_id}: {e}")
            return False

    # Run in parallel
    max_workers = max(1, config.threads // 4)
    with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers) as executor:
        results = list(executor.map(process_sra, sra_files))

    extracted_count = sum(results)
    logger.info(f"Fallback extraction completed: {extracted_count}/{len(sra_files)} samples extracted")

    # Clean up temp dir
    try:
        shutil.rmtree(tmp_dir)
    except (OSError, PermissionError, FileNotFoundError):
        pass

    return extracted_count


def manual_integration_fallback(config: AmalgkitWorkflowConfig) -> bool:
    """Fallback for when 'amalgkit integrate' fails (e.g. path detection bugs).

    Manually exposes FASTQ files from getfastq subdirectories to the location
    expected by 'amalgkit quant'.
    """
    logger.info("Attempting Manual Integration Fallback...")

    steps_config = config.extra_config.get("steps", {})
    getfastq_dir_raw = steps_config.get("getfastq", {}).get("out_dir", config.work_dir / "fastq" / "getfastq")
    getfastq_dir = Path(getfastq_dir_raw)

    # Also check standard amalgkit structure if configured path empty
    if not getfastq_dir.exists():
        getfastq_dir = config.work_dir / "fastq" / "getfastq"

    if not getfastq_dir.exists():
        logger.error(f"Cannot perform manual integration: getfastq dir {getfastq_dir} not found")
        return False

    # Target directory for quant (expects out_dir/getfastq/sample/sample_1.fastq.gz)
    # If quant.out_dir is work_dir, then it looks in work_dir/getfastq/
    steps_config = config.extra_config.get("steps", {})
    quant_out_dir = Path(steps_config.get("quant", {}).get("out_dir", config.work_dir))
    quant_input_dir = quant_out_dir / "getfastq"
    quant_input_dir.mkdir(parents=True, exist_ok=True)

    found_any = False

    # Search for FASTQ files in getfastq_dir recursively
    for fastq_file in getfastq_dir.glob("**/*.fastq*"):
        if fastq_file.is_file():
            # Identify sample ID from parent folder or filename
            # Standard amalgkit: getfastq/SRR12345/SRR12345_1.fastq.gz
            # SRA fallback: getfastq/SRR12345_1.fastq.gz

            sample_id = None
            if fastq_file.parent.name.startswith(("SRR", "ERR", "DRR")):
                sample_id = fastq_file.parent.name
            else:
                # Extract from filename
                import re

                match = re.search(r"(SRR\d+|ERR\d+|DRR\d+)", fastq_file.name)
                if match:
                    sample_id = match.group(1)

            if sample_id:
                sample_dest_dir = quant_input_dir / sample_id
                sample_dest_dir.mkdir(parents=True, exist_ok=True)
                dest_path = sample_dest_dir / fastq_file.name

                if not dest_path.exists():
                    try:
                        dest_path.symlink_to(fastq_file.resolve())
                        found_any = True
                        logger.info(f"Linked {fastq_file.name} -> {sample_id}")
                    except Exception:
                        try:
                            shutil.copy2(fastq_file, dest_path)
                            found_any = True
                        except Exception as e:
                            logger.warning(f"Failed to link {fastq_file}: {e}")
                else:
                    found_any = True

    if found_any:
        # Create a dummy metadata if amalgkit integrate was bypassed?
        # Actually, quant might need the integrated metadata.tsv.
        # If it doesn't exist, we might need to copy the selected one.
        metadata_tsv = config.work_dir / "metadata" / "metadata.tsv"
        selected_metadata = config.work_dir / "metadata" / "metadata_selected.tsv"
        if not metadata_tsv.exists() and selected_metadata.exists():
            shutil.copy2(selected_metadata, metadata_tsv)
            logger.info("Created metadata.tsv from metadata_selected.tsv for quant step")

        logger.info("Manual integration completed: valid FASTQ files expose for quant.")
        return True

    return False


def get_quantified_samples(config: AmalgkitWorkflowConfig) -> set[str]:
    """Get set of sample IDs that already have successful quantification results.

    Args:
        config: Workflow configuration

    Returns:
        Set of sample IDs (e.g., 'SRR12345') that have abundance.tsv files
    """
    quantified = set()

    # Get quant output directory from config
    steps_config = config.extra_config.get("steps", {})
    quant_dir_raw = steps_config.get("quant", {}).get("out_dir", config.work_dir / "quant")
    quant_dir = Path(quant_dir_raw)

    if not quant_dir.exists():
        return quantified

    # Find all abundance.tsv or *_abundance.tsv files and extract sample IDs
    patterns = ["**/abundance.tsv", "**/*_abundance.tsv"]
    for pattern in patterns:
        for abundance_file in quant_dir.glob(pattern):
            # Sample ID is the parent directory name (e.g., quant/SRR12345/abundance.tsv)
            sample_id = abundance_file.parent.name
            if sample_id.startswith(("SRR", "ERR", "DRR")):
                # Verify the file has content (not empty/corrupt)
                if abundance_file.stat().st_size > 100:  # Has more than just header
                    quantified.add(sample_id)

    if quantified:
        logger.info(f"Found {len(quantified)} samples already quantified")

    return quantified


def cleanup_after_quant(config: AmalgkitWorkflowConfig, dry_run: bool = False) -> dict[str, Any]:
    """Delete FASTQ and SRA files for samples with successful quantification.

    This function implements the cleanup sequence from workflow_resilience.md:
    "Only remove raw SRA and intermediate FASTQ data once the quant step
    has successfully produced abundance files, confirming the integrity
    of the upstream extraction."

    Args:
        config: Workflow configuration
        dry_run: If True, only report what would be deleted without deleting

    Returns:
        Dictionary with cleanup statistics:
        {
            'samples_cleaned': int,
            'fastq_files_deleted': int,
            'sra_files_deleted': int,
            'bytes_freed': int,
            'errors': list[str]
        }
    """
    result = {"samples_cleaned": 0, "fastq_files_deleted": 0, "sra_files_deleted": 0, "bytes_freed": 0, "errors": []}

    # Get quantified samples
    quantified_samples = get_quantified_samples(config)
    if not quantified_samples:
        logger.info("No quantified samples found - skipping cleanup")
        return result

    # Get FASTQ directory
    steps_config = config.extra_config.get("steps", {})
    fastq_dir_raw = steps_config.get("getfastq", {}).get("out_dir", config.work_dir / "fastq")
    fastq_dir = Path(fastq_dir_raw)

    # Check for getfastq subdirectory (amalgkit structure)
    if fastq_dir.name != "getfastq":
        getfastq_subdir = fastq_dir / "getfastq"
        if getfastq_subdir.exists():
            fastq_dir = getfastq_subdir

    if not fastq_dir.exists():
        logger.info(f"FASTQ directory does not exist: {fastq_dir}")
        return result

    logger.info(f"Cleaning up FASTQ/SRA files for {len(quantified_samples)} quantified samples")
    if dry_run:
        logger.info("DRY RUN - no files will be deleted")

    for sample_id in quantified_samples:
        sample_dir = fastq_dir / sample_id
        if not sample_dir.exists():
            continue

        files_to_delete = []

        # Find FASTQ files
        for pattern in ["*.fastq.gz", "*.fastq", "*.fq.gz", "*.fq", "*.amalgkit.fastq.gz"]:
            files_to_delete.extend(sample_dir.glob(pattern))

        # Find SRA files
        for pattern in ["*.sra", "*.sra.part"]:
            files_to_delete.extend(sample_dir.glob(pattern))

        if not files_to_delete:
            continue

        sample_bytes_freed = 0
        for file_path in files_to_delete:
            try:
                file_size = file_path.stat().st_size
                if not dry_run:
                    file_path.unlink()
                sample_bytes_freed += file_size

                if file_path.suffix in (".sra", ".part"):
                    result["sra_files_deleted"] += 1
                else:
                    result["fastq_files_deleted"] += 1

                logger.debug(
                    f"{'Would delete' if dry_run else 'Deleted'}: {file_path.name} ({file_size / 1024**2:.1f} MB)"
                )
            except Exception as e:
                error_msg = f"Failed to delete {file_path}: {e}"
                result["errors"].append(error_msg)
                logger.warning(error_msg)

        if sample_bytes_freed > 0:
            result["samples_cleaned"] += 1
            result["bytes_freed"] += sample_bytes_freed
            logger.info(f"{'Would clean' if dry_run else 'Cleaned'} {sample_id}: {sample_bytes_freed / 1024**3:.2f} GB")

    # Summary
    total_gb = result["bytes_freed"] / 1024**3
    action = "Would free" if dry_run else "Freed"
    logger.info(
        f"Cleanup complete: {result['samples_cleaned']} samples, "
        f"{result['fastq_files_deleted']} FASTQ files, "
        f"{result['sra_files_deleted']} SRA files, "
        f"{action} {total_gb:.2f} GB"
    )

    if result["errors"]:
        logger.warning(f"{len(result['errors'])} errors during cleanup")

    return result


def filter_metadata_for_unquantified(
    config: AmalgkitWorkflowConfig, source_metadata: Path, output_metadata: Path
) -> int:
    """Filter metadata to include only samples not yet quantified.

    Args:
        config: Workflow configuration
        source_metadata: Path to source metadata TSV
        output_metadata: Path to write filtered metadata

    Returns:
        Number of samples remaining after filtering
    """
    import csv

    # Get already quantified samples
    quantified = get_quantified_samples(config)

    if not quantified:
        # No samples quantified yet - copy metadata as-is
        if source_metadata.exists():
            shutil.copy2(source_metadata, output_metadata)
            # Count rows
            with open(source_metadata, "r") as f:
                return sum(1 for _ in f) - 1  # Subtract header
        return 0

    # Read and filter metadata
    try:
        with open(source_metadata, "r", newline="") as f:
            reader = csv.DictReader(f, delimiter="\t")
            fieldnames = reader.fieldnames
            rows = list(reader)
    except Exception as e:
        logger.error(f"Could not read metadata: {e}")
        return 0

    # Filter out quantified samples
    filtered_rows = []
    for row in rows:
        run_id = row.get("run", "")
        if run_id not in quantified:
            filtered_rows.append(row)
        else:
            logger.debug(f"Skipping already-quantified sample: {run_id}")

    skipped_count = len(rows) - len(filtered_rows)
    logger.info(f"Filtered metadata: {len(filtered_rows)} remaining, {skipped_count} already quantified")

    if not filtered_rows:
        logger.info("All samples already quantified - nothing to process")
        return 0

    # Write filtered metadata
    try:
        with open(output_metadata, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
            writer.writeheader()
            writer.writerows(filtered_rows)
        logger.info(f"Created filtered metadata with {len(filtered_rows)} unprocessed samples: {output_metadata}")
    except Exception as e:
        logger.error(f"Could not write filtered metadata: {e}")
        return 0

    return len(filtered_rows)


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
