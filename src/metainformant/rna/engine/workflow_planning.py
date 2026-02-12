"""RNA-seq workflow planning, step defaults, and pre-execution preparation.

This module handles workflow step planning, parameter resolution, path logic,
reference genome preparation, and pre-flight checks for RNA-seq workflows.
"""

from __future__ import annotations

import csv
import os
import shutil
import time as time_mod
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, Union

from metainformant.core.utils import logging
from metainformant.rna.engine.workflow_cleanup import (
    check_disk_space,
    check_disk_space_or_fail,
    cleanup_incorrectly_placed_sra_files,
    cleanup_temp_files,
)
from metainformant.rna.engine.workflow_core import AmalgkitWorkflowConfig

logger = logging.get_logger(__name__)

# Private aliases for backward compatibility with internal callers
_cleanup_incorrectly_placed_sra_files = cleanup_incorrectly_placed_sra_files
_cleanup_temp_files = cleanup_temp_files
_check_disk_space = check_disk_space
_check_disk_space_or_fail = check_disk_space_or_fail


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
            bool(per_step_dict.get("cstmm", {}).get("orthogroup_table"))
            or bool(per_step_dict.get("cstmm", {}).get("dir_busco"))
        )
        or (
            bool(per_step_dict.get("csca", {}).get("orthogroup_table"))
            or bool(per_step_dict.get("csca", {}).get("dir_busco"))
        )
        or bool(config.extra_config.get("orthogroup_table"))
        or bool(config.extra_config.get("dir_busco"))
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
            # amalgkit merge typically creates a 'merge' subdirectory, but sometimes outputs directly
            # Check for standard subdirectory first
            metadata_in_subdir = merge_dir / "merge" / "metadata.tsv"
            if metadata_in_subdir.exists():
                step_params["metadata"] = str(metadata_in_subdir)
            else:
                # Fallback to parent directory (seen in pbarbatus)
                step_params["metadata"] = str(merge_dir / "metadata.tsv")

        # PATH RESOLUTION: Auto-adjust integrate fastq_dir to include getfastq subdirectory
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
    skipped_count = output_text.count("Previously-downloaded sra file was detected")
    downloaded_count = output_text.count("Previously-downloaded sra file was not detected")
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
            summary_parts.append(f"{processed_count} file(s) processed")

        if summary_parts:
            summary = ", ".join(summary_parts)
            if total_samples:
                logger.info(f"Step getfastq summary: {summary} (total: {total_samples} samples)")
            else:
                logger.info(f"Step getfastq summary: {summary}")


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

    Args:
        getfastq_dir: Path to the getfastq output directory
        source_metadata: Path to the source metadata TSV file
        output_metadata: Path where filtered metadata will be written

    Returns:
        Number of samples with available SRA files
    """
    sra_dir = getfastq_dir / "sra"
    if not sra_dir.exists():
        logger.warning(f"No SRA directory found at {sra_dir}")
        return 0

    # Find all available SRA files
    available_sra = {}
    for sra_file in sra_dir.glob("*.sra"):
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
            index_dir.mkdir(parents=True, exist_ok=True)
            import shutil

            logger.info(f"Copying shared index to: {index_file}")
            shutil.copy2(shared_index_file, index_file)
            return True

        logger.info(f"Reference index missing at {index_file}. Preparing reference genome...")

        # 2. Get download URL
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

        kallisto_cmd = ["kallisto", "index", "-i", str(index_file), str(dest_file)]

        import subprocess

        result = subprocess.run(kallisto_cmd, capture_output=True, text=True)

        if result.returncode == 0:
            logger.info("Genome index built successfully")
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
            logger.info("SRA files detected. 'redo: no' should trigger extraction.")
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
        return False, None

    elif step_name == "integrate":
        integrated_meta = work_dir / "integration" / "integrated_metadata.json"
        if integrated_meta.exists():
            return True, str(integrated_meta)

        metadata_tsv = work_dir / "metadata" / "metadata.tsv"
        if metadata_tsv.exists():
            integration_dir = work_dir / "integration"
            if integration_dir.exists() and any(integration_dir.iterdir()):
                return True, f"Integration directory exists: {integration_dir}"
        return False, None

    elif step_name == "quant":
        return False, None

    elif step_name == "merge":
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
