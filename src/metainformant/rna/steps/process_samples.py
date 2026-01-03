"""Unified sample processing for RNA-seq workflows.

This module provides unified functions for downloading and quantifying RNA-seq samples,
supporting both sequential and parallel processing modes.
"""

from __future__ import annotations

import concurrent.futures
import os
import time
from pathlib import Path
from typing import Any

from metainformant.core import logging
from metainformant.rna.workflow import AmalgkitWorkflowConfig

logger = logging.get_logger(__name__)


def run_download_quant_workflow(config: AmalgkitWorkflowConfig, *,
                               num_workers: int = 1, progress: bool = True) -> list[int]:
    """Run unified download and quantification workflow.

    Args:
        config: Workflow configuration
        num_workers: Number of parallel workers (1 = sequential, >1 = parallel)
        progress: Whether to show progress

    Returns:
        List of return codes (0 = success)
    """
    if num_workers == 1:
        # Sequential processing
        return _run_sequential_workflow(config, progress=progress)
    else:
        # Parallel processing
        return _run_parallel_workflow(config, num_workers=num_workers, progress=progress)


def _run_sequential_workflow(config: AmalgkitWorkflowConfig, *, progress: bool = True) -> list[int]:
    """Run workflow in sequential mode (one sample at a time)."""
    from metainformant.rna.workflow import plan_workflow, execute_workflow

    # Plan the workflow
    steps = plan_workflow(config)

    # Execute steps sequentially
    return_codes = []
    for step_name, step_params in steps:
        if progress:
            logger.info(f"Executing step: {step_name}")

        # Execute single step
        try:
            # Import the appropriate step module
            step_module = __import__(f"metainformant.rna.steps.{step_name}",
                                    fromlist=[f"run_{step_name}"])
            run_func = getattr(step_module, f"run_{step_name}")

            # Execute step
            result = run_func(step_params)
            return_code = result.returncode if hasattr(result, 'returncode') else 0
            return_codes.append(return_code)

        except Exception as e:
            logger.error(f"Step {step_name} failed: {e}")
            return_codes.append(1)

    return return_codes


def _run_parallel_workflow(config: AmalgkitWorkflowConfig, *,
                          num_workers: int = 4, progress: bool = True) -> list[int]:
    """Run workflow in parallel mode."""
    from metainformant.rna.workflow import plan_workflow

    # Plan the workflow
    steps = plan_workflow(config)

    # Execute steps in parallel using ThreadPoolExecutor
    return_codes = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=num_workers) as executor:
        # Submit all steps
        futures = []
        for step_name, step_params in steps:
            future = executor.submit(_execute_step, step_name, step_params, progress)
            futures.append(future)

        # Collect results
        for future in concurrent.futures.as_completed(futures):
            try:
                return_code = future.result()
                return_codes.append(return_code)
            except Exception as e:
                logger.error(f"Step execution failed: {e}")
                return_codes.append(1)

    return return_codes


def _execute_step(step_name: str, step_params: Any, progress: bool) -> int:
    """Execute a single workflow step."""
    if progress:
        logger.info(f"Executing step: {step_name}")

    try:
        # Import the appropriate step module
        step_module = __import__(f"metainformant.rna.steps.{step_name}",
                                fromlist=[f"run_{step_name}"])
        run_func = getattr(step_module, f"run_{step_name}")

        # Execute step
        result = run_func(step_params)
        return_code = result.returncode if hasattr(result, 'returncode') else 0
        return return_code

    except Exception as e:
        logger.error(f"Step {step_name} failed: {e}")
        return 1


def _download_worker(sample_id: str, config: AmalgkitWorkflowConfig) -> bool:
    """Worker function for downloading a single sample.

    Args:
        sample_id: Sample identifier
        config: Workflow configuration

    Returns:
        True if download successful
    """
    try:
        # Import required modules
        from metainformant.rna.steps import getfastq

        # Create step parameters for this sample
        params = type('Params', (), {
            'work_dir': config.work_dir,
            'threads': config.threads,
            'sample_id': sample_id,
        })()

        # Execute download
        result = getfastq.run_getfastq(params)
        return result.returncode == 0

    except Exception as e:
        logger.error(f"Download failed for sample {sample_id}: {e}")
        return False


def _wait_for_fastq_files(sample_id: str, fastq_dir: str | Path,
                         timeout: int = 300, check_interval: int = 5) -> bool:
    """Wait for FASTQ files to be downloaded for a sample.

    Args:
        sample_id: Sample identifier
        fastq_dir: Directory containing FASTQ files
        timeout: Maximum time to wait in seconds
        check_interval: How often to check for files

    Returns:
        True if files found within timeout
    """
    fastq_path = Path(fastq_dir)
    start_time = time.time()

    while time.time() - start_time < timeout:
        # Check for FASTQ files
        fastq_files = list(fastq_path.glob(f"*{sample_id}*.fastq*"))
        if fastq_files:
            # Check that files are not empty and not being written
            valid_files = []
            for fq_file in fastq_files:
                if fq_file.stat().st_size > 0:
                    # Check if file is still being written (size changing)
                    time.sleep(1)
                    if fq_file.stat().st_size == fq_file.stat().st_size:  # Size stable
                        valid_files.append(fq_file)

            if len(valid_files) >= 1:  # At least one valid FASTQ file
                logger.info(f"FASTQ files ready for sample {sample_id}: {valid_files}")
                return True

        time.sleep(check_interval)

    logger.warning(f"Timeout waiting for FASTQ files for sample {sample_id}")
    return False


def _delete_fastq_for_sample(sample_id: str, fastq_dir: str | Path) -> bool:
    """Delete FASTQ files for a sample to free disk space.

    Args:
        sample_id: Sample identifier
        fastq_dir: Directory containing FASTQ files

    Returns:
        True if deletion successful
    """
    try:
        fastq_path = Path(fastq_dir)

        # Find and delete FASTQ files for this sample
        fastq_files = list(fastq_path.glob(f"*{sample_id}*.fastq*"))

        deleted_count = 0
        for fq_file in fastq_files:
            try:
                fq_file.unlink()
                deleted_count += 1
                logger.debug(f"Deleted FASTQ file: {fq_file}")
            except Exception as e:
                logger.warning(f"Failed to delete {fq_file}: {e}")

        if deleted_count > 0:
            logger.info(f"Deleted {deleted_count} FASTQ files for sample {sample_id}")
        else:
            logger.warning(f"No FASTQ files found to delete for sample {sample_id}")

        return True

    except Exception as e:
        logger.error(f"Error deleting FASTQ files for sample {sample_id}: {e}")
        return False


def _sample_already_quantified(sample_id: str, quant_dir: str | Path) -> bool:
    """Check if a sample has already been quantified.

    Args:
        sample_id: Sample identifier
        quant_dir: Quantification directory

    Returns:
        True if sample has abundance.tsv file (already quantified)
    """
    quant_path = Path(quant_dir)
    sample_quant_dir = quant_path / sample_id
    abundance_file = sample_quant_dir / "abundance.tsv"

    return abundance_file.exists() and abundance_file.stat().st_size > 0
