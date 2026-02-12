"""RNA-seq workflow orchestration utilities.

This module provides high-level orchestration functions for managing
complex RNA-seq workflows across multiple species and samples.

Functions:
    run_workflow_for_species: Execute workflow for a single species.
    cleanup_unquantified_samples: Remove failed sample artifacts.
    monitor_workflows: Track progress across all workflows.
    discover_species_configs: Find available species configuration files.
    run_parallel_workflows: Execute multiple species workflows concurrently.
    validate_and_execute: Validate environment then execute workflow.
    retry_failed_steps: Re-run only failed steps with exponential backoff.
    get_workflow_status: Query real-time workflow state from the filesystem.
    estimate_workflow_resources: Predict disk, memory, and time requirements.
"""

from __future__ import annotations

import concurrent.futures
import os
import time
from pathlib import Path
from typing import Any, Dict, List, Optional, Set

from metainformant.core.utils import logging

logger = logging.get_logger(__name__)


def run_workflow_for_species(
    config_path: str | Path, species: Optional[str] = None, steps: Optional[List[str]] = None, **kwargs: Any
) -> Dict[str, Any]:
    """Run RNA-seq workflow for a specific species.

    Args:
        config_path: Path to workflow configuration file
        species: Species name (if None, inferred from config)
        steps: Specific steps to run (if None, run all)
        **kwargs: Additional workflow options

    Returns:
        Workflow execution results
    """
    from metainformant.rna.engine.workflow import execute_workflow, load_workflow_config

    logger.info(f"Running workflow for species: {species}")

    # Load configuration
    config = load_workflow_config(config_path)

    # Override species if specified
    if species:
        config.species_list = [species]

    # Execute workflow
    workflow_result = execute_workflow(config, steps=steps, **kwargs)

    results = {
        "species": species or config.species_list[0],
        "config_path": str(config_path),
        "workflow_result": workflow_result,
        "return_codes": workflow_result.return_codes,  # Backward compatibility
        "success": workflow_result.success,
        "steps_executed": workflow_result.total_steps,
        "successful_steps": workflow_result.successful_steps,
        "failed_steps": workflow_result.failed_steps,
        "step_details": [step.__dict__ for step in workflow_result.steps_executed],
        # Backward compatibility for existing scripts
        "completed": [step.step_name for step in workflow_result.steps_executed if step.success],
        "failed": [step.step_name for step in workflow_result.steps_executed if not step.success],
    }

    logger.info(f"Workflow completed for {results['species']}: {'SUCCESS' if results['success'] else 'FAILED'}")

    return results


def cleanup_unquantified_samples(config_path: str | Path) -> tuple[int, int]:
    """Clean up samples that failed quantification.

    Args:
        config_path: Path to workflow configuration file

    Returns:
        Tuple of (quantified_count, failed_count)
    """
    from metainformant.rna.core.cleanup import cleanup_unquantified_samples as cleanup_func
    from metainformant.rna.engine.workflow import load_workflow_config

    config = load_workflow_config(config_path)

    logger.info(f"Cleaning up unquantified samples in {config.work_dir}")

    cleaned_samples = cleanup_func(config.work_dir)

    # Return (quantified=0, failed=len(cleaned)) since we're cleaning up failed samples
    quantified = 0
    failed = len(cleaned_samples)

    logger.info(f"Cleaned up {failed} unquantified samples")

    return (quantified, failed)


def monitor_workflows(work_dir: Path) -> Dict[str, Any]:
    """Monitor the status of all workflows in a directory.

    Args:
        work_dir: Main workflow directory

    Returns:
        Monitoring results for all workflows
    """
    from metainformant.rna.engine.monitoring import assess_all_species_progress

    logger.info(f"Monitoring workflows in {work_dir}")

    assessment = assess_all_species_progress(work_dir)

    # Add summary logging
    overall = assessment.get("overall_status", "unknown")
    total = assessment.get("species_summary", {}).get("overall", {}).get("total_species", 0)
    completed = assessment.get("species_summary", {}).get("overall", {}).get("completed_species", 0)

    logger.info(f"Workflow status: {overall} ({completed}/{total} species completed)")

    return assessment


def discover_species_configs(config_dir: Path) -> List[Dict[str, Any]]:
    """Discover available species configuration files.

    Args:
        config_dir: Directory containing configuration files

    Returns:
        List of species configuration information
    """
    configs = []

    if not config_dir.exists():
        return configs

    # Look for YAML files with amalgkit prefix
    for config_file in config_dir.glob("amalgkit_*.yaml"):
        try:
            species_name = config_file.stem.replace("amalgkit_", "")

            # Skip template and test configs
            if "template" in species_name.lower() or "test" in species_name.lower():
                continue

            config_info = {"species": species_name, "config_file": str(config_file), "path": config_file}

            configs.append(config_info)

        except Exception as e:
            logger.warning(f"Error processing config {config_file}: {e}")

    logger.info(f"Discovered {len(configs)} species configuration files")

    return configs


# ---------------------------------------------------------------------------
# New orchestration capabilities
# ---------------------------------------------------------------------------


def _execute_single_species(
    config: "AmalgkitWorkflowConfig",
    species_name: str,
    **kwargs: Any,
) -> tuple[str, "WorkflowExecutionResult"]:
    """Execute workflow for one species, returning ``(species_name, result)``.

    This is a top-level helper intended for use with :pyclass:`ProcessPoolExecutor`.
    It is deliberately kept at module scope so that it can be pickled.

    Args:
        config: Workflow configuration (must already have ``species_list`` set).
        species_name: Species identifier used as the result key.
        **kwargs: Forwarded to :pyfunc:`execute_workflow`.

    Returns:
        Tuple of ``(species_name, WorkflowExecutionResult)``.
    """
    from metainformant.rna.engine.workflow import WorkflowExecutionResult, WorkflowStepResult, execute_workflow

    try:
        result = execute_workflow(config, **kwargs)
        return species_name, result
    except Exception as exc:
        logger.error(f"Workflow for {species_name} raised an exception: {exc}")
        error_step = WorkflowStepResult(
            step_name="execute_workflow",
            return_code=1,
            success=False,
            error_message=str(exc),
        )
        return species_name, WorkflowExecutionResult(
            steps_executed=[error_step],
            success=False,
            total_steps=1,
            successful_steps=0,
            failed_steps=1,
        )


def run_parallel_workflows(
    configs: List["AmalgkitWorkflowConfig"],
    max_workers: int = 2,
    **kwargs: Any,
) -> Dict[str, "WorkflowExecutionResult"]:
    """Run multiple species workflows in parallel using a process pool.

    Each configuration in *configs* is assumed to target a distinct species.
    Workflows are submitted to a :pyclass:`concurrent.futures.ProcessPoolExecutor`
    so that CPU-bound quantification steps can proceed concurrently.

    A failure in one species does **not** cancel the remaining workflows.
    Instead, the failed species receives a ``WorkflowExecutionResult`` with
    ``success=False`` and an error step describing the exception.

    Args:
        configs: List of workflow configurations, one per species.
        max_workers: Maximum number of concurrent worker processes.
            Defaults to ``2`` to avoid saturating disk I/O on typical
            bioinformatics workstations.
        **kwargs: Additional keyword arguments forwarded to
            :pyfunc:`execute_workflow` for every species.

    Returns:
        Dictionary mapping species name to its ``WorkflowExecutionResult``.
        The species name is taken from ``config.species_list[0]`` (or
        ``"unknown"`` when the list is empty).

    Example::

        results = run_parallel_workflows([cfg_apis, cfg_bombus], max_workers=4)
        for species, result in results.items():
            print(f"{species}: {'OK' if result.success else 'FAIL'}")
    """
    from metainformant.rna.engine.workflow import WorkflowExecutionResult, WorkflowStepResult

    if not configs:
        logger.warning("run_parallel_workflows called with empty config list")
        return {}

    results: Dict[str, WorkflowExecutionResult] = {}
    total = len(configs)
    logger.info(f"Starting parallel workflows for {total} species (max_workers={max_workers})")

    # Build a mapping of species name -> config, guarding against duplicates.
    species_configs: List[tuple[str, "AmalgkitWorkflowConfig"]] = []
    for idx, cfg in enumerate(configs):
        name = cfg.species_list[0] if cfg.species_list else f"unknown_{idx}"
        species_configs.append((name, cfg))

    with concurrent.futures.ProcessPoolExecutor(max_workers=max_workers) as executor:
        future_to_species: Dict[concurrent.futures.Future, str] = {}
        for species_name, cfg in species_configs:
            logger.info(f"Submitting workflow for {species_name}")
            future = executor.submit(_execute_single_species, cfg, species_name, **kwargs)
            future_to_species[future] = species_name

        for future in concurrent.futures.as_completed(future_to_species):
            species_name = future_to_species[future]
            try:
                _name, result = future.result()
                results[_name] = result
                status = "SUCCESS" if result.success else "FAILED"
                completed = len(results)
                logger.info(f"[{completed}/{total}] Workflow for {_name}: {status}")
            except Exception as exc:
                logger.error(f"Unhandled error for {species_name}: {exc}")
                error_step = WorkflowStepResult(
                    step_name="parallel_dispatch",
                    return_code=1,
                    success=False,
                    error_message=str(exc),
                )
                results[species_name] = WorkflowExecutionResult(
                    steps_executed=[error_step],
                    success=False,
                    total_steps=1,
                    successful_steps=0,
                    failed_steps=1,
                )

    succeeded = sum(1 for r in results.values() if r.success)
    logger.info(f"Parallel workflows complete: {succeeded}/{total} succeeded")
    return results


def validate_and_execute(
    config: "AmalgkitWorkflowConfig",
    min_disk_gb: float = 10.0,
    **kwargs: Any,
) -> "WorkflowExecutionResult":
    """Validate environment and configuration, then execute the workflow.

    Pre-flight checks performed before any work begins:

    1. **Configuration validity** -- species list, threads, work_dir.
    2. **Disk space** -- at least *min_disk_gb* free on the work_dir volume.
    3. **External dependencies** -- amalgkit availability, SRA toolkit, kallisto.
    4. **Genome availability** -- genome config present when quant is planned.

    If any check fails the function raises :pyclass:`RuntimeError` with an
    actionable message describing exactly what is wrong and how to fix it.

    Args:
        config: Workflow configuration to validate then run.
        min_disk_gb: Minimum free disk space (GB) required to proceed.
        **kwargs: Forwarded to :pyfunc:`execute_workflow`.

    Returns:
        ``WorkflowExecutionResult`` from :pyfunc:`execute_workflow`.

    Raises:
        RuntimeError: When any pre-flight check fails.
    """
    from metainformant.rna.engine.workflow import execute_workflow, validate_workflow_config
    from metainformant.rna.engine.workflow_cleanup import check_disk_space

    logger.info("Running pre-flight validation...")

    # 1. Configuration validation
    is_valid, errors = validate_workflow_config(config)
    if not is_valid:
        msg = "Configuration validation failed:\n" + "\n".join(f"  - {e}" for e in errors)
        logger.error(msg)
        raise RuntimeError(msg)

    # 2. Disk space check
    work_dir = config.work_dir
    work_dir.mkdir(parents=True, exist_ok=True)
    sufficient, free_gb = check_disk_space(work_dir, min_free_gb=min_disk_gb)
    if not sufficient:
        raise RuntimeError(
            f"Insufficient disk space: {free_gb:.2f} GB free, need at least {min_disk_gb} GB. "
            f"Free up space on the volume containing {work_dir} before retrying."
        )
    logger.info(f"Disk space OK: {free_gb:.2f} GB free (minimum {min_disk_gb} GB)")

    # 3. External dependency checks
    _check_external_dependencies()

    # 4. Genome availability (needed for quant/merge steps)
    steps_requested = kwargs.get("steps")
    needs_genome = steps_requested is None or any(s in ("quant", "merge") for s in (steps_requested or []))
    if needs_genome and not config.genome:
        raise RuntimeError(
            "Genome configuration is required for quant/merge steps but none was provided. "
            "Set 'genome' in your config with at least an 'accession' field."
        )

    logger.info("All pre-flight checks passed -- starting workflow execution")
    return execute_workflow(config, **kwargs)


def _check_external_dependencies() -> None:
    """Verify that critical external tools are available on ``$PATH``.

    Raises:
        RuntimeError: If a required tool is missing, with installation
            instructions.
    """
    import shutil

    tools = {
        "kallisto": "Install kallisto: conda install -c bioconda kallisto",
        "prefetch": "Install SRA Toolkit: conda install -c bioconda sra-tools",
        "fasterq-dump": "Install SRA Toolkit: conda install -c bioconda sra-tools",
    }
    missing: List[str] = []
    for tool, install_hint in tools.items():
        if shutil.which(tool) is None:
            missing.append(f"  - {tool}: {install_hint}")

    if missing:
        msg = "Missing required external tools:\n" + "\n".join(missing)
        logger.warning(msg)
        # This is a warning, not a hard failure, because amalgkit may handle
        # installation automatically when auto_install_amalgkit is True.
        logger.info("Proceeding anyway -- amalgkit may install missing tools automatically")


def retry_failed_steps(
    config: "AmalgkitWorkflowConfig",
    previous_result: "WorkflowExecutionResult",
    max_retries: int = 3,
    base_delay: float = 5.0,
    **kwargs: Any,
) -> "WorkflowExecutionResult":
    """Re-run only the steps that failed in a previous execution.

    Uses exponential back-off between retry attempts (``base_delay * 2**attempt``
    seconds).  Steps that already succeeded are not re-executed; their original
    :pyclass:`WorkflowStepResult` entries are preserved in the combined result.

    Args:
        config: Workflow configuration (same one used for the original run).
        previous_result: The ``WorkflowExecutionResult`` from the prior run.
        max_retries: Maximum number of retry attempts per failed step group.
            Defaults to ``3``.
        base_delay: Initial delay in seconds before the first retry.
            Defaults to ``5.0``.
        **kwargs: Forwarded to :pyfunc:`execute_workflow`.

    Returns:
        Combined ``WorkflowExecutionResult`` merging the originally-successful
        steps with the results of the retried steps.

    Example::

        result = execute_workflow(config)
        if not result.success:
            result = retry_failed_steps(config, result, max_retries=2)
    """
    from metainformant.rna.engine.workflow import WorkflowExecutionResult, WorkflowStepResult, execute_workflow

    failed_step_names: List[str] = [step.step_name for step in previous_result.steps_executed if not step.success]

    if not failed_step_names:
        logger.info("No failed steps to retry -- returning previous result as-is")
        return previous_result

    logger.info(f"Retrying {len(failed_step_names)} failed steps: {failed_step_names}")

    # Collect originally-successful steps (these will be preserved)
    successful_originals: List[WorkflowStepResult] = [step for step in previous_result.steps_executed if step.success]

    last_retry_result: Optional["WorkflowExecutionResult"] = None

    for attempt in range(1, max_retries + 1):
        delay = base_delay * (2 ** (attempt - 1))
        logger.info(f"Retry attempt {attempt}/{max_retries} for steps {failed_step_names} (delay={delay:.1f}s)")
        time.sleep(delay)

        try:
            last_retry_result = execute_workflow(config, steps=failed_step_names, **kwargs)
        except Exception as exc:
            logger.error(f"Retry attempt {attempt} raised exception: {exc}")
            error_step = WorkflowStepResult(
                step_name=",".join(failed_step_names),
                return_code=1,
                success=False,
                error_message=f"Retry attempt {attempt} exception: {exc}",
            )
            last_retry_result = WorkflowExecutionResult(
                steps_executed=[error_step],
                success=False,
                total_steps=1,
                successful_steps=0,
                failed_steps=1,
            )

        # Check which steps still fail
        still_failing = [step.step_name for step in last_retry_result.steps_executed if not step.success]

        if not still_failing:
            logger.info(f"All failed steps recovered on attempt {attempt}")
            break

        logger.warning(f"After attempt {attempt}: still failing: {still_failing}")
        failed_step_names = still_failing

    # Merge: keep original successes + latest retry results
    assert last_retry_result is not None  # noqa: S101 -- guaranteed by loop
    merged_steps = successful_originals + list(last_retry_result.steps_executed)
    total = len(merged_steps)
    succeeded = sum(1 for s in merged_steps if s.success)
    failed = total - succeeded

    combined = WorkflowExecutionResult(
        steps_executed=merged_steps,
        success=(failed == 0),
        total_steps=total,
        successful_steps=succeeded,
        failed_steps=failed,
    )

    logger.info(f"Retry complete: {succeeded}/{total} steps succeeded")
    return combined


def get_workflow_status(config: "AmalgkitWorkflowConfig") -> Dict[str, Any]:
    """Query real-time workflow status by inspecting the filesystem.

    Examines the work directory for evidence of each pipeline step's
    completion, including manifest files, metadata, FASTQ downloads,
    quantification results, and merged output.

    Args:
        config: Workflow configuration whose ``work_dir`` will be inspected.

    Returns:
        Dictionary with comprehensive status information::

            {
                "work_dir": "/path/to/work",
                "species": ["Apis_mellifera"],
                "manifest_exists": True,
                "steps": {
                    "metadata": {"completed": True, "files": [...]},
                    "getfastq": {"completed": True, "sample_count": 42},
                    ...
                },
                "samples": {
                    "quantified": ["SRR123", ...],
                    "quantified_count": 42,
                    "with_fastq": ["SRR456", ...],
                    "fastq_count": 3,
                },
                "disk_free_gb": 123.4,
                "overall_status": "in_progress",
            }
    """
    from metainformant.rna.engine.workflow_cleanup import check_disk_space, get_quantified_samples

    work_dir = config.work_dir
    status: Dict[str, Any] = {
        "work_dir": str(work_dir),
        "species": list(config.species_list),
        "manifest_exists": config.manifest_path.exists(),
        "steps": {},
        "samples": {},
        "disk_free_gb": 0.0,
        "overall_status": "unknown",
    }

    # Disk space
    _, free_gb = check_disk_space(work_dir)
    status["disk_free_gb"] = round(free_gb, 2)

    # --- Per-step status ------------------------------------------------

    # Metadata step
    metadata_dir = work_dir / "metadata"
    metadata_files = list(metadata_dir.glob("*.tsv")) if metadata_dir.exists() else []
    status["steps"]["metadata"] = {
        "completed": len(metadata_files) > 0,
        "files": [str(f.name) for f in metadata_files],
    }

    # Config step
    config_file = work_dir / "config"
    config_files = list(config_file.glob("*")) if config_file.exists() else []
    status["steps"]["config"] = {
        "completed": len(config_files) > 0,
        "files": [str(f.name) for f in config_files],
    }

    # Select step
    select_dir = work_dir / "select"
    select_files = list(select_dir.glob("*.tsv")) if select_dir.exists() else []
    status["steps"]["select"] = {
        "completed": len(select_files) > 0,
        "files": [str(f.name) for f in select_files],
    }

    # Getfastq step -- count sample directories with FASTQ files
    fastq_dirs = [
        work_dir / "fastq" / "getfastq",
        work_dir / "fastq",
        work_dir / "getfastq",
    ]
    fastq_samples: Set[str] = set()
    for fdir in fastq_dirs:
        if not fdir.exists():
            continue
        for child in fdir.iterdir():
            if child.is_dir() and child.name.startswith(("SRR", "ERR", "DRR")):
                fq_files = list(child.glob("*.fastq*")) + list(child.glob("*.fq*"))
                if fq_files:
                    fastq_samples.add(child.name)
    status["steps"]["getfastq"] = {
        "completed": len(fastq_samples) > 0,
        "sample_count": len(fastq_samples),
    }

    # Quant step
    quantified: Set[str] = set()
    try:
        quantified = get_quantified_samples(config)
    except Exception as exc:
        logger.debug(f"Could not determine quantified samples: {exc}")
    status["steps"]["quant"] = {
        "completed": len(quantified) > 0,
        "sample_count": len(quantified),
    }

    # Merge step
    merge_dir = work_dir / "merge"
    merge_files = list(merge_dir.glob("*")) if merge_dir.exists() else []
    status["steps"]["merge"] = {
        "completed": len(merge_files) > 0,
        "files": [str(f.name) for f in merge_files[:20]],  # Cap to avoid huge lists
    }

    # Curate step
    curate_dir = work_dir / "curate"
    curate_files = list(curate_dir.glob("*")) if curate_dir.exists() else []
    status["steps"]["curate"] = {
        "completed": len(curate_files) > 0,
        "files": [str(f.name) for f in curate_files[:20]],
    }

    # --- Sample-level summary -------------------------------------------
    status["samples"] = {
        "quantified": sorted(quantified),
        "quantified_count": len(quantified),
        "with_fastq": sorted(fastq_samples),
        "fastq_count": len(fastq_samples),
    }

    # --- Overall status determination -----------------------------------
    if len(curate_files) > 0 and len(quantified) > 0:
        status["overall_status"] = "completed"
    elif len(quantified) > 0:
        status["overall_status"] = "in_progress"
    elif len(fastq_samples) > 0:
        status["overall_status"] = "downloading"
    elif len(metadata_files) > 0:
        status["overall_status"] = "metadata_ready"
    elif work_dir.exists():
        status["overall_status"] = "initialized"
    else:
        status["overall_status"] = "not_started"

    logger.info(
        f"Workflow status for {status['species']}: {status['overall_status']} "
        f"(quantified={len(quantified)}, fastq={len(fastq_samples)})"
    )
    return status


def estimate_workflow_resources(
    config: "AmalgkitWorkflowConfig",
    n_samples: int,
) -> Dict[str, Any]:
    """Estimate the disk, memory, and time needed for a workflow run.

    Estimates are based on empirical observations from RNA-seq workflows:

    * **Disk**: ~2 GB per sample for raw FASTQ + ~200 MB per sample for
      quantification artifacts, plus genome index overhead.
    * **Memory**: kallisto quantification typically needs 4-8 GB; the
      genome index is the dominant factor.
    * **Time**: Download bandwidth (~1 sample/min on average), plus
      ~2 min/sample for quantification (highly dependent on read depth).

    The returned estimates are intentionally *conservative* (rounded up)
    so callers can plan for worst-case resource requirements.

    Args:
        config: Workflow configuration (used to infer genome size and
            thread count).
        n_samples: Expected number of samples to process.

    Returns:
        Dictionary with resource estimates::

            {
                "n_samples": 100,
                "threads": 8,
                "estimated_disk_gb": 250.0,
                "estimated_memory_gb": 8.0,
                "estimated_hours": 6.5,
                "breakdown": {
                    "fastq_disk_gb": 200.0,
                    "quant_disk_gb": 20.0,
                    "index_disk_gb": 5.0,
                    "download_hours": 1.7,
                    "quant_hours": 3.3,
                    "merge_hours": 0.5,
                },
            }
    """
    # --- Constants (empirical) -------------------------------------------
    FASTQ_GB_PER_SAMPLE = 2.0  # compressed paired-end FASTQ
    QUANT_GB_PER_SAMPLE = 0.2  # abundance + aux files
    INDEX_BASE_GB = 5.0  # default genome index size estimate
    MEMORY_BASE_GB = 4.0  # baseline memory for kallisto
    MEMORY_PER_THREAD_GB = 0.5  # additional memory per thread
    DOWNLOAD_MINUTES_PER_SAMPLE = 1.0  # average download time
    QUANT_MINUTES_PER_SAMPLE = 2.0  # average quantification time
    MERGE_BASE_MINUTES = 10.0  # base time for merge step
    MERGE_MINUTES_PER_SAMPLE = 0.3  # additional merge time per sample

    threads = max(config.threads, 1)

    # Genome index size heuristic: if genome config specifies known info
    # we can refine, otherwise use default.
    genome_info = config.genome or {}
    index_gb = INDEX_BASE_GB
    # Some configs include an approximate genome_size_mb field
    genome_size_mb = genome_info.get("genome_size_mb", 0)
    if genome_size_mb:
        # Transcriptome index is roughly 1/10 genome size
        index_gb = max(INDEX_BASE_GB, genome_size_mb / 1024.0 * 0.1)

    # Disk
    fastq_disk = n_samples * FASTQ_GB_PER_SAMPLE
    quant_disk = n_samples * QUANT_GB_PER_SAMPLE
    total_disk = fastq_disk + quant_disk + index_gb
    # Round up to nearest integer
    total_disk = float(int(total_disk) + (1 if total_disk % 1 > 0 else 0))

    # Memory
    memory = MEMORY_BASE_GB + threads * MEMORY_PER_THREAD_GB
    memory = float(int(memory) + (1 if memory % 1 > 0 else 0))

    # Time (hours)
    download_minutes = n_samples * DOWNLOAD_MINUTES_PER_SAMPLE
    # Quantification parallelism: limited by thread count
    effective_parallelism = max(1, threads // 2)  # conservative
    quant_minutes = (n_samples * QUANT_MINUTES_PER_SAMPLE) / effective_parallelism
    merge_minutes = MERGE_BASE_MINUTES + n_samples * MERGE_MINUTES_PER_SAMPLE

    download_hours = download_minutes / 60.0
    quant_hours = quant_minutes / 60.0
    merge_hours = merge_minutes / 60.0
    total_hours = download_hours + quant_hours + merge_hours
    # Round up to one decimal
    total_hours = round(total_hours + 0.05, 1)

    estimate: Dict[str, Any] = {
        "n_samples": n_samples,
        "threads": threads,
        "estimated_disk_gb": total_disk,
        "estimated_memory_gb": memory,
        "estimated_hours": total_hours,
        "breakdown": {
            "fastq_disk_gb": round(fastq_disk, 1),
            "quant_disk_gb": round(quant_disk, 1),
            "index_disk_gb": round(index_gb, 1),
            "download_hours": round(download_hours, 1),
            "quant_hours": round(quant_hours, 1),
            "merge_hours": round(merge_hours, 1),
        },
    }

    logger.info(
        f"Resource estimate for {n_samples} samples ({threads} threads): "
        f"~{total_disk:.0f} GB disk, ~{memory:.0f} GB RAM, ~{total_hours:.1f} hours"
    )

    return estimate
