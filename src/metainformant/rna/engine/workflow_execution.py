"""RNA-seq workflow execution engine.

This module provides the main workflow execution functions including the
execute_workflow function that runs the complete amalgkit pipeline, streaming
mode for disk-constrained environments, and backward-compatible helpers.
"""

from __future__ import annotations

import csv
import json
import math
import shutil
import time as time_mod
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, Union

from metainformant.core.utils import logging
from metainformant.rna.amalgkit.metadata_filter import filter_selected_metadata
from metainformant.rna.amalgkit.metadata_utils import deduplicate_metadata
from metainformant.rna.amalgkit.metadata_utils import deduplicate_metadata
from metainformant.rna.engine.sra_extraction import manual_integration_fallback
from metainformant.rna.engine.workflow_cleanup import (
    check_disk_space,
    check_disk_space_or_fail,
    cleanup_fastqs,
    filter_metadata_for_unquantified,
)
from metainformant.rna.engine.workflow_core import AmalgkitWorkflowConfig, WorkflowExecutionResult, WorkflowStepResult

logger = logging.get_logger(__name__)

# Private aliases for backward compatibility
_check_disk_space = check_disk_space
_check_disk_space_or_fail = check_disk_space_or_fail


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
                    # Break chunk loop to prevents cascading failures
                    break
            except Exception as e:
                logger.error(f"  Exception in {step_name} chunk {chunk_idx}: {e}")
                all_step_results.append(WorkflowStepResult(f"{step_name}_chunk{chunk_idx}", 1, False, str(e)))
                success = False

        # Aggressive cleanup after chunk
        logger.info(f"  Cleaning up chunk {chunk_idx} FASTQ files...")
        chunk_ids = [row.get("run", "") for row in chunk_samples]
        cleanup_fastqs(config, chunk_ids)

    # After all chunks, run post-process steps (merge, etc)
    if success or not check:
        logger.info("Streaming chunks completed. processing post-chunk steps...")
        for step_name, step_params in post_process_steps:
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
        config as amalgkit_config,
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
    from metainformant.rna.engine.workflow_planning import (
        _is_step_completed,
        plan_workflow,
        prepare_reference_genome,
        sanitize_params_for_cli,
        verify_getfastq_prerequisites,
    )
    from metainformant.rna.engine.workflow_steps import (
        check_step_completion_status,
        handle_post_step_actions,
        log_workflow_summary,
        setup_vdb_config,
        validate_step_prerequisites,
    )

    # Check if we need to prepare reference genome (for quant step)
    has_quant_or_merge = (not steps) or any(s in ("quant", "merge") for s in steps)
    if not steps or has_quant_or_merge:
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
                return_code=126,
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
                            config_dir.symlink_to(config_base_dir.resolve())
                            logger.info(
                                f"Created symlink: {config_dir} -> {config_base_dir.resolve()} (for select step compatibility)"
                            )
                        except (OSError, FileExistsError) as e:
                            if config_dir.exists() or config_dir.is_symlink():
                                logger.debug(f"Config symlink/directory already exists: {config_dir}")
                            else:
                                logger.warning(f"Could not create config symlink: {e}")
                                try:
                                    config_dir.mkdir(parents=True, exist_ok=True)
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
            # Late-binding metadata path
            if step_name in ("getfastq", "integrate", "merge", "quant", "curate", "sanity"):
                selected_metadata = (config.work_dir / "metadata" / "metadata_selected.tsv").absolute()
                planned_metadata = step_params.get("metadata")
                if selected_metadata.exists() and (
                    not planned_metadata or "metadata_selected.tsv" not in str(planned_metadata)
                ):
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

            # Execute step
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
                        try:
                            config_dir.mkdir(parents=True, exist_ok=True)
                            for config_file in config_base_dir.glob("*.config"):
                                shutil.copy2(config_file, config_dir / config_file.name)
                            logger.info(f"Copied config files as fallback")
                        except Exception as e2:
                            logger.warning(f"Could not copy config files: {e2}")
                if (config_dir.exists() or config_dir.is_symlink()) and config_base_dir.exists():
                    step_params["config_dir"] = str(config_dir.absolute())
                    logger.debug(f"Using config directory (symlink): {config_dir.absolute()}")
                elif config_base_dir.exists():
                    step_params["config_dir"] = str(config_base_dir.absolute())
                    logger.debug(f"Using config_base directory (fallback): {config_base_dir.absolute()}")

            # STREAMING MODE DETECTION
            if step_name == "getfastq" and any(s[0] == "quant" for s in steps_planned):
                unquantified_metadata = config.work_dir / "metadata" / "metadata_unquantified.tsv"
                if not unquantified_metadata.exists():
                    selected_metadata = config.work_dir / "metadata" / "metadata_selected.tsv"
                    if selected_metadata.exists():
                        filter_metadata_for_unquantified(config, selected_metadata, unquantified_metadata)

                if unquantified_metadata.exists():
                    with open(unquantified_metadata, "r") as f:
                        reader = csv.DictReader(f, delimiter="\t")
                        remaining_count = sum(1 for _ in reader)

                    ok, free_gb = _check_disk_space(config.work_dir, min_free_gb=10.0)
                    estimated_need = remaining_count * 3.0

                    logger.debug(
                        f"Streaming detection - Remaining: {remaining_count}, Free: {free_gb:.2f}GB, Need: {estimated_need:.2f}GB"
                    )

                    if remaining_count > 5 and estimated_need > (free_gb * 0.8):
                        logger.warning(
                            f"Insufficient disk space for all-at-once processing: "
                            f"Need ~{estimated_need:.1f}GB, have {free_gb:.1f}GB. "
                            f"Switching to STREAMING MODE (chunk size 5)."
                        )

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
            result = step_func(
                sanitized_params,
                show_progress=progress,
                heartbeat_interval=5,
                check=check,
            )

            # Create step result
            error_message = None
            if result.returncode != 0:
                is_completed, completion_indicator = _is_step_completed(step_name, step_params, config)
                if is_completed:
                    logger.warning(
                        f"Step {step_name} reported error (return code {result.returncode}) but outputs exist: {completion_indicator}"
                    )
                    logger.warning(f"Continuing workflow - step appears to have completed successfully despite error")
                    result.returncode = 0
                    logger.info(f"Step {step_name} marked as successful based on output validation")
                elif step_name == "integrate":
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
                import subprocess

                raise subprocess.CalledProcessError(
                    result.returncode,
                    command_str or step_name,
                    output=getattr(result, "stdout", ""),
                    stderr=getattr(result, "stderr", ""),
                )

        except Exception as e:
            error_msg = f"Exception during execution: {e}"
            logger.error(f"Error executing step {step_name}: {e}")

            import time

            if "tqdm" in str(e).lower() or "refresh" in str(e).lower():
                logger.warning(f"tqdm error detected - waiting for subprocess to complete (up to 30 seconds)")
                for i in range(30):
                    time.sleep(1)
                    is_completed, completion_indicator = _is_step_completed(step_name, step_params, config)
                    if is_completed:
                        break
                    if i % 5 == 0:
                        logger.debug(f"Waiting for {step_name} outputs... ({i+1}/30 seconds)")
            else:
                time.sleep(2)

            is_completed, completion_indicator = _is_step_completed(step_name, step_params, config)
            if is_completed:
                logger.warning(f"Step {step_name} raised exception but outputs exist: {completion_indicator}")
                logger.warning(f"Continuing workflow - step appears to have completed successfully despite exception")
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
                continue
            else:
                step_results.append(
                    WorkflowStepResult(step_name=step_name, return_code=1, success=False, error_message=error_msg)
                )
                if check:
                    raise e
                continue

    successful_steps = sum(1 for sr in step_results if sr.success)
    total_steps = len(step_results)
    failed_steps = total_steps - successful_steps

    # Generate workflow summary with remediation steps
    if progress:
        log_workflow_summary(step_results, config)

    # After loop, ensure we record results in manifest
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


# Helper functions for backwards compatibility
def run_config_based_workflow(config_path: Union[str, Path], **kwargs: Any) -> Dict[str, Any]:
    """Run workflow from configuration file.

    Args:
        config_path: Path to configuration file
        **kwargs: Additional workflow options

    Returns:
        Workflow results dictionary
    """
    from metainformant.rna.engine.workflow_core import load_workflow_config

    config = load_workflow_config(config_path)

    # Apply any overrides
    for key, value in kwargs.items():
        if hasattr(config, key):
            setattr(config, key, value)

    return_codes = execute_workflow(config, **kwargs)

    return {"config": config.to_dict(), "return_codes": return_codes, "success": all(rc == 0 for rc in return_codes)}
