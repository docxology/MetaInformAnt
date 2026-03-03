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

    # Identify which steps are per-chunk (getfastq, quant) and which are post-process (merge, etc)
    # NOTE: 'integrate' is SKIPPED in streaming mode because amalgkit's integrate.py does a flat
    # os.listdir() scan for .fastq files, but getfastq places files in per-sample subdirectories
    # (e.g., getfastq/SRR123/SRR123.amalgkit.fastq.gz). The integrate step is designed for
    # enriching metadata with private FASTQ files — it's not needed when getfastq already manages
    # the download-extract-name pipeline and metadata is already complete.
    chunk_steps = []
    post_process_steps = []

    for name, params in steps_remaining:
        if name == "integrate":
            logger.info("Skipping 'integrate' step in streaming mode (incompatible with per-sample subdirectory layout)")
            continue
        elif name in ("getfastq", "quant"):
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
    
    # Symlink getfastq output dir to work dir so `quant` can find it
    fastq_src_dir = config.work_dir.parent / "fastq" / "getfastq"
    fastq_link_dir = config.work_dir / "getfastq"
    if fastq_src_dir.exists() and not fastq_link_dir.exists():
        try:
            fastq_link_dir.symlink_to(fastq_src_dir.resolve(), target_is_directory=True)
            logger.info("  Created symlink: work/getfastq -> fastq/getfastq")
        except Exception as e:
            logger.warning(f"  Failed to create getfastq symlink, quant may not find files: {e}")

    # We define a helper that processes exactly ONE sample at a time.
    def process_single_sample(sample_row: Dict[str, str], current_chunk_idx: int, sample_idx_in_chunk: int) -> list[WorkflowStepResult]:
        sample_results = []
        sample_id = sample_row.get(u"run", str(sample_idx_in_chunk))
        
        # Create single-sample metadata file
        single_meta_file = config.work_dir / "metadata" / f"metadata_chunk_{current_chunk_idx}_sample_{sample_id}.tsv"
        with open(single_meta_file, "w", newline="") as sf:
            swriter = csv.DictWriter(sf, fieldnames=fieldnames, delimiter="\t")
            swriter.writeheader()
            swriter.writerow(sample_row)

        sample_success = True
        for step_name, step_params in chunk_steps:
            logger.info(f"  > [Sample {sample_id}] Step: {step_name}")
            
            # Use the single-sample metadata, but scale down threads/jobs
            # so 16 concurrent samples don't spawn 16x16=256 threads.
            single_params = step_params.copy()
            single_params["metadata"] = str(single_meta_file)
            
            # Dynamically throttle cores to prevent system choking
            # If total threads=16, and chunk_size=16, give each sample 1 thread.
            base_threads = int(single_params.get("threads", 1))
            alloc_threads = max(1, base_threads // chunk_size)
            single_params["threads"] = alloc_threads
            if "jobs" in single_params:
                single_params["jobs"] = alloc_threads
                
            if walk:
                sample_results.append(
                    WorkflowStepResult(step_name=f"{step_name}_{sample_id}", return_code=0, success=True)
                )
                continue

            try:
                if step_name == "getfastq":
                    _check_disk_space_or_fail(
                        config.work_dir, min_free_gb=5.0, step_name=f"{step_name} (sample {sample_id})"
                    )

                step_func = step_functions.get(step_name)
                # Ensure the tool output differentiates from other parallel tasks
                result = step_func(single_params)

                step_res = WorkflowStepResult(
                    step_name=f"{step_name}_{sample_id}",
                    return_code=result.returncode,
                    success=result.returncode == 0,
                    error_message=result.stderr if result.returncode != 0 else None,
                )
                sample_results.append(step_res)

                if result.returncode != 0:
                    logger.error(f"  [Sample {sample_id}] Step {step_name} failed with code {result.returncode}")
                    if result.stderr:
                        logger.error(f"  [Sample {sample_id}] STDERR:\n{result.stderr.strip()}")
                    if result.stdout and not result.stderr:
                        logger.error(f"  [Sample {sample_id}] STDOUT (no stderr):\n{result.stdout.strip()}")

                    sample_success = False
                    break  # Stop processing this specific sample
                    
            except Exception as e:
                logger.error(f"  [Sample {sample_id}] Exception in {step_name}: {e}")
                sample_results.append(WorkflowStepResult(f"{step_name}_{sample_id}", 1, False, str(e)))
                sample_success = False
                break

        # Only clean up if this specific sample succeeded completely
        if sample_success:
            logger.info(f"  [Sample {sample_id}] Success. Cleaning up FASTQ files...")
            cleanup_fastqs(config, [sample_id])
        else:
            logger.warning(f"  [Sample {sample_id}] Failed. Skipping cleanup to preserve FASTQ files for debugging.")
            
        return sample_results

    # Continuous pool — as soon as one sample finishes, the next starts immediately.
    # No chunk barriers; max_workers limits concurrency to chunk_size.
    import concurrent.futures
    logger.info(f"Processing ALL {total_samples} samples with max {chunk_size} concurrent workers (no chunk barriers)")
    
    with concurrent.futures.ThreadPoolExecutor(max_workers=chunk_size) as executor:
        futures = {}
        for i, sample in enumerate(samples):
            fut = executor.submit(process_single_sample, sample, 0, i)
            futures[fut] = sample.get("run", str(i))
            
        # Process results as they complete — each finished slot immediately picks up next sample
        for fut in concurrent.futures.as_completed(futures):
            sample_id = futures[fut]
            try:
                res_list = fut.result()
                all_step_results.extend(res_list)
                completed = sum(1 for s in all_step_results if s.success and s.step_name.startswith("quant_"))
                logger.info(f"  [Pool] {sample_id} done. Progress: {completed}/{total_samples} quantified")
            except Exception as exc:
                logger.error(f"  [Pool] {sample_id} raised an unhandled exception: {exc}")
                all_step_results.append(WorkflowStepResult("async_sample_task", 1, False, str(exc)))
                
        # Check if failures should trigger early return
        if check and any(not s.success for s in all_step_results):
            return WorkflowExecutionResult(all_step_results, False, len(all_step_results), 0, 1)

    # After all chunks, run post-process steps (merge, etc)
    overall_success = all(s.success for s in all_step_results) if all_step_results else True
    if overall_success or not check:
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
        success=all(s.return_code == 0 for s in all_step_results) if all_step_results else True,
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
    has_quant_or_merge = (not steps) or any(s in ("quant", "merge", "index") for s in steps)
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

            # Note: config_dir is purposefully injected in workflow_planning.py to point 
            # to the global repository config/amalgkit directory and should not be mutated here.

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

                    if (remaining_count > 5 and estimated_need > (free_gb * 0.8)) or kwargs.get("stream", False):
                        chunk_size = kwargs.get("chunk_size", 5)
                        if kwargs.get("stream", False):
                            logger.info(f"Streaming mode forced with chunk size {chunk_size}")
                        else:
                            logger.warning(
                                f"Insufficient disk space for all-at-once processing: "
                                f"Need ~{estimated_need:.1f}GB, have {free_gb:.1f}GB. "
                                f"Switching to STREAMING MODE (chunk size {chunk_size})."
                            )

                        streaming_result = _execute_streaming_mode(
                            config,
                            steps_planned[i - 1 :],
                            unquantified_metadata,
                            chunk_size=chunk_size,
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
