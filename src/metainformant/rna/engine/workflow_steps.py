"""Workflow step helpers - prerequisite validation, post-step actions, and setup.

Extracted from workflow.py to reduce its complexity. These functions handle
the vdb-config setup, step prerequisite validation, post-step processing,
and workflow summary reporting.
"""

from __future__ import annotations

import csv
import json
import shutil
import subprocess
import time as time_mod
from pathlib import Path
from typing import TYPE_CHECKING, Any, Dict, List, Optional, Tuple

from metainformant.core.utils import logging
from metainformant.rna.engine.sra_extraction import extract_sra_directly, manual_integration_fallback

if TYPE_CHECKING:
    from metainformant.rna.engine.workflow import AmalgkitWorkflowConfig, WorkflowStepResult

logger = logging.get_logger(__name__)


def setup_vdb_config(
    config: AmalgkitWorkflowConfig,
    steps_planned: List[Tuple[str, Dict[str, Any]]],
) -> List[Tuple[str, Dict[str, Any]]]:
    """Configure vdb-config repository path and prepare getfastq step.

    Sets up SRA Toolkit's vdb-config to download to the correct location,
    cleans up misplaced SRA files, filters metadata for unquantified samples,
    and prepares extraction metadata.

    Args:
        config: Workflow configuration
        steps_planned: List of (step_name, params) tuples. May be mutated.

    Returns:
        Updated steps_planned (getfastq may be removed if all samples quantified)
    """
    try:
        getfastq_params = next((params for step, params in steps_planned if step == "getfastq"), None)
        if getfastq_params is None:
            logger.warning("Could not find getfastq step parameters, skipping vdb-config setup")
            return steps_planned

        fastq_out_dir = getfastq_params.get("out_dir", str(config.work_dir / "fastq"))
        fastq_dir = Path(fastq_out_dir)
        getfastq_dir = fastq_dir / "getfastq" if fastq_dir.name != "getfastq" else fastq_dir
        getfastq_dir.mkdir(parents=True, exist_ok=True)

        # Check disk space
        free_gb = check_disk_space_or_fail(getfastq_dir, min_free_gb=5.0, step_name="getfastq")
        if free_gb < 20.0:
            logger.warning(f"Low disk space ({free_gb:.1f}GB) - workflow may fail during downloads")

        # Clean up misplaced SRA files
        cleanup_incorrectly_placed_sra_files(getfastq_dir)

        # Clean up temp files
        repo_root = Path(__file__).resolve().parent.parent.parent.parent
        tmp_dir = repo_root / ".tmp" / "fasterq-dump"
        cleanup_temp_files(tmp_dir, max_size_gb=50.0)

        # Configure vdb-config repository path
        result = subprocess.run(
            ["vdb-config", "-s", f"/repository/user/main/public/root={getfastq_dir}"],
            capture_output=True,
            text=True,
            timeout=5,
        )
        if result.returncode == 0:
            logger.info(f"Configured vdb-config repository path for prefetch: {getfastq_dir}")
        else:
            logger.warning(f"Could not set vdb-config repository path (may require interactive): {result.stderr[:100]}")
            logger.info("Will rely on TMPDIR and VDB_CONFIG environment variables")

        # Prepare extraction directories
        from metainformant.rna.engine.workflow import prepare_extraction_directories

        cleaned = prepare_extraction_directories(getfastq_dir)
        if cleaned > 0:
            logger.info(f"Prepared {cleaned} sample directories for extraction")

        # Filter metadata for unquantified samples
        source_metadata = config.work_dir / "metadata" / "metadata_selected.tsv"
        unquantified_metadata = config.work_dir / "metadata" / "metadata_unquantified.tsv"

        if source_metadata.exists():
            remaining = filter_metadata_for_unquantified(config, source_metadata, unquantified_metadata)
            if remaining == 0:
                logger.info("All samples already quantified - skipping getfastq step")
                steps_planned = [(name, params) for name, params in steps_planned if name != "getfastq"]
            elif remaining < sum(1 for _ in open(source_metadata)) - 1:
                source_metadata = unquantified_metadata
                logger.info(f"Using unquantified metadata ({remaining} samples remaining)")

        # Create extraction metadata
        if source_metadata.exists() and any(s[0] == "getfastq" for s in steps_planned):
            from metainformant.rna.engine.workflow import create_extraction_metadata

            extraction_metadata = config.work_dir / "metadata" / "metadata_extraction.tsv"
            num_samples = create_extraction_metadata(getfastq_dir, source_metadata, extraction_metadata)
            if num_samples > 0:
                for i, (step_name, step_params) in enumerate(steps_planned):
                    if step_name == "getfastq":
                        steps_planned[i] = (step_name, {**step_params, "metadata": str(extraction_metadata)})
                        logger.info(f"Using extraction metadata ({num_samples} samples) for getfastq step")

    except Exception as e:
        logger.warning(f"Could not configure vdb-config: {e}")
        logger.info("Will rely on TMPDIR and VDB_CONFIG environment variables")

    return steps_planned


def check_step_completion_status(
    steps_planned: List[Tuple[str, Dict[str, Any]]],
    config: AmalgkitWorkflowConfig,
    **kwargs: Any,
) -> Tuple[List[Tuple[str, str]], List[str]]:
    """Check which planned steps are already completed.

    Args:
        steps_planned: List of (step_name, params) tuples
        config: Workflow configuration
        **kwargs: Extra options, including 'redo' override

    Returns:
        Tuple of (completed_steps, steps_to_run) where completed_steps
        is list of (name, indicator) tuples and steps_to_run is list of names
    """
    from metainformant.rna.engine.workflow import _is_step_completed

    completed_steps: List[Tuple[str, str]] = []
    steps_to_run: List[str] = []

    for step_name, step_params in steps_planned:
        is_completed, completion_indicator = _is_step_completed(step_name, step_params, config)

        redo_value = kwargs.get("redo", step_params.get("redo", "no"))
        if isinstance(redo_value, bool):
            force_redo = redo_value
        else:
            force_redo = str(redo_value).lower() in ("yes", "true", "1")

        if is_completed and not force_redo:
            completed_steps.append((step_name, completion_indicator))
        else:
            steps_to_run.append(step_name)

    logger.info("Workflow status summary:")
    logger.info(f"  Total steps planned: {len(steps_planned)}")
    if completed_steps:
        logger.info(f"  Steps already completed ({len(completed_steps)}): {', '.join([s[0] for s in completed_steps])}")
        for step_name, indicator in completed_steps:
            logger.info(f"    - {step_name}: {indicator}")
    if steps_to_run:
        logger.info(f"  Steps to run ({len(steps_to_run)}): {', '.join(steps_to_run)}")
    else:
        logger.info("  All steps already completed - nothing to run")

    return completed_steps, steps_to_run


def validate_step_prerequisites(
    step_name: str,
    step_params: Dict[str, Any],
    config: AmalgkitWorkflowConfig,
    steps_planned: List[Tuple[str, Dict[str, Any]]],
    steps_config: Dict[str, Any],
) -> Optional[str]:
    """Validate prerequisites for a workflow step before execution.

    Checks that required inputs exist (FASTQ files, metadata, tools, etc.)
    before executing the step to fail fast with actionable error messages.

    Args:
        step_name: Name of the step to validate
        step_params: Parameters for the step
        config: Workflow configuration
        steps_planned: All planned steps
        steps_config: Steps configuration dict

    Returns:
        Error message if validation fails, None if OK
    """
    if step_name == "integrate":
        return _validate_integrate_prerequisites(config, steps_config)
    elif step_name == "quant":
        return _validate_quant_prerequisites(config, steps_config)
    elif step_name == "merge":
        return _validate_merge_prerequisites(config, step_params, steps_config)
    elif step_name in ("curate", "cstmm"):
        return _validate_r_step_prerequisites(step_name, step_params, config, steps_config)
    return None


def _find_fastq_files(config: AmalgkitWorkflowConfig, steps_config: Dict[str, Any]) -> Tuple[List[Path], List[str]]:
    """Find FASTQ files across all possible locations.

    Returns:
        Tuple of (found_files, checked_directories)
    """
    fastq_locations = [
        Path(steps_config.get("getfastq", {}).get("out_dir", config.work_dir / "fastq")),
        config.work_dir.parent / "fastq",
        config.work_dir,
        config.work_dir / "fastq",
    ]

    fastq_files: List[Path] = []
    checked_dirs: List[str] = []

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

    return fastq_files, checked_dirs


def _validate_integrate_prerequisites(config: AmalgkitWorkflowConfig, steps_config: Dict[str, Any]) -> Optional[str]:
    """Validate prerequisites for integrate step."""
    fastq_files, checked_dirs = _find_fastq_files(config, steps_config)
    if not fastq_files:
        return (
            f"PREREQUISITE CHECK FAILED: No FASTQ files found before integrate step.\n"
            f"  - Checked locations: {', '.join(checked_dirs)}\n"
            f"  - Checked for: *.fastq, *.fastq.gz, *.fq, *.fq.gz\n\n"
            f"REMEDIATION:\n"
            f"  1. Ensure getfastq step completed successfully\n"
            f"  2. Check validation report: {config.work_dir / 'validation' / 'getfastq_validation.json'}\n"
            f"  3. Re-run getfastq step if FASTQ files are missing\n"
            f"  4. Verify amalgkit getfastq extracted files correctly"
        )
    return None


def _validate_quant_prerequisites(config: AmalgkitWorkflowConfig, steps_config: Dict[str, Any]) -> Optional[str]:
    """Validate prerequisites for quant step."""
    # Check integrate output
    integrated_metadata = config.work_dir / "integration" / "integrated_metadata.json"
    metadata_tsv = config.work_dir / "metadata" / "metadata.tsv"
    if not integrated_metadata.exists() and not metadata_tsv.exists():
        return (
            f"PREREQUISITE CHECK FAILED: Integrate step must complete before quant.\n"
            f"  - Expected: {integrated_metadata} or {metadata_tsv}\n\n"
            f"REMEDIATION:\n"
            f"  1. Ensure integrate step completed successfully\n"
            f"  2. Re-run integrate step if needed\n"
        )

    # Check FASTQ files
    fastq_files, checked_dirs = _find_fastq_files(config, steps_config)
    if not fastq_files:
        return (
            f"PREREQUISITE CHECK FAILED: No FASTQ files found before quant step.\n"
            f"  - Checked locations: {', '.join(checked_dirs)}\n\n"
            f"REMEDIATION:\n"
            f"  1. Ensure getfastq and integrate steps completed successfully\n"
            f"  2. Check validation reports in: {config.work_dir / 'validation'}\n"
            f"  3. Re-run getfastq step if FASTQ files are missing"
        )

    # Check quantification tools
    try:
        from metainformant.rna.core.deps import check_quantification_tools

        quant_tools = check_quantification_tools()
        available_tools = [tool for tool, (avail, _) in quant_tools.items() if avail]
        if not available_tools:
            return (
                f"PREREQUISITE CHECK FAILED: No quantification tools available.\n"
                f"  - Checked: kallisto, salmon\n"
                f"  - Status: {quant_tools}\n\n"
                f"REMEDIATION:\n"
                f"  1. Install kallisto: conda install -c bioconda kallisto\n"
                f"  2. Install salmon: conda install -c bioconda salmon\n"
                f"  3. Verify tools are in PATH: which kallisto / which salmon"
            )
    except Exception as e:
        logger.warning(f"Could not check quantification tools: {e}")

    return None


def _validate_merge_prerequisites(
    config: AmalgkitWorkflowConfig,
    step_params: Dict[str, Any],
    steps_config: Dict[str, Any],
) -> Optional[str]:
    """Validate prerequisites for merge step."""
    quant_dir_raw = steps_config.get("quant", {}).get("out_dir", config.work_dir / "quant")
    quant_dir = Path(quant_dir_raw)
    if quant_dir.name != "quant":
        actual_quant_dir = quant_dir / "quant"
        if actual_quant_dir.exists():
            quant_dir = actual_quant_dir

    quant_files = (
        list(quant_dir.glob("**/abundance.tsv"))
        + list(quant_dir.glob("**/*_abundance.tsv"))
        + list(quant_dir.glob("**/quant.sf"))
    )
    if not quant_files:
        return (
            f"PREREQUISITE CHECK FAILED: No quantification files found before merge step.\n"
            f"  - Expected location: {quant_dir}\n"
            f"  - Checked for: abundance.tsv, quant.sf\n\n"
            f"REMEDIATION:\n"
            f"  1. Ensure quant step completed successfully\n"
            f"  2. Re-run quant step if quantification files are missing"
        )

    # Check R/Rscript
    from metainformant.rna.core.environment import check_rscript

    r_available, r_message = check_rscript()
    if not r_available:
        return (
            f"PREREQUISITE CHECK FAILED: R/Rscript not available for merge step.\n"
            f"  Status: {r_message}\n\n"
            f"REMEDIATION:\n"
            f"  1. Install R: brew install r (on Mac) or apt-get install r-base-core\n"
        )

    # Check ggplot2
    try:
        check_ggplot2 = subprocess.run(
            ["Rscript", "-e", "library(ggplot2)"], capture_output=True, text=True, timeout=10
        )
        if check_ggplot2.returncode != 0:
            return (
                f"PREREQUISITE CHECK FAILED: R package 'ggplot2' not available for merge step.\n"
                f"  Installation: Rscript -e \"install.packages('ggplot2', repos='https://cloud.r-project.org')\"\n"
            )
    except Exception as e:
        logger.warning(f"Failed to check R packages: {e}")

    # Bridge quant results for merge
    merge_out_dir_path = Path(step_params.get("out_dir", config.work_dir / "merge"))
    merge_quant_link = merge_out_dir_path / "quant"
    if quant_dir.resolve() != merge_quant_link.resolve():
        if not merge_quant_link.exists() or not list(merge_quant_link.glob("**/*.tsv")):
            logger.info(f"Bridging quantification results for merge: {quant_dir} -> {merge_quant_link}")
            try:
                merge_quant_link.parent.mkdir(parents=True, exist_ok=True)
                if merge_quant_link.exists() and not merge_quant_link.is_symlink():
                    if merge_quant_link.is_dir() and not any(merge_quant_link.iterdir()):
                        merge_quant_link.rmdir()
                if not merge_quant_link.exists():
                    merge_quant_link.symlink_to(quant_dir.resolve())
            except Exception as e:
                logger.warning(f"Could not bridge quant results for merge: {e}")

    return None


def _validate_r_step_prerequisites(
    step_name: str,
    step_params: Dict[str, Any],
    config: AmalgkitWorkflowConfig,
    steps_config: Dict[str, Any],
) -> Optional[str]:
    """Validate prerequisites for R-dependent steps (curate, cstmm)."""
    if not shutil.which("Rscript"):
        return f"PREREQUISITE CHECK FAILED: Rscript not found for {step_name} step."

    # Bridge merge results for curate
    if step_name in ("curate", "cstmm"):
        merge_params = steps_config.get("merge", {})
        merge_out_dir = Path(merge_params.get("out_dir", config.work_dir / "merged"))
        merge_results_dir = merge_out_dir / "merge"
        curate_out_dir = Path(step_params.get("out_dir", config.work_dir / "curate"))
        curate_merge_link = curate_out_dir / "merge"

        if merge_results_dir.exists() and merge_results_dir.resolve() != curate_merge_link.resolve():
            if not curate_merge_link.exists():
                logger.info(f"Bridging merge results for {step_name}: {merge_results_dir} -> {curate_merge_link}")
                try:
                    curate_merge_link.parent.mkdir(parents=True, exist_ok=True)
                    curate_merge_link.symlink_to(merge_results_dir.resolve())
                except Exception as e:
                    logger.warning(f"Could not bridge merge results for {step_name}: {e}")

    return None


def handle_post_step_actions(
    step_name: str,
    step_params: Dict[str, Any],
    result: Any,
    config: AmalgkitWorkflowConfig,
    steps_planned: List[Tuple[str, Dict[str, Any]]],
    steps_config: Dict[str, Any],
    check: bool = False,
    step_results: Optional[List] = None,
) -> None:
    """Handle post-step actions like validation, cleanup, and metadata filtering.

    Args:
        step_name: Name of the step that just completed
        step_params: Parameters for the step
        result: Step execution result (has .returncode, .stdout, .stderr)
        config: Workflow configuration
        steps_planned: All planned steps
        steps_config: Steps configuration dict
        check: Whether to stop on first failure
        step_results: Mutable list to append additional step results to
    """
    from metainformant.rna.amalgkit.metadata_filter import filter_selected_metadata
    from metainformant.rna.amalgkit.metadata_utils import deduplicate_metadata

    if step_results is None:
        step_results = []

    # After integrate: deduplicate metadata
    if step_name == "integrate" and result.returncode == 0:
        metadata_dir = config.work_dir / "metadata"
        for m_file in [
            metadata_dir / "metadata_updated_for_private_fastq.tsv",
            metadata_dir / "metadata.tsv",
        ]:
            if m_file.exists():
                logger.info(f"Deduplicating metadata after integrate: {m_file}")
                deduplicate_metadata(m_file)

    # After config: create symlink for select step
    if step_name == "config":
        _ensure_config_symlink(config)

    # After select: create filtered metadata
    if step_name == "select" and result.returncode == 0:
        _create_filtered_metadata(config, check, step_results)

    # After getfastq: validate and run fallback extraction if needed
    if step_name == "getfastq" and result.returncode == 0:
        output_text = (result.stdout or "") + (result.stderr or "")
        if output_text:
            from metainformant.rna.engine.workflow import _log_getfastq_summary

            _log_getfastq_summary(output_text, logger)

        _validate_getfastq_results(config, steps_config, check, step_results)

    # After quant: validate and cleanup
    if step_name == "quant" and result.returncode == 0:
        _validate_quant_results(config, check, step_results)


def _ensure_config_symlink(config: AmalgkitWorkflowConfig) -> None:
    """Create config -> config_base symlink for select step compatibility."""
    config_base_dir = config.work_dir / "config_base"
    config_dir = config.work_dir / "config"
    if config_base_dir.exists() and not config_dir.exists():
        try:
            config_dir.symlink_to(config_base_dir.resolve())
            logger.info(f"Created symlink: {config_dir} -> {config_base_dir.resolve()}")
        except Exception as e:
            logger.warning(f"Could not create config symlink: {e}")
            try:
                config_dir.mkdir(parents=True, exist_ok=True)
                for config_file in config_base_dir.glob("*.config"):
                    shutil.copy2(config_file, config_dir / config_file.name)
                logger.info(f"Copied config files from {config_base_dir} to {config_dir} as fallback")
            except Exception as e2:
                logger.warning(f"Could not copy config files as fallback: {e2}")


def _create_filtered_metadata(config: AmalgkitWorkflowConfig, check: bool, step_results: List) -> None:
    """Create filtered metadata after select step."""
    from metainformant.rna.amalgkit.metadata_filter import filter_selected_metadata
    from metainformant.rna.engine.workflow import WorkflowStepResult

    try:
        try:
            filter_selected_metadata(
                config.work_dir / "metadata" / "metadata.tsv",
                config.work_dir / "metadata" / "metadata_selected.tsv",
                exclude_lite_files=True,
            )
        except ValueError as e:
            if "No samples meet the filtering criteria" in str(e):
                logger.warning(
                    "All selected samples are LITE files. Creating metadata_selected.tsv without LITE filtering."
                )
                filter_selected_metadata(
                    config.work_dir / "metadata" / "metadata.tsv",
                    config.work_dir / "metadata" / "metadata_selected.tsv",
                    exclude_lite_files=False,
                )
            else:
                raise
    except Exception as e:
        logger.error(f"Failed to create filtered metadata: {e}")
        if check:
            step_results.append(
                WorkflowStepResult(step_name="metadata_filtering", return_code=1, success=False, error_message=str(e))
            )


def _validate_getfastq_results(
    config: AmalgkitWorkflowConfig, steps_config: Dict[str, Any], check: bool, step_results: List
) -> None:
    """Validate getfastq results and run fallback extraction if needed."""
    from metainformant.rna.engine.workflow import WorkflowStepResult

    try:
        from metainformant.rna.analysis.validation import save_validation_report, validate_all_samples

        validation_result = validate_all_samples(config, stage="extraction")
        validation_dir = config.work_dir / "validation"
        validation_dir.mkdir(parents=True, exist_ok=True)
        save_validation_report(validation_result, validation_dir / "getfastq_validation.json")

        total = validation_result.get("total_samples", 0)
        validated = validation_result.get("validated", 0)
        failed = validation_result.get("failed", 0)
        logger.info(f"Validation after getfastq: {validated}/{total} samples have FASTQ files extracted")

        if validated == 0 and failed > 0 and total > 0:
            # Try fallback extraction from SRA files
            _attempt_fallback_extraction(config, steps_config, total, step_results, check)
        elif failed > 0:
            logger.warning(f"{failed} samples missing FASTQ files after getfastq step")
    except Exception as e:
        logger.warning(f"Validation after getfastq failed: {e}")


def _attempt_fallback_extraction(
    config: AmalgkitWorkflowConfig,
    steps_config: Dict[str, Any],
    total: int,
    step_results: List,
    check: bool,
) -> None:
    """Attempt fallback SRA extraction when getfastq produces no FASTQ files."""
    from metainformant.rna.engine.workflow import WorkflowStepResult

    sra_dir = config.work_dir / "fastq" / "getfastq"
    if "getfastq" in steps_config and "out_dir" in steps_config["getfastq"]:
        step_out = Path(steps_config["getfastq"]["out_dir"])
        sra_dir = step_out / "getfastq" if step_out.name != "getfastq" else step_out
    elif not sra_dir.exists():
        sra_dir = config.work_dir.parent / "fastq" / "getfastq"

    if sra_dir.exists():
        sra_files_exist = list(sra_dir.glob("*.sra")) + list(sra_dir.glob("sra/*.sra")) + list(sra_dir.glob("*/*.sra"))
        if sra_files_exist:
            logger.warning(f"Validation failed (0 FASTQ) but found {len(sra_files_exist)} SRA files.")
            logger.info("Triggering direct fallback extraction...")
            try:
                extracted = extract_sra_directly(config, sra_dir, sra_dir)
                if extracted > 0:
                    logger.info(f"Fallback extraction recovered {extracted} samples")
                    return
            except Exception as e:
                logger.error(f"Fallback extraction failed: {e}")

    error_msg = (
        f"CRITICAL: getfastq step failed to extract FASTQ files for any samples.\n"
        f"  - Total samples: {total}\n\n"
        f"REMEDIATION STEPS:\n"
        f"  1. Check amalgkit getfastq logs: {config.work_dir / 'logs' / 'getfastq.log'}\n"
        f"  2. Verify SRA files were downloaded\n"
        f"  3. Check if fasterq-dump tools are available\n"
        f"  4. Try re-running with redo: yes if files may be corrupted"
    )
    logger.error(error_msg)
    from metainformant.rna.engine.workflow import WorkflowStepResult

    step_results.append(
        WorkflowStepResult(step_name="getfastq_validation", return_code=1, success=False, error_message=error_msg)
    )


def _validate_quant_results(config: AmalgkitWorkflowConfig, check: bool, step_results: List) -> None:
    """Validate quant results and run post-quant cleanup."""
    try:
        from metainformant.rna.analysis.validation import save_validation_report, validate_all_samples

        validation_result = validate_all_samples(config, stage="quantification")
        validation_dir = config.work_dir / "validation"
        validation_dir.mkdir(parents=True, exist_ok=True)
        save_validation_report(validation_result, validation_dir / "quant_validation.json")

        total = validation_result.get("total_samples", 0)
        validated = validation_result.get("validated", 0)
        failed = validation_result.get("failed", 0)
        logger.info(f"Validation after quant: {validated}/{total} samples quantified")
        if failed > 0:
            logger.warning(f"{failed} samples missing quantification files after quant step")

        # Cleanup FASTQ/SRA files for quantified samples
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


def log_workflow_summary(
    step_results: List[WorkflowStepResult],
    config: AmalgkitWorkflowConfig,
) -> None:
    """Log workflow summary with remediation steps for failed steps.

    Args:
        step_results: List of WorkflowStepResult from execution
        config: Workflow configuration
    """
    successful_steps = sum(1 for sr in step_results if sr.success)
    total_steps = len(step_results)
    failed_steps = total_steps - successful_steps

    logger.info(f"Workflow completed: {successful_steps}/{total_steps} steps successful")

    if failed_steps == 0:
        return

    logger.info("=" * 80)
    logger.info("WORKFLOW SUMMARY")
    logger.info("=" * 80)

    failed_step_results = [sr for sr in step_results if not sr.success]
    logger.info(f"\nFailed steps: {len(failed_step_results)}/{total_steps}")

    for failed_step in failed_step_results:
        logger.info(f"\n  ❌ {failed_step.step_name}")
        if failed_step.error_message:
            for line in failed_step.error_message.split("\n"):
                logger.info(f"     {line}")

    logger.info("\n" + "=" * 80)
    logger.info("REMEDIATION STEPS")
    logger.info("=" * 80)

    _REMEDIATION_TEMPLATES = {
        "getfastq": (
            "    1. Check logs: {log_dir}/getfastq.log\n"
            "    2. Check validation: {work_dir}/validation/getfastq_validation.json\n"
            "    3. Verify SRA files downloaded\n"
            "    4. Check tool availability: which fasterq-dump which fastp"
        ),
        "getfastq_validation": None,  # Same as getfastq
        "integrate": (
            "    1. Ensure getfastq completed successfully\n"
            "    2. Check FASTQ files exist\n"
            "    3. Re-run: --steps getfastq integrate"
        ),
        "quant": (
            "    1. Ensure getfastq and integrate completed successfully\n"
            "    2. Check quantification tools: which kallisto which salmon\n"
            "    3. Check validation: {work_dir}/validation/quant_validation.json"
        ),
        "merge": (
            "    1. Ensure quant completed successfully\n"
            "    2. Check quant files exist\n"
            "    3. Install ggplot2: Rscript -e \"install.packages('ggplot2')\""
        ),
        "curate": (
            "    1. Ensure merge completed successfully\n" "    2. Install R if missing: apt-get install r-base-core"
        ),
    }

    for failed_step in failed_step_results:
        sn = failed_step.step_name
        template = _REMEDIATION_TEMPLATES.get(sn)
        if template is None and sn == "getfastq_validation":
            template = _REMEDIATION_TEMPLATES["getfastq"]
        if template:
            logger.info(f"\n  For {sn}:\n" + template.format(log_dir=config.log_dir, work_dir=config.work_dir))
        else:
            logger.info(
                f"\n  For {sn}:\n"
                f"    1. Check logs: {config.log_dir / f'{sn}.log'}\n"
                f"    2. Review error message above for specific guidance"
            )

    logger.info("\n" + "=" * 80)
