"""Workflow orchestration functions for RNA-seq analysis.

This module provides functions to orchestrate complete amalgkit workflows
for single species, including status checking, cleanup, and step execution.
"""

from __future__ import annotations

import logging
import subprocess
import time
from collections.abc import Sequence
from datetime import datetime
from pathlib import Path
from typing import Any

from ..core.io import read_delimited
from ..core.logging import get_logger
from .monitoring import (
    analyze_species_status,
    check_active_downloads,
    check_workflow_progress,
    count_quantified_samples,
    find_unquantified_samples,
)
from .steps import delete_sample_fastqs, quantify_sample
from .workflow import AmalgkitWorkflowConfig, execute_workflow, load_workflow_config

logger = get_logger(__name__)

# Valid amalgkit steps in execution order
VALID_STEPS = [
    "metadata",
    "integrate",
    "config",
    "select",
    "getfastq",
    "quant",
    "merge",
    "cstmm",
    "curate",
    "csca",
    "sanity",
]


def discover_species_configs(config_dir: Path) -> dict[str, dict[str, Any]]:
    """Discover all species configs and extract metadata.

    Args:
        config_dir: Directory containing amalgkit config files

    Returns:
        Dictionary mapping species_id -> {'name': str, 'total': int}
    """
    import yaml

    species_info = {}
    for config_file in sorted(config_dir.glob("amalgkit_*.yaml")):
        if "template" in config_file.stem.lower() or "test" in config_file.stem.lower():
            continue

        species_id = config_file.stem.replace("amalgkit_", "")

        try:
            config = yaml.safe_load(config_file.read_text())
            species_list = config.get("species_list", [])
            if species_list:
                name = species_list[0].replace("_", " ")
            else:
                name = species_id.replace("_", " ").title()

            # Try to get total from metadata if available
            work_dir = Path(config.get("work_dir", ""))
            if work_dir.exists():
                metadata_file = work_dir / "metadata" / "metadata.tsv"
                if not metadata_file.exists():
                    metadata_file = work_dir / "metadata" / "metadata.filtered.tissue.tsv"

                if metadata_file.exists():
                    rows = list(read_delimited(metadata_file, delimiter="\t"))
                    total = len([r for r in rows if r.get("run")])
                else:
                    total = 0
            else:
                total = 0

            species_info[species_id] = {"name": name, "total": total}
        except Exception:
            species_info[species_id] = {"name": species_id.replace("_", " ").title(), "total": 0}

    return species_info


def run_workflow_for_species(
    config_path: Path,
    steps: Sequence[str] | None = None,
    *,
    check: bool = False,
) -> dict[str, Any]:
    """Run workflow steps for a single species.

    Args:
        config_path: Path to species workflow config file
        steps: List of steps to run (default: all steps)
        check: If True, stop on first failure

    Returns:
        Dictionary with 'success', 'completed', 'failed' keys
    """
    logger.info(f"Running workflow for {config_path.name}")

    # Load config
    cfg = load_workflow_config(config_path)

    # Determine steps to run
    if steps is None:
        # Run all steps via execute_workflow
        return_codes = execute_workflow(cfg, check=check)
        steps_run = VALID_STEPS[: len(return_codes)]
        success = all(rc == 0 or rc == 204 for rc in return_codes)  # 204 = skipped
        completed = [step for step, rc in zip(steps_run, return_codes) if rc == 0]
        failed = [step for step, rc in zip(steps_run, return_codes) if rc != 0 and rc != 204]

        return {
            "success": success,
            "completed": completed,
            "failed": failed,
            "return_codes": return_codes,
        }
    else:
        # Run specific steps
        results = {"success": True, "completed": [], "failed": []}

        for step in steps:
            if step not in VALID_STEPS:
                logger.warning(f"Invalid step: {step}")
                continue

            success = _run_single_step(cfg, step)
            if success:
                results["completed"].append(step)
            else:
                results["failed"].append(step)
                results["success"] = False
                if check:
                    break

        return results


def _log_metadata_sample_count(work_dir: Path, logger: logging.Logger) -> None:
    """Log the number of samples found in metadata files after metadata step.
    
    Args:
        work_dir: Work directory where metadata files are created
        logger: Logger instance for output
    """
    # Check for metadata files in order of preference
    metadata_paths = [
        work_dir / "metadata" / "metadata.tsv",
        work_dir / "metadata" / "metadata.filtered.tissue.tsv",
        work_dir / "metadata" / "metadata.filtered.clean.tsv",
    ]
    
    for metadata_file in metadata_paths:
        if not metadata_file.exists():
            continue
        
        try:
            rows = list(read_delimited(metadata_file, delimiter="\t"))
            if not rows:
                continue
            
            # Check if this file has a 'run' column (sample-level metadata)
            if "run" in rows[0]:
                sample_count = len([r for r in rows if r.get("run")])
                logger.info(f"📊 Found {sample_count} samples in {metadata_file.name}")
                return
            
            # If no 'run' column, this might be a pivot table - log that
            logger.debug(f"Metadata file {metadata_file.name} exists but has no 'run' column")
        except Exception as e:
            logger.debug(f"Error reading {metadata_file.name}: {e}")
            continue
    
    # If we get here, no suitable metadata file was found
    logger.warning("⚠️  Metadata step completed but no sample metadata file found")


def _run_single_step(config: AmalgkitWorkflowConfig, step: str) -> bool:
    """Run a single amalgkit step.

    Args:
        config: Workflow configuration
        step: Step name to run

    Returns:
        True if step completed successfully
    """
    from metainformant.rna.workflow import _apply_step_defaults
    from metainformant.rna import steps as _steps_mod
    
    logger.info(f"Running {step}...")

    # Apply step defaults first (this modifies config.per_step in place)
    _apply_step_defaults(config)

    # Get step parameters (now with defaults applied)
    step_params = config.per_step.get(step, {})
    work_dir = config.work_dir
    log_dir = config.log_dir or (work_dir.parent / "logs")

    # Inject metadata for steps that need it if not provided
    if step in {"select", "getfastq", "integrate", "quant", "merge"} and not step_params.get("metadata"):
        # All steps need row-per-sample format (NOT pivot tables)
        # Pivot tables (pivot_selected.tsv, pivot_qualified.tsv) are summary tables, not suitable for processing
        candidates = [
            work_dir / "metadata" / "metadata.filtered.tissue.tsv",
            work_dir / "metadata" / "metadata.tsv",
        ]
        for candidate in candidates:
            if candidate.exists():
                step_params["metadata"] = str(candidate)
                break
        # If none found, use default
        if "metadata" not in step_params:
            step_params["metadata"] = str(work_dir / "metadata" / "metadata.tsv")
    
    # Inject config_dir for select step if not provided
    if step == "select" and not step_params.get("config_dir"):
        preferred = work_dir / "config_base"
        fallback = work_dir / "config"
        step_params["config_dir"] = str(preferred if preferred.exists() else fallback)

    # Use step runner if available (preferred - handles all parameter normalization)
    runner = _steps_mod.STEP_RUNNERS.get(step)
    if runner is not None:
        try:
            result = runner(
                step_params,
                work_dir=work_dir,
                log_dir=log_dir,
                check=False,
            )
            if result.returncode == 0:
                logger.info(f"✅ {step} completed successfully")
                return True
            else:
                logger.warning(f"⚠️  {step} failed with code {result.returncode}")
                return False
        except Exception as e:
            logger.error(f"Error running {step} via step runner: {e}", exc_info=True)
            # Fall through to manual command building

    # Fallback: Build command manually (for steps without runners or if runner fails)
    cmd = ["amalgkit", step]

    # Add common parameters
    # Only add --threads for steps that support it
    steps_with_threads = {"getfastq", "quant", "merge", "cstmm", "csca", "curate"}
    if config.threads and step in steps_with_threads:
        cmd.extend(["--threads", str(config.threads)])

    # Add step-specific parameters
    # Special handling for boolean parameters that need "yes"/"no" values
    boolean_yes_no_params = {"redo", "aws", "gcp", "ncbi", "pfd", "fastp", "keep_fastq", "build_index"}
    
    # Parameters that need absolute paths (amalgkit runs from work_dir, so relative paths fail)
    path_params = {"out_dir", "fastq_dir", "index_dir", "index-dir", "out", "dest_dir", "config_dir", "metadata"}
    
    for key, value in step_params.items():
        if value is None:
            continue
        
        # Handle boolean False
        if value is False:
            continue
        
        # Convert relative paths to absolute for path parameters
        if key in path_params and isinstance(value, (str, Path)):
            path_value = Path(value)
            if not path_value.is_absolute():
                # Resolve relative to repo root (assume paths are relative to repo root)
                try:
                    from metainformant.core.paths import get_repo_root
                    repo_root = get_repo_root()
                    path_value = (repo_root / path_value).resolve()
                except Exception:
                    # Fallback: resolve relative to current working directory
                    path_value = Path(value).resolve()
            cmd.extend([f"--{key}", str(path_value)])
        # Handle boolean True or string "yes" for yes/no params
        elif key in boolean_yes_no_params:
            if value is True or (isinstance(value, str) and value.lower() in ("yes", "true", "1")):
                cmd.extend([f"--{key}", "yes"])
            elif isinstance(value, str) and value.lower() in ("no", "false", "0"):
                cmd.extend([f"--{key}", "no"])
            else:
                # Already a string value like "yes" or "no"
                cmd.extend([f"--{key}", str(value)])
        elif value is True:
            # Regular boolean flag
            cmd.append(f"--{key}")
        else:
            # Regular value
            cmd.extend([f"--{key}", str(value)])

    # Special handling for integrate step: skip if no FASTQ files exist yet
    if step == "integrate":
        fastq_dir_param = step_params.get("fastq_dir")
        if fastq_dir_param:
            fastq_dir = Path(fastq_dir_param)
            if not fastq_dir.is_absolute():
                try:
                    from metainformant.core.paths import get_repo_root
                    repo_root = get_repo_root()
                    fastq_dir = (repo_root / fastq_dir).resolve()
                except Exception:
                    fastq_dir = Path(fastq_dir_param).resolve()
            
            # Check if fastq directory exists and has files
            if not fastq_dir.exists() or not any(fastq_dir.rglob("*.fastq*")):
                logger.info(f"⚠️  {step}: No FASTQ files found yet, skipping (will run after getfastq)")
                return True  # Skip gracefully

    # Run
    log_dir.mkdir(parents=True, exist_ok=True)
    log_file = log_dir / f'{step}_{datetime.now().strftime("%Y%m%d_%H%M%S")}.log'

    try:
        with open(log_file, "w") as log_f:
            result = subprocess.run(
                cmd,
                stdout=log_f,
                stderr=subprocess.STDOUT,
                cwd=str(work_dir),
                timeout=7200,  # 2 hour timeout
            )

        if result.returncode == 0:
            logger.info(f"✅ {step} completed successfully")
            
            # Special handling for metadata step: count samples found
            if step == "metadata":
                _log_metadata_sample_count(work_dir, logger)
            
            return True
        else:
            logger.warning(f"⚠️  {step} failed with code {result.returncode}")
            logger.warning(f"Log: {log_file}")
            return False
    except subprocess.TimeoutExpired:
        logger.error(f"❌ {step} timed out after 2 hours")
        return False
    except Exception as e:
        logger.error(f"❌ {step} error: {e}")
        return False


def check_workflow_status(
    config_path: Path,
    detailed: bool = False,
) -> dict[str, Any]:
    """Check workflow status for a species.

    Args:
        config_path: Path to species workflow config file
        detailed: If True, include detailed sample categories

    Returns:
        Dictionary with status information
    """
    if detailed:
        return analyze_species_status(config_path)
    else:
        return check_workflow_progress(config_path)


def cleanup_unquantified_samples(
    config_path: Path,
    *,
    log_dir: Path | None = None,
) -> tuple[int, int]:
    """Quantify downloaded samples and cleanup FASTQs.

    Args:
        config_path: Path to species workflow config file
        log_dir: Optional log directory

    Returns:
        Tuple of (quantified_count, failed_count)
    """
    cfg = load_workflow_config(config_path)

    # Find unquantified samples
    unquantified = find_unquantified_samples(config_path)
    if not unquantified:
        logger.info("No downloaded unquantified samples")
        return 0, 0

    logger.info(f"Found {len(unquantified)} downloaded but unquantified samples")

    # Get paths
    fastq_dir = Path(cfg.per_step.get("getfastq", {}).get("out_dir", cfg.work_dir / "fastq"))
    quant_dir = Path(cfg.per_step.get("quant", {}).get("out_dir", cfg.work_dir / "quant"))

    # Read metadata
    metadata_file = cfg.work_dir / "metadata" / "metadata.tsv"
    if not metadata_file.exists():
        metadata_file = cfg.work_dir / "metadata" / "metadata.filtered.tissue.tsv"

    if not metadata_file.exists():
        logger.error("No metadata file found")
        return 0, len(unquantified)

    rows = list(read_delimited(metadata_file, delimiter="\t"))

    # Get quant params
    quant_params = dict(cfg.per_step.get("quant", {}))
    quant_params["out_dir"] = str(quant_dir.absolute())
    quant_params["threads"] = cfg.threads or 12

    # Inject index_dir if needed
    if "index_dir" not in quant_params and "index-dir" not in quant_params:
        index_dir = quant_dir.parent / "work" / "index"
        if index_dir.exists():
            quant_params["index_dir"] = str(index_dir.absolute())

    quantified_count = 0
    failed_count = 0

    # Process each sample
    for sample_id in unquantified:
        try:
            sample_rows = [row for row in rows if row.get("run") == sample_id]

            if not sample_rows:
                logger.warning(f"⚠️  {sample_id} not in metadata, skipping")
                failed_count += 1
                continue

            # Quantify
            success, message, abundance_file = quantify_sample(
                sample_id=sample_id,
                metadata_rows=sample_rows,
                quant_params=quant_params,
                log_dir=log_dir or (cfg.log_dir or (cfg.work_dir / "logs")),
                step_name=f"quant_{sample_id}",
            )

            if success:
                logger.info(f"✅ Quantified: {sample_id}")
                quantified_count += 1
            else:
                logger.warning(f"❌ Failed: {sample_id} - {message}")
                failed_count += 1

            # Always cleanup FASTQs
            delete_sample_fastqs(sample_id, fastq_dir)
            logger.info(f"🗑️  Cleaned: {sample_id}")

        except Exception as e:
            logger.warning(f"❌ Error: {sample_id} - {e}")
            failed_count += 1
            # Attempt cleanup
            try:
                delete_sample_fastqs(sample_id, fastq_dir)
            except Exception:
                pass

    logger.info(f"Summary: {quantified_count} quantified, {failed_count} failed")
    return quantified_count, failed_count


def monitor_workflows(
    species_configs: dict[str, Path],
    watch_interval: int = 60,
) -> None:
    """Real-time monitoring of workflows.

    Args:
        species_configs: Dictionary mapping species_id -> config_path
        watch_interval: Update interval in seconds
    """
    try:
        while True:
            # Clear screen
            print("\033[2J\033[H", end="")

            print("\n" + "=" * 80)
            print("  MULTI-SPECIES RNA-SEQ WORKFLOW MONITOR")
            print("=" * 80)
            print(f"\n⏰ {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")

            total_samples = 0
            total_quantified = 0
            active_downloads = check_active_downloads()

            for species_id, config_path in species_configs.items():
                if not config_path.exists():
                    continue

                progress = check_workflow_progress(config_path)
                total_samples += progress["total"]
                total_quantified += progress["quantified"]

                print(f"📊 {species_id}")
                print(f"   Progress: {progress['quantified']}/{progress['total']} ({progress['percentage']:.1f}%)")
                print()

            print("=" * 80)
            overall_percent = (total_quantified * 100) // total_samples if total_samples > 0 else 0
            print(f"🎯 OVERALL: {total_quantified}/{total_samples} samples ({overall_percent}%)")
            print("=" * 80)

            remaining = total_samples - total_quantified
            if remaining > 0:
                estimated_hours = (remaining * 7.5) / 60
                print(f"\n⏳ Estimated remaining: {estimated_hours:.1f} hours ({remaining} samples)")

            print(f"\nPress Ctrl+C to stop monitoring (updates every {watch_interval}s)")

            time.sleep(watch_interval)
    except KeyboardInterrupt:
        print("\n\nMonitoring stopped.")

