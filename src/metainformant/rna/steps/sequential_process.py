"""Sequential sample processing: download → quantify → delete FASTQ → repeat.

This module implements disk-space-aware orchestration for RNA-seq workflows.
Instead of downloading all samples then quantifying all, we process one sample
at a time to prevent disk exhaustion.

Workflow:
1. Read metadata table to get sample list
2. For each sample (run ID):
   a. Download FASTQ files (getfastq for single sample)
   b. Quantify that sample (quant for single sample)
   c. Delete FASTQ files for that sample
   d. Skip if quantification already exists (resume capability)
3. Continue to next sample

This approach ensures disk usage remains bounded regardless of cohort size.
"""

from __future__ import annotations

import logging
import shutil
import subprocess
from collections.abc import Mapping
from pathlib import Path
from typing import Any

from ...core.io import read_delimited
from ..amalgkit import run_amalgkit
from .download_progress import DownloadProgressMonitor

logger = logging.getLogger(__name__)


def _get_sample_list(metadata_path: Path) -> list[str]:
    """Extract list of run IDs from metadata table.
    
    Args:
        metadata_path: Path to metadata TSV file (should contain 'run' column)
        
    Returns:
        List of run IDs (SRA accessions)
    """
    try:
        # read_delimited returns a generator, convert to list, explicitly use TSV
        rows = list(read_delimited(metadata_path, delimiter="\t"))
        
        if not rows:
            raise ValueError(f"Metadata file {metadata_path} is empty")
        
        # Check if 'run' column exists
        if "run" not in rows[0]:
            available_cols = list(rows[0].keys()) if rows else []
            raise ValueError(f"Metadata file {metadata_path} missing 'run' column. Available columns: {available_cols}")
        
        # Extract run IDs
        runs = [row.get("run", "") for row in rows]
        runs = [r.strip() for r in runs if r and str(r).strip()]
        
        logger.info(f"Found {len(runs)} samples to process from {metadata_path}")
        return runs
        
    except Exception as e:
        logger.error(f"Failed to read metadata from {metadata_path}: {e}")
        raise


def _sample_already_quantified(run_id: str, quant_dir: Path) -> bool:
    """Check if a sample has already been quantified.
    
    Args:
        run_id: SRA run accession
        quant_dir: Directory where quantification outputs are stored
        
    Returns:
        True if abundance.tsv exists for this sample
    """
    # Amalgkit typically creates run_id/abundance.tsv
    abundance_file = quant_dir / run_id / "abundance.tsv"
    exists = abundance_file.exists()
    
    if exists:
        logger.info(f"Sample {run_id} already quantified, skipping")
    
    return exists


def _delete_fastq_for_sample(run_id: str, fastq_dir: Path) -> None:
    """Delete FASTQ files for a specific sample to free disk space.
    
    Args:
        run_id: SRA run accession
        fastq_dir: Directory where FASTQ files are stored
    """
    # Amalgkit typically creates run_id/ subdirectory or run_id_*.fastq.gz files
    sample_dir = fastq_dir / run_id
    
    if sample_dir.exists() and sample_dir.is_dir():
        try:
            shutil.rmtree(sample_dir)
            logger.info(f"Deleted FASTQ directory for {run_id}: {sample_dir}")
        except Exception as e:
            logger.warning(f"Failed to delete {sample_dir}: {e}")
    
    # Also check for loose FASTQ files
    for pattern in [f"{run_id}_*.fastq*", f"{run_id}.fastq*"]:
        for fastq_file in fastq_dir.glob(pattern):
            try:
                fastq_file.unlink()
                logger.info(f"Deleted FASTQ file: {fastq_file}")
            except Exception as e:
                logger.warning(f"Failed to delete {fastq_file}: {e}")


def run_sequential_download_quant(
    metadata_path: str | Path,
    getfastq_params: Mapping[str, Any] | None = None,
    quant_params: Mapping[str, Any] | None = None,
    *,
    work_dir: str | Path | None = None,
    log_dir: str | Path | None = None,
    max_samples: int | None = None,
    skip_completed: bool = True,
) -> dict[str, Any]:
    """Process samples sequentially: download → quantify → delete FASTQ.
    
    Args:
        metadata_path: Path to metadata TSV with sample list
        getfastq_params: Parameters for amalgkit getfastq step
        quant_params: Parameters for amalgkit quant step
        work_dir: Working directory for amalgkit commands
        log_dir: Directory for step logs
        max_samples: Optional limit on number of samples to process (for testing)
        skip_completed: If True, skip samples that are already quantified
        
    Returns:
        Dictionary with processing statistics:
        - total_samples: Total number of samples
        - processed: Number of samples processed
        - skipped: Number of samples skipped (already done)
        - failed: Number of samples that failed
        - failed_runs: List of run IDs that failed
    """
    metadata_path = Path(metadata_path)
    
    if not metadata_path.exists():
        raise FileNotFoundError(f"Metadata file not found: {metadata_path}")
    
    # Get sample list
    run_ids = _get_sample_list(metadata_path)
    
    if max_samples is not None:
        run_ids = run_ids[:max_samples]
        logger.info(f"Limited to first {max_samples} samples")
    
    # Extract output directories from params and sanitize
    # Remove workflow-level params that shouldn't be passed to amalgkit steps
    workflow_only_params = {"species_list", "species-list", "accelerate", "genome_dir", "keep_fastq"}
    
    getfastq_params_dict = {
        k: v for k, v in (getfastq_params or {}).items() 
        if k not in workflow_only_params
    }
    quant_params_dict = {
        k: v for k, v in (quant_params or {}).items()
        if k not in workflow_only_params
    }
    
    fastq_dir = Path(getfastq_params_dict.get("out_dir", "output/amalgkit/fastq"))
    quant_dir = Path(quant_params_dict.get("out_dir", "output/amalgkit/quant"))
    
    fastq_dir.mkdir(parents=True, exist_ok=True)
    quant_dir.mkdir(parents=True, exist_ok=True)
    
    # Statistics
    stats = {
        "total_samples": len(run_ids),
        "processed": 0,
        "skipped": 0,
        "failed": 0,
        "failed_runs": [],
    }
    
    logger.info(f"Starting sequential processing of {len(run_ids)} samples")
    logger.info(f"FASTQ dir: {fastq_dir}")
    logger.info(f"Quant dir: {quant_dir}")
    
    # Initialize progress monitor if enabled
    progress_monitor: DownloadProgressMonitor | None = None
    show_progress = getfastq_params_dict.get("show_progress", True)
    if show_progress:
        update_interval = float(getfastq_params_dict.get("progress_update_interval", 2.0))
        use_bars = getfastq_params_dict.get("progress_style", "bar") == "bar"
        progress_monitor = DownloadProgressMonitor(
            out_dir=fastq_dir,
            update_interval=update_interval,
            use_progress_bars=use_bars,
            show_summary=not use_bars,
        )
        progress_monitor.start_monitoring()
        logger.info("📊 Progress tracking enabled for sequential downloads")
    
    for idx, run_id in enumerate(run_ids, 1):
        logger.info(f"[{idx}/{len(run_ids)}] Processing sample: {run_id}")
        
        # Check if already quantified
        if skip_completed and _sample_already_quantified(run_id, quant_dir):
            stats["skipped"] += 1
            continue
        
        try:
            # Step 1: Download FASTQ for this sample only
            # Create a temporary single-sample metadata file
            logger.info(f"  → Downloading FASTQ for {run_id} ({idx}/{len(run_ids)})")
            
            # Register with progress monitor
            if progress_monitor:
                progress_monitor.register_thread(1, run_id)  # Use thread_id=1 for sequential
            
            # Read full metadata and filter to this run only
            rows = list(read_delimited(metadata_path, delimiter="\t"))
            single_row = [row for row in rows if row.get("run") == run_id]
            
            if len(single_row) == 0:
                logger.warning(f"  ⚠ Run {run_id} not found in metadata, skipping")
                if progress_monitor:
                    progress_monitor.unregister_thread(1, success=False)
                stats["failed"] += 1
                stats["failed_runs"].append(run_id)
                continue
            
            # Write single-sample metadata
            temp_metadata = metadata_path.parent / f"metadata.single.{run_id}.tsv"
            from ...core.io import write_delimited
            write_delimited(single_row, temp_metadata, delimiter="\t")
            
            getfastq_single = getfastq_params_dict.copy()
            getfastq_single["metadata"] = str(temp_metadata)
            
            result_download = run_amalgkit(
                "getfastq",
                getfastq_single,
                work_dir=work_dir,
                log_dir=log_dir,
                step_name=f"getfastq_{run_id}",
                check=False,
            )
            
            # Clean up temp metadata
            try:
                temp_metadata.unlink()
            except Exception:
                pass
            
            # Unregister from progress monitor
            success = result_download.returncode == 0
            if progress_monitor:
                progress_monitor.unregister_thread(1, success=success)
            
            if not success:
                logger.error(f"  ✗ Download failed for {run_id} (code {result_download.returncode})")
                stats["failed"] += 1
                stats["failed_runs"].append(run_id)
                continue
            
            logger.info(f"  ✓ Downloaded FASTQ for {run_id}")
            
            # Step 2: Quantify this sample
            logger.info(f"  → Quantifying {run_id}")
            
            # Create temp metadata for quantification
            temp_metadata_quant = metadata_path.parent / f"metadata.quant.{run_id}.tsv"
            write_delimited(single_row, temp_metadata_quant, delimiter="\t")
            
            quant_single = quant_params_dict.copy()
            quant_single["metadata"] = str(temp_metadata_quant)
            
            result_quant = run_amalgkit(
                "quant",
                quant_single,
                work_dir=work_dir,
                log_dir=log_dir,
                step_name=f"quant_{run_id}",
                check=False,
            )
            
            # Clean up temp metadata
            try:
                temp_metadata_quant.unlink()
            except Exception:
                pass
            
            if result_quant.returncode != 0:
                logger.error(f"  ✗ Quantification failed for {run_id} (code {result_quant.returncode})")
                stats["failed"] += 1
                stats["failed_runs"].append(run_id)
                # Still delete FASTQ to free space
                _delete_fastq_for_sample(run_id, fastq_dir)
                continue
            
            logger.info(f"  ✓ Quantified {run_id}")
            
            # Step 3: Delete FASTQ files to free disk space
            logger.info(f"  → Deleting FASTQ for {run_id}")
            _delete_fastq_for_sample(run_id, fastq_dir)
            logger.info(f"  ✓ Freed disk space for {run_id}")
            
            stats["processed"] += 1
            logger.info(f"[{idx}/{len(run_ids)}] ✓ Completed {run_id}")
            
        except Exception as e:
            logger.error(f"  ✗ Unexpected error processing {run_id}: {e}", exc_info=True)
            if progress_monitor:
                try:
                    progress_monitor.unregister_thread(1, success=False)
                except Exception:
                    pass
            stats["failed"] += 1
            stats["failed_runs"].append(run_id)
            # Attempt cleanup
            try:
                _delete_fastq_for_sample(run_id, fastq_dir)
            except Exception:
                pass
    
    # Stop progress monitoring
    if progress_monitor:
        progress_monitor.stop_monitoring()
    
    # Final summary
    logger.info("=" * 80)
    logger.info("Sequential processing completed")
    logger.info(f"Total samples: {stats['total_samples']}")
    logger.info(f"Successfully processed: {stats['processed']}")
    logger.info(f"Skipped (already done): {stats['skipped']}")
    logger.info(f"Failed: {stats['failed']}")
    
    if stats["failed_runs"]:
        logger.warning(f"Failed runs: {', '.join(stats['failed_runs'][:10])}" + 
                      (f" ... and {len(stats['failed_runs']) - 10} more" if len(stats['failed_runs']) > 10 else ""))
    
    return stats


__all__ = ["run_sequential_download_quant"]
