"""Batched download and quantification for disk space management.

This module implements a simpler batching approach:
1. Download a batch of N samples (in parallel via amalgkit's internal parallelism)
2. Quantify all samples in the batch (using amalgkit's batch processing)
3. Delete all FASTQ files from the batch
4. Repeat with next batch

This is simpler than threading but still provides parallelism through amalgkit's
internal capabilities.
"""

from __future__ import annotations

import logging
import shutil
from collections.abc import Mapping
from pathlib import Path
from typing import Any

from ...core.io import read_delimited, write_delimited
from ..amalgkit import run_amalgkit
from .download_progress import DownloadProgressMonitor

logger = logging.getLogger(__name__)


def _sample_already_quantified(run_id: str, quant_dir: Path) -> bool:
    """Check if a sample has already been quantified.
    
    Args:
        run_id: SRA Run ID
        quant_dir: Base quantification output directory (output/amalgkit/species/quant/)
        
    Returns:
        True if abundance.tsv exists for this sample, False otherwise.
    """
    abundance_file = quant_dir / run_id / "abundance.tsv"
    return abundance_file.exists()


def _delete_fastq_batch(run_ids: list[str], fastq_dir: Path) -> None:
    """Delete FASTQ files for a batch of samples."""
    for run_id in run_ids:
        sample_dir = fastq_dir / run_id
        
        if sample_dir.exists() and sample_dir.is_dir():
            try:
                shutil.rmtree(sample_dir)
                logger.info(f"  🗑️  Deleted FASTQ: {run_id}")
            except Exception as e:
                logger.warning(f"Failed to delete {sample_dir}: {e}")
        
        # Also check for loose FASTQ files  
        for pattern in [f"{run_id}_*.fastq*", f"{run_id}.fastq*"]:
            for fastq_file in fastq_dir.glob(pattern):
                try:
                    fastq_file.unlink()
                except Exception as e:
                    logger.warning(f"Failed to delete {fastq_file}: {e}")


def run_batched_download_quant(
    metadata_path: str | Path,
    getfastq_params: Mapping[str, Any] | None = None,
    quant_params: Mapping[str, Any] | None = None,
    *,
    work_dir: str | Path | None = None,
    log_dir: str | Path | None = None,
    batch_size: int = 8,
    max_samples: int | None = None,
) -> dict[str, Any]:
    """Process samples in batches: download batch → quantify batch → delete FASTQs.
    
    This approach uses amalgkit's internal parallelism for downloading/quantifying
    multiple samples at once, then deletes FASTQs before the next batch.
    
    Args:
        metadata_path: Path to metadata TSV with sample list
        getfastq_params: Parameters for amalgkit getfastq step
        quant_params: Parameters for amalgkit quant step
        work_dir: Working directory for amalgkit commands
        log_dir: Directory for step logs
        batch_size: Number of samples to process per batch (default: 8)
        max_samples: Optional limit on total samples to process
        
    Returns:
        Dictionary with processing statistics
    """
    metadata_path = Path(metadata_path)
    
    if not metadata_path.exists():
        raise FileNotFoundError(f"Metadata file not found: {metadata_path}")
    
    # Sanitize params
    workflow_only_params = {"species_list", "species-list", "accelerate", "genome_dir", "keep_fastq"}
    
    getfastq_params_dict = {
        k: v for k, v in (getfastq_params or {}).items() 
        if k not in workflow_only_params
    }
    quant_params_dict = {
        k: v for k, v in (quant_params or {}).items()
        if k not in workflow_only_params
    }
    
    fastq_dir = Path(getfastq_params_dict.get("out_dir", "output/amalgkit/fastq")).absolute()
    quant_dir = Path(quant_params_dict.get("out_dir", "output/amalgkit/quant")).absolute()
    
    # Update params to use absolute paths
    getfastq_params_dict["out_dir"] = str(fastq_dir)
    quant_params_dict["out_dir"] = str(quant_dir)
    
    fastq_dir.mkdir(parents=True, exist_ok=True)
    quant_dir.mkdir(parents=True, exist_ok=True)
    
    # Get full sample list
    rows = list(read_delimited(metadata_path, delimiter="\t"))
    if not rows or "run" not in rows[0]:
        raise ValueError(f"Metadata file missing 'run' column: {metadata_path}")
    
    all_run_ids = [row.get("run", "").strip() for row in rows]
    all_run_ids = [r for r in all_run_ids if r]
    
    if max_samples:
        all_run_ids = all_run_ids[:max_samples]
    
    # Filter out already-quantified samples BEFORE batching
    logger.info("🔍 Checking for already-quantified samples...")
    samples_to_process = []
    already_quantified = []
    for run_id in all_run_ids:
        if _sample_already_quantified(run_id, quant_dir):
            already_quantified.append(run_id)
        else:
            samples_to_process.append(run_id)
    
    # Statistics
    stats = {
        "total_samples": len(all_run_ids),
        "processed": 0,
        "skipped": len(already_quantified),
        "failed": 0,
        "failed_runs": [],
        "batches": 0,
    }
    
    if already_quantified:
        logger.info(f"⏭️  Skipping {len(already_quantified)} already-quantified samples: {', '.join(already_quantified[:10])}{'...' if len(already_quantified) > 10 else ''}")
    
    if not samples_to_process:
        logger.info("✅ All samples already quantified! Nothing to process.")
        return stats
    
    # Split into batches (only samples that need processing)
    batches = [samples_to_process[i:i + batch_size] for i in range(0, len(samples_to_process), batch_size)]
    
    logger.info("=" * 80)
    logger.info(f"🚀 BATCHED DOWNLOAD-QUANT-DELETE WORKFLOW")
    logger.info(f"   Total samples: {len(all_run_ids)}")
    logger.info(f"   Already quantified (skipped): {len(already_quantified)}")
    logger.info(f"   Samples to process: {len(samples_to_process)}")
    
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
        logger.info("📊 Progress tracking enabled for batched downloads")
    logger.info(f"   Batch size: {batch_size}")
    logger.info(f"   Number of batches: {len(batches)}")
    logger.info("=" * 80)
    
    for batch_num, batch_run_ids in enumerate(batches, 1):
        logger.info("")
        logger.info("=" * 80)
        logger.info(f"📦 BATCH {batch_num}/{len(batches)}: {len(batch_run_ids)} samples")
        logger.info("=" * 80)
        
        # All samples in this batch need processing (already filtered)
        samples_to_process = batch_run_ids
        
        logger.info(f"  📋 Processing {len(samples_to_process)} samples (total skipped: {stats['skipped']})")
        
        # Create batch metadata file (only for samples that need processing)
        batch_rows = [row for row in rows if row.get("run") in samples_to_process]
        batch_metadata = metadata_path.parent / f"metadata.batch{batch_num}.tsv"
        write_delimited(batch_rows, batch_metadata, delimiter="\t")
        
        try:
            # Step 1: Download this batch
            logger.info(f"  ⬇️  Downloading batch {batch_num} ({len(batch_run_ids)} samples)...")
            
            # Register batch samples with progress monitor
            if progress_monitor:
                for idx, run_id in enumerate(batch_run_ids, 1):
                    progress_monitor.register_thread(batch_num * 1000 + idx, run_id)  # Use batch-based thread IDs
            
            download_params = getfastq_params_dict.copy()
            download_params["metadata"] = str(batch_metadata.absolute())  # Use absolute path
            
            download_result = run_amalgkit(
                "getfastq",
                download_params,
                work_dir=None,  # Don't set work_dir when using absolute paths
                log_dir=log_dir,
                step_name=f"getfastq_batch{batch_num}",
                check=False,
            )
            
            # Unregister batch samples from progress monitor
            if progress_monitor:
                for idx, run_id in enumerate(batch_run_ids, 1):
                    # Check if sample was downloaded
                    sample_dir = fastq_dir / run_id
                    success = sample_dir.exists() or any(fastq_dir.glob(f"{run_id}*.fastq*"))
                    progress_monitor.unregister_thread(batch_num * 1000 + idx, success=success)
            
            if download_result.returncode != 0:
                logger.warning(f"  ⚠️  Download batch {batch_num} had errors (code {download_result.returncode})")
                # Continue anyway - some samples might have succeeded
            else:
                logger.info(f"  ✅ Download batch {batch_num} completed")
            
            # Check which samples actually downloaded
            downloaded = []
            for run_id in batch_run_ids:
                # Check if FASTQ exists
                sample_dir = fastq_dir / run_id
                if sample_dir.exists() or any(fastq_dir.glob(f"{run_id}*.fastq*")):
                    downloaded.append(run_id)
            
            logger.info(f"  📊 {len(downloaded)}/{len(batch_run_ids)} samples downloaded successfully")
            
            if not downloaded:
                logger.warning(f"  ⚠️  No samples downloaded in batch {batch_num}, skipping quantification")
                stats["failed"] += len(batch_run_ids)
                stats["failed_runs"].extend(batch_run_ids)
                # Clean up batch metadata
                try:
                    batch_metadata.unlink()
                except Exception:
                    pass
                continue
            
            # Step 2: Quantify downloaded samples
            logger.info(f"  🧬 Quantifying batch {batch_num} ({len(downloaded)} samples)...")
            
            quant_params_copy = quant_params_dict.copy()
            quant_params_copy["metadata"] = str(batch_metadata.absolute())  # Use absolute path
            
            quant_result = run_amalgkit(
                "quant",
                quant_params_copy,
                work_dir=None,  # Don't set work_dir when using absolute paths
                log_dir=log_dir,
                step_name=f"quant_batch{batch_num}",
                check=False,
            )
            
            if quant_result.returncode != 0:
                logger.warning(f"  ⚠️  Quantification batch {batch_num} had errors (code {quant_result.returncode})")
            else:
                logger.info(f"  ✅ Quantification batch {batch_num} completed")
            
            # Check which samples were actually quantified
            quantified = []
            for run_id in downloaded:
                if _sample_already_quantified(run_id, quant_dir):
                    quantified.append(run_id)
                else:
                    stats["failed"] += 1
                    stats["failed_runs"].append(run_id)
            
            logger.info(f"  📊 {len(quantified)}/{len(downloaded)} samples quantified successfully")
            stats["processed"] += len(quantified)
            
            # Step 3: Delete FASTQ files for this batch
            logger.info(f"  🗑️  Deleting FASTQ files from batch {batch_num}...")
            _delete_fastq_batch(batch_run_ids, fastq_dir)
            logger.info(f"  ✅ Cleaned up batch {batch_num}")
            
            stats["batches"] += 1
            
        except Exception as e:
            logger.error(f"  ❌ Error processing batch {batch_num}: {e}", exc_info=True)
            stats["failed"] += len(batch_run_ids)
            stats["failed_runs"].extend(batch_run_ids)
        
        finally:
            # Clean up batch metadata
            try:
                batch_metadata.unlink()
            except Exception:
                pass
    
    # Stop progress monitoring
    if progress_monitor:
        progress_monitor.stop_monitoring()
    
    # Final summary
    logger.info("")
    logger.info("=" * 80)
    logger.info("🎉 BATCHED WORKFLOW COMPLETED")
    logger.info(f"   Total samples: {stats['total_samples']}")
    logger.info(f"   Successfully processed: {stats['processed']}")
    logger.info(f"   Failed: {stats['failed']}")
    logger.info(f"   Batches completed: {stats['batches']}/{len(batches)}")
    logger.info("=" * 80)
    
    if stats["failed_runs"]:
        logger.warning(f"Failed runs ({len(stats['failed_runs'])}): {', '.join(stats['failed_runs'][:10])}" + 
                      (f"... and {len(stats['failed_runs']) - 10} more" if len(stats['failed_runs']) > 10 else ""))
    
    return stats


__all__ = ["run_batched_download_quant"]

