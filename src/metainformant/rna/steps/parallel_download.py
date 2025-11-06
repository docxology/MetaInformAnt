"""Parallel download with sequential quantification for disk space management.

This module implements a producer-consumer pattern:
- Multiple download threads (producers) fetch FASTQ files in parallel
- Single quantification thread (consumer) processes them sequentially
- FASTQ files are deleted immediately after quantification

This maximizes download throughput while preventing disk exhaustion.
"""

from __future__ import annotations

import logging
import queue
import shutil
import threading
from collections.abc import Mapping
from pathlib import Path
from typing import Any

from ...core.io import read_delimited, write_delimited
from ..amalgkit import run_amalgkit
from .download_progress import DownloadProgressMonitor

logger = logging.getLogger(__name__)


def _delete_fastq_for_sample(run_id: str, fastq_dir: Path) -> None:
    """Delete FASTQ files for a specific sample to free disk space."""
    sample_dir = fastq_dir / run_id
    
    if sample_dir.exists() and sample_dir.is_dir():
        try:
            shutil.rmtree(sample_dir)
            logger.info(f"  🗑️  Deleted FASTQ directory: {sample_dir.name}")
        except Exception as e:
            logger.warning(f"Failed to delete {sample_dir}: {e}")
    
    # Also check for loose FASTQ files
    for pattern in [f"{run_id}_*.fastq*", f"{run_id}.fastq*"]:
        for fastq_file in fastq_dir.glob(pattern):
            try:
                fastq_file.unlink()
                logger.info(f"  🗑️  Deleted FASTQ file: {fastq_file.name}")
            except Exception as e:
                logger.warning(f"Failed to delete {fastq_file}: {e}")


def _download_worker(
    download_queue: queue.Queue,
    completion_queue: queue.Queue,
    metadata_path: Path,
    getfastq_params: dict[str, Any],
    work_dir: Path | None,
    log_dir: Path | None,
    worker_id: int,
    progress_monitor: DownloadProgressMonitor | None = None,
):
    """Worker thread that downloads FASTQ files.
    
    Args:
        download_queue: Queue of run IDs to download
        completion_queue: Queue to put completed downloads
        metadata_path: Path to metadata TSV
        getfastq_params: Parameters for amalgkit getfastq
        work_dir: Working directory
        log_dir: Log directory
        worker_id: Unique worker identifier
        progress_monitor: Optional progress monitor for tracking downloads
    """
    logger.info(f"📥 Download worker {worker_id} started")
    
    while True:
        try:
            run_id = download_queue.get(timeout=1)
            if run_id is None:  # Poison pill
                logger.info(f"📥 Download worker {worker_id} shutting down")
                break
            
            logger.info(f"  [{worker_id}] ⬇️  Downloading {run_id}")
            
            # Register with progress monitor if available
            if progress_monitor:
                progress_monitor.register_thread(worker_id, run_id)
            
            # Read metadata and filter to this run
            rows = list(read_delimited(metadata_path, delimiter="\t"))
            single_row = [row for row in rows if row.get("run") == run_id]
            
            if len(single_row) == 0:
                logger.warning(f"  [{worker_id}] ⚠️  {run_id} not in metadata, skipping")
                if progress_monitor:
                    progress_monitor.unregister_thread(worker_id, success=False)
                download_queue.task_done()
                continue
            
            # Create temp metadata for this sample
            temp_metadata = metadata_path.parent / f"metadata.download.{run_id}.tsv"
            write_delimited(single_row, temp_metadata, delimiter="\t")
            
            # Download
            params = getfastq_params.copy()
            params["metadata"] = str(temp_metadata.absolute())  # Use absolute path
            
            result = run_amalgkit(
                "getfastq",
                params,
                work_dir=None,  # Don't set work_dir - using absolute paths
                log_dir=log_dir,
                step_name=f"getfastq_{run_id}",
                check=False,
            )
            
            # Cleanup temp metadata
            try:
                temp_metadata.unlink()
            except Exception:
                pass
            
            success = result.returncode == 0
            
            # Unregister from progress monitor
            if progress_monitor:
                progress_monitor.unregister_thread(worker_id, success=success)
            
            if success:
                logger.info(f"  [{worker_id}] ✅ Downloaded {run_id}")
                # Add to quantification queue
                completion_queue.put(("success", run_id))
            else:
                logger.error(f"  [{worker_id}] ❌ Download failed for {run_id} (code {result.returncode})")
                completion_queue.put(("failed", run_id))
            
            download_queue.task_done()
            
        except queue.Empty:
            continue
        except Exception as e:
            logger.error(f"  [{worker_id}] 💥 Error in download worker: {e}", exc_info=True)
            if progress_monitor:
                try:
                    progress_monitor.unregister_thread(worker_id, success=False)
                except Exception:
                    pass
            download_queue.task_done()


def _quantification_worker(
    completion_queue: queue.Queue,
    metadata_path: Path,
    quant_params: dict[str, Any],
    fastq_dir: Path,
    quant_dir: Path,
    work_dir: Path | None,
    log_dir: Path | None,
    stats: dict[str, Any],
):
    """Worker thread that quantifies samples sequentially."""
    logger.info("🧬 Quantification worker started")
    
    while True:
        try:
            item = completion_queue.get(timeout=1)
            if item is None:  # Poison pill
                logger.info("🧬 Quantification worker shutting down")
                break
            
            status, run_id = item
            
            if status != "success":
                stats["failed"] += 1
                stats["failed_runs"].append(run_id)
                completion_queue.task_done()
                continue
            
            # Check if already quantified
            abundance_file = quant_dir / run_id / "abundance.tsv"
            if abundance_file.exists():
                logger.info(f"  🔬 {run_id} already quantified, skipping quant but deleting FASTQ")
                _delete_fastq_for_sample(run_id, fastq_dir)
                stats["skipped"] += 1
                completion_queue.task_done()
                continue
            
            logger.info(f"  🔬 Quantifying {run_id}")
            
            # Read metadata
            rows = list(read_delimited(metadata_path, delimiter="\t"))
            single_row = [row for row in rows if row.get("run") == run_id]
            
            if len(single_row) == 0:
                logger.warning(f"  ⚠️  {run_id} not in metadata")
                stats["failed"] += 1
                completion_queue.task_done()
                continue
            
            # Create temp metadata
            temp_metadata = metadata_path.parent / f"metadata.quant.{run_id}.tsv"
            write_delimited(single_row, temp_metadata, delimiter="\t")
            
            # Quantify
            params = quant_params.copy()
            params["metadata"] = str(temp_metadata.absolute())  # Use absolute path
            
            result = run_amalgkit(
                "quant",
                params,
                work_dir=None,  # Don't set work_dir - using absolute paths
                log_dir=log_dir,
                step_name=f"quant_{run_id}",
                check=False,
            )
            
            # Cleanup temp metadata
            try:
                temp_metadata.unlink()
            except Exception:
                pass
            
            if result.returncode == 0:
                logger.info(f"  ✅ Quantified {run_id}")
                stats["processed"] += 1
            else:
                logger.error(f"  ❌ Quantification failed for {run_id} (code {result.returncode})")
                stats["failed"] += 1
                stats["failed_runs"].append(run_id)
            
            # Always delete FASTQ to free space
            _delete_fastq_for_sample(run_id, fastq_dir)
            
            completion_queue.task_done()
            
        except queue.Empty:
            continue
        except Exception as e:
            logger.error(f"💥 Error in quantification worker: {e}", exc_info=True)
            completion_queue.task_done()


def run_parallel_download_sequential_quant(
    metadata_path: str | Path,
    getfastq_params: Mapping[str, Any] | None = None,
    quant_params: Mapping[str, Any] | None = None,
    *,
    work_dir: str | Path | None = None,
    log_dir: str | Path | None = None,
    num_download_workers: int = 4,
    max_samples: int | None = None,
    progress_monitor: DownloadProgressMonitor | None = None,
) -> dict[str, Any]:
    """Process samples with parallel downloads and sequential quantification.
    
    Architecture:
    - N download worker threads fetch FASTQ files in parallel
    - 1 quantification worker processes them sequentially
    - FASTQ files are deleted immediately after quantification
    
    Args:
        metadata_path: Path to metadata TSV with sample list
        getfastq_params: Parameters for amalgkit getfastq step
        quant_params: Parameters for amalgkit quant step
        work_dir: Working directory for amalgkit commands
        log_dir: Directory for step logs
        num_download_workers: Number of parallel download threads (default: 4)
        max_samples: Optional limit on number of samples to process
        
    Returns:
        Dictionary with processing statistics
    """
    metadata_path = Path(metadata_path)
    
    if not metadata_path.exists():
        raise FileNotFoundError(f"Metadata file not found: {metadata_path}")
    
    # Sanitize params (remove workflow-only params)
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
    
    # Ensure params use absolute paths
    getfastq_params_dict["out_dir"] = str(fastq_dir)
    quant_params_dict["out_dir"] = str(quant_dir)
    
    fastq_dir.mkdir(parents=True, exist_ok=True)
    quant_dir.mkdir(parents=True, exist_ok=True)
    
    # Get sample list
    rows = list(read_delimited(metadata_path, delimiter="\t"))
    if not rows or "run" not in rows[0]:
        raise ValueError(f"Metadata file missing 'run' column: {metadata_path}")
    
    run_ids = [row.get("run", "").strip() for row in rows]
    run_ids = [r for r in run_ids if r]
    
    if max_samples:
        run_ids = run_ids[:max_samples]
    
    logger.info(f"🚀 Starting parallel download workflow")
    logger.info(f"   Samples: {len(run_ids)}")
    logger.info(f"   Download workers: {num_download_workers}")
    logger.info(f"   Quantification: sequential (1 worker)")
    
    # Initialize progress monitor if not provided
    if progress_monitor is None:
        # Check if progress tracking is enabled in params
        show_progress = getfastq_params_dict.get("show_progress", True)
        if show_progress:
            update_interval = float(getfastq_params_dict.get("progress_update_interval", 2.0))
            use_bars = getfastq_params_dict.get("progress_style", "bar") == "bar"
            progress_monitor = DownloadProgressMonitor(
                out_dir=fastq_dir,
                update_interval=update_interval,
                use_progress_bars=use_bars,
                show_summary=not use_bars,  # Show text summary if not using bars
            )
            progress_monitor.start_monitoring()
    
    # Create queues
    download_queue = queue.Queue()
    completion_queue = queue.Queue()
    
    # Statistics
    stats = {
        "total_samples": len(run_ids),
        "processed": 0,
        "skipped": 0,
        "failed": 0,
        "failed_runs": [],
    }
    
    # Populate download queue
    for run_id in run_ids:
        download_queue.put(run_id)
    
    # Start download workers
    download_threads = []
    for i in range(num_download_workers):
        t = threading.Thread(
            target=_download_worker,
            args=(
                download_queue,
                completion_queue,
                metadata_path,
                getfastq_params_dict,
                work_dir,
                log_dir,
                i + 1,
                progress_monitor,
            ),
            daemon=True,
        )
        t.start()
        download_threads.append(t)
    
    # Start quantification worker
    quant_thread = threading.Thread(
        target=_quantification_worker,
        args=(
            completion_queue,
            metadata_path,
            quant_params_dict,
            fastq_dir,
            quant_dir,
            work_dir,
            log_dir,
            stats,
        ),
        daemon=True,
    )
    quant_thread.start()
    
    # Wait for all downloads to complete
    download_queue.join()
    logger.info("📥 All downloads queued/completed")
    
    # Send poison pills to download workers
    for _ in range(num_download_workers):
        download_queue.put(None)
    
    # Wait for download workers to finish
    for t in download_threads:
        t.join()
    logger.info("📥 All download workers finished")
    
    # Stop progress monitoring
    if progress_monitor:
        progress_monitor.stop_monitoring()
    
    # Wait for all quantifications to complete
    completion_queue.join()
    logger.info("🧬 All quantifications completed")
    
    # Send poison pill to quantification worker
    completion_queue.put(None)
    quant_thread.join()
    logger.info("🧬 Quantification worker finished")
    
    # Final summary
    logger.info("=" * 80)
    logger.info("PARALLEL DOWNLOAD + SEQUENTIAL QUANT COMPLETED")
    logger.info(f"Total samples: {stats['total_samples']}")
    logger.info(f"Successfully processed: {stats['processed']}")
    logger.info(f"Skipped (already done): {stats['skipped']}")
    logger.info(f"Failed: {stats['failed']}")
    logger.info("=" * 80)
    
    return stats


__all__ = ["run_parallel_download_sequential_quant"]

