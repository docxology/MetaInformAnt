"""Unified sample processing: download → quantify → delete FASTQ.

This module provides a single configurable function that handles both sequential
and parallel processing modes for RNA-seq workflows.

Processing Modes:
- Sequential (num_workers=1): Process one sample at a time
  - Download → Quantify → Delete FASTQ → Next sample
  - Maximum disk efficiency: only one sample's FASTQs exist at a time
  
- Parallel (num_workers>1): Multiple downloads, sequential quantification
  - N download workers fetch FASTQ files in parallel
  - 1 quantification worker processes them sequentially
  - FASTQ files deleted immediately after quantification
  - Maximizes throughput while preventing disk exhaustion

Both modes ensure disk usage remains bounded regardless of cohort size.
"""

from __future__ import annotations

import os
import queue
import shutil
import threading
import time
from collections.abc import Mapping
from pathlib import Path
from typing import Any

from ...core.io import read_delimited, write_delimited
from ...core.logging import get_logger
from ..amalgkit import run_amalgkit
from .download_progress import DownloadProgressMonitor

logger = get_logger(__name__)


def _get_sample_list(metadata_path: Path) -> list[str]:
    """Extract list of run IDs from metadata table.
    
    Args:
        metadata_path: Path to metadata TSV file (should contain 'run' column)
        
    Returns:
        List of run IDs (SRA accessions)
        
    Raises:
        ValueError: If metadata file is empty or missing 'run' column
    """
    try:
        rows = list(read_delimited(metadata_path, delimiter="\t"))
        
        if not rows:
            raise ValueError(f"Metadata file {metadata_path} is empty")
        
        if "run" not in rows[0]:
            available_cols = list(rows[0].keys()) if rows else []
            raise ValueError(
                f"Metadata file {metadata_path} missing 'run' column. "
                f"Available columns: {available_cols}"
            )
        
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
    abundance_file = quant_dir / run_id / "abundance.tsv"
    exists = abundance_file.exists()
    
    if exists:
        logger.info(f"Sample {run_id} already quantified, skipping")
    
    return exists


def _delete_fastq_for_sample(run_id: str, fastq_dir: Path) -> None:
    """Delete FASTQ files for a specific sample to free disk space.
    
    amalgkit getfastq places files in fastq_dir/getfastq/<SRR>/ by default.
    This function checks both locations for compatibility.
    
    Args:
        run_id: SRA run accession
        fastq_dir: Directory where FASTQ files are stored
    """
    # Check standard amalgkit location: fastq_dir/getfastq/<SRR>/
    sample_dir_getfastq = fastq_dir / "getfastq" / run_id
    if sample_dir_getfastq.exists() and sample_dir_getfastq.is_dir():
        try:
            shutil.rmtree(sample_dir_getfastq)
            logger.info(f"  🗑️  Deleted FASTQ directory: {sample_dir_getfastq.name}")
        except Exception as e:
            logger.warning(f"Failed to delete {sample_dir_getfastq}: {e}")
    
    # Also check direct location: fastq_dir/<SRR>/
    sample_dir = fastq_dir / run_id
    if sample_dir.exists() and sample_dir.is_dir():
        try:
            shutil.rmtree(sample_dir)
            logger.info(f"  🗑️  Deleted FASTQ directory: {sample_dir.name}")
        except Exception as e:
            logger.warning(f"Failed to delete {sample_dir}: {e}")
    
    # Also check for loose FASTQ files in getfastq subdirectory
    getfastq_subdir = fastq_dir / "getfastq"
    if getfastq_subdir.exists():
        for pattern in [f"{run_id}_*.fastq*", f"{run_id}.fastq*"]:
            for fastq_file in getfastq_subdir.glob(pattern):
                try:
                    fastq_file.unlink()
                    logger.info(f"  🗑️  Deleted FASTQ file: {fastq_file.name}")
                except Exception as e:
                    logger.warning(f"Failed to delete {fastq_file}: {e}")
    
    # Check for loose FASTQ files in root fastq_dir
    for pattern in [f"{run_id}_*.fastq*", f"{run_id}.fastq*"]:
        for fastq_file in fastq_dir.glob(pattern):
            try:
                fastq_file.unlink()
                logger.info(f"  🗑️  Deleted FASTQ file: {fastq_file.name}")
            except Exception as e:
                logger.warning(f"Failed to delete {fastq_file}: {e}")


def _wait_for_fastq_files(run_id: str, fastq_dir: Path, max_wait_minutes: int = 60) -> bool:
    """Wait for FASTQ files to exist for a sample.
    
    amalgkit getfastq may return success before fasterq-dump completes conversion.
    This function actively polls for FASTQ files to ensure they exist before
    attempting quantification.
    
    Args:
        run_id: Sample ID
        fastq_dir: Directory where FASTQ files should be
        max_wait_minutes: Maximum time to wait in minutes
        
    Returns:
        True if FASTQ files found, False if timeout
    """
    sample_dir_getfastq = fastq_dir / "getfastq" / run_id
    sample_dir_direct = fastq_dir / run_id
    
    patterns = [
        f"{run_id}_1.fastq.gz",
        f"{run_id}_2.fastq.gz",
        f"{run_id}_1.fastq",
        f"{run_id}_2.fastq",
        f"{run_id}.fastq.gz",
        f"{run_id}.fastq",
    ]
    
    start_time = time.time()
    max_wait_seconds = max_wait_minutes * 60
    check_interval = 10  # Check every 10 seconds
    
    while time.time() - start_time < max_wait_seconds:
        # Check getfastq subdirectory first (standard location)
        if sample_dir_getfastq.exists():
            for pattern in patterns:
                if (sample_dir_getfastq / pattern).exists():
                    return True
        
        # Check direct location
        if sample_dir_direct.exists():
            for pattern in patterns:
                if (sample_dir_direct / pattern).exists():
                    return True
        
        # Also check for any FASTQ files in the directories
        if sample_dir_getfastq.exists():
            fastq_files = list(sample_dir_getfastq.glob("*.fastq*"))
            if fastq_files:
                return True
        
        if sample_dir_direct.exists():
            fastq_files = list(sample_dir_direct.glob("*.fastq*"))
            if fastq_files:
                return True
        
        time.sleep(check_interval)
    
    return False


def _download_worker(
    download_queue: queue.Queue,
    completion_queue: queue.Queue,
    metadata_path: Path,
    getfastq_params: dict[str, Any],
    work_dir: Path | None,
    log_dir: Path | None,
    worker_id: int,
    progress_monitor: DownloadProgressMonitor | None = None,
) -> None:
    """Worker thread that downloads FASTQ files (parallel mode only).
    
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
            
            # Set PATH to include wrapper for fasterq-dump with --size-check off
            # This prevents "disk-limit exceeded" errors during extraction
            fastq_dir = Path(getfastq_params.get("out_dir", work_dir / "fastq" if work_dir else Path("fastq")))
            if not fastq_dir.is_absolute():
                fastq_dir = fastq_dir.resolve()
            wrapper_dir = fastq_dir / "temp"
            wrapper_path = wrapper_dir / "fasterq-dump"
            
            # Set up environment variables for external drive usage
            env = os.environ.copy()
            
            # Set temp directories to use external drive (not /tmp)
            if work_dir:
                # Find repo root by walking up from work_dir
                work_dir_path = Path(work_dir).resolve()
                repo_root = work_dir_path
                max_levels = 8  # Allow more levels for deep work_dir structures
                for _ in range(max_levels):
                    # Check for .git or pyproject.toml first (more reliable markers)
                    if (repo_root / ".git").exists() or (repo_root / "pyproject.toml").exists():
                        break
                    # Only use .cursorrules if we haven't found better markers
                    if (repo_root / ".cursorrules").exists() and "output" not in str(repo_root):
                        break
                    if repo_root == repo_root.parent:
                        # Reached filesystem root, use system temp directory as fallback
                        # This handles test environments where repo detection fails
                        import tempfile
                        repo_root = Path(tempfile.gettempdir()) / "metainformant_repo"
                        break
                    repo_root = repo_root.parent
                
                general_temp_dir = repo_root / "output" / ".tmp"
                sra_temp_dir = fastq_dir / "temp" / "sra"
                
                general_temp_dir.mkdir(parents=True, exist_ok=True)
                sra_temp_dir.mkdir(parents=True, exist_ok=True)
                
                env["TMPDIR"] = str(general_temp_dir)
                env["TEMP"] = str(general_temp_dir)
                env["TMP"] = str(general_temp_dir)
                env["NCBI_SRA_REPOSITORY"] = str(sra_temp_dir)
            
            if wrapper_path.exists() and wrapper_path.is_file():
                # Prepend wrapper directory to PATH so amalgkit finds our wrapper
                env["PATH"] = f"{wrapper_dir}:{env.get('PATH', '')}"
            
            result = run_amalgkit(
                "getfastq",
                params,
                work_dir=None,  # Don't set work_dir - using absolute paths
                log_dir=log_dir,
                step_name=f"getfastq_{run_id}",
                env=env,
                check=False,
            )
            
            # Cleanup temp metadata
            try:
                temp_metadata.unlink()
            except Exception:
                pass
            
            success = result.returncode == 0
            
            # Validate that FASTQ files actually exist before reporting success
            # getfastq may return success even when it produces 0 reads
            if success:
                # Extract fastq_dir from getfastq_params
                out_dir_str = getfastq_params.get("out_dir", "")
                if out_dir_str:
                    fastq_dir = Path(out_dir_str)
                    if not fastq_dir.is_absolute():
                        fastq_dir = fastq_dir.resolve()
                elif work_dir:
                    fastq_dir = (work_dir / "fastq").resolve()
                else:
                    # Fallback: use params dict which should have absolute path
                    fastq_dir = Path(params.get("out_dir", "fastq")).resolve()
                
                # Quick check for FASTQ files (don't wait long - conversion should be done)
                sample_dir_getfastq = fastq_dir / "getfastq" / run_id
                sample_dir_direct = fastq_dir / run_id
                has_fastq = False
                
                if sample_dir_getfastq.exists():
                    has_fastq = any(sample_dir_getfastq.glob("*.fastq*"))
                if not has_fastq and sample_dir_direct.exists():
                    has_fastq = any(sample_dir_direct.glob("*.fastq*"))
                
                if not has_fastq:
                    # Check if SRA exists (conversion may have failed)
                    has_sra = False
                    sra_file = None
                    if sample_dir_getfastq.exists():
                        sra_files = list(sample_dir_getfastq.glob("*.sra"))
                        has_sra = len(sra_files) > 0
                        if sra_files:
                            sra_file = sra_files[0]
                    if not has_sra and sample_dir_direct.exists():
                        sra_files = list(sample_dir_direct.glob("*.sra"))
                        has_sra = len(sra_files) > 0
                        if sra_files:
                            sra_file = sra_files[0]
                    
                    if has_sra and sra_file:
                        # Attempt automatic conversion
                        logger.warning(f"  [{worker_id}] ⚠️  {run_id}: SRA exists but no FASTQ - attempting automatic conversion")
                        from metainformant.rna.steps.getfastq import convert_sra_to_fastq
                        
                        conversion_success, conversion_msg, fastq_files = convert_sra_to_fastq(
                            run_id,
                            sra_file,
                            sra_file.parent,
                            threads=getfastq_params.get("threads", 24),
                            log_dir=log_dir,
                        )
                        
                        if conversion_success and fastq_files:
                            logger.info(f"  [{worker_id}] ✅ {run_id}: Successfully converted SRA to FASTQ ({len(fastq_files)} files)")
                            success = True
                        else:
                            logger.warning(f"  [{worker_id}] ⚠️  {run_id}: Automatic conversion failed: {conversion_msg}")
                            success = False
                    else:
                        logger.warning(f"  [{worker_id}] ⚠️  {run_id}: getfastq succeeded but no files found (may have produced 0 reads)")
                        success = False
            
            # Unregister from progress monitor
            if progress_monitor:
                progress_monitor.unregister_thread(worker_id, success=success)
            
            if success:
                logger.info(f"  [{worker_id}] ✅ Downloaded {run_id}")
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
) -> None:
    """Worker thread that quantifies samples sequentially.
    
    Used in parallel mode to process downloaded samples one at a time.
    
    Args:
        completion_queue: Queue of completed downloads to quantify
        metadata_path: Path to metadata TSV
        quant_params: Parameters for amalgkit quant
        fastq_dir: Directory containing FASTQ files
        quant_dir: Directory for quantification outputs
        work_dir: Working directory
        log_dir: Log directory
        stats: Statistics dictionary to update
    """
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
            
            # Wait for FASTQ files to exist (amalgkit getfastq may return success before conversion completes)
            logger.info(f"  🔬 Waiting for FASTQ files for {run_id}...")
            if not _wait_for_fastq_files(run_id, fastq_dir, max_wait_minutes=60):
                logger.error(f"  ❌ FASTQ files not found for {run_id} after waiting (conversion may have failed)")
                stats["failed"] += 1
                stats["failed_runs"].append(run_id)
                completion_queue.task_done()
                continue
            
            logger.info(f"  🔬 FASTQ files found, quantifying {run_id}")
            
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
            # Note: amalgkit quant automatically finds FASTQ files based on metadata and standard directory structure
            # It looks in fastq_dir/getfastq/<SRR>/ or fastq_dir/<SRR>/ relative to the work_dir
            
            # Inject index_dir if not already set and index exists
            if "index_dir" not in params and "index-dir" not in params:
                if work_dir:
                    index_dir = Path(work_dir) / "index"
                else:
                    index_dir = quant_dir.parent / "work" / "index"
                if index_dir.exists():
                    params["index_dir"] = str(index_dir.absolute())
                    logger.debug(f"Injected index_dir: {params['index_dir']}")
            
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


def run_download_quant_workflow(
    metadata_path: str | Path,
    getfastq_params: Mapping[str, Any] | None = None,
    quant_params: Mapping[str, Any] | None = None,
    *,
    work_dir: str | Path | None = None,
    log_dir: str | Path | None = None,
    num_workers: int = 1,
    max_samples: int | None = None,
    skip_completed: bool = True,
    progress_monitor: DownloadProgressMonitor | None = None,
) -> dict[str, Any]:
    """Process samples: download → quantify → delete FASTQ.
    
    Unified function that handles both sequential and parallel processing modes.
    
    This function automatically configures environment variables to use external
    drive for temporary files (not /tmp) to avoid disk space issues.
    
    Processing Modes:
    - Sequential (num_workers=1): Process one sample at a time
      - Download → Quantify → Delete FASTQ → Next sample
      - Maximum disk efficiency: only one sample's FASTQs exist at a time
      
    - Parallel (num_workers>1): Multiple downloads, sequential quantification
      - N download workers fetch FASTQ files in parallel
      - 1 quantification worker processes them sequentially
      - FASTQ files deleted immediately after quantification
      - Maximizes throughput while preventing disk exhaustion
    
    Args:
        metadata_path: Path to metadata TSV with sample list
        getfastq_params: Parameters for amalgkit getfastq step
        quant_params: Parameters for amalgkit quant step
        work_dir: Working directory for amalgkit commands
        log_dir: Directory for step logs
        num_workers: Number of parallel download workers (default: 1 for sequential)
        max_samples: Optional limit on number of samples to process
        skip_completed: If True, skip samples that are already quantified (sequential mode only)
        progress_monitor: Optional progress monitor (if None, will be created if enabled in params)
        
    Returns:
        Dictionary with processing statistics:
        - total_samples: Total number of samples
        - processed: Number of samples successfully processed
        - skipped: Number of samples skipped (already done)
        - failed: Number of samples that failed
        - failed_runs: List of run IDs that failed
        
    Raises:
        FileNotFoundError: If metadata file not found
        ValueError: If metadata file missing 'run' column
    """
    metadata_path = Path(metadata_path)
    
    if not metadata_path.exists():
        raise FileNotFoundError(f"Metadata file not found: {metadata_path}")
    
    # Set up environment variables for external drive usage BEFORE processing
    if work_dir:
        work_dir_path = Path(work_dir).resolve()
        # Find repo root using same logic as workflow.py
        # Prefer .git and pyproject.toml over .cursorrules (which might exist in output/)
        repo_root = work_dir_path
        max_levels = 20  # Allow more levels to reach filesystem root
        for _ in range(max_levels):
            # Check for .git or pyproject.toml first (more reliable markers)
            if (repo_root / ".git").exists() or (repo_root / "pyproject.toml").exists():
                break
            # Only use .cursorrules if we haven't found better markers
            if (repo_root / ".cursorrules").exists() and "output" not in str(repo_root):
                break
            if repo_root == repo_root.parent:
                # Reached filesystem root, use system temp directory as fallback
                # This handles test environments where repo detection fails
                import tempfile
                repo_root = Path(tempfile.gettempdir()) / "metainformant_repo"
                break
            repo_root = repo_root.parent
        else:
            # If we exhausted max_levels without finding markers or reaching root
            import tempfile
            repo_root = Path(tempfile.gettempdir()) / "metainformant_repo"

        general_temp_dir = repo_root / "output" / ".tmp"
        fastq_dir = Path(getfastq_params.get("out_dir", work_dir_path / "fastq") if getfastq_params else work_dir_path / "fastq")
        if not fastq_dir.is_absolute():
            fastq_dir = fastq_dir.resolve()
        sra_temp_dir = fastq_dir / "temp" / "sra"
        
        general_temp_dir.mkdir(parents=True, exist_ok=True)
        sra_temp_dir.mkdir(parents=True, exist_ok=True)
        
        # Set environment variables for this process and all subprocesses
        os.environ["TMPDIR"] = str(general_temp_dir)
        os.environ["TEMP"] = str(general_temp_dir)
        os.environ["TMP"] = str(general_temp_dir)
        os.environ["NCBI_SRA_REPOSITORY"] = str(sra_temp_dir)
        
        logger.info(f"🔧 Configured temp directories for external drive:")
        logger.info(f"   TMPDIR: {general_temp_dir}")
        logger.info(f"   NCBI_SRA_REPOSITORY: {sra_temp_dir}")
    
    # Sanitize params (remove workflow-only params)
    workflow_only_params = {
        "species_list",
        "species-list",
        "accelerate",
        "genome_dir",
        "keep_fastq",
        "num_download_workers",
        "parallel_workers",
    }
    
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
    
    # Statistics
    stats = {
        "total_samples": len(run_ids),
        "processed": 0,
        "skipped": 0,
        "failed": 0,
        "failed_runs": [],
    }
    
    # Choose processing mode based on num_workers
    if num_workers > 1:
        # PARALLEL MODE: Multiple download workers, sequential quantification
        logger.info(f"🚀 Starting parallel download workflow")
        logger.info(f"   Samples: {len(run_ids)}")
        logger.info(f"   Download workers: {num_workers}")
        logger.info(f"   Quantification: sequential (1 worker)")
        logger.info(f"   {num_workers} samples download in parallel → quantify sequentially → delete FASTQs")
        
        # Initialize progress monitor if not provided
        if progress_monitor is None:
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
        
        # Create queues
        download_queue = queue.Queue()
        completion_queue = queue.Queue()
        
        # Populate download queue
        for run_id in run_ids:
            download_queue.put(run_id)
        
        # Start download workers
        download_threads = []
        for i in range(num_workers):
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
        for _ in range(num_workers):
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
        
    else:
        # SEQUENTIAL MODE: One sample at a time
        logger.info(f"🚀 Starting sequential processing")
        logger.info(f"   Samples: {len(run_ids)}")
        logger.info(f"   Each sample: download → immediately quantify → immediately delete FASTQs → next sample")
        logger.info(f"   Maximum disk efficiency: only one sample's FASTQs exist at a time")
        
        # Initialize progress monitor if enabled
        if progress_monitor is None:
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
                logger.info(f"  ✓ {run_id} already quantified, skipping")
                stats["skipped"] += 1
                continue
            
            # Check if FASTQ files already exist (from previous download)
            sample_dir_getfastq = fastq_dir / "getfastq" / run_id
            sample_dir_direct = fastq_dir / run_id
            has_fastq = False
            has_sra = False
            
            if sample_dir_getfastq.exists():
                has_fastq = any(sample_dir_getfastq.glob("*.fastq*"))
                has_sra = any(sample_dir_getfastq.glob("*.sra"))
            if not has_fastq and sample_dir_direct.exists():
                has_fastq = has_fastq or any(sample_dir_direct.glob("*.fastq*"))
                has_sra = has_sra or any(sample_dir_direct.glob("*.sra"))
            
            skip_download = has_fastq or (has_sra and _wait_for_fastq_files(run_id, fastq_dir, max_wait_minutes=0.1))
            
            try:
                if skip_download:
                    logger.info(f"  ✓ FASTQ files already exist for {run_id}, skipping download")
                else:
                    # Step 1: Download FASTQ for this sample only
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
                    write_delimited(single_row, temp_metadata, delimiter="\t")
                    
                    getfastq_single = getfastq_params_dict.copy()
                    getfastq_single["metadata"] = str(temp_metadata.absolute())
                    
                    result_download = run_amalgkit(
                        "getfastq",
                        getfastq_single,
                        work_dir=None,  # Don't set work_dir - using absolute paths
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
                    
                    # Check for 0-read LITE file issue in output
                    # amalgkit may return success even when it produces 0 reads (LITE files)
                    output_text = (result_download.stdout or "") + (result_download.stderr or "")
                    if "fastp input reads: 0 bp" in output_text or "fastp output reads: 0 bp" in output_text:
                        if "Individual fastp input reads (bp): 0" in output_text or "Individual fastp output reads (bp): 0" in output_text:
                            logger.error(f"  ❌ {run_id}: Download produced 0 reads (LITE file - metadata only, no sequence data)")
                            stats["failed"] += 1
                            stats["failed_runs"].append(run_id)
                            _delete_fastq_for_sample(run_id, fastq_dir)
                            continue
                
                # Validate that FASTQ files actually exist before proceeding
                # getfastq may return success even when it produces 0 reads
                logger.info(f"  ✓ Checking FASTQ files for {run_id}...")
                
                # Wait for FASTQ files to exist (with shorter timeout for validation)
                if not _wait_for_fastq_files(run_id, fastq_dir, max_wait_minutes=5):
                    # Check if SRA file exists but no FASTQ (conversion failure)
                    sample_dir_getfastq = fastq_dir / "getfastq" / run_id
                    sample_dir_direct = fastq_dir / run_id
                    has_sra = False
                    if sample_dir_getfastq.exists():
                        has_sra = any(sample_dir_getfastq.glob("*.sra"))
                    if not has_sra and sample_dir_direct.exists():
                        has_sra = any(sample_dir_direct.glob("*.sra"))
                    
                    if has_sra:
                        # Attempt automatic conversion
                        logger.warning(f"  ⚠️  {run_id}: SRA file exists but FASTQ conversion failed - attempting automatic conversion")
                        from metainformant.rna.steps.getfastq import convert_sra_to_fastq
                        
                        sra_file = None
                        if sample_dir_getfastq.exists():
                            sra_files = list(sample_dir_getfastq.glob("*.sra"))
                            if sra_files:
                                sra_file = sra_files[0]
                        if not sra_file and sample_dir_direct.exists():
                            sra_files = list(sample_dir_direct.glob("*.sra"))
                            if sra_files:
                                sra_file = sra_files[0]
                        
                        if sra_file:
                            conversion_success, conversion_msg, fastq_files = convert_sra_to_fastq(
                                run_id,
                                sra_file,
                                sra_file.parent,
                                threads=quant_params_dict.get("threads", 24),
                                log_dir=log_dir,
                            )
                            
                            if conversion_success and fastq_files:
                                logger.info(f"  ✅ {run_id}: Successfully converted SRA to FASTQ ({len(fastq_files)} files)")
                                # Continue with quantification below
                            else:
                                logger.error(f"  ❌ {run_id}: Automatic conversion failed: {conversion_msg}")
                                stats["failed"] += 1
                                stats["failed_runs"].append(run_id)
                                _delete_fastq_for_sample(run_id, fastq_dir)
                                continue
                        else:
                            logger.error(f"  ❌ {run_id}: SRA file exists but FASTQ conversion failed (possibly 0 reads or conversion error)")
                            stats["failed"] += 1
                            stats["failed_runs"].append(run_id)
                            _delete_fastq_for_sample(run_id, fastq_dir)
                            continue
                    else:
                        logger.error(f"  ❌ {run_id}: No FASTQ files found and no SRA file (download may have produced 0 reads)")
                        stats["failed"] += 1
                        stats["failed_runs"].append(run_id)
                        _delete_fastq_for_sample(run_id, fastq_dir)
                        continue
                
                logger.info(f"  ✓ Downloaded FASTQ for {run_id}")
                
                # Step 2: Quantify this sample
                logger.info(f"  → Quantifying {run_id}")
                
                # Create temp metadata for quantification
                temp_metadata_quant = metadata_path.parent / f"metadata.quant.{run_id}.tsv"
                write_delimited(single_row, temp_metadata_quant, delimiter="\t")
                
                quant_single = quant_params_dict.copy()
                quant_single["metadata"] = str(temp_metadata_quant.absolute())
                
                # Inject index_dir if not already set and index exists
                if "index_dir" not in quant_single and "index-dir" not in quant_single:
                    if work_dir:
                        index_dir = Path(work_dir) / "index"
                    else:
                        index_dir = quant_dir.parent / "work" / "index"
                    if index_dir.exists():
                        quant_single["index_dir"] = str(index_dir.absolute())
                        logger.debug(f"Injected index_dir: {quant_single['index_dir']}")
                
                result_quant = run_amalgkit(
                    "quant",
                    quant_single,
                    work_dir=None,  # Don't set work_dir - using absolute paths
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
    if num_workers > 1:
        logger.info("PARALLEL DOWNLOAD + SEQUENTIAL QUANT COMPLETED")
    else:
        logger.info("SEQUENTIAL PROCESSING COMPLETED")
    logger.info(f"Total samples: {stats['total_samples']}")
    logger.info(f"Successfully processed: {stats['processed']}")
    logger.info(f"Skipped (already done): {stats['skipped']}")
    logger.info(f"Failed: {stats['failed']}")
    logger.info("=" * 80)
    
    if stats["failed_runs"]:
        logger.warning(
            f"Failed runs: {', '.join(stats['failed_runs'][:10])}" + 
            (f" ... and {len(stats['failed_runs']) - 10} more" if len(stats['failed_runs']) > 10 else "")
        )
    
    return stats


__all__ = ["run_download_quant_workflow"]

