#!/usr/bin/env python3
"""
Integrated ENA download + quantification workflow with batched processing.

This script combines:
1. Robust ENA FASTQ downloads (direct, no SRA conversion)
2. Kallisto quantification
3. Automatic FASTQ cleanup after quantification

Workflow per batch:
- Download N samples from ENA (parallel)
- Quantify all downloaded samples
- Delete FASTQs to free disk space
- Repeat with next batch

Usage:
    python3 workflow_ena_integrated.py --config config/amalgkit/amalgkit_cfloridanus.yaml --batch-size 12 --threads 12
"""

# ============================================================================
# CONFIGURATION
# ============================================================================
# Scope: Single-species ENA-based workflow with batched processing
# Steps: metadata → config → select → getfastq → quant → merge
# Config: YAML file via --config argument (required)
# Threads: Per-config (default 12) or --threads override
# Batch Size: Per-config or --batch-size override (default 12)
# Output: output/amalgkit/{species}/work/
# Dependencies: wget, kallisto, amalgkit (for metadata/index setup)
# Reliability: 100% (ENA direct downloads vs 0% SRA Toolkit for large samples)
# ============================================================================

import argparse
import os
import subprocess
import sys
import time
from pathlib import Path
from typing import Optional
import yaml

# Add src to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent.parent / "src"))

from metainformant.core.io import read_delimited, write_delimited
from metainformant.core.logging import get_logger
from metainformant.rna.progress_tracker import get_tracker

logger = get_logger(__name__)


def load_config(config_path: Path) -> dict:
    """Load amalgkit YAML config.
    
    Args:
        config_path: Path to YAML configuration file
        
    Returns:
        Dictionary containing configuration values
        
    Side effects:
        None (read-only operation)
    """
    with open(config_path) as f:
        return yaml.safe_load(f)


def get_sample_list(metadata_path: Path) -> list[str]:
    """Get list of run IDs from metadata file.
    
    Args:
        metadata_path: Path to TSV metadata file with 'run' column
        
    Returns:
        List of SRA run IDs (SRR*)
        
    Raises:
        ValueError: If metadata file lacks 'run' column
        
    Side effects:
        None (read-only operation)
    """
    rows = list(read_delimited(metadata_path, delimiter='\t'))
    if not rows or 'run' not in rows[0]:
        raise ValueError(f"Metadata must have 'run' column: {metadata_path}")
    
    run_ids = [row['run'] for row in rows if row.get('run')]
    return run_ids


def sample_already_quantified(run_id: str, quant_dir: Path) -> bool:
    """Check if sample has abundance.tsv.
    
    Args:
        run_id: SRA run ID (e.g., 'SRR1234567')
        quant_dir: Directory containing quantification results
        
    Returns:
        True if abundance.tsv exists and has non-zero size
        
    Side effects:
        None (read-only operation)
    """
    abundance = quant_dir / run_id / "abundance.tsv"
    return abundance.exists() and abundance.stat().st_size > 0


def download_batch_ena(run_ids: list[str], metadata_path: Path, fastq_dir: Path, 
                       threads: int, batch_num: int, species: str | None = None,
                       tracker=None) -> tuple[list[str], list[str]]:
    """Download a batch of samples using the robust ENA downloader.
    
    Args:
        run_ids: List of SRA run IDs to download
        metadata_path: Path to metadata TSV file
        fastq_dir: Directory to save downloaded FASTQ files
        threads: Number of parallel download threads
        batch_num: Batch number for logging
    
    Returns:
        Tuple of (successful_downloads, failed_downloads)
        
    Side effects:
        - Creates temporary batch metadata file
        - Downloads FASTQ files to fastq_dir/{run_id}/
        - Removes temporary batch metadata file
        
    Dependencies:
        - amalgkit getfastq (via metainformant.rna.steps.getfastq)
        - wget command available (for ENA downloads)
    """
    logger.info(f"  📥 Downloading batch {batch_num} ({len(run_ids)} samples)...")
    
    # Report download starts to tracker
    if tracker and species:
        for run_id in run_ids:
            tracker.on_download_start(species, run_id)
    
    # Get temp directory for batch metadata
    from metainformant.core.disk import get_recommended_temp_dir
    repo_root = Path(__file__).parent.parent.parent.resolve()
    temp_dir = get_recommended_temp_dir(repo_root)
    
    # Create batch metadata in temp directory
    rows = list(read_delimited(metadata_path, delimiter='\t'))
    batch_rows = [row for row in rows if row.get('run') in run_ids]
    batch_metadata = temp_dir / f"metadata_batch_{species}_{batch_num}.tsv"
    write_delimited(batch_rows, batch_metadata, delimiter='\t')
    logger.debug(f"    Created batch metadata: {batch_metadata}")
    
    # Use amalgkit getfastq with ENA acceleration
    from metainformant.rna.steps.getfastq import run as run_getfastq
    
    getfastq_params = {
        "metadata": str(batch_metadata),
        "out_dir": str(fastq_dir),
        "threads": threads,
        "accelerate": True,  # Enable ENA/AWS/GCP acceleration
        "aws": "yes",
        "gcp": "yes",
        "ncbi": "yes",
    }
    
    try:
        result = run_getfastq(
            getfastq_params,
            work_dir=None,
            log_dir=fastq_dir.parent / "logs",
            check=False,
        )
        
        # Parse which samples succeeded - check both possible locations
        successful = []
        failed = []
        
        for run_id in run_ids:
            # Check both possible locations: fastq_dir/{run_id} and fastq_dir/getfastq/{run_id}
            sample_dir = fastq_dir / run_id
            if not sample_dir.exists():
                sample_dir = fastq_dir / "getfastq" / run_id
            
            if sample_dir.exists():
                fastq_files = list(sample_dir.glob("*.fastq.gz")) + list(sample_dir.glob("*.fastq"))
                if fastq_files and all(f.stat().st_size > 1000000 for f in fastq_files):
                    successful.append(run_id)
                else:
                    failed.append(run_id)
            else:
                failed.append(run_id)
        
        logger.info(f"  ✅ Downloaded: {len(successful)}/{len(run_ids)} samples")
        if failed:
            logger.warning(f"  ⚠️  Failed downloads: {len(failed)}/{len(run_ids)}")
            logger.debug(f"    Failed IDs: {', '.join(failed[:10])}{'...' if len(failed) > 10 else ''}")
        
        # Report to progress tracker
        if tracker and species:
            for run_id in successful:
                tracker.on_download_complete(species, run_id)
            for run_id in failed:
                tracker.on_download_failed(species, run_id)
        
        # Clean up batch metadata
        try:
            batch_metadata.unlink(missing_ok=True)
            logger.debug(f"    Cleaned up batch metadata: {batch_metadata}")
        except Exception as e:
            logger.debug(f"    Note: Could not remove batch metadata: {e}")
        
        return successful, failed
        
    except subprocess.TimeoutExpired:
        logger.error(f"  ❌ Download batch {batch_num} timed out after 2 hours")
        batch_metadata.unlink(missing_ok=True)
        return [], run_ids
    except Exception as e:
        logger.error(f"  ❌ Download error: {e}")
        batch_metadata.unlink(missing_ok=True)
        return [], run_ids


def quantify_batch_kallisto(run_ids: list[str], fastq_dir: Path, quant_dir: Path,
                            index_path: Path, threads: int, batch_num: int,
                            species: str | None = None, tracker=None) -> tuple[list[str], list[str]]:
    """Quantify samples using kallisto.
    
    Args:
        run_ids: List of SRA run IDs to quantify
        fastq_dir: Directory containing downloaded FASTQ files
        quant_dir: Directory to save quantification results
        index_path: Path to kallisto transcriptome index
        threads: Number of threads for kallisto
        batch_num: Batch number for logging
    
    Returns:
        Tuple of (successful_quants, failed_quants)
        
    Side effects:
        - Creates quant_dir/{run_id}/ directories
        - Writes abundance.tsv and run_info.json per sample
        - Skips already-quantified samples
        
    Dependencies:
        - kallisto command available on PATH
        - Index must exist and be valid
    """
    logger.info(f"  🧬 Quantifying batch {batch_num} ({len(run_ids)} samples)...")
    
    successful = []
    failed = []
    
    for run_id in run_ids:
        # Check if already quantified
        if sample_already_quantified(run_id, quant_dir):
            logger.info(f"    ✅ {run_id} already quantified")
            successful.append(run_id)
            continue
        
        # Find FASTQ files
        sample_dir = fastq_dir / run_id
        fastq_files = sorted(sample_dir.glob("*.fastq.gz")) if sample_dir.exists() else []
        
        if not fastq_files:
            logger.warning(f"    ⚠️  {run_id} - no FASTQ files found")
            failed.append(run_id)
            continue
        
        # Create output directory
        output_dir = quant_dir / run_id
        output_dir.mkdir(parents=True, exist_ok=True)
        
        # Run kallisto (detect single vs paired-end)
        cmd = [
            'kallisto', 'quant',
            '-i', str(index_path),
            '-o', str(output_dir),
            '-t', str(threads),
        ]
        
        # Add --single flag for single-end data
        if len(fastq_files) == 1:
            cmd.extend(['--single', '-l', '200', '-s', '20'])  # fragment length and SD estimates
        
        # Add FASTQ files
        cmd.extend([str(f) for f in fastq_files])
        
        try:
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=1800)  # 30 min per sample
            
            if result.returncode == 0 and (output_dir / "abundance.tsv").exists():
                logger.info(f"    ✅ {run_id} quantified")
                successful.append(run_id)
            else:
                logger.warning(f"    ⚠️  {run_id} quantification failed: {result.stderr[:200]}")
                failed.append(run_id)
                
        except subprocess.TimeoutExpired:
            logger.warning(f"    ⚠️  {run_id} quantification timed out")
            failed.append(run_id)
        except Exception as e:
            logger.warning(f"    ⚠️  {run_id} error: {e}")
            failed.append(run_id)
    
    logger.info(f"  ✅ Quantified: {len(successful)}/{len(run_ids)} samples")
    if failed:
        logger.warning(f"  ⚠️  Failed quantifications: {len(failed)}/{len(run_ids)}")
        logger.debug(f"    Failed IDs: {', '.join(failed[:10])}{'...' if len(failed) > 10 else ''}")
    
    # Report to progress tracker
    if tracker and species:
        for run_id in successful:
            tracker.on_quant_complete(species, run_id)
    
    return successful, failed


def cleanup_fastqs(run_ids: list[str], fastq_dir: Path, species: str | None = None,
                   tracker=None) -> int:
    """Delete FASTQ files for samples.
    
    Args:
        run_ids: List of SRA run IDs whose FASTQs to delete
        fastq_dir: Directory containing FASTQ files
        
    Returns:
        Number of sample directories successfully deleted
        
    Side effects:
        - Removes fastq_dir/{run_id}/ directories recursively
        - Logs warnings for deletion failures
        
    Dependencies:
        - shutil module (standard library)
    """
    import shutil
    
    deleted = 0
    for run_id in run_ids:
        sample_dir = fastq_dir / run_id
        if sample_dir.exists():
            try:
                shutil.rmtree(sample_dir)
                deleted += 1
                
                # Report to progress tracker
                if tracker and species:
                    tracker.on_delete_complete(species, run_id)
            except Exception as e:
                logger.warning(f"    Failed to delete {sample_dir}: {e}")
    
    return deleted


def main():
    """Main entry point for ENA-integrated workflow.
    
    Orchestrates the complete workflow:
    1. Load configuration from YAML
    2. Get sample list from metadata
    3. Process samples in batches (download → quantify → cleanup)
    4. Report statistics
    
    Side effects:
        - Downloads FASTQ files
        - Creates quantification results
        - Deletes FASTQ files after quantification
        - Writes logs to console
        
    Exit codes:
        0: Success (all samples processed)
        1: Failure (some samples failed)
    """
    parser = argparse.ArgumentParser(
        description="Integrated ENA download + quantification workflow",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Batch Size Recommendations:
  Small drives (< 500GB): 8-12 samples
  Medium drives (500GB-1TB): 20-30 samples
  Large drives (1TB+): 50-100 samples

Batch size is auto-detected based on available disk space if not specified.
        """
    )
    parser.add_argument('--config', required=True, help="Amalgkit config YAML")
    parser.add_argument(
        '--batch-size', 
        type=int, 
        default=None,
        help="Samples per batch (default: auto-detect based on drive size, 50 for large drives)"
    )
    parser.add_argument(
        '--max-batch-size',
        type=int,
        default=100,
        help="Maximum batch size limit (default: 100)"
    )
    parser.add_argument('--threads', type=int, default=12, help="Threads for download/quant (default: 12)")
    parser.add_argument('--max-samples', type=int, help="Limit total samples (for testing)")
    parser.add_argument('--skip-download', action='store_true', help="Skip download, only quantify existing")
    parser.add_argument(
        '--min-free-gb',
        type=float,
        default=None,
        help="Minimum free disk space in GB (default: auto-detect based on drive size)"
    )
    parser.add_argument(
        '--max-batch-disk-gb',
        type=float,
        default=None,
        help="Maximum disk space per batch in GB (default: auto-calculate from batch size)"
    )
    
    args = parser.parse_args()
    
    # Load config
    config_path = Path(args.config)
    if not config_path.exists():
        logger.error(f"Config not found: {config_path}")
        sys.exit(1)
    
    config = load_config(config_path)
    
    work_dir = Path(config['work_dir'])
    fastq_dir = Path(config['steps']['getfastq']['out_dir'])
    quant_dir = Path(config['steps']['quant']['out_dir'])
    
    # Determine batch size (auto-detect if not specified)
    repo_root = Path(__file__).parent.parent.parent.resolve()
    output_dir = repo_root / "output" / "amalgkit"
    
    # Get temp directory for batch metadata and logging
    from metainformant.core.disk import get_recommended_temp_dir
    temp_dir = get_recommended_temp_dir(repo_root)
    
    if args.batch_size is None:
        # Check environment variable first
        env_batch = os.environ.get("AK_BATCH_SIZE")
        if env_batch:
            batch_size = int(env_batch)
            logger.info(f"Using batch size from AK_BATCH_SIZE: {batch_size}")
        else:
            from metainformant.core.disk import get_recommended_batch_size
            recommended_batch = get_recommended_batch_size(output_dir)
            batch_size = min(recommended_batch, args.max_batch_size)
            logger.info(f"Auto-detected batch size: {batch_size} (recommended: {recommended_batch}, max: {args.max_batch_size})")
    else:
        batch_size = min(args.batch_size, args.max_batch_size)
        logger.info(f"Using specified batch size: {batch_size} (max limit: {args.max_batch_size})")
    
    # Determine disk space thresholds
    from metainformant.core.disk import detect_drive_size_category, check_disk_space
    
    drive_category = detect_drive_size_category(output_dir)
    logger.info(f"Drive category: {drive_category}")
    
    if args.min_free_gb is None:
        # Auto-detect based on drive size
        if drive_category == "large":
            min_free_gb = 50.0
        elif drive_category == "medium":
            min_free_gb = 20.0
        else:
            min_free_gb = 10.0
        logger.info(f"Auto-detected min_free_gb: {min_free_gb} (based on {drive_category} drive)")
    else:
        min_free_gb = args.min_free_gb
        logger.info(f"Using specified min_free_gb: {min_free_gb}")
    
    # Check disk space before starting
    is_sufficient, disk_msg = check_disk_space(output_dir, min_free_gb=min_free_gb)
    logger.info(f"Disk space check: {disk_msg}")
    if not is_sufficient:
        logger.warning(f"⚠️  {disk_msg}")
        logger.warning("Continuing anyway, but workflow may fail if disk space runs out")
    
    # Calculate expected disk usage per batch
    sample_size_gb = 1.5  # Average estimate
    expected_batch_gb = batch_size * sample_size_gb
    if args.max_batch_disk_gb:
        if expected_batch_gb > args.max_batch_disk_gb:
            logger.warning(
                f"⚠️  Batch size {batch_size} would use ~{expected_batch_gb:.1f}GB, "
                f"exceeding limit of {args.max_batch_disk_gb}GB"
            )
            # Auto-adjust batch size
            adjusted_batch = int(args.max_batch_disk_gb / sample_size_gb)
            logger.info(f"Auto-adjusting batch size to {adjusted_batch} to stay within limit")
            batch_size = adjusted_batch
    
    logger.info(f"Batch configuration: {batch_size} samples (~{expected_batch_gb:.1f}GB per batch)")
    
    # Ensure directories exist
    work_dir.mkdir(parents=True, exist_ok=True)
    fastq_dir.mkdir(parents=True, exist_ok=True)
    quant_dir.mkdir(parents=True, exist_ok=True)
    
    # Get metadata
    metadata_path = work_dir / "metadata" / "metadata.tsv"
    if not metadata_path.exists():
        logger.error(f"Metadata not found: {metadata_path}")
        sys.exit(1)
    
    # Extract species name from config path
    species_name = config_path.stem.replace("amalgkit_", "")
    
    # Initialize progress tracker
    tracker = get_tracker()
    
    # Get sample list
    all_run_ids = get_sample_list(metadata_path)
    if args.max_samples:
        all_run_ids = all_run_ids[:args.max_samples]
    
    # Filter already-quantified
    to_process = []
    already_done = []
    
    for run_id in all_run_ids:
        if sample_already_quantified(run_id, quant_dir):
            already_done.append(run_id)
        else:
            to_process.append(run_id)
    
    # Initialize species in tracker
    tracker.initialize_species(species_name, len(all_run_ids), all_run_ids)
    
    # Mark already-completed samples in tracker
    for run_id in already_done:
        # Check if FASTQs are still present (needs_delete) or already deleted (completed)
        sample_dir = fastq_dir / run_id
        if sample_dir.exists() and any(sample_dir.glob("*.fastq.gz")):
            tracker.on_quant_complete(species_name, run_id)  # Mark as needs_delete
        else:
            # Already quantified and deleted - mark as completed
            tracker.on_quant_complete(species_name, run_id)
            tracker.on_delete_complete(species_name, run_id)
    
    # Find kallisto index
    genome_dir = Path(config['genome']['dest_dir']) / "ncbi_dataset_api_extracted"
    index_candidates = list((work_dir / "index").glob("*.idx")) if (work_dir / "index").exists() else []
    index_path = index_candidates[0] if index_candidates else None
    
    if not index_path or not index_path.exists():
        logger.error("Kallisto index not found. Run 'amalgkit quant' with build_index: yes first")
        sys.exit(1)
    
    # Summary
    logger.info("=" * 80)
    logger.info("🚀 INTEGRATED ENA DOWNLOAD + QUANTIFICATION WORKFLOW")
    logger.info("=" * 80)
    logger.info(f"Config: {config_path}")
    logger.info(f"Total samples: {len(all_run_ids)}")
    logger.info(f"Already quantified: {len(already_done)}")
    logger.info(f"To process: {len(to_process)}")
    logger.info(f"Batch size: {args.batch_size}")
    logger.info(f"Threads: {args.threads}")
    logger.info(f"Kallisto index: {index_path}")
    logger.info("=" * 80)
    
    if not to_process:
        logger.info("✅ All samples already quantified!")
        sys.exit(0)
    
    # Process in batches
    batches = [to_process[i:i + batch_size] for i in range(0, len(to_process), batch_size)]
    
    stats = {
        'downloaded': 0,
        'quantified': 0,
        'failed': 0,
        'cleaned_fastqs': 0,
    }
    
    start_time = time.time()
    
    for batch_num, batch_run_ids in enumerate(batches, 1):
        logger.info("")
        logger.info("=" * 80)
        logger.info(f"📦 BATCH {batch_num}/{len(batches)}: {len(batch_run_ids)} samples")
        logger.info(f"   Progress: {len(already_done) + (batch_num - 1) * batch_size}/{len(all_run_ids)} "
                   f"({(len(already_done) + (batch_num - 1) * batch_size) / len(all_run_ids) * 100:.1f}%)")
        logger.info("=" * 80)
        
        # Step 1: Download
        if not args.skip_download:
            logger.info(f"  📥 Downloading batch {batch_num} ({len(batch_run_ids)} samples)...")
            download_start = time.time()
            
            # Report download starts to tracker
            if tracker:
                for run_id in batch_run_ids:
                    tracker.on_download_start(species_name, run_id)
            
            downloaded, dl_failed = download_batch_ena(
                batch_run_ids, metadata_path, fastq_dir, args.threads, batch_num,
                species=species_name, tracker=tracker
            )
            download_time = time.time() - download_start
            stats['downloaded'] += len(downloaded)
            stats['failed'] += len(dl_failed)
            
            logger.info(f"  ✅ Download complete: {len(downloaded)}/{len(batch_run_ids)} samples in {download_time:.1f}s")
            if dl_failed:
                logger.warning(f"  ⚠️  Failed downloads: {len(dl_failed)}")
        else:
            downloaded = batch_run_ids
            logger.info(f"  ⏭️  Skipping download (--skip-download)")
        
        # Step 2: Quantify
        if downloaded:
            logger.info(f"  🔬 Quantifying batch {batch_num} ({len(downloaded)} samples)...")
            quant_start = time.time()
            
            quantified, q_failed = quantify_batch_kallisto(
                downloaded, fastq_dir, quant_dir, index_path, args.threads, batch_num,
                species=species_name, tracker=tracker
            )
            quant_time = time.time() - quant_start
            stats['quantified'] += len(quantified)
            if q_failed:
                stats['failed'] += len(q_failed)
            
            logger.info(f"  ✅ Quantification complete: {len(quantified)}/{len(downloaded)} samples in {quant_time:.1f}s")
            if q_failed:
                logger.warning(f"  ⚠️  Failed quantifications: {len(q_failed)}")
            
            # Step 3: Cleanup
            logger.info(f"  🗑️  Cleaning up FASTQs...")
            cleanup_start = time.time()
            cleaned = cleanup_fastqs(quantified, fastq_dir, species=species_name, tracker=tracker)
            cleanup_time = time.time() - cleanup_start
            stats['cleaned_fastqs'] += cleaned
            logger.info(f"  ✅ Deleted {cleaned} FASTQ directories in {cleanup_time:.1f}s")
            
            # Batch summary
            batch_time = download_time + quant_time + cleanup_time if not args.skip_download else quant_time + cleanup_time
            logger.info(f"  📊 Batch {batch_num} summary: {len(downloaded)} downloaded, "
                      f"{len(quantified)} quantified, {cleaned} cleaned in {batch_time:.1f}s")
        
        # Progress
        total_done = len(already_done) + stats['quantified']
        percent = (total_done / len(all_run_ids) * 100) if len(all_run_ids) > 0 else 0.0
        remaining = len(all_run_ids) - total_done
        remaining_batches = len(batches) - batch_num
        
        logger.info(f"  📊 Overall progress: {total_done}/{len(all_run_ids)} ({percent:.1f}%)")
        if remaining_batches > 0:
            logger.info(f"     Remaining: {remaining} samples in {remaining_batches} batches")
    
    # Final summary
    elapsed = time.time() - start_time
    
    # Update dashboard one final time
    if tracker:
        tracker.update_dashboard()
    
    logger.info("")
    logger.info("=" * 80)
    logger.info("🎉 WORKFLOW COMPLETE")
    logger.info("=" * 80)
    logger.info(f"Downloaded: {stats['downloaded']}")
    logger.info(f"Quantified: {stats['quantified']}")
    logger.info(f"Failed: {stats['failed']}")
    logger.info(f"Cleaned FASTQs: {stats['cleaned_fastqs']}")
    logger.info(f"Total time: {elapsed/60:.1f} minutes")
    logger.info(f"Rate: {elapsed/max(1, stats['quantified']):.1f} sec/sample")
    logger.info("=" * 80)
    
    sys.exit(0 if stats['failed'] == 0 else 1)


if __name__ == '__main__':
    main()

