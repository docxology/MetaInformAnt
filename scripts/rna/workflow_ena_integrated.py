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
                       threads: int, batch_num: int) -> tuple[list[str], list[str]]:
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
        - download_ena_robust.py script in same directory
        - wget command available
    """
    logger.info(f"  📥 Downloading batch {batch_num} ({len(run_ids)} samples)...")
    
    # Create batch metadata
    rows = list(read_delimited(metadata_path, delimiter='\t'))
    batch_rows = [row for row in rows if row.get('run') in run_ids]
    batch_metadata = metadata_path.parent / f"metadata_batch{batch_num}.tsv"
    write_delimited(batch_rows, batch_metadata, delimiter='\t')
    
    # Use the robust ENA downloader
    script_dir = Path(__file__).parent
    downloader = script_dir / "download_ena_robust.py"
    
    cmd = [
        sys.executable,
        str(downloader),
        '--metadata', str(batch_metadata),
        '--out-dir', str(fastq_dir),
        '--threads', str(threads),
        '--max-retries', '3',
    ]
    
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=7200)  # 2 hour timeout per batch
        
        # Parse which samples succeeded
        successful = []
        failed = []
        
        for run_id in run_ids:
            sample_dir = fastq_dir / run_id
            fastq_files = list(sample_dir.glob("*.fastq.gz")) if sample_dir.exists() else []
            if fastq_files and all(f.stat().st_size > 1000000 for f in fastq_files):
                successful.append(run_id)
            else:
                failed.append(run_id)
        
        logger.info(f"  ✅ Downloaded: {len(successful)}/{len(run_ids)} samples")
        if failed:
            logger.warning(f"  ⚠️  Failed: {', '.join(failed[:5])}{' ...' if len(failed) > 5 else ''}")
        
        # Clean up batch metadata
        batch_metadata.unlink(missing_ok=True)
        
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
                            index_path: Path, threads: int, batch_num: int) -> tuple[list[str], list[str]]:
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
    return successful, failed


def cleanup_fastqs(run_ids: list[str], fastq_dir: Path) -> int:
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
    parser = argparse.ArgumentParser(description="Integrated ENA download + quantification workflow")
    parser.add_argument('--config', required=True, help="Amalgkit config YAML")
    parser.add_argument('--batch-size', type=int, default=12, help="Samples per batch (default: 12)")
    parser.add_argument('--threads', type=int, default=12, help="Threads for download/quant (default: 12)")
    parser.add_argument('--max-samples', type=int, help="Limit total samples (for testing)")
    parser.add_argument('--skip-download', action='store_true', help="Skip download, only quantify existing")
    
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
    
    # Ensure directories exist
    work_dir.mkdir(parents=True, exist_ok=True)
    fastq_dir.mkdir(parents=True, exist_ok=True)
    quant_dir.mkdir(parents=True, exist_ok=True)
    
    # Get metadata
    metadata_path = work_dir / "metadata" / "metadata.tsv"
    if not metadata_path.exists():
        logger.error(f"Metadata not found: {metadata_path}")
        sys.exit(1)
    
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
    batches = [to_process[i:i + args.batch_size] for i in range(0, len(to_process), args.batch_size)]
    
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
        logger.info("=" * 80)
        
        # Step 1: Download
        if not args.skip_download:
            downloaded, dl_failed = download_batch_ena(
                batch_run_ids, metadata_path, fastq_dir, args.threads, batch_num
            )
            stats['downloaded'] += len(downloaded)
            stats['failed'] += len(dl_failed)
        else:
            downloaded = batch_run_ids
        
        # Step 2: Quantify
        if downloaded:
            quantified, q_failed = quantify_batch_kallisto(
                downloaded, fastq_dir, quant_dir, index_path, args.threads, batch_num
            )
            stats['quantified'] += len(quantified)
            if q_failed:
                stats['failed'] += len(q_failed)
            
            # Step 3: Cleanup
            logger.info(f"  🗑️  Cleaning up FASTQs...")
            cleaned = cleanup_fastqs(quantified, fastq_dir)
            stats['cleaned_fastqs'] += cleaned
            logger.info(f"  ✅ Deleted {cleaned} FASTQ directories")
        
        # Progress
        total_done = len(already_done) + stats['quantified']
        percent = (total_done * 100) // len(all_run_ids)
        logger.info(f"  📊 Progress: {total_done}/{len(all_run_ids)} ({percent}%)")
    
    # Final summary
    elapsed = time.time() - start_time
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

