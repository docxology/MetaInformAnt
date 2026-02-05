#!/usr/bin/env python3
"""Parallel processor for Apis mellifera samples.

This script processes multiple samples concurrently using a worker pool:
1. Download SRA file using prefetch
2. Extract FASTQ using fasterq-dump  
3. Quantify using kallisto
4. Clean up FASTQ files
5. Move to next sample

Uses a configurable number of parallel workers (default: 3) to significantly
speed up processing compared to the sequential version.
"""

import subprocess
import csv
import sys
from pathlib import Path
import shutil
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
import multiprocessing
import logging
from datetime import datetime, timedelta

# Configuration
WORK_DIR = Path("output/amalgkit/apis_mellifera_all/work")
FASTQ_DIR = Path("/Volumes/blue/data/apis_mellifera")  # External drive for downloads
QUANT_DIR = WORK_DIR / "quant"
INDEX_FILE = WORK_DIR / "index/Apis_mellifera_transcripts.idx"
METADATA_FILE = WORK_DIR / "metadata/metadata_selected.tsv"
PROGRESS_FILE = WORK_DIR / "parallel_progress.json"
LOG_FILE = WORK_DIR / "parallel_processing.log"

# Parallel processing settings
NUM_WORKERS = 10  # Number of parallel workers (maximized - external drive has 1.8TB)
THREADS_PER_SAMPLE = 2  # Threads for fasterq-dump and kallisto per sample
DISK_CHECK_INTERVAL = 10  # Check disk every N samples
MIN_FREE_GB = 50  # Minimum free disk space in GB
MIN_SRA_BYTES = 500 * 1024  # 500KB min - tiny samples fail kallisto
MAX_SRA_BYTES = 3 * 1024**3  # 3GB max - start with smaller samples to understand scaling

# Skip problematic sample prefixes that have abnormally slow extraction
# SRR340* = PRJNA1277852, SRR260* = PRJNA1017469
SKIP_PREFIXES = ("SRR340", "SRR260", "SRR345", "SRR351", "SRR353")

# Ensure directories exist
FASTQ_DIR.mkdir(parents=True, exist_ok=True)
QUANT_DIR.mkdir(parents=True, exist_ok=True)

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s [%(levelname)s] %(message)s',
    handlers=[
        logging.FileHandler(LOG_FILE),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)


def get_free_disk_gb() -> float:
    """Get free disk space in GB on the FASTQ download directory."""
    import os
    stat = os.statvfs(FASTQ_DIR)
    return (stat.f_frsize * stat.f_bavail) / (1024**3)


def get_processed_samples() -> set:
    """Get set of already quantified samples."""
    processed = set()
    if not QUANT_DIR.exists():
        return processed
    for sample_dir in QUANT_DIR.iterdir():
        if sample_dir.is_dir():
            abundance_file = sample_dir / "abundance.tsv"
            if abundance_file.exists():
                processed.add(sample_dir.name)
    return processed


def download_sample(sample_id: str, sample_dir: Path) -> tuple[bool, str]:
    """Download SRA file using prefetch."""
    sra_file = sample_dir / f"{sample_id}.sra"
    if sra_file.exists():
        return True, "SRA exists"
    
    cmd = [
        "prefetch",
        "--force", "no",
        "--max-size", "100G",
        "-O", str(FASTQ_DIR),
        sample_id
    ]
    
    result = subprocess.run(
        cmd, 
        capture_output=True, 
        text=True,
        timeout=7200  # 2 hour timeout for large files
    )
    
    if result.returncode != 0:
        return False, f"prefetch failed: {result.stderr[:200]}"
    
    return sra_file.exists(), "Downloaded"


def extract_fastq(sample_id: str, sample_dir: Path) -> tuple[bool, str]:
    """Extract FASTQ from SRA file."""
    sra_file = sample_dir / f"{sample_id}.sra"
    
    # Check if already extracted
    fastq_files = list(sample_dir.glob("*.fastq.gz"))
    if fastq_files:
        return True, f"FASTQ exists ({len(fastq_files)} files)"
    
    if not sra_file.exists():
        return False, "SRA file not found"
    
    cmd = [
        "fasterq-dump",
        "--outdir", str(sample_dir),
        "--temp", str(sample_dir),
        "--threads", str(THREADS_PER_SAMPLE),
        "--split-3",
        str(sra_file)
    ]
    
    result = subprocess.run(
        cmd,
        capture_output=True,
        text=True,
        timeout=7200  # 2 hour timeout for large files
    )
    
    if result.returncode != 0:
        # Check for disk limit errors
        if "disk-limit" in result.stderr.lower():
            return False, "disk-limit exceeded"
        return False, f"fasterq-dump failed: {result.stderr[:200]}"
    
    # Compress with gzip
    for fq in sample_dir.glob("*.fastq"):
        subprocess.run(
            ["gzip", "-f", str(fq)],
            capture_output=True,
            timeout=600
        )
    
    # Delete the .sra file to save space
    if sra_file.exists():
        sra_file.unlink()
    
    return True, "Extracted"


def quantify_sample(sample_id: str, sample_dir: Path, index_file: Path) -> tuple[bool, str]:
    """Run kallisto quantification."""
    quant_output = QUANT_DIR / sample_id
    abundance_file = quant_output / "abundance.tsv"
    
    if abundance_file.exists():
        return True, "Already quantified"
    
    quant_output.mkdir(parents=True, exist_ok=True)
    
    fastq_files = sorted(sample_dir.glob("*.fastq.gz"))
    if not fastq_files:
        return False, "No FASTQ files found"
    
    cmd = [
        "kallisto", "quant",
        "-i", str(index_file),
        "-o", str(quant_output),
        "-t", str(THREADS_PER_SAMPLE),
    ]
    
    # Paired or single end
    if len(fastq_files) >= 2:
        cmd.extend([str(fastq_files[0]), str(fastq_files[1])])
    else:
        cmd.extend(["--single", "-l", "200", "-s", "30", str(fastq_files[0])])
    
    result = subprocess.run(
        cmd,
        capture_output=True,
        text=True,
        timeout=1800  # 30 min timeout
    )
    
    if result.returncode != 0:
        return False, f"kallisto failed: {result.stderr[:200]}"
    
    return abundance_file.exists(), "Quantified"


def cleanup_sample(sample_dir: Path) -> None:
    """Delete FASTQ files after quantification."""
    if sample_dir.exists():
        shutil.rmtree(sample_dir, ignore_errors=True)


def process_single_sample(args: tuple) -> dict:
    """Process a single sample end-to-end. Returns result dict with per-stage timing."""
    sample_id, index_file = args
    sample_dir = FASTQ_DIR / sample_id
    sample_dir.mkdir(parents=True, exist_ok=True)
    
    result = {
        "sample_id": sample_id,
        "success": False,
        "steps": {},
        "timing": {},  # Per-stage timing in seconds
        "error": None,
        "duration": 0
    }
    
    start_time = time.time()
    
    try:
        # Step 1: Download
        t0 = time.time()
        ok, msg = download_sample(sample_id, sample_dir)
        result["timing"]["download"] = time.time() - t0
        result["steps"]["download"] = msg
        if not ok:
            result["error"] = f"Download: {msg}"
            cleanup_sample(sample_dir)
            return result
        
        # Step 2: Extract
        t0 = time.time()
        ok, msg = extract_fastq(sample_id, sample_dir)
        result["timing"]["extract"] = time.time() - t0
        result["steps"]["extract"] = msg
        if not ok:
            result["error"] = f"Extract: {msg}"
            cleanup_sample(sample_dir)
            return result
        
        # Step 3: Quantify
        t0 = time.time()
        ok, msg = quantify_sample(sample_id, sample_dir, index_file)
        result["timing"]["quantify"] = time.time() - t0
        result["steps"]["quantify"] = msg
        if not ok:
            result["error"] = f"Quantify: {msg}"
            cleanup_sample(sample_dir)
            return result
        
        # Step 4: Cleanup
        cleanup_sample(sample_dir)
        result["steps"]["cleanup"] = "Cleaned"
        result["success"] = True
        
    except subprocess.TimeoutExpired:
        result["error"] = "Timeout"
        cleanup_sample(sample_dir)
    except Exception as e:
        result["error"] = str(e)
        cleanup_sample(sample_dir)
    
    result["duration"] = time.time() - start_time
    return result


def save_failed_samples(failed: list[dict], path: Path) -> None:
    """Save failed samples to a file for retry later."""
    import json
    with open(path, "w") as f:
        json.dump(failed, f, indent=2)


def main():
    """Main parallel processing loop."""
    logger.info("="*60)
    logger.info("Apis mellifera Parallel Processor")
    logger.info("="*60)
    
    # Verify prerequisites
    if not INDEX_FILE.exists():
        logger.error(f"Kallisto index not found: {INDEX_FILE}")
        sys.exit(1)
    
    if not METADATA_FILE.exists():
        logger.error(f"Metadata file not found: {METADATA_FILE}")
        sys.exit(1)
    
    # Check disk space
    free_gb = get_free_disk_gb()
    logger.info(f"Free disk space: {free_gb:.1f} GB")
    if free_gb < MIN_FREE_GB:
        logger.error(f"Insufficient disk space. Need at least {MIN_FREE_GB} GB.")
        sys.exit(1)
    
    # Load samples
    with open(METADATA_FILE, "r") as f:
        reader = csv.DictReader(f, delimiter="\t")
        
        # Filter by size
        all_samples = []  # List of (sample_id, size_bytes) tuples
        skipped_large = 0
        
        for row in reader:
            sample_id = row.get("run")
            if not sample_id:
                continue
            
            # Skip problematic prefixes
            if sample_id.startswith(SKIP_PREFIXES):
                continue
                
            # Check size
            try:
                size_str = row.get("size", "0")
                size_bytes = float(size_str) if size_str else 0
                
                # Skip too small (fail kallisto)
                if size_bytes < MIN_SRA_BYTES:
                    continue
                
                if size_bytes > MAX_SRA_BYTES:
                    skipped_large += 1
                    logger.warning(f"Skipping {sample_id}: size {size_bytes/1024**3:.1f}GB > {MAX_SRA_BYTES/1024**3:.1f}GB limit")
                    continue
                    
                all_samples.append((sample_id, size_bytes))
            except (ValueError, TypeError):
                # If size can't be parsed, add with size 0
                all_samples.append((sample_id, 0))
    
    # Sort by size (smallest first) for timing analysis
    all_samples.sort(key=lambda x: x[1])
    logger.info(f"Loaded {len(all_samples)} samples, sorted by size (smallest first)")
    logger.info(f"Skipped {skipped_large} large samples > {MAX_SRA_BYTES/1024**3:.1f}GB")
    
    # Get already processed
    processed = get_processed_samples()
    logger.info(f"Already processed: {len(processed)} samples")
    
    # Filter to unprocessed (all_samples is list of (sample_id, size_bytes) tuples)
    to_process = [(s, sz) for s, sz in all_samples if s not in processed]
    logger.info(f"Remaining to process: {len(to_process)} samples")
    if to_process:
        logger.info(f"Size range: {to_process[0][1]/1024**2:.1f}MB - {to_process[-1][1]/1024**2:.1f}MB")
    
    if not to_process:
        logger.info("All samples already processed!")
        return
    
    # Stats tracking
    success_count = 0
    fail_count = 0
    failed_samples = []
    start_time = time.time()
    
    logger.info(f"Starting parallel processing with {NUM_WORKERS} workers")
    logger.info(f"Threads per sample: {THREADS_PER_SAMPLE}")
    
    # Create work items: (sample_id, index_file) - size kept for logging
    sample_sizes = {s: sz for s, sz in to_process}  # Map for size lookup
    work_items = [(s, INDEX_FILE) for s, sz in to_process]
    
    try:
        with ProcessPoolExecutor(max_workers=NUM_WORKERS) as executor:
            futures = {
                executor.submit(process_single_sample, item): item[0]
                for item in work_items
            }
            
            for i, future in enumerate(as_completed(futures), 1):
                sample_id = futures[future]
                
                try:
                    result = future.result(timeout=7200)  # 2 hour max per sample
                    
                    if result["success"]:
                        success_count += 1
                        size_mb = sample_sizes.get(sample_id, 0) / 1024**2
                        t = result.get("timing", {})
                        dl_t = t.get("download", 0)
                        ex_t = t.get("extract", 0)
                        qt_t = t.get("quantify", 0)
                        logger.info(
                            f"[{i}/{len(to_process)}] ✓ {sample_id} "
                            f"({result['duration']:.1f}s, {size_mb:.1f}MB) "
                            f"[dl:{dl_t:.0f}s ex:{ex_t:.0f}s qt:{qt_t:.0f}s]"
                        )
                        # Flush log handlers for real-time visibility
                        for handler in logger.handlers:
                            handler.flush()
                    else:
                        fail_count += 1
                        failed_samples.append(result)
                        logger.warning(
                            f"[{i}/{len(to_process)}] ✗ {sample_id}: "
                            f"{result['error']}"
                        )
                        # Flush log handlers for real-time visibility
                        for handler in logger.handlers:
                            handler.flush()
                    
                    # Progress update every 10 samples for better visibility
                    if i % 10 == 0:
                        elapsed = time.time() - start_time
                        rate = success_count / (elapsed / 3600) if elapsed > 0 else 0
                        remaining = len(to_process) - i
                        eta = remaining / rate if rate > 0 else 0
                        
                        logger.info(
                            f"Progress: {i}/{len(to_process)} "
                            f"({success_count} ok, {fail_count} failed) "
                            f"Rate: {rate:.1f}/hr, ETA: {eta/24:.1f} days"
                        )
                        
                        # Check disk space
                        free_gb = get_free_disk_gb()
                        if free_gb < MIN_FREE_GB:
                            logger.error(
                                f"Low disk space: {free_gb:.1f} GB. Stopping."
                            )
                            break
                
                except Exception as e:
                    fail_count += 1
                    logger.error(f"[{i}/{len(to_process)}] {sample_id}: {e}")
                    failed_samples.append({
                        "sample_id": sample_id,
                        "error": str(e)
                    })
    
    except KeyboardInterrupt:
        logger.info("Interrupted by user")
    
    # Final summary
    elapsed = time.time() - start_time
    elapsed_str = str(timedelta(seconds=int(elapsed)))
    
    logger.info("="*60)
    logger.info("SUMMARY")
    logger.info("="*60)
    logger.info(f"Elapsed time: {elapsed_str}")
    logger.info(f"Processed: {success_count + fail_count}")
    logger.info(f"Successful: {success_count}")
    logger.info(f"Failed: {fail_count}")
    logger.info(f"Total quantified: {len(get_processed_samples())}/{len(all_samples)}")
    
    if success_count > 0:
        rate = success_count / (elapsed / 3600)
        logger.info(f"Processing rate: {rate:.1f} samples/hour")
    
    # Save failed samples for retry
    if failed_samples:
        failed_file = WORK_DIR / "failed_samples.json"
        save_failed_samples(failed_samples, failed_file)
        logger.info(f"Failed samples saved to: {failed_file}")


if __name__ == "__main__":
    main()
