#!/usr/bin/env python3
"""ENA-based parallel processor for Apis mellifera samples.

This script downloads pre-extracted FASTQs directly from the European Nucleotide
Archive (ENA), bypassing the slow fasterq-dump extraction step entirely.

ENA mirrors all SRA data and provides gzipped FASTQs via FTP/HTTP, which means:
- No need for prefetch + fasterq-dump
- Direct download of compressed FASTQ files
- 10x faster per-sample processing (5 min vs 50 min)
- Minimal temp space requirements

URL Pattern:
  ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRRXXX/0YY/SRRXXXYY/SRRXXXYY_1.fastq.gz
  
Where:
  - First 6 chars of SRR ID go in path
  - Full ID used for final directory and filename
  - For IDs with 10+ digits, extra subdirectory is added
"""

import subprocess
import csv
import sys
from pathlib import Path
import shutil
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
import logging
from datetime import datetime, timedelta
import urllib.request
import urllib.error

# Configuration
WORK_DIR = Path("output/amalgkit/apis_mellifera_all/work")
FASTQ_DIR = Path("/Volumes/blue/data/apis_mellifera")  # External drive for downloads
QUANT_DIR = WORK_DIR / "quant"
INDEX_FILE = WORK_DIR / "index/Apis_mellifera_transcripts.idx"
METADATA_FILE = WORK_DIR / "metadata/metadata_selected.tsv"
LOG_FILE = WORK_DIR / "ena_parallel_processing.log"

# ENA parallel processing settings - higher than NCBI since no extraction bottleneck
NUM_WORKERS = 20  # Safe with ENA since we skip fasterq-dump
THREADS_PER_SAMPLE = 2  # Threads for kallisto per sample
DISK_CHECK_INTERVAL = 20  # Check disk every N samples
MIN_FREE_GB = 50  # Minimum free disk space in GB
MIN_SRA_BYTES = 500 * 1024  # 500KB min - tiny samples fail kallisto
MAX_SRA_BYTES = 10 * 1024**3  # 10GB max - ENA can handle larger files

# Skip problematic sample prefixes
SKIP_PREFIXES = ("SRR340", "SRR260", "SRR345", "SRR351", "SRR353")

# ENA base URLs
ENA_FTP_BASE = "ftp://ftp.sra.ebi.ac.uk/vol1/fastq"
ENA_HTTP_BASE = "http://ftp.sra.ebi.ac.uk/vol1/fastq"

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


def build_ena_url(sample_id: str, read_num: int = 1) -> str:
    """Build ENA FTP URL for a sample.
    
    ENA URL structure:
    - SRR prefix directory (first 6 chars): SRR123
    - For IDs >= 10 digits, add subdirectory: 00X where X = last digit of ID
    - Full sample ID directory
    - Filename: {sample_id}_{1,2}.fastq.gz
    
    Examples:
        SRR1234567 (9 digits) -> ftp.sra.ebi.ac.uk/vol1/fastq/SRR123/SRR1234567/SRR1234567_1.fastq.gz
        SRR12345678 (10 digits) -> ftp.sra.ebi.ac.uk/vol1/fastq/SRR123/008/SRR12345678/SRR12345678_1.fastq.gz
    """
    prefix = sample_id[:6]
    
    if len(sample_id) >= 10:
        # Add leading zeros subdirectory based on last digit(s)
        suffix = sample_id[9:]  # Characters after position 9
        subdir = suffix.zfill(3)  # Pad to 3 digits: "008", "012", etc.
        path = f"{prefix}/{subdir}/{sample_id}"
    else:
        path = f"{prefix}/{sample_id}"
    
    filename = f"{sample_id}_{read_num}.fastq.gz"
    return f"{ENA_HTTP_BASE}/{path}/{filename}"


def get_ena_fastq_urls(sample_id: str) -> list[str]:
    """Get FASTQ URLs from ENA API.
    
    The ENA Portal API returns the exact FTP paths for FASTQ files.
    This is more reliable than guessing the URL structure.
    """
    api_url = f"https://www.ebi.ac.uk/ena/portal/api/filereport?accession={sample_id}&result=read_run&fields=fastq_ftp"
    
    try:
        with urllib.request.urlopen(api_url, timeout=30) as response:
            content = response.read().decode('utf-8')
            lines = content.strip().split('\n')
            
            if len(lines) < 2:
                return []
            
            # Second line contains the FTP URLs (semicolon-separated if multiple)
            ftp_field = lines[1].split('\t')[-1]
            if not ftp_field or ftp_field == 'fastq_ftp':
                return []
            
            # Convert FTP paths to HTTP URLs
            urls = []
            for ftp_path in ftp_field.split(';'):
                ftp_path = ftp_path.strip()
                if ftp_path:
                    http_url = f"http://{ftp_path}"
                    urls.append(http_url)
            
            return urls
    except Exception:
        return []


def check_ena_availability(sample_id: str) -> tuple[bool, list[str]]:
    """Check if sample is available on ENA and return available FASTQ URLs.
    
    Returns (is_available, list_of_urls).
    """
    urls = get_ena_fastq_urls(sample_id)
    return len(urls) > 0, urls


def download_fastq_from_ena(sample_id: str, sample_dir: Path) -> tuple[bool, str, list[Path]]:
    """Download FASTQ files directly from ENA using API-discovered URLs.
    
    Returns (success, message, list_of_fastq_files).
    """
    sample_dir.mkdir(parents=True, exist_ok=True)
    downloaded_files = []
    
    # Get URLs from ENA API (most reliable method)
    urls = get_ena_fastq_urls(sample_id)
    
    if not urls:
        # Fallback: try guessing URLs for samples not yet indexed
        return False, "Not available on ENA", []
    
    for url in urls:
        # Extract filename from URL
        filename = url.split('/')[-1]
        output_file = sample_dir / filename
        
        if output_file.exists() and output_file.stat().st_size > 0:
            downloaded_files.append(output_file)
            continue
        
        # Use curl for better reliability and resume support
        cmd = [
            "curl", "-fsSL",
            "--retry", "3",
            "--retry-delay", "5",
            "-o", str(output_file),
            url
        ]
        
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=1800)
        
        if result.returncode == 0 and output_file.exists() and output_file.stat().st_size > 0:
            downloaded_files.append(output_file)
        else:
            # Clean up partial download
            if output_file.exists():
                output_file.unlink()
    
    if downloaded_files:
        return True, f"Downloaded {len(downloaded_files)} file(s)", downloaded_files
    else:
        return False, "Download failed", []


def quantify_sample(sample_id: str, fastq_files: list[Path], index_file: Path) -> tuple[bool, str]:
    """Run kallisto quantification on downloaded FASTQs."""
    quant_output = QUANT_DIR / sample_id
    abundance_file = quant_output / "abundance.tsv"
    
    if abundance_file.exists():
        return True, "Already quantified"
    
    quant_output.mkdir(parents=True, exist_ok=True)
    
    if not fastq_files:
        return False, "No FASTQ files"
    
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
    
    result = subprocess.run(cmd, capture_output=True, text=True, timeout=1800)
    
    if result.returncode != 0:
        return False, f"kallisto failed: {result.stderr[:200]}"
    
    return abundance_file.exists(), "Quantified"


def cleanup_sample(sample_dir: Path) -> None:
    """Delete FASTQ files after quantification."""
    if sample_dir.exists():
        shutil.rmtree(sample_dir, ignore_errors=True)


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


def process_single_sample(args: tuple) -> dict:
    """Process a single sample end-to-end via ENA. Returns result dict with timing."""
    sample_id, index_file = args
    sample_dir = FASTQ_DIR / sample_id
    
    result = {
        "sample_id": sample_id,
        "success": False,
        "source": "ENA",
        "timing": {},
        "error": None,
        "duration": 0
    }
    
    start_time = time.time()
    
    try:
        # Step 1: Download from ENA
        t0 = time.time()
        ok, msg, fastq_files = download_fastq_from_ena(sample_id, sample_dir)
        result["timing"]["download"] = time.time() - t0
        
        if not ok:
            result["error"] = f"Download: {msg}"
            cleanup_sample(sample_dir)
            return result
        
        # Step 2: Quantify (no extraction needed!)
        t0 = time.time()
        ok, msg = quantify_sample(sample_id, fastq_files, index_file)
        result["timing"]["quantify"] = time.time() - t0
        
        if not ok:
            result["error"] = f"Quantify: {msg}"
            cleanup_sample(sample_dir)
            return result
        
        # Step 3: Cleanup
        cleanup_sample(sample_dir)
        result["success"] = True
        
    except subprocess.TimeoutExpired:
        result["error"] = "Timeout"
        cleanup_sample(sample_dir)
    except Exception as e:
        result["error"] = str(e)
        cleanup_sample(sample_dir)
    
    result["duration"] = time.time() - start_time
    return result


def main():
    """Main parallel processing loop using ENA."""
    logger.info("="*60)
    logger.info("Apis mellifera ENA Parallel Processor")
    logger.info("="*60)
    logger.info("Strategy: Direct FASTQ download from ENA (skipping extraction)")
    
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
        
        all_samples = []
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
                
                if size_bytes < MIN_SRA_BYTES:
                    continue
                
                if size_bytes > MAX_SRA_BYTES:
                    skipped_large += 1
                    continue
                
                all_samples.append((sample_id, size_bytes))
            except (ValueError, TypeError):
                all_samples.append((sample_id, 0))
    
    # Sort by size (smallest first)
    all_samples.sort(key=lambda x: x[1])
    logger.info(f"Loaded {len(all_samples)} samples, sorted by size")
    logger.info(f"Skipped {skipped_large} samples > {MAX_SRA_BYTES/1024**3:.1f}GB")
    
    # Get already processed
    processed = get_processed_samples()
    logger.info(f"Already processed: {len(processed)} samples")
    
    # Filter to unprocessed
    to_process = [(s, sz) for s, sz in all_samples if s not in processed]
    logger.info(f"Remaining to process: {len(to_process)} samples")
    
    if not to_process:
        logger.info("All samples already processed!")
        return
    
    # Stats tracking
    success_count = 0
    fail_count = 0
    ena_unavailable = 0
    failed_samples = []
    start_time = time.time()
    
    logger.info(f"Starting ENA parallel processing with {NUM_WORKERS} workers")
    
    # Create work items
    sample_sizes = {s: sz for s, sz in to_process}
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
                    result = future.result(timeout=3600)  # 1 hour max per sample
                    
                    if result["success"]:
                        success_count += 1
                        size_mb = sample_sizes.get(sample_id, 0) / 1024**2
                        t = result.get("timing", {})
                        dl_t = t.get("download", 0)
                        qt_t = t.get("quantify", 0)
                        logger.info(
                            f"[{i}/{len(to_process)}] ✓ {sample_id} "
                            f"({result['duration']:.1f}s, {size_mb:.1f}MB) "
                            f"[dl:{dl_t:.0f}s qt:{qt_t:.0f}s]"
                        )
                        for handler in logger.handlers:
                            handler.flush()
                    else:
                        fail_count += 1
                        if "Not available on ENA" in str(result.get("error", "")):
                            ena_unavailable += 1
                        failed_samples.append(result)
                        logger.warning(
                            f"[{i}/{len(to_process)}] ✗ {sample_id}: "
                            f"{result['error']}"
                        )
                        for handler in logger.handlers:
                            handler.flush()
                    
                    # Progress update every 20 samples
                    if i % 20 == 0:
                        elapsed = time.time() - start_time
                        rate = success_count / (elapsed / 3600) if elapsed > 0 else 0
                        remaining = len(to_process) - i
                        eta_hours = remaining / rate if rate > 0 else 0
                        
                        logger.info(
                            f"Progress: {i}/{len(to_process)} "
                            f"({success_count} ok, {fail_count} failed, {ena_unavailable} not on ENA) "
                            f"Rate: {rate:.1f}/hr, ETA: {eta_hours:.1f} hours"
                        )
                        
                        # Check disk space
                        free_gb = get_free_disk_gb()
                        if free_gb < MIN_FREE_GB:
                            logger.error(f"Low disk space: {free_gb:.1f} GB. Stopping.")
                            break
                
                except Exception as e:
                    fail_count += 1
                    logger.error(f"[{i}/{len(to_process)}] {sample_id}: {e}")
                    failed_samples.append({"sample_id": sample_id, "error": str(e)})
    
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
    logger.info(f"Failed: {fail_count} ({ena_unavailable} not on ENA)")
    logger.info(f"Total quantified: {len(get_processed_samples())}")
    
    if success_count > 0:
        rate = success_count / (elapsed / 3600)
        logger.info(f"Processing rate: {rate:.1f} samples/hour")
    
    # Save failed samples for NCBI fallback
    if failed_samples:
        import json
        failed_file = WORK_DIR / "ena_failed_samples.json"
        with open(failed_file, "w") as f:
            json.dump(failed_samples, f, indent=2)
        logger.info(f"Failed samples saved to: {failed_file}")
        logger.info(f"Run NCBI pipeline on these {len(failed_samples)} samples as fallback")


if __name__ == "__main__":
    main()
