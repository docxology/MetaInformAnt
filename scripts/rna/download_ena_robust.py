#!/usr/bin/env python3
"""
Robust ENA FASTQ downloader with retry logic and resume capability.

Downloads FASTQ files directly from ENA (European Nucleotide Archive) which is:
- Much faster than SRA toolkit (100-500x)
- More reliable (direct FASTQ, no conversion)
- Has multiple mirror locations
- Supports resume for interrupted downloads

Usage:
    python3 download_ena_robust.py --metadata metadata.tsv --out-dir fastq/ --threads 12
"""

import argparse
import concurrent.futures
import logging
import subprocess
import sys
import time
from pathlib import Path
from typing import Optional

# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


def read_metadata(metadata_file: Path) -> list[dict]:
    """Read metadata TSV and extract run IDs."""
    import csv
    
    with open(metadata_file) as f:
        reader = csv.DictReader(f, delimiter='\t')
        rows = list(reader)
    
    if not rows or 'run' not in rows[0]:
        raise ValueError(f"Metadata file must have 'run' column: {metadata_file}")
    
    return rows


def get_ena_fastq_urls(run_id: str) -> list[str]:
    """
    Get ENA FASTQ URLs for a run ID using the ENA API.
    
    This is more reliable than guessing the directory structure.
    """
    import urllib.request
    import urllib.error
    
    try:
        # Query ENA API for FASTQ FTP URLs
        api_url = f"https://www.ebi.ac.uk/ena/portal/api/filereport?accession={run_id}&result=read_run&fields=fastq_ftp"
        
        with urllib.request.urlopen(api_url, timeout=30) as response:
            content = response.read().decode('utf-8')
        
        lines = content.strip().split('\n')
        if len(lines) < 2:
            logger.warning(f"  No FASTQ URLs found for {run_id}")
            return []
        
        # Parse the FTP URLs from the second line
        ftp_urls = lines[1].split('\t')[1] if len(lines[1].split('\t')) > 1 else ""
        
        if not ftp_urls:
            logger.warning(f"  No FASTQ files available for {run_id}")
            return []
        
        # Split multiple URLs (semicolon-separated for paired-end)
        urls = [f"ftp://{url}" for url in ftp_urls.split(';') if url]
        
        return urls
    
    except urllib.error.URLError as e:
        logger.warning(f"  Failed to query ENA API for {run_id}: {e}")
        return []
    except Exception as e:
        logger.warning(f"  Error getting URLs for {run_id}: {e}")
        return []


def download_file_robust(url: str, output_path: Path, max_retries: int = 5) -> bool:
    """
    Download a file with retry logic and resume capability.
    
    Uses wget with:
    - Retry logic (up to max_retries attempts)
    - Resume capability (-c flag)
    - Timeout handling
    - Connection retry
    """
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    for attempt in range(1, max_retries + 1):
        try:
            # Use wget for robust downloading with resume
            cmd = [
                'wget',
                '--continue',  # Resume partial downloads
                '--timeout=60',  # 60 second timeout
                '--tries=3',  # 3 tries per wget invocation
                '--no-verbose',
                '-O', str(output_path),
                url
            ]
            
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=3600)
            
            if result.returncode == 0:
                # Verify file exists and has content
                if output_path.exists() and output_path.stat().st_size > 0:
                    return True
                else:
                    logger.warning(f"  Downloaded file is empty: {output_path}")
            else:
                if "404" in result.stderr or "not found" in result.stderr.lower():
                    # File doesn't exist on server, don't retry
                    return False
                
                logger.warning(f"  Attempt {attempt}/{max_retries} failed: {result.stderr[:200]}")
        
        except subprocess.TimeoutExpired:
            logger.warning(f"  Attempt {attempt}/{max_retries} timed out after 1 hour")
        except Exception as e:
            logger.warning(f"  Attempt {attempt}/{max_retries} error: {e}")
        
        if attempt < max_retries:
            wait_time = min(60, 5 * attempt)  # Exponential backoff, max 60s
            logger.info(f"  Waiting {wait_time}s before retry...")
            time.sleep(wait_time)
    
    return False


def download_sample(run_id: str, out_dir: Path, max_retries: int = 5) -> tuple[str, bool, list[Path]]:
    """
    Download FASTQ files for a sample from ENA.
    
    Returns:
        (run_id, success, list_of_downloaded_files)
    """
    logger.info(f"📥 Downloading {run_id}")
    
    sample_dir = out_dir / run_id
    sample_dir.mkdir(parents=True, exist_ok=True)
    
    # Check if already downloaded
    existing_fastq = list(sample_dir.glob("*.fastq.gz"))
    if existing_fastq and all(f.stat().st_size > 1000000 for f in existing_fastq):  # >1MB
        logger.info(f"  ✅ Already downloaded: {run_id} ({len(existing_fastq)} files)")
        return (run_id, True, existing_fastq)
    
    # Get URLs from ENA API
    urls = get_ena_fastq_urls(run_id)
    
    if not urls:
        logger.error(f"  ❌ Failed: {run_id} (no URLs from ENA)")
        return (run_id, False, [])
    
    downloaded_files = []
    
    # Download each file
    for url in urls:
        filename = Path(url).name
        output_path = sample_dir / filename
        
        logger.info(f"  Downloading: {filename}")
        success = download_file_robust(url, output_path, max_retries=max_retries)
        
        if success:
            size_mb = output_path.stat().st_size / 1024 / 1024
            logger.info(f"  ✅ Downloaded: {filename} ({size_mb:.1f} MB)")
            downloaded_files.append(output_path)
        else:
            logger.warning(f"  ⚠️  Failed to download: {filename}")
    
    if downloaded_files:
        logger.info(f"  ✅ Completed: {run_id} ({len(downloaded_files)} files)")
        return (run_id, True, downloaded_files)
    else:
        logger.error(f"  ❌ Failed: {run_id} (no files downloaded)")
        return (run_id, False, [])


def main():
    parser = argparse.ArgumentParser(description="Robust ENA FASTQ downloader")
    parser.add_argument('--metadata', required=True, help="Metadata TSV with 'run' column")
    parser.add_argument('--out-dir', required=True, help="Output directory for FASTQ files")
    parser.add_argument('--threads', type=int, default=12, help="Number of parallel downloads (default: 12)")
    parser.add_argument('--max-retries', type=int, default=5, help="Max retries per file (default: 5)")
    parser.add_argument('--max-samples', type=int, help="Limit number of samples to download")
    
    args = parser.parse_args()
    
    metadata_file = Path(args.metadata)
    out_dir = Path(args.out_dir)
    
    if not metadata_file.exists():
        logger.error(f"Metadata file not found: {metadata_file}")
        sys.exit(1)
    
    out_dir.mkdir(parents=True, exist_ok=True)
    
    # Read metadata
    logger.info("="*80)
    logger.info("🚀 ENA ROBUST DOWNLOADER")
    logger.info("="*80)
    logger.info(f"Metadata: {metadata_file}")
    logger.info(f"Output: {out_dir}")
    logger.info(f"Threads: {args.threads}")
    logger.info(f"Max retries: {args.max_retries}")
    
    rows = read_metadata(metadata_file)
    run_ids = [row['run'] for row in rows if row.get('run')]
    
    if args.max_samples:
        run_ids = run_ids[:args.max_samples]
    
    logger.info(f"Samples to download: {len(run_ids)}")
    logger.info("="*80)
    
    # Download in parallel
    start_time = time.time()
    successful = []
    failed = []
    
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.threads) as executor:
        futures = {
            executor.submit(download_sample, run_id, out_dir, args.max_retries): run_id
            for run_id in run_ids
        }
        
        for future in concurrent.futures.as_completed(futures):
            run_id, success, files = future.result()
            if success:
                successful.append(run_id)
            else:
                failed.append(run_id)
    
    # Summary
    elapsed = time.time() - start_time
    logger.info("")
    logger.info("="*80)
    logger.info("📊 DOWNLOAD COMPLETE")
    logger.info("="*80)
    logger.info(f"Total samples: {len(run_ids)}")
    logger.info(f"Successful: {len(successful)} ({len(successful)*100//len(run_ids)}%)")
    logger.info(f"Failed: {len(failed)} ({len(failed)*100//len(run_ids) if len(run_ids) > 0 else 0}%)")
    logger.info(f"Time: {elapsed/60:.1f} minutes ({elapsed/len(run_ids):.1f} sec/sample)")
    logger.info("="*80)
    
    if failed:
        logger.warning(f"Failed samples: {', '.join(failed[:20])}")
        if len(failed) > 20:
            logger.warning(f"  ... and {len(failed)-20} more")
    
    sys.exit(0 if len(successful) > 0 else 1)


if __name__ == '__main__':
    main()

