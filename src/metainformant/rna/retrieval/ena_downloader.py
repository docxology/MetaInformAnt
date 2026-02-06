"""ENA (European Nucleotide Archive) FASTQ downloader with best practices.

This module provides functionality to download FASTQ files from the European
Nucleotide Archive (ENA) for RNA-seq samples. It handles URL resolution,
MD5 verification, resume-capable downloads, and implements best practices:

- Size-ordered downloads (smallest first for quick wins)
- Configurable size limits and timeouts
- Multi-source fallback (ENA FTP → ENA HTTP → NCBI SRA → AWS S3)

Example usage:
    >>> from metainformant.rna.retrieval.ena_downloader import download_sra_samples
    >>> download_sra_samples(["SRR12345"], Path("output/fastq"))
    
    # With best practices options:
    >>> download_sra_samples(
    ...     ["SRR12345", "SRR67890"], 
    ...     Path("output/fastq"),
    ...     max_size_bytes=5 * 1024**3,  # 5GB
    ...     timeout=1800,  # 30 min
    ...     sort_by_size=True,
    ...     use_fallback=True
    ... )
"""

from __future__ import annotations

import argparse
import hashlib
import os
import subprocess
import sys
import time
from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional, Tuple

import requests

from metainformant.core import logging

logger = logging.get_logger(__name__)

# Size constants
KB = 1024
MB = 1024 * KB
GB = 1024 * MB

# Default limits
DEFAULT_TIMEOUT = 600  # 10 minutes
DEFAULT_MAX_SIZE = 10 * GB  # 10 GB


@dataclass
class SampleInfo:
    """Information about an SRA sample including size and download URLs."""
    sra_id: str
    size_bytes: int
    ena_ftp_urls: List[Tuple[str, str]]  # (url, md5) tuples
    ena_http_urls: List[Tuple[str, str]]
    ncbi_urls: List[str]
    aws_urls: List[str]


def parse_size(size_str: str) -> int:
    """Parse human-readable size string to bytes.
    
    Args:
        size_str: Size like "500MB", "5GB", "1024KB"
        
    Returns:
        Size in bytes
    """
    size_str = size_str.strip().upper()
    if size_str.endswith("GB"):
        return int(float(size_str[:-2]) * GB)
    elif size_str.endswith("MB"):
        return int(float(size_str[:-2]) * MB)
    elif size_str.endswith("KB"):
        return int(float(size_str[:-2]) * KB)
    elif size_str.endswith("B"):
        return int(size_str[:-1])
    else:
        return int(size_str)


def format_size(size_bytes: int) -> str:
    """Format bytes to human-readable string."""
    if size_bytes >= GB:
        return f"{size_bytes / GB:.2f}GB"
    elif size_bytes >= MB:
        return f"{size_bytes / MB:.2f}MB"
    elif size_bytes >= KB:
        return f"{size_bytes / KB:.2f}KB"
    else:
        return f"{size_bytes}B"


def get_ena_sample_info(sra_id: str) -> Optional[SampleInfo]:
    """Query ENA to get complete sample information including sizes and URLs.

    Args:
        sra_id: SRA accession ID (e.g., SRR12345, ERR12345, DRR12345)

    Returns:
        SampleInfo object with all download options, or None on error.
    """
    # Request more fields including size
    fields = "fastq_ftp,fastq_md5,fastq_bytes,submitted_ftp,sra_ftp"
    url = (
        f"https://www.ebi.ac.uk/ena/portal/api/filereport?accession={sra_id}"
        f"&result=read_run&fields={fields}&format=tsv"
    )

    try:
        response = requests.get(url, timeout=15)
        response.raise_for_status()

        lines = response.text.strip().split("\n")
        if len(lines) < 2:
            logger.warning(f"[{sra_id}] No data found in ENA.")
            return None

        header = lines[0].split("\t")
        data = lines[1].split("\t")
        row = dict(zip(header, data))

        # Parse size
        size_bytes = 0
        if "fastq_bytes" in row and row["fastq_bytes"]:
            # Multiple files separated by ;
            sizes = row["fastq_bytes"].split(";")
            size_bytes = sum(int(s) for s in sizes if s.strip())

        # Parse FTP URLs
        ena_ftp_urls: List[Tuple[str, str]] = []
        ena_http_urls: List[Tuple[str, str]] = []
        
        if "fastq_ftp" in row and row["fastq_ftp"]:
            ftps = row["fastq_ftp"].split(";")
            md5s = row.get("fastq_md5", "").split(";") if row.get("fastq_md5") else [""] * len(ftps)
            
            for ftp, md5 in zip(ftps, md5s):
                if not ftp.startswith(("ftp://", "http://", "https://")):
                    ftp_url = f"ftp://{ftp}"
                    http_url = f"http://{ftp}"
                else:
                    ftp_url = ftp
                    http_url = ftp.replace("ftp://", "http://")
                ena_ftp_urls.append((ftp_url, md5))
                ena_http_urls.append((http_url, md5))

        # Build NCBI/AWS fallback URLs
        ncbi_urls = [f"https://sra-pub-run-odp.s3.amazonaws.com/sra/{sra_id}/{sra_id}"]
        aws_urls = [f"s3://sra-pub-run-odp/sra/{sra_id}/{sra_id}"]

        return SampleInfo(
            sra_id=sra_id,
            size_bytes=size_bytes,
            ena_ftp_urls=ena_ftp_urls,
            ena_http_urls=ena_http_urls,
            ncbi_urls=ncbi_urls,
            aws_urls=aws_urls,
        )

    except requests.exceptions.RequestException as e:
        logger.error(f"[{sra_id}] API Request failed: {e}")
        return None


def get_ena_links(sra_id: str) -> List[Tuple[str, str]]:
    """Query ENA to get FASTQ download links and MD5 checksums.
    
    Legacy function for backwards compatibility.

    Args:
        sra_id: SRA accession ID (e.g., SRR12345, ERR12345, DRR12345)

    Returns:
        List of (url, md5) tuples for each FASTQ file associated with the sample.
        Empty list if no data found or on error.
    """
    info = get_ena_sample_info(sra_id)
    if info and info.ena_ftp_urls:
        return info.ena_ftp_urls
    return []


def calculate_md5(file_path: Path) -> str:
    """Calculate MD5 checksum of a file efficiently.

    Args:
        file_path: Path to the file to checksum

    Returns:
        Hexadecimal MD5 digest string
    """
    hash_md5 = hashlib.md5()
    with open(file_path, "rb") as f:
        for chunk in iter(lambda: f.read(8192), b""):
            hash_md5.update(chunk)
    return hash_md5.hexdigest()


def clean_stagnant_file(file_path: Path) -> None:
    """Delete a file if it exists (cleanup for corrupted/incomplete downloads).

    Args:
        file_path: Path to the file to delete
    """
    if file_path.exists():
        try:
            os.remove(file_path)
            logger.info(f"Deleted incomplete/stagnant file: {file_path}")
        except OSError as e:
            logger.warning(f"Error deleting {file_path}: {e}")


def download_file(
    url: str, 
    dest_path: Path, 
    expected_md5: str = "", 
    retries: int = 3,
    timeout: int = DEFAULT_TIMEOUT
) -> bool:
    """Download file using curl with retries and optional MD5 verification.

    Args:
        url: URL to download from
        dest_path: Local path to save the file
        expected_md5: Expected MD5 checksum for verification (empty to skip)
        retries: Number of retry attempts
        timeout: Max time in seconds for the download

    Returns:
        True if download succeeded (and MD5 verified if provided), False otherwise
    """
    logger.info(f"Downloading {url} -> {dest_path}")

    for attempt in range(1, retries + 1):
        try:
            cmd = [
                "curl",
                "-L",
                "-C", "-",  # Resume
                "-o", str(dest_path),
                "--retry", "2",
                "--retry-delay", "5",
                "--connect-timeout", "60",
                "--max-time", str(timeout),
                "--keepalive-time", "60",
                "-f",  # Fail on HTTP errors
                url,
            ]

            subprocess.run(cmd, check=True, capture_output=True)

            if dest_path.exists():
                if expected_md5:
                    logger.info(f"Verifying MD5 for {dest_path.name}...")
                    current_md5 = calculate_md5(dest_path)
                    if current_md5 == expected_md5:
                        logger.info(f"[OK] {dest_path.name} verified.")
                        return True
                    else:
                        logger.warning(
                            f"[FAIL] MD5 mismatch (Attempt {attempt}/{retries}). "
                            f"Expected {expected_md5}, got {current_md5}"
                        )
                        clean_stagnant_file(dest_path)
                else:
                    # No MD5 to verify, assume success
                    logger.info(f"[OK] {dest_path.name} downloaded (no MD5 verification).")
                    return True

        except subprocess.CalledProcessError as e:
            logger.warning(f"[FAIL] Download failed (Attempt {attempt}/{retries}). Exit: {e.returncode}")
            if attempt < retries:
                time.sleep(5 * attempt)

    return False


def download_with_fallback(
    sample_info: SampleInfo,
    sample_dir: Path,
    timeout: int = DEFAULT_TIMEOUT
) -> bool:
    """Download sample FASTQ files with multi-source fallback.
    
    Tries sources in order: ENA FTP → ENA HTTP → NCBI SRA
    
    Args:
        sample_info: SampleInfo with all URL options
        sample_dir: Directory to save files
        timeout: Timeout per file in seconds
        
    Returns:
        True if all files downloaded successfully
    """
    sra_id = sample_info.sra_id
    sample_dir.mkdir(parents=True, exist_ok=True)
    
    # Strategy 1: ENA FTP (fastest, with MD5)
    if sample_info.ena_ftp_urls:
        logger.info(f"[{sra_id}] Trying ENA FTP...")
        success = True
        for url, md5 in sample_info.ena_ftp_urls:
            fname = Path(url).name
            dest = sample_dir / fname
            if not download_file(url, dest, md5, retries=2, timeout=timeout):
                success = False
                break
        if success:
            return True
        logger.warning(f"[{sra_id}] ENA FTP failed, trying HTTP...")
    
    # Strategy 2: ENA HTTP (backup)
    if sample_info.ena_http_urls:
        logger.info(f"[{sra_id}] Trying ENA HTTP...")
        success = True
        for url, md5 in sample_info.ena_http_urls:
            fname = Path(url).name
            dest = sample_dir / fname
            if not download_file(url, dest, md5, retries=2, timeout=timeout):
                success = False
                break
        if success:
            return True
        logger.warning(f"[{sra_id}] ENA HTTP failed, trying NCBI...")
    
    # Strategy 3: NCBI SRA (no MD5, requires fasterq-dump)
    if sample_info.ncbi_urls:
        logger.info(f"[{sra_id}] Trying NCBI SRA with fasterq-dump...")
        try:
            cmd = [
                "fasterq-dump",
                "--outdir", str(sample_dir),
                "--temp", str(sample_dir),
                "--threads", "2",
                "--progress",
                sra_id
            ]
            result = subprocess.run(cmd, capture_output=True, timeout=timeout)
            if result.returncode == 0:
                # Compress the output
                for fq in sample_dir.glob("*.fastq"):
                    subprocess.run(["gzip", str(fq)], check=True)
                return True
        except (subprocess.CalledProcessError, subprocess.TimeoutExpired) as e:
            logger.warning(f"[{sra_id}] NCBI fasterq-dump failed: {e}")
    
    return False


def sort_samples_by_size(
    sra_ids: List[str],
    max_size_bytes: Optional[int] = None
) -> List[SampleInfo]:
    """Query sizes and sort samples from smallest to largest.
    
    Args:
        sra_ids: List of SRA IDs to process
        max_size_bytes: Maximum file size to include (None = no limit)
        
    Returns:
        List of SampleInfo sorted by size, filtered by max_size
    """
    logger.info(f"Querying sizes for {len(sra_ids)} samples...")
    
    samples: List[SampleInfo] = []
    skipped_large = 0
    skipped_error = 0
    
    for i, sra_id in enumerate(sra_ids):
        if (i + 1) % 50 == 0:
            logger.info(f"  Progress: {i + 1}/{len(sra_ids)}")
        
        info = get_ena_sample_info(sra_id)
        if info is None:
            skipped_error += 1
            continue
            
        if max_size_bytes and info.size_bytes > max_size_bytes:
            logger.debug(f"[{sra_id}] Skipping: {format_size(info.size_bytes)} > {format_size(max_size_bytes)}")
            skipped_large += 1
            continue
            
        samples.append(info)
    
    # Sort by size
    samples.sort(key=lambda x: x.size_bytes)
    
    logger.info(f"Sample summary: {len(samples)} to download, {skipped_large} too large, {skipped_error} errors")
    if samples:
        logger.info(f"Size range: {format_size(samples[0].size_bytes)} - {format_size(samples[-1].size_bytes)}")
    
    return samples


def download_sra_samples(
    sra_ids: List[str], 
    base_out_dir: Path,
    max_size_bytes: Optional[int] = None,
    timeout: int = DEFAULT_TIMEOUT,
    sort_by_size: bool = False,
    use_fallback: bool = False
) -> Tuple[int, int]:
    """Download FASTQ files for a list of SRA samples from ENA.

    Creates amalgkit-compatible directory structure:
    base_out_dir/getfastq/SRA_ID/SRA_ID_1.amalgkit.fastq.gz

    Args:
        sra_ids: List of SRA accession IDs to download
        base_out_dir: Base output directory
        max_size_bytes: Skip samples larger than this (None = no limit)
        timeout: Timeout per file in seconds
        sort_by_size: If True, process smallest samples first
        use_fallback: If True, try NCBI/AWS when ENA fails

    Returns:
        Tuple of (success_count, fail_count)
    """
    base_out_dir = Path(base_out_dir)
    getfastq_dir = base_out_dir / "getfastq"
    getfastq_dir.mkdir(parents=True, exist_ok=True)

    # Get sample info and optionally sort
    if sort_by_size or max_size_bytes:
        samples = sort_samples_by_size(sra_ids, max_size_bytes)
    else:
        # Just get basic info without sorting
        samples = []
        for sra_id in sra_ids:
            info = get_ena_sample_info(sra_id)
            if info:
                samples.append(info)

    success_count = 0
    fail_count = 0
    
    total_size = sum(s.size_bytes for s in samples)
    logger.info(f"Starting download of {len(samples)} samples ({format_size(total_size)} total)")

    for i, sample in enumerate(samples):
        sra_id = sample.sra_id
        logger.info(f"[{i+1}/{len(samples)}] Processing {sra_id} ({format_size(sample.size_bytes)})...")

        sample_dir = getfastq_dir / sra_id
        
        if use_fallback:
            sample_success = download_with_fallback(sample, sample_dir, timeout)
        else:
            # Original ENA-only approach
            if not sample.ena_ftp_urls:
                logger.warning(f"[{sra_id}] No ENA links, skipping.")
                fail_count += 1
                continue
                
            sample_dir.mkdir(exist_ok=True)
            sample_success = True
            
            for url, md5 in sample.ena_ftp_urls:
                fname = Path(url).name
                dest = sample_dir / fname
                if not download_file(url, dest, md5, timeout=timeout):
                    sample_success = False
                    break

        if sample_success:
            # Rename to amalgkit format
            logger.info(f"Renaming files to .amalgkit.fastq.gz format...")
            for f in sample_dir.glob("*.fastq.gz"):
                if ".amalgkit." not in f.name:
                    new_name = f.name.replace(".fastq.gz", ".amalgkit.fastq.gz")
                    f.rename(sample_dir / new_name)
                    logger.debug(f"Renamed {f.name} -> {new_name}")

            success_count += 1
            logger.info(f"[SUCCESS] {sra_id} fully downloaded.")
        else:
            fail_count += 1
            logger.error(f"[FAILURE] {sra_id} failed to download.")

    logger.info(f"Summary: {success_count} Success, {fail_count} Failed.")
    return success_count, fail_count


def main() -> int:
    """CLI entry point for ENA downloader."""
    parser = argparse.ArgumentParser(
        description="ENA FASTQ Downloader with Best Practices",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic download
  python -m metainformant.rna.retrieval.ena_downloader --ids SRR12345 --out output/

  # Download with size limit, smallest first
  python -m metainformant.rna.retrieval.ena_downloader \\
    --file samples.txt --out output/ --max-size 5GB --sort-by-size

  # Full best practices (size-ordered, fallback enabled)  
  python -m metainformant.rna.retrieval.ena_downloader \\
    --file samples.txt --out output/ --max-size 5GB --timeout 1800 \\
    --sort-by-size --fallback
        """
    )
    parser.add_argument("--ids", nargs="+", help="List of SRA IDs to download")
    parser.add_argument("--file", help="File containing SRA IDs (one per line)")
    parser.add_argument("--out", required=True, help="Output base directory")
    parser.add_argument(
        "--max-size", 
        help="Maximum file size (e.g., 500MB, 5GB). Larger samples skipped."
    )
    parser.add_argument(
        "--timeout",
        type=int,
        default=DEFAULT_TIMEOUT,
        help=f"Download timeout per file in seconds (default: {DEFAULT_TIMEOUT})"
    )
    parser.add_argument(
        "--sort-by-size",
        action="store_true",
        help="Process smallest samples first"
    )
    parser.add_argument(
        "--fallback",
        action="store_true", 
        help="Enable fallback to NCBI/AWS when ENA fails"
    )

    args = parser.parse_args()

    id_list: List[str] = []
    if args.ids:
        id_list.extend(args.ids)
    if args.file:
        try:
            with open(args.file) as f:
                id_list.extend([line.strip() for line in f if line.strip()])
        except FileNotFoundError:
            logger.error(f"File {args.file} not found.")
            return 1

    if not id_list:
        logger.error("No IDs provided.")
        return 1

    max_size_bytes = parse_size(args.max_size) if args.max_size else None

    success, failed = download_sra_samples(
        id_list, 
        Path(args.out),
        max_size_bytes=max_size_bytes,
        timeout=args.timeout,
        sort_by_size=args.sort_by_size,
        use_fallback=args.fallback
    )
    return 0 if failed == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
