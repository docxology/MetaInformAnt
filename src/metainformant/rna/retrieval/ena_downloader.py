"""
ENA Downloader Module.

Handles direct FASTQ downloads from the European Nucleotide Archive (ENA).
This bypasses the NCBI SRA Toolkit (prefetch/fasterq-dump) bottleneck by
fetching pre-extracted .fastq.gz files directly via HTTP/FTP.

This module provides a high-performance alternative to traditional SRA
download methods, particularly useful for:
- High-throughput RNA-seq projects (1000+ samples)
- Environments where SRA Toolkit installation is problematic
- Projects requiring specific ENA samples not yet synced to SRA
- Parallel batch downloads with resume support

Example:
    >>> from metainformant.rna.retrieval.ena_downloader import ENADownloader
    >>> downloader = ENADownloader(timeout=1800, retries=3)
    >>> success, msg, files = downloader.download_run("SRR1234567", Path("output/fastq"))
    >>> print(f"Download {msg}")
    Download downloaded 2 files
"""

import gzip
import hashlib
import logging
import shutil
import subprocess
import urllib.error
import urllib.request
from pathlib import Path
from typing import List, Optional, Tuple

from metainformant.core.utils.logging import get_logger

logger = get_logger(__name__)


def calculate_md5(file_path: Path, chunk_size: int = 4096) -> str:
    """Calculate MD5 checksum of a file.

    Args:
        file_path: Path to the file.
        chunk_size: Size of chunks for reading.

    Returns:
        Hex digest of the MD5 checksum.
    """
    md5_hash = hashlib.md5()
    with open(file_path, "rb") as f:
        for chunk in iter(lambda: f.read(chunk_size), b""):
            md5_hash.update(chunk)
    return md5_hash.hexdigest()


def clean_stagnant_file(file_path: Path) -> None:
    """Remove a stagnant/incomplete download file if it exists.

    Args:
        file_path: Path to the file to remove.
    """
    file_path = Path(file_path)
    if file_path.exists():
        file_path.unlink()
        logger.info(f"Cleaned stagnant file: {file_path}")


def verify_gzip_integrity(file_path: Path) -> bool:
    """Verify the integrity of a gzip file.

    Args:
        file_path: Path to the file to verify.

    Returns:
        True if the file is valid gzip or not a .gz file, False if corrupted.
    """
    file_path = Path(file_path)
    if not file_path.suffix == ".gz":
        return True  # Non-gz files are assumed valid

    try:
        with gzip.open(file_path, "rb") as f:
            while True:
                chunk = f.read(8192)
                if not chunk:
                    break
        return True
    except (gzip.BadGzipFile, OSError, EOFError):
        return False

class ENADownloader:
    """
    Downloads FASTQ files directly from ENA.
    """
    
    ENA_HTTP_BASE = "http://ftp.sra.ebi.ac.uk/vol1/fastq"

    def __init__(self, timeout: int = 1800, retries: int = 3):
        """
        Initialize the downloader.

        Args:
            timeout: Maximum download time in seconds (default: 1800/30mins).
            retries: Number of curl retries (default: 3).
        """
        self.timeout = timeout
        self.retries = retries

    def get_fastq_urls(self, sample_id: str) -> List[str]:
        """
        Discover FASTQ URLs for a sample using the ENA Portal API.

        Queries the ENA Portal API to discover available FASTQ files for a
        given SRA run accession. Converts FTP paths to HTTP URLs for use
        with curl.

        Args:
            sample_id: SRA run accession (e.g., SRR1234567, ERR1234567).
                Must be a valid run accession in the ENA database.

        Returns:
            List of HTTP URLs for the FASTQ files. Returns empty list if
            no FASTQ files are found or the sample doesn't exist.

        Raises:
            urllib.error.URLError: If the API request fails

        Example:
            >>> downloader = ENADownloader()
            >>> urls = downloader.get_fastq_urls("SRR1234567")
            >>> print(f"Found {len(urls)} FASTQ files")
            Found 2 FASTQ files

        Note:
            The ENA Portal API returns FTP paths which are converted to
            HTTP URLs for curl compatibility. Single-end runs return
            one URL; paired-end runs return two URLs (forward and reverse).
        """
        api_url = (
            f"https://www.ebi.ac.uk/ena/portal/api/filereport?"
            f"accession={sample_id}&result=read_run&fields=fastq_ftp"
        )

        try:
            with urllib.request.urlopen(api_url, timeout=30) as response:
                content = response.read().decode("utf-8")
                lines = content.strip().split("\n")

                if len(lines) < 2:
                    return []

                # Second line contains the FTP URLs (semicolon-separated)
                # Header line is: run_accession	fastq_ftp
                parts = lines[1].split("\t")
                if len(parts) < 2:
                    return []
                
                ftp_field = parts[1].strip()
                
                # Check if field is valid/not empty/not just header
                if not ftp_field or ftp_field == "fastq_ftp":
                    return []

                # Convert FTP paths to HTTP URLs for curl
                urls = []
                for ftp_path in ftp_field.split(";"):
                    ftp_path = ftp_path.strip()
                    if ftp_path:
                        # ENA API returns "ftp.sra...", so prepend http://
                        http_url = f"http://{ftp_path}"
                        urls.append(http_url)

                return urls
        except Exception as e:
            logger.warning(f"ENA API query failed for {sample_id}: {e}")
            return []

    def download_run(self, sample_id: str, output_dir: Path) -> Tuple[bool, str, List[Path]]:
        """
        Download all FASTQ files for a run.

        Downloads all FASTQ files associated with a run accession using
        curl with automatic retry. Checks for existing files to enable
        resume functionality.

        Args:
            sample_id: Run accession (e.g., SRR1234567, ERR1234567).
                Must be a valid SRA/ENA run accession.
            output_dir: Directory to save downloaded files. Created if
                it doesn't exist.

        Returns:
            Tuple of (success, message, list_of_downloaded_files):
            - success: True if all files downloaded successfully
            - message: Status message (e.g., "Downloaded 2 files", error description)
            - list_of_downloaded_files: List of Path objects for downloaded files

        Raises:
            subprocess.TimeoutExpired: If download exceeds timeout

        Example:
            >>> from pathlib import Path
            >>> from metainformant.rna.retrieval.ena_downloader import ENADownloader
            >>>
            >>> downloader = ENADownloader(timeout=3600, retries=3)
            >>> success, msg, files = downloader.download_run(
            ...     "SRR1234567",
            ...     Path("output/fastq/SRR1234567")
            ... )
            >>> if success:
            ...     print(f"Downloaded: {[f.name for f in files]}")

        Note:
            - Existing non-empty files are skipped (resume support)
            - Uses curl with retry logic for reliability
            - Cleans up partial downloads on failure
        """
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        downloaded_files = []

        # 1. Discover URLs
        urls = self.get_fastq_urls(sample_id)
        if not urls:
            return False, "Not found on ENA", []

        # 2. Download each file
        for url in urls:
            filename = url.split("/")[-1]
            output_file = output_dir / filename

            # Check if already exists and non-empty
            if output_file.exists() and output_file.stat().st_size > 0:
                downloaded_files.append(output_file)
                continue

            # Use curl for reliability
            cmd = [
                "curl",
                "-fsSL",
                "--retry", str(self.retries),
                "--retry-delay", "10",
                "--retry-connrefused",
                "--retry-all-errors",
                "--connect-timeout", "30",
                "-o", str(output_file),
                url
            ]

            try:
                result = subprocess.run(
                    cmd, 
                    capture_output=True, 
                    text=True, 
                    timeout=self.timeout
                )

                if result.returncode == 0 and output_file.exists() and output_file.stat().st_size > 0:
                    downloaded_files.append(output_file)
                else:
                    # Clean up partial
                    if output_file.exists():
                        output_file.unlink()
                    return False, f"Download failed for {filename}: {result.stderr}", []
            
            except subprocess.TimeoutExpired:
                if output_file.exists():
                    output_file.unlink()
                return False, "Download timed out", []
            except Exception as e:
                return False, f"Download error: {str(e)}", []

        return True, f"Downloaded {len(downloaded_files)} files", downloaded_files
