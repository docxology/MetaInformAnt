"""
ENA Downloader Module.

Handles direct FASTQ downloads from the European Nucleotide Archive (ENA).
This bypasses the NCBI SRA Toolkit (prefetch/fasterq-dump) bottleneck by
fetching pre-extracted .fastq.gz files directly via HTTP/FTP.
"""

import logging
import shutil
import subprocess
import urllib.error
import urllib.request
from pathlib import Path
from typing import List, Optional, Tuple

from metainformant.core.utils.logging import get_logger

logger = get_logger(__name__)

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

        Args:
            sample_id: SRA run accession (e.g., SRR1234567).

        Returns:
            List of HTTP URLs for the FASTQ files.
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
                ftp_field = lines[1].split("\t")[-1]
                
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

        Args:
            sample_id: Run accession (e.g., SRR1234567).
            output_dir: Directory to save files.

        Returns:
            Tuple of (success, message, list_of_downloaded_files).
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
                "--retry-delay", "5",
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
