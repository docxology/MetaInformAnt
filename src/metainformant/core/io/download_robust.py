"""Robust file download utilities for MetaInformAnt.

This module provides reliable download functions using available system tools
(curl, wget) with built-in retries, timeouts, and validation to handle
unstable network conditions or large files (e.g., SRA/FASTQ).
"""

import shutil
import subprocess
import time
from pathlib import Path
from typing import Dict, Optional

from metainformant.core.utils.logging import get_logger

logger = get_logger(__name__)


def get_remote_file_size(url: str) -> int:
    """Get the size of a remote file in bytes using curl/wget.

    Returns:
        Size in bytes, or 0 if unknown/failed.
    """
    # Try curl first
    if shutil.which("curl"):
        try:
            # -I/--head fetches header only, -L follows redirects
            cmd = ["curl", "-s", "-I", "-L", url]
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=10)
            if result.returncode == 0:
                for line in result.stdout.splitlines():
                    if line.lower().startswith("content-length:"):
                        try:
                            return int(line.split(":")[1].strip())
                        except ValueError:
                            pass
        except Exception as e:
            logger.debug(f"Failed to get size with curl for {url}: {e}")

    # Fallback to wget
    if shutil.which("wget"):
        try:
            # --spider checks file existence and size
            cmd = ["wget", "--spider", "--server-response", url]
            # Capture stderr since wget writes info there
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=10)
            # wget output parsing is messier, look for Content-Length in stderr
            for line in result.stderr.splitlines():
                if "Content-Length:" in line:
                    try:
                        # Format often: "  Content-Length: 12345 (12K)"
                        # Just take the first number
                        parts = line.split(":")[1].strip().split()
                        if parts:
                            return int(parts[0])
                    except ValueError:
                        pass
        except Exception:
            pass

    return 0


def robust_download_url(
    url: str, dest_path: Path, max_retries: int = 5, retry_delay: int = 10, timeout: int = 300
) -> bool:
    """Download a file from a URL robustly using curl or wget.

    Args:
        url: URL to download from
        dest_path: Destination file path
        max_retries: Number of retries
        retry_delay: Delay between retries in seconds
        timeout: Connection timeout in seconds

    Returns:
        True if successful, False otherwise
    """
    dest_path = Path(dest_path)
    dest_path.parent.mkdir(parents=True, exist_ok=True)
    temp_path = dest_path.with_suffix(dest_path.suffix + ".part")

    # Determined preferred tool
    has_curl = shutil.which("curl") is not None
    has_wget = shutil.which("wget") is not None

    if not has_curl and not has_wget:
        logger.error("No download tool (curl/wget) found!")
        return False

    for attempt in range(1, max_retries + 1):
        try:
            logger.info(f"Downloading {url} to {dest_path} (Attempt {attempt}/{max_retries})...")

            if has_curl:
                # curl: -L (follow redirects), -C - (resume), --retry (built-in retry), --fail (fail on 404)
                cmd = [
                    "curl",
                    "-L",
                    "-C",
                    "-",
                    "--fail",
                    "--retry",
                    str(3),
                    "--retry-delay",
                    str(5),
                    "--connect-timeout",
                    str(timeout),
                    "-o",
                    str(temp_path),
                    url,
                ]
            else:  # wget
                # wget: -c (continue), -t (tries), -T (timeout)
                cmd = ["wget", "-c", "-t", str(3), "-T", str(timeout), "-O", str(temp_path), url]

            result = subprocess.run(cmd, capture_output=True, text=True)

            if result.returncode == 0:
                # Success
                if temp_path.exists():
                    shutil.move(str(temp_path), str(dest_path))
                    logger.info(f"Download successful: {dest_path}")
                    return True
            else:
                logger.warning(f"Download failed (code {result.returncode}): {result.stderr}")

        except Exception as e:
            logger.warning(f"Download exception: {e}")

        if attempt < max_retries:
            logger.info(f"Retrying in {retry_delay} seconds...")
            time.sleep(retry_delay)

    return False


def download_sra_files_from_metadata(
    metadata_path: Path, output_base_dir: Path, max_count: Optional[int] = None
) -> Dict[str, bool]:
    """Download SRA files listed in metadata TSV using parallel TUI manager.

    Parses amalgkit metadata to find S3/FTP URLs and downloads them
    to {output_base_dir}/{sample_id}/{sample_id}.sra

    Args:
        metadata_path: Path to metadata.tsv
        output_base_dir: Base directory (e.g. output/fastq/getfastq)
        max_count: Limit number of downloads (for testing)

    Returns:
        Dictionary of sample_id -> success status
    """
    import csv

    from metainformant.core.io.download_manager import DownloadManager

    results = {}
    if not metadata_path.exists():
        logger.error(f"Metadata file not found: {metadata_path}")
        return results

    manager = DownloadManager(max_threads=5)
    queued_count = 0

    with open(metadata_path, "r") as f:
        reader = csv.DictReader(f, delimiter="\t")

        for row in reader:
            if max_count and queued_count >= max_count:
                break

            sample_id = row.get("run") or row.get("Run")
            if not sample_id:
                continue

            # Construct target path
            sample_dir = output_base_dir / sample_id
            target_file = sample_dir / f"{sample_id}.sra"

            # Also check backup directory usually used by workflow manager
            backup_file = output_base_dir / "sra" / f"{sample_id}.sra"

            if (target_file.exists() and target_file.stat().st_size > 1000) or (
                backup_file.exists() and backup_file.stat().st_size > 1000
            ):
                logger.info(f"File already exists (original or backup): {sample_id}")
                results[sample_id] = True
                continue

            # Find best URL
            # 1. Trust AWS_Link if present (best for amalgkit metadata)
            # 2. Trust sra_url if present
            # 3. Fallback to constructed S3 URL

            url = ""
            if row.get("AWS_Link"):
                url = row["AWS_Link"]
            elif row.get("sra_url"):
                overrides = row["sra_url"].split(";")
                if overrides:
                    # Filter out sralite
                    valid = [u for u in overrides if "sralite" not in u]
                    if valid:
                        url = valid[0]

            if not url:
                url = f"https://sra-pub-run-odp.s3.amazonaws.com/sra/{sample_id}/{sample_id}"

            manager.add_download(url, target_file, label=sample_id)
            queued_count += 1

    # Execute all
    if queued_count > 0:
        logger.info(f"Starting parallel download of {queued_count} files...")
        batch_results = manager.start()
        results.update(batch_results)

    return results
