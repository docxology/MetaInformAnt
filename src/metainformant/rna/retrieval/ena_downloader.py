"""ENA (European Nucleotide Archive) FASTQ downloader.

This module provides functionality to download FASTQ files from the European
Nucleotide Archive (ENA) for RNA-seq samples. It handles URL resolution,
MD5 verification, and resume-capable downloads.

Example usage:
    >>> from metainformant.rna.retrieval.ena_downloader import download_sra_samples
    >>> download_sra_samples(["SRR12345"], Path("output/fastq"))
"""

from __future__ import annotations

import argparse
import hashlib
import os
import subprocess
import sys
import time
from pathlib import Path
from typing import List, Tuple

import requests

from metainformant.core import logging

logger = logging.get_logger(__name__)


def get_ena_links(sra_id: str) -> List[Tuple[str, str]]:
    """Query ENA to get FASTQ download links and MD5 checksums.

    Args:
        sra_id: SRA accession ID (e.g., SRR12345, ERR12345, DRR12345)

    Returns:
        List of (url, md5) tuples for each FASTQ file associated with the sample.
        Empty list if no data found or on error.
    """
    url = (
        f"https://www.ebi.ac.uk/ena/portal/api/filereport?accession={sra_id}"
        f"&result=read_run&fields=fastq_ftp,fastq_md5&format=tsv"
    )

    try:
        response = requests.get(url, timeout=10)
        response.raise_for_status()

        lines = response.text.strip().split("\n")
        if len(lines) < 2:
            logger.warning(f"[{sra_id}] No data found in ENA.")
            return []

        header = lines[0].split("\t")
        data = lines[1].split("\t")

        # Parse TSV map
        row = dict(zip(header, data))

        if "fastq_ftp" not in row or not row["fastq_ftp"]:
            logger.warning(f"[{sra_id}] No FASTQ FTP links returned.")
            return []

        ftps = row["fastq_ftp"].split(";")
        md5s = row["fastq_md5"].split(";")

        if len(ftps) != len(md5s):
            logger.warning(f"[{sra_id}] Mismatch between FTP links and MD5s.")
            return []

        results = []
        for ftp, md5 in zip(ftps, md5s):
            # ENA links often lack protocol prefix (just "ftp.sra.ebi.ac.uk/...")
            # Check for actual protocol prefixes (ftp://, http://, https://)
            if not ftp.startswith(("ftp://", "http://", "https://")):
                # Prefer FTP as HTTP seems to hang on port 80 for this host
                full_url = f"ftp://{ftp}"
            else:
                full_url = ftp
            results.append((full_url, md5))

        return results

    except requests.exceptions.RequestException as e:
        logger.error(f"[{sra_id}] API Request failed: {e}")
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
        for chunk in iter(lambda: f.read(4096), b""):
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


def download_file(url: str, dest_path: Path, expected_md5: str, retries: int = 5) -> bool:
    """Download file using curl with retries and MD5 verification.

    Args:
        url: URL to download from
        dest_path: Local path to save the file
        expected_md5: Expected MD5 checksum for verification
        retries: Number of retry attempts

    Returns:
        True if download succeeded and MD5 verified, False otherwise
    """
    logger.info(f"Downloading {url} -> {dest_path}")

    for attempt in range(1, retries + 1):
        try:
            # Use curl for robustness
            # -C - : Continue/Resume
            # --connect-timeout 60 : Don't wait forever for connection
            # --retry 3 : Internal curl retries
            cmd = [
                "curl",
                "-L",
                "-C",
                "-",
                "-o",
                str(dest_path),
                "--retry",
                "3",
                "--retry-delay",
                "5",
                "--connect-timeout",
                "60",
                "--keepalive-time",
                "60",
                url,
            ]

            # Allow stderr to be seen for debugging
            subprocess.run(cmd, check=True)

            if dest_path.exists():
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
                    # If MD5 fails, the file is likely corrupted (wrong offset or bad network).
                    # Delete it to start fresh next time.
                    clean_stagnant_file(dest_path)

        except subprocess.CalledProcessError as e:
            logger.warning(f"[FAIL] Download command failed (Attempt {attempt}/{retries}). Exit code: {e.returncode}")

        time.sleep(5 * attempt)

    return False


def download_sra_samples(sra_ids: List[str], base_out_dir: Path) -> Tuple[int, int]:
    """Download FASTQ files for a list of SRA samples from ENA.

    Creates amalgkit-compatible directory structure:
    base_out_dir/getfastq/SRA_ID/SRA_ID_1.amalgkit.fastq.gz

    Args:
        sra_ids: List of SRA accession IDs to download
        base_out_dir: Base output directory

    Returns:
        Tuple of (success_count, fail_count)
    """
    base_out_dir = Path(base_out_dir)
    getfastq_dir = base_out_dir / "getfastq"
    getfastq_dir.mkdir(parents=True, exist_ok=True)

    success_count = 0
    fail_count = 0

    for sra_id in sra_ids:
        logger.info(f"Processing {sra_id}...")

        links = get_ena_links(sra_id)
        if not links:
            logger.warning(f"Skipping {sra_id} (No links found).")
            fail_count += 1
            continue

        # Prepare sample directory
        sample_dir = getfastq_dir / sra_id
        sample_dir.mkdir(exist_ok=True)

        sample_success = True

        for i, (url, md5) in enumerate(links):
            # Amalgkit naming convention: SRR..._1.fastq.gz, SRR..._2.fastq.gz
            # ENA files usually end in _1.fastq.gz or _2.fastq.gz
            fname = Path(url).name

            # Map ENA filenames to Amalgkit expectation if needed
            # ENA: SRR..._1.fastq.gz -> Amalgkit: SRR..._1.fastq.gz (matches)
            # ENA: SRR....fastq.gz -> Amalgkit: SRR....fastq.gz (single)

            dest = sample_dir / fname

            # Amalgkit final file often has .amalgkit.fastq.gz
            # But the 'getfastq' step often produces just the raw fastq first.
            # We will download the raw files.

            if not download_file(url, dest, md5):
                sample_success = False
                break

        if sample_success:
            # Create a marker file or symlink to .amalgkit.fastq.gz to trick amalgkit?
            # Ideally, we just rename them to what amalgkit expects for the 'integrate' step.
            # Amalgkit integrate looks for *.amalgkit.fastq.gz

            logger.info(f"Renaming files to .amalgkit.fastq.gz format for compatibility...")
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
    parser = argparse.ArgumentParser(description="Modular ENA Downloader")
    parser.add_argument("--ids", nargs="+", help="List of SRA IDs to download")
    parser.add_argument("--file", help="File containing SRA IDs (one per line)")
    parser.add_argument(
        "--out", required=True, help="Output base directory (e.g. output/amalgkit/pbarbatus_all/fastq)"
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

    success, failed = download_sra_samples(id_list, Path(args.out))
    return 0 if failed == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
