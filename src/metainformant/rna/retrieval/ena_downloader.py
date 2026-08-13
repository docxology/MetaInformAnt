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

import csv
import gzip
import hashlib
import os
import signal
import subprocess
import time
import urllib.error
import urllib.request
from pathlib import Path
from typing import List, Mapping, Tuple

from metainformant.core.utils.logging import get_logger

logger = get_logger(__name__)

_RETRYABLE_CURL_EXIT_CODES = frozenset({6, 7, 18, 28, 35, 52, 55, 56, 92})
_RETRYABLE_HTTP_STATUS_CODES = frozenset({403, 408, 425, 429, 500, 502, 503, 504})


def _run_command_in_process_group(
    command: list[str],
    timeout: int,
    *,
    env: Mapping[str, str] | None = None,
) -> subprocess.CompletedProcess[str]:
    """Run a transfer command and terminate its descendants on timeout.

    ``curl`` normally has no descendants, but keeping all acquisition commands
    in their own process group makes timeout recovery uniform and prevents a
    future wrapper or retry helper from surviving as an orphan. Partial files
    are intentionally left to the caller so a later run can resume them.
    """

    process = subprocess.Popen(
        command,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        start_new_session=True,
        env=dict(env) if env is not None else None,
    )
    try:
        stdout, stderr = process.communicate(timeout=timeout)
    except subprocess.TimeoutExpired as exc:
        try:
            os.killpg(process.pid, signal.SIGTERM)
        except (OSError, ProcessLookupError):
            process.terminate()
        try:
            stdout, stderr = process.communicate(timeout=10)
        except subprocess.TimeoutExpired:
            try:
                os.killpg(process.pid, signal.SIGKILL)
            except (OSError, ProcessLookupError):
                process.kill()
            stdout, stderr = process.communicate()
        raise subprocess.TimeoutExpired(
            command,
            timeout,
            output=stdout if stdout else exc.output,
            stderr=stderr if stderr else exc.stderr,
        ) from exc
    return subprocess.CompletedProcess(command, process.returncode, stdout, stderr)


def _curl_http_status(result: subprocess.CompletedProcess[str]) -> int | None:
    """Return the HTTP status emitted by curl's ``--write-out`` option."""

    output = result.stdout.strip()
    if len(output) >= 3 and output[-3:].isdigit():
        return int(output[-3:])
    return None


def _is_retryable_transfer_failure(result: subprocess.CompletedProcess[str]) -> bool:
    """Classify failures that can safely resume the retained partial file."""

    if result.returncode in _RETRYABLE_CURL_EXIT_CODES:
        return True
    if result.returncode == 22:
        return _curl_http_status(result) in _RETRYABLE_HTTP_STATUS_CODES
    return False


def calculate_md5(file_path: Path, chunk_size: int = 4096) -> str:
    """Calculate MD5 checksum of a file.

    Args:
        file_path: Path to the file.
        chunk_size: Size of chunks for reading.

    Returns:
        Hex digest of the MD5 checksum.
    """
    # MD5 is used only for transfer-integrity comparison, never authentication.
    md5_hash = hashlib.md5(usedforsecurity=False)
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


def preserve_invalid_file(file_path: Path) -> Path:
    """Move an invalid transfer aside without overwriting prior evidence."""
    file_path = Path(file_path)
    candidate = file_path.with_name(f"{file_path.name}.invalid")
    counter = 1
    while candidate.exists():
        candidate = file_path.with_name(f"{file_path.name}.invalid.{counter}")
        counter += 1
    file_path.replace(candidate)
    logger.warning("Preserved invalid transfer as %s", candidate)
    return candidate


def record_invalid_transfer(
    file_path: Path,
    *,
    reason: str,
    max_witness_bytes: int = 16 * 1024 * 1024,
) -> Path | None:
    """Record a corrupt payload and retain at most one full witness per file.

    Repeated ENA retries can otherwise consume more space than the valid
    dataset. The first invalid payload is retained; later payloads are hashed,
    recorded, and removed. The manifest is append-only and local to the sample
    directory so it remains available after the raw FASTQ is reclaimed.
    """

    file_path = Path(file_path)
    if not file_path.is_file():
        return None
    try:
        size_bytes = file_path.stat().st_size
        digest = calculate_md5(file_path)
    except OSError as exc:
        logger.warning("Unable to fingerprint invalid transfer %s: %s", file_path, exc)
        size_bytes = 0
        digest = ""

    prior_witnesses = list(file_path.parent.glob(f"{file_path.name}.invalid*"))
    retained_path: Path | None = None
    if not prior_witnesses and size_bytes <= max_witness_bytes:
        retained_path = preserve_invalid_file(file_path)
    else:
        file_path.unlink(missing_ok=True)
        logger.warning("Discarded duplicate invalid transfer after fingerprinting: %s", file_path)

    manifest = file_path.parent / ".transfer_integrity_failures.tsv"
    needs_header = not manifest.exists()
    with manifest.open("a", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["source_name", "size_bytes", "md5", "reason", "retained_witness", "recorded_at_utc"],
            delimiter="\t",
        )
        if needs_header:
            writer.writeheader()
        writer.writerow(
            {
                "source_name": file_path.name,
                "size_bytes": size_bytes,
                "md5": digest,
                "reason": reason,
                "retained_witness": retained_path.name if retained_path else "",
                "recorded_at_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
            }
        )
        handle.flush()
        os.fsync(handle.fileno())
    return retained_path


def verify_gzip_integrity(file_path: Path) -> bool:
    """Verify the integrity of a gzip file.

    Args:
        file_path: Path to the file to verify.

    Returns:
        True if the file is valid gzip or not a .gz file, False if corrupted.
    """
    file_path = Path(file_path)
    if not (file_path.name.endswith(".gz") or file_path.name.endswith(".gz.part")):
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

    ENA_HTTP_BASE = "https://ftp.sra.ebi.ac.uk/vol1/fastq"

    def __init__(
        self,
        timeout: int = 1800,
        retries: int = 5,
        integrity_retries: int = 1,
        invalid_witness_max_bytes: int = 16 * 1024 * 1024,
        retry_delay_seconds: int = 5,
        speed_limit_bytes: int = 1024,
        speed_time_seconds: int = 600,
        api_retries: int = 2,
        api_retry_delay_seconds: int = 2,
        api_timeout_seconds: int = 30,
    ):
        """
        Initialize the downloader.

        Args:
            timeout: Maximum download time in seconds (default: 1800/30mins).
            retries: Maximum consecutive resumable attempts that make no byte
                progress (default: 5).
            integrity_retries: Maximum fresh-transfer retries after a completed
                payload fails gzip integrity (default: 1). This is separate
                from resumable transport retries so corrupt endpoints fall
                through to NCBI promptly.
            invalid_witness_max_bytes: Maximum size of one retained invalid
                payload witness per FASTQ; larger payloads are fingerprinted
                and removed (default: 16 MiB).
                Productive premature closes reset the no-progress budget and
                continue within ``timeout``. Each retry starts a new curl
                invocation so ``--continue-at -`` advances from the retained
                ``.part`` file.
            retry_delay_seconds: Initial retry delay in seconds (default: 5).
                The delay doubles between attempts up to 60 seconds.
            speed_limit_bytes: Minimum sustained transfer rate before curl
                considers a transfer stalled (default: 1024 bytes/second).
            speed_time_seconds: Seconds below ``speed_limit_bytes`` before
                curl aborts the transfer and preserves its ``.part`` file
                (default: 600 seconds).
            api_retries: Number of retries for transient ENA API failures
                (default: 2).
            api_retry_delay_seconds: Seconds between ENA API retries
                (default: 2).
            api_timeout_seconds: Timeout for each ENA metadata request
                (default: 30). Keep this shorter than an enclosing test or
                workflow watchdog when operating in bounded environments.
        """
        if timeout <= 0:
            raise ValueError("timeout must be positive")
        if retries < 0:
            raise ValueError("retries must be non-negative")
        if integrity_retries < 0:
            raise ValueError("integrity_retries must be non-negative")
        if invalid_witness_max_bytes < 0:
            raise ValueError("invalid_witness_max_bytes must be non-negative")
        if retry_delay_seconds <= 0:
            raise ValueError("retry_delay_seconds must be positive")
        if speed_limit_bytes <= 0:
            raise ValueError("speed_limit_bytes must be positive")
        if speed_time_seconds <= 0:
            raise ValueError("speed_time_seconds must be positive")
        if api_retries < 0:
            raise ValueError("api_retries must be non-negative")
        if api_retry_delay_seconds <= 0:
            raise ValueError("api_retry_delay_seconds must be positive")
        if api_timeout_seconds <= 0:
            raise ValueError("api_timeout_seconds must be positive")
        self.timeout = timeout
        self.retries = retries
        self.integrity_retries = integrity_retries
        self.invalid_witness_max_bytes = invalid_witness_max_bytes
        self.retry_delay_seconds = retry_delay_seconds
        self.speed_limit_bytes = speed_limit_bytes
        self.speed_time_seconds = speed_time_seconds
        self.api_retries = api_retries
        self.api_retry_delay_seconds = api_retry_delay_seconds
        self.api_timeout_seconds = api_timeout_seconds

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

        for attempt in range(self.api_retries + 1):
            try:
                with urllib.request.urlopen(api_url, timeout=self.api_timeout_seconds) as response:  # nosec B310
                    content = response.read().decode("utf-8")
                    lines = content.strip().split("\n")

                    if len(lines) < 2:
                        return []

                    # Second line contains the FTP URLs (semicolon-separated)
                    # Header line is: run_accession\tfastq_ftp
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
                            # ENA API returns "ftp.sra...".  Use HTTPS so that
                            # large transfers are not downgraded to a less auditable
                            # clear-text endpoint.
                            urls.append(f"https://{ftp_path}")

                    return urls
            except Exception as exc:
                if attempt >= self.api_retries:
                    logger.warning(
                        "ENA API query failed for %s after %d retries: %s",
                        sample_id,
                        self.api_retries,
                        exc,
                    )
                    return []
                logger.info(
                    "ENA API query failed for %s (%s); retrying in %ss (%d/%d)",
                    sample_id,
                    exc,
                    self.api_retry_delay_seconds,
                    attempt + 1,
                    self.api_retries,
                )
                time.sleep(self.api_retry_delay_seconds)
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

            partial_file = output_file.with_name(f"{output_file.name}.part")

            # Check if already exists and non-empty. Resume support must still
            # reject corrupt gzip files from interrupted prior downloads.
            if output_file.exists() and output_file.stat().st_size > 0:
                if verify_gzip_integrity(output_file):
                    downloaded_files.append(output_file)
                    continue
                logger.warning(
                    "Existing FASTQ failed gzip integrity check; preserving it and starting a fresh transfer: %s",
                    output_file,
                )
                record_invalid_transfer(
                    output_file,
                    reason="existing_invalid_gzip",
                    max_witness_bytes=self.invalid_witness_max_bytes,
                )

            # A previous interrupted curl is retained as ``.part``.  ENA
            # advertises byte ranges, so curl can continue without repeating
            # the already acquired portion of a multi-gigabyte FASTQ.

            # Run each retry as a separate curl invocation. Curl's
            # ``--retry-all-errors`` restarts a failed transfer within the
            # same invocation and can discard bytes received during that
            # attempt. Separate invocations allow ``--continue-at -`` to
            # advance from the retained partial after exit 18/52 failures.
            cmd = [
                "curl",
                "-fsSL",
                "--connect-timeout",
                "30",
                "--speed-limit",
                str(self.speed_limit_bytes),
                "--speed-time",
                str(self.speed_time_seconds),
                "--continue-at",
                "-",
                "-o",
                str(partial_file),
                "--write-out",
                "%{http_code}",
                url,
            ]

            deadline = time.monotonic() + self.timeout
            consecutive_no_progress = 0
            retry_events = 0
            gzip_integrity_retries = 0
            while True:
                remaining_seconds = deadline - time.monotonic()
                if remaining_seconds <= 0:
                    return False, f"Download timed out; retained partial {partial_file.name}", []

                size_before = partial_file.stat().st_size if partial_file.exists() else 0
                try:
                    result = _run_command_in_process_group(
                        cmd,
                        timeout=max(1, int(remaining_seconds)),
                    )
                except subprocess.TimeoutExpired:
                    return False, f"Download timed out; retained partial {partial_file.name}", []
                except Exception as e:
                    return False, f"Download error: {str(e)}", []

                if (
                    result.returncode == 0
                    and partial_file.exists()
                    and partial_file.stat().st_size > 0
                    and verify_gzip_integrity(partial_file)
                ):
                    partial_file.replace(output_file)
                    downloaded_files.append(output_file)
                    break

                if result.returncode == 0:
                    if partial_file.exists():
                        record_invalid_transfer(
                            partial_file,
                            reason="completed_transfer_invalid_gzip",
                            max_witness_bytes=self.invalid_witness_max_bytes,
                        )
                    if output_file.exists():
                        record_invalid_transfer(
                            output_file,
                            reason="completed_output_invalid_gzip",
                            max_witness_bytes=self.invalid_witness_max_bytes,
                        )
                    gzip_integrity_retries += 1
                    if gzip_integrity_retries <= self.integrity_retries:
                        logger.warning(
                            "ENA transfer for %s failed gzip integrity; starting fresh retry "
                            "%d/%d while preserving the invalid payload",
                            filename,
                            gzip_integrity_retries,
                            self.integrity_retries,
                        )
                        consecutive_no_progress = 0
                        continue
                    return False, f"Download failed gzip integrity check for {filename}", []

                if not _is_retryable_transfer_failure(result):
                    return False, f"Download failed for {filename}: {result.stderr}", []

                retained_bytes = partial_file.stat().st_size if partial_file.exists() else 0
                gained_bytes = max(0, retained_bytes - size_before)
                if gained_bytes > 0:
                    consecutive_no_progress = 0
                else:
                    consecutive_no_progress += 1
                if consecutive_no_progress > self.retries:
                    return (
                        False,
                        f"Download retry budget exhausted after {self.retries} "
                        f"consecutive no-progress failures; retained partial {partial_file.name}",
                        [],
                    )

                retry_events += 1
                delay_seconds = min(
                    self.retry_delay_seconds * (2 ** min(max(consecutive_no_progress - 1, 0), 4)),
                    60,
                )
                remaining_seconds = deadline - time.monotonic()
                if remaining_seconds <= delay_seconds:
                    return (
                        False,
                        f"Download retry budget exhausted; retained partial {partial_file.name}",
                        [],
                    )
                logger.info(
                    "Retrying resumable ENA transfer for %s after curl %d / HTTP %s "
                    "(retry event %d; consecutive no-progress failures %d/%d; "
                    "retained %.2f MiB, gained %.2f MiB, backoff %ds)",
                    filename,
                    result.returncode,
                    _curl_http_status(result) or "unknown",
                    retry_events,
                    consecutive_no_progress,
                    self.retries,
                    retained_bytes / 1024 / 1024,
                    gained_bytes / 1024 / 1024,
                    delay_seconds,
                )
                time.sleep(delay_seconds)

        return True, f"Downloaded {len(downloaded_files)} files", downloaded_files


def download_sra_samples(
    sra_ids: List[str],
    base_out_dir: Path,
    *,
    sort_by_size: bool = True,
    use_fallback: bool = True,
    downloader: ENADownloader | None = None,
) -> Tuple[int, int]:
    """Download a batch of SRA run accessions into amalgkit getfastq layout.

    Args:
        sra_ids: Run accessions to download.
        base_out_dir: Amalgkit work directory. Files are written under
            ``base_out_dir/getfastq/{run_id}/``.
        sort_by_size: Reserved for compatibility with callers that may pass
            size-sorted IDs. The current ENA API path does not pre-fetch sizes.
        use_fallback: Reserved for compatibility; direct ENA download is the
            only implementation in this helper.
        downloader: Optional downloader instance for tests or custom timeouts.

    Returns:
        Tuple of ``(success_count, fail_count)``.
    """
    _ = sort_by_size, use_fallback
    base_out_dir = Path(base_out_dir)
    getfastq_dir = base_out_dir / "getfastq"
    getfastq_dir.mkdir(parents=True, exist_ok=True)
    downloader = downloader or ENADownloader()

    success_count = 0
    fail_count = 0
    for sra_id in sra_ids:
        sample_dir = getfastq_dir / sra_id
        success, message, _files = downloader.download_run(sra_id, sample_dir)
        if success:
            success_count += 1
        else:
            fail_count += 1
            logger.warning(f"Failed to download {sra_id}: {message}")

    return success_count, fail_count
