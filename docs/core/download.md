# Core: Downloads with Progress and Heartbeat

Centralized, robust download utilities for METAINFORMANT with real-time progress tracking, heartbeat monitoring, retry logic, and resume capability. The `download` module provides both high-level functions and low-level handlers for HTTP(S), FTP, and local file transfers.

## Purpose

Bioinformatics workflows frequently download large datasets:
- **Reference genomes**: Gigabytes to terabytes (SRA, ENA, EBI)
- **Public datasets**: SRA run accessions, GEO series, dbSNP variants
- **External tools**: Binaries, indexes (BWA, STAR, kallisto)
- **Intermediate results**: Pre-computed alignments, annotations

These downloads need to be:
- **Observable**: Long-running (hours for terabyte-scale transfers) need progress indication
- **Robust**: Network failures should not require restart from scratch
- **Monitorable**: External tools (e.g., `amalgkit getfastq`) may handle downloads internally but still need heartbeat files for orchestrators
- **Resumable**: Support HTTP Range requests to continue interrupted transfers
- **Rate-limited**: Some servers require throttling

The `download` module centralizes this logic to ensure consistent behavior across all METAINFORMANT domains.

## Design Principles

### 1. **Heartbeat Files**
Every download writes a machine-readable JSON heartbeat file every `heartbeat_interval` seconds (default 5) to `dest.parent/.downloads/<filename>.heartbeat.json`. The heartbeat contains:
- Current progress (bytes downloaded, percent complete)
- Transfer speed (MB/s) and ETA
- Current status (`starting`, `downloading`, `completed`, `failed`)
- Error accumulation for troubleshooting

Orchestrators (workflow engines, TUI) poll heartbeat files to update UI without needing to wrap download calls.

### 2. **Progress Bars Integrated**
If `tqdm` is available via `metainformant.core.progress`, a progress bar is shown by default. Downloads can disable progress for background/cron runs with `show_progress=False`.

### 3. **Retry with Exponential Backoff**
Transient network errors trigger automatic retry (default 3 attempts) with delays of 5s, 10s, 20s. Errors are accumulated in the heartbeat file for post-mortem analysis.

### 4. **Best-Effort Resume**
HTTP(S) downloads support resume via `Range` headers. If the destination file already exists and `resume=True` (default), the download continues from the existing size. Non-resumable protocols (FTP, file) restart from beginning.

### 5. **Protocol Abstraction**
`DownloadHandler` Protocol defines a common interface. Built-in handlers:
- `HTTPDownloadHandler` — requests-based with Range support
- `FTPDownloadHandler` — urllib-based
- `FileDownloadHandler` — local file:// copies

### 6. **Subprocess Monitoring**
External tools (e.g., SRA Toolkit's `prefetch`) that manage their own downloads can be monitored via:
- `monitor_subprocess_directory_growth()` — Watch bytes written in a directory
- `monitor_subprocess_file_count()` — Count expected output files
- `monitor_subprocess_sample_progress()` — Track sample completion by glob pattern

These monitor functions write heartbeat files just like direct downloads, unifying the monitoring interface.

### 7. **Atomic Comms**
Heartbeat files are written atomically via tmp-rename to avoid partial reads by consumers. Directory size calculations are cached (TTL 1s) to avoid excessive `stat()` calls during fast-growing downloads.

### 8. **No Global State**
All functions accept explicit parameters; no module-level configuration.

## Module Organization

The `download` module is in `src/metainformant/core/io/download.py`. It's part of the `io` subpackage.

**Key exports**:
- `download_with_progress()` — Main high-level download function
- `monitor_subprocess_*()` — External tool monitoring functions
- `DownloadResult` — Dataclass with download outcome
- `DownloadHeartbeatState` — Mutable state written to heartbeat
- `DownloadHandler` — Protocol for custom handlers
- `HTTPDownloadHandler`, `FTPDownloadHandler`, `FileDownloadHandler` — Built-in handlers
- `get_download_handler()` — Factory for protocol-specific handler

**Helper functions** (internal):
- `_download_http()` — HTTP implementation
- `_download_ftp()` — FTP implementation
- `_download_file_url()` — file:// implementation
- `_heartbeat_path()` — Compute heartbeat file location
- `_atomic_write_json()` — Atomic JSON writes
- `_compute_speed_mbps()`, `_compute_eta_seconds()`, `_compute_progress_percent()` — Metrics

## Heartbeat Format

Heartbeat files are JSON with this schema:
```json
{
  "url": "https://example.com/file.fastq.gz",
  "destination": "/path/to/output/file.fastq.gz",
  "started_at": "2026-01-15T10:30:00Z",
  "last_update": "2026-01-15T10:35:00Z",
  "bytes_downloaded": 524288000,
  "total_bytes": 1073741824,
  "progress_percent": 48.82,
  "speed_mbps": 12.5,
  "eta_seconds": 42.3,
  "status": "downloading",  // or: starting, completed, failed
  "step": null,
  "progress": null,
  "errors": []
}
```

Fields:
- `status` transitions: `starting` → `downloading` → (`completed` | `failed`)
- `errors` array accumulates all retry failures
- `progress` field (optional) includes structured progress for subprocess monitors (type: `directory_size`, `file_count`, `sample_count`)

## API Reference

### Dataclasses

#### `DownloadResult`

Immutable result of a download attempt.

**Fields**:
- `url: str` — Source URL
- `dest_path: Path` — Final destination path
- `success: bool` — Whether download succeeded
- `bytes_downloaded: int = 0` — Total bytes written to disk
- `total_bytes: int | None = None` — Expected total from Content-Length (if known)
- `elapsed_seconds: float = 0.0` — Time taken for final attempt
- `resumed: bool = False` — Whether download resumed partial file
- `error: str | None = None` — Error message on failure
- `attempts: int = 1` — Number of attempts made

**Example**:
```python
result = download_with_progress(
    "https://example.com/large.bam",
    "output/large.bam"
)
if result.success:
    print(f"Downloaded {result.bytes_downloaded} bytes in {result.elapsed_seconds:.1f}s")
    print(f"Speed: {result.bytes_downloaded / result.elapsed_seconds / 1e6:.2f} MB/s")
else:
    print(f"Failed after {result.attempts} attempts: {result.error}")
```

#### `DownloadHeartbeatState`

Mutable state written to heartbeat file.

**Fields**:
- `url: str`
- `destination: str`
- `started_at: str` — ISO format UTC timestamp
- `last_update: str` — ISO format UTC timestamp
- `bytes_downloaded: int`
- `total_bytes: int | None`
- `progress_percent: float | None`
- `speed_mbps: float | None`
- `eta_seconds: float | None`
- `status: str` — One of: `starting`, `downloading`, `completed`, `failed`
- `step: str | None` — Optional step description (for subprocess monitoring)
- `progress: dict[str, Any] | None` — Structured progress object
- `errors: list[str]` — Accumulated error messages

**Methods**:
- `to_json() -> dict[str, Any]` — Serialize to JSON-compatible dict

### Main Function

#### `download_with_progress()`

Download a file with full observability.

**Signature**:
```text
def download_with_progress(
    url: str,
    dest_path: str | Path,
    *,
    protocol: str = "auto",
    heartbeat_interval: int = 5,
    show_progress: bool = True,
    resume: bool = True,
    max_retries: int = 3,
    retry_delay: float = 5.0,
    chunk_size: int = 1024 * 1024,
    timeout: int = 300,
) -> DownloadResult
```

**Parameters**:

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `url` | `str` | required | Source URL (http, https, ftp, file) |
| `dest_path` | `str \| Path` | required | Destination file path or directory |
| `protocol` | `str` | `"auto"` | Force protocol: `"http"`, `"https"`, `"ftp"`, `"file"` |
| `heartbeat_interval` | `int` | `5` | Seconds between heartbeat writes (0 disables) |
| `show_progress` | `bool` | `True` | Show tqdm progress bar if available |
| `resume` | `bool` | `True` | Attempt HTTP Range resume if partial file exists |
| `max_retries` | `int` | `3` | Number of retry attempts on failure |
| `retry_delay` | `float` | `5.0` | Base delay for exponential backoff (seconds) |
| `chunk_size` | `int` | `1048576` (1 MiB) | Read/write chunk size in bytes |
| `timeout` | `int` | `300` | Request timeout in seconds |

**Returns**: `DownloadResult` dataclass

**Raises**:
- `ValueError` — Unsupported protocol
- `ImportError` — If `requests` unavailable for HTTP (should be present)

**Behavior notes**:
- If `dest_path` is a directory, filename is inferred from URL
- Heartbeat file automatically created in `dest.parent/.downloads/`
- Progress bar uses `tqdm` if available; otherwise silent unless `show_progress=False` explicitly
- Retries use exponential backoff: delays = `retry_delay * (2 ^ attempt)` (5s, 10s, 20s, ...)
- After final failure, heartbeat file remains with `"status": "failed"` and error details

**Example**:
```python
result = download_with_progress(
    "https://ftp.ncbi.nlm.nih.gov/sra/sra-instant/read/ByRun/sra/SRR/SRR001/SRR001000/SRR001000_1.fastq.gz",
    "output/sra/",
    heartbeat_interval=10,
    show_progress=True,
    max_retries=5,
    resume=True
)
if not result.success:
    raise RuntimeError(f"Download failed: {result.error}")
```

### Subprocess Monitoring Functions

These functions are designed for external tools (e.g., `subprocess.Popen`) that perform their own downloads but need heartbeat integration.

#### `monitor_subprocess_directory_growth()`

Monitor a subprocess by watching byte growth in a directory.

**Signature**:
```text
def monitor_subprocess_directory_growth(
    *,
    process: Any,
    watch_dir: str | Path,
    heartbeat_path: str | Path,
    total_bytes: int | None = None,
    heartbeat_interval: int = 5,
    show_progress: bool = True,
    desc: str | None = None,
    errors: list[str] | None = None,
) -> tuple[int, int]
```

**Parameters**:
- `process` — `subprocess.Popen` object (must not have stdout/stderr piped to avoid deadlock)
- `watch_dir` — Directory to monitor for growing file sizes
- `heartbeat_path` — Where to write heartbeat JSON
- `total_bytes` — Expected total size (for progress bar and ETA); if None, progress indeterminate
- `heartbeat_interval` — Write heartbeat every N seconds
- `show_progress` — Show tqdm progress bar
- `desc` — Progress bar description
- `errors` — Mutable list to accumulate error strings (passed by caller)

**Returns**: `(returncode, bytes_written)` tuple when process terminates

**Example**:
```python
import subprocess
from metainformant.core.io.download import monitor_subprocess_directory_growth

proc = subprocess.Popen(
    ["prefetch", "--max-size", "100G", "SRR001000"],
    stdout=subprocess.PIPE,  # OK if you consume it
    stderr=subprocess.PIPE,
    cwd="/data/sra"
)

rc, bytes_done = monitor_subprocess_directory_growth(
    process=proc,
    watch_dir="/data/sra",
    heartbeat_path="output/.downloads/prefetch.heartbeat.json",
    total_bytes=50_000_000_000,  # 50 GB expected
    heartbeat_interval=30,
    show_progress=True,
    desc="Downloading SRA run"
)
if rc != 0:
    raise RuntimeError(f"prefetch failed with code {rc}")
```

#### `monitor_subprocess_file_count()`

Monitor by counting expected output files.

**Signature**:
```text
def monitor_subprocess_file_count(
    *,
    process: Any,
    watch_dir: str | Path,
    heartbeat_path: str | Path,
    expected_files: Iterable[str],
    heartbeat_interval: int = 5,
    show_progress: bool = True,
    desc: str | None = None,
    errors: list[str] | None = None,
) -> int
```

**Parameters**:
- `expected_files` — Iterable of relative file paths (from `watch_dir`) expected to be created

**Returns**: Process return code

**Example**:
```python
expected = [f"sample_{i}.fastq.gz" for i in range(1, 101)]
proc = subprocess.Popen(["parallel-fastq-dump", "--split-3", ...])

rc = monitor_subprocess_file_count(
    process=proc,
    watch_dir="output/fastq",
    heartbeat_path="output/.downloads/fastq_dump.heartbeat.json",
    expected_files=expected,
)
```

#### `monitor_subprocess_sample_progress()`

Monitor by counting files matching a glob pattern (for unknown exact filenames ahead of time).

**Signature**:
```text
def monitor_subprocess_sample_progress(
    *,
    process: Any,
    watch_dir: str | Path,
    heartbeat_path: str | Path,
    completion_glob: str,
    total_samples: int | None = None,
    heartbeat_interval: int = 5,
    show_progress: bool = True,
    desc: str | None = None,
    errors: list[str] | None = None,
) -> int
```

**Parameters**:
- `completion_glob` — Glob pattern (relative to `watch_dir`) matching completed files (e.g., `"*.fastq.gz"` or `"sample_*/finished.txt"`)
- `total_samples` — Expected final count (for percentage); if `None`, progress bar is indeterminate

**Returns**: Process return code

**Example**:
```python
proc = subprocess.Popen(
    ["snakemake", "--cores", "8", "all_quantifications"]
)

rc = monitor_subprocess_sample_progress(
    process=proc,
    watch_dir="output/quant",
    heartbeat_path="output/.downloads/quant.heartbeat.json",
    completion_glob="**/*.tsv",  # All TSV files count as completed samples
    total_samples=1000,
    desc="Quantifying samples"
)
```

### Handler Classes

#### `HTTPDownloadHandler`

Handler for HTTP/HTTPS URLs. Supports resume via Range headers.

**Methods**:
```text
get_file_size(url: str, *, timeout: int = 60) -> int | None
download(url: str, dest_path: str | Path, **kwargs) -> DownloadResult
```

**Example**:
```python
handler = HTTPDownloadHandler()
size = handler.get_file_size("https://example.com/file.bin", timeout=30)
result = handler.download(
    "https://example.com/file.bin",
    "local.bin",
    heartbeat_interval=5,
    show_progress=True,
)
```

#### `FTPDownloadHandler`

Handler for FTP URLs (no resume support).

**Example**:
```python
handler = FTPDownloadHandler()
result = handler.download(
    "ftp://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna_sm.toplevel.fa.gz",
    "reference/genome.fa.gz",
    show_progress=False,
)
```

#### `FileDownloadHandler`

Handler for local `file://` URLs (fast copy).

**Example**:
```python
handler = FileDownloadHandler()
result = handler.download(
    "file:///data/local/sample.fastq",
    "output/sample.fastq",
)
```

#### `get_download_handler(url: str) -> DownloadHandler`

Factory returning appropriate handler for URL scheme.

**Raises**: `ValueError` if no handler for scheme

**Example**:
```python
handler = get_download_handler("https://example.com/data.tar.gz")
size = handler.get_file_size()
result = handler.download(...)
```

## Usage Examples

### Basic File Download

```python
from metainformant.core.io.download import download_with_progress

result = download_with_progress(
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE1nnn/GSE1000/suppl/GSE1000_RAW.tar",
    "data/raw/",
    show_progress=True,
    max_retries=3,
)
if result.success:
    print(f"Downloaded {result.bytes_downloaded / 1e9:.2f} GB")
else:
    print(f"Failed: {result.error}")
```

### Download with Heartbeat Polling

```python
import json
import time
from pathlib import Path
from metainformant.core.io.download import download_with_progress

dest = Path("output/large_dataset.bam")
heartbeat_path = dest.parent / ".downloads" / f"{dest.name}.heartbeat.json"

# Start download in background thread or separate process
# Meanwhile, poll heartbeat:
def monitor_heartbeat(heartbeat_path: Path, timeout: int = 3600):
    start = time.time()
    while time.time() - start < timeout:
        if heartbeat_path.exists():
            with open(heartbeat_path) as fh:
                hb = json.load(fh)
            print(f"Status: {hb['status']}, "
                  f"Progress: {hb['progress_percent']:.1f}%, "
                  f"Speed: {hb['speed_mbps']:.2f} MB/s, "
                  f"ETA: {hb['eta_seconds'] / 60:.1f} min")
            if hb['status'] in ('completed', 'failed'):
                break
        time.sleep(5)

# In real usage, this monitoring runs in separate thread/process
```

### Batch Download with Parallelism

```python
from metainformant.core.execution import parallel
from metainformant.core.io.download import download_with_progress

def download_one(url_dest_pair):
    url, dest = url_dest_pair
    return download_with_progress(
        url,
        dest,
        show_progress=False,  # Don't clutter output with many bars
        heartbeat_interval=0,  # Disable heartbeat for batch mode
        max_retries=2,
    )

# List of (url, destination) pairs
downloads = [
    ("https://example.com/data1.bam", "output/sample1.bam"),
    ("https://example.com/data2.bam", "output/sample2.bam"),
    # ... 100s more
]

# Use thread_map for I/O-bound parallelism
results = parallel.thread_map(
    download_one,
    downloads,
    max_workers=8,
    timeout=600,
)

# Check results
successes = [r for r in results if r.success]
failures = [r for r in results if not r.success]
print(f"Downloaded {len(successes)}/{len(results)} files")
for fail in failures:
    print(f"  FAILED: {fail.dest_path} – {fail.error}")
```

### Download with Custom Retry Policy

```python
from metainformant.core.io.download import download_with_progress
import time

result = download_with_progress(
    "https://unreliable.server.edu/dataset.tar.gz",
    "output/dataset.tar.gz",
    max_retries=5,
    retry_delay=30.0,  # Wait 30s, 60s, 120s, 240s, 480s
    chunk_size=512 * 1024,  # 512 KiB chunks
    timeout=60,  # 60-second request timeout
)

if not result.success:
    print(f"All retries exhausted. Last error: {result.error}")
    print(f"Partial file size: {result.bytes_downloaded} bytes")
    # Decide whether to keep partial or delete
    if result.bytes_downloaded > 0:
        print("Partial download preserved—resume may be possible")
```

### Handling Non-Resumable Protocols

```python
# FTP does not support resume—be cautious
result = download_with_progress(
    "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR001/001/SRR001000/SRR001000_1.fastq.gz",
    "output/srr.fastq.gz",
    resume=False,  # Explicit (FTP handler ignores anyway)
    max_retries=3,
)

if result.success:
    print("FTP download complete")
```

### Local File Copy

```python
# Useful for testing or post-processing within pipeline
result = download_with_progress(
    "file:///data/interim/processed.bam",
    "output/final/processed.bam",
    show_progress=False,
)
```

### Subprocess Monitoring: External Download Tool

```python
import subprocess
from metainformant.core.io.download import monitor_subprocess_directory_growth

# Use fastq-dump or sra-tools which handles SRA complexities
proc = subprocess.Popen(
    ["fastq-dump", "--split-3", "--gzip", "SRR001000"],
    cwd="/data/sra",
    stdout=subprocess.PIPE,  # Capture to avoid blocking
    stderr=subprocess.PIPE,
)

# Monitor size of output directory
rc, bytes_written = monitor_subprocess_directory_growth(
    process=proc,
    watch_dir="/data/sra",
    heartbeat_path="output/.downloads/fastq-dump.heartbeat.json",
    total_bytes=2_500_000_000,  # 2.5 GB expected
    heartbeat_interval=10,
    show_progress=True,
    desc="fastq-dump SRR001000",
)

if rc == 0:
    print(f"Downloaded {bytes_written} bytes successfully")
else:
    print(f"fastq-dump failed with code {rc}")
```

### Subprocess Monitoring: Counting Output Files

```python
import subprocess
from metainformant.core.io.download import monitor_subprocess_file_count

# Tool that produces 96 FASTQ files (48 samples × 2 read directions)
expected = []
for i in range(1, 49):
    expected.append(f"sample_{i}_R1.fastq.gz")
    expected.append(f"sample_{i}_R2.fastq.gz")

proc = subprocess.Popen(["parallel-fastq-dump", "--split-3", ...])

rc = monitor_subprocess_file_count(
    process=proc,
    watch_dir="output/fastq",
    heartbeat_path="output/.downloads/parallel_fastq_dump.heartbeat.json",
    expected_files=expected,
    heartbeat_interval=15,
    desc="Downloading 48 samples",
)
```

### Custom Download Handler

```python-snippet
from metainformant.core.io.download import DownloadHandler, DownloadResult
import requests

class S3DownloadHandler:
    """Custom handler for S3 URLs (s3://bucket/key)."""

    def get_file_size(self, url: str, *, timeout: int = 60) -> int | None:
        import boto3
        parsed = url.replace("s3://", "").split("/", 1)
        bucket, key = parsed[0], parsed[1]
        s3 = boto3.client("s3")
        resp = s3.head_object(Bucket=bucket, Key=key)
        return resp["ContentLength"]

    def download(
        self,
        url: str,
        dest_path: str | Path,
        *,
        heartbeat_interval: int = 5,
        show_progress: bool = True,
        **kwargs
    ) -> DownloadResult:
        import boto3
        from botocore.exceptions import ClientError

        dest = Path(dest_path)
        dest.parent.mkdir(parents=True, exist_ok=True)

        parsed = url.replace("s3://", "").split("/", 1)
        bucket, key = parsed[0], parsed[1]

        s3 = boto3.client("s3")
        start = time.time()
        bytes_written = 0

        try:
            with open(dest, "wb") as fh:
                s3.download_fileobj(bucket, key, fh)
            return DownloadResult(
                url=url,
                dest_path=dest,
                success=True,
                bytes_downloaded=dest.stat().st_size,
                elapsed_seconds=time.time() - start,
            )
        except ClientError as e:
            return DownloadResult(
                url=url,
                dest_path=dest,
                success=False,
                error=str(e),
            )

# Register custom handler by wrapping download_with_progress
def download_s3(url, dest, **kwargs):
    handler = S3DownloadHandler()
    return handler.download(url, dest, **kwargs)

result = download_s3("s3://my-bucket/genome.fa.gz", "reference/genome.fa.gz")
```

### Reading Heartbeat Files

```python
import json
from pathlib import Path
import time

def read_heartbeat(heartbeat_path: str | Path) -> dict | None:
    path = Path(heartbeat_path)
    if not path.exists():
        return None
    try:
        with open(path) as fh:
            return json.load(fh)
    except json.JSONDecodeError:
        return None

def monitor_multiple_downloads(heartbeat_dir: Path):
    """Monitor all active downloads in a directory."""
    while True:
        hb_files = list(heartbeat_dir.glob("*.heartbeat.json"))
        for hb_path in hb_files:
            hb = read_heartbeat(hb_path)
            if hb and hb["status"] == "downloading":
                print(f"{hb['destination']}: {hb['progress_percent']:.1f}% "
                      f"at {hb['speed_mbps']:.2f} MB/s")
        time.sleep(10)

# Usage
monitor_multiple_downloads(Path("output/.downloads"))
```

## Error Handling

### Network Failures

**Symptom**: `requests.exceptions.ConnectionError`, `Timeout`, all retries exhausted → `DownloadResult.success=False`.

**Common causes**:
1. **Server down** — Wait and retry later; check service status
2. **DNS resolution failure** — Check network, DNS server
3. **Connection timeout** — Server slow, increase `timeout` parameter
4. **Rate limiting (HTTP 429)** — Increase `retry_delay`, add jitter

**Solutions**:
```python
result = download_with_progress(
    url,
    dest,
    max_retries=5,
    retry_delay=60.0,  # Longer delays for rate-limited servers
    timeout=120,  # Longer timeout
)
```

**Backoff jitter** to avoid synchronized retries (if you're downloading many files from same server):
```python
import random

def download_with_jitter(url, dest):
    # Add randomness to retry_delay to avoid thundering herd
    base_delay = 30.0
    jitter = random.uniform(0.8, 1.2)  # ±20%
    return download_with_progress(
        url, dest,
        max_retries=5,
        retry_delay=base_delay * jitter,
    )
```

### HTTP Errors

**Symptom**: `requests.exceptions.HTTPError: 404 Client Error` or other status codes.

**Handling**:
- `404 Not Found` — URL incorrect or file removed; verify URL
- `403 Forbidden` — Permission denied; check credentials or referer
- `429 Too Many Requests` — Rate limited; retry with longer delays
- `500/502/503` — Server error; retry with backoff

```python
result = download_with_progress(url, dest, max_retries=3)
if not result.success and "404" in result.error:
    print(f"File not found: {url}")
elif not result.success:
    # Retry later with more attempts
    pass
```

### Disk Space Exhaustion

**Symptom**: `OSError: [Errno 28] No space left on device` during download.

**Pre-check**:
```python
import shutil
from metainformant.core.io.download import _http_head_size

def safe_download(url, dest, reserve_gb=10):
    """Download only if enough disk space available."""
    total_size = _http_head_size(url, timeout=30) or 0
    free_gb = shutil.disk_usage(Path(dest.parent)).free / (1024**3)
    needed_gb = total_size / (1024**3) + reserve_gb

    if free_gb < needed_gb:
        raise OSError(
            f"Insufficient disk space: {free_gb:.1f} GB free, "
            f"need {needed_gb:.1f} GB"
        )

    return download_with_progress(url, dest)
```

### Corrupted Partial Files

**Symptom**: Download "completes" but file is corrupted (checksum mismatch).

**Fix**: Enable checksum validation:
```python
import hashlib

def download_and_verify(url, dest, expected_md5):
    result = download_with_progress(url, dest, resume=True)
    if not result.success:
        raise RuntimeError(result.error)

    # Verify checksum
    md5 = hashlib.md5()
    with open(dest, "rb") as fh:
        while chunk := fh.read(1024 * 1024):
            md5.update(chunk)

    if md5.hexdigest() != expected_md5:
        dest.unlink(missing_ok=True)  # Delete corrupted file
        raise ValueError(
            f"Checksum mismatch: got {md5.hexdigest()}, "
            f"expected {expected_md5}"
        )
    return result
```

### Stalled Downloads

**Symptom**: Progress bar frozen, no bytes written for extended period.

**Causes**:
- Server not sending data (connection hung)
- Extremely slow network (may be valid)
- Progress display bug

**Diagnosis**:
```python
# Check heartbeat file for last update timestamp
import json, time

hb = json.load(open("file.heartbeat.json"))
last_update = time.mktime(time.strptime(hb["last_update"], "%Y-%m-%dT%H:%M:%SZ"))
if time.time() - last_update > 300:  # 5 minutes
    print("Download stalled—terminating")
```

### Handle Unsupported Protocols

**Symptom**: `ValueError: Unsupported protocol/scheme: s3`.

**Fix**: Implement custom `DownloadHandler` (see Custom Handler example above) or use specialized library (boto3 for S3, gsutil for GCS).

### Broken Pipe with Subprocess Monitoring

**Symptom**: Subprocess hangs when using `monitor_subprocess_directory_growth`.

**Cause**: You piped `stdout` or `stderr` of subprocess; child blocks on full pipe buffer.

**Fix**: Either consume subprocess output or redirect to file:
```python
# Option 1: Consume pipes (requires threading or async)
proc = subprocess.Popen(
    cmd,
    stdout=subprocess.PIPE,
    stderr=subprocess.PIPE,
)
# Must read stdout/stderr in separate threads to avoid deadlock

# Option 2: Redirect to files (recommended for monitoring)
with open("stdout.log", "w") as out, open("stderr.log", "w") as err:
    proc = subprocess.Popen(cmd, stdout=out, stderr=err)
    # monitor_subprocess...() now works fine
```

## Performance Considerations

### Chunk Size Tuning

Default chunk size is 1 MiB (1,048,576 bytes). Adjust based on storage medium:

| Storage | Recommended chunk size | Rationale |
|---------|----------------------|-----------|
| HDD | 1–4 MiB | Larger chunks reduce seek overhead |
| SSD/NVMe | 512 KiB–1 MiB | Smaller chunks keep I/O pipelines full |
| Network | 1–4 MiB | Balance latency vs throughput |
| Memory filesystem | 64–256 KiB | Minimal I/O overhead; small chunks fine |

**Benchmarking**:
```python
import time

def benchmark_chunk(url, dest, chunk_size):
    start = time.time()
    download_with_progress(
        url, dest,
        chunk_size=chunk_size,
        show_progress=False,
    )
    elapsed = time.time() - start
    size = Path(dest).stat().st_size
    print(f"Chunk {chunk_size}: {size / elapsed / 1e6:.2f} MB/s")

for cs in [64*1024, 256*1024, 1024*1024, 4*1024*1024]:
    benchmark_chunk(url, f"test_{cs}.bin", cs)
```

### Concurrent Download Limits

Too many parallel downloads can:
- Overwhelm local disk I/O (especially HDD)
- Trigger server-side rate limiting
- Exhaust file descriptor limits

**Rule**: Total parallel downloads ≤ CPU cores × 2, and ≤ 8–12 for most systems. Monitor system load (`iostat`, `iotop`) and network utilization.

### HTTP Keep-Alive

Requests session reuse is handled internally by `requests` library. For many small files to same host, connection pooling reduces handshake overhead. Custom handlers can maintain session:

```python
import requests
from functools import lru_cache

@lru_cache(maxsize=1)
def get_requests_session():
    session = requests.Session()
    adapter = requests.adapters.HTTPAdapter(
        pool_connections=10,
        pool_maxsize=20,
        max_retries=3
    )
    session.mount("http://", adapter)
    session.mount("https://", adapter)
    return session
```

### Resume Efficiency

Resuming is most efficient when:
- Server supports HTTP `Range` header (most do)
- Partial file is contiguous (no fragmentation)
- Chunk size aligns with server's preferred range granularity

Check for resume support:
```python
resp = requests.head(url)
accepts_ranges = resp.headers.get("Accept-Ranges", "none")
supports_resume = "bytes" in accepts_ranges.lower()
```

### Heartbeat Overhead

Heartbeat writes are small (<2KB) and infrequent (default every 5s). The atomic write (tmp + rename) is very fast on modern filesystems (<1ms). For 10,000 downloads, heartbeat overhead is <5 seconds total.

### Directory Monitoring Performance

`monitor_subprocess_directory_growth()` uses a cached directory size calculator (TTL 1s) to avoid repeatedly walking the entire tree. For directories with millions of files, this caching is critical:

```python
# Without caching, 1-second sleep + rglob over 1M files = expensive
# With 1s TTL cache, same result returned from memory until TTL expires
```

## Security Notes

### URL Validation

Untrusted URLs should be validated:
```python
from urllib.parse import urlparse

def is_safe_url(url: str) -> bool:
    parsed = urlparse(url)
    # Only allow HTTP(S) and trusted protocols
    if parsed.scheme not in ("http", "https", "ftp"):
        return False
    # Block private IP ranges if downloading from internet
    # (use ipaddress module to check)
    return True
```

### Path Traversal Protection

`download_with_progress()` does not sanitize destination paths—you must ensure dest_path is within allowed output directories:

```python
from metainformant.core.io import paths

dest = "output/" + user_provided_filename
if not paths.is_safe_path(dest):
    raise ValueError(f"Unsafe path: {dest}")

result = download_with_progress(url, dest)
```

### Temporary File Security

Download creates temp files with default umask (usually 0o022) → readable by all users on shared system. For sensitive data:
```python
import os
import tempfile

# Create temp file with restricted permissions
temp_fd, temp_path = tempfile.mkstemp(
    dir=dest.parent,
    prefix=".dl_",
    suffix=".part"
)
os.close(temp_fd)
os.chmod(temp_path, 0o600)  # Owner read/write only
```

### Content-Type Validation

For security-critical downloads, validate Content-Type matches expected:
```python
import requests

resp = requests.head(url, allow_redirects=True)
content_type = resp.headers.get("Content-Type", "")
if "application/gzip" not in content_type and "application/x-gzip" not in content_type:
    raise ValueError(f"Unexpected content type: {content_type}")
```

## Integration Examples

### With Parallel Processing

```python
from metainformant.core.execution import parallel
from metainformant.core.io.download import download_with_progress

urls = [
    ("https://example.com/a.fastq.gz", "output/a.fastq.gz"),
    ("https://example.com/b.fastq.gz", "output/b.fastq.gz"),
    # ... many more
]

def download_task(pair):
    url, dest = pair
    return download_with_progress(
        url, dest,
        show_progress=False,  # Manage progress at batch level
        heartbeat_interval=0,  # Disable per-file heartbeats
    )

results = parallel.thread_map(
    download_task,
    urls,
    max_workers=8,
    timeout=900,
)

successful = [r for r in results if r.success]
print(f"Downloaded {len(successful)} / {len(urls)} files")
```

### With I/O Module

`download_with_progress()` is itself part of `metainformant.core.io` package. Use with other I/O utilities:

```python
from metainformant.core.io import download, io, cache

def download_and_cache(url, cache_dir):
    # Check cache first
    cache_key = cache.key_from_url(url)
    cached = cache.get(cache_key, cache_dir)
    if cached:
        return cached

    # Download to temp then move to cache atomically
    dest = Path(cache_dir) / cache_key
    result = download.download_with_progress(
        url,
        dest,
        show_progress=True,
    )
    if result.success:
        cache.set(cache_key, dest, cache_dir)
    return dest
```

### With Paths Module

```python
from metainformant.core.io import download
from metainformant.core.io import paths

def safe_download(url, user_proposed_name, output_root="output"):
    # Sanitize and resolve destination
    safe_name = paths.sanitize_filename(user_proposed_name)
    dest = paths.expand_and_resolve(Path(output_root) / "downloads" / safe_name)

    # Ensure directory exists
    paths.prepare_file_path(dest)

    # Verify destination is within output_root
    if not paths.is_within(dest, paths.expand_and_resolve(output_root)):
        raise ValueError(f"Destination outside output root: {dest}")

    return download.download_with_progress(url, dest)
```

### With Configuration

```python
from metainformant.core.io import download
from metainformant.core.utils import config

cfg = config.load_mapping_from_file("config/downloads.yaml")

for source in cfg["sources"]:
    result = download.download_with_progress(
        url=source["url"],
        dest_path=source["destination"],
        timeout=source.get("timeout", 300),
        max_retries=source.get("max_retries", 3),
    )
    if not result.success:
        print(f"FAILED: {source['name']}: {result.error}")
```

### In a Snakemake Rule

```text
# In a Snakemake workflow, use heartbeat for resume across runs
rule download_genome:
    output:
        genome="reference/genome.fa.gz",
        heartbeat="reference/.downloads/genome.fa.gz.heartbeat.json"
    params:
        url="https://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna_sm.toplevel.fa.gz"
    shell:
        """
        python -c "
from metainformant.core.io.download import download_with_progress
result = download_with_progress(
    '{params.url}',
    '{output.genome}',
    heartbeat_interval=10,
    show_progress=False,
)
assert result.success, result.error
"
        """
```

## Testing Guidelines

### Test with Temporary Files

```python
from pathlib import Path
from metainformant.core.io.download import download_with_progress

def test_download_to_temp(tmp_path: Path):
    dest = tmp_path / "download.bin"
    result = download_with_progress(
        "file:///etc/hosts",  # Small local file for test
        dest,
        show_progress=False,
        heartbeat_interval=0,
    )
    assert result.success
    assert dest.exists()
    assert dest.stat().st_size > 0
```

### Test Heartbeat Creation

```python
def test_heartbeat_created(tmp_path: Path):
    dest = tmp_path / "file.dat"
    hb_dir = tmp_path / ".downloads"
    hb_file = hb_dir / "file.dat.heartbeat.json"

    result = download_with_progress(
        "file:///etc/hosts",
        dest,
        heartbeat_interval=1,
        show_progress=False,
    )

    assert hb_file.exists()
    hb = json.loads(hb_file.read_text())
    assert hb["status"] == "completed"
    assert "bytes_downloaded" in hb
    assert "speed_mbps" in hb
```

### Test Retry Logic

Retry behavior should be tested with a small local HTTP server that returns
deterministic failure responses before serving the final file. That keeps the
download stack real while avoiding external network dependence.

### Test with Real HTTP Server

```python
import http.server
import threading

def test_download_from_local_server(tmp_path):
    # Start local HTTP server
    class Handler(http.server.BaseHTTPRequestHandler):
        def do_GET(self):
            self.send_response(200)
            self.send_header("Content-Length", "1024")
            self.end_headers()
            self.wfile.write(b"x" * 1024)

    server = http.server.HTTPServer(("127.0.0.1", 0), Handler)
    port = server.server_address[1]
    thread = threading.Thread(target=server.serve_forever)
    thread.daemon = True
    thread.start()

    try:
        result = download_with_progress(
            f"http://127.0.0.1:{port}/file.bin",
            tmp_path / "file.bin",
            show_progress=False,
        )
        assert result.success
        assert result.bytes_downloaded == 1024
    finally:
        server.shutdown()
```

### Monitor Subprocess Tests

```python
def test_monitor_directory_growth_subprocess(tmp_path):
    import subprocess, sys, time

    out_dir = tmp_path / "out"
    out_dir.mkdir()
    hb = tmp_path / "hb.json"

    # Child process writes 1MB every 0.5s for 3 iterations
    code = f"""
import time, pathlib
p = pathlib.Path(r'{out_dir}') / 'data.bin'
for _ in range(3):
    p.write_bytes(b'x' * 1024*1024)
    pathlib.Path(p).touch()
    time.sleep(0.5)
"""
    proc = subprocess.Popen([sys.executable, "-c", code])

    from metainformant.core.io.download import monitor_subprocess_directory_growth
    rc, bytes_written = monitor_subprocess_directory_growth(
        process=proc,
        watch_dir=out_dir,
        heartbeat_path=hb,
        total_bytes=3*1024*1024,
        heartbeat_interval=1,
        show_progress=False,
    )

    assert rc == 0
    assert bytes_written >= 3*1024*1024
    assert hb.exists()
```

### Skip Slow Network Tests

```python
import pytest

@pytest.mark.network
def test_download_real_sra_file():
    """Integration test requiring network and ~100MB download."""
    result = download_with_progress(
        "https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR001/001/SRR001000/SRR001000_1.fastq.gz",
        tmp_path / "sra.fastq.gz",
        show_progress=False,
        timeout=120,
    )
    assert result.success
    assert result.bytes_downloaded > 100_000_000
```

## Dependencies

### Required
- `requests` — HTTP library (already used elsewhere in METAINFORMANT)
- Standard library: `json`, `os`, `time`, `pathlib`, `urllib.parse`, `dataclasses`, `typing`

### Optional
- `tqdm` — Progress bars (imported via `metainformant.core.progress`)

### External Tools
None—pure Python implementation.

## Troubleshooting

| Symptom | Likely Cause | Fix |
|---------|--------------|-----|
| Heartbeat file not updating | `heartbeat_interval=0` or process hung | Check `download_with_progress()` args; verify process still running |
| Progress bar flickering | tqdm + multi-thread collision | Use `show_progress=False` in parallel downloads |
| `Pool exhausted` errors | Too many parallel downloads | Reduce `max_workers` in `parallel.thread_map()` |
| Partial file never resumes | Server lacks Range support | Set `resume=False` or pre-delete partial |
| `[Errno 28]` No space | Disk full | Check `df -h`; clean up; use `--max-size` filter |
| Download hangs indefinitely | Server not sending data, timeout too high | Increase timeout or set moderately (e.g., 60s) and let retry logic handle |

## Further Reading

- HTTP Range Requests: RFC 7233
- Exponential backoff: https://aws.amazon.com/blogs/architecture/exponential-backoff-and-jitter/
- tqdm documentation: https://tqdm.github.io/
