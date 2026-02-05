"""Amalgkit CLI integration for RNA-seq workflow orchestration.

This module provides Python wrappers around the amalgkit command-line interface
for comprehensive RNA-seq analysis workflows.

PATH RESOLUTION NOTES:
- getfastq: Creates FASTQ files in {out_dir}/getfastq/{sample_id}/ (automatically creates getfastq subdirectory)
- integrate: Expects fastq_dir to point to {fastq_dir}/getfastq/ (the getfastq subdirectory)
- quant: Should use out_dir = work_dir so it can find getfastq output in {out_dir}/getfastq/
- merge: Looks for abundance files in {out_dir}/quant/{sample_id}/{sample_id}_abundance.tsv

See docs/rna/amalgkit/PATH_RESOLUTION.md for complete path resolution documentation.
"""

from __future__ import annotations

import concurrent.futures
import csv
import math
import os
import shutil
import subprocess
import sys
import tempfile
import time
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, Union

from metainformant.core import logging
from metainformant.core.io.download import (
    monitor_subprocess_directory_growth,
    monitor_subprocess_file_count,
    monitor_subprocess_sample_progress,
)

logger = logging.get_logger(__name__)


class AmalgkitParams:
    """Parameters for amalgkit commands."""

    def __init__(
        self, work_dir: Union[str, Path], threads: int = 8, species_list: Optional[List[str]] = None, **kwargs
    ):
        """Initialize amalgkit parameters.

        Args:
            work_dir: Working directory for amalgkit operations. Output will be structured here.
            threads: Number of CPU threads to use for parallel operations. Defaults to 8.
            species_list: List of scientific names (e.g. ['Pogonomyrmex barbatus']) to process.
            **kwargs: Additional parameters passed directly to the amalgkit CLI.
                     Common keys: 'redo', 'pfd', 'aws', 'jobs', etc.
        """
        self.work_dir = Path(work_dir)
        self.threads = threads
        self.species_list = species_list or []
        self.extra_params = kwargs

    def to_dict(self) -> Dict[str, Any]:
        """Convert parameters to dictionary."""
        return {
            "work_dir": str(self.work_dir),
            "threads": self.threads,
            "species_list": self.species_list,
            **self.extra_params,
        }


def build_cli_args(
    params: AmalgkitParams | Dict[str, Any] | None, subcommand: str | None = None, *, for_cli: bool = False
) -> List[str]:
    """Build command-line arguments for amalgkit.

    Converts an AmalgkitParams object or dictionary into a list of command-line
    arguments suitable for passing to subprocess.run() or similar.

    Args:
        params: Amalgkit parameters. Can be:
            - AmalgkitParams: Structured parameter object with work_dir, threads, species_list
            - Dict[str, Any]: Dictionary with parameter key-value pairs
            - None: Returns empty list
        subcommand: The amalgkit subcommand being run (e.g., 'getfastq', 'quant').
            Used to filter out unsupported arguments for specific subcommands.
        for_cli: Whether to format for direct CLI usage. Currently unused but
            reserved for future formatting options.

    Returns:
        List of command-line argument strings, e.g. ['--out_dir', '/path', '--threads', '8']

    Example:
        >>> params = {'out_dir': '/work', 'threads': 4, 'redo': True}
        >>> build_cli_args(params, subcommand='getfastq')
        ['--out_dir', '/work', '--threads', '4', '--redo', 'yes']
    """
    args = []

    if params is None:
        return args

    # Handle dict parameters
    if isinstance(params, dict):
        # amalgkit uses --out_dir for the main working directory in all subcommands
        if "out_dir" in params:
            args.extend(["--out_dir", str(params["out_dir"])])
        elif "work_dir" in params:
            args.extend(["--out_dir", str(params["work_dir"])])

        if "threads" in params:
            args.extend(["--threads", str(params["threads"])])

        # skip_keys should check both underscores and hyphens
        skip_keys = {
            "work_dir",
            "out_dir",
            "threads",
            "extra_params",
            "genome",
            "filters",
            "taxon_id",
            "auto_install_amalgkit",
            "log_dir",
            "num_download_workers",
        }

        # Flags that expect explicit yes/no instead of just being present/absent
        yes_no_flags = {
            "redo",
            "pfd",
            "resolve_names",
            "mark_redundant_biosamples",
            "aws",
            "ncbi",
            "gcp",
            "fastp",
            "remove_sra",
            "remove_tmp",
            "pfd_print",
            "fastp_print",
            "build_index",
            "overwrite",
        }

        for key, value in params.items():
            if key in skip_keys or value is None:
                continue

            # Skip arguments not supported by specific subcommands
            if subcommand in ("merge", "sanity", "select", "curate") and key in ("out", "threads", "priority", "redo"):
                continue

            # Amalgkit 0.12.20+ generally uses underscores for its CLI flags.
            # We'll normalize to underscores to ensure compatibility.
            cli_key = key.replace("-", "_")

            if isinstance(value, bool):
                if cli_key in yes_no_flags:
                    args.extend([f"--{cli_key}", "yes" if value else "no"])
                elif value:
                    args.append(f"--{cli_key}")
            elif isinstance(value, (int, float)):
                if key == "threads" and subcommand in ("merge", "metadata", "config", "sanity", "select", "curate"):
                    continue
                args.extend([f"--{cli_key}", str(value)])
            elif isinstance(value, list):
                # Handle lists (multiple arguments)
                if key == "species-list" or key == "species_list" or key == "species":
                    if subcommand not in (
                        "integrate",
                        "config",
                        "help",
                        "merge",
                        "metadata",
                        "getfastq",
                        "quant",
                        "sanity",
                        "select",
                    ):
                        for val in value:
                            args.extend(["--species", str(val)])
                else:
                    for val in value:
                        args.extend([f"--{cli_key}", str(val)])
            elif key == "species-list" or key == "species_list" or key == "species":
                if subcommand not in ("integrate", "config", "help", "merge", "metadata", "getfastq", "quant"):
                    args.extend(["--species", str(value)])
            else:
                args.extend([f"--{cli_key}", str(value)])
        return args

    # Handle AmalgkitParams object
    args.extend(["--out_dir", str(params.work_dir)])
    if subcommand not in ("merge", "metadata", "config"):
        args.extend(["--threads", str(params.threads)])

    if params.species_list and subcommand not in ("integrate", "config", "merge", "metadata", "getfastq", "quant"):
        for species in params.species_list:
            args.extend(["--species", species])

    # Add extra parameters
    for key, value in params.extra_params.items():
        if isinstance(value, bool):
            if value:
                # Keep underscores for amalgkit compatibility
                args.extend([f"--{key}", "yes"])
        elif isinstance(value, (int, float)):
            args.extend([f"--{key}", str(value)])
        elif isinstance(value, str):
            args.extend([f"--{key}", value])
        elif isinstance(value, list):
            for item in value:
                args.extend([f"--{key}", str(item)])

    return args


def build_amalgkit_command(subcommand: str, params: AmalgkitParams | Dict[str, Any] | None = None) -> List[str]:
    """Build full amalgkit command.

    Args:
        subcommand: Amalgkit subcommand
        params: Parameters for the command (AmalgkitParams object or dict)

    Returns:
        Complete command as list of strings
    """
    command = ["amalgkit", subcommand]

    if params:
        command.extend(build_cli_args(params, subcommand=subcommand))

    return command


def check_cli_available() -> Tuple[bool, str]:
    """Check if amalgkit CLI is available.

    Returns:
        Tuple of (available, message)
    """
    try:
        result = subprocess.run(["amalgkit", "--help"], capture_output=True, text=True, timeout=10)
        if result.returncode == 0:
            return True, "amalgkit CLI is available"
        else:
            return False, f"amalgkit CLI returned error: {result.stderr}"
    except FileNotFoundError:
        return False, "amalgkit CLI not found in PATH"
    except subprocess.TimeoutExpired:
        return False, "amalgkit CLI check timed out"
    except (OSError, PermissionError) as e:
        return False, f"Error checking amalgkit CLI: {e}"


def ensure_cli_available(*, auto_install: bool = False) -> Tuple[bool, str, Dict | None]:
    """Ensure amalgkit CLI is available, optionally installing it.

    Args:
        auto_install: Whether to attempt automatic installation

    Returns:
        Tuple of (success, message, version_info)
    """
    available, message = check_cli_available()

    if available:
        # Try to get version info
        try:
            result = subprocess.run(["amalgkit", "--version"], capture_output=True, text=True, timeout=5)
            version_info = {"version": result.stdout.strip()} if result.returncode == 0 else None
        except (subprocess.SubprocessError, OSError, TimeoutError):
            version_info = None

        return True, message, version_info

    if not auto_install:
        return False, message, None

    # Attempt automatic installation using UV package manager
    logger.info("Attempting to install amalgkit via UV package manager...")
    try:
        # Use UV package manager (per METAINFORMANT policy)
        # Note: subprocess is imported at module level
        install_result = subprocess.run(
            ["uv", "pip", "install", "amalgkit"], capture_output=True, text=True, timeout=300
        )

        if install_result.returncode == 0:
            logger.info("Successfully installed amalgkit via UV")
            # Verify installation after install
            return ensure_cli_available(auto_install=False)
        else:
            error_msg = install_result.stderr if install_result.stderr else "Installation failed with no error message"
            logger.error(f"UV installation of amalgkit failed: {error_msg}")
            return False, f"Installation failed: {error_msg}", None
    except subprocess.TimeoutExpired:
        return False, "Installation timed out (timeout: 300s)", None
    except FileNotFoundError:
        return False, "UV package manager not found. Install UV to use automatic installation.", None
    except Exception as e:
        logger.error(f"Unexpected error during amalgkit installation: {e}")
        return False, f"Installation failed with unexpected error: {e}", None


def run_amalgkit(
    subcommand: str, params: AmalgkitParams | Dict[str, Any] | None = None, **kwargs: Any
) -> subprocess.CompletedProcess[str]:
    """Run amalgkit command.

    Args:
        subcommand: Amalgkit subcommand
        params: Parameters for the command
        **kwargs: Additional arguments passed to subprocess.run

    Returns:
        CompletedProcess instance
    """
    command = build_amalgkit_command(subcommand, params)

    logger.info(f"Running amalgkit command: {' '.join(command)}")

    # Extract custom keys to avoid passing them to subprocess.run/Popen twice
    cwd = kwargs.pop("cwd", None) or kwargs.pop("work_dir", None)
    log_dir = kwargs.pop("log_dir", None)
    step_name = kwargs.pop("step_name", None)
    monitor = kwargs.pop("monitor", None)
    heartbeat_interval = int(kwargs.pop("heartbeat_interval", 5))
    show_progress = bool(kwargs.pop("show_progress", True))

    if log_dir:
        import time

        ts = int(time.time())
        log_dir_path = Path(log_dir)
        log_dir_path.mkdir(parents=True, exist_ok=True)
        stdout_log = log_dir_path / f"{ts}.{subcommand}.stdout.log"
        stderr_log = log_dir_path / f"{ts}.{subcommand}.stderr.log"
        stdout_fh = open(stdout_log, "w")
        stderr_fh = open(stderr_log, "w")

    try:
        # Default behavior stays the same unless monitoring is enabled.
        # Enable monitoring for long-running steps by default.
        if monitor is None:
            monitor = subcommand in {"metadata", "integrate", "getfastq", "quant", "merge", "cstmm", "curate", "csca"}

        if monitor and isinstance(params, dict):
            out_dir_val = params.get("out_dir")
            metadata_val = params.get("metadata")
            if out_dir_val:
                out_dir = Path(str(out_dir_val))
                hb_path = out_dir / ".downloads" / f"amalgkit-{subcommand}.heartbeat.json"

                watch_dir = out_dir
                if subcommand == "getfastq":
                    watch_dir = out_dir / "getfastq"
                elif subcommand == "quant":
                    watch_dir = out_dir / "quant"
                elif subcommand == "metadata":
                    watch_dir = out_dir / "metadata"
                elif subcommand == "merge":
                    watch_dir = out_dir / "merge"
                elif subcommand == "cstmm":
                    watch_dir = out_dir / "cstmm"
                elif subcommand == "curate":
                    watch_dir = out_dir / "curate"
                elif subcommand == "csca":
                    watch_dir = out_dir / "csca"

                # Best-effort total estimations
                meta_path: Path | None = None
                if metadata_val:
                    meta_path = Path(str(metadata_val))
                else:
                    meta_path = _infer_metadata_path_from_out_dir(out_dir)

                total_bytes: int | None = None
                total_runs: int | None = None
                if meta_path is not None:
                    total_runs = _estimate_total_runs_from_metadata(meta_path)
                    if subcommand == "getfastq":
                        total_bytes = _estimate_total_bytes_from_metadata(meta_path)

                env = os.environ.copy()
                # Check tqdm compatibility
                try:
                    import tqdm

                    if tuple(map(int, tqdm.__version__.split("."))) < (4, 60, 0):
                        env["AMALGKIT_PROGRESS"] = "false"
                except (ImportError, ValueError, AttributeError):
                    pass

                if subcommand == "getfastq":
                    repo_root = Path(__file__).resolve().parent.parent.parent.parent
                    tmp_dir = repo_root / ".tmp" / "fasterq-dump"
                    tmp_dir.mkdir(parents=True, exist_ok=True)
                    env["TMPDIR"] = str(tmp_dir)
                    vdb_cache = repo_root / ".tmp" / "vdb"
                    vdb_cache.mkdir(parents=True, exist_ok=True)
                    env["VDB_CONFIG"] = str(vdb_cache)

                proc_kwargs = {"env": env}
                if cwd:
                    proc_kwargs["cwd"] = cwd

                if log_dir:
                    proc_kwargs["stdout"] = stdout_fh
                    proc_kwargs["stderr"] = stderr_fh

                proc = subprocess.Popen(command, **proc_kwargs)

                if subcommand == "quant":
                    rc = monitor_subprocess_sample_progress(
                        process=proc,
                        watch_dir=watch_dir,
                        heartbeat_path=hb_path,
                        completion_glob="*/abundance.tsv",
                        total_samples=total_runs,
                        heartbeat_interval=heartbeat_interval,
                        show_progress=show_progress,
                        desc=f"amalgkit {subcommand}",
                    )
                elif subcommand == "metadata":
                    expected = ["metadata.tsv", "metadata_original.tsv", "pivot_selected.tsv", "pivot_qualified.tsv"]
                    rc = monitor_subprocess_file_count(
                        process=proc,
                        watch_dir=watch_dir,
                        heartbeat_path=hb_path,
                        expected_files=expected,
                        heartbeat_interval=heartbeat_interval,
                        show_progress=show_progress,
                        desc=f"amalgkit {subcommand}",
                    )
                elif subcommand in {"merge", "cstmm", "curate", "csca"}:
                    expected = (
                        ["merged_abundance.tsv"]
                        if subcommand == "merge"
                        else (
                            ["cstmm.tsv"]
                            if subcommand == "cstmm"
                            else ["tables"] if subcommand == "curate" else ["csca.tsv"]
                        )
                    )
                    rc = monitor_subprocess_file_count(
                        process=proc,
                        watch_dir=watch_dir,
                        heartbeat_path=hb_path,
                        expected_files=expected,
                        heartbeat_interval=heartbeat_interval,
                        show_progress=show_progress,
                        desc=f"amalgkit {subcommand}",
                    )
                else:
                    rc, _ = monitor_subprocess_directory_growth(
                        process=proc,
                        watch_dir=watch_dir,
                        heartbeat_path=hb_path,
                        total_bytes=total_bytes,
                        heartbeat_interval=heartbeat_interval,
                        show_progress=show_progress,
                        desc=f"amalgkit {subcommand}",
                    )
        # Non-monitored path
        proc_kwargs = kwargs.copy()
        if cwd:
            proc_kwargs["cwd"] = cwd

        if log_dir:
            proc_kwargs["stdout"] = stdout_fh
            proc_kwargs["stderr"] = stderr_fh
        else:
            if "capture_output" not in proc_kwargs:
                proc_kwargs["capture_output"] = True
                proc_kwargs["text"] = True

        result = subprocess.run(command, **proc_kwargs)
        if log_dir:
            stdout_fh.close()
            stderr_fh.close()
        return result

    except Exception as e:
        logger.error(f"Error running amalgkit {subcommand}: {e}")
        # Return failed process if check=False, otherwise re-raise if it was a budget error
        if kwargs.get("check"):
            raise
        return subprocess.CompletedProcess(args=command, returncode=1, stdout="", stderr=str(e))


def _estimate_total_bytes_from_metadata(metadata_path: Path) -> int | None:
    """Estimate total download size from amalgkit metadata.tsv.

    We sum the `size` column when present (typically bytes of the run payload).
    Returns None if size is unavailable.
    """
    if not metadata_path.exists():
        return None

    try:
        with open(metadata_path, "r", encoding="utf-8") as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            if reader.fieldnames is None or "size" not in reader.fieldnames:
                return None
            total = 0
            for row in reader:
                val = (row.get("size") or "").strip()
                if not val:
                    continue
                try:
                    total += int(float(val))
                except ValueError:
                    continue
            return total if total > 0 else None
    except (OSError, IOError, csv.Error) as e:
        logger.debug(f"Could not estimate bytes from metadata: {e}")
        return None


def _estimate_total_runs_from_metadata(metadata_path: Path) -> int | None:
    """Estimate number of runs/samples from metadata TSV (row-per-run)."""
    if not metadata_path.exists():
        return None
    try:
        with open(metadata_path, "r", encoding="utf-8") as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            if reader.fieldnames is None or "run" not in reader.fieldnames:
                return None
            n = 0
            for _ in reader:
                n += 1
            return n if n > 0 else None
    except (OSError, IOError, csv.Error) as e:
        logger.debug(f"Could not estimate runs from metadata: {e}")
        return None


def _infer_metadata_path_from_out_dir(out_dir: Path) -> Path | None:
    """Best-effort: locate metadata.tsv under a typical amalgkit work dir."""
    candidates = [
        out_dir / "metadata" / "metadata_selected.tsv",
        out_dir / "metadata" / "metadata.tsv",
        out_dir / "work" / "metadata" / "metadata_selected.tsv",
        out_dir / "work" / "metadata" / "metadata.tsv",
    ]
    for c in candidates:
        if c.exists():
            return c
    return None


def metadata(params: AmalgkitParams | Dict[str, Any] | None = None, **kwargs: Any) -> subprocess.CompletedProcess[str]:
    """Run amalgkit metadata command.

    Args:
        params: Amalgkit parameters (AmalgkitParams object or dict)
        **kwargs: Additional arguments

    Returns:
        CompletedProcess instance
    """
    return run_amalgkit("metadata", params, **kwargs)


def integrate(params: AmalgkitParams | Dict[str, Any] | None = None, **kwargs: Any) -> subprocess.CompletedProcess[str]:
    """Run amalgkit integrate command.

    Args:
        params: Amalgkit parameters
        **kwargs: Additional arguments

    Returns:
        CompletedProcess instance
    """
    return run_amalgkit("integrate", params, **kwargs)


def config(params: AmalgkitParams | Dict[str, Any] | None = None, **kwargs: Any) -> subprocess.CompletedProcess[str]:
    """Run amalgkit config command.

    Args:
        params: Amalgkit parameters
        **kwargs: Additional arguments

    Returns:
        CompletedProcess instance
    """
    return run_amalgkit("config", params, **kwargs)


def select(params: AmalgkitParams | Dict[str, Any] | None = None, **kwargs: Any) -> subprocess.CompletedProcess[str]:
    """Run amalgkit select command.

    Args:
        params: Amalgkit parameters
        **kwargs: Additional arguments

    Returns:
        CompletedProcess instance
    """
    return run_amalgkit("select", params, **kwargs)


def getfastq(params: AmalgkitParams | Dict[str, Any] | None = None, **kwargs: Any) -> subprocess.CompletedProcess[str]:
    """Run amalgkit getfastq command with retry logic.

    Args:
        params: Amalgkit parameters
        **kwargs: Additional arguments

    Returns:
        CompletedProcess instance
    """
    # Check if parallel execution is requested
    jobs = 1
    if isinstance(params, dict):
        jobs = int(params.get("jobs", 1))
    elif isinstance(params, AmalgkitParams):
        jobs = int(params.extra_params.get("jobs", 1))

    # Retry configuration
    max_retries = 3
    retry_delay = 5  # seconds, will be multiplied by attempt number

    last_result = None
    for attempt in range(1, max_retries + 1):
        if attempt > 1:
            logger.info(f"Retry attempt {attempt}/{max_retries} for getfastq...")
            import time

            time.sleep(retry_delay * attempt)  # Exponential backoff

        try:
            if jobs > 1:
                result = _run_parallel_getfastq(jobs, params, **kwargs)
            else:
                result = run_amalgkit("getfastq", params, **kwargs)

            last_result = result

            # Check if successful (return code 0)
            if result.returncode == 0:
                return result

            # Check if 'check' was True and command failed
            check = kwargs.get("check", False)
            if result.returncode != 0 and check:
                raise subprocess.CalledProcessError(
                    result.returncode, command, output=result.stdout, stderr=result.stderr
                )

            # Check if partial success (some files downloaded)
            # Don't retry if the command ran but just had no work to do
            stderr_lower = (result.stderr or "").lower()
            if "no samples" in stderr_lower or "nothing to download" in stderr_lower:
                logger.info("getfastq completed with no samples to process")
                return result

            # Log retry reason
            logger.warning(f"getfastq failed with return code {result.returncode} (attempt {attempt}/{max_retries})")
            if result.stderr:
                logger.warning(f"Error: {result.stderr[:500]}")

        except Exception as e:
            logger.error(f"getfastq exception on attempt {attempt}/{max_retries}: {e}")
            last_result = subprocess.CompletedProcess(
                args=["amalgkit", "getfastq"], returncode=1, stdout="", stderr=str(e)
            )

    # Return last result after all retries exhausted
    logger.error(f"getfastq failed after {max_retries} attempts")
    return last_result or subprocess.CompletedProcess(
        args=["amalgkit", "getfastq"], returncode=1, stdout="", stderr="All retry attempts failed"
    )


def _split_metadata_by_worker(metadata_path: Path, n_workers: int) -> List[Path]:
    """Split metadata file into chunks for parallel processing.

    Args:
        metadata_path: Path to original metadata file
        n_workers: Number of workers/chunks

    Returns:
        List of paths to temporary metadata chunk files
    """
    if not metadata_path.exists():
        raise FileNotFoundError(f"Metadata file not found: {metadata_path}")

    # Read all rows
    rows = []
    header = []
    with open(metadata_path, "r", newline="") as f:
        reader = csv.reader(f, delimiter="\t")
        try:
            header = next(reader)
            rows = list(reader)
        except StopIteration:
            pass

    if not rows:
        return [metadata_path]

    # Calculate chunk size
    n_samples = len(rows)
    chunk_size = math.ceil(n_samples / n_workers)

    chunk_paths = []
    temp_dir = Path(tempfile.mkdtemp(prefix="amalgkit_parallel_"))

    for i in range(n_workers):
        start_idx = i * chunk_size
        end_idx = min((i + 1) * chunk_size, n_samples)

        if start_idx >= n_samples:
            break

        chunk_rows = rows[start_idx:end_idx]
        if not chunk_rows:
            continue

        chunk_path = temp_dir / f"metadata_chunk_{i}.tsv"
        with open(chunk_path, "w", newline="") as f:
            writer = csv.writer(f, delimiter="\t")
            writer.writerow(header)
            writer.writerows(chunk_rows)

        chunk_paths.append(chunk_path)

    return chunk_paths


def _run_parallel_getfastq(
    jobs: int, params: AmalgkitParams | Dict[str, Any] | None, **kwargs: Any
) -> subprocess.CompletedProcess[str]:
    """Execute getfastq in parallel using multiple workers.

    Args:
        jobs: Number of parallel jobs
        params: Amalgkit parameters
        **kwargs: Additional arguments

    Returns:
        Combined CompletedProcess
    """
    # Extract metadata path
    metadata_path = None
    if isinstance(params, dict):
        metadata_path = Path(params.get("metadata", ""))
    else:
        metadata_path = Path(params.extra_params.get("metadata", ""))

    if not metadata_path.exists():
        # Fallback to single process if metadata not found
        logger.warning(f"Metadata not found at {metadata_path}, falling back to sequential execution")
        return run_amalgkit("getfastq", params, **kwargs)

    # Split metadata
    try:
        chunk_paths = _split_metadata_by_worker(metadata_path, jobs)
    except Exception as e:
        logger.warning(f"Failed to split metadata: {e}, falling back to sequential execution")
        return run_amalgkit("getfastq", params, **kwargs)

    if len(chunk_paths) <= 1:
        # Only one chunk needed (small sample size), run normally
        # MUST strip 'jobs' from params because CLI doesn't support it
        if isinstance(params, dict):
            params.pop("jobs", None)
        elif isinstance(params, AmalgkitParams):
            params.extra_params.pop("jobs", None)

        return run_amalgkit("getfastq", params, **kwargs)

    logger.info(f"Parallelizing getfastq with {len(chunk_paths)} workers")

    # Create worker params
    worker_futures = []
    results = []

    with concurrent.futures.ThreadPoolExecutor(max_workers=len(chunk_paths)) as executor:
        for i, chunk_path in enumerate(chunk_paths):
            # Create specific params for this worker
            if isinstance(params, dict):
                worker_params = params.copy()
                worker_params["metadata"] = str(chunk_path)
                # Ensure jobs param is removed/reset to prevent recursion if passed down
                worker_params.pop("jobs", None)
            else:
                # Clone object
                worker_params = AmalgkitParams(
                    work_dir=params.work_dir,
                    threads=params.threads,  # Each worker gets full requested threads? Or divide?
                    # Usually better to set threads=1 per worker if bandwidth limited,
                    # OR keep threads high if fasterq-dump needs it.
                    # We'll rely on workflow.py to set appropriate threads per job.
                    species_list=params.species_list,
                    **params.extra_params.copy(),
                )
                worker_params.extra_params["metadata"] = str(chunk_path)
                worker_params.extra_params.pop("jobs", None)

            # Submit job
            # We disable individual monitoring for workers to avoid console fighting
            # Instead we monitor the main directory globally in the main thread?
            # Or just let them run. Since `run_amalgkit` runs a subprocess, blocking here is fine.
            # We pass `monitor=False` to workers to prevent them from starting their own monitoring loop?
            # Yes, `monitor=False` is safer.
            # But we want to monitor overall progress.
            # We can start a separate monitor thread for the main directory.

            worker_kwargs = kwargs.copy()
            worker_kwargs["monitor"] = False  # Disable individual monitoring

            future = executor.submit(run_amalgkit, "getfastq", worker_params, **worker_kwargs)
            worker_futures.append(future)

        # Optional: Start a global monitor here if requested
        # For now, we rely on the fact that individual workers log to stdout/stderr.
        # But we disabled their monitor loop (heartbeat), not the logging.
        # run_amalgkit logs "Running amalgkit command..."

        # Wait for all with heartbeat
        start_time = time.time()
        while True:
            # Check for completion
            done, not_done = concurrent.futures.wait(
                worker_futures, timeout=60, return_when=concurrent.futures.FIRST_COMPLETED
            )

            # Log heartbeat
            # NOTE: Local import to avoid circular dependency with workflow.py
            from metainformant.rna.engine.workflow import _log_heartbeat

            _log_heartbeat(f"getfastq parallel ({len(done)}/{len(chunk_paths)} chunks done)", start_time=start_time)

            # Process completed
            for future in done:
                if future in worker_futures:
                    worker_futures.remove(future)
                    try:
                        res = future.result()
                        results.append(res)
                        logger.debug(f"Worker completed with rc={res.returncode}")
                        if res.stdout:
                            logger.debug(f"Worker stdout:\n{res.stdout}")
                        if res.stderr:
                            logger.debug(f"Worker stderr:\n{res.stderr}")
                    except Exception as e:
                        logger.error(f"Worker failed: {e}")
                        results.append(subprocess.CompletedProcess(args=[], returncode=1, stderr=str(e)))

            if not worker_futures:
                break

    # Cleanup temp chunks
    try:
        shutil.rmtree(chunk_paths[0].parent)
    except (OSError, PermissionError, FileNotFoundError) as e:
        logger.debug(f"Failed to cleanup temp chunk directory: {e}")

    # Aggregate results for return
    # Combine stdout and stderr from all workers so summary parsing works for all
    combined_stdout = "\n".join([r.stdout for r in results if r.stdout])
    combined_stderr = "\n".join([r.stderr for r in results if r.stderr])

    # Check for failures
    failed = [r for r in results if r.returncode != 0]
    if failed:
        # Return a failed process with combined output
        return subprocess.CompletedProcess(
            args=failed[0].args, returncode=failed[0].returncode, stdout=combined_stdout, stderr=combined_stderr
        )

    # Return success with combined output
    return subprocess.CompletedProcess(
        args=results[0].args if results else [], returncode=0, stdout=combined_stdout, stderr=combined_stderr
    )


def quant(params: AmalgkitParams | Dict[str, Any] | None = None, **kwargs: Any) -> subprocess.CompletedProcess[str]:
    """Run amalgkit quant command.

    Args:
        params: Amalgkit parameters
        **kwargs: Additional arguments

    Returns:
        CompletedProcess instance
    """
    return run_amalgkit("quant", params, **kwargs)


def merge(params: AmalgkitParams | Dict[str, Any] | None = None, **kwargs: Any) -> subprocess.CompletedProcess[str]:
    """Run amalgkit merge command.

    Args:
        params: Amalgkit parameters
        **kwargs: Additional arguments

    Returns:
        CompletedProcess instance
    """
    return run_amalgkit("merge", params, **kwargs)


def cstmm(params: AmalgkitParams | Dict[str, Any] | None = None, **kwargs: Any) -> subprocess.CompletedProcess[str]:
    """Run amalgkit cstmm command.

    Args:
        params: Amalgkit parameters
        **kwargs: Additional arguments

    Returns:
        CompletedProcess instance
    """
    return run_amalgkit("cstmm", params, **kwargs)


def curate(params: AmalgkitParams | Dict[str, Any] | None = None, **kwargs: Any) -> subprocess.CompletedProcess[str]:
    """Run amalgkit curate command.

    Args:
        params: Amalgkit parameters
        **kwargs: Additional arguments

    Returns:
        CompletedProcess instance
    """
    return run_amalgkit("curate", params, **kwargs)


def csca(params: AmalgkitParams | Dict[str, Any] | None = None, **kwargs: Any) -> subprocess.CompletedProcess[str]:
    """Run amalgkit csca command.

    Args:
        params: Amalgkit parameters
        **kwargs: Additional arguments

    Returns:
        CompletedProcess instance
    """
    return run_amalgkit("csca", params, **kwargs)


def sanity(params: AmalgkitParams | Dict[str, Any] | None = None, **kwargs: Any) -> subprocess.CompletedProcess[str]:
    """Run amalgkit sanity command.

    Args:
        params: Amalgkit parameters
        **kwargs: Additional arguments

    Returns:
        CompletedProcess instance
    """
    return run_amalgkit("sanity", params, **kwargs)
