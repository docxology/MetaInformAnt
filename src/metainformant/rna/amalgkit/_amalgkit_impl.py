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
import importlib.metadata
import importlib.util
import json
import math
import os
import re
import shutil
import subprocess
import sys
import tempfile
import threading
import time
import uuid
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, Union

from metainformant.core.io.download import (
    monitor_subprocess_directory_growth,
    monitor_subprocess_file_count,
    monitor_subprocess_sample_progress,
)
from metainformant.core.utils import logging
from metainformant.rna.amalgkit.commands import (
    METAINFORMANT_ONLY_STEP_KEYS,
    NO_THREAD_STEPS,
    SUPPORTED_CLI_OPTIONS,
    YES_NO_CLI_OPTIONS,
)
from metainformant.rna.amalgkit.sra_environment import build_sra_environment
from metainformant.rna.core.sample_utils import read_sample_ids_from_metadata

logger = logging.get_logger(__name__)

# The MetaInformAnt wrapper is validated against exactly one upstream Amalgkit
# build. Upstream tags (and bioconda) only publish x.y.0 releases, so a patch
# build is pinned by its immutable commit, never by a tag or a moving branch.
# The same commit appears in pyproject.toml (the ``rna`` and ``all`` extras)
# and uv.lock; tests/rna/test_amalgkit_pin.py fails when the three disagree.
# The version and revision are the only declared values; everything else below
# is derived from them.
AMALGKIT_REPOSITORY_URL = "https://github.com/kfuku52/amalgkit"
REQUIRED_AMALGKIT_VERSION = "0.16.108"
AMALGKIT_SOURCE_REVISION = "1d1566638abe32e79ee262193bb88769c3482e76"
# Source label recorded in provenance sidecars. It is not a git tag: no
# ``v0.16.108`` tag exists upstream.
AMALGKIT_RELEASE_TAG = f"master@{AMALGKIT_SOURCE_REVISION[:12]}"
AMALGKIT_INSTALL_SPEC = f"git+{AMALGKIT_REPOSITORY_URL}@{AMALGKIT_SOURCE_REVISION}"
# The only supported way to obtain the pinned build: install it from uv.lock.
# MetaInformAnt never installs or upgrades packages at runtime.
AMALGKIT_SYNC_HINT = "uv sync --frozen --extra rna"
YES_NO_CLI_FLAGS = {
    "redo",
    "resolve_names",
    "mark_redundant_biosamples",
    "aws",
    "ncbi",
    "gcp",
    "ena",
    "ddbj",
    "fastp",
    "fasterq_size_check",
    "rrna_filter",
    "contam_filter",
    "remove_sra",
    "remove_tmp",
    "dump_print",
    "fastp_print",
    "build_index",
    "clean_fastq",
    "overwrite",
    "clip_negative",
    "maintain_zero",
}
_ALL_AMALGKIT_OPTIONS = frozenset().union(*SUPPORTED_CLI_OPTIONS.values())
_warned_dropped_keys: set[tuple[str, str]] = set()
_warning_lock = threading.Lock()


def _yes_no_flags(subcommand: str | None) -> frozenset[str] | set[str]:
    """Return the yes/no-valued options for ``subcommand`` (all known ones when unknown)."""

    if subcommand is not None and subcommand in YES_NO_CLI_OPTIONS:
        return YES_NO_CLI_OPTIONS[subcommand]
    return YES_NO_CLI_FLAGS


def _report_dropped_key(subcommand: str, key: str) -> None:
    """Log a key that ``build_cli_args`` leaves out of an Amalgkit command line.

    MetaInformAnt-only keys and options that belong to a different Amalgkit
    command are expected in shared parameter maps and are logged at debug level.
    A key that no command of the pinned build accepts is a typo or an option
    this build lacks; it is reported once per (command, key) as a warning so a
    silently ignored configuration value cannot change results unnoticed.
    """

    if key in METAINFORMANT_ONLY_STEP_KEYS or key in _ALL_AMALGKIT_OPTIONS:
        logger.debug("amalgkit %s: leaving out %r (not an option of this command)", subcommand, key)
        return
    with _warning_lock:
        if (subcommand, key) in _warned_dropped_keys:
            return
        _warned_dropped_keys.add((subcommand, key))
    logger.warning(
        "amalgkit %s: dropping configuration key %r: no command of the pinned Amalgkit %s accepts it",
        subcommand,
        key,
        REQUIRED_AMALGKIT_VERSION,
    )


class AmalgkitParams:
    """Parameters for amalgkit commands."""

    def __init__(
        self, work_dir: Union[str, Path], threads: int = 8, species_list: Optional[List[str]] = None, **kwargs: Any
    ):
        """Initialize amalgkit parameters.

        Args:
            work_dir: Working directory for amalgkit operations. Output will be structured here.
            threads: Number of CPU threads to use for parallel operations. Defaults to 8.
            species_list: List of scientific names (e.g. ['Pogonomyrmex barbatus']) to process.
            **kwargs: Additional parameters passed directly to the amalgkit CLI.
                     Common keys: 'redo', 'fastp', 'aws', 'internal_jobs', etc.
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
    args: list[str] = []

    if params is None:
        return args
    yes_no_flags = _yes_no_flags(subcommand)

    # Handle dict parameters
    if isinstance(params, dict):
        # amalgkit uses --out_dir for the main working directory in all subcommands
        if "out_dir" in params:
            args.extend(["--out_dir", str(params["out_dir"])])
        elif "work_dir" in params:
            args.extend(["--out_dir", str(params["work_dir"])])

        if "threads" in params and subcommand not in NO_THREAD_STEPS:
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
            "log_dir",
            "num_download_workers",
        }

        for key, value in params.items():
            if key in skip_keys or value is None:
                continue

            # Skip arguments not supported by the pinned command. The planner
            # carries a few metadata-only values for orchestration; they must
            # never be emitted as raw CLI flags. Unknown keys are reported.
            normalized_key = key.replace("-", "_")
            supported = SUPPORTED_CLI_OPTIONS.get(subcommand or "")
            if supported is not None and normalized_key not in supported:
                _report_dropped_key(subcommand or "", normalized_key)
                continue

            # Skip arguments with command-specific incompatibilities that are
            # not represented by a simple allow-list (kept for clarity when a
            # shared parameter map is passed by a caller).
            unsupported_by_step = {
                "cstmm": {"out", "threads", "jobs"},
                "merge": {"out", "redo", "jobs"},
                "sanity": {"out", "redo", "jobs"},
            }
            if key in unsupported_by_step.get(subcommand or "", set()):
                continue
            if subcommand in ("getfastq", "metadata", "quant") and key == "jobs":
                continue

            # Amalgkit spells every CLI flag with underscores.
            cli_key = normalized_key

            if isinstance(value, bool):
                if cli_key in yes_no_flags:
                    args.extend([f"--{cli_key}", "yes" if value else "no"])
                elif value:
                    args.append(f"--{cli_key}")
            elif isinstance(value, (int, float)):
                if key == "threads" and subcommand in NO_THREAD_STEPS:
                    continue
                args.extend([f"--{cli_key}", str(value)])
            elif isinstance(value, list):
                # Handle lists (multiple arguments)
                if key == "species-list" or key == "species_list" or key == "species":
                    if subcommand not in (
                        "integrate",
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
                if subcommand not in ("integrate", "help", "merge", "metadata", "getfastq", "quant"):
                    args.extend(["--species", str(value)])
            else:
                args.extend([f"--{cli_key}", str(value)])
        return args

    # Handle AmalgkitParams object
    args.extend(["--out_dir", str(params.work_dir)])
    if subcommand not in NO_THREAD_STEPS:
        args.extend(["--threads", str(params.threads)])

    supported = SUPPORTED_CLI_OPTIONS.get(subcommand or "")
    if params.species_list and (supported is None or "species" in supported):
        for species in params.species_list:
            args.extend(["--species", species])

    # Add extra parameters
    for key, value in params.extra_params.items():
        cli_key = key.replace("-", "_")
        if supported is not None and cli_key not in supported:
            _report_dropped_key(subcommand or "", cli_key)
            continue
        if isinstance(value, bool):
            if cli_key in yes_no_flags:
                args.extend([f"--{cli_key}", "yes" if value else "no"])
            elif value:
                args.append(f"--{cli_key}")
        elif isinstance(value, (int, float)):
            args.extend([f"--{cli_key}", str(value)])
        elif isinstance(value, str):
            args.extend([f"--{cli_key}", value])
        elif isinstance(value, list):
            for item in value:
                args.extend([f"--{cli_key}", str(item)])

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


def _resolve_amalgkit_command() -> list[str] | None:
    """Resolve the project-owned Amalgkit command without mutating PATH.

    Hosted runners can invoke pytest through ``uv`` while the package's CLI is
    installed in the project environment.  Relying only on the inherited PATH
    makes capability checks disagree with the subprocess that the workflow
    actually owns.  Resolution remains read-only and bounded to the active
    interpreter environment, PATH, and this checkout's ``.venv``.

    The interpreter's own console script comes first: it belongs to the same
    environment whose installed distribution :func:`verify_installed_amalgkit`
    checks, so the command that runs is the build that was verified rather than
    whichever ``amalgkit`` happens to be earlier on PATH.
    """
    candidates: list[str] = [str(Path(sys.executable).with_name("amalgkit"))]
    virtual_env = os.environ.get("VIRTUAL_ENV")
    if virtual_env:
        candidates.append(str(Path(virtual_env) / "bin" / "amalgkit"))
    discovered = shutil.which("amalgkit")
    if discovered:
        candidates.append(discovered)

    repository_root = Path(__file__).resolve().parents[4]
    candidates.append(str(repository_root / ".venv" / "bin" / "amalgkit"))

    for candidate in candidates:
        path = Path(candidate).expanduser()
        if not path.is_absolute():
            path = (Path.cwd() / path).resolve()
        if path.is_file() and os.access(path, os.X_OK):
            return [str(path)]

    # Some Amalgkit source revisions install the Python package and its
    # ``__main__`` module but omit a console-script entry point.  Treat the
    # module runner as the equivalent project-owned CLI capability.
    try:
        module_runner_available = importlib.util.find_spec("amalgkit.__main__") is not None
    except ModuleNotFoundError:
        module_runner_available = False
    if module_runner_available:
        return [str(Path(sys.executable).resolve()), "-m", "amalgkit"]
    return None


def check_cli_available() -> Tuple[bool, str]:
    """Check if amalgkit CLI is available.

    Returns:
        Tuple of (available, message)
    """
    command_prefix = _resolve_amalgkit_command()
    if command_prefix is None:
        return False, "amalgkit CLI not found in PATH or the project environment"

    try:
        result = subprocess.run([*command_prefix, "--help"], capture_output=True, text=True, timeout=10)
        if result.returncode == 0:
            return True, f"amalgkit CLI is available: {' '.join(command_prefix)}"
        else:
            return False, f"amalgkit CLI returned error: {result.stderr}"
    except FileNotFoundError:
        return False, f"amalgkit CLI disappeared after resolution: {' '.join(command_prefix)}"
    except subprocess.TimeoutExpired:
        return False, "amalgkit CLI check timed out"
    except (OSError, PermissionError) as e:
        return False, f"Error checking amalgkit CLI: {e}"


def parse_and_check_version(
    version_output: str,
    min_version: str = REQUIRED_AMALGKIT_VERSION,
    *,
    exact: bool = False,
) -> Tuple[bool, str]:
    """Parse version string and check a minimum or exact requirement.

    Args:
        version_output: Output from 'amalgkit --version' (e.g. "amalgkit 0.4.1")
        min_version: Minimum required version, or the expected version when
            ``exact`` is true.
        exact: Require the parsed version to equal ``min_version``.

    Returns:
        Tuple of (valid, message)
    """
    output = version_output.strip()
    if not output:
        return False, "Empty version output"

    # Amalgkit emits ``amalgkit version X.Y.Z``.
    version_match = re.search(r"(?<!\d)(\d+(?:\.\d+){1,2})(?!\d)", output)
    if version_match is None:
        return False, f"Unexpected version format: {output}"
    current_version = version_match.group(1)

    # Simple semantic version comparison
    # We assume X.Y.Z format
    try:
        current_parts = [int(x) for x in current_version.split(".")]
        min_parts = [int(x) for x in min_version.split(".")]

        # Pad with zeros if needed
        while len(current_parts) < 3:
            current_parts.append(0)
        while len(min_parts) < 3:
            min_parts.append(0)

        current_tuple = tuple(current_parts)
        required_tuple = tuple(min_parts)
        if exact and current_tuple == required_tuple:
            return True, f"Version {current_version} matches required {min_version}"
        if not exact and current_tuple >= required_tuple:
            return True, f"Version {current_version} meets requirement >={min_version}"
        if exact:
            return False, f"Version {current_version} does not match required {min_version}"
        else:
            return False, f"Version {current_version} is older than required {min_version}"

    except ValueError:
        # Exact-version validation must fail closed when the CLI emits an
        # unrecognized value; accepting it would defeat the pinned contract.
        logger.warning(f"Could not parse version string: {current_version}")
        return False, f"Could not parse version {current_version}"


@dataclass(frozen=True)
class InstalledAmalgkit:
    """Identity of the Amalgkit distribution installed in this interpreter.

    ``commit_id`` and ``url`` come from the PEP 610 ``direct_url.json`` record
    that uv writes for a git-sourced install; they are ``None`` for an install
    whose origin was not recorded (for example a bioconda or index package).
    """

    version: str
    url: str | None
    commit_id: str | None
    requested_revision: str | None
    editable: bool


def _normalize_repository_url(url: str) -> str:
    normalized = url.strip().removeprefix("git+").rstrip("/")
    return normalized.removesuffix(".git").lower()


def installed_amalgkit_distribution() -> InstalledAmalgkit | None:
    """Read the installed Amalgkit distribution's version and VCS origin.

    Returns:
        The distribution identity, or ``None`` when this interpreter has no
        ``amalgkit`` distribution.
    """

    try:
        distribution = importlib.metadata.distribution("amalgkit")
    except importlib.metadata.PackageNotFoundError:
        return None
    payload: dict[str, Any] = {}
    raw = distribution.read_text("direct_url.json")
    if raw:
        try:
            loaded = json.loads(raw)
        except json.JSONDecodeError:
            loaded = {}
        if isinstance(loaded, dict):
            payload = loaded
    raw_vcs, raw_dir = payload.get("vcs_info"), payload.get("dir_info")
    vcs_info: dict[str, Any] = raw_vcs if isinstance(raw_vcs, dict) else {}
    dir_info: dict[str, Any] = raw_dir if isinstance(raw_dir, dict) else {}
    url = payload.get("url")
    commit_id = vcs_info.get("commit_id")
    requested_revision = vcs_info.get("requested_revision")
    return InstalledAmalgkit(
        version=str(distribution.version),
        url=str(url) if isinstance(url, str) else None,
        commit_id=str(commit_id) if isinstance(commit_id, str) else None,
        requested_revision=str(requested_revision) if isinstance(requested_revision, str) else None,
        editable=bool(dir_info.get("editable", False)),
    )


def verify_installed_amalgkit(installed: InstalledAmalgkit | None = None) -> Tuple[bool, str]:
    """Fail closed unless this interpreter holds exactly the pinned Amalgkit build.

    The installed distribution must report :data:`REQUIRED_AMALGKIT_VERSION`
    and its ``direct_url.json`` must name :data:`AMALGKIT_REPOSITORY_URL` at
    commit :data:`AMALGKIT_SOURCE_REVISION`. A version string alone is not
    accepted: upstream does not tag patch builds, so two different commits can
    report the same version.

    Args:
        installed: Distribution identity to check; read from this interpreter
            when omitted.

    Returns:
        Tuple of (matches, message). A mismatch message ends with the
        lock-based provisioning hint.
    """

    if installed is None:
        installed = installed_amalgkit_distribution()
    if installed is None:
        return False, f"amalgkit is not installed for {sys.executable}. Run `{AMALGKIT_SYNC_HINT}`."
    problems: list[str] = []
    if installed.version != REQUIRED_AMALGKIT_VERSION:
        problems.append(f"version {installed.version} is not {REQUIRED_AMALGKIT_VERSION}")
    if installed.commit_id is None:
        problems.append("no VCS commit is recorded in direct_url.json, so the build cannot be identified")
    elif installed.commit_id != AMALGKIT_SOURCE_REVISION:
        problems.append(f"commit {installed.commit_id} is not {AMALGKIT_SOURCE_REVISION}")
    if installed.url is not None and _normalize_repository_url(installed.url) != _normalize_repository_url(
        AMALGKIT_REPOSITORY_URL
    ):
        problems.append(f"source {installed.url} is not {AMALGKIT_REPOSITORY_URL}")
    if installed.editable:
        problems.append("it is an editable checkout, whose bytes can change after installation")
    if problems:
        return False, f"Installed amalgkit is not the pinned build ({'; '.join(problems)}). Run `{AMALGKIT_SYNC_HINT}`."
    return True, f"amalgkit {installed.version} at commit {AMALGKIT_SOURCE_REVISION[:12]} matches the pin"


def validate_amalgkit_version(
    min_version: str = REQUIRED_AMALGKIT_VERSION,
    *,
    exact: bool = True,
) -> Tuple[bool, str]:
    """Validate the installed Amalgkit CLI against the project contract.

    For the pinned contract (``exact`` with the required version) the CLI's
    reported version must match and the interpreter's installed distribution
    must come from the pinned commit (:func:`verify_installed_amalgkit`).

    Args:
        min_version: Required version string.
        exact: Require the pinned version by default; set false only for a
            deliberate compatibility check against a lower bound.

    Returns:
        Tuple of (valid, message)
    """
    command_prefix = _resolve_amalgkit_command()
    if command_prefix is None:
        return False, f"amalgkit CLI not found in PATH or the project environment. Run `{AMALGKIT_SYNC_HINT}`."

    try:
        result = subprocess.run([*command_prefix, "--version"], capture_output=True, text=True, timeout=5)
    except (subprocess.SubprocessError, OSError, TimeoutError) as e:
        return False, f"Version check error: {e}"
    if result.returncode != 0:
        return False, f"Version check failed: {result.stderr}"

    valid, message = parse_and_check_version(result.stdout, min_version, exact=exact)
    if not valid or not exact or min_version != REQUIRED_AMALGKIT_VERSION:
        return valid, message
    return verify_installed_amalgkit()


def ensure_cli_available(
    *, min_version: str = REQUIRED_AMALGKIT_VERSION, exact: bool = True
) -> Tuple[bool, str, Dict | None]:
    """Ensure the Amalgkit CLI is available and is the pinned build.

    This check never installs, upgrades or downgrades anything: runtime
    installs bypass ``uv.lock`` and have silently replaced Amalgkit and its
    transitive dependencies in the past. A missing or different build fails
    closed with the provisioning hint :data:`AMALGKIT_SYNC_HINT`.

    Args:
        min_version: Required version
        exact: Require the pinned version by default.

    Returns:
        Tuple of (success, message, version_info). ``version_info`` reports the
        observed version, whether it is valid, ``attempted: False`` and the
        provisioning hint.
    """
    installed = installed_amalgkit_distribution()
    version_info: Dict[str, Any] = {
        "version": installed.version if installed is not None else None,
        "commit_id": installed.commit_id if installed is not None else None,
        "valid": False,
        "attempted": False,
        "install_hint": AMALGKIT_SYNC_HINT,
    }
    available, message = check_cli_available()
    if not available:
        return False, f"{message}. Run `{AMALGKIT_SYNC_HINT}`.", version_info

    valid_version, version_msg = validate_amalgkit_version(min_version, exact=exact)
    version_info["valid"] = valid_version
    if not valid_version:
        logger.warning(f"Amalgkit is not usable: {version_msg}")
        hint = "" if AMALGKIT_SYNC_HINT in version_msg else f" Run `{AMALGKIT_SYNC_HINT}`."
        return False, f"Amalgkit available but not the required build: {version_msg}{hint}", version_info
    return True, message, version_info


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

    The standard ``timeout`` subprocess keyword also bounds monitored steps.
    A monitored timeout terminates the child, writes a timeout heartbeat, and
    returns exit code 124. Use ``monitor_timeout`` when the subprocess timeout
    and the monitoring wall-clock limit need to differ.
    """
    command = build_amalgkit_command(subcommand, params)
    command_prefix = _resolve_amalgkit_command()
    if command_prefix is not None:
        command = [*command_prefix, *command[1:]]

    logger.info(f"Running amalgkit command: {' '.join(command)}")

    # Extract custom keys
    cwd = kwargs.pop("cwd", None) or kwargs.pop("work_dir", None)
    log_dir = kwargs.pop("log_dir", None)
    params_monitor = kwargs.pop("monitor", None)  # Rename to avoid conflict if any
    monitor_timeout = kwargs.pop("monitor_timeout", None)
    if monitor_timeout is None:
        monitor_timeout = kwargs.get("timeout")
    heartbeat_interval = int(kwargs.pop("heartbeat_interval", 5))
    show_progress = bool(kwargs.pop("show_progress", True))

    # Determine monitoring strategy
    monitor = False
    if params_monitor is None:
        # Default monitoring for long-running steps
        monitor = subcommand in {
            "metadata",
            "integrate",
            "getfastq",
            "quant",
            "merge",
            "wsfilter",
            "cstmm",
            "csfilter",
            "finalize",
        }
    else:
        monitor = bool(params_monitor)

    # Setup file handles if logging requested
    stdout_fh = None
    stderr_fh = None
    if log_dir:
        ts = int(time.time())
        log_dir_path = Path(log_dir)
        log_dir_path.mkdir(parents=True, exist_ok=True)
        stdout_log = log_dir_path / f"{ts}.{subcommand}.stdout.log"
        stderr_log = log_dir_path / f"{ts}.{subcommand}.stderr.log"
        stdout_fh = open(stdout_log, "w")
        stderr_fh = open(stderr_log, "w")

    try:
        # Environmental setup
        env = os.environ.copy()
        if subcommand == "getfastq":
            # The installed Amalgkit process inherits this environment when it
            # invokes fasterq-dump. Keep its cache and temporary files inside
            # the selected campaign root without mutating site-packages or
            # user-global NCBI configuration.
            env = build_sra_environment(base_environment=env)

        # Monitor Branch
        if monitor and isinstance(params, dict):
            out_dir_val = params.get("out_dir")
            metadata_val = params.get("metadata")
            if out_dir_val:
                out_dir = Path(str(out_dir_val))
                hb_path = out_dir / ".downloads" / f"amalgkit-{subcommand}-{uuid.uuid4().hex[:8]}.heartbeat.json"

                # Determine watch directory
                watch_dir = out_dir
                if subcommand in {
                    "getfastq",
                    "quant",
                    "metadata",
                    "merge",
                    "wsfilter",
                    "cstmm",
                    "csfilter",
                    "finalize",
                }:
                    watch_dir = out_dir / subcommand

                # Metadata estimation
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

                # Special env vars for getfastq
                if subcommand == "getfastq":
                    # Check tqdm compatibility
                    try:
                        import tqdm

                        if tuple(map(int, tqdm.__version__.split("."))) < (4, 60, 0):
                            env["AMALGKIT_PROGRESS"] = "false"
                    except ImportError, ValueError, AttributeError:
                        pass

                proc_kwargs: dict[str, Any] = {"env": env}
                if cwd:
                    proc_kwargs["cwd"] = cwd

                if log_dir:
                    proc_kwargs["stdout"] = stdout_fh
                    proc_kwargs["stderr"] = stderr_fh
                else:
                    # If monitoring but not logging to file, we might capture or let it stream
                    # For monitoring functions to work, they often expect to manage the process
                    # But if we don't pass capturing, output goes to console
                    pass

                proc = subprocess.Popen(command, **proc_kwargs)

                # Call appropriate monitoring function
                monitor_returncode: int | None = None
                if subcommand == "quant":
                    monitor_returncode = monitor_subprocess_sample_progress(
                        process=proc,
                        watch_dir=watch_dir,
                        heartbeat_path=hb_path,
                        completion_glob="*/abundance.tsv",
                        total_samples=total_runs,
                        heartbeat_interval=heartbeat_interval,
                        show_progress=show_progress,
                        desc=f"amalgkit {subcommand}",
                        timeout_seconds=monitor_timeout,
                    )
                elif subcommand == "metadata":
                    expected = ["metadata.tsv", "metadata_original.tsv", "pivot_selected.tsv", "pivot_qualified.tsv"]
                    monitor_returncode = monitor_subprocess_file_count(
                        process=proc,
                        watch_dir=watch_dir,
                        heartbeat_path=hb_path,
                        expected_files=expected,
                        heartbeat_interval=heartbeat_interval,
                        show_progress=show_progress,
                        desc=f"amalgkit {subcommand}",
                        timeout_seconds=monitor_timeout,
                    )
                elif subcommand in {"merge", "wsfilter", "cstmm", "csfilter", "finalize"}:
                    expected = {
                        "merge": ["merged_abundance.tsv"],
                        "wsfilter": ["metadata.tsv"],
                        "cstmm": ["cstmm.tsv"],
                        "csfilter": ["metadata.tsv"],
                        "finalize": ["metadata.tsv"],
                    }[subcommand]
                    monitor_returncode = monitor_subprocess_file_count(
                        process=proc,
                        watch_dir=watch_dir,
                        heartbeat_path=hb_path,
                        expected_files=expected,
                        heartbeat_interval=heartbeat_interval,
                        show_progress=show_progress,
                        desc=f"amalgkit {subcommand}",
                        timeout_seconds=monitor_timeout,
                    )
                else:
                    # Default: Directory growth (getfastq, integrate, etc)
                    monitor_returncode, _ = monitor_subprocess_directory_growth(
                        process=proc,
                        watch_dir=watch_dir,
                        heartbeat_path=hb_path,
                        total_bytes=total_bytes,
                        heartbeat_interval=heartbeat_interval,
                        show_progress=show_progress,
                        desc=f"amalgkit {subcommand}",
                        timeout_seconds=monitor_timeout,
                    )

                # Wait for process to complete if monitor function didn't (most do)
                if proc.poll() is None:
                    proc.wait()

                # Create a CompletedProcess to return
                # We can't get stdout/stderr if they were streamed/logged, but that's expected
                return subprocess.CompletedProcess(
                    args=command,
                    returncode=(124 if monitor_returncode == 124 else proc.returncode),
                    stdout="" if log_dir else "Output monitored/streamed",
                    stderr="" if log_dir else "Output monitored/streamed",
                )

        # Non-monitored execution
        proc_kwargs = kwargs.copy()
        proc_kwargs["env"] = env  # Use safe env
        if cwd:
            proc_kwargs["cwd"] = cwd

        if log_dir:
            proc_kwargs["stdout"] = stdout_fh
            proc_kwargs["stderr"] = stderr_fh
        else:
            if "capture_output" not in proc_kwargs and "stdout" not in proc_kwargs:
                proc_kwargs["capture_output"] = True
                proc_kwargs["text"] = True

        result: subprocess.CompletedProcess[str] = subprocess.run(command, **proc_kwargs)

        # Log failure info if captured
        if result.returncode != 0:
            logger.warning(f"amalgkit {subcommand} failed with code {result.returncode}")
            if result.stderr:
                # Log last few lines of stderr
                lines = result.stderr.strip().splitlines()
                start_line = max(0, len(lines) - 20)
                for line in lines[start_line:]:
                    logger.warning(f"  [stderr] {line}")
            if result.stdout:
                # Log last few lines of stdout if potentially relevant
                lines = result.stdout.strip().splitlines()
                start_line = max(0, len(lines) - 10)
                for line in lines[start_line:]:
                    logger.info(f"  [stdout] {line}")

        return result

    except Exception as e:
        logger.error(f"Error running amalgkit {subcommand}: {e}")
        # Return failed process if check=False, otherwise re-raise
        if kwargs.get("check"):
            raise
        return subprocess.CompletedProcess(args=command, returncode=1, stdout="", stderr=str(e))
    finally:
        if stdout_fh:
            stdout_fh.close()
        if stderr_fh:
            stderr_fh.close()


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


def select(params: AmalgkitParams | Dict[str, Any] | None = None, **kwargs: Any) -> subprocess.CompletedProcess[str]:
    """Run amalgkit select command.

    Args:
        params: Amalgkit parameters
        **kwargs: Additional arguments

    Returns:
        CompletedProcess instance
    """
    return run_amalgkit("select", params, **kwargs)


def _run_getfastq_via_metainformant(
    params: AmalgkitParams | Dict[str, Any] | None, **kwargs: Any
) -> subprocess.CompletedProcess[str]:
    """Run getfastq using internal ena_downloader (metainformant backend)."""
    from metainformant.rna.retrieval import ena_downloader

    # Extract parameters
    p_dict = params.to_dict() if isinstance(params, AmalgkitParams) else (params or {})
    work_dir = Path(kwargs.get("work_dir") or p_dict.get("work_dir") or ".")

    # Locate metadata
    metadata_path = None
    if p_dict.get("metadata"):
        metadata_path = Path(p_dict["metadata"])
    else:
        metadata_path = _infer_metadata_path_from_out_dir(work_dir)

    if not metadata_path or not metadata_path.exists():
        return subprocess.CompletedProcess(
            args=["internal_ena_downloader"], returncode=1, stderr="Metadata file not found"
        )

    try:
        run_ids = read_sample_ids_from_metadata(metadata_path)
    except Exception as e:
        return subprocess.CompletedProcess(
            args=["internal_ena_downloader"], returncode=1, stderr=f"Error reading metadata: {e}"
        )

    if not run_ids:
        return subprocess.CompletedProcess(
            args=["internal_ena_downloader"],
            returncode=1,
            stderr="No run IDs found in metadata. Expected one of: run, Run, SRA_Run, sample_id, accession",
        )

    getfastq_dir = work_dir / "getfastq"
    getfastq_dir.mkdir(parents=True, exist_ok=True)

    logger.info(f"Downloading {len(run_ids)} samples using metainformant backend to {getfastq_dir}")

    # We pass work_dir as base_out_dir so it creates proper structure
    success_count, fail_count = ena_downloader.download_sra_samples(sra_ids=run_ids, base_out_dir=work_dir)

    if fail_count > 0:
        logger.warning(f"ENA/SRA download completed with {fail_count} failures")
        if success_count == 0:
            return subprocess.CompletedProcess(
                args=["internal_ena_downloader"], returncode=1, stderr="All downloads failed"
            )

    return subprocess.CompletedProcess(
        args=["internal_ena_downloader"],
        returncode=0,
        stdout=f"Downloaded {success_count} samples, {fail_count} failed",
        stderr="",
    )


def getfastq(params: AmalgkitParams | Dict[str, Any] | None = None, **kwargs: Any) -> subprocess.CompletedProcess[str]:
    """Run amalgkit getfastq command with retry logic.

    Args:
        params: Amalgkit parameters
        **kwargs: Additional arguments

    Returns:
        CompletedProcess instance
    """
    # Check backend
    p_dict = params.to_dict() if isinstance(params, AmalgkitParams) else (params or {})
    backend = p_dict.get("backend", "amalgkit")

    if backend == "metainformant":
        return _run_getfastq_via_metainformant(params, **kwargs)

    # Check if parallel execution is requested
    jobs = 1
    if isinstance(params, dict):
        jobs = int(params.get("jobs", 1))
    elif isinstance(params, AmalgkitParams):
        jobs = int(params.extra_params.get("jobs", 1))

    # Retry configuration.  An empty parameter map is a capability/dispatch
    # probe, not a real acquisition request; retrying a deliberately invalid
    # command only wastes the external-test timeout.  Real configured runs
    # retain the resilient three-attempt default, with explicit overrides for
    # callers that need a different recovery policy.
    default_retries = 1 if not p_dict else 3
    default_delay = 0 if not p_dict else 5
    max_retries = int(kwargs.pop("max_retries", default_retries))
    retry_delay = float(kwargs.pop("retry_delay", default_delay))
    if max_retries < 1:
        raise ValueError("max_retries must be at least 1")

    last_result = None
    for attempt in range(1, max_retries + 1):
        if attempt > 1:
            logger.info(f"Retry attempt {attempt}/{max_retries} for getfastq...")

            time.sleep(retry_delay * attempt)  # Linear backoff between attempts

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
                    result.returncode, result.args, output=result.stdout, stderr=result.stderr
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

        except subprocess.CalledProcessError:
            # check=True is a caller-requested hard-failure contract: re-raise
            # instead of downgrading the failure to a retryable CompletedProcess.
            raise
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

    try:
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
    except BaseException:
        # Never leak a chunk temp dir when splitting fails mid-write.
        shutil.rmtree(temp_dir, ignore_errors=True)
        raise

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
    metadata_path: Path | None = None
    if isinstance(params, dict):
        metadata_path = Path(params.get("metadata", ""))
    elif isinstance(params, AmalgkitParams):
        metadata_path = Path(params.extra_params.get("metadata", ""))

    if metadata_path is None or not metadata_path.exists():
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
            single_dict = params.copy()
            single_dict.pop("jobs", None)
            single_params: AmalgkitParams | Dict[str, Any] | None = single_dict
        elif isinstance(params, AmalgkitParams):
            extra_params = params.extra_params.copy()
            extra_params.pop("jobs", None)
            single_params = AmalgkitParams(
                params.work_dir,
                threads=params.threads,
                species_list=list(params.species_list),
                **extra_params,
            )
        else:
            single_params = params

        try:
            # A lone chunk may still live in the split temp dir; reclaim it
            # after the run. The [metadata_path] no-row placeholder must not
            # be deleted.
            return run_amalgkit("getfastq", single_params, **kwargs)
        finally:
            if chunk_paths[0] != metadata_path:
                shutil.rmtree(chunk_paths[0].parent, ignore_errors=True)

    logger.info(f"Parallelizing getfastq with {len(chunk_paths)} workers")

    # Create worker params
    worker_futures = []
    results = []

    with concurrent.futures.ThreadPoolExecutor(max_workers=len(chunk_paths)) as executor:
        for i, chunk_path in enumerate(chunk_paths):
            # Create specific params for this worker
            if isinstance(params, dict):
                worker_dict = params.copy()
                worker_dict["metadata"] = str(chunk_path)
                # Ensure jobs param is removed/reset to prevent recursion if passed down
                worker_dict.pop("jobs", None)
                worker_params: AmalgkitParams | Dict[str, Any] | None = worker_dict
            elif isinstance(params, AmalgkitParams):
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
            from metainformant.rna.engine.workflow_planning import _log_heartbeat

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


def busco(params: AmalgkitParams | Dict[str, Any] | None = None, **kwargs: Any) -> subprocess.CompletedProcess[str]:
    """Run the current Amalgkit BUSCO/compleasm support command."""
    return run_amalgkit("busco", params, **kwargs)


def cstmm(params: AmalgkitParams | Dict[str, Any] | None = None, **kwargs: Any) -> subprocess.CompletedProcess[str]:
    """Run amalgkit cstmm command.

    Args:
        params: Amalgkit parameters
        **kwargs: Additional arguments

    Returns:
        CompletedProcess instance
    """
    return run_amalgkit("cstmm", params, **kwargs)


def finalize(params: AmalgkitParams | Dict[str, Any] | None = None, **kwargs: Any) -> subprocess.CompletedProcess[str]:
    """Run the current amalgkit finalize command directly."""
    return run_amalgkit("finalize", params, **kwargs)


def wsfilter(params: AmalgkitParams | Dict[str, Any] | None = None, **kwargs: Any) -> subprocess.CompletedProcess[str]:
    """Run the current Amalgkit within-species filter command.

    Args:
        params: Amalgkit parameters
        **kwargs: Additional arguments

    Returns:
        CompletedProcess instance
    """
    return run_amalgkit("wsfilter", params, **kwargs)


def csfilter(params: AmalgkitParams | Dict[str, Any] | None = None, **kwargs: Any) -> subprocess.CompletedProcess[str]:
    """Run the current Amalgkit cross-species filter command.

    Args:
        params: Amalgkit parameters
        **kwargs: Additional arguments

    Returns:
        CompletedProcess instance
    """
    return run_amalgkit("csfilter", params, **kwargs)


def sanity(params: AmalgkitParams | Dict[str, Any] | None = None, **kwargs: Any) -> subprocess.CompletedProcess[str]:
    """Run amalgkit sanity command.

    Args:
        params: Amalgkit parameters
        **kwargs: Additional arguments

    Returns:
        CompletedProcess instance
    """
    return run_amalgkit("sanity", params, **kwargs)


def rerun(params: AmalgkitParams | Dict[str, Any] | None = None, **kwargs: Any) -> subprocess.CompletedProcess[str]:
    """Run the current Amalgkit failed-target recovery command."""
    return run_amalgkit("rerun", params, **kwargs)


def dataset(params: AmalgkitParams | Dict[str, Any] | None = None, **kwargs: Any) -> subprocess.CompletedProcess[str]:
    """Run the current Amalgkit workspace/rule-set export command."""
    return run_amalgkit("dataset", params, **kwargs)
