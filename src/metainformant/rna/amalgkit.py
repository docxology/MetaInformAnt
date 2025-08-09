from __future__ import annotations

# Thin, typed utilities to invoke the external `amalgkit` RNA-seq toolkit.
#
# Provides:
# - Parameter normalization: turn a mapping of Python values into CLI flags
# - Command builders for each amalgkit subcommand
# - A single safe runner that handles working directories and optional logging
#
# Repository policy:
# - Do not persist outside `output/` by default (callers pass work_dir/log_dir)
# - Return `subprocess.CompletedProcess[str]` for testability
# - Avoid side effects beyond the requested execution and optional log writes


import os
import shutil
import subprocess
from collections.abc import Mapping
from datetime import datetime
from pathlib import Path
from typing import Any


AmalgkitParams = Mapping[str, Any]


def _normalize_key_to_flag(key: str) -> str:
    """Normalize a mapping key to a CLI flag name.

    Rules:
    - Convert underscores to hyphens
    - Ensure leading "--"
    """
    key = key.strip().replace("_", "-")
    if not key.startswith("--"):
        key = f"--{key}"
    return key


def _ensure_str(value: Any) -> str:
    """Stringify values, preserving paths as filesystem strings."""
    if isinstance(value, Path):
        return str(value)
    return str(value)


def build_cli_args(params: AmalgkitParams | None) -> list[str]:
    """Convert a params mapping into a flat list of CLI args for `amalgkit`.

    Rules:
    - None values are skipped
    - bool True adds a flag (e.g., {dry_run: True} → ["--dry-run"]) ; False is skipped
    - list/tuple produces repeated pairs (e.g., {species: [a,b]} → [--species a --species b])
    - Path values are stringified
    - other scalars are appended as flag + value
    """
    if not params:
        return []

    args: list[str] = []
    for raw_key, value in params.items():
        if value is None:
            continue

        flag = _normalize_key_to_flag(raw_key)

        if isinstance(value, bool):
            if value:
                args.append(flag)
            continue

        if isinstance(value, list | tuple):
            for item in value:
                args.append(flag)
                args.append(_ensure_str(item))
            continue

        args.append(flag)
        args.append(_ensure_str(value))

    return args


def build_amalgkit_command(subcommand: str, params: AmalgkitParams | None = None) -> list[str]:
    """Prepare the `amalgkit` CLI command as a token list.

    Example: build_amalgkit_command("metadata", {"threads": 8}) →
        ["amalgkit", "metadata", "--threads", "8"]
    """
    return ["amalgkit", subcommand] + build_cli_args(params)


def check_cli_available() -> tuple[bool, str]:
    """Check if `amalgkit` is available on PATH.

    Returns a tuple: (is_available, help_or_version_text_or_error_message)
    """
    if shutil.which("amalgkit") is None:
        return False, "amalgkit not found on PATH"

    try:
        proc = subprocess.run(["amalgkit", "-h"], capture_output=True, text=True, check=False)
        output = proc.stdout or proc.stderr
        return proc.returncode == 0, output
    except Exception as exc:  # pragma: no cover - defensive
        return False, f"error invoking amalgkit: {exc}"


def run_amalgkit(
    subcommand: str,
    params: AmalgkitParams | None = None,
    *,
    work_dir: str | Path | None = None,
    env: Mapping[str, str] | None = None,
    check: bool = False,
    capture_output: bool = True,
    log_dir: str | Path | None = None,
    step_name: str | None = None,
) -> subprocess.CompletedProcess[str]:
    """Execute an `amalgkit` subcommand with optional logging.

    Parameters
    - subcommand: e.g., "metadata", "integrate", "quant", ...
    - params: mapping of flags/values following `build_cli_args` rules
    - work_dir: working directory to run in (created if missing)
    - env: additional environment variables to merge with current env
    - check: if True, raise CalledProcessError on non-zero exit
    - capture_output: capture stdout/stderr as text
    - log_dir: if provided, write timestamped stdout/stderr logs per step
    - step_name: optional label to prefix log filenames
    """
    cmd = build_amalgkit_command(subcommand, params)
    run_env = os.environ.copy()
    if env:
        run_env.update(env)

    # Ensure working directory exists if provided
    if work_dir is not None:
        Path(work_dir).mkdir(parents=True, exist_ok=True)

    result = subprocess.run(
        cmd,
        cwd=str(work_dir) if work_dir is not None else None,
        env=run_env,
        capture_output=capture_output,
        text=True,
        check=check,
    )
    # Optionally write logs per step
    if log_dir is not None:
        log_path = Path(log_dir)
        log_path.mkdir(parents=True, exist_ok=True)
        ts = datetime.utcnow().strftime("%Y%m%dT%H%M%SZ")
        base = step_name or subcommand
        (log_path / f"{ts}.{base}.stdout.log").write_text(result.stdout or "")
        (log_path / f"{ts}.{base}.stderr.log").write_text(result.stderr or "")
    return result


"""Convenience wrappers for each `amalgkit` subcommand."""


def metadata(params: AmalgkitParams | None = None, **kwargs: Any) -> subprocess.CompletedProcess[str]:
    """Run `amalgkit metadata` (fetch SRA/ENA metadata)."""
    return run_amalgkit("metadata", params, **kwargs)


def integrate(params: AmalgkitParams | None = None, **kwargs: Any) -> subprocess.CompletedProcess[str]:
    """Run `amalgkit integrate` (combine local FASTQs with metadata)."""
    return run_amalgkit("integrate", params, **kwargs)


def config(params: AmalgkitParams | None = None, **kwargs: Any) -> subprocess.CompletedProcess[str]:
    """Run `amalgkit config` (generate tool configs for downstream steps)."""
    return run_amalgkit("config", params, **kwargs)


def select(params: AmalgkitParams | None = None, **kwargs: Any) -> subprocess.CompletedProcess[str]:
    """Run `amalgkit select` (filter samples based on criteria)."""
    return run_amalgkit("select", params, **kwargs)


def getfastq(params: AmalgkitParams | None = None, **kwargs: Any) -> subprocess.CompletedProcess[str]:
    """Run `amalgkit getfastq` (download raw FASTQ files)."""
    return run_amalgkit("getfastq", params, **kwargs)


def quant(params: AmalgkitParams | None = None, **kwargs: Any) -> subprocess.CompletedProcess[str]:
    """Run `amalgkit quant` (quantify transcript abundances)."""
    return run_amalgkit("quant", params, **kwargs)


def merge(params: AmalgkitParams | None = None, **kwargs: Any) -> subprocess.CompletedProcess[str]:
    """Run `amalgkit merge` (combine per-sample quantifications)."""
    return run_amalgkit("merge", params, **kwargs)


def cstmm(params: AmalgkitParams | None = None, **kwargs: Any) -> subprocess.CompletedProcess[str]:
    """Run `amalgkit cstmm` (cross-species TMM normalization)."""
    return run_amalgkit("cstmm", params, **kwargs)


def curate(params: AmalgkitParams | None = None, **kwargs: Any) -> subprocess.CompletedProcess[str]:
    """Run `amalgkit curate` (detect outliers and unwanted biases)."""
    return run_amalgkit("curate", params, **kwargs)


def csca(params: AmalgkitParams | None = None, **kwargs: Any) -> subprocess.CompletedProcess[str]:
    """Run `amalgkit csca` (cross-species correlation analysis and plots)."""
    return run_amalgkit("csca", params, **kwargs)


def sanity(params: AmalgkitParams | None = None, **kwargs: Any) -> subprocess.CompletedProcess[str]:
    """Run `amalgkit sanity` (final integrity checks)."""
    return run_amalgkit("sanity", params, **kwargs)


__all__ = [
    "AmalgkitParams",
    "build_cli_args",
    "build_amalgkit_command",
    "check_cli_available",
    "run_amalgkit",
    # subcommands
    "metadata",
    "integrate",
    "config",
    "select",
    "getfastq",
    "quant",
    "merge",
    "cstmm",
    "curate",
    "csca",
    "sanity",
]


