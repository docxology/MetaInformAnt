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
import sys
from collections.abc import Mapping
from datetime import datetime
from pathlib import Path
from typing import Any


AmalgkitParams = Mapping[str, Any]


def _normalize_key_to_flag(key: str, *, for_cli: bool) -> str:
    """Normalize a mapping key to a CLI flag token.

    - When for_cli is True: keep snake_case (e.g., --out_dir) for compatibility
      with the external amalgkit CLI which commonly uses underscores.
    - When for_cli is False: prefer kebab-case (e.g., --out-dir) for display/tests.
    """
    name = key.strip()
    if not for_cli:
        name = name.replace("_", "-")
    if not name.startswith("--"):
        name = f"--{name}"
    return name


_KEY_ALIASES = {
    # Common hyphenated -> underscored variants
    "out-dir": "out_dir",
    "genome-dir": "genome_dir",
    "fastq-dir": "fastq_dir",
    "quant-dir": "quant_dir",
    "search-string": "search_string",
    "entrez-email": "entrez_email",
    "redo": "redo",  # passthrough
}


def _normalize_param_keys(params: AmalgkitParams | None) -> dict[str, Any]:
    if not params:
        return {}
    out: dict[str, Any] = {}
    for k, v in params.items():
        kk = str(k).strip()
        if kk in _KEY_ALIASES:
            kk = _KEY_ALIASES[kk]
        out[kk] = v
    return out


def _ensure_str(value: Any) -> str:
    """Stringify values, preserving paths as filesystem strings."""
    if isinstance(value, Path):
        return str(value)
    return str(value)


_BOOL_VALUE_FLAGS: set[str] = {"redo"}


def build_cli_args(params: AmalgkitParams | None, *, for_cli: bool = False) -> list[str]:
    """Convert a params mapping into a flat list of CLI args for `amalgkit`.

    Rules:
    - None values are skipped
    - For most flags: bool True adds a flag (e.g., {dry_run: True} → ["--dry-run"]) ; False is skipped
    - For certain flags that require explicit values (e.g., {redo: True|False}),
      emit yes/no as a value pair (e.g., ["--redo", "yes"]).
    - list/tuple produces repeated pairs (e.g., {species: [a,b]} → [--species a --species b])
    - Path values are stringified
    - other scalars are appended as flag + value
    """
    if not params:
        return []

    params = _normalize_param_keys(params)

    args: list[str] = []
    for raw_key, value in params.items():
        if value is None:
            continue

        key_name = str(raw_key)
        flag = _normalize_key_to_flag(key_name, for_cli=for_cli)

        if isinstance(value, bool):
            # Some flags (e.g., --redo) expect an explicit yes/no argument
            if key_name in _BOOL_VALUE_FLAGS:
                args.append(flag)
                args.append("yes" if value else "no")
            elif value:
                args.append(flag)
            continue

        if isinstance(value, (list, tuple)):
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
    # Use CLI-compatible flag style (snake_case) when building the real command
    return ["amalgkit", subcommand] + build_cli_args(params, for_cli=True)


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


def ensure_cli_available(*, auto_install: bool = False) -> tuple[bool, str, dict | None]:
    """Ensure `amalgkit` CLI is available; optionally attempt auto-install.

    Returns (ok, message, install_record_dict_or_none).
    install_record contains keys: {"attempted": bool, "return_code": int, "stdout": str, "stderr": str}
    """
    ok, msg = check_cli_available()
    if ok or not auto_install:
        return ok, msg, None

    # Attempt installation via pip
    cmd = [sys.executable, "-m", "pip", "install", "--no-input", "--no-warn-script-location", "git+https://github.com/kfuku52/amalgkit"]
    try:
        proc = subprocess.run(cmd, capture_output=True, text=True, check=False)
        ok2, msg2 = check_cli_available()
        install_rec = {
            "attempted": True,
            "return_code": proc.returncode,
            "stdout": proc.stdout or "",
            "stderr": proc.stderr or "",
            "command": " ".join(cmd),
        }
        return ok2, (msg2 if ok2 else msg), install_rec
    except Exception as exc:  # pragma: no cover - defensive
        return False, f"auto-install failed: {exc}", {"attempted": True, "return_code": -1, "stdout": "", "stderr": str(exc), "command": " ".join(cmd)}


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

    # Streaming mode (for long-running steps) when MI_STREAM_AMALGKIT_LOGS=1 and log_dir is set
    stream = os.getenv("MI_STREAM_AMALGKIT_LOGS") == "1" and log_dir is not None
    ts = datetime.utcnow().strftime("%Y%m%dT%H%M%SZ")
    base = step_name or subcommand

    if stream:
        log_path = Path(log_dir)
        log_path.mkdir(parents=True, exist_ok=True)
        stdout_file = log_path / f"{ts}.{base}.stdout.log"
        stderr_file = log_path / f"{ts}.{base}.stderr.log"
        with open(stdout_file, "w", encoding="utf-8") as out_fh, open(stderr_file, "w", encoding="utf-8") as err_fh:
            proc = subprocess.Popen(
                cmd,
                cwd=str(work_dir) if work_dir is not None else None,
                env=run_env,
                stdout=out_fh,
                stderr=err_fh,
                text=True,
            )
            rc = proc.wait()
        # Build a CompletedProcess with empty captured text (logs are in files)
        result = subprocess.CompletedProcess(cmd, rc, stdout="", stderr="")
        if check and rc != 0:
            raise subprocess.CalledProcessError(rc, cmd)
        return result

    # Default: capture then write logs
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
    """Run `amalgkit config` (generate tool configs for downstream tools)."""
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
    "ensure_cli_available",
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


