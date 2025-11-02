from __future__ import annotations

import os
import shutil
import subprocess
import sys
import threading
from collections.abc import Mapping
from datetime import datetime
from pathlib import Path
from typing import Any

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


_BOOL_VALUE_FLAGS: set[str] = {
    "redo", "pfd", "fastp", "remove_sra", "remove_tmp", "pfd_print", "fastp_print",
    "ncbi", "aws", "gcp", "cleanup_raw", "build_index"
}


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
    cmd = [
        sys.executable,
        "-m",
        "pip",
        "install",
        "--no-input",
        "--no-warn-script-location",
        "git+https://github.com/kfuku52/amalgkit",
    ]
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
        return (
            False,
            f"auto-install failed: {exc}",
            {"attempted": True, "return_code": -1, "stdout": "", "stderr": str(exc), "command": " ".join(cmd)},
        )


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
    # Prepare safe, mutable params and ensure important directories exist
    safe_params: dict[str, Any] = {}
    if params:
        safe_params.update({str(k): v for k, v in params.items()})

    # Normalize and pre-create output directories used by amalgkit
    def _mkdir_parent(path_str: str) -> str:
        p = Path(path_str).expanduser()
        # If relative, resolve relative to repository root
        try:
            p = p.resolve()
        except Exception:
            p = p
        p.parent.mkdir(parents=True, exist_ok=True)
        # Some subcommands expect the directory itself to exist (metadata)
        if subcommand in {"metadata", "getfastq", "quant", "cstmm", "curate", "csca"}:
            p.mkdir(parents=True, exist_ok=True)
        return str(p)

    # out_dir style
    for key in ("out_dir", "out-dir"):
        if key in safe_params and isinstance(safe_params[key], (str, Path)):
            safe_params["out_dir"] = _mkdir_parent(str(safe_params[key]))
            break
    # merge 'out' file path
    if subcommand == "merge" and isinstance(safe_params.get("out"), (str, Path)):
        out_path = Path(str(safe_params["out"]))
        out_path = out_path.expanduser()
        try:
            out_path = out_path.resolve()
        except Exception as e:
            # Log error but continue - path resolution issues shouldn't block execution
            import logging
            logging.getLogger(__name__).debug(f"Could not resolve output path: {e}")
        out_path.parent.mkdir(parents=True, exist_ok=True)
        safe_params["out"] = str(out_path)

    cmd = build_amalgkit_command(subcommand, safe_params)
    run_env = os.environ.copy()
    if env:
        run_env.update(env)

    # Ensure working directory exists if provided
    # NOTE: work_dir is the CWD where amalgkit runs, NOT used in command args
    if work_dir is not None:
        Path(work_dir).mkdir(parents=True, exist_ok=True)
        # CRITICAL: Do NOT change to work_dir - let amalgkit use absolute paths in command
        work_dir = None  # Run from current directory to avoid path resolution issues

    # Streaming mode (for long-running steps) when MI_STREAM_AMALGKIT_LOGS=1 and log_dir is set
    stream = os.getenv("MI_STREAM_AMALGKIT_LOGS") == "1" and log_dir is not None
    ts = datetime.utcnow().strftime("%Y%m%dT%H%M%SZ")
    base = step_name or subcommand

    if stream:
        # Announce start
        try:
            start_hhmm = datetime.utcnow().strftime("%H:%M:%S")
            sys.stdout.write(f"[{start_hhmm}] starting step '{base}' -> {' '.join(cmd)}\n")
            sys.stdout.flush()
        except Exception as e:
            # Log error but continue - path resolution issues shouldn't block execution
            import logging
            logging.getLogger(__name__).debug(f"Could not resolve output path: {e}")
        log_path = Path(log_dir)
        log_path.mkdir(parents=True, exist_ok=True)
        stdout_file = log_path / f"{ts}.{base}.stdout.log"
        stderr_file = log_path / f"{ts}.{base}.stderr.log"

        proc = subprocess.Popen(
            cmd,
            cwd=str(work_dir) if work_dir is not None else None,
            env=run_env,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            bufsize=1,
        )

        def _tee(stream_in, stream_out, file_path: Path):
            with open(file_path, "w", encoding="utf-8") as fh:
                for line in iter(stream_in.readline, ""):
                    fh.write(line)
                    fh.flush()
                    try:
                        stream_out.write(line)
                        stream_out.flush()
                    except Exception as e:
                        # Log error but continue - streaming issues shouldn't block execution
                        import logging
                        logging.getLogger(__name__).debug(f"Streaming output error: {e}")

        threads: list[threading.Thread] = []
        if proc.stdout is not None:
            t_out = threading.Thread(target=_tee, args=(proc.stdout, sys.stdout, stdout_file), daemon=True)
            t_out.start()
            threads.append(t_out)
        if proc.stderr is not None:
            t_err = threading.Thread(target=_tee, args=(proc.stderr, sys.stderr, stderr_file), daemon=True)
            t_err.start()
            threads.append(t_err)

        # Heartbeat reporter for long-running quiet periods
        stop_heartbeat = threading.Event()

        def _heartbeat():
            while not stop_heartbeat.is_set():
                # Every 30s, emit a heartbeat if still running
                stop_heartbeat.wait(30.0)
                if stop_heartbeat.is_set():
                    break
                if proc.poll() is None:
                    try:
                        now_hhmm = datetime.utcnow().strftime("%H:%M:%S")
                        sys.stdout.write(f"[{now_hhmm}] still running step '{base}' (pid={proc.pid})...\n")
                        sys.stdout.flush()
                    except Exception as e:
                        # Log error but continue - streaming issues shouldn't block execution
                        import logging
                        logging.getLogger(__name__).debug(f"Streaming output error: {e}")

        hb = threading.Thread(target=_heartbeat, daemon=True)
        hb.start()

        rc = proc.wait()
        stop_heartbeat.set()
        for t in threads:
            t.join(timeout=1.0)
        try:
            end_hhmm = datetime.utcnow().strftime("%H:%M:%S")
            sys.stdout.write(f"[{end_hhmm}] finished step '{base}' with code {rc}\n")
            sys.stdout.flush()
        except Exception as e:
            # Log error but continue - path resolution issues shouldn't block execution
            import logging
            logging.getLogger(__name__).debug(f"Could not resolve output path: {e}")

        # Build a CompletedProcess with empty captured text (logs are in files and console)
        result = subprocess.CompletedProcess(cmd, rc, stdout="", stderr="")
        if check and rc != 0:
            raise subprocess.CalledProcessError(rc, cmd)
        return result

    # Default: capture then write logs
    try:
        start_hhmm = datetime.utcnow().strftime("%H:%M:%S")
        sys.stdout.write(f"[{start_hhmm}] starting step '{base}' -> {' '.join(cmd)}\n")
        sys.stdout.flush()
    except Exception:
        pass
    result = subprocess.run(
        cmd,
        cwd=str(work_dir) if work_dir is not None else None,
        env=run_env,
        capture_output=capture_output,
        text=True,
        check=check,
    )
    try:
        end_hhmm = datetime.utcnow().strftime("%H:%M:%S")
        sys.stdout.write(f"[{end_hhmm}] finished step '{base}' with code {result.returncode}\n")
        sys.stdout.flush()
    except Exception:
        pass
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
