from __future__ import annotations

import os
import shutil
import subprocess
import sys
import threading
import time
from collections.abc import Mapping
from datetime import datetime
from pathlib import Path
from typing import Any

from ..core.logging import get_logger

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
    "ncbi", "aws", "gcp", "cleanup_raw", "build_index", "resolve_names", "mark_redundant_biosamples"
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
    # Check standard PATH first
    amalgkit_path = shutil.which("amalgkit")
    
    # Also check user local bin (common for --user installs)
    if amalgkit_path is None:
        user_bin = Path.home() / ".local" / "bin" / "amalgkit"
        if user_bin.exists() and user_bin.is_file():
            amalgkit_path = str(user_bin)
    
    if amalgkit_path is None:
        return False, "amalgkit not found on PATH"

    try:
        # Use the found path or default to "amalgkit" (which will use PATH)
        cmd = [amalgkit_path if amalgkit_path else "amalgkit", "-h"]
        proc = subprocess.run(cmd, capture_output=True, text=True, check=False)
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
    
    # config_dir style (for select step) - ensure absolute path
    if "config_dir" in safe_params and isinstance(safe_params["config_dir"], (str, Path)):
        cfg_path = Path(str(safe_params["config_dir"])).expanduser()
        try:
            cfg_path = cfg_path.resolve()
        except Exception:
            # If resolution fails, keep original but log
            get_logger(__name__).debug(f"Could not resolve config_dir path: {safe_params['config_dir']}")
            cfg_path = Path(str(safe_params["config_dir"]))
        # Ensure directory exists (amalgkit expects it to exist)
        if not cfg_path.exists():
            try:
                cfg_path.mkdir(parents=True, exist_ok=True)
            except Exception:
                # If we can't create it, log but continue - amalgkit will handle the error
                get_logger(__name__).warning(f"Could not create config_dir: {cfg_path}")
        safe_params["config_dir"] = str(cfg_path)
    
    # merge 'out' file path
    if subcommand == "merge" and isinstance(safe_params.get("out"), (str, Path)):
        out_path = Path(str(safe_params["out"]))
        out_path = out_path.expanduser()
        try:
            out_path = out_path.resolve()
        except Exception as e:
            # Log error but continue - path resolution issues shouldn't block execution
            get_logger(__name__).debug(f"Could not resolve output path: {e}")
        out_path.parent.mkdir(parents=True, exist_ok=True)
        safe_params["out"] = str(out_path)

    cmd = build_amalgkit_command(subcommand, safe_params)
    run_env = os.environ.copy()
    if env:
        run_env.update(env)

    # For getfastq, set working directory to out_dir so prefetch downloads to correct location
    # This prevents prefetch from using --output-directory ./ which would download to repo root
    effective_work_dir = None
    if work_dir is not None:
        Path(work_dir).mkdir(parents=True, exist_ok=True)
        effective_work_dir = str(work_dir)
    elif subcommand == "getfastq" and "out_dir" in safe_params:
        # For getfastq, use out_dir as working directory to ensure prefetch downloads there
        out_dir_path = Path(str(safe_params["out_dir"])).expanduser().resolve()
        out_dir_path.mkdir(parents=True, exist_ok=True)
        effective_work_dir = str(out_dir_path)

    # Streaming mode (for long-running steps) when log_dir is set
    # Enable by default for long-running steps, or if MI_STREAM_AMALGKIT_LOGS=1
    long_running_steps = {"getfastq", "quant", "merge"}
    stream = (
        log_dir is not None
        and (
            os.getenv("MI_STREAM_AMALGKIT_LOGS") == "1"
            or subcommand in long_running_steps
        )
    )
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
            get_logger(__name__).debug(f"Could not resolve output path: {e}")
        log_path = Path(log_dir)
        log_path.mkdir(parents=True, exist_ok=True)
        stdout_file = log_path / f"{ts}.{base}.stdout.log"
        stderr_file = log_path / f"{ts}.{base}.stderr.log"

        proc = subprocess.Popen(
            cmd,
            cwd=effective_work_dir,
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
                        get_logger(__name__).debug(f"Streaming output error: {e}")

        threads: list[threading.Thread] = []
        if proc.stdout is not None:
            t_out = threading.Thread(target=_tee, args=(proc.stdout, sys.stdout, stdout_file), daemon=True)
            t_out.start()
            threads.append(t_out)
        if proc.stderr is not None:
            t_err = threading.Thread(target=_tee, args=(proc.stderr, sys.stderr, stderr_file), daemon=True)
            t_err.start()
            threads.append(t_err)

        # Heartbeat reporter for long-running quiet periods with progress tracking
        stop_heartbeat = threading.Event()
        
        # Extract output directory for progress monitoring (for getfastq/quant steps)
        progress_monitor_dir: Path | None = None
        if subcommand == "getfastq":
            # Try to find out_dir from params
            out_dir_val = safe_params.get("out_dir") or safe_params.get("out-dir")
            if out_dir_val:
                progress_monitor_dir = Path(str(out_dir_val)).expanduser()
                try:
                    progress_monitor_dir = progress_monitor_dir.resolve()
                except Exception:
                    pass
                # Check for getfastq subdirectory
                getfastq_subdir = progress_monitor_dir / "getfastq"
                if getfastq_subdir.exists():
                    progress_monitor_dir = getfastq_subdir
        elif subcommand == "quant":
            # For quant, monitor the quant output directory
            quant_dir_val = safe_params.get("out_dir") or safe_params.get("out-dir") or safe_params.get("quant_dir") or safe_params.get("quant-dir")
            if quant_dir_val:
                progress_monitor_dir = Path(str(quant_dir_val)).expanduser()
                try:
                    progress_monitor_dir = progress_monitor_dir.resolve()
                except Exception:
                    pass

        def _heartbeat():
            """Emit periodic heartbeat messages with progress information."""
            heartbeat_count = 0
            last_size = 0
            last_check_time = time.time()
            
            while not stop_heartbeat.is_set():
                # Wait 30 seconds or until event is set
                event_set = stop_heartbeat.wait(30.0)
                if event_set:
                    # Event was set, process finished
                    break
                # Timeout occurred - check if process is still running
                if proc.poll() is None:
                    heartbeat_count += 1
                    try:
                        now_hhmm = datetime.utcnow().strftime("%H:%M:%S")
                        now_time = time.time()
                        
                        # Build progress message
                        progress_info = ""
                        if progress_monitor_dir and progress_monitor_dir.exists():
                            try:
                                # Efficiently calculate total size and find active samples
                                total_size = 0
                                file_count = 0
                                recent_samples: list[tuple[str, float]] = []  # (sample_id, mtime)
                                
                                # For getfastq, check sample directories directly (more efficient)
                                if subcommand == "getfastq":
                                    # Check sample directories (SRR*)
                                    for sample_dir in progress_monitor_dir.iterdir():
                                        if sample_dir.is_dir() and sample_dir.name.startswith("SRR"):
                                            sample_mtime = 0
                                            sample_size = 0
                                            sample_files = 0
                                            
                                            # Check files in this sample directory
                                            for item in sample_dir.rglob("*"):
                                                if item.is_file():
                                                    try:
                                                        stat = item.stat()
                                                        total_size += stat.st_size
                                                        file_count += 1
                                                        sample_size += stat.st_size
                                                        sample_files += 1
                                                        # Track most recent modification
                                                        if stat.st_mtime > sample_mtime:
                                                            sample_mtime = stat.st_mtime
                                                    except (OSError, FileNotFoundError):
                                                        pass
                                            
                                            # If sample was modified recently (last 5 min), track it
                                            if sample_mtime > 0 and (now_time - sample_mtime) < 300:
                                                recent_samples.append((sample_dir.name, sample_mtime))
                                else:
                                    # For other steps (quant, etc.), scan all files
                                    for item in progress_monitor_dir.rglob("*"):
                                        if item.is_file():
                                            try:
                                                stat = item.stat()
                                                total_size += stat.st_size
                                                file_count += 1
                                            except (OSError, FileNotFoundError):
                                                pass
                                
                                # Format size
                                if total_size > 1024 * 1024 * 1024:  # GB
                                    size_str = f"{total_size / (1024**3):.2f}GB"
                                elif total_size > 1024 * 1024:  # MB
                                    size_str = f"{total_size / (1024**2):.1f}MB"
                                elif total_size > 1024:  # KB
                                    size_str = f"{total_size / 1024:.1f}KB"
                                else:
                                    size_str = f"{total_size}B"
                                
                                # Calculate download rate
                                elapsed = now_time - last_check_time
                                size_delta = total_size - last_size
                                rate_str = ""
                                if elapsed > 0 and size_delta > 0:
                                    rate_mbps = (size_delta / (1024 * 1024)) / elapsed
                                    rate_str = f" @ {rate_mbps:.2f}MB/s"
                                
                                # Find most recently active sample
                                current_sample = ""
                                if recent_samples:
                                    recent_samples.sort(key=lambda x: x[1], reverse=True)
                                    current_sample = f" | Current: {recent_samples[0][0]}"
                                
                                progress_info = f" | {size_str} ({file_count} files){rate_str}{current_sample}"
                                last_size = total_size
                                last_check_time = now_time
                            except Exception as e:
                                # If monitoring fails, continue without progress info
                                get_logger(__name__).debug(f"Progress monitoring error: {e}")
                        
                        sys.stdout.write(
                            f"[{now_hhmm}] still running step '{base}' (pid={proc.pid}, "
                            f"heartbeat #{heartbeat_count}){progress_info}...\n"
                        )
                        sys.stdout.flush()
                    except Exception as e:
                        # Log error but continue - streaming issues shouldn't block execution
                        get_logger(__name__).debug(f"Streaming output error: {e}")
                else:
                    # Process finished, exit loop
                    break

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
            get_logger(__name__).debug(f"Could not resolve output path: {e}")

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
        cwd=effective_work_dir,
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
