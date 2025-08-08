from __future__ import annotations

import os
import shutil
import subprocess
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence


AmalgkitParams = Mapping[str, Any]


def _normalize_key_to_flag(key: str) -> str:
    key = key.strip().replace("_", "-")
    if not key.startswith("--"):
        key = f"--{key}"
    return key


def _ensure_str(value: Any) -> str:
    if isinstance(value, Path):
        return str(value)
    return str(value)


def build_cli_args(params: AmalgkitParams | None) -> list[str]:
    """Convert a params dict into a flat list of CLI args for amalgkit.

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

        if isinstance(value, (list, tuple)):
            for item in value:
                args.append(flag)
                args.append(_ensure_str(item))
            continue

        args.append(flag)
        args.append(_ensure_str(value))

    return args


def build_amalgkit_command(subcommand: str, params: AmalgkitParams | None = None) -> list[str]:
    """Prepare the amalgkit CLI command as a token list."""
    return ["amalgkit", subcommand] + build_cli_args(params)


def check_cli_available() -> tuple[bool, str]:
    """Return (available, help_or_version_text)."""
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
) -> subprocess.CompletedProcess[str]:
    """Execute an amalgkit subcommand.

    Parameters
    - subcommand: e.g., "metadata", "integrate", "quant", ...
    - params: dict of flags/values following build_cli_args rules
    - work_dir: working directory to run in
    - env: additional environment variables to merge with current env
    - check: if True, raise CalledProcessError on non-zero exit
    - capture_output: capture stdout/stderr as text
    """
    cmd = build_amalgkit_command(subcommand, params)
    run_env = os.environ.copy()
    if env:
        run_env.update(env)

    return subprocess.run(
        cmd,
        cwd=str(work_dir) if work_dir is not None else None,
        env=run_env,
        capture_output=capture_output,
        text=True,
        check=check,
    )


# Convenience wrappers for each documented amalgkit module command
def metadata(params: AmalgkitParams | None = None, **kwargs: Any) -> subprocess.CompletedProcess[str]:
    return run_amalgkit("metadata", params, **kwargs)


def integrate(params: AmalgkitParams | None = None, **kwargs: Any) -> subprocess.CompletedProcess[str]:
    return run_amalgkit("integrate", params, **kwargs)


def config(params: AmalgkitParams | None = None, **kwargs: Any) -> subprocess.CompletedProcess[str]:
    return run_amalgkit("config", params, **kwargs)


def select(params: AmalgkitParams | None = None, **kwargs: Any) -> subprocess.CompletedProcess[str]:
    return run_amalgkit("select", params, **kwargs)


def getfastq(params: AmalgkitParams | None = None, **kwargs: Any) -> subprocess.CompletedProcess[str]:
    return run_amalgkit("getfastq", params, **kwargs)


def quant(params: AmalgkitParams | None = None, **kwargs: Any) -> subprocess.CompletedProcess[str]:
    return run_amalgkit("quant", params, **kwargs)


def merge(params: AmalgkitParams | None = None, **kwargs: Any) -> subprocess.CompletedProcess[str]:
    return run_amalgkit("merge", params, **kwargs)


def cstmm(params: AmalgkitParams | None = None, **kwargs: Any) -> subprocess.CompletedProcess[str]:
    return run_amalgkit("cstmm", params, **kwargs)


def curate(params: AmalgkitParams | None = None, **kwargs: Any) -> subprocess.CompletedProcess[str]:
    return run_amalgkit("curate", params, **kwargs)


def csca(params: AmalgkitParams | None = None, **kwargs: Any) -> subprocess.CompletedProcess[str]:
    return run_amalgkit("csca", params, **kwargs)


def sanity(params: AmalgkitParams | None = None, **kwargs: Any) -> subprocess.CompletedProcess[str]:
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


