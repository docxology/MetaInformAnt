"""Amalgkit CLI integration for RNA-seq workflow orchestration.

This module provides Python wrappers around the amalgkit command-line interface
for comprehensive RNA-seq analysis workflows.
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, Union

from metainformant.core import logging

logger = logging.get_logger(__name__)


class AmalgkitParams:
    """Parameters for amalgkit commands."""

    def __init__(self,
                 work_dir: Union[str, Path],
                 threads: int = 8,
                 species_list: Optional[List[str]] = None,
                 **kwargs):
        """Initialize amalgkit parameters.

        Args:
            work_dir: Working directory for amalgkit
            threads: Number of threads to use
            species_list: List of species to process
            **kwargs: Additional parameters
        """
        self.work_dir = Path(work_dir)
        self.threads = threads
        self.species_list = species_list or []
        self.extra_params = kwargs

    def to_dict(self) -> Dict[str, Any]:
        """Convert parameters to dictionary."""
        return {
            'work_dir': str(self.work_dir),
            'threads': self.threads,
            'species_list': self.species_list,
            **self.extra_params
        }


def build_cli_args(params: AmalgkitParams | None, *, for_cli: bool = False) -> List[str]:
    """Build command-line arguments for amalgkit.

    Args:
        params: Amalgkit parameters
        for_cli: Whether to format for CLI usage

    Returns:
        List of command-line arguments
    """
    args = []

    if params:
        args.extend(['--work-dir', str(params.work_dir)])
        args.extend(['--threads', str(params.threads)])

        if params.species_list:
            for species in params.species_list:
                args.extend(['--species', species])

        # Add extra parameters
        for key, value in params.extra_params.items():
            if isinstance(value, bool):
                if value:
                    args.append(f'--{key.replace("_", "-")}')
            elif isinstance(value, (int, float)):
                args.extend([f'--{key.replace("_", "-")}', str(value)])
            elif isinstance(value, str):
                args.extend([f'--{key.replace("_", "-")}', value])
            elif isinstance(value, list):
                for item in value:
                    args.extend([f'--{key.replace("_", "-")}', str(item)])

    return args


def build_amalgkit_command(subcommand: str, params: AmalgkitParams | None = None) -> List[str]:
    """Build full amalgkit command.

    Args:
        subcommand: Amalgkit subcommand
        params: Parameters for the command

    Returns:
        Complete command as list of strings
    """
    command = ['amalgkit', subcommand]

    if params:
        command.extend(build_cli_args(params))

    return command


def check_cli_available() -> Tuple[bool, str]:
    """Check if amalgkit CLI is available.

    Returns:
        Tuple of (available, message)
    """
    try:
        result = subprocess.run(
            ['amalgkit', '--help'],
            capture_output=True,
            text=True,
            timeout=10
        )
        if result.returncode == 0:
            return True, "amalgkit CLI is available"
        else:
            return False, f"amalgkit CLI returned error: {result.stderr}"
    except FileNotFoundError:
        return False, "amalgkit CLI not found in PATH"
    except subprocess.TimeoutExpired:
        return False, "amalgkit CLI check timed out"
    except Exception as e:
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
            result = subprocess.run(
                ['amalgkit', '--version'],
                capture_output=True,
                text=True,
                timeout=5
            )
            version_info = {'version': result.stdout.strip()} if result.returncode == 0 else None
        except:
            version_info = None

        return True, message, version_info

    if not auto_install:
        return False, message, None

    # Attempt installation (placeholder - would need actual installation logic)
    logger.info("Attempting to install amalgkit...")
    try:
        # This would be the actual installation command
        # For now, just return failure
        return False, "Automatic installation not implemented", None
    except Exception as e:
        return False, f"Installation failed: {e}", None


def run_amalgkit(subcommand: str, params: AmalgkitParams | None = None, **kwargs: Any) -> subprocess.CompletedProcess[str]:
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

    # Default kwargs
    run_kwargs = {
        'capture_output': True,
        'text': True,
        'check': False,  # Don't raise on non-zero exit
        **kwargs
    }

    try:
        result = subprocess.run(command, **run_kwargs)
        logger.debug(f"amalgkit {subcommand} completed with return code {result.returncode}")
        return result
    except Exception as e:
        logger.error(f"Error running amalgkit {subcommand}: {e}")
        raise


def metadata(params: AmalgkitParams | None = None, **kwargs: Any) -> subprocess.CompletedProcess[str]:
    """Run amalgkit metadata command.

    Args:
        params: Amalgkit parameters
        **kwargs: Additional arguments

    Returns:
        CompletedProcess instance
    """
    return run_amalgkit('metadata', params, **kwargs)


def integrate(params: AmalgkitParams | None = None, **kwargs: Any) -> subprocess.CompletedProcess[str]:
    """Run amalgkit integrate command.

    Args:
        params: Amalgkit parameters
        **kwargs: Additional arguments

    Returns:
        CompletedProcess instance
    """
    return run_amalgkit('integrate', params, **kwargs)


def config(params: AmalgkitParams | None = None, **kwargs: Any) -> subprocess.CompletedProcess[str]:
    """Run amalgkit config command.

    Args:
        params: Amalgkit parameters
        **kwargs: Additional arguments

    Returns:
        CompletedProcess instance
    """
    return run_amalgkit('config', params, **kwargs)


def select(params: AmalgkitParams | None = None, **kwargs: Any) -> subprocess.CompletedProcess[str]:
    """Run amalgkit select command.

    Args:
        params: Amalgkit parameters
        **kwargs: Additional arguments

    Returns:
        CompletedProcess instance
    """
    return run_amalgkit('select', params, **kwargs)


def getfastq(params: AmalgkitParams | None = None, **kwargs: Any) -> subprocess.CompletedProcess[str]:
    """Run amalgkit getfastq command.

    Args:
        params: Amalgkit parameters
        **kwargs: Additional arguments

    Returns:
        CompletedProcess instance
    """
    return run_amalgkit('getfastq', params, **kwargs)


def quant(params: AmalgkitParams | None = None, **kwargs: Any) -> subprocess.CompletedProcess[str]:
    """Run amalgkit quant command.

    Args:
        params: Amalgkit parameters
        **kwargs: Additional arguments

    Returns:
        CompletedProcess instance
    """
    return run_amalgkit('quant', params, **kwargs)


def merge(params: AmalgkitParams | None = None, **kwargs: Any) -> subprocess.CompletedProcess[str]:
    """Run amalgkit merge command.

    Args:
        params: Amalgkit parameters
        **kwargs: Additional arguments

    Returns:
        CompletedProcess instance
    """
    return run_amalgkit('merge', params, **kwargs)


def cstmm(params: AmalgkitParams | None = None, **kwargs: Any) -> subprocess.CompletedProcess[str]:
    """Run amalgkit cstmm command.

    Args:
        params: Amalgkit parameters
        **kwargs: Additional arguments

    Returns:
        CompletedProcess instance
    """
    return run_amalgkit('cstmm', params, **kwargs)


def curate(params: AmalgkitParams | None = None, **kwargs: Any) -> subprocess.CompletedProcess[str]:
    """Run amalgkit curate command.

    Args:
        params: Amalgkit parameters
        **kwargs: Additional arguments

    Returns:
        CompletedProcess instance
    """
    return run_amalgkit('curate', params, **kwargs)


def csca(params: AmalgkitParams | None = None, **kwargs: Any) -> subprocess.CompletedProcess[str]:
    """Run amalgkit csca command.

    Args:
        params: Amalgkit parameters
        **kwargs: Additional arguments

    Returns:
        CompletedProcess instance
    """
    return run_amalgkit('csca', params, **kwargs)


def sanity(params: AmalgkitParams | None = None, **kwargs: Any) -> subprocess.CompletedProcess[str]:
    """Run amalgkit sanity command.

    Args:
        params: Amalgkit parameters
        **kwargs: Additional arguments

    Returns:
        CompletedProcess instance
    """
    return run_amalgkit('sanity', params, **kwargs)



