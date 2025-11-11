"""Step runner for `amalgkit sanity` (final integrity checks)."""

from __future__ import annotations

import subprocess
from collections.abc import Mapping
from pathlib import Path
from typing import Any

from ..amalgkit import sanity as _sanity


def run(
    params: Mapping[str, Any] | None = None,
    *,
    work_dir: str | Path | None = None,
    log_dir: str | Path | None = None,
    check: bool = False,
) -> subprocess.CompletedProcess[str]:
    """Run `amalgkit sanity` (final integrity checks and validation).
    
    Performs final validation checks on workflow outputs to ensure data integrity
    and completeness. This is typically the last step in a workflow.
    
    Args:
        params: Parameters for amalgkit sanity step. Common parameters:
            - out_dir: Directory containing workflow outputs to validate
            - all: If True, check all output types (default: False)
        work_dir: Working directory for amalgkit commands
        log_dir: Directory for step logs
        check: If True, raise CalledProcessError on non-zero exit
        
    Returns:
        subprocess.CompletedProcess with return code and output
        
    Raises:
        subprocess.CalledProcessError: If check=True and step fails
    """
    return _sanity(params, work_dir=work_dir, log_dir=log_dir, step_name="sanity", check=check)
