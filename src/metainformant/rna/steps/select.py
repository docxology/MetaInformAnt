"""Step runner for `amalgkit select` (filter samples based on criteria)."""

from __future__ import annotations

import subprocess
from collections.abc import Mapping
from pathlib import Path
from typing import Any

from ..amalgkit import select as _select


def run(
    params: Mapping[str, Any] | None = None,
    *,
    work_dir: str | Path | None = None,
    log_dir: str | Path | None = None,
    check: bool = False,
) -> subprocess.CompletedProcess[str]:
    """Run `amalgkit select` (filter samples based on quality criteria).
    
    Filters RNA-seq samples based on quality criteria defined in configuration
    files. This step produces a qualified sample list for downstream processing.
    
    Args:
        params: Parameters for amalgkit select step. Common parameters:
            - metadata: Path to input metadata TSV file
            - config_dir: Directory containing .config files with quality criteria
            - out_dir: Output directory for qualified sample list
        work_dir: Working directory for amalgkit commands
        log_dir: Directory for step logs
        check: If True, raise CalledProcessError on non-zero exit
        
    Returns:
        subprocess.CompletedProcess with return code and output
        
    Raises:
        subprocess.CalledProcessError: If check=True and step fails
    """
    return _select(params, work_dir=work_dir, log_dir=log_dir, step_name="select", check=check)
