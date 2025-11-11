"""Step runner for `amalgkit config`.

Provides a stable `run(...)` surface used by the RNA workflow orchestrator.
"""

from __future__ import annotations

import subprocess
from collections.abc import Mapping
from pathlib import Path
from typing import Any

from ..amalgkit import config as _config


def run(
    params: Mapping[str, Any] | None = None,
    *,
    work_dir: str | Path | None = None,
    log_dir: str | Path | None = None,
    check: bool = False,
) -> subprocess.CompletedProcess[str]:
    """Run `amalgkit config` (generate configuration files for downstream tools).
    
    Generates configuration files (.config) for downstream analysis tools based
    on metadata. These configs define quality criteria for sample selection.
    
    Args:
        params: Parameters for amalgkit config step. Common parameters:
            - out_dir: Output directory for .config files
        work_dir: Working directory for amalgkit commands
        log_dir: Directory for step logs
        check: If True, raise CalledProcessError on non-zero exit
        
    Returns:
        subprocess.CompletedProcess with return code and output
        
    Raises:
        subprocess.CalledProcessError: If check=True and step fails
        
    Note:
        This step typically runs after metadata retrieval and before sample selection.
    """
    return _config(params, work_dir=work_dir, log_dir=log_dir, step_name="config", check=check)
