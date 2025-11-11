"""Step runner for `amalgkit csca` (correlation analysis and plots)."""

from __future__ import annotations

import subprocess
from collections.abc import Mapping
from pathlib import Path
from typing import Any

from ..amalgkit import csca as _csca


def run(
    params: Mapping[str, Any] | None = None,
    *,
    work_dir: str | Path | None = None,
    log_dir: str | Path | None = None,
    check: bool = False,
) -> subprocess.CompletedProcess[str]:
    """Run `amalgkit csca` (cross-species correlation analysis and plots).
    
    Performs correlation analysis across multiple species to assess expression
    conservation and generate comparative plots.
    
    Args:
        params: Parameters for amalgkit csca step. Common parameters:
            - metadata: Path to metadata TSV file
            - curate_dir: Directory containing curated expression matrices
            - out_dir: Output directory for correlation results and plots
        work_dir: Working directory for amalgkit commands
        log_dir: Directory for step logs
        check: If True, raise CalledProcessError on non-zero exit
        
    Returns:
        subprocess.CompletedProcess with return code and output
        
    Raises:
        subprocess.CalledProcessError: If check=True and step fails
        
    Note:
        This step requires R (Rscript) for statistical analysis and plotting.
    """
    return _csca(params, work_dir=work_dir, log_dir=log_dir, step_name="csca", check=check)
