"""Step runner for `amalgkit curate` (outlier removal and bias checks)."""

from __future__ import annotations

import subprocess
from collections.abc import Mapping
from pathlib import Path
from typing import Any

from ..amalgkit import curate as _curate


def run(
    params: Mapping[str, Any] | None = None,
    *,
    work_dir: str | Path | None = None,
    log_dir: str | Path | None = None,
    check: bool = False,
) -> subprocess.CompletedProcess[str]:
    """Run `amalgkit curate` (outlier removal and bias correction).
    
    Performs quality control, outlier detection, and batch effect correction
    on expression matrices. Generates QC plots and corrected expression data.
    
    Args:
        params: Parameters for amalgkit curate step. Common parameters:
            - metadata: Path to metadata TSV file
            - merge_dir: Directory containing merged expression matrix
            - out_dir: Output directory for curated matrices and QC plots
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
    return _curate(params, work_dir=work_dir, log_dir=log_dir, step_name="curate", check=check)
