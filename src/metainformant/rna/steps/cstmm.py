"""Step runner for `amalgkit cstmm` (cross-species TMM normalization)."""

from __future__ import annotations

import subprocess
from collections.abc import Mapping
from pathlib import Path
from typing import Any

from ..amalgkit import cstmm as _cstmm


def run(
    params: Mapping[str, Any] | None = None,
    *,
    work_dir: str | Path | None = None,
    log_dir: str | Path | None = None,
    check: bool = False,
) -> subprocess.CompletedProcess[str]:
    """Run `amalgkit cstmm` (cross-species TMM normalization).
    
    Performs Trimmed Mean of M-values (TMM) normalization across multiple
    species to enable cross-species expression comparisons.
    
    Args:
        params: Parameters for amalgkit cstmm step. Common parameters:
            - metadata: Path to metadata TSV file
            - merge_dir: Directory containing merged expression matrices
            - out_dir: Output directory for normalized matrices
        work_dir: Working directory for amalgkit commands
        log_dir: Directory for step logs
        check: If True, raise CalledProcessError on non-zero exit
        
    Returns:
        subprocess.CompletedProcess with return code and output
        
    Raises:
        subprocess.CalledProcessError: If check=True and step fails
        
    Note:
        This step requires merged expression matrices from multiple species.
    """
    return _cstmm(params, work_dir=work_dir, log_dir=log_dir, step_name="cstmm", check=check)
