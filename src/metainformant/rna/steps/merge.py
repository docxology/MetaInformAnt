from __future__ import annotations

"""Step runner for `amalgkit merge` (combine per-sample quantifications)."""

import subprocess
from collections.abc import Mapping
from pathlib import Path
from typing import Any

from ..amalgkit import merge as _merge


def run(
    params: Mapping[str, Any] | None = None,
    *,
    work_dir: str | Path | None = None,
    log_dir: str | Path | None = None,
    check: bool = False,
) -> subprocess.CompletedProcess[str]:
    """Run `amalgkit merge` (combine per-sample quantifications into expression matrix).
    
    Merges per-sample abundance files (abundance.tsv) into a single expression
    matrix with samples as columns and transcripts as rows.
    
    Args:
        params: Parameters for amalgkit merge step. Common parameters:
            - metadata: Path to metadata TSV file with sample list
            - quant_dir: Directory containing per-sample quantification results
            - out: Output file path for merged expression matrix (TSV format)
        work_dir: Working directory for amalgkit commands
        log_dir: Directory for step logs
        check: If True, raise CalledProcessError on non-zero exit
        
    Returns:
        subprocess.CompletedProcess with return code and output
        
    Raises:
        subprocess.CalledProcessError: If check=True and step fails
    """
    return _merge(params, work_dir=work_dir, log_dir=log_dir, step_name="merge", check=check)
