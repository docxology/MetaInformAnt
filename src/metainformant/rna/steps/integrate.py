"""Step runner for `amalgkit integrate` (augment metadata with local FASTQs)."""

from __future__ import annotations

import subprocess
from collections.abc import Mapping
from pathlib import Path
from typing import Any

from ..amalgkit import integrate as _integrate


def run(
    params: Mapping[str, Any] | None = None,
    *,
    work_dir: str | Path | None = None,
    log_dir: str | Path | None = None,
    check: bool = False,
) -> subprocess.CompletedProcess[str]:
    """Run `amalgkit integrate` (augment metadata with local FASTQ paths).
    
    Integrates local FASTQ file paths into metadata tables. This step discovers
    FASTQ files in the filesystem and adds their paths to the metadata.
    
    Args:
        params: Parameters for amalgkit integrate step. Common parameters:
            - metadata: Path to input metadata TSV file
            - fastq_dir: Directory containing FASTQ files to discover
            - out_dir: Output directory for updated metadata
        work_dir: Working directory for amalgkit commands
        log_dir: Directory for step logs
        check: If True, raise CalledProcessError on non-zero exit
        
    Returns:
        subprocess.CompletedProcess with return code and output
        
    Raises:
        subprocess.CalledProcessError: If check=True and step fails
        
    Note:
        This step gracefully skips if no FASTQ files are found yet (expected
        before getfastq step completes).
    """
    return _integrate(params, work_dir=work_dir, log_dir=log_dir, step_name="integrate", check=check)
