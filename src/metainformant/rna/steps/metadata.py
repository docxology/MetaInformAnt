"""Step runner for `amalgkit metadata` (fetch SRA/ENA metadata)."""

from __future__ import annotations

import subprocess
from collections.abc import Mapping
from pathlib import Path
from typing import Any

from ..amalgkit import metadata as _metadata


def run(
    params: Mapping[str, Any] | None = None,
    *,
    work_dir: str | Path | None = None,
    log_dir: str | Path | None = None,
    check: bool = False,
) -> subprocess.CompletedProcess[str]:
    """Run `amalgkit metadata` (fetch SRA/ENA metadata).
    
    Retrieves RNA-seq sample metadata from NCBI SRA/ENA databases based on
    search criteria. This is typically the first step in an RNA-seq workflow.
    
    Args:
        params: Parameters for amalgkit metadata step. Common parameters:
            - search_string: NCBI Entrez search query (e.g., "species[Organism] AND RNA-Seq[Strategy]")
            - out_dir: Output directory for metadata files
            - entrez_email: Email for NCBI Entrez API (or set NCBI_EMAIL env var)
        work_dir: Working directory for amalgkit commands
        log_dir: Directory for step logs
        check: If True, raise CalledProcessError on non-zero exit
        
    Returns:
        subprocess.CompletedProcess with return code and output
        
    Raises:
        subprocess.CalledProcessError: If check=True and step fails
        
    Examples:
        >>> from metainformant.rna.steps import run_metadata
        >>> result = run_metadata({
        ...     "search_string": "Pogonomyrmex barbatus[Organism] AND RNA-Seq[Strategy]",
        ...     "out_dir": "output/amalgkit/pbarbatus/work"
        ... })
    """
    return _metadata(params, work_dir=work_dir, log_dir=log_dir, step_name="metadata", check=check)
