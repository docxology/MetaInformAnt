"""Step runner for `amalgkit getfastq` (download raw FASTQ files)."""

from __future__ import annotations

from collections.abc import Mapping
from pathlib import Path
from typing import Any

from ..amalgkit import getfastq as _getfastq


def run(
    params: Mapping[str, Any] | None = None,
    *,
    work_dir: str | Path | None = None,
    log_dir: str | Path | None = None,
    check: bool = False,
):
    """Run `amalgkit getfastq` (downloading FASTQ files)."""
    return _getfastq(params, work_dir=work_dir, log_dir=log_dir, step_name="getfastq", check=check)


