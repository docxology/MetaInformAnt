from __future__ import annotations

"""Step runner for `amalgkit merge` (combine per-sample quantifications)."""

from pathlib import Path
from typing import Any, Mapping

from ..amalgkit import merge as _merge


def run(
    params: Mapping[str, Any] | None = None,
    *,
    work_dir: str | Path | None = None,
    log_dir: str | Path | None = None,
    check: bool = False,
):
    """Run `amalgkit merge` (aggregating quantifications)."""
    return _merge(params, work_dir=work_dir, log_dir=log_dir, step_name="merge", check=check)
