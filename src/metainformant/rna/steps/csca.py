"""Step runner for `amalgkit csca` (correlation analysis and plots)."""

from __future__ import annotations

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
):
    """Run `amalgkit csca` (sample clustering/assessment)."""
    return _csca(params, work_dir=work_dir, log_dir=log_dir, step_name="csca", check=check)


