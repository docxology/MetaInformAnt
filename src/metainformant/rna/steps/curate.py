"""Step runner for `amalgkit curate` (outlier removal and bias checks)."""

from __future__ import annotations

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
):
    """Run `amalgkit curate`."""
    return _curate(params, work_dir=work_dir, log_dir=log_dir, step_name="curate", check=check)
