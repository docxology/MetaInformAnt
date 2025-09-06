"""Step runner for `amalgkit select` (filter samples based on criteria)."""

from __future__ import annotations

from collections.abc import Mapping
from pathlib import Path
from typing import Any

from ..amalgkit import select as _select


def run(
    params: Mapping[str, Any] | None = None,
    *,
    work_dir: str | Path | None = None,
    log_dir: str | Path | None = None,
    check: bool = False,
):
    """Run `amalgkit select`."""
    return _select(params, work_dir=work_dir, log_dir=log_dir, step_name="select", check=check)
