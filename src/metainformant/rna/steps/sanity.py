"""Step runner for `amalgkit sanity` (final integrity checks)."""

from __future__ import annotations

from collections.abc import Mapping
from pathlib import Path
from typing import Any

from ..amalgkit import sanity as _sanity


def run(
    params: Mapping[str, Any] | None = None,
    *,
    work_dir: str | Path | None = None,
    log_dir: str | Path | None = None,
    check: bool = False,
):
    """Run `amalgkit sanity` (final checks)."""
    return _sanity(params, work_dir=work_dir, log_dir=log_dir, step_name="sanity", check=check)


