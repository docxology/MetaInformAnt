"""Step runner for `amalgkit metadata` (fetch SRA/ENA metadata)."""

from __future__ import annotations

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
):
    """Run `amalgkit metadata` (fetch SRA/ENA metadata)."""
    return _metadata(params, work_dir=work_dir, log_dir=log_dir, step_name="metadata", check=check)
