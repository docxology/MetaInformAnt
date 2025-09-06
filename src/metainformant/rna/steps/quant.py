"""Step runner for `amalgkit quant` (quantify transcript abundances)."""

from __future__ import annotations

from collections.abc import Mapping
from pathlib import Path
from typing import Any

from ..amalgkit import quant as _quant


def run(
    params: Mapping[str, Any] | None = None,
    *,
    work_dir: str | Path | None = None,
    log_dir: str | Path | None = None,
    check: bool = False,
):
    """Run `amalgkit quant` (quantification, e.g., Salmon)."""
    return _quant(params, work_dir=work_dir, log_dir=log_dir, step_name="quant", check=check)
