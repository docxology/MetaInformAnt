"""Step runner for `amalgkit integrate` (augment metadata with local FASTQs)."""

from __future__ import annotations

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
):
    """Run `amalgkit integrate`."""
    return _integrate(params, work_dir=work_dir, log_dir=log_dir, step_name="integrate", check=check)


