from __future__ import annotations

from typing import Mapping, Any
from pathlib import Path

from ..amalgkit import integrate as _integrate


def run(
    params: Mapping[str, Any] | None = None,
    *,
    work_dir: str | Path | None = None,
    log_dir: str | Path | None = None,
    check: bool = False,
):
    return _integrate(params, work_dir=work_dir, log_dir=log_dir, step_name="integrate", check=check)


