"""Step runner for `amalgkit config`.

Provides a stable `run(...)` surface used by the RNA workflow orchestrator.
"""

from __future__ import annotations

from collections.abc import Mapping
from pathlib import Path
from typing import Any

from ..amalgkit import config as _config


def run(
    params: Mapping[str, Any] | None = None,
    *,
    work_dir: str | Path | None = None,
    log_dir: str | Path | None = None,
    check: bool = False,
):
    """Run `amalgkit config`.

    Parameters follow the common wrapper contract; see `metainformant.rna.amalgkit.run_amalgkit`.
    """
    return _config(params, work_dir=work_dir, log_dir=log_dir, step_name="config", check=check)


