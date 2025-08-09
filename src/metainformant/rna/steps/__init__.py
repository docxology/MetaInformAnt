from __future__ import annotations

from collections.abc import Callable

# Expose a uniform runner function per amalgkit step. Each module keeps
# any future step-specific logic isolated from the main orchestrator.

from .config import run as run_config
from .csca import run as run_csca
from .cstmm import run as run_cstmm
from .curate import run as run_curate
from .getfastq import run as run_getfastq
from .integrate import run as run_integrate
from .merge import run as run_merge
from .metadata import run as run_metadata
from .quant import run as run_quant
from .select import run as run_select
from .sanity import run as run_sanity


STEP_RUNNERS: dict[str, Callable[..., object]] = {
    "metadata": run_metadata,
    "integrate": run_integrate,
    "config": run_config,
    "select": run_select,
    "getfastq": run_getfastq,
    "quant": run_quant,
    "merge": run_merge,
    "cstmm": run_cstmm,
    "curate": run_curate,
    "csca": run_csca,
    "sanity": run_sanity,
}

__all__ = [
    "STEP_RUNNERS",
    "run_metadata",
    "run_integrate",
    "run_config",
    "run_select",
    "run_getfastq",
    "run_quant",
    "run_merge",
    "run_cstmm",
    "run_curate",
    "run_csca",
    "run_sanity",
]


