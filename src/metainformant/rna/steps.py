"""Public registry for RNA workflow step runners."""

from __future__ import annotations

from metainformant.rna.amalgkit.amalgkit import (
    csfilter as run_csfilter,
    cstmm as run_cstmm,
    finalize as run_finalize,
    getfastq as run_getfastq,
    integrate as run_integrate,
    merge as run_merge,
    metadata as run_metadata,
    quant as run_quant,
    sanity as run_sanity,
    select as run_select,
    wsfilter as run_wsfilter,
)

STEP_RUNNERS = {
    "metadata": run_metadata,
    "integrate": run_integrate,
    "select": run_select,
    "getfastq": run_getfastq,
    "quant": run_quant,
    "merge": run_merge,
    "wsfilter": run_wsfilter,
    "cstmm": run_cstmm,
    "csfilter": run_csfilter,
    "finalize": run_finalize,
    "sanity": run_sanity,
}

__all__ = [
    "STEP_RUNNERS",
    "run_metadata",
    "run_integrate",
    "run_select",
    "run_getfastq",
    "run_quant",
    "run_merge",
    "run_wsfilter",
    "run_cstmm",
    "run_csfilter",
    "run_finalize",
    "run_sanity",
]
