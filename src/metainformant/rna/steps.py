"""Public registry for RNA workflow step runners."""

from __future__ import annotations

from metainformant.rna.amalgkit.amalgkit import csfilter as run_csfilter
from metainformant.rna.amalgkit.amalgkit import cstmm as run_cstmm
from metainformant.rna.amalgkit.amalgkit import finalize as run_finalize
from metainformant.rna.amalgkit.amalgkit import getfastq as run_getfastq
from metainformant.rna.amalgkit.amalgkit import integrate as run_integrate
from metainformant.rna.amalgkit.amalgkit import merge as run_merge
from metainformant.rna.amalgkit.amalgkit import metadata as run_metadata
from metainformant.rna.amalgkit.amalgkit import quant as run_quant
from metainformant.rna.amalgkit.amalgkit import sanity as run_sanity
from metainformant.rna.amalgkit.amalgkit import select as run_select
from metainformant.rna.amalgkit.amalgkit import wsfilter as run_wsfilter

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
