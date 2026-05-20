"""Public registry for RNA workflow step runners."""

from __future__ import annotations

from metainformant.rna.amalgkit.amalgkit import config as run_config
from metainformant.rna.amalgkit.amalgkit import csca as run_csca
from metainformant.rna.amalgkit.amalgkit import cstmm as run_cstmm
from metainformant.rna.amalgkit.amalgkit import curate as run_curate
from metainformant.rna.amalgkit.amalgkit import getfastq as run_getfastq
from metainformant.rna.amalgkit.amalgkit import integrate as run_integrate
from metainformant.rna.amalgkit.amalgkit import merge as run_merge
from metainformant.rna.amalgkit.amalgkit import metadata as run_metadata
from metainformant.rna.amalgkit.amalgkit import quant as run_quant
from metainformant.rna.amalgkit.amalgkit import sanity as run_sanity
from metainformant.rna.amalgkit.amalgkit import select as run_select

STEP_RUNNERS = {
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
