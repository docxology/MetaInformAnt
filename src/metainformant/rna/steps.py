"""Shim module for RNA workflow steps to maintain backward compatibility.

This module exposes the STEP_RUNNERS dictionary which maps step names to their
respective runner functions.
"""

from __future__ import annotations

from metainformant.rna.amalgkit.amalgkit import (
    config,
    csca,
    cstmm,
    curate,
    getfastq,
    integrate,
    merge,
    metadata,
    quant,
    sanity,
    select,
)

STEP_RUNNERS = {
    "metadata": metadata,
    "integrate": integrate,
    "config": config,
    "select": select,
    "getfastq": getfastq,
    "quant": quant,
    "merge": merge,
    "cstmm": cstmm,
    "curate": curate,
    "csca": csca,
    "sanity": sanity,
}
