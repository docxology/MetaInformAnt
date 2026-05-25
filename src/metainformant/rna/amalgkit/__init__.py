"""Amalgkit RNA-seq workflow module exports."""

from __future__ import annotations

from . import amalgkit, genome_prep, index_prep, metadata_filter, metadata_utils, tissue_normalizer

# Expose core functions from amalgkit.py at package level for convenience
from .amalgkit import (
    AMALGKIT_INSTALL_SPEC,
    MIN_AMALGKIT_VERSION,
    AmalgkitParams,
    build_amalgkit_command,
    build_cli_args,
    check_cli_available,
    config,
    csca,
    cstmm,
    curate,
    ensure_cli_available,
    getfastq,
    integrate,
    merge,
    metadata,
    quant,
    run_amalgkit,
    sanity,
    select,
)

__all__ = [
    "amalgkit",
    "genome_prep",
    "index_prep",
    "metadata_filter",
    "metadata_utils",
    "tissue_normalizer",
    "AMALGKIT_INSTALL_SPEC",
    "AmalgkitParams",
    "MIN_AMALGKIT_VERSION",
    "build_amalgkit_command",
    "build_cli_args",
    "check_cli_available",
    "config",
    "csca",
    "cstmm",
    "curate",
    "ensure_cli_available",
    "getfastq",
    "integrate",
    "merge",
    "metadata",
    "quant",
    "run_amalgkit",
    "sanity",
    "select",
]
