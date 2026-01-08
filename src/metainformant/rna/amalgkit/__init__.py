"""Amalgkit RNA-seq workflow module exports."""

from .amalgkit import (
    AmalgkitParams,
    build_cli_args,
    build_amalgkit_command,
    check_cli_available,
    ensure_cli_available,
    run_amalgkit,
    metadata,
    integrate,
    config,
    select,
    getfastq,
    quant,
)
from . import genome_prep
from . import metadata_filter

__all__ = [
    "AmalgkitParams",
    "build_cli_args",
    "build_amalgkit_command",
    "check_cli_available",
    "ensure_cli_available",
    "run_amalgkit",
    "metadata",
    "integrate",
    "config",
    "select",
    "getfastq",
    "quant",
    "genome_prep",
    "metadata_filter",
]
